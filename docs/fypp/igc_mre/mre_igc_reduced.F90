! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
!
! Reduced reproducer: wrong results of an OpenMP target region on an Intel GPU (ifx -fopenmp-targets=spir64, IGC
! with FP64 emulation). Standalone, single precision only, no external library.
!
! Each work-item runs a return-mapping loop on one point. Lane 3 of every aligned block of 16 work-items is
! "plastic" (it iterates through yield_df: two cyclic Jacobi decompositions and SCALE); the other 15 are elastic
! (the loop exits at i = 1, iters = 0). Expected output: "bad points: 0 of 256, elastic 240".
! Observed (-O2, -O3): "bad points: 128": the 8 elastic lanes of the SIMD8 thread that shares the block with the
! plastic lane's thread (the other half of the block) return iters = 50 instead of 0; the stress is right.
!
! Needed: SCALE on real32 in the plastic branch (ifx lowers it to emulated FP64 calls; x*2.0**(-e) gives 0); the
! rotation indices read from PARAMETER arrays (literal indices give 0); -O2 or above (-O0 and, for this file, -O1
! give 0). IGC_EmulationFunctionControl=1 gives 0. The failure is sensitive to register allocation: replacing a
! run-time value by a constant (e.g. the elastic constants or the initial plastic strain) also makes it disappear.
!
! Build and run (ifx 2026.2.0, Intel UHD Graphics 770; IGC 2.28.4 / NEO 26.05 and IGC 2.42.6 / NEO 26.40):
!   ifx -O2 -fiopenmp -fopenmp-targets=spir64 -fpp mre_igc_reduced.F90 -o mre_red
!   OMP_TARGET_OFFLOAD=MANDATORY OverrideDefaultFP64Settings=1 IGC_EnableDPEmulation=1 ./mre_red   ! 128 bad
!   OMP_TARGET_OFFLOAD=DISABLED ./mre_red                                                          ! 0 bad
module mre
    use, intrinsic :: iso_fortran_env, only: real32, real64
    implicit none
    integer, parameter :: sp = real32
    type :: sym
        real(sp) :: vals(6)
    end type sym
    type :: yield_t
        real(sp) :: Lp(3,3)
        real(sp) :: cp_shear(3)
        real(sp) :: Lpp(3,3)
        real(sp) :: cpp_shear(3)
        real(sp) :: a
    end type yield_t
    type :: mat_t
        real(sp) :: lame(2) = [121153.84_sp, 80769.23_sp]
        type(yield_t) :: yield
    end type mat_t
contains
    pure function add(a, b) result(r)
        !$omp declare target
        type(sym), intent(in) :: a, b
        type(sym) :: r
        r%vals = a%vals + b%vals
    end function add
    pure function ddot(a, b) result(r)
        !$omp declare target
        type(sym), intent(in) :: a, b
        real(sp) :: r
        r = a%vals(1)*b%vals(1) + a%vals(2)*b%vals(2) + a%vals(3)*b%vals(3) &
            + 2*(a%vals(4)*b%vals(4) + a%vals(5)*b%vals(5) + a%vals(6)*b%vals(6))
    end function ddot
    pure function grad_eigen_sp(V, dphi) result(res)
        implicit none
!$omp declare target
        real(sp), intent(in) :: V(3,3)
        real(sp), intent(in) :: dphi(3)
        type(sym) :: res
        res%vals(1) = dphi(1)*V(1,1)*V(1,1) + dphi(2)*V(1,2)*V(1,2) + dphi(3)*V(1,3)*V(1,3)
        res%vals(2) = dphi(1)*V(2,1)*V(2,1) + dphi(2)*V(2,2)*V(2,2) + dphi(3)*V(2,3)*V(2,3)
        res%vals(3) = dphi(1)*V(3,1)*V(3,1) + dphi(2)*V(3,2)*V(3,2) + dphi(3)*V(3,3)*V(3,3)
        res%vals(4) = dphi(1)*V(1,1)*V(2,1) + dphi(2)*V(1,2)*V(2,2) + dphi(3)*V(1,3)*V(2,3)
        res%vals(5) = dphi(1)*V(2,1)*V(3,1) + dphi(2)*V(2,2)*V(3,2) + dphi(3)*V(2,3)*V(3,3)
        res%vals(6) = dphi(1)*V(1,1)*V(3,1) + dphi(2)*V(1,2)*V(3,2) + dphi(3)*V(1,3)*V(3,3)
    end function grad_eigen_sp
    pure subroutine jacobi_tangent(app, aqq, apq, t)
        implicit none
!$omp declare target
        real(sp), intent(in) :: app, aqq, apq
        real(sp), intent(out) :: t
        real(sp) :: diff, theta
        diff = aqq - app
        if (abs(2.0_sp*apq) <= 1.0e-15_sp*abs(diff)) then
            t = apq/diff
        else
            theta = diff/(2.0_sp*apq)
            t = sign(1.0_sp, theta)/(abs(theta) + sqrt(theta*theta + 1.0_sp))
        end if
    end subroutine jacobi_tangent
    pure subroutine jacobi_core(a, v, max_sweeps)
        implicit none
!$omp declare target
        real(sp), intent(inout) :: a(3, 3), v(3, 3)
        integer, intent(in) :: max_sweeps
        integer, parameter :: P_OF(3) = [1, 1, 2], Q_OF(3) = [2, 3, 3], R_OF(3) = [3, 2, 1]
        real(sp) :: apq, t, c, sn, tau, g, h
        integer :: sweep, i, j, p, q, r
        do sweep = 1, max_sweeps
            do i = 1, 3
                p = P_OF(i)
                q = Q_OF(i)
                r = R_OF(i)
                apq = a(p, q)
                call jacobi_tangent(a(p, p), a(q, q), apq, t)
                c = 1.0_sp/sqrt(t*t + 1.0_sp)
                sn = t*c
                tau = sn/(1.0_sp + c)
                a(p, p) = a(p, p) - t*apq
                a(q, q) = a(q, q) + t*apq
                a(p, q) = 0.0_sp
                a(q, p) = 0.0_sp
                g = a(r, p)
                h = a(r, q)
                a(r, p) = g - sn*(h + g*tau)
                a(p, r) = a(r, p)
                a(r, q) = h + sn*(g - h*tau)
                a(q, r) = a(r, q)
                do j = 1, 3
                    g = v(j, p)
                    h = v(j, q)
                    v(j, p) = g - sn*(h + g*tau)
                    v(j, q) = h + sn*(g - h*tau)
                end do
            end do
        end do
    end subroutine jacobi_core
    pure subroutine jacobi(s, lam, v)
        implicit none
!$omp declare target
        type(sym), intent(in) :: s
        real(sp), intent(out) :: lam(3), v(3, 3)
        real(sp) :: a(3, 3)
        integer :: i
        a(1, 1) = s%vals(1)
        a(2, 2) = s%vals(2)
        a(3, 3) = s%vals(3)
        a(1, 2) = s%vals(4)
        a(2, 1) = s%vals(4)
        a(2, 3) = s%vals(5)
        a(3, 2) = s%vals(5)
        a(1, 3) = s%vals(6)
        a(3, 1) = s%vals(6)
        v = 0.0_sp
        v(1, 1) = 1.0_sp
        v(2, 2) = 1.0_sp
        v(3, 3) = 1.0_sp
        call jacobi_core(a, v, 4)
        do i = 1, 3
            lam(i) = a(i, i)
        end do
    end subroutine jacobi
    pure function yield_f(prm, sigma) result(res)
        implicit none
!$omp declare target
        type(yield_t), intent(in) :: prm
        type(sym), intent(in) :: sigma
        real(sp) :: res
        real(sp) :: lp(3), lpp(3)
        type(sym) :: tp, tpp
        tp = linear_transform(prm%Lp, prm%cp_shear, sigma)
        tpp = linear_transform(prm%Lpp, prm%cpp_shear, sigma)
        lp = tp%vals(1:3)
        lpp = tpp%vals(1:3)
        res = maxval(abs(lp)) + maxval(abs(lpp))
    end function yield_f
    pure function yield_df(prm, sigma) result(res)
        implicit none
!$omp declare target
        type(yield_t), intent(in) :: prm
        type(sym), intent(in) :: sigma
        type(sym) :: res
        real(sp) :: lp(3), lpp(3), Vp(3,3), Vpp(3,3), dfp(3), dfpp(3), x(3), y(3)
        integer :: e
        call jacobi(linear_transform(prm%Lp, prm%cp_shear, sigma), lp, Vp)
        call jacobi(linear_transform(prm%Lpp, prm%cpp_shear, sigma), lpp, Vpp)
        e = 8
        x = scale(lp, -e)
        y = scale(lpp, -e)
        dfp = sign(1.0_sp, x)
        dfpp = -sign(1.0_sp, y)
        res = add(grad_eigen_sp(Vp, dfp), grad_eigen_sp(Vpp, dfpp))
    end function yield_df
    pure function linear_transform(L, c_shear, stress) result(res)
        implicit none
!$omp declare target
        real(sp), intent(in) :: L(3,3)
        real(sp), intent(in) :: c_shear(3)
        type(sym), intent(in) :: stress
        type(sym) :: res
        res%vals(1) = L(1,1)*stress%vals(1) + L(1,2)*stress%vals(2) + L(1,3)*stress%vals(3)
        res%vals(2) = L(2,1)*stress%vals(1) + L(2,2)*stress%vals(2) + L(2,3)*stress%vals(3)
        res%vals(3) = L(3,1)*stress%vals(1) + L(3,2)*stress%vals(2) + L(3,3)*stress%vals(3)
        res%vals(4) = c_shear(1)*stress%vals(4)
        res%vals(5) = c_shear(2)*stress%vals(5)
        res%vals(6) = c_shear(3)*stress%vals(6)
    end function linear_transform
    pure subroutine solve(mat, lame, strain, sp0, stress, iters)
        !$omp declare target
        type(mat_t), intent(in) :: mat
        real(sp), intent(in) :: lame(2)
        type(sym), intent(in) :: strain
        real(sp), intent(in) :: sp0(7)
        type(sym), intent(out) :: stress
        integer, intent(out) :: iters
        type(sym) :: strain_p, df, e
        real(sp) :: strain_pf, hard, dhard, f, ddgamma, tr
        integer :: i
        strain_p%vals = sp0(1:6)
        strain_pf = sp0(7)
        do i = 1, 50
            e%vals = strain%vals - strain_p%vals
            tr = e%vals(1) + e%vals(2) + e%vals(3)
            stress%vals(1:3) = lame(1)*tr + 2.0_sp*lame(2)*e%vals(1:3)
            stress%vals(4:6) = 2.0_sp*lame(2)*e%vals(4:6)
            hard = 200.0_sp + 500.0_sp*strain_pf
            f = yield_f(mat%yield, stress) - hard
            if (f < 1.0e-5_sp*hard) exit
            df = yield_df(mat%yield, stress)
            dhard = 500.0_sp
            tr = df%vals(1) + df%vals(2) + df%vals(3)
            ddgamma = f/(lame(1)*tr*tr + 2.0_sp*lame(2)*ddot(df, df) + dhard)
            strain_pf = strain_pf + ddgamma
            strain_p%vals = strain_p%vals + ddgamma*df%vals
        end do
        iters = i - 1
    end subroutine solve
    pure subroutine point(mat, p, npts, eps, hsv, sig, its)
        !$omp declare target
        type(mat_t), intent(in) :: mat
        integer, intent(in) :: p, npts
        real(sp), intent(in) :: eps(npts, 6), hsv(npts, 7)
        real(sp), intent(inout) :: sig(npts, 6)
        integer, intent(inout) :: its(npts)
        type(sym) :: e, s
        e%vals = eps(p, :)
        call solve(mat, mat%lame, e, hsv(p, :), s, its(p))
        sig(p, :) = s%vals
    end subroutine point
    subroutine batch(mat, npts, eps, hsv, sig, its)
        type(mat_t), intent(in) :: mat
        integer, intent(in) :: npts
        real(sp), intent(in) :: eps(npts, 6), hsv(npts, 7)
        real(sp), intent(inout) :: sig(npts, 6)
        integer, intent(inout) :: its(npts)
        integer :: i
        !$omp target teams distribute parallel do map(present, alloc: mat, eps, hsv, sig, its)
        do i = 1, npts
            call point(mat, i, npts, eps, hsv, sig, its)
        end do
    end subroutine batch
end module mre
program kern
    use mre
    implicit none
    integer, parameter :: N = 256
    real(sp), parameter :: A_Y = 200.0_sp/210000.0_sp
    type(mat_t) :: mat
    real(sp) :: hsv(N, 7), eps(N, 6), sg(N, 6), sc(N, 6), d1(6), es
    integer :: itg(N), itc(N), i, nbad
        mat%yield%Lp = reshape([2.8884000e-01_sp, 2.8159133e-01_sp, 1.0456600e-01_sp, &
            3.5872800e-01_sp, 3.6073433e-01_sp, -7.3387300e-01_sp, &
            -6.4756800e-01_sp, -6.4232567e-01_sp, 6.2930700e-01_sp], [3, 3])
        mat%yield%Lpp = reshape([4.8597067e-01_sp, -9.4601667e-02_sp, -7.8977133e-01_sp, &
            -4.9520033e-01_sp, 4.8071433e-01_sp, 4.3453267e-01_sp, &
            9.2296667e-03_sp, -3.8611267e-01_sp, 3.5523867e-01_sp], [3, 3])
        mat%yield%cp_shear = [0.954322_sp, 1.023770_sp, 1.069060_sp]
        mat%yield%cpp_shear = [1.404620_sp, 1.051660_sp, 1.147100_sp]
        mat%yield%a = 8.0_sp
    d1 = [1.0_sp, -0.3_sp, -0.3_sp, 0.2_sp, 0.05_sp, 0.1_sp]
    do i = 1, N
        ! one plastic work-item (lane 3) in every aligned block of 16, the others elastic
        if (mod(i - 1, 16) == 3) then
            eps(i, :) = 3.0_sp*A_Y*d1
        else
            eps(i, :) = 0.5_sp*A_Y*d1
        end if
    end do
    hsv = 0.0_sp
    sg = 0.0_sp
    itg = -1
    !$omp target data map(to: mat, eps, hsv) map(tofrom: sg, itg)
    call batch(mat, N, eps, hsv, sg, itg)
    !$omp end target data
    do i = 1, N
        call point(mat, i, N, eps, hsv, sc, itc)
    end do
    nbad = 0
    do i = 1, N
        es = maxval(abs(sg(i, :) - sc(i, :)))/maxval(abs(sc(i, :)))
        if (itg(i) /= itc(i) .or. .not. es <= 1.0e-3_sp) nbad = nbad + 1
    end do
    print '(a, i0, a, i0, a, i0)', 'bad points: ', nbad, ' of ', N, ', elastic ', count(itc == 0)
end program kern
