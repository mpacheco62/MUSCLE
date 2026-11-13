! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>
!
! Reproducer: wrong results of an OpenMP target region on an Intel GPU (ifx -fopenmp-targets=spir64, IGC with
! FP64 emulation). Standalone, single precision only, no external library.
!
! Each work-item runs an elastic-plastic return mapping on one point (a reduction of a material kernel). With this
! input 233 of the 256 points are elastic (the loop exits at i = 1, `iters` = 0) and 23 are plastic. The target
! region is compared with the same routine on the host. Expected output: "bad points: 0 of 256".
! Observed on the Intel GPU at -O1, -O2, -O3: 48 elastic points come back with the right stress but iters = 50: the
! 8 work-items of a SIMD8 thread with no plastic point whose neighbouring SIMD8 thread (same aligned block of 16)
! has one. -DPATTERN=1..7 (patterns file) chooses the plastic work-items of every block of 16.
!
! Ingredients (each one is needed):
! - the intrinsic SCALE on real32 values (yield_f, scaled): ifx compiles it with double precision operations, so
!   the device image needs FP64, emulated by IGC on this GPU (OverrideDefaultFP64Settings=1
!   IGC_EnableDPEmulation=1). Replacing SCALE(x, -e) by x*2.0**(-e) gives 0 bad points;
! - the cyclic Jacobi rotations of jacobi_core, which index the local arrays a and v with p, q, r read from
!   PARAMETER arrays (P_OF(i), ...). The same three rotations written with literal indices give 0 bad points;
! - divergence between two adjacent SIMD8 threads: one executes the plastic branch (yield_df, two Jacobi
!   decompositions), the other does not; the 6x6 tangent arguments of the kernel are needed as well;
! - -O1 or above (-O0 is correct). IGC_EmulationFunctionControl=1 (emulation functions inlined) gives 0 bad points.
!
! Build and run (ifx 2026.2.0, Intel UHD Graphics 770, IGC 2.28.4, intel-opencl-icd 26.05.37020):
!   ifx -O2 -fiopenmp -fopenmp-targets=spir64 -fpp mre_igc_scale_fp64emu.F90 -o mre
!   OMP_TARGET_OFFLOAD=MANDATORY OverrideDefaultFP64Settings=1 IGC_EnableDPEmulation=1 ./mre
!   OMP_TARGET_OFFLOAD=DISABLED ./mre      ! reference: bad points: 0
module mre
    use, intrinsic :: iso_fortran_env, only: real32, real64
    implicit none
    integer, parameter :: sp = real32
    type :: sym
        real(sp) :: vals(6)
    end type sym
    type :: iden_t
        integer :: dummy = 0
    end type iden_t
    type :: iden_s
        real(sp) :: val
    end type iden_s
    interface operator(+)
        module procedure add
    end interface
    interface operator(-)
        module procedure sub, sub_iden
    end interface
    interface operator(*)
        module procedure mul, mul_iden
    end interface
    interface operator(.ddot.)
        module procedure ddot
    end interface
    integer, parameter :: CAND_FIXED_SWEEPS = 4
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
        real(sp) :: tol = 1.0e-5_sp
        integer :: iter_nw = 50
    end type mat_t
contains
    pure function square(a) result(r)
        !$omp declare target
        type(sym), intent(in) :: a
        type(sym) :: r
        r%vals(1) = a%vals(1)*a%vals(1) + a%vals(4)*a%vals(4) + a%vals(6)*a%vals(6)
        r%vals(2) = a%vals(4)*a%vals(4) + a%vals(2)*a%vals(2) + a%vals(5)*a%vals(5)
        r%vals(3) = a%vals(6)*a%vals(6) + a%vals(5)*a%vals(5) + a%vals(3)*a%vals(3)
        r%vals(4) = a%vals(1)*a%vals(4) + a%vals(4)*a%vals(2) + a%vals(6)*a%vals(5)
        r%vals(5) = a%vals(4)*a%vals(6) + a%vals(2)*a%vals(5) + a%vals(5)*a%vals(3)
        r%vals(6) = a%vals(1)*a%vals(6) + a%vals(4)*a%vals(5) + a%vals(6)*a%vals(3)
    end function square
    pure function to_matrix(a) result(m)
        !$omp declare target
        type(sym), intent(in) :: a
        real(sp) :: m(3, 3)
        m(1, 1) = a%vals(1)
        m(2, 2) = a%vals(2)
        m(3, 3) = a%vals(3)
        m(1, 2) = a%vals(4)
        m(2, 1) = a%vals(4)
        m(2, 3) = a%vals(5)
        m(3, 2) = a%vals(5)
        m(1, 3) = a%vals(6)
        m(3, 1) = a%vals(6)
    end function to_matrix
    pure function add(a, b) result(r)
        !$omp declare target
        type(sym), intent(in) :: a, b
        type(sym) :: r
        r%vals = a%vals + b%vals
    end function add
    pure function sub(a, b) result(r)
        !$omp declare target
        type(sym), intent(in) :: a, b
        type(sym) :: r
        r%vals = a%vals - b%vals
    end function sub
    pure function sub_iden(a, b) result(r)
        !$omp declare target
        type(sym), intent(in) :: a
        type(iden_s), intent(in) :: b
        type(sym) :: r
        r%vals = a%vals
        r%vals(1:3) = a%vals(1:3) - b%val
    end function sub_iden
    pure function mul(c, a) result(r)
        !$omp declare target
        real(sp), intent(in) :: c
        type(sym), intent(in) :: a
        type(sym) :: r
        r%vals = c*a%vals
    end function mul
    pure function mul_iden(c, a) result(r)
        !$omp declare target
        real(sp), intent(in) :: c
        type(iden_t), intent(in) :: a
        type(iden_s) :: r
        r%val = c
    end function mul_iden
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
    pure function eigenvals_sp(mat) result(res)
        implicit none
!$omp declare target
        type(sym), intent(in) :: mat
        real(sp) :: res(3)
        real(sp), parameter :: EPS = 1.0e-9_sp
        real(sp) :: I1, J2, s, d, alpha, n1, n2
        real(sp) :: J2_sqrt
        real(sp) :: cd, lam_a, lam_b, lam_c, sd
        integer :: sj
        type(sym) :: Sm, T, d1, d2
        type(iden_t) :: Iden
        I1 = sum(mat%vals(1:3))
        J2 = ((mat%vals(1)-mat%vals(2))**2.0_sp + (mat%vals(2)-mat%vals(3))**2.0_sp &
             + (mat%vals(3)-mat%vals(1))**2.0_sp)/6.0_sp &
             + (mat%vals(4)**2 + mat%vals(5)**2 + mat%vals(6)**2)
        s = (J2/3.0_sp)**0.5_sp
        if (abs(s) <= EPS) then
            res = I1/3.0_sp
            return
        end if
        Sm = mat - (I1/3.0_sp)*Iden
        T = square(Sm) - (2.0_sp/3.0_sp)*J2*Iden
        d1 = T - s*Sm
        d2 = T + s*Sm
        n1 = d1 .ddot. d1
        n2 = d2 .ddot. d2
        J2_sqrt = J2**0.5_sp
        if (n2 <= n1/huge(n1)) then
            sj = -1
            alpha = 0.0_sp
        else
            d = (n1 / n2)**0.5_sp
            if (abs(1.0_sp - d) <= EPS) then
                res(1) = J2_sqrt
                res(2) = 0.0_sp
                res(3) = -res(1)
                res = res + I1/3.0_sp
                return
            end if
            sj = int(sign(1.0_sp,(1.0_sp-d)))
            alpha = 2.0_sp/3.0_sp * atan(d**sj)
        end if
        cd = sj*s*cos(alpha)
        lam_a = 2.0_sp*cd
        sd = J2_sqrt*sin(alpha)
        lam_b = -cd + sd
        lam_c = -cd - sd
        res(1) = max(lam_a, max(lam_b, lam_c))
        res(3) = min(lam_a, min(lam_b, lam_c))
        res(2) = lam_a + lam_b + lam_c - res(1) - res(3)
        res = res + I1/3.0_sp
    end function eigenvals_sp
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
        a = to_matrix(s)
        v = 0.0_sp
        v(1, 1) = 1.0_sp
        v(2, 2) = 1.0_sp
        v(3, 3) = 1.0_sp
        call jacobi_core(a, v, CAND_FIXED_SWEEPS)
        do i = 1, 3
            lam(i) = a(i, i)
        end do
    end subroutine jacobi
    pure subroutine scaled(lp, lpp, x, y, e)
        implicit none
!$omp declare target
        real(sp), intent(in) :: lp(3)
        real(sp), intent(in) :: lpp(3)
        real(sp), intent(out) :: x(3)
        real(sp), intent(out) :: y(3)
        integer, intent(out) :: e
        e = 8
        x = scale(lp, -e)
        y = scale(lpp, -e)
    end subroutine scaled
    pure subroutine principal(prm, x, y, e, phi, dfp, dfpp)
        implicit none
!$omp declare target
        type(yield_t), intent(in) :: prm
        real(sp), intent(in) :: x(3), y(3)
        integer, intent(in) :: e
        real(sp), intent(out) :: phi, dfp(3), dfpp(3)
        real(sp) :: d(3,3), dh(3,3), c
        integer :: i, j
        do j = 1, 3
            do i = 1, 3
                d(i,j) = x(i) - y(j)
            end do
        end do
        phi = sum(abs(d)**prm%a)
        dfp = 0.0_sp
        dfpp = 0.0_sp
        if (phi == 0.0_sp) return
        c = (0.25_sp*phi)**(1.0_sp/prm%a)/(prm%a*phi)
        dh = prm%a*abs(d)**(prm%a - 2.0_sp)*d
        do i = 1, 3
            dfp(i) = c*sum(dh(i,:))
            dfpp(i) = -c*sum(dh(:,i))
        end do
    end subroutine principal
    pure function yield_f(prm, sigma) result(res)
        implicit none
!$omp declare target
        type(yield_t), intent(in) :: prm
        type(sym), intent(in) :: sigma
        real(sp) :: res
        real(sp) :: lp(3), lpp(3), x(3), y(3)
        integer :: e
        lp = eigenvals_sp(linear_transform(prm%Lp, prm%cp_shear, sigma))
        lpp = eigenvals_sp(linear_transform(prm%Lpp, prm%cpp_shear, sigma))
        call scaled(lp, lpp, x, y, e)
        res = (0.25_sp*(abs(x(1) - y(1))**prm%a + abs(x(1) - y(2))**prm%a &
                         + abs(x(1) - y(3))**prm%a                         &
                         + abs(x(2) - y(1))**prm%a + abs(x(2) - y(2))**prm%a &
                         + abs(x(2) - y(3))**prm%a                         &
                         + abs(x(3) - y(1))**prm%a + abs(x(3) - y(2))**prm%a &
                         + abs(x(3) - y(3))**prm%a))**(1.0_sp/prm%a)
        res = scale(res, e)
    end function yield_f
    pure function yield_df(prm, sigma) result(res)
        implicit none
!$omp declare target
        type(yield_t), intent(in) :: prm
        type(sym), intent(in) :: sigma
        type(sym) :: res
        real(sp) :: lp(3), lpp(3), Vp(3,3), Vpp(3,3), dfp(3), dfpp(3), phi, x(3), y(3)
        integer :: e
        call jacobi(linear_transform(prm%Lp, prm%cp_shear, sigma), lp, Vp)
        call jacobi(linear_transform(prm%Lpp, prm%cpp_shear, sigma), lpp, Vpp)
        call scaled(lp, lpp, x, y, e)
        call principal(prm, x, y, e, phi, dfp, dfpp)
        if (phi == 0.0_sp) then
            res%vals = 0.0_sp
            return
        end if
        res = yield_gradient(prm, Vp, Vpp, dfp, dfpp)
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
    pure function yield_gradient(prm, Vp, Vpp, dfp, dfpp) result(res)
        implicit none
!$omp declare target
        type(yield_t), intent(in) :: prm
        real(sp), intent(in) :: Vp(3,3), Vpp(3,3)
        real(sp), intent(in) :: dfp(3), dfpp(3)
        type(sym) :: res
        res = grad_eigen_sp(Vp, dfp) + grad_eigen_sp(Vpp, dfpp)
    end function yield_gradient
    pure subroutine solve(mat, lame, strain, sp0, stress, tangent, iters)
        !$omp declare target
        type(mat_t), intent(in) :: mat
        real(sp), intent(in) :: lame(2)
        type(sym), intent(in) :: strain
        real(sp), intent(in) :: sp0(7)
        type(sym), intent(out) :: stress
        real(sp), intent(out) :: tangent(6, 6)
        integer, intent(out) :: iters
        type(sym) :: strain_p, df, e
        real(sp) :: strain_pf, hard, dhard, f, ddgamma, tr
        integer :: i, j
        strain_p%vals = sp0(1:6)
        strain_pf = sp0(7)
        do i = 1, 50
            e = strain - strain_p
            tr = e%vals(1) + e%vals(2) + e%vals(3)
            stress%vals(1:3) = lame(1)*tr + 2.0_sp*lame(2)*e%vals(1:3)
            stress%vals(4:6) = 2.0_sp*lame(2)*e%vals(4:6)
            hard = 200.0_sp + 500.0_sp*strain_pf
            f = yield_f(mat%yield, stress) - hard
            if (f < 1.0e-5_sp*hard) exit
            df = yield_df(mat%yield, stress)
            dhard = 500.0_sp
            tr = df%vals(1) + df%vals(2) + df%vals(3)
            ddgamma = f/(lame(1)*tr*tr + 2.0_sp*lame(2)*(df .ddot. df) + dhard)
            strain_pf = strain_pf + ddgamma
            strain_p = strain_p + ddgamma*df
        end do
        iters = i - 1
        tangent = 0.0_sp
    end subroutine solve
    pure subroutine point(mat, p, npts, eps, hsv, sig, tan, its)
        !$omp declare target
        type(mat_t), intent(in) :: mat
        integer, intent(in) :: p, npts
        real(sp), intent(in) :: eps(npts, 6), hsv(npts, 7)
        real(sp), intent(inout) :: sig(npts, 6), tan(npts, 6, 6)
        integer, intent(inout) :: its(npts)
        type(sym) :: e, s
        real(sp) :: t(6, 6)
        e%vals = eps(p, :)
        call solve(mat, mat%lame, e, hsv(p, :), s, t, its(p))
        sig(p, :) = s%vals
    end subroutine point
    subroutine batch(mat, npts, eps, hsv, sig, tan, its)
        type(mat_t), intent(in) :: mat
        integer, intent(in) :: npts
        real(sp), intent(in) :: eps(npts, 6), hsv(npts, 7)
        real(sp), intent(inout) :: sig(npts, 6), tan(npts, 6, 6)
        integer, intent(inout) :: its(npts)
        integer :: i
        !$omp target teams distribute parallel do map(present, alloc: mat, eps, hsv, sig, tan, its)
        do i = 1, npts
            call point(mat, i, npts, eps, hsv, sig, tan, its)
        end do
    end subroutine batch
end module mre
program kern
    use mre
    implicit none
    integer, parameter :: N = 256, I0 = 21000
    real(sp), parameter :: A_Y = 200.0_sp/210000.0_sp
    type(mat_t) :: mat
    real(sp) :: hsv(N, 7), eps(N, 6), sg(N, 6), tg(N, 6, 6), sc(N, 6), tc(N, 6, 6), d1(6), d2(6), f, es, et
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
    d2 = [0.1_sp, 0.2_sp, -0.3_sp, 0.15_sp, -0.05_sp, 0.1_sp]
    do i = 1, N
        f = 4.0_sp*(real(i + I0, sp) - 0.5_sp)/100000.0_sp
        eps(i, :) = f*A_Y*d1 + 0.3_sp*A_Y*sin(real(i + I0, sp))*d2
    end do
    hsv = 0.0_sp
    sg = 0.0_sp
    tg = 0.0_sp
    itg = -1
    !$omp target data map(to: mat, eps, hsv) map(tofrom: sg, tg, itg)
    call batch(mat, N, eps, hsv, sg, tg, itg)
    !$omp end target data
    do i = 1, N
        call point(mat, i, N, eps, hsv, sc, tc, itc)
    end do
    nbad = 0
    do i = 1, N
        es = maxval(abs(sg(i, :) - sc(i, :)))/maxval(abs(sc(i, :)))
        et = 0.0_sp
        if (itg(i) /= itc(i) .or. .not. (es <= 1.0e-3_sp .and. et <= 1.0e-3_sp)) nbad = nbad + 1
    end do
    print '(a, i0, a, i0, a, i0)', 'bad points: ', nbad, ' of ', N, ', elastic ', count(itc == 0)
end program kern
