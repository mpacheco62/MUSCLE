! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

program test_muscle_yield_yld2004
    implicit none
    logical :: passed

    call test_Yld2004_stresseq_1(passed)
    if (.not. passed) STOP 1

    call test_Yld2004_gradient_2(passed)
    if (.not. passed) STOP 2

    call test_Yld2004_gradient_numeric_3(passed)
    if (.not. passed) STOP 3

    print*, "Passed!", passed
    STOP 0
end program test_muscle_yield_yld2004

subroutine Yld2004_test_materials(yld)
    ! yld(1): isotropic (all 18 coefficients 1, so S' = S'' = s), a = 2: von Mises.
    ! yld(2): AA2090-T3, a = 8 (Barlat et al., 2005, Table 2). The table lists the shears as
    ! (yz, zx, xy); here (c44, c55, c66) = (xy, yz, xz) = (c66, c44, c55) of the table.
    use muscle_yield_yld2004
    implicit none
    type(Yld2004), intent(out) :: yld(2)

    call yld(1)%init(cp12=1D0, cp13=1D0, cp21=1D0, cp23=1D0, cp31=1D0, cp32=1D0, &
                     cp44=1D0, cp55=1D0, cp66=1D0, &
                     cpp12=1D0, cpp13=1D0, cpp21=1D0, cpp23=1D0, cpp31=1D0, cpp32=1D0, &
                     cpp44=1D0, cpp55=1D0, cpp66=1D0, a=2D0)
    call yld(2)%init(cp12=-0.069888D0, cp13=0.936408D0, cp21=0.079143D0, cp23=1.003060D0, &
                     cp31=0.524741D0, cp32=1.363180D0, cp44=0.954322D0, cp55=1.023770D0, cp66=1.069060D0, &
                     cpp12=0.981171D0, cpp13=0.476741D0, cpp21=0.575316D0, cpp23=0.866827D0, &
                     cpp31=1.145010D0, cpp32=-0.079294D0, cpp44=1.404620D0, cpp55=1.051660D0, &
                     cpp66=1.147100D0, a=8D0)
end subroutine Yld2004_test_materials

subroutine test_Yld2004_stresseq_1(passed)
    ! Equivalent stress against independent references. The AA2090-T3 values come from a
    ! Python script with mpmath at 90 digits that writes phi with traces of powers of S' and
    ! S'' (valid for even a), without eigenvalues.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_yld2004
    use muscle_yield_vonmises
    implicit none
    logical, intent(out) :: passed

    ! Measured errors are below 4e-16; 1e-12 leaves room for other compilers and still separates
    ! the shear assignments, which differ by more than 1e-3
    real(real64), parameter :: TOL = 1.0D-12
    type(Yld2004) :: yld(2)
    type(VonMises) :: vm
    type(ten_3D2Osym) :: s, s_aa(4)
    real(real64) :: expected(4), err
    integer :: i

    call Yld2004_test_materials(yld)

    ! Isotropic limit, general 3D state: 4 f**2 = sum_ij (s_i - s_j)**2, i.e. von Mises
    ! (three distinct eigenvalues with all shears present)
    call s%init((/180D0, -40D0, 25D0, 60D0, -30D0, 45D0/))
    err = abs(yld(1)%stress_eq(s) - vm%stress_eq(s))/vm%stress_eq(s)
    passed = err <= TOL
    if (.not. passed) print*, "Yld2004 isotropic vs von Mises, rel. error =", err
    if (.not. passed) return

    ! Uniaxial RD: S' and S'' diagonal, pairs lam'_i - lam''_j of both signs
    call s_aa(1)%init((/100D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    ! Uniaxial TD: the same with the other rows of C' and C''
    call s_aa(2)%init((/0D0, 100D0, 0D0, 0D0, 0D0, 0D0/))
    ! Uniaxial at 45 deg: the xy shear picks c44 (c66 in the paper)
    call s_aa(3)%init((/50D0, 50D0, 0D0, 50D0, 0D0, 0D0/))
    ! Rotated 3D state: the yz and xz shears pick c55 and c66
    call s_aa(4)%init((/120.5D0, -35.25D0, 60.75D0, 42.5D0, -28.0D0, 33.25D0/))
    expected = (/99.932153929799531134D0, 110.32696825195229388D0, &
                 122.30436942533559492D0, 187.28212180970985109D0/)
    do i = 1, 4
        err = abs(yld(2)%stress_eq(s_aa(i)) - expected(i))/expected(i)
        passed = err <= TOL
        if (.not. passed) print*, "Yld2004 AA2090-T3 state", i, "rel. error =", err
        if (.not. passed) return
    end do

    ! Homogeneity of degree 1 and pressure insensitivity on the rotated state
    s = 2.5D0*s_aa(4)
    err = abs(yld(2)%stress_eq(s) - 2.5D0*expected(4))/expected(4)
    passed = err <= TOL
    if (.not. passed) print*, "Yld2004 homogeneity, rel. error =", err
    if (.not. passed) return
    call s%init(s_aa(4)%vals + (/300D0, 300D0, 300D0, 0D0, 0D0, 0D0/))
    err = abs(yld(2)%stress_eq(s) - expected(4))/expected(4)
    passed = err <= TOL
    if (.not. passed) print*, "Yld2004 pressure insensitivity, rel. error =", err
end subroutine test_Yld2004_stresseq_1

subroutine test_Yld2004_gradient_2(passed)
    ! yld%dstressEq_dstress against analytical values: von Mises in the isotropic limit and an
    ! independent high-precision reference for AA2090-T3. Both states have normal and shear
    ! components at once.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_yld2004
    use muscle_yield_vonmises
    implicit none
    logical, intent(out) :: passed

    ! The base-class finite differences agree to about 1e-10 on these states (measured 9.8e-11
    ! and 3.9e-11); 1e-8 keeps a margin and still separates any wrong coefficient (> 1e-3)
    real(real64), parameter :: TOL = 1.0D-8
    type(Yld2004) :: yld(2)
    type(VonMises) :: vm
    type(ten_3D2Osym) :: s, g, expected

    call Yld2004_test_materials(yld)

    ! Isotropic limit, general 3D state: three distinct eigenvalues with all shears present
    call s%init((/180D0, -40D0, 25D0, 60D0, -30D0, 45D0/))
    g = yld(1)%dstressEq_dstress(s)
    passed = g%is_approx(vm%dstressEq_dstress(s), tol=TOL)
    if (.not. passed) print*, "Yld2004 gradient vs von Mises:", g%vals
    if (.not. passed) return

    ! AA2090-T3, rotated 3D state: reference from mpmath at 90 digits (phi from traces of powers
    ! of S' and S'', central differences, no eigenvalues)
    call s%init((/120.5D0, -35.25D0, 60.75D0, 42.5D0, -28.0D0, 33.25D0/))
    call expected%init((/0.46131648419506646702D0, -0.81171119189451749695D0, &
                         0.35039470769945102992D0, 0.31754139064342329668D0, &
                         -0.45754809859404116078D0, 0.43880414633419044754D0/))
    g = yld(2)%dstressEq_dstress(s)
    passed = g%is_approx(expected, tol=TOL)
    if (.not. passed) print*, "Yld2004 gradient vs reference:", g%vals
end subroutine test_Yld2004_gradient_2

subroutine Yld2004_derivative_materials(yld)
    ! Materials of the derivative tests 3 and 4: yld(1:2) of Yld2004_test_materials and
    ! yld(3): C' isotropic and C'' of AA2090-T3 with c''12 = c''13 = 1, a = 2: the first rows
    ! of C' and C'' match, so states in the material axes give a zero difference
    ! d_ij = lam'_i - lam''_j in an anisotropic material, where |d|**(a-2) = 0**0.
    ! yld(4): the coefficients of AA2090-T3 (as yld(2)) with a = 6, the exponent for BCC metals
    ! (Barlat et al., 2005, p. 1017).
    use muscle_yield_yld2004
    implicit none
    type(Yld2004), intent(out) :: yld(4)

    call Yld2004_test_materials(yld(1:2))
    call yld(3)%init(cp12=1D0, cp13=1D0, cp21=1D0, cp23=1D0, cp31=1D0, cp32=1D0, &
                     cp44=1D0, cp55=1D0, cp66=1D0, &
                     cpp12=1D0, cpp13=1D0, cpp21=0.575316D0, cpp23=0.866827D0, &
                     cpp31=1.145010D0, cpp32=-0.079294D0, cpp44=1.404620D0, cpp55=1.051660D0, &
                     cpp66=1.147100D0, a=2D0)
    call yld(4)%init(cp12=-0.069888D0, cp13=0.936408D0, cp21=0.079143D0, cp23=1.003060D0, &
                     cp31=0.524741D0, cp32=1.363180D0, cp44=0.954322D0, cp55=1.023770D0, cp66=1.069060D0, &
                     cpp12=0.981171D0, cpp13=0.476741D0, cpp21=0.575316D0, cpp23=0.866827D0, &
                     cpp31=1.145010D0, cpp32=-0.079294D0, cpp44=1.404620D0, cpp55=1.051660D0, &
                     cpp66=1.147100D0, a=6D0)
end subroutine Yld2004_derivative_materials

subroutine Yld2004_test_states(s)
    ! States shared by the derivative tests 3 and 4 (Voigt xx, yy, zz, xy, yz, xz).
    use muscle_tensors
    implicit none
    type(ten_3D2Osym), intent(out) :: s(5)

    ! Rotated 3D state: three distinct eigenvalues in S' and S'', all shears present
    call s(1)%init((/120.5D0, -35.25D0, 60.75D0, 42.5D0, -28.0D0, 33.25D0/))
    ! Uniaxial RD: repeated pair in S' = S'' (isotropic); in yld(3), repeated pair in S' and
    ! d_11 = 1e-14 (eigenvals round-off)
    call s(2)%init((/100D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    ! Equibiaxial: repeated pair in S' = S'' (isotropic, the other side of the deviatoric
    ! plane); in yld(3), repeated pair in S' and d_22 = 0 exactly
    call s(3)%init((/100D0, 100D0, 0D0, 0D0, 0D0, 0D0/))
    ! Pure shear xy: a zero eigenvalue in S' and S''
    call s(4)%init((/0D0, 0D0, 0D0, 50D0, 0D0, 0D0/))
    ! Uniaxial along (1, 2, 2)/3: normal and shear components at once with a repeated pair in
    ! S' (isotropic and yld(3)) up to round-off
    call s(5)%init((100D0/9D0)*(/1D0, 4D0, 4D0, 2D0, 4D0, 2D0/))
end subroutine Yld2004_test_states

subroutine test_Yld2004_gradient_numeric_3(passed)
    ! yld%dstressEq_dstress against the finite differences of the base class (four materials,
    ! five states), and against von Mises in the isotropic limit and the AA2090-T3 reference of
    ! test 2 with a tolerance at round-off level.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_yield_yld2004
    use muscle_yield_vonmises
    implicit none
    logical, intent(out) :: passed

    ! The base-class finite differences agree to 2e-10 or better on these states; von Mises and
    ! the reference hold to round-off (measured 7e-16): 1e-12 leaves room for other compilers
    real(real64), parameter :: TOL_NUM = 1.0D-8, TOL = 1.0D-12
    type(Yld2004) :: yld(4)
    type(VonMises) :: vm
    type(ten_3D2Osym) :: s(5), g, expected
    integer :: m, i

    call Yld2004_derivative_materials(yld)
    call Yld2004_test_states(s)
    do m = 1, 4
        do i = 1, 5
            g = yld(m)%dstressEq_dstress(s(i))
            passed = g%is_approx(yld(m)%dstressEq_dstress_numeric(s(i)), tol=TOL_NUM)
            if (.not. passed) print*, "Yld2004 gradient vs numeric, material", m, "state", i, g%vals
            if (.not. passed) return
            if (m == 1) then
                passed = g%is_approx(vm%dstressEq_dstress(s(i)), tol=TOL)
                if (.not. passed) print*, "Yld2004 gradient vs von Mises, state", i, g%vals
                if (.not. passed) return
            end if
        end do
    end do

    ! AA2090-T3, rotated 3D state: the reference of test 2 (mpmath at 90 digits)
    call expected%init((/0.46131648419506646702D0, -0.81171119189451749695D0, &
                         0.35039470769945102992D0, 0.31754139064342329668D0, &
                         -0.45754809859404116078D0, 0.43880414633419044754D0/))
    g = yld(2)%dstressEq_dstress(s(1))
    passed = g%is_approx(expected, tol=TOL)
    if (.not. passed) print*, "Yld2004 gradient vs reference:", g%vals
end subroutine test_Yld2004_gradient_numeric_3
