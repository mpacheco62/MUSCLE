! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module test_2d_vs_3d_helpers
    !! Embedding/projection helpers between the 2D and 3D tensor types.
    !!
    !! A 2D tensor is exactly a 3D tensor whose out-of-plane coupling components
    !! (13, 23, 31, 32) vanish. Every 2D operator must therefore give the same result as
    !! the (independently implemented and tested) 3D operator applied to the embedded
    !! operands. This is the reference solution used by `test_operator_2d_vs_3d`.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none

    real(real64), parameter :: TOL = 1D-12
        !! Relative tolerance: the 2D and 3D kernels evaluate the same sums of a few
        !! products (or one LU inversion of a well-conditioned matrix), so they agree
        !! to a few ulp (~1e-15); 1e-12 leaves a wide margin.
    logical :: ok = .true.
        !! Accumulated result of all checks (see `muscle_test_utils`).

contains

    function e2s(a) result(r)
        !! Embeds a `ten_2D2Osym` into a `ten_3D2Osym` (yz = xz = 0).
        type(ten_2D2Osym), intent(in) :: a
            !! 2D symmetric tensor.
        type(ten_3D2Osym) :: r
            !! Embedded 3D tensor.
        call r%init(xx=a%vals(1), yy=a%vals(2), zz=a%vals(3), xy=a%vals(4), yz=0D0, xz=0D0)
    end function e2s

    function e2g(a) result(r)
        !! Embeds a `ten_2D2O` into a `ten_3D2O`.
        type(ten_2D2O), intent(in) :: a
            !! 2D general tensor.
        type(ten_3D2O) :: r
            !! Embedded 3D tensor.
        call r%init(xx=a%vals(1), xy=a%vals(3), xz=0D0, &
                    yx=a%vals(2), yy=a%vals(4), yz=0D0, &
                    zx=0D0,       zy=0D0,       zz=a%vals(5))
    end function e2g

    function e42(c, fill) result(r)
        !! Embeds a `ten_2D4O2sym` into a `ten_3D4O2sym`. With `fill`, the out-of-plane
        !! shear diagonal (yzyz, xzxz) is set to 1 so that the embedded tensor is invertible
        !! (block-diagonal, the 2D block is unaffected).
        type(ten_2D4O2sym), intent(in) :: c
            !! 2D minor-symmetric tensor.
        logical, intent(in) :: fill
            !! Fill the out-of-plane shear diagonal with 1.
        type(ten_3D4O2sym) :: r
            !! Embedded 3D tensor.
        real(real64) :: m(6,6)
        m = 0D0
        m(1:4,1:4) = c%vals
        if (fill) then
            m(5,5) = 1D0
            m(6,6) = 1D0
        end if
        call r%init(m)
    end function e42

    function e43(c, fill) result(r)
        !! Embeds a `ten_2D4O3sym` into a `ten_3D4O3sym` (see `e42`).
        type(ten_2D4O3sym), intent(in) :: c
            !! 2D fully symmetric tensor.
        logical, intent(in) :: fill
            !! Fill the out-of-plane shear diagonal with 1.
        type(ten_3D4O3sym) :: r
            !! Embedded 3D tensor.
        type(ten_2D4O2sym) :: full
        type(ten_3D4O2sym) :: r2
        full = c
        r2 = e42(full, fill)
        r = r2%convert_3sym()
    end function e43

    function p42(c) result(r)
        !! Projects a `ten_3D4O2sym` onto its in-plane 4x4 Voigt block.
        type(ten_3D4O2sym), intent(in) :: c
            !! 3D minor-symmetric tensor.
        type(ten_2D4O2sym) :: r
            !! 2D block.
        call r%init(c%vals(1:4,1:4))
    end function p42

    function p43(c) result(r)
        !! Projects a `ten_3D4O3sym` onto its in-plane block (fully symmetric).
        type(ten_3D4O3sym), intent(in) :: c
            !! 3D fully symmetric tensor.
        type(ten_2D4O3sym) :: r
            !! 2D block.
        type(ten_3D4O2sym) :: full
        type(ten_2D4O2sym) :: blk
        full = c
        blk = p42(full)
        r = blk%convert_3sym()
    end function p43

end module test_2d_vs_3d_helpers


program test_operator_2d_vs_3d
    !! Checks every 2D operator against the 3D operator applied to the embedded operands.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    use test_2d_vs_3d_helpers
    implicit none

    type(ten_2D2Osym) :: s1, s2
    type(ten_2D2O) :: g1, g2
    type(ten_2D4O3sym) :: C1, C2
    type(ten_2D4O2sym) :: D1
    type(iden_2O) :: I2
    type(iden_2OS) :: I2S
    type(iden_4O3T) :: I3
    type(iden_4O3TS) :: I3S
    type(iden_4O4T) :: I4
    type(iden_4O4TS) :: I4S

    type(ten_3D2Osym) :: S3a, S3b
    type(ten_3D2O) :: G3a, G3b
    type(ten_3D4O3sym) :: C3a, C3b
    type(ten_3D4O2sym) :: D3

    ! Generic (non-special) operands
    call s1%init(xx=1.5D0, yy=-2D0, zz=0.75D0, xy=0.6D0)
    call s2%init(xx=-0.4D0, yy=3D0, zz=1.25D0, xy=-1.1D0)
    call g1%init(xx=1.2D0, xy=0.3D0, yx=-0.7D0, yy=0.9D0, zz=1.1D0)
    call g2%init(xx=-0.5D0, xy=2D0, yx=1.5D0, yy=0.25D0, zz=-3D0)
    call C1%init((/12D0, 9D0, 7D0, 3D0, 2D0, 1.5D0, -0.5D0, 1D0, 0.25D0, 0.4D0/))
    call C2%init((/-1D0, 0.5D0, 2D0, 1D0, 3D0, -2D0, 0.1D0, 1D0, 4D0, -3D0/))
    call D1%init(reshape((/ 5D0, -1D0, 0.5D0, 1D0, &
                            1D0,  6D0, 2D0, -0.5D0, &
                            2D0,  1D0, 7D0, 0.25D0, &
                            0.5D0, 1D0, -1D0, 3D0 /), (/4,4/)))
    call I2S%init(1.7D0)
    I3S%val = -0.8D0
    I4S%val = 2.3D0

    S3a = e2s(s1); S3b = e2s(s2)
    G3a = e2g(g1); G3b = e2g(g2)
    C3a = e43(C1, .false.); C3b = e43(C2, .false.)
    D3 = e42(D1, .false.)

    call second_order()
    call fourth_order()

    if (.not. ok) then
        print *, "Some 2D vs 3D operator checks failed (see FAIL messages above)."
        flush(output_unit)
        stop 1
    end if
    print *, "All 2D vs 3D operator checks passed."
    stop 0

contains

    subroutine second_order()
        !! Second-order operators, methods and assignments.
        type(ten_2D2Osym) :: rs
        type(ten_2D2O) :: rg
        type(ten_3D2Osym) :: r3s
        type(ten_3D2O) :: r3g

        ! --- methods / unary operators ---
        call check_approx(ok, "dev 2D2Osym", e2s(.dev. s1), (.dev. S3a))
        call check_approx(ok, "dev 2D2O", e2g(.dev. g1), (.dev. G3a))
        call check_approx(ok, "square 2D2Osym", e2s(s1%square()), S3a%square())
        call check_approx(ok, "det 2D2Osym", s1%det(), S3a%det(), tol=TOL)
        call check_approx(ok, "det 2D2O", g1%det(), G3a%det(), tol=TOL)
        call check_approx(ok, "inv 2D2Osym", e2s(.inv. s1), (.inv. S3a))
        call check_approx(ok, "transpose 2D2O", e2g(transpose(g1)), transpose(G3a))
        call check_approx(ok, "norm 2D2Osym", s1%norm(), S3a%norm(), tol=TOL)
        call check_approx(ok, "norm 2D2O", g1%norm(), G3a%norm(), tol=TOL)

        ! --- double contraction (scalar) ---
        call check_approx(ok, "2D2Osym : 2D2Osym", s1 .ddot. s2, S3a .ddot. S3b, tol=TOL)
        call check_approx(ok, "2D2O : 2D2O", g1 .ddot. g2, G3a .ddot. G3b, tol=TOL)
        call check_approx(ok, "2D2O : 2D2Osym", g1 .ddot. s1, G3a .ddot. S3a, tol=TOL)
        call check_approx(ok, "2D2Osym : 2D2O", s1 .ddot. g1, S3a .ddot. G3a, tol=TOL)
        call check_approx(ok, "I : 2D2O", I2 .ddot. g1, I2 .ddot. G3a, tol=TOL)
        call check_approx(ok, "2D2O : cI", g1 .ddot. I2S, G3a .ddot. I2S, tol=TOL)

        ! --- 4th-order identities : 2nd order ---
        rs = I4 .ddot. g1;   r3s = I4 .ddot. G3a;  call check_approx(ok, "I^S : 2D2O", e2s(rs), r3s)
        rs = g1 .ddot. I4;   r3s = G3a .ddot. I4;  call check_approx(ok, "2D2O : I^S", e2s(rs), r3s)
        rs = I4S .ddot. g1;  r3s = I4S .ddot. G3a; call check_approx(ok, "cI^S : 2D2O", e2s(rs), r3s)
        rs = g1 .ddot. I4S;  r3s = G3a .ddot. I4S; call check_approx(ok, "2D2O : cI^S", e2s(rs), r3s)
        rs = I4S .ddot. s1;  r3s = I4S .ddot. S3a; call check_approx(ok, "cI^S : 2D2Osym", e2s(rs), r3s)
        rs = s1 .ddot. I4S;  r3s = S3a .ddot. I4S; call check_approx(ok, "2D2Osym : cI^S", e2s(rs), r3s)

        ! --- single contraction ---
        rg = s1 * s2;  r3g = S3a * S3b;  call check_approx(ok, "2D2Osym * 2D2Osym", e2g(rg), r3g)
        rg = g1 * s1;  r3g = G3a * S3a;  call check_approx(ok, "2D2O * 2D2Osym", e2g(rg), r3g)
        rg = s1 * g1;  r3g = S3a * G3a;  call check_approx(ok, "2D2Osym * 2D2O", e2g(rg), r3g)
        rg = g1 * g2;  r3g = G3a * G3b;  call check_approx(ok, "2D2O * 2D2O", e2g(rg), r3g)
        rg = I2S * g1; r3g = I2S * G3a;  call check_approx(ok, "cI * 2D2O", e2g(rg), r3g)
        rg = g1 * I2S; r3g = G3a * I2S;  call check_approx(ok, "2D2O * cI", e2g(rg), r3g)

        ! --- addition / subtraction ---
        rg = g1 + s1;  r3g = G3a + S3a;  call check_approx(ok, "2D2O + 2D2Osym", e2g(rg), r3g)
        rg = s1 + g1;  r3g = S3a + G3a;  call check_approx(ok, "2D2Osym + 2D2O", e2g(rg), r3g)
        rg = g1 - s1;  r3g = G3a - S3a;  call check_approx(ok, "2D2O - 2D2Osym", e2g(rg), r3g)
        rg = s1 - g1;  r3g = S3a - G3a;  call check_approx(ok, "2D2Osym - 2D2O", e2g(rg), r3g)
        rg = I2 + g1;  r3g = I2 + G3a;   call check_approx(ok, "I + 2D2O", e2g(rg), r3g)
        rg = g1 - I2;  r3g = G3a - I2;   call check_approx(ok, "2D2O - I", e2g(rg), r3g)
        rg = I2 - g1;  r3g = I2 - G3a;   call check_approx(ok, "I - 2D2O", e2g(rg), r3g)
        rg = I2S + g1; r3g = I2S + G3a;  call check_approx(ok, "cI + 2D2O", e2g(rg), r3g)
        rg = g1 - I2S; r3g = G3a - I2S;  call check_approx(ok, "2D2O - cI", e2g(rg), r3g)
        rg = I2S - g1; r3g = I2S - G3a;  call check_approx(ok, "cI - 2D2O", e2g(rg), r3g)

        ! --- assignments ---
        rg = s1;  r3g = S3a;  call check_approx(ok, "2D2O = 2D2Osym", e2g(rg), r3g)
        rs = I2;  r3s = I2;   call check_approx(ok, "2D2Osym = I", e2s(rs), r3s)
        rs = I2S; r3s = I2S;  call check_approx(ok, "2D2Osym = cI", e2s(rs), r3s)

        ! --- congruence transformation ---
        rs = g1 .transform. s1; r3s = G3a .transform. S3a
        call check_approx(ok, "2D2O .transform. 2D2Osym", e2s(rs), r3s)
    end subroutine second_order

    subroutine fourth_order()
        !! Fourth-order operators.
        type(ten_2D2Osym) :: rs
        type(ten_3D2Osym) :: r3s
        type(ten_2D4O3sym) :: r43
        type(ten_2D4O2sym) :: r42
        type(ten_3D4O3sym) :: r343
        type(ten_3D4O2sym) :: r342

        ! --- 4th order : 2nd order ---
        rs = C1 .ddot. s1;  r3s = C3a .ddot. S3a; call check_approx(ok, "2D4O3sym : 2D2Osym", e2s(rs), r3s)
        rs = s1 .ddot. C1;  r3s = S3a .ddot. C3a; call check_approx(ok, "2D2Osym : 2D4O3sym", e2s(rs), r3s)
        rs = D1 .ddot. s1;  r3s = D3 .ddot. S3a;  call check_approx(ok, "2D4O2sym : 2D2Osym", e2s(rs), r3s)
        rs = s1 .ddot. D1;  r3s = S3a .ddot. D3;  call check_approx(ok, "2D2Osym : 2D4O2sym", e2s(rs), r3s)

        ! --- 4th order : 4th order ---
        r42 = C1 .ddot. C2; r342 = C3a .ddot. C3b
        call check_approx(ok, "2D4O3sym : 2D4O3sym", r42, p42(r342))

        ! --- dyadic products ---
        r42 = s1 .tdot. s2;  r342 = S3a .tdot. S3b; call check_approx(ok, "2D2Osym .tdot. 2D2Osym", r42, p42(r342))
        r42 = s1 .tdot. I2;  r342 = S3a .tdot. I2;  call check_approx(ok, "2D2Osym .tdot. I", r42, p42(r342))
        r42 = I2 .tdot. s1;  r342 = I2 .tdot. S3a;  call check_approx(ok, "I .tdot. 2D2Osym", r42, p42(r342))
        r43 = s1 .tdotsym. s2; r343 = S3a .tdotsym. S3b
        call check_approx(ok, "2D2Osym .tdotsym. 2D2Osym", r43, p43(r343))
        r43 = .tdotsym. s1;    r343 = .tdotsym. S3a
        call check_approx(ok, ".tdotsym. 2D2Osym", r43, p43(r343))
        r43 = s1 .tdotsym. I2; r343 = S3a .tdotsym. I2
        call check_approx(ok, "2D2Osym .tdotsym. I", r43, p43(r343))
        r43 = I2 .tdotsym. s1; r343 = I2 .tdotsym. S3a
        call check_approx(ok, "I .tdotsym. 2D2Osym", r43, p43(r343))

        ! --- inverse (embedded with unit out-of-plane shear so that it is invertible) ---
        r43 = .inv. C1; r343 = .inv. e43(C1, .true.)
        call check_approx(ok, ".inv. 2D4O3sym", r43, p43(r343), tol=TOL)

        ! --- identities: +, - and assignment ---
        r43 = C1 + I3;  r343 = C3a + I3;  call check_approx(ok, "2D4O3sym + I3T", r43, p43(r343))
        r43 = I3S - C1; r343 = I3S - C3a; call check_approx(ok, "cI3T - 2D4O3sym", r43, p43(r343))
        r43 = C1 - I4;  r343 = C3a - I4;  call check_approx(ok, "2D4O3sym - I4T", r43, p43(r343))
        r43 = I4S + C1; r343 = I4S + C3a; call check_approx(ok, "cI4T + 2D4O3sym", r43, p43(r343))
        r42 = D1 + I3S; r342 = D3 + I3S;  call check_approx(ok, "2D4O2sym + cI3T", r42, p42(r342))
        r42 = I3 - D1;  r342 = I3 - D3;   call check_approx(ok, "I3T - 2D4O2sym", r42, p42(r342))
        r42 = D1 - I4S; r342 = D3 - I4S;  call check_approx(ok, "2D4O2sym - cI4T", r42, p42(r342))
        r42 = I4 + D1;  r342 = I4 + D3;   call check_approx(ok, "I4T + 2D4O2sym", r42, p42(r342))
        r43 = I3S;      r343 = I3S;       call check_approx(ok, "2D4O3sym = cI3T", r43, p43(r343))
        r43 = I4S;      r343 = I4S;       call check_approx(ok, "2D4O3sym = cI4T", r43, p43(r343))
        r42 = I4;       r342 = I4;        call check_approx(ok, "2D4O2sym = I4T", r42, p42(r342))
        r42 = I4S;      r342 = I4S;       call check_approx(ok, "2D4O2sym = cI4T", r42, p42(r342))
        r42 = C1;       r342 = C3a;       call check_approx(ok, "2D4O2sym = 2D4O3sym", r42, p42(r342))

        ! --- fourth-order congruence transformation ---
        r42 = g1 .transform. D1; r342 = G3a .transform. D3
        call check_approx(ok, "2D2O .transform. 2D4O2sym", r42, p42(r342), tol=TOL)
        ! Second general (non-symmetric, non-diagonal) A with a non-unit out-of-plane stretch
        r42 = g2 .transform. D1; r342 = G3b .transform. D3
        call check_approx(ok, "2D2O .transform. 2D4O2sym (A = g2)", r42, p42(r342), tol=TOL)

        call transform_rotation_analytic()
    end subroutine fourth_order

    subroutine transform_rotation_analytic()
        !! Analytic case: an in-plane rotation by 90 degrees (e1 -> e2, e2 -> -e1, e3 -> e3)
        !! maps C_1111 <-> C_2222, C_1133 <-> C_2233, C_1112 -> -C_2221 = -C_2212, etc.
        !! For R = [[0,-1,0],[1,0,0],[0,0,1]]: res_ijkl = R_iI R_jJ R_kK R_lL C_IJKL, so
        !! res_1111 = C_2222, res_2222 = C_1111, res_1122 = C_2211, res_3333 = C_3333,
        !! res_1133 = C_2233, res_2233 = C_1133, res_1212 = C_2121 = C_1212,
        !! res_1112 = -C_2221 = -C_2212, res_2212 = -C_1112, res_3312 = -C_3312.
        !! Entries are products of 0/+-1, so the result is exact.
        type(ten_2D2O) :: Rot
        type(ten_2D4O2sym) :: r42, expected
        real(real64) :: c(4,4)

        call Rot%init(xx=0D0, xy=-1D0, yx=1D0, yy=0D0, zz=1D0)
        c = D1%vals
        expected%vals(1,:) = (/ c(2,2), c(2,1), c(2,3), -c(2,4) /)
        expected%vals(2,:) = (/ c(1,2), c(1,1), c(1,3), -c(1,4) /)
        expected%vals(3,:) = (/ c(3,2), c(3,1), c(3,3), -c(3,4) /)
        expected%vals(4,:) = (/ -c(4,2), -c(4,1), -c(4,3), c(4,4) /)
        r42 = Rot .transform. D1
        call check_approx(ok, "2D2O .transform. 2D4O2sym (90 deg rotation)", r42, expected, tol=TOL)
    end subroutine transform_rotation_analytic


end program test_operator_2d_vs_3d
