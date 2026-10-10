program test_muscle_tensor_2d4o3sym
    use muscle_tensors
    implicit none
    
    logical :: passed

    call test_ten_2D4O3sym_approx(passed)
    if (.not. passed) STOP 1

    call test_ten_2D4O3sym_sum(passed)
    if (.not. passed) STOP 2

    call test_ten_2D4O3sym_sub(passed)
    if (.not. passed) STOP 3

    call test_ten_2D4O3sym_mul(passed)
    if (.not. passed) STOP 4

    call test_ten_2D4O3sym_div(passed)
    if (.not. passed) STOP 5

    call test_ten_2D4O3sym_ddot(passed)
    if (.not. passed) STOP 6

    call test_ten_2D4O3sym_tdot(passed)
    if (.not. passed) STOP 7

    call test_ten_2D4O3sym_inv_isotropic(passed)
    if (.not. passed) STOP 8

    call test_ten_2D4O3sym_set_norm_assign(passed)
    if (.not. passed) STOP 9

    call test_ten_2D4O3sym_identities(passed)
    if (.not. passed) STOP 10

    call test_ten_2D4O3sym_ddot_4O(passed)
    if (.not. passed) STOP 11

    STOP 0
end program test_muscle_tensor_2d4o3sym

subroutine test_ten_2D4O3sym_approx(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    
    logical, intent(out) :: passed

    type(ten_2D4O3sym) :: to_test1, to_test2

    passed = .true.
    call to_test1%init(xxxx=1D0, yyyy=2D0, zzzz=3D0, xyxy=4D0, &
                       xxyy=5D0, yyzz=6D0, zzxy=7D0,           &
                       xxzz=8D0, yyxy=9D0,                     &
                       xxxy=10D0                               &
                       )
    call to_test2%init(xxxx=1D0, yyyy=2D0, zzzz=3D0, xyxy=4D0, &
                       xxyy=5D0, yyzz=6D0, zzxy=7D0,           &
                       xxzz=8D0, yyxy=9D0,                     &
                       xxxy=10D0                               &
                       )
    call check_true(passed, ".approx. of two identical tensors", to_test1 .approx. to_test2)

    call to_test2%init(xxxx=1D0, yyyy=2D0, zzzz=3D0, xyxy=4D0, &
                       xxyy=5D0, yyzz=6D0, zzxy=7D0,           &
                       xxzz=8D0, yyxy=9D0,                     &
                       xxxy=10.001D0                           &
                       )
    call check_true(passed, ".approx. rejects xxxy = 10 vs 10.001", .not. (to_test1 .approx. to_test2))

    to_test1 = 0D0
    to_test2 = 0D0
    call check_true(passed, ".approx. of two zero tensors", to_test1 .approx. to_test2)
end subroutine

subroutine test_ten_2D4O3sym_sum(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    
    logical, intent(out) :: passed

    type(ten_2D4O3sym) :: to_test1, to_test2
    type(ten_2D4O3sym) :: expected_result

    passed = .true.
    to_test1 = 0D0
    call check_approx(passed, "0 + 0", to_test1 + to_test1, to_test1)
    
    call to_test2%init(xxxx=1D0, yyyy=2D0, zzzz=3D0, xyxy=4D0, &
                       xxyy=5D0, yyzz=6D0, zzxy=7D0,           &
                       xxzz=8D0, yyxy=9D0,                     &
                       xxxy=10D0                               &
                       )
    call check_approx(passed, "0 + b", to_test1 + to_test2, to_test2)
    call check_approx(passed, "b + 0", to_test2 + to_test1, to_test2)

    call to_test1%init(xxxx=1D1, yyyy=2D1, zzzz=3D1, xyxy=4D1, &
                       xxyy=5D1, yyzz=6D1, zzxy=7D1,           &
                       xxzz=8D1, yyxy=9D1,                     &
                       xxxy=10D1                               &
                       )
    call expected_result%init(xxxx=11D0, yyyy=22D0, zzzz=33D0, xyxy=44D0, &
                              xxyy=55D0, yyzz=66D0, zzxy=77D0,            &
                              xxzz=88D0, yyxy=99D0,                       &
                              xxxy=110D0                                  &
                              )
    call check_approx(passed, "b + a", to_test2 + to_test1, expected_result)
    call check_approx(passed, "a + b", to_test1 + to_test2, expected_result)
end subroutine

subroutine test_ten_2D4O3sym_sub(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    
    logical, intent(out) :: passed

    type(ten_2D4O3sym) :: to_test1, to_test2
    type(ten_2D4O3sym) :: expected_result

    passed = .true.
    to_test1 = 0D0
    call check_approx(passed, "0 - 0", to_test1 - to_test1, to_test1)
    
    call to_test2%init(xxxx=1D0, yyyy=2D0, zzzz=3D0, xyxy=4D0, &
                       xxyy=5D0, yyzz=6D0, zzxy=7D0,           &
                       xxzz=8D0, yyxy=9D0,                     &
                       xxxy=10D0                               &
                       )
    call check_approx(passed, "b - 0", to_test2 - to_test1, to_test2)
    call check_approx(passed, "0 - b", to_test1 - to_test2, -to_test2)

    call to_test1%init(xxxx=11D0, yyyy=12D0, zzzz=13D0, xyxy=14D0, &
                       xxyy=15D0, yyzz=16D0, zzxy=17D0,            &
                       xxzz=18D0, yyxy=19D0,                       &
                       xxxy=20D0                                   &
                       )
    expected_result = 10D0
    call check_approx(passed, "a - b", to_test1 - to_test2, expected_result)
    call check_approx(passed, "b - a", to_test2 - to_test1, -expected_result)
end subroutine

subroutine test_ten_2D4O3sym_mul(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    
    logical, intent(out) :: passed

    type(ten_2D4O3sym) :: to_test1
    type(ten_2D4O3sym) :: expected_result

    passed = .true.
    to_test1 = 0D0
    call check_approx(passed, "10*0", 10D0*to_test1, to_test1)
    call check_approx(passed, "0*10", to_test1*10D0, to_test1)

    call to_test1%init(xxxx=1D0, yyyy=2D0, zzzz=3D0, xyxy=4D0, &
                       xxyy=5D0, yyzz=6D0, zzxy=7D0,           &
                       xxzz=8D0, yyxy=9D0,                     &
                       xxxy=10D0                               &
                       )
    expected_result = 0D0
    call check_approx(passed, "0*a", 0D0*to_test1, expected_result)
    call check_approx(passed, "a*0", to_test1*0D0, expected_result)

    call expected_result%init(xxxx=2D0, yyyy=4D0, zzzz=6D0, xyxy=8D0, &
                              xxyy=10D0, yyzz=12D0, zzxy=14D0,        &
                              xxzz=16D0, yyxy=18D0,                   &
                              xxxy=20D0                               &
                              )
    call check_approx(passed, "2*a", 2D0*to_test1, expected_result)
    call check_approx(passed, "a*2", to_test1*2D0, expected_result)
end subroutine

subroutine test_ten_2D4O3sym_div(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    
    logical, intent(out) :: passed

    type(ten_2D4O3sym) :: to_test1
    type(ten_2D4O3sym) :: expected_result

    passed = .true.
    to_test1 = 0D0
    call check_approx(passed, "0/2", to_test1/2D0, to_test1)

    call to_test1%init(xxxx=2D0, yyyy=4D0, zzzz=6D0, xyxy=8D0, &
                       xxyy=10D0, yyzz=12D0, zzxy=14D0,        &
                       xxzz=16D0, yyxy=18D0,                   &
                       xxxy=20D0                               &
                       )
    call expected_result%init(xxxx=1D0, yyyy=2D0, zzzz=3D0, xyxy=4D0, &
                              xxyy=5D0, yyzz=6D0, zzxy=7D0,           &
                              xxzz=8D0, yyxy=9D0,                     &
                              xxxy=10D0                               &
                              )
    call check_approx(passed, "a/2", to_test1/2D0, expected_result)
end subroutine

subroutine test_ten_2D4O3sym_ddot(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    
    logical, intent(out) :: passed

    type(ten_2D4O3sym) :: to_test1
    type(ten_2D2Osym) :: to_test2, expected_result
    
    passed = .true.
    to_test1 = 0D0
    call to_test2%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0)
    expected_result = 0D0
    call check_approx(passed, "0 .ddot. x", to_test1 .ddot. to_test2, expected_result)
    call check_approx(passed, "x .ddot. 0", to_test2 .ddot. to_test1, expected_result)

    call to_test1%init(xxxx=1D0, yyyy=2D0, zzzz=3D0, xyxy=4D0, &
                       xxyy=5D0, yyzz=6D0, zzxy=7D0,           &
                       xxzz=8D0, yyxy=9D0,                     &
                       xxxy=10D0                                &
                       )
    call expected_result%init(xx=115D0, yy=99D0, zz=85D0, xy=81D0)
    call check_approx(passed, "C .ddot. x", to_test1 .ddot. to_test2, expected_result)
    call check_approx(passed, "x .ddot. C", to_test2 .ddot. to_test1, expected_result)
end subroutine


subroutine test_ten_2D4O3sym_tdot(passed)
    !! Dyadic products with a = (1,2,3,4), b = (5,6,7,8) in Voigt (xx,yy,zz,xy):
    !! - a .tdot. b is the 4x4 matrix C_IJ = a_I b_J (ten_2D4O2sym, no major symmetry);
    !! - a .tdotsym. b = (a(x)b + b(x)a)/2 and .tdotsym. a = a(x)a (ten_2D4O3sym);
    !! - a .tdotsym. I  = (a(x)I + I(x)a)/2.
    !! All entries are exact products of small integers.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D2Osym) :: a, b
    type(ten_2D4O2sym) :: ab, ab_sym_full
    type(ten_2D4O3sym) :: expected
    type(iden_2O) :: I2
    real(real64) :: m(4,4), expected_m(4,4)
    integer :: i, j

    passed = .true.
    call a%init(xx=1D0, yy=2D0, zz=3D0, xy=4D0)
    call b%init(xx=5D0, yy=6D0, zz=7D0, xy=8D0)

    ab = a .tdot. b
    do j = 1, 4
        do i = 1, 4
            m(i,j) = a%vals(i)*b%vals(j)
        end do
    end do
    call check_approx(passed, "a .tdot. b", ab%vals, m, tol=1D-15)

    ab_sym_full = a .tdotsym. b
    call check_approx(passed, "a .tdotsym. b", ab_sym_full%vals, 0.5D0*(m + transpose(m)), tol=1D-15)

    ! a(x)a: (11,22,33,44,12,23,34,13,24,14)
    call expected%init((/1D0, 4D0, 9D0, 16D0, 2D0, 6D0, 12D0, 3D0, 8D0, 4D0/))
    call check_approx(passed, ".tdotsym. a", .tdotsym. a, expected)
    call check_approx(passed, "a .tdotsym. a", a .tdotsym. a, expected)

    ! (a(x)I + I(x)a)/2: normal-normal block (a_I + a_J)/2, normal-shear a_xy/2, shear-shear 0
    call expected%init((/1D0, 2D0, 3D0, 0D0, 1.5D0, 2.5D0, 2D0, 2D0, 2D0, 2D0/))
    call check_approx(passed, "a .tdotsym. I", a .tdotsym. I2, expected)
    call check_approx(passed, "I .tdotsym. a", I2 .tdotsym. a, expected)

    ! a .tdot. I fills the three normal columns with a; I .tdot. a is its transpose
    expected_m = 0D0
    expected_m(:,1:3) = spread(a%vals, 2, 3)
    ab = a .tdot. I2
    call check_approx(passed, "a .tdot. I", ab%vals, expected_m, tol=0D0)
    ab = I2 .tdot. a
    call check_approx(passed, "I .tdot. a", ab%vals, transpose(expected_m), tol=0D0)
end subroutine test_ten_2D4O3sym_tdot


subroutine test_ten_2D4O3sym_inv_isotropic(passed)
    !! Analytic solution: the inverse of the isotropic stiffness
    !! C = lambda I(x)I + 2 mu I^S is the isotropic compliance
    !! S_1111 = 1/E, S_1122 = -nu/E, S_1212 = (1+nu)/(2E)  (tensorial components).
    !! Built here with E = 210, nu = 0.3 through the identity-assignment/sum operators.
    !! Tolerance: one 4x4 LU inversion of a well-conditioned matrix (cond ~ 10),
    !! error ~ 1e-15 relative, so 1e-12 is a safe margin.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D4O3sym) :: C, S, expected
    type(ten_2D4O2sym) :: SC, Isym
    type(iden_4O3TS) :: I3S
    type(iden_4O4TS) :: I4S
    type(iden_4O4T) :: I4
    real(real64), parameter :: E = 210D0, NU = 0.3D0
    real(real64) :: lambda, mu

    passed = .true.
    lambda = E*NU/((1D0 + NU)*(1D0 - 2D0*NU))
    mu = E/(2D0*(1D0 + NU))
    I4S%val = 2D0*mu
    I3S%val = lambda
    C = I4S
    C = C + I3S

    call expected%init(xxxx=1D0/E, yyyy=1D0/E, zzzz=1D0/E, xyxy=(1D0+NU)/(2D0*E), &
                       xxyy=-NU/E, yyzz=-NU/E, zzxy=0D0, xxzz=-NU/E, yyxy=0D0, xxxy=0D0)
    S = .inv. C
    call check_approx(passed, ".inv. C (isotropic compliance)", S, expected, tol=1D-12)

    ! S : C = I^S (minor-symmetric result)
    SC = S .ddot. C
    Isym = I4
    call check_approx(passed, "S .ddot. C = I^S", SC, Isym, tol=1D-12)
end subroutine test_ten_2D4O3sym_inv_isotropic


subroutine test_ten_2D4O3sym_set_norm_assign(passed)
    !! set modifies only the named components; norm = sum|diag| + 2 sum|offdiag|.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D4O3sym) :: a, expected

    passed = .true.
    a = 0D0
    call a%set(xyxy=3D0, xxzz=-2D0)
    call expected%init((/0D0, 0D0, 0D0, 3D0, 0D0, 0D0, 0D0, -2D0, 0D0, 0D0/))
    call check_approx(passed, "set(xyxy=3, xxzz=-2)", a, expected)
    call check_approx(passed, "norm = 3 + 2*2", a%norm(), 7D0, tol=1D-14)
    a = 1.5D0
    call check_approx(passed, "a = 1.5 (scalar assignment)", a%vals, spread(1.5D0, 1, 10), tol=0D0)
end subroutine test_ten_2D4O3sym_set_norm_assign


subroutine test_ten_2D4O3sym_identities(passed)
    !! Identity assignments and +/- with I(x)I and I^S, on both 4th-order 2D types.
    !! I(x)I adds 1 to the 3x3 normal block; I^S adds 1 to the normal diagonal and 1/2 to xyxy.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D4O3sym) :: a, r, expected
    type(ten_2D4O2sym) :: b, rb, expb
    type(iden_4O3T) :: I3
    type(iden_4O3TS) :: I3S
    type(iden_4O4T) :: I4
    type(iden_4O4TS) :: I4S

    passed = .true.
    I3S%val = 2D0
    I4S%val = 2D0

    a = I3S
    call expected%init((/2D0, 2D0, 2D0, 0D0, 2D0, 2D0, 0D0, 2D0, 0D0, 0D0/))
    call check_approx(passed, "2D4O3sym = 2 I(x)I", a, expected)
    a = I4S
    call expected%init((/2D0, 2D0, 2D0, 1D0, 0D0, 0D0, 0D0, 0D0, 0D0, 0D0/))
    call check_approx(passed, "2D4O3sym = 2 I^S", a, expected)

    a = 1D0
    call expected%init((/2D0, 2D0, 2D0, 1D0, 2D0, 2D0, 1D0, 2D0, 1D0, 1D0/))
    call check_approx(passed, "a + I(x)I", a + I3, expected)
    call check_approx(passed, "I(x)I + a", I3 + a, expected)
    call check_approx(passed, "I(x)I - a", I3 - a, expected - 2D0*a)
    call check_approx(passed, "a - I(x)I", a - I3, 2D0*a - expected)
    call expected%init((/3D0, 3D0, 3D0, 1D0, 3D0, 3D0, 1D0, 3D0, 1D0, 1D0/))
    call check_approx(passed, "a + 2I(x)I", a + I3S, expected)
    call check_approx(passed, "2I(x)I + a", I3S + a, expected)
    call check_approx(passed, "a - 2I(x)I", a - I3S, 2D0*a - expected)
    call check_approx(passed, "2I(x)I - a", I3S - a, expected - 2D0*a)
    call expected%init((/2D0, 2D0, 2D0, 1.5D0, 1D0, 1D0, 1D0, 1D0, 1D0, 1D0/))
    call check_approx(passed, "a + I^S", a + I4, expected)
    call check_approx(passed, "I^S + a", I4 + a, expected)
    call check_approx(passed, "a - I^S", a - I4, 2D0*a - expected)
    call check_approx(passed, "I^S - a", I4 - a, expected - 2D0*a)
    call expected%init((/3D0, 3D0, 3D0, 2D0, 1D0, 1D0, 1D0, 1D0, 1D0, 1D0/))
    call check_approx(passed, "a + 2I^S", a + I4S, expected)
    call check_approx(passed, "2I^S + a", I4S + a, expected)
    call check_approx(passed, "a - 2I^S", a - I4S, 2D0*a - expected)
    call check_approx(passed, "2I^S - a", I4S - a, expected - 2D0*a)

    ! Same operations on the minor-symmetric type, checked against the expanded 2D4O3sym
    b%vals = 1D0
    r = a + I3;  expb = r; rb = b + I3;  call check_approx(passed, "2D4O2sym + I(x)I", rb, expb)
                           rb = I3 + b;  call check_approx(passed, "I(x)I + 2D4O2sym", rb, expb)
    r = a + I3S; expb = r; rb = b + I3S; call check_approx(passed, "2D4O2sym + 2I(x)I", rb, expb)
                           rb = I3S + b; call check_approx(passed, "2I(x)I + 2D4O2sym", rb, expb)
    r = a + I4;  expb = r; rb = b + I4;  call check_approx(passed, "2D4O2sym + I^S", rb, expb)
                           rb = I4 + b;  call check_approx(passed, "I^S + 2D4O2sym", rb, expb)
    r = a + I4S; expb = r; rb = b + I4S; call check_approx(passed, "2D4O2sym + 2I^S", rb, expb)
                           rb = I4S + b; call check_approx(passed, "2I^S + 2D4O2sym", rb, expb)
    r = a - I3;  expb = r; rb = b - I3;  call check_approx(passed, "2D4O2sym - I(x)I", rb, expb)
    r = I3 - a;  expb = r; rb = I3 - b;  call check_approx(passed, "I(x)I - 2D4O2sym", rb, expb)
    r = a - I3S; expb = r; rb = b - I3S; call check_approx(passed, "2D4O2sym - 2I(x)I", rb, expb)
    r = I3S - a; expb = r; rb = I3S - b; call check_approx(passed, "2I(x)I - 2D4O2sym", rb, expb)
    r = a - I4;  expb = r; rb = b - I4;  call check_approx(passed, "2D4O2sym - I^S", rb, expb)
    r = I4 - a;  expb = r; rb = I4 - b;  call check_approx(passed, "I^S - 2D4O2sym", rb, expb)
    r = a - I4S; expb = r; rb = b - I4S; call check_approx(passed, "2D4O2sym - 2I^S", rb, expb)
    r = I4S - a; expb = r; rb = I4S - b; call check_approx(passed, "2I^S - 2D4O2sym", rb, expb)
    rb = I4;     expb = I4;              call check_approx(passed, "2D4O2sym = I^S", rb%vals, expb%vals, tol=0D0)
    r = I4S;     expb = r; rb = I4S;     call check_approx(passed, "2D4O2sym = 2I^S", rb, expb)
end subroutine test_ten_2D4O3sym_identities


subroutine test_ten_2D4O3sym_ddot_4O(passed)
    !! A:B for fully symmetric 2D tensors equals the Voigt product A W B, W = diag(1,1,1,2),
    !! and the minor-symmetric : second-order contractions match the 2D4O3sym ones.
    use, intrinsic :: iso_fortran_env
    use muscle_tensors
    use muscle_test_utils
    implicit none
    logical, intent(out) :: passed
    type(ten_2D4O3sym) :: a, b
    type(ten_2D4O2sym) :: ab, a_full, b_full
    type(ten_2D2Osym) :: x
    real(real64) :: w(4,4)

    passed = .true.
    call a%init((/1D0, 2D0, 3D0, 4D0, 5D0, 6D0, 7D0, 8D0, 9D0, 10D0/))
    call b%init((/-1D0, 0.5D0, 2D0, 1D0, 3D0, -2D0, 0D0, 1D0, 4D0, -3D0/))
    a_full = a
    b_full = b
    w = 0D0
    w(1,1) = 1D0; w(2,2) = 1D0; w(3,3) = 1D0; w(4,4) = 2D0
    ab = a .ddot. b
    call check_approx(passed, "2D4O3sym .ddot. 2D4O3sym = A W B", ab%vals, &
                      matmul(a_full%vals, matmul(w, b_full%vals)), tol=1D-14)

    call x%init(xx=1D0, yy=-2D0, zz=0.5D0, xy=3D0)
    call check_approx(passed, "2D4O2sym .ddot. 2D2Osym", a_full .ddot. x, a .ddot. x)
    call check_approx(passed, "2D2Osym .ddot. 2D4O2sym", x .ddot. a_full, x .ddot. a)
end subroutine test_ten_2D4O3sym_ddot_4O
