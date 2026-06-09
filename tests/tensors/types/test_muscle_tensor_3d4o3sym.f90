program test_muscle_tensor_3d4o3sym
    use muscle_tensor_3d4o3sym
    implicit none
    logical :: passed

    print*, "Starting ten_3D4O3sym (Elasticity Tensors) tests..."

    call test_3D4O3sym_approx(passed)
    if (.not. passed) stop 1

    call test_3D4O3sym_arithmetic(passed)
    if (.not. passed) stop 2

    call test_3D4O3sym_inv(passed)
    if (.not. passed) stop 3

    print*, "All ten_3D4O3sym tests passed successfully!"
end program test_muscle_tensor_3d4o3sym

! ==========================================
! Test: Approximation and Initialization
! ==========================================
subroutine test_3D4O3sym_approx(passed)
    use muscle_tensor_3d4o3sym
    implicit none
    logical, intent(out) :: passed
    type(ten_3D4O3sym) :: t1, t2

    ! 1.- Identity check
    call t1%init(1D0, 2D0, 3D0, 4D0, 5D0, 6D0, 7D0, 8D0, 9D0, 10D0, &
                 11D0, 12D0, 13D0, 14D0, 15D0, 16D0, 17D0, 18D0, 19D0, 20D0, 21D0)
    t2 = t1
    passed = t1 .approx. t2
    if (.not. passed) then
        print*, "1.- Error: Approx failed for identical 4th order tensors"
        return
    end if

    ! 2.- Difference detection
    t2%vals(21) = 21.001D0
    passed = .not. (t1 .approx. t2)
    if (.not. passed) then
        print*, "2.- Error: Approx failed to detect difference in component 21"
        return
    end if
end subroutine

! ==========================================
! Test: Arithmetic (+, -, *, /)
! ==========================================
subroutine test_3D4O3sym_arithmetic(passed)
    use muscle_tensor_3d4o3sym
    implicit none
    logical, intent(out) :: passed
    type(ten_3D4O3sym) :: t1, t2, res, exp

    t1 = 1D0
    t2 = 2D0

    ! 1.- Sum
    res = t1 + t1
    passed = res .approx. t2
    if (.not. passed) then
        print*, "1.- Error: Sum of 4th order tensors failed"
        return
    end if

    ! 2.- Scalar Multiplication
    res = t1 * 2.0D0
    passed = res .approx. t2
    if (.not. passed) then
        print*, "2.- Error: Scalar multiplication failed"
        return
    end if

    ! 3.- Unary Minus
    res = -t1
    exp = -1D0
    passed = res .approx. exp
    if (.not. passed) then
        print*, "3.- Error: Unary minus failed"
        return
    end if
end subroutine

! ==========================================
! Test: Tensor Inverse (.inv.)
! ==========================================
subroutine test_3D4O3sym_inv(passed)
    use, intrinsic :: iso_fortran_env
    use muscle_tensor_3d4o3sym
    implicit none
    logical, intent(out) :: passed
    type(ten_3D4O3sym) :: C_iso, S_compliance
    real(real64) :: lambda, mu, val_xxxx, val_xxyy
    real(real64) :: tmp

    ! Define an isotropic elasticity tensor
    ! C_ijkl = lambda*d_ij*d_kl + mu*(d_ik*d_jl + d_il*d_jk)
    lambda = 100.0D0
    mu = 50.0D0

    ! In Voigt notation for isotropic:
    ! C11 = lambda + 2*mu = 200
    ! C12 = lambda        = 100
    ! C44 = mu            = 50 (this is the xyxy component)
    
    C_iso = 0D0
    ! Diagonals (xxxx, yyyy, zzzz)
    ! Shear mu (xyxy, yzyz, xzxz)
    ! Off-diagonals lambda (xxyy, yyzz, xxzz)
    tmp = lambda + 2.0D0 * mu
    call C_iso%set(xxxx=tmp, yyyy=tmp, zzzz=tmp, &
                   xyxy=mu, yzyz=mu, xzxz=mu, &
                   xxyy=lambda, yyzz=lambda, xxzz=lambda)

    ! Compute Inverse (Compliance tensor S)
    S_compliance = .inv. C_iso

    ! Verification: S_1111 should be (lambda+mu)/(mu*(3*lambda+2*mu))
    ! For lambda=100, mu=50 -> S11 = 150 / (50 * 400) = 150 / 20000 = 0.0075
    val_xxxx = (lambda + mu) / (mu * (3.0D0 * lambda + 2.0D0 * mu))
    
    passed = abs(S_compliance%vals(1) - val_xxxx) < 1D-8
    if (.not. passed) then
        print*, "1.- Error: Inverse component S_1111 incorrect"
        print*, "Expected:", val_xxxx, " Actual:", S_compliance%vals(1)
        return
    end if

    ! Verification: S_1122 should be -lambda / (2*mu*(3*lambda+2*mu))
    ! S12 = -100 / (100 * 400) = -0.0025
    val_xxyy = -lambda / (2.0D0 * mu * (3.0D0 * lambda + 2.0D0 * mu))
    passed = abs(S_compliance%vals(7) - val_xxyy) < 1D-8
    if (.not. passed) then
        print*, "2.- Error: Inverse component S_1122 incorrect"
        print*, "Expected:", val_xxyy, " Actual:", S_compliance%vals(7)
        return
    end if
end subroutine