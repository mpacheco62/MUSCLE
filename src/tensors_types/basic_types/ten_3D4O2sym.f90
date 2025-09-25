module mod_ten_3D4O2sym
    !! Module mod_ten_3D4O2sym
    !! =======================
    !!
    !! Defines the type for 3D fourth-order tensors with minor symmetries and associated operations.
    !!
    !! This module provides the derived type `ten_3D4O2sym` to represent a fourth-order
    !! tensor in three dimensions possessing minor symmetries, i.e., \(C_{ijkl} = C_{jikl} = C_{ijlk}\).
    !! Such tensors commonly appear in constitutive models (e.g., elasticity tensors before
    !! assuming major symmetry).
    !!
    !! The tensor is stored internally using a 6x6 matrix based on the Voigt index mapping
    !! (11->1, 22->2, 33->3, 12->4, 23->5, 13->6). Note that due to only having minor
    !! symmetries, this 6x6 matrix is generally **not** symmetric (\(C_{IJ} \neq C_{JI}\)).
    !!
    !! The module overloads standard arithmetic operators (+, -, *, /), a custom
    !! equality comparison operator (.isequal.), and provides an inverse function (`inv_3D4O2sym`).
    !! It also provides methods for initialization.
    !!
    !! Public Entities
    !! ---------------
    !!
    !! ### Derived Type:
    !!
    !! - `ten_3D4O2sym`: Represents a 3D fourth-order tensor with minor symmetries.
    !!     - Component: `vals(6,6) :: real(real64)` - Stores the 36 independent components
    !!       in a 6x6 matrix using Voigt index mapping.
    !!     - Generic Procedure: `init` - Initializes the tensor either from a
    !!       6x6 array or from 36 individual components.
    !!
    !! ### Operators:
    !!
    !! - `.isequal.`: Compares two `ten_3D4O2sym` tensors for approximate equality using an L1 norm.
    !! - `+`: Adds two `ten_3D4O2sym` tensors.
    !! - `-`: Subtracts two `ten_3D4O2sym` tensors (binary) or computes the unary negation.
    !! - `*`: Multiplies a `ten_3D4O2sym` tensor by a `real(real64)` scalar (or vice-versa).
    !! - `/`: Divides a `ten_3D4O2sym` tensor by a `real(real64)` scalar.
    !! - `.inv.`: (Interface commented out, but function `inv_3D4O2sym` exists) Computes the inverse of the tensor based on its 6x6 matrix representation.
    !!
    !! Usage
    !! -----
    !!
    !! ```fortran
    !! program example_ten_3d4o2sym_usage
    !!   use mod_ten_3D4O2sym
    !!   use iso_fortran_env, only: real64
    !!   implicit none
    !!
    !!   type(ten_3D4O2sym) :: tensor_c, tensor_d, tensor_sum
    !!   real(real64) :: C_matrix(6,6)
    !!   logical :: are_equal
    !!
    !!   ! Initialize using a 6x6 array (Voigt notation)
    !!   ! (Example: Isotropic elasticity tensor - note it also has major symmetry)
    !!   real(real64) :: lambda, mu
    !!   lambda = 100.0D0
    !!   mu = 50.0D0
    !!   C_matrix = 0.0D0
    !!   C_matrix(1:3, 1:3) = lambda  ! Lambda in upper-left 3x3 block
    !!   C_matrix(1,1) = lambda + 2.0*mu
    !!   C_matrix(2,2) = lambda + 2.0*mu
    !!   C_matrix(3,3) = lambda + 2.0*mu
    !!   C_matrix(4,4) = mu
    !!   C_matrix(5,5) = mu
    !!   C_matrix(6,6) = mu
    !!   call tensor_c%init(C_matrix)
    !!
    !!   ! Initialize tensor_d to zero
    !!   tensor_d%vals = 0.0D0
    !!   ! Set a component
    !!   tensor_d%vals(1, 4) = 10.0 ! C_1112 = 10.0
    !!
    !!   ! Operations
    !!   tensor_sum = tensor_c + tensor_d * 2.0D0
    !!
    !!   ! Comparison
    !!   are_equal = (tensor_c .isequal. tensor_sum)
    !!
    !!   print *, "Tensor C(1,1) (C_1111):", tensor_c%vals(1,1)
    !!   print *, "Tensor C(4,4) (C_1212):", tensor_c%vals(4,4)
    !!   print *, "Tensor Sum(1,4) (C_1112):", tensor_sum%vals(1,4)
    !!   print *, "Are C and Sum equal?", are_equal
    !!
    !!   ! Inverse (using the function directly)
    !!   ! type(ten_3D4O2sym) :: tensor_c_inv
    !!   ! tensor_c_inv = inv_3D4O2sym(tensor_c)
    !!   ! print *, "Inv C(1,1):", tensor_c_inv%vals(1,1)
    !!
    !! end program example_ten_3d4o2sym_usage
    !! ```
    !!
    !! For more information see [[tensors_types]]

    use, intrinsic :: iso_fortran_env
    implicit none
    private
 
    type, public :: ten_3D4O2sym
        !! 3D Fourth-Order Tensor with Minor Symmetries (6x6 Voigt Storage)
        !! ================================================================
        !!
        !! Represents a fourth-order tensor \(C_{ijkl}\) in three dimensions that possesses
        !! minor symmetries:
        !! \[ C_{ijkl} = C_{jikl} = C_{ijlk} \]
        !! This implies that the tensor is symmetric with respect to the first two indices
        !! and the last two indices separately. It does **not** imply major symmetry
        !! (\(C_{ijkl} = C_{klij}\)).
        !!
        !! Storage:
        !! --------
        !! The 36 independent components resulting from minor symmetries are stored
        !! internally in a 6x6 `real(real64)` matrix `vals`. This uses the Voigt index
        !! mapping convention:
        !! - 1 <-> (1,1) or xx
        !! - 2 <-> (2,2) or yy
        !! - 3 <-> (3,3) or zz
        !! - 4 <-> (1,2) or (2,1) or xy
        !! - 5 <-> (2,3) or (3,2) or yz
        !! - 6 <-> (1,3) or (3,1) or xz
        !!
        !! The component \(C_{ijkl}\) is mapped to `vals(I, J)`, where \(I\) corresponds to the
        !! pair \((ij)\) and \(J\) corresponds to the pair \((kl)\) according to the Voigt mapping.
        !! For example, \(C_{1122}\) is stored in `vals(1, 2)`, and \(C_{1223}\) is stored in `vals(4, 5)`.
        !!
        !! **Important**: Due to the lack of major symmetry, the 6x6 matrix `vals` is generally
        !! **not** symmetric, i.e., `vals(I, J)` is not necessarily equal to `vals(J, I)`.
        !!
        !! Access:
        !! -------
        !! Components are typically accessed directly via the `vals` array using the
        !! appropriate Voigt indices (I, J).
        !!
        !! Initialization:
        !! ---------------
        !! Use the generic `init` procedure to initialize either from a 6x6 `real(real64)`
        !! array (corresponding directly to the `vals` matrix) or by providing the 36 components
        !! individually (see `init2_ten_3D4O2sym` for the required input order).
        !!
        !! For more information see [[tensors_types]]
        real(real64), dimension(6,6) :: vals
            !! Stores the 36 independent components as a 6x6 matrix using Voigt index mapping.
        contains
            generic, public :: init => init_ten_3D4O2sym, init2_ten_3D4O2sym
                !! Generic interface for initialization.
            procedure, private :: init_ten_3D4O2sym, init2_ten_3D4O2sym 
            procedure, public :: norm => norm_3D4O2sym
    end type ten_3D4O2sym

    public :: operator(.isequal.)
    interface operator (.isequal.)
        module procedure isequal_3D4O2sym
    end interface

    ! public :: operator(.inv.)
    ! interface operator (.inv.)
    !     module procedure inv_3D4O2sym
    ! end interface

    public :: operator(+)
    interface operator (+)
        module procedure sum_3D4O2sym
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_3D4O2sym
        module procedure subU_3D4O2sym
    end interface

    public :: operator(*)
    interface operator (*)
        module procedure mul_real64_3D4O2sym
        module procedure mul_3D4O2sym_real64
    end interface

    public :: operator( / )
    interface operator ( / )
        module procedure div_3D4O2sym_real64
    end interface
contains
    
    pure subroutine init_ten_3D4O2sym(self, vals)
        !! Initializes a ten_3D4O2sym tensor from a 6x6 array (Voigt matrix).
        implicit none
        class(ten_3D4O2sym), intent(inout) :: self
        real(real64), intent(in) :: vals(6,6)
        self%vals = vals
    end subroutine init_ten_3D4O2sym

    pure subroutine init2_ten_3D4O2sym(self,                               & 
                                              xxxx, xxyy, xxzz, xxxy, xxyz, xxxz, &
                                              yyxx, yyyy, yyzz, yyxy, yyyz, yyxz, &
                                              zzxx, zzyy, zzzz, zzxy, zzyz, zzxz, &
                                              xyxx, xyyy, xyzz, xyxy, xyyz, xyxz, &
                                              yzxx, yzyy, yzzz, yzxy, yzyz, yzxz, &
                                              xzxx, xzyy, xzzz, xzxy, xzyz, xzxz  &
                                              )
        !! Initializes a ten_3D4O2sym tensor from its 36 individual components.
        !! Input arguments correspond row-wise to the 6x6 Voigt matrix components.
        !! E.g., xxxx=C(1,1), xxyy=C(1,2), ..., xxxz=C(1,6), yyxx=C(2,1), ... xzxz=C(6,6).
        !
        !  Voigt Matrix Layout (I, J):
        !  | (1,1) (1,2) (1,3) (1,4) (1,5) (1,6) |   <- xxxx, xxyy, xxzz, xxxy, xxyz, xxxz
        !  | (2,1) (2,2) (2,3) (2,4) (2,5) (2,6) |   <- yyxx, yyyy, yyzz, yyxy, yyyz, yyxz
        !  | (3,1) (3,2) (3,3) (3,4) (3,5) (3,6) |   <- zzxx, zzyy, zzzz, zzxy, zzyz, zzxz
        !  | (4,1) (4,2) (4,3) (4,4) (4,5) (4,6) |   <- xyxx, xyyy, xyzz, xyxy, xyyz, xyxz
        !  | (5,1) (5,2) (5,3) (5,4) (5,5) (5,6) |   <- yzxx, yzyy, yzzz, yzxy, yzyz, yzxz
        !  | (6,1) (6,2) (6,3) (6,4) (6,5) (6,6) |   <- xzxx, xzyy, xzzz, xzxy, xzyz, xzxz
        implicit none
        class(ten_3D4O2sym), intent(inout) :: self
        real(real64), intent(in) :: xxxx, xxyy, xxzz, xxxy, xxyz, xxxz
        real(real64), intent(in) :: yyxx, yyyy, yyzz, yyxy, yyyz, yyxz
        real(real64), intent(in) :: zzxx, zzyy, zzzz, zzxy, zzyz, zzxz
        real(real64), intent(in) :: xyxx, xyyy, xyzz, xyxy, xyyz, xyxz
        real(real64), intent(in) :: yzxx, yzyy, yzzz, yzxy, yzyz, yzxz
        real(real64), intent(in) :: xzxx, xzyy, xzzz, xzxy, xzyz, xzxz
        self%vals(:,1) = (/ xxxx, yyxx, zzxx, xyxx, yzxx, xzxx /)
        self%vals(:,2) = (/ xxyy, yyyy, zzyy, xyyy, yzyy, xzyy /)
        self%vals(:,3) = (/ xxzz, yyzz, zzzz, xyzz, yzzz, xzzz /)
        self%vals(:,4) = (/ xxxy, yyxy, zzxy, xyxy, yzxy, xzxy /)
        self%vals(:,5) = (/ xxyz, yyyz, zzyz, xyyz, yzyz, xzyz /)
        self%vals(:,6) = (/ xxxz, yyxz, zzxz, xyxz, yzxz, xzxz /)
    end subroutine

    pure function isequal_3D4O2sym(a, b) result(res)
        !! `.isequal.` Compares two ten_3D4O2sym tensors for approximate equality.
        !! Uses the L1 norm of the difference of the 6x6 Voigt matrices with relative
        !! and absolute tolerances (EPS, EPS_ABS).
        !! norm(A) = sum(|A_IJ|) for I,J=1..6
        !! Condition: norm(a-b) / max(norm(a), norm(b), EPS_ABS) <= EPS
        implicit none
        class(ten_3D4O2sym), intent(in) :: a, b
        type(ten_3D4O2sym) :: temp
        logical :: res
        real(real64), parameter :: EPS=1e-7, EPS_ABS=1e-30
        real(real64) :: norm_a, norm_b, norm_max, norm
        integer :: i, j

        norm_a = 0D0
        norm_b = 0D0
        do i=1,6
            do j=1,6
                norm_a = norm_a + abs(a%vals(i,j))
                norm_b = norm_b + abs(b%vals(i,j))
            end do
        end do

        norm_max = max(max(norm_a, norm_b), EPS_ABS)

        temp = a - b
        norm = 0D0
        do i=1,6
            do j=1,6
                norm = norm + abs(temp%vals(i,j))
            end do
        end do

        if (norm/norm_max .gt. EPS) res=.false.
        if (norm/norm_max .le. EPS) res=.true.
    end function isequal_3D4O2sym

    pure function norm_3D4O2sym(a) result(norm)
        !! Uses the L1 norm of the difference of the 6x6 Voigt matrices with relative
        !! norm(A) = sum(|A_IJ|) for I,J=1..6
        implicit none
        class(ten_3D4O2sym), intent(in) :: a
        real(real64) :: norm
        integer :: i, j

        norm = 0D0
        do i=1,6
            do j=1,6
                norm = norm + abs(a%vals(i,j))
            end do
        end do
    end function norm_3D4O2sym

    pure function norm_3D4O2sym(a) result(norm)
        !! Uses the L1 norm of the difference of the 6x6 Voigt matrices with relative
        !! norm(A) = sum(|A_IJ|) for I,J=1..6
        implicit none
        class(ten_3D4O2sym), intent(in) :: a
        real(real64) :: norm
        integer :: i, j

        norm = 0D0
        do i=1,6
            do j=1,6
                norm = norm + abs(a%vals(i,j))
            end do
        end do
    end function norm_3D4O2sym

    pure function sum_3D4O2sym(a, b) result(res)
        implicit none
        class(ten_3D4O2sym), intent(in) :: a, b
        type(ten_3D4O2sym) :: res
        res%vals = a%vals + b%vals
    end function sum_3D4O2sym

    pure function sub_3D4O2sym(a, b) result(res)
        implicit none
        class(ten_3D4O2sym), intent(in) :: a, b
        type(ten_3D4O2sym) :: res
        res%vals = a%vals - b%vals
    end function sub_3D4O2sym

    pure function subU_3D4O2sym(a) result(res)
        implicit none
        class(ten_3D4O2sym), intent(in) :: a
        type(ten_3D4O2sym) :: res
        res%vals = -a%vals
    end function subU_3D4O2sym

    pure function mul_real64_3D4O2sym(a, b) result(res)
        implicit none
        real(real64), intent(in) :: a
        class(ten_3D4O2sym), intent(in) :: b
        type(ten_3D4O2sym) :: res
        res%vals = a * b%vals
    end function mul_real64_3D4O2sym

    pure function mul_3D4O2sym_real64(b, a) result(res)
        implicit none
        real(real64), intent(in) :: a
        class(ten_3D4O2sym), intent(in) :: b
        type(ten_3D4O2sym) :: res
        res%vals = a * b%vals
    end function mul_3D4O2sym_real64

    pure function div_3D4O2sym_real64(b, a) result(res)
        implicit none
        real(real64), intent(in) :: a
        class(ten_3D4O2sym), intent(in) :: b
        type(ten_3D4O2sym) :: res
        res%vals = b%vals/a
    end function div_3D4O2sym_real64

    pure function inv_3D4O2sym(a) result(res)
        !! Computes the inverse of a ten_3D4O2sym tensor.
        !! This calculates the inverse of the 6x6 Voigt matrix representation.
        !! Note: Appropriate scaling factors (1/2, 1/4) are applied to the result
        !! to ensure consistency with tensor operations like the double dot product
        !! when using Voigt notation where shear strains/stresses might be represented
        !! differently (e.g., tensorial vs engineering).
        !! Requires an external `M66INV` routine for 6x6 matrix inversion
        use, intrinsic :: iso_fortran_env
        use inverses_mat
        implicit none
        class(ten_3D4O2sym), intent(in) :: a
        type(ten_3D4O2sym) :: res
        real(real64) :: mat_a(6,6), mat_b(6,6)
        logical :: ok
        mat_a = a%vals

        call M66INV(mat_a, mat_b, ok)
        ! TODO: Add error handling if 'ok' is .false.

        ! Apply scaling factors for consistency with Voigt operations
        ! This assumes the inverse is needed for operations where Voigt notation
        ! implicitly includes factors of 2 for shear terms in contractions.
        mat_b(1:3, 4:6) = mat_b(1:3, 4:6)/2D0
        mat_b(4:6, 1:3) = mat_b(4:6, 1:3)/2D0
        mat_b(4:6, 4:6) = mat_b(4:6, 4:6)/4D0
        
        call res%init(mat_b)

    end function inv_3D4O2sym

end module mod_ten_3D4O2sym