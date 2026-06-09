module muscle_tensor_3d4o3sym
    !! Module muscle_tensor_3d4o3sym
    !! =======================
    !!
    !! Defines the type for fully symmetric 3D fourth-order tensors and associated operations.
    !!
    !! This module provides the derived type `ten_3D4O3sym` to represent a fourth-order
    !! tensor in three dimensions possessing both major and minor symmetries:
    !! \[ C_{ijkl} = C_{jikl} = C_{ijlk} = C_{klij} \]
    !! Such tensors are common in linear elasticity (the elasticity tensor).
    !!
    !! The tensor is stored internally using a compressed format with 21 components,
    !! corresponding to the independent elements of the symmetric 6x6 Voigt matrix representation.
    !! The storage order follows a specific convention (see type description).
    !!
    !! The module overloads standard arithmetic operators (+, -, *, /), a custom
    !! equality comparison operator (.approx.), an inverse operator (.inv.),
    !! and provides methods for initialization.
    !!
    !! Public Entities
    !! ---------------
    !!
    !! ### Derived Type:
    !!
    !! - `ten_3D4O3sym`: Represents a 3D fourth-order tensor with major and minor symmetries.
    !!     - Component: `vals(21) :: real(real64)` - Stores the 21 independent components.
    !!     - Generic Procedure: `init` - Initializes the tensor either from a
    !!       21-element array or from 21 individual components.
    !!
    !! ### Operators:
    !!
    !! - `.approx.`: Compares two `ten_3D4O3sym` tensors for approximate equality using a modified L1 norm.
    !! - `.inv.`: Computes the inverse of the tensor based on its 6x6 symmetric Voigt matrix representation.
    !! - `+`: Adds two `ten_3D4O3sym` tensors.
    !! - `-`: Subtracts two `ten_3D4O3sym` tensors (binary) or computes the unary negation.
    !! - `*`: Multiplies a `ten_3D4O3sym` tensor by a `real(real64)` scalar (or vice-versa).
    !! - `/`: Divides a `ten_3D4O3sym` tensor by a `real(real64)` scalar.
    !!
    !! Usage
    !! -----
    !!
    !! ```fortran
    !! program example_ten_3d4o3sym_usage
    !!   use muscle_tensors ! Includes muscle_tensor_3d4o3sym and others
    !!   use iso_fortran_env, only: real64
    !!   implicit none
    !!
    !!   type(ten_3D4O3sym) :: C_iso, C_inv
    !!   real(real64) :: lambda, mu
    !!   real(real64) :: C_vals(21)
    !!   logical :: are_equal
    !!
    !!   ! Initialize isotropic elasticity tensor using individual components
    !!   lambda = 120.0D0 ! Lame's first parameter
    !!   mu = 80.0D0      ! Shear modulus (Lame's second parameter)
    !!
    !!   ! Components based on init2 order:
    !!   ! xxxx, yyyy, zzzz, xyxy, yzyz, xzxz, xxyy, yyzz, zzxy, xyyz, yzxz,
    !!   ! xxzz, yyxy, zzyz, xyxz, xxxy, yyyz, zzxz, xxyz, yyxz, xxxz
    !!   call C_iso%init( &
    !!       lambda + 2*mu, lambda + 2*mu, lambda + 2*mu, & ! xxxx, yyyy, zzzz
    !!       mu, mu, mu,                                   & ! xyxy, yzyz, xzxz
    !!       lambda, lambda,                               & ! xxyy, yyzz
    !!       0.0, 0.0, 0.0,                                & ! zzxy, xyyz, yzxz
    !!       lambda,                                       & ! xxzz
    !!       0.0, 0.0, 0.0,                                & ! yyxy, zzyz, xyxz
    !!       0.0, 0.0, 0.0,                                & ! xxxy, yyyz, zzxz
    !!       0.0, 0.0,                                     & ! xxyz, yyxz
    !!       0.0                                          ) ! xxxz
    !!
    !!   ! Calculate the inverse (compliance tensor S = C^-1)
    !!   C_inv = .inv. C_iso
    !!
    !!   ! Check a component of the inverse
    !!   ! S_1111 = C_inv%vals(1) should be (lambda+mu)/(mu*(3*lambda+2*mu))
    !!   print *, "C_iso(1) (C_1111):", C_iso%vals(1)
    !!   print *, "C_inv(1) (S_1111):", C_inv%vals(1)
    !!   print *, "Expected S_1111:", (lambda+mu)/(mu*(3*lambda+2*mu))
    !!
    !!   ! Comparison
    !!   are_equal = (C_iso .approx. C_iso)
    !!   print *, "Is C_iso equal to itself?", are_equal
    !!
    !! end program example_ten_3d4o3sym_usage
    !! ```
    !!
    !! For more information see [[muscle_tensors]]

    use, intrinsic :: iso_fortran_env
    implicit none
    private
 
    type, public :: ten_3D4O3sym
        !! Fully Symmetric 3D Fourth-Order Tensor (21-Component Storage)
        !! ==============================================================
        !!
        !! Represents a fourth-order tensor \(C_{ijkl}\) in three dimensions that possesses
        !! both minor and major symmetries:
        !! \[ C_{ijkl} = C_{jikl} = C_{ijlk} = C_{klij} \]
        !! This is the standard symmetry for linear elastic constitutive tensors.
        !!
        !! Storage:
        !! --------
        !! Due to the symmetries, only 21 independent components exist. These are stored
        !! internally in a 1D array `vals` of size 21. The storage order corresponds
        !! to the upper triangle of the 6x6 symmetric Voigt matrix representation,
        !! stored column-wise:
        !!
        !! Voigt Matrix (Indices IJ):
        !! ```
        !! | 11 12 13 14 15 16 |
        !! |    22 23 24 25 26 |
        !! |       33 34 35 36 |
        !! |          44 45 46 |
        !! |             55 56 |
        !! |                66 |
        !! ```
        !! Storage order in `vals(1:21)`:
        !! (11, 22, 33, 44, 55, 66, 12, 23, 34, 45, 56, 13, 24, 35, 46, 14, 25, 36, 15, 26, 16)
        !!   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21
        !!
        !! where the Voigt index mapping is:
        !! - 1 <-> (1,1) or xx
        !! - 2 <-> (2,2) or yy
        !! - 3 <-> (3,3) or zz
        !! - 4 <-> (1,2) or (2,1) or xy
        !! - 5 <-> (2,3) or (3,2) or yz
        !! - 6 <-> (1,3) or (3,1) or xz
        !!
        !! Access:
        !! -------
        !! Components are typically accessed directly via the `vals` array using the
        !! appropriate index based on the storage order.
        !!
        !! Initialization:
        !! ---------------
        !! Use the generic `init` procedure to initialize either from a 21-element `real(real64)`
        !! array (following the storage order above) or by providing the 21 components
        !! individually (see `init2_ten_3D4O3sym` for the required input order).
        !!
        !! For more information see [[muscle_tensors]]

        real(real64), dimension(21) :: vals
            !! Stores the 21 independent components following the compressed Voigt storage order.
        contains
            generic, public :: init => init_ten_3D4O3sym, init2_ten_3D4O3sym
                !! Generic interface for initialization.
            procedure, private :: init_ten_3D4O3sym, init2_ten_3D4O3sym 
            procedure, public :: norm => norm_3D4O3sym
            procedure, public :: set => set_ten3D4O3sym
            procedure, public :: is_approx => is_approx_3D4O3sym
    end type ten_3D4O3sym

    public :: operator(.approx.)
    interface operator (.approx.)
        module procedure approx_3D4O3sym
    end interface

    public :: operator(.inv.)
    interface operator (.inv.)
        module procedure inv_3D4O3sym
    end interface

    public :: operator(+)
    interface operator (+)
        module procedure sum_3D4O3sym
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_3D4O3sym
        module procedure subU_3D4O3sym
    end interface

    public :: operator(*)
    interface operator (*)
        module procedure mul_real64_3D4O3sym
        module procedure mul_3D4O3sym_real64
    end interface

    public :: operator( / )
    interface operator ( / )
        module procedure div_3D4O3sym_real64
    end interface

    public :: assignment (=)
    interface assignment (=)
        module procedure assign_ten_3D4O3sym_real64
    end interface
contains
    
    pure subroutine init_ten_3D4O3sym(self, vals)
        !! Initializes a ten_3D4O3sym tensor from a 21-element array (compressed Voigt order).
        implicit none
        class(ten_3D4O3sym), intent(inout) :: self
        real(real64), intent(in) :: vals(21)
        self%vals = vals
    end subroutine init_ten_3D4O3sym

    pure subroutine init2_ten_3D4O3sym(self,             & 
                                              xxxx, yyyy, zzzz, &
                                              xyxy, yzyz, xzxz, &
                                              xxyy, yyzz,       &
                                              zzxy, xyyz, yzxz, &
                                              xxzz,             &
                                              yyxy, zzyz, xyxz, &
                                              xxxy, yyyz, zzxz, &
                                              xxyz, yyxz, xxxz  &
                                              )
        !! Initializes a ten_3D4O3sym tensor from its 21 individual components.
        !! Input arguments correspond to the independent Voigt matrix components C(I,J)
        !! in a specific order (see implementation and type documentation).
        !
        !  Voigt Matrix Layout (I, J):
        !  | (1,1) (1,2) (1,3) (1,4) (1,5) (1,6) |   <- xxxx, xxyy, xxzz, xxxy, xxyz, xxxz
        !  | (2,1) (2,2) (2,3) (2,4) (2,5) (2,6) |   <- yyxx, yyyy, yyzz, yyxy, yyyz, yyxz
        !  | (3,1) (3,2) (3,3) (3,4) (3,5) (3,6) |   <- zzxx, zzyy, zzzz, zzxy, zzyz, zzxz
        !  | (4,1) (4,2) (4,3) (4,4) (4,5) (4,6) |   <- xyxx, xyyy, xyzz, xyxy, xyyz, xyxz
        !  | (5,1) (5,2) (5,3) (5,4) (5,5) (5,6) |   <- yzxx, yzyy, yzzz, yzxy, yzyz, yzxz
        !  | (6,1) (6,2) (6,3) (6,4) (6,5) (6,6) |   <- xzxx, xzyy, xzzz, xzxy, xzyz, xzxz
        implicit none
        class(ten_3D4O3sym), intent(inout) :: self
        real(real64), intent(in) :: xxxx, xxyy, xxzz, xxxy, xxyz, xxxz
        real(real64), intent(in) :: yyyy, yyzz, yyxy, yyyz, yyxz
        real(real64), intent(in) :: zzzz, zzxy, zzyz, zzxz
        real(real64), intent(in) :: xyxy, xyyz, xyxz
        real(real64), intent(in) :: yzyz, yzxz
        real(real64), intent(in) :: xzxz
        self%vals = (/xxxx, yyyy, zzzz, xyxy, yzyz, xzxz, &
                      xxyy, yyzz, zzxy, xyyz, yzxz, &
                      xxzz, yyxy, zzyz, xyxz, &
                      xxxy, yyyz, zzxz, &
                      xxyz, yyxz, &
                      xxxz &
                      /)
    end subroutine

    pure function is_approx_3D4O3sym(a, b, tol) result(res)
        implicit none
        class(ten_3D4O3sym), intent(in) :: a, b
        real(real64), optional, intent(in) :: tol
        logical :: res
        real(real64), parameter :: EPS_ABS=1e-30
        real(real64) :: eps_check
        real(real64) :: tol2
        real(real64) :: max_val

        tol2 = 1E-12
        if (present(tol)) tol2 = tol
        max_val = max(maxval(abs(a%vals)), maxval(abs(b%vals)), EPS_ABS)
        eps_check = tol2 * max_val

        res = all(abs(a%vals - b%vals) .le. eps_check)
        return
    end function is_approx_3D4O3sym


    pure function approx_3D4O3sym(a, b) result(res)
        !! `.approx.` Compares two ten_3D4O3sym tensors for approximate equality.
        !! Uses a modified L1 norm based on the 21 stored components, where components
        !! corresponding to off-diagonal Voigt matrix entries are weighted by 2.
        !! norm(A) = sum(|A_diag|) + 2*sum(|A_offdiag|) based on 6x6 Voigt matrix.
        !! Condition: norm(a-b) / max(norm(a), norm(b), EPS_ABS) <= EPS
        implicit none
        type(ten_3D4O3sym), intent(in) :: a, b
        logical :: res
 
        res = a%is_approx(b)

        return
    end function approx_3D4O3sym

    pure function norm_3D4O3sym(a) result(norms)
        !! Uses a modified L1 norm based on the 21 stored components, where components
        !! norm(A) = sum(|A_diag|) + 2*sum(|A_offdiag|) based on 6x6 Voigt matrix.
        implicit none
        class(ten_3D4O3sym), intent(in) :: a
        real(real64) :: norms
        norms = sum(abs(a%vals(1:6))) + 2.0D0 * sum(abs(a%vals(7:21)))
    end function norm_3D4O3sym


    pure function sum_3D4O3sym(a, b) result(res)
        implicit none
        type(ten_3D4O3sym), intent(in) :: a, b
        type(ten_3D4O3sym) :: res
        res%vals = a%vals + b%vals
    end function sum_3D4O3sym

    pure function sub_3D4O3sym(a, b) result(res)
        implicit none
        type(ten_3D4O3sym), intent(in) :: a, b
        type(ten_3D4O3sym) :: res
        res%vals = a%vals - b%vals
    end function sub_3D4O3sym

    pure function subU_3D4O3sym(a) result(res)
        implicit none
        type(ten_3D4O3sym), intent(in) :: a
        type(ten_3D4O3sym) :: res
        res%vals = -a%vals
    end function subU_3D4O3sym

    pure function mul_real64_3D4O3sym(a, b) result(res)
        implicit none
        real(real64), intent(in) :: a
        type(ten_3D4O3sym), intent(in) :: b
        type(ten_3D4O3sym) :: res
        res%vals = a * b%vals
    end function mul_real64_3D4O3sym

    pure function mul_3D4O3sym_real64(b, a) result(res)
        implicit none
        real(real64), intent(in) :: a
        type(ten_3D4O3sym), intent(in) :: b
        type(ten_3D4O3sym) :: res
        res%vals = a * b%vals
    end function mul_3D4O3sym_real64

    pure function div_3D4O3sym_real64(b, a) result(res)
        implicit none
        real(real64), intent(in) :: a
        type(ten_3D4O3sym), intent(in) :: b
        type(ten_3D4O3sym) :: res
        res%vals = b%vals/a
    end function div_3D4O3sym_real64

    pure function inv_3D4O3sym(a) result(res)
        !! `.inv.` Computes the inverse of a ten_3D4O3sym tensor.
        !! Reconstructs the 6x6 symmetric Voigt matrix, inverts it using `M66INV`,
        !! and extracts the 21 independent components of the inverse, applying
        !! scaling factors (1/2, 1/4) to the off-diagonal Voigt components of the
        !! resulting inverse matrix before storing them. This ensures consistency
        !! for subsequent Voigt-based tensor operations.
        use, intrinsic :: iso_fortran_env
        use muscle_math_inverses
        implicit none
        type(ten_3D4O3sym), intent(in) :: a
        type(ten_3D4O3sym) :: res
        real(real64) :: mat_a(6,6), mat_b(6,6), v(21)
        logical :: ok
        v = a%vals
        mat_a = reshape((/  v(1),  v(7), v(12), v(16), v(19), v(21), &
                            v(7),  v(2),  v(8), v(13), v(17), v(20), &
                           v(12),  v(8),  v(3),  v(9), v(14), v(18), & 
                           v(16), v(13),  v(9),  v(4), v(10), v(15), &
                           v(19), v(17), v(14), v(10),  v(5), v(11), &
                           v(21), v(20), v(18), v(15), v(11),  v(6)  &
                        /), (/6,6/))

        call M66INV2(mat_a, mat_b, ok)


        ! call FINDInv(mat_a, mat_b, 6, iok)

        ! v = (/ mat_b(1,1), mat_b(2,2), mat_b(3,3), mat_b(4,4), mat_b(5,5), mat_b(6,6),  &
        !        mat_b(1,2), mat_b(2,3), mat_b(3,4), mat_b(4,5), mat_b(5,6),              &
        !        mat_b(1,3), mat_b(2,4), mat_b(3,5), mat_b(4,6),                          &
        !        mat_b(1,4), mat_b(2,5), mat_b(3,6),                                      &
        !        mat_b(1,5), mat_b(2,6),                                                  &
        !        mat_b(1,6)                                                               &
        !       /)

        ! TODO: Add error handling if 'ok' is .false. or iok indicates failure

        ! Extract the 21 components of the inverse matrix (mat_b)
        ! Apply scaling factors to off-diagonal Voigt components of the inverse
        ! This scaling is necessary for the inverse tensor stored in 21-component
        ! format to work correctly in standard Voigt tensor contractions like C:E
        ! where factors of 2 or 4 might be implicitly assumed depending on the
        ! definition of the contraction and the storage format.

        v = (/ mat_b(1,1), mat_b(2,2), mat_b(3,3), mat_b(4,4)/4D0, mat_b(5,5)/4D0, mat_b(6,6)/4D0,  &
               mat_b(1,2), mat_b(2,3), mat_b(3,4)/2D0, mat_b(4,5)/4D0, mat_b(5,6)/4D0,              &
               mat_b(1,3), mat_b(2,4)/2D0, mat_b(3,5)/2D0, mat_b(4,6)/4D0,                          &
               mat_b(1,4)/2D0, mat_b(2,5)/2D0, mat_b(3,6)/2D0,                                      &
               mat_b(1,5)/2D0, mat_b(2,6)/2D0,                                                  &
               mat_b(1,6)/2D0                                                               &
              /)


        ! mat_a = reshape((/  v(1),  v(7), v(12), v(16), v(19), v(21), &
        !                     v(7),  v(2),  v(8), v(13), v(17), v(20), &
        !                     v(12),  v(8),  v(3),  v(9), v(14), v(18), & 
        !                     v(16), v(13),  v(9),  v(4), v(10), v(15), &
        !                     v(19), v(17), v(14), v(10),  v(5), v(11), &
        !                     v(21), v(20), v(18), v(15), v(11),  v(6)  &
        !                     /), (/6,6/))
        ! call M66INV(mat_a, mat_b, ok)
        ! call FINDInv(mat_a, mat_b, 6, iok)
        
        call res%init(v)
        
            !
            !  | ( 1:1111) ( 7:1122) (12:1133) (16:1112) (19:1123) (21:1113) |
            !  | ( 7:2211) ( 2:2222) ( 8:2233) (13:2212) (17:2223) (20:2213) |
            !  | (12:3311) ( 8:3322) ( 3:3333) ( 9:3312) (14:3323) (18:3313) |
            !  | (16:1211) (13:1222) ( 9:1233) ( 4:1212) (10:1223) (15:1213) |
            !  | (19:2311) (17:2322) (14:2333) (10:2312) ( 5:2323) (11:2313) |
            !  | (21:1311) (20:1322) (18:1333) (15:1312) (11:1323) ( 6:1313) |

    end function inv_3D4O3sym

    pure subroutine assign_ten_3D4O3sym_real64(a, b)
        implicit none
        type(ten_3D4O3sym), intent(out) :: a
        real(real64), intent(in) :: b
        a%vals = b
    end subroutine assign_ten_3D4O3sym_real64

    pure subroutine set_ten3D4O3sym(a, &
                                    xxxx, yyyy, zzzz, xyxy, yzyz, xzxz, &
                                    xxyy, yyzz, zzxy, xyyz, yzxz, xxzz, &
                                    yyxy, zzyz, xyxz, xxxy, yyyz, zzxz, &
                                    xxyz, yyxz, xxxz  &
                                    )
        implicit none
        class(ten_3D4O3sym), intent(inout) :: a
        real(real64), optional, intent(in) :: xxxx, yyyy, zzzz, xyxy, yzyz, xzxz, &
                                              xxyy, yyzz, zzxy, xyyz, yzxz, xxzz, &
                                              yyxy, zzyz, xyxz, xxxy, yyyz, zzxz, &
                                              xxyz, yyxz, xxxz

        if (present(xxxx)) a%vals(1) = xxxx
        if (present(yyyy)) a%vals(2) = yyyy
        if (present(zzzz)) a%vals(3) = zzzz
        if (present(xyxy)) a%vals(4) = xyxy
        if (present(yzyz)) a%vals(5) = yzyz
        if (present(xzxz)) a%vals(6) = xzxz
        if (present(xxyy)) a%vals(7) = xxyy
        if (present(yyzz)) a%vals(8) = yyzz
        if (present(zzxy)) a%vals(9) = zzxy
        if (present(xyyz)) a%vals(10) = xyyz
        if (present(yzxz)) a%vals(11) = yzxz
        if (present(xxzz)) a%vals(12) = xxzz
        if (present(yyxy)) a%vals(13) = yyxy
        if (present(zzyz)) a%vals(14) = zzyz
        if (present(xyxz)) a%vals(15) = xyxz
        if (present(xxxy)) a%vals(16) = xxxy
        if (present(yyyz)) a%vals(17) = yyyz
        if (present(zzxz)) a%vals(18) = zzxz
        if (present(xxyz)) a%vals(19) = xxyz
        if (present(yyxz)) a%vals(20) = yyxz
        if (present(xxxz)) a%vals(21) = xxxz
    end subroutine

end module muscle_tensor_3d4o3sym