module mod_ten_3D2O
    !! Module mod_ten_3D2O
    !! ====================
    !!
    !! Defines the type for general (non-symmetric) 3D second-order tensors and associated operations.
    !!
    !! This module provides the derived type `ten_3D2O` to represent a general
    !! second-order tensor in three dimensions. Unlike `ten_3D2Osym`, no symmetry
    !! is assumed, and all 9 components are stored explicitly.
    !! The tensor is stored internally in a 1D array representing the components
    !! in column-major order: (11, 21, 31, 12, 22, 32, 13, 23, 33).
    !!
    !! The module overloads standard arithmetic operators (+, -, *, /), a custom
    !! equality comparison operator (.approx.), the deviatoric operator (.dev.),
    !! and the double dot product operator (.ddot.) for this tensor type.
    !! It also provides methods for initialization and assignment from a scalar.
    !!
    !! Public Entities
    !! ---------------
    !!
    !! ### Derived Type:
    !!
    !! - `ten_3D2O`: Represents a general 3D second-order tensor.
    !!     - Component: `vals(9) :: real(real64)` - Stores the 9 components in
    !!       column-major order (11, 21, 31, 12, 22, 32, 13, 23, 33).
    !!     - Generic Procedure: `init` - Initializes the tensor either from a
    !!       9-element array or from individual xx, xy, xz, yx, yy, yz, zx, zy, zz components.
    !!
    !! ### Operators:
    !!
    !! - `.approx.`: Compares two `ten_3D2O` tensors for approximate equality using an L1 norm.
    !! - `+`: Adds two `ten_3D2O` tensors.
    !! - `-`: Subtracts two `ten_3D2O` tensors (binary) or computes the unary negation.
    !! - `*`: Multiplies a `ten_3D2O` tensor by a `real(real64)` scalar (or vice-versa).
    !! - `/`: Divides a `ten_3D2O` tensor by a `real(real64)` scalar.
    !! - `.dev.`: Computes the deviatoric part of a `ten_3D2O` tensor.
    !! - `.ddot.`: Computes the double dot product (Frobenius inner product, scalar result) of two `ten_3D2O` tensors.
    !!
    !! ### Assignment:
    !!
    !! - `=`: Allows assigning a single `real(real64)` scalar value to all components
    !!        of a `ten_3D2O` tensor.
    !!
    !! Usage
    !! -----
    !!
    !! ```fortran
    !! program example_ten_3d2o_usage
    !!   use mod_ten_3D2O
    !!   use iso_fortran_env, only: real64
    !!   implicit none
    !!
    !!   type(ten_3D2O) :: tensor_a, tensor_b, tensor_c, tensor_dev
    !!   real(real64) :: dot_product
    !!   logical :: are_equal
    !!
    !!   ! Initialize using individual components (xx, xy, xz, yx, yy, yz, zx, zy, zz)
    !!   call tensor_a%init(1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0)
    !!
    !!   ! Initialize tensor_b to zero using scalar assignment
    !!   tensor_b = 0.0D0
    !!   ! Set some components (using column-major indices)
    !!   tensor_b%vals(1) = 1.0 ! (1,1)
    !!   tensor_b%vals(5) = 1.0 ! (2,2)
    !!   tensor_b%vals(9) = 1.0 ! (3,3)
    !!
    !!   ! Operations
    !!   tensor_c = tensor_a + tensor_b * 2.0D0
    !!   tensor_dev = .dev. tensor_a
    !!   dot_product = tensor_a .ddot. tensor_b
    !!
    !!   ! Comparison
    !!   are_equal = (tensor_b .approx. (tensor_a - tensor_c + tensor_b * 2.0D0))
    !!
    !!   print *, "Tensor A(1,1):", tensor_a%vals(1)
    !!   print *, "Tensor A(2,1):", tensor_a%vals(2)
    !!   print *, "Tensor A(1,2):", tensor_a%vals(4)
    !!   print *, "Deviatoric A(1,1):", tensor_dev%vals(1)
    !!   print *, "A : B =", dot_product
    !!   print *, "Are equal?", are_equal
    !!
    !! end program example_ten_3d2o_usage
    !! ```
    !!
    !! For more information see [[tensors_types]]

    use, intrinsic :: iso_fortran_env
    implicit none
    private

    type, public :: ten_3D2O
        !! General 3D Second-Order Tensor (Column-Major Storage)
        !! ======================================================
        !!
        !! Represents a general (potentially non-symmetric) second-order tensor in
        !! three dimensions, \(A_{ij}\).
        !!
        !! Storage:
        !! --------
        !! The tensor components are stored internally in a 1D array `vals` of size 9
        !! using **column-major** ordering, consistent with Fortran's default array storage:
        !!
        !! - `vals(1)`: Component (1,1) or xx
        !! - `vals(2)`: Component (2,1) or yx
        !! - `vals(3)`: Component (3,1) or zx
        !! - `vals(4)`: Component (1,2) or xy
        !! - `vals(5)`: Component (2,2) or yy
        !! - `vals(6)`: Component (3,2) or zy
        !! - `vals(7)`: Component (1,3) or xz
        !! - `vals(8)`: Component (2,3) or yz
        !! - `vals(9)`: Component (3,3) or zz
        !!
        !! Access:
        !! -------
        !! Components are typically accessed directly via the `vals` array using the
        !! appropriate index based on the column-major storage.
        !!
        !! Initialization:
        !! ---------------
        !! Use the generic `init` procedure to initialize either from a 9-element `real(real64)`
        !! array (assuming column-major order) or by providing the 9 components
        !! individually in the order: (xx, xy, xz, yx, yy, yz, zx, zy, zz). Note that the
        !! `init2` procedure internally rearranges these into the column-major `vals` array.
        !!
        !! For more information see [[tensors_types]]

        real(real64), dimension(9) :: vals
            !! Stores the 9 components in column-major order: (11, 21, 31, 12, 22, 32, 13, 23, 33).
    contains
        generic, public :: init => init_ten_3D2O, init2_ten_3D2O
            !! Generic interface for initialization.
        procedure, private :: init_ten_3D2O, init2_ten_3D2O
        procedure, public :: norm => norm_3D2O
            !! Computes the norm of the tensor.
        procedure, public :: is_approx => is_approx_3D2O
            !! Compares two tensors for approximate equality.
    end type ten_3D2O

    public :: operator(.approx.)
    interface operator (.approx.)
        module procedure approx_3D2O
    end interface

    public :: operator(+)
    interface operator (+)
        module procedure sum_3D2O
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_3D2O
        module procedure subU_3D2O
    end interface

    public :: operator(*)
    interface operator (*)
        module procedure mul_real64_3D2O
        module procedure mul_3D2O_real64
    end interface

    public :: operator( / )
    interface operator ( / )
        module procedure div_3D2O_real64
    end interface

    public :: operator(.dev.)
    interface operator (.dev.)
        module procedure dev_3D2O
    end interface

    public :: operator(.ddot.)
    interface operator (.ddot.)
        module procedure ddot_3D2O_3D2O
    end interface

    public :: assignment (=)
    interface assignment (=)
        module procedure ten_3D2O_real64_assign
    end interface

    public :: write(formatted)
    interface write(formatted)
        module procedure print_ten_3D2O
    end interface
contains

    pure subroutine ten_3D2O_real64_assign(a, b)
        implicit none
        type(ten_3D2O), intent(out) :: a
        real(real64), intent(in) :: b
        a%vals = b
    end subroutine ten_3D2O_real64_assign

    subroutine init_ten_3D2O(self, vals)
        !! Initializes a ten_3D2O tensor from a 9-element array (assumed column-major).
        implicit none
        class(ten_3D2O), intent(inout) :: self
        real(real64), intent(in) :: vals(9)
        self%vals = vals
    end subroutine init_ten_3D2O

    subroutine init2_ten_3D2O(self, xx, xy, xz, yx, yy, yz, zx, zy, zz)
        !! Initializes a ten_3D2O tensor from its 9 individual components.
        !! Input order is (xx, xy, xz, yx, yy, yz, zx, zy, zz).
        !! Internal storage is column-major: (xx, yx, zx, xy, yy, zy, xz, yz, zz).
        implicit none
        class(ten_3D2O), intent(inout) :: self
        real(real64), intent(in) :: xx, xy, xz, yx, yy, yz, zx, zy, zz
        self%vals = (/xx, yx, zx, xy, yy, zy, xz, yz, zz/)
    end subroutine init2_ten_3D2O

    pure function norm_3D2O(a) result(res)
        !! Computes the norm of a ten_3D2O tensor.
        implicit none
        class(ten_3D2O), intent(in) :: a
        real(real64) :: res

        res =   abs(a%vals(1)) + abs(a%vals(2)) + abs(a%vals(3)) &
              + abs(a%vals(4)) + abs(a%vals(5)) + abs(a%vals(6)) &
              + abs(a%vals(7)) + abs(a%vals(8)) + abs(a%vals(9))
    end function norm_3D2O

    pure function is_approx_3D2O(a, b, tol) result(res)
        implicit none
        class(ten_3D2O), intent(in) :: a, b
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
    end function is_approx_3D2O

    pure function approx_3D2O(a, b) result(res)
        !! `.approx.` Compares two ten_3D2O tensors for approximate equality.
        !! Uses the standard L1 norm with relative and absolute tolerances (EPS, EPS_ABS).
        !! norm(a) = sum(|a_ij|) for all i,j
        !! Condition: norm(a-b) / max(norm(a), norm(b), EPS_ABS) <= EPS
        implicit none
        class(ten_3D2O), intent(in) :: a, b
        logical :: res

        res = a%is_approx(b)
        
    end function approx_3D2O

    pure function sum_3D2O(a, b) result(res)
        implicit none
        class(ten_3D2O), intent(in) :: a, b
        type(ten_3D2O) :: res

        res%vals = a%vals + b%vals
    end function sum_3D2O

    pure function sub_3D2O(a, b) result(res)
        implicit none
        class(ten_3D2O), intent(in) :: a, b
        type(ten_3D2O) :: res

        res%vals = a%vals - b%vals
    end function sub_3D2O

    pure function subU_3D2O(a) result(res)
        implicit none
        class(ten_3D2O), intent(in) :: a
        type(ten_3D2O) :: res

        res%vals = -a%vals
    end function subU_3D2O

    pure function mul_real64_3D2O(a, b) result(res)
        implicit none
        real(real64), intent(in) :: a
        class(ten_3D2O), intent(in) :: b
        type(ten_3D2O) :: res

        res%vals = a * b%vals
    end function mul_real64_3D2O

    pure function mul_3D2O_real64(a, b) result(res)
        implicit none
        real(real64), intent(in) :: b
        class(ten_3D2O), intent(in) :: a
        type(ten_3D2O) :: res

        res%vals =  a%vals * b
    end function mul_3D2O_real64

    pure function div_3D2O_real64(a, b) result(res)
        implicit none
        real(real64), intent(in) :: b
        class(ten_3D2O), intent(in) :: a
        type(ten_3D2O) :: res

        res%vals = a%vals/b
    end function div_3D2O_real64

    pure function ddot_3D2O_3D2O(a, b) result(res)
        implicit none
        class(ten_3D2O), intent(in) :: a, b
        real(real64) :: res
    
        res =   a%vals(1)*b%vals(1) &
              + a%vals(2)*b%vals(2) &
              + a%vals(3)*b%vals(3) &
              + a%vals(4)*b%vals(4) &
              + a%vals(5)*b%vals(5) &
              + a%vals(6)*b%vals(6) &
              + a%vals(7)*b%vals(7) &
              + a%vals(8)*b%vals(8) &
              + a%vals(9)*b%vals(9)
    
    end function ddot_3D2O_3D2O

    pure function dev_3D2O(a) result(res)
        implicit none
        class(ten_3D2O), intent(in) :: a
        type(ten_3D2O) :: res

        real(real64) :: hydro
        hydro = (a%vals(1) + a%vals(5) + a%vals(9))/3D0
        res%vals(1) = a%vals(1) - hydro
        res%vals(5) = a%vals(5) - hydro
        res%vals(9) = a%vals(9) - hydro
        res%vals(2:4) = a%vals(2:4)
        res%vals(6:8) = a%vals(6:8)
    end function dev_3D2O


    subroutine print_ten_3D2O(dtv, unit, iotype, v_list, iostat, iomsg)
        !! Custom I/O formatting for ten_3D2Osym.
        class(ten_3D2O), intent(in) :: dtv
        integer, intent(in)            :: unit
        character(len=*), intent(in)   :: iotype
        integer, intent(in)            :: v_list(:)
        integer, intent(out)           :: iostat
        character(len=*), intent(inout) :: iomsg
        
        character(len=100) :: fmt_string
        integer :: w, d

        w = 11
        d = 4  
        iostat = 0
        
        if (size(v_list) >= 1) w = v_list(1)
        if (size(v_list) >= 2) d = v_list(2)

        if (iotype == "DTLIST" .or. iotype == "LIST") then
            ! Modo lista
            write(fmt_string, "('(A, 9(ES', I0, '.', I0, ', 1X), A)')") w, d
            write(unit, fmt_string, iostat=iostat) &
                "[", dtv%vals(1), dtv%vals(2), dtv%vals(3), dtv%vals(4), dtv%vals(5), & 
                     dtv%vals(6), dtv%vals(7), dtv%vals(8), dtv%vals(9), "]"
        else
            ! MODO MATRIZ: Un solo formato con saltos de línea incorporados
            ! Usamos el descriptor '/' para obligar el salto de línea entre filas.
            ! El primer '/' asegura que empiece en una línea nueva.
            write(fmt_string, "('(/, 3(ES', I0, '.', I0, ', 2X), /, 3(F', I0, '.', I0, ', 2X), /, 3(F', I0, '.', I0, ', 2X))')") &
                w, d, w, d, w, d
            
            write(unit, fmt_string, iostat=iostat) &
                dtv%vals(1), dtv%vals(4), dtv%vals(7), & ! Fila 1
                dtv%vals(2), dtv%vals(5), dtv%vals(8), & ! Fila 2
                dtv%vals(3), dtv%vals(6), dtv%vals(9)    ! Fila 3
        end if

        if (iostat /= 0) then
            iomsg = "Error in print_ten_3D2Osym: Failed to write to the specified unit."
        end if
    end subroutine print_ten_3D2O

    ! ! !TODO HACERLE UN TESTTTSSS!!!!!
    ! ! module procedure tdot_3D2Osym_3D2Osym
    ! !     implicit none
    ! !     !
    ! !     !  | ( 1:1111) ( 7:1122) (12:1133) (16:1112) (19:1123) (21:1113) |
    ! !     !  | ( 7:2211) ( 2:2222) ( 8:2233) (13:2212) (17:2223) (20:2213) |
    ! !     !  | (12:3311) ( 8:3322) ( 3:3333) ( 9:3312) (14:3323) (18:3313) |
    ! !     !  | (16:1211) (13:1222) ( 9:1233) ( 4:1212) (10:1223) (15:1213) |
    ! !     !  | (19:2311) (17:2322) (14:2333) (10:2312) ( 5:2323) (11:2313) |
    ! !     !  | (21:1311) (20:1322) (18:1333) (15:1312) (11:1323) ( 6:1313) |
    ! !     res%vals( 1) = a%vals(1)*b%vals(1)
    ! !     res%vals( 7) = a%vals(1)*b%vals(2)
    ! !     res%vals(12) = a%vals(1)*b%vals(3)
    ! !     res%vals(16) = a%vals(1)*b%vals(4)
    ! !     res%vals(19) = a%vals(1)*b%vals(5)
    ! !     res%vals(21) = a%vals(1)*b%vals(6)

    ! !     res%vals( 2) = a%vals(2)*b%vals(2)
    ! !     res%vals( 8) = a%vals(2)*b%vals(3)
    ! !     res%vals(13) = a%vals(2)*b%vals(4)
    ! !     res%vals(17) = a%vals(2)*b%vals(5)
    ! !     res%vals(20) = a%vals(2)*b%vals(6)

    ! !     res%vals( 3) = a%vals(3)*b%vals(3)
    ! !     res%vals( 9) = a%vals(3)*b%vals(4)
    ! !     res%vals(14) = a%vals(3)*b%vals(5)
    ! !     res%vals(18) = a%vals(3)*b%vals(6)

    ! !     res%vals( 4) = a%vals(4)*b%vals(4)
    ! !     res%vals(10) = a%vals(4)*b%vals(5)
    ! !     res%vals(15) = a%vals(4)*b%vals(6)

    ! !     res%vals( 5) = a%vals(5)*b%vals(5)
    ! !     res%vals(11) = a%vals(5)*b%vals(6)
        
    ! !     res%vals( 6) = a%vals(6)*b%vals(6)
    ! ! end procedure tdot_3D2Osym_3D2Osym

end module mod_ten_3D2O