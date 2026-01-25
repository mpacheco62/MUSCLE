module mod_ten_3D2Osym
    !! Module mod_ten_3D2Osym
    !! =======================
    !!
    !! Defines the type for symmetric 3D second-order tensors and associated operations.
    !!
    !! This module provides the derived type `ten_3D2Osym` to represent a symmetric
    !! second-order tensor in three dimensions (like stress or strain tensors).
    !! The tensor is stored internally using Voigt notation with 6 components
    !! in the order (11, 22, 33, 12, 23, 13).
    !!
    !! The module overloads standard arithmetic operators (+, -, *, /), a custom
    !! equality comparison operator (.approx.), the deviatoric operator (.dev.),
    !! and the double dot product operator (.ddot.) for this tensor type.
    !! It also provides methods for initialization and accessing individual components.
    !!
    !! Public Entities
    !! ---------------
    !!
    !! ### Derived Type:
    !!
    !! - `ten_3D2Osym`: Represents a symmetric 3D second-order tensor.
    !!     - Component: `vals(6) :: real(real64)` - Stores the 6 Voigt components
    !!       (xx, yy, zz, xy, yz, xz).
    !!     - Generic Procedure: `init` - Initializes the tensor either from a
    !!       6-element array or from individual xx, yy, zz, xy, yz, xz components.
    !!     - Procedures: `xx`, `yy`, `zz`, `xy`, `yz`, `xz` - Accessor functions
    !!       for individual tensor components.
    !!
    !! ### Operators:
    !!
    !! - `.approx.`: Compares two `ten_3D2Osym` tensors for approximate equality.
    !! - `+`: Adds two `ten_3D2Osym` tensors.
    !! - `-`: Subtracts two `ten_3D2Osym` tensors (binary) or computes the unary negation.
    !! - `*`: Multiplies a `ten_3D2Osym` tensor by a `real(real64)` scalar (or vice-versa).
    !! - `/`: Divides a `ten_3D2Osym` tensor by a `real(real64)` scalar.
    !! - `.dev.`: Computes the deviatoric part of a `ten_3D2Osym` tensor.
    !! - `.ddot.`: Computes the double dot product (scalar result) of two `ten_3D2Osym` tensors.
    !!
    !! ### Assignment:
    !!
    !! - `=`: Allows assigning a single `real(real64)` scalar value to all components
    !!        of a `ten_3D2Osym` tensor.
    !!
    !! Usage
    !! -----
    !!
    !! ```fortran
    !! program example_ten_3d2osym_usage
    !!   use mod_ten_3D2Osym
    !!   use iso_fortran_env, only: real64
    !!   implicit none
    !!
    !!   type(ten_3D2Osym) :: stress, strain, stress_dev
    !!   real(real64) :: trace_stress, dot_product
    !!   logical :: are_equal
    !!
    !!   ! Initialize using individual components (Voigt: 11, 22, 33, 12, 23, 13)
    !!   call stress%init(100.0D0, 50.0D0, 20.0D0, 10.0D0, 5.0D0, -2.0D0)
    !!
    !!   ! Initialize strain to zero using scalar assignment
    !!   strain = 0.0D0
    !!   ! Set some strain components
    !!   strain%vals(1) = 0.001 ! e_xx
    !!   strain%vals(4) = 0.002 ! e_xy
    !!
    !!   ! Calculate deviatoric stress
    !!   stress_dev = .dev. stress
    !!
    !!   ! Calculate double dot product
    !!   dot_product = stress .ddot. strain
    !!
    !!   ! Access a component
    !!   print *, "Stress XX component:", stress%xx() ! Note the parentheses for the function call
    !!   print *, "Stress XY component:", stress%xy()
    !!   print *, "Deviatoric Stress ZZ:", stress_dev%zz()
    !!   print *, "Stress : Strain =", dot_product
    !!
    !!   ! Comparison
    !!   are_equal = (stress .approx. stress_dev)
    !!   print *, "Is stress equal to its deviatoric part?", are_equal
    !!
    !! end program example_ten_3d2osym_usage
    !! ```
    !!
    !! For more information see [[tensors_types]]

    use, intrinsic :: iso_fortran_env
    implicit none
    private

    type, public :: ten_3D2Osym
        !! Symmetric 3D Second-Order Tensor (Voigt Notation)
        !! ==================================================
        !!
        !! Represents a symmetric second-order tensor in three dimensions, such as
        !! stress (\(\sigma_{ij}\)) or strain (\(\epsilon_{ij}\)), where \(\sigma_{ij} = \sigma_{ji}\).
        !! Due to symmetry, only 6 independent components are needed.
        !!
        !! Storage:
        !! --------
        !! The tensor components are stored internally in a 1D array `vals` of size 6
        !! using Voigt notation with the following mapping:
        !! - `vals(1)`: Component (1,1) or xx
        !! - `vals(2)`: Component (2,2) or yy
        !! - `vals(3)`: Component (3,3) or zz
        !! - `vals(4)`: Component (1,2) or xy (Note: This is the tensorial shear component, not the engineering one)
        !! - `vals(5)`: Component (2,3) or yz
        !! - `vals(6)`: Component (1,3) or xz
        !!
        !! Access:
        !! -------
        !! Components can be accessed directly via the `vals` array or more conveniently
        !! using the type-bound procedures `xx()`, `yy()`, `zz()`, `xy()`, `yz()`, `xz()`.
        !!
        !! Initialization:
        !! ---------------
        !! Use the generic `init` procedure to initialize either from a 6-element `real(real64)`
        !! array (following the Voigt order above) or by providing the 6 components
        !! individually (xx, yy, zz, xy, yz, xz).
        !!
        !! For more information see [[tensors_types]]
        real(real64), dimension(6) :: vals
            !! Stores the 6 independent components in Voigt notation: (xx, yy, zz, xy, yz, xz).
        contains
            generic, public :: init => init_ten_3D2Osym, init2_ten_3D2Osym
                !! Generic interface for initialization.
            procedure, private :: init_ten_3D2Osym, init2_ten_3D2Osym
            procedure, public :: xx 
                !! Accessor for the xx (1,1) component.
            procedure, public :: yy 
                !! Accessor for the yy (2,2) component.
            procedure, public :: zz
                !! Accessor for the zz (3,3) component.
            procedure, public :: xy
                !! Accessor for the xy (1,2) component.
            procedure, public :: yz
                !! Accessor for the yz (2,3) component.
            procedure, public :: xz
                !! Accessor for the xz (1,3) component.
            procedure, public :: square => square_3D2Osym
                !! Computes the square of the tensor.
            procedure, public :: norm => norm_3D2Osym
                !! Computes the norm of the tensor.
    end type ten_3D2Osym

    public :: operator(.approx.)
    interface operator (.approx.)
        module procedure approx_3D2Osym
    end interface

    public :: operator(+)
    interface operator (+)
        module procedure sum_3D2Osym
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_3D2Osym
        module procedure subU_3D2Osym
    end interface

    public :: operator(*)
    interface operator (*)
        module procedure mul_real64_3D2Osym
        module procedure mul_3D2Osym_real64
    end interface

    public :: operator( / )
    interface operator ( / )
        module procedure div_3D2Osym_real64
    end interface

    public :: operator(.dev.)
    interface operator (.dev.)
        module procedure dev_3D2Osym
    end interface

    public :: operator(.ddot.)
    interface operator (.ddot.)
        module procedure ddot_3D2Osym_3D2Osym
    end interface

    ! interface operator (.tdot.)
    !     module procedure tdot_3D2Osym_3D2Osym
    ! end interface

    public :: assignment (=)
    interface assignment (=)
        module procedure ten_3D2Osym_real64_assign
    end interface

    public :: write(formatted)
    interface write(formatted)
        module procedure print_ten_3D2Osym
    end interface

contains

    pure subroutine ten_3D2Osym_real64_assign(a, b)
        implicit none
        type(ten_3D2Osym), intent(out) :: a
        real(real64), intent(in) :: b
        a%vals = b
    end subroutine

    pure subroutine init_ten_3D2Osym(self, vals)
        !! Initializes a ten_3D2Osym tensor from a 6-element array (Voigt order).
        !!
        !! voigt notation used: 11, 22, 33, 12, 23, 13
        implicit none
        class(ten_3D2Osym), intent(inout) :: self
        real(real64), intent(in) :: vals(6)
        self%vals = vals
    end subroutine

    pure subroutine init2_ten_3D2Osym(self, xx, yy, zz, xy, yz, xz)
        !! Initializes a ten_3D2Osym tensor from its 6 individual components.
        implicit none
        class(ten_3D2Osym), intent(inout) :: self
        real(real64), intent(in) :: xx, yy, zz, xy, yz, xz
        self%vals = (/xx, yy, zz, xy, yz, xz/)
    end subroutine

    pure function norm_3D2Osym(a) result(res)
        !! Computes the norm of a ten_3D2Osym tensor.
        !! Uses a modified L1 norm (shear components weighted by 2) with relative
        !! and absolute tolerances (EPS, EPS_ABS).
        !! norm(a) = |a_11| + |a_22| + |a_33| + 2|a_12| + 2|a_23| + 2|a_13|
        implicit none
        class(ten_3D2Osym), intent(in) :: a
        real(real64) :: res

        res =   abs(a%vals(1)) + abs(a%vals(2)) + abs(a%vals(3)) &
              + 2*abs(a%vals(4)) + 2*abs(a%vals(5)) + 2*abs(a%vals(6))
    end function norm_3D2Osym

    pure function approx_3D2Osym(a, b) result(res)
        !! `.approx.` Compares two ten_3D2Osym tensors for approximate equality.
        !! Uses a modified L1 norm (shear components weighted by 2) with relative
        !! and absolute tolerances (EPS, EPS_ABS).
        !! norm(a) = |a_11| + |a_22| + |a_33| + 2|a_12| + 2|a_23| + 2|a_13|
        !! Condition: norm(a-b) / max(norm(a), norm(b), EPS_ABS) <= EPS
        implicit none
        class(ten_3D2Osym), intent(in) :: a, b
        logical :: res

        real(real64), parameter :: EPS=1e-7, EPS_ABS=1e-30
        real(real64) :: norm_a, norm_b, norm_max, norm

        norm_a =   abs(a%vals(1)) + abs(a%vals(2)) + abs(a%vals(3)) &
                 + 2*abs(a%vals(4)) + 2*abs(a%vals(5)) + 2*abs(a%vals(6))
        norm_b =   abs(b%vals(1)) + abs(b%vals(2)) + abs(b%vals(3)) &
                 + 2*abs(b%vals(4)) + 2*abs(b%vals(5)) + 2*abs(b%vals(6))

        norm_max = max(max(norm_a, norm_b), EPS_ABS)

        norm =     abs(a%vals(1)-b%vals(1)) +   abs(a%vals(2)-b%vals(2)) +   abs(a%vals(3)-b%vals(3)) &
               + 2*abs(a%vals(4)-b%vals(4)) + 2*abs(a%vals(5)-b%vals(5)) + 2*abs(a%vals(6)-b%vals(6))

        
        if (norm/norm_max .gt. EPS) res=.false.
        if (norm/norm_max .le. EPS) res=.true.
        
    end function approx_3D2Osym

    pure function sum_3D2Osym(a, b) result(res)
        implicit none
        class(ten_3D2Osym), intent(in) :: a, b
        type(ten_3D2Osym) :: res
        res%vals = a%vals + b%vals
    end function sum_3D2Osym

    pure function sub_3D2Osym(a, b) result(res)
        implicit none
        class(ten_3D2Osym), intent(in) :: a, b
        type(ten_3D2Osym) :: res
        res%vals = a%vals - b%vals
    end function sub_3D2Osym

    pure function subU_3D2Osym(a) result(res)
        implicit none
        class(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals = -a%vals
    end function subU_3D2Osym

    pure function mul_real64_3D2Osym(a, b) result(res)
        implicit none
        real(real64), intent(in) :: a
        class(ten_3D2Osym), intent(in) :: b
        type(ten_3D2Osym) :: res
        res%vals = a * b%vals
    end function mul_real64_3D2Osym

    pure function mul_3D2Osym_real64(a, b) result(res)
        implicit none
        real(real64), intent(in) :: b
        class(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals =  a%vals * b
    end function mul_3D2Osym_real64

    pure function div_3D2Osym_real64(a, b) result(res)
        implicit none
        real(real64), intent(in) :: b
        class(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals = a%vals/b
    end function div_3D2Osym_real64

    pure function ddot_3D2Osym_3D2Osym(a, b) result(res)
        implicit none
        class(ten_3D2Osym), intent(in) :: a, b
        real(real64) :: res
        res =   a%vals(1)*b%vals(1) &
                             + a%vals(2)*b%vals(2) &
                             + a%vals(3)*b%vals(3) &
                             + 2*a%vals(4)*b%vals(4) &
                             + 2*a%vals(5)*b%vals(5) &
                             + 2*a%vals(6)*b%vals(6)
    end function ddot_3D2Osym_3D2Osym

    pure function dev_3D2Osym(a) result(res)
        implicit none
        class(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym) :: res
        real(real64) :: hydro
        hydro = (a%vals(1) + a%vals(2) + a%vals(3))/3D0
        res%vals(1:3) = a%vals(1:3) - hydro
        res%vals(4:6) = a%vals(4:6)
    end function dev_3D2Osym

    pure function square_3D2Osym(a) result(res)
        implicit none
        class(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals(1) = a%vals(1)**2 + a%vals(4)**2 + a%vals(6)**2
        res%vals(2) = a%vals(2)**2 + a%vals(4)**2 + a%vals(5)**2
        res%vals(3) = a%vals(3)**2 + a%vals(5)**2 + a%vals(6)**2
        res%vals(4) = a%vals(4)*(a%vals(1) + a%vals(2)) + a%vals(5)*a%vals(6)
        res%vals(5) = a%vals(5)*(a%vals(2) + a%vals(3)) + a%vals(4)*a%vals(6)
        res%vals(6) = a%vals(6)*(a%vals(1) + a%vals(3)) + a%vals(4)*a%vals(5)
    end function square_3D2Osym
    
    pure function xx(a) result(res)
        !! Accessor function for the xx (11) component (vals(1)).
        implicit none
        class(ten_3D2Osym), intent(in) :: a
        real(real64) :: res
        res = a%vals(1)
    end function xx
    
    pure function yy(a) result(res)
        !! Accessor function for the yy (22) component (vals(2)).
        implicit none
        class(ten_3D2Osym), intent(in) :: a
        real(real64) :: res
        res = a%vals(2)
    end function yy

    pure function zz(a) result(res)
        !! Accessor function for the zz (33) component (vals(3)).
        implicit none
        class(ten_3D2Osym), intent(in) :: a
        real(real64) :: res
        res = a%vals(3)
    end function zz

    pure function xy(a) result(res)
        !! Accessor function for the xy (12) component (vals(4)).
        implicit none
        class(ten_3D2Osym), intent(in) :: a
        real(real64) :: res
        res = a%vals(4)
    end function xy

    pure function yz(a) result(res)
        !! Accessor function for the yz (23) component (vals(5)).
        implicit none
        class(ten_3D2Osym), intent(in) :: a
        real(real64) :: res
        res = a%vals(5)
    end function yz

    pure function xz(a) result(res)
        !! Accessor function for the xz (13) component (vals(6)).
        implicit none
        class(ten_3D2Osym), intent(in) :: a
        real(real64) :: res
        res = a%vals(6)
    end function xz


    subroutine print_ten_3D2Osym(dtv, unit, iotype, v_list, iostat, iomsg)
        !! Custom I/O formatting for ten_3D2Osym.
        class(ten_3D2Osym), intent(in) :: dtv
        integer, intent(in)            :: unit
        character(len=*), intent(in)   :: iotype
        integer, intent(in)            :: v_list(:)
        integer, intent(out)           :: iostat
        character(len=*), intent(inout) :: iomsg
        
        character(len=100) :: fmt_string
        integer :: w, d

        w = 10 
        d = 4  
        iostat = 0
        
        if (size(v_list) >= 1) w = v_list(1)
        if (size(v_list) >= 2) d = v_list(2)

        if (iotype == "DTLIST" .or. iotype == "LIST") then
            ! Modo lista
            write(fmt_string, "('(A, 6(F', I0, '.', I0, ', 1X), A)')") w, d
            write(unit, fmt_string, iostat=iostat) &
                "[", dtv%xx(), dtv%yy(), dtv%zz(), dtv%xy(), dtv%xz(), dtv%yz(), "]"
        else
            ! MODO MATRIZ: Un solo formato con saltos de línea incorporados
            ! Usamos el descriptor '/' para obligar el salto de línea entre filas.
            ! El primer '/' asegura que empiece en una línea nueva.
            write(fmt_string, "('(/, 3(F', I0, '.', I0, ', 2X), /, 3(F', I0, '.', I0, ', 2X), /, 3(F', I0, '.', I0, ', 2X))')") &
                w, d, w, d, w, d
            
            write(unit, fmt_string, iostat=iostat) &
                dtv%xx(), dtv%xy(), dtv%xz(), & ! Fila 1
                dtv%xy(), dtv%yy(), dtv%yz(), & ! Fila 2
                dtv%xz(), dtv%yz(), dtv%zz()    ! Fila 3
        end if

        if (iostat /= 0) then
            iomsg = "Error in print_ten_3D2Osym: Failed to write to the specified unit."
        end if
    end subroutine print_ten_3D2Osym
end module mod_ten_3D2Osym