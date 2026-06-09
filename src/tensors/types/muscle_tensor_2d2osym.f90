! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_tensor_2d2osym
    !! Module muscle_tensor_2d2osym
    !! =======================
    !!
    !! Defines the type for symmetric 3D second-order tensors and associated operations.
    !!
    !! This module provides the derived type `ten_2D2Osym` to represent a symmetric
    !! second-order tensor in three dimensions (like stress or strain tensors).
    !! The tensor is stored internally using Voigt notation with 6 components
    !! in the order (11, 22, 33, 12).
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
    !! - `ten_2D2Osym`: Represents a symmetric 3D second-order tensor.
    !!     - Component: `vals(4) :: real(real64)` - Stores the 4 Voigt components
    !!       (xx, yy, zz, xy).
    !!     - Generic Procedure: `init` - Initializes the tensor either from a
    !!       4-element array or from individual xx, yy, zz, xy components.
    !!     - Procedures: `xx`, `yy`, `zz`, `xy` - Accessor functions
    !!       for individual tensor components.
    !!
    !! ### Operators:
    !!
    !! - `.approx.`: Compares two `ten_2D2Osym` tensors for approximate equality.
    !! - `+`: Adds two `ten_2D2Osym` tensors.
    !! - `-`: Subtracts two `ten_2D2Osym` tensors (binary) or computes the unary negation.
    !! - `*`: Multiplies a `ten_2D2Osym` tensor by a `real(real64)` scalar (or vice-versa).
    !! - `/`: Divides a `ten_2D2Osym` tensor by a `real(real64)` scalar.
    !! - `.dev.`: Computes the deviatoric part of a `ten_2D2Osym` tensor.
    !! - `.ddot.`: Computes the double dot product (scalar result) of two `ten_2D2Osym` tensors.
    !!
    !! ### Assignment:
    !!
    !! - `=`: Allows assigning a single `real(real64)` scalar value to all components
    !!        of a `ten_2D2Osym` tensor.
    !!
    !! Usage
    !! -----
    !!
    !! ```fortran
    !! program example_ten_2D2Osym_usage
    !!   use muscle_tensor_2d2osym
    !!   use iso_fortran_env, only: real64
    !!   implicit none
    !!
    !!   type(ten_2D2Osym) :: stress, strain, stress_dev
    !!   real(real64) :: trace_stress, dot_product
    !!   logical :: are_equal
    !!
    !!   ! Initialize using individual components (Voigt: 11, 22, 33, 12)
    !!   call stress%init(100.0D0, 50.0D0, 20.0D0, 10.0D0)
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
    !! end program example_ten_2D2Osym_usage
    !! ```
    !!
    !! For more information see [[muscle_tensors]]

    use, intrinsic :: iso_fortran_env
    implicit none
    private

    type, public :: ten_2D2Osym
        !! Symmetric 2D Second-Order Tensor (Voigt Notation)
        !! ==================================================
        !!
        !! Represents a symmetric second-order tensor in three dimensions, such as
        !! stress (\(\sigma_{ij}\)) or strain (\(\epsilon_{ij}\)), where \(\sigma_{ij} = \sigma_{ji}\).
        !! Due to symmetry, only 4 independent components are needed.
        !!
        !! Storage:
        !! --------
        !! The tensor components are stored internally in a 1D array `vals` of size 4
        !! using Voigt notation with the following mapping:
        !! - `vals(1)`: Component (1,1) or xx
        !! - `vals(2)`: Component (2,2) or yy
        !! - `vals(3)`: Component (3,3) or zz
        !! - `vals(4)`: Component (1,2) or xy (Note: This is the tensorial shear component, not the engineering one)
        !!
        !! Access:
        !! -------
        !! Components can be accessed directly via the `vals` array or more conveniently
        !! using the type-bound procedures `xx()`, `yy()`, `zz()`, `xy()`.
        !!
        !! Initialization:
        !! ---------------
        !! Use the generic `init` procedure to initialize either from a 4-element `real(real64)`
        !! array (following the Voigt order above) or by providing the 4 components
        !! individually (xx, yy, zz, xy).
        !!
        !! For more information see [[muscle_tensors]]
        real(real64), dimension(4) :: vals
            !! Stores the 4 independent components in Voigt notation: (xx, yy, zz, xy).
        contains
            generic, public :: init => init_ten_2D2Osym, init2_ten_2D2Osym
                !! Generic interface for initialization.
            procedure, private :: init_ten_2D2Osym, init2_ten_2D2Osym
            procedure, public :: xx 
                !! Accessor for the xx (1,1) component.
            procedure, public :: yy 
                !! Accessor for the yy (2,2) component.
            procedure, public :: zz
                !! Accessor for the zz (3,3) component.
            procedure, public :: xy
                !! Accessor for the xy (1,2) component.
            procedure, public :: square => square_2D2Osym
    end type ten_2D2Osym

    public :: operator(.approx.)
    interface operator (.approx.)
        module procedure approx_2D2Osym
    end interface

    public :: operator(+)
    interface operator (+)
        module procedure sum_2D2Osym
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_2D2Osym
        module procedure subU_2D2Osym
    end interface

    public :: operator(*)
    interface operator (*)
        module procedure mul_real64_2D2Osym
        module procedure mul_2D2Osym_real64
    end interface

    public :: operator( / )
    interface operator ( / )
        module procedure div_2D2Osym_real64
    end interface

    public :: operator(.dev.)
    interface operator (.dev.)
        module procedure dev_2D2Osym
    end interface

    public :: operator(.ddot.)
    interface operator (.ddot.)
        module procedure ddot_2D2Osym_2D2Osym
    end interface

    ! interface operator (.tdot.)
    !     module procedure tdot_2D2Osym_2D2Osym
    ! end interface

    public :: assignment (=)
    interface assignment (=)
        module procedure ten_2D2Osym_real64_assign
    end interface

contains

    pure subroutine ten_2D2Osym_real64_assign(a, b)
        implicit none
        type(ten_2D2Osym), intent(out) :: a
        real(real64), intent(in) :: b
        a%vals = b
    end subroutine

    pure subroutine init_ten_2D2Osym(self, vals)
        !! Initializes a ten_2D2Osym tensor from a 6-element array (Voigt order).
        !!
        !! voigt notation used: 11, 22, 33, 12, 23, 13
        implicit none
        class(ten_2D2Osym), intent(inout) :: self
        real(real64), intent(in) :: vals(4)
        self%vals = vals
    end subroutine

    pure subroutine init2_ten_2D2Osym(self, xx, yy, zz, xy)
        !! Initializes a ten_2D2Osym tensor from its 6 individual components.
        implicit none
        class(ten_2D2Osym), intent(inout) :: self
        real(real64), intent(in) :: xx, yy, zz, xy
        self%vals = (/xx, yy, zz, xy/)
    end subroutine

    pure function approx_2D2Osym(a, b) result(res)
        !! `.approx.` Compares two ten_2D2Osym tensors for approximate equality.
        !! Uses a modified L1 norm (shear components weighted by 2) with relative
        !! and absolute tolerances (EPS, EPS_ABS).
        !! norm(a) = |a_11| + |a_22| + |a_33| + 2|a_12|
        !! Condition: norm(a-b) / max(norm(a), norm(b), EPS_ABS) <= EPS
        implicit none
        type(ten_2D2Osym), intent(in) :: a, b
        logical :: res

        real(real64), parameter :: EPS=1e-7, EPS_ABS=1e-30
        real(real64) :: norm_a, norm_b, norm_max, norm

        norm_a =   abs(a%vals(1)) + abs(a%vals(2)) + abs(a%vals(3)) &
                 + 2*abs(a%vals(4))
        norm_b =   abs(b%vals(1)) + abs(b%vals(2)) + abs(b%vals(3)) &
                 + 2*abs(b%vals(4))

        norm_max = max(max(norm_a, norm_b), EPS_ABS)

        norm =     abs(a%vals(1)-b%vals(1)) +   abs(a%vals(2)-b%vals(2)) +   abs(a%vals(3)-b%vals(3)) &
               + 2*abs(a%vals(4)-b%vals(4))

        
        if (norm/norm_max .gt. EPS) res=.false.
        if (norm/norm_max .le. EPS) res=.true.
        
    end function approx_2D2Osym

    pure function sum_2D2Osym(a, b) result(res)
        implicit none
        type(ten_2D2Osym), intent(in) :: a, b
        type(ten_2D2Osym) :: res
        res%vals = a%vals + b%vals
    end function sum_2D2Osym

    pure function sub_2D2Osym(a, b) result(res)
        implicit none
        type(ten_2D2Osym), intent(in) :: a, b
        type(ten_2D2Osym) :: res
        res%vals = a%vals - b%vals
    end function sub_2D2Osym

    pure function subU_2D2Osym(a) result(res)
        implicit none
        type(ten_2D2Osym), intent(in) :: a
        type(ten_2D2Osym) :: res
        res%vals = -a%vals
    end function subU_2D2Osym

    pure function mul_real64_2D2Osym(a, b) result(res)
        implicit none
        real(real64), intent(in) :: a
        type(ten_2D2Osym), intent(in) :: b
        type(ten_2D2Osym) :: res
        res%vals = a * b%vals
    end function mul_real64_2D2Osym

    pure function mul_2D2Osym_real64(a, b) result(res)
        implicit none
        real(real64), intent(in) :: b
        type(ten_2D2Osym), intent(in) :: a
        type(ten_2D2Osym) :: res
        res%vals =  a%vals * b
    end function mul_2D2Osym_real64

    pure function div_2D2Osym_real64(a, b) result(res)
        implicit none
        real(real64), intent(in) :: b
        type(ten_2D2Osym), intent(in) :: a
        type(ten_2D2Osym) :: res
        res%vals = a%vals/b
    end function div_2D2Osym_real64

    pure function ddot_2D2Osym_2D2Osym(a, b) result(res)
        implicit none
        type(ten_2D2Osym), intent(in) :: a, b
        real(real64) :: res
        res =   a%vals(1)*b%vals(1) &
                             + a%vals(2)*b%vals(2) &
                             + a%vals(3)*b%vals(3) &
                             + 2*a%vals(4)*b%vals(4)
    end function ddot_2D2Osym_2D2Osym

    pure function dev_2D2Osym(a) result(res)
        implicit none
        type(ten_2D2Osym), intent(in) :: a
        type(ten_2D2Osym) :: res
        real(real64) :: hydro
        hydro = (a%vals(1) + a%vals(2) + a%vals(3))/3D0
        res%vals(1:3) = a%vals(1:3) - hydro
        res%vals(4) = a%vals(4)
    end function dev_2D2Osym

    pure function square_2D2Osym(a) result(res)
        implicit none
        class(ten_2D2Osym), intent(in) :: a
        type(ten_2D2Osym) :: res
        res%vals(1) = a%vals(1)**2 + a%vals(4)**2
        res%vals(2) = a%vals(2)**2 + a%vals(4)**2
        res%vals(3) = a%vals(3)**2
        res%vals(4) = a%vals(4)*(a%vals(1) + a%vals(2))
    end function square_2D2Osym
    
    pure function xx(a) result(res)
        !! Accessor function for the xx (11) component (vals(1)).
        implicit none
        class(ten_2D2Osym), intent(in) :: a
        real(real64) :: res
        res = a%vals(1)
    end function xx
    
    pure function yy(a) result(res)
        !! Accessor function for the yy (22) component (vals(2)).
        implicit none
        class(ten_2D2Osym), intent(in) :: a
        real(real64) :: res
        res = a%vals(2)
    end function yy

    pure function zz(a) result(res)
        !! Accessor function for the zz (33) component (vals(3)).
        implicit none
        class(ten_2D2Osym), intent(in) :: a
        real(real64) :: res
        res = a%vals(3)
    end function zz

    pure function xy(a) result(res)
        !! Accessor function for the xy (12) component (vals(4)).
        implicit none
        class(ten_2D2Osym), intent(in) :: a
        real(real64) :: res
        res = a%vals(4)
    end function xy

end module muscle_tensor_2d2osym