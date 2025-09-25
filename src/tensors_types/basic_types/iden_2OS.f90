module mod_iden_2OS
    !! author: MPacheco
    !! version: 1.0 - Initial documentation
    !!
    !! Module mod_iden_2OS
    !! ======================
    !!
    !! This module defines the scaled 3D second-order identity \(c \cdot \delta_{ij}\).
    !!
    !! Denoted mathematically as \(c \delta_{ij}\), where \(c\) is a scalar value and
    !! \(\delta_{ij}\) is the Kronecker delta.
    !!
    !! Unlike its non-scaled counterpart `iden_2O` (defined in `mod_iden_2O`),
    !! the `iden_2OS` type explicitly stores the scaling factor \(c\) in its `val` component.
    !! This module provides the derived type definition and overloaded operators for
    !! arithmetic operations involving these scaled identity tensors and scalars.
    !! All arithmetic operations defined within this module result in another `iden_2OS` object.
    !!
    !! Public Entities
    !! ---------------
    !!
    !! ### Derived Type:
    !!
    !! - `iden_2OS`: Represents a scaled 3D second-order identity tensor \(c \delta_{ij}\).
    !!     - Component: `val :: real(real64)` - Stores the scaling factor \(c\).
    !!     - Generic Procedure: `init` - Used to initialize or set the `val` component.
    !!
    !! ### Operators:
    !!
    !! The following operators are overloaded for `iden_2OS` objects:
    !!
    !! - `+`: Addition:
    !!     - `iden_2OS + iden_2OS`: Adds two scaled identities. Result: `iden_2OS` with `val = I2a%val + I2b%val`.
    !! - `-`: Subtraction:
    !!     - `- iden_2OS`: Unary negation. Result: `iden_2OS` with `val = -I2a%val`.
    !!     - `iden_2OS - iden_2OS`: Subtracts two scaled identities. Result: `iden_2OS` with `val = I2a%val - I2b%val`.
    !! - `*`: Multiplication:
    !!     - `iden_2OS * real(real64)`: Scalar multiplication. Result: `iden_2OS` with `val = IMod%val * a`.
    !!     - `real(real64) * iden_2OS`: Scalar multiplication. Result: `iden_2OS` with `val = a * IMod%val`.
    !!     - `iden_2OS * iden_2OS`: Multiplies two scaled identities. Result: `iden_2OS` with `val = I2a%val * I2b%val`.
    !! - `/`: Division:
    !!     - `iden_2OS / real(real64)`: Scalar division. Result: `iden_2OS` with `val = IMod%val / a`.
    !!
    !! Usage
    !! -----
    !!
    !! ```fortran
    !! program example_iden_2Os_usage
    !!   use mod_iden_2OS
    !!   use iso_fortran_env, only: real64
    !!   implicit none
    !!
    !!   type(iden_2OS) :: scaled_id1, scaled_id2, result_id
    !!   real(real64) :: factor = 3.0D0
    !!
    !!   ! Initialize scaled identities
    !!   call scaled_id1%init(5.0D0)  ! scaled_id1%val is 5.0
    !!   scaled_id2%val = -2.0D0      ! scaled_id2%val is -2.0 (Direct assignment also works)
    !!
    !!   ! Perform operations
    !!   result_id = scaled_id1 + scaled_id2  ! result_id%val = 5.0 + (-2.0) = 3.0
    !!   print *, "Sum:", result_id%val
    !!
    !!   result_id = scaled_id1 * factor      ! result_id%val = 5.0 * 3.0 = 15.0
    !!   print *, "Scalar Multiplication:", result_id%val
    !!
    !!   result_id = -scaled_id2              ! result_id%val = -(-2.0) = 2.0
    !!   print *, "Unary Negation:", result_id%val
    !!
    !!   result_id = scaled_id1 / 2.0D0       ! result_id%val = 5.0 / 2.0 = 2.5
    !!   print *, "Scalar Division:", result_id%val
    !!
    !!   result_id = scaled_id1 * scaled_id2  ! result_id%val = 5.0 * (-2.0) = -10.0
    !!   print *, "Identity Multiplication:", result_id%val
    !!
    !! end program example_iden_2Os_usage
    !! ```
    !!
    !! for more information see [[tensors_types]]
    
    use, intrinsic :: iso_fortran_env
    implicit none
    private

    type, public :: iden_2OS  ! val*\delta_ij
        !! Scaled 3D second-order identity tensor \(c \cdot \delta_{ij}\).
        !! ===============================================================
        !!
        !! Type for the scaled 3D second-order identity tensor \(c \delta_{ij}\).
        !!
        !! Represents the scaled 3D second-order identity tensor, mathematically \(c \delta_{ij}\),
        !! where \(c\) is a scalar scaling factor and \(\delta_{ij}\) is the Kronecker delta.
        !! The scaling factor \(c\) is stored explicitly in the `val` component.
        !! This contrasts with the symbolic, non-scaled `iden_2O` type.
        !!
        !! for more information see [[tensors_types]]
        real(real64) :: val
            !! The scaling factor \(c\) for the identity tensor \(c \delta_{ij}\).
        contains
            generic, public :: init => init_iden_2OS
                !! Generic interface for initializing or setting the scaling factor `val`.
            procedure, private :: init_iden_2OS
    end type iden_2OS


    public :: operator(+)
    interface operator (+)
        module procedure sum_I2OS_I2OS
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure subU_I2OS
        module procedure sub_I2OS_I2OS
    end interface

    public :: operator(*)
    interface operator (*)
        module procedure mul_I2OS_real64
        module procedure mul_real64_I2OS
        module procedure mul_I2OS_I2OS
    end interface

    public :: operator( / )
    interface operator ( / )
        module procedure div_I2OS_real64
    end interface

    contains

    subroutine init_iden_2OS(self, val)
        !! Initializes or sets the scaling factor of an `iden_2OS` object.
        implicit none
        class(iden_2OS), intent(inout) :: self
            !! self The `iden_2OS` object to initialize/modify.
        real(real64), intent(in) :: val
            !! val The `real(real64)` value to assign as the scaling factor.
        self%val = val
    end subroutine init_iden_2OS

    pure function sum_I2OS_I2OS(I2a, I2b) result(res)
        implicit none
        class(iden_2OS), intent(in) :: I2a, I2b
        type(iden_2OS) :: res
        res%val = I2a%val + I2b%val
    end function sum_I2OS_I2OS

    pure function subU_I2OS(I2a) result(res)
        implicit none
        class(iden_2OS), intent(in) :: I2a
        type(iden_2OS) :: res
        res%val = -I2a%val
    end function subU_I2OS

    pure function sub_I2OS_I2OS(I2a, I2b) result(res)
        implicit none
        class(iden_2OS), intent(in) :: I2a, I2b
        type(iden_2OS) :: res
        res%val = I2a%val - I2b%val
    end function sub_I2OS_I2OS

    pure function mul_I2OS_I2OS(I2a, I2b) result(res)
        implicit none
        class(iden_2OS), intent(in) :: I2a, I2b
        type(iden_2OS) :: res
        res%val = I2a%val*I2b%val
    end function mul_I2OS_I2OS

    pure function mul_I2OS_real64(IMod, a) result(res)
        implicit none
        class(iden_2OS), intent(in) :: IMod
        real(real64), intent(in) :: a
        type(iden_2OS) :: res
        res%val = IMod%val * a 
    end function mul_I2OS_real64

    pure function mul_real64_I2OS(a, IMod) result(res)
        implicit none
        class(iden_2OS), intent(in) :: IMod
        real(real64), intent(in) :: a
        type(iden_2OS) :: res
        res%val = IMod%val * a 
    end function mul_real64_I2OS

    pure function div_I2OS_real64(IMod, a) result(res)
        implicit none
        class(iden_2OS), intent(in) :: IMod
        real(real64), intent(in) :: a
        type(iden_2OS) :: res
        res%val = IMod%val/a 
    end function div_I2OS_real64

end module mod_iden_2OS