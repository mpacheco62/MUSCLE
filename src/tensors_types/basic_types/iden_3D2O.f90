module mod_iden_3D2O
    !! author: MPacheco
    !! version: 1.0 - Initial documentation
    !!
    !! Module mod_iden_3D2O
    !! =====================
    !!
    !! This module defines the standard (non-scaled) 3D second-order identity \(\delta_{ij}\).
    !!
    !! This module defines the standard (non-scaled) 3D second-order identity tensor type,
    !! denoted mathematically as \(\delta_{ij}\). It represents the identity matrix.
    !! The module also provides overloaded operators for arithmetic operations involving
    !! this identity tensor, often interacting with its scaled counterpart (`iden_3D2OS`)
    !! defined in `mod_iden_3D2OS`.
    !!
    !! The `iden_3D2O` type itself does not store any numerical value; it acts as a
    !! symbolic representation of the identity tensor where the implicit scaling factor is 1.0.
    !! Operations involving `iden_3D2O` typically result in a scaled identity tensor
    !! (`iden_3D2OS`), which explicitly stores the scaling factor.
    !!
    !! Public Entities
    !! ---------------
    !!
    !! ### Derived Type:
    !!
    !! - `iden_3D2O`: Represents the standard 3D second-order identity tensor \(\delta_{ij}\).
    !!                It has no components and implicitly represents the identity with a
    !!                scaling factor of 1.0.
    !!
    !! ### Operators:
    !!
    !! The following operators are overloaded. Note that most operations involving `iden_3D2O`
    !! return a scaled identity tensor `iden_3D2OS`.
    !!
    !! - `+`: Addition:
    !!     - `iden_3D2O + iden_3D2O`: Results in `iden_3D2OS` with value `2.0`.
    !!     - `iden_3D2O + iden_3D2OS`: Results in `iden_3D2OS` (value `1.0 + iden_3D2OS%val`).
    !!     - `iden_3D2OS + iden_3D2O`: Results in `iden_3D2OS` (value `iden_3D2OS%val + 1.0`).
    !! - `-`: Subtraction:
    !!     - `- iden_3D2O`: Unary negation. Results in `iden_3D2OS` with value `-1.0`.
    !!     - `iden_3D2O - iden_3D2O`: Results in `iden_3D2OS` with value `0.0`.
    !!     - `iden_3D2O - iden_3D2OS`: Results in `iden_3D2OS` (value `1.0 - iden_3D2OS%val`).
    !!     - `iden_3D2OS - iden_3D2O`: Results in `iden_3D2OS` (value `iden_3D2OS%val - 1.0`).
    !! - `*`: Multiplication:
    !!     - `iden_3D2O * real(real64)`: Scalar multiplication. Results in `iden_3D2OS`.
    !!     - `real(real64) * iden_3D2O`: Scalar multiplication. Results in `iden_3D2OS`.
    !!     - `iden_3D2O * iden_3D2O`: Identity tensor multiplication (conceptually \(I \cdot I = I\)). Results in `iden_3D2O`.
    !! - `/`: Division:
    !!     - `iden_3D2O / real(real64)`: Scalar division. Results in `iden_3D2OS` (value `1.0 / scalar`).
    !!
    !! Usage
    !! -----
    !!
    !! ```fortran
    !! program example_iden_3d2o_usage
    !!   use mod_iden_3D2O
    !!   use mod_iden_3D2OS
    !!   use iso_fortran_env, only: real64
    !!   implicit none
    !!
    !!   type(iden_3D2O) :: id_tensor       ! Represents delta_ij
    !!   type(iden_3D2OS) :: scaled_id1, scaled_id2
    !!   real(real64) :: factor = 5.0D0
    !!
    !!   ! Create scaled identities from the standard identity
    !!   scaled_id1 = id_tensor * factor      ! scaled_id1%val = 5.0
    !!   scaled_id2 = id_tensor / 2.0D0       ! scaled_id2%val = 0.5
    !!
    !!   ! Operations involving iden_3D2O
    !!   scaled_id1 = id_tensor + id_tensor   ! scaled_id1%val = 2.0
    !!   scaled_id2 = -id_tensor              ! scaled_id2%val = -1.0
    !!   scaled_id1 = scaled_id1 + id_tensor  ! scaled_id1%val = 2.0 + 1.0 = 3.0
    !!
    !!   print *, "Scaled ID 1 value:", scaled_id1%val ! Expected: 3.0
    !!   print *, "Scaled ID 2 value:", scaled_id2%val ! Expected: -1.0
    !!
    !! end program example_iden_3d2o_usage
    !! ```
    !! For more information see [[tensors_types]]
    use, intrinsic :: iso_fortran_env
    use mod_iden_3D2OS
    implicit none
    private

    type, public :: iden_3D2O
        !! 3D second-order identity tensor \(\delta_{ij}\)
        !! ===============================================
        !! 
        !! Type for the standard (non-scaled) 3D second-order identity tensor \(\delta_{ij}\).
        !!
        !! Represents the standard 3D second-order identity tensor, mathematically the Kronecker delta \(\delta_{ij}\).
        !! This type acts as a symbolic marker and contains **no data components**. The implicit
        !! scaling factor is always considered to be `1.0D0`. It is primarily used in operations
        !! (like addition, subtraction, or scalar multiplication/division) that typically result
        !! in a scaled identity tensor (`iden_3D2OS`), which explicitly stores a scaling factor.
        !! Multiplying two `iden_3D2O` types results in another `iden_3D2O` type (\(I \cdot I = I\)).
        !! 
        !! For more information see [[tensors_types]]
    end type iden_3D2O


    public :: operator(+)
    interface operator (+)
        module procedure sum_I3D2O_I3D2O
        module procedure sum_I3D2O_I3D2OS
        module procedure sum_I3D2OS_I3D2O
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure subU_I3D2O
        module procedure sub_I3D2O_I3D2O
        module procedure sub_I3D2O_I3D2OS
        module procedure sub_I3D2OS_I3D2O
    end interface

    public :: operator(*)
    interface operator (*)
        module procedure mul_I3D2O_real64
        module procedure mul_real64_I3D2O
        module procedure mul_I3D2O_I3D2O
    end interface

    public :: operator( / )
    interface operator ( / )
        module procedure div_I3D2O_real64
    end interface

contains

    pure function sum_I3D2O_I3D2O(I2a, I2b) result(res)
        !! `+` Adds two standard identity tensors. Result is a scaled identity with value 2.0.
        implicit none
        class(iden_3D2O), intent(in) :: I2a, I2b
        type(iden_3D2OS) :: res
        res%val = 2D0
    end function sum_I3D2O_I3D2O

    pure function sum_I3D2O_I3D2OS(I2a, I2b) result(res)
        !! `+` Adds a standard identity tensor and a scaled identity tensor.
        implicit none
        class(iden_3D2O), intent(in) :: I2a
        class(iden_3D2OS), intent(in) :: I2b
        type(iden_3D2OS) :: res
        res%val = 1D0 + I2b%val
    end function sum_I3D2O_I3D2OS

    pure function sum_I3D2OS_I3D2O(I2a, I2b) result(res)
        !! `+` Adds a scaled identity tensor and a standard identity tensor.
        implicit none
        class(iden_3D2OS), intent(in) :: I2a
        class(iden_3D2O), intent(in) :: I2b
        type(iden_3D2OS) :: res
        res%val = I2a%val + 1D0
    end function sum_I3D2OS_I3D2O

    pure function subU_I3D2O(I2a) result(res)
        !! `-` Unary negation of the standard identity tensor. Result is a scaled identity with value -1.0.
        implicit none
        class(iden_3D2O), intent(in) :: I2a
        type(iden_3D2OS) :: res
        res%val = -1D0
    end function subU_I3D2O

    pure function sub_I3D2O_I3D2O(I2a, I2b) result(res)
        !! `-` Subtracts one standard identity tensor from another. Result is a scaled identity with value 0.0.
        implicit none
        class(iden_3D2O), intent(in) :: I2a, I2b
        type(iden_3D2OS) :: res
        res%val = 0D0
    end function sub_I3D2O_I3D2O

    pure function sub_I3D2O_I3D2OS(I2a, I2b) result(res)
        !! `-` Subtracts a scaled identity tensor from a standard identity tensor.
        implicit none
        class(iden_3D2O), intent(in) :: I2a
        class(iden_3D2OS), intent(in) :: I2b
        type(iden_3D2OS) :: res
        res%val = 1D0 - I2b%val
    end function sub_I3D2O_I3D2OS

    pure function sub_I3D2OS_I3D2O(I2a, I2b) result(res)
        !! `-` Subtracts a standard identity tensor from a scaled identity tensor.
        implicit none
        class(iden_3D2OS), intent(in) :: I2a
        class(iden_3D2O), intent(in) :: I2b
        type(iden_3D2OS) :: res
        res%val = I2a%val - 1D0
    end function sub_I3D2OS_I3D2O

    pure function mul_I3D2O_real64(I2, a) result(res)
        !! `*` Multiplies the standard identity tensor by a scalar. Result is a scaled identity.
        implicit none
        class(iden_3D2O), intent(in) :: I2
        real(real64), intent(in) :: a
        type(iden_3D2OS) :: res
        res%val = a 
    end function mul_I3D2O_real64

    pure function mul_real64_I3D2O(a, I2) result(res)
        !! `*` Multiplies a scalar by the standard identity tensor. Result is a scaled identity.
        implicit none
        class(iden_3D2O), intent(in) :: I2
        real(real64), intent(in) :: a
        type(iden_3D2OS) :: res
        res%val = a 
    end function mul_real64_I3D2O

    pure module function mul_I3D2O_I3D2O(I2a, I2b) result(res)
        !! `\cdot` Multiplies two standard identity tensors. Conceptually I*I = I. Result is the standard identity
        implicit none
        class(iden_3D2O), intent(in) :: I2a, I2b
        type(iden_3D2O) :: res
    end function mul_I3D2O_I3D2O

    pure module function div_I3D2O_real64(I2, a) result(res)
        !! Divides the standard identity tensor by a scalar. Result is a scaled identity.
        implicit none
        class(iden_3D2O), intent(in) :: I2
        real(real64), intent(in) :: a
        type(iden_3D2OS) :: res
        res%val = 1D0/a
    end function div_I3D2O_real64

end module mod_iden_3D2O