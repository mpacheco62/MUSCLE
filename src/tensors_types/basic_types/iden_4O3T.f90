module mod_iden_4O3T
    !! Module mod_iden_4O3T
    !! ======================
    !!
    !! Defines the standard (non-scaled) 3D fourth-order identity tensor \(\delta_{ij}\delta_{kl}\).
    !!
    !! This module provides the definition for the standard 3D fourth-order identity tensor
    !! corresponding to the tensor product of two second-order identity tensors, mathematically
    !! represented as \(\delta_{ij}\delta_{kl}\). This is referred to as the "type 3" identity
    !! in the library's nomenclature.
    !!
    !! Similar to the second-order identity `iden_2O`, the `iden_4O3T` type defined here
    !! is symbolic. It does not store any numerical data and implicitly represents the
    !! \(\delta_{ij}\delta_{kl}\) structure with a scaling factor of 1.0. Operations involving
    !! this type and scalars (multiplication or division) typically yield an instance of the
    !! corresponding scaled identity tensor type, `iden_4O3TS` (defined in `mod_iden_4O3TS`),
    !! which explicitly stores the resulting scaling factor.
    !!
    !! Public Entities
    !! ---------------
    !!
    !! ### Derived Type:
    !!
    !! - `iden_4O3T`: Represents the standard 3D fourth-order identity tensor of type 3, \(\delta_{ij}\delta_{kl}\).
    !!                It contains no data components and has an implicit scale factor of 1.0.
    !!
    !! ### Operators:
    !!
    !! The following operators are overloaded for interactions involving `iden_4O3T`:
    !!
    !! - `*`: Multiplication:
    !!     - `iden_4O3T * real(real64)`: Multiplies the identity by a scalar. Returns `iden_4O3TS` with `val = scalar`.
    !!     - `real(real64) * iden_4O3T`: Multiplies a scalar by the identity. Returns `iden_4O3TS` with `val = scalar`.
    !! - `/`: Division:
    !!     - `iden_4O3T / real(real64)`: Divides the identity by a scalar. Returns `iden_4O3TS` with `val = 1.0D0 / scalar`.
    !!
    !! Usage
    !! -----
    !!
    !! ```fortran
    !! program example_iden_4O3T_usage
    !!   use mod_iden_4O3T
    !!   use mod_iden_4O3TS
    !!   use iso_fortran_env, only: real64
    !!   implicit none
    !!
    !!   type(iden_4O3T) :: id4_tensor3     ! Represents delta_ij * delta_kl
    !!   type(iden_4O3TS) :: scaled_id4_3t
    !!   real(real64) :: factor = 10.0D0
    !!
    !!   ! Create scaled identity from the standard type 3 identity
    !!   scaled_id4_3t = id4_tensor3 * factor      ! scaled_id4_3t%val becomes 10.0
    !!   print *, "Scaled ID4 Type 3 (Mult):", scaled_id4_3t%val
    !!
    !!   scaled_id4_3t = id4_tensor3 / 4.0D0       ! scaled_id4_3t%val becomes 0.25
    !!   print *, "Scaled ID4 Type 3 (Div):", scaled_id4_3t%val
    !!
    !! end program example_iden_4O3T_usage
    !! ```
    !!
    !! For more information see [[tensors_types]]
    use, intrinsic :: iso_fortran_env
    use mod_iden_4O3TS
    implicit none
    private

    type, public :: iden_4O3T
        !! Standard 3D fourth-order identity tensor, type 3 (\(\delta_{ij}\delta_{kl}\)).
        !! ==========================================================================
        !!
        !! Type for the standard (non-scaled) 3D fourth-order identity tensor of type 3.
        !!
        !! Represents the standard 3D fourth-order identity tensor formed by the tensor product
        !! of two second-order identities: \(\mathbf{I} \otimes \mathbf{I}\), which corresponds to
        !! components \(\delta_{ij}\delta_{kl}\). This is referred to as the "type 3" identity
        !! in the library's nomenclature.
        !!
        !! This type acts as a symbolic marker and contains **no data components**. The implicit
        !! scaling factor is always considered to be `1.0D0`. It is primarily used in operations
        !! (like scalar multiplication/division) that typically result
        !! in a scaled identity tensor of type 3 (`iden_4O3TS`), which explicitly stores a scaling factor.
        !!
        !! For more information see [[tensors_types]]
    end type iden_4O3T


    public :: operator(*)
    interface operator (*)
        module procedure mul_I4O3T_real64
        module procedure mul_real64_I4O3T
    end interface

    public :: operator( / )
    interface operator ( / )
        module procedure div_I4O3T_real64
    end interface

contains

    pure module function mul_I4O3T_real64(I2, a) result(res)
        implicit none
        class(iden_4O3T), intent(in) :: I2
        real(real64), intent(in) :: a
        type(iden_4O3TS) :: res
        res%val = a 
    end function mul_I4O3T_real64

    pure module function mul_real64_I4O3T(a, I2) result(res)
        implicit none
        class(iden_4O3T), intent(in) :: I2
        real(real64), intent(in) :: a
        type(iden_4O3TS) :: res
        res%val = a 
    end function mul_real64_I4O3T

    pure module function div_I4O3T_real64(I2, a) result(res)
        implicit none
        class(iden_4O3T), intent(in) :: I2
        real(real64), intent(in) :: a
        type(iden_4O3TS) :: res
        res%val = 1D0/a
    end function div_I4O3T_real64

end module mod_iden_4O3T