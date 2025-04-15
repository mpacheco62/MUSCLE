module mod_iden_3D4O3TS
    !! Module mod_iden_3D4O3TS
    !! ======================
    !!
    !! Defines the scaled 3D fourth-order identity tensor of type 3.
    !!
    !! This module provides the definition for the scaled 3D fourth-order identity tensor
    !! corresponding to the tensor product of two second-order identity tensors, mathematically
    !! represented as \(c \cdot (\delta_{ij}\delta_{kl})\), where \(c\) is a scalar scaling factor.
    !! This is referred to as the "scaled type 3" identity in the library's nomenclature.
    !!
    !! Similar to `iden_3D2OS`, the `iden_3D4O3TS` type explicitly stores the scaling factor \(c\)
    !! in its `val` component. This contrasts with the non-scaled symbolic type `iden_3D4O3T`.
    !! This module defines the derived type and overloads operators for scalar multiplication
    !! and division involving this scaled identity tensor. Operations defined here result
    !! in another `iden_3D4O3TS` object.
    !!
    !! Public Entities
    !! ---------------
    !!
    !! ### Derived Type:
    !!
    !! - `iden_3D4O3TS`: Represents the scaled 3D fourth-order identity tensor of type 3, \(c \cdot (\delta_{ij}\delta_{kl})\).
    !!     - Component: `val :: real(real64)` - Stores the scaling factor \(c\).
    !!
    !! ### Operators:
    !!
    !! The following operators are overloaded for interactions involving `iden_3D4O3TS`:
    !!
    !! - `*`: Multiplication:
    !!     - `iden_3D4O3TS * real(real64)`: Multiplies the scaled identity by a scalar. Returns `iden_3D4O3TS` with `val = IMod%val * a`.
    !!     - `real(real64) * iden_3D4O3TS`: Multiplies a scalar by the scaled identity. Returns `iden_3D4O3TS` with `val = a * IMod%val`.
    !! - `/`: Division:
    !!     - `iden_3D4O3TS / real(real64)`: Divides the scaled identity by a scalar. Returns `iden_3D4O3TS` with `val = IMod%val / a`.
    !!
    !! Usage
    !! -----
    !!
    !! ```fortran
    !! program example_iden_3d4o3ts_usage
    !!   use mod_iden_3D4O3TS
    !!   use iso_fortran_env, only: real64
    !!   implicit none
    !!
    !!   type(iden_3D4O3TS) :: scaled_id4_t3_1, scaled_id4_t3_2
    !!   real(real64) :: factor = -2.0D0
    !!
    !!   ! Initialize scaled type 3 identity
    !!   scaled_id4_t3_1%val = 10.0D0     ! scaled_id4_t3_1%val is 10.0
    !!
    !!   ! Perform operations
    !!   scaled_id4_t3_2 = scaled_id4_t3_1 * factor  ! scaled_id4_t3_2%val = 10.0 * (-2.0) = -20.0
    !!   print *, "Scaled ID4 Type 3 (Mult):", scaled_id4_t3_2%val
    !!
    !!   scaled_id4_t3_2 = scaled_id4_t3_1 / 5.0D0   ! scaled_id4_t3_2%val = 10.0 / 5.0 = 2.0
    !!   print *, "Scaled ID4 Type 3 (Div):", scaled_id4_t3_2%val
    !!
    !! end program example_iden_3d4o3ts_usage
    !! ```
    !!
    !! For more information see [[tensors_types]]
    use, intrinsic :: iso_fortran_env
    implicit none
    private

    type, public :: iden_3D4O3TS
        !! Scaled 3D fourth-order identity tensor, type 3 (\(c \cdot \delta_{ij}\delta_{kl}\)).
        !! ================================================================================
        !!
        !! Type for the scaled 3D fourth-order identity tensor of type 3.
        !!
        !! Represents the scaled 3D fourth-order identity tensor formed by scaling the tensor product
        !! of two second-order identities: \(c \cdot (\mathbf{I} \otimes \mathbf{I})\), which corresponds to
        !! components \(c \cdot (\delta_{ij}\delta_{kl})\), where \(c\) is a scalar scaling factor.
        !! The scaling factor \(c\) is stored explicitly in the `val` component.
        !! This contrasts with the symbolic, non-scaled `iden_3D4O3T` type.
        !!
        !! For more information see [[tensors_types]]
        real(real64) :: val
            !! The scaling factor \(c\) for the type 3 identity tensor \(c \cdot (\delta_{ij}\delta_{kl})\).
    end type iden_3D4O3TS

    public :: operator(*)
    interface operator (*)
        module procedure mul_I3D4O3TS_real64
        module procedure mul_real64_I3D4O3TS
    end interface

    public :: operator( / )
    interface operator ( / )
        module procedure div_I3D4O3TS_real64
    end interface

    contains

    pure module function mul_I3D4O3TS_real64(IMod, a) result(res)
        implicit none
        class(iden_3D4O3TS), intent(in) :: IMod
        real(real64), intent(in) :: a
        type(iden_3D4O3TS) :: res
        res%val = IMod%val * a 
    end function mul_I3D4O3TS_real64

    pure module function mul_real64_I3D4O3TS(a, IMod) result(res)
        implicit none
        class(iden_3D4O3TS), intent(in) :: IMod
        real(real64), intent(in) :: a
        type(iden_3D4O3TS) :: res
        res%val = IMod%val * a 
    end function mul_real64_I3D4O3TS

    pure module function div_I3D4O3TS_real64(IMod, a) result(res)
        implicit none
        class(iden_3D4O3TS), intent(in) :: IMod
        real(real64), intent(in) :: a
        type(iden_3D4O3TS) :: res
        res%val = IMod%val/a 
    end function div_I3D4O3TS_real64

end module mod_iden_3D4O3TS