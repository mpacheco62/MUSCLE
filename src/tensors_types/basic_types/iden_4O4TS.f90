module mod_iden_4O4TS
    !! Module mod_iden_4O4TS
    !! ======================
    !!
    !! Defines the scaled 3D fourth-order symmetric identity tensor \(c \cdot \frac{1}{2}(\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk})\).
    !!
    !! This module provides the definition for the scaled 3D fourth-order symmetric identity tensor,
    !! mathematically represented as \(c \cdot \mathbf{I}^{sym} = c \cdot \frac{1}{2}(\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk})\),
    !! where \(c\) is a scalar scaling factor. This is referred to as the "scaled type 4" identity
    !! in the library's nomenclature.
    !!
    !! Similar to `iden_2OS`, the `iden_4O4TS` type explicitly stores the scaling factor \(c\)
    !! in its `val` component. This contrasts with the non-scaled symbolic type `iden_4O4T`.
    !! This module defines the derived type and overloads operators for scalar multiplication
    !! and division involving this scaled symmetric identity tensor. Operations defined here
    !! result in another `iden_4O4TS` object.
    !!
    !! Public Entities
    !! ---------------
    !!
    !! ### Derived Type:
    !!
    !! - `iden_4O4TS`: Represents the scaled 3D fourth-order symmetric identity tensor \(c \cdot \mathbf{I}^{sym}\).
    !!     - Component: `val :: real(real64)` - Stores the scaling factor \(c\).
    !!
    !! ### Operators:
    !!
    !! The following operators are overloaded for interactions involving `iden_4O4TS`:
    !!
    !! - `*`: Multiplication:
    !!     - `iden_4O4TS * real(real64)`: Multiplies the scaled identity by a scalar. Returns `iden_4O4TS` with `val = IMod%val * a`.
    !!     - `real(real64) * iden_4O4TS`: Multiplies a scalar by the scaled identity. Returns `iden_4O4TS` with `val = a * IMod%val`.
    !! - `/`: Division:
    !!     - `iden_4O4TS / real(real64)`: Divides the scaled identity by a scalar. Returns `iden_4O4TS` with `val = IMod%val / a`.
    !!
    !! Usage
    !! -----
    !!
    !! ```fortran
    !! program example_iden_4O4Ts_usage
    !!   use mod_iden_4O4TS
    !!   use iso_fortran_env, only: real64
    !!   implicit none
    !!
    !!   type(iden_4O4TS) :: scaled_id4_sym1, scaled_id4_sym2
    !!   real(real64) :: factor = 0.5D0
    !!
    !!   ! Initialize scaled symmetric identity
    !!   scaled_id4_sym1%val = 5.0D0      ! scaled_id4_sym1%val is 5.0
    !!
    !!   ! Perform operations
    !!   scaled_id4_sym2 = scaled_id4_sym1 * factor  ! scaled_id4_sym2%val = 5.0 * 0.5 = 2.5
    !!   print *, "Scaled ID4 Sym (Mult):", scaled_id4_sym2%val
    !!
    !!   scaled_id4_sym2 = scaled_id4_sym1 / 10.0D0  ! scaled_id4_sym2%val = 5.0 / 10.0 = 0.5
    !!   print *, "Scaled ID4 Sym (Div):", scaled_id4_sym2%val
    !!
    !! end program example_iden_4O4Ts_usage
    !! ```
    !!
    !! For more information see [[tensors_types]]
    use, intrinsic :: iso_fortran_env
    implicit none
    private

    type, public :: iden_4O4TS
        !! Scaled 3D fourth-order symmetric identity tensor \(c \cdot \frac{1}{2}(\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk})\).
        !! ============================================================================
        !!
        !! Type for the scaled 3D fourth-order symmetric identity tensor \(c \cdot \frac{1}{2}(\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk})\).
        !!
        !! Represents the scaled 3D fourth-order symmetric identity tensor, mathematically
        !! \(c \cdot \mathbf{I}^{sym} = c \cdot \frac{1}{2}(\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk})\),
        !! where \(c\) is a scalar scaling factor. This is often used as the identity element
        !! for operations involving symmetric fourth-order tensors.
        !! The scaling factor \(c\) is stored explicitly in the `val` component.
        !! This contrasts with the symbolic, non-scaled `iden_4O4T` type.
        !!
        !! For more information see [[tensors_types]]
        real(real64) :: val
            !! The scaling factor \(c\) for the symmetric identity tensor \(c \cdot \mathbf{I}^{sym}\).
    end type iden_4O4TS

    public :: operator(*)
    interface operator (*)
        module procedure mul_I3D4O4TS_real64
        module procedure mul_real64_I3D4O4TS
    end interface

    public :: operator( / )
    interface operator ( / )
        module procedure div_I3D4O4TS_real64
    end interface

    contains

    pure module function mul_I3D4O4TS_real64(IMod, a) result(res)
        implicit none
        class(iden_4O4TS), intent(in) :: IMod
        real(real64), intent(in) :: a
        type(iden_4O4TS) :: res
        res%val = IMod%val * a 
    end function mul_I3D4O4TS_real64

    pure module function mul_real64_I3D4O4TS(a, IMod) result(res)
        implicit none
        class(iden_4O4TS), intent(in) :: IMod
        real(real64), intent(in) :: a
        type(iden_4O4TS) :: res
        res%val = IMod%val * a 
    end function mul_real64_I3D4O4TS

    pure module function div_I3D4O4TS_real64(IMod, a) result(res)
        implicit none
        class(iden_4O4TS), intent(in) :: IMod
        real(real64), intent(in) :: a
        type(iden_4O4TS) :: res
        res%val = IMod%val/a 
    end function div_I3D4O4TS_real64

end module mod_iden_4O4TS