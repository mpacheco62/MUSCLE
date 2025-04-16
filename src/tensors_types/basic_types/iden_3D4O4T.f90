module mod_iden_3D4O4T
    !! Module mod_iden_3D4O4T
    !! ======================
    !!
    !! Defines the standard (non-scaled) 3D fourth-order symmetric identity tensor \(\frac{1}{2}(\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk})\).
    !!
    !! This module provides the definition for the standard 3D fourth-order symmetric identity tensor,
    !! mathematically represented as \(\mathbf{I}^{sym} = \frac{1}{2}(\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk})\).
    !! This is referred to as the "type 4" identity in the library's nomenclature and is often used
    !! as the identity element for operations involving symmetric fourth-order tensors.
    !!
    !! Similar to other non-scaled identity types (`iden_3D2O`, `iden_3D4O3T`), the `iden_3D4O4T`
    !! type defined here is symbolic. It does not store any numerical data and implicitly represents the
    !! symmetric identity structure with a scaling factor of 1.0. Operations involving
    !! this type and scalars (multiplication or division) typically yield an instance of the
    !! corresponding scaled identity tensor type, `iden_3D4O4TS` (defined in `mod_iden_3D4O4TS`),
    !! which explicitly stores the resulting scaling factor.
    !!
    !! Public Entities
    !! ---------------
    !!
    !! ### Derived Type:
    !!
    !! - `iden_3D4O4T`: Represents the standard 3D fourth-order symmetric identity tensor \(\mathbf{I}^{sym}\).
    !!                It contains no data components and has an implicit scale factor of 1.0.
    !!
    !! ### Operators:
    !!
    !! The following operators are overloaded for interactions involving `iden_3D4O4T`:
    !!
    !! - `*`: Multiplication:
    !!     - `iden_3D4O4T * real(real64)`: Multiplies the identity by a scalar. Returns `iden_3D4O4TS` with `val = scalar`.
    !!     - `real(real64) * iden_3D4O4T`: Multiplies a scalar by the identity. Returns `iden_3D4O4TS` with `val = scalar`.
    !! - `/`: Division:
    !!     - `iden_3D4O4T / real(real64)`: Divides the identity by a scalar. Returns `iden_3D4O4TS` with `val = 1.0D0 / scalar`.
    !!
    !! Usage
    !! -----
    !!
    !! ```fortran
    !! program example_iden_3d4o4t_usage
    !!   use mod_iden_3D4O4T
    !!   use mod_iden_3D4O4TS
    !!   use iso_fortran_env, only: real64
    !!   implicit none
    !!
    !!   type(iden_3D4O4T) :: id4_sym_tensor   ! Represents I^sym
    !!   type(iden_3D4O4TS) :: scaled_id4_sym
    !!   real(real64) :: factor = 7.5D0
    !!
    !!   ! Create scaled symmetric identity from the standard type 4 identity
    !!   scaled_id4_sym = id4_sym_tensor * factor      ! scaled_id4_sym%val becomes 7.5
    !!   print *, "Scaled ID4 Sym (Mult):", scaled_id4_sym%val
    !!
    !!   scaled_id4_sym = id4_sym_tensor / 2.0D0       ! scaled_id4_sym%val becomes 0.5
    !!   print *, "Scaled ID4 Sym (Div):", scaled_id4_sym%val
    !!
    !! end program example_iden_3d4o4t_usage
    !! ```
    !!
    !! For more information see [[tensors_types]]
    use, intrinsic :: iso_fortran_env
    use mod_iden_3D4O4TS
    implicit none
    private

    type, public :: iden_3D4O4T
        !! Standard 3D fourth-order symmetric identity tensor \(\mathbf{I}^{sym}\).
        !! =====================================================================
        !!
        !! Type for the standard (non-scaled) 3D fourth-order symmetric identity tensor (type 4).
        !!
        !! Represents the standard 3D fourth-order symmetric identity tensor, mathematically
        !! \(\mathbf{I}^{sym} = \frac{1}{2}(\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk})\).
        !! This is often used as the identity element for operations involving symmetric
        !! fourth-order tensors.
        !!
        !! This type acts as a symbolic marker and contains **no data components**. The implicit
        !! scaling factor is always considered to be `1.0D0`. It is primarily used in operations
        !! (like scalar multiplication/division) that typically result
        !! in a scaled symmetric identity tensor (`iden_3D4O4TS`), which explicitly stores a scaling factor.
        !!
        !! For more information see [[tensors_types]]
    end type iden_3D4O4T


    public :: operator(*)
    interface operator (*)
        module procedure mul_I3D4O4T_real64
        module procedure mul_real64_I3D4O4T
    end interface

    public :: operator( / )
    interface operator ( / )
        module procedure div_I3D4O4T_real64
    end interface

contains

    pure module function mul_I3D4O4T_real64(I2, a) result(res)
        implicit none
        class(iden_3D4O4T), intent(in) :: I2
        real(real64), intent(in) :: a
        type(iden_3D4O4TS) :: res
        res%val = a 
    end function mul_I3D4O4T_real64

    pure module function mul_real64_I3D4O4T(a, I2) result(res)
        implicit none
        class(iden_3D4O4T), intent(in) :: I2
        real(real64), intent(in) :: a
        type(iden_3D4O4TS) :: res
        res%val = a 
    end function mul_real64_I3D4O4T

    pure module function div_I3D4O4T_real64(I2, a) result(res)
        implicit none
        class(iden_3D4O4T), intent(in) :: I2
        real(real64), intent(in) :: a
        type(iden_3D4O4TS) :: res
        res%val = 1D0/a
    end function div_I3D4O4T_real64

end module mod_iden_3D4O4T