module tensors_types
    !! author: MPacheco
    !! version: 1.0 - Initial documentation
    !!
    !! Tensors Types Module
    !! ====================
    !!
    !! This module serves as the central access point for various tensor types.
    !! 
    !! This module serves as the central access point for various tensor types and their associated
    !! operations within the UmatLib library. It aggregates tensor definitions and operator
    !! implementations from specialized submodules, providing a unified interface for users.
    !!
    !! Instead of defining the tensor types and operations directly, this module uses other
    !! modules (`mod_ten_*`, `mod_iden_*`, `mod_operator_*`) and makes specific types and
    !! operators publicly available. This promotes modularity and organization within the library.
    !!
    !! Nomenclature
    !! ------------
    !!
    !! The module exposes two main categories of tensor types: Regular Tensors and Identity Tensors.
    !!
    !! ### Regular Tensors:
    !!
    !! Nomenclature `ten_xyz`:
    !!
    !! - `x`: Dimensions of the tensor (e.g., 3D)
    !! - `y`: Order of the tensor (e.g., 2O, 4O)
    !! - `z`: Number of symmetries of the tensor (e.g., sym, 2sym, 3sym)
    !!
    !! Examples:
    !!
    !! - `ten_3D2O`: 3D tensor, second order, no specific symmetry assumed in the type name (likely stored as 3x3).
    !! - `ten_3D2Osym`: 3D tensor, second order, symmetric (likely stored efficiently, e.g., 6 components).
    !! - `ten_3D4O2sym`: 3D tensor, fourth order, with minor symmetries (\(C_{ijkl} = C_{jikl} = C_{ijlk}\)).
    !! - `ten_3D4O3sym`: 3D tensor, fourth order, with major and minor symmetries (\(C_{ijkl} = C_{jikl} = C_{ijlk} = C_{klij}\)).
    !!
    !! ### Identity Tensors:
    !!
    !! Nomenclature `iden_xyz` or `iden_xyzS`:
    !!
    !! - `x`: Dimensions of the tensor (e.g., 3D)
    !! - `y`: Order of the tensor (e.g., 2O, 4O)
    !! - `z`: Type of identity (e.g., None, 1T, 2T, 3T, 4T)
    !!     * `None`: Standard second-order identity tensor \(\delta_{ij}\). (`iden_2O`)
    !!     * `1T`: \(\delta_{ik}\delta_{jl}\), only valid for a fourth order tensor.
    !!     * `2T`: \(\delta_{il}\delta_{jk}\), only valid for a fourth order tensor.
    !!     * `3T`: \(\delta_{ij}\delta_{kl}\), only valid for a fourth order tensor. (`iden_4O3T`)
    !!     * `4T`: \(\frac{\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk}}{2}\), symmetric fourth-order identity. (`iden_4O4T`)
    !! - `S`: Suffix indicating a *scaled* identity tensor (e.g., `val * identity`). (`iden_2OS`, `iden_4O3TS`, `iden_4O4TS`)
    !!
    !! Examples:
    !!
    !! - `iden_2O`: 3D identity tensor, second order (\(\delta_{ij}\)).
    !! - `iden_2OS`: Scaled 3D identity tensor, second order (\(c \delta_{ij}\)).
    !! - `iden_4O3T`: 3D identity tensor, fourth order, type 3 (\(\delta_{ij}\delta_{kl}\)).
    !! - `iden_4O4T`: 3D identity tensor, fourth order, type 4 (symmetric identity).
    !! - `iden_4O4TS`: Scaled 3D identity tensor, fourth order, type 4 (\(c \cdot (\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk})/2\)).
    !!
    !! Public Entities
    !! ---------------
    !!
    !! This module makes the following entities publicly available:
    !!
    !! ### Derived Types:
    !!
    !! - [[ten_3D2O]]:       General 3D second-order tensor (from [[mod_ten_3D2O]]).
    !! - [[ten_3D2Osym]]:    Symmetric 3D second-order tensor (from [[mod_ten_3D2Osym]]).
    !! - [[ten_3D4O2sym]]:   3D fourth-order tensor with minor symmetries (from [[mod_ten_3D4O2sym]]).
    !! - [[ten_3D4O3sym]]:   3D fourth-order tensor with major and minor symmetries (from [[mod_ten_3D4O3sym]]).
    !! - [[iden_2O]]:      Standard 3D second-order identity tensor (from [[mod_iden_2O]]).
    !! - [[iden_2OS]]:   Scaled 3D second-order identity tensor (from [[mod_iden_2OS]]).
    !! - [[iden_4O3T]]:    Type 3 3D fourth-order identity tensor (from [[mod_iden_4O3T]]).
    !! - [[iden_4O3TS]]: Scaled Type 3 3D fourth-order identity tensor (from [[mod_iden_4O3TS]]).
    !! - [[iden_4O4T]]:    Type 4 (symmetric) 3D fourth-order identity tensor (from [[mod_iden_4O4T]]).
    !! - [[iden_4O4TS]]: Scaled Type 4 3D fourth-order identity tensor (from [[mod_iden_4O4TS]]).
    !!
    !! ### Operators:
    !!
    !! The following operators are overloaded for various combinations of the public tensor types
    !! and intrinsic types (like `real(real64)`). The implementations are provided by the
    !! `mod_operator_*` modules.
    !!
    !! - `+`: Addition (e.g., `tensor + tensor`).
    !! - `-`: Subtraction (e.g., `tensor - tensor`).
    !! - `*`: Multiplication (scalar * tensor, tensor * scalar, potentially tensor * tensor depending on definitions).
    !! - `/`: Division (tensor / scalar).
    !! - `.isequal.`: Custom equality comparison for tensors (checks for approximate equality of components).
    !! - `.dev.`: Deviatoric part of a second-order tensor.
    !! - `.ddot.`: Double dot product (e.g., `tensor4 : tensor2`, `tensor2 : tensor2`).
    !! - `.tdot.`: Tensor dot product (specific definition depends on implementation, often \( (A \otimes B)_{ijkl} = A_{ij} B_{kl} \) for second order).
    !! - `.inv.`: Inverse of a 3D fourth-order tensor, with major and minor symmetries.
    !!
    !! ### Assignment:
    !!
    !! - `=`: Overloaded assignment allows copying between compatible tensor types.
    !!
    !! Usage
    !! -----
    !!
    !! To use the tensor types and operations defined here, simply add `use tensors_types`
    !! to your Fortran code.
    !!
    !! ```fortran
    !! program example_usage
    !!   use tensors_types
    !!   implicit none
    !!
    !!   type(ten_3D2Osym) :: stress, strain, stress_dev
    !!   type(ten_3D4O3sym) :: C_elastic
    !!   real(real64) :: scalar_val
    !!
    !!   ! ... initialize tensors ...
    !!
    !!   stress = C_elastic .ddot. strain  ! Double dot product via overloaded operator
    !!   stress_dev = .dev. stress   ! Deviatoric part
    !!   stress = stress + stress_dev * scalar_val ! Addition, multiplication
    !!
    !!   ! ... etc ...
    !!
    !! end program example_usage
    !! ```

    use, intrinsic :: iso_fortran_env
    ! use mod_ten_3D2O, only : ten_3D2O
    use mod_ten_3D2O
    use mod_ten_3D2Osym
    use mod_ten_3D4O3sym
    use mod_ten_3D4O2sym
    use mod_iden_2OS
    use mod_iden_2O
    use mod_iden_4O3TS
    use mod_iden_4O3T
    use mod_iden_4O4TS
    use mod_iden_4O4T
    use mod_operator_I2O_3D2Osym
    use mod_operator_I2O_3D2O
    use mod_operator_I2OS_3D2Osym
    use mod_operator_I2OS_3D2O
    use mod_operator_I4O3T_3D4O3sym
    use mod_operator_I4O3TS_3D4O3sym
    use mod_operator_I4O4T_3D4O3sym
    use mod_operator_I4O4TS_3D4O3sym
    use mod_operator_3D2Osym_3D2O
    use mod_operator_3D2Osym_3D4O3sym
    use mod_operator_3D4O2sym_3D4O3sym
    use mod_operator_3D2Osym_3D4O2sym
    use mod_operator_I4O3T_3D4O2sym
    use mod_operator_I4O4T_3D4O2sym

    implicit None

    public :: ten_3D2O
    public :: ten_3D2Osym
    public :: ten_3D4O3sym
    public :: ten_3D4O2sym
    public :: iden_2OS
    public :: iden_2O
    public :: iden_4O3TS
    public :: iden_4O3T
    public :: iden_4O4TS
    public :: iden_4O4T

    private 
    ! INCLUDE 'ten_3D4O3sym_3D2Osym/ten_3D4O3sym_3D2Osym.inc'

    ! iden_4O3T ! \delta_ij\delta_kl

    ! iden_4O3TS ! val*\delta_ij\delta_kl

    ! type, public :: iden_4O4T ! (\delta_ik\delta_jl + \delta_il\delta_jk)/2
    ! end type iden_4O4T

    ! iden_4O4TS ! val*(\delta_ik\delta_jl + \delta_il\delta_jk)/2
    

    public :: operator(.isequal.)
    public :: operator(+)
    public :: operator(-)
    public :: operator(*)
    public :: operator( / )
    
    public :: operator(.dev.)

    public :: operator(.ddot.)
    ! interface operator (.ddot.)
    !     module procedure ddot_3D4O3sym_3D2Osym
    !     module procedure ddot_3D2Osym_3D4O3sym
    ! end interface

    public :: operator(.tdot.)
    public :: operator(.inv.)
    

    public :: assignment (=)

    ! interface ! iden_4O4T
    !     module pure function sum_I4O4T_3D4O3sym(I4, a) result(res)
    !         implicit none
    !         class(iden_4O4T), intent(in) :: I4
    !         type(ten_3D4O3sym), intent(in) :: a
    !         type(ten_3D4O3sym) :: res
    !     end function sum_I4O4T_3D4O3sym

    !     module pure function sum_3D4O3sym_I4O4T(a, I4) result(res)
    !         implicit none
    !         class(iden_4O4T), intent(in) :: I4
    !         type(ten_3D4O3sym), intent(in) :: a
    !         type(ten_3D4O3sym) :: res
    !     end function sum_3D4O3sym_I4O4T

    !     module pure function sub_I4O4T_3D4O3sym(I4, a) result(res)
    !         implicit none
    !         class(iden_4O4T), intent(in) :: I4
    !         type(ten_3D4O3sym), intent(in) :: a
    !         type(ten_3D4O3sym) :: res
    !     end function sub_I4O4T_3D4O3sym

    !     module pure function sub_3D4O3sym_I4O4T(a, I4) result(res)
    !         implicit none
    !         class(iden_4O4T), intent(in) :: I4
    !         type(ten_3D4O3sym), intent(in) :: a
    !         type(ten_3D4O3sym) :: res
    !     end function sub_3D4O3sym_I4O4T

    !     module pure function mul_I4O4T_real64(I4, a) result(res)
    !         implicit none
    !         class(iden_4O4T), intent(in) :: I4
    !         real(real64), intent(in) :: a
    !         type(ten_3D4O3sym) :: res
    !     end function mul_I4O4T_real64

    !     module pure function mul_real64_I4O4T(a, I4) result(res)
    !         implicit none
    !         class(iden_4O4T), intent(in) :: I4
    !         real(real64), intent(in) :: a
    !         type(ten_3D4O3sym) :: res
    !     end function mul_real64_I4O4T

    !     module pure function div_I4O4T_real64(I4, a) result(res)
    !         implicit none
    !         class(iden_4O4T), intent(in) :: I4
    !         real(real64), intent(in) :: a
    !         type(ten_3D4O3sym) :: res
    !     end function div_I4O4T_real64
    ! end interface

    ! interface  !O4 3D 3Sym

    !     module pure function ddot_3D4O3sym_3D2Osym(a, b) result(res)
    !         !
    !         !  | ( 1:1111) ( 7:1122) (12:1133) (16:1112) (19:1123) (21:1113) |
    !         !  | ( 7:2211) ( 2:2222) ( 8:2233) (13:2212) (17:2223) (20:2213) |
    !         !  | (12:3311) ( 8:3322) ( 3:3333) ( 9:3312) (14:3323) (18:3313) |
    !         !  | (16:1211) (13:1222) ( 9:1233) ( 4:1212) (10:1223) (15:1213) |
    !         !  | (19:2311) (17:2322) (14:2333) (10:2312) ( 5:2323) (11:2313) |
    !         !  | (21:1311) (20:1322) (18:1333) (15:1312) (11:1323) ( 6:1313) |
    !         implicit none
    !         class(ten_3D4O3sym), intent(in) :: a
    !         class(ten_3D2Osym), intent(in) :: b
    !         type(ten_3D2Osym) :: res
    !     end function ddot_3D4O3sym_3D2Osym

    !     module pure function ddot_3D2Osym_3D4O3sym(b, a) result(res)
    !         !
    !         !  | ( 1:1111) ( 7:1122) (12:1133) (16:1112) (19:1123) (21:1113) |
    !         !  | ( 7:2211) ( 2:2222) ( 8:2233) (13:2212) (17:2223) (20:2213) |
    !         !  | (12:3311) ( 8:3322) ( 3:3333) ( 9:3312) (14:3323) (18:3313) |
    !         !  | (16:1211) (13:1222) ( 9:1233) ( 4:1212) (10:1223) (15:1213) |
    !         !  | (19:2311) (17:2322) (14:2333) (10:2312) ( 5:2323) (11:2313) |
    !         !  | (21:1311) (20:1322) (18:1333) (15:1312) (11:1323) ( 6:1313) |
    !         implicit none
    !         class(ten_3D4O3sym), intent(in) :: a
    !         class(ten_3D2Osym), intent(in) :: b
    !         type(ten_3D2Osym) :: res            
    !     end function ddot_3D2Osym_3D4O3sym
    ! end interface

end module tensors_types