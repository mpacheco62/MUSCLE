module tensors_types
    !! Tensors types
    !! =============
    !!
    !! This module contains the definitions of tensors. So far, there are two tensors implemented.
    !! 
    !! Regular Tensors:
    !! ----------------
    !!
    !! Nomenclature `ten_xyz`:
    !!
    !! - `x`: Dimensions of the tensor (e.g., 3D)
    !! - `y`: Order of the tensor (e.g., 2O, 4O)
    !! - `z`: Number of symmetries of the tensor (e.g., sym, 3sym)
    !!
    !! Examples:
    !!
    !! - `ten_3D2Osym`: 3D tensor, second order, with symmetries.
    !! - `ten_3D4O3sym`: 3D tensor, fourth order, with major and minor, symmetries.
    !!
    !! Identity Tensors:
    !! ----------------
    !!
    !! Nomenclature `iden_xyz`:
    !!
    !! - `x`: Dimensions of the tensor (e.g., 3D)
    !! - `y`: Order of the tensor (e.g., 2O, 4O)
    !! - `z`: Type of identity (e.g., None, 1T, 2T, 3T, 4T)
    !!     * None: Only valid for a second order tensor
    !!     * `1`: \(\delta_{ik}\delta_{jl}\), only valid for a fourth order tensor
    !!     * `2`: \(\delta_{il}\delta_{jk}\), only valid for a fourth order tensor
    !!     * `3`: \(\delta_{ij}\delta_{kl}\), only valid for a fourth order tensor
    !!     * `4`: \(\frac{\delta_{ik}\delta_{jl} + \delta_{il}\delta_{jk}}{2}\), only valid for a fourth order tensor
    !!
    !! Examples:
    !!
    !! - `iden_3D2O`: 3D indentity type, second order.
    !! - `iden_3D4O4T`: 3D identity type, fourth order, fourth type.
    use, intrinsic :: iso_fortran_env
    ! use mod_ten_3D2O, only : ten_3D2O
    use mod_ten_3D2O
    use mod_ten_3D2Osym
    use mod_ten_3D4O3sym
    use mod_iden_3D2OMod
    use mod_iden_3D2O
    use mod_iden_3D4O3TMod
    use mod_iden_3D4O3T
    use mod_iden_3D4O4TMod
    use mod_iden_3D4O4T
    use mod_operator_I3D2O_3D2Osym
    use mod_operator_I3D2O_3D2O
    use mod_operator_I3D2OMod_3D2Osym
    use mod_operator_I3D2OMod_3D2O
    use mod_operator_I3D4O3T_3D4O3sym
    use mod_operator_I3D4O3TMod_3D4O3sym
    use mod_operator_I3D4O4T_3D4O3sym
    use mod_operator_I3D4O4TMod_3D4O3sym
    use mod_operator_3D2Osym_3D4O3sym

    implicit None

    public :: ten_3D2O
    public :: ten_3D2Osym
    public :: ten_3D4O3sym
    public :: iden_3D2Omod
    public :: iden_3D2O
    public :: iden_3D4O3TMod
    public :: iden_3D4O3T
    public :: iden_3D4O4TMod
    public :: iden_3D4O4T

    private 
    ! INCLUDE 'ten_3D4O3sym_3D2Osym/ten_3D4O3sym_3D2Osym.inc'

    ! iden_3D4O3T ! \delta_ij\delta_kl

    ! iden_3D4O3TMod ! val*\delta_ij\delta_kl

    ! type, public :: iden_3D4O4T ! (\delta_ik\delta_jl + \delta_il\delta_jk)/2
    ! end type iden_3D4O4T

    ! iden_3D4O4TMod ! val*(\delta_ik\delta_jl + \delta_il\delta_jk)/2
    

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
    
    public :: assignment (=)

    ! interface ! iden_3D4O4T
    !     module pure function sum_I3D4O4T_3D4O3sym(I4, a) result(res)
    !         implicit none
    !         class(iden_3D4O4T), intent(in) :: I4
    !         type(ten_3D4O3sym), intent(in) :: a
    !         type(ten_3D4O3sym) :: res
    !     end function sum_I3D4O4T_3D4O3sym

    !     module pure function sum_3D4O3sym_I3D4O4T(a, I4) result(res)
    !         implicit none
    !         class(iden_3D4O4T), intent(in) :: I4
    !         type(ten_3D4O3sym), intent(in) :: a
    !         type(ten_3D4O3sym) :: res
    !     end function sum_3D4O3sym_I3D4O4T

    !     module pure function sub_I3D4O4T_3D4O3sym(I4, a) result(res)
    !         implicit none
    !         class(iden_3D4O4T), intent(in) :: I4
    !         type(ten_3D4O3sym), intent(in) :: a
    !         type(ten_3D4O3sym) :: res
    !     end function sub_I3D4O4T_3D4O3sym

    !     module pure function sub_3D4O3sym_I3D4O4T(a, I4) result(res)
    !         implicit none
    !         class(iden_3D4O4T), intent(in) :: I4
    !         type(ten_3D4O3sym), intent(in) :: a
    !         type(ten_3D4O3sym) :: res
    !     end function sub_3D4O3sym_I3D4O4T

    !     module pure function mul_I3D4O4T_real64(I4, a) result(res)
    !         implicit none
    !         class(iden_3D4O4T), intent(in) :: I4
    !         real(real64), intent(in) :: a
    !         type(ten_3D4O3sym) :: res
    !     end function mul_I3D4O4T_real64

    !     module pure function mul_real64_I3D4O4T(a, I4) result(res)
    !         implicit none
    !         class(iden_3D4O4T), intent(in) :: I4
    !         real(real64), intent(in) :: a
    !         type(ten_3D4O3sym) :: res
    !     end function mul_real64_I3D4O4T

    !     module pure function div_I3D4O4T_real64(I4, a) result(res)
    !         implicit none
    !         class(iden_3D4O4T), intent(in) :: I4
    !         real(real64), intent(in) :: a
    !         type(ten_3D4O3sym) :: res
    !     end function div_I3D4O4T_real64
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