module mod_operator_I4O4T_3D2Osym
    use, intrinsic :: iso_fortran_env
    use mod_iden_4O4T
    use mod_ten_3D2Osym
    implicit none
    private
    
    ! NOTE: The .ddot. operator it is not necessary for the operations in I:Tensor
    ! public :: operator(.ddot.)
    ! interface operator (.ddot.)
    !     module procedure ddot_3D2Osym_I4O4T
    !     module procedure ddot_I4O4T_3D2Osym
    ! end interface

    contains


    ! pure function ddot_3D2Osym_I4O4T(I2, a) result(res)
    !     implicit none
    !     class(iden_4O4T), intent(in) :: I2
    !     type(ten_3D2Osym), intent(in) :: a
    !     type(ten_3D2Osym) :: res
    !     res%vals = a%vals
    ! end function ddot_3D2Osym_I4O4T

    ! pure function ddot_I4O4T_3D2Osym(a, I2) result(res)
    !     implicit none
    !     class(iden_4O4T), intent(in) :: I2
    !     type(ten_3D2Osym), intent(in) :: a
    !     type(ten_3D2Osym) :: res
    !     res%vals = a%vals
    ! end function ddot_I4O4T_3D2Osym

end module mod_operator_I4O4T_3D2Osym