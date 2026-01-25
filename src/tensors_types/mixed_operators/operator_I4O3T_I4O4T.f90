module mod_operator_I4O3T_I4O4T
    use, intrinsic :: iso_fortran_env
    use mod_iden_4O3T
    use mod_iden_4O4T
    use mod_ten_3D4O3sym
    implicit none
    private
    
    public :: operator(+)
    interface operator (+)
        module procedure sum_I4O3T_I4O4T
        module procedure sum_I4O4T_I4O3T
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_I4O3T_I4O4T
        module procedure sub_I4O4T_I4O3T
    end interface
    
    contains

    pure function sum_I4O3T_I4O4T(I3, I4) result(res)
        implicit none
        class(iden_4O3T), intent(in) :: I3
        class(iden_4O4T), intent(in) :: I4
        type(ten_3D4O3sym) :: res
        res%vals = (/2D0, 2D0, 2D0, 0.5D0, 0.5D0, 0.5D0, &
                     1D0, 1D0, 0D0,   0D0,   0D0,   1D0, &
                     0D0, 0D0, 0D0,   0D0,   0D0,   0D0, &
                     0D0, 0D0, 0D0/)
        !              1    2    3      4      5      6    
        !              7    8    9     10     11     12
        !             13   14   15     16     17     18
        !             19   20   21
    end function sum_I4O3T_I4O4T

    pure function sum_I4O4T_I4O3T(I4, I3) result(res)
        implicit none
        class(iden_4O3T), intent(in) :: I3
        class(iden_4O4T), intent(in) :: I4
        type(ten_3D4O3sym) :: res
        res%vals = (/2D0, 2D0, 2D0, 0.5D0, 0.5D0, 0.5D0, &
                     1D0, 1D0, 0D0,   0D0,   0D0,   1D0, &
                     0D0, 0D0, 0D0,   0D0,   0D0,   0D0, &
                     0D0, 0D0, 0D0/)
        !              1    2    3      4      5      6    
        !              7    8    9     10     11     12
        !             13   14   15     16     17     18
        !             19   20   21
    end function sum_I4O4T_I4O3T

    pure function sub_I4O3T_I4O4T(I3, I4) result(res)
        implicit none
        class(iden_4O3T), intent(in) :: I3
        class(iden_4O4T), intent(in) :: I4
        type(ten_3D4O3sym) :: res
        res%vals = (/0D0, 0D0, 0D0,-0.5D0,-0.5D0,-0.5D0, &
                     1D0, 1D0, 0D0,   0D0,   0D0,   1D0, &
                     0D0, 0D0, 0D0,   0D0,   0D0,   0D0, &
                     0D0, 0D0, 0D0/)
        !              1    2    3      4      5      6    
        !              7    8    9     10     11     12
        !             13   14   15     16     17     18
        !             19   20   21
    end function sub_I4O3T_I4O4T

    pure function sub_I4O4T_I4O3T(I4, I3) result(res)
        implicit none
        class(iden_4O3T), intent(in) :: I3
        class(iden_4O4T), intent(in) :: I4
        type(ten_3D4O3sym) :: res
        res%vals = (/ 0D0,  0D0, 0D0, 0.5D0, 0.5D0,  0.5D0, &
                     -1D0, -1D0, 0D0,   0D0,   0D0,   -1D0, &
                      0D0,  0D0, 0D0,   0D0,   0D0,    0D0, &
                      0D0,  0D0, 0D0/)
        !               1     2    3      4      5      6    
        !               7     8    9     10     11     12
        !              13    14   15     16     17     18
        !              19    20   21
    end function sub_I4O4T_I4O3T
end module mod_operator_I4O3T_I4O4T