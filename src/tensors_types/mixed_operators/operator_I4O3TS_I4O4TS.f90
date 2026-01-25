module mod_operator_I4O3TS_I4O4TS
    use, intrinsic :: iso_fortran_env
    use mod_iden_4O3TS
    use mod_iden_4O4TS
    use mod_ten_3D4O3sym
    implicit none
    private
    
    public :: operator(+)
    interface operator (+)
        module procedure sum_I4O3TS_I4O4TS
        module procedure sum_I4O4TS_I4O3TS
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_I4O3TS_I4O4TS
        module procedure sub_I4O4TS_I4O3TS
    end interface
    
    contains

    pure function sum_I4O3TS_I4O4TS(I3, I4) result(res)
        implicit none
        class(iden_4O3TS), intent(in) :: I3
        class(iden_4O4TS), intent(in) :: I4
        type(ten_3D4O3sym) :: res
        res%vals = (/I4%val+I3%val, I4%val+I3%val, I4%val+I3%val, 0.5D0*I4%val, 0.5D0*I4%val, 0.5D0*I4%val, &
                     I3%val, I3%val, 0D0, 0D0, 0D0, I3%val, &
                     0D0, 0D0, 0D0, 0D0, 0D0, 0D0, &
                     0D0, 0D0, 0D0/)
        !              1    2    3      4      5      6    
        !              7    8    9     10     11     12
        !             13   14   15     16     17     18
        !             19   20   21
    end function sum_I4O3TS_I4O4TS

    pure function sum_I4O4TS_I4O3TS(I4, I3) result(res)
        implicit none
        class(iden_4O3TS), intent(in) :: I3
        class(iden_4O4TS), intent(in) :: I4
        type(ten_3D4O3sym) :: res
        res%vals = (/I4%val+I3%val, I4%val+I3%val, I4%val+I3%val, 0.5D0*I4%val, 0.5D0*I4%val, 0.5D0*I4%val, &
                     I3%val, I3%val, 0D0, 0D0, 0D0, I3%val, &
                     0D0, 0D0, 0D0, 0D0, 0D0, 0D0, &
                     0D0, 0D0, 0D0/)
        !              1    2    3      4      5      6    
        !              7    8    9     10     11     12
        !             13   14   15     16     17     18
        !             19   20   21
    end function sum_I4O4TS_I4O3TS

    pure function sub_I4O3TS_I4O4TS(I3, I4) result(res)
        implicit none
        class(iden_4O3TS), intent(in) :: I3
        class(iden_4O4TS), intent(in) :: I4
        type(ten_3D4O3sym) :: res
        res%vals = (/I3%val-I4%val, I3%val-I4%val, I3%val-I4%val, -0.5D0*I4%val, -0.5D0*I4%val, -0.5D0*I4%val, &
                     I3%val, I3%val, 0D0, 0D0, 0D0, I3%val, &
                     0D0, 0D0, 0D0, 0D0, 0D0, 0D0, &
                     0D0, 0D0, 0D0/)
        !              1    2    3      4      5      6    
        !              7    8    9     10     11     12
        !             13   14   15     16     17     18
        !             19   20   21
    end function sub_I4O3TS_I4O4TS

    pure function sub_I4O4TS_I4O3TS(I4, I3) result(res)
        implicit none
        class(iden_4O3TS), intent(in) :: I3
        class(iden_4O4TS), intent(in) :: I4
        type(ten_3D4O3sym) :: res
        res%vals = (/-I3%val+I4%val, -I3%val+I4%val, -I3%val+I4%val, 0.5D0*I4%val, 0.5D0*I4%val, 0.5D0*I4%val, &
                     -I3%val, -I3%val, 0D0, 0D0, 0D0, -I3%val, &
                     0D0, 0D0, 0D0, 0D0, 0D0, 0D0, &
                     0D0, 0D0, 0D0/)
        !               1     2    3      4      5      6    
        !               7     8    9     10     11     12
        !              13    14   15     16     17     18
        !              19    20   21
    end function sub_I4O4TS_I4O3TS
end module mod_operator_I4O3TS_I4O4TS