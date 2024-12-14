! my_tests.f90
module my_tests_mod
    use iso_fortran_env
    use iso_c_binding
    implicit none
    contains
    
    integer function my_test_success() result(ret) bind(C)
      ret = 0
    end function my_test_success
    
    integer function my_test_fail() result(ret) bind(C)
      ret = 1
    end function my_test_fail
    
end module my_tests_mod