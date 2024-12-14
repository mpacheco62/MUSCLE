module basic_operations
  use, intrinsic :: iso_fortran_env
  
  public


  contains
      ! pure function compare_mat(mat1, mat2, eps)
      !     logical :: compare_mat
      !     real(real64), intent(in), dimension(:, :) :: mat1, mat2
      !     real(real64), intent(in) :: eps
      !     integer, dimension(2) :: dim
      !     integer :: i, j
      !     real(real64) :: error
          
      !     error = 0.0D0
      !     dim = shape(mat1)

      !     do i=1,dim(1)
      !       do j=1,dim(2)
      !         error = error + (mat1(i,j) + mat2(i,j))**2
      !       end do
      !     end do

      !     if (error > eps) then
      !       compare_mat = .FALSE.
      !     else
      !       compare_mat = .TRUE.
      !     end if
      ! end function

      ! pure function compare_vec(vec1, vec2, eps)
      !     logical :: compare_vec
      !     real(real64), intent(in), dimension(:) :: vec1, vec2
      !     real(real64), intent(in) :: eps
      !     integer :: dim, i
      !     real(real64) :: error
          
      !     error = 0.0D0
      !     dim = size(vec1)

      !     do i=1,dim
      !       error = error + (vec1(i) + vec2(i))**2
      !     end do

      !     if (error > eps) then
      !       compare_vec = .FALSE.
      !     else
      !       compare_vec = .TRUE.
      !     end if
      ! end function



        ! function derivative_escalar_tensor3x3(func, type, mat)
        !   implicit none
        !   interface
        !     function f_scalar_3x3(self, x)
        !       use, intrinsic :: iso_fortran_env
        !       implicit none
        !       class(*) :: self
        !       real(real64), intent(in), dimension(3,3) :: x
        !       real(real64) :: f_scalar_3x3
        !     end function f_scalar_3x3
        !   end interface
  
        !   procedure(f_scalar_3x3) :: func
        !   class(*) :: type
        !   real(real64), intent(in) :: mat(3,3)
        !   real(real64) :: derivative_escalar_tensor3x3(3,3)
          
        !   real(real64) :: mat_var(3,3), val_func, val_func_var
        !   real(real64), parameter :: EPS=1.0D-8
        !   integer :: i,j
  
  
        !   derivative_escalar_tensor3x3 = 0.0D0
        !   val_func = func(type, mat)
        !   do i=1,3
        !       do j=1,3
        !           mat_var = mat
        !           mat_var(i,j) = mat_var(i,j) + EPS
        !           val_func_var = func(type, mat_var)
        !           derivative_escalar_tensor3x3(i,j) = (val_func_var - val_func)/EPS
        !       end do
        !   end do
        ! end function

      pure function derivative_escalar_tensor3x3(func, mat)
        implicit none
        interface
          pure function f_scalar_3x3(x)
            use, intrinsic :: iso_fortran_env
            implicit none
            real(real64), intent(in), dimension(3,3) :: x
            real(real64) :: f_scalar_3x3
          end function f_scalar_3x3
        end interface

        procedure(f_scalar_3x3) :: func
        real(real64), intent(in) :: mat(3,3)
        real(real64) :: derivative_escalar_tensor3x3(3,3)
        
        real(real64) :: mat_var(3,3), val_func, val_func_var
        real(real64), parameter :: EPS=1.0D-8
        integer :: i,j


        derivative_escalar_tensor3x3 = 0.0D0
        val_func = func(mat)
        do i=1,3
            do j=1,3
                mat_var = mat
                mat_var(i,j) = mat_var(i,j) + EPS
                val_func_var = func(mat_var)
                derivative_escalar_tensor3x3(i,j) = (val_func_var - val_func)/EPS
            end do
        end do
      end function derivative_escalar_tensor3x3

      pure function derivative_escalar_tensor3x3sym(func, mat)
        implicit none
        interface
          pure function f_scalar_3x3(x)
            use, intrinsic :: iso_fortran_env
            implicit none
            real(real64), intent(in), dimension(3,3) :: x
            real(real64) :: f_scalar_3x3
          end function f_scalar_3x3
        end interface

        procedure(f_scalar_3x3) :: func
        real(real64), intent(in) :: mat(3,3)
        real(real64) :: derivative_escalar_tensor3x3sym(3,3)
        
        real(real64) :: mat_var(3,3), val_func, val_func_var
        real(real64), parameter :: EPS=1.0D-8
        integer :: i,j


        derivative_escalar_tensor3x3sym = 0.0D0
        val_func = func(mat)
        do i=1,3
            do j=1,3
                mat_var = mat
                mat_var(i,j) = mat_var(i,j) + EPS/2.0D0
                mat_var(j,i) = mat_var(j,i) + EPS/2.0D0
                val_func_var = func(mat_var)
                derivative_escalar_tensor3x3sym(i,j) = (val_func_var - val_func)/EPS
            end do
        end do
      end function derivative_escalar_tensor3x3sym

end module