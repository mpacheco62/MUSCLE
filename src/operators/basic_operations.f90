module basic_operations
  use, intrinsic :: iso_fortran_env
  
  public
  private :: eigenvals_3x3sym
  interface eigenvals
    module procedure eigenvals_3x3sym
  end interface


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

      pure function eigenvals_3x3sym(mat)
        ! based on https://doi.org/10.1002/nme.7153
        use tensors_types
        implicit none
        type(ten_3D2Osym), intent(in) :: mat
        real(real64) :: eigenvals_3x3sym(3)

        real(real64), parameter :: EPS=1e-10
        real(real64) :: I1, J2, s, d, alpha
        real(real64) :: J2_sqrt
        real(real64) :: cd, lam_a, lam_b, lam_c, sd
        integer :: sj
        type(ten_3D2Osym) :: Sm, T, d1, d2
        type(iden_3D2O) :: Iden



        I1 = sum(mat%vals(1:3))
        J2 = ((mat%vals(1)-mat%vals(2))**2D0 + (mat%vals(2)-mat%vals(3))**2D0 + (mat%vals(3)-mat%vals(1))**2D0)/6D0 &
             + (mat%vals(4)**2 + mat%vals(5)**2 + mat%vals(6)**2)
        s = (J2/3D0)**0.5D0

        if (abs(s).le.EPS) then
          eigenvals_3x3sym = I1/3D0
          return
        end if
        
        Sm = mat - (I1/3D0)*Iden
        T = Sm%square() - (2D0/3D0)*J2*Iden
        d1 = T-s*Sm
        d2 = T+s*Sm
        d = ((d1 .ddot. d1) / (d2 .ddot. d2))**0.5D0

        J2_sqrt = J2**0.5D0
        if ((1D0-d) .le. EPS) then
          eigenvals_3x3sym(1) = J2_sqrt
          eigenvals_3x3sym(2) = 0
          eigenvals_3x3sym(3) = -eigenvals_3x3sym(1)
          eigenvals_3x3sym = eigenvals_3x3sym + I1/3D0
          return
        end if

        sj = int((1D0-d)/abs(1D0-d))
        alpha = 2D0/3D0 * datan(d**sj)
        cd = sj*s*cos(alpha)
        lam_a = 2D0*cd

        sd = J2_sqrt*sin(alpha)
        lam_b = -cd + sd
        lam_c = -cd - sd

        eigenvals_3x3sym(1) = max(lam_a, max(lam_b, lam_c))
        eigenvals_3x3sym(3) = min(lam_a, min(lam_b, lam_c))
        eigenvals_3x3sym(2) = lam_a + lam_b + lam_c - eigenvals_3x3sym(1) - eigenvals_3x3sym(3)
        eigenvals_3x3sym = eigenvals_3x3sym + I1/3D0

      end function eigenvals_3x3sym
end module