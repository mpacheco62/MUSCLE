module mod_yield_criteria
    !! Module mod_yield_criteria
    !! =========================
    !!
    !! Defines the abstract base type for yield criteria used in plasticity models.
    !!
    !! This module provides an abstract derived type, `Base_yield_critera`, which serves as
    !! a blueprint for implementing various yield criteria (e.g., Von Mises, Tresca, Drucker-Prager).
    !! It defines the essential interface that any concrete yield criterion must provide:
    !! a procedure to calculate the equivalent stress.
    !!
    !! Additionally, it provides default implementations for calculating the first and second
    !! derivatives of the equivalent stress with respect to the stress tensor components using
    !! numerical differentiation (finite differences). Concrete implementations can override
    !! these derivative procedures with analytical solutions if available for better performance
    !! and accuracy.
    !!
    !! Public Entities
    !! ---------------
    !!
    !! ### Abstract Derived Type:
    !!
    !! - `Base_yield_critera`: Abstract base type for yield criteria.
    !!     - Deferred Procedure: `stress_eq(stress)` - Interface for a function that computes
    !!       the equivalent stress `res` (`real(real64)`) from a given second-order symmetric
    !!       stress tensor `stress` (`ten_3D2Osym`). **Must be implemented by concrete subtypes.**
    !!     - Procedure: `dstressEq_dstress(stress)` - Computes the first derivative of the
    !!       equivalent stress with respect to the stress tensor (\( \frac{\partial \sigma_{eq}}{\partial \sigma} \)),
    !!       returning a second-order symmetric tensor (`ten_3D2Osym`). By default, binds to
    !!       `dstressEq_dstress_numeric` which uses finite differences. Can be overridden by
    !!       concrete subtypes with an analytical implementation.
    !!     - Procedure: `ddstressEq_ddstress(stress)` - Computes the second derivative of the
    !!       equivalent stress with respect to the stress tensor (\( \frac{\partial^2 \sigma_{eq}}{\partial \sigma \partial \sigma} \)),
    !!       returning a fourth-order fully symmetric tensor (`ten_3D4O3sym`). By default, binds to
    !!       `ddstressEq_ddstress_numeric` which uses finite differences. Can be overridden by
    !!       concrete subtypes with an analytical implementation.
    !!
    !! Usage (Conceptual)
    !! ------------------
    !!
    !! Concrete yield criteria will extend `Base_yield_critera`:
    !!
    !! ```fortran
    !! module mod_von_mises_yield
    !!   use mod_yield_criteria
    !!   use tensors_types
    !!   use iso_fortran_env, only: real64
    !!   implicit none
    !!
    !!   type, extends(Base_yield_critera) :: VonMisesYield
    !!     ! No specific parameters needed for standard Von Mises
    !!   contains
    !!     procedure :: stress_eq => calculate_von_mises_eq_stress
    !!     ! Optionally override derivatives for analytical versions:
    !!     ! procedure :: dstressEq_dstress => calculate_von_mises_dstressEq_dstress
    !!     ! procedure :: ddstressEq_ddstress => calculate_von_mises_ddstressEq_ddstress
    !!   end type VonMisesYield
    !!
    !! contains
    !!
    !!   pure function calculate_von_mises_eq_stress(self, stress) result(res)
    !!     class(VonMisesYield), intent(in) :: self
    !!     class(ten_3D2Osym), intent(in) :: stress
    !!     real(real64) :: res
    !!     type(ten_3D2Osym) :: s_dev
    !!     s_dev = .dev. stress
    !!     res = sqrt( (3.0D0/2.0D0) * (s_dev .ddot. s_dev) )
    !!   end function calculate_von_mises_eq_stress
    !!
    !!   ! ... (Optional analytical derivative implementations) ...
    !!
    !! end module mod_von_mises_yield
    !!
    !! program use_yield_criterion
    !!   use mod_von_mises_yield
    !!   use tensors_types
    !!   use iso_fortran_env, only: real64
    !!   implicit none
    !!
    !!   type(VonMisesYield) :: material_yield
    !!   class(Base_yield_critera), allocatable :: yield_ptr
    !!   type(ten_3D2Osym) :: stress_tensor, first_derivative
    !!   type(ten_3D4O3sym) :: second_derivative
    !!   real(real64) :: eq_stress
    !!
    !!   ! Point the abstract pointer to the concrete implementation
    !!   allocate(yield_ptr, source=material_yield)
    !!
    !!   ! Define a stress state
    !!   call stress_tensor%init(xx=100.0, yy=50.0, zz=0.0, xy=20.0, yz=0.0, xz=0.0)
    !!
    !!   ! Use the abstract interface to call the concrete/default methods
    !!   eq_stress = yield_ptr%stress_eq(stress_tensor)
    !!   first_derivative = yield_ptr%dstressEq_dstress(stress_tensor) ! Uses numerical default if not overridden
    !!   second_derivative = yield_ptr%ddstressEq_ddstress(stress_tensor) ! Uses numerical default if not overridden
    !!
    !!   print *, "Stress State (xx, yy, zz, xy):", stress_tensor%vals(1:4)
    !!   print *, "Equivalent Stress (Von Mises):", eq_stress
    !!   print *, "d(EqStress)/d(Stress) (xx component):", first_derivative%xx()
    !!   print *, "d2(EqStress)/d(Stress)2 (1111 component):", second_derivative%vals(1)
    !!
    !!   deallocate(yield_ptr)
    !!
    !! end program use_yield_criterion
    !! ```
    !!
    !! For tensor type definitions see [[tensors_types]].

    use, intrinsic :: iso_fortran_env
    use tensors_types
    implicit None
    PRIVATE

    public :: Base_yield_critera
    type, abstract :: Base_yield_critera
        !! Abstract Base Type for Yield Criteria
        !! =====================================
        !!
        !! Serves as the fundamental interface for all yield criteria within the library.
        !! Any concrete yield criterion (e.g., Von Mises, Tresca) must extend this type
        !! and provide an implementation for the deferred procedure `stress_eq`.
        !! Default numerical implementations for the first and second derivatives are provided
        !! via `dstressEq_dstress_numeric` and `ddstressEq_ddstress_numeric`, which can be
        !! overridden by analytical versions in concrete subtypes for efficiency and accuracy.
        contains
            procedure(stress_eq_interface), deferred :: stress_eq
                !! Computes the equivalent stress for a given stress tensor. Must be implemented by concrete types.
            procedure :: dstressEq_dstress => dstressEq_dstress_numeric
                !! Computes d(stress_eq)/d(stress). Default uses numerical differentiation.
            procedure :: ddstressEq_ddstress => ddstressEq_ddstress_numeric
                !! Computes d^2(stress_eq)/d(stress)^2. Default uses numerical differentiation.
            procedure :: dstressEq_dstress_numeric
                !! Numerical implementation of d(stress_eq)/d(stress) using finite differences.
            procedure :: ddstressEq_ddstress_numeric
                !! Numerical implementation of d^2(stress_eq)/d(stress)^2 using finite differences.
            procedure, private :: ddf_dxdy_numeric_vals
                !! Helper subroutine for numerical second derivative calculation.
    end type Base_yield_critera

    interface
        pure function stress_eq_interface(self, stress) result(res)
            !! Interface required for the `stress_eq` procedure.
            use, intrinsic :: iso_fortran_env
            use tensors_types
            import Base_yield_critera
            class(Base_yield_critera), intent(in) :: self  !! The yield criterion object.
            class(ten_3D2Osym), intent(in) :: stress       !! Input stress tensor (`ten_3D2Osym`).
            real(real64) :: res                            !! Output equivalent stress (`real(real64)`).
        end function stress_eq_interface

        pure function dstressEq_dstress_interface(self, stress) result(res)
            !! Interface required for the `dstressEq_dstress` procedure (analytical version).
            use, intrinsic :: iso_fortran_env
            use tensors_types
            import Base_yield_critera
            class(Base_yield_critera), intent(in) :: self  !! The yield criterion object.
            class(ten_3D2Osym), intent(in) :: stress       !! Input stress tensor (`ten_3D2Osym`) at which the derivative is evaluated.
            type(ten_3D2Osym) :: res                       !! Output first derivative tensor (`ten_3D2Osym`).
        end function dstressEq_dstress_interface

        pure function ddstressEq_ddstress_interface(self, stress) result(res)
            !! Interface required for the `dstressEq_dstress` procedure (analytical version).
            use, intrinsic :: iso_fortran_env
            use tensors_types
            import Base_yield_critera
            class(Base_yield_critera), intent(in) :: self  !! The yield criterion object.
            class(ten_3D2Osym), intent(in) :: stress       !! Input stress tensor (`ten_3D2Osym`) at which the second derivative is evaluated.  
            type(ten_3D4O3sym) :: res                      !! Output second derivative tensor (`ten_3D4O3sym`).
        end function ddstressEq_ddstress_interface


    end interface


    contains
        pure function dstressEq_dstress_numeric(self, stress) result(res)
            !! Computes d(stress_eq)/d(stress) numerically using central finite differences.
            !! This is the default implementation bound to the `dstressEq_dstress` procedure.
            use, intrinsic :: iso_fortran_env
            class(Base_yield_critera), intent(in) :: self  !! The yield criterion object.
            class(ten_3D2Osym), intent(in) :: stress       !! Input stress tensor (`ten_3D2Osym`) at which the derivative is evaluated.
            type(ten_3D2Osym) :: res                       !! Output first derivative tensor (`ten_3D2Osym`).

            type(ten_3D2Osym) :: stress_var1, stress_var2
            real(real64) :: stress_eq
            real(real64), parameter :: DIVEPS = 1D-7, MAX_EPS=1D-40  !! Relative and minimum absolute step size
            real(real64) :: eps

            res = 0.0D0
            stress_eq = self%stress_eq(stress)

            ! Use a relative step size, but ensure it's not too small
            eps = max(abs(stress_eq*DIVEPS), MAX_EPS)

            ! d/d(sigma_11)
            stress_var1 = stress; stress_var1%vals(1) = stress_var1%vals(1) + eps
            stress_var2 = stress; stress_var2%vals(1) = stress_var2%vals(1) - eps
            res%vals(1) = (self%stress_eq(stress_var1) - self%stress_eq(stress_var2))/(2D0*eps)

            ! d/d(sigma_22)
            stress_var1%vals(1) = stress_var1%vals(1) - eps; stress_var1%vals(2) = stress_var1%vals(2) + eps
            stress_var2%vals(1) = stress_var2%vals(1) + eps; stress_var2%vals(2) = stress_var2%vals(2) - eps
            res%vals(2) = (self%stress_eq(stress_var1) - self%stress_eq(stress_var2))/(2D0*eps)

            ! d/d(sigma_33)
            stress_var1%vals(2) = stress_var1%vals(2) - eps; stress_var1%vals(3) = stress_var1%vals(3) + eps
            stress_var2%vals(2) = stress_var2%vals(2) + eps; stress_var2%vals(3) = stress_var2%vals(3) - eps
            res%vals(3) = (self%stress_eq(stress_var1) - self%stress_eq(stress_var2))/(2D0*eps)

            ! d/d(sigma_12)
            stress_var1%vals(3) = stress_var1%vals(3) - eps; stress_var1%vals(4) = stress_var1%vals(4) + eps/2D0
            stress_var2%vals(3) = stress_var2%vals(3) + eps; stress_var2%vals(4) = stress_var2%vals(4) - eps/2D0
            res%vals(4) = (self%stress_eq(stress_var1) - self%stress_eq(stress_var2))/(2D0*eps)

            ! d/d(sigma_23)
            stress_var1%vals(4) = stress_var1%vals(4) - eps/2D0; stress_var1%vals(5) = stress_var1%vals(5) + eps/2D0
            stress_var2%vals(4) = stress_var2%vals(4) + eps/2D0; stress_var2%vals(5) = stress_var2%vals(5) - eps/2D0
            res%vals(5) = (self%stress_eq(stress_var1) - self%stress_eq(stress_var2))/(2D0*eps)

            ! d/d(sigma_13)
            stress_var1%vals(5) = stress_var1%vals(5) - eps/2D0; stress_var1%vals(6) = stress_var1%vals(6) + eps/2D0
            stress_var2%vals(5) = stress_var2%vals(5) + eps/2D0; stress_var2%vals(6) = stress_var2%vals(6) - eps/2D0
            res%vals(6) = (self%stress_eq(stress_var1) - self%stress_eq(stress_var2))/(2D0*eps)

        end function dstressEq_dstress_numeric

        pure function ddstressEq_ddstress_numeric(self, stress) result(res)
            !! Computes d^2(stress_eq)/d(stress)^2 numerically using central finite differences.
            !! This is the default implementation bound to the `ddstressEq_ddstress` procedure.
            use, intrinsic :: iso_fortran_env
            implicit none
            class(Base_yield_critera), intent(in) :: self  !! The yield criterion object.
            class(ten_3D2Osym), intent(in) :: stress       !! Input stress tensor (`ten_3D2Osym`) at which the second derivative is evaluated.
            type(ten_3D4O3sym) :: res                      !! Output second derivative tensor (`ten_3D4O3sym`).

            type(ten_3D2Osym) :: stress_var1
            real(real64) :: stress_eq
            real(real64), parameter :: DIVEPS = 1D-4, MAX_EPS=1D-40  !! Relative and minimum absolute step size for 2nd derivative
            real(real64) :: eps

            real(real64) :: f11p_11p, f11m_11m
            real(real64) :: f22p_22p, f22m_22m
            real(real64) :: f33p_33p, f33m_33m
            real(real64) :: f12p_12p, f12m_12m
            real(real64) :: f23p_23p, f23m_23m
            real(real64) :: f13p_13p, f13m_13m

            real(real64) :: f11p_22p, f11p_22m, f11m_22p, f11m_22m
            real(real64) :: f11p_33p, f11p_33m, f11m_33p, f11m_33m
            real(real64) :: f11p_12p, f11p_12m, f11m_12p, f11m_12m
            real(real64) :: f11p_23p, f11p_23m, f11m_23p, f11m_23m
            real(real64) :: f11p_13p, f11p_13m, f11m_13p, f11m_13m

            real(real64) :: f22p_11p, f22p_11m, f22m_11p, f22m_11m
            real(real64) :: f22p_33p, f22p_33m, f22m_33p, f22m_33m
            real(real64) :: f22p_12p, f22p_12m, f22m_12p, f22m_12m
            real(real64) :: f22p_23p, f22p_23m, f22m_23p, f22m_23m
            real(real64) :: f22p_13p, f22p_13m, f22m_13p, f22m_13m

            real(real64) :: f33p_11p, f33p_11m, f33m_11p, f33m_11m
            real(real64) :: f33p_22p, f33p_22m, f33m_22p, f33m_22m
            real(real64) :: f33p_12p, f33p_12m, f33m_12p, f33m_12m
            real(real64) :: f33p_23p, f33p_23m, f33m_23p, f33m_23m
            real(real64) :: f33p_13p, f33p_13m, f33m_13p, f33m_13m

            real(real64) :: f12p_11p, f12p_11m, f12m_11p, f12m_11m
            real(real64) :: f12p_22p, f12p_22m, f12m_22p, f12m_22m
            real(real64) :: f12p_33p, f12p_33m, f12m_33p, f12m_33m
            real(real64) :: f12p_23p, f12p_23m, f12m_23p, f12m_23m
            real(real64) :: f12p_13p, f12p_13m, f12m_13p, f12m_13m

            real(real64) :: f23p_11p, f23p_11m, f23m_11p, f23m_11m
            real(real64) :: f23p_22p, f23p_22m, f23m_22p, f23m_22m
            real(real64) :: f23p_33p, f23p_33m, f23m_33p, f23m_33m
            real(real64) :: f23p_12p, f23p_12m, f23m_12p, f23m_12m
            real(real64) :: f23p_13p, f23p_13m, f23m_13p, f23m_13m

            real(real64) :: f13p_11p, f13p_11m, f13m_11p, f13m_11m
            real(real64) :: f13p_22p, f13p_22m, f13m_22p, f13m_22m
            real(real64) :: f13p_33p, f13p_33m, f13m_33p, f13m_33m
            real(real64) :: f13p_12p, f13p_12m, f13m_12p, f13m_12m
            real(real64) :: f13p_23p, f13p_23m, f13m_23p, f13m_23m

            ! res = 0.0D0
            stress_eq = self%stress_eq(stress)


            eps = max(abs(stress_eq*DIVEPS), MAX_EPS)

            stress_var1 = stress

            stress_var1%vals(1) = stress%vals(1) + 2*eps
            f11p_11p = self%stress_eq(stress_var1)
            stress_var1%vals(1) = stress%vals(1) - 2*eps
            f11m_11m = self%stress_eq(stress_var1)
            stress_var1%vals(1) = stress%vals(1)

            stress_var1%vals(2) = stress%vals(2) + 2*eps
            f22p_22p = self%stress_eq(stress_var1)
            stress_var1%vals(2) = stress%vals(2) - 2*eps
            f22m_22m = self%stress_eq(stress_var1)
            stress_var1%vals(2) = stress%vals(2)

            stress_var1%vals(3) = stress%vals(3) + 2*eps
            f33p_33p = self%stress_eq(stress_var1)
            stress_var1%vals(3) = stress%vals(3) - 2*eps
            f33m_33m = self%stress_eq(stress_var1)
            stress_var1%vals(3) = stress%vals(3)

            stress_var1%vals(4) = stress%vals(4) + eps
            f12p_12p = self%stress_eq(stress_var1)
            stress_var1%vals(4) = stress%vals(4) - eps
            f12m_12m = self%stress_eq(stress_var1)
            stress_var1%vals(4) = stress%vals(4)

            stress_var1%vals(5) = stress%vals(5) + eps
            f23p_23p = self%stress_eq(stress_var1)
            stress_var1%vals(5) = stress%vals(5) - eps
            f23m_23m = self%stress_eq(stress_var1)
            stress_var1%vals(5) = stress%vals(5)

            stress_var1%vals(6) = stress%vals(6) + eps
            f13p_13p = self%stress_eq(stress_var1)
            stress_var1%vals(6) = stress%vals(6) - eps
            f13m_13m = self%stress_eq(stress_var1)
            stress_var1%vals(6) = stress%vals(6)

            call self%ddf_dxdy_numeric_vals(1, 2,     eps,     eps, f11p_22p, f11m_22p, f11m_22m, f11p_22m, stress, stress_var1)
            call self%ddf_dxdy_numeric_vals(1, 3,     eps,     eps, f11p_33p, f11m_33p, f11m_33m, f11p_33m, stress, stress_var1)
            call self%ddf_dxdy_numeric_vals(1, 4,     eps, eps/2D0, f11p_12p, f11m_12p, f11m_12m, f11p_12m, stress, stress_var1)
            call self%ddf_dxdy_numeric_vals(1, 5,     eps, eps/2D0, f11p_23p, f11m_23p, f11m_23m, f11p_23m, stress, stress_var1)
            call self%ddf_dxdy_numeric_vals(1, 6,     eps, eps/2D0, f11p_13p, f11m_13p, f11m_13m, f11p_13m, stress, stress_var1)

            call self%ddf_dxdy_numeric_vals(2, 1,     eps,     eps, f22p_11p, f22m_11p, f22m_11m, f22p_11m, stress, stress_var1)
            call self%ddf_dxdy_numeric_vals(2, 3,     eps,     eps, f22p_33p, f22m_33p, f22m_33m, f22p_33m, stress, stress_var1)
            call self%ddf_dxdy_numeric_vals(2, 4,     eps, eps/2D0, f22p_12p, f22m_12p, f22m_12m, f22p_12m, stress, stress_var1)
            call self%ddf_dxdy_numeric_vals(2, 5,     eps, eps/2D0, f22p_23p, f22m_23p, f22m_23m, f22p_23m, stress, stress_var1)
            call self%ddf_dxdy_numeric_vals(2, 6,     eps, eps/2D0, f22p_13p, f22m_13p, f22m_13m, f22p_13m, stress, stress_var1)

            call self%ddf_dxdy_numeric_vals(3, 1,     eps,     eps, f33p_11p, f33m_11p, f33m_11m, f33p_11m, stress, stress_var1)
            call self%ddf_dxdy_numeric_vals(3, 2,     eps,     eps, f33p_22p, f33m_22p, f33m_22m, f33p_22m, stress, stress_var1)
            call self%ddf_dxdy_numeric_vals(3, 4,     eps, eps/2D0, f33p_12p, f33m_12p, f33m_12m, f33p_12m, stress, stress_var1)
            call self%ddf_dxdy_numeric_vals(3, 5,     eps, eps/2D0, f33p_23p, f33m_23p, f33m_23m, f33p_23m, stress, stress_var1)
            call self%ddf_dxdy_numeric_vals(3, 6,     eps, eps/2D0, f33p_13p, f33m_13p, f33m_13m, f33p_13m, stress, stress_var1)

            call self%ddf_dxdy_numeric_vals(4, 1, eps/2D0,     eps, f12p_11p, f12m_11p, f12m_11m, f12p_11m, stress, stress_var1)
            call self%ddf_dxdy_numeric_vals(4, 2, eps/2D0,     eps, f12p_22p, f12m_22p, f12m_22m, f12p_22m, stress, stress_var1)
            call self%ddf_dxdy_numeric_vals(4, 3, eps/2D0,     eps, f12p_33p, f12m_33p, f12m_33m, f12p_33m, stress, stress_var1)
            call self%ddf_dxdy_numeric_vals(4, 5, eps/2D0, eps/2D0, f12p_23p, f12m_23p, f12m_23m, f12p_23m, stress, stress_var1)
            call self%ddf_dxdy_numeric_vals(4, 6, eps/2D0, eps/2D0, f12p_13p, f12m_13p, f12m_13m, f12p_13m, stress, stress_var1)

            call self%ddf_dxdy_numeric_vals(5, 1, eps/2D0,     eps, f23p_11p, f23m_11p, f23m_11m, f23p_11m, stress, stress_var1)
            call self%ddf_dxdy_numeric_vals(5, 2, eps/2D0,     eps, f23p_22p, f23m_22p, f23m_22m, f23p_22m, stress, stress_var1)
            call self%ddf_dxdy_numeric_vals(5, 3, eps/2D0,     eps, f23p_33p, f23m_33p, f23m_33m, f23p_33m, stress, stress_var1)
            call self%ddf_dxdy_numeric_vals(5, 4, eps/2D0, eps/2D0, f23p_12p, f23m_12p, f23m_12m, f23p_12m, stress, stress_var1)
            call self%ddf_dxdy_numeric_vals(5, 6, eps/2D0, eps/2D0, f23p_13p, f23m_13p, f23m_13m, f23p_13m, stress, stress_var1)

            call self%ddf_dxdy_numeric_vals(6, 1, eps/2D0,     eps, f13p_11p, f13m_11p, f13m_11m, f13p_11m, stress, stress_var1)
            call self%ddf_dxdy_numeric_vals(6, 2, eps/2D0,     eps, f13p_22p, f13m_22p, f13m_22m, f13p_22m, stress, stress_var1)
            call self%ddf_dxdy_numeric_vals(6, 3, eps/2D0,     eps, f13p_33p, f13m_33p, f13m_33m, f13p_33m, stress, stress_var1)
            call self%ddf_dxdy_numeric_vals(6, 4, eps/2D0, eps/2D0, f13p_12p, f13m_12p, f13m_12m, f13p_12m, stress, stress_var1)
            call self%ddf_dxdy_numeric_vals(6, 5, eps/2D0, eps/2D0, f13p_23p, f13m_23p, f13m_23m, f13p_23m, stress, stress_var1)

            ! temp = f11p_11p - 2D0*stress_eq + f11m_11m
            res%vals(1) = (f11p_11p - 2D0*stress_eq + f11m_11m)/(4D0*eps**2)
            res%vals(2) = (f22p_22p - 2D0*stress_eq + f22m_22m)/(4D0*eps**2)
            res%vals(3) = (f33p_33p - 2D0*stress_eq + f33m_33m)/(4D0*eps**2)
            res%vals(4) = (f12p_12p - 2D0*stress_eq + f12m_12m)/(4D0*eps**2)
            res%vals(5) = (f23p_23p - 2D0*stress_eq + f23m_23m)/(4D0*eps**2)
            res%vals(6) = (f13p_13p - 2D0*stress_eq + f13m_13m)/(4D0*eps**2)

            !
            !  | ( 1:1111) ( 7:1122) (12:1133) (16:1112) (19:1123) (21:1113) |
            !  | ( 7:2211) ( 2:2222) ( 8:2233) (13:2212) (17:2223) (20:2213) |
            !  | (12:3311) ( 8:3322) ( 3:3333) ( 9:3312) (14:3323) (18:3313) |
            !  | (16:1211) (13:1222) ( 9:1233) ( 4:1212) (10:1223) (15:1213) |
            !  | (19:2311) (17:2322) (14:2333) (10:2312) ( 5:2323) (11:2313) |
            !  | (21:1311) (20:1322) (18:1333) (15:1312) (11:1323) ( 6:1313) |

            res%vals( 7) = (f11p_22p - f11m_22p -f11p_22m + f11m_22m)/(4D0*eps**2)
            res%vals( 8) = (f22p_33p - f22m_33p -f22p_33m + f22m_33m)/(4D0*eps**2)
            res%vals( 9) = (f33p_12p - f33m_12p -f33p_12m + f33m_12m)/(4D0*eps**2)
            res%vals(10) = (f12p_23p - f12m_23p -f12p_23m + f12m_23m)/(4D0*eps**2)
            res%vals(11) = (f23p_13p - f23m_13p -f23p_13m + f23m_13m)/(4D0*eps**2)
            res%vals(12) = (f11p_33p - f11m_33p -f11p_33m + f11m_33m)/(4D0*eps**2)
            res%vals(13) = (f22p_12p - f22m_12p -f22p_12m + f22m_12m)/(4D0*eps**2)
            res%vals(14) = (f33p_23p - f33m_23p -f33p_23m + f33m_23m)/(4D0*eps**2)
            res%vals(15) = (f12p_13p - f12m_13p -f12p_13m + f12m_13m)/(4D0*eps**2)
            res%vals(16) = (f11p_12p - f11m_12p -f11p_12m + f11m_12m)/(4D0*eps**2)
            res%vals(17) = (f22p_23p - f22m_23p -f22p_23m + f22m_23m)/(4D0*eps**2)
            res%vals(18) = (f33p_13p - f33m_13p -f33p_13m + f33m_13m)/(4D0*eps**2)
            res%vals(19) = (f11p_23p - f11m_23p -f11p_23m + f11m_23m)/(4D0*eps**2)
            res%vals(20) = (f22p_13p - f22m_13p -f22p_13m + f22m_13m)/(4D0*eps**2)
            res%vals(21) = (f11p_13p - f11m_13p -f11p_13m + f11m_13m)/(4D0*eps**2)

        end function ddstressEq_ddstress_numeric

        pure subroutine ddf_dxdy_numeric_vals(self, x_idx, y_idx, epsx, epsy, xp_yp, xm_yp, xm_ym, xp_ym, stress, stress_temp)
            !! Helper subroutine to calculate the four function evaluations needed for the central difference
            !! approximation of the mixed partial derivative d^2f / dx dy.
            use, intrinsic :: iso_fortran_env
            implicit none
            class(Base_yield_critera), intent(in) :: self            !! The yield criterion object.
            integer, intent(in) :: x_idx, y_idx                      !! Indices (1-6) corresponding to the stress components in Voigt notation.
            real(real64), intent(in) :: epsx, epsy                   !! Step sizes for the x and y directions.
            real(real64), intent(out) :: xp_yp, xm_yp, xm_ym, xp_ym  !! Function evaluations at f(x+epsx, y+epsy), f(x-epsx, y+epsy), f(x-epsx, y-epsy), f(x+epsx, y-epsy).
            class(ten_3D2Osym), intent(in) :: stress                 !! The original stress tensor at which the derivative is evaluated.
            class(ten_3D2Osym), intent(inout) :: stress_temp         !! A temporary stress tensor used for perturbations.
    
    
            ! f(x+h, y+k)
            stress_temp%vals(x_idx) = stress%vals(x_idx) + epsx
            stress_temp%vals(y_idx) = stress%vals(y_idx) + epsy
            xp_yp = self%stress_eq(stress_temp)
    
            ! f(x-h, y+k)  
            stress_temp%vals(x_idx) = stress%vals(x_idx) - epsx
            xm_yp = self%stress_eq(stress_temp)
    
            ! f(x-h, y-k)
            stress_temp%vals(y_idx) = stress%vals(y_idx) - epsy
            xm_ym = self%stress_eq(stress_temp)
    
            ! f(x+h, y-k)
            stress_temp%vals(x_idx) = stress%vals(x_idx) + epsx
            xp_ym = self%stress_eq(stress_temp)
            
            ! Reset temporary tensor to original state (optional, but good practice)
            stress_temp%vals(x_idx) = stress%vals(x_idx)
            stress_temp%vals(y_idx) = stress%vals(y_idx)
        end subroutine ddf_dxdy_numeric_vals
end module