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
            use :: muscle_math_derivatives
            class(Base_yield_critera), intent(in) :: self  !! The yield criterion object.
            class(ten_3D2Osym), intent(in) :: stress       !! Input stress tensor (`ten_3D2Osym`) at which the derivative is evaluated.
            type(ten_3D2Osym) :: res                       !! Output first derivative tensor (`ten_3D2Osym`).

            res = derivative(wrapper, stress)

            contains
                pure function wrapper(x1) result(res1)
                    use, intrinsic :: iso_fortran_env
                    use tensors_types, only : ten_3D2Osym
                    implicit none
                    type(ten_3D2Osym), intent(in) :: x1
                    real(real64) :: res1
                    res1 = self%stress_eq(x1)
                end function wrapper

        end function dstressEq_dstress_numeric

        pure function ddstressEq_ddstress_numeric(self, stress) result(res)
            !! Computes d^2(stress_eq)/d(stress)^2 numerically using central finite differences.
            !! This is the default implementation bound to the `ddstressEq_ddstress` procedure.
            use, intrinsic :: iso_fortran_env
            use :: muscle_math_derivatives
            implicit none
            class(Base_yield_critera), intent(in) :: self  !! The yield criterion object.
            class(ten_3D2Osym), intent(in) :: stress       !! Input stress tensor (`ten_3D2Osym`) at which the second derivative is evaluated.
            type(ten_3D4O3sym) :: res                      !! Output second derivative tensor (`ten_3D4O2sym`).
  
            res = derivative2O(wrapper, stress)
            contains
                pure function wrapper(x1) result(res1)
                    use, intrinsic :: iso_fortran_env
                    use tensors_types, only : ten_3D2Osym
                    implicit none
                    type(ten_3D2Osym), intent(in) :: x1
                    real(real64) :: res1
                    res1 = self%stress_eq(x1)
                end function wrapper
        end function ddstressEq_ddstress_numeric
end module