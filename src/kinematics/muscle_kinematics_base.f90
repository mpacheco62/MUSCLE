! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_kinematics_base
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    implicit none
    private

    public :: Base_kinematics

    type, abstract :: Base_kinematics
        !! Abstract Base Type for Kinematic Strategies
    contains
        ! --- 1. Non-deferred update procedures with fallback implementations ---
        procedure :: update_F           => base_update_F_noop
        procedure :: update_strain      => base_update_strain_noop
        procedure :: update_corotational => base_update_corotational_noop

        ! --- 2. Kinematic Quantities ---
        procedure(jacobian_interface), deferred :: jacobian
        procedure(is_finite_interface), deferred :: is_finite_strain

        ! --- 3. Strain Measures ---
        procedure(strain_interface), deferred :: get_strain

        ! --- 4. Stress Transformations ---
        procedure(push_stress_interface), deferred :: push_forward_stress
        procedure(pull_stress_interface), deferred :: pull_back_stress

        ! --- 5. Tangent Transformations ---
        procedure(push_tangent_interface), deferred :: push_forward_tangent
        procedure(pull_tangent_interface), deferred :: pull_back_tangent

        ! --- 6. Corotational Frame Rotations ---
        procedure :: rotate_to_corotational   => base_rotate_to_corotational
        procedure :: rotate_from_corotational => base_rotate_from_corotational

        ! --- Generic Update Interface (Unambiguous) ---
        generic :: update => update_F, update_strain, update_corotational

        procedure(cauchy_interface), deferred :: to_Cauchy
            !! Context-aware converter: Converts the constitutive flow stress 
            !! to spatial Cauchy stress according to the formulation rules.
    end type Base_kinematics

    abstract interface
        pure function jacobian_interface(self) result(J)
            use, intrinsic :: iso_fortran_env, only : real64
            import Base_kinematics
            class(Base_kinematics), intent(in) :: self
            real(real64)                        :: J
        end function jacobian_interface

        pure function is_finite_interface(self) result(res)
            import Base_kinematics
            class(Base_kinematics), intent(in) :: self
            logical                            :: res
        end function is_finite_interface

        pure function strain_interface(self) result(strain)
            use muscle_tensors, only : ten_3D2Osym
            import Base_kinematics
            class(Base_kinematics), intent(in) :: self
            type(ten_3D2Osym)                  :: strain
        end function strain_interface

        pure function push_stress_interface(self, S_material) result(sigma_spatial)
            use muscle_tensors, only : ten_3D2Osym
            import Base_kinematics
            class(Base_kinematics), intent(in) :: self
            type(ten_3D2Osym), intent(in)      :: S_material
            type(ten_3D2Osym)                  :: sigma_spatial
        end function push_stress_interface

        pure function pull_stress_interface(self, sigma_spatial) result(S_material)
            use muscle_tensors, only : ten_3D2Osym
            import Base_kinematics
            class(Base_kinematics), intent(in) :: self
            type(ten_3D2Osym), intent(in)      :: sigma_spatial
            type(ten_3D2Osym)                  :: S_material
        end function pull_stress_interface

        pure function push_tangent_interface(self, C_material, stress_spatial) result(c_spatial)
            use muscle_tensors, only : ten_3D4O2sym, ten_3D2Osym
            import Base_kinematics
            class(Base_kinematics), intent(in) :: self
            type(ten_3D4O2sym), intent(in)     :: C_material
            type(ten_3D2Osym), intent(in), optional :: stress_spatial
            type(ten_3D4O2sym)                 :: c_spatial
        end function push_tangent_interface

        pure function pull_tangent_interface(self, c_spatial, stress_material) result(C_material)
            use muscle_tensors, only : ten_3D4O2sym, ten_3D2Osym
            import Base_kinematics
            class(Base_kinematics), intent(in) :: self
            type(ten_3D4O2sym), intent(in)     :: c_spatial
            type(ten_3D2Osym), intent(in), optional :: stress_material
            type(ten_3D4O2sym)                 :: C_material
        end function pull_tangent_interface

        pure function cauchy_interface(self, constitutive_stress) result(sigma_cauchy)
            use muscle_tensors, only : ten_3D2Osym
            import Base_kinematics
            class(Base_kinematics), intent(in) :: self
            type(ten_3D2Osym), intent(in)      :: constitutive_stress
                !! Stress outputted by the return mapping algorithm.
            type(ten_3D2Osym)                  :: sigma_cauchy
                !! True spatial Cauchy stress tensor required by FEA solvers.
        end function cauchy_interface
    end interface

contains

    pure subroutine base_update_F_noop(self, F)
        class(Base_kinematics), intent(inout) :: self
        type(ten_3D2O), intent(in)            :: F
    end subroutine base_update_F_noop

    pure subroutine base_update_strain_noop(self, strain)
        class(Base_kinematics), intent(inout) :: self
        type(ten_3D2Osym), intent(in)         :: strain
    end subroutine base_update_strain_noop

    pure subroutine base_update_corotational_noop(self, strain, R)
        class(Base_kinematics), intent(inout) :: self
        type(ten_3D2Osym), intent(in)         :: strain
        type(ten_3D2O), intent(in)            :: R
    end subroutine base_update_corotational_noop

    pure function base_rotate_to_corotational(self, T, R) result(res)
        class(Base_kinematics), intent(in) :: self
        type(ten_3D2Osym), intent(in)      :: T
        type(ten_3D2O), intent(in), optional :: R
        type(ten_3D2Osym)                  :: res

        if (present(R)) then
            res = R%transpose() .transform. T
        else
            res = T
        end if
    end function base_rotate_to_corotational

    pure function base_rotate_from_corotational(self, T_corot, R) result(res)
        class(Base_kinematics), intent(in) :: self
        type(ten_3D2Osym), intent(in)      :: T_corot
        type(ten_3D2O), intent(in), optional :: R
        type(ten_3D2Osym)                  :: res

        if (present(R)) then
            res = R .transform. T_corot
        else
            res = T_corot
        end if
    end function base_rotate_from_corotational

end module muscle_kinematics_base