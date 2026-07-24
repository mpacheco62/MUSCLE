! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_kin_material_logarithmic
    !! Module muscle_kin_material_logarithmic
    !! ======================================
    !! Kinematic strategy using Material Hencky strain H = ln(U) = 0.5*ln(C).
    !! Conjugate stress is Rotated Kirchhoff stress tau_R.

    use muscle_tensors
    use muscle_kin_finite_base
    implicit none
    private

    public :: Material_logarithmic_kinematics

    type, extends(Base_F_kinematics) :: Material_logarithmic_kinematics
        type(ten_3D2Osym) :: strain !! Cached Hencky strain H = ln(U)
    contains
        procedure :: update_F   => mat_log_update_F
        procedure :: get_strain => mat_log_get_strain
        procedure :: to_Cauchy  => mat_log_to_Cauchy
        procedure :: to_spatial_tangent => mat_log_to_spatial_tangent
    end type Material_logarithmic_kinematics

contains

    pure subroutine mat_log_update_F(self, F)
        class(Material_logarithmic_kinematics), intent(inout) :: self
        type(ten_3D2O), intent(in)                            :: F

        call self%update_F_base(F)
        self%strain = self%get_Material_Logarithmic()
    end subroutine mat_log_update_F

    pure function mat_log_get_strain(self) result(strain)
        class(Material_logarithmic_kinematics), intent(in) :: self
        type(ten_3D2Osym)                                  :: strain

        strain = self%strain
    end function mat_log_get_strain

    pure function mat_log_to_Cauchy(self, constitutive_stress) result(sigma_cauchy)
        !! Converts Rotated Kirchhoff stress (tau_R) to spatial Cauchy stress:
        !! sigma = (1/J) * R * tau_R * R^T
        class(Material_logarithmic_kinematics), intent(in) :: self
        type(ten_3D2Osym), intent(in)                       :: constitutive_stress
        type(ten_3D2Osym)                                   :: sigma_cauchy

        sigma_cauchy = (self%R_rot .transform. constitutive_stress) / self%J
    end function mat_log_to_Cauchy

    pure function mat_log_to_spatial_tangent(self, C_constitutive, stress_spatial) result(c_spatial)
        !! Converts material logarithmic tangent d(tau_R)/dH to spatial tangent stiffness:
        !! c_spat = (1/J) * R_iI * R_jJ * R_kK * R_lL * C_IJKL
        class(Material_logarithmic_kinematics), intent(in) :: self
        type(ten_3D4O2sym), intent(in)                     :: C_constitutive
        type(ten_3D2Osym), intent(in), optional            :: stress_spatial
        type(ten_3D4O2sym)                                 :: c_spatial

        c_spatial = (self%R_rot .transform. C_constitutive) / self%J
    end function mat_log_to_spatial_tangent

end module muscle_kin_material_logarithmic