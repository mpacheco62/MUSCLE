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

end module muscle_kin_material_logarithmic