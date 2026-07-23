! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_kin_spatial_logarithmic
    !! Module muscle_kin_spatial_logarithmic
    !! =====================================
    !! Kinematic strategy using Spatial Hencky strain h = ln(V) = 0.5*ln(b).
    !! Conjugate stress is spatial Kirchhoff stress tau.

    use muscle_tensors
    use muscle_kin_finite_base
    implicit none
    private

    public :: Spatial_logarithmic_kinematics

    type, extends(Base_F_kinematics) :: Spatial_logarithmic_kinematics
        type(ten_3D2Osym) :: strain !! Cached spatial Hencky strain h = ln(V)
    contains
        procedure :: update_F   => spat_log_update_F
        procedure :: get_strain => spat_log_get_strain
        procedure :: to_Cauchy  => spat_log_to_Cauchy
    end type Spatial_logarithmic_kinematics

contains

    pure subroutine spat_log_update_F(self, F)
        class(Spatial_logarithmic_kinematics), intent(inout) :: self
        type(ten_3D2O), intent(in)                           :: F

        call self%update_F_base(F)
        self%strain = self%get_Spatial_Logarithmic()
    end subroutine spat_log_update_F

    pure function spat_log_get_strain(self) result(strain)
        class(Spatial_logarithmic_kinematics), intent(in) :: self
        type(ten_3D2Osym)                                 :: strain

        strain = self%strain
    end function spat_log_get_strain

    pure function spat_log_to_Cauchy(self, constitutive_stress) result(sigma_cauchy)
        !! Converts spatial Kirchhoff stress (tau) to Cauchy stress:
        !! sigma = (1/J) * tau
        class(Spatial_logarithmic_kinematics), intent(in) :: self
        type(ten_3D2Osym), intent(in)                     :: constitutive_stress
        type(ten_3D2Osym)                                 :: sigma_cauchy

        sigma_cauchy = constitutive_stress / self%J
    end function spat_log_to_Cauchy

end module muscle_kin_spatial_logarithmic