! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_kin_corotational
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_kinematics_base
    implicit none
    private

    public :: Corotational_kinematics

    type, extends(Base_kinematics) :: Corotational_kinematics
        type(ten_3D2Osym) :: strain             !! Strain in unrotated frame
        type(ten_3D2O)    :: R_rot              !! Rotation matrix R
        real(real64)      :: J = 1.0D0          !! Volume ratio det(F) = exp(tr(eps))
        logical           :: has_rotation = .FALSE.
    contains
        procedure :: update_corotational => corotational_update
        procedure :: jacobian            => corotational_jacobian
        procedure :: is_finite_strain    => corotational_is_finite
        procedure :: get_strain          => corotational_get_strain
        procedure :: push_forward_stress => corotational_push_stress
        procedure :: pull_back_stress     => corotational_pull_stress
        procedure :: push_forward_tangent => corotational_push_tangent
        procedure :: pull_back_tangent   => corotational_pull_tangent
        procedure :: to_Cauchy => corotational_to_Cauchy
    end type Corotational_kinematics

contains

    pure subroutine corotational_update(self, strain, R)
        class(Corotational_kinematics), intent(inout) :: self
        type(ten_3D2Osym), intent(in)                 :: strain
        type(ten_3D2O), intent(in)                    :: R
        real(real64) :: tr_eps

        self%strain = strain
        tr_eps = strain%xx() + strain%yy() + strain%zz()
        self%J = exp(tr_eps)
        self%R_rot = R
        self%has_rotation = .TRUE.
    end subroutine corotational_update

    pure function corotational_jacobian(self) result(J)
        class(Corotational_kinematics), intent(in) :: self
        real(real64)                               :: J
        J = self%J
    end function corotational_jacobian

    pure function corotational_is_finite(self) result(res)
        class(Corotational_kinematics), intent(in) :: self
        logical                                    :: res
        res = .TRUE.
    end function corotational_is_finite

    pure function corotational_get_strain(self) result(strain)
        class(Corotational_kinematics), intent(in) :: self
        type(ten_3D2Osym)                          :: strain
        strain = self%strain
    end function corotational_get_strain

    pure function corotational_push_stress(self, S_material) result(sigma_spatial)
        class(Corotational_kinematics), intent(in) :: self
        type(ten_3D2Osym), intent(in)              :: S_material
        type(ten_3D2Osym)                          :: sigma_spatial

        if (self%has_rotation) then
            sigma_spatial = self%R_rot .transform. S_material
        else
            sigma_spatial = S_material
        end if
    end function corotational_push_stress

    pure function corotational_pull_stress(self, sigma_spatial) result(S_material)
        class(Corotational_kinematics), intent(in) :: self
        type(ten_3D2Osym), intent(in)              :: sigma_spatial
        type(ten_3D2Osym)                          :: S_material

        if (self%has_rotation) then
            S_material = self%R_rot%transpose() .transform. sigma_spatial
        else
            S_material = sigma_spatial
        end if
    end function corotational_pull_stress

    pure function corotational_push_tangent(self, C_material, stress_spatial) result(c_spatial)
        class(Corotational_kinematics), intent(in) :: self
        type(ten_3D4O2sym), intent(in)             :: C_material
        type(ten_3D2Osym), intent(in), optional    :: stress_spatial
        type(ten_3D4O2sym)                         :: c_spatial
        c_spatial = C_material
    end function corotational_push_tangent

    pure function corotational_pull_tangent(self, c_spatial, stress_material) result(C_material)
        class(Corotational_kinematics), intent(in) :: self
        type(ten_3D4O2sym), intent(in)             :: c_spatial
        type(ten_3D2Osym), intent(in), optional    :: stress_material
        type(ten_3D4O2sym)                         :: C_material
        C_material = c_spatial
    end function corotational_pull_tangent

    pure function corotational_to_Cauchy(self, constitutive_stress) result(sigma_cauchy)
        class(Corotational_kinematics), intent(in) :: self
        type(ten_3D2Osym), intent(in)              :: constitutive_stress
        type(ten_3D2Osym)                          :: sigma_cauchy
        if (self%has_rotation) then
            sigma_cauchy = self%R_rot .transform. constitutive_stress
        else
            sigma_cauchy = constitutive_stress
        end if
    end function corotational_to_Cauchy

end module muscle_kin_corotational