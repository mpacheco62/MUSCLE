! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_kin_small_strain
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_kinematics_base
    implicit none
    private

    public :: Small_strain_kinematics

    type, extends(Base_kinematics) :: Small_strain_kinematics
        type(ten_3D2Osym) :: strain !! Current infinitesimal strain tensor
    contains
    ! Specific update procedures for small strain
        procedure :: update_F      => small_update_F
        procedure :: update_strain => small_update_strain
        generic   :: update        => update_F, update_strain

        procedure :: jacobian           => small_jacobian
        procedure :: is_finite_strain   => small_is_finite_strain
        procedure :: get_strain         => small_get_strain
        procedure :: push_forward_stress => small_push_stress
        procedure :: pull_back_stress    => small_pull_stress
        procedure :: push_forward_tangent=> small_push_tangent
        procedure :: pull_back_tangent  => small_pull_tangent
        procedure :: to_Cauchy => small_to_Cauchy
        procedure :: to_spatial_tangent => small_to_spatial_tangent
    end type Small_strain_kinematics

contains

    pure subroutine small_update_F(self, F)
        class(Small_strain_kinematics), intent(inout) :: self
        type(ten_3D2O), intent(in)                   :: F
        type(ten_3D2O) :: F_sym
        F_sym = 0.5D0 * (F + F%transpose())
        self%strain%vals(1) = F_sym%vals(1) - 1.0D0 ! xx
        self%strain%vals(2) = F_sym%vals(5) - 1.0D0 ! yy
        self%strain%vals(3) = F_sym%vals(9) - 1.0D0 ! zz
        self%strain%vals(4) = F_sym%vals(2)         ! xy
        self%strain%vals(5) = F_sym%vals(6)         ! yz
        self%strain%vals(6) = F_sym%vals(3)         ! xz
    end subroutine small_update_F

    pure subroutine small_update_strain(self, strain)
        class(Small_strain_kinematics), intent(inout) :: self
        type(ten_3D2Osym), intent(in)                 :: strain
        self%strain = strain
    end subroutine small_update_strain

    pure function small_jacobian(self) result(J)
        class(Small_strain_kinematics), intent(in) :: self
        real(real64)                               :: J
        J = 1.0D0
    end function small_jacobian

    pure function small_is_finite_strain(self) result(res)
        class(Small_strain_kinematics), intent(in) :: self
        logical                                    :: res
        res = .FALSE.
    end function small_is_finite_strain

    pure function small_get_strain(self) result(strain)
        class(Small_strain_kinematics), intent(in) :: self
        type(ten_3D2Osym)                          :: strain
        strain = self%strain
    end function small_get_strain

    pure function small_push_stress(self, S_material) result(sigma_spatial)
        class(Small_strain_kinematics), intent(in) :: self
        type(ten_3D2Osym), intent(in)              :: S_material
        type(ten_3D2Osym)                          :: sigma_spatial
        sigma_spatial = S_material
    end function small_push_stress

    pure function small_pull_stress(self, sigma_spatial) result(S_material)
        class(Small_strain_kinematics), intent(in) :: self
        type(ten_3D2Osym), intent(in)              :: sigma_spatial
        type(ten_3D2Osym)                          :: S_material
        S_material = sigma_spatial
    end function small_pull_stress

    pure function small_push_tangent(self, C_material, stress_spatial) result(c_spatial)
        class(Small_strain_kinematics), intent(in) :: self
        type(ten_3D4O2sym), intent(in)             :: C_material
        type(ten_3D2Osym), intent(in), optional    :: stress_spatial
        type(ten_3D4O2sym)                         :: c_spatial
        c_spatial = C_material
    end function small_push_tangent

    pure function small_pull_tangent(self, c_spatial, stress_material) result(C_material)
        class(Small_strain_kinematics), intent(in) :: self
        type(ten_3D4O2sym), intent(in)             :: c_spatial
        type(ten_3D2Osym), intent(in), optional    :: stress_material
        type(ten_3D4O2sym)                         :: C_material
        C_material = c_spatial
    end function small_pull_tangent

    pure function small_to_Cauchy(self, constitutive_stress) result(sigma_cauchy)
        class(Small_strain_kinematics), intent(in) :: self
        type(ten_3D2Osym), intent(in)              :: constitutive_stress
        type(ten_3D2Osym)                          :: sigma_cauchy
        sigma_cauchy = constitutive_stress
    end function small_to_Cauchy

    pure function small_to_spatial_tangent(self, C_constitutive, stress_spatial) result(c_spatial)
        class(Small_strain_kinematics), intent(in) :: self
        type(ten_3D4O2sym), intent(in)             :: C_constitutive
        type(ten_3D2Osym), intent(in), optional    :: stress_spatial
        type(ten_3D4O2sym)                         :: c_spatial

        c_spatial = C_constitutive
    end function small_to_spatial_tangent

end module muscle_kin_small_strain