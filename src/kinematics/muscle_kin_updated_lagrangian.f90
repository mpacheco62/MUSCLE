module muscle_kin_updated_lagrangian
    use muscle_tensors
    use muscle_kin_finite_base
    implicit none
    private
    public :: Updated_lagrangian_kinematics

    type, extends(Base_F_kinematics) :: Updated_lagrangian_kinematics
        type(ten_3D2Osym) :: strain !! Cached Euler-Almansi strain e
    contains
        procedure :: update_F   => updated_lagrangian_update_F
        procedure :: get_strain => updated_lagrangian_get_strain
        procedure :: to_Cauchy => updated_lagrangian_to_Cauchy
        procedure :: to_spatial_tangent => updated_to_spatial_tangent
    end type Updated_lagrangian_kinematics

contains
    pure subroutine updated_lagrangian_update_F(self, F)
        class(Updated_lagrangian_kinematics), intent(inout) :: self
        type(ten_3D2O), intent(in)                          :: F

        call self%update_F_base(F)
        self%strain = self%get_EulerAlmansi()
    end subroutine updated_lagrangian_update_F

    pure function updated_lagrangian_get_strain(self) result(strain)
        class(Updated_lagrangian_kinematics), intent(in) :: self
        type(ten_3D2Osym)                                :: strain
        strain = self%strain
    end function updated_lagrangian_get_strain

    pure function updated_lagrangian_to_Cauchy(self, constitutive_stress) result(sigma_cauchy)
        class(Updated_lagrangian_kinematics), intent(in) :: self
        type(ten_3D2Osym), intent(in)                    :: constitutive_stress
        type(ten_3D2Osym)                                :: sigma_cauchy
        ! Stress is ALREADY spatial Cauchy stress in Updated Lagrangian!
        sigma_cauchy = constitutive_stress
    end function updated_lagrangian_to_Cauchy

    pure function updated_to_spatial_tangent(self, C_constitutive, stress_spatial) result(c_spatial)
        class(Updated_lagrangian_kinematics), intent(in) :: self
        type(ten_3D4O2sym), intent(in)                   :: C_constitutive
        type(ten_3D2Osym), intent(in), optional          :: stress_spatial
        type(ten_3D4O2sym)                               :: c_spatial

        c_spatial = C_constitutive ! Tangent is ALREADY spatial in Updated Lagrangian
    end function updated_to_spatial_tangent
end module muscle_kin_updated_lagrangian