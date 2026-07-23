module muscle_kin_total_lagrangian
    use muscle_tensors
    use muscle_kin_finite_base
    implicit none
    private
    public :: Total_lagrangian_kinematics

    type, extends(Base_F_kinematics) :: Total_lagrangian_kinematics
        type(ten_3D2Osym) :: strain !! Cached Green-Lagrange strain E
    contains
        procedure :: update_F   => total_lagrangian_update_F
        procedure :: get_strain => total_lagrangian_get_strain
        procedure :: to_Cauchy => total_lagrangian_to_Cauchy
    end type Total_lagrangian_kinematics

contains
    pure subroutine total_lagrangian_update_F(self, F)
        class(Total_lagrangian_kinematics), intent(inout) :: self
        type(ten_3D2O), intent(in)                        :: F

        call self%update_F_base(F)
        self%strain = self%get_GreenLagrange()
    end subroutine total_lagrangian_update_F

    pure function total_lagrangian_get_strain(self) result(strain)
        class(Total_lagrangian_kinematics), intent(in) :: self
        type(ten_3D2Osym)                             :: strain
        strain = self%strain
    end function total_lagrangian_get_strain

    pure function total_lagrangian_to_Cauchy(self, constitutive_stress) result(sigma_cauchy)
        class(Total_lagrangian_kinematics), intent(in) :: self
        type(ten_3D2Osym), intent(in)                  :: constitutive_stress
        type(ten_3D2Osym)                              :: sigma_cauchy
        ! Push-forward PK2 -> Cauchy
        sigma_cauchy = (self%F .transform. constitutive_stress) / self%J
    end function total_lagrangian_to_Cauchy
end module muscle_kin_total_lagrangian