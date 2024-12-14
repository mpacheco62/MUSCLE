module mod_base_elasticity
    use, intrinsic :: iso_fortran_env
    use tensors_types
    implicit None

    type, abstract :: Base_elasticity
        contains
            procedure(stress_interface), deferred :: stress
            procedure(dstress_dstrain_interface), deferred :: dstress_dstrain
    end type Base_elasticity

    interface
        pure function stress_interface(self, strain) result(res)
            use, intrinsic :: iso_fortran_env
            use tensors_types
            import Base_elasticity
            class(Base_elasticity), intent(in) :: self
            class(ten_3D2Osym), intent(in) :: strain
            type(ten_3D2Osym) :: res
        end function stress_interface

        function dstress_dstrain_interface(self, strain) result(res)
            use, intrinsic :: iso_fortran_env
            use tensors_types
            import Base_elasticity
            class(Base_elasticity), intent(inout) :: self
            class(ten_3D2Osym), intent(in) :: strain
            type(ten_3D4O3sym) :: res
        end function dstress_dstrain_interface

    end interface

end module