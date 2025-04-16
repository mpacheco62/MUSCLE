module mod_hardening_law
    use, intrinsic :: iso_fortran_env
    implicit None
    PRIVATE
    
    PUBLIC :: Base_hardening_law
    type, abstract :: Base_hardening_law
        contains
            procedure(stress_interface), deferred :: stress
            procedure(dstress_dep_interface), deferred :: dstress_dep
            procedure(ddstress_ddep_interface), deferred :: ddstress_ddep
    end type Base_hardening_law

    interface
        pure function stress_interface(self, ep) result(res)
            use, intrinsic :: iso_fortran_env
            import Base_hardening_law
            class(Base_hardening_law), intent(in) :: self
            real(real64), intent(in) :: ep
            real(real64) :: res
        end function stress_interface

        pure function dstress_dep_interface(self, ep) result(res)
            use, intrinsic :: iso_fortran_env
            import Base_hardening_law
            class(Base_hardening_law), intent(in) :: self
            real(real64), intent(in) :: ep
            real(real64) :: res
        end function dstress_dep_interface

        pure function ddstress_ddep_interface(self, ep) result(res)
            use, intrinsic :: iso_fortran_env
            import Base_hardening_law
            class(Base_hardening_law), intent(in) :: self
            real(real64), intent(in) :: ep
            real(real64) :: res
        end function ddstress_ddep_interface

    end interface
end module