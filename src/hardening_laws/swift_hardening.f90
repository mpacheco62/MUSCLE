module mod_swift_hardening
    use, intrinsic :: iso_fortran_env
    use mod_hardening_law
    implicit none
    PRIVATE

    PUBLIC :: Swift_hardening
    type, extends(Base_hardening_law) :: Swift_hardening
        real(real64) :: k, e0, n
    contains
        procedure :: stress => stress_swift
        procedure :: dstress_dep => dstress_dep_swift
        procedure :: ddstress_ddep => ddstress_ddep_swift
    end type Swift_hardening

    contains
    pure function stress_swift(self, ep) result(res)
        implicit none
        Class(Swift_hardening), intent(in) :: self
        real(real64), intent(in) :: ep
        real(real64) :: res
        res = self%K*(self%e0 + ep)**self%n
    end function stress_swift

    pure function dstress_dep_swift(self, ep) result(res)
        implicit none
        Class(Swift_hardening), intent(in) :: self  ! parameters
        real(real64), intent(in) :: ep
        real(real64) :: res
        res = self%K*self%n*(self%e0 + ep)**(self%n-1D0)
    end function dstress_dep_swift

    pure function ddstress_ddep_swift(self, ep) result(res)
        implicit none
        Class(Swift_hardening), intent(in) :: self  ! parameters
        real(real64), intent(in) :: ep
        real(real64) :: res
        res = self%K*self%n*(self%n-1D0)*(self%e0 + ep)**(self%n-2D0)
    end function ddstress_ddep_swift
end module