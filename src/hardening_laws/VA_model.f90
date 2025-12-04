module mod_VA_viscoplastic
    use, intrinsic :: iso_fortran_env
    use mod_basis_viscoplastic_law
    implicit none
    private
    public :: VA_viscoplastic

    type, extends(basis_viscoplastic_law) :: VA_viscoplastic
        real(real64) :: B    
        real(real64) :: n    
        real(real64) :: B1   
        real(real64) :: B2   
        real(real64) :: m    
        real(real64) :: sig_u
    contains
        procedure :: flow_stress => flow_VA
    end type VA_viscoplastic

contains

    pure function flow_VA(self, ep, epd) result(res)
        class(VA_viscoplastic), intent(in) :: self
        real(real64), intent(in) :: ep, epd
        real(real64) :: res
        real(real64) :: factor_rate

        if (self%m /= 0.0) then
            factor_rate = 1.0 - self%B1 * epd**(1.0/self%m) - self%B2
        else
            factor_rate = 1.0 - self%B2
        end if

        res = self%B * ep**(self%n) * factor_rate + self%sig_u
    end function flow_VA

end module mod_VA_viscoplastic
