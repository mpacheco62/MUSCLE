module mod_JC_viscoplastic
    use, intrinsic :: iso_fortran_env
    use mod_basis_viscoplastic_law
    implicit none
    private
    public :: JC_viscoplastic

    type, extends(basis_viscoplastic_law) :: JC_viscoplastic
        real(real64) :: C     
        real(real64) :: epdmax
    contains
        procedure :: flow_stress => flow_JC
    end type JC_viscoplastic

contains

    pure function flow_JC(self, ep, epd) result(res)
        class(JC_viscoplastic), intent(in) :: self
        real(real64), intent(in) :: ep, epd
        real(real64) :: res
        real(real64) :: sigma0, factor_rate

        if (epd <= 0.0) then
            res = 0.0
            return
        end if

        if (.not. associated(self%hard_law)) then
            res = 0.0
            return
        end if

        sigma0 = self%hard_law%stress(ep)
        factor_rate = 1.0 + self%C * log(epd / self%epdmax)

        res = sigma0 * factor_rate
    end function flow_JC

end module mod_JC_viscoplastic
