module mod_MRK_viscoplastic
    use, intrinsic :: iso_fortran_env
    use mod_basis_viscoplastic_law
    implicit none
    private
    public :: MRK_viscoplastic

    type, extends(basis_viscoplastic_law) :: MRK_viscoplastic
        real(real64) :: B01    
        real(real64) :: B02    
        real(real64) :: epdmax 
        real(real64) :: nu1    
        real(real64) :: nu2    
        real(real64) :: n0     
        real(real64) :: D2     
        real(real64) :: epdmin 
        real(real64) :: chi1  
        real(real64) :: chi2  
        real(real64) :: sig_u 
    contains
        procedure :: flow_stress => flow_MRK
    end type MRK_viscoplastic

contains

    pure function flow_MRK(self, ep, epd) result(res)
        class(MRK_viscoplastic), intent(in) :: self
        real(real64), intent(in) :: ep, epd
        real(real64) :: res
        real(real64) :: n_val, B1_val, B2_val
        real(real64) :: sig0_a, sig_a
        real(real64) :: factor_n, factor_rate, log_arg

        if (epd <= 0.0) then
            res = self%sig_u
            return
        end if

        factor_n = 1.0 - self%D2 * log10(epd / self%epdmin)
        if (factor_n > 0.0) then
            n_val = self%n0 * factor_n
        else
            n_val = 0.0
        end if

        log_arg = log10(self%epdmax / epd)
        B1_val  = self%B01 * log_arg**(-self%nu1)
        B2_val  = self%B02 * log_arg**(-self%nu2)

        sig0_a = B1_val*ep + B2_val*(1.0 - exp(-ep*n_val))

        factor_rate = 1.0 - self%chi1 * log10(self%epdmax / epd)
        if (factor_rate > 0.0) then
            sig_a = factor_rate**(1.0/self%chi2) * sig0_a
        else
            sig_a = 0.0
        end if

        res = sig_a + self%sig_u
    end function flow_MRK

end module mod_MRK_viscoplastic
