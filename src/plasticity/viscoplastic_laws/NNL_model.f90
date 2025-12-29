module mod_NNL_viscoplastic
    use, intrinsic :: iso_fortran_env
    use mod_basis_viscoplastic_law
    implicit none
    private
    public :: NNL_viscoplastic

    type, extends(basis_viscoplastic_law) :: NNL_viscoplastic
        real(real64) :: sig_a 
        real(real64) :: n1    
        real(real64) :: sig_0 
        real(real64) :: KG0   
        real(real64) :: epd0  
        real(real64) :: at    
        real(real64) :: n0    
        real(real64) :: q     
        real(real64) :: p     

        real(real64) :: eps_log  = 1.0e-12
        real(real64) :: eps_base = 1.0e-12
    contains
        procedure :: flow_stress => flow_NNL
    end type NNL_viscoplastic

contains

    pure function safe_log_scalar(x, eps) result(res)
        real(real64), intent(in) :: x, eps
        real(real64) :: res
        if (x > eps) then
            res = log(x)
        else
            res = log(eps)
        end if
    end function safe_log_scalar

    pure function sig_u_fn(self, ep) result(res)
        class(NNL_viscoplastic), intent(in) :: self
        real(real64), intent(in) :: ep
        real(real64) :: res
        res = self%sig_a * ep**(self%n1)
    end function sig_u_fn

    pure function sig_as_fn(self, ep, epd) result(res)
        class(NNL_viscoplastic), intent(in) :: self
        real(real64), intent(in) :: ep, epd
        real(real64) :: res
        real(real64) :: term1, term2, S, X, X_pow, base2
        real(real64) :: epsl, epsb

        epsl = self%eps_log
        epsb = self%eps_base

        term1 = safe_log_scalar(epd / self%epd0, epsl)
        term2 = safe_log_scalar(1.0 + self%at * ep**(self%n0), epsl)
        S = term1 + term2

        X = -self%KG0 * S
        if (X < 0.0) X = 0.0

        if (self%q /= 0.0) then
            X_pow = X**(1.0 / self%q)
        else
            X_pow = 0.0
        end if

        base2 = 1.0 - X_pow
        if (base2 < epsb) base2 = epsb

        if (self%p /= 0.0) then
            res = self%sig_0 * base2**(1.0 / self%p) * &
                  (1.0 + self%at * ep**(self%n0))
        else
            res = self%sig_0 * (1.0 + self%at * ep**(self%n0))
        end if
    end function sig_as_fn

    pure function flow_NNL(self, ep, epd) result(res)
        class(NNL_viscoplastic), intent(in) :: self
        real(real64), intent(in) :: ep, epd
        real(real64) :: res
        res = sig_u_fn(self, ep) + sig_as_fn(self, ep, epd)
    end function flow_NNL

end module mod_NNL_viscoplastic
