module muscle_vp_nnl
    !! Module muscle_vp_nnl
    !! ===========================
    !! Implementación robusta del modelo Nemat-Nasser-Li (NNL).
    !! Incluye protecciones térmicas y cinemáticas para asegurar convergencia.

    use, intrinsic :: iso_fortran_env
    use muscle_vp_base
    implicit none
    private
    public :: NNL_viscoplastic

    type, extends(Base_viscoplastic_law) :: NNL_viscoplastic
        ! Parámetros del modelo (Tabla 3.4 de la Tesis)
        real(real64) :: sig_a, n1, sig_0, KG0, epd0, at, n0, q, p 
        
        ! Constantes de seguridad numérica
        real(real64) :: eps_log  = 1.0d-12
        real(real64) :: eps_base = 1.0d-12
    contains
        procedure :: flow_stress => flow_NNL
        procedure :: dstress_dep => dstress_dep_NNL
    end type NNL_viscoplastic

contains

    pure function safe_log_scalar(x, eps) result(res)
        real(real64), intent(in) :: x, eps
        real(real64) :: res
        res = log(max(x, eps))
    end function safe_log_scalar

    pure function sig_u_fn(self, ep) result(res)
        class(NNL_viscoplastic), intent(in) :: self
        real(real64), intent(in) :: ep
        real(real64) :: res
        ! Protección para ep=0 con n1 pequeño (evita derivada infinita)
        res = self%sig_a * (max(ep, 1.0d-10)**self%n1)
    end function sig_u_fn

    pure function sig_as_fn(self, ep, epd) result(res)
        class(NNL_viscoplastic), intent(in) :: self
        real(real64), intent(in) :: ep, epd
        real(real64) :: res
        real(real64) :: S, X, X_pow, base2, h_ep, ep_eff, epd_eff

        ! 1. Estabilización de variables de entrada
        ep_eff  = max(ep, 1.0d-10)
        epd_eff = max(epd, self%eps_log)

        ! 2. Término de evolución microestructural
        h_ep = 1.0d0 + self%at * (ep_eff**self%n0)
        
        ! 3. Cálculo de la energía de activación normalizada (S)
        S = safe_log_scalar(epd_eff / self%epd0, self%eps_log) + safe_log_scalar(h_ep, self%eps_log)

        ! 4. Cálculo de X con BLINDAJE: X debe estar estrictamente en [0, 1)
        X = -self%KG0 * S
        X = max(0.0d0, min(0.999d0, X)) 

        ! 5. Perfil de la barrera de energía (p, q)
        X_pow = X**(1.0d0 / self%q)
        base2 = max(self%eps_base, 1.0d0 - X_pow)

        ! 6. Esfuerzo térmico final [cite: 196]
        res = self%sig_0 * (base2**(1.0d0 / self%p)) * h_ep
    end function sig_as_fn

    pure function flow_NNL(self, ep, epd) result(res)
        class(NNL_viscoplastic), intent(in) :: self
        real(real64), intent(in) :: ep, epd
        real(real64) :: res
        ! Suma de componentes atérmica y térmica [cite: 187-188]
        res = sig_u_fn(self, ep) + sig_as_fn(self, ep, epd)
    end function flow_NNL

    pure function dstress_dep_NNL(self, ep, epd, dt) result(res)
        class(NNL_viscoplastic), intent(in) :: self
        real(real64), intent(in) :: ep, epd, dt
        real(real64) :: res
        real(real64) :: ep_eff, epd_eff, S, X, X_pow, base2, factor_h
        real(real64) :: dsig_u_dep, dsig_as_dep, dsig_as_depd
        real(real64) :: dS_dep, dX_dep, dXpow_dep, dbase2_dep
        real(real64) :: dX_depd, dXpow_depd

        ep_eff = max(ep, 1.0d-7) 
        epd_eff = max(epd, self%eps_log)

        dsig_u_dep = self%sig_a * self%n1 * (ep_eff**(self%n1 - 1.0d0))

        factor_h = 1.0d0 + self%at * (ep_eff**self%n0)
        S = log(epd_eff / self%epd0) + log(factor_h)
        X = max(0.0d0, min(0.999d0, -self%KG0 * S))
        
        X_pow = X**(1.0d0 / self%q)
        base2 = max(self%eps_base, 1.0d0 - X_pow)

        dS_dep = (self%at * self%n0 * ep_eff**(self%n0 - 1.0d0)) / factor_h
        dX_dep = -self%KG0 * dS_dep
        
        dXpow_dep = 0.0d0
        if (X > 1.0d-10) then
            dXpow_dep = (1.0d0/self%q) * (X**(1.0d0/self%q - 1.0d0)) * dX_dep
        end if
        dbase2_dep = -dXpow_dep
        
        dsig_as_dep = (self%sig_0 * (1.0d0/self%p) * base2**(1.0d0/self%p - 1.0d0) * dbase2_dep) * factor_h + &
                      (self%sig_0 * base2**(1.0d0/self%p)) * (self%at * self%n0 * ep_eff**(self%n0 - 1.0d0))

        dX_depd = -self%KG0 / epd_eff
        dXpow_depd = 0.0d0
        if (X > 1.0d-10) then
            dXpow_depd = (1.0d0/self%q) * (X**(1.0d0/self%q - 1.0d0)) * dX_depd
        end if
        dsig_as_depd = (self%sig_0 * (1.0d0/self%p) * base2**(1.0d0/self%p - 1.0d0) * (-dXpow_depd)) * factor_h

        res = (dsig_u_dep + dsig_as_dep) + (dsig_as_depd / dt)
    end function dstress_dep_NNL

end module muscle_vp_nnl
