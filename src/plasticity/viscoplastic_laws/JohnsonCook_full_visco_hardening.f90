!======================================================================
! File: JohnsonCook_full_visco_hardening.f90
!======================================================================
! Module mod_JohnsonCook_full_visco_hardening
! ===========================================
!
! Implementa el modelo Johnson–Cook completo usando la base viscoplástica:
!
!   sigma_y(ep, epdot, T)
!     = sigma_hard(ep)
!       * [ 1 + C * log( epdot / epdot0 ) ]
!       * [ 1 - (Tstar)^m ]
!
! con:
!
!   Tstar = (T - Troom) / (Tmelt - Troom)   limitado a [0,1].
!
! Nota:
! -----
!  - Este modelo extiende Base_viscoplastic_law (como JC_model).
!  - La parte de endurecimiento sigma_hard(ep) se obtiene desde hard_law,
!    o, si no está asociada, desde los parámetros A, B, n.
!  - La temperatura actual se almacena en T_current y debe actualizarse
!    desde el algoritmo local antes de llamar a flow_stress().
!
!======================================================================

module mod_JohnsonCook_full_visco_hardening
    use, intrinsic :: iso_fortran_env, only : real64
    use mod_viscoplastic_law, only : Base_viscoplastic_law
    implicit none
    private

    public :: JohnsonCook_full_visco_hardening

    type, extends(Base_viscoplastic_law) :: JohnsonCook_full_visco_hardening
        !! Johnson–Cook full viscoplastic law with thermal softening:
        !!
        !!   sigma_y(ep, epdot, T)
        !!     = sigma_hard(ep)
        !!       * [ 1 + C * log( epdot / epdot0 ) ]
        !!       * [ 1 - (Tstar)^m ]
        !!
        !! Parámetros de endurecimiento (fallback si no hay hard_law):
        real(real64) :: A      = 0.0d0    !! Límite elástico cuasiestático
        real(real64) :: B      = 0.0d0    !! Coeficiente de endurecimiento
        real(real64) :: n      = 0.0d0    !! Exponente de endurecimiento
        !! Parámetros viscoplásticos y térmicos:
        real(real64) :: C      = 0.0d0    !! Sensibilidad a la tasa
        real(real64) :: epdmax = 1.0d0    !! Tasa de referencia
        real(real64) :: m      = 0.0d0    !! Exponente ablandamiento térmico
        real(real64) :: Troom  = 293.15d0 !! Temp. de referencia (K)
        real(real64) :: Tmelt  = 1.0d3    !! Temp. de fusión aproximada (K)
        real(real64) :: T_current = 0.0d0 !! Temp. actual (K)
    contains
        procedure :: flow_stress => flow_JC_full
    end type JohnsonCook_full_visco_hardening

contains

    pure function flow_JC_full(self, ep, epd) result(res)
        class(JohnsonCook_full_visco_hardening), intent(in) :: self
        real(real64), intent(in) :: ep
        real(real64), intent(in) :: epd
        real(real64) :: res
        real(real64) :: sigma0, factor_rate, factor_temp
        real(real64) :: epd_eff, Tstar, denomT

        if (epd <= 0.0d0) then
            res = 0.0d0
            return
        end if

        if (associated(self%hard_law)) then
            sigma0 = self%hard_law%stress(ep)
        else
            sigma0 = self%A + self%B * ep**self%n
        end if

        epd_eff = max(epd, tiny(1.0d0))
        factor_rate = 1.0d0 + self%C * log(epd_eff / self%epdmax)

        denomT = self%Tmelt - self%Troom
        if (denomT > 0.0d0) then
            Tstar = (self%T_current - self%Troom) / denomT
        else
            Tstar = 0.0d0
        end if
        Tstar = max(0.0d0, min(1.0d0, Tstar))

        if (self%m /= 0.0d0) then
            factor_temp = 1.0d0 - Tstar**self%m
        else
            factor_temp = 1.0d0
        end if

        res = sigma0 * factor_rate * factor_temp
    end function flow_JC_full

end module mod_JohnsonCook_full_visco_hardening
