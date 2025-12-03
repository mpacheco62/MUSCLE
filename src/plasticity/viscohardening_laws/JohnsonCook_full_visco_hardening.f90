!======================================================================
! File: JohnsonCook_full_visco_hardening.f90
!======================================================================
! Module mod_JohnsonCook_full_visco_hardening
! ===========================================
!
! Implementa la ley de endurecimiento tipo Johnson–Cook completa:
!
!   sigma_y(ep, epdot, T)
!     = (A + B * ep^n)
!       * [ 1 + C * log( epdot / epdot0 ) ]
!       * [ 1 - (Tstar)^m ]
!
! con:
!
!   Tstar = (T - Troom) / (Tmelt - Troom)   limitado a [0,1].
!
! La interfaz de la clase base sólo pasa ep, por lo que la tasa de
! deformación plástica y la temperatura se almacenan como campos
! internos del tipo (epdot_current, T_current), que deben ser
! actualizados por el algoritmo de integración local antes de llamar
! a stress(), dstress_dep(), ddstress_ddep().
!
!======================================================================

module mod_JohnsonCook_full_visco_hardening
    use, intrinsic :: iso_fortran_env, only: real64
    use mod_visco_hardening_law
    implicit none
    private

    !------------------------------------------------------------------
    ! Tipo público
    !------------------------------------------------------------------
    public :: JohnsonCook_full_visco_hardening

    type, extends(Base_visco_hardening_law) :: JohnsonCook_full_visco_hardening
        !! Johnson–Cook full isotropic visco-hardening law:
        !!
        !!   sigma_y(ep, epdot, T)
        !!     = (A + B * ep^n)
        !!       * [ 1 + C * log( epdot / epdot0 ) ]
        !!       * [ 1 - (Tstar)^m ]
        !!
        !! donde:
        !!   - ep           : deformación plástica equivalente
        !!   - epdot        : tasa de deformación plástica equivalente
        !!   - T            : temperatura
        !!   - Tstar        : temperatura adimensional en [0,1]
        !!
        !! Parámetros del modelo:
        real(real64) :: A      = 0.0d0  !! Límite elástico cuasiestático (ref. Troom, epdot0)
        real(real64) :: B      = 0.0d0  !! Coeficiente de endurecimiento por deformación
        real(real64) :: n      = 0.0d0  !! Exponente de endurecimiento
        real(real64) :: C      = 0.0d0  !! Parámetro de sensibilidad a la tasa
        real(real64) :: epdot0 = 1.0d0  !! Tasa de deformación de referencia
        real(real64) :: m      = 0.0d0  !! Exponente de ablandamiento térmico
        real(real64) :: Troom  = 293.15d0 !! Temperatura de referencia (K)
        real(real64) :: Tmelt  = 1.0d3    !! Temperatura de fusión aproximada (K)
        !!
        !! Estado viscoplástico "externo" (se actualiza desde el modelo local):
        real(real64) :: epdot_current = 1.0d0   !! Tasa de deformación plástica equivalente actual
        real(real64) :: T_current     = 293.15d0 !! Temperatura actual
    contains
        procedure :: stress        => stress_JC_full
        procedure :: dstress_dep   => dstress_dep_JC_full
        procedure :: ddstress_ddep => ddstress_ddep_JC_full
    end type JohnsonCook_full_visco_hardening

contains

    !------------------------------------------------------------------
    ! sigma_y(ep) = F(ep) * G(epdot_current) * H(T_current)
    !  F(ep) = A + B ep^n
    !  G(epdot) = 1 + C log(epdot / epdot0)
    !  H(T) = 1 - (Tstar)^m
    !------------------------------------------------------------------
    pure function stress_JC_full(self, ep) result(res)
        class(JohnsonCook_full_visco_hardening), intent(in) :: self
        real(real64), intent(in) :: ep
        real(real64) :: res
        real(real64) :: F, G, H
        real(real64) :: epdot_eff, Tstar, denomT

        ! Parte de endurecimiento por deformación
        F = self%A + self%B * ep**self%n

        ! Parte de tasa de deformación: proteger contra epdot <= 0
        epdot_eff = max(self%epdot_current, tiny(1.0d0))
        G = 1.0d0 + self%C * log(epdot_eff / self%epdot0)

        ! Parte térmica: Tstar en [0,1]
        denomT = self%Tmelt - self%Troom
        if (denomT > 0.0d0) then
            Tstar = (self%T_current - self%Troom) / denomT
        else
            Tstar = 0.0d0
        end if
        Tstar = max(0.0d0, min(1.0d0, Tstar))

        if (self%m /= 0.0d0) then
            H = 1.0d0 - Tstar**self%m
        else
            H = 1.0d0
        end if

        res = F * G * H
    end function stress_JC_full

    !------------------------------------------------------------------
    ! d sigma_y / d ep = F'(ep) * G * H
    !  F'(ep) = B n ep^{n-1}
    !------------------------------------------------------------------
    pure function dstress_dep_JC_full(self, ep) result(res)
        class(JohnsonCook_full_visco_hardening), intent(in) :: self
        real(real64), intent(in) :: ep
        real(real64) :: res
        real(real64) :: Fprime, G, H
        real(real64) :: ep_eff, epdot_eff, Tstar, denomT

        ! Para evitar problemas numéricos con ep muy pequeño
        ep_eff = max(ep, 1.0d-16)

        ! Derivada de la parte de endurecimiento
        Fprime = self%B * self%n * ep_eff**(self%n - 1.0d0)

        ! Mismo G y H que en stress()
        epdot_eff = max(self%epdot_current, tiny(1.0d0))
        G = 1.0d0 + self%C * log(epdot_eff / self%epdot0)

        denomT = self%Tmelt - self%Troom
        if (denomT > 0.0d0) then
            Tstar = (self%T_current - self%Troom) / denomT
        else
            Tstar = 0.0d0
        end if
        Tstar = max(0.0d0, min(1.0d0, Tstar))

        if (self%m /= 0.0d0) then
            H = 1.0d0 - Tstar**self%m
        else
            H = 1.0d0
        end if

        res = Fprime * G * H
    end function dstress_dep_JC_full

    !------------------------------------------------------------------
    ! d^2 sigma_y / d ep^2 = F''(ep) * G * H
    !  F''(ep) = B n (n-1) ep^{n-2}
    !------------------------------------------------------------------
    pure function ddstress_ddep_JC_full(self, ep) result(res)
        class(JohnsonCook_full_visco_hardening), intent(in) :: self
        real(real64), intent(in) :: ep
        real(real64) :: res
        real(real64) :: Fsecond, G, H
        real(real64) :: ep_eff, epdot_eff, Tstar, denomT

        ! Proteger contra ep muy pequeño
        ep_eff = max(ep, 1.0d-16)

        ! Segunda derivada de la parte de endurecimiento
        Fsecond = self%B * self%n * (self%n - 1.0d0) * ep_eff**(self%n - 2.0d0)

        ! Mismo G y H que en stress()
        epdot_eff = max(self%epdot_current, tiny(1.0d0))
        G = 1.0d0 + self%C * log(epdot_eff / self%epdot0)

        denomT = self%Tmelt - self%Troom
        if (denomT > 0.0d0) then
            Tstar = (self%T_current - self%Troom) / denomT
        else
            Tstar = 0.0d0
        end if
        Tstar = max(0.0d0, min(1.0d0, Tstar))

        if (self%m /= 0.0d0) then
            H = 1.0d0 - Tstar**self%m
        else
            H = 1.0d0
        end if

        res = Fsecond * G * H
    end function ddstress_ddep_JC_full

end module mod_JohnsonCook_full_visco_hardening
