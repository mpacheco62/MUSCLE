!======================================================================
! File: JohnsonCook_visco_hardening.f90
!======================================================================
! Module mod_JohnsonCook_visco_hardening
! ======================================
!
! Define una implementación concreta de la ley de endurecimiento
! tipo potencia:
!
!   sigma_flow(ep) = K * (e0 + ep)^n
!
! que corresponde a la parte de endurecimiento por deformación de
! un modelo tipo Johnson–Cook / Hollomon.
!
! Esta ley depende sólo de la deformación plástica equivalente
! acumulada ep. La dependencia viscoplástica (tasa y temperatura)
! se puede introducir asignando parámetros K, e0, n adecuados
! en cada estado de carga desde el algoritmo de integración local.
!
!======================================================================

module mod_JohnsonCook_visco_hardening
    use, intrinsic :: iso_fortran_env, only: real64
    use mod_visco_hardening_law
    implicit none
    private

    !------------------------------------------------------------------
    ! Tipo público
    !------------------------------------------------------------------
    public :: JohnsonCook_visco_hardening

    type, extends(Base_visco_hardening_law) :: JohnsonCook_visco_hardening
        !! Concrete type for a Johnson–Cook/Hollomon-like isotropic
        !! hardening law:
        !!
        !!   sigma_flow = K * (e0 + ep)^n
        !!
        !! donde:
        !!   - K  : coeficiente de resistencia
        !!   - e0 : offset de deformación (para evitar ep = 0)
        !!   - n  : exponente de endurecimiento
        real(real64) :: k  = 0.0d0  !! Strength coefficient (K)
        real(real64) :: e0 = 0.0d0  !! Initial strain offset (epsilon_0)
        real(real64) :: n  = 0.0d0  !! Hardening exponent (n)
    contains
        procedure :: stress        => stress_JohnsonCook
        procedure :: dstress_dep   => dstress_dep_JohnsonCook
        procedure :: ddstress_ddep => ddstress_ddep_JohnsonCook
    end type JohnsonCook_visco_hardening

contains

    !------------------------------------------------------------------
    ! sigma_flow(ep) = K * (e0 + ep)^n
    !------------------------------------------------------------------
    pure function stress_JohnsonCook(self, ep) result(res)
        !! Calcula la tensión de fluencia según:
        !!   sigma_flow = K * (e0 + ep)^n.
        class(JohnsonCook_visco_hardening), intent(in) :: self
        real(real64), intent(in) :: ep
        real(real64) :: res

        res = self%k * (self%e0 + ep)**self%n
    end function stress_JohnsonCook

    !------------------------------------------------------------------
    ! d sigma_flow / d ep = K * n * (e0 + ep)^{n - 1}
    !------------------------------------------------------------------
    pure function dstress_dep_JohnsonCook(self, ep) result(res)
        !! Calcula el módulo de endurecimiento:
        !!   H = d(sigma_flow)/d(ep) = K * n * (e0 + ep)^(n - 1).
        class(JohnsonCook_visco_hardening), intent(in) :: self
        real(real64), intent(in) :: ep
        real(real64) :: res

        res = self%k * self%n * (self%e0 + ep)**(self%n - 1.0d0)
    end function dstress_dep_JohnsonCook

    !------------------------------------------------------------------
    ! d^2 sigma_flow / d ep^2 = K * n * (n - 1) * (e0 + ep)^{n - 2}
    !------------------------------------------------------------------
    pure function ddstress_ddep_JohnsonCook(self, ep) result(res)
        !! Calcula la segunda derivada:
        !!   d^2(sigma_flow)/d(ep)^2
        !!   = K * n * (n - 1) * (e0 + ep)^(n - 2).
        class(JohnsonCook_visco_hardening), intent(in) :: self
        real(real64), intent(in) :: ep
        real(real64) :: res

        res = self%k * self%n * (self%n - 1.0d0) * (self%e0 + ep)**(self%n - 2.0d0)
    end function ddstress_ddep_JohnsonCook

end module mod_JohnsonCook_visco_hardening
