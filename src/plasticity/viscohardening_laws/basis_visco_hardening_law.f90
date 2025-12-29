!======================================================================
! File: basis_visco_hardening_law.f90
!======================================================================
! Module mod_visco_hardening_law
! ==============================
!
! Define el tipo abstracto base para leyes de endurecimiento que se usan
! en modelos viscoplásticos.
!
! La idea es tener una interfaz común para distintas leyes de endurecimiento
! isotrópico (lineal, ley de potencia, saturación, Johnson–Cook, etc.),
! de forma que el código del modelo constitutivo pueda trabajar de manera
! polimórfica con cualquier ley concreta.
!
! Cada ley concreta debe extender el tipo abstracto
!     type, extends(Base_visco_hardening_law) :: MiLey
! y proveer implementaciones para:
!   - stress(ep)         : tensión de fluencia (flow stress)
!   - dstress_dep(ep)    : d(stress)/d(ep)          (módulo de endurecimiento)
!   - ddstress_ddep(ep)  : d²(stress)/d(ep)²        (segunda derivada)
!
! Nota:
! -----
!  - El argumento "ep" es la deformación plástica equivalente acumulada.
!  - La dependencia con tasa de deformación y temperatura (viscoplasticidad)
!    puede entrar a través de parámetros almacenados en el objeto (self),
!    que se actualizan desde el algoritmo de integración local antes de
!    llamar a estas funciones.
!
!======================================================================

module mod_visco_hardening_law
    use, intrinsic :: iso_fortran_env, only: real64
    use mod_hardening_law, only: Base_hardening_law
    implicit none
    private

    !------------------------------------------------------------------
    ! Tipo abstracto público
    !------------------------------------------------------------------
    public :: Base_visco_hardening_law

    type, abstract, extends(Base_hardening_law) :: Base_visco_hardening_law
        !! Abstract Base Type for Isotropic (Visco-)Hardening Laws
        !!
        !! Sirve como interfaz fundamental para todas las leyes de
        !! endurecimiento isotrópico dentro de la librería.
        !!
        !! Cualquier modelo concreto (lineal, ley de potencia, saturación,
        !! Johnson–Cook, etc.) debe extender este tipo y proporcionar
        !! implementaciones para:
        !!   - stress(ep)
        !!   - dstress_dep(ep)
        !!   - ddstress_ddep(ep)
        !!
        !! La variable independiente principal es la deformación plástica
        !! equivalente acumulada "ep". Dependencias adicionales en tasa
        !! de deformación y temperatura se manejan como parámetros internos
        !! del propio objeto (self).

         !! Aquí puedes agregar cosas propias de "visco", ej:
        real(real64) :: epdot_current = 0.0d0
        real(real64) :: T_current     = 0.0d0
        !! Los procedimientos diferidos (stress, dstress_dep, ddstress_ddep)
        !! se heredan de Base_hardening_law; no hace falta redeclararlos aquí.
    end type Base_visco_hardening_law

end module mod_visco_hardening_law
