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
    implicit none
    private

    !------------------------------------------------------------------
    ! Tipo abstracto público
    !------------------------------------------------------------------
    public :: Base_visco_hardening_law

    type, abstract :: Base_visco_hardening_law
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
    contains
        procedure(stress_interface),      deferred :: stress
        !! Devuelve la tensión de fluencia (flow stress) para una ep dada.
        procedure(dstress_dep_interface), deferred :: dstress_dep
        !! Devuelve el módulo de endurecimiento d(stress)/d(ep).
        procedure(ddstress_ddep_interface), deferred :: ddstress_ddep
        !! Devuelve la segunda derivada d²(stress)/d(ep)².
    end type Base_visco_hardening_law

    !==================================================================
    ! Interfaces de los procedimientos diferidos
    !==================================================================
    interface
        pure function stress_interface(self, ep) result(res)
            !! Interface for the `stress` procedure.
            !! Debe ser implementada por subtipos concretos.
            use, intrinsic :: iso_fortran_env, only: real64
            import :: Base_visco_hardening_law
            class(Base_visco_hardening_law), intent(in) :: self
            !! Objeto ley de endurecimiento.
            real(real64), intent(in) :: ep
            !! Deformación plástica equivalente.
            real(real64) :: res
            !! Tensión de fluencia (flow stress) correspondiente.
        end function stress_interface

        pure function dstress_dep_interface(self, ep) result(res)
            !! Interface for the `dstress_dep` procedure.
            !! Debe ser implementada por subtipos concretos.
            use, intrinsic :: iso_fortran_env, only: real64
            import :: Base_visco_hardening_law
            class(Base_visco_hardening_law), intent(in) :: self
            !! Objeto ley de endurecimiento.
            real(real64), intent(in) :: ep
            !! Deformación plástica equivalente en la que se evalúa la derivada.
            real(real64) :: res
            !! Módulo de endurecimiento d(stress)/d(ep).
        end function dstress_dep_interface

        pure function ddstress_ddep_interface(self, ep) result(res)
            !! Interface for the `ddstress_ddep` procedure.
            !! Debe ser implementada por subtipos concretos.
            use, intrinsic :: iso_fortran_env, only: real64
            import :: Base_visco_hardening_law
            class(Base_visco_hardening_law), intent(in) :: self
            !! Objeto ley de endurecimiento.
            real(real64), intent(in) :: ep
            !! Deformación plástica equivalente donde se evalúa la segunda derivada.
            real(real64) :: res
            !! Segunda derivada d²(stress)/d(ep)².
        end function ddstress_ddep_interface

    end interface

end module mod_visco_hardening_law
