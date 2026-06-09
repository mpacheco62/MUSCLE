module mod_cutting_plane
    use muscle_tensors
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_hard_base, only : Base_hardening_laws
    use muscle_yield_base, only : Base_yield_critera
    use muscle_elasticity_base, only : Base_elasticity
    contains
        subroutine cutting_plane(strain, elasticity, hardening, yield, strain_pf, strain_p, error, stress)
            implicit none
            class(ten_3D2Osym), intent(in) :: strain
            class(Base_elasticity), intent(in) :: elasticity
            class(Base_hardening_laws), intent(in) :: hardening
            class(Base_yield_critera), intent(in) :: yield
            real(real64), intent(inout) :: strain_pf
            type(ten_3D2Osym), intent(inout) :: strain_p
            type(ten_3D2Osym), intent(out) :: stress
            logical, intent(out) :: error
            real(real64), parameter :: TOL=1D-5
            type(ten_3D2Osym) :: df
            type(ten_3D4O3sym) :: elas_tan
            real(real64) :: ddgamma, hard, f, dhard
            integer :: i

            error = .False.

            elas_tan = elasticity%dstress_dstrain(strain-strain_p)  ! constant
            do i=1,100
                stress = elasticity%stress(strain-strain_p)
                hard = hardening%stress(strain_pf)
                f = yield%stress_eq(stress) - hard

                write(*,*) "f:", f
                if (f .lt. tol) return  ! Elastic case non varing

                df = yield%dstressEq_dstress(stress)
                dhard = hardening%dstress_dep(strain_pf)

                ddgamma = f/(((df .ddot. elas_tan) .ddot. df) + dhard)
                strain_pf = strain_pf + ddgamma
                strain_p = strain_p + ddgamma*df
            end do

            error = .True.

        end subroutine cutting_plane

end module