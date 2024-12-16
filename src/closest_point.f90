module mod_closest_point
    use tensors_types
    use, intrinsic :: iso_fortran_env, only : real64
    use mod_hardening_law, only : Base_hardening_law
    use mod_yield_criteria, only : Base_yield_critera
    use mod_base_elasticity, only : Base_elasticity
    contains
        subroutine closest_point(strain, elasticity, hardening, yield, strain_pf, strain_p, error, stress)
            implicit none
            class(ten_3D2Osym), intent(in) :: strain
            class(Base_elasticity), intent(in) :: elasticity
            class(Base_hardening_law), intent(in) :: hardening
            class(Base_yield_critera), intent(in) :: yield
            real(real64), intent(inout) :: strain_pf
            type(ten_3D2Osym), intent(inout) :: strain_p
            type(ten_3D2Osym), intent(out) :: stress
            logical, intent(out) :: error
            real(real64), parameter :: TOL=1D-5
            type(ten_3D2Osym) :: df, strain_p_init, residual1
            type(ten_3D4O3sym) :: elas_tan, hess, ddf
            real(real64) :: dgamma, ddgamma, hard, dhard, f, strain_pf_init, norm_res
            real(real64) :: residual(7), residual2
            integer :: i

            error = .False.
            dgamma = 0d0

            strain_p_init = strain_p
            strain_pf_init = strain_pf
            elas_tan = elasticity%dstress_dstrain(strain-strain_p)  ! constant
            do i=1,100
                stress = elasticity%stress(strain-strain_p)
                hard = hardening%stress(strain_pf)
                dhard = hardening%dstress_dep(strain_pf)
                f = yield%stress_eq(stress) - hard
                df = yield%dstressEq_dstress(stress)

                residual1 = strain_p_init - strain_p + (dgamma*df)
                residual2 = strain_pf_init - strain_pf + dgamma

                residual(1:6) = residual1%vals**2
                norm_res = (residual(1)+residual(2)+residual(3)+residual(4)+residual(5)+residual(6))**0.5

                ! write(*,*) "f:", f, norm_res, residual1, dgamma
                write(*,*) "f:", f,  dgamma
                if ((abs(f) .lt. tol) .and. (norm_res .lt. tol)) return  ! Elastic case non varing

                ddf = yield%ddstressEq_ddstress(stress)
                hess = .inv. ((.inv. elas_tan) + dgamma*ddf)
                ! print*, f, residual1 .ddot. hess .ddot. df
                ddgamma = (f-(df .ddot. hess .ddot. residual1) - dhard*residual2)/((df .ddot. hess .ddot. df) + dhard)
                ! ddgamma = (f-(residual1 .ddot. hess .ddot. df))/((df .ddot. hess .ddot. df) + dhard)
                dgamma = dgamma + ddgamma
                strain_pf = strain_pf_init + dgamma
                strain_p = strain_p_init + dgamma*df
                ! write(*,*) "f:", f,  ddgamma, dgamma 

            end do

            error = .True.

        end subroutine closest_point

end module