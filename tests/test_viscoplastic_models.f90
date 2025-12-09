
module test_viscoplastic_mod
    use mod_basis_viscoplastic_law
    use mod_voce_m_hardening
    use mod_JC_viscoplastic
    use mod_RK_viscoplastic
    use mod_MRK_viscoplastic
    use mod_NNL_viscoplastic
    use mod_VA_viscoplastic
    implicit none
    private

    public :: test_viscoplastic

    double precision, parameter :: tolerance = 0.05D0

    ! ================================
    !   DATOS EXPERIMENTALES
    ! ================================
    double precision, parameter :: ep_exp(15) = (/ &
        1.95456507D-03, 1.01665618D-02, 1.82321656D-02, 2.70017864D-02, &
        3.71373728D-02, 4.60355078D-02, 5.56854902D-02, 6.50158637D-02, &
        7.49742416D-02, 8.57663638D-02, 9.39227091D-02, 1.03962761D-01, &
        1.16123074D-01, 1.30946821D-01, 1.48177371D-01 /)

    double precision, parameter :: stress_exp(15) = (/ &
        2.00064799D+02, 2.56080625D+02, 2.65324635D+02, 2.68667987D+02, &
        2.71711623D+02, 2.74414834D+02, 2.77372186D+02, 2.80012858D+02, &
        2.82808165D+02, 2.85900979D+02, 2.88143097D+02, 2.90925994D+02, &
        2.94237227D+02, 2.98234489D+02, 3.02465123D+02 /)

    double precision, parameter :: epd_const = 6.10D-05

contains

    ! =========================================================
    !   NRMSD
    ! =========================================================
    pure function nrmsd(exp, pred) result(res)
        double precision, intent(in) :: exp(:), pred(:)
        double precision :: res
        double precision :: rmse, denom
        integer :: n

        n = size(exp)
        rmse = sqrt(sum((pred - exp)**2) / dble(n))
        denom = maxval(exp) - minval(exp)

        if (denom .gt. 0.D0) then
            res = rmse / denom
        else
            res = 0.D0
        end if
    end function nrmsd

    subroutine test_viscoplastic(passed)
        logical, intent(out) :: passed

        type(Voce_modified_hardening) :: voce
        type(JC_viscoplastic) :: jc
        type(RK_viscoplastic) :: rk
        type(MRK_viscoplastic) :: mrk
        type(NNL_viscoplastic) :: nnl
        type(VA_viscoplastic) :: va
        class(Base_hardening_law), pointer :: p_hard

        double precision :: pred_JC(15), pred_RK(15), pred_MRK(15)
        double precision :: pred_NNL(15), pred_VA(15)
        double precision :: err_JC, err_RK, err_MRK, err_NNL, err_VA
        logical :: ok_JC, ok_RK, ok_MRK, ok_NNL, ok_VA
        integer :: i

        ! ---- JC ----
        voce%k = 111.497480D0  ! A
        voce%q = 249.653725D0  ! B
        voce%n = 0.08512057D0  ! ns

        p_hard => voce
        jc%hard_law => p_hard
        jc%C      = 0.01006816D0
        jc%epdmax = 0.32D0

        ! ---- RK ----
        rk%B0     = 18.5805384D0
        rk%epdmax = 0.01D0
        rk%nu     = 0.02014578D0
        rk%ep0    = 2.36959316D0
        rk%n0     = 3.05066530D0
        rk%D2     = 1.5302D-08
        rk%D1     = 0.41425354D0
        rk%epdmin = 1.0D-06
        rk%sig0   = 1.1693D-05
        rk%m      = 0.94914871D0

        ! ---- MRK ----
        mrk%B01    = 100.020701D0
        mrk%B02    = 100.0D0
        mrk%epdmax = 0.01D0
        mrk%nu1    = 0.01439459D0
        mrk%nu2    = 0.14313556D0
        mrk%n0     = 3.95295384D0
        mrk%D2     = 5.6375D-12
        mrk%epdmin = 1.0D-05
        mrk%sig_u  = 253.772550D0
        mrk%chi1   = 9.1503D-04
        mrk%chi2   = 0.01942062D0

        ! ---- NNL ----
        nnl%sig_a = 334.776175D0
        nnl%sig_0 = 5.88318003D0
        nnl%KG0   = 0.21906996D0
        nnl%n1    = 0.06301105D0
        nnl%epd0  = 0.01425422D0
        nnl%at    = 22.8752986D0
        nnl%n0    = 0.97427303D0
        nnl%q     = 2.D0
        nnl%p     = 0.6666667D0
        nnl%eps_log  = 1.D-12
        nnl%eps_base = 1.D-12

        ! ---- VA ----
        va%B     = 374.525856D0
        va%B1    = 0.00331794D0
        va%B2    = 0.50280316D0
        va%n     = 0.26412620D0
        va%m     = 0.03766668D0
        va%sig_u = 195.233128D0

        do i = 1, 15
            pred_JC(i)  = jc%flow_stress(ep_exp(i),  epd_const)
            pred_RK(i)  = rk%flow_stress(ep_exp(i),  epd_const)
            pred_MRK(i) = mrk%flow_stress(ep_exp(i), epd_const)
            pred_NNL(i) = nnl%flow_stress(ep_exp(i), epd_const)
            pred_VA(i)  = va%flow_stress(ep_exp(i),  epd_const)
        end do

        ! ============================
        ! 3) NRMSD
        ! ============================
        err_JC  = nrmsd(stress_exp, pred_JC)
        err_RK  = nrmsd(stress_exp, pred_RK)
        err_MRK = nrmsd(stress_exp, pred_MRK)
        err_NNL = nrmsd(stress_exp, pred_NNL)
        err_VA  = nrmsd(stress_exp, pred_VA)

        ok_JC  = (err_JC  < tolerance)
        ok_RK  = (err_RK  < tolerance)
        ok_MRK = (err_MRK < tolerance)
        ok_NNL = (err_NNL < tolerance)
        ok_VA  = (err_VA  < tolerance)

        ! ============================
        ! 4) Results
        ! ============================
        print *, "===================================================="
        print *, "         TEST Viscoplastic Models"
        print *, "NRMSD < ", tolerance
        print *, "===================================================="

        if (ok_JC) then
            print *, "JC   PASSED   NRMSD = ", err_JC
        else
            print *, "JC   FAILED   NRMSD = ", err_JC
        end if

        if (ok_RK) then
            print *, "RK   PASSED   NRMSD = ", err_RK
        else
            print *, "RK   FAILED   NRMSD = ", err_RK
        end if

        if (ok_MRK) then
            print *, "MRK  PASSED   NRMSD = ", err_MRK
        else
            print *, "MRK  FAILED   NRMSD = ", err_MRK
        end if

        if (ok_NNL) then
            print *, "NNL  PASSED   NRMSD = ", err_NNL
        else
            print *, "NNL  FAILED   NRMSD = ", err_NNL
        end if

        if (ok_VA) then
            print *, "VA   PASSED   NRMSD = ", err_VA
        else
            print *, "VA   FAILED   NRMSD = ", err_VA
        end if

        passed = ok_JC .and. ok_RK .and. ok_MRK .and. ok_NNL .and. ok_VA

        print *, "===================================================="
        if (passed) then
            print *, " PASSED"
        else
            print *, " FAILED"
        end if
        print *, "===================================================="

    end subroutine test_viscoplastic

end module test_viscoplastic_mod


! ================================================================
!   Program
! ================================================================
program test_viscoplastic_main
    use test_viscoplastic_mod
    implicit none

    logical :: passed

    call test_viscoplastic(passed)

    if (.not. passed) stop 1
    stop 0
end program test_viscoplastic_main


