module test_viscoplastic_mod
    use mod_voce_m_hardening
    use mod_JC_viscoplastic
    use mod_RK_viscoplastic
    use mod_MRK_viscoplastic
    use mod_NNL_viscoplastic
    use mod_VA_viscoplastic
    implicit none
    private
    public :: test_viscoplastic

    integer, parameter :: nRates = 5
    integer, parameter :: nPts   = 3
    double precision, parameter :: tolerance = 1e-5

    ! ----------------------------
    ! Tasa de deformación 
    ! ----------------------------
    double precision, parameter :: epd_rates(nRates) = (/ &
        6.10D-05, 1.90D-04, 7.20D-04, 1.80D-03, 7.00D-03 /)

    character(len=8), parameter :: rate_name(nRates) = (/ &
        "0.125   ", "0.5     ", "2.0     ", "8.0     ", "32.0    " /)

    ! ==========================================================
    ! Datos experimentales

    double precision, parameter :: ep_test(nPts, nRates) = reshape((/ &
        ! 0.125
        1.0167D-02, 6.5016D-02, 1.48177D-01, &
        ! 0.5
        8.3060D-03, 7.0821D-02, 1.41501D-01, &
        ! 2.0
        1.0442D-02, 8.5531D-02, 1.86225D-01, &
        ! 8.0
        1.6275D-02, 9.3995D-02, 1.68051D-01, &
        ! 32.0
        2.1193D-02, 7.9007D-02, 1.76218D-01  &
    /), (/nPts, nRates/))

    double precision, parameter :: ref_JC(nPts, nRates) = reshape((/ &
        2.54122901D+02, 2.75180654D+02, 3.04011499D+02, &
        2.56341724D+02, 2.80554401D+02, 3.05356309D+02, &
        2.60693918D+02, 2.89809269D+02, 3.24093795D+02, &
        2.65493492D+02, 2.95586794D+02, 3.21206442D+02, &
        2.71116773D+02, 2.93969226D+02, 3.28239369D+02  &
    /), (/nPts, nRates/))

    double precision, parameter :: ref_RK(nPts, nRates) = reshape((/ &
        2.57500594D+02, 2.76037431D+02, 3.05821337D+02, &
        2.58194546D+02, 2.79465716D+02, 3.04898323D+02, &
        2.61048679D+02, 2.86995567D+02, 3.24435956D+02, &
        2.65283137D+02, 2.92536674D+02, 3.20198734D+02, &
        2.75534653D+02, 2.96369472D+02, 3.33745319D+02  &
    /), (/nPts, nRates/))

    double precision, parameter :: ref_MRK(nPts, nRates) = reshape((/ &
        2.57844849D+02, 2.77783408D+02, 3.02609442D+02, &
        2.57287918D+02, 2.81083725D+02, 3.03268359D+02, &
        2.58519646D+02, 2.88532576D+02, 3.19826380D+02, &
        2.61625905D+02, 2.94143186D+02, 3.18951988D+02, &
        2.66352743D+02, 2.96608134D+02, 3.36791644D+02  &
    /), (/nPts, nRates/))

    double precision, parameter :: ref_NNL(nPts, nRates) = reshape((/ &
        2.50717404D+02, 2.81822545D+02, 2.97339220D+02, &
        2.47622744D+02, 2.84254012D+02, 2.98493682D+02, &
        2.51932894D+02, 2.90674594D+02, 3.11406703D+02, &
        2.60253715D+02, 2.96555422D+02, 3.15925555D+02, &
        2.68476571D+02, 3.02529156D+02, 3.30768917D+02  &
    /), (/nPts, nRates/))

    double precision, parameter :: ref_VA(nPts, nRates) = reshape((/ &
        2.50651355D+02, 2.85701644D+02, 3.07691536D+02, &
        2.47769724D+02, 2.87768484D+02, 3.06330367D+02, &
        2.51044682D+02, 2.92498114D+02, 3.14689185D+02, &
        2.57984507D+02, 2.94952719D+02, 3.11492696D+02, &
        2.62517373D+02, 2.90480900D+02, 3.12959149D+02  &
    /), (/nPts, nRates/))

contains

    pure function nrmsd(ref, pred) result(res)
        double precision, intent(in) :: ref(:), pred(:)
        double precision :: res, rmse, denom
        integer :: n

        n = size(ref)
        rmse  = sqrt( sum( (pred - ref)**2 ) / dble(n) )
        denom = maxval(ref) - minval(ref)

        if (denom .gt. 0.D0) then
            res = rmse / denom
        else
            res = 0.D0
        end if
    end function nrmsd


    subroutine test_viscoplastic(passed)
        logical, intent(out) :: passed

        type(Voce_modified_hardening) :: voce
        type(JC_viscoplastic)  :: jc
        type(RK_viscoplastic)  :: rk
        type(MRK_viscoplastic) :: mrk
        type(NNL_viscoplastic) :: nnl
        type(VA_viscoplastic)  :: va

        double precision :: pred_JC(nPts), pred_RK(nPts), pred_MRK(nPts)
        double precision :: pred_NNL(nPts), pred_VA(nPts)
        double precision :: errJC, errRK, errMRK, errNNL, errVA
        logical :: okJC, okRK, okMRK, okNNL, okVA
        logical :: all_ok
        integer :: r, i

        character(len=7)  :: stJC, stRK, stMRK, stNNL, stVA
        character(len=30) :: stFINAL

        ! ---- JC
        voce = Voce_modified_hardening(sy=272.075784D0, k=0.04230655D0, q=297.603923D0, n=1.48375781D0)
        jc = JC_viscoplastic(hard_law=voce, C=0.00946101D0, epdmax=0.32D0)

        ! ---- RK ----
        rk = RK_viscoplastic(B0 = 18.5805384D0, epdmax = 0.01D0, nu = 0.02014578D0, &
                             ep0 = 2.36959316D0, n0 = 3.05066530D0, D2 = 1.5302D-08, D1 = 0.41425354D0, &
                             epdmin = 1.0D-06, sig0 = 1.1693D-05, m = 0.94914871D0)

        ! ---- MRK ----
        mrk = MRK_viscoplastic(B01=100.020701D0, B02=100.0D0, epdmax=0.01D0, &
                               nu1=0.01439459D0, nu2=0.14313556D0, n0=3.95295384D0, D2=5.6375D-12, &
                               epdmin=1.0D-05, sig_u=253.772550D0, chi1=9.1503D-04, chi2=0.01942062D0)

        ! ---- NNL ----
        nnl = NNL_viscoplastic(sig_a=334.776175D0, sig_0=5.88318003D0, KG0=0.21906996D0, &
                               n1=0.06301105D0, epd0=0.01425422D0, at=22.8752986D0, &
                               n0=0.97427303D0, q=2.D0, p=0.6666667D0)

        ! ---- VA ----
        va = VA_viscoplastic(B=374.525856D0, B1=0.00331794D0, B2=0.50280316D0, &
                             n=0.26412620D0, m=0.03766668D0, sig_u=195.233128D0)
        va%B     = 374.525856D0
        va%B1    = 0.00331794D0
        va%B2    = 0.50280316D0
        va%n     = 0.26412620D0
        va%m     = 0.03766668D0
        va%sig_u = 195.233128D0

        all_ok = .true.

        print *, "===================================================="
        print *, "   TEST: Fortran models vs Python reference values"
        print *, "   NRMSD < ", tolerance
        print *, "===================================================="

        do r = 1, nRates

            do i = 1, nPts
                pred_JC(i)  = jc%flow_stress(ep_test(i,r),  epd_rates(r))
                pred_RK(i)  = rk%flow_stress(ep_test(i,r),  epd_rates(r))
                pred_MRK(i) = mrk%flow_stress(ep_test(i,r), epd_rates(r))
                pred_NNL(i) = nnl%flow_stress(ep_test(i,r), epd_rates(r))
                pred_VA(i)  = va%flow_stress(ep_test(i,r),  epd_rates(r))
            end do

            errJC  = nrmsd(ref_JC(:,r),  pred_JC)
            errRK  = nrmsd(ref_RK(:,r),  pred_RK)
            errMRK = nrmsd(ref_MRK(:,r), pred_MRK)
            errNNL = nrmsd(ref_NNL(:,r), pred_NNL)
            errVA  = nrmsd(ref_VA(:,r),  pred_VA)

            okJC  = (errJC  < tolerance)
            okRK  = (errRK  < tolerance)
            okMRK = (errMRK < tolerance)
            okNNL = (errNNL < tolerance)
            okVA  = (errVA  < tolerance)

            all_ok = all_ok .and. okJC .and. okRK .and. okMRK .and. okNNL .and. okVA

            if (okJC)  then; stJC  = "JC OK  ";  else; stJC  = "JC FAIL"; end if
            if (okRK)  then; stRK  = "RK OK  ";  else; stRK  = "RK FAIL"; end if
            if (okMRK) then; stMRK = "MRK OK ";  else; stMRK = "MRKFAIL"; end if
            if (okNNL) then; stNNL = "NNL OK ";  else; stNNL = "NNLFAIL"; end if
            if (okVA)  then; stVA  = "VA OK  ";  else; stVA  = "VA FAIL"; end if

            print *
            write(*,'(A,1X,A,1X,A,ES12.4)') "Rate", trim(rate_name(r)), "epd=", epd_rates(r)
            write(*,'(A,1X,A,2X,A,2X,A,2X,A,2X,A)') "Status:", stJC, stRK, stMRK, stNNL, stVA
            write(*,'(A,1X,ES12.4,2X,A,1X,ES12.4,2X,A,1X,ES12.4,2X,A,1X,ES12.4,2X,A,1X,ES12.4)') &
                "NRMSD:", errJC, "RK", errRK, "MRK", errMRK, "NNL", errNNL, "VA", errVA

            if (.not.(okJC .and. okRK .and. okMRK .and. okNNL .and. okVA)) then
                print *
                print *, "---------------------------------------------------------------------------------------------------------"
                print *, " JC_ref  JC_F    RK_ref  RK_F    MRK_ref MRK_F   NNL_ref NNL_F   VA_ref VA_F"
                print *, "---------------------------------------------------------------------------------------------------------"
                do i = 1, nPts
                    write(*,'(F8.2,1X,F8.2,1X,F8.2,1X,F8.2,1X,F8.2,1X,F8.2,1X,F8.2,1X,F8.2,1X,F8.2,1X,F8.2)') &
                        ref_JC(i,r), pred_JC(i), &
                        ref_RK(i,r), pred_RK(i), &
                        ref_MRK(i,r), pred_MRK(i), &
                        ref_NNL(i,r), pred_NNL(i), &
                        ref_VA(i,r), pred_VA(i)
                end do
                print *, "---------------------------------------------------------------------------------------------------------"
            end if

        end do

        if (all_ok) then
            stFINAL = "PASSED      "
        else
            stFINAL = "FAILED "
        end if

        print *
        print *, "==============================================================="
        write(*,'(A,1X,A)') "FINAL:", trim(stFINAL)
        print *, "==============================================================="

        passed = all_ok

    end subroutine test_viscoplastic

end module test_viscoplastic_mod


program test_viscoplastic_main
    use test_viscoplastic_mod
    implicit none
    logical :: passed

    call test_viscoplastic(passed)

    if (.not. passed) stop 1
    stop 0
end program test_viscoplastic_main


