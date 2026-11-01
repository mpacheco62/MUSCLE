! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module plastic_history_golden_cases
    !! Cases of the golden test of muscle_plastic_history: the defaults of `Plastic_state` and
    !! `Plastic_material_history`, every combination of the optional arguments of `init` (on a history
    !! whose committed and candidate states differ), `commit`, `rollback`, `unpack_from_fea` and
    !! `pack_to_fea` (arrays of exactly 13 and of more components, the extra ones untouched), and
    !! multi-step transactional sequences (unpack, candidate updates, pack, commit or rollback).
    use muscle_kinds, only : dp
    use muscle_tensors, only : ten_3D2Osym
    use muscle_plastic_history, only : Plastic_state, Plastic_material_history
    use muscle_test_golden, only : golden_check
    implicit none
    private

    public :: tensor_of, state_of, record_history

contains

    function tensor_of(seed) result(t)
        !! Deterministic tensor with six distinct, non-round components.
        real(dp), intent(in) :: seed
        type(ten_3D2Osym) :: t
        integer :: i
        t%vals = [(seed*sin(1.3_dp*i + seed) + 0.1_dp*i, i = 1, 6)]
    end function tensor_of

    function state_of(seed) result(s)
        !! Deterministic plastic state.
        real(dp), intent(in) :: seed
        type(Plastic_state) :: s
        s%stress = tensor_of(100.0_dp*seed)
        s%strain_p = tensor_of(1.0e-2_dp*seed)
        s%strain_pf = 3.0e-2_dp*seed + 1.0_dp/7.0_dp
    end function state_of

    subroutine record_history(key, h)
        !! Records both states of a history.
        character(len=*), intent(in) :: key
        type(Plastic_material_history), intent(in) :: h
        call golden_check(key // ' | state_n stress, strain_p, strain_pf; state_np1 stress, strain_p, strain_pf', &
            [h%state_n%stress%vals, h%state_n%strain_p%vals, h%state_n%strain_pf, &
             h%state_np1%stress%vals, h%state_np1%strain_p%vals, h%state_np1%strain_pf])
    end subroutine record_history

end module plastic_history_golden_cases


program test_plastic_history_golden
    !! Golden (characterization) test of muscle_plastic_history, double precision, bit for bit
    !! against tests/reference/plasticity/plastic_history_dp.txt (captured from the original OO code;
    !! see regenerate.sh there). Exit codes: those of `golden_finish` (1-4). The `error stop` of
    !! `unpack_from_fea`/`pack_to_fea` with fewer than 13 components cannot be recorded here.
    use muscle_kinds, only : dp
    use muscle_tensors, only : ten_3D2Osym, operator(+)
    use muscle_plastic_history, only : Plastic_state, Plastic_material_history
    use muscle_test_golden, only : golden_init, golden_check, golden_finish
    use plastic_history_golden_cases
    implicit none

    type(Plastic_material_history) :: h
    real(dp) :: hsv13(13), hsv16(16)
    integer :: mask, istep
    character(len=64) :: buf

    call golden_init('plasticity/plastic_history_dp.txt')

    ! --- Defaults ---
    block
        type(Plastic_state) :: s0
        type(Plastic_material_history) :: h0
        call golden_check('defaults | Plastic_state stress, strain_p, strain_pf', &
            [s0%stress%vals, s0%strain_p%vals, s0%strain_pf])
        call record_history('defaults | history', h0)
    end block

    ! --- init: every combination of the optional arguments, on differing committed/candidate states ---
    do mask = 0, 7
        block
            type(Plastic_material_history) :: hi
            hi%state_n = state_of(1.0_dp)
            hi%state_np1 = state_of(2.0_dp)
            select case (mask)
            case (0); call hi%init()
            case (1); call hi%init(strain_p=tensor_of(-0.03_dp))
            case (2); call hi%init(strain_pf=0.0625_dp)
            case (3); call hi%init(strain_p=tensor_of(-0.03_dp), strain_pf=0.0625_dp)
            case (4); call hi%init(stress=tensor_of(-250.0_dp))
            case (5); call hi%init(strain_p=tensor_of(-0.03_dp), stress=tensor_of(-250.0_dp))
            case (6); call hi%init(strain_pf=0.0625_dp, stress=tensor_of(-250.0_dp))
            case default; call hi%init(tensor_of(-0.03_dp), 0.0625_dp, tensor_of(-250.0_dp))
            end select
            write(buf, '(a,i0)') 'init mask ', mask
            call record_history(trim(buf), hi)
        end block
    end do
    block
        type(Plastic_material_history) :: hi
        call hi%init(strain_pf=0.5_dp)
        call record_history('init on default history', hi)
    end block

    ! --- commit / rollback on differing states ---
    block
        type(Plastic_material_history) :: hc
        hc%state_n = state_of(3.0_dp)
        hc%state_np1 = state_of(4.0_dp)
        call hc%commit()
        call record_history('commit', hc)
        hc%state_n = state_of(3.0_dp)
        hc%state_np1 = state_of(4.0_dp)
        call hc%rollback()
        call record_history('rollback', hc)
    end block

    ! --- unpack (13 and 16 components) and pack (extra components untouched) ---
    hsv13 = [(0.37_dp*istep - 1.1_dp*sin(real(istep, dp)), istep = 1, 13)]
    hsv16 = [(-0.53_dp*istep + 2.3_dp*cos(real(istep, dp)), istep = 1, 16)]
    h%state_n = state_of(5.0_dp)
    h%state_np1 = state_of(6.0_dp)
    call h%unpack_from_fea(hsv13)
    call record_history('unpack 13', h)
    h%state_np1 = state_of(7.0_dp)
    call h%unpack_from_fea(hsv16)
    call record_history('unpack 16', h)
    h%state_n = state_of(8.0_dp)
    h%state_np1 = state_of(9.0_dp)
    hsv13 = -1.0_dp
    call h%pack_to_fea(hsv13)
    call golden_check('pack 13 | hsv', hsv13)
    hsv16 = [(real(istep, dp)/3.0_dp, istep = 1, 16)]
    call h%pack_to_fea(hsv16)
    call golden_check('pack 16 | hsv', hsv16)
    call record_history('pack leaves history', h)

    ! --- Multi-step transactional sequence (as an FE driver: unpack, iterate, pack, commit/rollback) ---
    call h%init(tensor_of(0.0_dp), 0.0_dp, tensor_of(0.0_dp))
    hsv16 = 0.0_dp
    do istep = 1, 6
        write(buf, '(a,i0)') 'sequence step ', istep
        call h%unpack_from_fea(hsv16)
        call record_history(trim(buf) // ' unpack', h)
        h%state_np1%stress = h%state_n%stress + tensor_of(10.0_dp*istep)
        h%state_np1%strain_p = h%state_n%strain_p + tensor_of(1.0e-3_dp*istep)
        h%state_np1%strain_pf = h%state_n%strain_pf + 1.0e-3_dp*istep
        if (mod(istep, 3) == 0) then
            call h%rollback()
            call record_history(trim(buf) // ' rollback', h)
            h%state_np1%strain_pf = h%state_n%strain_pf + 5.0e-4_dp*istep
        end if
        call h%pack_to_fea(hsv16)
        call golden_check(trim(buf) // ' pack | hsv', hsv16)
        call h%commit()
        call record_history(trim(buf) // ' commit', h)
    end do

    call golden_finish()
end program test_plastic_history_golden
