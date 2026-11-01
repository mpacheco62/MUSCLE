! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_plastic_history
    !! Module muscle_plastic_history
    !! =============================
    !! Provides persistent, transaction-aware material state management for plasticity (historical
    !! OO API, host only; removed in Phase 7 with the OO adapters, D-017.1).
    !!
    !! The states are the kernel type `plastic_state_dp` of muscle_plastic_state (historical name
    !! `Plastic_state`, same components and defaults); `pack_to_fea`/`unpack_from_fea` call its kernel
    !! procedures `plastic_state_pack`/`plastic_state_unpack` (FEA history array layout there) after
    !! checking the size of the array. Kernels keep the two states as separate arguments instead of
    !! this container.

    use muscle_kinds, only : dp
    use muscle_tensors, only : ten_3D2Osym
    use muscle_plastic_state, only : PLASTIC_STATE_SIZE, Plastic_state => plastic_state_dp, &
        plastic_state_pack, plastic_state_unpack
    implicit none
    private

    public :: Plastic_state
    public :: Plastic_material_history

    type :: Plastic_material_history
        !! Transactional history container holding converged state (t_n) and trial state (t_n+1).
        type(Plastic_state) :: state_n   !! Converged state at beginning of step t_n (Read-Only during NR)
        type(Plastic_state) :: state_np1 !! Trial/Candidate state at end of step t_n+1 (Updated during NR)
    contains
        procedure :: init            => history_init
        procedure :: commit          => history_commit
        procedure :: rollback        => history_rollback
        procedure :: unpack_from_fea => history_unpack
        procedure :: pack_to_fea     => history_pack
    end type Plastic_material_history

contains

    pure subroutine history_init(self, strain_p, strain_pf, stress)
        !! Initializes material history state for both t_n and t_n+1.
        class(Plastic_material_history), intent(inout) :: self
        type(ten_3D2Osym), intent(in), optional        :: strain_p
        real(dp), intent(in), optional                 :: strain_pf
        type(ten_3D2Osym), intent(in), optional        :: stress

        if (present(strain_p)) then
            self%state_n%strain_p   = strain_p
            self%state_np1%strain_p = strain_p
        end if

        if (present(strain_pf)) then
            self%state_n%strain_pf   = strain_pf
            self%state_np1%strain_pf = strain_pf
        end if

        if (present(stress)) then
            self%state_n%stress   = stress
            self%state_np1%stress = stress
        end if
    end subroutine history_init

    pure subroutine history_commit(self)
        !! Commits trial state_np1 to converged state_n upon global FEA convergence.
        class(Plastic_material_history), intent(inout) :: self
        self%state_n = self%state_np1
    end subroutine history_commit

    pure subroutine history_rollback(self)
        !! Discards candidate state_np1 and resets it to state_n if global FEA step fails.
        class(Plastic_material_history), intent(inout) :: self
        self%state_np1 = self%state_n
    end subroutine history_rollback

    pure subroutine history_unpack(self, hsv)
        !! Unpacks state_n from a flat 1D FEA history array (e.g. ANSYS / LS-DYNA statev/hsv array).
        !! Standard Mapping:
        !! - hsv(1:6)  : Constitutive Stress (xx, yy, zz, xy, yz, xz)
        !! - hsv(7:12) : Plastic Strain (xx, yy, zz, xy, yz, xz)
        !! - hsv(13)   : Equivalent Plastic Strain
        class(Plastic_material_history), intent(inout) :: self
        real(dp), intent(in)                           :: hsv(:)

        if (size(hsv) < PLASTIC_STATE_SIZE) then
            error stop "ERROR FATAL [Plastic_material_history%unpack_from_fea]: " // &
                       "El arreglo 'hsv' (statev) debe tener al menos 13 componentes."
        end if

        call plastic_state_unpack(hsv(1:PLASTIC_STATE_SIZE), self%state_n)

        ! Initialize candidate state_np1 to match state_n at step start
        self%state_np1 = self%state_n
    end subroutine history_unpack

    pure subroutine history_pack(self, hsv)
        !! Packs trial state_np1 into a flat 1D FEA history array for solver output.
        class(Plastic_material_history), intent(in) :: self
        real(dp), intent(inout)                     :: hsv(:)

        if (size(hsv) < PLASTIC_STATE_SIZE) then
            error stop "ERROR FATAL [Plastic_material_history%pack_to_fea]: " // &
                       "El arreglo 'hsv' (statev) debe tener al menos 13 componentes."
        end if

        call plastic_state_pack(self%state_np1, hsv(1:PLASTIC_STATE_SIZE))
    end subroutine history_pack

end module muscle_plastic_history