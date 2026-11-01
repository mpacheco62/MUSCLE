! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_visco_state
    !! Module muscle_visco_state
    !! =========================
    !! Persistent, transaction-aware state of a viscoelastic material point (historical OO API,
    !! host only; removed in Phase 7 with the OO adapters, D-017.1).
    !!
    !! The states are the kernel type `visco_state_dp` of muscle_visco_types (historical name
    !! `Viscoelastic_state`, same components and defaults): fixed size (`MAX_VISCO_BRANCHES` history
    !! tensors), so copying a history (commit, rollback, private copies in OpenMP regions) never
    !! allocates memory. The meaning of the history tensors depends on the model: for the
    !! generalized Maxwell model they are the deviatoric stresses of each Maxwell branch. Kernels
    !! keep the two states as separate arguments (`visco_update(prm, coefs, state_n, strain_e,
    !! state_np1)`) instead of this container.
    !!
    !! FEA history array layout (`pack_to_fea` / `unpack_from_fea`):
    !!
    !! - hsv(1:6)              : Elastic strain at the end of the step (xx, yy, zz, xy, yz, xz)
    !! - hsv(6*k+1 : 6*k+6)    : History tensor of branch k (k = 1..n_branches)
    !!
    !! `unpack_from_fea` needs `n_branches` set before (`init`, D-017.13c).

    use muscle_kinds, only : dp
    use muscle_tensors, only : ten_3D2Osym
    use muscle_status, only : STATUS_OK
    use muscle_visco_types, only : MAX_VISCO_BRANCHES, Viscoelastic_state => visco_state_dp, visco_state_init
    implicit none
    private

    public :: MAX_VISCO_BRANCHES
    public :: Viscoelastic_state
    public :: Viscoelastic_history

    type :: Viscoelastic_history
        !! Transactional history container holding converged state (t_n) and trial state (t_n+1).
        type(Viscoelastic_state) :: state_n   !! Converged state at beginning of step t_n (Read-Only during NR)
        type(Viscoelastic_state) :: state_np1 !! Trial/Candidate state at end of step t_n+1
    contains
        procedure :: init            => history_init
        procedure :: commit          => history_commit
        procedure :: rollback        => history_rollback
        procedure :: unpack_from_fea => history_unpack
        procedure :: pack_to_fea     => history_pack
    end type Viscoelastic_history

contains

    pure subroutine history_init(self, n_branches, strain_e)
        !! Initializes both states with `n_branches` zero history tensors.
        class(Viscoelastic_history), intent(inout) :: self
        integer, intent(in)                        :: n_branches
        type(ten_3D2Osym), intent(in), optional    :: strain_e
        integer :: status

        call visco_state_init(self%state_n, n_branches, status, strain_e)
        if (status /= STATUS_OK) then
            error stop "ERROR FATAL [Viscoelastic_history%init]: " // &
                       "n_branches must be between 0 and MAX_VISCO_BRANCHES."
        end if
        self%state_np1 = self%state_n
    end subroutine history_init

    pure subroutine history_commit(self)
        !! Commits trial state_np1 to converged state_n upon global FEA convergence.
        class(Viscoelastic_history), intent(inout) :: self
        self%state_n = self%state_np1
    end subroutine history_commit

    pure subroutine history_rollback(self)
        !! Discards candidate state_np1 and resets it to state_n if global FEA step fails.
        class(Viscoelastic_history), intent(inout) :: self
        self%state_np1 = self%state_n
    end subroutine history_rollback

    pure subroutine history_unpack(self, hsv)
        !! Unpacks state_n from a flat 1D FEA history array. `n_branches` must be set (`init`).
        class(Viscoelastic_history), intent(inout) :: self
        real(dp), intent(in)                       :: hsv(:)
        integer :: n

        n = self%state_n%n_branches
        if (size(hsv) < 6*(n + 1)) then
            error stop "ERROR FATAL [Viscoelastic_history%unpack_from_fea]: " // &
                       "El arreglo 'hsv' debe tener al menos 6*(n_branches + 1) componentes."
        end if

        self%state_n%strain_e%vals = hsv(1:6)
        self%state_n%h(:, 1:n) = reshape(hsv(7:6*(n + 1)), [6, n])

        ! Initialize candidate state_np1 to match state_n at step start
        self%state_np1 = self%state_n
    end subroutine history_unpack

    pure subroutine history_pack(self, hsv)
        !! Packs trial state_np1 into a flat 1D FEA history array.
        class(Viscoelastic_history), intent(in) :: self
        real(dp), intent(inout)                 :: hsv(:)
        integer :: n

        n = self%state_np1%n_branches
        if (size(hsv) < 6*(n + 1)) then
            error stop "ERROR FATAL [Viscoelastic_history%pack_to_fea]: " // &
                       "El arreglo 'hsv' debe tener al menos 6*(n_branches + 1) componentes."
        end if

        hsv(1:6) = self%state_np1%strain_e%vals
        hsv(7:6*(n + 1)) = reshape(self%state_np1%h(:, 1:n), [6*n])
    end subroutine history_pack

end module muscle_visco_state
