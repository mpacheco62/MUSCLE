! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_visco_state
    !! Module muscle_visco_state
    !! =========================
    !! Persistent, transaction-aware state of a viscoelastic material point.
    !!
    !! The state has a fixed size (`MAX_VISCO_BRANCHES` history tensors), so copying a history
    !! (commit, rollback, private copies in OpenMP regions) never allocates memory. The meaning
    !! of the history tensors depends on the model: for the generalized Maxwell model they are
    !! the deviatoric stresses of each Maxwell branch.
    !!
    !! FEA history array layout (`pack_to_fea` / `unpack_from_fea`):
    !!
    !! - n_branches            : Number of history tensors in use (-1: state not initialized)
    !! - hsv(1:6)              : Elastic strain at the end of the step (xx, yy, zz, xy, yz, xz)
    !! - hsv(6*k+1 : 6*k+6)    : History tensor of branch k (k = 1..n_branches)

    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    implicit none
    private

    integer, parameter, public :: MAX_VISCO_BRANCHES = 10
        !! Maximum number of branches (history tensors) of a viscoelastic state

    public :: Viscoelastic_state
    public :: Viscoelastic_history

    type :: Viscoelastic_state
        !! Viscoelastic state of a material point at a specific instant
        integer           :: n_branches = -1                                               !! Number of active branches
        type(ten_3D2Osym) :: strain_e = ten_3D2Osym((/0d0, 0d0, 0d0, 0d0, 0d0, 0d0/))  !! Elastic strain
        real(real64)      :: h(6, MAX_VISCO_BRANCHES) = 0D0                                !! History tensors (Voigt)
    end type Viscoelastic_state

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

        if (n_branches < 0 .or. n_branches > MAX_VISCO_BRANCHES) then
            error stop "ERROR FATAL [Viscoelastic_history%init]: " // &
                       "n_branches debe estar entre 0 y MAX_VISCO_BRANCHES."
        end if

        self%state_n%n_branches = n_branches
        self%state_n%h = 0D0
        self%state_n%strain_e%vals = 0D0
        if (present(strain_e)) self%state_n%strain_e = strain_e
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
        !! Unpacks state_n from a flat 1D FEA history array. `n_branches` must be set (`init`);
        !! a state that was never initialized stops with an error instead of reading no history.
        class(Viscoelastic_history), intent(inout) :: self
        real(real64), intent(in)                   :: hsv(:)
        integer :: n

        n = self%state_n%n_branches
        if (n < 0) then
            error stop "ERROR FATAL [Viscoelastic_history%unpack_from_fea]: " // &
                       "El estado no fue iniciado: llame a 'init' antes de 'unpack_from_fea'."
        end if
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
        real(real64), intent(inout)             :: hsv(:)
        integer :: n

        n = self%state_np1%n_branches
        if (n < 0) then
            error stop "ERROR FATAL [Viscoelastic_history%pack_to_fea]: " // &
                       "El estado no fue iniciado: llame a 'init' antes de 'pack_to_fea'."
        end if
        if (size(hsv) < 6*(n + 1)) then
            error stop "ERROR FATAL [Viscoelastic_history%pack_to_fea]: " // &
                       "El arreglo 'hsv' debe tener al menos 6*(n_branches + 1) componentes."
        end if

        hsv(1:6) = self%state_np1%strain_e%vals
        hsv(7:6*(n + 1)) = reshape(self%state_np1%h(:, 1:n), [6*n])
    end subroutine history_pack

end module muscle_visco_state
