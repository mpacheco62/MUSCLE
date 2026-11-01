! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_kin_finite_base
    !! Historical OO base of the kinematics driven by the total deformation gradient F (host only,
    !! removed in Phase 7 with Base_kinematics, D-017.1). The kinematic state is cached in the
    !! historical components (`F`, `F_inv`, `R_rot`, `J`, `C_mat`, `b_spat`, `b_inv`) and every
    !! binding calls the dp kernel procedures of muscle_kin_finite (D-017.14), ignoring their
    !! `status`: same results as the historical code, also for an invalid F (D-017.14a).
    !!
    !! `finite_state_of(self)` and `finite_set_state(self, state)` copy the cached state to and
    !! from a `kin_state_dp` (for the adapters of the concrete kinematics). Documented, not changed
    !! (D-017.14b): `update_incremental` builds F_{n+1} = (drot + dstrain) F_n (additive) and ignores
    !! `dt`; `update_F_base` does not update the strain of the concrete type; the tangent
    !! transformations ignore their optional stress argument.
    use muscle_kinds, only : dp
    use muscle_tensors
    use muscle_kinematics_base
    use muscle_kin_finite, only : kin_state_dp, kin_finite_update_dp, kin_green_lagrange_dp, &
        kin_euler_almansi_dp, kin_material_log_dp, kin_spatial_log_dp, kin_pk2_to_cauchy_dp, &
        kin_cauchy_to_pk2_dp, kin_pk2_to_kirchhoff_dp, kin_kirchhoff_to_pk2_dp, kin_cauchy_to_kirchhoff_dp, &
        kin_kirchhoff_to_cauchy_dp, kin_pk2_to_pk1_dp, kin_pk1_to_pk2_dp, kin_cauchy_to_pk1_dp, &
        kin_pk1_to_cauchy_dp, kin_push_forward_tangent_dp, kin_pull_back_tangent_dp
    implicit none
    private

    public :: Base_F_kinematics, finite_state_of, finite_set_state

    type, abstract, extends(Base_kinematics) :: Base_F_kinematics
        !! Abstract Base Class for Kinematics driven by total Deformation Gradient F
        type(ten_3D2O)    :: F          !! Total deformation gradient F
        type(ten_3D2O)    :: F_inv      !! Inverse deformation gradient F^-1 (Cached)
        type(ten_3D2O)    :: R_rot      !! Rotation tensor R from polar decomp F = R*U (Cached)
        real(dp)          :: J = 1.0_dp !! Jacobian det(F) (Cached)
        type(ten_3D2Osym) :: C_mat      !! Right Cauchy-Green tensor C = F^T * F (Cached)
        type(ten_3D2Osym) :: b_spat     !! Left Cauchy-Green tensor b = F * F^T (Cached)
        type(ten_3D2Osym) :: b_inv      !! Inverse Left Cauchy-Green tensor b^-1 = F^-T * F^-1 (Cached)
    contains
        procedure :: update_F_base      => finite_update_F
        procedure :: update_F           => finite_update_F
        procedure :: update_incremental => finite_update_incremental
        generic   :: update             => update_F, update_incremental

        procedure :: jacobian           => finite_jacobian
        procedure :: is_finite_strain   => finite_is_finite_strain

        ! Default Base_kinematics Push/Pull Overrides
        procedure :: push_forward_stress => PK2_to_Cauchy
        procedure :: pull_back_stress    => Cauchy_to_PK2
        procedure :: push_forward_tangent=> finite_push_tangent
        procedure :: pull_back_tangent  => finite_pull_back_tangent

        ! --- Explicit Strain Measure Getters ---
        procedure :: get_RightCauchyGreen   => finite_get_C
        procedure :: get_LeftCauchyGreen_inv=> finite_get_binv
        procedure :: get_GreenLagrange      => finite_get_E_GL
        procedure :: get_EulerAlmansi       => finite_get_e_EA
        procedure :: get_Material_Logarithmic => finite_get_H_mat
        procedure :: get_Spatial_Logarithmic  => finite_get_h_spat
        procedure :: get_Rotation             => finite_get_R

        ! --- Direct Compile-Time Bound Stress Conversions ---
        procedure :: convert_PK2_to_Cauchy        => PK2_to_Cauchy
        procedure :: convert_Cauchy_to_PK2        => Cauchy_to_PK2
        procedure :: convert_PK2_to_Kirchhoff     => PK2_to_Kirchhoff
        procedure :: convert_Kirchhoff_to_PK2     => Kirchhoff_to_PK2
        procedure :: convert_Cauchy_to_Kirchhoff  => Cauchy_to_Kirchhoff
        procedure :: convert_Kirchhoff_to_Cauchy  => Kirchhoff_to_Cauchy
        procedure :: convert_PK2_to_PK1           => PK2_to_PK1
        procedure :: convert_PK1_to_PK2           => PK1_to_PK2
        procedure :: convert_Cauchy_to_PK1        => Cauchy_to_PK1
        procedure :: convert_PK1_to_Cauchy        => PK1_to_Cauchy
    end type Base_F_kinematics

contains

    pure function finite_state_of(self) result(state)
        !! Kernel state (`kin_state_dp`) of the cached components; zero strain.
        class(Base_F_kinematics), intent(in) :: self
            !! Kinematics.
        type(kin_state_dp) :: state
            !! Kinematic state.

        state%F = self%F
        state%F_inv = self%F_inv
        state%R = self%R_rot
        state%J = self%J
        state%C = self%C_mat
        state%b = self%b_spat
        state%b_inv = self%b_inv
    end function finite_state_of

    pure subroutine finite_set_state(self, state)
        !! Caches a kernel state (`kin_state_dp`) in the historical components (not the strain).
        class(Base_F_kinematics), intent(inout) :: self
            !! Kinematics.
        type(kin_state_dp), intent(in) :: state
            !! Kinematic state.

        self%F = state%F
        self%F_inv = state%F_inv
        self%R_rot = state%R
        self%J = state%J
        self%C_mat = state%C
        self%b_spat = state%b
        self%b_inv = state%b_inv
    end subroutine finite_set_state

    pure subroutine finite_update_F(self, F)
        !! Caches the kinematic state of F (`kin_finite_update_dp`; its status is ignored, D-017.14a).
        class(Base_F_kinematics), intent(inout) :: self
        type(ten_3D2O), intent(in)              :: F
        type(kin_state_dp) :: state
        integer :: status

        call kin_finite_update_dp(F, state, status)
        call finite_set_state(self, state)
    end subroutine finite_update_F

    pure subroutine finite_update_incremental(self, dstrain, drot, dt)
        class(Base_F_kinematics), intent(inout) :: self
        type(ten_3D2Osym), intent(in)            :: dstrain
        type(ten_3D2O), intent(in)                 :: drot
        real(dp), intent(in), optional             :: dt
        type(ten_3D2O) :: F_inc
        F_inc = drot + dstrain
        call self%update_F(F_inc * self%F)
    end subroutine finite_update_incremental

    pure function finite_jacobian(self) result(J)
        class(Base_F_kinematics), intent(in) :: self
        real(dp)                             :: J
        J = self%J
    end function finite_jacobian

    pure function finite_is_finite_strain(self) result(res)
        class(Base_F_kinematics), intent(in) :: self
        logical                              :: res
        res = .TRUE.
    end function finite_is_finite_strain

    pure function finite_get_C(self) result(C_tensor)
        class(Base_F_kinematics), intent(in) :: self
        type(ten_3D2Osym)                     :: C_tensor
        C_tensor = self%C_mat
    end function finite_get_C

    pure function finite_get_binv(self) result(binv_tensor)
        class(Base_F_kinematics), intent(in) :: self
        type(ten_3D2Osym)                     :: binv_tensor
        binv_tensor = self%b_inv
    end function finite_get_binv

    pure function finite_get_E_GL(self) result(E_gl)
        class(Base_F_kinematics), intent(in) :: self
        type(ten_3D2Osym)                     :: E_gl
        E_gl = kin_green_lagrange_dp(finite_state_of(self))
    end function finite_get_E_GL

    pure function finite_get_e_EA(self) result(e_ea)
        class(Base_F_kinematics), intent(in) :: self
        type(ten_3D2Osym)                     :: e_ea
        e_ea = kin_euler_almansi_dp(finite_state_of(self))
    end function finite_get_e_EA

    pure function finite_get_E_log(self) result(E_log)
        !! Unused private procedure, kept (D-017.14b): the material Hencky strain.
        class(Base_F_kinematics), intent(in) :: self
        type(ten_3D2Osym)                     :: E_log
        E_log = kin_material_log_dp(finite_state_of(self))
    end function finite_get_E_log

    pure function finite_get_b(self) result(b_tensor)
        !! Unused private procedure, kept (D-017.14b): the left Cauchy-Green tensor.
        class(Base_F_kinematics), intent(in) :: self
        type(ten_3D2Osym)                     :: b_tensor
        b_tensor = self%b_spat
    end function finite_get_b

    pure function finite_get_R(self) result(R_tensor)
        class(Base_F_kinematics), intent(in) :: self
        type(ten_3D2O)                        :: R_tensor
        R_tensor = self%R_rot
    end function finite_get_R

    pure function finite_get_H_mat(self) result(H_mat)
        !! Material Hencky strain H = 0.5 * ln(C) = ln(U)
        class(Base_F_kinematics), intent(in) :: self
        type(ten_3D2Osym)                     :: H_mat
        H_mat = kin_material_log_dp(finite_state_of(self))
    end function finite_get_H_mat

    pure function finite_get_h_spat(self) result(h_spat)
        !! Spatial Hencky strain h = 0.5 * ln(b) = ln(V)
        class(Base_F_kinematics), intent(in) :: self
        type(ten_3D2Osym)                     :: h_spat
        h_spat = kin_spatial_log_dp(finite_state_of(self))
    end function finite_get_h_spat

    pure function PK2_to_Cauchy(self, S_material) result(sigma)
        class(Base_F_kinematics), intent(in) :: self
        type(ten_3D2Osym), intent(in)        :: S_material
        type(ten_3D2Osym)                    :: sigma
        sigma = kin_pk2_to_cauchy_dp(finite_state_of(self), S_material)
    end function PK2_to_Cauchy

    pure function Cauchy_to_PK2(self, sigma_spatial) result(S)
        class(Base_F_kinematics), intent(in) :: self
        type(ten_3D2Osym), intent(in)        :: sigma_spatial
        type(ten_3D2Osym)                    :: S
        S = kin_cauchy_to_pk2_dp(finite_state_of(self), sigma_spatial)
    end function Cauchy_to_PK2

    pure function PK2_to_Kirchhoff(self, S) result(tau)
        class(Base_F_kinematics), intent(in) :: self
        type(ten_3D2Osym), intent(in)        :: S
        type(ten_3D2Osym)                    :: tau
        tau = kin_pk2_to_kirchhoff_dp(finite_state_of(self), S)
    end function PK2_to_Kirchhoff

    pure function Kirchhoff_to_PK2(self, tau) result(S)
        class(Base_F_kinematics), intent(in) :: self
        type(ten_3D2Osym), intent(in)        :: tau
        type(ten_3D2Osym)                    :: S
        S = kin_kirchhoff_to_pk2_dp(finite_state_of(self), tau)
    end function Kirchhoff_to_PK2

    pure function Cauchy_to_Kirchhoff(self, sigma) result(tau)
        class(Base_F_kinematics), intent(in) :: self
        type(ten_3D2Osym), intent(in)        :: sigma
        type(ten_3D2Osym)                    :: tau
        tau = kin_cauchy_to_kirchhoff_dp(finite_state_of(self), sigma)
    end function Cauchy_to_Kirchhoff

    pure function Kirchhoff_to_Cauchy(self, tau) result(sigma)
        class(Base_F_kinematics), intent(in) :: self
        type(ten_3D2Osym), intent(in)        :: tau
        type(ten_3D2Osym)                    :: sigma
        sigma = kin_kirchhoff_to_cauchy_dp(finite_state_of(self), tau)
    end function Kirchhoff_to_Cauchy

    pure function PK2_to_PK1(self, S) result(P)
        class(Base_F_kinematics), intent(in) :: self
        type(ten_3D2Osym), intent(in)        :: S
        type(ten_3D2O)                       :: P
        P = kin_pk2_to_pk1_dp(finite_state_of(self), S)
    end function PK2_to_PK1

    pure function PK1_to_PK2(self, P) result(S)
        class(Base_F_kinematics), intent(in) :: self
        type(ten_3D2O), intent(in)           :: P
        type(ten_3D2Osym)                    :: S
        S = kin_pk1_to_pk2_dp(finite_state_of(self), P)
    end function PK1_to_PK2

    pure function Cauchy_to_PK1(self, sigma) result(P)
        class(Base_F_kinematics), intent(in) :: self
        type(ten_3D2Osym), intent(in)        :: sigma
        type(ten_3D2O)                       :: P
        P = kin_cauchy_to_pk1_dp(finite_state_of(self), sigma)
    end function Cauchy_to_PK1

    pure function PK1_to_Cauchy(self, P) result(sigma)
        class(Base_F_kinematics), intent(in) :: self
        type(ten_3D2O), intent(in)           :: P
        type(ten_3D2Osym)                    :: sigma
        sigma = kin_pk1_to_cauchy_dp(finite_state_of(self), P)
    end function PK1_to_Cauchy

    pure function finite_push_tangent(self, C_material, stress_spatial) result(c_spatial)
        !! True 4th-order Push-Forward: c_ijkl = (1/J) * F_iI * F_jJ * F_kK * F_lL * C_IJKL
        class(Base_F_kinematics), intent(in)     :: self
        type(ten_3D4O2sym), intent(in)           :: C_material
        type(ten_3D2Osym), intent(in), optional  :: stress_spatial
        type(ten_3D4O2sym)                       :: c_spatial

        c_spatial = kin_push_forward_tangent_dp(finite_state_of(self), C_material)
    end function finite_push_tangent

    pure function finite_pull_back_tangent(self, c_spatial, stress_material) result(C_material)
        !! True 4th-order Pull-Back: C_IJKL = J * F^-1_Ii * F^-1_Jj * F^-1_Kk * F^-1_Ll * c_ijkl
        class(Base_F_kinematics), intent(in)     :: self
        type(ten_3D4O2sym), intent(in)           :: c_spatial
        type(ten_3D2Osym), intent(in), optional  :: stress_material
        type(ten_3D4O2sym)                       :: C_material

        C_material = kin_pull_back_tangent_dp(finite_state_of(self), c_spatial)
    end function finite_pull_back_tangent

end module muscle_kin_finite_base
