! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_kin_finite_base
    use, intrinsic :: iso_fortran_env, only : real64
    use muscle_tensors
    use muscle_kinematics_base
    use muscle_math_operations, only : eigenvals
    use muscle_math_spectral_derivs, only : dEigenvalues_dTensor
    implicit none
    private

    public :: Base_F_kinematics

    type, abstract, extends(Base_kinematics) :: Base_F_kinematics
        !! Abstract Base Class for Kinematics driven by total Deformation Gradient F
        type(ten_3D2O)    :: F          !! Total deformation gradient F
        type(ten_3D2O)    :: F_inv      !! Inverse deformation gradient F^-1 (Cached)
        type(ten_3D2O)    :: R_rot      !! Rotation tensor R from polar decomp F = R*U (Cached)
        real(real64)      :: J = 1.0D0  !! Jacobian det(F) (Cached)
        type(ten_3D2Osym) :: C_mat      !! Right Cauchy-Green tensor C = F^T * F (Cached)
        type(ten_3D2Osym) :: b_spat     !! Left Cauchy-Green tensor b = F * F^T (Cached)
        type(ten_3D2Osym) :: b_inv      !! Inverse Left Cauchy-Green tensor b^-1 = F^-T * F^-1 (Cached)
    contains
        procedure :: update_F_base      => finite_update_F
        procedure :: update_F           => finite_update_F
        procedure :: update_incremental => finite_update_incremental
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

    pure subroutine finite_update_F(self, F)
        class(Base_F_kinematics), intent(inout) :: self
        type(ten_3D2O), intent(in)              :: F
        real(real64) :: F11, F21, F31, F12, F22, F32, F13, F23, F33
        real(real64) :: G11, G21, G31, G12, G22, G32, G13, G23, G33
        real(real64) :: inv_J, c11, c22, c33, c12, c23, c13
        real(real64) :: b11, b22, b33, b12, b23, b13
        real(real64) :: binv11, binv22, binv33, binv12, binv23, binv13
        real(real64) :: eigC(3)
        type(ten_3D2Osym) :: H_projectors(3), U_inv
        integer :: a

        self%F = F
        self%J = F%det()
        inv_J  = 1.0D0 / self%J

        ! 1. Unpack F
        F11 = F%vals(1); F21 = F%vals(2); F31 = F%vals(3)
        F12 = F%vals(4); F22 = F%vals(5); F32 = F%vals(6)
        F13 = F%vals(7); F23 = F%vals(8); F33 = F%vals(9)

        ! 2. Compute F_inv
        G11 = (F22*F33 - F23*F32)*inv_J; G12 = (F13*F32 - F12*F33)*inv_J; G13 = (F12*F23 - F13*F22)*inv_J
        G21 = (F23*F31 - F21*F33)*inv_J; G22 = (F11*F33 - F13*F31)*inv_J; G23 = (F13*F21 - F11*F23)*inv_J
        G31 = (F21*F32 - F22*F31)*inv_J; G32 = (F12*F31 - F11*F32)*inv_J; G33 = (F11*F22 - F12*F21)*inv_J

        call self%F_inv%init( &
            xx = G11, xy = G12, xz = G13, &
            yx = G21, yy = G22, yz = G23, &
            zx = G31, zy = G32, zz = G33  &
        )

        ! 3. Compute C = F^T * F
        c11 = F11*F11 + F21*F21 + F31*F31; c22 = F12*F12 + F22*F22 + F32*F32; c33 = F13*F13 + F23*F23 + F33*F33
        c12 = F11*F12 + F21*F22 + F31*F32; c23 = F12*F13 + F22*F23 + F32*F33; c13 = F11*F13 + F21*F23 + F31*F33
        call self%C_mat%init(xx=c11, yy=c22, zz=c33, xy=c12, yz=c23, xz=c13)

        ! 4. Compute b = F * F^T
        b11 = F11*F11 + F12*F12 + F13*F13; b22 = F21*F21 + F22*F22 + F23*F23; b33 = F31*F31 + F32*F32 + F33*F33
        b12 = F11*F21 + F12*F22 + F13*F23; b23 = F21*F31 + F22*F32 + F23*F33; b13 = F11*F31 + F12*F32 + F13*F33
        call self%b_spat%init(xx=b11, yy=b22, zz=b33, xy=b12, yz=b23, xz=b13)

        ! 5. Compute b^-1 = F^-T * F^-1
        binv11 = G11*G11 + G21*G21 + G31*G31; binv22 = G12*G12 + G22*G22 + G32*G32; binv33 = G13*G13 + G23*G23 + G33*G33
        binv12 = G11*G12 + G21*G22 + G31*G32; binv23 = G12*G13 + G22*G23 + G32*G33; binv13 = G11*G13 + G21*G23 + G31*G33
        call self%b_inv%init(xx=binv11, yy=binv22, zz=binv33, xy=binv12, yz=binv23, xz=binv13)
        
        ! 6. Analytic Rotation R = F * U^-1 via spectral decomposition of C
        eigC = eigenvals(self%C_mat)
        H_projectors = dEigenvalues_dTensor(self%C_mat, eigC)
        U_inv = 0.0D0
        do a = 1, 3
            if (eigC(a) > 1.0D-12) then
                U_inv = U_inv + (1.0D0 / sqrt(eigC(a))) * H_projectors(a)
            end if
        end do
        self%R_rot = F * U_inv
    end subroutine finite_update_F

    pure subroutine finite_update_incremental(self, dstrain, drot, dt)
        class(Base_F_kinematics), intent(inout) :: self
        type(ten_3D2Osym), intent(in)            :: dstrain
        type(ten_3D2O), intent(in)                 :: drot
        real(real64), intent(in), optional         :: dt
        type(ten_3D2O) :: F_inc
        F_inc = drot + dstrain
        call self%update_F(F_inc * self%F)
    end subroutine finite_update_incremental

    pure function finite_jacobian(self) result(J)
        class(Base_F_kinematics), intent(in) :: self
        real(real64)                         :: J
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

        call E_gl%init( &
            xx = 0.5D0 * (self%C_mat%vals(1) - 1.0D0), &
            yy = 0.5D0 * (self%C_mat%vals(2) - 1.0D0), &
            zz = 0.5D0 * (self%C_mat%vals(3) - 1.0D0), &
            xy = 0.5D0 * self%C_mat%vals(4),           &
            yz = 0.5D0 * self%C_mat%vals(5),           &
            xz = 0.5D0 * self%C_mat%vals(6)            &
        )
    end function finite_get_E_GL

    pure function finite_get_e_EA(self) result(e_ea)
        class(Base_F_kinematics), intent(in) :: self
        type(ten_3D2Osym)                     :: e_ea

        call e_ea%init( &
            xx = 0.5D0 * (1.0D0 - self%b_inv%vals(1)), &
            yy = 0.5D0 * (1.0D0 - self%b_inv%vals(2)), &
            zz = 0.5D0 * (1.0D0 - self%b_inv%vals(3)), &
            xy = -0.5D0 * self%b_inv%vals(4),          &
            yz = -0.5D0 * self%b_inv%vals(5),          &
            xz = -0.5D0 * self%b_inv%vals(6)           &
        )
    end function finite_get_e_EA

    pure function finite_get_E_log(self) result(E_log)
        class(Base_F_kinematics), intent(in) :: self
        type(ten_3D2Osym)                     :: E_log
        real(real64)      :: eigC(3)
        type(ten_3D2Osym) :: H_projectors(3)
        integer           :: a

        eigC = eigenvals(self%C_mat)
        H_projectors = dEigenvalues_dTensor(self%C_mat, eigC)

        E_log = 0.0D0
        do a = 1, 3
            if (eigC(a) > 1.0D-12) then
                E_log = E_log + (0.5D0 * log(eigC(a))) * H_projectors(a)
            end if
        end do
    end function finite_get_E_log

    pure function finite_get_b(self) result(b_tensor)
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
        real(real64)      :: eigC(3)
        type(ten_3D2Osym) :: H_projectors(3)
        integer           :: a

        eigC = eigenvals(self%C_mat)
        H_projectors = dEigenvalues_dTensor(self%C_mat, eigC)

        H_mat = 0.0D0
        do a = 1, 3
            if (eigC(a) > 1.0D-12) then
                H_mat = H_mat + (0.5D0 * log(eigC(a))) * H_projectors(a)
            end if
        end do
    end function finite_get_H_mat

    pure function finite_get_h_spat(self) result(h_spat)
        !! Spatial Hencky strain h = 0.5 * ln(b) = ln(V)
        class(Base_F_kinematics), intent(in) :: self
        type(ten_3D2Osym)                     :: h_spat
        real(real64)      :: eigb(3)
        type(ten_3D2Osym) :: H_projectors(3)
        integer           :: a

        eigb = eigenvals(self%b_spat)
        H_projectors = dEigenvalues_dTensor(self%b_spat, eigb)

        h_spat = 0.0D0
        do a = 1, 3
            if (eigb(a) > 1.0D-12) then
                h_spat = h_spat + (0.5D0 * log(eigb(a))) * H_projectors(a)
            end if
        end do
    end function finite_get_h_spat

    pure function PK2_to_Cauchy(self, S_material) result(sigma)
        class(Base_F_kinematics), intent(in) :: self
        type(ten_3D2Osym), intent(in)        :: S_material
        type(ten_3D2Osym)                    :: sigma
        sigma = (self%F .transform. S_material) / self%J
    end function PK2_to_Cauchy

    pure function Cauchy_to_PK2(self, sigma_spatial) result(S)
        class(Base_F_kinematics), intent(in) :: self
        type(ten_3D2Osym), intent(in)        :: sigma_spatial
        type(ten_3D2Osym)                    :: S
        S = (self%F_inv .transform. sigma_spatial) * self%J
    end function Cauchy_to_PK2

    pure function PK2_to_Kirchhoff(self, S) result(tau)
        class(Base_F_kinematics), intent(in) :: self
        type(ten_3D2Osym), intent(in)        :: S
        type(ten_3D2Osym)                    :: tau
        tau = self%F .transform. S
    end function PK2_to_Kirchhoff

    pure function Kirchhoff_to_PK2(self, tau) result(S)
        class(Base_F_kinematics), intent(in) :: self
        type(ten_3D2Osym), intent(in)        :: tau
        type(ten_3D2Osym)                    :: S
        S = self%F_inv .transform. tau
    end function Kirchhoff_to_PK2

    pure function Cauchy_to_Kirchhoff(self, sigma) result(tau)
        class(Base_F_kinematics), intent(in) :: self
        type(ten_3D2Osym), intent(in)        :: sigma
        type(ten_3D2Osym)                    :: tau
        tau = sigma * self%J
    end function Cauchy_to_Kirchhoff

    pure function Kirchhoff_to_Cauchy(self, tau) result(sigma)
        class(Base_F_kinematics), intent(in) :: self
        type(ten_3D2Osym), intent(in)        :: tau
        type(ten_3D2Osym)                    :: sigma
        sigma = tau / self%J
    end function Kirchhoff_to_Cauchy

    pure function PK2_to_PK1(self, S) result(P)
        class(Base_F_kinematics), intent(in) :: self
        type(ten_3D2Osym), intent(in)        :: S
        type(ten_3D2O)                       :: P
        P = self%F * S
    end function PK2_to_PK1

    pure function PK1_to_PK2(self, P) result(S)
        class(Base_F_kinematics), intent(in) :: self
        type(ten_3D2O), intent(in)           :: P
        type(ten_3D2Osym)                    :: S
        type(ten_3D2O)                       :: S_gen
        S_gen = self%F_inv * P
        call S%init(xx=S_gen%vals(1), yy=S_gen%vals(5), zz=S_gen%vals(9), &
                    xy=S_gen%vals(4), yz=S_gen%vals(8), xz=S_gen%vals(7))
    end function PK1_to_PK2

    pure function Cauchy_to_PK1(self, sigma) result(P)
        class(Base_F_kinematics), intent(in) :: self
        type(ten_3D2Osym), intent(in)        :: sigma
        type(ten_3D2O)                       :: P
        P = (sigma * self%F_inv%transpose()) * self%J
    end function Cauchy_to_PK1

    pure function PK1_to_Cauchy(self, P) result(sigma)
        class(Base_F_kinematics), intent(in) :: self
        type(ten_3D2O), intent(in)           :: P
        type(ten_3D2Osym)                    :: sigma
        type(ten_3D2O)                       :: sig_gen
        sig_gen = (P * self%F%transpose()) / self%J
        call sigma%init(xx=sig_gen%vals(1), yy=sig_gen%vals(5), zz=sig_gen%vals(9), &
                         xy=sig_gen%vals(4), yz=sig_gen%vals(8), xz=sig_gen%vals(7))
    end function PK1_to_Cauchy

    pure function finite_push_tangent(self, C_material, stress_spatial) result(c_spatial)
        !! True 4th-order Push-Forward: c_ijkl = (1/J) * F_iI * F_jJ * F_kK * F_lL * C_IJKL
        class(Base_F_kinematics), intent(in)     :: self
        type(ten_3D4O2sym), intent(in)           :: C_material
        type(ten_3D2Osym), intent(in), optional  :: stress_spatial
        type(ten_3D4O2sym)                       :: c_spatial

        c_spatial = (self%F .transform. C_material) / self%J
    end function finite_push_tangent

    pure function finite_pull_back_tangent(self, c_spatial, stress_material) result(C_material)
        !! True 4th-order Pull-Back: C_IJKL = J * F^-1_Ii * F^-1_Jj * F^-1_Kk * F^-1_Ll * c_ijkl
        class(Base_F_kinematics), intent(in)     :: self
        type(ten_3D4O2sym), intent(in)           :: c_spatial
        type(ten_3D2Osym), intent(in), optional  :: stress_material
        type(ten_3D4O2sym)                       :: C_material

        C_material = (self%F_inv .transform. c_spatial) * self%J
    end function finite_pull_back_tangent

end module muscle_kin_finite_base