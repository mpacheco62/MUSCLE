! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2025 Matias Pacheco-Alarcon <matias.pacheco.a@gmail.com>

module muscle_tensor_2d2osym
    !! Module muscle_tensor_2d2osym
    !! ============================
    !!
    !! Defines the type for symmetric 2D (plane / axisymmetric) second-order tensors
    !! and associated operations.
    !!
    !! This module provides the derived type `ten_2D2Osym` to represent a symmetric
    !! second-order tensor in a two-dimensional setting (plane strain, plane stress or
    !! axisymmetry), where the out-of-plane shear components vanish
    !! (\(A_{13} = A_{23} = 0\)) but the out-of-plane normal component \(A_{33}\) is kept.
    !! The tensor is stored internally using Voigt notation with 4 components
    !! in the order (11, 22, 33, 12).
    !!
    !! The module overloads standard arithmetic operators (+, -, *, /), a custom
    !! equality comparison operator (.approx.), the deviatoric operator (.dev.),
    !! the double dot product operator (.ddot.) and the inverse operator (.inv.)
    !! for this tensor type. It also provides methods for initialization and
    !! accessing individual components.
    !!
    !! Public Entities
    !! ---------------
    !!
    !! ### Derived Type:
    !!
    !! - `ten_2D2Osym`: Represents a symmetric 2D second-order tensor.
    !!     - Component: `vals(4) :: real(real64)` - Stores the 4 Voigt components
    !!       (xx, yy, zz, xy).
    !!     - Generic Procedure: `init` - Initializes the tensor either from a
    !!       4-element array or from individual xx, yy, zz, xy components.
    !!     - Procedures: `xx`, `yy`, `zz`, `xy` - Accessor functions
    !!       for individual tensor components.
    !!     - Procedures: `square`, `norm`, `is_approx`, `det`, `inv`.
    !!
    !! ### Operators:
    !!
    !! - `.approx.`: Compares two `ten_2D2Osym` tensors for approximate equality.
    !! - `+`: Adds two `ten_2D2Osym` tensors.
    !! - `-`: Subtracts two `ten_2D2Osym` tensors (binary) or computes the unary negation.
    !! - `*`: Multiplies a `ten_2D2Osym` tensor by a `real(real64)` scalar (or vice-versa).
    !! - `/`: Divides a `ten_2D2Osym` tensor by a `real(real64)` scalar.
    !! - `.dev.`: Computes the deviatoric part of a `ten_2D2Osym` tensor.
    !! - `.ddot.`: Computes the double dot product (scalar result) of two `ten_2D2Osym` tensors.
    !! - `.inv.`: Computes the analytical inverse of a `ten_2D2Osym` tensor.
    !!
    !! ### Assignment:
    !!
    !! - `=`: Allows assigning a single `real(real64)` scalar value to all components
    !!        of a `ten_2D2Osym` tensor.
    !!
    !! Usage
    !! -----
    !!
    !! ```fortran
    !! program example_ten_2D2Osym_usage
    !!   use muscle_tensor_2d2osym
    !!   use iso_fortran_env, only: real64
    !!   implicit none
    !!
    !!   type(ten_2D2Osym) :: stress, strain, stress_dev
    !!   real(real64) :: dot_product
    !!   logical :: are_equal
    !!
    !!   ! Initialize using individual components (Voigt: 11, 22, 33, 12)
    !!   call stress%init(100.0D0, 50.0D0, 20.0D0, 10.0D0)
    !!
    !!   ! Initialize strain to zero using scalar assignment
    !!   strain = 0.0D0
    !!   ! Set some strain components
    !!   strain%vals(1) = 0.001 ! e_xx
    !!   strain%vals(4) = 0.002 ! e_xy
    !!
    !!   ! Calculate deviatoric stress
    !!   stress_dev = .dev. stress
    !!
    !!   ! Calculate double dot product
    !!   dot_product = stress .ddot. strain
    !!
    !!   ! Access a component
    !!   print *, "Stress XX component:", stress%xx() ! Note the parentheses for the function call
    !!   print *, "Stress XY component:", stress%xy()
    !!   print *, "Deviatoric Stress ZZ:", stress_dev%zz()
    !!   print *, "Stress : Strain =", dot_product
    !!
    !!   ! Comparison
    !!   are_equal = (stress .approx. stress_dev)
    !!   print *, "Is stress equal to its deviatoric part?", are_equal
    !!
    !! end program example_ten_2D2Osym_usage
    !! ```
    !!
    !! For more information see [[muscle_tensors]]

    use, intrinsic :: iso_fortran_env
    implicit none
    private

    type, public :: ten_2D2Osym
        !! Symmetric 2D Second-Order Tensor (Voigt Notation)
        !! ==================================================
        !!
        !! Represents a symmetric second-order tensor in a 2D setting (plane or
        !! axisymmetric), such as stress (\(\sigma_{ij}\)) or strain (\(\epsilon_{ij}\)),
        !! where \(\sigma_{ij} = \sigma_{ji}\) and \(\sigma_{13} = \sigma_{23} = 0\).
        !! Only 4 independent components are needed.
        !!
        !! Storage:
        !! --------
        !! The tensor components are stored internally in a 1D array `vals` of size 4
        !! using Voigt notation with the following mapping:
        !! - `vals(1)`: Component (1,1) or xx
        !! - `vals(2)`: Component (2,2) or yy
        !! - `vals(3)`: Component (3,3) or zz
        !! - `vals(4)`: Component (1,2) or xy (Note: This is the tensorial shear component, not the engineering one)
        !!
        !! Access:
        !! -------
        !! Components can be accessed directly via the `vals` array or more conveniently
        !! using the type-bound procedures `xx()`, `yy()`, `zz()`, `xy()`.
        !!
        !! Initialization:
        !! ---------------
        !! Use the generic `init` procedure to initialize either from a 4-element `real(real64)`
        !! array (following the Voigt order above) or by providing the 4 components
        !! individually (xx, yy, zz, xy).
        !!
        !! For more information see [[muscle_tensors]]
        real(real64), dimension(4) :: vals
            !! Stores the 4 independent components in Voigt notation: (xx, yy, zz, xy).
        contains
            generic, public :: init => init_ten_2D2Osym, init2_ten_2D2Osym
                !! Generic interface for initialization.
            procedure, private :: init_ten_2D2Osym, init2_ten_2D2Osym
            procedure, public :: xx
                !! Accessor for the xx (1,1) component.
            procedure, public :: yy
                !! Accessor for the yy (2,2) component.
            procedure, public :: zz
                !! Accessor for the zz (3,3) component.
            procedure, public :: xy
                !! Accessor for the xy (1,2) component.
            procedure, public :: square => square_2D2Osym
                !! Computes the square of the tensor.
            procedure, public :: norm => norm_2D2Osym
                !! Computes the norm of the tensor.
            procedure, public :: is_approx => is_approx_2D2Osym
                !! Compares two tensors for approximate equality.
            procedure, public :: det => det_2D2Osym
                !! Computes the determinant of the symmetric 2nd-order tensor.
            procedure, public :: inv => inv_2D2Osym
                !! Computes the analytical inverse tensor (A^-1).
    end type ten_2D2Osym

    public :: operator(.approx.)
    interface operator (.approx.)
        module procedure approx_2D2Osym
    end interface

    public :: operator(+)
    interface operator (+)
        module procedure sum_2D2Osym
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_2D2Osym
        module procedure subU_2D2Osym
    end interface

    public :: operator(*)
    interface operator (*)
        module procedure mul_real64_2D2Osym
        module procedure mul_2D2Osym_real64
    end interface

    public :: operator( / )
    interface operator ( / )
        module procedure div_2D2Osym_real64
    end interface

    public :: operator(.dev.)
    interface operator (.dev.)
        module procedure dev_2D2Osym
    end interface

    public :: operator(.ddot.)
    interface operator (.ddot.)
        module procedure ddot_2D2Osym_2D2Osym
    end interface

    public :: assignment (=)
    interface assignment (=)
        module procedure ten_2D2Osym_real64_assign
    end interface

    public :: operator(.inv.)
    interface operator (.inv.)
        module procedure inv_2D2Osym
    end interface

#ifdef ENABLE_UDTIO
    public :: write(formatted)
    interface write(formatted)
        module procedure print_ten_2D2Osym
    end interface
#endif

contains

    pure subroutine ten_2D2Osym_real64_assign(a, b)
        !! Assigns the scalar `b` to every stored component of `a`.
        implicit none
        type(ten_2D2Osym), intent(out) :: a
            !! Target tensor.
        real(real64), intent(in) :: b
            !! Scalar value copied to all 4 components.
        a%vals = b
    end subroutine

    pure subroutine init_ten_2D2Osym(self, vals)
        !! Initializes a ten_2D2Osym tensor from a 4-element array (Voigt order).
        !!
        !! voigt notation used: 11, 22, 33, 12
        implicit none
        class(ten_2D2Osym), intent(inout) :: self
            !! Tensor to initialize.
        real(real64), intent(in) :: vals(4)
            !! Components in Voigt order (xx, yy, zz, xy).
        self%vals = vals
    end subroutine

    pure subroutine init2_ten_2D2Osym(self, xx, yy, zz, xy)
        !! Initializes a ten_2D2Osym tensor from its 4 individual components.
        implicit none
        class(ten_2D2Osym), intent(inout) :: self
            !! Tensor to initialize.
        real(real64), intent(in) :: xx
            !! Component (1,1).
        real(real64), intent(in) :: yy
            !! Component (2,2).
        real(real64), intent(in) :: zz
            !! Component (3,3).
        real(real64), intent(in) :: xy
            !! Component (1,2) = (2,1) (tensorial shear).
        self%vals = (/xx, yy, zz, xy/)
    end subroutine

    pure function norm_2D2Osym(a) result(res)
        !! Computes the norm of a ten_2D2Osym tensor.
        !! Uses a modified L1 norm (shear component weighted by 2):
        !! norm(a) = |a_11| + |a_22| + |a_33| + 2|a_12|
        implicit none
        class(ten_2D2Osym), intent(in) :: a
            !! Tensor whose norm is computed.
        real(real64) :: res
            !! Modified L1 norm.

        res =   abs(a%vals(1)) + abs(a%vals(2)) + abs(a%vals(3)) &
              + 2*abs(a%vals(4))
    end function norm_2D2Osym

    pure function is_approx_2D2Osym(a, b, tol) result(res)
        !! Component-wise approximate equality with a relative tolerance:
        !! \( |a_I - b_I| \le tol \cdot \max(\max|a|, \max|b|, 10^{-30}) \) for all I.
        implicit none
        class(ten_2D2Osym), intent(in) :: a
            !! First tensor.
        class(ten_2D2Osym), intent(in) :: b
            !! Second tensor.
        real(real64), optional, intent(in) :: tol
            !! Relative tolerance (default `1E-12`).
        logical :: res
            !! `.TRUE.` if both tensors are approximately equal.
        real(real64), parameter :: EPS_ABS=1e-30
        real(real64) :: eps_check
        real(real64) :: tol2
        real(real64) :: max_val

        tol2 = 1E-12
        if (present(tol)) tol2 = tol
        max_val = max(maxval(abs(a%vals)), maxval(abs(b%vals)), EPS_ABS)
        eps_check = tol2 * max_val

        res = all(abs(a%vals - b%vals) .le. eps_check)
    end function is_approx_2D2Osym

    pure function approx_2D2Osym(a, b) result(res)
        !! `.approx.` Compares two ten_2D2Osym tensors for approximate equality.
        !! Delegates to `is_approx` with its default tolerance.
        implicit none
        type(ten_2D2Osym), intent(in) :: a
            !! First tensor.
        type(ten_2D2Osym), intent(in) :: b
            !! Second tensor.
        logical :: res
            !! `.TRUE.` if both tensors are approximately equal.

        res = a%is_approx(b)
    end function approx_2D2Osym

    pure function sum_2D2Osym(a, b) result(res)
        !! Component-wise sum \(\mathbf{a} + \mathbf{b}\).
        implicit none
        type(ten_2D2Osym), intent(in) :: a
            !! First operand.
        type(ten_2D2Osym), intent(in) :: b
            !! Second operand.
        type(ten_2D2Osym) :: res
            !! Sum.
        res%vals = a%vals + b%vals
    end function sum_2D2Osym

    pure function sub_2D2Osym(a, b) result(res)
        !! Component-wise difference \(\mathbf{a} - \mathbf{b}\).
        implicit none
        type(ten_2D2Osym), intent(in) :: a
            !! Minuend.
        type(ten_2D2Osym), intent(in) :: b
            !! Subtrahend.
        type(ten_2D2Osym) :: res
            !! Difference.
        res%vals = a%vals - b%vals
    end function sub_2D2Osym

    pure function subU_2D2Osym(a) result(res)
        !! Unary negation \(-\mathbf{a}\).
        implicit none
        type(ten_2D2Osym), intent(in) :: a
            !! Operand.
        type(ten_2D2Osym) :: res
            !! Negated tensor.
        res%vals = -a%vals
    end function subU_2D2Osym

    pure function mul_real64_2D2Osym(a, b) result(res)
        !! Scalar product \(a\,\mathbf{b}\).
        implicit none
        real(real64), intent(in) :: a
            !! Scalar factor.
        type(ten_2D2Osym), intent(in) :: b
            !! Tensor.
        type(ten_2D2Osym) :: res
            !! Scaled tensor.
        res%vals = a * b%vals
    end function mul_real64_2D2Osym

    pure function mul_2D2Osym_real64(a, b) result(res)
        !! Scalar product \(\mathbf{a}\,b\).
        implicit none
        type(ten_2D2Osym), intent(in) :: a
            !! Tensor.
        real(real64), intent(in) :: b
            !! Scalar factor.
        type(ten_2D2Osym) :: res
            !! Scaled tensor.
        res%vals =  a%vals * b
    end function mul_2D2Osym_real64

    pure function div_2D2Osym_real64(a, b) result(res)
        !! Scalar division \(\mathbf{a}/b\).
        implicit none
        type(ten_2D2Osym), intent(in) :: a
            !! Tensor.
        real(real64), intent(in) :: b
            !! Scalar divisor (must be non-zero).
        type(ten_2D2Osym) :: res
            !! Scaled tensor.
        res%vals = a%vals/b
    end function div_2D2Osym_real64

    pure function ddot_2D2Osym_2D2Osym(a, b) result(res)
        !! Double contraction \(\mathbf{a} : \mathbf{b} = a_{ij} b_{ij}\)
        !! (the shear product is counted twice).
        implicit none
        type(ten_2D2Osym), intent(in) :: a
            !! First operand.
        type(ten_2D2Osym), intent(in) :: b
            !! Second operand.
        real(real64) :: res
            !! Scalar result.
        res =   a%vals(1)*b%vals(1) &
              + a%vals(2)*b%vals(2) &
              + a%vals(3)*b%vals(3) &
              + 2*a%vals(4)*b%vals(4)
    end function ddot_2D2Osym_2D2Osym

    pure function dev_2D2Osym(a) result(res)
        !! Deviatoric part \(\mathbf{a} - \frac{1}{3}\mathrm{tr}(\mathbf{a})\mathbf{I}\).
        implicit none
        type(ten_2D2Osym), intent(in) :: a
            !! Operand.
        type(ten_2D2Osym) :: res
            !! Deviatoric tensor.
        real(real64) :: hydro
        hydro = (a%vals(1) + a%vals(2) + a%vals(3))/3D0
        res%vals(1:3) = a%vals(1:3) - hydro
        res%vals(4) = a%vals(4)
    end function dev_2D2Osym

    pure function square_2D2Osym(a) result(res)
        !! Square \(\mathbf{a}\cdot\mathbf{a}\) (symmetric for a symmetric tensor).
        implicit none
        class(ten_2D2Osym), intent(in) :: a
            !! Operand.
        type(ten_2D2Osym) :: res
            !! Squared tensor.
        res%vals(1) = a%vals(1)**2 + a%vals(4)**2
        res%vals(2) = a%vals(2)**2 + a%vals(4)**2
        res%vals(3) = a%vals(3)**2
        res%vals(4) = a%vals(4)*(a%vals(1) + a%vals(2))
    end function square_2D2Osym

    pure function xx(a) result(res)
        !! Accessor function for the xx (11) component (vals(1)).
        implicit none
        class(ten_2D2Osym), intent(in) :: a
            !! Tensor.
        real(real64) :: res
            !! Component (1,1).
        res = a%vals(1)
    end function xx

    pure function yy(a) result(res)
        !! Accessor function for the yy (22) component (vals(2)).
        implicit none
        class(ten_2D2Osym), intent(in) :: a
            !! Tensor.
        real(real64) :: res
            !! Component (2,2).
        res = a%vals(2)
    end function yy

    pure function zz(a) result(res)
        !! Accessor function for the zz (33) component (vals(3)).
        implicit none
        class(ten_2D2Osym), intent(in) :: a
            !! Tensor.
        real(real64) :: res
            !! Component (3,3).
        res = a%vals(3)
    end function zz

    pure function xy(a) result(res)
        !! Accessor function for the xy (12) component (vals(4)).
        implicit none
        class(ten_2D2Osym), intent(in) :: a
            !! Tensor.
        real(real64) :: res
            !! Component (1,2).
        res = a%vals(4)
    end function xy

#ifdef ENABLE_UDTIO
    subroutine print_ten_2D2Osym(dtv, unit, iotype, v_list, iostat, iomsg)
        !! Custom I/O formatting for ten_2D2Osym.
        !! List mode prints the 4 Voigt components; matrix mode prints the full 3x3
        !! matrix with the vanishing out-of-plane shear terms.
        class(ten_2D2Osym), intent(in) :: dtv
            !! Tensor to print.
        integer, intent(in)            :: unit
            !! Output unit.
        character(len=*), intent(in)   :: iotype
            !! Edit descriptor type (`LISTDIRECTED`, `DT...`).
        integer, intent(in)            :: v_list(:)
            !! Optional (width, decimals).
        integer, intent(out)           :: iostat
            !! I/O status.
        character(len=*), intent(inout) :: iomsg
            !! I/O error message.

        character(len=100) :: fmt_string
        integer :: w, d

        w = 11
        d = 4
        iostat = 0

        if (size(v_list) >= 1) w = v_list(1)
        if (size(v_list) >= 2) d = v_list(2)

        if (iotype == "DTLIST" .or. iotype == "LIST") then
            write(fmt_string, "('(A, 4(ES', I0, '.', I0, ', 1X), A)')") w, d
            write(unit, fmt_string, iostat=iostat) &
                "[", dtv%xx(), dtv%yy(), dtv%zz(), dtv%xy(), "]"
        else
            write(fmt_string, "('(/, 3(ES', I0, '.', I0, ', 2X), /, 3(ES', I0, '.', I0, ', 2X), /, 3(ES', I0, '.', I0, ', 2X))')") &
                w, d, w, d, w, d

            write(unit, fmt_string, iostat=iostat) &
                dtv%xx(), dtv%xy(), 0D0,      & ! Row 1
                dtv%xy(), dtv%yy(), 0D0,      & ! Row 2
                0D0,      0D0,      dtv%zz()    ! Row 3
        end if

        if (iostat /= 0) then
            iomsg = "Error in print_ten_2D2Osym: Failed to write to the specified unit."
        end if
    end subroutine print_ten_2D2Osym
#endif

    pure function det_2D2Osym(self) result(res)
        !! Determinant of a 2D symmetric 2nd-order tensor:
        !! \( \det\mathbf{A} = (A_{11}A_{22} - A_{12}^2)\,A_{33} \).
        implicit none
        class(ten_2D2Osym), intent(in) :: self
            !! Tensor.
        real(real64)                   :: res
            !! Determinant.

        res = (self%vals(1)*self%vals(2) - self%vals(4)**2) * self%vals(3)
    end function det_2D2Osym

    pure function inv_2D2Osym(self) result(res)
        !! Analytical inverse of a 2D symmetric 2nd-order tensor. The in-plane 2x2
        !! block and the out-of-plane component are inverted independently:
        !! \[ \mathbf{A}^{-1} = \frac{1}{A_{11}A_{22}-A_{12}^2}
        !!    \begin{bmatrix} A_{22} & -A_{12} \\ -A_{12} & A_{11} \end{bmatrix}
        !!    \oplus \frac{1}{A_{33}} \]
        !! If \(\det\mathbf{A} \approx 0\) (\(|\det| < 10^{-30}\)), a zero tensor is returned.
        implicit none
        class(ten_2D2Osym), intent(in) :: self
            !! Tensor to invert.
        type(ten_2D2Osym)             :: res
            !! Inverse tensor (zero if singular).

        real(real64) :: det2, det_val, inv_det2
        real(real64), parameter :: EPS_DET = 1.0D-30

        det2 = self%vals(1)*self%vals(2) - self%vals(4)**2
        det_val = det2 * self%vals(3)

        if (abs(det_val) < EPS_DET) then
            res%vals = 0.0D0
            return
        end if

        inv_det2 = 1.0D0 / det2
        res%vals(1) =  self%vals(2) * inv_det2 ! xx
        res%vals(2) =  self%vals(1) * inv_det2 ! yy
        res%vals(3) =  1.0D0 / self%vals(3)    ! zz
        res%vals(4) = -self%vals(4) * inv_det2 ! xy
    end function inv_2D2Osym

end module muscle_tensor_2d2osym
