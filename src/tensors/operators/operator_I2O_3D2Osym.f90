module mod_operator_I2O_3D2Osym
    !! Module mod_operator_I2O_3D2Osym
    !! ====================================
    !!
    !! Defines mixed algebraic operations between the standard 3D identity tensor
    !! (`iden_2O`) and 3D symmetric second-order tensors (`ten_3D2Osym`).
    !!
    !! This module handles the 6-component Voigt storage format where diagonal 
    !! elements occupy the first three positions. It also includes the dyadic 
    !! product with the identity, resulting in 4th-order tensors.
    !!
    !! Overloaded Operators
    !! --------------------
    !!
    !! - `+` : Addition of a tensor and the identity.
    !! - `-` : Subtraction of a tensor and the identity.
    !! - `*` : Single contraction (matrix product) with the identity.
    !! - `.ddot.` : Double contraction, computing the trace of the tensor.
    !! - `.tdot.` : Dyadic product (\(\mathbf{a} \otimes \mathbf{I}\) or \(\mathbf{I} \otimes \mathbf{a}\)).
    !!
    !! For tensor type definitions, see [[tensors_types]].

    use, intrinsic :: iso_fortran_env
    use mod_iden_2O
    use mod_ten_3D2Osym
    use mod_ten_3D4O3sym
    implicit none
    private
    
    public :: operator(*)
    interface operator (*)
        module procedure mul_3D2Osym_I2O
        module procedure mul_I2O_3D2Osym
    end interface

    public :: operator(.ddot.)
    interface operator (.ddot.)
        module procedure ddot_I2O_3D2Osym
        module procedure ddot_3D2Osym_I2O
    end interface

    public :: operator(.tdot.)
    interface operator (.tdot.)
        module procedure tdot_3D2Osym_I2O
        module procedure tdot_I2O_3D2Osym
    end interface

    public :: operator(.tdotsym.)
    interface operator (.tdotsym.)
        module procedure tdotsym_3D2Osym_I2O
        module procedure tdotsym_I2O_3D2Osym
    end interface

    public :: assignment (=)
    interface assignment (=)
        module procedure assign_3D2Osym_I2O
    end interface

contains

    pure function ddot_I2O_3D2Osym(I2, b) result(res)
        !! Computes the double contraction \(\text{tr}(\mathbf{B}) = \mathbf{I} : \mathbf{B}\).
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: b
        real(real64) :: res
        res = b%vals(1) + b%vals(2) + b%vals(3)
    end function ddot_I2O_3D2Osym

    pure function ddot_3D2Osym_I2O(b, I2) result(res)
        !! Computes the double contraction \(\text{tr}(\mathbf{B}) = \mathbf{B} : \mathbf{I}\).
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: b
        real(real64) :: res
        res = b%vals(1) + b%vals(2) + b%vals(3)
    end function ddot_3D2Osym_I2O

    pure function mul_I2O_3D2Osym(I2, a) result(res)
        !! Computes the single contraction \(\mathbf{res} = \mathbf{I} \cdot \mathbf{A} = \mathbf{A}\).
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals = a%vals
    end function mul_I2O_3D2Osym

    pure function mul_3D2Osym_I2O(a, I2) result(res)
        !! Computes the single contraction \(\mathbf{res} = \mathbf{A} \cdot \mathbf{I} = \mathbf{A}\).
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym) :: res
        res%vals = a%vals
    end function mul_3D2Osym_I2O

    pure function tdot_3D2Osym_I2O(a, I2) result(res)
        !! Computes the dyadic product \(\mathbb{C} = \mathbf{a} \otimes \mathbf{I}\).
        !! In 6x6 Voigt notation, this fills the first three columns with the vector \(\mathbf{a}\).
        use mod_ten_3D4O2sym
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: a
        type(ten_3D4O2sym) :: res
        res%vals = 0.0D0
        res%vals(:,1) = a%vals
        res%vals(:,2) = a%vals
        res%vals(:,3) = a%vals
    end function tdot_3D2Osym_I2O

    pure function tdot_I2O_3D2Osym(I2, a) result(res)
        !! Computes the dyadic product \(\mathbb{C} = \mathbf{I} \otimes \mathbf{a}\).
        !! In 6x6 Voigt notation, this fills the first three rows with the vector \(\mathbf{a}\).
        use mod_ten_3D4O2sym
        implicit none
        type(iden_2O), intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: a
        type(ten_3D4O2sym) :: res
        res%vals = 0.0D0
        res%vals(1,:) = a%vals
        res%vals(2,:) = a%vals
        res%vals(3,:) = a%vals
    end function tdot_I2O_3D2Osym

    pure function tdotsym_3D2Osym_I2O(a, I2) result(res)
        !! Computes the symmetrized dyadic product \(\mathbb{C} = \frac{1}{2}(\mathbf{a} \otimes \mathbf{I} + \mathbf{I} \otimes \mathbf{a})\).
        !!
        !! Unlike the simple dyadic product, the symmetrized version guarantees major symmetry 
        !! (\(\mathbb{C}_{ijkl} = \mathbb{C}_{klij}\)). Therefore, the result is returned as a fully 
        !! symmetric fourth-order tensor (`ten_3D4O3sym` with 21 compressed components).
        !!
        !! In 6x6 Voigt notation, this operation symmetrically distributes the components of \(\mathbf{a}\) 
        !! across the rows and columns associated with the normal components of the identity.
        implicit none
        type(ten_3D2Osym), intent(in) :: a   !! Symmetric second-order tensor \(\mathbf{a}\)
        type(iden_2O),     intent(in) :: I2  !! Standard identity tensor \(\mathbf{I}\)
        type(ten_3D4O3sym)            :: res !! Fully symmetric fourth-order tensor
        
        real(real64) :: xx, yy, zz, xy, yz, xz
        
        ! Extract components for faster access
        xx = a%vals(1)
        yy = a%vals(2)
        zz = a%vals(3)
        xy = a%vals(4)
        yz = a%vals(5)
        xz = a%vals(6)

        ! first row
        res%vals(1) = xx
        res%vals(7) = 0.5*(xx + yy)
        res%vals(12) = 0.5*(xx + zz)
        res%vals(16) = 0.5D0*xy
        res%vals(19) = 0.5D0*yz
        res%vals(21) = 0.5D0*xz

        ! second row
        res%vals(2) = yy
        res%vals(8) = 0.5*(yy+zz)
        res%vals(13) = 0.5D0*xy
        res%vals(17) = 0.5D0*yz
        res%vals(20) = 0.5D0*xz

        ! third row
        res%vals(3) = zz
        res%vals(9) = 0.5D0*xy
        res%vals(14) = 0.5D0*yz
        res%vals(18) = 0.5D0*xz

        ! fourth row
        res%vals(4) = 0
        res%vals(10) = 0
        res%vals(15) = 0

        ! fifth row
        res%vals(5) = 0
        res%vals(11) = 0

        ! sixth row
        res%vals(6) = 0

    end function tdotsym_3D2Osym_I2O

    pure function tdotsym_I2O_3D2Osym(I2, a) result(res)
        !! Computes the symmetrized dyadic product \(\mathbb{C} = \frac{1}{2}(\mathbf{I} \otimes \mathbf{a} + \mathbf{a} \otimes \mathbf{I})\).
        !!
        !! Due to the commutativity of the addition, this is mathematically identical to `a .tdotsym. I2`.
        implicit none
        type(iden_2O),     intent(in) :: I2
        type(ten_3D2Osym), intent(in) :: a
        type(ten_3D4O3sym)            :: res
        
        ! Delegate to the primary implementation
        res = tdotsym_3D2Osym_I2O(a, I2)
    end function tdotsym_I2O_3D2Osym
    
    pure subroutine assign_3D2Osym_I2O(a, b)
        !! Explicit assignment from a second order identity to a 3D symmetric tensor.
        !! 
        !! This is necessary to allow statements like `a = I2` where `a` is a `ten_3D2Osym` and `I2` is an `iden_2O`.
        !! The resulting tensor `a` will have its diagonal components set to 1 and off-diagonal components set to 0.
        implicit none
        type(ten_3D2Osym), intent(out) :: a
            !! The target symmetric tensor to be overwritten.
        type(iden_2O), intent(in) :: b
            !! The source second-order identity tensor.

        a%vals(1) = 1.0D0 ! xx
        a%vals(2) = 1.0D0 ! yy
        a%vals(3) = 1.0D0 ! zz
        a%vals(4) = 0.0D0 ! xy
        a%vals(5) = 0.0D0 ! yz
        a%vals(6) = 0.0D0 ! xz
    end subroutine assign_3D2Osym_I2O
end module mod_operator_I2O_3D2Osym