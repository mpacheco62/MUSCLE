module mod_ops_dyadic
    !! This module defines the overloaded dyadic (outer) product operators `.tdot.`
    !! and `.tdotsym.` for mixed second-order symmetric and identity tensors.
    !!
    !! The standard dyadic product `a .tdot. b` computes the tensor product:
    !! \[ \mathbb{C} = \mathbf{a} \otimes \mathbf{b} \quad \implies \quad C_{ijkl} = a_{ij} b_{kl} \]
    !! The result is a fourth-order tensor with minor symmetries (`ten_3D4O2sym`).
    !!
    !! The symmetrized dyadic product `a .tdotsym. b` computes the symmetric part:
    !! \[ \mathbb{C} = \frac{1}{2} (\mathbf{a} \otimes \mathbf{b} + \mathbf{b} \otimes \mathbf{a}) \]
    !! This formulation preserves major symmetry, yielding a fully symmetric fourth-order tensor (`ten_3D4O3sym`).
    !! It is also overloaded for the unary-like self-product `.tdotsym. a` which computes \(\mathbf{a} \otimes \mathbf{a}\).

    use, intrinsic :: iso_fortran_env, only : real64
    
    ! Import necessary types
    use mod_iden_2O
    use mod_ten_2D2Osym
    use mod_ten_2D4O3sym
    use mod_ten_3D2Osym
    use mod_ten_3D4O2sym
    use mod_ten_3D4O3sym

    implicit none
    private

    ! =========================================================================
    ! PUBLIC INTERFACES
    ! =========================================================================
    
    public :: operator(.tdot.)
    interface operator (.tdot.)
        module procedure tdot_2D2Osym_2D2Osym
        module procedure tdot_3D2Osym_3D2Osym
        module procedure tdot_3D2Osym_I2O
        module procedure tdot_I2O_3D2Osym
    end interface

    public :: operator(.tdotsym.)
    interface operator (.tdotsym.)
        module procedure tdotsym_3D2Osym_3D2Osym
        module procedure tdotsym_3D2Osym
        module procedure tdotsym_3D2Osym_I2O
        module procedure tdotsym_I2O_3D2Osym
    end interface

contains

    ! =========================================================================
    ! 1. STANDARD DYADIC PRODUCT (.tdot.)
    ! =========================================================================

    pure function tdot_2D2Osym_2D2Osym(a, b) result(res)
        !! Computes the dyadic (outer) tensor product between two 2D symmetric 
        !! second-order tensors.
        !!
        !! Mathematically: \( \mathbb{C}_{ijkl} = a_{ij} b_{kl} \).
        !!
        !! The result is returned as a fully symmetric fourth-order tensor (`ten_2D4O3sym`).
        !! Note: This inherently enforces the major symmetry structure of the resulting 
        !! Voigt matrix, so this operator is most accurately used when \(\mathbf{a}\) 
        !! and \(\mathbf{b}\) are physically equivalent or proportional (e.g., \(\mathbf{I} \otimes \mathbf{I}\)).
        !!
        !! Voigt Mapping Reference:
        !! ```
        !!  | ( 1:1111) ( 5:1122) ( 8:1133) (10:1112) |
        !!  | ( 5:2211) ( 2:2222) ( 6:2233) ( 9:2212) |
        !!  | ( 8:3311) ( 6:3322) ( 3:3333) ( 7:3312) |
        !!  | (10:1211) ( 9:1222) ( 7:1233) ( 4:1212) |
        !! ```
        implicit none
        type(ten_2D2Osym), intent(in) :: a
            !! First second-order symmetric tensor \(\mathbf{a}\)
        type(ten_2D2Osym), intent(in) :: b
            !! Second second-order symmetric tensor \(\mathbf{b}\)
        type(ten_2D4O3sym) :: res
            !! Resulting fourth-order fully symmetric tensor \(\mathbb{C}\)
            
        res%vals( 1) = a%vals(1)*b%vals(1)
        res%vals( 5) = a%vals(1)*b%vals(2)
        res%vals( 8) = a%vals(1)*b%vals(3)
        res%vals(10) = a%vals(1)*b%vals(4)

        res%vals( 2) = a%vals(2)*b%vals(2)
        res%vals( 6) = a%vals(2)*b%vals(3)
        res%vals( 9) = a%vals(2)*b%vals(4)

        res%vals( 3) = a%vals(3)*b%vals(3)
        res%vals( 7) = a%vals(3)*b%vals(4)

        res%vals( 4) = a%vals(4)*b%vals(4)
    end function tdot_2D2Osym_2D2Osym

    pure function tdot_3D2Osym_3D2Osym(a, b) result(res)
        !! Computes the dyadic (outer) tensor product: \(\mathbb{C} = \mathbf{a} \otimes \mathbf{b}\).
        !!
        !! Mathematically: \( C_{ijkl} = a_{ij} b_{kl} \).
        !!
        !! The result is a fourth-order tensor that only has minor symmetries, even if
        !! \(\mathbf{a}\) and \(\mathbf{b}\) are symmetric. Major symmetry (\(C_{ijkl} = C_{klij}\))
        !! is only guaranteed if \(\mathbf{a}\) and \(\mathbf{b}\) are proportional.
        !! The result is therefore correctly returned as a `ten_3D4O2sym` type.
        !!
        !! In Voigt notation, this operation constructs a 6x6 matrix by taking the
        !! outer product of the 6-component Voigt vectors of \(\mathbf{a}\) and \(\mathbf{b}\).
        !!
        implicit none
        type(ten_3D2Osym), intent(in) :: a  
            !! The first second-order symmetric tensor \(\mathbf{a}\).
        type(ten_3D2Osym), intent(in) :: b  
            !! The second second-order symmetric tensor \(\mathbf{b}\).
        type(ten_3D4O2sym) :: res         
            !! The resulting fourth-order tensor with minor symmetries \(\mathbb{C}\).

        ! Build the 6x6 Voigt matrix C_IJ = a_I * b_J by constructing each row.
        res%vals(1,:) = a%vals(1) * b%vals(:)
        res%vals(2,:) = a%vals(2) * b%vals(:)
        res%vals(3,:) = a%vals(3) * b%vals(:)
        res%vals(4,:) = a%vals(4) * b%vals(:)
        res%vals(5,:) = a%vals(5) * b%vals(:)
        res%vals(6,:) = a%vals(6) * b%vals(:)

    end function tdot_3D2Osym_3D2Osym

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

    ! =========================================================================
    ! 2. SYMMETRIZED DYADIC PRODUCT (.tdotsym.)
    ! =========================================================================

    pure function tdotsym_3D2Osym_3D2Osym(a, b) result(res)
        !! Computes the symmetrized dyadic product: \(\mathbb{C} = \frac{1}{2} (\mathbf{a} \otimes \mathbf{b} + \mathbf{b} \otimes \mathbf{a})\).
        !! The result is always a fully symmetric 4th-order tensor (`ten_3D4O3sym`).
        !! When called as `a .tdotsym. a`, a modern compiler will optimize this to `a ⊗ a`.
        implicit none
        type(ten_3D2Osym), intent(in) :: a
        type(ten_3D2Osym), intent(in) :: b
        type(ten_3D4O3sym) :: res
        
        res%vals(1)  = 0.5D0 * (a%vals(1)*b%vals(1) + b%vals(1)*a%vals(1))
        res%vals(2)  = 0.5D0 * (a%vals(2)*b%vals(2) + b%vals(2)*a%vals(2))
        res%vals(3)  = 0.5D0 * (a%vals(3)*b%vals(3) + b%vals(3)*a%vals(3))
        res%vals(4)  = 0.5D0 * (a%vals(4)*b%vals(4) + b%vals(4)*a%vals(4))
        res%vals(5)  = 0.5D0 * (a%vals(5)*b%vals(5) + b%vals(5)*a%vals(5))
        res%vals(6)  = 0.5D0 * (a%vals(6)*b%vals(6) + b%vals(6)*a%vals(6))
        res%vals(7)  = 0.5D0 * (a%vals(1)*b%vals(2) + b%vals(1)*a%vals(2))
        res%vals(8)  = 0.5D0 * (a%vals(2)*b%vals(3) + b%vals(2)*a%vals(3))
        res%vals(9)  = 0.5D0 * (a%vals(3)*b%vals(4) + b%vals(3)*a%vals(4))
        res%vals(10) = 0.5D0 * (a%vals(4)*b%vals(5) + b%vals(4)*a%vals(5))
        res%vals(11) = 0.5D0 * (a%vals(5)*b%vals(6) + b%vals(5)*a%vals(6))
        res%vals(12) = 0.5D0 * (a%vals(1)*b%vals(3) + b%vals(1)*a%vals(3))
        res%vals(13) = 0.5D0 * (a%vals(2)*b%vals(4) + b%vals(2)*a%vals(4))
        res%vals(14) = 0.5D0 * (a%vals(3)*b%vals(5) + b%vals(3)*a%vals(5))
        res%vals(15) = 0.5D0 * (a%vals(4)*b%vals(6) + b%vals(4)*a%vals(6))
        res%vals(16) = 0.5D0 * (a%vals(1)*b%vals(4) + b%vals(1)*a%vals(4))
        res%vals(17) = 0.5D0 * (a%vals(2)*b%vals(5) + b%vals(2)*a%vals(5))
        res%vals(18) = 0.5D0 * (a%vals(3)*b%vals(6) + b%vals(3)*a%vals(6))
        res%vals(19) = 0.5D0 * (a%vals(1)*b%vals(5) + b%vals(1)*a%vals(5))
        res%vals(20) = 0.5D0 * (a%vals(2)*b%vals(6) + b%vals(2)*a%vals(6))
        res%vals(21) = 0.5D0 * (a%vals(1)*b%vals(6) + b%vals(1)*a%vals(6))

    end function tdotsym_3D2Osym_3D2Osym

    pure function tdotsym_3D2Osym(a) result(res)
        !! Computes the direct dyadic product of a tensor with itself: \(\mathbb{C} = \mathbf{a} \otimes \mathbf{a}\).
        !! This function provides a syntactically "unary" way to perform the dyadic product.
        !! Mathematically: \( C_{ijkl} = a_{ij} a_{kl} \).
        implicit none
        type(ten_3D2Osym), intent(in) :: a
            !! The symmetric 2nd-order tensor to be multiplied by itself.
        type(ten_3D4O3sym) :: res
            !! The resulting fully symmetric 4th-order tensor.
        
        ! Compute C_IJ = a_I * a_J for the 21 unique components
        res%vals(1)  = a%vals(1)*a%vals(1)
        res%vals(2)  = a%vals(2)*a%vals(2)
        res%vals(3)  = a%vals(3)*a%vals(3)
        res%vals(4)  = a%vals(4)*a%vals(4)
        res%vals(5)  = a%vals(5)*a%vals(5)
        res%vals(6)  = a%vals(6)*a%vals(6)
        res%vals(7)  = a%vals(1)*a%vals(2)
        res%vals(8)  = a%vals(2)*a%vals(3)
        res%vals(9)  = a%vals(3)*a%vals(4)
        res%vals(10) = a%vals(4)*a%vals(5)
        res%vals(11) = a%vals(5)*a%vals(6)
        res%vals(12) = a%vals(1)*a%vals(3)
        res%vals(13) = a%vals(2)*a%vals(4)
        res%vals(14) = a%vals(3)*a%vals(5)
        res%vals(15) = a%vals(4)*a%vals(6)
        res%vals(16) = a%vals(1)*a%vals(4)
        res%vals(17) = a%vals(2)*a%vals(5)
        res%vals(18) = a%vals(3)*a%vals(6)
        res%vals(19) = a%vals(1)*a%vals(5)
        res%vals(20) = a%vals(2)*a%vals(6)
        res%vals(21) = a%vals(1)*a%vals(6)

    end function tdotsym_3D2Osym

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

end module mod_ops_dyadic