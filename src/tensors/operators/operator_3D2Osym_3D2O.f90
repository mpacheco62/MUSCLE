module mod_operator_3D2Osym_3D2O
    !! Module mod_operator_3D2Osym_3D2O
    !! ================================
    !!
    !! Defines mixed algebraic operations involving 3D symmetric second-order 
    !! tensors (`ten_3D2Osym`) and general (non-symmetric) 3D second-order 
    !! tensors (`ten_3D2O`).
    !!
    !! Since operations between a symmetric and a non-symmetric tensor generally 
    !! break symmetry, the resulting tensors from these operations are mostly 
    !! cast into the general `ten_3D2O` format (9 components, column-major).
    !!
    !! Overloaded Operators
    !! --------------------
    !!
    !! - `+` : Addition.
    !!   - `3D2O + 3D2Osym` -> `3D2O`
    !!   - `3D2Osym + 3D2O` -> `3D2O`
    !!
    !! - `-` : Subtraction.
    !!   - `3D2O - 3D2Osym` -> `3D2O`
    !!   - `3D2Osym - 3D2O` -> `3D2O`
    !!
    !! - `*` : Single tensor contraction (Matrix dot product).
    !!   - \(\mathbf{C} = \mathbf{A} \cdot \mathbf{B}\)
    !!   - `3D2O * 3D2Osym` -> `3D2O`
    !!   - `3D2Osym * 3D2O` -> `3D2O`
    !!   - `3D2Osym * 3D2Osym` -> `3D2O` (Since the product of two symmetric tensors is generally non-symmetric).
    !!
    !! - `.ddot.` : Double tensor contraction (Frobenius inner product).
    !!   - \(\alpha = \mathbf{A} : \mathbf{B}\) -> `real(real64)`
    !!
    !! - `=` : Assignment.
    !!   - Implicit casting from a symmetric tensor to a general tensor format.

    use, intrinsic :: iso_fortran_env
    use mod_ten_3D2Osym
    use mod_ten_3D2O
    implicit none
    private
    
    public :: operator(*)
    interface operator (*)
        module procedure dot_3D2Osym_3D2Osym
        module procedure dot_3D2O_3D2Osym
        module procedure dot_3D2Osym_3D2O
    end interface

    public :: operator(+)
    interface operator (+)
        module procedure sum_3D2O_3D2Osym
        module procedure sum_3D2Osym_3D2O
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure sub_3D2O_3D2Osym
        module procedure sub_3D2Osym_3D2O
    end interface

    public :: operator(.ddot.)
    interface operator (.ddot.)
        module procedure ddot_3D2O_3D2Osym
        module procedure ddot_3D2Osym_3D2O
    end interface

    public :: assignment (=)
    interface assignment (=)
        module procedure assign_3D2O_3D2Osym
    end interface

contains

    ! =========================================================================
    ! SINGLE CONTRACTION (MATRIX PRODUCT)
    ! =========================================================================

    pure function dot_3D2Osym_3D2Osym(a, b) result(res)
        !! Computes the single tensor contraction (matrix product) between two 
        !! symmetric second-order tensors.
        !!
        !! Mathematically: \( res_{ij} = a_{ik} b_{kj} \)
        !!
        !! Note: The product of two symmetric tensors is only symmetric if they 
        !! commute. Thus, the result is safely returned as a general `ten_3D2O` tensor.
        !!
        !! Performance: Manual loop unrolling mapping the 6-component Voigt 
        !! arrays to a 9-component column-major array.
        implicit none
        type(ten_3D2Osym), intent(in) :: a
            !! First symmetric second-order tensor \(\mathbf{A}\) (6 components)
        type(ten_3D2Osym), intent(in) :: b
            !! Second symmetric second-order tensor \(\mathbf{B}\) (6 components)
        type(ten_3D2O) :: res
            !! Resulting general second-order tensor \(\mathbf{res}\) (9 components)
        
        ! Column 1 (xx, yx, zx)
        res%vals(1) = a%vals(1)*b%vals(1) + a%vals(4)*b%vals(4) + a%vals(6)*b%vals(6)
        res%vals(2) = a%vals(4)*b%vals(1) + a%vals(2)*b%vals(4) + a%vals(5)*b%vals(6)
        res%vals(3) = a%vals(6)*b%vals(1) + a%vals(5)*b%vals(4) + a%vals(3)*b%vals(6)
        
        ! Column 2 (xy, yy, zy)
        res%vals(4) = a%vals(1)*b%vals(4) + a%vals(4)*b%vals(2) + a%vals(6)*b%vals(5)
        res%vals(5) = a%vals(4)*b%vals(4) + a%vals(2)*b%vals(2) + a%vals(5)*b%vals(5)
        res%vals(6) = a%vals(6)*b%vals(4) + a%vals(5)*b%vals(2) + a%vals(3)*b%vals(5)
        
        ! Column 3 (xz, yz, zz)
        res%vals(7) = a%vals(1)*b%vals(6) + a%vals(4)*b%vals(5) + a%vals(6)*b%vals(3)
        res%vals(8) = a%vals(4)*b%vals(6) + a%vals(2)*b%vals(5) + a%vals(5)*b%vals(3)
        res%vals(9) = a%vals(6)*b%vals(6) + a%vals(5)*b%vals(5) + a%vals(3)*b%vals(3)
    end function dot_3D2Osym_3D2Osym

    pure function dot_3D2O_3D2Osym(a, b) result(res)
        !! Computes the single tensor contraction between a general and a 
        !! symmetric second-order tensor.
        !!
        !! Mathematically: \( res_{ij} = a_{ik} b_{kj} \)
        !!
        !! Performance: Manual loop unrolling mapping the 9-component array 
        !! and the 6-component array directly.
        implicit none
        type(ten_3D2O), intent(in) :: a     
            !! General second-order tensor \(\mathbf{A}\) (9 components)
        type(ten_3D2Osym), intent(in) :: b  
            !! Symmetric second-order tensor \(\mathbf{B}\) (6 components)
        type(ten_3D2O) :: res
            !! Resulting general second-order tensor \(\mathbf{res}\) (9 components)
        
        ! Column 1
        res%vals(1) = a%vals(1)*b%vals(1) + a%vals(4)*b%vals(4) + a%vals(7)*b%vals(6)
        res%vals(2) = a%vals(2)*b%vals(1) + a%vals(5)*b%vals(4) + a%vals(8)*b%vals(6)
        res%vals(3) = a%vals(3)*b%vals(1) + a%vals(6)*b%vals(4) + a%vals(9)*b%vals(6)
        
        ! Column 2
        res%vals(4) = a%vals(1)*b%vals(4) + a%vals(4)*b%vals(2) + a%vals(7)*b%vals(5)
        res%vals(5) = a%vals(2)*b%vals(4) + a%vals(5)*b%vals(2) + a%vals(8)*b%vals(5)
        res%vals(6) = a%vals(3)*b%vals(4) + a%vals(6)*b%vals(2) + a%vals(9)*b%vals(5)
        
        ! Column 3
        res%vals(7) = a%vals(1)*b%vals(6) + a%vals(4)*b%vals(5) + a%vals(7)*b%vals(3)
        res%vals(8) = a%vals(2)*b%vals(6) + a%vals(5)*b%vals(5) + a%vals(8)*b%vals(3)
        res%vals(9) = a%vals(3)*b%vals(6) + a%vals(6)*b%vals(5) + a%vals(9)*b%vals(3)
    end function dot_3D2O_3D2Osym

    pure function dot_3D2Osym_3D2O(a, b) result(res)
        !! Computes the single tensor contraction between a symmetric and a 
        !! general second-order tensor.
        !!
        !! Mathematically: \( res_{ij} = a_{ik} b_{kj} \)
        !!
        !! Performance: Manual loop unrolling.
        implicit none
        type(ten_3D2Osym), intent(in) :: a  
            !! Symmetric second-order tensor \(\mathbf{A}\) (6 components)
        type(ten_3D2O), intent(in) :: b     
            !! General second-order tensor \(\mathbf{B}\) (9 components)
        type(ten_3D2O) :: res
            !! Resulting general second-order tensor \(\mathbf{res}\) (9 components)
        
        ! Column 1
        res%vals(1) = a%vals(1)*b%vals(1) + a%vals(4)*b%vals(2) + a%vals(6)*b%vals(3)
        res%vals(2) = a%vals(4)*b%vals(1) + a%vals(2)*b%vals(2) + a%vals(5)*b%vals(3)
        res%vals(3) = a%vals(6)*b%vals(1) + a%vals(5)*b%vals(2) + a%vals(3)*b%vals(3)
        
        ! Column 2
        res%vals(4) = a%vals(1)*b%vals(4) + a%vals(4)*b%vals(5) + a%vals(6)*b%vals(6)
        res%vals(5) = a%vals(4)*b%vals(4) + a%vals(2)*b%vals(5) + a%vals(5)*b%vals(6)
        res%vals(6) = a%vals(6)*b%vals(4) + a%vals(5)*b%vals(5) + a%vals(3)*b%vals(6)
        
        ! Column 3
        res%vals(7) = a%vals(1)*b%vals(7) + a%vals(4)*b%vals(8) + a%vals(6)*b%vals(9)
        res%vals(8) = a%vals(4)*b%vals(7) + a%vals(2)*b%vals(8) + a%vals(5)*b%vals(9)
        res%vals(9) = a%vals(6)*b%vals(7) + a%vals(5)*b%vals(8) + a%vals(3)*b%vals(9)
    end function dot_3D2Osym_3D2O


    ! =========================================================================
    ! ADDITION
    ! =========================================================================

    pure function sum_3D2O_3D2Osym(a, b) result(res)
        !! Computes the addition of a general and a symmetric 2nd-order tensor.
        !!
        !! Mathematically: \( res_{ij} = A_{ij} + B_{ij} \)
        !! The symmetric tensor components are expanded to match the column-major 
        !! layout of the general tensor.
        implicit none
        type(ten_3D2O), intent(in) :: a
            !! General second-order tensor \(\mathbf{A}\)
        type(ten_3D2Osym), intent(in) :: b
            !! Symmetric second-order tensor \(\mathbf{B}\)
        type(ten_3D2O) :: res
            !! Resulting general second-order tensor \(\mathbf{res}\)
        
        res%vals(1) = a%vals(1) + b%vals(1) ! xx
        res%vals(2) = a%vals(2) + b%vals(4) ! yx (sym: xy)
        res%vals(3) = a%vals(3) + b%vals(6) ! zx (sym: xz)
        res%vals(4) = a%vals(4) + b%vals(4) ! xy
        res%vals(5) = a%vals(5) + b%vals(2) ! yy
        res%vals(6) = a%vals(6) + b%vals(5) ! zy (sym: yz)
        res%vals(7) = a%vals(7) + b%vals(6) ! xz
        res%vals(8) = a%vals(8) + b%vals(5) ! yz
        res%vals(9) = a%vals(9) + b%vals(3) ! zz
    end function sum_3D2O_3D2Osym

    pure function sum_3D2Osym_3D2O(a, b) result(res)
        !! Computes the addition of a symmetric and a general 2nd-order tensor.
        !!
        !! Mathematically: \( res_{ij} = A_{ij} + B_{ij} \)
        !! The symmetric tensor components are expanded to match the column-major 
        !! layout of the general tensor.
        implicit none
        type(ten_3D2Osym), intent(in) :: a
            !! Symmetric second-order tensor \(\mathbf{A}\)
        type(ten_3D2O), intent(in) :: b
            !! General second-order tensor \(\mathbf{B}\)
        type(ten_3D2O) :: res
            !! Resulting general second-order tensor \(\mathbf{res}\)
        
        res%vals(1) = a%vals(1) + b%vals(1) ! xx
        res%vals(2) = a%vals(4) + b%vals(2) ! yx (sym: xy)
        res%vals(3) = a%vals(6) + b%vals(3) ! zx (sym: xz)
        res%vals(4) = a%vals(4) + b%vals(4) ! xy
        res%vals(5) = a%vals(2) + b%vals(5) ! yy
        res%vals(6) = a%vals(5) + b%vals(6) ! zy (sym: yz)
        res%vals(7) = a%vals(6) + b%vals(7) ! xz
        res%vals(8) = a%vals(5) + b%vals(8) ! yz
        res%vals(9) = a%vals(3) + b%vals(9) ! zz
    end function sum_3D2Osym_3D2O


    ! =========================================================================
    ! SUBTRACTION
    ! =========================================================================

    pure function sub_3D2O_3D2Osym(a, b) result(res)
        !! Computes the subtraction of a symmetric tensor from a general tensor.
        !!
        !! Mathematically: \( res_{ij} = A_{ij} - B_{ij} \)
        implicit none
        type(ten_3D2O), intent(in) :: a     
            !! General second-order tensor \(\mathbf{A}\)
        type(ten_3D2Osym), intent(in) :: b  
            !! Symmetric second-order tensor \(\mathbf{B}\)
        type(ten_3D2O) :: res
            !! Resulting general second-order tensor \(\mathbf{res}\)
        
        res%vals(1) = a%vals(1) - b%vals(1) ! xx
        res%vals(2) = a%vals(2) - b%vals(4) ! yx - xy
        res%vals(3) = a%vals(3) - b%vals(6) ! zx - xz
        res%vals(4) = a%vals(4) - b%vals(4) ! xy
        res%vals(5) = a%vals(5) - b%vals(2) ! yy
        res%vals(6) = a%vals(6) - b%vals(5) ! zy - yz
        res%vals(7) = a%vals(7) - b%vals(6) ! xz
        res%vals(8) = a%vals(8) - b%vals(5) ! yz
        res%vals(9) = a%vals(9) - b%vals(3) ! zz
    end function sub_3D2O_3D2Osym

    pure function sub_3D2Osym_3D2O(a, b) result(res)
        !! Computes the subtraction of a general tensor from a symmetric tensor.
        !!
        !! Mathematically: \( res_{ij} = A_{ij} - B_{ij} \)
        implicit none
        type(ten_3D2Osym), intent(in) :: a  
            !! Symmetric second-order tensor \(\mathbf{A}\)
        type(ten_3D2O), intent(in) :: b     
            !! General second-order tensor \(\mathbf{B}\)
        type(ten_3D2O) :: res
            !! Resulting general second-order tensor \(\mathbf{res}\)
        
        res%vals(1) = a%vals(1) - b%vals(1) ! xx
        res%vals(2) = a%vals(4) - b%vals(2) ! yx (sym: xy)
        res%vals(3) = a%vals(6) - b%vals(3) ! zx (sym: xz)
        res%vals(4) = a%vals(4) - b%vals(4) ! xy
        res%vals(5) = a%vals(2) - b%vals(5) ! yy
        res%vals(6) = a%vals(5) - b%vals(6) ! zy (sym: yz)
        res%vals(7) = a%vals(6) - b%vals(7) ! xz
        res%vals(8) = a%vals(5) - b%vals(8) ! yz
        res%vals(9) = a%vals(3) - b%vals(9) ! zz
    end function sub_3D2Osym_3D2O


    ! =========================================================================
    ! DOUBLE CONTRACTION (FROBENIUS INNER PRODUCT)
    ! =========================================================================

    pure function ddot_3D2O_3D2Osym(a, b) result(res)
        !! Computes the double contraction between a general and a symmetric tensor.
        !!
        !! Mathematically: \( \alpha = A_{ij} B_{ij} \)
        !! Result is a scalar `real(real64)`.
        !!
        !! Since B is symmetric (\(B_{ij} = B_{ji}\)), we can group terms:
        !! \(\mathbf{A}:\mathbf{B} = A_{11} B_{11} + A_{22} B_{22} + A_{33} B_{33} 
        !!     + (A_{12} + A_{21}) B_{12} + (A_{23} + A_{32}) B_{23} + (A_{13} + A_{31}) B_{13}\)
        implicit none
        type(ten_3D2O), intent(in) :: a     
            !! General second-order tensor \(\mathbf{A}\)
        type(ten_3D2Osym), intent(in) :: b  
            !! Symmetric second-order tensor \(\mathbf{B}\)
        real(real64) :: res
            !! Scalar result \(\alpha\)
        
        res = a%vals(1)*b%vals(1) + a%vals(5)*b%vals(2) + a%vals(9)*b%vals(3) &
            + (a%vals(4) + a%vals(2))*b%vals(4) &  ! (xy + yx) * xy_sym
            + (a%vals(8) + a%vals(6))*b%vals(5) &  ! (yz + zy) * yz_sym
            + (a%vals(7) + a%vals(3))*b%vals(6)    ! (xz + zx) * xz_sym
    end function ddot_3D2O_3D2Osym

    pure function ddot_3D2Osym_3D2O(a, b) result(res)
        !! Computes the double contraction between a symmetric and a general tensor.
        !!
        !! Mathematically: \( \alpha = A_{ij} B_{ij} \)
        !! Due to the commutative property of the inner product, this is identical
        !! to the computation in `ddot_3D2O_3D2Osym`.
        implicit none
        type(ten_3D2Osym), intent(in) :: a  
            !! Symmetric second-order tensor \(\mathbf{A}\)
        type(ten_3D2O), intent(in) :: b     
            !! General second-order tensor \(\mathbf{B}\)
        real(real64) :: res
            !! Scalar result \(\alpha\)
        
        res = a%vals(1)*b%vals(1) + a%vals(2)*b%vals(5) + a%vals(3)*b%vals(9) &
            + a%vals(4)*(b%vals(4) + b%vals(2)) &  ! xy_sym * (xy + yx)
            + a%vals(5)*(b%vals(8) + b%vals(6)) &  ! yz_sym * (yz + zy)
            + a%vals(6)*(b%vals(7) + b%vals(3))    ! xz_sym * (xz + zx)
    end function ddot_3D2Osym_3D2O


    ! =========================================================================
    ! ASSIGNMENT
    ! =========================================================================

    pure subroutine assign_3D2O_3D2Osym(a, b)
        !! Explicit assignment from a symmetric tensor to a general tensor.
        !! 
        !! Unpacks the 6-component array of the symmetric tensor into the 
        !! 9-component column-major array of the general tensor.
        implicit none
        type(ten_3D2O), intent(out) :: a
            !! The target general tensor to be overwritten.
        type(ten_3D2Osym), intent(in) :: b  
            !! The source symmetric tensor.
        
        a%vals(1) = b%vals(1) ! xx
        a%vals(2) = b%vals(4) ! yx <- xy
        a%vals(3) = b%vals(6) ! zx <- xz
        a%vals(4) = b%vals(4) ! xy
        a%vals(5) = b%vals(2) ! yy
        a%vals(6) = b%vals(5) ! zy <- yz
        a%vals(7) = b%vals(6) ! xz
        a%vals(8) = b%vals(5) ! yz
        a%vals(9) = b%vals(3) ! zz
    end subroutine assign_3D2O_3D2Osym

end module mod_operator_3D2Osym_3D2O