module mod_ops_assignment
    !! This module defines the overloaded assignment operator `=` for mixed-type conversions
    !! and safe casting between different second-order, fourth-order, and identity tensors.
    !!
    !! Overloading assignment allows operations like `tensor_general = tensor_symmetric` 
    !! or `tensor_symmetric = identity_scaled` to be handled transparently.
    !! It automatically maps the storage formats and handles any necessary scaling
    !! (e.g., setting the Voigt diagonal components while zeroing out shear terms).

    use, intrinsic :: iso_fortran_env, only : real64
    
    ! Import necessary types
    use mod_iden_2O
    use mod_iden_2OS
    use mod_iden_4O3TS
    use mod_iden_4O4T
    use mod_iden_4O4TS
    use mod_ten_3D2O
    use mod_ten_3D2Osym
    use mod_ten_3D4O2sym
    use mod_ten_3D4O3sym

    implicit none
    private

    ! =========================================================================
    ! PUBLIC INTERFACES
    ! =========================================================================
    
    public :: assignment(=)
    interface assignment (=)
        ! --- 2nd Order Assignments (Castings) ---
        module procedure assign_3D2O_3D2Osym
        module procedure assign_3D2Osym_I2O
        module procedure assign_3D2Osym_I2OS

        ! --- 4th Order Assignments (Castings) ---
        module procedure assign_3D4O2sym_3D4O3sym
        module procedure assign_3D4O2sym_I4O4T
        module procedure assign_3D4O2sym_I4O4TS
        module procedure assign_3D4O3sym_I4O3TS
        module procedure assign_3D4O3sym_I4O4TS
    end interface

contains

    ! =========================================================================
    ! 1. 2ND ORDER ASSIGNMENTS
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

    pure subroutine assign_3D2Osym_I2OS(a, b)
        !! Explicit assignment from a scaled second order identity to a 3D symmetric tensor.
        !! 
        !! This is necessary to allow statements like `a = I2` where `a` is a `ten_3D2Osym` and `I2` is an `iden_2OS`.
        !! The resulting tensor `a` will have its diagonal components set to the scale factor and off-diagonal components set to 0.
        implicit none
        type(ten_3D2Osym), intent(out) :: a
            !! The target symmetric tensor to be overwritten.
        type(iden_2OS), intent(in) :: b
            !! The source scaled second-order identity tensor.

        a%vals(1) = b%val ! xx
        a%vals(2) = b%val ! yy
        a%vals(3) = b%val ! zz
        a%vals(4) = 0.0D0 ! xy
        a%vals(5) = 0.0D0 ! yz
        a%vals(6) = 0.0D0 ! xz
    end subroutine assign_3D2Osym_I2OS

    ! =========================================================================
    ! 2. 4TH ORDER ASSIGNMENTS
    ! =========================================================================

    pure subroutine assign_3D4O2sym_3D4O3sym(self, b)
        !! Assigns a fully symmetric tensor (21 components) to a minor-symmetric
        !! tensor (36 components) by expanding it to a full 6x6 Voigt matrix.
        implicit none
        class(ten_3D4O2sym), intent(inout) :: self
            !! The target minor-symmetric tensor.
        class(ten_3D4O3sym), intent(in) :: b
            !! The source fully symmetric tensor.
        self%vals(1,1)=b%vals(1) ; self%vals(1,2)=b%vals(7) ;  self%vals(1,3)=b%vals(12)
        self%vals(1,4)=b%vals(16); self%vals(1,5)=b%vals(19);  self%vals(1,6)=b%vals(21)
        self%vals(2,1)=b%vals(7) ; self%vals(2,2)=b%vals(2) ;  self%vals(2,3)=b%vals(8)
        self%vals(2,4)=b%vals(13); self%vals(2,5)=b%vals(17);  self%vals(2,6)=b%vals(20)
        self%vals(3,1)=b%vals(12); self%vals(3,2)=b%vals(8) ;  self%vals(3,3)=b%vals(3)
        self%vals(3,4)=b%vals(9) ; self%vals(3,5)=b%vals(14);  self%vals(3,6)=b%vals(18)
        self%vals(4,1)=b%vals(16); self%vals(4,2)=b%vals(13);  self%vals(4,3)=b%vals(9) 
        self%vals(4,4)=b%vals(4) ; self%vals(4,5)=b%vals(10);  self%vals(4,6)=b%vals(15)
        self%vals(5,1)=b%vals(19); self%vals(5,2)=b%vals(17);  self%vals(5,3)=b%vals(14)
        self%vals(5,4)=b%vals(10); self%vals(5,5)=b%vals(5) ;  self%vals(5,6)=b%vals(11)
        self%vals(6,1)=b%vals(21); self%vals(6,2)=b%vals(20);  self%vals(6,3)=b%vals(18)
        self%vals(6,4)=b%vals(15); self%vals(6,5)=b%vals(11);  self%vals(6,6)=b%vals(6)
    end subroutine assign_3D4O2sym_3D4O3sym

    pure subroutine assign_3D4O2sym_I4O4T(a, I4)
        !! Assigns the symmetric identity tensor directly to a fourth-order tensor: \(\mathbf{A} = \mathbb{I}^S\).
        implicit none
        type(ten_3D4O2sym), intent(out) :: a
        type(iden_4O4T),   intent(in)  :: I4
        
        a%vals = 0.0D0
        a%vals(1,1) = 1.0D0
        a%vals(2,2) = 1.0D0
        a%vals(3,3) = 1.0D0
        a%vals(4,4) = 0.5D0
        a%vals(5,5) = 0.5D0
        a%vals(6,6) = 0.5D0
    end subroutine assign_3D4O2sym_I4O4T

    pure subroutine assign_3D4O2sym_I4O4TS(a, I4S)
        !! Directly assigns a scaled symmetric identity to a 4th-order tensor: \(\mathbb{A} = c\mathbb{I}^S\).
        implicit none
        type(ten_3D4O2sym), intent(out) :: a
        type(iden_4O4TS),  intent(in)  :: I4S
        
        real(real64) :: c, half_c

        c = I4S%val
        half_c = 0.5D0 * c

        a%vals = 0.0D0
        a%vals(1,1) = c
        a%vals(2,2) = c
        a%vals(3,3) = c
        a%vals(4,4) = half_c
        a%vals(5,5) = half_c
        a%vals(6,6) = half_c
    end subroutine assign_3D4O2sym_I4O4TS

    pure subroutine assign_3D4O3sym_I4O3TS(a, b)
        implicit none
        type(ten_3D4O3sym), intent(out) :: a
            !! The target general tensor to be overwritten.
        type(iden_4O3TS), intent(in) :: b  
            !! The source symmetric tensor.
        
        a%vals = 0D0
        a%vals(1:3) = b%val
        a%vals(7:8) = b%val
        a%vals(12)  = b%val
    end subroutine assign_3D4O3sym_I4O3TS

    pure subroutine assign_3D4O3sym_I4O4TS(a, b)
        implicit none
        type(ten_3D4O3sym), intent(out) :: a
            !! The target general tensor to be overwritten.
        type(iden_4O4TS), intent(in) :: b  
            !! The source symmetric tensor.
        real(real64) :: c, half_c

        c = b%val
        half_c = 0.5D0 * c

        a%vals = 0D0
        ! Add scaled identity components to the Voigt diagonal
        a%vals(1:3) = c
        a%vals(4:6) = half_c
    end subroutine assign_3D4O3sym_I4O4TS

end module mod_ops_assignment