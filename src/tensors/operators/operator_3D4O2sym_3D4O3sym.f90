module mod_operator_3D4O2sym_3D4O3sym
        !! Module mod_operator_3D4O2sym_3D4O3sym
    !! ====================================
    !!
    !! Defines mixed algebraic operations involving 3D fully symmetric fourth-order
    !! tensors (`ten_3D4O3sym`) and 3D fourth-order tensors with minor symmetries (`ten_3D4O2sym`).
    !!
    !! This module provides:
    !! - Double contraction (`.ddot.`) between two fully symmetric tensors.
    !! - Assignment (`=`) to cast a fully symmetric tensor into a minor-symmetric one.
    !!
    !! Overloaded Operators & Assignments
    !! ----------------------------------
    !!
    !! - `.ddot.` : Double tensor contraction between two `ten_3D4O3sym` tensors.
    !!   - \(\mathbb{C} = \mathbb{A} : \mathbb{B}\) (`3D4O3sym .ddot. 3D4O3sym` -> `3D4O2sym`)
    !!     *Note: This computes the standard matrix product of the 6x6 Voigt representations.*
    !!
    !! - `=` : Assignment (Type casting).
    !!   - `ten_3D4O2sym = ten_3D4O3sym`
    !!
    !! For tensor type definitions, see [[tensors_types]].
    use, intrinsic :: iso_fortran_env
    use mod_ten_3D4O3sym
    use mod_ten_3D4O2sym
    implicit none
    private


    public :: operator(.ddot.)
    interface operator (.ddot.)
        module procedure ddot_3D4O3sym_3D4O3sym
    end interface

     public :: assignment (=)
     interface assignment (=)
         module procedure assign_3D4O2sym_3D4O3sym
     end interface

    contains

    pure function ddot_3D4O3sym_3D4O3sym(aT, bT) result(res)
        !! Computes the double contraction of two fully symmetric 4th-order tensors.
        !! \([\mathbf{res}]_{ijmn} = [\mathbf{A}]_{ijkl} [\mathbf{B}]_{klmn}\).
        !!
        !! This implementation is fully unrolled for maximum performance, avoiding temporary arrays
        !! and leveraging compiler optimizations like register allocation and SIMD.
        !!
        implicit none
        class(ten_3D4O3sym), intent(in) :: aT
        class(ten_3D4O3sym), intent(in) :: bT
        type(ten_3D4O2sym) :: res
        
        real(real64) :: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, x17, x18, x19
        real(real64) :: x20, x21, x22, x23, x24, x25, x26, x27, x28, x29

        x0 = aT%vals(7)*bT%vals(7)
        x1 = aT%vals(12)*bT%vals(12)
        x2 = aT%vals(16)*bT%vals(16)
        x3 = aT%vals(21)*bT%vals(21)
        x4 = aT%vals(19)*bT%vals(19)
        x5 = 2*aT%vals(16)
        x6 = 2*aT%vals(21)
        x7 = 2*aT%vals(19)
        x8 = 2*aT%vals(13)
        x9 = 2*aT%vals(20)
        x10 = 2*aT%vals(17)
        x11 = aT%vals(8)*bT%vals(8)
        x12 = aT%vals(13)*bT%vals(13)
        x13 = aT%vals(20)*bT%vals(20)
        x14 = aT%vals(17)*bT%vals(17)
        x15 = 2*aT%vals(9)
        x16 = 2*aT%vals(18)
        x17 = 2*aT%vals(14)
        x18 = aT%vals(9)*bT%vals(9)
        x19 = aT%vals(18)*bT%vals(18)
        x20 = aT%vals(14)*bT%vals(14)
        x21 = 2*aT%vals(4)
        x22 = 2*aT%vals(15)
        x23 = 2*aT%vals(10)
        x24 = bT%vals(15)*x22
        x25 = bT%vals(10)*x23
        x26 = 2*aT%vals(11)
        x27 = 2*aT%vals(5)
        x28 = bT%vals(11)*x26
        x29 = 2*aT%vals(6)
        res%vals(1, 1) = aT%vals(1)*bT%vals(1) + x0 + x1 + 2*x2 + 2*x3 + 2*x4
        res%vals(2, 1) = aT%vals(7)*bT%vals(1) + aT%vals(2)*bT%vals(7) + aT%vals(8)*bT%vals(12) + bT%vals(16)*x8 + &
                         bT%vals(21)*x9 + bT%vals(19)*x10
        res%vals(3, 1) = aT%vals(12)*bT%vals(1) + aT%vals(8)*bT%vals(7) + aT%vals(3)*bT%vals(12) + bT%vals(16)*x15 + &
                         bT%vals(21)*x16 + bT%vals(19)*x17
        res%vals(4, 1) = aT%vals(16)*bT%vals(1) + aT%vals(13)*bT%vals(7) + aT%vals(9)*bT%vals(12) + bT%vals(16)*x21 + &
                         bT%vals(21)*x22 + bT%vals(19)*x23
        res%vals(5, 1) = aT%vals(19)*bT%vals(1) + aT%vals(17)*bT%vals(7) + aT%vals(14)*bT%vals(12) + bT%vals(16)*x23 + &
                         bT%vals(21)*x26 + bT%vals(19)*x27
        res%vals(6, 1) = aT%vals(21)*bT%vals(1) + aT%vals(20)*bT%vals(7) + aT%vals(18)*bT%vals(12) + bT%vals(16)*x22 + &
                         bT%vals(21)*x29 + bT%vals(19)*x26
        res%vals(1, 2) = aT%vals(1)*bT%vals(7) + aT%vals(7)*bT%vals(2) + aT%vals(12)*bT%vals(8) + bT%vals(13)*x5 + &
                         bT%vals(20)*x6 + bT%vals(17)*x7
        res%vals(2, 2) = aT%vals(2)*bT%vals(2) + x0 + x11 + 2*x12 + 2*x13 + 2*x14
        res%vals(3, 2) = aT%vals(12)*bT%vals(7) + aT%vals(8)*bT%vals(2) + aT%vals(3)*bT%vals(8) + bT%vals(13)*x15 + &
                         bT%vals(20)*x16 + bT%vals(17)*x17
        res%vals(4, 2) = aT%vals(16)*bT%vals(7) + aT%vals(13)*bT%vals(2) + aT%vals(9)*bT%vals(8) + bT%vals(13)*x21 + &
                         bT%vals(20)*x22 + bT%vals(17)*x23
        res%vals(5, 2) = aT%vals(19)*bT%vals(7) + aT%vals(17)*bT%vals(2) + aT%vals(14)*bT%vals(8) + bT%vals(13)*x23 + &
                         bT%vals(20)*x26 + bT%vals(17)*x27
        res%vals(6, 2) = aT%vals(21)*bT%vals(7) + aT%vals(20)*bT%vals(2) + aT%vals(18)*bT%vals(8) + bT%vals(13)*x22 + &
                         bT%vals(20)*x29 + bT%vals(17)*x26
        res%vals(1, 3) = aT%vals(1)*bT%vals(12) + aT%vals(7)*bT%vals(8) + aT%vals(12)*bT%vals(3) + bT%vals(9)*x5 + &
                         bT%vals(18)*x6 + bT%vals(14)*x7
        res%vals(2, 3) = aT%vals(7)*bT%vals(12) + aT%vals(2)*bT%vals(8) + aT%vals(8)*bT%vals(3) + bT%vals(9)*x8 + &
                         bT%vals(18)*x9 + bT%vals(14)*x10
        res%vals(3, 3) = aT%vals(3)*bT%vals(3) + x1 + x11 + 2*x18 + 2*x19 + 2*x20
        res%vals(4, 3) = aT%vals(16)*bT%vals(12) + aT%vals(13)*bT%vals(8) + aT%vals(9)*bT%vals(3) + bT%vals(9)*x21 + &
                         bT%vals(18)*x22 + bT%vals(14)*x23
        res%vals(5, 3) = aT%vals(19)*bT%vals(12) + aT%vals(17)*bT%vals(8) + aT%vals(14)*bT%vals(3) + bT%vals(9)*x23 + &
                         bT%vals(18)*x26 + bT%vals(14)*x27
        res%vals(6, 3) = aT%vals(21)*bT%vals(12) + aT%vals(20)*bT%vals(8) + aT%vals(18)*bT%vals(3) + bT%vals(9)*x22 + &
                         bT%vals(18)*x29 + bT%vals(14)*x26
        res%vals(1, 4) = aT%vals(1)*bT%vals(16) + aT%vals(7)*bT%vals(13) + aT%vals(12)*bT%vals(9) + bT%vals(4)*x5 + &
                         bT%vals(15)*x6 + bT%vals(10)*x7
        res%vals(2, 4) = aT%vals(7)*bT%vals(16) + aT%vals(2)*bT%vals(13) + aT%vals(8)*bT%vals(9) + bT%vals(4)*x8 + &
                         bT%vals(15)*x9 + bT%vals(10)*x10
        res%vals(3, 4) = aT%vals(12)*bT%vals(16) + aT%vals(8)*bT%vals(13) + aT%vals(3)*bT%vals(9) + bT%vals(4)*x15 + &
                         bT%vals(15)*x16 + bT%vals(10)*x17
        res%vals(4, 4) = bT%vals(4)*x21 + x12 + x18 + x2 + x24 + x25
        res%vals(5, 4) = aT%vals(19)*bT%vals(16) + aT%vals(17)*bT%vals(13) + aT%vals(14)*bT%vals(9) + bT%vals(4)*x23 + &
                         bT%vals(15)*x26 + bT%vals(10)*x27
        res%vals(6, 4) = aT%vals(21)*bT%vals(16) + aT%vals(20)*bT%vals(13) + aT%vals(18)*bT%vals(9) + bT%vals(4)*x22 + &
                         bT%vals(15)*x29 + bT%vals(10)*x26
        res%vals(1, 5) = aT%vals(1)*bT%vals(19) + aT%vals(7)*bT%vals(17) + aT%vals(12)*bT%vals(14) + bT%vals(10)*x5 + &
                         bT%vals(11)*x6 + bT%vals(5)*x7
        res%vals(2, 5) = aT%vals(7)*bT%vals(19) + aT%vals(2)*bT%vals(17) + aT%vals(8)*bT%vals(14) + bT%vals(10)*x8 + &
                         bT%vals(11)*x9 + bT%vals(5)*x10
        res%vals(3, 5) = aT%vals(12)*bT%vals(19) + aT%vals(8)*bT%vals(17) + aT%vals(3)*bT%vals(14) + bT%vals(10)*x15 + &
                         bT%vals(11)*x16 + bT%vals(5)*x17
        res%vals(4, 5) = aT%vals(16)*bT%vals(19) + aT%vals(13)*bT%vals(17) + aT%vals(9)*bT%vals(14) + bT%vals(10)*x21 + &
                         bT%vals(11)*x22 + bT%vals(5)*x23
        res%vals(5, 5) = bT%vals(5)*x27 + x14 + x20 + x25 + x28 + x4
        res%vals(6, 5) = aT%vals(21)*bT%vals(19) + aT%vals(20)*bT%vals(17) + aT%vals(18)*bT%vals(14) + bT%vals(10)*x22 + &
                         bT%vals(11)*x29 + bT%vals(5)*x26
        res%vals(1, 6) = aT%vals(1)*bT%vals(21) + aT%vals(7)*bT%vals(20) + aT%vals(12)*bT%vals(18) + bT%vals(15)*x5 + &
                         bT%vals(6)*x6 + bT%vals(11)*x7
        res%vals(2, 6) = aT%vals(7)*bT%vals(21) + aT%vals(2)*bT%vals(20) + aT%vals(8)*bT%vals(18) + bT%vals(15)*x8 + &
                         bT%vals(6)*x9 + bT%vals(11)*x10
        res%vals(3, 6) = aT%vals(12)*bT%vals(21) + aT%vals(8)*bT%vals(20) + aT%vals(3)*bT%vals(18) + bT%vals(15)*x15 + &
                         bT%vals(6)*x16 + bT%vals(11)*x17
        res%vals(4, 6) = aT%vals(16)*bT%vals(21) + aT%vals(13)*bT%vals(20) + aT%vals(9)*bT%vals(18) + bT%vals(15)*x21 + &
                         bT%vals(6)*x22 + bT%vals(11)*x23
        res%vals(5, 6) = aT%vals(19)*bT%vals(21) + aT%vals(17)*bT%vals(20) + aT%vals(14)*bT%vals(18) + bT%vals(15)*x23 + &
                         bT%vals(6)*x26 + bT%vals(11)*x27
        res%vals(6, 6) = bT%vals(6)*x29 + x13 + x19 + x24 + x28 + x3

    end function ddot_3D4O3sym_3D4O3sym


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
    end subroutine

end module mod_operator_3D4O2sym_3D4O3sym