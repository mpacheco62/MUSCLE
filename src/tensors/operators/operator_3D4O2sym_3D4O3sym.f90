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

     public :: assignment (=)
     interface assignment (=)
         module procedure assign_3D4O2sym_3D4O3sym
     end interface

    contains

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