module mod_operator_2D2Osym_2D4O3sym
    !! Module mod_operator_2D2Osym_2D4O3sym
    !! ====================================
    !!
    !! Defines mixed algebraic operations between 2D symmetric second-order
    !! tensors (`ten_2D2Osym`) and 2D fully symmetric fourth-order tensors
    !! (`ten_2D4O3sym`).
    !!
    !! This module overloads Fortran intrinsic operators to allow direct and
    !! readable computation of constitutive laws in 2D mechanics (e.g., plane
    !! strain or plane stress). It internally manages the Voigt notation mapping,
    !! including the necessary factor of 2 for shear components during tensor
    !! contractions.
    !!
    !! Overloaded Operators
    !! --------------------
    !!
    !! - `.ddot.` : Double tensor contraction (Inner product).
    !!   - \(\boldsymbol{\sigma} = \mathbb{C} : \boldsymbol{\varepsilon}\) (`2D4O3sym .ddot. 2D2Osym` -> `2D2Osym`)
    !!   - \(\boldsymbol{\varepsilon} = \boldsymbol{\sigma} : \mathbb{S}\) (`2D2Osym .ddot. 2D4O3sym` -> `2D2Osym`)
    !!
    !! - `.tdot.` : Dyadic tensor product (Outer product).
    !!   - \(\mathbb{P} = \mathbf{a} \otimes \mathbf{b}\) (`2D2Osym .tdot. 2D2Osym` -> `2D4O3sym`)
    !!     *Note: Results in a fourth-order tensor assuming major symmetry.*

    use, intrinsic :: iso_fortran_env
    use mod_ten_2D2Osym
    use mod_ten_2D4O3sym
    implicit none
    private

    public :: operator(.tdot.)
    interface operator (.tdot.)
        module procedure tdot_2D2Osym_2D2Osym
    end interface

contains

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
        class(ten_2D2Osym), intent(in) :: a
            !! First second-order symmetric tensor \(\mathbf{a}\)
        class(ten_2D2Osym), intent(in) :: b
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

end module mod_operator_2D2Osym_2D4O3sym