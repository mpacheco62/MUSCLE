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

    public :: operator(.ddot.)
    interface operator (.ddot.)
        module procedure ddot_2D4O3sym_2D2Osym
        module procedure ddot_2D2Osym_2D4O3sym
    end interface

    public :: operator(.tdot.)
    interface operator (.tdot.)
        module procedure tdot_2D2Osym_2D2Osym
    end interface

contains

    pure function ddot_2D4O3sym_2D2Osym(a, b) result(res)
        !! Computes the double contraction product: \(\mathbf{res} = \mathbb{A} : \mathbf{b}\).
        !!
        !! Mathematically: \( res_{ij} = A_{ijkl} b_{kl} \).
        !!
        !! Performance Optimization:
        !! -------------------------
        !! Performs manual loop unrolling of the contraction using Voigt notation 
        !! (4x4 matrix contracted with a 4x1 vector). Shear terms are implicitly 
        !! multiplied by 2 due to the nature of symmetric tensor contractions.
        !!
        !! Voigt Mapping Reference:
        !! ```
        !!  | ( 1:1111) ( 5:1122) ( 8:1133) (10:1112) |
        !!  | ( 5:2211) ( 2:2222) ( 6:2233) ( 9:2212) |
        !!  | ( 8:3311) ( 6:3322) ( 3:3333) ( 7:3312) |
        !!  | (10:1211) ( 9:1222) ( 7:1233) ( 4:1212) |
        !! ```
        implicit none
        class(ten_2D4O3sym), intent(in) :: a
            !! Fourth-order fully symmetric tensor \(\mathbb{A}\)
        class(ten_2D2Osym), intent(in) :: b
            !! Second-order symmetric tensor \(\mathbf{b}\)
        type(ten_2D2Osym) :: res
            !! Resulting second-order symmetric tensor \(\mathbf{res}\)
        
        real(real64) :: a1111b11, a1122b22, a1133b33, a1112b12
        real(real64) :: a2211b11, a2222b22, a2233b33, a2212b12
        real(real64) :: a3311b11, a3322b22, a3333b33, a3312b12
        real(real64) :: a1112b11, a2212b22, a3312b33, a1212b12

        a1111b11 = a%vals( 1)*b%vals(1)
        a1122b22 = a%vals( 5)*b%vals(2)
        a1133b33 = a%vals( 8)*b%vals(3)
        a1112b12 = a%vals(10)*b%vals(4)

        a2211b11 = a%vals( 5)*b%vals(1)
        a2222b22 = a%vals( 2)*b%vals(2)
        a2233b33 = a%vals( 6)*b%vals(3)
        a2212b12 = a%vals( 9)*b%vals(4)

        a3311b11 = a%vals( 8)*b%vals(1)
        a3322b22 = a%vals( 6)*b%vals(2)
        a3333b33 = a%vals( 3)*b%vals(3)
        a3312b12 = a%vals( 7)*b%vals(4)

        a1112b11 = a%vals(10)*b%vals(1)
        a2212b22 = a%vals( 9)*b%vals(2)
        a3312b33 = a%vals( 7)*b%vals(3)
        a1212b12 = a%vals( 4)*b%vals(4)
        
        res%vals(1) = a1111b11 + a1122b22 + a1133b33 + 2.0D0 * a1112b12
        res%vals(2) = a2211b11 + a2222b22 + a2233b33 + 2.0D0 * a2212b12
        res%vals(3) = a3311b11 + a3322b22 + a3333b33 + 2.0D0 * a3312b12
        res%vals(4) = a1112b11 + a2212b22 + a3312b33 + 2.0D0 * a1212b12
        
    end function ddot_2D4O3sym_2D2Osym

    pure function ddot_2D2Osym_2D4O3sym(b, a) result(res)
        !! Computes the double contraction product: \(\mathbf{res} = \mathbf{b} : \mathbb{A}\).
        !!
        !! Mathematically: \( res_{ij} = b_{kl} A_{klij} \).
        !! Due to the major symmetry of \(\mathbb{A}\), this is equivalent to \(\mathbb{A} : \mathbf{b}\).
        !!
        !! Performance Optimization:
        !! -------------------------
        !! Performs manual loop unrolling of the contraction using Voigt notation 
        !! (1x4 vector contracted with a 4x4 matrix).
        !!
        !! Voigt Mapping Reference:
        !! ```
        !!  | ( 1:1111) ( 5:1122) ( 8:1133) (10:1112) |
        !!  | ( 5:2211) ( 2:2222) ( 6:2233) ( 9:2212) |
        !!  | ( 8:3311) ( 6:3322) ( 3:3333) ( 7:3312) |
        !!  | (10:1211) ( 9:1222) ( 7:1233) ( 4:1212) |
        !! ```
        implicit none
        class(ten_2D2Osym), intent(in) :: b
            !! Second-order symmetric tensor \(\mathbf{b}\)
        class(ten_2D4O3sym), intent(in) :: a
            !! Fourth-order fully symmetric tensor \(\mathbb{A}\)
        type(ten_2D2Osym) :: res
            !! Resulting second-order symmetric tensor \(\mathbf{res}\)
        
        real(real64) :: a1111b11, a2211b22, a3311b33, a1211b12
        real(real64) :: a1122b11, a2222b22, a3322b33, a1222b12
        real(real64) :: a1133b11, a2233b22, a3333b33, a1233b12
        real(real64) :: a1112b11, a2212b22, a3312b33, a1212b12

        a1111b11 = a%vals( 1)*b%vals(1)
        a2211b22 = a%vals( 5)*b%vals(2)
        a3311b33 = a%vals( 8)*b%vals(3)
        a1211b12 = a%vals(10)*b%vals(4)

        a1122b11 = a%vals( 5)*b%vals(1)
        a2222b22 = a%vals( 2)*b%vals(2)
        a3322b33 = a%vals( 6)*b%vals(3)
        a1222b12 = a%vals( 9)*b%vals(4)

        a1133b11 = a%vals( 8)*b%vals(1)
        a2233b22 = a%vals( 6)*b%vals(2)
        a3333b33 = a%vals( 3)*b%vals(3)
        a1233b12 = a%vals( 7)*b%vals(4)

        a1112b11 = a%vals(10)*b%vals(1)
        a2212b22 = a%vals( 9)*b%vals(2)
        a3312b33 = a%vals( 7)*b%vals(3)
        a1212b12 = a%vals( 4)*b%vals(4)
        
        res%vals(1) = a1111b11 + a2211b22 + a3311b33 + 2.0D0 * a1211b12
        res%vals(2) = a1122b11 + a2222b22 + a3322b33 + 2.0D0 * a1222b12
        res%vals(3) = a1133b11 + a2233b22 + a3333b33 + 2.0D0 * a1233b12
        res%vals(4) = a1112b11 + a2212b22 + a3312b33 + 2.0D0 * a1212b12
   
    end function ddot_2D2Osym_2D4O3sym

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