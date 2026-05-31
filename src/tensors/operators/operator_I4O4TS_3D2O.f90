module mod_operator_I4O4TS_3D2O
    !! Module mod_operator_I4O4TS_3D2O
    !! ===============================
    !!
    !! Defines mixed algebraic operations between the scaled 3D fourth-order 
    !! symmetric identity tensor (`iden_4O4TS`, \(c\mathbb{I}^S\)) and 
    !! general 3D second-order tensors (`ten_3D2O`).
    !!
    !! The double contraction (\(c\mathbb{I}^S : \mathbf{A}\)) extracts the 
    !! symmetric part of the general tensor and scales it by the factor \(c\):
    !! \[ \mathbf{res} = \frac{c}{2}(\mathbf{A} + \mathbf{A}^T) \]

    use, intrinsic :: iso_fortran_env
    use mod_iden_4O4TS
    use mod_ten_3D2Osym
    use mod_ten_3D2O
    implicit none
    private
    
    public :: operator(.ddot.)
    interface operator (.ddot.)
        module procedure ddot_I4O4TS_3D2O
        module procedure ddot_3D2O_I4O4TS
    end interface

contains

    pure function ddot_I4O4TS_3D2O(I4S, a) result(res)
        !! Computes the scaled double contraction product: \(\mathbf{res} = c\mathbb{I}^S : \mathbf{A}\).
        !!
        !! Mathematically, this yields the scaled symmetric part of tensor \(\mathbf{A}\).
        !!
        !! Mapping logic from 9-component (column-major) to 6-component (Voigt):
        !! - res(1) = c * A(11)
        !! - res(2) = c * A(22)
        !! - res(3) = c * A(33)
        !! - res(4) = 0.5 * c * (A(12) + A(21))
        !! - res(5) = 0.5 * c * (A(23) + A(32))
        !! - res(6) = 0.5 * c * (A(13) + A(31))
        implicit none
        type(iden_4O4TS), intent(in) :: I4S
            !! Scaled fourth-order symmetric identity tensor (\(c\mathbb{I}^S\)).
        type(ten_3D2O),   intent(in) :: a
            !! General second-order tensor (9 components, column-major).
        type(ten_3D2Osym)           :: res
            !! Resulting symmetric second-order tensor (6 components, Voigt).
        
        real(real64) :: c, half_c

        c = I4S%val
        half_c = 0.5D0 * c

        ! Diagonal components: c * A_ii
        res%vals(1) = a%vals(1) * c ! xx
        res%vals(2) = a%vals(5) * c ! yy
        res%vals(3) = a%vals(9) * c ! zz
        
        ! Symmetric off-diagonal components: 0.5 * c * (A_ij + A_ji)
        ! Note indices for ten_3D2O: (11:1, 21:2, 31:3, 12:4, 22:5, 32:6, 13:7, 23:8, 33:9)
        res%vals(4) = (a%vals(2) + a%vals(4)) * half_c ! xy part
        res%vals(5) = (a%vals(6) + a%vals(8)) * half_c ! yz part
        res%vals(6) = (a%vals(3) + a%vals(7)) * half_c ! xz part
        
    end function ddot_I4O4TS_3D2O

    pure function ddot_3D2O_I4O4TS(a, I4S) result(res)
        !! Computes the scaled double contraction product: \(\mathbf{res} = \mathbf{A} : c\mathbb{I}^S\).
        !!
        !! Equivalent to \(c\mathbb{I}^S : \mathbf{A}\) due to symmetry of the contraction.
        implicit none
        type(ten_3D2O),   intent(in) :: a
        type(iden_4O4TS), intent(in) :: I4S
        type(ten_3D2Osym)           :: res
        
        ! Delegate to the primary implementation to ensure consistency
        res = ddot_I4O4TS_3D2O(I4S, a)
        
    end function ddot_3D2O_I4O4TS

end module mod_operator_I4O4TS_3D2O