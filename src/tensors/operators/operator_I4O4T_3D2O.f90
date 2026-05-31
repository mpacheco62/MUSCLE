module mod_operator_I4O4T_3D2O
    !! Module mod_operator_I4O4T_3D2O
    !! ==============================
    !!
    !! Defines mixed algebraic operations between the standard 3D fourth-order 
    !! symmetric identity tensor (`iden_4O4T`) and general (non-symmetric) 
    !! 3D second-order tensors (`ten_3D2O`).
    !!
    !! The primary operation is the double contraction (\(\mathbb{I}^S : \mathbf{A}\)),
    !! which mathematically extracts the symmetric part of a general tensor:
    !! \[ \mathbf{res} = \frac{1}{2}(\mathbf{A} + \mathbf{A}^T) \]
    !!
    !! This is particularly useful in mechanics when projecting non-symmetric 
    !! gradients (like the velocity gradient) into symmetric measures (like the 
    !! rate of deformation tensor).

    use, intrinsic :: iso_fortran_env
    use mod_iden_4O4T
    use mod_ten_3D2Osym
    use mod_ten_3D2O
    implicit none
    private
    
    public :: operator(.ddot.)
    interface operator (.ddot.)
        module procedure ddot_I4O4T_3D2O
        module procedure ddot_3D2O_I4O4T
    end interface

contains

    pure function ddot_I4O4T_3D2O(I4, a) result(res)
        !! Computes the double contraction product: \(\mathbf{res} = \mathbb{I}^S : \mathbf{A}\).
        !!
        !! This results in the symmetric part of tensor \(\mathbf{A}\).
        !!
        !! Mapping logic from 9-component (column-major) to 6-component (Voigt):
        !! - res(11) = A(11)
        !! - res(22) = A(22)
        !! - res(33) = A(33)
        !! - res(12) = 0.5 * (A(12) + A(21))
        !! - res(23) = 0.5 * (A(23) + A(32))
        !! - res(13) = 0.5 * (A(13) + A(31))
        implicit none
        type(iden_4O4T), intent(in) :: I4
            !! Standard fourth-order symmetric identity tensor.
        type(ten_3D2O),  intent(in) :: a
            !! General second-order tensor (9 components, column-major).
        type(ten_3D2Osym)           :: res
            !! Resulting symmetric second-order tensor (6 components, Voigt).

        ! Diagonal components
        res%vals(1) = a%vals(1) ! xx
        res%vals(2) = a%vals(5) ! yy
        res%vals(3) = a%vals(9) ! zz
        
        ! Symmetric off-diagonal components: 0.5 * (A_ij + A_ji)
        ! Note: ten_3D2O storage is (11, 21, 31, 12, 22, 32, 13, 23, 33)
        res%vals(4) = 0.5D0 * (a%vals(2) + a%vals(4)) ! xy part: (21 + 12)/2
        res%vals(5) = 0.5D0 * (a%vals(6) + a%vals(8)) ! yz part: (32 + 23)/2
        res%vals(6) = 0.5D0 * (a%vals(3) + a%vals(7)) ! xz part: (31 + 13)/2
        
    end function ddot_I4O4T_3D2O

    pure function ddot_3D2O_I4O4T(a, I4) result(res)
        !! Computes the double contraction product: \(\mathbf{res} = \mathbf{A} : \mathbb{I}^S\).
        !!
        !! Equivalent to \(\mathbb{I}^S : \mathbf{A}\) due to the symmetry of the operator.
        implicit none
        type(ten_3D2O),  intent(in) :: a
        type(iden_4O4T), intent(in) :: I4
        type(ten_3D2Osym)           :: res
        
        ! Delegate to the primary implementation
        res = ddot_I4O4T_3D2O(I4, a)
        
    end function ddot_3D2O_I4O4T

end module mod_operator_I4O4T_3D2O