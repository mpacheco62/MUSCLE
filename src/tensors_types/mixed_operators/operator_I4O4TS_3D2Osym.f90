module mod_operator_I4O4TS_3D2Osym
    !! Module mod_operator_I4O4TS_3D2Osym
    !! =====================================
    !!
    !! Defines mixed algebraic operations between the scaled 3D fourth-order 
    !! symmetric identity tensor (`iden_4O4TS`, \(c\mathbb{I}^S\)) and 
    !! 3D symmetric second-order tensors (`ten_3D2Osym`).
    !!
    !! Since the second-order tensor is already symmetric, the double 
    !! contraction with the symmetric identity results in a simple 
    !! scalar-tensor multiplication:
    !! \[ \mathbf{res} = c\mathbb{I}^S : \mathbf{A} = c\mathbf{A} \]

    use, intrinsic :: iso_fortran_env
    use mod_iden_4O4TS
    use mod_ten_3D2Osym
    implicit none
    private
    
    public :: operator(.ddot.)
    interface operator (.ddot.)
        module procedure ddot_I4O4TS_3D2Osym
        module procedure ddot_3D2Osym_I4O4TS
    end interface

contains

    pure function ddot_I4O4TS_3D2Osym(I4S, a) result(res)
        !! Computes the scaled double contraction product: \(\mathbf{res} = c\mathbb{I}^S : \mathbf{A}\).
        !!
        !! For a symmetric tensor \(\mathbf{A}\), this is equivalent to \(c \cdot \mathbf{A}\).
        implicit none
        type(iden_4O4TS),  intent(in) :: I4S
            !! Scaled fourth-order symmetric identity tensor (\(c\mathbb{I}^S\)).
        type(ten_3D2Osym), intent(in) :: a
            !! Symmetric second-order tensor (6 components, Voigt).
        type(ten_3D2Osym)             :: res
            !! Resulting symmetric second-order tensor.

        res%vals = I4S%val * a%vals
    end function ddot_I4O4TS_3D2Osym

    pure function ddot_3D2Osym_I4O4TS(a, I4S) result(res)
        !! Computes the scaled double contraction product: \(\mathbf{res} = \mathbf{A} : c\mathbb{I}^S\).
        !!
        !! Equivalent to \(c\mathbb{I}^S : \mathbf{A}\).
        implicit none
        type(ten_3D2Osym), intent(in) :: a
        type(iden_4O4TS),  intent(in) :: I4S
        type(ten_3D2Osym)             :: res
        
        ! Delegate to the primary implementation
        res = ddot_I4O4TS_3D2Osym(I4S, a)
    end function ddot_3D2Osym_I4O4TS

end module mod_operator_I4O4TS_3D2Osym