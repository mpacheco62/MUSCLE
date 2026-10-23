! Standalone compiler probes, not a replacement constitutive implementation.
! invert6 includes the unchanged GPL-3.0-or-later MUSCLE LU kernel.
module probe_core
    use iso_fortran_env, only: real64
    implicit none
    private
    public :: material_parameters, point_kernel, operator(.contract.)
    type :: material_parameters
        real(real64) :: lambda=115384.61538461538d0, mu=76923.07692307692d0
        real(real64) :: k=100d0, e0=1d0, exponent=0.1d0
    end type
    type :: tensor6
        real(real64) :: vals(6)
    end type
    interface operator(.contract.)
        module procedure contract_tensor
    end interface
    !$omp declare target(point_kernel, invert6, contract_tensor)
contains
    pure function contract_tensor(a,b) result(value)
        type(tensor6), intent(in) :: a,b
        real(real64) :: value
        value=sum(a%vals(1:3)*b%vals(1:3))+2d0*sum(a%vals(4:6)*b%vals(4:6))
    end function

    pure subroutine invert6(A,AINV,OK_FLAG)
        integer, parameter :: N=6
        real(real64), intent(in) :: A(N,N)
        real(real64), intent(out) :: AINV(N,N)
        logical, intent(out) :: OK_FLAG
        include 'muscle_math_lu_inverse.inc'
    end subroutine

    pure subroutine point_kernel(input,parameters,index,stress,vm,hard,inverse_residual,status)
        real(real64), intent(in) :: input(6)
        type(material_parameters), intent(in) :: parameters
        integer, intent(in) :: index
        real(real64), intent(out) :: stress(6),vm,hard,inverse_residual
        integer, intent(out) :: status
        type(tensor6) :: deviator
        real(real64) :: trace,A(6,6),inverse(6,6),identity_error(6,6)
        integer :: row
        logical :: ok
        trace=sum(input(1:3))
        stress=2d0*parameters%mu*input
        stress(1:3)=stress(1:3)+parameters%lambda*trace
        deviator%vals=stress
        deviator%vals(1:3)=stress(1:3)-sum(stress(1:3))/3d0
        vm=sqrt(max(0d0,1.5d0*(deviator .contract. deviator)))
        hard=parameters%k*(parameters%e0+abs(input(1)))**parameters%exponent
        A=0.02d0
        do row=1,6
            A(row,row)=real(row+1,real64)+abs(input(row))
        end do
        ! Use an exactly zero row: the upstream kernel detects exact zero
        ! pivots, not arbitrary near-singularity after floating-point rounding.
        if (mod(index,17)==0) A(6,:)=0d0
        call invert6(A,inverse,ok)
        if (ok) then
            identity_error=matmul(A,inverse)
            do row=1,6
                identity_error(row,row)=identity_error(row,row)-1d0
            end do
            inverse_residual=maxval(abs(identity_error))
            status=0
        else
            inverse_residual=maxval(abs(inverse))
            status=1
        end if
    end subroutine
end module
