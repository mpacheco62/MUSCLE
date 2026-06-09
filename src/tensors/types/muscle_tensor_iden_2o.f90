module muscle_tensor_iden_2o
    !! author: MPacheco
    !! version: 1.0 - Initial documentation
    !!
    !! Module muscle_tensor_iden_2o
    !! =====================
    !!
    !! This module defines the standard (non-scaled) 3D second-order identity \(\delta_{ij}\).
    !!
    !! This module defines the standard (non-scaled) 3D second-order identity tensor type,
    !! denoted mathematically as \(\delta_{ij}\). It represents the identity matrix.
    !! The module also provides overloaded operators for arithmetic operations involving
    !! this identity tensor, often interacting with its scaled counterpart (`iden_2OS`)
    !! defined in `muscle_tensor_iden_2os`.
    !!
    !! The `iden_2O` type itself does not store any numerical value; it acts as a
    !! symbolic representation of the identity tensor where the implicit scaling factor is 1.0.
    !! Operations involving `iden_2O` typically result in a scaled identity tensor
    !! (`iden_2OS`), which explicitly stores the scaling factor.
    !!
    !! Public Entities
    !! ---------------
    !!
    !! ### Derived Type:
    !!
    !! - `iden_2O`: Represents the standard 3D second-order identity tensor \(\delta_{ij}\).
    !!                It has no components and implicitly represents the identity with a
    !!                scaling factor of 1.0.
    !!
    !! ### Operators:
    !!
    !! The following operators are overloaded. Note that most operations involving `iden_2O`
    !! return a scaled identity tensor `iden_2OS`.
    !!
    !! - `+`: Addition:
    !!     - `iden_2O + iden_2O`: Results in `iden_2OS` with value `2.0`.
    !!     - `iden_2O + iden_2OS`: Results in `iden_2OS` (value `1.0 + iden_2OS%val`).
    !!     - `iden_2OS + iden_2O`: Results in `iden_2OS` (value `iden_2OS%val + 1.0`).
    !! - `-`: Subtraction:
    !!     - `- iden_2O`: Unary negation. Results in `iden_2OS` with value `-1.0`.
    !!     - `iden_2O - iden_2O`: Results in `iden_2OS` with value `0.0`.
    !!     - `iden_2O - iden_2OS`: Results in `iden_2OS` (value `1.0 - iden_2OS%val`).
    !!     - `iden_2OS - iden_2O`: Results in `iden_2OS` (value `iden_2OS%val - 1.0`).
    !! - `*`: Multiplication:
    !!     - `iden_2O * real(real64)`: Scalar multiplication. Results in `iden_2OS`.
    !!     - `real(real64) * iden_2O`: Scalar multiplication. Results in `iden_2OS`.
    !!     - `iden_2O * iden_2O`: Identity tensor multiplication (conceptually \(I \cdot I = I\)). Results in `iden_2O`.
    !! - `/`: Division:
    !!     - `iden_2O / real(real64)`: Scalar division. Results in `iden_2OS` (value `1.0 / scalar`).
    !!
    !! Usage
    !! -----
    !!
    !! ```fortran
    !! program example_iden_2O_usage
    !!   use muscle_tensor_iden_2o
    !!   use muscle_tensor_iden_2os
    !!   use iso_fortran_env, only: real64
    !!   implicit none
    !!
    !!   type(iden_2O) :: id_tensor       ! Represents delta_ij
    !!   type(iden_2OS) :: scaled_id1, scaled_id2
    !!   real(real64) :: factor = 5.0D0
    !!
    !!   ! Create scaled identities from the standard identity
    !!   scaled_id1 = id_tensor * factor      ! scaled_id1%val = 5.0
    !!   scaled_id2 = id_tensor / 2.0D0       ! scaled_id2%val = 0.5
    !!
    !!   ! Operations involving iden_2O
    !!   scaled_id1 = id_tensor + id_tensor   ! scaled_id1%val = 2.0
    !!   scaled_id2 = -id_tensor              ! scaled_id2%val = -1.0
    !!   scaled_id1 = scaled_id1 + id_tensor  ! scaled_id1%val = 2.0 + 1.0 = 3.0
    !!
    !!   print *, "Scaled ID 1 value:", scaled_id1%val ! Expected: 3.0
    !!   print *, "Scaled ID 2 value:", scaled_id2%val ! Expected: -1.0
    !!
    !! end program example_iden_2O_usage
    !! ```
    !! For more information see [[muscle_tensors]]
    use, intrinsic :: iso_fortran_env
    use muscle_tensor_iden_2os
    implicit none
    private

    type, public :: iden_2O
        !! 3D second-order identity tensor \(\delta_{ij}\)
        !! ===============================================
        !! 
        !! Type for the standard (non-scaled) 3D second-order identity tensor \(\delta_{ij}\).
        !!
        !! Represents the standard 3D second-order identity tensor, mathematically the Kronecker delta \(\delta_{ij}\).
        !! This type acts as a symbolic marker and contains **no data components**. The implicit
        !! scaling factor is always considered to be `1.0D0`. It is primarily used in operations
        !! (like addition, subtraction, or scalar multiplication/division) that typically result
        !! in a scaled identity tensor (`iden_2OS`), which explicitly stores a scaling factor.
        !! Multiplying two `iden_2O` types results in another `iden_2O` type (\(I \cdot I = I\)).
        !! 
        !! For more information see [[muscle_tensors]]
    end type iden_2O


    public :: operator(+)
    interface operator (+)
        module procedure sum_I2O_I2O
        module procedure sum_I2O_I2OS
        module procedure sum_I2OS_I2O
    end interface

    public :: operator(-)
    interface operator (-)
        module procedure subU_I2O
        module procedure sub_I2O_I2O
        module procedure sub_I2O_I2OS
        module procedure sub_I2OS_I2O
    end interface

    public :: operator(*)
    interface operator (*)
        module procedure mul_I2O_real64
        module procedure mul_real64_I2O
        module procedure mul_I2O_I2O
    end interface

    public :: operator( / )
    interface operator ( / )
        module procedure div_I2O_real64
    end interface

    public :: write(formatted)
    interface write(formatted)
        module procedure print_ten_I2O
    end interface

contains

    pure function sum_I2O_I2O(I2a, I2b) result(res)
        !! `+` Adds two standard identity tensors. Result is a scaled identity with value 2.0.
        implicit none
        type(iden_2O), intent(in) :: I2a, I2b
        type(iden_2OS) :: res
        res%val = 2D0
    end function sum_I2O_I2O

    pure function sum_I2O_I2OS(I2a, I2b) result(res)
        !! `+` Adds a standard identity tensor and a scaled identity tensor.
        implicit none
        type(iden_2O), intent(in) :: I2a
        type(iden_2OS), intent(in) :: I2b
        type(iden_2OS) :: res
        res%val = 1D0 + I2b%val
    end function sum_I2O_I2OS

    pure function sum_I2OS_I2O(I2a, I2b) result(res)
        !! `+` Adds a scaled identity tensor and a standard identity tensor.
        implicit none
        type(iden_2OS), intent(in) :: I2a
        type(iden_2O), intent(in) :: I2b
        type(iden_2OS) :: res
        res%val = I2a%val + 1D0
    end function sum_I2OS_I2O

    pure function subU_I2O(I2a) result(res)
        !! `-` Unary negation of the standard identity tensor. Result is a scaled identity with value -1.0.
        implicit none
        type(iden_2O), intent(in) :: I2a
        type(iden_2OS) :: res
        res%val = -1D0
    end function subU_I2O

    pure function sub_I2O_I2O(I2a, I2b) result(res)
        !! `-` Subtracts one standard identity tensor from another. Result is a scaled identity with value 0.0.
        implicit none
        type(iden_2O), intent(in) :: I2a, I2b
        type(iden_2OS) :: res
        res%val = 0D0
    end function sub_I2O_I2O

    pure function sub_I2O_I2OS(I2a, I2b) result(res)
        !! `-` Subtracts a scaled identity tensor from a standard identity tensor.
        implicit none
        type(iden_2O), intent(in) :: I2a
        type(iden_2OS), intent(in) :: I2b
        type(iden_2OS) :: res
        res%val = 1D0 - I2b%val
    end function sub_I2O_I2OS

    pure function sub_I2OS_I2O(I2a, I2b) result(res)
        !! `-` Subtracts a standard identity tensor from a scaled identity tensor.
        implicit none
        type(iden_2OS), intent(in) :: I2a
        type(iden_2O), intent(in) :: I2b
        type(iden_2OS) :: res
        res%val = I2a%val - 1D0
    end function sub_I2OS_I2O

    pure function mul_I2O_real64(I2, a) result(res)
        !! `*` Multiplies the standard identity tensor by a scalar. Result is a scaled identity.
        implicit none
        type(iden_2O), intent(in) :: I2
        real(real64), intent(in) :: a
        type(iden_2OS) :: res
        res%val = a 
    end function mul_I2O_real64

    pure function mul_real64_I2O(a, I2) result(res)
        !! `*` Multiplies a scalar by the standard identity tensor. Result is a scaled identity.
        implicit none
        type(iden_2O), intent(in) :: I2
        real(real64), intent(in) :: a
        type(iden_2OS) :: res
        res%val = a 
    end function mul_real64_I2O

    pure function mul_I2O_I2O(I2a, I2b) result(res)
        !! `\cdot` Multiplies two standard identity tensors. Conceptually I*I = I. Result is the standard identity
        implicit none
        type(iden_2O), intent(in) :: I2a, I2b
        type(iden_2O) :: res
    end function mul_I2O_I2O

    pure function div_I2O_real64(I2, a) result(res)
        !! Divides the standard identity tensor by a scalar. Result is a scaled identity.
        implicit none
        type(iden_2O), intent(in) :: I2
        real(real64), intent(in) :: a
        type(iden_2OS) :: res
        res%val = 1D0/a
    end function div_I2O_real64


    subroutine print_ten_I2O(dtv, unit, iotype, v_list, iostat, iomsg)
        !! Custom I/O formatting for iden_2O.
        class(iden_2O), intent(in)      :: dtv
        integer, intent(in)             :: unit
        character(len=*), intent(in)    :: iotype
        integer, intent(in)             :: v_list(:)
        integer, intent(out)            :: iostat
        character(len=*), intent(inout) :: iomsg
        
        character(len=100) :: fmt_string
        integer :: w, d

        w = 10 
        d = 4  
        iostat = 0
        
        if (size(v_list) >= 1) w = v_list(1)
        if (size(v_list) >= 2) d = v_list(2)

        if (iotype == "DTLIST" .or. iotype == "LIST") then
            ! Modo lista
            write(fmt_string, "('(A, 1(F', I0, '.', I0, ', 1X), A)')") w, d
            write(unit, fmt_string, iostat=iostat) &
                "[", 1D0, "]"
        else
            ! MODO MATRIZ: Un solo formato con saltos de línea incorporados
            ! Usamos el descriptor '/' para obligar el salto de línea entre filas.
            ! El primer '/' asegura que empiece en una línea nueva.
            write(fmt_string, "('(/, 3(F', I0, '.', I0, ', 2X), /, 3(F', I0, '.', I0, ', 2X), /, 3(F', I0, '.', I0, ', 2X))')") &
                w, d, w, d, w, d
            
            write(unit, fmt_string, iostat=iostat) &
                1D0, 0D0, 0D0, & ! Fila 1
                0D0, 1D0, 0D0, & ! Fila 2
                0D0, 0D0, 1D0    ! Fila 3
        end if

        if (iostat /= 0) then
            iomsg = "Error in print_ten_3D2Osym: Failed to write to the specified unit."
        end if
    end subroutine print_ten_I2O
end module muscle_tensor_iden_2o