! module mod_operator_3D2Osym_3D2O
!     use, intrinsic :: iso_fortran_env
!     use mod_ten_3D2Osym
!     use mod_ten_3D2O
!     implicit none
!     private
    
!     public :: operator(+)
!     interface operator (+)
!         module procedure sum_3D2O_3D2Osym
!         module procedure sum_3D2Osym_3D2O
!     end interface

!     public :: operator(-)
!     interface operator (-)
!         module procedure sub_3D2O_3D2Osym
!         module procedure sub_3D2Osym_3D2O
!     end interface

!     public :: operator(.ddot.)
!     interface operator (.ddot.)
!         module procedure ddot_3D2O_3D2Osym
!         module procedure ddot_3D2Osym_3D2O
!     end interface

!     public :: assignment (=)
!     interface assignment (=)
!         module procedure assign_3D2O_3D2Osym
!     end interface

!     contains

! end module mod_operator_3D2Osym_3D2O