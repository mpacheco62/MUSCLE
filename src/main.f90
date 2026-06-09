Program test
    use, intrinsic :: iso_fortran_env
    implicit None
    ! use muscle_yield_vonmises
    ! use basic_operations
    integer :: i
    real*8 :: eps22Init_t = 0.001216D0, eps22Init_c=0.002731D0, E22 = 4.76D0
    real*8 :: El=0.01D0, Gm=2.00000D-6, fmt, wmt, wmc
    
    do i=0,900
      fmt = 1D0 + i/90D0
      wmt  = 1D0 - 1D0/fmt * exp(-E22*eps22Init_t**2 * (fmt - 1D0)*El/Gm)
      wmc  = 1D0 - 1D0/fmt * exp(-E22*eps22Init_c**2 * (fmt - 1D0)*El/Gm)
      write(900, *) fmt, wmt, wmc
    end do

end program