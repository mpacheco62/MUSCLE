program probe_host
    use iso_fortran_env, only: real64
    use ieee_arithmetic, only: ieee_is_finite
    use probe_core
    !$ use omp_lib, only: omp_get_num_threads
    implicit none
    integer, parameter :: n=257
    type(material_parameters) :: parameters
    real(real64) :: input(6,n),reference(9,n),actual(9,n),error
    integer :: statuses(n),expected(n),i,j,threads
    character(len=32) :: mode
    call get_command_argument(1,mode)
    do i=1,n
        do j=1,6
            input(j,i)=real(mod(i*j,19)-9,real64)*1d-4
        end do
        call point_kernel(input(:,i),parameters,i,reference(1:6,i),reference(7,i), &
                          reference(8,i),reference(9,i),expected(i))
    end do
    threads=1
    !$omp parallel shared(threads)
    !$omp single
    !$ threads=omp_get_num_threads()
    !$omp end single
    !$omp end parallel
    !$omp parallel do private(i)
    do i=1,n
        call point_kernel(input(:,i),parameters,i,actual(1:6,i),actual(7,i), &
                          actual(8,i),actual(9,i),statuses(i))
    end do
    !$omp end parallel do
    if (.not.all(ieee_is_finite(actual))) error stop 'Non-finite result.'
    if (any(statuses/=expected)) error stop 'Status mismatch.'
    if (count(statuses==1)/=n/17) error stop 'Singular matrices not detected.'
    if (maxval(actual(9,:))>1d-12) error stop 'Inverse residual too large.'
    error=maxval(abs(actual-reference)/max(1d0,abs(reference)))
    if (error>1d-12) error stop 'Serial/OpenMP mismatch.'
    if (trim(mode)=='--require-parallel' .and. threads<2) error stop 'No parallel team.'
    print '(a,i0)', 'threads=',threads
    print '(a,i0)', 'samples=',n
    print '(a,es24.16)', 'normalized_error=',error
    print '(a)', 'PASS_HOST'
end program
