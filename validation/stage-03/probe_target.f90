program probe_target
    use iso_fortran_env, only: real64
    use ieee_arithmetic, only: ieee_is_finite
    use omp_lib, only: omp_get_num_devices,omp_is_initial_device
    use probe_core
    implicit none
    integer, parameter :: n=257
    type(material_parameters) :: parameters
    real(real64) :: input(6,n),reference(9,n),actual(9,n),error,checksum,expected_checksum
    integer :: statuses(n),expected(n),i,j,devices
    logical :: initial_device,require_device
    character(len=32) :: mode
    call get_command_argument(1,mode)
    require_device=trim(mode)=='--require-device'
    devices=omp_get_num_devices()
    print '(a,i0)', 'omp_num_devices=',devices
    initial_device=.true.
    !$omp target map(from:initial_device)
        initial_device=omp_is_initial_device()
    !$omp end target
    print '(a,l1)', 'executed_on_initial_device=',initial_device
    if (require_device .and. initial_device) error stop 'GPU required; host fallback rejected.'
    if (require_device .and. devices==0) error stop 'GPU required; no OpenMP devices.'
    do i=1,n
        do j=1,6
            input(j,i)=real(mod(i*j,19)-9,real64)*1d-4
        end do
        call point_kernel(input(:,i),parameters,i,reference(1:6,i),reference(7,i), &
                          reference(8,i),reference(9,i),expected(i))
    end do
    actual=0d0
    statuses=-1
    checksum=0d0
    !$omp target data map(to:input,parameters) map(alloc:actual,statuses)
        !$omp target teams distribute parallel do private(i)
        do i=1,n
            call point_kernel(input(:,i),parameters,i,actual(1:6,i),actual(7,i), &
                              actual(8,i),actual(9,i),statuses(i))
        end do
        !$omp end target teams distribute parallel do
        ! Reuse resident results in a second region before copying them back.
        !$omp target teams distribute parallel do reduction(+:checksum) map(tofrom:checksum)
        do i=1,n
            checksum=checksum+actual(7,i)+actual(8,i)
        end do
        !$omp end target teams distribute parallel do
        !$omp target update from(actual,statuses)
    !$omp end target data
    if (.not.all(ieee_is_finite(actual))) error stop 'Non-finite device result.'
    if (.not.ieee_is_finite(checksum)) error stop 'Non-finite reduction.'
    if (any(statuses/=expected)) error stop 'Device status mismatch.'
    if (count(statuses==1)/=n/17) error stop 'Singular matrices not detected.'
    if (maxval(actual(9,:))>1d-12) error stop 'Device inverse residual too large.'
    error=maxval(abs(actual-reference)/max(1d0,abs(reference)))
    if (error>1d-11) error stop 'Device/host mismatch.'
    expected_checksum=sum(reference(7,:)+reference(8,:))
    if (abs(checksum-expected_checksum)>1d-11*max(1d0,abs(expected_checksum))) &
        error stop 'Reduction mismatch.'
    print '(a,i0)', 'samples=',n
    print '(a,es24.16)', 'normalized_error=',error
    if (initial_device) then
        print '(a)', 'PASS_HOST_TARGET_ONLY'
    else
        print '(a)', 'PASS_DEVICE_TARGET'
    end if
end program
