! SPDX-License-Identifier: GPL-3.0-or-later
! Test the generated library kernels against the original CPU API and singular inverses.
program test_generated_kernels
    use iso_fortran_env, only: real64
    use ieee_arithmetic, only: ieee_is_finite
    use omp_lib, only: omp_is_initial_device, omp_get_num_devices
    use muscle_hard_swift, only: Swift_hardening
    use muscle_math_inverses, only: M66INV
    use muscle_kernels_cpu, only: cpu_parameters => swift_parameters, &
        cpu_stress => swift_stress, cpu_dstress => swift_dstress, &
        cpu_ddstress => swift_ddstress, cpu_inverse => invert6
    use muscle_kernels_target, only: target_parameters => swift_parameters, &
        target_stress => swift_stress, target_dstress => swift_dstress, &
        target_ddstress => swift_ddstress, target_inverse => invert6
    implicit none
    integer, parameter :: npoints=257
    type(Swift_hardening) :: original
    type(cpu_parameters) :: cpu
    type(target_parameters) :: target
    real(real64) :: matrices(6,6,npoints), strains(npoints), reference(3,npoints), actual(3,npoints)
    real(real64) :: original_inverse(6,6,npoints), inverses(6,6,npoints), local_inverse(6,6)
    real(real64) :: local_values(3), max_error
    logical :: original_ok(npoints), actual_ok(npoints), local_ok, initial(npoints), require_device
    integer :: point, row
    character(len=32) :: argument

    call get_command_argument(1, argument)
    if (trim(argument)/='--host' .and. trim(argument)/='--require-device') error stop 'Expected execution policy'
    require_device=trim(argument)=='--require-device'
    original=Swift_hardening(k=100d0, e0=1d0, n=0.1d0)
    cpu=cpu_parameters(k=original%k, e0=original%e0, n=original%n)
    target=target_parameters(k=original%k, e0=original%e0, n=original%n)
    do point=1,npoints
        strains(point)=real(point-129,real64)*1d-3
        reference(:,point)=[original%stress(strains(point)), original%dstress_dep(strains(point)), &
                            original%ddstress_ddep(strains(point))]
        matrices(:,:,point)=0.02d0
        do row=1,6
            matrices(row,row,point)=real(row+1,real64)+abs(strains(point))
        end do
        if (mod(point,17)==0) matrices(6,:,point)=0d0
        call M66INV(matrices(:,:,point),original_inverse(:,:,point),original_ok(point))
        local_values=[cpu_stress(cpu,strains(point)),cpu_dstress(cpu,strains(point)), &
                      cpu_ddstress(cpu,strains(point))]
        call check_values(local_values,reference(:,point),'CPU Swift')
        call cpu_inverse(matrices(:,:,point),local_inverse,local_ok)
        if (local_ok .neqv. original_ok(point)) error stop 'CPU inversion status'
        if (maxval(abs(local_inverse-original_inverse(:,:,point)))>1d-12) error stop 'CPU inverse'
    end do
    actual=0d0
    actual_ok=.false.
    inverses=0d0
    initial=.true.
    print '(a,i0)', 'omp_num_devices=',omp_get_num_devices()
    !$omp target teams distribute parallel do map(to:target,strains,matrices) &
    !$omp& map(from:actual,inverses,actual_ok,initial)
    do point=1,npoints
        initial(point)=omp_is_initial_device()
        actual(:,point)=[target_stress(target,strains(point)),target_dstress(target,strains(point)), &
                         target_ddstress(target,strains(point))]
        call target_inverse(matrices(:,:,point),inverses(:,:,point),actual_ok(point))
    end do
    !$omp end target teams distribute parallel do
    print '(a,i0)', 'calculation_host_points=',count(initial)
    if (require_device .and. (any(initial) .or. omp_get_num_devices()==0)) &
        error stop 'OFFLOAD_ERROR: real device execution is required'
    if (.not.require_device .and. .not.all(initial)) error stop 'Expected explicit host execution'
    max_error=0d0
    do point=1,npoints
        call check_values(actual(:,point),reference(:,point),'Target Swift')
        if (actual_ok(point) .neqv. original_ok(point)) error stop 'Target inversion status'
        if (.not.all(ieee_is_finite(inverses(:,:,point)))) error stop 'Nonfinite target inverse'
        if (maxval(abs(inverses(:,:,point)-original_inverse(:,:,point)))>1d-12) error stop 'Target inverse'
        max_error=max(max_error,maxval(abs(actual(:,point)-reference(:,point))))
    end do
    print '(a,i0)', 'samples=',npoints
    print '(a,es24.16)', 'swift_max_absolute_difference=',max_error
    if (require_device) then
        print '(a)', 'PASS_GENERATED_KERNELS_DEVICE'
    else
        print '(a)', 'PASS_GENERATED_KERNELS_HOST'
    end if
contains
    subroutine check_values(values,expected,label)
        real(real64), intent(in) :: values(:),expected(:)
        character(len=*), intent(in) :: label
        if (.not.all(ieee_is_finite(values))) then
            print *, label
            error stop 'Nonfinite Swift output'
        end if
        if (any(abs(values-expected)>1d-12+1d-11*abs(expected))) then
            print *, label, values, expected
            error stop 'Swift difference exceeds fixed tolerance'
        end if
    end subroutine
end program
