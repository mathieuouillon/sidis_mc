!c------------------------------------------------------------------------------
!c Unit test for calculate_vz_position() subroutine
!c Verifies vertex z-position ranges for different targets
!c------------------------------------------------------------------------------
program test_vz
    implicit none

    real :: vz
    integer :: i, failures, total
    real :: vz_min, vz_max
    logical :: in_range

    failures = 0
    total = 0

    ! Seed intrinsic RNG for reproducibility
    call srand(12345)

    ! =========================================================================
    ! Test 1: Sn target (type 10)
    ! target_pos = -3.5, target_length = 0.018
    ! vz should be in [-3.5 - 0.009, -3.5 + 0.009] = [-3.509, -3.491]
    ! =========================================================================
    total = total + 1
    in_range = .true.
    do i = 1, 1000
        call calculate_vz_position(10, vz)
        if (vz < -3.510 .or. vz > -3.490) then
            in_range = .false.
            write(*,*) 'FAIL Test 1: Sn vz out of range:', vz
            exit
        end if
    end do
    if (in_range) then
        write(*,*) 'PASS Test 1: Sn target vz in expected range'
    else
        failures = failures + 1
    end if

    ! =========================================================================
    ! Test 2: Cu target (type 15)
    ! target_pos = -8.5, target_length = 0.009
    ! vz should be in [-8.5 - 0.0045, -8.5 + 0.0045] = [-8.5045, -8.4955]
    ! =========================================================================
    total = total + 1
    in_range = .true.
    do i = 1, 1000
        call calculate_vz_position(15, vz)
        if (vz < -8.506 .or. vz > -8.494) then
            in_range = .false.
            write(*,*) 'FAIL Test 2: Cu vz out of range:', vz
            exit
        end if
    end do
    if (in_range) then
        write(*,*) 'PASS Test 2: Cu target vz in expected range'
    else
        failures = failures + 1
    end if

    ! =========================================================================
    ! Test 3: Carbon target (type 7)
    ! Randomly chooses target_pos = -3.5 or -8.5, target_length = 0.2
    ! vz should be near -3.5 or -8.5, within ±0.1
    ! =========================================================================
    total = total + 1
    in_range = .true.
    do i = 1, 1000
        call calculate_vz_position(7, vz)
        if (.not. ((vz >= -3.61 .and. vz <= -3.39) .or. &
                   (vz >= -8.61 .and. vz <= -8.39))) then
            in_range = .false.
            write(*,*) 'FAIL Test 3: Carbon vz out of range:', vz
            exit
        end if
    end do
    if (in_range) then
        write(*,*) 'PASS Test 3: Carbon target vz in expected range'
    else
        failures = failures + 1
    end if

    ! =========================================================================
    ! Test 4: LD2 target (type 1)
    ! target_pos = -5.0, target_length = 5.0
    ! vz should be in [-7.5, -2.5]
    ! =========================================================================
    total = total + 1
    in_range = .true.
    do i = 1, 1000
        call calculate_vz_position(1, vz)
        if (vz < -7.51 .or. vz > -2.49) then
            in_range = .false.
            write(*,*) 'FAIL Test 4: LD2 vz out of range:', vz
            exit
        end if
    end do
    if (in_range) then
        write(*,*) 'PASS Test 4: LD2 target vz in expected range'
    else
        failures = failures + 1
    end if

    ! =========================================================================
    ! Test 5: Default target (type 0, proton) should return vz = 0.0
    ! =========================================================================
    total = total + 1
    call calculate_vz_position(0, vz)
    if (vz /= 0.0) then
        write(*,*) 'FAIL Test 5: Proton vz should be 0.0, got:', vz
        failures = failures + 1
    else
        write(*,*) 'PASS Test 5: Proton target returns vz = 0.0'
    end if

    ! =========================================================================
    ! Summary
    ! =========================================================================
    write(*,*) ''
    write(*,*) '================================'
    write(*,*) 'Total tests: ', total
    write(*,*) 'Passed:      ', total - failures
    write(*,*) 'Failed:      ', failures
    write(*,*) '================================'

    if (failures > 0) then
        stop 1
    end if

end program test_vz
