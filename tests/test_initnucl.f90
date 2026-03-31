!c------------------------------------------------------------------------------
!c Unit test for InitNucl() subroutine
!c Verifies iZ, iA, rFM for each target type (0-15)
!c------------------------------------------------------------------------------
program test_initnucl
    use fermi_motion_module
    use config_module
    ! Need to use all modules that InitNucl uses
    use kinematics_module
    use file_names_module
    use event_info_module
    use particles_module
    use density_module
    use interaction_module
    use quenching_module
    use misc_module
    use pythia_commons
    implicit none

    integer :: failures, total
    integer :: expected_iZ, expected_iA
    real :: expected_rFM

    failures = 0
    total = 0

    ! Target 0: proton
    call check_target(0, 1, 1, 0.0, failures, total)
    ! Target 1: deuterium
    call check_target(1, 1, 2, 0.07, failures, total)
    ! Target 2: tritium
    call check_target(2, 1, 3, 0.1, failures, total)
    ! Target 3: 3He
    call check_target(3, 2, 3, 0.1, failures, total)
    ! Target 4: 4He
    call check_target(4, 2, 4, 0.120, failures, total)
    ! Target 5: 6Li
    call check_target(5, 3, 6, 0.17, failures, total)
    ! Target 6: 7Li
    call check_target(6, 3, 7, 0.17, failures, total)
    ! Target 7: Carbon
    call check_target(7, 6, 12, 0.221, failures, total)
    ! Target 8: Aluminum
    call check_target(8, 13, 27, 0.235, failures, total)
    ! Target 9: Iron
    call check_target(9, 26, 56, 0.260, failures, total)
    ! Target 10: Tin
    call check_target(10, 50, 120, 0.260, failures, total)
    ! Target 11: Lead
    call check_target(11, 82, 208, 0.265, failures, total)
    ! Target 12: Neon
    call check_target(12, 10, 20, 0.228, failures, total)
    ! Target 13: Krypton
    call check_target(13, 36, 84, 0.260, failures, total)
    ! Target 14: Xenon
    call check_target(14, 54, 132, 0.260, failures, total)
    ! Target 15: Cu
    call check_target(15, 29, 63, 0.260, failures, total)

    write(*,*) ''
    write(*,*) '================================'
    write(*,*) 'Total tests: ', total
    write(*,*) 'Passed:      ', total - failures
    write(*,*) 'Failed:      ', failures
    write(*,*) '================================'

    if (failures > 0) then
        stop 1
    end if

contains

    subroutine check_target(tg, exp_Z, exp_A, exp_rFM, failures, total)
        integer, intent(in) :: tg, exp_Z, exp_A
        real, intent(in) :: exp_rFM
        integer, intent(inout) :: failures, total

        total = total + 1
        iTg = tg
        call InitNucl()

        if (iZ /= exp_Z .or. iA /= exp_A .or. rFM /= exp_rFM) then
            write(*,*) 'FAIL Target ', tg, ':'
            if (iZ /= exp_Z) write(*,*) '  iZ: expected', exp_Z, ' got', iZ
            if (iA /= exp_A) write(*,*) '  iA: expected', exp_A, ' got', iA
            if (rFM /= exp_rFM) write(*,*) '  rFM: expected', exp_rFM, ' got', rFM
            failures = failures + 1
        else
            write(*,*) 'PASS Target ', tg, ': iZ=', iZ, ' iA=', iA, ' rFM=', rFM
        end if
    end subroutine check_target

end program test_initnucl
