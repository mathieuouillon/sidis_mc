!c------------------------------------------------------------------------------
!c Unit test for the Lorentz transformation subroutine TL()
!c Tests: analytical boost, roundtrip, invariant mass preservation
!c------------------------------------------------------------------------------
program test_tl
    implicit none

    real :: E, PP, Px, Py, Pz
    real :: BB, Bx, By, Bz
    real :: E_save, Px_save, Py_save, Pz_save
    real :: mass_before, mass_after
    real :: gamma, expected_E, expected_Pz
    real :: tol
    integer :: failures, total

    tol = 1.0e-5
    failures = 0
    total = 0

    ! =========================================================================
    ! Test 1: Boost along z with beta = 0.5 on a particle at rest (mass=1 GeV)
    ! Analytical: E' = gamma*E, Pz' = -gamma*beta*E, Px'=Py'=0
    ! =========================================================================
    total = total + 1
    E = 1.0; PP = 0.0; Px = 0.0; Py = 0.0; Pz = 0.0
    BB = 0.5; Bx = 0.0; By = 0.0; Bz = 0.5
    gamma = 1.0 / sqrt(1.0 - BB**2)
    expected_E = gamma * E
    expected_Pz = -gamma * BB * E

    call TL(E, PP, Px, Py, Pz, BB, Bx, By, Bz)

    if (abs(E - expected_E) / expected_E > tol) then
        write(*,*) 'FAIL Test 1a: E after z-boost'
        write(*,*) '  Expected:', expected_E, ' Got:', E
        failures = failures + 1
    else if (abs(Px) > tol) then
        write(*,*) 'FAIL Test 1b: Px should be 0'
        write(*,*) '  Got:', Px
        failures = failures + 1
    else if (abs(Py) > tol) then
        write(*,*) 'FAIL Test 1c: Py should be 0'
        write(*,*) '  Got:', Py
        failures = failures + 1
    else if (abs(Pz - expected_Pz) / abs(expected_Pz) > tol) then
        write(*,*) 'FAIL Test 1d: Pz after z-boost'
        write(*,*) '  Expected:', expected_Pz, ' Got:', Pz
        failures = failures + 1
    else
        write(*,*) 'PASS Test 1: Boost along z on particle at rest'
    end if

    ! =========================================================================
    ! Test 2: Invariant mass preservation for a moving particle
    ! m^2 = E^2 - p^2 should be unchanged after boost
    ! BB must equal sqrt(Bx^2 + By^2 + Bz^2) exactly
    ! =========================================================================
    total = total + 1
    E = 2.0; Px = 0.3; Py = -0.5; Pz = 1.2
    PP = sqrt(Px**2 + Py**2 + Pz**2)
    mass_before = E**2 - Px**2 - Py**2 - Pz**2

    Bx = 0.2; By = -0.3; Bz = 0.4
    BB = sqrt(Bx**2 + By**2 + Bz**2)

    call TL(E, PP, Px, Py, Pz, BB, Bx, By, Bz)
    mass_after = E**2 - Px**2 - Py**2 - Pz**2

    if (abs(mass_before - mass_after) / abs(mass_before) > tol) then
        write(*,*) 'FAIL Test 2: Invariant mass not preserved'
        write(*,*) '  Before:', mass_before, ' After:', mass_after
        failures = failures + 1
    else
        write(*,*) 'PASS Test 2: Invariant mass preserved'
    end if

    ! =========================================================================
    ! Test 3: Boost then inverse boost (roundtrip)
    ! Apply boost with beta, then with -beta, should recover original
    ! =========================================================================
    total = total + 1
    E = 3.5; Px = 0.5; Py = -1.0; Pz = 2.0
    PP = sqrt(Px**2 + Py**2 + Pz**2)
    E_save = E; Px_save = Px; Py_save = Py; Pz_save = Pz

    ! Forward boost (BB = sqrt(0.1^2 + 0.2^2 + 0.3^2))
    Bx = 0.1; By = 0.2; Bz = 0.3
    BB = sqrt(Bx**2 + By**2 + Bz**2)
    call TL(E, PP, Px, Py, Pz, BB, Bx, By, Bz)

    ! Inverse boost (negate beta components)
    call TL(E, PP, Px, Py, Pz, BB, -Bx, -By, -Bz)

    if (abs(E - E_save) / E_save > tol) then
        write(*,*) 'FAIL Test 3a: E not recovered after roundtrip'
        write(*,*) '  Expected:', E_save, ' Got:', E
        failures = failures + 1
    else if (abs(Px - Px_save) / (abs(Px_save) + 1e-10) > tol) then
        write(*,*) 'FAIL Test 3b: Px not recovered'
        write(*,*) '  Expected:', Px_save, ' Got:', Px
        failures = failures + 1
    else if (abs(Py - Py_save) / (abs(Py_save) + 1e-10) > tol) then
        write(*,*) 'FAIL Test 3c: Py not recovered'
        write(*,*) '  Expected:', Py_save, ' Got:', Py
        failures = failures + 1
    else if (abs(Pz - Pz_save) / (abs(Pz_save) + 1e-10) > tol) then
        write(*,*) 'FAIL Test 3d: Pz not recovered'
        write(*,*) '  Expected:', Pz_save, ' Got:', Pz
        failures = failures + 1
    else
        write(*,*) 'PASS Test 3: Roundtrip boost recovers original 4-vector'
    end if

    ! =========================================================================
    ! Test 4: Boost along x with beta = 0.8 on particle at rest
    ! =========================================================================
    total = total + 1
    E = 0.938; PP = 0.0; Px = 0.0; Py = 0.0; Pz = 0.0
    BB = 0.8; Bx = 0.8; By = 0.0; Bz = 0.0
    gamma = 1.0 / sqrt(1.0 - BB**2)
    expected_E = gamma * E

    call TL(E, PP, Px, Py, Pz, BB, Bx, By, Bz)

    if (abs(E - expected_E) / expected_E > tol) then
        write(*,*) 'FAIL Test 4a: E after x-boost'
        write(*,*) '  Expected:', expected_E, ' Got:', E
        failures = failures + 1
    else if (abs(Px - (-gamma * 0.8 * 0.938)) / abs(gamma * 0.8 * 0.938) > tol) then
        write(*,*) 'FAIL Test 4b: Px after x-boost'
        failures = failures + 1
    else if (abs(Py) > tol .or. abs(Pz) > tol) then
        write(*,*) 'FAIL Test 4c: Py/Pz should be 0'
        failures = failures + 1
    else
        write(*,*) 'PASS Test 4: Boost along x on proton at rest'
    end if

    ! =========================================================================
    ! Test 5: Small beta boost (numerical stability)
    ! =========================================================================
    total = total + 1
    E = 1.5; Px = 0.1; Py = 0.2; Pz = 0.3
    PP = sqrt(Px**2 + Py**2 + Pz**2)
    mass_before = E**2 - Px**2 - Py**2 - Pz**2

    BB = 0.001; Bx = 0.0; By = 0.0; Bz = 0.001

    call TL(E, PP, Px, Py, Pz, BB, Bx, By, Bz)
    mass_after = E**2 - Px**2 - Py**2 - Pz**2

    if (abs(mass_before - mass_after) / abs(mass_before) > tol) then
        write(*,*) 'FAIL Test 5: Invariant mass not preserved at small beta'
        write(*,*) '  Before:', mass_before, ' After:', mass_after
        failures = failures + 1
    else
        write(*,*) 'PASS Test 5: Small beta boost preserves invariant mass'
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

end program test_tl
