!==============================================================================
! MODULE: simulation_helpers_module
! Purpose: Helper subroutines for the main simulation loop
!==============================================================================
module simulation_helpers_module
    implicit none

contains

! ------------------------------------------------------------------------------
! Initialize random number generator
! ------------------------------------------------------------------------------
subroutine init_random()
    use config_module, only: user_seed
    use pythia_commons, only: MRPY, RRPY
    implicit none

    integer :: i                          ! Loop counter
    real :: test                          ! Dummy variable
    real :: ranf                          ! Random number generator from CERNLIB
    integer :: initrm1, initrm2          ! To initialize ranf
    real(kind=8) :: PYR
    integer :: today(3), now(3)
    integer :: iseed                      ! Effective seed value

    if (user_seed >= 0) then
        ! Deterministic seeding for reproducibility
        iseed = user_seed
        call ranset(iseed)

        MRPY(2) = 0
        MRPY(3) = mod(iseed, 85635)
        MRPY(4) = mod(iseed, 67)
        MRPY(5) = mod(iseed, 56)

        do i = 1, 100
            RRPY(i) = ranf(0)
        end do

        ! Fixed warm-up count for reproducibility
        do i = 1, 1000
            test = ranf(0)
            test = PYR(0)
        end do

        ! Seed the Fortran intrinsic RNG (used by calculate_vz_position)
        call srand(iseed)
    else
        ! Original time-based seeding
        call idate(today)
        call itime(now)

        test = real(now(1)*now(2))/real(now(3))

        call datime(initrm1, initrm2)
        call ranset(initrm1*initrm2*int(test))

        MRPY(2) = 0
        MRPY(3) = mod(initrm1, 85635)
        MRPY(4) = mod(initrm1, 67)
        MRPY(5) = mod(initrm2, 56)

        do i = 1, 100
            RRPY(i) = ranf(0)
        end do

        do i = 1, initrm2*int(test)
            test = ranf(0)
            test = PYR(0)
        end do
    end if

end subroutine init_random

! ------------------------------------------------------------------------------
! Initialize kinematics values
! ------------------------------------------------------------------------------
subroutine init_kin()
    use kinematics_module, only: EEe, PPe, Pex, Pey, Pez, EEn, PPn, Pnx, Pny, Pnz, &
        ThFM, PhiFM, Kf, ele_ene, ele_the, ele_phi, nuc_mom, nuc_the, nuc_phi
    use config_module, only: E0
    implicit none

    real, parameter :: electron_mass = 0.000511
    real, parameter :: nucleon_mass = 0.938

    PPe = E0
    EEe = sqrt(PPe**2 + electron_mass**2)
    Pex = 0.0
    Pey = 0.0
    Pez = E0
    PPn = Kf
    EEn = sqrt(PPn**2 + nucleon_mass**2)
    Pnx = sin(ThFM)*cos(PhiFM)*PPn
    Pny = sin(ThFM)*sin(PhiFM)*PPn
    Pnz = cos(ThFM)*PPn

    ele_ene = EEe
    ele_the = 0.0
    ele_phi = 0.0
    nuc_mom = Kf
    nuc_the = ThFM
    nuc_phi = PhiFM

end subroutine init_kin

! ------------------------------------------------------------------------------
! Create spectators
! ------------------------------------------------------------------------------
subroutine create_spec()
    use kinematics_module, only: nuc_the, nuc_phi, nuc_mom
    use fermi_motion_module, only: FMintact
    use config_module, only: iTg, nucleon, specId
    use pythia_commons, only: N, K, P
    implicit none

    ! Decide if there is a spectator (RW model only)
    if (rand(0) > FMintact) return

    ! Decide the kind of spectator
    if (iTg == 1 .and. nucleon == 2112) then
        specId = 2212
    else if (iTg == 1 .and. nucleon == 2212) then
        specId = 2112
    else if (iTg == 2 .and. nucleon == 2112) then
        specId = 1000010020
    else if (iTg == 2 .and. nucleon == 2212) then
        specId = 1000000020
    else if (iTg == 3 .and. nucleon == 2112) then
        specId = 1000020020
    else if (iTg == 3 .and. nucleon == 2212) then
        specId = 1000010020
    else if (iTg == 4 .and. nucleon == 2112) then
        specId = 1000020030
    else if (iTg == 4 .and. nucleon == 2212) then
        specId = 1000010030
    end if

    ! Add Spectator
    N = N + 1

    k(N, 1) = 1
    k(N, 2) = specId
    k(N, 3) = 2

    p(N, 1) = -sin(nuc_the)*cos(nuc_phi)*nuc_mom
    p(N, 2) = -sin(nuc_the)*sin(nuc_phi)*nuc_mom
    p(N, 3) = -cos(nuc_the)*nuc_mom

    select case (k(N, 2))
    case (2212)
        p(N, 5) = 0.938272
    case (2112)
        p(N, 5) = 0.939566
    case (1000000020)
        p(N, 5) = 1.87913
    case (1000010020)
        p(N, 5) = 1.876124
    case (1000020020)
        p(N, 5) = 1.87654
    case (1000020030)
        p(N, 5) = 2.809356
    case (1000010030)
        p(N, 5) = 2.809356
    end select

    P(N, 4) = sqrt(P(N, 1)**2 + P(N, 2)**2 + P(N, 3)**2 + P(N, 5)**2)

end subroutine create_spec

! ------------------------------------------------------------------------------
! Lorentz transformation for Fermi motion
! ------------------------------------------------------------------------------
subroutine lorentz_fm(transform_type)
    use kinematics_module, only: EEe, PPe, Pex, Pey, Pez, EEn, PPn, Pnx, Pny, Pnz, &
        BB1, B1x, B1y, B1z, Thi, Phi
    implicit none

    integer, intent(in) :: transform_type

    if (transform_type == 1) then
        BB1 = PPn/EEn
        B1x = Pnx/EEn
        B1y = Pny/EEn
        B1z = Pnz/EEn

        call TL(EEe, PPe, Pex, Pey, Pez, BB1, B1x, B1y, B1z)
        call TL(EEn, PPn, Pnx, Pny, Pnz, BB1, B1x, B1y, B1z)

        ! Rotate around y
        Thi = -atan2(Pex, Pez)
        call InitRotY(Thi)

        ! Rotate around z
        Phi = -atan2(Pey, Pez)
        call InitRotZ(Phi)
    end if

end subroutine lorentz_fm

! ------------------------------------------------------------------------------
! Inverse Lorentz transformation for Fermi motion
! ------------------------------------------------------------------------------
subroutine lorentz_fm_back(transform_type)
    use kinematics_module, only: BB1, B1x, B1y, B1z, Thi, Phi
    use pythia_commons, only: N, P
    implicit none

    integer, intent(in) :: transform_type

    real :: mom1, mom2, mom3, mom4, ppp
    integer :: ip

    ! Rotate around z
    if (transform_type == 1) call FinalRotZ(-Phi)

    ! Rotate around y
    if (transform_type == 1) call FinalRotY(-Thi)

    ! Lorentz boost of all the particles
    do ip = 1, N
        mom1 = p(ip, 1)
        mom2 = p(ip, 2)
        mom3 = p(ip, 3)
        mom4 = p(ip, 4)

        if (transform_type == 1) then
            call TL(mom4, ppp, mom1, mom2, mom3, -BB1, -B1x, -B1y, -B1z)
        end if

        p(ip, 1) = mom1
        p(ip, 2) = mom2
        p(ip, 3) = mom3
        p(ip, 4) = mom4
    end do

end subroutine lorentz_fm_back

end module simulation_helpers_module
