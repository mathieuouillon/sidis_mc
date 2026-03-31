program monte_carlo_simulation
    use kinematics_module
    use file_names_module
    use event_info_module
    use particles_module
    use fermi_motion_module
    use density_module
    use interaction_module
    use quenching_module
    use config_module
    use acceptance_module
    use misc_module
    use pythia_commons
    implicit none

    ! ------------------------------------------------------------------------------
    ! TO DO LIST:
    !   Implement some radiative effect
    !   Have alpha s changing with Q2 in QW
    ! ------------------------------------------------------------------------------

    ! Local variables
    real(kind=8) :: t1, t2               ! For time of computation
    real(kind=8) :: beam_energy          ! Input value for pythia
    integer(kind=8) :: i, l              ! Loop counters
    character(len=100) :: lund_file      ! Lund output file name

    ! Command line argument variables
    integer :: num_args, iostat, i_arg
    character(len=100) :: arg_str, next_arg
    character(len=10) :: target_symbol

    !logical :: nevent_set, output_set, target_set, nkin_set, e0_set

    ! Initialize flags
    nevent = 0
    lund_file = ''
    iTg = -1
    nkin = 0
    e0 = 0.0
    iColl = 0
    eColl = 0.0
    iFM = 5

    ! Parse command line arguments
    num_args = command_argument_count()
    call parse_command_line(num_args, nevent, lund_file, iTg, nkin, e0, iColl, eColl, iFM)
    call get_symbol_name(iTg, target_symbol)

    write (*, *) 'Monte Carlo Simulation Parameters:'
    write (*, *) '  Number of events: ', nevent
    write (*, *) '  Lund output file: ', trim(lund_file)
    write (*, *) '  Target type:      ', trim(target_symbol)
    write (*, *) '  Nkin value:       ', nkin
    write (*, *) '  Electron energy:  ', e0, ' GeV'
    write (*, *) '  Fermi motion:     ', iFM
    if (iColl == 1) then
        write (*, *) '  Collider mode:    Enabled'
        write (*, *) '  Collider energy:  ', eColl, ' GeV'
    else
        write (*, *) '  Collider mode:    Disabled'
    end if
    write (*, *) ''

    call timex(t1)

    ! Beginning of the simulation
    FMlimit = 0.5

    ! Isospin symmetry respected (0) or split at half (1)
    iIso = 0

    ! Nuclear spectator: 0 = no nuclear spectator, 1 = nuclear spectator
    ! This option is only for 2H and 4He targets
    iNS = 0

    ! Acceptance flags
    iAccept = 0                      ! CLAS12 Acceptance put 1
    iLund = 1                        ! Lund File

    ! Initialize quenching weights
    iQuenching = 1                   ! 0 deactivate Quenching
    iqw = 1                          ! 1 SW, 2 Arleo
    alphas = 1.0d0/3.0d0
    scor = 1
    ncor = 0
    sfthrd = 1
    qhat = 0.34                      ! Transport coefficient (GeV^2.fm^-1)
    ehat = 0.0                       ! Drag coefficient
    iDens = 1                        ! 0= hard sphere, 1= Wood Saxon param
    iqg = 1                          ! 1 -> quark and gluons are quenched other -> only q
    iEg = 0                          ! 1 -> a gluon is added to satisfy energy conservation
    iPtF = 3                         ! 0 -> no Pt; 1 -> from qhat; 2 -> BDMPS; 3 -> gluon angle from SW
    SupFac = qhat/(qhat + ehat)    ! Suppression factor

    iSim = 1                         ! 0 -> Turn off Pythia for tests

    ! To save time with useless Pythia initialization
    i = 0
    if (iTg == 0 .or. iFM == 0) then
        nkin = nevent + 1
    end if

    ! Initialize
    ievent = 0
    call InitNucl()
    if (rFM /= 0) then
        call InitFM()
        call GenNucDens()
    end if
    call init_random()
    if (iAccept == 1) call init_recoil()
    if (iLund == 1) open (unit=59, file=lund_file)

    ! Main Loop
    do while (ievent < nevent)

        ! Initialize the simulation
        if (ievent == 0 .or. i == nkin .or. ievent == int(nevent*iZ/iA) .or. &
            mod(ievent, 500000) == 0) then
            i = 0

            restart_kinematics: do
                ! Determine target: 1 for proton, 0 neutron
                if (ievent < (nevent*iZ/iA)) then
                    nucleon = 2212
                else
                    nucleon = 2112
                end if

                ! Randomize Theta Phi and Kf
                if (rFM /= 0 .and. iFM /= 0) call FMParam()

                ! Initialize the kinematics
                call init_kin()

                ! Going in the nucleon rest frame
                if (iColl /= 0) call lorentz_fm(2)
                if (iFM /= 0) call lorentz_fm(1)

                ! Security for low energies
                if (PPe >= 4) exit restart_kinematics
            end do restart_kinematics

            ! Parameters for Pythia
            call PythiaConfigDIS()

            ! Stop fragmentation for QW to be applied
            if (iQuenching /= 0) MSTJ(1) = 0

            ! PYTHIA init
            beam_energy = PPe
            if ((iIso == 0 .and. ievent < (nevent*iZ/iA)) .or. &
                (iIso == 1 .and. ievent < (nevent*0.5))) then
                if (iSim /= 0) call pyinit('FIXT', 'gamma/e-', 'p+', beam_energy)
            else
                if (iSim /= 0) write (*, *) 'X sec 99 = ', XSEC(99, 1)
                if (iSim /= 0) call pyinit('FIXT', 'gamma/e-', 'n0', beam_energy)
            end if
        end if

        ! Counter
        if (mod(ievent, 10000) == 0) write (*, *) ievent, ' events processed'

        ! Initialization of the kinematic variables
        call InitKin2Book()

        ! Event generation
        if (iSim /= 0) call pyevnt()

        ! Come back in target frame
        if (iFM /= 0) call lorentz_fm_back(1)

        ! Energy loss of the partons
        if (iQuenching /= 0 .and. iTg > 1) then
            call InterPos()              ! Pick the position of the interaction in the nuclei
            call ApplyQW()               ! Compute QW
        end if

        ! Fragmentation (if needed)
        if (iQuenching /= 0) then
            MSTJ(1) = 1
            if (iSim /= 0) call pyexec()
            MSTJ(1) = 0
        end if

        if (iNS == 1 .and. (iTg >= 1 .or. iTg <= 4)) call create_spec()

        ! Go back in lab frame (for collider mode)
        if (iColl /= 0) call lorentz_fm_back(2)

        ! Compute physical values for output
        call ComputV()

        ! Book the ntuple
        call calculate_vz_position(iTg, vz)

        if (iLund == 1) then
            write (59, '(I12,I12,I12,I12,I12,I12,F12.7,I17,I12,F12.8)') Nb_part, iA, iZ, 0, 0, 11, E0, nucleon, 1, 1.0
            do l = 1, Nb_part
                write (59, '(I12,I12,I12,I12,I12,I12,E16.8,E16.8,E16.8,E16.8,E16.8,I15,I12,F12.8)') l, &
                    ch_part(l), 1, id_part(l), id_mother(l), 0, &
                    px_part(l), py_part(l), pz_part(l), E_part(l), m_part(l), 0, 0, vz
            end do
        end if

        ievent = ievent + 1
        i = i + 1
    end do

    ! Close the file
    if (iLund == 1) close (59)

    write (*, *) 'X sec 99 = ', XSEC(99, 1)
    write (*, *) 'q hat = ', QW_qhat/QW_nb

    call timex(t2)
    write (*, *) nevent, ' events in ', t2 - t1, ' s'

end program monte_carlo_simulation

! ------------------------------------------------------------------------------
! Set target vertex z for Lund (RG-D; 2024)
! ------------------------------------------------------------------------------
subroutine calculate_vz_position(target_type, vz)
    implicit none

    integer, intent(in) :: target_type    ! Input target type
    real, intent(out) :: vz               ! Output vertex z position

    ! Local variables
    real, parameter :: first_pos = -2.5
    real, parameter :: second_pos = -7.5
    real :: target_pos
    real :: target_length
    real :: random_number

    random_number = rand()                 ! Choose a value in [0,1]
    vz = 0.0                              ! Default value for vz

    select case (target_type)
    case (10)
        ! Case for Sn
        target_pos = -3.5
        target_length = 0.018

    case (15)
        ! Case for Cu
        target_pos = -8.5
        target_length = 0.009

    case (7)
        ! Case for Carbon
        target_length = 0.2
        if (random_number < 0.5) then
            target_pos = -3.5
        else
            target_pos = -8.5
        end if

    case (1)
        ! Case for LD2
        target_pos = -5.0
        target_length = 5.0

    case default
        ! Do nothing for unspecified targets
        return
    end select

    vz = target_pos + target_length*(rand() - 0.5)

end subroutine calculate_vz_position

! ------------------------------------------------------------------------------
! Initialize random number generator
! ------------------------------------------------------------------------------
subroutine init_random()
    implicit none

    integer :: i                          ! Loop counter
    real :: test                          ! Dummy variable
    real :: ranf                          ! Random number generator from CERNLIB
    integer :: initrm1, initrm2          ! To initialize ranf
    integer :: MRPY(6)
    real(kind=8) :: RRPY(100), PYR
    integer :: today(3), now(3)

    common/PYDATR/MRPY, RRPY
    save/PYDATR/

    call idate(today)                     ! today(1)=day, (2)=month, (3)=year
    call itime(now)                       ! now(1)=hour, (2)=minute, (3)=second

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

end subroutine init_random

! ------------------------------------------------------------------------------
! Initialize kinematics values
! ------------------------------------------------------------------------------
subroutine init_kin()
    use kinematics_module
    use file_names_module
    use event_info_module
    use particles_module
    use fermi_motion_module
    use density_module
    use interaction_module
    use quenching_module
    use config_module
    use acceptance_module
    use misc_module
    use pythia_commons
    implicit none

    real, parameter :: electron_mass = 0.000511
    real, parameter :: nucleon_mass = 0.938

    ! include 'common.f90'

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
    use kinematics_module
    use file_names_module
    use event_info_module
    use particles_module
    use fermi_motion_module
    use density_module
    use interaction_module
    use quenching_module
    use config_module
    use acceptance_module
    use misc_module
    use pythia_commons
    implicit none

    !include 'common.f90'

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
    use kinematics_module
    use file_names_module
    use event_info_module
    use particles_module
    use fermi_motion_module
    use density_module
    use interaction_module
    use quenching_module
    use config_module
    use acceptance_module
    use misc_module
    use pythia_commons
    implicit none

    integer, intent(in) :: transform_type
    ! include 'common.f90'

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

    else if (transform_type == 2) then
        BB2 = -sqrt(EColl**2 - 0.939**2)/EColl
        B2x = 0.0
        B2y = 0.0
        B2z = BB2

        call TL(EEe, PPe, Pex, Pey, Pez, BB2, B2x, B2y, B2z)
    end if

end subroutine lorentz_fm

! ------------------------------------------------------------------------------
! Inverse Lorentz transformation for Fermi motion
! ------------------------------------------------------------------------------
subroutine lorentz_fm_back(transform_type)
    use kinematics_module
    use file_names_module
    use event_info_module
    use particles_module
    use fermi_motion_module
    use density_module
    use interaction_module
    use quenching_module
    use config_module
    use acceptance_module
    use misc_module
    use pythia_commons
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
        else if (transform_type == 2) then
            call TL(mom4, ppp, mom1, mom2, mom3, -BB2, -B2x, -B2y, -B2z)
        end if

        p(ip, 1) = mom1
        p(ip, 2) = mom2
        p(ip, 3) = mom3
        p(ip, 4) = mom4
    end do

end subroutine lorentz_fm_back

subroutine parse_command_line(num_args, nevent, lund_file, iTg, nkin, e0, iColl, eColl, iFM)
    implicit none

    integer, intent(in) :: num_args
    integer, intent(out) :: nevent
    character(len=100), intent(out) :: lund_file
    integer, intent(out) :: iTg, nkin, iColl, iFM
    real(kind=4), intent(out) :: e0, eColl

    ! Local variables
    integer :: i_arg, iostat
    character(len=100) :: arg_str, next_arg
    logical :: nevent_set, output_set, target_set, nkin_set, e0_set, iColl_set, eColl_set, iFM_set

    ! Initialize flags
    nevent_set = .false.
    output_set = .false.
    target_set = .false.
    nkin_set = .false.
    e0_set = .false.
    iColl_set = .false.
    eColl_set = .false.
    iFM_set = .false.
    nevent = 0
    lund_file = ''
    iTg = -1
    nkin = 0
    e0 = 0.0
    iColl = 0
    eColl = 0.0
    iFM = 5

    if (num_args == 0) then
        write (*, *) 'Usage: ./monte_carlo_simulation --nevent <number_of_events> --output <lund_output_file> &
           --target <target_type> --nkin <nkin_value> --e0 <electron_energy> --iColl <collider_option> --eColl <collider_energy> &
           --iFM <fermi_motion_option>'
        write (*, *) 'Example: ./monte_carlo_simulation --nevent 10000 --output 120k_D.txt --target 1 --nkin 20000 --e0 10.5 &
         --iColl 1 --eColl 0.5 --iFM 5'
        write (*, *) 'Use --help for detailed information about all options'
        stop 1
    end if

    ! Parse named arguments
    i_arg = 1
    do while (i_arg <= num_args)
        call get_command_argument(i_arg, arg_str)

        select case (trim(arg_str))
        case ('--nevent', '-n')
            if (i_arg == num_args) then
                write (*, *) 'Error: --nevent requires a value'
                stop 1
            end if
            call get_command_argument(i_arg + 1, next_arg)
            read (next_arg, *, iostat=iostat) nevent
            if (iostat /= 0 .or. nevent <= 0) then
                write (*, *) 'Error: Invalid number of events. Must be a positive integer.'
                write (*, *) 'Provided: ', trim(next_arg)
                stop 1
            end if
            nevent_set = .true.
            i_arg = i_arg + 2

        case ('--output', '-o')
            if (i_arg == num_args) then
                write (*, *) 'Error: --output requires a value'
                stop 1
            end if
            call get_command_argument(i_arg + 1, lund_file)
            if (len_trim(lund_file) == 0) then
                write (*, *) 'Error: Lund file name cannot be empty.'
                stop 1
            end if
            output_set = .true.
            i_arg = i_arg + 2

        case ('--target', '-t')
            if (i_arg == num_args) then
                write (*, *) 'Error: --target requires a value'
                stop 1
            end if
            call get_command_argument(i_arg + 1, next_arg)
            read (next_arg, *, iostat=iostat) iTg
            if (iostat /= 0 .or. iTg < 0 .or. iTg > 15) then
                write (*, *) 'Error: Invalid target type. Must be an integer between 0 and 15.'
                write (*, *) 'Target types: 0->p, 1->2H, 2->3H, 3->3He, 4->4He, 5->6Li,'
                write (*, *) '              6->7Li, 7->C, 8->Al, 9->Fe, 10->Sn, 11->Pb,'
                write (*, *) '              12->Ne, 13->Kr, 14->Xe, 15->Cu'
                write (*, *) 'Provided: ', trim(next_arg)
                stop 1
            end if
            target_set = .true.
            i_arg = i_arg + 2

        case ('--nkin')
            if (i_arg == num_args) then
                write (*, *) 'Error: --nkin requires a value'
                stop 1
            end if
            call get_command_argument(i_arg + 1, next_arg)
            read (next_arg, *, iostat=iostat) nkin
            if (iostat /= 0 .or. nkin <= 0) then
                write (*, *) 'Error: Invalid nkin value. Must be a positive integer.'
                write (*, *) 'Provided: ', trim(next_arg)
                stop 1
            end if
            nkin_set = .true.
            i_arg = i_arg + 2

        case ('--e0')
            if (i_arg == num_args) then
                write (*, *) 'Error: --e0 requires a value'
                stop 1
            end if
            call get_command_argument(i_arg + 1, next_arg)
            read (next_arg, *, iostat=iostat) e0
            if (iostat /= 0 .or. e0 <= 0.0) then
                write (*, *) 'Error: Invalid electron energy. Must be a positive number.'
                write (*, *) 'Provided: ', trim(next_arg)
                stop 1
            end if
            e0_set = .true.
            i_arg = i_arg + 2

        case ('--iColl')
            if (i_arg == num_args) then
                write (*, *) 'Error: --iColl requires a value'
                stop 1
            end if
            call get_command_argument(i_arg + 1, next_arg)
            read (next_arg, *, iostat=iostat) iColl
            if (iostat /= 0 .or. (iColl /= 0 .and. iColl /= 1)) then
                write (*, *) 'Error: Invalid collider option. Must be 0 (no collider) or 1 (collider).'
                write (*, *) 'Provided: ', trim(next_arg)
                stop 1
            end if
            iColl_set = .true.
            i_arg = i_arg + 2

        case ('--eColl')
            if (i_arg == num_args) then
                write (*, *) 'Error: --eColl requires a value'
                stop 1
            end if
            call get_command_argument(i_arg + 1, next_arg)
            read (next_arg, *, iostat=iostat) eColl
            if (iostat /= 0 .or. eColl < 0.0) then
                write (*, *) 'Error: Invalid collider energy. Must be a non-negative number.'
                write (*, *) 'Provided: ', trim(next_arg)
                stop 1
            end if
            eColl_set = .true.
            i_arg = i_arg + 2

        case ('--iFM')
            if (i_arg == num_args) then
                write (*, *) 'Error: --iFM requires a value'
                stop 1
            end if
            call get_command_argument(i_arg + 1, next_arg)
            read (next_arg, *, iostat=iostat) iFM
            if (iostat /= 0 .or. iFM < 0 .or. iFM > 5) then
                write (*, *) 'Error: Invalid Fermi motion option. Must be an integer between 0 and 5.'
                write (*, *) 'Provided: ', trim(next_arg)
                stop 1
            end if
            iFM_set = .true.
            i_arg = i_arg + 2

        case ('--help', '-h')
            write (*, *) 'Monte Carlo Simulation Program'
            write (*, *) ''
            write (*, *) 'Usage: ./monte_carlo_simulation [OPTIONS]'
            write (*, *) ''
            write (*, *) 'Required Options:'
            write (*, *) '  --nevent, -n    Number of events to generate (positive integer)'
            write (*, *) '  --output, -o    Output file name for Lund format data'
            write (*, *) '  --target, -t    Target type (integer 0-15)'
            write (*, *) '  --nkin          Number of kinematic iterations (positive integer)'
            write (*, *) '  --e0            Electron energy in GeV (positive number)'
            write (*, *) ''
            write (*, *) 'Optional:'
            write (*, *) '  --help, -h      Show this help message'
            write (*, *) '  --iColl         Collider option: 0 (no collider), 1 (collider) [default: 0]'
            write (*, *) '  --eColl        Collider energy in GeV (non-negative number) [default: 0.0]'
            write (*, *) '  --iFM          Fermi motion option (integer 0-5) [default: 5]'
            write (*, *) ''
            write (*, *) 'Target Types:'
            write (*, *) '  0  -> p      (proton)'
            write (*, *) '  1  -> 2H     (deuterium)'
            write (*, *) '  2  -> 3H     (tritium)'
            write (*, *) '  3  -> 3He    (helium-3)'
            write (*, *) '  4  -> 4He    (helium-4)'
            write (*, *) '  5  -> 6Li    (lithium-6)'
            write (*, *) '  6  -> 7Li    (lithium-7)'
            write (*, *) '  7  -> C      (carbon)'
            write (*, *) '  8  -> Al     (aluminum)'
            write (*, *) '  9  -> Fe     (iron)'
            write (*, *) '  10 -> Sn     (tin)'
            write (*, *) '  11 -> Pb     (lead)'
            write (*, *) '  12 -> Ne     (neon)'
            write (*, *) '  13 -> Kr     (krypton)'
            write (*, *) '  14 -> Xe     (xenon)'
            write (*, *) '  15 -> Cu     (copper)'
            write (*, *) ''
            write (*, *) 'Fermi Motion Options:'
            write (*, *) '  0 -> No Fermi motion'
            write (*, *) '  1 -> like 4 plus a tail from [2] (deut, C, Al, Fe, Sn, Pb, 4He)'
            write (*, *) '  2 -> Accardi SVG (7Li, C, O, Ne, Al, Ar, Ca, Ni, Cu, Zr, Sn, Pb)'
            write (*, *) '  3 -> Accardi CS (2H, 3He, 4He, C, O, Ca, Fe, Pb)'
            write (*, *) '  4 -> Hard sphere with values from [1]'
            write (*, *) '  5 -> R. Wiringa et al. PRC 89, 024305 (2014)'
            write (*, *) ''    
            write (*, *) 'Notes:'
            write (*, *) '  - Fermi motion distributions are limited to 1 GeV nucleons.'
            write (*, *) '  - [1] E. J. Moniz et al. PRL 26, 445 (1971)'
            write (*, *) '  - [2] A. Bodek and J. L. Ritchie PRD 23, 1070 (1981)'
            write (*, *) ''
            write (*, *) 'Examples:'
            write (*, *) '  ./monte_carlo_simulation --nevent 10000 --output 120k_D.txt --target 1 --nkin 20000 --e0 10.5 &
            --iColl 1 --eColl 0.5 --iFM 5'
            write (*, *) '  ./monte_carlo_simulation -n 50000 -o results.txt -t 7 --nkin 15000 --e0 12.0'
            write (*, *) '  ./monte_carlo_simulation --target 0 --nevent 25000 --output proton.txt --nkin 10000 --e0 8.5'
            stop 0

        case default
            write (*, *) 'Error: Unknown argument: ', trim(arg_str)
            write (*, *) 'Use --help for usage information'
            stop 1
        end select
    end do

    ! Check if required arguments were provided
    if (.not. nevent_set) then
        write (*, *) 'Error: --nevent argument is required'
        write (*, *) 'Use --help for usage information'
        stop 1
    end if

    if (.not. output_set) then
        write (*, *) 'Error: --output argument is required'
        write (*, *) 'Use --help for usage information'
        stop 1
    end if

    if (.not. target_set) then
        write (*, *) 'Error: --target argument is required'
        write (*, *) 'Use --help for usage information'
        stop 1
    end if

    if (.not. nkin_set) then
        write (*, *) 'Error: --nkin argument is required'
        write (*, *) 'Use --help for usage information'
        stop 1
    end if

    if (.not. e0_set) then
        write (*, *) 'Error: --e0 argument is required'
        write (*, *) 'Use --help for usage information'
        stop 1
    end if
end subroutine parse_command_line


subroutine get_symbol_name(iTg, symbol_name)
    implicit none
    integer, intent(in) :: iTg
    character(len=10), intent(out) :: symbol_name

    select case (iTg)
    case (0)
        symbol_name = 'p'
    case (1)
        symbol_name = '2H'
    case (2)
        symbol_name = '3H'
    case (3)
        symbol_name = 'He3'
    case (4)
        symbol_name = 'He4'
    case (5)
        symbol_name = 'Li6'
    case (6)
        symbol_name = 'Li7'
    case (7)
        symbol_name = 'C'
    case (8)
        symbol_name = 'Al'
    case (9)
        symbol_name = 'Fe'
    case (10)
        symbol_name = 'Sn'
    case (11)
        symbol_name = 'Pb'
    case (12)
        symbol_name = 'Ne'
    case (13)
        symbol_name = 'Kr'
    case (14)
        symbol_name = 'Xe'
    case (15)
        symbol_name = 'Cu'
    case default
        symbol_name = 'Unknown'
    end select
end subroutine get_symbol_name
