program monte_carlo_simulation
    use kinematics_module, only: ievent, PPe
    use particles_module, only: Nb_part, ch_part, id_part, id_mother, &
        px_part, py_part, pz_part, E_part, m_part
    use fermi_motion_module, only: iZ, iA
    use quenching_module, only: QW_nb, QW_qhat, alphas, iqw, scor, ncor, sfthrd
    use config_module, only: nevent, nkin, iTg, iFM, E0, user_seed, &
        iIso, iNS, iLund, iQuenching, iSim, nucleon, qhat, ehat, &
        iDens, iqg, iEg, iPtF, SupFac, FMlimit, rFM
    use misc_module, only: vz
    use pythia_commons, only: MSTJ, XSEC
    use simulation_helpers_module
    use cli_module
    implicit none

    ! Local variables
    real(kind=8) :: t1, t2               ! For time of computation
    real(kind=8) :: beam_energy          ! Input value for pythia
    integer(kind=8) :: i, l              ! Loop counters
    character(len=100) :: lund_file      ! Lund output file name
    integer :: io_err

    ! Command line argument variables
    integer :: num_args
    character(len=10) :: target_symbol

    ! Initialize flags
    nevent = 0
    lund_file = ''
    iTg = -1
    nkin = 0
    e0 = 0.0
    iFM = 5

    ! Parse command line arguments
    num_args = command_argument_count()
    call parse_command_line(num_args, nevent, lund_file, iTg, nkin, e0, iFM, user_seed)
    call get_symbol_name(iTg, target_symbol)

    write (*, *) 'Monte Carlo Simulation Parameters:'
    write (*, *) '  Number of events: ', nevent
    write (*, *) '  Lund output file: ', trim(lund_file)
    write (*, *) '  Target type:      ', trim(target_symbol)
    write (*, *) '  Nkin value:       ', nkin
    write (*, *) '  Electron energy:  ', e0, ' GeV'
    write (*, *) '  Fermi motion:     ', iFM
    if (user_seed >= 0) then
        write (*, *) '  Random seed:      ', user_seed
    else
        write (*, *) '  Random seed:       time-based'
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
    if (iLund == 1) then
        open (unit=59, file=lund_file, iostat=io_err)
        if (io_err /= 0) then
            write(*,*) 'ERROR: Cannot open output file: ', trim(lund_file)
            stop 1
        end if
    end if

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
