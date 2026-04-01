!==============================================================================
! MODULE: c_interface
! Purpose: ISO_C_BINDING wrappers exposing Fortran physics routines to C/C++
!==============================================================================
module c_interface
    use iso_c_binding
    use simulation_helpers_module, only: init_random, init_kin, create_spec, &
        lorentz_fm, lorentz_fm_back
    implicit none

contains

    ! =========================================================================
    ! Configuration setters
    ! =========================================================================

    subroutine farm_set_config(c_nevent, c_iTg, c_iFM, c_E0, c_nkin, c_seed) &
            bind(C, name="farm_set_config")
        use config_module, only: nevent, iTg, iFM, E0, nkin, user_seed
        integer(c_int), value, intent(in) :: c_nevent, c_iTg, c_iFM, c_nkin, c_seed
        real(c_float), value, intent(in) :: c_E0

        nevent = c_nevent
        iTg = c_iTg
        iFM = c_iFM
        E0 = c_E0
        nkin = c_nkin
        user_seed = c_seed
    end subroutine farm_set_config

    subroutine farm_set_physics_params(c_FMlimit, c_iIso, c_iNS, c_iLund, &
            c_iQuenching, c_iSim, c_iqw, c_qhat, c_ehat, c_iDens, &
            c_iqg, c_iEg, c_iPtF) bind(C, name="farm_set_physics_params")
        use config_module, only: FMlimit, iIso, iNS, iLund, iQuenching, iSim, &
            iqg, iEg, iPtF, qhat, ehat, iDens, SupFac
        use quenching_module, only: alphas, iqw, scor, ncor, sfthrd
        integer(c_int), value, intent(in) :: c_iIso, c_iNS, c_iLund, &
            c_iQuenching, c_iSim, c_iqw, c_iDens, c_iqg, c_iEg, c_iPtF
        real(c_float), value, intent(in) :: c_FMlimit, c_qhat, c_ehat

        FMlimit = c_FMlimit
        iIso = c_iIso
        iNS = c_iNS
        iLund = c_iLund
        iQuenching = c_iQuenching
        iSim = c_iSim
        iqw = c_iqw
        alphas = 1.0d0 / 3.0d0
        scor = 1
        ncor = 0
        sfthrd = 1
        qhat = c_qhat
        ehat = c_ehat
        iDens = c_iDens
        iqg = c_iqg
        iEg = c_iEg
        iPtF = c_iPtF
        SupFac = qhat / (qhat + ehat)
    end subroutine farm_set_physics_params

    subroutine farm_set_nucleon(c_nucleon) bind(C, name="farm_set_nucleon")
        use config_module, only: nucleon
        integer(c_int), value, intent(in) :: c_nucleon
        nucleon = c_nucleon
    end subroutine farm_set_nucleon

    subroutine farm_set_nkin(c_nkin) bind(C, name="farm_set_nkin")
        use config_module, only: nkin
        integer(c_int), value, intent(in) :: c_nkin
        nkin = c_nkin
    end subroutine farm_set_nkin

    ! =========================================================================
    ! Configuration getters
    ! =========================================================================

    subroutine farm_get_state(c_ievent, c_PPe, c_rFM, c_iZ, c_iA, c_nucleon) &
            bind(C, name="farm_get_state")
        use kinematics_module, only: ievent, PPe
        use fermi_motion_module, only: iZ, iA
        use config_module, only: rFM, nucleon
        integer(c_int), intent(out) :: c_ievent, c_iZ, c_iA, c_nucleon
        real(c_float), intent(out) :: c_PPe, c_rFM

        c_ievent = ievent
        c_PPe = PPe
        c_rFM = rFM
        c_iZ = iZ
        c_iA = iA
        c_nucleon = nucleon
    end subroutine farm_get_state

    subroutine farm_set_ievent(c_ievent) bind(C, name="farm_set_ievent")
        use kinematics_module, only: ievent
        integer(c_int), value, intent(in) :: c_ievent
        ievent = c_ievent
    end subroutine farm_set_ievent

    subroutine farm_get_xsec99(c_xsec) bind(C, name="farm_get_xsec99")
        use pythia_commons, only: XSEC
        real(c_double), intent(out) :: c_xsec
        c_xsec = XSEC(99, 1)
    end subroutine farm_get_xsec99

    subroutine farm_get_qw_stats(c_qhat, c_nb) bind(C, name="farm_get_qw_stats")
        use quenching_module, only: QW_qhat, QW_nb
        real(c_float), intent(out) :: c_qhat
        integer(c_int), intent(out) :: c_nb
        c_qhat = QW_qhat
        c_nb = QW_nb
    end subroutine farm_get_qw_stats

    ! =========================================================================
    ! Event data access
    ! =========================================================================

    subroutine farm_get_event_data(c_nb_part, c_ids, c_charges, c_mothers, &
            c_px, c_py, c_pz, c_E, c_m, max_part) bind(C, name="farm_get_event_data")
        use particles_module, only: Nb_part, id_part, ch_part, id_mother, &
            px_part, py_part, pz_part, E_part, m_part, MAX_PARTICLES
        integer(c_int), intent(out) :: c_nb_part
        integer(c_int), intent(in), value :: max_part
        integer(c_int), intent(out) :: c_ids(max_part), c_charges(max_part), c_mothers(max_part)
        real(c_float), intent(out) :: c_px(max_part), c_py(max_part), c_pz(max_part)
        real(c_float), intent(out) :: c_E(max_part), c_m(max_part)
        integer :: j, n

        c_nb_part = Nb_part
        n = min(Nb_part, max_part)
        do j = 1, n
            c_ids(j) = id_part(j)
            c_charges(j) = ch_part(j)
            c_mothers(j) = id_mother(j)
            c_px(j) = px_part(j)
            c_py(j) = py_part(j)
            c_pz(j) = pz_part(j)
            c_E(j) = E_part(j)
            c_m(j) = m_part(j)
        end do
    end subroutine farm_get_event_data

    ! =========================================================================
    ! PYTHIA control
    ! =========================================================================

    subroutine farm_set_mstj1(c_val) bind(C, name="farm_set_mstj1")
        use pythia_commons, only: MSTJ
        integer(c_int), value, intent(in) :: c_val
        MSTJ(1) = c_val
    end subroutine farm_set_mstj1

    subroutine farm_pyinit_proton(c_energy) bind(C, name="farm_pyinit_proton")
        real(c_double), value, intent(in) :: c_energy
        real(kind=8) :: energy
        energy = c_energy
        call pyinit('FIXT', 'gamma/e-', 'p+', energy)
    end subroutine farm_pyinit_proton

    subroutine farm_pyinit_neutron(c_energy) bind(C, name="farm_pyinit_neutron")
        real(c_double), value, intent(in) :: c_energy
        real(kind=8) :: energy
        energy = c_energy
        call pyinit('FIXT', 'gamma/e-', 'n0', energy)
    end subroutine farm_pyinit_neutron

    subroutine farm_pyevnt() bind(C, name="farm_pyevnt")
        call pyevnt()
    end subroutine farm_pyevnt

    subroutine farm_pyexec() bind(C, name="farm_pyexec")
        call pyexec()
    end subroutine farm_pyexec

    ! =========================================================================
    ! Physics initialization
    ! =========================================================================

    subroutine farm_init_nucl() bind(C, name="farm_init_nucl")
        call InitNucl()
    end subroutine farm_init_nucl

    subroutine farm_init_fm() bind(C, name="farm_init_fm")
        call InitFM()
    end subroutine farm_init_fm

    subroutine farm_gen_nuc_dens() bind(C, name="farm_gen_nuc_dens")
        call GenNucDens()
    end subroutine farm_gen_nuc_dens

    subroutine farm_init_random() bind(C, name="farm_init_random")
        call init_random()
    end subroutine farm_init_random

    ! =========================================================================
    ! Physics event generation
    ! =========================================================================

    subroutine farm_fm_param() bind(C, name="farm_fm_param")
        call FMParam()
    end subroutine farm_fm_param

    subroutine farm_init_kin() bind(C, name="farm_init_kin")
        call init_kin()
    end subroutine farm_init_kin

    subroutine farm_lorentz_fm(c_type) bind(C, name="farm_lorentz_fm")
        integer(c_int), value, intent(in) :: c_type
        call lorentz_fm(c_type)
    end subroutine farm_lorentz_fm

    subroutine farm_lorentz_fm_back(c_type) bind(C, name="farm_lorentz_fm_back")
        integer(c_int), value, intent(in) :: c_type
        call lorentz_fm_back(c_type)
    end subroutine farm_lorentz_fm_back

    subroutine farm_pythia_config_dis() bind(C, name="farm_pythia_config_dis")
        call PythiaConfigDIS()
    end subroutine farm_pythia_config_dis

    subroutine farm_init_kin2book() bind(C, name="farm_init_kin2book")
        call InitKin2Book()
    end subroutine farm_init_kin2book

    subroutine farm_compute_v() bind(C, name="farm_compute_v")
        call ComputV()
    end subroutine farm_compute_v

    subroutine farm_inter_pos() bind(C, name="farm_inter_pos")
        call InterPos()
    end subroutine farm_inter_pos

    subroutine farm_apply_qw() bind(C, name="farm_apply_qw")
        call ApplyQW()
    end subroutine farm_apply_qw

    subroutine farm_create_spec() bind(C, name="farm_create_spec")
        call create_spec()
    end subroutine farm_create_spec

    subroutine farm_calc_vz(c_target_type, c_vz) bind(C, name="farm_calc_vz")
        use misc_module, only: vz
        integer(c_int), value, intent(in) :: c_target_type
        real(c_float), intent(out) :: c_vz
        call calculate_vz_position(c_target_type, vz)
        c_vz = vz
    end subroutine farm_calc_vz

    ! =========================================================================
    ! Combined initialization (avoids C/Fortran I/O interleaving issues)
    ! =========================================================================

    subroutine farm_initialize(c_rFM_out, c_iZ_out, c_iA_out) &
            bind(C, name="farm_initialize")
        use config_module, only: rFM
        use fermi_motion_module, only: iZ, iA
        real(c_float), intent(out) :: c_rFM_out
        integer(c_int), intent(out) :: c_iZ_out, c_iA_out

        call InitNucl()
        if (rFM /= 0) then
            call InitFM()
            call GenNucDens()
        end if
        call init_random()

        c_rFM_out = rFM
        c_iZ_out = iZ
        c_iA_out = iA
    end subroutine farm_initialize

    ! =========================================================================
    ! High-level event loop wrappers
    ! =========================================================================

    !> Setup kinematics: Fermi motion sampling, Lorentz boost, PYTHIA init.
    !> Retries until beam energy >= 4 GeV. Handles proton/neutron selection.
    subroutine farm_setup_kinematics(c_ievent, c_nevent) &
            bind(C, name="farm_setup_kinematics")
        use kinematics_module, only: PPe, Kf, ThFM, PhiFM
        use fermi_motion_module, only: iZ, iA
        use config_module, only: rFM, iFM, iQuenching, iSim, iIso, nucleon, nevent
        use pythia_commons, only: MSTJ, XSEC
        integer(c_int), value, intent(in) :: c_ievent, c_nevent

        real(kind=8) :: beam_energy

        ! Restart kinematics loop until PPe >= 4
        do
            ! Set nucleon type based on Z/A fraction
            if (c_ievent < (c_nevent * iZ / iA)) then
                nucleon = 2212
            else
                nucleon = 2112
            end if

            ! Fermi motion sampling
            if (rFM /= 0 .and. iFM /= 0) call FMParam()

            ! Initialize electron-nucleon kinematics
            call init_kin()

            ! Lorentz boost to nucleon rest frame
            if (iFM /= 0) call lorentz_fm(1)

            ! Exit if beam energy is sufficient
            if (PPe >= 4.0) exit
        end do

        ! Configure PYTHIA
        call PythiaConfigDIS()

        ! Stop fragmentation for quenching
        if (iQuenching /= 0) MSTJ(1) = 0

        ! PYTHIA initialization
        beam_energy = PPe
        if ((iIso == 0 .and. c_ievent < (c_nevent * iZ / iA)) .or. &
            (iIso == 1 .and. c_ievent < (c_nevent / 2))) then
            if (iSim /= 0) call pyinit('FIXT', 'gamma/e-', 'p+', beam_energy)
        else
            if (iSim /= 0) then
                write (*, *) 'X sec 99 = ', XSEC(99, 1)
                call pyinit('FIXT', 'gamma/e-', 'n0', beam_energy)
            end if
        end if
    end subroutine farm_setup_kinematics

    !> Generate and process one complete physics event:
    !> PYTHIA generation, Lorentz boost back, quenching, fragmentation,
    !> spectators, and output variable computation.
    subroutine farm_generate_event() bind(C, name="farm_generate_event")
        use config_module, only: iFM, iQuenching, iSim, iNS, iTg
        use pythia_commons, only: MSTJ
        implicit none

        ! Initialize kinematic variables for bookkeeping
        call InitKin2Book()

        ! PYTHIA event generation
        if (iSim /= 0) call pyevnt()

        ! Inverse Lorentz boost back to lab frame
        if (iFM /= 0) call lorentz_fm_back(1)

        ! Nuclear energy loss (quenching weights)
        if (iQuenching /= 0 .and. iTg > 1) then
            call InterPos()
            call ApplyQW()
        end if

        ! Fragmentation
        if (iQuenching /= 0) then
            MSTJ(1) = 1
            if (iSim /= 0) call pyexec()
            MSTJ(1) = 0
        end if

        ! Nuclear spectators
        if (iNS == 1 .and. (iTg >= 1 .or. iTg <= 4)) call create_spec()

        ! Compute output variables (DIS kinematics, hadron variables)
        call ComputV()
    end subroutine farm_generate_event

    !> Check if PYTHIA needs re-initialization for this event.
    subroutine farm_needs_reinit(c_ievent, c_nkin_counter, c_nevent, c_result) &
            bind(C, name="farm_needs_reinit")
        use fermi_motion_module, only: iZ, iA
        integer(c_int), value, intent(in) :: c_ievent, c_nkin_counter, c_nevent
        integer(c_int), intent(out) :: c_result

        if (c_ievent == 0 .or. c_nkin_counter == 0 .or. &
            c_ievent == (c_nevent * iZ / iA) .or. &
            mod(c_ievent, 500000) == 0) then
            c_result = 1
        else
            c_result = 0
        end if
    end subroutine farm_needs_reinit

    ! =========================================================================
    ! Timing (CERNLIB)
    ! =========================================================================

    subroutine farm_timex(c_t) bind(C, name="farm_timex")
        real(c_double), intent(out) :: c_t
        real(kind=8) :: t
        call timex(t)
        c_t = t
    end subroutine farm_timex

    ! =========================================================================
    ! LUND file I/O (Fortran-formatted for golden file compatibility)
    ! =========================================================================

    subroutine farm_open_lund(c_filename, c_len, c_unit, c_err) bind(C, name="farm_open_lund")
        character(c_char), intent(in) :: c_filename(*)
        integer(c_int), value, intent(in) :: c_len
        integer(c_int), intent(out) :: c_unit, c_err
        character(len=256) :: filename
        integer :: i

        filename = ''
        do i = 1, min(c_len, 256)
            filename(i:i) = c_filename(i)
        end do
        open(newunit=c_unit, file=trim(filename), iostat=c_err)
    end subroutine farm_open_lund

    subroutine farm_close_lund(c_unit) bind(C, name="farm_close_lund")
        integer(c_int), value, intent(in) :: c_unit
        close(c_unit)
    end subroutine farm_close_lund

    subroutine farm_write_event(c_unit) bind(C, name="farm_write_event")
        use particles_module, only: Nb_part, ch_part, id_part, id_mother, &
            px_part, py_part, pz_part, E_part, m_part
        use fermi_motion_module, only: iZ, iA
        use config_module, only: E0, nucleon
        use misc_module, only: vz
        integer(c_int), value, intent(in) :: c_unit
        integer :: l

        write (c_unit, '(I12,I12,I12,I12,I12,I12,F12.7,I17,I12,F12.8)') &
            Nb_part, iA, iZ, 0, 0, 11, E0, nucleon, 1, 1.0
        do l = 1, Nb_part
            write (c_unit, '(I12,I12,I12,I12,I12,I12,E16.8,E16.8,E16.8,E16.8,E16.8,I15,I12,F12.8)') &
                l, ch_part(l), 1, id_part(l), id_mother(l), 0, &
                px_part(l), py_part(l), pz_part(l), E_part(l), m_part(l), 0, 0, vz
        end do
    end subroutine farm_write_event

end module c_interface
