module kinematics_module
    implicit none

    ! Initial Kinematic variables
    integer :: ievent
    real :: EEn, PPn, Pnx, Pny, Pnz           ! nucleon kinematics
    real :: EEe, PPe, Pex, Pey, Pez           ! electron kinematics
    real :: ele_ene, ele_the, ele_phi         ! values to book
    real :: nuc_mom, nuc_the, nuc_phi

    ! Transformation variables
    real :: ThFM, PhiFM, Kf                   ! Theta, Phi and K Fermi
    real :: BB1, B1x, B1y, B1z                ! Lorentz kinematics
    real :: BB2, B2x, B2y, B2z                ! Lorentz kinematics
    real :: Thi, Phi                          ! angles for rotation
    real :: Beta, ECoM                        ! For calculation of the CoM energy

end module kinematics_module

!==============================================================================
! MODULE: file_names_module
! Purpose: Output file names
!==============================================================================
module file_names_module
    implicit none
    character(len=60) :: bosout, hbookout
end module file_names_module

!==============================================================================
! MODULE: event_info_module
! Purpose: Event kinematic information
!==============================================================================
module event_info_module
    implicit none
    real :: Q22, W, Nu, XBj, y_ele
end module event_info_module

!==============================================================================
! MODULE: particles_module
! Purpose: Recoiled particles data
!==============================================================================
module particles_module
    implicit none

    integer, parameter :: MAX_PARTICLES = 100

    integer :: TrkGS
    integer :: Nb_part
    integer :: id_part(MAX_PARTICLES)
    integer :: id_mother(MAX_PARTICLES)
    integer :: acc_part(MAX_PARTICLES)
    integer :: ch_part(MAX_PARTICLES)

    real :: p_part(MAX_PARTICLES)
    real :: px_part(MAX_PARTICLES)
    real :: py_part(MAX_PARTICLES)
    real :: pz_part(MAX_PARTICLES)
    real :: E_part(MAX_PARTICLES)
    real :: m_part(MAX_PARTICLES)
    real :: z_part(MAX_PARTICLES)
    real :: th_part(MAX_PARTICLES)
    real :: tt_part(MAX_PARTICLES)
    real :: Pts_part(MAX_PARTICLES)
    real :: phih_part(MAX_PARTICLES)
    real :: phi_part(MAX_PARTICLES)
    real :: vx_part(MAX_PARTICLES)
    real :: Xf_part(MAX_PARTICLES)
    real :: Als(MAX_PARTICLES)

end module particles_module

!==============================================================================
! MODULE: fermi_motion_module
! Purpose: Fermi motion distribution tables
!==============================================================================
module fermi_motion_module
    implicit none

    integer, parameter :: FM_TABLE_SIZE = 1000

    real :: step_size_FM
    real :: FMintact
    real :: FM_table(FM_TABLE_SIZE)
    real :: FM_n(FM_TABLE_SIZE)
    real :: FM_p(FM_TABLE_SIZE)
    real :: FM_i(FM_TABLE_SIZE)

    integer :: iZ, iA  ! target Z and A
    integer :: FMnb

end module fermi_motion_module

!==============================================================================
! MODULE: density_module
! Purpose: Nuclear density distribution
!==============================================================================
module density_module
    implicit none

    integer, parameter :: DENSITY_TABLE_SIZE = 2000

    real :: density_table(DENSITY_TABLE_SIZE)
    real :: step_size_dens
    real :: quantity_table(DENSITY_TABLE_SIZE)
    real :: init_dens

end module density_module

!==============================================================================
! MODULE: interaction_module
! Purpose: Interaction position in nucleus
!==============================================================================
module interaction_module
    implicit none

    real :: x_inter, y_inter, z_inter
    real :: pos_radius, pos_theta, pos_phi

end module interaction_module

!==============================================================================
! MODULE: quenching_module
! Purpose: Quenching weight variables
!==============================================================================
module quenching_module
    implicit none

    ! Quenching Weight variables
    integer :: QW_nb
    real :: QW_wc, QW_R, QW_L, QW_w
    real :: QW_qhat, QW_chi, QW_th

    ! AA routine variables
    real(kind=8) :: alphas
    integer :: iqw, scor, ncor, sfthrd, irw

end module quenching_module

!==============================================================================
! MODULE: config_module
! Purpose: Simulation flags and configuration values
!==============================================================================
module config_module
    implicit none

    ! Target and physics flags
    integer :: iTg          ! Target type
    integer :: iFM          ! Fermi motion flag
    integer :: iDens        ! Density model
    integer :: iQuenching   ! Quenching flag
    integer :: iSim         ! Simulation flag
    integer :: iNS          ! Nuclear spectator
    integer :: iAccept      ! Acceptance flag
    integer :: iLund        ! Lund output flag
    integer :: iAlert       ! ALERT acceptance
    integer :: iIso         ! Isospin flag

    ! Event and kinematic parameters
    integer :: nevent       ! Number of events
    integer :: nkin         ! Kinematic iterations
    integer :: nucleon      ! Nucleon type (2212=p, 2112=n)
    integer :: specId       ! Spectator ID

    ! Collider parameters
    integer :: iColl        ! Collider flag
    real :: EColl           ! Collider energy

    ! Quenching parameters
    integer :: iqg          ! Quark/gluon quenching
    integer :: iEg          ! Energy conservation gluon
    integer :: iPtF         ! Pt broadening model

    ! Energy and momentum parameters
    real :: rFM             ! Fermi motion parameter
    real :: E0              ! Beam energy
    real :: qhat            ! Transport coefficient
    real :: ehat            ! Drag coefficient
    real :: FMlimit         ! Fermi momentum limit
    real :: SupFac          ! Suppression factor

    ! Seed for reproducibility (-1 = use time-based seed)
    integer :: user_seed = -1

end module config_module

!==============================================================================
! MODULE: acceptance_module
! Purpose: Recoil particle acceptance tables
!==============================================================================
module acceptance_module
    implicit none

    integer, parameter :: ACC_SIZE_1 = 25
    integer, parameter :: ACC_SIZE_2 = 25

    real :: pro_acc(ACC_SIZE_1, ACC_SIZE_2)
    real :: deu_acc(ACC_SIZE_1, ACC_SIZE_2)
    real :: tri_acc(ACC_SIZE_1, ACC_SIZE_2)
    real :: he3_acc(ACC_SIZE_1, ACC_SIZE_2)
    real :: he4_acc(ACC_SIZE_1, ACC_SIZE_2)

end module acceptance_module

!==============================================================================
! MODULE: misc_module
! Purpose: Miscellaneous variables and constants
!==============================================================================
module misc_module
    implicit none

    real :: vxz             ! Vertex position for GSIM
    real :: vz              ! Vertex z position
    real, parameter :: pi = 3.1415926535

    interface
        real function ranf(dummy)
            integer, intent(in) :: dummy
        end function ranf
    end interface
end module misc_module

!==============================================================================
! MODULE: pythia_commons
! Purpose: PYTHIA 6.4 common blocks (must remain as common blocks)
! Note: These CANNOT be converted to modules - they interface with PYTHIA
!==============================================================================
module pythia_commons
    implicit none

    ! PYPARS common block
    real(kind=8) :: PARP(200), PARI(200)
    integer :: MSTP(200), MSTI(200)
    COMMON/PYPARS/MSTP, PARP, MSTI, PARI

    ! PYDAT1 common block
    real(kind=8) :: PARU(200), PARJ(200)
    integer :: MSTU(200), MSTJ(200)
    COMMON/PYDAT1/MSTU, PARU, MSTJ, PARJ

    ! PYJETS common block
    real(kind=8) :: P(4000, 5), V(4000, 5)
    integer :: N, NPAD, K(4000, 5)
    COMMON/PYJETS/N, NPAD, K, P, V

    ! PYSUBS common block
    integer :: MSEL, MSELPD, MSUB(500), KFIN(2, -40:40)
    real(kind=8) :: CKIN(200)
    COMMON/PYSUBS/MSEL, MSELPD, MSUB, KFIN, CKIN

    ! PYINT5 common block
    integer :: NGENPD, NGEN(0:500, 3)
    real(kind=8) :: XSEC(0:500, 3)
    COMMON/PYINT5/NGENPD, NGEN, XSEC

end module pythia_commons
