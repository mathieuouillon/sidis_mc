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
