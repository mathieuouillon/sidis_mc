!c------------------------------------------------------------------------------
!c Initialize the config values of LEPTO
!c------------------------------------------------------------------------------
      subroutine PythiaConfigAll
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

!ccccc Include all the common blocks
          ! include 'common.f90'

          MSEL = 2
!c Kind of possible multiple interaction (needed to avoid bugs) D = 4 H =1
          MSTP(82) = 1

!c Lowest CM energy D = 10 H = 3
          PARP(2) = 2 ! Modify because of the FM

!c remaining energy below witch the fragmentation is stopped D = 0.8
          PARJ(33) = 0.3

!c Lower limit for sqrt(s)
          CKIN(1) = 1.

!c To avoid random crash of init:
!c     CKIN(1) = 2
!c     CKIN(3) = 0
!c Q2 Limits
          CKIN(65) = 1.
          CKIN(66) = -1
!c W limit
          CKIN(77) = 1.90
          CKIN(78) = -1.

!c HERMES Params

! c     PARJ(1) = 0.02
! c     PARJ(2) = 0.25
! c     PARJ(11) = 0.51
! c     PARJ(12) = 0.57
! c     PARJ(21) = 0.42
! c     PARJ(33) = 0.47
! c     PARJ(41) = 0.68
! c     PARJ(42) = 0.35
! c     PARJ(45) = 0.74
!
! c     PARJ(21) = 0.44
! c     PARJ(23) = 0.01
! c     PARJ(24) = 2.0

      end

      subroutine PythiaConfigDIS
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

!ccccc Include all the common blocks
          ! include 'common.f90'
          !integer i

!c DIS only
          MSTP(14) = 26

!c Kind of possible multiple interaction (needed to avoid bugs) D = 4 H =1
          MSTP(82) = 0

!c Lowest CM energy D = 10 H = 3
          PARP(2) = 2

!c remaining energy below witch the fragmentation is stopped D = 0.8
          PARJ(33) = 0.3

!c Lower limit for sqrt(s)
          CKIN(1) = 2.

!c To avoid crash of init during long simulation:
          CKIN(1) = 2
          CKIN(3) = 0
!c Q2 Limits
          CKIN(65) = 1
          CKIN(66) = -1
!c W limit
          CKIN(77) = 2.
          CKIN(78) = -1.
      end

!c------------------------------------------------------------------------------
!c Initialize the config values of LEPTO
!c------------------------------------------------------------------------------
      subroutine PythiaConfigOWN
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

!ccccc Include all the common blocks
          ! include 'common.f90'

!ccc Deeply Inelastic Scattering and γ ∗ γ ∗ physics
!ccc MSEL = 1, 2, 35, 36, 37, 38
!c MSEL = 0 Tout les process doivent etre activer individuellement
!c MSEL = 1 in our case call MSTP(14)
!c MSEL = 2 include 1 + elastic (hadronic) + diffractive (hadronic)
!c                    + low Pt process
!c     MSEL = 2

! ccc MSUB(ISUB) = 0 : the subprocess is excluded.
!ccc MSUB(ISUB) = 1 : the subprocess is included.

! C ... Set kinematic cuts
! c min invariant mass = sqrt(s)
          CKIN(1) = 1.d0             ! lower lim
          CKIN(65) = 0.85d0           ! min Q2
          CKIN(77) = 1.8d0            ! min W

!c DIS only
          MSTP(14) = 26
!c Suppression of VMD (0:none)
          MSTP(20) = 0
!c Q2 definition override
          MSTP(22) = 4
!c retain x and Q2 of original scat.
          MSTP(23) = 0
!c Nb of massless quark
          MSTP(38) = 3
!c Master swich for decays
!c      MSTP(41)  = 1
!c Choice of the proton parton distrib.
          MSTP(51) = 6
!c Number of quark flavour in the PDF
          MSTP(58) = 3
!c QED and QCD  radiative effects
          MSTP(61) = 2
!c Master switch for multiple interaction
          MSTP(81) = 0
!c structure of multiple interaction
          MSTP(82) = 1
!c structure of diffractive system
          MSTP(101) = 1

!c ... Allow low c.m. energies
          PARP(2) = 2.d0    ! min(Ecm) = 2 GeV, default 10 GeV
          PARP(62) = .3d0    ! reduce the min. mass of time-like parton
          PARP(65) = .3d0    ! reduce the min. mass of time-like parton

!c Allow diffractive excitement above M+PARP(102)
          PARP(102) = .4d0
!c Minimum energy for elastic and difractive effects
          PARP(104) = .3d0

          PARP(111) = 0.d0    ! reduce the min remnant invariant mass, default 2 GeV

!c Minor change about vector meson production
          PARP(161) = 2.69d0
          PARP(162) = 24.6d0
          PARP(163) = 18.8d0
          PARP(165) = .33d0  ! reduce the min. mass of time-like parton

!c Stuff about the number of flavour
!c     MSTJ(45) = 3
!c     MSTU(112) = 3
!c     MSTU(113) = 2
!c     MSTU(114) = 3

!c remaining energy below witch the fragmentation is stopped D = 0.8
          PARJ(33) = 0.47

      end

!c------------------------------------------------------------------------------
!c Initialize the config values of LEPTO
!c------------------------------------------------------------------------------
      subroutine PythiaConfigHayk
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

!ccccc Include all the common blocks
          ! include 'common.f90'

!c Maximum number of generations D = 3 / H = 2
          MSTP(1) = 2
!c Calculation of alpha s D = H = 1
!c Comments in thesis let think he want to put 0...
!c      MSTP(2) = 1
!c Not in Pythia manual ...
          MSTP(5) = 0
!c Fix the range over wich electron emit photons D = 1 H = 2
!c text in manual seem to push for 1
          MSTP(13) = 2
!c structure of incoming photon beam D = H = 30
!c      MSTP(14) = 30
!c Choice of definition of the fractional part taken by photon D = H = 1
!c lets try 0
!c      MSTP(16) = 1
!c Possibility of a extra factor for processes involving resolved virtual photons
!c D = 4 H = 6   6 is not in the manual !
          MSTP(17) = 6
!c suppression of resolved (VMD or GVMD) cross sections D = 3 H = 0
          MSTP(20) = 0
!c Keep the final electron giving the right x and q2 D = H = 1
!c      MSTP(23) = 1
!c Handling of quark loops, number of allowed quarks D = 5 H = 4
          MSTP(38) = 4
!c master switch for decay  D = 2 H = 1
          MSTP(41) = 1
!c Choice of the PDF D = H = 7
          MSTP(51) = 7
!c Max number of quarks in pdf D = 5 H = 4
          MSTP(58) = 4

!c General switch for initial state radiations D = 2 H = 0
!c try to turn on 1 or 2
          MSTP(61) = 0
!c normaly no effect from 6x processes
          MSTP(62) = 3
          MSTP(63) = 2
          MSTP(64) = 2
          MSTP(65) = 1
          MSTP(66) = 5
          MSTP(67) = 2
          MSTP(68) = 1
          MSTP(69) = 0
!c Master switch for final state radiations D = 1 H = 0
          MSTP(71) = 0
!c Master switch for multiple interactions D = 1 H = 0
          MSTP(81) = 0
!c normaly no effect from 8x processes
          MSTP(82) = 1
          MSTP(83) = 100
          MSTP(86) = 2
!c
!c Change the energy partitionning between remnant D = 3 H = 4
          MSTP(92) = 4
!c Structure of diffractive system D = 3 H = 1
          MSTP(101) = 1
!c Calculation of kinematic coefficients D = 0 H = 1
          MSTP(121) = 1

!c Lowest CM energy D = 10 H = 3
          PARP(2) = 2 ! Modify because of the FM
!c Scale for GVMD process D = 0.4 H = 0.17
          PARP(18) = 0.17
!c effective Q or transverse k cut off for parton shower D = 1 H = 0.5
          PARP(62) = 0.5
!c cut off for energy (cm) D = 2 H = 0.5
          PARP(65) = 0.5
!c modify Q2 scale D = 4 H = 1
          PARP(67) = 1.
!c width of primordial gaussian inside hadron D = 2 H = 0.44
          PARP(91) = 0.44
!c upper cut off for kt inside hadron D = 5 H = 2
          PARP(93) = 2.
!c With of primordial kt in photon D = 1 H = .44
          PARP(99) = 0.44
!c upper cut off for kt inside photon D = 5 H = 2
          PARP(100) = 2
!c Mass spectrum of diffractive state D = 0.28 H = 0.5
          PARP(102) = 0.5
!c Mass cut for isotropic decay of diffractive state D = 1 H = 0.5
          PARP(103) = 0.5
!c Energy cut off for hadron hadron D = 0.8 H = 0.3
!c manual say it cannot go below .8
          PARP(104) = 0.3
!c Minimum invariant mass of remnant D = 2 H = 0
          PARP(111) = 0.
!c Some coupling constants for vector mesons
          PARP(161) = 2.69  ! D = 2.20
          PARP(162) = 24.6  ! D = 23.6
          PARP(163) = 18.8  ! D = 18.4
          PARP(164) = 11.5  ! D = 11.5
!c simple factor to supress some transverse resolved photons D = 0.5 H = 0.33
          PARP(165) = 0.33

!c Suppression factor for some quark processes
          PARJ(1) = 0.025 ! D = 0.1
          PARJ(2) = 0.120 ! D = 0.3
          PARJ(3) = 0.25  ! D = 0.4
!c Factors to determine spin of mesons
!c Probability for spin 1 for light mesons
          PARJ(11) = 0.25 ! D =0.5
!c Probability for spin 1 for strange mesons
          PARJ(12) = 0.3  !D = 0.6
!c Width of gaussian for transverse momentum of primary hadrons D = 0.36 H = 0.63
          PARJ(21) = 0.63
!c Make a non gaussian tail to the transverse momentum (HUGE MODIF)
!c D = 0.01 H = 0.3
          PARJ(23) = 0.3
!c D = 2 H = 5
          PARJ(24) = 5.

!c remaining energy below witch the fragmentation is stopped D = 0.8
          PARJ(33) = 0.6
!c Parameters for Lund fragmentation
          PARJ(41) = 1.13 ! D = 0.3
          PARJ(42) = 0.37 ! D = 0.58

!c Parameter for lund fragmentation D = 0.5
          PARJ(45) = 0.8

!c Choice of fragmentation scheme D = 1 H = 1
!c      MSTJ(1) = 1
!c Choice of gluon fragmentation scheme D = 3 H = 3
!c      MSTJ(2) = 3
!c Choice of Baryon production model D = 2 H = 1
          MSTJ(12) = 1 !Here I made most important change 1->2 ??
!c Maximum flavour produce by gluon fragmentation D = 5 H = 4
          MSTJ(45) = 4
!c Number min or max of flavor for different process
          MSTU(112) = 4 !D = 5
          MSTU(113) = 4 !D = 5
          MSTU(114) = 4 !D = 5

!c Lower limit for sqrt(s)
          CKIN(1) = 1.
!c Max Pt for hard process
          CKIN(4) = 2.
!c Bunch of range of pseudo rapidity
          CKIN(9) = -10.
          CKIN(10) = 10.
          CKIN(11) = -10.
          CKIN(12) = 10.
          CKIN(13) = -10.
          CKIN(14) = 10.
          CKIN(15) = -10.
          CKIN(16) = 10.

!c Q2 Limits
          CKIN(65) = .9
          CKIN(66) = 4.
!c W limit
          CKIN(77) = 1.90
          CKIN(78) = -1.

      end

!c------------------------------------------------------------------------------
!c Initialize the config values of LEPTO
!c------------------------------------------------------------------------------
      subroutine PythiaConfigBrahim
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

!ccccc Include all the common blocks
          ! include 'common.f90'

          MSEL = 2
!c      MSUB(29) = 1
!c
!c      MSEL = 0
!c      MSUB(99)=1

!C ... Set kinematic cuts
          CKIN(1) = 1.d0             ! min invariant mass = sqrt(s)
!c      CKIN(21) = .5d0                      ! min x1
!c      CKIN(23) = .5d0                      ! min x2
!c      CKIN(35) = .1d0                      ! min t
!c      CKIN(61) = .5d0                      ! min x=nu/E0 (beam)
!c      CKIN(63) = .5d0                      ! min x=..    (target)
          CKIN(65) = 1.d0             ! min Q2
          CKIN(77) = 2.d0             ! min W

!c ... Keep only photon exchange
          MSTP(1) = 2
          MSTP(13) = 2
!c      MSTP(14) = 30
!c      MSTP(15) = 0
!c      MSTP(16) = 1
          MSTP(17) = 3
          MSTP(18) = 2
!c      MSTP(19) = 4
          MSTP(20) = 0
          MSTP(21) = 2                      ! neutral g* only
!c      MSTP(22)  = 0                      ! Q2 definition override
!c      MSTP(23)  = 1                      ! retain x and Q2 of original scat.
!c      MSTP(32) = 8                       ! Q2 definition
          MSTP(38) = 4                       ! quark loop masses
!c
          MSTP(51) = 11
          MSTP(52) = 1
!c
          MSTP(58) = 4
!c
          MSTP(61) = 0
          MSTP(62) = 3
          MSTP(63) = 2
          MSTP(64) = 2
          MSTP(65) = 1
          MSTP(66) = 5
          MSTP(67) = 2
          MSTP(68) = 1
          MSTP(69) = 0
!c
          MSTP(71) = 0
!c
          MSTP(81) = 0
          MSTP(82) = 1
          MSTP(83) = 100
          MSTP(86) = 2
!c
          MSTP(91) = 1
          MSTP(92) = 4
          MSTP(93) = 1
          MSTP(94) = 3
!c
          MSTP(101) = 1

!c ... Allow low c.m. energies
          PARP(2) = 1.d0    ! min(Ecm) = 3 GeV, default 10 GeV
!c      PARP(2)   = 1.d0    ! min(Ecm) = 3 GeV, default 10 GeV
!c      PARP(18)  = .17d0    ! min(Ecm) = 3 GeV, default 10 GeV

!c      PARP(42)  = 0.d0    ! reduce the min. mass of Z/W

!c      PARP(62)  = .5d0    ! reduce the min. mass of time-like parton
!c      PARP(65)  = .5d0    ! reduce the min. mass of time-like parton
!c      PARP(67)  = 1.d0    ! reduce the min. mass of time-like parton

!c      PARP(91)  = 0.44d0  ! reduce the min. mass of time-like parton
!c      PARP(93)  = 2.d0    ! reduce the min. mass of time-like parton
!c      PARP(94)  = 1.d0    ! reduce the min. mass of time-like parton

!c      PARP(99)  = 0.44d0  ! reduce the min. mass of time-like parton
!c      PARP(100) = 2

!c      PARP(102) = .6d0
!c      PARP(103) = .5d0

!c      PARP(111) = 0.d0    ! reduce the min remnant invariant mass, default 2 GeV

!c      PARP(165) = .33d0  ! reduce the min. mass of time-like parton

      end

!c------------------------------------------------------------------------------
!c Initialize the config values of LEPTO
!c------------------------------------------------------------------------------
      subroutine PythiaConfigHERMES
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

!ccccc Include all the common blocks
          ! include 'common.f90'

!C****************************************************************************
!C Latest pythia tune with fragmentation tune 2004_C
!C*****************************************************************************
!C NEW: neglect error messages
          MSTU(21) = 1
!C
!C NON_DEFAULT VALUES OR VALUES WHICH CAN BE CHANGED FOR STUDIES
!C*****************************************************************************
!C Selects the typ of processes Pythia uses to generate events
!C MSEL 1)=2
!C
!C Q2 range over which electrons are assumed to radiate photons Pythia6
!C default 1 MSTP 13)=2 better for photo-production
          MSTP(13) = 2
!C
!C photon structure choice
          MSTP(14) = 30
!c DIS only
          MSTP(14) = 26
!C
!C MSTP(15) regulates the pt_min treatment for anaomolous in respect to VMD
!C Default MSTP 15)=0
          MSTP(15) = 0
!C
!C Variable to tell pythia in which variable the generation happens
!C MSTP(16)=1 --> generation in y
!C MSTP(16)=0 --> generation in pythia-x
          MSTP(16) = 1
!C
!C MSTP(17)=6 is the R-rho measured as by hermes MSTP 17=4 is default
!C check also PARP165/PARP166 if MSTP 17)=6
          MSTP(17) = 6
!C
!C scale for pt_min VMD/GVMD  MSTP 18)=3 by W^2
!C scale for pt_min VMD+GVMD to direct processes MSTP 18)=2  pt_min = PARP(15)
          MSTP(18) = 3
!C
!C choice of partonic cross section in process 99 and dampening factor
          MSTP(19) = 4
!C
!C An additional suppression of resolved (VMD or GVMD) c.s.
!C (W^2/(W^2 + Q_1^2 + Q_2^2))^MSTP(20) ---> MSTP 20)=3
!C if MSTP 20)=0 than the q^2 -slope in pyxtot.F is changed to 2.575
          MSTP(20) = 4
!C
!C Q^2 definition in hard scattering for 2 --> 2 processes
!C Q2 )= pT**2 + (P1**2 + P2**2 +m3**2 + m4**2)/2.
          MSTP(32) = 8
!C
!C handling of masses in quark loops
          MSTP(38) = 4
!C
!C proton parton distribution: genSet_PartonSet
!c     MSTP(51)=4046
!c     MSTP(52)=2
!C MSTP(51)=11
!C MSTP(52)=1
!C
!C pion parton distribution set
          MSTP(53) = 3
!C MSTP54)=1 internal ones from pythia MSTP54=2 pfdf-lib
!c     MSTP(54)=1
!C
!C photon parton distributions default MSTP 55)=5 SaS1D
!C MSTP(55)=5003  GRV
!c     MSTP(55)=5
!C choice of photon parton distributions according to pythia intern MSTP 56)=1
!C MSTP(56)=2  use PDFLIB
!c     MSTP(56)=1
!C
!C Q2 dependence in PDFs default Pythia6 1, but 0 means PDFs are made Q2
!C independent and used at lower cut-off value Q^2_o
          MSTP(57) = 1
!C max # of quark flavours in PDFs (GMC: 58)=genSet_PyMaxFl)
          MSTP(58) = 4
!C
!C extension of the SaS real-photon distributions to off-shell photons,
!C especially for the anomalous component
          MSTP(60) = 7
!C
!C master switch QCD and QED ISR (D: 1)
          MSTP(61) = 0
!C master switch QCD and QED FSR (D: 1)
          MSTP(71) = 0
!C
!C master switch for multiple interactions
!C Default: 81)=1, 82=1
          MSTP(81) = 0
          MSTP(82) = 1
!C
!C shape of primordial k_t in hadron
          MSTP(91) = 1
!C energy sharing between two coloured beam remnants
          MSTP(92) = 4
!C shape of primordial k_t in photon
          MSTP(93) = 1
!C
!C structure of diffrative system Default:3
          MSTP(101) = 1
!C rho^0 decays according to angular distribution
          MSTP(102) = 1
!C Switch to turn off fragmentation completely than MSTP(111)=0
          MSTP(111) = 1
!C
!C Initialization of maxima,  Default MSTP(121))=0, MST(121)=1: multiply by PARP(121)
          MSTP(121) = 1
!C
!C************************************************************PARP settings
!C min CMS energy allowed for event as a whole (set for gamma-p automatically)
          PARP(2) = 7.
!C
!C maximum scale for photoproduction if using MSTP(13))=2
          PARP(13) = 1
!C
!C Suppression factor for GVMD compared to VMD
          PARP(18) = 0.17
!C
!C MSTP(18))=3 pTmin used there is parameterized as
!C PARP(81)*(W/PARP(89))**PARP(90)
!C default: PARP(81))=1.9GeV, PARP(89)=1000GeV, PARP(90)=0.16
          PARP(81) = 1.9
          PARP(89) = 1000
          PARP(90) = 0.16
!C
!C intrinsic kT of initial state partons in hadron default pythia6 91)=1.
          PARP(91) = 0.40
!C upper cut on primordial kT spectrum default pythia6 93)=5 pythi5 93=2
          PARP(93) = 2.
!C
!C DEFAULT in Pythia5 PARP 99 )= 0.44 in Pythia6 99=1
          PARP(99) = 0.40
!C upper cut on primordial kT photon spectrum default pythia6 100)=5
!C pythia5 100)=2
          PARP(100) = 5
!C
!C mass above the VM (rho, omega, phi) mass, where the single- and double
!C diffractive spectrum starts (default: 0.28)
          PARP(102) = 0.5
! C mass above the VM (rho, omega, phi) mass, where the fragmentation of
! C diffractive states starts (below a two-body decay is forced (default: 1.)
          PARP(103) = 0.5
!C
!C minimum energy above threshold, where hadronic cross sections are defined
          PARP(104) = 0.3
!C
!C hadronic beam remnant has an energy of at least PARP(111) in the rest frame
!C of the event default pythia5 111)=2
          PARP(111) = 0.
!C
!C Factor to multiply cross section maxima with
          PARP(121) = 2.
!C
!C PARP 161-164 are the coupling constant photon-VM for VMD model
!C (default: 2.2 (rho), 23.6(omega), 18.4(phi), 11.5(J/psi))
          PARP(161) = 3.00
          PARP(162) = 24.6
          PARP(163) = 18.8
          PARP(164) = 11.5
!C
!C PARP 165/166 are linked to MSTP 17 as R_rho of HERMES is used
!C Fit to world data
          PARP(165) = 0.47679
          PARP(166) = 0.67597
!C values from the old He-3 paper for 2< W < 7GeV
!C PARP(165)=0.33
!C PARP(166)=0.61
!C values for only longitudinal rho's
!C PARP(165)=1000.
!C PARP(166)=1.
!C values for only transverse rho's
!C PARP(165)=0.
!C PARP(166)=0.67597
!C
!C******************************************************** PYTHIA/JETSET control switches
!C
!C  diquark suppression P(qq)/P(q)  default:0.1
          PARJ(1) = 0.029
!C
!C s quark suppression P(s)/P(u) default:0.3
          PARJ(2) = 0.283
!C
!C extra suppression of strange diquarks default:0.4 tuned by P. Kravtsov
          PARJ(3) = 1.2
!C
!C suppression of spin-1 diquarks over spin-0 default:0.05
          PARJ(4) = 0.05
!C
!C suppression of BMBbar to BBbar in popcorn model default: 0.5/0.5/0.5
          PARJ(5) = 0.5
          PARJ(6) = 0.5
          PARJ(7) = 0.5
!C
!C Vectormeson to pseudoscaler suppression default:0.5
          PARJ(11) = 0.5
!C
!C Vectormeson to pseudoscaler suppression for strangeness default:0.6
          PARJ(12) = 0.6
!C
!C width for px, py transverse momentum distributions for primary hadrons default:0.36
          PARJ(21) = 0.400
!C
!C fraction of non gaussian tails to the pt distribution times a factor to increase PARJ 21
!C default:0.01 / 2
          PARJ(23) = 0.01
          PARJ(24) = 2.0
!C
!C minimum allowable energy for color-singlet jet system default:1.0
          PARJ(32) = 1.0
!C
!C (* D: 0.8) minimum energy used to stop fragmenting a jet default:0.8
          PARJ(33) = 0.800
!C
!C a parameter for symmetric Lund fragmentation function default:0.3
          PARJ(41) = 1.94
!C
!C b parameter for symmetric Lund fragmentation function default:0.58
          PARJ(42) = 0.544
!C
!C a parameter for the symmetric Lund fragmentation function for diquarks
!C default:0.5
          PARJ(45) = 1.05
!C
!C Select the fragmentation modell 0: no fragmentation, 1: string, 2: independent
          MSTJ(1) = 1
!C
!C choice of the baryon production model: MSTJ 12)=2 Popcorn scheme
          MSTJ(12) = 1
!C
!C parton showering is turned of
!C MSTJ(41)=0
!C
!C maximum flavour that can be produced in shower by g->qq
          MSTJ(45) = 4
!C
!C nominal number of flavours assumed in alpha_s expression
          MSTU(112) = 4
!C minimum number of flavours that may be assumed in alpha_s expression
          MSTU(113) = 4
!C maximum number of flavours that may be assumed in alpha_s expression
          MSTU(114) = 4
!C
!C******************************************************** kinematical cuts
!C
!C range of m^ )= sqrt(s^)
!C (D: 2., -1.)
          CKIN(1) = 1.
          CKIN(2) = -1.
!C range of p^_t (DIS )= Q^2)
!C (*) (D: 0., -1.) (GMC: 0, -1)
          CKIN(3) = 0.
!C CKIN( 3)=1.
          CKIN(4) = -1.
!C p^_t,min for singular processes in limit p^_t -> 0
          CKIN(5) = 1.00
!C m_0 hard 2->2 process is classified singular for p^_t -> 0 if the mass
!C of at least one of the two outgoing partons is below m_0
          CKIN(6) = 1.00
!C range for W^2 in DIS processes (W^2 )= Q^2(1-x)/x, neglecting M^2 and ISR!)
!C (*) (D: 4.,-1) before (10. -1)
          CKIN(39) = 4.
          CKIN(40) = -1.
!C range for Q^2
          CKIN(65) = 1.
          CKIN(66) = 100.
!C CKIN(65)=0.1
!C CKIN(66)=0.01
!C allowed range for W, i.e. either the photon-hadron or photon-photon invariant mass
          CKIN(77) = 2.0
          CKIN(78) = -1.

      end

