c------------------------------------------------------------------------------
c Initialize the config values of LEPTO
c------------------------------------------------------------------------------
      subroutine PythiaConfigCLAS
      implicit none

ccccc Include all the common blocks
      include 'common.f'

c Kind of possible multiple interaction (needed to avoid bugs) D = 4 H =1
      MSTP(82) = 1 

c Lowest CM energy D = 10 H = 3
      PARP(2) = 2 ! Modify because of the FM

c remaining energy below witch the fragmentation is stopped D = 0.8
      PARJ(33) = 0.3 ! This parameter have huge effect on z distribution

c Lower limit for sqrt(s)
      CKIN(1) = 1. 

c Q2 Limits
      CKIN(65) = .9 
      CKIN(66) = 4. 
c W limit
      CKIN(77) = 1.90 
      CKIN(78) = -1. 


      end

c------------------------------------------------------------------------------
c Initialize the config values of LEPTO
c------------------------------------------------------------------------------
      subroutine PythiaConfigOWN
      implicit none

ccccc Include all the common blocks
      include 'common.f'

ccc Deeply Inelastic Scattering and γ ∗ γ ∗ physics
ccc MSEL = 1, 2, 35, 36, 37, 38
c MSEL = 0 Tout les process doivent etre activer individuellement
c MSEL = 1 in our case call MSTP(14)
c MSEL = 2 include 1 + elastic? + diffractive? + low Pt process
      MSEL = 2

ccc MSUB(ISUB) = 0 : the subprocess is excluded.
ccc MSUB(ISUB) = 1 : the subprocess is included.

      
C ... Set kinematic cuts
c min invariant mass = sqrt(s)
      CKIN(1)  = 1.d0             ! lower lim 
c      CKIN(2)  = -1.d0            ! no higher limit
c range of Pt
c      CKIN(3)  = 0.d0  
c      CKIN(4)  = 2.d0  
c      CKIN(21) = .5d0                      ! min x1
c      CKIN(23) = .5d0                      ! min x2
c      CKIN(35) = .1d0                      ! min t
c      CKIN(61) = .5d0                      ! min x=nu/E0 (beam)
c      CKIN(63) = .5d0                      ! min x=..    (target)
      CKIN(65) = 0.85d0           ! min Q2
c      CKIN(74) = .9d0             ! max y
      CKIN(77) = 1.8d0            ! min W

c Suppression of VMD (0:none)
      MSTP(20) = 0
c Q2 definition override
      MSTP(22)  = 4                      
c retain x and Q2 of original scat.
      MSTP(23)  = 0                      
c Nb of massless quark
      MSTP(38) = 3
c Master swich for decays
c      MSTP(41)  = 1
c Choice of the proton parton distrib.
      MSTP(51)  = 6
c Number of quark flavour in the PDF 
      MSTP(58)  = 3
c QED and QCD  radiative effects
      MSTP(61)  = 2 
c Master switch for multiple interaction
      MSTP(81)  = 0
c structure of multiple interaction
      MSTP(82)  = 1
c structure of diffractive system      
      MSTP(101)  = 1

c ... Allow low c.m. energies
      PARP(2)   = 2.d0    ! min(Ecm) = 2 GeV, default 10 GeV
      PARP(62)  = .3d0    ! reduce the min. mass of time-like parton
      PARP(65)  = .3d0    ! reduce the min. mass of time-like parton

c Allow diffractive excitement above M+PARP(102)
      PARP(102) = .4d0   
c Minimum energy for elastic and difractive effects
      PARP(104) = .3d0   

      PARP(111) = 0.d0    ! reduce the min remnant invariant mass, default 2 GeV 

c Minor change about vector meson production
      PARP(161) = 2.69d0 
      PARP(162) = 24.6d0 
      PARP(163) = 18.8d0 
      PARP(165) = .33d0  ! reduce the min. mass of time-like parton

c Stuff about the number of flavour
      MSTJ(45) = 3
      MSTU(112) = 3
      MSTU(113) = 2
      MSTU(114) = 3


      end

c------------------------------------------------------------------------------
c Initialize the config values of LEPTO
c------------------------------------------------------------------------------
      subroutine PythiaConfigHayk
      implicit none

ccccc Include all the common blocks
      include 'common.f'

c Maximum number of generations D = 3 / H = 2
      MSTP(1) = 2 
c Calculation of alpha s D = H = 1
c Comments in thesis let think he want to put 0...
c      MSTP(2) = 1 
c Not in Pythia manual ...
       MSTP(5) = 0 
c Fix the range over wich electron emit photons D = 1 H = 2
c text in manual seem to push for 1 
      MSTP(13) = 2 
c structure of incoming photon beam D = H = 30
c      MSTP(14) = 30 
c Choice of definition of the fractional part taken by photon D = H = 1
c lets try 0
c      MSTP(16) = 1 
c Possibility of a extra factor for processes involving resolved virtual photons
c D = 4 H = 6   6 is not in the manual !
      MSTP(17) = 6 
c suppression of resolved (VMD or GVMD) cross sections D = 3 H = 0
      MSTP(20) = 0 
c Keep the final electron giving the right x and q2 D = H = 1
c      MSTP(23) = 1 
c Handling of quark loops, number of allowed quarks D = 5 H = 4
      MSTP(38) = 4 
c master switch for decay  D = 2 H = 1
      MSTP(41) = 1 
c Choice of the PDF D = H = 7
c      MSTP(51) = 7 
c Max number of quarks in pdf D = 5 H = 4
      MSTP(58) = 4 

c General switch for initial state radiations D = 2 H = 0
c try to turn on 1 or 2
      MSTP(61) = 0 
c normaly no effect from 6x processes
      MSTP(62) = 3 
      MSTP(63) = 2 
      MSTP(64) = 2 
      MSTP(65) = 1 
      MSTP(66) = 5 
      MSTP(67) = 2 
      MSTP(68) = 1 
      MSTP(69) = 0 
c Master switch for final state radiations D = 1 H = 0
      MSTP(71) = 0 
c Master switch for multiple interactions D = 1 H = 0
      MSTP(81) = 0 
c normaly no effect from 8x processes
      MSTP(82) = 1 
      MSTP(83) = 100 
      MSTP(86) = 2 
c 
c Change the energy partitionning between remnant D = 3 H = 4
      MSTP(92) = 4 
c Structure of diffractive system D = 3 H = 1
      MSTP(101) = 1 
c Calculation of kinematic coefficients D = 0 H = 1
      MSTP(121) = 1 


c Lowest CM energy D = 10 H = 3
      PARP(2) = 2 ! Modify because of the FM
c Scale for GVMD process D = 0.4 H = 0.17
      PARP(18) = 0.17 
c effective Q or transverse k cut off for parton shower D = 1 H = 0.5
      PARP(62) = 0.5 
c cut off for energy (cm) D = 2 H = 0.5 
      PARP(65) = 0.5 
c modify Q2 scale D = 4 H = 1
      PARP(67) = 1. 
c width of primordial gaussian inside hadron D = 2 H = 0.44
      PARP(91) = 0.44 
c upper cut off for kt inside hadron D = 5 H = 2
      PARP(93) = 2. 
c With of primordial kt in photon D = 1 H = .44
      PARP(99) = 0.44 
c upper cut off for kt inside photon D = 5 H = 2
      PARP(100) = 2 
c Mass spectrum of diffractive state D = 0.28 H = 0.5
      PARP(102) = 0.5 
c Mass cut for isotropic decay of diffractive state D = 1 H = 0.5
      PARP(103) = 0.5 
c Energy cut off for hadron hadron D = 0.8 H = 0.3
c manual say it cannot go below .8
      PARP(104) = 0.3 
c Minimum invariant mass of remnant D = 2 H = 0
      PARP(111) = 0. 
c Some coupling constants for vector mesons
      PARP(161) = 2.69  ! D = 2.20
      PARP(162) = 24.6  ! D = 23.6
      PARP(163) = 18.8  ! D = 18.4
      PARP(164) = 11.5  ! D = 11.5
c simple factor to supress some transverse resolved photons D = 0.5 H = 0.33
      PARP(165) = 0.33 



c Suppression factor for some quark processes
      PARJ(1) = 0.025 ! D = 0.1
      PARJ(2) = 0.120 ! D = 0.3
      PARJ(3) = 0.25  ! D = 0.4
c Factors to determine spin of mesons
c Probability for spin 1 for light mesons
      PARJ(11) = 0.25 ! D =0.5
c Probability for spin 1 for strange mesons
      PARJ(12) = 0.3  !D = 0.6
c Width of gaussian for transverse momentum of primary hadrons D = 0.36 H = 0.63
      PARJ(21) = 0.63 
c Make a non gaussian tail to the transverse momentum (HUGE MODIF)
c D = 0.01 H = 0.3
      PARJ(23) = 0.3 
c D = 2 H = 5
      PARJ(24) = 5. 

c remaining energy below witch the fragmentation is stopped D = 0.8
      PARJ(33) = 0.6 
c Parameters for Lund fragmentation
      PARJ(41) = 1.13 ! D = 0.3
      PARJ(42) = 0.37 ! D = 0.58

c Parameter for lund fragmentation D = 0.5
      PARJ(45) = 0.8 

c Choice of fragmentation scheme D = 1 H = 1 
c      MSTJ(1) = 1 
c Choice of gluon fragmentation scheme D = 3 H = 3 
c      MSTJ(2) = 3 
c Choice of Baryon production model D = 2 H = 1
      MSTJ(12) = 1 !Here I made most important change 1->2 ??
c Maximum flavour produce by gluon fragmentation D = 5 H = 4
      MSTJ(45) = 4 
c Number min or max of flavor for different process
      MSTU(112) = 4 !D = 5
      MSTU(113) = 4 !D = 5 
      MSTU(114) = 4 !D = 5 

c Lower limit for sqrt(s)
      CKIN(1) = 1. 
c Max Pt for hard process
      CKIN(4) = 2. 
c Bunch of range of pseudo rapidity
      CKIN(9) = -10. 
      CKIN(10) = 10. 
      CKIN(11) = -10. 
      CKIN(12) = 10. 
      CKIN(13) = -10. 
      CKIN(14) = 10. 
      CKIN(15) = -10. 
      CKIN(16) = 10. 


c Q2 Limits
      CKIN(65) = .9 
      CKIN(66) = 4. 
c W limit
      CKIN(77) = 1.90 
      CKIN(78) = -1. 


      end

c------------------------------------------------------------------------------
c Initialize the config values of LEPTO
c------------------------------------------------------------------------------
      subroutine PythiaConfigBrahim
      implicit none

ccccc Include all the common blocks
      include 'common.f'


      MSEL = 2
c      MSUB(29) = 1
c
c      MSEL = 0
c      MSUB(99)=1    
      
C ... Set kinematic cuts
      CKIN(1)  = 1.d0             ! min invariant mass = sqrt(s)
c      CKIN(21) = .5d0                      ! min x1
c      CKIN(23) = .5d0                      ! min x2
c      CKIN(35) = .1d0                      ! min t
c      CKIN(61) = .5d0                      ! min x=nu/E0 (beam)
c      CKIN(63) = .5d0                      ! min x=..    (target)
      CKIN(65) = 1.d0             ! min Q2
      CKIN(77) = 2.d0             ! min W

c ... Keep only photon exchange
      MSTP(1)  = 2
      MSTP(13) = 2
c      MSTP(14) = 30
c      MSTP(15) = 0
c      MSTP(16) = 1
      MSTP(17) = 3
      MSTP(18) = 2
c      MSTP(19) = 4
      MSTP(20) = 0
      MSTP(21)  = 2                      ! neutral g* only
c      MSTP(22)  = 0                      ! Q2 definition override
c      MSTP(23)  = 1                      ! retain x and Q2 of original scat.
c      MSTP(32) = 8                       ! Q2 definition
      MSTP(38) = 4                       ! quark loop masses
c
      MSTP(51)  = 11
      MSTP(52)  = 1
c      
      MSTP(58)  = 4
c
      MSTP(61)  = 0
      MSTP(62)  = 3
      MSTP(63)  = 2
      MSTP(64)  = 2
      MSTP(65)  = 1
      MSTP(66)  = 5
      MSTP(67)  = 2
      MSTP(68)  = 1
      MSTP(69)  = 0
c      
      MSTP(71)  = 0
c
      MSTP(81)  = 0
      MSTP(82)  = 1
      MSTP(83)  = 100
      MSTP(86)  = 2
c
      MSTP(91)  = 1
      MSTP(92)  = 4
      MSTP(93)  = 1
      MSTP(94)  = 3
c      
      MSTP(101)  = 1

c ... Allow low c.m. energies
      PARP(2)   = 1.d0    ! min(Ecm) = 3 GeV, default 10 GeV
c      PARP(2)   = 1.d0    ! min(Ecm) = 3 GeV, default 10 GeV
c      PARP(18)  = .17d0    ! min(Ecm) = 3 GeV, default 10 GeV


c      PARP(42)  = 0.d0    ! reduce the min. mass of Z/W
            
c      PARP(62)  = .5d0    ! reduce the min. mass of time-like parton
c      PARP(65)  = .5d0    ! reduce the min. mass of time-like parton
c      PARP(67)  = 1.d0    ! reduce the min. mass of time-like parton

c      PARP(91)  = 0.44d0  ! reduce the min. mass of time-like parton
c      PARP(93)  = 2.d0    ! reduce the min. mass of time-like parton
c      PARP(94)  = 1.d0    ! reduce the min. mass of time-like parton

c      PARP(99)  = 0.44d0  ! reduce the min. mass of time-like parton
c      PARP(100) = 2    

c      PARP(102) = .6d0   
c      PARP(103) = .5d0   

c      PARP(111) = 0.d0    ! reduce the min remnant invariant mass, default 2 GeV 

c      PARP(165) = .33d0  ! reduce the min. mass of time-like parton


      end




