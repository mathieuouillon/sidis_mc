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

      MSTP(1) = 2 
      MSTP(2) = 1 
      MSTP(3) = 2 
      MSTP(4) = 0 
      MSTP(5) = 0 
      MSTP(7) = 0 
      MSTP(8) = 0 
      MSTP(9) = 0 
      MSTP(11) = 1 
      MSTP(12) = 0 
      MSTP(13) = 2 
      MSTP(14) = 30 
      MSTP(15) = 0 
      MSTP(16) = 1 
      MSTP(17) = 6 
      MSTP(18) = 3 
      MSTP(19) = 4 
      MSTP(20) = 0 
      MSTP(21) = 1 
      MSTP(22) = 0 
      MSTP(23) = 1 
      MSTP(31) = 1 
      MSTP(32) = 8 
      MSTP(33) = 0 
      MSTP(34) = 1 
      MSTP(35) = 0 
      MSTP(36) = 2 
      MSTP(37) = 1 
      MSTP(38) = 4 
      MSTP(39) = 2 
      MSTP(40) = 0 
      MSTP(41) = 1 
      MSTP(42) = 1 
      MSTP(43) = 3 
      MSTP(44) = 7 
      MSTP(45) = 3 
      MSTP(46) = 1 
      MSTP(47) = 1 
      MSTP(48) = 0 
      MSTP(49) = 1 
      MSTP(50) = 0 
      MSTP(51) = 7 
      MSTP(52) = 1 
      MSTP(53) = 3 
      MSTP(54) = 1 
      MSTP(55) = 5 
      MSTP(56) = 1 
      MSTP(57) = 1 
      MSTP(58) = 4 
      MSTP(59) = 1 
      MSTP(60) = 7 
      MSTP(61) = 0 
      MSTP(61) = 0 
      MSTP(62) = 3 
      MSTP(63) = 2 
      MSTP(64) = 2 
      MSTP(65) = 1 
      MSTP(66) = 5 
      MSTP(67) = 2 
      MSTP(68) = 1 
      MSTP(69) = 0 
      MSTP(71) = 0 
      MSTP(81) = 0 
      MSTP(82) = 1 
      MSTP(83) = 100 
      MSTP(86) = 2 
      MSTP(91) = 1 
      MSTP(92) = 4 
      MSTP(93) = 1 
      MSTP(94) = 3 
      MSTP(101) = 1 
      MSTP(102) = 1 
      MSTP(111) = 1 
      MSTP(121) = 1 
      MSTP(131) = 0 
      MSTP(171) = 0 
      MSTP(172) = 2 
      MSTP(173) = 0 

      PARP(1) = 0.25 
      PARP(2) = 3 
      PARP(13) = 1 
      PARP(14) = 0.01 
      PARP(15) = 0.5 
      PARP(16) = 1 
      PARP(17) = 1 
      PARP(18) = 0.17 
      PARP(61) = 0.25 
      PARP(62) = 0.5 
      PARP(63) = 0.25 
      PARP(64) = 1 
      PARP(65) = 0.5 
      PARP(66) = 0.001 
      PARP(67) = 1. 
      PARP(68) = 0.001 
      PARP(71) = 4. 
      PARP(72) = 0.25 
      PARP(91) = 0.44 
      PARP(93) = 2. 
      PARP(94) = 1. 
      PARP(95) = 0. 
      PARP(96) = 3. 
      PARP(97) = 1. 
      PARP(98) = 0.75 
      PARP(99) = 0.44 
      PARP(100) = 2 
      PARP(102) = 0.5 
      PARP(103) = 0.5 
      PARP(104) = 0.3 
      PARP(111) = 0. 
      PARP(161) = 2.69 
      PARP(162) = 24.6 
      PARP(163) = 18.8 
      PARP(164) = 11.5 
      PARP(165) = 0.33 
      PARJ(1) = 0.025 
      PARJ(2) = 0.120 
      PARJ(3) = 0.25 
      PARJ(4) = 0.05 
      PARJ(5) = 0.5 
      PARJ(6) = 0.5 
      PARJ(7) = 0.5 
      PARJ(11) = 0.25 
      PARJ(12) = 0.3 
      PARJ(21) = 0.63 
      PARJ(23) = 0.3 
      PARJ(24) = 5. 
      PARJ(32) = 1.0 
      PARJ(33) = 0.6 
      PARJ(41) = 1.13 
      PARJ(42) = 0.37 
      PARJ(45) = 0.8 

      MSTJ(1) = 1 
      MSTJ(2) = 3 
      MSTJ(3) = 0 
      MSTJ(12) = 1 !Here I made most important change 1->2
      MSTJ(40) = 0 
      MSTJ(45) = 4 
      MSTU(112) = 4 
      MSTU(113) = 4 
      MSTU(114) = 4 

      CKIN(1) = 1. 
      CKIN(2) = -1. 
      CKIN(3) = 0. 
      CKIN(4) = 2. 
      CKIN(5) = 1.00 
      CKIN(6) = 1.00 
      CKIN(7) = -10. 
      CKIN(8) = 10. 
      CKIN(9) = -10. 
      CKIN(10) = 10. 
      CKIN(11) = -10. 
      CKIN(12) = 10. 
      CKIN(13) = -10. 
      CKIN(14) = 10. 
      CKIN(15) = -10. 
      CKIN(16) = 10. 
      CKIN(17) = -1. 
      CKIN(18) = 1. 
      CKIN(19) = -1. 
      CKIN(20) = 1. 
      CKIN(21) = 0. 
      CKIN(22) = 1. 
      CKIN(23) = 0. 
      CKIN(24) = 1. 
      CKIN(25) = -1. 
      CKIN(26) = 1. 
      CKIN(27) = -1. 
      CKIN(28) = 1. 
      CKIN(31) = 2. 
      CKIN(32) = -1. 
      CKIN(35) = 0. 
      CKIN(36) = -1 
      CKIN(37) = 0. 
      CKIN(38) = -1. 
      CKIN(39) = 4. 
      CKIN(40) = -1. 
      CKIN(65) = 1. 
      CKIN(66) = 4. 
      CKIN(77) = 2.0 
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
      PARP(2)   = 3.d0    ! min(Ecm) = 3 GeV, default 10 GeV
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




