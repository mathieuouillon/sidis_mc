      integer function simulation()
      implicit none
c------------------------------------------------------------------------------
c TO DO LIST:
c   Implement some radiative effect
c   Have alpha s changing with Q2 in QW
c------------------------------------------------------------------------------

ccccc Include all the common blocks
      include 'common.f'
c     include 'include/names.inc'

ccccc Miscellaneous
      integer icycle ! for hbook
      real T1,T2 ! For time of computation
      double precision BeamE !Input value for pythia
      integer ip ! For do
      integer i
 

      call TIMEX(T1)
CCCCCC Begining of the simulation
ccc Integer j,nkin ! number of kinematics (??? is it still the case?)
      nkin = 10000
ccc Number of events per kinematics 
      nevent = 10000000
ccc Electron energy (GeV)
      E0 = 27.5
ccc Target type ! 0-> p, 1-> 2H, 2-> 3H, 3-> 3He, 4-> 4He, 5-> 6Li, 
ccc               6-> 7Li, 7-> C, 8-> Al, 9-> Fe, 10-> Sn, 11-> Pb
ccc               12-> Ne, 13-> Kr, 14-> Xe
      iTg = 14

ccc Collider options
c     integer iColl !1 = activate collider kinematic
      iColl = 0
c     real EColl ! energy of the nuclei (GeV/nucleon)
      EColl = 38.

ccc Fermimotion flag 
c     0 = no FM
c     1 = like 4 plus a tail from [2] (deut, C, Al, Fe, Sn, Pb, 4He)
c     2 = Accardi SVG (7Li, C, O, Ne, Al, Ar, Ca, Ni, Cu, Zr, Sn, Pb)
c     3 = Accardi CS (2H, 3He, 4He, C, O, Ca, Fe, Pb)
c     4 = hard sphere with values from [1], 
c     5 = R. Wiringa private communication
c     All FM distributions are limited to 1 GeV nucleons
c     [1] E. J. Moniz et al. PRL 26, 445 (1971)
c     [2] A. Bodek and J. L. Ritchie PRD 23, 1070 (1981)
      iFM = 3
      FMlimit = 1

ccc Init for the quenching weights
c     integer iQuenching ! 0 desactivate Quenching
      iQuenching = 1
c     integer iqw 1 SW, 2 Arleo
      iqw = 1
      alphas = 1d0/3d0
      scor = 1
      ncor = 0
      sfthrd = 1
c     real qhat !Transport coefficient (GeV^2.fm^-1)
      qhat = 0.40
c     drag coefficient
      ehat = 0.0
c     integer iDens !0= hard sphere, 1= Wood Saxon param
      iDens = 1
c     iqg = 1 -> quark and gluons are quenched other -> only q
      iqg = 1
c     iEg = 1 -> a gluon is added to satisfy energy conservation
      iEg = 0
c     iPtF = 0 -> no Pt; 1 -> from qhat; 2 -> BDMPS; 3 -> gluon angle from SW
      iPtF = 3
c     Suppretion factor
      SupFac= qhat /(qhat+ehat)
c      SupFac= 1.

c     integer iSim ! 0 -> Turn off Pythia for tests
      iSim = 1

ccc To save time with useless Pythia initialization
      i = 0
      if (iTg .eq. 0 .or. iFM.eq.0) then
        nkin = nevent +1
      endif
ccc Initialize
      ievent = 0
      call InitNucl
      if(rFM.ne.0) then
        call InitFM
        call GenNucDens
      endif
      call InitRandom
c      call CLASBOSINIT('MCEVENT')

ccc Main Loop
      do while (ievent.lt.nevent)

ccc Initialize the simulation
        if(iSim.ne.0.and.
     &   (ievent.eq.0.or.i.eq.nkin.or.ievent.eq.int(nevent*iZ/iA)
     &    .or. mod(ievent,500000).eq.0)) then
          i = 0
 100      continue
ccc Determine target: 1 for proton, 0 neutron
          if(ievent.lt.(nevent*iZ/iA)) then
            nucleon = 1
          else
            nucleon = 0
          endif

ccc Randomize Theta Phi and Kf
          if(rFM.ne.0.and.iFM.ne.0) call FMParam

ccc Initialize the kinematics
          call InitKin

ccc Going in the nucleon rest frame
          if(iColl.ne.0) call LorentzFM(2)
          if(iFM.ne.0) call LorentzFM(1)

ccc Security for low energies
          if (PPe .lt. 4) goto 100
 200      continue
ccc Parameters for Pythia
          call PythiaConfigDIS

ccc Stop fragmentation for QW to be applied
          if (iQuenching.ne.0) MSTJ(1) =0

ccc PYTHIA init
          BeamE = PPe
          if(ievent.lt.(nevent*iZ/iA)) then
            call pyinit('FIXT','gamma/e-','p+',BeamE)
          else
            write(*,*) 'X sec 99 = ', XSEC(99,1)
            call pyinit('FIXT','gamma/e-','n0',BeamE)
          endif
          if (XSEC(99,1).eq.0) then
            call pyrest
            goto 200
          endif
        endif

ccc Counter
        if (MOD(ievent,10000).eq.0) write(*,*) ievent,'events proceded'

ccc Initialization of the kinematic variables
        call InitKin2Book

ccc Event generation
        if(iSim.ne.0) CALL pyevnt
c        if(iSim.ne.0) CALL pylist(1)

ccc Come back in target frame
        if(iFM.ne.0) call LorentzFMBack(1)

ccc Energy loss of the partons
        if (iQuenching.ne.0.and.iTg.gt.1) then
          call InterPos ! Pick the position of the interaction in the nuclei
          call ApplyQW  ! Compute QW
        endif

ccc Fragmentation (if needed)
        if (iQuenching.ne.0) then
          MSTJ(1) =1
          if(iSim.ne.0) call PYEXEC
c          if(iSim.ne.0) CALL pylist(1)
          MSTJ(1) =0
        endif

ccc Go back in lab frame (for collider mode)
        if(iColl.ne.0) call LorentzFMBack(2)

ccc Compute of physical values for output
        call ComputV

ccc Book the ntuple
        call fillroot()
        ievent = ievent+ 1
        i = i+ 1
c          call CLASBOSFILL(iTg)
      enddo

ccc Close the file
c      call CLASBOSEND('MCEVENT')

      write(*,*) 'X sec 99 = ', XSEC(99,1)
      write(*,*) 'q hat = ', QW_qhat/QW_nb

      call TIMEX(T2)

      write(*,*) nevent,'events in ',T2-T1,'s'
   
      end

c------------------------------------------------------------------------------
c Initialize random number
c------------------------------------------------------------------------------
      subroutine InitRandom
      implicit none

      integer i ! for do loop
      real test ! dummy
      real ranf ! random number generator from CERNLIB
      integer initrm1,initrm2 ! to initialize ranf
      integer MRPY(6)
      double precision RRPY(100),PYR
      COMMON/PYDATR/MRPY,RRPY
      SAVE /PYDATR/
      integer*4 today(3), now(3)

      call idate(today)   ! today(1)=day, (2)=month, (3)=year
      call itime(now)     ! now(1)=hour, (2)=minute, (3)=second

      test = now(1)*now(2)/now(3)

      call datime(initrm1,initrm2)
      call ranset(initrm1*initrm2*int(test))

      MRPY(2) = 0
      MRPY(3) = mod(initrm1,85635)
      MRPY(4) = mod(initrm1,67)
      MRPY(5) = mod(initrm2,56)

      do i=1,100
       RRPY(i) = ranf(0)
      enddo

      do i=1,initrm2*int(test)
        test = ranf(0)
        test = PYR(0)
      enddo

      end

c------------------------------------------------------------------------------
c Initialize kinematics values
c------------------------------------------------------------------------------
      subroutine InitKin
      implicit none

      real me
      data me/0.000511/
      real mn
      data mn/0.938/
ccccc Include all the common blocks
      include 'common.f'

      PPe = E0
      EEe = sqrt(PPe**2 + me**2)
      Pex = 0
      Pey = 0
      Pez = E0
      PPn = Kf
      EEn = sqrt(PPn**2+mn**2)
      Pnx = sin(ThFM)*cos(PhiFM)*PPn
      Pny = sin(ThFM)*sin(PhiFM)*PPn
      Pnz = cos(ThFM)*PPn

c     if (iColl.eq.1) then
c       Pnz = Pnz - EColl
c       PPn = sqrt(Pnx**2 + Pny**2 + Pnz**2)
c       EEn = sqrt(PPn**2+mn**2)
c     endif

      ele_ene = EEe
      ele_the = 0
      ele_phi = 0
      nuc_mom = Kf
      nuc_the = ThFM
      nuc_phi = PhiFM

      end

c------------------------------------------------------------------------------
c Create spectators
c------------------------------------------------------------------------------
      subroutine CreateSpec()
      implicit none

      include 'common.f'

      N = N+1
      if (iTg .eq. 1 .and. nucleon .eq. 0) then
        specId = 2212
      else if (iTg .eq. 1 .and. nucleon .eq. 1) then
        specId = 2112
      else if (iTg .eq. 2 .and. nucleon .eq. 0) then
        specId = 10102
      else if (iTg .eq. 2 .and. nucleon .eq. 1) then
        specId = 10002
      else if (iTg .eq. 3 .and. nucleon .eq. 0) then
        specId = 10202
      else if (iTg .eq. 3 .and. nucleon .eq. 1) then
        specId = 10102
      else if (iTg .eq. 4 .and. nucleon .eq. 0) then
        specId = 10203
      else if (iTg .eq. 4 .and. nucleon .eq. 1) then
        specId = 10103
      endif
 
      k(N,1) = 1
      k(N,2) = specId
      k(N,3) = 2

      p(N,1) = -sin(nuc_the)*cos(nuc_phi)*nuc_mom
      p(N,2) = -sin(nuc_the)*sin(nuc_phi)*nuc_mom
      p(N,3) = -cos(nuc_the)*nuc_mom
      if (k(N,2) .eq. 2212) then
        p(N,5) = .938272
      else if (k(N,2) .eq. 2112) then
        p(N,5) = .939566
      else if (k(N,2) .eq. 10002) then
        p(N,5) = 1.87913
      else if (k(N,2) .eq. 10102) then
        p(N,5) = 1.876124
      else if (k(N,2) .eq. 10202) then
        p(N,5) = 1.87654
      else if (k(N,2) .eq. 10203) then
        p(N,5) = 2.809356
      else if (k(N,2) .eq. 10103) then
        p(N,5) = 2.809356
      endif
      P(N,4) = sqrt(P(N,1)**2+P(N,2)**2+P(N,3)**2+P(N,5)**2)

c      write(*,*) K(N,1),K(N,2),K(N,3)
c      write(*,*) P(N,1),P(N,2),P(N,3),P(N,4),P(N,5)

      end

      subroutine LorentzFM(i)
      implicit none

      include 'common.f'
      integer i

      if(i.eq.1) then
        BB1 = PPn/EEn
        B1x = Pnx/EEn
        B1y = Pny/EEn
        B1z = Pnz/EEn

        call TL(EEe,PPe,Pex,Pey,Pez,BB1,B1x,B1y,B1z)
        call TL(EEn,PPn,Pnx,Pny,Pnz,BB1,B1x,B1y,B1z)
      
ccc Rotate around y
        Thi = -atan2(Pex,Pez)
        call InitRotY(Thi)
ccc Rotate around z
        Phi = -atan2(Pey,Pez)
        call InitRotZ(Phi)
      else if(i.eq.2) then
        BB2 = -sqrt(EColl**2-.939**2)/EColl
        B2x = 0
        B2y = 0
        B2z = BB2

        call TL(EEe,PPe,Pex,Pey,Pez,BB2,B2x,B2y,B2z)
      endif

      end

      subroutine LorentzFMBack(i)
      implicit none

      include 'common.f'
      real Mom1,Mom2,Mom3,Mom4 
      real PPP 
      integer ip,i

ccc Rotate around z
      if (i.eq.1) call FinalRotZ(-Phi)
ccc Rotate around y
      if (i.eq.1) call FinalRotY(-Thi)

ccc Lorentz boost of all the particles
      do ip=1,N 
        Mom1 = p(ip,1)
        Mom2 = p(ip,2)
        Mom3 = p(ip,3)
        Mom4 = p(ip,4)
        if (i.eq.1) call TL(Mom4,PPP,Mom1,Mom2,Mom3,
     &                           -BB1,-B1x,-B1y,-B1z)
        if (i.eq.2) call TL(Mom4,PPP,Mom1,Mom2,Mom3,
     &                           -BB2,-B2x,-B2y,-B2z)
        p(ip,1) = Mom1
        p(ip,2) = Mom2
        p(ip,3) = Mom3
        p(ip,4) = Mom4
      enddo

      end

