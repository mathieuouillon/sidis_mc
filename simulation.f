      integer function simulation()
      implicit none
c------------------------------------------------------------------------------
c TO DO LIST:
c   Implement some radiative effect
c   Check if there is flag for every options
c   CoM output energy broken
c------------------------------------------------------------------------------

ccccc Include all the common blocks
      include 'common.f'
      include 'names.inc'

ccccc Miscellaneous
      integer icycle ! for hbook
      real T1,T2 ! For time of computation
      real PPP ! Dummy value
      real Mom1,Mom2,Mom3,Mom4 ! Dummy value
      double precision BeamE !Input value for pythia
      integer ip ! For do
      real ipx,ipy,ipz,E !dummy only for test
      real ipl,ipt
      integer flag
      real iplx,iply,iplz
      real iptx,ipty,iptz
      real ipix,ipiy,ipiz
      real th,ph
      integer i,j
 

      call TIMEX(T1)
CCCCCC Begining of the simulation
ccc Integer j,nkin ! number of kinematics
      nkin = 2000
ccc Number of events per kinematics
      nevent = 500
ccc Electron energy (GeV)
      E0 = 5.014
ccc Target type ! 0 proton, 1 deut, 2 C, 3 Al, 4 Fe, 5 Sn, 6 Pb, 7 He4
      iTg = 4

ccc Collider options
c     integer iColl !1 = activate collider kinematic
      iColl = 0
c     real EColl ! energy of the nuclei (GeV/nucleon)
      EColl = 5.0

ccc Fermimotion flag 
c     0 = hard sphere with values from [1], 
c     1 = like 0 plus a tail from [2],
c     2 = Accardi SVG (put list of available nuclei)
c     3 = Accardi CS (put list of available nuclei)
c     All FM distributions are limited to 1 GeV nucleons
c     [1] E. J. Moniz et al. PRL 26, 445 (1971)
c     [2] A. Bodek and J. L. Ritchie PRD 23, 1070 (1981)
      iFM = 3

ccc Integer iNS ! 0 = no nuclear spectator, 1 = nuclear spectator
c                 ! this option is only for 2H and 4He targets
      iNS = 0

ccc CLAS12 Acceptance put 1 
      iAccept = 0

ccc dummy
      bosout = 'test.A00'
      hbookout = 'helium.hbook'

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
      qhat =0.6
c     integer iDens !0= hard sphere, 1= Wood Saxon param
      iDens = 1

c     integer iSim ! 0 = Turn off Pythia for tests
      iSim = 1

ccc To save time with useless Pythia initialization
      if (iTg .eq. 0) then
        nevent = nkin*nevent
        nkin = 1
      endif
ccc Initialize
      ievent = 0
      call InitNucl
      if(rFM.ne.0) then
        call InitFM
        call GenNucDens
      endif
      call InitRandom
c      call InitHbook
      if (iAccept.eq.1) call readtables
c      call CLASBOSINIT('MCEVENT')

      do j=1,nkin

 100    continue
ccc Randomize Theta Phi and Kf
      if(rFM.ne.0) call FMParam

ccc Initialize the kinematics
        call InitKin

ccc Going in the nucleon rest frame
        if(iTg .ne. 0 .or. iColl.ne.0) then
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

ccc Center of mass energy calculation
          Beta = (EEe - 0.939) / PPe
          ECoM = 2*EEe*(1-Beta)/sqrt(1-Beta**2)

        else
          BB1 = 0
          B1x = 0
          B1y = 0
          B1z = 0
          Thi = 0 
          Phi = 0 
          Beta = (EEe - 0.939) / PPe 
          ECoM = 2*EEe*(1-Beta)/sqrt(1-Beta**2)
        endif

        if (PPe .lt. 2.5) goto 100

ccc Parameters for Pythia
c        call PythiaConfigBrahim
c        call PythiaConfigOWN
c        call PythiaConfigHayk
        call PythiaConfigCLAS

ccc Block fragmentation if QW will be applied
        if (iQuenching.ne.0) MSTJ(1) =0


ccc Initialize the simulation
        if(iSim.ne.0) write(*,*) 'Momentum of the electron: ',PPe
          BeamE = PPe
        if(iSim.ne.0) then
          if(j.lt.(nkin*iZ/iA)) then
            call pyinit('FIXT','gamma/e-','p+',BeamE)
            nucleon = 1
          else
            call pyinit('FIXT','gamma/e-','n0',BeamE)
            nucleon = 0
          endif
        endif

ccc Loop over events of a given kinematic
        do i=1,nevent

ccc Counter
          ievent = ievent + 1
          if (MOD(i,10000) .eq. 0) write(*,*) ievent,'events proceded'

ccc Some initialization
          call InitKin2Book

ccc Event generation
          if(iSim.ne.0) CALL pyevnt
c          if(iSim.ne.0) CALL PYLIST(1)

ccc Come back in Lab frame
          if(iTg .ne. 0 .or. iColl.ne.0) then
ccc Rotate around z
            call FinalRotZ(-Phi)
ccc Rotate around y
            call FinalRotY(-Thi)

ccc Lorentz boost of all the particles
            do ip=1,N 
              Mom1 = p(ip,1)
              Mom2 = p(ip,2)
              Mom3 = p(ip,3)
              Mom4 = p(ip,4)
              call TL(Mom4,PPP,Mom1,Mom2,Mom3,
     &                                 -BB1,-B1x,-B1y,-B1z)
              p(ip,1) = Mom1
              p(ip,2) = Mom2
              p(ip,3) = Mom3
              p(ip,4) = Mom4
            enddo
          endif

ccc Energy loss of the partons
          if (iQuenching.ne.0.and.iTg.gt.1) then
            call InterPos
            flag = 0
            do ip =1,N
c              if((abs(K(ip,2)).lt.6.or.K(ip,2).eq.21).and.K(ip,1).lt.9) then
              if(abs(K(ip,2)).lt.6.and.K(ip,1).lt.9) then
 101            continue
                call QWComput(P(ip,1),P(ip,2),P(ip,3),P(ip,4),K(ip,2))
                if (QW_w .gt. 0.) then
                  if (QW_w.lt.sqrt(P(ip,4)**2 -P(ip,5)**2)-.25) then
                    ipt = sqrt(8*QW_w/3/alphas/QW_L)
                    ipl = dsqrt(P(ip,1)**2+P(ip,2)**2+P(ip,3)**2)
                    iplx = P(ip,1)/ipl
                    iply = P(ip,2)/ipl
                    iplz = P(ip,3)/ipl
                    iptz = 0.
                    if (iply .ne. 0) then
c                      ipty = sqrt(iplx**2/(iply**2*((iplx/iply)**2+1)))
c                      iptx = sqrt(1-ipty**2)
c                      if (iplx*iply.gt.0) ipty = -ipty
                       ipiz = ranf(0)*4*asin(1.)
                       ipix = cos(iptz)
                       ipiy = sin(iptz)
                       ipiz = 0
                       iptx = -ipix*iplx*iplz/sqrt(1-iplz**2)
     &                        -ipiy*iply/sqrt(1-iplz**2)
                       ipty = -ipix*iply*iplz/sqrt(1-iplz**2)
     &                        -ipiy*iplx/sqrt(1-iplz**2)
                       iptz = ipix*sqrt(1-iplz**2)
                    else
                      ipty = 1.
                      iptx = 0.
                    endif
c                   write(*,*) 'Pl:',ipl,iplx,iply,iplz
c                   write(*,*) 'Pt:',ipt,iptx,ipty,iptz
c                   write(*,*) 'QW:',QW_w
c                   write(*,*) 'test',iplx*iptx+iply*ipty+iplz*iptz

                    ipl = ipl - QW_w
                    if (ipt.ge.ipl) then
                      ipx = ipl*iptx
                      ipy = ipl*ipty
                      ipz = ipl*iptz
                    else
                      ipl = sqrt(ipl**2 - ipt**2)
                      ipx = ipt*iptx+ipl*iplx
                      ipy = ipt*ipty+ipl*iply
                      ipz = ipt*iptz+ipl*iplz
                    endif

                    P(ip,1) = ipx
                    P(ip,2) = ipy
                    P(ip,3) = ipz
                    P(ip,4) = sqrt(P(ip,5)**2+ipx**2+ipy**2+ipz**2)
                  else
                    th = acos(2*ranf(0)-1)
                    ph = 2*pi*ranf(0)
                    P(ip,1) = sin(th)*cos(ph)*.25
                    P(ip,2) = sin(th)*sin(ph)*.25
                    P(ip,3) = cos(th)*.25
                    P(ip,4) = sqrt(P(ip,5)**2+ipx**2+ipy**2+ipz**2)
                  endif
                endif
              endif
            enddo
          endif

ccc Fragmentation
          if (iQuenching.ne.0) then
            MSTJ(1) =1
            if(iSim.ne.0) call PYEXEC
            MSTJ(1) =0
c            CALL PYLIST(1)
          endif
ccc Output to check Lorentz transforamtions
c          write(*,*) 'Momentum of the electron: ',BeamE
c          write(*,*) 'Momentum of the electron: ',p(1,4)
c           write(*,*) PPn,Pnx,Pny,Pnz
c           write(*,*) nuc_mom,nuc_the,nuc_phi 

          if (iNS .eq. 1 .and. (iTg .eq. 1 .or. iTg .eq. 7)) then
            N = N+1
            if (iTg .eq. 1 .and. nucleon .eq. 0) then
              specId = 2212
            else if (iTg .eq. 1 .and. nucleon .eq. 1) then
              specId = 2112
            else if (iTg .eq. 7 .and. nucleon .eq. 0) then
              specId = 10203
            else if (iTg .eq. 7 .and. nucleon .eq. 1) then
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
            else if (k(N,2) .eq. 10203) then
              p(N,5) = 2.809356
            else if (k(N,2) .eq. 10103) then
              p(N,5) = 2.809356
            endif
            P(N,4) = sqrt(P(N,1)**2+P(N,2)**2+P(N,3)**2+P(N,5)**2)

c            write(*,*) K(N,1),K(N,2),K(N,3)
c            write(*,*) P(N,1),P(N,2),P(N,3),P(N,4),P(N,5)
          endif

ccc Compute of physical values for the hbook
          call ComputV

ccc Book the ntuple
          call fillroot()
c          call hfnt(33)
c          call CLASBOSFILL(iTg)
        enddo
      enddo

ccc Close the hbook file
c      call hrout(33,icycle,' ')
c      call hrend('out')
c      call CLASBOSEND('MCEVENT')

      call TIMEX(T2)

      write(*,*) nevent*nkin,'events in ',T2-T1,'s'
   
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

      call datime(initrm1,initrm2)
      call ranset(initrm1*initrm2)

      MRPY(2) = 0
      MRPY(3) = mod(initrm1,85635)
      MRPY(4) = mod(initrm1,67)
      MRPY(5) = mod(initrm2,56)

      do i=1,100
       RRPY(i) = ranf(0)
      enddo

      do i=1,initrm2
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

      if (iColl.eq.1) then
        Pnz = Pnz - EColl
        PPn = sqrt(Pnx**2 + Pny**2 + Pnz**2)
        EEn = sqrt(PPn**2+mn**2)
      endif

      ele_ene = EEe
      ele_the = 0
      ele_phi = 0
      nuc_mom = Kf
      nuc_the = ThFM
      nuc_phi = PhiFM

      end
