      program simulation
      implicit none
c------------------------------------------------------------------------------
c TO DO LIST:
c   Implement position in the nuclei of the interaction (book it in the ntuple)
c   Implement Quenching 
c   Implement some radiative effect
c   Separate FM from nuclei calculcation
c   Check if there is flag for every options
c   Gluon or quark in QW to add
c   Smearing of the form Delta E = 3/8 alphas Delta Pts L 
c                              L = R/omega
c   Check if Pt of quarks are coherant
c
c------------------------------------------------------------------------------

ccccc Values for the simulation
      integer iTg ! = 0 no FM, = 1 deut, = 2 C, = 3 Al, = 4 Fe, = 5 Sn, = 6 Pb
      real rFM ! Fermi momentum  in the target (GeV)
      integer iFM ! 
c 0 = hard sphere with values from [1], 
c 1 = like 0 plus a tail from [2],
c 2 = Accardi SVG or Deuterium from Taya
c 3 = Accardi CS or Deuterium from Taya, 
c All FM distributions are limited to 1 GeV nucleons
c [1] E. J. Moniz et al. PRL 26, 445 (1971)
c [2] A. Bodek and J. L. Ritchie PRD 23, 1070 (1981)
      integer iDens !0= hard sphere, 1= Wood Saxon param
      integer iQuenching ! 0 desactivate Quenching
c     integer iqw 1 SW, 2 Arleo
      integer iSim ! 0 = Turn off Pythia
      real E0 ! beam energy (GeV)
      integer i,nevent ! number of events
      integer j,nkin ! number of kinematics
      real qhat !Transport coefficient (GeV^2.fm^-1)

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
 

ccc Begining of the simulation
      call TIMEX(T1)
      E0 = 5.014
      iTg = 1
      iFM = 3
      iDens = 1
      iSim = 1
      nkin = 1000
      nevent = 200
      ievent = 0
      bosout = 'test.A00'

ccc Init for the quenching weights
      iQuenching = 1
      iqw = 1
      alphas = 1d0/3d0
      scor = 1
      ncor = 0
      sfthrd = 1
      qhat =0.6

      if (iTg .eq. 0) then
        nevent = nkin*nevent
        nkin = 1
      endif
ccc Initialize
      call InitFM(iTg,rFM,iFM)
      call GenNucDens(iDens)
      call InitRandom
      call InitHbook
c      call CLASBOSINIT('MCEVENT')

      do j=1,nkin

 100    continue
ccc Randomize Theta Phi and Kf
        call FMParam(rFM,iFM,iTg)

ccc Initialize the kinematics
        call InitKin(E0,rFM)

ccc Going in the nucleon rest frame
        if(iTg .ne. 0) then
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
        call PythiaConfigHayk
        MSTJ(1) =0

ccc Initialize the simulation
        if(iSim.ne.0) write(*,*) 'Momentum of the electron: ',PPe
          BeamE = PPe
        if(iSim.ne.0) then
          if(j.lt.(nkin*iZ/iA)) then
            call pyinit('FIXT','gamma/e-','p+',BeamE)
          else
            call pyinit('FIXT','gamma/e-','n0',BeamE)
          endif
        endif

ccc Loop over events of a given kinematic
        do i=1,nevent

ccc Counter
          ievent = ievent + 1
          if (MOD(i+(nevent*(j-1)),nevent*nkin/20) .eq. 0) 
     &          write(*,*) i+(nevent*(j-1)),'events proceded'

ccc Some initialization
          call InitKin2Book

ccc Event generation
          if(iSim.ne.0) CALL pyevnt
c          if(iSim.ne.0) CALL PYLIST(1)

ccc Come back in Lab frame
          if(iTg .ne. 0) then
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
              if(K(ip,1).eq.2) then
 101            continue
                call QWComput(qhat,
     &                       P(ip,1),P(ip,2),P(ip,3),P(ip,4),K(ip,2))
                if (QW_w .gt. 0.) then
                  if (QW_w.lt.sqrt(P(ip,4)**2 -P(ip,5)**2)) then
                    ipt = sqrt(8*QW_w/3/alphas/QW_L)
                    ipl = dsqrt(P(ip,1)**2+P(ip,2)**2+P(ip,3)**2)
                    iplx = P(ip,1)/ipl
                    iply = P(ip,2)/ipl
                    iplz = P(ip,3)/ipl
                    iptz = 0.
                    if (iply .ne. 0) then
                      ipty = sqrt(iplx**2/(iply**2*((iplx/iply)**2+1)))
                      iptx = sqrt(1-ipty**2)
                      if (iplx*iply.gt.0) ipty = -ipty
                    else
                      ipty = 1.
                      iptx = 0.
                    endif
c                   write(*,*) 'Pl:',ipl,iplx,iply,iplz
c                   write(*,*) 'Pt:',ipt,iptx,ipty,iptz
c                   write(*,*) 'QW:',QW_w

                    ipl = ipl - QW_w
                    if (ipt.ge.ipl) then
                      ipx = ipl*iptx
                      ipy = ipl*ipty
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
                    goto 101
                  endif
                endif
              endif
            enddo
          endif

ccc Fragmentation
          MSTJ(1) =1
          if(iSim.ne.0) call PYEXEC
          MSTJ(1) =0
c          CALL PYLIST(1)

ccc Output to check Lorentz transforamtions
c          write(*,*) 'Momentum of the electron: ',BeamE
c          write(*,*) 'Momentum of the electron: ',p(1,4)

ccc Compute of physical values for the hbook
          call ComputV

ccc Book the ntuple
          call hfnt(33)
c          call CLASBOSFILL()
        enddo
      enddo

ccc Close the hbook file
      call hrout(33,icycle,' ')
      call hrend('out')
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
      subroutine InitKin(E0,rFM)
      implicit none

      real E0 ! beam energy (GeV)
      real rFM ! Fermi momentum  in the target (GeV)
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

      ele_ene = EEe
      ele_the = 0
      ele_phi = 0
      nuc_mom = Kf
      nuc_the = ThFM
      nuc_phi = PhiFM

      end
