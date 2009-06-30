      program simulation
      implicit none
c------------------------------------------------------------------------------
c TO DO LIST:
c 
c
c------------------------------------------------------------------------------
c RD Mar 20 2009: Integration of Pythia 6.4
c------------------------------------------------------------------------------
c RD Feb 19 2009: Integration of LEPTO
c------------------------------------------------------------------------------
c RD Feb 17 2009: Add the tail of fermi motion distribution
c------------------------------------------------------------------------------
c RD Feb 16 2009: hbook booking
c------------------------------------------------------------------------------
c RD Feb 12 2009: fermi motion in nuclei interaction
c------------------------------------------------------------------------------

ccccc Values for the simulation
      integer iTg ! = 0 no FM, = 1 deut, = 2 C, = 3 Al, = 4 Fe, = 5 Sn, = 6 Pb
      real rFM ! Fermi momentum  in the target (GeV)
      real E0 ! beam energy (GeV)
      integer i,nevent ! number of events
      integer j,nkin ! number of kinematics

ccccc Include all the common blocks
      include 'common.f'
      include 'names.inc'

ccccc Miscellaneous
      real ranf ! random number generator from CERNLIB
      integer icycle ! for hbook
      real T1,T2 ! For time of computation
      real PPP ! Dummy value
      real Mom1,Mom2,Mom3,Mom4 ! Dummy value
      double precision BeamE !Input value for pythia
      integer ip ! For do

ccc Begining of the simulation
      call TIMEX(T1)
      E0 = 5.014
      iTg = 4
      nkin = 100
      nevent = 200
      ievent = 0
      bosout = 'test.A00'

      if (iTg .eq. 0) then
        nevent = nkin*nevent
        nkin = 1
      endif
ccc Initialize
      call InitFM(iTg,rFM)
      call InitRandom
      call InitHbook
c      call CLASBOSINIT('MCEVENT')

      do j=1,nkin

 100    continue
ccc Randomize Theta Phi and Kf
        call FMParam(rFM)

ccc Initialize the kinematics
        call InitKin(E0,rFM)

ccc Radiative effects in targets

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

ccc Simulation
c        call PythiaConfigBrahim
        call PythiaConfigOWN
c        call PythiaConfigHayk
        MSTJ(1) =0
c        MSTP(143) =1

        write(*,*) 'Momentum of the electron: ',PPe
        BeamE = PPe
        write(*,*) 'Momentum of the electron: ',BeamE
        call pyinit('FIXT','gamma/e-','p+',BeamE)
c        call pyinit('FIXT','gamma/e-','n0',BeamE)


        do i=1,nevent
          ievent = ievent + 1
          if (MOD(i+(nevent*(j-1)),nevent*nkin/20) .eq. 0) 
     &          write(*,*) i+(nevent*(j-1)),'events proceded'

          call InitKin2Book

c see MSTP and PARP 171 for variable beam energy
          CALL pyevnt
ccc New stuff, that can be try
c          CALL pyevnw

          CALL PYLIST(1)
c         do ip=1,N
c           if(K(ip,1).lt.10.and.K(ip,2).lt.4.and.K(ip,2).gt.-4) then
c             write (*,*) P(ip,1)
c             p(ip,1) = .5 * p(ip,1)
c             p(ip,4) =sqrt(p(ip,1)**2+p(ip,2)**2+p(ip,3)**2+p(ip,5)**2)
c           endif
c         enddo
          MSTJ(1) =1
c          call PYSHOW
          call PYEXEC
          MSTJ(1) =0
          CALL PYLIST(1)

ccc Radiative correction


          if(iTg .ne. 0) then
ccc Rotate around z
            call FinalRotZ(-Phi)
ccc Rotate around y
            call FinalRotY(-Thi)

ccc Lorentz boost of all the particles
            do ip=1,N 
c              write(*,*) 'part', ip
c              write(*,*) p(ip,1),p(ip,2),p(ip,3),p(ip,4)
              Mom1 = p(ip,1)
              Mom2 = p(ip,2)
              Mom3 = p(ip,3)
              Mom4 = p(ip,4)
            
              call TL(Mom4,PPP,Mom1,Mom2,Mom3,
     &                                 -BB1,-B1x,-B1y,-B1z)

c              call TL(p(ip,4),PPP,p(ip,1),p(ip,2),p(ip,3),
c     &                                 -BB1,-B1x,-B1y,-B1z)

              p(ip,1) = Mom1
              p(ip,2) = Mom2
              p(ip,3) = Mom3
              p(ip,4) = Mom4
c              write(*,*) p(ip,1),p(ip,2),p(ip,3),p(ip,4)
            enddo
          endif

c          write(*,*) 'Momentum of the electron: ',BeamE
c          write(*,*) 'Momentum of the electron: ',p(1,4)
c          write(*,*) -BB1,-B1x,-B1y,-B1z
          call ComputV

ccc Radiative effects of the scattered particles

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
