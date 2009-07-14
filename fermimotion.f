c------------------------------------------------------------------------------
c Initialize fermi momentum
c------------------------------------------------------------------------------
      subroutine  InitFM(iTg,rFM,iFM)
      implicit none

      integer iTg ! = 0 no FM, = 1 deut, = 2 C, = 3 Al, = 4 Fe, = 5 Sn, = 6 Pb
      integer iZ, iA ! target Z and A
      real rFM ! Fermi momentum  in the target (GeV)
      
        select case(iTg)
          case (0) 
            rFM = 0
            iZ = 0
            iA = 0
          case (1) 
            rFM = 0.07
            iZ = 1
            iA = 2
          case (2) 
            rFM = 0.221
            iZ = 6
            iA = 12
          case (3) 
            rFM = 0.235
            iZ = 13
            iA = 27
          case (4) 
            rFM = 0.260
            iZ = 26
            iA = 56
          case (5) 
            rFM = 0.260
            iZ = 50
            iA = 120
          case (6) 
            rFM = 0.265
            iZ = 82
            iA = 208
          case default 
            rFM = 0
            iZ = 0
            iA = 0
        end select
c      if (iFM.eq.2 .and. iTg.gt.1) then 
        
     
c      endif

      end


c------------------------------------------------------------------------------
c Initialize kinematics values with fermi motion
c------------------------------------------------------------------------------
      subroutine FMParam(rFM,iFM)
      implicit none

      real E0 ! beam energy (GeV)
      real rFM ! Fermi momentum  in the target (GeV)
      integer iFM ! Flag for the kind of FM
      real pi
      data pi/3.1415926535/
      real a
      data a/2./
      real Plim
      data Plim/1./
      real R
      real ranf ! random number generator from CERNLIB
ccccc Important variables for simulation
      real Ps,Thr,C !variables for Kf computation
ccccc Include all the common blocks
      include 'common.f'

ccc Generate Theta and Phi of the particle's fermi momentum
      ThFM = acos(2*ranf(0)-1)
      PhiFM = 2*pi*ranf(0)

ccc Selector for the FM
      if(iFM.eq.0) then
ccc Generation of FM in a Fermi sphere
        Kf = rFM*(ranf(0))**(1./3.)
      else if(iFM.eq.1) then
c        write(*,*) 'Begining'
c        write(*,*) 'Kf',Kf
ccc Thresholds
        Thr = 1 - 6*(rFM*a/pi)**2
c        write(*,*) 'Thr',Thr
        Ps = 5
ccc Constant calculation
        C = 4./3.*pi*rFM**3.
c        write(*,*) 'C',C

ccc Remove events with Pf > 4 GeV/c or negative values
        do while (Ps.gt.Plim .or. Ps.lt.0)
ccc Generation of random number
          Kf = ranf(0)
ccc Apply the threshold and produce the tail
          if (Kf .le. Thr) then
            Kf = (3.*C*Kf/4./pi/Thr)**(1./3.)
c        write(*,*) 'Kf',Kf
          else
            R = 1. / (1.-rFM/4.)
            Kf = - rFM / ( (Kf-Thr)*pi*C/8./rFM**5/a**2 -1. )
c        write(*,*) 'Kf',Kf
          endif
          Ps = Kf
        enddo
c        write(*,*) 'Kf',Kf

      else
        write(*,*) 'this iFM is not implemented'
      endif

      end
