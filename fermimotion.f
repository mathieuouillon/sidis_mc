c------------------------------------------------------------------------------
c Initialize fermi momentum
c------------------------------------------------------------------------------
      subroutine  InitFM(iTg,rFM)
      implicit none

      integer iTg ! = 0 no FM, = 1 deut, = 2 C, = 3 Al, = 4 Fe, = 5 Sn, = 6 Pb
      real rFM ! Fermi momentum  in the target (GeV)

      select case(iTg)
        case (0) 
          rFM = 0
        case (1) 
          rFM = 0.07
        case (2) 
          rFM = 0.221
        case (3) 
          rFM = 0.235
        case (4) 
          rFM = 0.260
        case (5) 
          rFM = 0.260
        case (6) 
          rFM = 0.265
        case default 
          rFM = 0
      end select

      end


c------------------------------------------------------------------------------
c Initialize kinematics values with fermi motion
c------------------------------------------------------------------------------
      subroutine FMParam(rFM)
      implicit none

      real E0 ! beam energy (GeV)
      real rFM ! Fermi momentum  in the target (GeV)
      real pi
      data pi/3.1415926535/
      real a
      data a/2./
      real Plim
      data Plim/4./
      real R
      real ranf ! random number generator from CERNLIB
ccccc Important variables for simulation
      real Ps,Thr,C !variables for Kf computation
ccccc Include all the common blocks
      include 'common.f'

ccc Generate Theta and Phi of the particle's fermi momentum
      ThFM = acos(2*ranf(0)-1)
      PhiFM = 2*pi*ranf(0)

ccc Thresholds
      Thr = 1 - 6*(rFM*a/pi)**2
      Ps = 5

ccc Constant calculation
      C = 4./3.*pi*rFM**3.

ccc Remove events with Pf > 4 GeV/c
c      do while (Ps.gt.Plim)

ccc Generate the random number
c        Kf = ranf(0)
        Kf = rFM*(ranf(0))**(1./3.)

ccc Apply the threshold and compute the Kf value
c        if (Kf .le. Thr) then
c          Kf = (3.*C*Kf/4./pi/Thr)**(1./3.)
c        else
c          R = 1. / (1.-rFM/4.)
c          Kf = - rFM / ( (Kf-Thr)*pi*C/8./rFM**5/a**2 -1. )
c        endif
c        Ps = Kf
c      enddo

      if (Kf.gt.4) write(*,*) 'HEYYYYYYYYYYYYYYYY'

c      write(*,*) ThFM,PhiFM,Kf

      end


