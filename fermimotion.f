c------------------------------------------------------------------------------
c Initialize fermi momentum
c------------------------------------------------------------------------------
      subroutine  InitNucl
      implicit none

ccccc Include all the common blocks
      include 'common.f'
      
ccc Fill rFM, iZ and iA in function of the Target
        select case(iTg)
          case (0) 
            rFM = 0
            iZ = 1
            iA = 1
!I(RD) extrapolate rFM values between 6Li and 2H
          case (1) 
            rFM = 0.07
            iZ = 1
            iA = 2
          case (2) 
            rFM = 0.1
            iZ = 1
            iA = 3
          case (3) 
            rFM = 0.1
            iZ = 2
            iA = 3
          case (4) 
            rFM = 0.120 
            iZ = 2
            iA = 4
          case (5) 
            rFM = 0.17
            iZ = 3
            iA = 6
          case (6) 
            rFM = 0.17
            iZ = 3
            iA = 7
          case (7) 
            rFM = 0.221
            iZ = 6
            iA = 12
          case (8) 
            rFM = 0.235
            iZ = 13
            iA = 27
          case (9) 
            rFM = 0.260
            iZ = 26
            iA = 56
          case (10) 
            rFM = 0.260
            iZ = 50
            iA = 120
          case (11) 
            rFM = 0.265
            iZ = 82
            iA = 208
          case default 
            rFM = 0
            iZ = 1
            iA = 1
        end select

      end

      subroutine  InitFM
      implicit none

ccccc Include all the common blocks
      include 'common.f'
      integer irho !dummy
      

ccc Produce the table for FM generation (CS)
      if (iFM.eq.2 .or. iFM.eq.3) then 
        irho = iFM - 1
        call GenFMtable(irho)  
      endif

      end


c------------------------------------------------------------------------------
c Initialize kinematics values with fermi motion
c------------------------------------------------------------------------------
      subroutine FMParam
      implicit none

      real a
      data a/2./
      real R

ccccc Important variables for simulation
      real Ps,Thr,C,Rd !variables for Kf computation
      integer i
ccccc Include all the common blocks
      include 'common.f'

ccc Generate Theta and Phi of the particle's fermi momentum
      ThFM = acos(2*ranf(0)-1)
      PhiFM = 2*pi*ranf(0)

ccc Selector for the FM
      if(iFM.eq.4) then
ccc Generation of FM in a Fermi sphere
        Kf = rFM*(ranf(0))**(1./3.)
      else if(iFM.eq.1) then
ccc Thresholds
        Thr = 1 - 6*(rFM*a/pi)**2
        Ps = 5
ccc Constant calculation
        C = 4./3.*pi*rFM**3.

ccc Remove events with Pf > 4 GeV/c or negative values
        do while (Ps.gt.FMlimit .or. Ps.lt.0)
ccc Generation of random number
          Kf = ranf(0)
ccc Apply the threshold and produce the tail
          if (Kf .le. Thr) then
            Kf = (3.*C*Kf/4./pi/Thr)**(1./3.)
          else
            R = 1. / (1.-rFM/4.)
            Kf = - rFM / ( (Kf-Thr)*pi*C/8./rFM**5/a**2 -1. )
          endif
          Ps = Kf
        enddo

ccc Fermi Momentum from Accardi routines
      else if (iFM.eq.2 .or. iFM.eq.3) then
        Rd = ranf(0)
        i = 1
        do while (FM_table(i).lt.Rd)
          i= i + 1
        enddo

        Kf =  (i - rand(0)) * step_size_FM

      else
        write(*,*) 'this iFM is not implemented'
        stop
      endif

      end
