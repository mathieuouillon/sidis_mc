c------------------------------------------------------------------------------
c Initialize fermi momentum
c------------------------------------------------------------------------------
      subroutine  InitFM(iTg,rFM,iFM)
      implicit none

ccccc Include all the common blocks
      include 'common.f'
      integer iTg ! = 0 no FM, = 1 deut, = 2 C, = 3 Al, = 4 Fe, = 5 Sn, = 6 Pb
      integer iFM ! flag for the kind of FM
      real rFM ! Fermi momentum  in the target (GeV)
      integer irho !dummy
      
ccc Fill rFM, iZ and iA in function of the Target
        select case(iTg)
          case (0) 
            rFM = 0
            iZ = 1
            iA = 2
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
          case (7) 
            rFM = 0.120 !no source for this, i extrapolate between 6Li and 2H
            iZ = 2
            iA = 4
          case default 
            rFM = 0
            iZ = 0
            iA = 0
        end select

ccc Produce the table for FM generation (CS)
      if (iFM.eq.2 .or. iFM.eq.3) then 
        irho = iFM - 1
        call GenFMtable(irho)  
      endif

      end


c------------------------------------------------------------------------------
c Initialize kinematics values with fermi motion
c------------------------------------------------------------------------------
      subroutine FMParam(rFM,iFM,iTg)
      implicit none

      real E0 ! beam energy (GeV)
      real rFM ! Fermi momentum  in the target (GeV)
      integer iFM,iTg ! Flag for the kind of FM
      real a
      data a/2./
      real Plim
      data Plim/1./
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
      if(iFM.eq.0) then
ccc Generation of FM in a Fermi sphere
        Kf = rFM*(ranf(0))**(1./3.)
      else if(iFM.eq.1) then
ccc Thresholds
        Thr = 1 - 6*(rFM*a/pi)**2
        Ps = 5
ccc Constant calculation
        C = 4./3.*pi*rFM**3.

ccc Remove events with Pf > 4 GeV/c or negative values
        do while (Ps.gt.Plim .or. Ps.lt.0)
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
