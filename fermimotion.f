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
ccc Produce the table for FM generation (RW)
      else if (iFM.eq.5) then
        call GenRWtable()
      endif

      end

c Read table from RW
      subroutine GenRWtable()
      implicit none

      include 'common.f'

      integer ERROR
      real a,b,c,d,e
      real sum_n,sum_p
      integer i,j
      character*100 str

      write(*,*) 'Enter GenRWtable'

      ERROR = 0
      select case(iA)
        case(2)
          if(iZ.eq.1) then
            open(8,file='datafiles/fmrw/h2.momentum',status='old')
          else 
            ERROR = 1
          endif
        case(3)
          if(iZ.eq.1) then
            open(8,file='datafiles/fmrw/h3.momentum',status='old')
          else if(iZ.eq.2) then
            open(8,file='datafiles/fmrw/he3.momentum',status='old')
          else 
            ERROR = 1
          endif
        case(4)
          if(iZ.eq.2) then
            open(8,file='datafiles/fmrw/he4.momentum',status='old')
          else 
            ERROR = 1
          endif
        case(6)
          if(iZ.eq.3) then
            open(8,file='datafiles/fmrw/lad.momentum',status='old')
          else
            ERROR = 1
          endif
        case(7)
          if(iZ.eq.3) then
            open(8,file='datafiles/fmrw/lat.momentum',status='old')
          else 
            ERROR = 1
          endif
        case default
          ERROR = 1
      end select

      if (ERROR.eq.1) then
        write(*,*) 'ERROR'
        close (8)
        stop
      endif


      
      do while (str(2:2).ne.'*'.or.str(3:3).ne.'*')
        read(8,*) str
        write(*,*) str
        write(*,*) str(2:2)
        write(*,*) str(3:3)
      enddo            

      i = 0
      sum_n = 0
      sum_p = 0
      if(iZ.eq.1 .and. iA.eq.2) then
        do while (a.lt.FMlimit/0.1973269602)
          i = i+1
          read(8,*) a,b,c,d
          write(*,*) a,b,c,d
          FM_n(i) = d
          FM_p(i) = d
          sum_n = sum_n + d
          sum_p = sum_p + d
        enddo
      else if((iZ.eq.1 .and. iA.eq.3) .or. (iZ.eq.2 .and. iA.eq.3)
     &            .or. (iZ.eq.3 .and. iA.eq.7)) then
        do while (a.lt.FMlimit/0.1973269602)
          i = i+1
          read(8,*) a,b,c,d,e
          write(*,*) a,b,c,d,e
          FM_n(i) = d
          FM_p(i) = b
          sum_n = sum_n + d
          sum_p = sum_p + b
        enddo
      else if((iZ.eq.2 .and. iA.eq.4) .or. (iZ.eq.3 .and. iA.eq.6)) then
        do while (a.lt.FMlimit/0.1973269602)
          i = i+1
          read(8,*) a,b,c
          write(*,*) a,b,c
          FM_n(i) = b
          FM_p(i) = b
          sum_n = sum_n + b
          sum_p = sum_p + b
        enddo
      else
        write(*,*) 'ERROR'
        close (8)
        stop
      endif
              
      step_size_FM = FMlimit / (i-1)
      write(*,*) 'step ',step_size_FM

      do j=1,i
        FM_n(j) = FM_n(j) / sum_n
        FM_p(j) = FM_p(j) / sum_p
      enddo
                       
      do j=2,i
        FM_n(j) = FM_n(j-1) + FM_n(j)
        FM_p(j) = FM_p(j-1) + FM_p(j)
        write(*,*) j,FM_n(j),FM_p(j)
      enddo

      write(*,*) 'Exit GenRWtable'
      close (8)
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

ccc Fermi Momentum from R. Wiringa
      else if (iFM.eq.5) then
        Rd = ranf(0)
        i = 1
        if(nucleon.eq.0) then
          do while (FM_n(i).lt.Rd)
            i= i + 1
          enddo
        else if(nucleon.eq.1) then
          do while (FM_p(i).lt.Rd)
            i= i + 1
          enddo
        endif

        Kf =  (i - rand(0)) * step_size_FM
        if (Kf.gt.FMlimit) then
          write(*,*) 'warning ',i,step_size_FM
          stop
        endif

      else
        write(*,*) 'this iFM is not implemented'
        stop
      endif

      end
