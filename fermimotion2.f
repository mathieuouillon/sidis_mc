      subroutine GenFMtable(iZ,iA,irho)
      implicit none

      integer i,iZ,iA,irho
      double precision rhofermi,mom,proba,ptot,step
      include 'common.f'

      step_size = 0.001
      FMnb = 1000
      mom = 0.
      proba = 0.
      ptot = 0.

      do i=1,FMnb
       proba = 3.14159265*mom*mom*rhofermi(iZ,iA,mom,irho)*step_size
       ptot = ptot + proba
       table(i) = ptot
       mom = mom + step_size
      enddo

      do i=1,FMnb
        table(i) = table(i)/ptot
      enddo

      end
************************************************************************
*  Subroutines for nucleon Fermi motion                                *
*  programmer: Alberto Accardi                                         *
*  created: 12 Jun 2005                                                *
*                                                                      *
*  SUMMARY:                                                            *
*                                                                      *
*   I. Single nucleon Fermi distributions - main front-end             *
*  II. SVG distributions (Sandel-Vary-Garpman)                         *
* III. CS distributions (Ciofi-Simula)                                 *
*                                                                      *
*  HISTORY:                                                            *
*                                                                      *
*  * fermimotion   : the _first_ one                                   * 
*    (12 Jun 05)    --- WORKING VERSION ---                            *
*  * fermimotion2  : Extracts soft CS (normalized and unnormalized     * 
*    (31 Aug 06)     and the unnormalized soft SVG                     *
*                                                                      *
*  NEEDS:                                                              *
*   - fermimotion2.SVG.tbl                                             *
*   - fermimotion2.CS.tbl                                              *
*                                                                      *
*  TO DO LIST                                                          *
*  ----------                                                          *
*                                                                      *
************************************************************************


************************************************************************
*                                                                      *
*   I. Single nucleon Fermi momentum distributions                     * 
*                                                                      *
************************************************************************

***********************************************************************
*     Nucleon Fermi momentum distribution [GeV^-3]
      function rhofermi(iZ,iA,k,irho)
*     programmer: Alberto Accardi
*     date: 25 May 2005
*
*  A. COMMENTARY
*
*     Returns the nucleon Fermi momentum distribution in a 
*     nucleus (iZ,iA) normalized to 1.
*
*       rhofermi [GeV^-3] = (dp) Fermi momentum distribution 
*                                normalized to 1
*       iZ,iA             = (i)  Atomic and mass numbers
*       k [GeV]           = (dp) nucleon momentum
*       irho              = (i)  distribution model
*                                1=SVG [1]  2=CS [3]
*
*     SVG model
*     ---------
*
*     Short range correations (hard part) fitted to computations in [2].
*     Long range correlations (soft part) taken from [1], and 
*     normalized so that Int(rho)=1. 
*
*       rho = N*rho_soft + rho_hard
*
*     Where rho_soft = rhosoftSVG(p) and rho_hard = rhohardSVG(p) 
*     have [fm^3] dimension and argument p with [fm^-1] dimension        
*
*     It allows any A>=7. However: 
*     1) The soft part is parametrized only for      
*        a subset of nuclei (7Li, 12C, 160, 20Ne, 27Al, 40Ar, 40Ca, 58Ni, 
*        63Cu, 90Zr, 120Sn, 208Pb). For other nuclei the one in the above 
*        list which is closest in A is used (a warning is printed on screen).
*     2) The hard tail (short range correlations) is taken from the 
*        parametrization in Ref. [2]. The slope is fixed to 0.22 fm^2. 
*        The height presented in [2] for a few nuclei (2H 3He 12C 16O 
*        40Ca 56Fe an 208Pb - see coments below for 8Be) has been fitted 
*        to a functional dependence on A. 
*        The fit is accurate to the +-3% level, except for 4He, which 
*        is underestimated by ~30%. In this case, the original 
*        value from [1] is used. Caution should be used when considering 
*        light nuclei, in general.
*
*     CS model
*     --------
*
*     It is a parametrization taken from [2] of theoretical computation 
*     for a few nuclei (12C 16O 40Ca 56Fe and 208Pb) made by many authors. 
*     Parameters for 8Be fitted to results in [3] by A.A. 
*     It allows use of the above listed nuclei only.
*
*     References: 
*     [1] M.Sandel, J.P.Vary and S.I.A.Garpman, Phys.Rev.C20(1979)744
*     [2] C.Ciofi degli Atti and S.Simula, Pys.Rev.C53(1996)1689
*     [3] R.B.Wiringa et al., Phys.Rev.C62(2001)014001
*
*  B. DECLARATIONS
*
      implicit none

      double precision rhofermi,k
      integer iZ,iA,irho

*    *** variables 

      double precision kk

*    *** functions

      double precision rhoSVGsh,rhoCSsh

*    *** constants & initial values

*     Conversion factors
      double precision FmGeV,Fm3GeV3,pi,eightpi
      parameter(FmGeV=0.1973269602,Fm3GeV3=FmGeV*FmGeV*FmGeV
     &     ,pi=3.1415926535898d0,eightpi=8d0*pi)

*
*  C. ACTION   
*
  
      kk = k/FmGeV

      if (irho.eq.1) then
        rhofermi = rhoSVGsh(iZ,iA,kk,1) / Fm3GeV3
      else if (irho.eq.2) then
         rhofermi = rhoCSsh(iZ,iA,kk,1) / Fm3GeV3
      end if

      return
      end


************************************************************************
*                                                                      *
*  II. SVG distributions                                               * 
*                                                                      *
************************************************************************

************************************************************************
*     SVG parametrization of Fermi momentum distribution [fm^3]
      function rhoSVGsh(iZ,iA,k,idist)
*     programmer: Alberto Accardi
*     date: 31 Aug 2006
*
*  A. COMMENTARY
*
*     Front-end which retutrns "soft" and "hard" SVG distributions.
*
*     Short range correations (hard part) fitted to computations in [2].
*     Long range correlations (soft part) taken from [1], and 
*     normalized so that Int(rho)=1. 
*
*       rho = N*rho_soft + rho_hard
*
*     Where rho_soft = rhosoftSVG(p) and rho_hard = rhohardSVG(p) 
*     have [fm^3] dimension and argument p with [fm^-1] dimension        
*
*     It allows any A>=7. However: 
*     1) The soft part is parametrized only for      
*        a subset of nuclei (7Li, 12C, 160, 20Ne, 27Al, 40Ar, 40Ca, 58Ni, 
*        63Cu, 90Zr, 120Sn, 208Pb). For other nuclei the one in the above 
*        list which is closest in A is used (a warning is printed on screen).
*     2) The hard tail (short range correlations) is taken from the 
*        parametrization in Ref. [2]. The slope is fixed to 0.22 fm^2. 
*        The height presented in [2] for a few nuclei (2H 3He 12C 16O 
*        40Ca 56Fe an 208Pb - see coments below for 8Be) has been fitted 
*        to a functional dependence on A. 
*        The fit is accurate to the +-3% level, except for 4He, which 
*        is underestimated by ~30%. In this case, the original 
*        value from [1] is used. Caution should be used when considering 
*        light nuclei, in general.
*
*     Normalization: \int(d3k rhoCSsh) = 1 (except for idist=2,4)
*
*       rhoSVGsh [fm3] = (dp) Fermi momentum distribution 
*                          normalized to 4*pi
*       iZ,iA          = (i)  Atomic and mass numbers
*       k [fm^-1]      = (dp) nucleon momentum
*       idist          = (i)  1 = soft+hard
*                             2 = soft (unnormalized)
*                             3 = soft (normalized to 1)
*                             4 = hard (unnormalized)
*                            (note: 1=2+4)
*
*     References: 
*     [1] M.Sandel, J.P.Vary and S.I.A.Garpman, Phys.Rev.C20(1979)744
*
*  B. DECLARATIONS
*
      implicit none

      double precision rhoSVGsh,k
      integer iZ,iA,idist

*    *** variables 

      double precision norm,rs,rh,height,slope

*    *** functions

      double precision rhosoftSVG,rhohardSVG

*    *** constants & initial values

*     Conversion factors
      double precision FmGeV,Fm3GeV3,pi,eightpi
      parameter(FmGeV=0.1973269602,Fm3GeV3=FmGeV*FmGeV*FmGeV
     &     ,pi=3.1415926535898d0,eightpi=8d0*pi)


*
*  C. ACTION   
*

      if (idist.eq.1) then
*       ... soft+hard
         rs = rhosoftSVG(iZ,iA,k)
         rh = rhohardSVG(iZ,iA,k,height,slope,norm)
         rhoSVGsh = (norm*rs + rh) 
      else if (idist.eq.2) then
*       ... soft (unnormalized)
         rs = rhosoftSVG(iZ,iA,k)
         rh = rhohardSVG(iZ,iA,k,height,slope,norm)
         rhoSVGsh = norm*rs 
      else if (idist.eq.3) then
*       ... soft (normalized)
         rhoSVGsh = rhosoftSVG(iZ,iA,k)
      else if (idist.eq.4) then
*       ... soft (unnormalized)
         rhoSVGsh = rhohardSVG(iZ,iA,k,height,slope,norm)
      end if
      return
      end

************************************************************************
*     Soft Nucleon Fermi momentum distribution [fm^3]
      function rhosoftSVG(iZ,iA,k)
*     programmer: Alberto Accardi
*     date: 25 May 2005
*
*  A. COMMENTARY
*
*     Returns the (absolute value) of the soft part of nucleon 
*     Fermi momentum distribution in a nucleus (iZ,iA), according 
*     to Ref.[1] 
*
*       rhosoftSVG [fm3] = (dp) Fermi momentum distribution 
*                               normalized to 1
*       iZ,iA            = (i)  Atomic and mass numbers
*       k [fm^-1]        = (dp) nucleon momentum
*
*     NOTE: absolute value of the soft Fermi distribution is taken 
*     because the parametrization chosen in [1] leads to negative 
*     values at very large momenta, well outside the region of 
*     validity of the parametrization itself. At these values of the 
*     nucleon momentum the distribution is nonetheless very small 
*     and completely subdominant compared to the hard tail.
*
*     Reference: 
*     [1] M.Sandel, J.P.Vary and S.I.A.Garpman, Phys.Rev.C20(1979)744
*
*  B. DECLARATIONS
*
      implicit none

      double precision rhosoftSVG,k
      integer iZ,iA

*    *** variables 

      integer nmax,nin,i,j,n
      save n

      parameter (nmax=20)
      double precision a(6,nmax),alpha(nmax),beta(nmax),nj
      integer ZZ(nmax),AA(nmax)
      save a,alpha,beta,ZZ,AA

      character string*100

*    *** functions

      integer nextunit
      double precision gaussnd

*    *** constants & initial values

*     Constants
      double precision pi
      parameter(pi=3.1415926535898d0)

*     initialization variables
      logical firsttime
      integer Zold,Aold,nnuke
      save Zold,Aold,firsttime,nnuke
 
      data Zold/0/ Aold/0/ firsttime/.true./
*
*  C. ACTION   
*

      if (iA.lt.7) then
         print*
         print*, 'ERROR (rhosoftSVG): called with A<7: ', iA
         print*, ' -- Try using Ciofi-Simula parametrization instead'
         print*
         stop
      end if

*    *** INITIALIZATION 
      if (firsttime) then
         firsttime = .false.
         nin = nextunit()
         open(nin,file='fermimotion2.SVG.tbl',status='old')
*       ... skips headers
         do i = 1, 7 
            read(nin,*) string
         end do
*       ... reads values
         j=0
 10      j=j+1
         read(nin,*) string,ZZ(j),AA(j),a(1,j),a(2,j),a(3,j),a(4,j)
     &        ,a(5,j),a(6,j),alpha(j),beta(j)
         if (ZZ(j).ne.0) goto 10 
         nnuke = j-1
         close(nin)
      end if

*       ... if new nucleus, locates nearest nucleus
      if ((iZ.ne.Zold).or.(iA.ne.Aold)) then
         Zold = iZ
         Aold = iA
         call ilocatetab(ZZ,nnuke,iZ,j)
         if (j.eq.0) then 
            n = 1
         else if ((ZZ(j).eq.iZ).and.(AA(j).eq.iA)) then
*          ... if iZ and iA are found in the table OK 
            n = j
         else 
*          ... otherwise chooses the nucleus with closest A
            if ((iA-AA(j)).le.(AA(j+1)-iA)) then
               n = j
            else
               n = j+1
            end if
         end if
      end if

*    *** COMPUTES the soft Fermi momentum distribution
      rhosoftSVG = 0d0
      do j = 1, 6
         nj = alpha(n)*beta(n)**j
         rhosoftSVG = rhosoftSVG + a(j,n)*(nj/pi)**1.5
     &        * dexp(-(nj*k*k))
      end do

      rhosoftSVG = dabs(rhosoftSVG)

*      rhosoftSVG = gaussnd(k,0d0,0.0351019399/3.8937966d-2/3d0,3)

      return
      end


************************************************************************
*     Hard nucleon Fermi momentum distribution [fm^3]
      function rhohardSVG(iZ,iA,k,height,slope,norm)
*     programmer: Alberto Accardi
*     date: 25 May 2005
*
*  A. COMMENTARY
*
*     Returns the "hard tai" of nucleon Fermi momentum distribution in a 
*     nucleus (iZ,iA), and its height, slope and normalization. 
*
*       rhofermi [fm^3]   = (dp) Fermi momentum distribution 
*       iZ,iA             = (i)  Atomic and mass numbers
*       k [fm^-1]         = (dp) nucleon momentum
*       height            = (dp) [OUT] height of exponential (hard) tail
*       slope             = (dp) [OUT] slope of exponential (hard) tail
*       norm              = (dp) [OUT] normalization = int(d3k*rhohard) 
*
*     Short range correlations (hard tail) are taken from the 
*     parametrization in Ref. [1]. The slope is fixed to 0.22 fm^2. 
*     The height presented in [1] has been fitted to a functional 
*     dependence on A. 
*     The fit included the following nuclei: 2H 3He 12C 16O 
*     40Ca 56Fe an 208Pb, taken from [1] and 8Be fitted by A.Accardi 
*     to the computations of [2]. The fit is accurate to the +-3% level, 
*     except for 4He, which is underestimated by ~30%. In this case, 
*     the original value from [1] is used. Caution should be used when 
*     considering light nuclei, in general.
*
*     Reference: 
*     [1] C.Ciofi degli Atti and S.Simula, Pys.Rev.C53(1996)1689
*     [2] R.B.Wiringa et al., Phys.Rev.C62(2001)014001
*
*     TO DO LIST
*     ----------
*
*  B. DECLARATIONS
*
      implicit none

      double precision rhohardSVG,k,height,slope,norm
      integer iZ,iA

*    *** variables 

      double precision AA

*    *** constants & initial values

      double precision pi,fourpi
      parameter(pi=3.1415926535898d0,fourpi=4d0*pi)

      logical firsttime
      save firsttime
      data firsttime/.true./

*
*  C. ACTION   
*

      if (firsttime) then 
         firsttime = .false.
         if (iA.lt.8) then
            print*
            print*, '**************************************************'
            print*, 'WARNING (rhohardSVG): called with 2<A<8: A=', iA
            print*, ' -- parametrization of hard tail''s height'
            print*, '    should be used with caution'
            print*, '**************************************************'
            print*
         end if
      end if

      if ((iZ.eq.2).and.(iA.eq.3)) then
         slope = 0.234d0
      else
         slope = 0.220d0
      end if

      if ((iZ.eq.2).and.(iA.eq.4)) then
         height = 0.0244/fourpi
      else
         AA = dble(iA)
         height = 0.00623/fourpi
     &     * (1d0 + 2.69d0*dlog(AA/2d0)**0.695d0/AA**0.133d0)
      end if

      norm = 1d0 - height * (pi/slope)**1.5 

      rhohardSVG = height * dexp(-slope*k*k)
      
      return
      end


************************************************************************
*                                                                      *
*  III. SVG distributions                                              * 
*                                                                      *
************************************************************************


************************************************************************
*     Ciofi-Simula parametrization of Fermi momentum distribution [fm^3]
      function rhoCS(iZ,iA,k)
*     programmer: Alberto Accardi
*     date: 31 Aug 2006
*
*  A. COMMENTARY
*
*     Front-end to "rhoCSsh" for compatibility with previous versions
*     and to simplify the syntax. It returns the full soft+hard CS 
*     parametrization taken from [1] of theoretical computation 
*     for a few nuclei (12C 16O 40Ca 56Fe and 208Pb) made by many authors. 
*     Parameters for 8Be fitted to results in [2] by A.A. 
*     It allows use of the above listed nuclei only.
*
*     Normalization: \int(d3k rhoCSsh) = 1
*
*       rhoCS [fm3] = (dp) Fermi momentum distribution 
*                          normalized to 4*pi
*       iZ,iA       = (i)  Atomic and mass numbers
*       k [fm^-1]   = (dp) nucleon momentum
*
*     References: 
*     [1] C.Ciofi degli Atti and S.Simula, Pys.Rev.C53(1996)1689
*     [2] R.B.Wiringa et al., Phys.Rev.C62(2001)014001
*
*
*  B. DECLARATIONS
*
      implicit none

      double precision rhoCS,rhoCSsh,k
      integer iZ,iA
*
*  C. ACTION   
*

      rhoCS = rhoCSsh(iZ,iA,k,1)

      return
      end


************************************************************************
*     Ciofi-Simula parametrization of Fermi momentum distribution [fm^3]
      function rhoCSsh(iZ,iA,k,idist)
*     programmer: Alberto Accardi
*     date: 25 May 2005
*           31 Aug 2006
*
*  A. COMMENTARY
*
*     It is a parametrization taken from [1] of theoretical computation 
*     for a few nuclei (12C 16O 40Ca 56Fe and 208Pb) made by many authors. 
*     Parameters for 8Be fitted to results in [2] by A.A. 
*     It allows use of the above listed nuclei only.
*
*     Normalization: \int(d3k rhoCSsh) = 1 (except for idist=2,4)
*
*       rhoCS [fm3] = (dp) Fermi momentum distribution 
*                          normalized to 4*pi
*       iZ,iA       = (i)  Atomic and mass numbers
*       k [fm^-1]   = (dp) nucleon momentum
*       idist       = (i)  1 = soft+hard
*                          2 = soft (unnormalized)
*                          3 = soft (normalized to 1)
*                          4 = hard (unnormalized)
*                          5 = n0, as defined in [1]
*                          6 = n1, as defined in [1]
*                          (note: 1=2+4)
*
*     References: 
*     [1] C.Ciofi degli Atti and S.Simula, Pys.Rev.C53(1996)1689
*     [2] R.B.Wiringa et al., Phys.Rev.C62(2001)014001
*
*
*  B. DECLARATIONS
*
      implicit none

      double precision rhoCSsh,k
      integer iZ,iA,idist

*    *** variables 

      integer nmax,nin,i,j,n
      save n

      parameter (nmax=20)
      double precision a0(nmax),b0(nmax),c0(nmax),d0(nmax),e0(nmax)
     &     ,f0(nmax),a1(nmax),b1(nmax),b2(nmax),c1(nmax),d1(nmax)
     &     ,n0,n1,soft,hard,norm,tmp1,k2,k4,k6,k8
      integer ZZ(nmax),AA(nmax)
      save a0,b0,c0,d0,e0,f0,a1,b1,b2,c1,d1,ZZ,AA

      character string*100

*    *** functions

      integer nextunit

*    *** constants & initial values

*     Constants
      double precision pi,fourpi,pith
      parameter(pi=3.1415926535898d0,fourpi=4*pi
     &     ,pith=pi*5.5683279968317d0)

      double precision normCS
      parameter(normCS=1d0/(4d0*pi))

*     initialization variables
      logical firsttime
      integer Zold,Aold,nnuke
      save Zold,Aold,firsttime,nnuke
 
      data Zold/0/ Aold/0/ firsttime/.true./
*
*  C. ACTION   
*

*    *** INITIALIZATION 
      if (firsttime) then
         firsttime = .false.
         nin = nextunit()
         open(nin,file='fermimotion2.CS.tbl',status='old')
*       ... skips headers
         do i = 1, 9
            read(nin,*) string
         end do
*       ... reads n0(k) parameters
         j=0
 10      j=j+1
         read(nin,*) string,ZZ(j),AA(j),a0(j),b0(j),c0(j),d0(j)
     &        ,e0(j),f0(j)
         if (ZZ(j).ne.0) goto 10 
         nnuke = j-1
*       ... skips text
         do i = 1, 4
            read(nin,*) string
         end do
*       ... reads n1(k) parameters
         j=0
 20      j=j+1
         read(nin,*) string,ZZ(j),AA(j),a1(j),b1(j),b2(j),c1(j),d1(j)
         if (ZZ(j).ne.0) goto 20 
         nnuke = j-1
         close(nin)
      end if

*    *** if new nucleus, locates it in the tables
      if ((iZ.ne.Zold).or.(iA.ne.Aold)) then
         Zold = iZ
         Aold = iA
         call ilocatetab(ZZ,nnuke,iZ,j)
*       ... if iZ and iA are found in the table OK 
         if ((j.gt.0).and.(j.le.nnuke)) then
            if ((ZZ(j).eq.iZ).and.(AA(j).eq.iA)) then 
               n = j
            else if ((ZZ(j-1).eq.iZ).and.(AA(j-1).eq.iA)) then 
               n = j-1
            else if ((ZZ(j+1).eq.iZ).and.(AA(j+1).eq.iA)) then 
               n = j+1
            else 
*             ... otherwise stops
               print*
               print*, 'ERROR (rhosoftCS): Z,A outside parameter table'
     &              ,iZ,iA
               print*
               stop
            end if
         else 
*          ... otherwise stops
            print*
            print*, 'ERROR (rhosoftCS): Z,A outside parameter table'
     &           ,iZ,iA
            print*
            stop
         end if
      end if


*    *** COMPUTES the Fermi momentum distribution

      k2 = k*k

*    ...  n0 (as defined in [1])
      if (iA.le.4) then
         n0 = a0(n)*dexp(-b0(n)*k2)/((1+c0(n)*k2)**2)
     &        + d0(n)*dexp(-e0(n)*k2)/((1+f0(n)*k2)**2)
      else
         k4 = k2*k2
         k6 = k4*k2
         k8 = k4*k4
         n0 = a0(n)*dexp(-b0(n)*k2) 
     &        * (1d0 + c0(n)*k2 + d0(n)*k4 + e0(n)*k6 + f0(n)*k8)
      end if
      
*    ... soft, hard, and n1 (as defined in [1])
      tmp1 = a1(n)*dexp(-b1(n)*k2)/((1+b2(n)*k2)**2)
      hard = c1(n)*dexp(-d1(n)*k2)
      soft = n0 + tmp1
      n1 = tmp1 + hard


*    *** final result
      if (idist.eq.1) then
*       ... soft+hard 
         rhoCSsh = normCS * (n0 + n1)
      else if (idist.eq.2) then
*       ... soft (unnormalized)
         rhoCSsh = normCS * soft
      else if (idist.eq.3) then
*       ... soft (normalized to 1)
         norm = 1d0 - (c1(n)/fourpi) * (pith/d1(n))
         rhoCSsh = normCS * (soft / norm)
      else if (idist.eq.4) then
*       ... hard (unnormalized)
         rhoCSsh = normCS * hard
      else if (idist.eq.5) then
*       ... n0 [1]
         rhoCSsh = normCS * n0
      else if (idist.eq.6) then
*       ... n1 [1]
         rhoCSsh = normCS * n1
      end if

      return
      end

************************************************************************
*     Search an ordered integer table
      SUBROUTINE ilocatetab(xx,n,x,j) 
*     from "Numerical REcipes in F77"
*
*  A. COMMENTARY
*
*     Given an array xx(1:n)and given a value x returns a value j such 
*     that x is between xx(j) and xx(j+1). xx(1:n) must be monotonic, 
*     either increasing r decreasing.j=0 or j=n is returned to indicate 
*     that x is out f range. 
*
*  B. DECLARATIONS
*
      implicit none

      INTEGER j,n 
      integer x,xx(n) 

*    *** variables

      INTEGER jl,jm,ju 
*
*  C. ACTION
*
      jl=0 
      ju = n+1 
 10   if (ju-jl.gt.1) then 
         jm=(ju+jl)/2 
         if ((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm))) then 
            jl=jm 
         else 
            ju=jm 
         endif 
         goto 10
      endif 
      if (x.eq.xx(1)) then 
         j=1 
      else if (x.eq.xx(n)) then 
         j=n-1 
      else 
         j=jl 
      endif 

      return 
      end

************************************************************************
* INPUT/OUTPUT subroutines                                             *
* v 1.0, May 2005                                                      *
* collected or written by A.Accardi                                    *
*                                                                      *
* Contents:                                                            *
*                                                                      *
*  - NextUnit        Returns unallocated i/o unit                      *
*                                                                      *
************************************************************************

************************************************************************
*     Finds available i/o unit
      integer function NextUnit()
*     Author: CTEQ collab. (taken from "Cteq6Pdf-2004.f")
*
*  A. COMMENTARY
*
*     Returns an unallocated FORTRAN i/o unit between 10 and 300.
*
*  B. DECLARATIONS
*
      implicit none

*    *** variables

      Logical EX
      integer N
*
*  C. ACTION   
*
      Do 10 N = 10, 300
         INQUIRE (UNIT=N, OPENED=EX)
         If (.NOT. EX) then
            NextUnit = N
            Return
         Endif
 10   Continue
      Stop ' There is no available I/O unit. '

      return
      End
      

