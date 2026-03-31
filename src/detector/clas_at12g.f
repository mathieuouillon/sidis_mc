c==================================================
      subroutine clas_at12g(id,p,thetad,phis,torcur,d_phi,weight,cfile)
*     
*     Version 1.1 / April 18, 2001
*     Simulates acceptance and smearing for CLAS12
*
*     Nov. 30, 2004 - major changes as I separate functions into
*     separate subroutines - Mestayer
*
* Inputs -
*     ic      Configuration number (read configuration from-> conf[i].dat)
*     id:     Particle ID according to PDG (11-electron, 2212-proton, 
*     211-pi+, -211-pi-, 22-photon ...)
*     p:      momentum, in GeV/c, on output will be smeared
*     phis:   phi in the sector in degrees (at mid plane phi=0.), 
*                 on output will be smeared
*     thetad: scattering angle in degrees, 
*                 on output will be smeared
*     torcur: torus current in Amps (e.g. -1500.)
* Outputs -
*     d_phi:  part of the phi acceptance, full acceptance should be 
*                 delta_phi=2*pi*d_phi
*     weight: =1. if accepted,  0. if not

c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ckm New version of clas_at12g with reading parameters from data file
c___________________________________________________________________
      implicit none
      include "includes/read_par_clas12g.inc"
      real p,thetad,phis,torcur,d_phi,weight
      real pout,thout,phiout,beta
      integer ic,id
c     
      real phie_mod,thcut,pnorm,expon
      logical pcutl
      real ec_r,ec_th,rl,rh
      data ec_r,ec_th,rl,rh/510.32,0.436332,562.22,215.67/
      real ylow,yhi,tgrho,dl,rc0
      data ylow,yhi,tgrho,dl,rc0/-182.974,189.956,1.95325,100.,120./
      real r,xh,yh,zh,xcut
      real uh,vh,wh,xi,yi
      real sinth,costh,sinph,cosph,dth
      real rcoil,hcoil,dcoil
c     
ckm      real dthe,dphi,tmax,rmass
      real dthe,dphi,rmass
      data dthe,dphi/0.0286,0.0573/
      real d2r,pi,r2d,dgap,rran(3)
      data d2r,pi,r2d,dgap/0.01745,3.1415926,57.299578,2./
      integer lrun
      data lrun/3/
      real charge
c  common block keeping the logical values of the track type;      
      logical inbend,outbend,neutral,fwd,central
      common/track_type/inbend,outbend,neutral,fwd,central
c
      logical lfirst/.true./
      logical linit/.true./
      character*30 cfile
c
c
      if(lfirst) then
         call init_fastmc(cfile)
         lfirst=.false.
         linit=.false.
c         pause                  !check parameters
      endif
c  initialize weight and d_phi     
      weight=0.
      d_phi=0.
c
c  1st decide if particle charge is +1, -1, 0 or unknown
c
      call charge_type_(id,charge,rmass)
      beta=1.
      if(rmass.ne.999..and.rmass.ne.0.)then
         beta=p/sqrt(p**2+rmass**2)
      endif
c      print *,'id, charge,beta',id,charge,beta
c  a charge type of 999. is an unknown
      if(charge.eq.999.) then
         print *,'charge is 999. - unknown particle'
         return
      endif
c
c  2nd, determine which track type we are concerned with;
c  the subroutine tracktype_ fills the logical variables in common
c  block, track_type
c


      call tracktype_(charge,torcur,thetad)
c      print *,'in,out,neu,fwd,cent',inbend,outbend,neutral,fwd,central
c
c  3rd, call acceptance function with charge, torus current,
c  and p,th,phi as inputs; return is weight (1, accepted; 0, not)
c  and d_phi which is the phi part of the acceptance.  Note, this
c  is the tracking acceptance only, it does NOT contain the additional
c  fiducial volume cuts of Cerenkov counters imposed on electrons.
c
c  get a 3-dim. random number
      call rnorml(rran,lrun)
c
      if(abs(charge).eq.1.or.charge.eq.0)then
       call accep_fun(id,torcur,p,thetad,phis,weight,d_phi)
c         print *,'accep - weight,dphi',weight,d_phi
      endif
c
c  4th, if particle is accepted (weight=1.) call smearing routine
c   with input momenta and output momenta
c
      if(weight.eq.1.) then
         call rnorml(rran,lrun)
c
         call smear_fun(rran,charge,beta,torcur,p,thetad,phis,pout,
     $        thout,phiout)
c  now set p,th,phi to smeared values
         p=pout
         thetad=thout
         phis=phiout
      endif
c     
      return
      end
C=======================================================================
      Subroutine ec_xyz_duvw(x,y,z,u,v,w,xi,yi)
C=======================================================================
      implicit none
      REAL x,y,z,u,v,w,xi,yi,zi
      REAL EC_the,ec_phi,phi,ylow,yhi,tgrho,sinrho,cosrho
      data EC_the/0.4363323/
      data ylow,yhi/-182.974,189.956/
      data tgrho,sinrho,cosrho/1.95325,0.8901256,0.455715/
      real rot(3,3)
c
      phi=atan2(y,x)*57.29578
      if(phi.lt.0.)phi=phi+360.
      phi=phi+30.
      if(phi.ge.360.)phi=phi-360. 	
      Ec_phi=int(phi/60.)*1.0471975
c     
      rot(1,1)=cos(Ec_the)*cos(Ec_phi)
      rot(1,2)=-sin(Ec_phi)
      rot(1,3)=sin(Ec_the)*cos(Ec_phi)
      rot(2,1)=cos(Ec_the)*sin(Ec_phi)
      rot(2,2)=cos(Ec_phi)
      rot(2,3)=sin(Ec_the)*sin(Ec_phi)
      rot(3,1)=-sin(Ec_the)
      rot(3,2)=0.
      rot(3,3)=cos(Ec_the)
c     
      yi=x*rot(1,1)+y*rot(2,1)+z*rot(3,1)
      xi=x*rot(1,2)+y*rot(2,2)+z*rot(3,2)
      zi=x*rot(1,3)+y*rot(2,3)+z*rot(3,3)
      zi=zi-510.32
      u=(yi-ylow)/sinrho
      v=(yhi-ylow)/tgrho-xi+(yhi-yi)/tgrho
      w=((yhi-ylow)/tgrho+xi+(yhi-yi)/tgrho)/2./cosrho
      end
C++++++++++++++++++++++++++++++++++++++++
      subroutine charge_type_(id,charge,rmass)
C++++++++++++++++++++++++++++++++++++++++
c subroutine which returns charge of a particle given its PDG id
      implicit none
      real amp,ampi,amka
      parameter (amp=0.93827,ampi=0.13957,amka=0.49368)
      real amn
      parameter (amn=0.93957)
      real charge,rmass
      integer id
c
      if(id.eq.11)then
         charge=-1.
         rmass=0.5e-03
      elseif(id.eq.-11) then
         charge=1.
         rmass=0.5e-03
      elseif(id.eq.-211) then
         charge=-1.
         rmass=ampi
      elseif(id.eq.-321) then
         charge=-1.
         rmass=amka
      elseif(id.eq.321) then
         charge=1.
         rmass=amka
      elseif(id.eq.2212)then
         charge=1.
         rmass=amp
	 elseif(id.eq.-2212)then
         charge=-1.
         rmass=amp         
      elseif(id.eq.211)then
         charge=1.
         rmass=ampi
      elseif(id.eq.45)then
         charge=1.
         rmass=0.5
      elseif(id.eq.22)then
         charge=0.
         rmass=0.
      elseif(id.eq.2112)then
         charge=0.
         rmass=amn
      else   
         print *,'unsupported ID ->',id
         charge=999.
         rmass=999.
      endif   
      return
      end
c+++++++++++++++++++++++++++++++++++++++++++++++
C++++++++++++++++++++++++++++++++++++++++
      subroutine tracktype_(charge,torcur,thetad)
*OK
c  subroutine which figures out the track type
c  common block keeping the logical values of the track type;      
      implicit none
      real charge,torcur,thetad
      logical inbend,outbend,neutral,fwd,central
      common/track_type/inbend,outbend,neutral,fwd,central
      include "includes/read_par_clas12g.inc"
c  initialize logicals to false
      inbend=.false.
      outbend=.false.
      neutral=.false.
      fwd=.false.
      central=.false.
c  6 functions: (inbend,outbend,neutral) track X (fwd,central) det.
c
      if(charge.eq.0) neutral=.true.
      if(abs(charge).eq.1)then
         if(charge*torcur.gt.0.)outbend=.true.
         if(charge*torcur.lt.0.)inbend=.true.
      endif
      if(thetad.lt.Thf_lim)fwd=.true.   !Thf_lim max forward angle
      if(thetad.ge.Thf_lim)central=.true.
c
      return
      end
