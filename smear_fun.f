ckm
      SUBROUTINE smear_fun(ran,charge,beta,torcur,p,thetad,phid,
     $     pout,thout,phiout)
c---- ----------------------------------------------------------
c      Smearing momentum and theta 
c      !common for fastmc and clasev packages! 
c      April 15 2005
c---- ----------------------------------------------------------
c      May 23 2005. K.Mikhailov: add new parameters to conf-files
c      for the momentum resolution depedence on theta
c---- ----------------------------------------------------------
c      May 28 2005. K.Mikhailov: add new parameters to conf-files
c      for the forward phi-resolution depedence on theta
c      Change the coge according Gail's phi-resolution study:
c sigma_1 = 0.0142/sin(theta) in units of mrad (where theta is the polar angle)
c sigma_1=sigma_1/sin(12^0)=sigma_1/0.2079=0.0142/0.2079/sin(theta)=
c         0.0683/sin(theta)[mrad]=0.683E-4/sin(theta)[rad]
c sigma_2 = a0 + a1*theta + a2*theta^2   where theta is in degrees
c       a0 = 530.6 mrad-MeV; a1= -5.38 mrad-MeV/deg; a2=0.0898 mrad-MeV/deg2
c       a0 = 0.5306 rad-MeV; a1= -0.00538 rad-MeV/deg; a2=0.0000898 rad-MeV/deg2
c       a0 = 0.0005306 rad-GeV; a1= -0.00000538 rad-GeV/deg; a2=0.0000000898 rad-GeV/deg2
c       a0 = 0.5306E-3 rad-GeV; a1= -0.538E-5 rad-GeV/deg; a2=0.898E-7 rad-GeV/deg2
c---- ----------------------------------------------------------
c      May 31 2005. K.Mikhailov: add new parameters to conf-files
c      The new results for the forward phi resolution (vertex+svt+dc)
c      sigma_1 = a0/theta + a1 + a2*theta where theta is in degrees
c	    a0 = 3.125; a1 = 0.0262; a2 = -0.00035
c	    a0 = 0.3125E-2 rad; a1 = 0.262E-4 rad; a2 = -0.35E-6 rad
c
c      sigma_2 = a0 + a1*theta + a2*theta^2 where theta is in degrees
c           a0 = 541.86; a1 = -5.307; a2 = 0.0886
c           a0 = 0.542E-3 rad-GeV; a1 = -0.531E-5 rad-GeV; a2 = 0.886E-7rad-Gev
c
c
c
c---- ----------------------------------------------------------

      implicit none
      real ran(3),charge,torcur,p,phid,pout,thout,phiout
      real thetad,thr,phir
      real sig1p,sig2p,sig1x,sig2x,sig1fi,sig2fi,beta
      real sig1th,sig2th
      real dmom,dfi,dtet,dpperp
      logical inbend,outbend,neutral,fwd,central
      real pnorm,pperp,pperpout
      real dthe,dphi
      integer lrun
      data lrun/3/
c common block containing track_type; defined in sub. tracktype_
      real d2r/0.0174532925199433/,r2d/57.29577951308232/
      common/track_type/inbend,outbend,neutral,fwd,central
      real theta_ic_max,theta_ic_min
      parameter (theta_ic_max=5.0,theta_ic_min=1.6)
      real erg, ic_enrg_res, ic_ang_res
      real p_ic_erg(4), p_ic_ang(4)
      data p_ic_erg / 0.02997, -0.007688, 0.001068, -4.636e-05 /
      data p_ic_ang / 0.0572505, -0.0144679, 0.00204576, -9.26138e-05  /
      include "read_par_clas12g.inc"
c
c  since our acceptance functions care about sagitta and not momentum, per se,
c  we will form a "normalized momentum", pnorm
cmac   this normalized momentum is only used for the fwd. tracker
      pnorm=p*Tmax/abs(torcur)
c
      thr=thetad*d2r
      phir=phid*d2r
c     select which of the 6 functions we want
      if(inbend.and.fwd)then
ckm 14-apr-2005 (fit on all thetas), theta, phi resolutions:
         sig1p=Tmax/torcur*
     &        (FoMT1S1+FoMT1S2*thetad+FoMT1S3*thetad*thetad)
         sig2p=Tmax/torcur*
     &        (FoMT2S1+FoMT2S2*thetad+FoMT2S3*thetad*thetad)
ckm   sigma_theta is at 5 degrees(as upper limit)
         sig1th=FoInThS1
         sig2th=FoInThS2
ckm May 28, 2005: Gail's phi-resolution as function of theta
ckm         sig1fi=2.*sig1th
ckm         sig2fi=2.*sig2th

ckm         sig1fi=FoIOFiS1/SIN(thr)
c         print*,'sig1fi,FoIOFiS1,thr',sig1fi,FoIOFiS1,thr*r2d
ckm         sig2fi=FoA0FiS2+FoA1FiS2*thr*r2d+FoA2FiS2*(thr*r2d)**2!theta in degrees

ckm May 31, 2005
         sig1fi=FoA0FiS1/(thr*r2d) + FoA1FiS1 + FoA2FiS1*(thr*r2d)!theta in degrees
         sig2fi=FoA0FiS2+FoA1FiS2*thr*r2d+FoA2FiS2*(thr*r2d)**2!theta in degrees

         pout=p+ran(1)*p*sqrt((sig1p*p)**2+(sig2p/beta)**2)
         thout=thr+ran(2)*sqrt(sig1th**2+(sig2th/p/beta)**2)
         phiout=phir+ran(3)*sqrt(sig1fi**2+(sig2fi/p/beta)**2)
c     
      elseif(central.AND..NOT.neutral)then
         sig1p=CeInMoS1
         sig2p=CeInMoS2
ckm May 31, 2005: sigma_2 = a0 + a1*theta + a2*theta^2
         sig1th=CeInThS1*(sin(thr))**2
         sig2th=CeInThS2*(sin(thr))**2
cmac         sig1th=CeInThS1*(sin(thr))**2
cmac         sig2th=CeA0ThS2 + CeA1ThS2*(thr*r2d) + CeA2ThS2*(thr*r2d)**2

         sig1fi=CeInFiS1
         sig2fi=CeInFiS2
         pperp=p*sin(thr)
cmac         dpperp=p*sqrt((CeInMoS1*p)**2+(CeInMoS2/beta)**2)
         dpperp=pperp*sqrt((CeInMoS1*p)**2+(CeInMoS2/beta)**2)
         pperpout=pperp+ran(1)*dpperp
         thout=thr+ran(2)*sqrt(sig1th**2+(sig2th/p/beta)**2)
         pout=pperpout/sin(thout)
         phiout=phir+ran(3)*sqrt(sig1fi**2+(sig2fi/p/beta)**2)
c     
      elseif(outbend.and.fwd)then
ckm 14-apr-2005 (fit on all thetas), theta, phi resolutions:
         sig1p=Tmax/torcur*
     &        (FoMT1S1+FoMT1S2*thetad+FoMT1S3*thetad*thetad)
         sig2p=Tmax/torcur*
     &        (FoMT2S1+FoMT2S2*thetad+FoMT2S3*thetad*thetad)
ckm   sigma_theta is at 5 degrees(as upper limit)
         sig1th=FoInThS1
         sig2th=FoInThS2
ckm May 28, 2005: Gail's phi-resolution as function of theta
ckm         sig1fi=2.*sig1th
ckm         sig2fi=2.*sig2th

ckm         sig1fi=FoIOFiS1/SIN(thr)
c         print*,'sig1fi,FoIOFiS1,thr',sig1fi,FoIOFiS1,thr*r2d
ckm         sig2fi=FoA0FiS2+FoA1FiS2*thr*r2d+FoA2FiS2*(thr*r2d)**2!theta in degrees
 
ckm May 31, 2005
         sig1fi=FoA0FiS1/(thr*r2d) + FoA1FiS1 + FoA2FiS1*(thr*r2d)!theta in degrees
         sig2fi=FoA0FiS2+FoA1FiS2*thr*r2d+FoA2FiS2*(thr*r2d)**2!theta in degrees

         pout=p+ran(1)*p*sqrt((sig1p*p)**2+(sig2p/beta)**2)
         thout=thr+ran(2)*sqrt(sig1th**2+(sig2th/p/beta)**2)
         phiout=phir+ran(3)*sqrt(sig1fi**2+(sig2fi/p/beta)**2)
c     
      elseif(neutral.and.fwd)then
         if( thetad .GT. theta_ic_max ) then 
c  Neutral in FWD EC            
            pout=p+ran(1)*sqrt(p)*FoNeMoS1
            thout=thr+ran(2)*FoNeThS1
            phiout=phir+ran(3)*FoNeFiS1
         elseif( thetad .GT. theta_ic_min .AND. 
     &           thetad .LT. theta_ic_max ) then 
c  Neutral in IC 
            if( beta .GT. 0 .AND. beta .LT. 1 ) then 
               erg = p / beta
            else
               erg = p
            endif
c Find IC energy esolution from Rustam parameterization

            ic_enrg_res = p_ic_erg(1)*erg +  p_ic_erg(2)*erg**2 + 
     +           p_ic_erg(3)*erg**3 + p_ic_erg(4)*erg**4
c Find IC angular esolution from Rustam parameterization 
c Assume that theta and phi resolutions are same 
            ic_ang_res  = p_ic_ang(1) + p_ic_ang(2)*erg + 
     +         p_ic_ang(3)*erg**2 + p_ic_ang(4)*erg**3  
c Smear energy and angles with a gaussian            
            pout = p + ran(1) * ic_enrg_res
            thout  = thr + ran(2) * ic_ang_res * d2r
            phiout = phir+ ran(3) * ic_ang_res * d2r
         endif
c     
      elseif(neutral.and.central)then
         pout=p+ran(1)*sqrt(p)*CeNeMoS1
         thout=thr+ran(2)*CeNeThS1
         phiout=phir+ran(3)*CeNeFiS1
      endif
c  convert from radians to degrees
      thout=thout*r2d
      phiout=phiout*r2d


c     
      return
      end
c     

