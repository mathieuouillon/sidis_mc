c+++++++++++++++++++++++++++++++++++++++++++++++
      subroutine accep_fun(id,torcur,p,thetad,phi0,weight,d_phi)
cOK
      implicit none
      include "read_par_clas12g.inc"
      integer id
      real torcur,p,thetad,phi0,weight,d_phi
      real phis
c
c  IC acceptance (to be added in the conf file)
c
      real theta_ic_max,theta_ic_min,coil_calo_width
      parameter (theta_ic_max=5.0,theta_ic_min=1.6)  ! 5 and 1.6 degree
      parameter (coil_calo_width=0)  ! coil calo width 0-default
c     parameter (coil_calo_width=6)  ! coil calo width 6-if included
c
c common block containing track_type; defined in sub. tracktype_
      logical inbend,outbend,neutral,fwd,central
      common/track_type/inbend,outbend,neutral,fwd,central

      real pnorm,thcut,expon,d2r,sinth,costh
      real sinph,cosph,r,ec_r,ec_th,dl,xh,yh,zh,uh,vh,wh
      real xi,yi,rcoil,rc0,hcoil,dcoil,ylow,yhi,xcut,tgrho,dth
c
      real theta_min_cut,p_pion_min,el_e_min
      data theta_min_cut,p_pion_min,el_e_min/5.0,0.2,0.5/
      data ec_r,ec_th /510.32,0.436332/
      data ylow,yhi,tgrho,dl,rc0/-182.974,189.956,1.95325,100.,120./
c

      real r2d,dgap
      logical pcutl
      data d2r,r2d,dgap/0.01745,57.299578,2./

      real deflPhi, c_spec
      data c_spec/ 0.299792458/
       
      real charge, rmass, Enrg, BdL
      real findBdL
c
c  Take into account phi deflection due to 
c  solenoid magnet
c
      call charge_type_( id, charge, rmass )
      Enrg = sqrt( p*p + rmass*rmass )
      if ( abs( p * sin( thetad*d2r ) ) .GT. 0  ) then 
         BdL =  findBdL( p, Enrg, thetad*d2r )
         deflPhi = r2d * c_spec * BdL / ( p * sin( thetad * d2r ) )
      else
         deflPhi = 0.0
      endif
      if( charge .GT. 0 ) then 
         phis = amod( ( phi0 + deflPhi ) + 30.0 + 60000.0, 60.0 ) - 30.0 
      elseif( charge .LT. 0 ) then 
         phis = amod( ( phi0 - deflPhi ) + 30.0 + 60000.0, 60.0 ) - 30.0
      else 
         phis = phi0
      endif
c      print *,'my version',charge,phi0,phis,phis-phi0

c
c  since our acceptance functions care about sagitta and not momentum, per se,
c  we will form a "normalized momentum", pnorm
c
      pnorm=p*Tmax/abs(torcur)
c
c  initialize values
      weight=0.
      d_phi=0.
c
c     hard cuts from old version
c
      if(thetad.lt.theta_min_cut.and.id.ne.22)return
      if(id.eq.11.and.p.lt.el_e_min)return
      if((abs(id).eq.211.or.abs(id).eq.321)
     6.and.p.lt.p_pion_min)return
c
c    negatives
c
      IF(id.eq.11.or.id.eq.-211.or.id.eq.-321.or.id.eq.-2212)THEN
*     Define acceptance cuts for electrons
         if(thetad.gt.theta_min_cut.and.fwd)then
            If(outbend)then
               thcut = Thacc01 - 25.*(1.-p*1500./abs(torcur)/5.)**24
               pnorm = p*tmax/abs(torcur)
               pcutl=pnorm.gt.(Thmom1+Thmom2*(thetad/Thmom3)**Thmom4)
     6                        *abs(torcur)/tmax
               if(thetad.gt.thcut.and.pcutl)then
            d_phi=Phimom1*sin((thetad-Phimom2)*d2r)**Phimom3* 
     *           (p*1500./abs(torcur)/5.)**Phimom5
c
                  if(abs(phis).lt.d_phi)then 
                     weight=1.
                  else
                     weight=0.
                  endif
               else
                  weight=0.
                  d_phi=0.
               endif
c
            ElseIf(inbend)then
               thcut =Thacc1 + Thacc2/((p+0.05)*tmax/abs(torcur))
               if(thetad.gt.thcut)then
                  expon = Thacc4*(pnorm)**Thacc5
                  d_phi = Thacc3*sin((thetad-thcut)*d2r)**expon 
                  if(id.eq.-211) d_phi=d_phi*Thacc6
                  if(abs(phis).lt.d_phi)then
                     weight=1.
                  else
                     weight=0.
                  endif
               else
                  d_phi= 0.
                  weight=0.
               endif
            Endif
c
         elseif(thetad.gt.Thf_lim+dgap
     6   .and.((id.eq.-211).or.(id.eq.-321).or.(id.eq.-2212)))then
            weight=1.
            d_phi=30.
         endif
c     
         d_phi=d_phi/30.
c
c     End of simulation for negatives
c
c     Positives come here
c
      ELSEIF(id.eq.2212.or.id.eq.211.or.id.eq.321.or.id.eq.45)THEN
*     minimum momentum cut
         if(p.gt.0.4.and.fwd.and.thetad.gt.theta_min_cut)then
            If(inbend)then
               thcut = 5.-17.*(1.-p*1500./abs(torcur)/5.)**24
               pnorm = p*tmax/abs(torcur)
               pcutl=pnorm.gt.(.6+4.*(thetad/45.)**6)*abs(torcur)/tmax
               if(thetad.gt.thcut.and.pcutl)then
                  d_phi=25.*sin((thetad-5.)*d2r)**0.15* 
     *                 (p*1500./abs(torcur)/5.)**(1./24.)
               else
                  d_phi=0.
               endif
            ElseIf(outbend)then
               d_phi = 28.*cos((thetad-37.)*d2r)*
     *              (p*tmax/abs(torcur)/5.)**(1./64.)
            EndIf
            if(abs(phis).lt.d_phi) weight=1.
c
         elseif(thetad.gt.Thf_lim+dgap.and.thetad.lt.100.
     &           .and.(p.gt.(0.25/sin(thetad*d2r))).and.p.lt.2.)then
            d_phi=30.
            weight=1.
         else
            d_phi=0.
            weight=0.
         endif 
         d_phi=d_phi/30.
c
c     End of simulation for positives
c
      ELSEIF(id.eq.22)THEN
c     Photon
         if(p.lt.0.1)return
         dth=atan(2./r/sqrt(p))*r2d
c
c        Central calorimeter 
         if(thetad.gt.Thf_lim+dgap.and.thetad.lt.100.)then
C        Comment out these two  line if you want central EC            
            d_phi = 0.0
            weight = 0.0
c        Uncomment these two lines if you want central EC    
c            d_phi=28.
c            weight=1.
c
c        Forward EC            
         elseif(thetad.lt.Thf_lim.and.thetad.gt.theta_ic_max)then
c     
            if(p.lt.0.2)return  ! check this
c     
            sinth=sin(thetad*d2r)
            costh=cos(thetad*d2r)
            sinph=sin(phis*d2r)
            cosph=cos(phis*d2r)
c     
            r=(ec_r*(sin(ec_th)**2+cos(ec_th))+dl)/
     /           (costh+sinth*cosph*sin(ec_th))
            xh=r*cosph*sinth
            yh=r*sinph*sinth
            zh=r*costh-dl
c     
            call ec_xyz_duvw(xh,yh,zh,uh,vh,wh,xi,yi)
c     
            rcoil=dl*costh+sqrt(rc0**2-(dl*sinth)**2)
            hcoil=rcoil*sinth*cosph
            dcoil=hcoil*sin(0.5236)/cos(0.5236)-abs(rcoil*sinth*sinph)
c     
            if(yi.gt.ylow+10..and.yi.lt.yhi-5.)then
               xcut=(yi-ylow)/tgrho
               if(abs(xi).lt.xcut-5.)then
                  weight=1
               elseif(abs(xi).gt.xcut.and.dcoil.lt.coil_calo_width)then
                  weight=1
               endif
               d_phi=1.-5./(abs(r*sinph*sinth)+8.)
            endif
c
c      only IC part
c
         elseif(thetad.gt.theta_ic_min.and.thetad.le.theta_ic_max) then
            weight=1.0
            d_phi=30.0
         else
            weight=0.
            d_phi=0.
       endif
      ELSE
c     
cc     No acceptance and smearing for this particle
         weight=0.
         d_phi=0.
      ENDIF
*     
      return
      end


C
C
C=================================================================================
C
C A function to calculate BdL of the solenoid 
C assuming uniform magnetic field of solenoid 
C specified in the pCDR of CLAS12
C  
      real function findBdL( P, E, Theta )
c     Input: momentum, energy and the polar angle 
c     of a charged particle in GeV and rad
c     Output: calculated BdL in T x m 
      implicit none

      real P, E, Theta
      real Time, TimeL, TimeR

c     B-field, diameter and legth of the solenoid from Vlassov plots
      real B, dBarrel, lBarrel 
      data B, dBarrel, lBarrel / 5.0, 0.78, 0.9 /
      real rBarrel

c     spped of light in m/ns 
      real cLight
      data cLight / 0.299792458 /

      real v, v_t, v_l

      real P_t, r_curv

c     Find the velocity and the components  
      if ( E .GT. 0 .AND. P .LE. E ) then  
         v = cLight * P / E
      else
         write( 6, * ) 'Error: Bad kinematics'
         v = 0.0
      endif
      v_l = v * cos( Theta )
      v_t = v * sin( Theta )
      
      P_t = P * sin( Theta )

c     find the radious of the curvature in the solenoid 
      if( abs( B ) .GT. 0 ) then 
         r_curv = P_t / ( cLight * B )
      else
         r_curv = 1.0e+15
      endif

c     Find the time it would take to leave the solenoid 
c     through the base of the cylinder
      if ( v_l .GT. 0 ) then 
         TimeL = lBarrel / 2.0 / v_l
      else
         TimeL = 1.0e+10
      endif
c     find te time it would take to leave the solenoid
c     through the walls of the cylinder
      rBarrel = dBarrel / 2.0
      if ( ( rBarrel .LT. 2. * r_curv ) .AND. ( v_t .GT. 0 ) ) then 
         TimeR = 2. * r_curv * asin( rBarrel / (2.*r_curv) ) / v_t
      else 
         TimeR = TimeL   
      endif
c     Take the minimum of the two
      Time = min( TimeL, TimeR )

c     Final value
      findBdL = B * v_t * Time
      
      return
      end
