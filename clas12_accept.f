c=================================================
c  Code to read the output of Pepsi and send
c  it to CLAS fastMC
c  By Koko Sept 3 2008 
c  PAC34 preparation
c  Modified to fit in PyQM
c=================================================


         integer function clas12_accept(nat2,k1,k2,k3)
         implicit none
         
         
         integer nat2
                  
         character*30 cl_conf

         real r2d,d2r
         real newphi11,PolarTheta,Azimuthalphi
         real thetad, phis, hp, wt, dphi
         real torcur, thetad_b, phis_b
         double precision k1, k2, k3
         real kin1,kin2,kin3

         parameter (torcur = 2250.)
         data r2d,d2r/57.2957795,0.017453293/
         real pi
         
         pi = acos(-1.0)

         kin1 = real(k1)
         kin2 = real(k2)
         kin3 = real(k3)

         cl_conf = 'datafiles/clasfm/conf5.dat'
         
         hp = sqrt(kin1**2 + kin2**2 + kin3**2)
         thetad = PolarTheta(kin1,kin2,kin3)
         thetad = thetad * r2d
         phis = Azimuthalphi(kin1,kin2)
         
         phis = phis*r2d
         phis = newphi11(phis)

         thetad_b = thetad
         phis_b = phis

         if(abs(nat2).eq.11.or.abs(nat2).eq.211.or.
     &         abs(nat2).eq.321.or.abs(nat2).eq.2212.or.
     &  nat2.eq.2112.or.nat2.eq.22) then


         call clas_at12g(nat2,hp,thetad,phis,torcur,dphi,wt,
     &  cl_conf)
         
         else
         
         thetad_b = 400.
         phis_b = 400.
         thetad = 400.
         phis = 400.         
         dphi = 400.
         wt = 400.
         
         endif

         clas12_accept = wt

         end
         
c========================================================
      real Function newphi11(phi)
      real phi,phinew
      if (phi.gt.330.) then
        phinew = phi-360.
      elseif (phi.ge.0.0.and.phi.le.30.) then
        phinew = phi
      elseif (phi.gt.30.0.and.phi.le.90.) then
        phinew = phi-60.
      elseif (phi.gt.90.0.and.phi.le.150.) then
        phinew = phi-120.
      elseif (phi.gt.150.0.and.phi.le.210.) then
        phinew = phi-180.
      elseif (phi.gt.210.0.and.phi.le.270.) then
        phinew = phi-240.
      elseif (phi.gt.270.0.and.phi.le.330.) then
        phinew = phi-300.
      endif
c       print *,phi,phinew
      newphi11=phinew
      end  
c========================================================
c
        real function  PolarTheta(vx,vy,vz)
        implicit none
        real vx,vy,vz,pmod,theta
        pmod=vx*vx+vy*vy+vz*vz
        if(pmod .gt. 0 ) then
         theta=acos(vz/sqrt(pmod))
        else
         theta=-100
        endif
        PolarTheta=theta
        return
        end
c=======================================================
        real function Azimuthalphi(vx,vy )
        implicit none
        real pi
        parameter ( PI=3.1415926)
        real vx,vy,pmod,phi,cosf
        pmod=vx*vx+vy*vy
        if(pmod .gt. 0 ) then
         pmod=sqrt(pmod)
         cosf=vx/pmod
        else
         cosf=1.0
        endif
        if(abs(cosf) .le. 1.0) phi=acos(cosf);
        if(vy .lt. 0.0) phi= 2*PI-phi;
        Azimuthalphi=phi
        return
        end
         
         
         
