!c------------------------------------------------------------------------------
!c Subroutine for lorentz transformation
!c------------------------------------------------------------------------------
      subroutine TL(EEj,PPj,Pjx,Pjy,Pjz,BB,Bx,By,Bz)
      implicit none

      real EEi,PPi,Pix,Piy,Piz ! intermediate kinematics
      real EEj,PPj,Pjx,Pjy,Pjz ! jntermedjate kjnematjcs
      real BB,Bx,By,Bz,GG ! Lorentz kinematics

      GG = 1/sqrt(1-BB**2)

      EEi = GG*(EEj-Bx*Pjx-By*Pjy-Bz*Pjz)
      Pix = -Bx*GG*EEj+Pjx+Pjx*(GG-1)*Bx**2/BB**2+Pjy*(GG-1)*By*Bx/BB**2+Pjz*(GG-1)*Bz*Bx/BB**2
      Piy = -By*GG*EEj+Pjx*(GG-1)*Bx*By/BB**2+Pjy+Pjy*(GG-1)*By*By/BB**2+Pjz*(GG-1)*Bz*By/BB**2
      Piz = -Bz*GG*EEj+Pjx*(GG-1)*Bx*Bz/BB**2+Pjz+Pjy*(GG-1)*By*Bz/BB**2+Pjz*(GG-1)*Bz*Bz/BB**2
      PPi = sqrt(Pix**2 + Piy**2 + Piz**2)
      
!c      write(*,*) 'Beta :',BB,Bx,By,Bz,GG
!c      write(*,*) 'Int :',EEi,PPi,Pix,Piy,Piz

      EEj = EEi
      PPj = PPi
      Pjx = Pix
      Pjy = Piy
      Pjz = Piz
      
      end

!c------------------------------------------------------------------------------
!c Subroutine for rotation around y
!c------------------------------------------------------------------------------
subroutine InitRotY(Theta)
    use kinematics_module
    use file_names_module
    use event_info_module
    use particles_module
    use fermi_motion_module
    implicit none

!ccccc Include all the common blocks
      include 'common.f90'
      real Pix,Piy,Piz ! intermediate kinematics
      real Theta !angle for rotation

      Pix = Pex*cos(Theta) + Pez*sin(Theta)
      Piy = Pey
      Piz = -Pex*sin(Theta) + Pez*cos(Theta)
      Pex = Pix
      Pey = Piy
      Pez = Piz

      Pix = Pnx*cos(Theta) + Pnz*sin(Theta)
      Piy = Pny
      Piz = -Pnx*sin(Theta) + Pnz*cos(Theta)
      Pnx = Pix
      Pny = Piy
      Pnz = Piz

      end

!c------------------------------------------------------------------------------
!c Subroutine for rotation around z
!c------------------------------------------------------------------------------
      subroutine InitRotZ(PhiAng)
    use kinematics_module
    use file_names_module
    use event_info_module
    use particles_module
    use fermi_motion_module
      implicit none

!ccccc Include all the common blocks
      include 'common.f90'
      real Pix,Piy,Piz ! intermediate kinematics
      real PhiAng !angle for rotation

      Pix = Pex
      Piy = Pey*cos(PhiAng) + Pez*sin(PhiAng)
      Piz = -Pey*sin(PhiAng) + Pez*cos(PhiAng)
      Pex = Pix
      Pey = Piy
      Pez = Piz

      Pix = Pnx
      Piy = Pny*cos(PhiAng) + Pnz*sin(PhiAng)
      Piz = -Pny*sin(PhiAng) + Pnz*cos(PhiAng)
      Pnx = Pix
      Pny = Piy
      Pnz = Piz

      end
!c------------------------------------------------------------------------------
!c Subroutine for rotation around y
!c------------------------------------------------------------------------------
subroutine FinalRotY(Theta)
    use kinematics_module
    use file_names_module
    use event_info_module
    use particles_module
    use fermi_motion_module
    implicit none

!ccccc Include all the common blocks
      include 'common.f90'
      real Pix,Piy,Piz ! intermediate kinematics
      real Theta !angle for rotation
      integer ip ! For do

      do ip=1,N

        Pix = p(ip,1)*cos(Theta) + p(ip,3)*sin(Theta)
        Piy = p(ip,2)
        Piz = -p(ip,1)*sin(Theta) + p(ip,3)*cos(Theta)
      
        p(ip,1) = Pix
        p(ip,2) = Piy
        p(ip,3) = Piz

      enddo      

      end

!c------------------------------------------------------------------------------
!c Subroutine for rotation around z
!c------------------------------------------------------------------------------
      subroutine FinalRotZ(PhiAng)
    use kinematics_module
    use file_names_module
    use event_info_module
    use particles_module
    use fermi_motion_module
      implicit none

!ccccc Include all the common blocks
      include 'common.f90'
      real Pix,Piy,Piz ! intermediate kinematics
      real PhiAng !angle for rotation
      integer ip ! For do

      do ip=1,N

        Pix = p(ip,1)
        Piy = p(ip,2)*cos(PhiAng) + p(ip,3)*sin(PhiAng)
        Piz = -p(ip,2)*sin(PhiAng) + p(ip,3)*cos(PhiAng)
     
        p(ip,1) = Pix
        p(ip,2) = Piy
        p(ip,3) = Piz

      enddo      

      end
