ccccc Generate the interaction position
      subroutine InterPos
      implicit none
      include 'common.f'

      real r
      integer i

ccccc Randomize the radius and angles
      pos_radius = ranf(0)*quantity_table(2000)
      pos_theta = acos(2*ranf(0)-1)
      pos_phi = 2*pi*ranf(0)

      r = 0.
      i = 1.

      do while (pos_radius.ge.quantity_table(i))
        r = r + step_size_dens
        i = i + 1
      enddo
      pos_radius = r + ranf(0)*step_size_dens

cccc Calculate position on cartesian axis
      x_inter = pos_radius * sin(pos_theta) * cos(pos_phi)
      y_inter = pos_radius * sin(pos_theta) * sin(pos_phi)
      z_inter = pos_radius * cos(pos_theta)

      end

ccccc Generate the table of density in function of r
      subroutine GenNucDens
      implicit none

      include 'common.f'

      integer i,idist,irho
      double precision nucdens,r
      real integral

ccccc Some init
      QW_nb = 0
      QW_qhat = 0.

      integral = 0.
      idist = iDens
      irho = 1
      step_size_dens = 0.01
      r = 0.005

      do i=1,2000
        density_table(i) = nucdens(r,iZ,iA,idist,irho)
        r = r + step_size_dens
ccccc calculate the quantity of mater in function of the radius
        integral = integral + density_table(i)*r**3
        quantity_table(i) = integral
      enddo

      end
