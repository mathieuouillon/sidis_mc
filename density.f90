!ccccc Generate the interaction position
      subroutine InterPos
          use kinematics_module
          use file_names_module
          use event_info_module
          use particles_module
          use fermi_motion_module
          use density_module
          use interaction_module
          use quenching_module
          use config_module
          use acceptance_module
          use misc_module
          use pythia_commons
          implicit none
          ! include 'common.f90'

          real r
          integer i

!ccccc Randomize the radius and angles
          pos_radius = ranf(0)*quantity_table(2000)
          pos_theta = acos(2*ranf(0) - 1)
          pos_phi = 2*pi*ranf(0)

          r = init_dens
          i = 1.

          do while (pos_radius .ge. quantity_table(i))
              r = r + step_size_dens
              i = i + 1
          end do
          pos_radius = r + (ranf(0) - 0.5)*step_size_dens

!cccc Calculate position on cartesian axis
          x_inter = pos_radius*sin(pos_theta)*cos(pos_phi)
          y_inter = pos_radius*sin(pos_theta)*sin(pos_phi)
          z_inter = pos_radius*cos(pos_theta)

      end

!ccccc Generate the table of density in function of r
      subroutine GenNucDens
          use kinematics_module
          use file_names_module
          use event_info_module
          use particles_module
          use fermi_motion_module
          use density_module
          use interaction_module
          use quenching_module
          use config_module
          use acceptance_module
          use misc_module
          use pythia_commons
          implicit none

          ! include 'common.f90'

          integer i, idist, irho
          double precision nucdens, r
          real integral

!ccccc Some init
          QW_nb = 0
          QW_qhat = 0.

          integral = 0.
          idist = iDens
          irho = 1
          step_size_dens = 0.01
          init_dens = 0.005
          r = init_dens

          do i = 1, 2000
              density_table(i) = nucdens(r, iZ, iA, idist, irho)
!ccccc calculate the quantity of mater in function of the radius
              integral = integral + 4*PI*density_table(i)*r**2*step_size_dens
              quantity_table(i) = integral
              r = r + step_size_dens
          end do

      end
