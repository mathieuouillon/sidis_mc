      subroutine QWComput(qhat,ipx,ipy,ipz)
      implicit none

      include 'common.f'

      real ipx,ipy,ipz
      real qhat 
      real radius
      real x,y,z
      real pp,px,py,pz
      real integral_step
      parameter(integral_step = 0.1)
      real d

      QW_wc = 0.
      QW_R  = 0.
      d = 0

ccc normalize momentum
      pp = sqrt(ipx**2+ipy**2+ipz**2)
      px = ipx/pp * integral_step
      py = ipy/pp * integral_step
      pz = ipz/pp * integral_step

      x = x_inter
      y = y_inter
      z = z_inter
      
      radius = sqrt(x**2+y**2+z**2)
      do while (radius.lt.20)
        d = d + integral_step
        x = x + px
        y = y + py
        z = z + pz
        radius = sqrt(x**2+y**2+z**2)
        QW_wc = QW_wc + integral_step * d *
     &          density_table(INT(radius/step_size_dens))
        QW_R = QW_R + integral_step *
     &          density_table(INT(radius/step_size_dens))
      enddo

      QW_wc = qhat/density_table(1) * QW_wc
      QW_R = 2 * density_table(1) * QW_wc**2 / QW_R / qhat

      end
