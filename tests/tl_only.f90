!c------------------------------------------------------------------------------
!c Standalone copy of TL() for unit testing (from transfo.f90)
!c------------------------------------------------------------------------------
      subroutine TL(EEj, PPj, Pjx, Pjy, Pjz, BB, Bx, By, Bz)
          implicit none

          real EEi, PPi, Pix, Piy, Piz
          real EEj, PPj, Pjx, Pjy, Pjz
          real BB, Bx, By, Bz, GG

          GG = 1/sqrt(1 - BB**2)

          EEi = GG*(EEj - Bx*Pjx - By*Pjy - Bz*Pjz)
          Pix = -Bx*GG*EEj + Pjx + Pjx*(GG - 1)*Bx**2/BB**2 + Pjy*(GG - 1)*By*Bx/BB**2 + Pjz*(GG - 1)*Bz*Bx/BB**2
          Piy = -By*GG*EEj + Pjx*(GG - 1)*Bx*By/BB**2 + Pjy + Pjy*(GG - 1)*By*By/BB**2 + Pjz*(GG - 1)*Bz*By/BB**2
          Piz = -Bz*GG*EEj + Pjx*(GG - 1)*Bx*Bz/BB**2 + Pjz + Pjy*(GG - 1)*By*Bz/BB**2 + Pjz*(GG - 1)*Bz*Bz/BB**2
          PPi = sqrt(Pix**2 + Piy**2 + Piz**2)

          EEj = EEi
          PPj = PPi
          Pjx = Pix
          Pjy = Piy
          Pjz = Piz

      end
