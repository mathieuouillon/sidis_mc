      function RTPC_accept(id,HeE,theta)

      implicit none

c Alpha energy
      real HeE
c Input id
      integer id
      real theta
c output accept
      integer RTPC_accept

c Parameters
      real TarGasL,TarGasD
      real TarWalL,TarWalD
      real HelBagL,HelBagD
      real AlMFoiL,AlMFoiD
      real NeoDmeL,NeoDMED
      real DriftGL,DriftGD
      real CurrDen,CurrPos
c 1 He / 2 Kapton / 3 Neon
      integer CurrMat

      real CurrSP,eloss
      real scoef

c For loops
      integer dist
      integer i

      include 'common.f'

      RTPC_accept = 0

c  Values are in cm and g
      TarGasL = .3
      TarGasD = 0.000166*(75./14.7+1.)
      TarWalL = 0.0030
      TarWalD = 1.42
      HelBagL = 1.7
      HelBagD = 0.000166
      AlMFoiL = 0.0006
      AlMFoiD = 1.42
      NeoDmeL = 1.
      NeoDmeD = 0.0009
      DriftGL = 3.
      DriftGD = 0.0009

      if (id.eq.10204) then
       scoef = 1.
      else if (id.eq.10203) then
       scoef = 9./16.
      else if (id.eq.10103) then
       scoef = 9./64.
      else 
       scoef = 0
       write(*,*) 'ID ERROR'
       goto 225
      endif

      if ((1-.3/abs(tan(theta))).lt.0) goto 225

      do dist=0,6000
        CurrPos = real(dist)/1000
        i = 1
        do while (HeE.gt.AlphaHe(1,i)) 
          i = i+1
        enddo
        if(abs(HeE-AlphaHe(1,i)).gt.abs(HeE-AlphaHe(1,i-1))) i = i-1
c        write(*,*) HeE,AlphaHe(1,i)

        if(CurrPos.lt.TarGasL) then
          CurrDen = TarGasD
          CurrMat = 1
        else if(CurrPos.lt.TarGasL+TarWalL) then
          CurrDen = TarWalD
          CurrMat = 2
        else if(CurrPos.lt.TarGasL+TarWalL+HelBagL) then
          CurrDen = HelBagD
          CurrMat = 1
        else if(CurrPos.lt.TarGasL+TarWalL+HelBagL+AlMFoiL) then
          CurrDen = AlMFoiD
          CurrMat = 2
        else if(CurrPos.lt.TarGasL+TarWalL+HelBagL
     &                                    +AlMFoiL+NeoDmeL) then
          CurrDen = NeoDmeD
          CurrMat = 3
        else if(CurrPos.lt.TarGasL+TarWalL+HelBagL
     &                                    +2*AlMFoiL+NeoDmeL) then
          CurrDen = AlMFoiD
          CurrMat = 2
        else if(CurrPos.lt.TarGasL+TarWalL
     &                      +HelBagL+2*AlMFoiL+NeoDmeL+driftGL) then
          CurrDen = DriftGD
          CurrMat = 3
        else 
          CurrDen = 0
        endif
        if (CurrMat .eq.1) then
          CurrSP = AlphaHe(2,i)
        else if (CurrMat .eq.2) then
          CurrSP = AlphaKa(2,i)
        else if (CurrMat .eq.3) then
          CurrSP = AlphaNe(2,i)
        else
          CurrSP = 0
        endif

        eloss =  CurrSP * CurrDen * 0.001 *scoef * 
     &             sqrt(1+cos(theta)**2)
        HeE = HeE - eloss
c        if (mod(dist,1000).eq.0)
c     &     write(*,*) CurrPos,CurrDen,CurrSP,eloss,HeE

      enddo

      if (HeE.gt.0) RTPC_accept = 1

 225  continue

      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine readtables()
      implicit none
      include 'common.f'

c Dummy
      character a1
      character a2
      character a3

      integer i

c Read tables
      open(8,file='AlphaInHe',status='old')
      do i=1,8 
         read(8,*)
      enddo

      do i=1,121
      read(8,"(ES9.3E3,A1,ES9.3E3,A1,ES9.3E3,A1,F6.4)") 
     &      AlphaHe(1,i),a1,AlphaHe(2,i),a2,AlphaHe(3,i),a3,AlphaHe(4,i)
      enddo

      close(8)

      open(8,file='AlphaInKapton',status='old')
      do i=1,8 
         read(8,*)
      enddo

      do i=1,121
      read(8,"(ES9.3E3,A1,ES9.3E3,A1,ES9.3E3,A1,F6.4)") 
     &      AlphaKa(1,i),a1,AlphaKa(2,i),a2,AlphaKa(3,i),a3,AlphaKa(4,i)
      enddo

      close(8)

      open(8,file='AlphaInNeon',status='old')
      do i=1,8 
         read(8,*)
      enddo

      do i=1,121
      read(8,"(ES9.3E3,A1,ES9.3E3,A1,ES9.3E3,A1,F6.4)") 
     &      AlphaNe(1,i),a1,AlphaNe(2,i),a2,AlphaNe(3,i),a3,AlphaNe(4,i)
      enddo

      close(8)

      end
