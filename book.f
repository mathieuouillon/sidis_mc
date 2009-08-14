c------------------------------------------------------------------------------
c Variables to book in the ntuple
c------------------------------------------------------------------------------
      subroutine InitHbook
      implicit none

ccccc Paw variables
      integer nwpawc
      parameter (nwpawc=900000)
      integer LREC,ISTAT
ccccc Include all the common blocks
      include 'common.f'


      LREC=  8190

      call hlimit(nwpawc)

      call hropen(10,'out',hbookout,'N',LREC,ISTAT)
      call HBNT(33,'out',' ')
    
      call HBNAME(33,'kinemati',ievent,'ievent')
      call HBNAME(33,'kinemati',EEe,'EEe')
      call HBNAME(33,'kinemati',PPn,'PPn')

      call HBNAME(33,'TransVar',PhiFM,'PhiFM')
      call HBNAME(33,'TransVar',ThFM ,'ThFM ')
      call HBNAME(33,'TransVar',Kf   ,'Pf  ')
      call HBNAME(33,'TransVar',ECoM ,'ECoM ')

      call HBNAME(33,'Position',x_inter ,'x_inter ')
      call HBNAME(33,'Position',y_inter ,'y_inter ')
      call HBNAME(33,'Position',z_inter ,'z_inter ')

      if (iQW .ne. 0) then
        call HBNAME(33,'QWeight',QW_wc ,'QW_wc')
        call HBNAME(33,'QWeight',QW_R  ,'QW_R ')
      endif

      call HBNAME(33,'EvntInfo',Q22  ,'Q2   ')
      call HBNAME(33,'EvntInfo',W    ,'W    ')
      call HBNAME(33,'EvntInfo',Nu   ,'GamNu')
      call HBNAME(33,'EvntInfo',XBj  ,'XBj  ')
      call HBNAME(33,'EvntInfo',y_ele,'y    ')
    
      call HBNAME(33,'PartInfo',Nb_part  ,'Nb_part[0,30]      ')
      call HBNAME(33,'PartInfo',id_part  ,'Npart_id(Nb_part)  ')
      call HBNAME(33,'PartInfo',id_mother,'Nmother_id(Nb_part)')
      call HBNAME(33,'PartInfo',p_part   ,'part_P(Nb_part)    ')
      call HBNAME(33,'PartInfo',px_part  ,'part_Px(Nb_part)   ')
      call HBNAME(33,'PartInfo',py_part  ,'part_Py(Nb_part)   ')
      call HBNAME(33,'PartInfo',pz_part  ,'part_Pz(Nb_part)   ')
      call HBNAME(33,'PartInfo',E_part   ,'part_E(Nb_part)    ')
      call HBNAME(33,'PartInfo',m_part   ,'part_m(Nb_part)    ')
      call HBNAME(33,'PartInfo',z_part   ,'part_z(Nb_part)    ')
      call HBNAME(33,'PartInfo',th_part  ,'part_th(Nb_part)   ')
      call HBNAME(33,'PartInfo',phi_part ,'part_phi(Nb_part)  ')
      call HBNAME(33,'PartInfo',phih_part,'part_phih(Nb_part) ')
      call HBNAME(33,'PartInfo',tt_part  ,'part_tt(Nb_part)   ')
      call HBNAME(33,'PartInfo',Pts_part ,'part_Pts(Nb_part)  ')

      end
c------------------------------------------------------------------------------
c Initialization of variables to book in the ntuple
c------------------------------------------------------------------------------
      subroutine InitKin2Book
      implicit none

ccccc Include all the common blocks
      include 'common.f'

      Nu = 0.
      Q22 = 0
      XBj = 0
      W = 0
      TrkGS = 0
      Nb_part = 0

      end
c------------------------------------------------------------------------------
c Computation of variables to book in the ntuple
c------------------------------------------------------------------------------
      subroutine ComputV
      implicit none

ccccc Include all the common blocks
      include 'common.f'

      integer ip
ccccc for calculation of Phih
      real A1,A2,A3,AA
      real B1,B2,B3,BB
      real phi_ele

C.. Booking the hbook
      do ip=1,N


        if (k(ip,2).eq.22 .and. k(ip,3).eq.1) then
          Nu = p(ip,4)
c          write(*,*) 'Nu =', Nu, p(ip,4)
          Q22 = P(ip,5)**2
          if (Nu .ne. 0) XBj = Q22 /2 /0.938 /Nu
          if (Nu .ne. 0) y_ele = Nu / p(1,4)
          W = 0.
          W = dsqrt((p(2,4)+p(ip,4))**2-(p(2,3)+p(ip,3))**2
     &                                 -p(ip,2)**2-p(ip,1)**2)
c          W = sqrt(W2)
          TrkGS = ip
          phi_ele = atan2(p(ip,2),p(ip,1))*57.2958 +210
        endif

        if (k(ip,1).eq.1 .or. k(ip,2).eq.111 .or .k(ip,2).eq.310
     &       .or. k(ip,2).eq.221 .or .k(ip,2).eq.333 .or.
     &        k(ip,2).eq.2114 .or .k(ip,2).eq.3122 ) then
          Nb_part           = Nb_part + 1
          id_part(Nb_part)  = k(ip,2)
          id_mother(Nb_part)= k(k(ip,3),2)
          px_part(Nb_part)  = p(ip,1)
          py_part(Nb_part)  = p(ip,2)
          pz_part(Nb_part)  = p(ip,3)
          p_part(Nb_part)   = sqrt(p(ip,4)**2-p(ip,5)**2)
          E_part(Nb_part)   = p(ip,4)
          m_part(Nb_part)   = p(ip,5)
          z_part(Nb_part)   = p(ip,4)/Nu
          th_part(Nb_part)  = 57.2957795*acos(pz_part(Nb_part)/
     &                sqrt(p(ip,1)**2+p(ip,2)**2+pz_part(Nb_part)**2))
          phi_part(Nb_part) = atan2(p(ip,2),p(ip,1))*57.2958 +30
          if(phi_part(Nb_part).lt.0) 
     &    phi_part(Nb_part) = phi_part(Nb_part) + 360

          A1 = sin(phi_ele)
          A2 = -cos(phi_ele)
          A3 = 0
          AA = A1**2 + A2**2 + A3**2
       
          B1 = p(TrkGS,2)*p(ip,3) - p(TrkGS,3)*p(ip,2)
          B2 = p(TrkGS,3)*p(ip,1) - p(TrkGS,1)*p(ip,3)
          B3 = p(TrkGS,1)*p(ip,2) - p(TrkGS,2)*p(ip,1)
          BB = B1**2 + B2**2 + B3**2
     
          phih_part(Nb_part) = 
     &         acos((A1*B1+A2*B2+A3*B3)/sqrt(AA*BB))*57.2958

          tt_part(Nb_part)  =
     &               (Nu-p(ip,4))**2 - (p(TrkGS,1)-p(ip,1))**2 - 
     &               (p(TrkGS,2)-p(ip,2))**2 - (p(TrkGS,3)-p(ip,3))**2
          Pts_part(Nb_part) = (p(ip,1)**2+p(ip,2)**2+p(ip,3)**2)
     &    -((p(TrkGS,1)*p(ip,1)+p(TrkGS,2)*p(ip,2)+p(TrkGS,3)*p(ip,3))
     &    / sqrt(p(TrkGS,1)**2+p(TrkGS,2)**2+p(TrkGS,3)**2))**2

        endif
      enddo      

      end


