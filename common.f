ccccc Initial Kinematic Block
      integer ievent
      real EEn,PPn,Pnx,Pny,Pnz ! nucleon kinematics
      real EEe,PPe,Pex,Pey,Pez ! electron kinematics
      real ele_ene,ele_the,ele_phi,nuc_mom,nuc_the,nuc_phi !value to book
      common/kinematics/ievent,EEn,PPn,Pnx,Pny,Pnz,EEe,PPe,Pex,Pey,Pez,
     &                  ele_ene,ele_the,ele_phi,nuc_mom,nuc_the,nuc_phi

ccccc Important variables for transformations
      real ThFM,PhiFM,Kf !Theta, Phi and K Fermi
      real BB1,B1x,B1y,B1z ! Lorentz kinematics
      real BB2,B2x,B2y,B2z ! Lorentz kinematics
      real Thi,Phi !angles for rotation
      real Beta,ECoM ! For calculation of the CoM energy
      common/TransfoVar/ThFM,PhiFM,Kf,BB1,B1x,B1y,B1z,BB2,B2x,B2y,B2z,
     &                  Thi,Phi,Beta,ECoM

ccccc Files names
      character*60 bosout,hbookout
      common /OUT_NAMES/ bosout,hbookout

C...Stuff for PYTHIA 6.4
      double precision PARP(200),PARI(200)
      integer MSTP(200),MSTI(200)
      COMMON/PYPARS/MSTP,PARP,MSTI,PARI
      double precision PARU(200),PARJ(200)
      integer MSTU(200),MSTJ(200)
      COMMON/PYDAT1/MSTU,PARU,MSTJ,PARJ
      double precision P(4000,5),V(4000,5)
      integer N,NPAD,K(4000,5)
      COMMON/PYJETS/N,NPAD,K,P,V
      integer MSEL,MSELPD,MSUB(500),KFIN(2,-40:40)
      double precision CKIN(200)
      COMMON/PYSUBS/MSEL,MSELPD,MSUB,KFIN,CKIN
      integer NGENPD,NGEN(0:500,3)
      double precision XSEC(0:500,3)
      COMMON/PYINT5/NGENPD,NGEN,XSEC

ccccc Events information
      real Q22,W,Nu,XBj,y_ele
      common/nt/Q22,W,Nu,XBj,y_ele

ccccc Recoiled Particles
      integer TrkGS
      integer Nb_part,id_part(100),id_mother(100),acc_part(100)
      integer ch_part(100)
      real p_part(100),px_part(100),py_part(100),pz_part(100)
      real E_part(100),m_part(100),z_part(100),th_part(100)
      real tt_part(100),Pts_part(100),phih_part(100),phi_part(100)
      real vx_part(100),Xf_part(100),Als(100)
      common/part/TrkGS,Nb_part,id_part,id_mother,acc_part,ch_part,
     &      p_part,px_part,py_part,pz_part,E_part,m_part,
     &      z_part,th_part,tt_part,Pts_part,phih_part,phi_part,
     &      vx_part,Xf_part,Als

ccccc Table for FM distribution
      real FM_table(1000),step_size_FM,FM_n(1000),FM_p(1000)
      real FM_i(1000), FMintact
      integer iZ, iA ! target Z and A
      integer FMnb
      common/FMvar/FM_table,step_size_FM,FM_n,FM_p,FM_i,FMintact,iZ,iA,
     &             FMnb

ccccc Table for density distribution
      real density_table(2000),step_size_dens
      real quantity_table(2000),init_dens
      common/Density/density_table,step_size_dens,quantity_table,
     &               init_dens

ccccc Interaction position
      real x_inter,y_inter,z_inter
      real pos_radius,pos_theta,pos_phi
      common/InteracPos/x_inter,y_inter,z_inter,
     &      pos_radius,pos_theta,pos_phi

ccccc Quenching Weight variables
      integer QW_nb
      real QW_wc,QW_R,QW_L,QW_w,QW_qhat,QW_chi,QW_th
      common/QuenWei/QW_nb,QW_wc,QW_R,QW_L,QW_w,QW_qhat,QW_chi,QW_th

ccccc AA routine variables
      double precision alphas
            integer iqw,scor,ncor,sfthrd,irw
      common/qw/alphas,iqw,scor,ncor,sfthrd,irw

ccccc Miscellanous
c random number generator from CERNLIB
      real ranf

c Position of the interaction in CLAS for GSIM
      real vxz
      common/vertex/vxz

      real pi
      data pi/3.1415926535/

c Flags and config values
      integer iTg,iFM,iDens,iQuenching,iSim,iNS,iAccept,nevent 
      integer nkin,nucleon,specId,iColl,iqg,iEg,iPtF,iIso,iAlert
      integer iLund
      real rFM,E0,qhat,ehat,EColl,FMlimit,SupFac
      common/flags/iTg,iFM,iDens,iQuenching,iSim,iNS,iAccept,iLund,
     &             nevent,
     &             nkin,nucleon,specId,rFM,E0,qhat,ehat,iColl,EColl,
     &             FMlimit,iqg,iEg,SupFac,iPtF,iIso,iAlert

c Recoil acceptance
      real pro_acc(25,25),deu_acc(25,25),tri_acc(25,25)
      real he3_acc(25,25),he4_acc(25,25)
      common/recoil/ pro_acc,deu_acc,tri_acc,he3_acc,he4_acc

      real vz