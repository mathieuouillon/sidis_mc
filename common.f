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
      real Thi,Phi !angles for rotation
      real Beta,ECoM ! For calculation of the CoM energy
      common/TransfoVar/ThFM,PhiFM,BB1,B1x,B1y,B1z,Thi,Phi,Beta,ECoM,Kf

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



ccccc Events information
      real Q22,W,Nu,XBj
      common/nt/Q22,W,Nu,XBj

ccccc Recoiled Particles
      integer TrkGS
      integer Nb_part,id_part(100),id_mother(100)
      real px_part(100),py_part(100),pz_part(100)
      real E_part(100),m_part(100),z_part(100),th_part(100)
      real tt_part(100),Pts_part(100)
      common/part/Nb_part,id_part,id_mother,
     &      px_part,py_part,pz_part,E_part,m_part,
     &      z_part,th_part,tt_part,Pts_part


