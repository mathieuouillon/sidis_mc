#ifndef MAIN_H
#define MAIN_H

extern "C" {

  int simulation_();

// Initial Kinematic Block
  extern struct {
      int ievent;
      float EEn,PPn,Pnx,Pny,Pnz;
      float EEe,PPe,Pex,Pey,Pez;
      float ele_ene,ele_the,ele_phi,nuc_mom,nuc_the,nuc_phi;
  }kinematics_;

// Important variables for transformations
  extern struct {
      float ThFM,PhiFM,Kf;
      float BB1,B1x,B1y,B1z;
      float BB2,B2x,B2y,B2z;
      float Thi,Phi;
      float Beta,ECoM;
  }transfovar_;

// Events information
  extern struct {
      float Q22,W,Nu,XBj,y_ele;
  }nt_;

// Recoiled Particles
  extern struct {
      int TrkGS;
      int Nb_part,id_part[100],id_mother[100],acc_part[100];
      float p_part[100],px_part[100],py_part[100],pz_part[100];
      float E_part[100],m_part[100],z_part[100],th_part[100];
      float tt_part[100],Pts_part[100],phih_part[100],phi_part[100];
      float vx_part[100],Xf_part[100];
  }part_;

// Interaction position
  extern struct {
      float x_inter,y_inter,z_inter;
      float pos_radius,pos_theta,pos_phi;
  }interacpos_;

// Quenching Weight variables
  extern struct {
      int QW_nb;
      float QW_wc,QW_R,QW_L,QW_w,QW_qhat;
  }quenwei_;
}

#endif
