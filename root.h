//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Feb 18 14:34:39 2010 by ROOT version 5.18/00b
// from TTree h33/out
// found on file: helium.root
//////////////////////////////////////////////////////////

#ifndef root_h
#define root_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// ****************** Variables that can be in the ROOT output
// Initial Kinematic Block
extern "C" {
  extern struct {
      int ievent;
      float EEn,PPn,Pnx,Pny,Pnz;
      float EEe,PPe,Pex,Pey,Pez;
      float ele_ene,ele_the,ele_phi,nuc_mom,nuc_the,nuc_phi;
  }kinematics_;
}
// Important variables for transformations
extern "C" {
  extern struct {
      float ThFM,PhiFM,Kf;
      float BB1,B1x,B1y,B1z;
      float Thi,Phi;
      float Beta,ECoM;
  }transfovar_;
}

// Events information
extern "C" {
  extern struct {
      float Q22,W,Nu,XBj,y_ele;
  }nt_;
}

// Recoiled Particles
extern "C" {
  extern struct {
      int TrkGS;
      int Nb_part,id_part[100],id_mother[100],acc_part[100];
      float p_part[100],px_part[100],py_part[100],pz_part[100];
      float E_part[100],m_part[100],z_part[100],th_part[100];
      float tt_part[100],Pts_part[100],phih_part[100],phi_part[100];
      float vx_part[100];
  }part_;
}

// Interaction position
extern "C" {
  extern struct {
      float x_inter,y_inter,z_inter;
      float pos_radius,pos_theta,pos_phi;
  }interacpos_;
}

// Quenching Weight variables
extern "C" {
  extern struct {
      float QW_wc,QW_R,QW_L,QW_w;
  }quenwei_;
}


class root {
public :
   TTree          *h33;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           ievent;
   Float_t         Eee;
   Float_t         Ppn;
   Float_t         Phifm;
   Float_t         Thfm;
   Float_t         Pf;
   Float_t         Ecom;
   Float_t         x_inter;
   Float_t         y_inter;
   Float_t         z_inter;
   Float_t         Qw_wc;
   Float_t         Qw_r;
   Float_t         Qw_l;
   Float_t         Qw_w;
   Float_t         Q2;
   Float_t         W;
   Float_t         Gamnu;
   Float_t         Xbj;
   Float_t         y;
   Int_t           Nb_part;
   Int_t           Npart_id[109];   //[Nb_part]
   Int_t           Nmother_id[109];   //[Nb_part]
   Int_t           Naccept[109];   //[Nb_part]
   Float_t         part_p[109];   //[Nb_part]
   Float_t         part_px[109];   //[Nb_part]
   Float_t         part_py[109];   //[Nb_part]
   Float_t         part_pz[109];   //[Nb_part]
   Float_t         part_e[109];   //[Nb_part]
   Float_t         part_m[109];   //[Nb_part]
   Float_t         part_z[109];   //[Nb_part]
   Float_t         part_th[109];   //[Nb_part]
   Float_t         part_phi[109];   //[Nb_part]
   Float_t         part_phih[109];   //[Nb_part]
   Float_t         part_tt[109];   //[Nb_part]
   Float_t         part_pts[109];   //[Nb_part]

   // List of branches
   TBranch        *b_ievent;   //!
   TBranch        *b_Eee;   //!
   TBranch        *b_Ppn;   //!
   TBranch        *b_Phifm;   //!
   TBranch        *b_Thfm;   //!
   TBranch        *b_Pf;   //!
   TBranch        *b_Ecom;   //!
   TBranch        *b_x_inter;   //!
   TBranch        *b_y_inter;   //!
   TBranch        *b_z_inter;   //!
   TBranch        *b_Qw_wc;   //!
   TBranch        *b_Qw_r;   //!
   TBranch        *b_Qw_l;   //!
   TBranch        *b_Qw_w;   //!
   TBranch        *b_Q2;   //!
   TBranch        *b_W;   //!
   TBranch        *b_Gamnu;   //!
   TBranch        *b_Xbj;   //!
   TBranch        *b_y;   //!
   TBranch        *b_Nb_part;   //!
   TBranch        *b_Npart_id;   //!
   TBranch        *b_Nmother_id;   //!
   TBranch        *b_Naccept;   //!
   TBranch        *b_part_p;   //!
   TBranch        *b_part_px;   //!
   TBranch        *b_part_py;   //!
   TBranch        *b_part_pz;   //!
   TBranch        *b_part_e;   //!
   TBranch        *b_part_m;   //!
   TBranch        *b_part_z;   //!
   TBranch        *b_part_th;   //!
   TBranch        *b_part_phi;   //!
   TBranch        *b_part_phih;   //!
   TBranch        *b_part_tt;   //!
   TBranch        *b_part_pts;   //!

   root();
   virtual ~root();
   void FillValues();
   void CloseFile();
};

#endif

#ifdef root_cxx
root::root()
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   TFile f("helium.root","recreate");
   TTree *h33 = new TTree("h33","ntuple");
  
   h33->Branch("ievent"    , &ievent               , "ievent/I");
// h33->Branch("Eee"       , &kinematics_.EEe      , "Eee                /F");
// h33->Branch("Ppn"       , &kinematics_.PPn      , "Ppn                /F");
// h33->Branch("Phifm"     , &transfovar_.PhiFM    , "Phifm              /F");
// h33->Branch("Thfm"      , &transfovar_.ThFM     , "Thfm               /F");
// h33->Branch("Pf"        , &transfovar_.Kf       , "Pf                 /F");
// h33->Branch("Ecom"      , &transfovar_.ECoM     , "Ecom               /F");
// h33->Branch("x_inter"   , &interacpos_.x_inter  , "x_inter            /F");
// h33->Branch("y_inter"   , &interacpos_.y_inter  , "y_inter            /F");
// h33->Branch("z_inter"   , &interacpos_.z_inter  , "z_inter            /F");
// h33->Branch("Qw_wc"     , &quenwei_.QW_wc       , "Qw_wc              /F");
// h33->Branch("Qw_r"      , &quenwei_.QW_R        , "Qw_r               /F");
// h33->Branch("Qw_l"      , &quenwei_.QW_L        , "Qw_l               /F");
// h33->Branch("Qw_w"      , &quenwei_.QW_w        , "Qw_w               /F");
// h33->Branch("Q2"        , &nt_.Q22              , "Q2                 /F");
// h33->Branch("W"         , &nt_.W                , "W                  /F");
// h33->Branch("Gamnu"     , &nt_.Nu               , "Gamnu              /F");
// h33->Branch("Xbj"       , &nt_.XBj              , "Xbj                /F");
// h33->Branch("y"         , &nt_.y_ele            , "y                  /F");
// h33->Branch("Nb_part"   , &Nb_part  , "Nb_part            /I");
// h33->Branch("Npart_id"  , Npart_id  , "Npart_id[Nb_part]  /I");
// h33->Branch("Nmother_id", Nmother_id, "Nmother_id[Nb_part]/I");
// h33->Branch("Naccept"   , Naccept   , "Naccept[Nb_part]   /I");
// h33->Branch("part_p"    , part_p    , "part_p[Nb_part]    /F");
// h33->Branch("part_px"   , part_px   , "part_px[Nb_part]   /F");
// h33->Branch("part_py"   , part_py   , "part_py[Nb_part]   /F");
// h33->Branch("part_pz"   , part_pz   , "part_pz[Nb_part]   /F");
// h33->Branch("part_e"    , part_e    , "part_e[Nb_part]    /F");
// h33->Branch("part_m"    , part_m    , "part_m[Nb_part]    /F");
// h33->Branch("part_z"    , part_z    , "part_z[Nb_part]    /F");
// h33->Branch("part_th"   , part_th   , "part_th[Nb_part]   /F");
// h33->Branch("part_phi"  , part_phi  , "part_phi[Nb_part]  /F");
// h33->Branch("part_phih" , part_phih , "part_phih[Nb_part] /F");
// h33->Branch("part_tt"   , part_tt   , "part_tt[Nb_part]   /F");
// h33->Branch("part_pts"  , part_pts  , "part_pts[Nb_part]  /F");

}

root::~root()
{
   if (!h33) return;
   delete h33->GetCurrentFile();
}

#endif // #ifdef root_cxx
