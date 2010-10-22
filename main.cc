#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

#include "main.h"

using namespace std;

TFile *f;
TTree *t2;
void branching();

int main(int argc, char * argv[])
{
        int i;
        if (argc > 1) f= new TFile(argv[1],"recreate");
        else f= new TFile("default.root","recreate");
        t2= new TTree("t2","Tree with simulated events");
        branching();
        i = simulation_();
        t2->Write();
        return 0;
}

extern "C" {void fillroot_() {
  //  printf("Nb part %i\n",part_.Nb_part);
  //  for( int ii = 0 ; ii<part_.Nb_part ; ii++)
  //  printf("part id %i / %i\n",ii,part_.id_part[ii]);

      t2->Fill();
      return; } }

void branching(){
        t2->Branch("ievent"    , &kinematics_.ievent   , "ievent/I"             );
        t2->Branch("Eee"       , &kinematics_.EEe      , "Eee/F"                );
        t2->Branch("Ppn"       , &kinematics_.PPn      , "Ppn/F"                );
        t2->Branch("Phifm"     , &transfovar_.PhiFM    , "Phifm/F"              );
        t2->Branch("Thfm"      , &transfovar_.ThFM     , "Thfm/F"               );
        t2->Branch("Pf"        , &transfovar_.Kf       , "Pf/F"                 );
        t2->Branch("Ecom"      , &transfovar_.ECoM     , "Ecom/F"               );
        t2->Branch("x_inter"   , &interacpos_.x_inter  , "x_inter/F"            );
        t2->Branch("y_inter"   , &interacpos_.y_inter  , "y_inter/F"            );
        t2->Branch("z_inter"   , &interacpos_.z_inter  , "z_inter/F"            );
        t2->Branch("Qw_wc"     , &quenwei_.QW_wc       , "Qw_wc/F"              );
        t2->Branch("Qw_r"      , &quenwei_.QW_R        , "Qw_r/F"               );
        t2->Branch("Qw_l"      , &quenwei_.QW_L        , "Qw_l/F"               );
        t2->Branch("Qw_w"      , &quenwei_.QW_w        , "Qw_w/F"               );
        t2->Branch("Qw_th"     , &quenwei_.QW_th       , "Qw_th/F"              );
        t2->Branch("Q2"        , &nt_.Q22              , "Q2/F"                 );
        t2->Branch("W"         , &nt_.W                , "W/F"                  );
        t2->Branch("Gamnu"     , &nt_.Nu               , "Gamnu/F"              );
        t2->Branch("Xbj"       , &nt_.XBj              , "Xbj/F"                );
        t2->Branch("y"         , &nt_.y_ele            , "y/F"                  );
        t2->Branch("Nb_part"   , &part_.Nb_part        , "Nb_part/I"            );
        t2->Branch("Npart_id"  , part_.id_part         , "Npart_id[Nb_part]/I"  );
        t2->Branch("Nmother_id", part_.id_mother       , "Nmother_id[Nb_part]/I");
        t2->Branch("Naccept"   , part_.acc_part        , "Naccept[Nb_part]/I"   );
        t2->Branch("part_p"    , part_.p_part          , "part_p[Nb_part]/F"    );
        t2->Branch("part_px"   , part_.px_part         , "part_px[Nb_part]/F"   );
        t2->Branch("part_py"   , part_.py_part         , "part_py[Nb_part]/F"   );
        t2->Branch("part_pz"   , part_.pz_part         , "part_pz[Nb_part]/F"   );
        t2->Branch("part_e"    , part_.E_part          , "part_e[Nb_part]/F"    );
        t2->Branch("part_m"    , part_.m_part          , "part_m[Nb_part]/F"    );
        t2->Branch("part_z"    , part_.z_part          , "part_z[Nb_part]/F"    );
        t2->Branch("part_th"   , part_.th_part         , "part_th[Nb_part]/F"   );
        t2->Branch("part_phi"  , part_.phi_part        , "part_phi[Nb_part]/F"  );
        t2->Branch("part_phih" , part_.phih_part       , "part_phih[Nb_part]/F" );
        t2->Branch("part_tt"   , part_.tt_part         , "part_tt[Nb_part]/F"   );
        t2->Branch("part_pts"  , part_.Pts_part        , "part_pts[Nb_part]/F"  );
        t2->Branch("part_xf"   , part_.Xf_part         , "part_xf[Nb_part]/F"   );

        return; }
