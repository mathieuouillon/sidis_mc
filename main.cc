#include <stdio.h>
#include <stdlib.h>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

#include "main.h"

int main(int argc, char * argv[])
{
        int i;
        printf("Hello from %s!\n", argv[0]);
        i = simulation_();
        return 0;
}

int InitROOT()
{
   TFile f("iron.root","recreate");
   TTree t1("t1","a simple Tree with simple variables");
   
   t1.Branch("ievent",&kinematics_.ievent,"ievent");
   t1.Branch("EEe   ",&kinematics_.EEe   ,"EEe"   );
   t1.Branch("PPn   ",&kinematics_.PPn   ,"PPn"   );

   return 0;
}

int FillROOT()
{
    
//     t1.Fill();

//    call HBNAME(33,'TransVar',PhiFM,'PhiFM')
//    call HBNAME(33,'TransVar',ThFM ,'ThFM ')
//    call HBNAME(33,'TransVar',Kf   ,'Pf  ')
//    call HBNAME(33,'TransVar',ECoM ,'ECoM ')

//    call HBNAME(33,'Position',x_inter ,'x_inter ')
//    call HBNAME(33,'Position',y_inter ,'y_inter ')
//    call HBNAME(33,'Position',z_inter ,'z_inter ')

//    call HBNAME(33,'QWeight',QW_wc ,'QW_wc')
//    call HBNAME(33,'QWeight',QW_R  ,'QW_R ')
//    call HBNAME(33,'QWeight',QW_L  ,'QW_L ')
//    call HBNAME(33,'QWeight',QW_w  ,'QW_w ')

//    call HBNAME(33,'EvntInfo',Q22  ,'Q2   ')
//    call HBNAME(33,'EvntInfo',W    ,'W    ')
//    call HBNAME(33,'EvntInfo',Nu   ,'GamNu')
//    call HBNAME(33,'EvntInfo',XBj  ,'XBj  ')
//    call HBNAME(33,'EvntInfo',y_ele,'y    ')
//  
//    call HBNAME(33,'PartInfo',Nb_part  ,'Nb_part[0,99]      ')
//    call HBNAME(33,'PartInfo',id_part  ,'Npart_id(Nb_part)  ')
//    call HBNAME(33,'PartInfo',id_mother,'Nmother_id(Nb_part)')
//    call HBNAME(33,'PartInfo',acc_part ,'Naccept(Nb_part)   ')
//    call HBNAME(33,'PartInfo',p_part   ,'part_P(Nb_part)    ')
//    call HBNAME(33,'PartInfo',px_part  ,'part_Px(Nb_part)   ')
//    call HBNAME(33,'PartInfo',py_part  ,'part_Py(Nb_part)   ')
//    call HBNAME(33,'PartInfo',pz_part  ,'part_Pz(Nb_part)   ')
//    call HBNAME(33,'PartInfo',E_part   ,'part_E(Nb_part)    ')
//    call HBNAME(33,'PartInfo',m_part   ,'part_m(Nb_part)    ')
//    call HBNAME(33,'PartInfo',z_part   ,'part_z(Nb_part)    ')
//    call HBNAME(33,'PartInfo',th_part  ,'part_th(Nb_part)   ')
//    call HBNAME(33,'PartInfo',phi_part ,'part_phi(Nb_part)  ')
//    call HBNAME(33,'PartInfo',phih_part,'part_phih(Nb_part) ')
//    call HBNAME(33,'PartInfo',tt_part  ,'part_tt(Nb_part)   ')
//    call HBNAME(33,'PartInfo',Pts_part ,'part_Pts(Nb_part)  ')


   return 0;
}
