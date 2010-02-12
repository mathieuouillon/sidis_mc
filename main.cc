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
// TFile f("tree1.root","recreate");
// TTree t1("t1","a simple Tree with simple variables");
// Float_t px, py, pz;
// Double_t random;
// Int_t ev;
// t1.Branch("px",&px,"px/F");
// t1.Branch("py",&py,"py/F");
// t1.Branch("pz",&pz,"pz/F");
// t1.Branch("random",&random,"random/D");
// t1.Branch("ev",&ev,"ev/I");
// 
// //fill the tree
// for (Int_t i=0;i<10000;i++) {
//   gRandom->Rannor(px,py);
//   pz = px*px + py*py;
//   random = gRandom->Rndm();
//   ev = i;
//   t1.Fill();
//}
   return 0;
}

int FillROOT()
{

   return 0;
}
