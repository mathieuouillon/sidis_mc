#define root_cxx
#include "root.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


void root::FillValues()
{
   printf("test\n");
   ievent = (Int_t)kinematics_.ievent; 
   h33->Fill();
   return;
}

void root::CloseFile()
{
   h33->Write();
   return;
}

