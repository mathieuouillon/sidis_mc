#include <stdio.h>
#include <stdlib.h>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"

#include "main.h"
#include "root.h"

int main(int argc, char * argv[])
{
        int i;
  //      root *routout;
        routout = new root();
        printf("Hello from %s!\n", argv[0]);
        i = simulation_();
        delete routout;
        return 0;
}

extern "C" {void fillroot_()
{
        routout->FillValues();
        return;
}

void closeroot_()
{
        routout->CloseFile();
        return;
}}
