//FastMC for A Low Energy Recoil Tracker 
//
//Author: Gabriel Charles

// An example on how to use the propagation function is provided in the main(), one can find below how to understand it

// An object MyParticle has been created, it has 8 arguments that should be set before propagating the particle
// and 5 arguments that are obtained by the propagation
// The 8 arguments that must be filled before propagation are
// Pid: the particle number following pdg (see below for the particles treated)
// Mom: the initial momentum of the particle
// VertexX: the X position of the vertex
// VertexY: the Y position of the vertex
// VertexZ: the Z position of the vertex
// VertexMomX: the projection of the initial direction of the particle over X
// VertexMomY: the projection of the initial direction of the particle over Y
// VertexMomZ: the projection of the initial direction of the particle over Z
// The following conditions must be respected TMath::Sqrt(vertex_px*vertex_px+vertex_py*vertex_py+vertex_pz*vertex_pz) = 1
// The momentum of the particle in MeV/c

// The 5 arguments of MyParticle filled during the propagation are:
// pT: the transverse momentum
// Theta
// Phi
// VertexZ_recon, which is the z vertex position obtained after reconstruction
// isRecon, that tells if the particle has been reconstructed or not

// To set an arguement just do:
// MyParticle pc; // to define a new particle
// pc.Mom = 100.; // to set the value of the initial momentum
// pc.vertex_px/py/pz/... // to set all the required parameters
// MyParticle propagated = Propagate(pc); // to propagate the particle and store propagate values in a new particle
// cout << propagated.pT << endl; // would display the transverse momentum reconstructed after propagation
// // Note that initial values are NOT passed to the propagated particle.

// Here is a list of the particles that can be input in the fastMC with their pdg code
// proton 2212
// alpha: 1000020040
// deuteron: 1000010020
// triton: 1000010030
// helium 3: 1000020030

// There is no particle identification for now, so the code is returning the same pid as the one given in input, like if the particle
// identification was 100% efficient

#include <iostream>   
#include <fstream>
#include <string>
#include <iomanip>
#include <locale>
#include <sstream>
#include <math.h>
#include <cstdlib>
#include <cstring>

#define PI 3.14159265358979323846
extern "C" {

   double mass_proton = 938.; // MeV/c2
   double mass_alpha = 3727.;
   double mass_tritium = 2809.;
   double mass_deuterium = 1876.;
   double mass_He3 = 2809.;
 
   double q_proton = 1.;
   double q_alpha = 2.;
   double q_tritium = 1.;
   double q_deuterium = 1.;
   double q_He3 = 2.;

using namespace std;

const int Nparticles = 5;
const int NVar = 5;
const int MaxNbin_Th = 100;
const int MaxNbin_E = 184;
double ResoPt[Nparticles][MaxNbin_E][MaxNbin_Th];
double ResoPhi[Nparticles][MaxNbin_E][MaxNbin_Th];
double ResoTheta[Nparticles][MaxNbin_E][MaxNbin_Th];
double ResoZ[Nparticles][MaxNbin_E][MaxNbin_Th];
double Acc[Nparticles][MaxNbin_E][MaxNbin_Th];

double GausRand(double, double);
bool IsReconFunc(double);
double funcEk(double, double);
double funcMom(double, double);

// Recoiled Particles
  extern struct {
      int TrkGS;
      int Nb_part,id_part[100],id_mother[100],acc_part[100];
      float p_part[100],px_part[100],py_part[100],pz_part[100];
      float E_part[100],m_part[100],z_part[100],th_part[100];
      float tt_part[100],Pts_part[100],phih_part[100],phi_part[100];
      float vx_part[100],Xf_part[100],Als[100];
  }part_;
void initalert_(){
 
   char Title[100];

   double ener, theta;   
   double val, err_val;
   double acc;
   
   int Nbin_E;
   double Eini, Efin;
   const int Nbin_th=100;
   const double Thini=0.;
   const double Thfin = PI;
   double stepTh, stepE;
   
   int ind_th = 0;
   int ind_th_prev = -1; // used to solve problems of roundings
   
   for(int pc=0;pc<Nparticles;pc++){
      for(int var=0;var<NVar;var++){
   
         strcpy(Title,"datafiles/alert/Bare_Mix_3atm_1atm_");

         if(pc==0){
            strcat(Title,"Protons");
            Nbin_E=38.;
            Eini=2;
            Efin=40.;
         }
         else if(pc==1){
            strcat(Title,"Alphas");
            Nbin_E=42.;
            Eini=4;
            Efin=25.;
         }
         else if(pc==2){
            strcat(Title,"Deutons");
            Nbin_E=46.;
            Eini=2;
            Efin=25.;
         }
         else if(pc==3){
            strcat(Title,"H3");
            Nbin_E=46.;
            Eini=2;
            Efin=25.;
         }
         else if(pc==4){
            strcat(Title,"He3");
            Nbin_E=46.;
            Eini=2;
            Efin=25.;
         } 
	 
	 stepTh = (Thfin-Thini)/(1.*Nbin_th);
         stepE = (Efin-Eini)/(1.*Nbin_E);
	 
	 if(var==0) strcat(Title,"_Pt.txt");
	 else if(var==1) strcat(Title,"_Phi.txt");
	 else if(var==2) strcat(Title,"_Th.txt");
	 else if(var==3) strcat(Title,"_Z.txt");
	 else if(var==4) strcat(Title,"_Acc.txt");
	    
         ifstream file(Title, ios::in);  
         if(file){ 
            while(!file.eof()){  
               file >> ener >> theta >> val >> err_val; 
	       ind_th = (int) ((theta)/stepTh);
	       if(ind_th==ind_th_prev) ind_th++; // because of problems due to roundings
               if(var==0)      ResoPt[pc][(int) ((ener-Eini)/stepE)][ind_th] = val;
	       else if(var==1) ResoPhi[pc][(int) ((ener-Eini)/stepE)][ind_th] = val;
	       else if(var==2) ResoTheta[pc][(int) ((ener-Eini)/stepE)][ind_th] = val;
	       else if(var==3) ResoZ[pc][(int) ((ener-Eini)/stepE)][ind_th] = val;
	       else if(var==4) Acc[pc][(int) ((ener-Eini)/stepE)][ind_th] = val;
	       
	       ind_th_prev=ind_th;
            }
            file.close();
            cout << Title << " loaded." << endl; 
         }
         else cerr << Title << " not found." << endl;
      
      }  
   }   
}


int alertaccept_(int &a, float &b, float &c){

   int iPid = a;
   float Mom = b*1000;
   float Th = c/180*PI;

   int pc;
   double mymass, myq;
   int Nbin_E;
   int isRecon;
   double Eini, Efin;
   const int Nbin_th=100;
   const double Thini=0.;
   const double Thfin = PI;
   
//________________________________________________________________________________________________
// ______________________________________ Identify pc type _______________________________________
//________________________________________________________________________________________________


   if(iPid==2212){ // The particle is a proton      
//      cout << "It is a proton." << endl;      
      pc=0;
      mymass = mass_proton;
      myq = q_proton;
            Nbin_E=38.;
            Eini=2;
            Efin=40.;
   }
   else if(iPid==10204){ // The particle is an alpha
//      cout << "It is an alpha." << endl;      
      pc=1;
      mymass = mass_alpha;
      myq = q_alpha;
            Nbin_E=42.;
            Eini=4;
            Efin=25.;
   }
   else if(iPid==10102){  
//      cout << "It is a deuteron." << endl;      
      pc=2;
      mymass = mass_deuterium;
      myq = q_deuterium;
            Nbin_E=46.;
            Eini=2;
            Efin=25.;
   }
   else if(iPid==10103){   
//      cout << "It is a triton." << endl;      
      pc=3;
      mymass = mass_tritium;
      myq = q_tritium;
            Nbin_E=46.;
            Eini=2;
            Efin=25.;
   }
   else if(iPid==10203){   
//      cout << "It is a helium 3." << endl;      
      pc=4;
      mymass = mass_He3;
      myq = q_He3;
            Nbin_E=46.;
            Eini=2;
            Efin=25.;
   }
   else {return 0;}
   
   double stepTh = (Thfin-Thini)/(1.*Nbin_th);
   double stepE = (Efin-Eini)/(1.*Nbin_E);
   
   double pT = Mom*sin(Th);
   double th_gen=Th;

   double Ene=funcEk(mymass,double(Mom));
   if(Ene<Eini) return 0;

   if(th_gen<0) th_gen+=PI;
   
   isRecon = IsReconFunc(Acc[pc][(int) ((Ene-Eini)/stepE)][(int) (double(Th)/stepTh)]);

   return isRecon; 
}



double GausRand(double mu, double sigma){

   double x1, x2, w, y1, y2;
   while(1){
      x1 = 2.0 * (double)rand()/RAND_MAX - 1.0;
      x2 = 2.0 * (double)rand()/RAND_MAX - 1.0;

      w = x1 * x1 + x2 * x2;
      if(w<1.0) break;
   }
 
   w = sqrt( (-2.0 * log( w ) ) / w );
   y1 = x1 * w;
   y2 = x2 * w;

   return y1 * sigma + mu; 
}

bool IsReconFunc(double Prob){
   bool isRec = false;
   if(((double)rand()/RAND_MAX)<Prob/100.) isRec = true;
   return isRec;
}

double funcEk(double mass, double momentum){

 return sqrt(momentum*momentum+mass*mass)-mass;

}

double funcMom(double mass, double Ek){

 return sqrt(Ek*Ek+2.*mass*Ek);

}

}
