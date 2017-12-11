/* ----------- Head file for C, C++ and Root ---------------- */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stddef.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <signal.h>
#include <TROOT.h>
#include <TFile.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include "JGenPhotonEnergy.h"   //JLA  including new jlab routine ...
#include "JGenPhaseSpace.h"     //  ... replacing TGenPhaseSpace
#include "JGenFermiMomentum.h"
using namespace std;
/* ----------- Head file for C, C++ and Root ---------------- */


/* ----------- Function declaration ---------------- */
vector<double> readdcs(string filename); //read unpolarized differential cross section for the first step
double calcdcs_pin(vector<double> vect, double w, double theta, int &status);
double calcdcs_kl(vector<double> vect, double w, double theta,int &status);
double calcdcs_ks(vector<double> vect, double w, double theta, int &status);
double calcdcs_piplusn(vector<double> vect, double w, double theta, int &status); // calculate unpolarized differential cross section for the first step

void readcxczpy_kl(string filename, double e_left[16][10], double e_right[16][10], double e_aver[16][10], double ctk_left[16][10], double ctk_right[16][10], double ctk_aver[16][10], double cx[16][10], double cz[16][10], double py[16][10]); // read circular polarization observables Cx, Cz and P for KN rescattering and QF of gamma + p -> kaon + Lambda
vector<double> calccxczpy_kl(double e, double ctk, int &obs_status, double e_left[16][10], double e_right[16][10], double e_aver[16][10], double ctk_left[16][10], double ctk_right[16][10], double ctk_aver[16][10], double cx[16][10], double cz[16][10], double py[16][10]); // calcualte Cx, Cz and P for KN rescattering and QF of gamma + p -> kaon + Lambda

vector<vector<double> > readSTOxOzPforGPKL(string filename); // read linear polarization observables Sigma, T, Ox, Oz and P for KN rescattering and QF of gamma + p -> kaon + Lambda
vector<double> calcSTOxOzPforGPKL(double w, double ct, int &obs_status, vector<vector<double> > vectObs); // calcualte Sigma, T, Ox, Oz and P for KN rescattering and QF of gamma + p -> kaon + Lambda

vector<double> reads(string filename); //read Sigma for the first step of two pion mediated mechanisms
double calcs_pin(vector<double> vect, double e, double theta,int &obs_status);
double calcs_piplusn(vector<double> vect, double e, double theta,int &obs_status);//calculate Sigma for the first step of two pion mediated mechanisms

double calcdcsKL_original(double Wv, double ctv, int &Sstatus); // calculate unpolarized differential cross section for the second step using oringial table
vector<double> readSdcs(string filename);
double calcdcsKL_adjust(vector<double> vect, double w, double ct, int &status);

void PrintUsage (char *processName);
/* ----------- Function declaration ---------------- */

/* ----------- main function ---------------- */
int main (int argc, char **argv){
  //// Define a object for the class of JGenPhotonEnergy & Deal with arguments & set lower and upper of random number for selected channels////  
  extern char *optarg; //save option argument
  char *ProgName = argv[0]; // save command
  int c; // save option order
  char *ROOTFILE = NULL; // save output rootfile name
  int WillBeRootOutput = 0; // condition for root file output
  int max = 1; // save number of events
  int Nchannel=9; // save number of channels with implementation
  char *channel="123456789"; // save channel IDs
  string cs1="ncs"; // save argument for differential cross section of the first step
  string cs2="ncs"; // save argument for differential cross section of the second step
  string rm="diff"; // save argument for the choice of maximum of random value for the first step. 0: Different channels have different maximum values dependent of their own different cross section table; 1: All channels have the same maximum value, which is maximum value of all tables for the case of unpolarized differential cross section and is maximum value of all tables multiplies a factor for the case of polairzed differential cross section.

  JGenPhotonEnergy *beamE=new JGenPhotonEnergy(0.95,2.6); //JLA  default plain distribution for beam energy

  //Options and Arguments
  if (argc == 1){
    PrintUsage (ProgName);
    exit (0);
  }

  while ((c = getopt (argc, argv, "hE:R:M:n:c:p:q:s:")) != -1){
    switch (c){
    case 'h':
      PrintUsage (ProgName); // help
      break;
    case 'R':
      WillBeRootOutput = 1; // condition of root file output
      ROOTFILE = optarg;
      break;
    case 'M':
      max = atoi (optarg); // number of events
      break;
    case 'E':          
      beamE = new JGenPhotonEnergy(optarg);   //JLA added option -E to generate energy distribution
      break;
    case 'n':          
      Nchannel = atoi(optarg);   //How many processed channels
      break;
    case 'c':          
      channel = optarg;   //Channel number --- For example: If you want to process channel 2, 4, 5, and 6, "-c 2456" should be conclude in the executive command. 
      break;
    case 'p':          
      cs1 = optarg;   //differential cross section for the first step
      break;
    case 'q':          
      cs2 = optarg;   //differential cross section for the second step
      break;
    case 's':          
      rm = optarg;   //differential cross section for the second step
      break;
    default:
      fprintf (stderr, "Unrecognized argument: [-%c]\n\n", c);
      PrintUsage (ProgName);
      exit (1);
      break;
    }
  }

  double channel_random_left[9]={0}, channel_random_right[9]={0}; 
  for(int i=0;i<Nchannel;i++){
    char temp[10]={channel[i]};
    channel_random_left[atoi(temp)-1]=i*1./Nchannel;
    channel_random_right[atoi(temp)-1]=(i+1)*1./Nchannel; // set lower and upper of random number for selected channels
  }

  //// Define a object for the class of JGenPhotonEnergy & Deal with arguments & set lower and upper of random number for selected channels//// 

  //// Variables defination ////
  double weight,Sweight;
  int label = 0; //indicates different reactions
  double rndnumber; //a ramdon number in order to decide which reaction to produce
  double cmpol; //polar angle in center of mass for the first step
  double tinW_kl=0.0, tinW_ks=0.0, tinW_pin=0.0, tinW_piplusn=0.0,tinW_ssk=0.0,tinW_ssminusk=0.0; // W for the first step

  double ran1=0, ran2=0, ran3=0, ran4l=0, ran4c=0, ran5=0; //random number: ran1 and ran2 are for extraction of Fermi momentum, ran3 is for comparison to unpolarized differential cross section of the first step, ran4l is for comparison to linearly polarized differential cross section of the first step, ran4c is for comparison to circularly polarized differential cross section of the first step, and ran5 is for comparsion to unpolarized differential cross section of the second step
  double fcos=0, fphi=0,fpx, fpy, fpz, fEp, fEn; //variables for Fermi Momentum
  double betax, betay, betaz; //variables for lorantz boost for the first step
  TLorentzVector temp1, temp2; //temp1: template of Lorentz vector of product 1 or product 2 of the first step; temp2: template of Lorentz vector of beam

  double Sbetax, Sbetay, Sbetaz; //variables for lorantz boost for the second step
  TLorentzVector Stemp1, Stemp2; //Stemp1: template of Lorentz vector of product 2 of the second step; Stemp2: template of Lorentz vector of beam

  double masses_lk[2],masses_pin[2],Smasses_ln[2],Smasses_kn[2],Smasses_lk[2],masses_sk[2],Smasses_sn[2],masses_piplusn[2]; // masses of two products of the first and second steps
  masses_lk[0] = 1.115683;
  masses_lk[1] = 0.493677;
  masses_pin[0] = 0.1349764;
  masses_pin[1] = 0.939565378;
  masses_sk[0] = 1.192642;
  masses_sk[1] = 0.493677;
  masses_piplusn[0] = 0.13975018;
  masses_piplusn[1] = 0.939565378;
  Smasses_ln[0] = 1.115683;
  Smasses_ln[1] = 0.939565378;
  Smasses_kn[0] = 0.493677;
  Smasses_kn[1] = 0.939565378;
  Smasses_lk[0] = 1.115683;
  Smasses_lk[1] = 0.493677;
  Smasses_sn[0] = 1.192642;
  Smasses_sn[1] =  0.939565378;

  double masses_ssk[2],masses_ssminusk[2],Smasses_lpi[2],Smasses_lpiminus[2]; // masses of two products of the first and second steps
  masses_ssk[0]=1.3837;
  masses_ssk[1]=0.493677;
  masses_ssminusk[0]=1.3872;
  masses_ssminusk[1]=0.493677;
  Smasses_lpi[0]=1.115683;
  Smasses_lpi[1]=0.1349764;
  Smasses_lpiminus[0]=1.115683;
  Smasses_lpiminus[1]=0.13975018;

  float X=-1000., Y=-1000., Z=-1000.; //Position

  int csstat=0; // Identify W is inside or outside the range of cross section table for the first step: 0 - outside; 1 - inside
  double csval=0.0,loc_wval=0.0,loc_theta=0.0; // W, theta and calculated cross section for the first step
  //// Variables defination ////

  //// Variable defination for polaried differential cross section of the first step ////
  //Polarized cross section are adopted by channels 1, 2 and 5. Channel 2 is for circular polarized dcs calculated by Cx, Cz, Py and unpolarized dcs, and channels 1 & 5 are for linear polarized differential dcs calculated by Sigma and unpolarided dcs
  double Dmasses_ppi[2];
  Dmasses_ppi[0]=0.938272046;
  Dmasses_ppi[1]=0.13957018;
  JGenPhaseSpace event_decayppi;

  TLorentzVector *pDp1 = NULL;
  TLorentzVector Dp1;
  pDp1 = &Dp1;	
  TLorentzVector *pDp2 = NULL;
  TLorentzVector Dp2;
  pDp2 = &Dp2;  // About availabes are only for channel 2 KN re-scattering which need to know information of decay proton

  int obsstat=0; // Identify beam energy and polar angle is inside or outside the range of observable tables
  double lcsval=0.0; // calculated linearly polarized cross section
  double ccsval=0.0; // calculated circularly polarized cross section
  double cx=0,cz=0,py=0; // calculated Cx, Cz and Py
  double costhetax,costhetay,costhetaz; //cosine theta angle of proton in the RF of Lambda
  int cirpol=0; // Heilicty of circular polarized beam
  double phi; //Used for observables Sigma
  double sigma=0,t=0,ox=0,oz=0; // observable Sigma, T, Ox, Oz;
  int linpol=0; // para or perp of linear polarized beam

  TLorentzVector Pphot,Pkaon,Pprot,Plamb; // Only for channel 2
  double beta2x,beta2y,beta2z; //Boost in RF of Lambda for channel 2 
  TVector3 Vphot,Vkaon,Vprot,Vx,Vy,Vz; // Only for channel 2
  vector<double> vect_cxczpy; // save Cx, Cz and Py calculated by the function "calccxczpy_KL"
  vector<double> vect_stoxozpy; // save Sigma, T, Ox, Oz and P calculated by the function "calcSTOxOzPforGPKL"

  int Nlabel2=0; //Used for setting cirpol for the channel 2
  int Nlabel8=0; //Used for setting cirpol for the channel 8
  int Nlabel1=0; //Used for setting linpol for the channel 1
  int Nlabel5=0; //Used for setting linpol for the channel 5
  //// Variable defination for polaried cross section of the first step ////

  //// Variable defination for the second step ////
  double Scmpol; //polar angle in center of mass for the second step
  double StinW_KL=0.0; // W for the second step
  int Scsstat=0; // Identify W is inside or outside the range of cross section table for the second step: 0 - outside; 1 - inside
  double Scsval_original=0.0, Scsval_adjust=0.0, loc_Swval=0.0,loc_Stheta=0.0; // W, theta and calculated cross section for the second step
  //// Variable defination for the second step ////

  //// TLorentzVector for beam target and final particles ////
  TLorentzVector *pbeam = NULL;
  TLorentzVector beam;
  pbeam = &beam;
  TLorentzVector *ptarget = NULL;
  TLorentzVector target;
  ptarget = &target;
  TLorentzVector *pW = NULL;
  TLorentzVector W;
  pW = &W;
  TLorentzVector *pQp1 = NULL;
  TLorentzVector Qp1;
  pQp1 = &Qp1;
  TLorentzVector *pQp2 = NULL;
  TLorentzVector Qp2;
  pQp2 = &Qp2;	
  TLorentzVector *pSbeam = NULL;
  TLorentzVector Sbeam;
  pSbeam = &Sbeam;
  TLorentzVector *pStarget = NULL;
  TLorentzVector Starget;
  pStarget = &Starget;	
  TLorentzVector *pSW =NULL;
  TLorentzVector SW;
  pSW = &SW;	
  TLorentzVector *pRp1 = NULL;
  TLorentzVector Rp1;
  pRp1 = &Rp1;	
  TLorentzVector *pRp2 = NULL;
  TLorentzVector Rp2;
  pRp2 = &Rp2; 
  //// TLorentzVector for beam target and final particles ////

  ////Some preparation////
  TRandom3 randomNum(0);  // randomNum is used in the main function
  TRandom3 *myRNG=new TRandom3(0);
  gRandom = myRNG;  //gRandom is used in the JGenPhotonEnergy and JGenPhaseSpace classed

  JGenPhaseSpace event,event_finalln,event_finalkn,event_pin,event_finalpip,event_sk,event_finalsn,event_piplusn,event_finalpiplusn; // JGenPhaseSpace
  JGenPhaseSpace event_ssk,event_ssminusk, event_finallpi,event_finallpiminus; //Phase space calculation of the first and second steps

  int Nevents = 0; //Serial number of events
  TLorentzVector deuteron;
  deuteron.SetXYZT(0,0,0,1.87561);

  double StinW_KL_recalc=0;
  TLorentzVector Sbeam_recalc, Starget_recalc; // Unpolarized different cross section for the second step of channels 1 and 5 is calcuated through another channel with application of Clebsch–Gordan coefficients for isospin multiplets of targets and beams. Masses of targets and beams may be different for different channel, so center-of-mass energy W should be recalculated.
  ////Some preparation////

  //// Root Stuff ////
  TFile RootOut (ROOTFILE, "recreate", "");	// This is ROOT output file
  TTree *mytree = new TTree ("mytree", "A TTree object");
  mytree->Branch ("beam", "TLorentzVector", &pbeam, 32000, 1); // Lorentz vector of beam of the first step
  mytree->Branch ("target", "TLorentzVector", &ptarget, 32000, 1); // Lorentz vector of target of the first step
  mytree->Branch ("W", "TLorentzVector", &pW, 32000, 1); // Lorentz vecotr of beam+target of the first step
  mytree->Branch ("Qp1", "TLorentzVector", &pQp1, 32000, 1); // Lorentz vector of product 1 of the first step
  mytree->Branch ("Qp2", "TLorentzVector", &pQp2, 32000, 1); // Lorentz vector of product 2 of the first step
  mytree->Branch ("Sbeam", "TLorentzVector", &pSbeam, 32000, 1); // Lorentz vector of beam of the second step
  mytree->Branch ("Starget", "TLorentzVector", &pStarget, 32000, 1); // Lorentz vecotr of target of the second step
  mytree->Branch ("SW", "TLorentzVector", &pSW, 32000, 1); // Lorentz vector of beam+target of the second step
  mytree->Branch ("Rp1", "TLorentzVector", &pRp1, 32000, 1); // Lorentz vector of product 1 of the second step
  mytree->Branch ("Rp2", "TLorentzVector", &pRp2, 32000, 1); // Lorentz vector of product 2 of the second step
  mytree->Branch ("Dp1", "TLorentzVector", &pDp1, 32000, 1); // Lorentz vector of product 1 of Lambda decay; NULL except channel 2
  mytree->Branch ("Dp2", "TLorentzVector", &pDp2, 32000, 1); // Lorentz vector of prodcut 2 of Lambda decay; NULL except channel 2
  mytree->Branch ("weight", &weight, "weight/D");
  mytree->Branch ("Sweight", &Sweight, "Sweight/D");
  mytree->Branch ("Nevents", &Nevents, "Nevents/I"); // Event ID
  mytree->Branch ("label", &label, "label/I"); // Channel ID
  mytree->Branch ("rnumber", &rndnumber, "rndnumber/D"); // Random number for random channel
  mytree->Branch ("X", &X, "X/F"); // Position X; default value; useless currently 
  mytree->Branch ("Y", &Y, "Y/F"); // Position Y; default value; useless currently
  mytree->Branch ("Z", &Z, "Z/F"); // Position Z; default value; useless currently
  mytree->Branch ("csstat", &csstat, "csstat/I"); // Identify W is outside or inside of differential cross section table of the first step: 0 - outside; 1 - inside
  mytree->Branch ("csval", &csval, "csval/D"); // Calculated value of unpolairzed differential cross section of the first step
  mytree->Branch ("wval", &loc_wval, "wval/D"); // W of the first step
  mytree->Branch ("theta", &loc_theta, "theta/D"); // theta of the first step in CMS

  mytree->Branch ("Scsstat", &Scsstat, "Scsstat/I"); // Identify W is outside or inside of differential cross section table of the second step: 0 - outside; 1 - inside
  mytree->Branch ("Scsval_original", &Scsval_original, "Scsval_original/D"); // Calculated value of unpolairzed differential cross section of the second step using original table
  mytree->Branch ("Scsval_adjust", &Scsval_adjust, "Scsval_adjust/D"); // Calculated value of unpolairzed differential cross section of the second step using adjusted table
  mytree->Branch ("Swval", &loc_Swval, "Swval/D"); // W of the second step
  mytree->Branch ("Stheta", &loc_Stheta, "Stheta/D"); // theta of the second step in CMS


  mytree->Branch ("obsstat", &obsstat, "obsstat/I"); // Identify E_gamma or W is outside or inside of polarization observable tables: 0 - outside; 1 - inside
  mytree->Branch ("lcsval", &lcsval, "lcsval/D"); // Calculated value of linearly polarized differential cross section of the first step
  mytree->Branch ("ccsval", &ccsval, "ccsval/D"); // Calculated value of circularly polarized differential cross section of the first step
  mytree->Branch ("costhetax", &costhetax, "costhetax/D"); // cosines of x direction of decay proton in Lambda RF
  mytree->Branch ("costhetay", &costhetay, "costhetay/D"); // cosines of y direction of decay proton in Lambda RF
  mytree->Branch ("costhetaz", &costhetaz, "costhetaz/D"); // cosines of z direction of decay proton in Lambda RF
  mytree->Branch ("cx", &cx, "cx/D"); // Calculated Cx for calculation polarized different cross section of the first step
  mytree->Branch ("cz", &cz, "cz/D"); // Calculated Cz for calculation polarized different cross section of the first step
  mytree->Branch ("py", &py, "py/D"); // Calculated Py for calculation polarized different cross section of the first step
  mytree->Branch ("cirpol", &cirpol, "cirpol/I"); // Helicity of polarized helicity beam
  mytree->Branch("phi",&phi,"phi/D"); // phi 
  mytree->Branch("sigma",&sigma,"sigma/D");  
  mytree->Branch("t",&t,"t/D");  
  mytree->Branch("ox",&ox,"ox/D");  
  mytree->Branch("oz",&oz,"oz/D");  
  mytree->Branch ("linpol", &linpol, "linpol/I");
  //// Root Stuff ////

  //// Read differential cross section for different table lists for the first step using the function "readdcs" ////
  vector<double> vect_pin=readdcs("dcs_pin.txt");
  vector<double> vect_kl=readdcs("dcs_kl.txt");
  vector<double> vect_ks=readdcs("dcs_ks.txt");
  vector<double> vect_piplusn=readdcs("dcs_piplusn.txt");
  //// Read differential cross section for different table lists for the first step using the function "readdcs" ////

  //// Read Sigma for different table lists for the first step using the function "reads" ////
  vector<double> vect_spin=reads("Spin.txt");
  vector<double> vect_spiplusn=reads("Spiplusn.txt");
  //// Read Sigma for different table lists for the first step using the function "reads" ////

  //// Read circular polarization observables for KN rescattering and QF of gamma + p -> kaon + Lambda
  double e_left[16][10], e_right[16][10], e_aver[16][10], ctk_left[16][10], ctk_right[16][10], ctk_aver[16][10], cxTable[16][10], czTable[16][10], pyTable[16][10];
  readcxczpy_kl("obs_kl.txt", e_left, e_right, e_aver, ctk_left, ctk_right, ctk_aver, cxTable, czTable, pyTable);
  //// Read circular polarization observables for KN rescattering and QF of gamma + p -> kaon + Lambda

  //// Read linear polarization observables for KN rescattering and QF of gamma + p -> kaon + Lambda
  vector<vector<double> > vectSTOxOzPforGPKL;
  vectSTOxOzPforGPKL=readSTOxOzPforGPKL("STOxOzPforGPKL");
  //// Read linear polarization observables for KN rescattering and QF of gamma + p -> kaon + Lambda



  //// Read differential cross section for table list for the second step using the function "readSdcs" ////
  vector<double> vect_SKL=readSdcs("csKL_adjust");
  //// Read differential cross section for table list for the second step using the function "readSdcs" ////

  //// Loop begin ////
  cout << "Generating the events..." << endl;
  for (Nevents = 0; Nevents < max;){

    //// Polarization Setup ////
    //For channels are not adopt circularly polarized cross section or linearly polarized cross section, 100 is set for default values.
    obsstat=100; lcsval=100; ccsval=100; cx=100; cz=100; py=100; costhetax=100; costhetay=100; costhetaz=100; cirpol=100; phi=100; sigma=100; t=100; ox=100; oz=100; linpol=100;
    //Just channel 2 includes decay. For channels wihtout decay, pDp1 and pDp2 are NULL
    pDp1 = NULL; pDp2 = NULL;
    //For quasi-free channel 8 and 9, only the first step is involved
    pRp1 = NULL; pRp2 = NULL; pSbeam=NULL; pSW=NULL; Sweight=0;
    //// Polarization Setup ////

    //// Cross section setup for the second step ////
    Scsstat=-100; Scsval_original=-100; Scsval_adjust=-100; loc_Swval=-100; loc_Stheta=-100;

    //// Cross section setup for the second step ////

    //// Print out a message every 100 events ////
    if (Nevents % 100 == 0){
      fprintf (stderr, "%d\r", Nevents);
      fflush (stderr);
    }
    //// Print out a message every 100 events ////

    //// Use reject-accept method to select a value of Fermi Momentum ////
    ran1 = 0.5*randomNum.Rndm();
    ran2 = 11.0*randomNum.Rndm();
    if (JGenFermiMomentum::Instance().Spectral(ran1)<ran2) continue;
    fcos = (2*randomNum.Rndm())-1;
    fphi = 2*3.141592653*randomNum.Rndm();
    fpx = ran1*sqrt(1-fcos*fcos)*cos(fphi);
    fpy = ran1*sqrt(1-fcos*fcos)*sin(fphi);
    fpz = ran1*fcos;
    fEp = sqrt(fpx*fpx+fpy*fpy+fpz*fpz+0.93827231*0.93827231);
    fEn = sqrt(fpx*fpx+fpy*fpy+fpz*fpz+0.939565378*0.939565378);
    //// Use reject-accept method to select a value of Fermi Momentum ////

    ////Use JGenPhotonEnergy class to set beam ////
    beamE -> Generate (); // ... to generate an event
    beam.SetXYZT (0.0, 0.0, beamE->GetE(),beamE->GetE());
    ////Use JGenPhotonEnergy class to set beam ////

    rndnumber=randomNum.Rndm (); // Random reaction 

    //-----------generate events for pi0 mediated reaction (label=1)---------//
    if (rndnumber>=channel_random_left[0] && rndnumber<=channel_random_right[0]) {
      Starget.SetXYZT (-fpx, -fpy, -fpz, fEp);
      target=deuteron-Starget;
      W = beam + target;
      tinW_pin = 1000*W.Mag();
      label = 1;

      if (event_pin.SetDecay (W, 2, masses_pin)){
	weight = event_pin.Generate ();
	pQp1 = event_pin.GetDecay (0);
	pQp2 = event_pin.GetDecay (1);

	// Differential cross section : Calculate w and theta firstly, and then obtain differential cross section //
	if(rm == "diff" || rm == "same"){
	  ran3 = 30.9*randomNum.Rndm(); // For unpolarized differential cross section, 30.9 is the maximum value of unpolarized differential cross section table for the first step of pi0 mediated channel, and also is maximum for channels 1, 2, 3, 4, 5, 8 and 9 
	  ran4l = 24*randomNum.Rndm(); // For linearly polarized differential cross section, 24 is the maximum value of linear polarized differential cross section for the first step of pi0 mediated channel, and also is maximum for channel 1,2,5 and 8
	}

	temp1 = *pQp1;
	temp2 = beam;
	betax = -(beam.Px()+target.Px())/(beam.Energy()+target.Energy());
	betay = -(beam.Py()+target.Py())/(beam.Energy()+target.Energy());
	betaz = -(beam.Pz()+target.Pz())/(beam.Energy()+target.Energy());
	temp1.Boost(betax, betay, betaz);
	temp2.Boost(betax, betay, betaz);
	cmpol = (180/3.141592653)*acos((temp1.Px()*temp2.Px()+temp1.Py()*temp2.Py()+temp1.Pz()*temp2.Pz())/(temp1.P()*temp2.P()));
	   
	csval=calcdcs_pin(vect_pin,tinW_pin,cmpol, csstat);

	loc_wval=tinW_pin;
	loc_theta=cmpol;
	// Differential cross section : Calculate w and theta firstly, and then obtain differential cross section 
	if(cs1=="ncs"){
	  Sbeam=*pQp1;
	  SW = Sbeam + Starget;
	  StinW_KL=1000*SW.Mag();

	  Sbeam_recalc.SetXYZM(Sbeam.Px(),Sbeam.Py(),Sbeam.Pz(),0.13975018);
	  StinW_KL_recalc=1000*(Sbeam_recalc+Starget).Mag();

	  if (event_finalpip.SetDecay (SW, 2, Smasses_lk)){		
	    Sweight = event_finalpip.Generate ();
	    pRp1 = event_finalpip.GetDecay (0);
	    pRp2 = event_finalpip.GetDecay (1);

	    ran5 = 113.3*randomNum.Rndm(); //113.3 is the maximum value of unpolarized differential cross section based on table of the second step
	    Stemp1 = *pRp2;
	    Stemp2 = Sbeam;
	    Sbetax = -(Sbeam.Px()+Starget.Px())/(Sbeam.Energy()+Starget.Energy());
	    Sbetay = -(Sbeam.Py()+Starget.Py())/(Sbeam.Energy()+Starget.Energy());
	    Sbetaz = -(Sbeam.Pz()+Starget.Pz())/(Sbeam.Energy()+Starget.Energy());
	    Stemp1.Boost(Sbetax, Sbetay, Sbetaz);
	    Stemp2.Boost(Sbetax, Sbetay, Sbetaz);
	    Scmpol = (180/3.141592653)*acos((Stemp1.Px()*Stemp2.Px()+Stemp1.Py()*Stemp2.Py()+Stemp1.Pz()*Stemp2.Pz())/(Stemp1.P()*Stemp2.P()));

	    Scsval_original=calcdcsKL_original(StinW_KL_recalc, cos(Scmpol/180*3.141592653), Scsstat)/2;
	    Scsval_adjust=calcdcsKL_adjust(vect_SKL, StinW_KL_recalc, cos(Scmpol/180*3.141592653), Scsstat)/2;

	    loc_Swval=StinW_KL;
	    loc_Stheta=Scmpol;

	    if(cs2=="ncs"){
	      if (event_decayppi.SetDecay (*pRp1, 2, Dmasses_ppi)){		 
		event_decayppi.Generate ();
		pDp1 = event_decayppi.GetDecay (0);
		pDp2 = event_decayppi.GetDecay (1);	      
		Nevents++;
	    
		// Root Stuff //
		if (WillBeRootOutput){
		  mytree->Fill ();
		}	
		// Root Stuff //
	      }
	    }

	    else if(cs2=="ocs"){
	      if(Scsval_original<ran5) continue;
	      if (event_decayppi.SetDecay (*pRp1, 2, Dmasses_ppi)){		 
		event_decayppi.Generate ();
		pDp1 = event_decayppi.GetDecay (0);
		pDp2 = event_decayppi.GetDecay (1);
		Nevents++;
	    
		// Root Stuff //
		if (WillBeRootOutput){
		  mytree->Fill ();
		}	
		// Root Stuff //
	      }
	    }

	    else if(cs2=="acs"){
	      if(Scsval_adjust<ran5) continue;
	      if (event_decayppi.SetDecay (*pRp1, 2, Dmasses_ppi)){		 
		event_decayppi.Generate ();
		pDp1 = event_decayppi.GetDecay (0);
		pDp2 = event_decayppi.GetDecay (1);
		Nevents++;
	    
		// Root Stuff //
		if (WillBeRootOutput){
		  mytree->Fill ();
		}	
		// Root Stuff //
	      }
	    }
	  }
	}

	if(cs1=="ucs" || cs1=="ccs"){
	  if(csval<ran3) continue;	 
	  Sbeam=*pQp1;
	  SW = Sbeam + Starget;
	  StinW_KL=1000*SW.Mag();

	  Sbeam_recalc.SetXYZM(Sbeam.Px(),Sbeam.Py(),Sbeam.Pz(),0.13975018);
	  StinW_KL_recalc=1000*(Sbeam_recalc+Starget).Mag();

	  if (event_finalpip.SetDecay (SW, 2, Smasses_lk)){		
	    Sweight = event_finalpip.Generate ();
	    pRp1 = event_finalpip.GetDecay (0);
	    pRp2 = event_finalpip.GetDecay (1);

	    ran5 = 113.3*randomNum.Rndm(); //113.3 is the maximum value of unpolarized differential cross section based on table of the second step
	    Stemp1 = *pRp2;
	    Stemp2 = Sbeam;
	    Sbetax = -(Sbeam.Px()+Starget.Px())/(Sbeam.Energy()+Starget.Energy());
	    Sbetay = -(Sbeam.Py()+Starget.Py())/(Sbeam.Energy()+Starget.Energy());
	    Sbetaz = -(Sbeam.Pz()+Starget.Pz())/(Sbeam.Energy()+Starget.Energy());
	    Stemp1.Boost(Sbetax, Sbetay, Sbetaz);
	    Stemp2.Boost(Sbetax, Sbetay, Sbetaz);
	    Scmpol = (180/3.141592653)*acos((Stemp1.Px()*Stemp2.Px()+Stemp1.Py()*Stemp2.Py()+Stemp1.Pz()*Stemp2.Pz())/(Stemp1.P()*Stemp2.P()));

	    Scsval_original=calcdcsKL_original(StinW_KL_recalc, cos(Scmpol/180*3.141592653), Scsstat)/2;
	    Scsval_adjust=calcdcsKL_adjust(vect_SKL, StinW_KL_recalc, cos(Scmpol/180*3.141592653), Scsstat)/2;

	    loc_Swval=StinW_KL;
	    loc_Stheta=Scmpol;


	    if(cs2=="ncs"){
	      if (event_decayppi.SetDecay (*pRp1, 2, Dmasses_ppi)){		 
		event_decayppi.Generate ();
		pDp1 = event_decayppi.GetDecay (0);
		pDp2 = event_decayppi.GetDecay (1);
		Nevents++;
	    
		// Root Stuff //
		if (WillBeRootOutput){
		  mytree->Fill ();
		}	
		// Root Stuff //
	      }
	    }

	    else if(cs2=="ocs"){
	      if(Scsval_original<ran5) continue;
	      if (event_decayppi.SetDecay (*pRp1, 2, Dmasses_ppi)){		 
		event_decayppi.Generate ();
		pDp1 = event_decayppi.GetDecay (0);
		pDp2 = event_decayppi.GetDecay (1);
		Nevents++;
	    
		// Root Stuff //
		if (WillBeRootOutput){
		  mytree->Fill ();
		}	
		// Root Stuff //
	      }
	    }

	    else if(cs2=="acs"){
	      if(Scsval_adjust<ran5) continue;
	      if (event_decayppi.SetDecay (*pRp1, 2, Dmasses_ppi)){		 
		event_decayppi.Generate ();
		pDp1 = event_decayppi.GetDecay (0);
		pDp2 = event_decayppi.GetDecay (1);
		Nevents++;
	    
		// Root Stuff //
		if (WillBeRootOutput){
		  mytree->Fill ();
		}	
		// Root Stuff //
	      }
	    }
	  }
	}

	if(cs1=="lcs"){
	  Sbeam=*pQp1;
	  SW = Sbeam + Starget;
	  StinW_KL=1000*SW.Mag();

	  Sbeam_recalc.SetXYZM(Sbeam.Px(),Sbeam.Py(),Sbeam.Pz(),0.13975018);
	  StinW_KL_recalc=1000*(Sbeam_recalc+Starget).Mag();

	  if (event_finalpip.SetDecay (SW, 2, Smasses_lk)){		
	    Sweight = event_finalpip.Generate ();
	    pRp1 = event_finalpip.GetDecay (0);
	    pRp2 = event_finalpip.GetDecay (1);

	    if(Nlabel1%2==0) linpol=1;
	    else linpol=0;

	    phi=pQp1->Phi();
	    sigma=calcs_pin(vect_spin, beam.E(), cmpol, obsstat); // phi is defined as azimuthal angle of pi in the lab frame
	    if(linpol==1) lcsval=csval*(1-sigma*cos(2*phi));
	    else lcsval=csval*(1+sigma*cos(2*phi));


	    if(lcsval<ran4l) continue;

	    ran5 = 113.3*randomNum.Rndm(); //113.3 is the maximum value of unpolarized differential cross section based on table of the second step
	    Stemp1 = *pRp2;
	    Stemp2 = Sbeam;
	    Sbetax = -(Sbeam.Px()+Starget.Px())/(Sbeam.Energy()+Starget.Energy());
	    Sbetay = -(Sbeam.Py()+Starget.Py())/(Sbeam.Energy()+Starget.Energy());
	    Sbetaz = -(Sbeam.Pz()+Starget.Pz())/(Sbeam.Energy()+Starget.Energy());
	    Stemp1.Boost(Sbetax, Sbetay, Sbetaz);
	    Stemp2.Boost(Sbetax, Sbetay, Sbetaz);
	    Scmpol = (180/3.141592653)*acos((Stemp1.Px()*Stemp2.Px()+Stemp1.Py()*Stemp2.Py()+Stemp1.Pz()*Stemp2.Pz())/(Stemp1.P()*Stemp2.P()));

	    Scsval_original=calcdcsKL_original(StinW_KL_recalc, cos(Scmpol/180*3.141592653), Scsstat)/2;
	    Scsval_adjust=calcdcsKL_adjust(vect_SKL, StinW_KL_recalc, cos(Scmpol/180*3.141592653), Scsstat)/2;

	    loc_Swval=StinW_KL;
	    loc_Stheta=Scmpol;


	    if(cs2=="ncs"){
	      if (event_decayppi.SetDecay (*pRp1, 2, Dmasses_ppi)){		 
		event_decayppi.Generate ();
		pDp1 = event_decayppi.GetDecay (0);
		pDp2 = event_decayppi.GetDecay (1);
		Nlabel1++;
		Nevents++;
	    
		// Root Stuff //
		if (WillBeRootOutput){
		  mytree->Fill ();
		}	
		// Root Stuff //
	      }
	    }

	    else if(cs2=="ocs"){
	      if(Scsval_original<ran5) continue;
	      if (event_decayppi.SetDecay (*pRp1, 2, Dmasses_ppi)){		 
		event_decayppi.Generate ();
		pDp1 = event_decayppi.GetDecay (0);
		pDp2 = event_decayppi.GetDecay (1);
		Nlabel1++;
		Nevents++;
	    
		// Root Stuff //
		if (WillBeRootOutput){
		  mytree->Fill ();
		}	
		// Root Stuff //
	      }
	    }

	    else if(cs2=="acs"){
	      if(Scsval_adjust<ran5) continue;
	      if (event_decayppi.SetDecay (*pRp1, 2, Dmasses_ppi)){		 
		event_decayppi.Generate ();
		pDp1 = event_decayppi.GetDecay (0);
		pDp2 = event_decayppi.GetDecay (1);
		Nlabel1++;
		Nevents++;
	    
		// Root Stuff //
		if (WillBeRootOutput){
		  mytree->Fill ();
		}	
		// Root Stuff //
	      }
	    }
	  }
	}
      } // Condition for the first step      
    } // Condition for channel

    //-----------generate events for pi0 mediated reaction (label=1)---------//


    //------------------------generate event for g+d-->KLn Kn rescattering (label 2)---------------------//  
    if (rndnumber>=channel_random_left[1] && rndnumber<=channel_random_right[1]) {
      Starget.SetXYZT (-fpx, -fpy, -fpz, fEn);
      target=deuteron-Starget;
      W = beam + target;
      tinW_kl = 1000*W.Mag();
      label = 2;

      if (event.SetDecay (W, 2, masses_lk)){
	weight = event.Generate ();
	pQp1 = event.GetDecay (0);
	pQp2 = event.GetDecay (1);
	    
	// Differential cross section : Calculate w and theta firstly, and then obtain differential cross section //
	if(rm == "diff"){
	  ran3 = 0.355*randomNum.Rndm(); // For unpolarized differetial cross section, 0.355 is the maximum value of unpolarized cross section table for the first step of qf, Kn and Ln rescattering
	  ran4l = 0.80*randomNum.Rndm(); // For linearly polarized differetial cross section, 0.80 is the maximum value of linearly polarized differential cross section for the first step of Kn rescattering
	}
	else if(rm == "same"){
	  ran3 = 30.9*randomNum.Rndm(); // For unpolarized differetial cross section, 30.9 is the maximum value of unpolarized cross section table for the first step of channels 1, 2, 3, 4, 5, 8 and 9
	  ran4l = 24*randomNum.Rndm(); // For linearly polarized differetial cross section, 24 is the maximum value of linearly polarized differential cross section for channels 1, 2, 5 and 8
	}

	ran4c = 0.58*randomNum.Rndm(); // For circularly polarized differetial cross section, 0.58 is the maximum value of circularly polarized diffrential cross section for the first step of Kn rescattering

	temp1 = *pQp2;
	temp2 = beam;
	betax = -(beam.Px()+target.Px())/(beam.Energy()+target.Energy());
	betay = -(beam.Py()+target.Py())/(beam.Energy()+target.Energy());
	betaz = -(beam.Pz()+target.Pz())/(beam.Energy()+target.Energy());
	temp1.Boost(betax, betay, betaz);
	temp2.Boost(betax, betay, betaz);
	cmpol = (180/3.141592653)*acos((temp1.Px()*temp2.Px()+temp1.Py()*temp2.Py()+temp1.Pz()*temp2.Pz())/(temp1.P()*temp2.P()));

	csval=calcdcs_kl(vect_kl,tinW_kl,cmpol, csstat);
	loc_wval=tinW_kl;
	loc_theta=cmpol;
	// Differential cross section : Calculate w and theta firstly, and then obtain differential cross section

	if(cs1=="ncs"){	
	  Sbeam=*pQp2;
	  SW = Sbeam + Starget;
	  if (event_finalkn.SetDecay (SW, 2, Smasses_kn)){
	    Sweight = event_finalkn.Generate ();
	    pRp1 = event_finalkn.GetDecay (0);
	    pRp2 = event_finalkn.GetDecay (1);

	    if (event_decayppi.SetDecay (*pQp1, 2, Dmasses_ppi)){		 
	      event_decayppi.Generate ();
	      pDp1 = event_decayppi.GetDecay (0);
	      pDp2 = event_decayppi.GetDecay (1);

	      Nevents++;
 
	      // Root Stuff //
	      if (WillBeRootOutput){
		mytree->Fill ();
	      }
	      // Root Stuff //
	    }
	  }
	}   

	if(cs1=="ucs"){
	  if(csval<ran3) continue;	
	  Sbeam=*pQp2;
	  SW = Sbeam + Starget;
	  if (event_finalkn.SetDecay (SW, 2, Smasses_kn)){
	    Sweight = event_finalkn.Generate ();
	    pRp1 = event_finalkn.GetDecay (0);
	    pRp2 = event_finalkn.GetDecay (1);

	    if (event_decayppi.SetDecay (*pQp1, 2, Dmasses_ppi)){		 
	      event_decayppi.Generate ();
	      pDp1 = event_decayppi.GetDecay (0);
	      pDp2 = event_decayppi.GetDecay (1);

	      Nevents++;
 
	      // Root Stuff //
	      if (WillBeRootOutput){
		mytree->Fill ();
	      }
	      // Root Stuff //
	    }
	  }
	}

	if(cs1=="ccs"){
	  Sbeam=*pQp2;
	  SW = Sbeam + Starget;
	  if (event_finalkn.SetDecay (SW, 2, Smasses_kn)){
	    Sweight = event_finalkn.Generate ();
	    pRp1 = event_finalkn.GetDecay (0);
	    pRp2 = event_finalkn.GetDecay (1);

	    if (event_decayppi.SetDecay (*pQp1, 2, Dmasses_ppi)){		 
	      event_decayppi.Generate ();
	      pDp1 = event_decayppi.GetDecay (0);
	      pDp2 = event_decayppi.GetDecay (1);

	      Pphot=beam;
	      Pphot.Boost(betax, betay, betaz);
	      Pkaon=temp1;
	      Vphot.SetXYZ(Pphot.Px(),Pphot.Py(),Pphot.Pz());
	      Vz=Vphot.Unit();
	      Vkaon.SetXYZ(Pkaon.Px(),Pkaon.Py(),Pkaon.Pz());
	      Vy=Vz.Cross(Vkaon);
	      Vy=Vy.Unit();
	      Vx=Vy.Cross(Vz);

	      Plamb=*pQp1;
	      beta2x=-Plamb.Px()/Plamb.E();    beta2y=-Plamb.Py()/Plamb.E();    beta2z=-Plamb.Pz()/Plamb.E();
	      Pprot=*pDp1;
	      Pprot.Boost(beta2x,beta2y,beta2z);
	      Vprot.SetXYZ(Pprot.Px(),Pprot.Py(),Pprot.Pz());

	      costhetax=Vprot*Vx/Vprot.Mag();
	      costhetay=Vprot*Vy/Vprot.Mag();
	      costhetaz=Vprot*Vz/Vprot.Mag();

	      vect_cxczpy.clear();
	      vect_cxczpy=calccxczpy_kl(beam.E(),cos(cmpol/180*3.141592653),obsstat,e_left, e_right, e_aver, ctk_left, ctk_right, ctk_aver, cxTable, czTable, pyTable);
	      cx=vect_cxczpy[0];
	      cz=vect_cxczpy[1];
	      py=vect_cxczpy[2];

	      if(Nlabel2%2==0) cirpol=1;
	      else cirpol=0;

	      if(cirpol==1) ccsval=csval*(1+0.642*costhetax*cx+0.642*costhetaz*cz+0.642*costhetay*py);
	      else ccsval=csval*(1-0.642*costhetax*cx-0.642*costhetaz*cz+0.642*costhetay*py);

	      if(ccsval<ran4c) continue;
	      Nlabel2++;

	      Nevents++; 
	      // Root Stuff //
	      if (WillBeRootOutput){
		mytree->Fill ();
	      }
	      // Root Stuff //
	    }	    
	  }
	}

	if(cs1=="lcs"){
	  Sbeam=*pQp2;
	  SW = Sbeam + Starget;
	  if (event_finalkn.SetDecay (SW, 2, Smasses_kn)){
	    Sweight = event_finalkn.Generate ();
	    pRp1 = event_finalkn.GetDecay (0);
	    pRp2 = event_finalkn.GetDecay (1);

	    if (event_decayppi.SetDecay (*pQp1, 2, Dmasses_ppi)){		 
	      event_decayppi.Generate ();
	      pDp1 = event_decayppi.GetDecay (0);
	      pDp2 = event_decayppi.GetDecay (1);

	      Pphot=beam;
	      Pphot.Boost(betax, betay, betaz);
	      Pkaon=temp1;
	      Vphot.SetXYZ(Pphot.Px(),Pphot.Py(),Pphot.Pz());
	      Vz=Vphot.Unit();
	      Vkaon.SetXYZ(Pkaon.Px(),Pkaon.Py(),Pkaon.Pz());
	      Vy=Vz.Cross(Vkaon);
	      Vy=Vy.Unit();
	      Vx=Vy.Cross(Vz);

	      Plamb=*pQp1;
	      beta2x=-Plamb.Px()/Plamb.E();    beta2y=-Plamb.Py()/Plamb.E();    beta2z=-Plamb.Pz()/Plamb.E();
	      Pprot=*pDp1;
	      Pprot.Boost(beta2x,beta2y,beta2z);
	      Vprot.SetXYZ(Pprot.Px(),Pprot.Py(),Pprot.Pz());

	      costhetax=Vprot*Vx/Vprot.Mag();
	      costhetay=Vprot*Vy/Vprot.Mag();
	      costhetaz=Vprot*Vz/Vprot.Mag();

	      vect_stoxozpy.clear();
	      vect_stoxozpy=calcSTOxOzPforGPKL(tinW_kl,cmpol,obsstat,vectSTOxOzPforGPKL);

	      sigma=vect_stoxozpy[0];
	      t=vect_stoxozpy[1];
	      ox=vect_stoxozpy[2];
	      oz=vect_stoxozpy[3];
	      py=vect_stoxozpy[4];

	      if(Nlabel2%2==0) linpol=1;
	      else linpol=0;

	      phi=pQp1->Phi();
	      if(linpol==1) lcsval=csval*(1-sigma*cos(2*phi)-0.642*costhetay*t*cos(2*phi)+0.642*costhetax*ox*sin(2*phi)+0.642*costhetaz*oz*sin(2*phi)+0.642*costhetay*py);
	      else lcsval=csval*(1+sigma*cos(2*phi)+0.642*costhetay*t*cos(2*phi)-0.642*costhetax*ox*sin(2*phi)-0.642*costhetaz*oz*sin(2*phi)+0.642*costhetay*py);

	      if(lcsval<ran4l) continue;
	      Nlabel2++;

	      Nevents++; 
	      // Root Stuff //
	      if (WillBeRootOutput){
		mytree->Fill ();
	      }
	      // Root Stuff //
	    }	    
	  }
	}
	
      } // Condition for the first step      
    } // Condition for channel     
    //------------------------generate event for g+d-->KLn Kn rescattering (label 2)---------------------//  
      
    //----------------------generate event g+d-->KLn Ln rescattering (label 3)----------------------//
    if (rndnumber>=channel_random_left[2] && rndnumber<=channel_random_right[2]) {
      Starget.SetXYZT (-fpx, -fpy, -fpz, fEn);
      target=deuteron-Starget;
      W = beam + target;
      tinW_kl = 1000*W.Mag();
      label = 3;

      if (event.SetDecay (W, 2, masses_lk)){
	weight = event.Generate ();
	pQp1 = event.GetDecay (0);
	pQp2 = event.GetDecay (1);
	// Differential cross section : Calculate w and theta firstly, and then obtain differential cross section //
	if(rm == "diff"){
	  ran3 = 0.355*randomNum.Rndm(); // For unpolarized differetial cross section, 0.355 is the maximum value of unpolarized cross section table for the first step of qf, Kn and Ln rescattering
	}
	else if(rm == "same"){
	  ran3 = 30.9*randomNum.Rndm(); // For unpolarized differetial cross section, 30.9 is the maximum value of unpolarized cross section table for the first step of all channels
	}

	temp1 = *pQp2;
	temp2 = beam;
	betax = -(beam.Px()+target.Px())/(beam.Energy()+target.Energy());
	betay = -(beam.Py()+target.Py())/(beam.Energy()+target.Energy());
	betaz = -(beam.Pz()+target.Pz())/(beam.Energy()+target.Energy());
	temp1.Boost(betax, betay, betaz);
	temp2.Boost(betax, betay, betaz);
	cmpol = (180/3.141592653)*acos((temp1.Px()*temp2.Px()+temp1.Py()*temp2.Py()+temp1.Pz()*temp2.Pz())/(temp1.P()*temp2.P()));

	csval=calcdcs_kl(vect_kl,tinW_kl,cmpol, csstat);
	loc_wval=tinW_kl;
	loc_theta=cmpol;          
	// Differential cross section : Calculate w and theta firstly, and then obtain differential cross section 

	if(cs1=="ncs"){
	  Sbeam=*pQp1;
	  SW = Sbeam + Starget;
	  if (event_finalln.SetDecay (SW, 2, Smasses_ln)){		 
	    Sweight = event_finalln.Generate ();
	    pRp1 = event_finalln.GetDecay (0);
	    pRp2 = event_finalln.GetDecay (1);
	    if (event_decayppi.SetDecay (*pRp1, 2, Dmasses_ppi)){		 
	      event_decayppi.Generate ();
	      pDp1 = event_decayppi.GetDecay (0);
	      pDp2 = event_decayppi.GetDecay (1);
	      Nevents++;
	    
	      // Root Stuff //
	      if (WillBeRootOutput){
		mytree->Fill ();
	      }
	      // Root Stuff //
	    } 
	  }
	}

	if(cs1=="ucs" || cs1=="ccs" || cs1=="lcs"){
	  if(csval<ran3) continue;
	  Sbeam=*pQp1;
	  SW = Sbeam + Starget;
	  if (event_finalln.SetDecay (SW, 2, Smasses_ln)){		 
	    Sweight = event_finalln.Generate ();
	    pRp1 = event_finalln.GetDecay (0);
	    pRp2 = event_finalln.GetDecay (1);
	    if (event_decayppi.SetDecay (*pRp1, 2, Dmasses_ppi)){		 
	      event_decayppi.Generate ();
	      pDp1 = event_decayppi.GetDecay (0);
	      pDp2 = event_decayppi.GetDecay (1);
	      Nevents++;
	    
	      // Root Stuff //
	      if (WillBeRootOutput){
		mytree->Fill ();
	      }
	      // Root Stuff //
	    } 
	  }
	}
      } // Condition for the first step      
    } // Condition for channel
    //----------------------generate event g+d-->KLn Ln rescattering (label 3)----------------------//

    //--------------generate event for g+d-->KSn Sn rescattering(label 4)----------------------//
    if (rndnumber>=channel_random_left[3] && rndnumber<=channel_random_right[3]) {
      Starget.SetXYZT (-fpx, -fpy, -fpz, fEn);
      target=deuteron-Starget;
      W = beam + target;
      tinW_ks = 1000*W.Mag();
      label = 4;

      if (event_sk.SetDecay (W, 2, masses_sk)){
	weight = event_sk.Generate ();
	pQp1 = event_sk.GetDecay (0);
	pQp2 = event_sk.GetDecay (1);

	// Differential cross section : Calculate w and theta firstly, and then obtain differential cross section //
	if(rm == "diff"){
	  ran3 = 0.265*randomNum.Rndm(); // For unpolarized differetial cross section, 0.265 is the maximum value of unpolarized cross section table for the first step of qf and Sn rescattering
	}
	else if(rm == "same"){
	  ran3 = 30.9*randomNum.Rndm(); // For unpolarized differetial cross section, 30.9 is the maximum value of unpolarized cross section table for the first step of all channels
	}

	temp1 = *pQp2;
	temp2 = beam;
	betax = -(beam.Px()+target.Px())/(beam.Energy()+target.Energy());
	betay = -(beam.Py()+target.Py())/(beam.Energy()+target.Energy());
	betaz = -(beam.Pz()+target.Pz())/(beam.Energy()+target.Energy());
	temp1.Boost(betax, betay, betaz);
	temp2.Boost(betax, betay, betaz);
	cmpol = (180/3.141592653)*acos((temp1.Px()*temp2.Px()+temp1.Py()*temp2.Py()+temp1.Pz()*temp2.Pz())/(temp1.P()*temp2.P()));

	csval=calcdcs_ks(vect_ks,tinW_ks,cmpol, csstat);
	loc_wval=tinW_ks;
	loc_theta=cmpol;          
	// Differential cross section : Calculate w and theta firstly, and then obtain differential cross section 
	
	if(cs1=="ncs"){
	  Sbeam=*pQp1;
	  SW = Sbeam+Starget;
	  if (event_finalsn.SetDecay (SW, 2, Smasses_sn)){	
	    Sweight = event_finalsn.Generate ();
	    pRp1 = event_finalsn.GetDecay (0);
	    pRp2 = event_finalsn.GetDecay (1);
	    Nevents++;
	    
	    // Root Stuff //
	    if (WillBeRootOutput){
	      mytree->Fill ();
	    }
	    // Root Stuff //	
	  }
	}

	if(cs1=="ucs" || cs1=="ccs" || cs1 == "lcs"){
	  if(csval<ran3) continue;
	  Sbeam=*pQp1;
	  SW = Sbeam+Starget;
	  if (event_finalsn.SetDecay (SW, 2, Smasses_sn)){	
	    Sweight = event_finalsn.Generate ();
	    pRp1 = event_finalsn.GetDecay (0);
	    pRp2 = event_finalsn.GetDecay (1);
	    Nevents++;
	    
	    // Root Stuff //
	    if (WillBeRootOutput){
	      mytree->Fill ();
	    }
	    // Root Stuff //	
	  }
	}
      } // Condition for the first step      
    } // Condition for channel
    //--------------generate event for g+d-->KSn Sn rescattering(label 4)----------------------//

    //---------------generate events for pi+ mediated reaction (label 5)------------------//
    if (rndnumber>=channel_random_left[4] && rndnumber<=channel_random_right[4]) {
      Starget.SetXYZT (-fpx, -fpy, -fpz, fEn);
      target=deuteron-Starget;
      W = beam + target;
      tinW_piplusn = 1000*W.Mag();
      label = 5;

      if (event_piplusn.SetDecay (W, 2, masses_piplusn)){
	weight = event_piplusn.Generate ();
	pQp1 = event_piplusn.GetDecay (0);
	pQp2 = event_piplusn.GetDecay (1);

	// Differential cross section : Calculate w and theta firstly, and then obtain differential cross section //
	if(rm == "diff"){
	  ran3 = 24.8*randomNum.Rndm(); // For unpolarized differetial cross section, 24.8 is the maximum value of unpolarized cross section table for the first step of pi+ mediated
	  ran4l = 23*randomNum.Rndm(); // For linearly polarized differetial cross section, 23 is the maximum value of linearly polarizied differential cross section table for the first step of pi+ mediated
	}
	else if(rm == "same"){
	  ran3 = 30.9*randomNum.Rndm(); // For unpolarized differetial cross section, 30.9 is the maximum value of unpolarized cross section table for the first step of all channels
	  ran4l = 24*randomNum.Rndm(); // For linearly polarized differetial cross section, 24 is the maximum value of linearly polarized differential cross section for channels 1, 2, 5 and 8
	}

	temp1 = *pQp1;
	temp2 = beam;
	betax = -(beam.Px()+target.Px())/(beam.Energy()+target.Energy());
	betay = -(beam.Py()+target.Py())/(beam.Energy()+target.Energy());
	betaz = -(beam.Pz()+target.Pz())/(beam.Energy()+target.Energy());
	temp1.Boost(betax, betay, betaz);
	temp2.Boost(betax, betay, betaz);
	cmpol = (180/3.141592653)*acos((temp1.Px()*temp2.Px()+temp1.Py()*temp2.Py()+temp1.Pz()*temp2.Pz())/(temp1.P()*temp2.P()));

	csval=calcdcs_piplusn(vect_piplusn,tinW_piplusn,cmpol, csstat);
	loc_wval=tinW_piplusn;
	loc_theta=cmpol;              
	// Differential cross section : Calculate w and theta firstly, and then obtain differential cross section

	if(cs1=="ncs"){
	  Sbeam=*pQp1;
	  SW = Sbeam+Starget;
	  StinW_KL=1000*SW.Mag();

	  Starget_recalc.SetXYZM(Starget.Px(),Starget.Py(),Starget.Pz(),0.938272046);
	  StinW_KL_recalc=1000*(Sbeam+Starget_recalc).Mag();

	  if (event_finalpiplusn.SetDecay (SW, 2, Smasses_lk)){
	    Sweight = event_finalpiplusn.Generate ();
	    pRp1 = event_finalpiplusn.GetDecay (0);
	    pRp2 = event_finalpiplusn.GetDecay (1);

	    ran5 = 113.3*randomNum.Rndm(); //113.3 is the maximum value of unpolarized differential cross section based on table of the second step
	    Stemp1 = *pRp2;
	    Stemp2 = Sbeam;
	    Sbetax = -(Sbeam.Px()+Starget.Px())/(Sbeam.Energy()+Starget.Energy());
	    Sbetay = -(Sbeam.Py()+Starget.Py())/(Sbeam.Energy()+Starget.Energy());
	    Sbetaz = -(Sbeam.Pz()+Starget.Pz())/(Sbeam.Energy()+Starget.Energy());
	    Stemp1.Boost(Sbetax, Sbetay, Sbetaz);
	    Stemp2.Boost(Sbetax, Sbetay, Sbetaz);
	    Scmpol = (180/3.141592653)*acos((Stemp1.Px()*Stemp2.Px()+Stemp1.Py()*Stemp2.Py()+Stemp1.Pz()*Stemp2.Pz())/(Stemp1.P()*Stemp2.P()));

	    Scsval_original=calcdcsKL_original(StinW_KL_recalc, cos(Scmpol/180*3.141592653), Scsstat);
	    Scsval_adjust=calcdcsKL_adjust(vect_SKL, StinW_KL_recalc, cos(Scmpol/180*3.141592653), Scsstat);

	    loc_Swval=StinW_KL;
	    loc_Stheta=Scmpol;


	    if(cs2=="ncs"){
	      if (event_decayppi.SetDecay (*pRp1, 2, Dmasses_ppi)){		 
		event_decayppi.Generate ();
		pDp1 = event_decayppi.GetDecay (0);
		pDp2 = event_decayppi.GetDecay (1);
		Nevents++;
	      
		// Root Stuff //
		if (WillBeRootOutput){
		  mytree->Fill ();
		}
		// Root Stuff //	
	      }
	    }

	    else if(cs2=="ocs"){
	      if(Scsval_original<ran5) continue;
	      if (event_decayppi.SetDecay (*pRp1, 2, Dmasses_ppi)){		 
		event_decayppi.Generate ();
		pDp1 = event_decayppi.GetDecay (0);
		pDp2 = event_decayppi.GetDecay (1);
		Nevents++;
	    
		// Root Stuff //
		if (WillBeRootOutput){
		  mytree->Fill ();
		}	
		// Root Stuff //
	      }
	    }

	    else if(cs2=="acs"){
	      if(Scsval_adjust<ran5) continue;
	      if (event_decayppi.SetDecay (*pRp1, 2, Dmasses_ppi)){		 
		event_decayppi.Generate ();
		pDp1 = event_decayppi.GetDecay (0);
		pDp2 = event_decayppi.GetDecay (1);
		Nevents++;
	    
		// Root Stuff //
		if (WillBeRootOutput){
		  mytree->Fill ();
		}	
		// Root Stuff //
	      }
	    }
	  }
	}

	if(cs1=="ucs" || cs1 == "ccs"){
	  if(csval<ran3) continue;
	  Sbeam=*pQp1;
	  SW = Sbeam+Starget;
	  StinW_KL=1000*SW.Mag();

	  Starget_recalc.SetXYZM(Starget.Px(),Starget.Py(),Starget.Pz(),0.938272046);
	  StinW_KL_recalc=1000*(Sbeam+Starget_recalc).Mag();

	  if (event_finalpiplusn.SetDecay (SW, 2, Smasses_lk)){
	    Sweight = event_finalpiplusn.Generate ();
	    pRp1 = event_finalpiplusn.GetDecay (0);
	    pRp2 = event_finalpiplusn.GetDecay (1);

	    ran5 = 113.3*randomNum.Rndm(); //113.3 is the maximum value of unpolarized differential cross section based on table of the second step
	    Stemp1 = *pRp2;
	    Stemp2 = Sbeam;
	    Sbetax = -(Sbeam.Px()+Starget.Px())/(Sbeam.Energy()+Starget.Energy());
	    Sbetay = -(Sbeam.Py()+Starget.Py())/(Sbeam.Energy()+Starget.Energy());
	    Sbetaz = -(Sbeam.Pz()+Starget.Pz())/(Sbeam.Energy()+Starget.Energy());
	    Stemp1.Boost(Sbetax, Sbetay, Sbetaz);
	    Stemp2.Boost(Sbetax, Sbetay, Sbetaz);
	    Scmpol = (180/3.141592653)*acos((Stemp1.Px()*Stemp2.Px()+Stemp1.Py()*Stemp2.Py()+Stemp1.Pz()*Stemp2.Pz())/(Stemp1.P()*Stemp2.P()));

	    Scsval_original=calcdcsKL_original(StinW_KL_recalc, cos(Scmpol/180*3.141592653), Scsstat);
	    Scsval_adjust=calcdcsKL_adjust(vect_SKL, StinW_KL_recalc, cos(Scmpol/180*3.141592653), Scsstat);

	    loc_Swval=StinW_KL;
	    loc_Stheta=Scmpol;

	    if(cs2=="ncs"){
	      if (event_decayppi.SetDecay (*pRp1, 2, Dmasses_ppi)){		 
		event_decayppi.Generate ();
		pDp1 = event_decayppi.GetDecay (0);
		pDp2 = event_decayppi.GetDecay (1);
		Nevents++;
	      
		// Root Stuff //
		if (WillBeRootOutput){
		  mytree->Fill ();
		}
		// Root Stuff //	
	      }
	    }

	    else if(cs2=="ocs"){
	      if(Scsval_original<ran5) continue;
	      if (event_decayppi.SetDecay (*pRp1, 2, Dmasses_ppi)){		 
		event_decayppi.Generate ();
		pDp1 = event_decayppi.GetDecay (0);
		pDp2 = event_decayppi.GetDecay (1);
		Nevents++;
	    
		// Root Stuff //
		if (WillBeRootOutput){
		  mytree->Fill ();
		}	
		// Root Stuff //
	      }
	    }

	    else if(cs2=="acs"){
	      if(Scsval_adjust<ran5) continue;
	      if (event_decayppi.SetDecay (*pRp1, 2, Dmasses_ppi)){		 
		event_decayppi.Generate ();
		pDp1 = event_decayppi.GetDecay (0);
		pDp2 = event_decayppi.GetDecay (1);
		Nevents++;
	    
		// Root Stuff //
		if (WillBeRootOutput){
		  mytree->Fill ();
		}	
		// Root Stuff //
	      }
	    }
	  }
	}
	 
	if(cs1=="lcs"){
	  Sbeam=*pQp1;
	  SW = Sbeam+Starget;
	  StinW_KL=1000*SW.Mag();

	  Starget_recalc.SetXYZM(Starget.Px(),Starget.Py(),Starget.Pz(),0.938272046);
	  StinW_KL_recalc=1000*(Sbeam+Starget_recalc).Mag();

	  if (event_finalpiplusn.SetDecay (SW, 2, Smasses_lk)){
	    Sweight = event_finalpiplusn.Generate ();
	    pRp1 = event_finalpiplusn.GetDecay (0);
	    pRp2 = event_finalpiplusn.GetDecay (1);
	      
	    if(Nlabel5%2==0) linpol=1;
	    else linpol=0;
	      
	    phi=pQp1->Phi();
	    sigma=calcs_piplusn(vect_spiplusn, beam.E(), cmpol, obsstat); // phi is defined as azimuthal angle of pi in the lab frame
	    if(linpol==1) lcsval=csval*(1-sigma*cos(2*phi));
	    else lcsval=csval*(1+sigma*cos(2*phi));
	      
	    if(lcsval<ran4l) continue;

	    ran5 = 113.3*randomNum.Rndm(); //113.3 is the maximum value of unpolarized differential cross section based on table of the second step
	    Stemp1 = *pRp2;
	    Stemp2 = Sbeam;
	    Sbetax = -(Sbeam.Px()+Starget.Px())/(Sbeam.Energy()+Starget.Energy());
	    Sbetay = -(Sbeam.Py()+Starget.Py())/(Sbeam.Energy()+Starget.Energy());
	    Sbetaz = -(Sbeam.Pz()+Starget.Pz())/(Sbeam.Energy()+Starget.Energy());
	    Stemp1.Boost(Sbetax, Sbetay, Sbetaz);
	    Stemp2.Boost(Sbetax, Sbetay, Sbetaz);
	    Scmpol = (180/3.141592653)*acos((Stemp1.Px()*Stemp2.Px()+Stemp1.Py()*Stemp2.Py()+Stemp1.Pz()*Stemp2.Pz())/(Stemp1.P()*Stemp2.P()));

	    Scsval_original=calcdcsKL_original(StinW_KL_recalc, cos(Scmpol/180*3.141592653), Scsstat);
	    Scsval_adjust=calcdcsKL_adjust(vect_SKL, StinW_KL_recalc, cos(Scmpol/180*3.141592653), Scsstat);

	    loc_Swval=StinW_KL;
	    loc_Stheta=Scmpol;

	    if(cs2=="ncs"){
	      if (event_decayppi.SetDecay (*pRp1, 2, Dmasses_ppi)){		 
		event_decayppi.Generate ();
		pDp1 = event_decayppi.GetDecay (0);
		pDp2 = event_decayppi.GetDecay (1);
		Nlabel5++;
		Nevents++; 

		// Root Stuff //
		if (WillBeRootOutput){
		  mytree->Fill ();
		}
	      }
	    }

	    else if(cs2=="ocs"){
	      if(Scsval_original<ran5) continue;
	      if (event_decayppi.SetDecay (*pRp1, 2, Dmasses_ppi)){		 
		event_decayppi.Generate ();
		pDp1 = event_decayppi.GetDecay (0);
		pDp2 = event_decayppi.GetDecay (1);
		Nlabel5++;
		Nevents++;
	    
		// Root Stuff //
		if (WillBeRootOutput){
		  mytree->Fill ();
		}	
		// Root Stuff //
	      }
	    }

	    else if(cs2=="acs"){
	      if(Scsval_adjust<ran5) continue;
	      if (event_decayppi.SetDecay (*pRp1, 2, Dmasses_ppi)){		 
		event_decayppi.Generate ();
		pDp1 = event_decayppi.GetDecay (0);
		pDp2 = event_decayppi.GetDecay (1);
		Nlabel5++;
		Nevents++;
	    
		// Root Stuff //
		if (WillBeRootOutput){
		  mytree->Fill ();
		}	
		// Root Stuff //
	      }
	    }
	  }
	}
      } // Condition for the first step      
    } // Condition for channel
    //---------------generate events for pi+ mediated reaction (label 5)------------------//

    //--------------generate event for g+d-->KSSn (label 6)----------------------//
    if (rndnumber>=channel_random_left[5] && rndnumber<=channel_random_right[5]) { 
      Starget.SetXYZT (-fpx, -fpy, -fpz, fEn);
      target=deuteron-Starget;
      W = beam + target;
      tinW_ssk = 1000*W.Mag();
      label = 6;

      if (event_ssk.SetDecay (W, 2, masses_ssk)){
	weight = event_ssk.Generate ();
	pQp1 = event_ssk.GetDecay (0);
	pQp2 = event_ssk.GetDecay (1);

	temp1 = *pQp2;
	temp2 = beam;
	betax = -(beam.Px()+target.Px())/(beam.Energy()+target.Energy());
	betay = -(beam.Py()+target.Py())/(beam.Energy()+target.Energy());
	betaz = -(beam.Pz()+target.Pz())/(beam.Energy()+target.Energy());
	temp1.Boost(betax, betay, betaz);
	temp2.Boost(betax, betay, betaz);
	cmpol = (180/3.141592653)*acos((temp1.Px()*temp2.Px()+temp1.Py()*temp2.Py()+temp1.Pz()*temp2.Pz())/(temp1.P()*temp2.P()));
	
	loc_wval=tinW_ssk;
	loc_theta=cmpol;   

	if(cs1=="ncs" || cs1=="ucs" || cs1=="ccs" || cs1 == "lcs"){
	  if (event_finallpi.SetDecay (*pQp1, 2, Smasses_lpi)){	
	    Sweight = event_finallpi.Generate ();
	    pDp1 = event_finallpi.GetDecay (0);
	    pDp2 = event_finallpi.GetDecay (1);
	    Nevents++;
	    
	    // Root Stuff //
	    if (WillBeRootOutput){
	      mytree->Fill ();
	    }
	    // Root Stuff //
	  }	
	}    
      } // Condition for the first step      
    } // Condition for channel
    //--------------generate event for g+d-->KSSn (label 6)----------------------//     

    //-----------generate events for g+d-->KSSminusp mediated reaction (label=7)---------//
    if (rndnumber>=channel_random_left[6] && rndnumber<=channel_random_right[6]) {
      Starget.SetXYZT (-fpx, -fpy, -fpz, fEp);
      target=deuteron-Starget;
      W = beam + target;
      tinW_ssminusk = 1000*W.Mag();
      label = 7;

      if (event_ssminusk.SetDecay (W, 2, masses_ssminusk)){
	weight = event_ssminusk.Generate ();
	pQp1 = event_ssminusk.GetDecay (0);
	pQp2 = event_ssminusk.GetDecay (1);

	temp1 = *pQp2;
	temp2 = beam;
	betax = -(beam.Px()+target.Px())/(beam.Energy()+target.Energy());
	betay = -(beam.Py()+target.Py())/(beam.Energy()+target.Energy());
	betaz = -(beam.Pz()+target.Pz())/(beam.Energy()+target.Energy());
	temp1.Boost(betax, betay, betaz);
	temp2.Boost(betax, betay, betaz);
	cmpol = (180/3.141592653)*acos((temp1.Px()*temp2.Px()+temp1.Py()*temp2.Py()+temp1.Pz()*temp2.Pz())/(temp1.P()*temp2.P()));

	loc_wval=tinW_ssminusk;
	loc_theta=cmpol;   
	
	if(cs1=="ncs" || cs1=="ucs" || cs1=="ccs" || cs1 == "lcs"){
	  if (event_finallpiminus.SetDecay (*pQp1, 2, Smasses_lpiminus)){	
	    Sweight = event_finallpiminus.Generate ();
	    pDp1 = event_finallpiminus.GetDecay (0);
	    pDp2 = event_finallpiminus.GetDecay (1);
	    Nevents++;
	    
	    // Root Stuff //
	    if (WillBeRootOutput){
	      mytree->Fill ();
	    }
	    // Root Stuff //
	  }
	}
      } // Condition for the first step      
    } // Condition for channel
    //-----------generate events for g+d-->KSSminusp mediated reaction (label=7)---------//

    //----------------------generate event g+d-->KLn qf mechanism (label 8)----------------------//
    if (rndnumber>=channel_random_left[7] && rndnumber<=channel_random_right[7]) {
      Starget.SetXYZT (-fpx, -fpy, -fpz, fEn);
      target=deuteron-Starget;
      W = beam + target;
      tinW_kl = 1000*W.Mag();
      label = 8;

      if (event.SetDecay (W, 2, masses_lk)){
	weight = event.Generate ();
	pQp1 = event.GetDecay (0);
	pQp2 = event.GetDecay (1);
	// Differential cross section : Calculate w and theta firstly, and then obtain differential cross section //
	if(rm == "diff"){
	  ran3 = 0.355*randomNum.Rndm(); // For unpolarized differetial cross section, 0.355 is the maximum value of unpolarized cross section table for the first step of qf, Kn and Ln rescattering
	  ran4l = 0.80*randomNum.Rndm(); // For linearly polarized differetial cross section, 0.80 is the maximum value of linearly polarized differential cross section for the first step of Kn rescattering
	}
	else if(rm == "same"){
	  ran3 = 30.9*randomNum.Rndm(); // For unpolarized differetial cross section, 30.9 is the maximum value of unpolarized cross section table for the first step of all channels
	  ran4l = 24*randomNum.Rndm(); // For linearly polarized differetial cross section, 24 is the maximum value of linearly polarized differential cross section for channels 1, 2, 5 and 8
	}


	ran4c = 0.58*randomNum.Rndm(); // For circularly polarized differetial cross section, 0.58 is the maximum value of circularly differential polarized cross section for the first step of Kn rescattering

	temp1 = *pQp2;
	temp2 = beam;
	betax = -(beam.Px()+target.Px())/(beam.Energy()+target.Energy());
	betay = -(beam.Py()+target.Py())/(beam.Energy()+target.Energy());
	betaz = -(beam.Pz()+target.Pz())/(beam.Energy()+target.Energy());
	temp1.Boost(betax, betay, betaz);
	temp2.Boost(betax, betay, betaz);
	cmpol = (180/3.141592653)*acos((temp1.Px()*temp2.Px()+temp1.Py()*temp2.Py()+temp1.Pz()*temp2.Pz())/(temp1.P()*temp2.P()));

	csval=calcdcs_kl(vect_kl,tinW_kl,cmpol, csstat);
	loc_wval=tinW_kl;
	loc_theta=cmpol;
	// Differential cross section : Calculate w and theta firstly, and then obtain differential cross section

	if(cs1=="ncs"){	
	  if (event_decayppi.SetDecay (*pQp1, 2, Dmasses_ppi)){		 
	    event_decayppi.Generate ();
	    pDp1 = event_decayppi.GetDecay (0);
	    pDp2 = event_decayppi.GetDecay (1);

	    Nevents++;
 
	    // Root Stuff //
	    if (WillBeRootOutput){
	      mytree->Fill ();
	    }
	    // Root Stuff //
	  }
	}   

	if(cs1=="ucs"){
	  if(csval<ran3) continue;	
	  if (event_decayppi.SetDecay (*pQp1, 2, Dmasses_ppi)){		 
	    event_decayppi.Generate ();
	    pDp1 = event_decayppi.GetDecay (0);
	    pDp2 = event_decayppi.GetDecay (1);

	    Nevents++;
 
	    // Root Stuff //
	    if (WillBeRootOutput){
	      mytree->Fill ();
	    }
	    // Root Stuff //
	  }
	}

	if(cs1=="ccs"){
	  if (event_decayppi.SetDecay (*pQp1, 2, Dmasses_ppi)){		 
	    event_decayppi.Generate ();
	    pDp1 = event_decayppi.GetDecay (0);
	    pDp2 = event_decayppi.GetDecay (1);

	    Pphot=beam;
	    Pphot.Boost(betax, betay, betaz);
	    Pkaon=temp1;
	    Vphot.SetXYZ(Pphot.Px(),Pphot.Py(),Pphot.Pz());
	    Vz=Vphot.Unit();
	    Vkaon.SetXYZ(Pkaon.Px(),Pkaon.Py(),Pkaon.Pz());
	    Vy=Vz.Cross(Vkaon);
	    Vy=Vy.Unit();
	    Vx=Vy.Cross(Vz);

	    Plamb=*pQp1;
	    beta2x=-Plamb.Px()/Plamb.E();    beta2y=-Plamb.Py()/Plamb.E();    beta2z=-Plamb.Pz()/Plamb.E();
	    Pprot=*pDp1;
	    Pprot.Boost(beta2x,beta2y,beta2z);
	    Vprot.SetXYZ(Pprot.Px(),Pprot.Py(),Pprot.Pz());

	    costhetax=Vprot*Vx/Vprot.Mag();
	    costhetay=Vprot*Vy/Vprot.Mag();
	    costhetaz=Vprot*Vz/Vprot.Mag();

	    vect_cxczpy.clear();
	    vect_cxczpy=calccxczpy_kl(beam.E(),cos(cmpol/180*3.141592653),obsstat,e_left, e_right, e_aver, ctk_left, ctk_right, ctk_aver, cxTable, czTable, pyTable);
	    cx=vect_cxczpy[0];
	    cz=vect_cxczpy[1];
	    py=vect_cxczpy[2];

	    if(Nlabel8%2==0) cirpol=1;
	    else cirpol=0;

	    if(cirpol==1) ccsval=csval*(1+0.642*costhetax*cx+0.642*costhetaz*cz+0.642*costhetay*py);
	    else ccsval=csval*(1-0.642*costhetax*cx-0.642*costhetaz*cz+0.642*costhetay*py);

	    if(ccsval<ran4c) continue;
	    Nlabel8++;

	    Nevents++; 
	    // Root Stuff //
	    if (WillBeRootOutput){
	      mytree->Fill ();
	    }
	    // Root Stuff //
	  }	    
	}


	if(cs1=="lcs"){
	  if (event_decayppi.SetDecay (*pQp1, 2, Dmasses_ppi)){		 
	    event_decayppi.Generate ();
	    pDp1 = event_decayppi.GetDecay (0);
	    pDp2 = event_decayppi.GetDecay (1);

	    Pphot=beam;
	    Pphot.Boost(betax, betay, betaz);
	    Pkaon=temp1;
	    Vphot.SetXYZ(Pphot.Px(),Pphot.Py(),Pphot.Pz());
	    Vz=Vphot.Unit();
	    Vkaon.SetXYZ(Pkaon.Px(),Pkaon.Py(),Pkaon.Pz());
	    Vy=Vz.Cross(Vkaon);
	    Vy=Vy.Unit();
	    Vx=Vy.Cross(Vz);

	    Plamb=*pQp1;
	    beta2x=-Plamb.Px()/Plamb.E();    beta2y=-Plamb.Py()/Plamb.E();    beta2z=-Plamb.Pz()/Plamb.E();
	    Pprot=*pDp1;
	    Pprot.Boost(beta2x,beta2y,beta2z);
	    Vprot.SetXYZ(Pprot.Px(),Pprot.Py(),Pprot.Pz());

	    costhetax=Vprot*Vx/Vprot.Mag();
	    costhetay=Vprot*Vy/Vprot.Mag();
	    costhetaz=Vprot*Vz/Vprot.Mag();

	    vect_stoxozpy.clear();
	    vect_stoxozpy=calcSTOxOzPforGPKL(tinW_kl,cmpol,obsstat,vectSTOxOzPforGPKL);

	    sigma=vect_stoxozpy[0];
	    t=vect_stoxozpy[1];
	    ox=vect_stoxozpy[2];
	    oz=vect_stoxozpy[3];
	    py=vect_stoxozpy[4];

	    if(Nlabel8%2==0) linpol=1;
	    else linpol=0;

	    phi=pQp1->Phi();
	    if(linpol==1) lcsval=csval*(1-sigma*cos(2*phi)-0.642*costhetay*t*cos(2*phi)+0.642*costhetax*ox*sin(2*phi)+0.642*costhetaz*oz*sin(2*phi)+0.642*costhetay*py);
	    else lcsval=csval*(1+sigma*cos(2*phi)+0.642*costhetay*t*cos(2*phi)-0.642*costhetax*ox*sin(2*phi)-0.642*costhetaz*oz*sin(2*phi)+0.642*costhetay*py);

	    if(lcsval<ran4l) continue;
	    Nlabel8++;

	    Nevents++; 
	    // Root Stuff //
	    if (WillBeRootOutput){
	      mytree->Fill ();
	    }
	    // Root Stuff //
	  }	    
	}	
      } // Condition for the first step      
    } // Condition for channel  
    //----------------------generate event g+d-->KLn qf mechanism (label 8)----------------------//

    //--------------generate event for g+d-->KSn qf mechanism (label 9)----------------------//
    if (rndnumber>=channel_random_left[8] && rndnumber<=channel_random_right[8]) {
      Starget.SetXYZT (-fpx, -fpy, -fpz, fEn);
      target=deuteron-Starget;
      W = beam + target;
      tinW_ks = 1000*W.Mag();
      label = 9;

      if (event_sk.SetDecay (W, 2, masses_sk)){
	weight = event_sk.Generate ();
	pQp1 = event_sk.GetDecay (0);
	pQp2 = event_sk.GetDecay (1);

	// Differential cross section : Calculate w and theta firstly, and then obtain differential cross section //
	if(rm == "diff"){
	  ran3 = 0.265*randomNum.Rndm(); // For unpolarized differetial cross section, 0.265 is the maximum value of unpolarized cross section table for the first step of qf and Sn rescattering
	}
	else if(rm == "same"){
	  ran3 = 30.9*randomNum.Rndm(); // For unpolarized differetial cross section, 30.9 is the maximum value of unpolarized cross section table for the first step of all channels
	}

	temp1 = *pQp2;
	temp2 = beam;
	betax = -(beam.Px()+target.Px())/(beam.Energy()+target.Energy());
	betay = -(beam.Py()+target.Py())/(beam.Energy()+target.Energy());
	betaz = -(beam.Pz()+target.Pz())/(beam.Energy()+target.Energy());
	temp1.Boost(betax, betay, betaz);
	temp2.Boost(betax, betay, betaz);
	cmpol = (180/3.141592653)*acos((temp1.Px()*temp2.Px()+temp1.Py()*temp2.Py()+temp1.Pz()*temp2.Pz())/(temp1.P()*temp2.P()));

	csval=calcdcs_ks(vect_ks,tinW_ks,cmpol, csstat);
	loc_wval=tinW_ks;
	loc_theta=cmpol;          
	// Differential cross section : Calculate w and theta firstly, and then obtain differential cross section 
	
	if(cs1=="ncs"){
	  Nevents++;
	    
	  // Root Stuff //
	  if (WillBeRootOutput){
	    mytree->Fill ();
	  }
	  // Root Stuff //	
	}

	if(cs1=="ucs" || cs1=="ccs" || cs1 == "lcs"){
	  if(csval<ran3) continue;
	  Nevents++;
	    
	  // Root Stuff //
	  if (WillBeRootOutput){
	    mytree->Fill ();
	  }
	  // Root Stuff //	
	}
      } // Condition for the first step      
    } // Condition for channel
    //--------------generate event for g+d-->KSn qf mechanism (label 9)----------------------//

  }       
  //// Loop end //// 

  // Root Stuff //
  if (WillBeRootOutput)
    {
      RootOut.Write ();
    }
  // Root Stuff //
    
  cout << "Number of events processed: " << Nevents << endl;
    
  return 0;  
}
/* ----------- main function ---------------- */

/* ----------- Calculate differential cross section ---------------- */
// The interval of w is 5, and the interval of theta is 2.
// For the outside of w range of table lists, we use dcs of lower w as the dcs if w is smaller than minimum w, and dcs of upper w as the dcs if w is larger than maximum w
//// Function "readdcs" to read differential cross section of different table lists ////
vector<double> readdcs(string filename){
  ifstream input(filename.c_str());

  vector<double> vect_dcs;
  vector<string> table;
  string line,tempstring,subline,subline1,subline2,subline3,subline4;
  int tindex1,tindex2,tindex3,tindex4;
  double dcs;

  while(getline(input,line)) {
    if (line.empty() || line.find("W")!=string::npos) continue;
    table.push_back(line);

    tindex1 = line.find_first_not_of(" ");
    subline=line.substr(tindex1);
    tindex2 = subline.find_first_of(" ");
    subline1 = subline.substr(0, tindex2);
     
    subline2 = subline.substr(tindex2);
    tindex3 = subline2.find_first_not_of(" ");
    subline3 = subline2.substr(tindex3);
    tindex4 = subline3.find_first_of(" ");
    subline4 = subline3.substr(0, tindex4);
    istringstream(subline4)>>dcs;
    vect_dcs.push_back(dcs);
  }
  return vect_dcs;
}
//// Function "readdcs" to read differential cross section of different table lists ////

//// Function "calcdcs_pin" to calculate cross section for gamma+n -> pi0+n ////
double calcdcs_pin(vector<double> vect, double w, double theta, int &status){
  double w_min=1075, w_max=2000;
  int w_order,theta_order;
  double wmin,thetamin;
  int index_wmin_thetamin,index_wmin_thetamax,index_wmax_thetamin,index_wmax_thetamax;
  double wmin_dcs,wmax_dcs,calcdcs;

  int index_thetamin,index_thetamax;
    
  if (w<w_min){
    theta_order=theta/2; //order from 0
    thetamin=theta_order*2;
    index_thetamin=theta_order;
    index_thetamax=theta_order+1;
        
    calcdcs=(vect[index_thetamax]-vect[index_thetamin])/2*(theta-thetamin)+vect[index_thetamin];
        
    status=0;
    return calcdcs;
  }
  else if(w>w_max){
            
    theta_order=theta/2; //order from 0
    thetamin=theta_order*2;
    index_thetamin=185*91+theta_order;
    index_thetamax=185*91+theta_order+1;
    calcdcs=(vect[index_thetamax]-vect[index_thetamin])/2*(theta-thetamin)+vect[index_thetamin];
    status=0;
    return calcdcs;
  }

  else{
    w_order=(w-w_min)/5; //order from 0
    theta_order=theta/2; //order from 0
    wmin=w_min+w_order*5;
    thetamin=theta_order*2;
    index_wmin_thetamin=w_order*91+theta_order;
    index_wmin_thetamax=w_order*91+theta_order+1;
    index_wmax_thetamin=(w_order+1)*91+theta_order;
    index_wmax_thetamax=(w_order+1)*91+theta_order+1;

    wmin_dcs=(vect[index_wmin_thetamax]-vect[index_wmin_thetamin])/2*(theta-thetamin)+vect[index_wmin_thetamin];
    wmax_dcs=(vect[index_wmax_thetamax]-vect[index_wmax_thetamin])/2*(theta-thetamin)+vect[index_wmax_thetamin];
    calcdcs = (wmax_dcs-wmin_dcs)/5*(w-wmin)+wmin_dcs;
    status=1;

    return calcdcs;
  }
}
//// Function "calcdcs_pin" to calculate cross section for gamma+n -> pi0+n ////

//// Function "calcdcs_kl" to calculate cross section for gamma+p -> kaon+Lambda ////
double calcdcs_kl(vector<double> vect, double w, double theta, int &status){
  double w_min=1610, w_max=2200;
  int w_order,theta_order;
  double wmin,thetamin;
  int index_wmin_thetamin,index_wmin_thetamax,index_wmax_thetamin,index_wmax_thetamax;
  double wmin_dcs,wmax_dcs,calcdcs;
  int index_thetamin,index_thetamax;
    
  if (w<w_min){
    theta_order=theta/2; //order from 0
    thetamin=theta_order*2;
    index_thetamin=theta_order;
    index_thetamax=theta_order+1;
        
    calcdcs=(vect[index_thetamax]-vect[index_thetamin])/2*(theta-thetamin)+vect[index_thetamin];
        
    status=0;
    return calcdcs;
  }
  else if(w>w_max){
            
    theta_order=theta/2; //order from 0
    thetamin=theta_order*2;
    index_thetamin=118*91+theta_order;
    index_thetamax=118*91+theta_order+1;
            
    calcdcs=(vect[index_thetamax]-vect[index_thetamin])/2*(theta-thetamin)+vect[index_thetamin];
    status=0;
    return calcdcs;
  }

        
  else{
    w_order=(w-w_min)/5; //order from 0
    theta_order=theta/2; //order from 0
    wmin=w_min+w_order*5;
    thetamin=theta_order*2;
    index_wmin_thetamin=w_order*91+theta_order;
    index_wmin_thetamax=w_order*91+theta_order+1;
    index_wmax_thetamin=(w_order+1)*91+theta_order;
    index_wmax_thetamax=(w_order+1)*91+theta_order+1;

    wmin_dcs=(vect[index_wmin_thetamax]-vect[index_wmin_thetamin])/2*(theta-thetamin)+vect[index_wmin_thetamin];
    wmax_dcs=(vect[index_wmax_thetamax]-vect[index_wmax_thetamin])/2*(theta-thetamin)+vect[index_wmax_thetamin];
    calcdcs = (wmax_dcs-wmin_dcs)/5*(w-wmin)+wmin_dcs;
    status=1;

    return calcdcs;
  }
}
//// Function "calcdcs_kl" to calculate cross section for gamma+p -> kaon+Lambda ////

//// Function "calcdcs_ks" to calculate cross section for gamma+p -> kaon+Sigma ////
// This case is a little bit special because the interval between first and second w is 3.
double calcdcs_ks(vector<double> vect, double w, double theta, int &status){
  double w_min=1690, w_max=2200;
  int w_order,theta_order;
  double wmin,thetamin;
  int index_wmin_thetamin,index_wmin_thetamax,index_wmax_thetamin,index_wmax_thetamax;
  double wmin_dcs,wmax_dcs,calcdcs;

  int index_thetamin,index_thetamax;
    
  if (w<1687){
    theta_order=theta/2; //order from 0
    thetamin=theta_order*2;
    index_thetamin=theta_order;
    index_thetamax=theta_order+1;
        
    calcdcs=(vect[index_thetamax]-vect[index_thetamin])/2*(theta-thetamin)+vect[index_thetamin];
        
    status=0;
    return calcdcs;
  }
  else if(w>w_max){
            
    theta_order=theta/2; //order from 0
    thetamin=theta_order*2;
    index_thetamin=102*91+theta_order;
    index_thetamax=102*91+theta_order+1;
            
    calcdcs=(vect[index_thetamax]-vect[index_thetamin])/2*(theta-thetamin)+vect[index_thetamin];
    status=0;
    return calcdcs;
  }
    
    
    
  else if(w>1687 && w<w_min){
    w_order=0; //order from 0
    theta_order=theta/2; //order from 0
    wmin=1687;

    thetamin=theta_order*2;
    index_wmin_thetamin=w_order*91+theta_order;
    index_wmin_thetamax=w_order*91+theta_order+1;
    index_wmax_thetamin=(w_order+1)*91+theta_order;
    index_wmax_thetamax=(w_order+1)*91+theta_order+1;

    wmin_dcs=(vect[index_wmin_thetamax]-vect[index_wmin_thetamin])/2*(theta-thetamin)+vect[index_wmin_thetamin];
    wmax_dcs=(vect[index_wmax_thetamax]-vect[index_wmax_thetamin])/2*(theta-thetamin)+vect[index_wmax_thetamin];
    calcdcs = (wmax_dcs-wmin_dcs)/3*(w-wmin)+wmin_dcs;
    status=1;

    return calcdcs;
  }

  else{
    w_order=(w-w_min)/5; //order from 0
    theta_order=theta/2; //order from 0
    wmin=w_min+w_order*5;

    thetamin=theta_order*2;

    index_wmin_thetamin=(w_order+1)*91+theta_order;
    index_wmin_thetamax=(w_order+1)*91+theta_order+1;
    index_wmax_thetamin=(w_order+2)*91+theta_order;
    index_wmax_thetamax=(w_order+2)*91+theta_order+1;

    wmin_dcs=(vect[index_wmin_thetamax]-vect[index_wmin_thetamin])/2*(theta-thetamin)+vect[index_wmin_thetamin];
    wmax_dcs=(vect[index_wmax_thetamax]-vect[index_wmax_thetamin])/2*(theta-thetamin)+vect[index_wmax_thetamin];
    calcdcs = (wmax_dcs-wmin_dcs)/5*(w-wmin)+wmin_dcs;

    return calcdcs;
  }
}
//// Function "calcdcs_ks" to calculate cross section for gamma+p -> kaon+Sigma ////

//// Function "calcdcs_piplusn" to calculate cross section for gamma+p -> piplus+n ////
double calcdcs_piplusn(vector<double> vect, double w, double theta, int &status){
  double w_min=1080, w_max=2000;
  int w_order,theta_order;
  double wmin,thetamin;
  int index_wmin_thetamin,index_wmin_thetamax,index_wmax_thetamin,index_wmax_thetamax;
  double wmin_dcs,wmax_dcs,calcdcs;

  int index_thetamin,index_thetamax;

  if (w<w_min){
    theta_order=theta/2; //order from 0
    thetamin=theta_order*2;

    index_thetamin=theta_order;
    index_thetamax=theta_order+1;
      
    calcdcs=(vect[index_thetamax]-vect[index_thetamin])/2*(theta-thetamin)+vect[index_thetamin];
      
    status=0;
    return calcdcs;
  }
  else if(w>w_max){
          
    theta_order=theta/2; //order from 0
    thetamin=theta_order*2;

    index_thetamin=184*91+theta_order;
    index_thetamax=184*91+theta_order+1;
          
    calcdcs=(vect[index_thetamax]-vect[index_thetamin])/2*(theta-thetamin)+vect[index_thetamin];
    status=0;
    return calcdcs;
  }
  else{
    w_order=(w-w_min)/5; //order from 0
    theta_order=theta/2; //order from 0
    wmin=w_min+w_order*5;

    thetamin=theta_order*2;

    index_wmin_thetamin=w_order*91+theta_order;
    index_wmin_thetamax=w_order*91+theta_order+1;
    index_wmax_thetamin=(w_order+1)*91+theta_order;
    index_wmax_thetamax=(w_order+1)*91+theta_order+1;

    wmin_dcs=(vect[index_wmin_thetamax]-vect[index_wmin_thetamin])/2*(theta-thetamin)+vect[index_wmin_thetamin];
    wmax_dcs=(vect[index_wmax_thetamax]-vect[index_wmax_thetamin])/2*(theta-thetamin)+vect[index_wmax_thetamin];
    calcdcs = (wmax_dcs-wmin_dcs)/5*(w-wmin)+wmin_dcs;
    status=1;
    return calcdcs;
  }
}
//// Function "calcdcs_piplusn" to calculate cross section for gamma+p -> piplus+n ////

//// Function "readcxczpy_kl" to read cx, cz and py from a circular polarization observable table for gamma+p -> kaon+Lambda ////
void readcxczpy_kl(string filename, double e_left[16][10], double e_right[16][10], double e_aver[16][10], double ctk_left[16][10], double ctk_right[16][10], double ctk_aver[16][10], double cx[16][10], double cz[16][10], double py[16][10]){
  ifstream infile(filename.c_str());
  const int ebinnum=16;
  const int ctkbinnum=10;
  for(int i=0;i<ebinnum;i++){
    for(int j=0;j<ctkbinnum;j++){
      infile>>e_left[i][j]>>e_right[i][j]>>e_aver[i][j]
	    >>ctk_left[i][j]>>ctk_right[i][j]>>ctk_aver[i][j]
	    >>cx[i][j]>>cz[i][j]>>py[i][j];
    }
  }
  infile.close();
}
//// Function "readcxczpy_kl" to read cx, cz and py from a circular polarization observable table for gamma+p -> kaon+Lambda ////

//// Function "calccxczpy_kl" to calculate cx, cz and py for gamma+p -> kaon+Lambda ////
vector<double> calccxczpy_kl(double e, double ctk, int &obs_status, double e_left[16][10], double e_right[16][10], double e_aver[16][10], double ctk_left[16][10], double ctk_right[16][10], double ctk_aver[16][10], double cx[16][10], double cz[16][10], double py[16][10]){
  //ifstream infile("obs_kl.txt");
  const int ebinnum=16;
  const int ctkbinnum=10;
  double calc_cx,calc_cz,calc_py;
  double e_min=0.9, e_max=2.6;
  double weight;
    
  if (e<=e_min){
    for(int j=0;j<ctkbinnum;j++){
      if(ctk>= ctk_left[0][j] &&  ctk<ctk_right[0][j]){
	calc_cx=cx[0][j];
	calc_cz=cz[0][j];
	calc_py=py[0][j];
	break;
      }
    }    
    obs_status=0;
    vector<double> vect;
    vect.push_back(calc_cx);
    vect.push_back(calc_cz);
    vect.push_back(calc_py);
    return vect;
  }

  else if (e>e_max){
    for(int j=0;j<ctkbinnum;j++){
      if(ctk>= ctk_left[ebinnum-1][j] &&  ctk<ctk_right[ebinnum-1][j]){
	calc_cx=cx[ebinnum-1][j];
	calc_cz=cz[ebinnum-1][j];
	calc_py=py[ebinnum-1][j];
	break;
      }
    }    
    obs_status=0;
    vector<double> vect;
    vect.push_back(calc_cx);
    vect.push_back(calc_cz);
    vect.push_back(calc_py);
    return vect;
  }

  else{
    for(int i=0;i<ebinnum;i++){
      if(e>=e_left[i][0] && e<e_right[i][0]){
	if(ctk<=ctk_aver[i][0]){
	  calc_cx=cx[i][0];
	  calc_cz=cz[i][0];
	  calc_py=py[i][0];
	  break;
	}
 
	else if(ctk>ctk_aver[i][ctkbinnum-1]){
	  calc_cx=cx[i][ctkbinnum-1];
	  calc_cz=cz[i][ctkbinnum-1];
	  calc_py=py[i][ctkbinnum-1];
	  break;
	}

	else{
	  for(int j=0;j<ctkbinnum;j++){
	    if(ctk>=ctk_aver[i][j] && ctk<ctk_aver[i][j+1]){
	      weight=(ctk-ctk_aver[i][j])/(ctk_aver[i][j+1]-ctk);
	      calc_cx=(weight*cx[i][j+1]+cx[i][j])/(weight+1);
	      calc_cz=(weight*cz[i][j+1]+cz[i][j])/(weight+1);
	      calc_py=(weight*py[i][j+1]+py[i][j])/(weight+1);
	      break;
	    }
	  }
	}
      }
    }

    obs_status=1;
    vector<double> vect;
    vect.push_back(calc_cx);
    vect.push_back(calc_cz);
    vect.push_back(calc_py);
    return vect;
  }
}
//// Function "calccxczpy_kl" to calculate cx, cz and py for gamma+p -> kaon+Lambda ////

vector<double> reads(string filename){
  ifstream input(filename.c_str());

  vector<double> vect_s;
  vector<string> table;
  string line,tempstring,subline,subline1,subline2,subline3,subline4;
  int tindex1,tindex2,tindex3,tindex4;
  double s;

  while(getline(input,line)) {
    if (line.empty() || line.find("E=")!=string::npos) continue;
    table.push_back(line);

    tindex1 = line.find_first_not_of(" ");
    subline=line.substr(tindex1);
    tindex2 = subline.find_first_of(" ");
    subline1 = subline.substr(0, tindex2);
    subline2 = subline.substr(tindex2);
    tindex3 = subline2.find_first_not_of(" ");
    subline3 = subline2.substr(tindex3);
    tindex4 = subline3.find_first_of(" ");
    subline4 = subline3.substr(0, tindex4);
    istringstream(subline4)>>s;
    vect_s.push_back(s);
  }
  return vect_s;
}
//// Function "reads" to read s of different table lists ////

//// Function "calcs_pin" to calculate cross section for gamma+n -> pi0+n ////
double calcs_pin(vector<double> vect, double e, double theta, int &obs_status){
  double e_min=900, e_max=2600;
  int e_order,theta_order;
  double emin,thetamin;
  int index_emin_thetamin,index_emin_thetamax,index_emax_thetamin,index_emax_thetamax;
  double emin_s,emax_s,calcs;

  int index_thetamin,index_thetamax;
    
  if (e<e_min){
    theta_order=theta/3; //order from 0
    thetamin=theta_order*3;
    index_thetamin=theta_order;
    index_thetamax=theta_order+1;
        
    calcs=(vect[index_thetamax]-vect[index_thetamin])/3*(theta-thetamin)+vect[index_thetamin];
        
    obs_status=0;
    return calcs;
  }
  else if(e>e_max){
            
    theta_order=theta/3; //order from 0
    thetamin=theta_order*3;
    index_thetamin=85*61+theta_order;
    index_thetamax=85*61+theta_order+1;
            
    calcs=(vect[index_thetamax]-vect[index_thetamin])/3*(theta-thetamin)+vect[index_thetamin];
    obs_status=0;
    return calcs;
  }

  else{
    e_order=(e-e_min)/20; //order from 0
    theta_order=theta/3; //order from 0
    emin=e_min+e_order*20;
    thetamin=theta_order*3;
    index_emin_thetamin=e_order*61+theta_order;
    index_emin_thetamax=e_order*61+theta_order+1;
    index_emax_thetamin=(e_order+1)*61+theta_order;
    index_emax_thetamax=(e_order+1)*61+theta_order+1;

    emin_s=(vect[index_emin_thetamax]-vect[index_emin_thetamin])/3*(theta-thetamin)+vect[index_emin_thetamin];
    emax_s=(vect[index_emax_thetamax]-vect[index_emax_thetamin])/3*(theta-thetamin)+vect[index_emax_thetamin];
    calcs = (emax_s-emin_s)/20*(e-emin)+emin_s;
    obs_status=1;

    return calcs;
  }
}
//// Function "calcs_pin" to calculate cross section for gamma+n -> pi+n ////

//// Function "calcs_piplusn" to calculate cross section for gamma+p -> piplus+n ////
double calcs_piplusn(vector<double> vect, double e, double theta, int &obs_status){
  double e_min=900, e_max=2600;
  int e_order,theta_order;
  double emin,thetamin;
  int index_emin_thetamin,index_emin_thetamax,index_emax_thetamin,index_emax_thetamax;
  double emin_s,emax_s,calcs;

  int index_thetamin,index_thetamax;
    
  if (e<e_min){
    theta_order=theta/3; //order from 0
    thetamin=theta_order*3;
    index_thetamin=theta_order;
    index_thetamax=theta_order+1;
        
    calcs=(vect[index_thetamax]-vect[index_thetamin])/3*(theta-thetamin)+vect[index_thetamin];
        
    obs_status=0;
    return calcs;
  }
  else if(e>e_max){
            
    theta_order=theta/3; //order from 0
    thetamin=theta_order*3;
    index_thetamin=85*61+theta_order;
    index_thetamax=85*61+theta_order+1;
            
    calcs=(vect[index_thetamax]-vect[index_thetamin])/3*(theta-thetamin)+vect[index_thetamin];
    obs_status=0;
    return calcs;
  }

  else{
    e_order=(e-e_min)/20; //order from 0
    theta_order=theta/3; //order from 0
    emin=e_min+e_order*20;
    thetamin=theta_order*3;
    index_emin_thetamin=e_order*61+theta_order;
    index_emin_thetamax=e_order*61+theta_order+1;
    index_emax_thetamin=(e_order+1)*61+theta_order;
    index_emax_thetamax=(e_order+1)*61+theta_order+1;

    emin_s=(vect[index_emin_thetamax]-vect[index_emin_thetamin])/3*(theta-thetamin)+vect[index_emin_thetamin];
    emax_s=(vect[index_emax_thetamax]-vect[index_emax_thetamin])/3*(theta-thetamin)+vect[index_emax_thetamin];
    calcs = (emax_s-emin_s)/20*(e-emin)+emin_s;
    obs_status=1;

    return calcs;
  }
}
//// Function "calcs_piplusn" to calculate cross section for gamma+p -> piplus+n ////


//// Function "calcdcsKL_original" to calculate dcs for piminus + p ->K + Lambda////
double calcdcsKL_original(double Wv, double ctv, int &Sstatus){
  ifstream infile("csKL_orignial");
  string str,strsub,str1,str2;
  int ch_first,ch_last;
  const int W_binnum=49;
  const int ct_binnum=50;
  double W[W_binnum],W_min=1626,W_max=2405;
  double ct[W_binnum][ct_binnum],dcs[W_binnum][ct_binnum];
  
  ostringstream sstr;
  int Wi=-1,ctj=0;
  while(getline(infile,str)) {
    if (str.empty()) continue;
    else if(str.find("W")!=string::npos){
      ch_first=str.find_last_of("=")+1;
      ch_last=str.length()-1;
      Wi++;
      ctj=0;
      W[Wi]=atof(str.substr(ch_first,ch_last).c_str());
      continue;
    }

    else{
      str1=str.substr(0,str.find_first_of(" "));
      str2=str.substr(str.find_last_of(" ")+1,str.length()-1);
      ct[Wi][ctj]=atof(str1.c_str());
      dcs[Wi][ctj]=atof(str2.c_str());
      ctj++;
      sstr.str("");
      continue;
    }  
  }

  double calc_dcs;
  double W1=0,W2=0,ct11=0,ct12=0,ct21=0,ct22=0,dcs11=0,dcs12=0,dcs21=0,dcs22=0,dcs1=0,dcs2=0;
  int ith;
    
  if (Wv<=W_min){
    if(ctv<=ct[0][0]){
      calc_dcs=dcs[0][0];
    }
    else{
      for(int j=0;j<ct_binnum;j++){
	if(ctv>= ct[0][j] &&  ctv<ct[0][j+1]){
	  calc_dcs=(ctv-ct[0][j])/(ct[0][j+1]-ct[0][j])*(dcs[0][j+1]-dcs[0][j])+dcs[0][j];
	  break;
	}
      }
    }    
    Sstatus=0;
    return calc_dcs;
  }

  else if (Wv>W_max){
    if(ctv<=ct[W_binnum-1][0]){
      calc_dcs=dcs[W_binnum-1][0];
    }
    else{
      for(int j=0;j<ct_binnum;j++){
	if(ctv>= ct[W_binnum-1][j] &&  ctv<ct[W_binnum-1][j+1]){
	  calc_dcs=(ctv-ct[W_binnum-1][j])/(ct[W_binnum-1][j+1]-ct[W_binnum-1][j])*(dcs[W_binnum-1][j+1]-dcs[W_binnum-1][j])+dcs[W_binnum-1][j];
	  break;
	}
      }   
    }
    Sstatus=0;
    return calc_dcs;
  }

  else{
    for(int i=0;i<W_binnum;i++){
      if(Wv>=W[i] && Wv<W[i+1]){
	  W1=W[i];
	  W2=W[i+1];
	  ith=i;
	  if(ctv<=ct[i][0]){
	    dcs1=dcs[i][0];
	    dcs2=dcs[i+1][0];
	    calc_dcs=(Wv-W1)/(W2-W1)*(dcs2-dcs1)+dcs1;
	    break;
	  }
	  else{
	    for(int j=0;j<ct_binnum;j++){
	      if(ctv>ct[ith][j] && ctv<=ct[ith][j+1]){
	      ct11=ct[ith][j];
	      ct12=ct[ith][j+1];
	      dcs11=dcs[ith][j];
	      dcs12=dcs[ith][j+1];
	      break;
	    }
	  }
	    for(int j=0;j<ct_binnum;j++){
	      if(ctv>ct[ith+1][j] && ctv<=ct[ith+1][j+1]){
	      ct21=ct[ith+1][j];
	      ct22=ct[ith+1][j+1];
	      dcs21=dcs[ith+1][j];
	      dcs22=dcs[ith+1][j+1];
	      break;
	    }
	  }
	  dcs1=(ctv-ct11)/(ct12-ct11)*(dcs12-dcs11)+dcs11;
	  dcs2=(ctv-ct21)/(ct22-ct21)*(dcs22-dcs21)+dcs21;
	  calc_dcs=(Wv-W1)/(W2-W1)*(dcs2-dcs1)+dcs1;
	  break;
	}
      }
    }

    Sstatus=1;
    return calc_dcs;
  }
}
//// Function "calcdcsKL_adjust" to calculate dcs for piminus + p ->K + Lambda////

//// Function "readSdcs" to read differential cross section from table of the second step ////
vector<double> readSdcs(string filename){
  ifstream input(filename.c_str());

  vector<double> vect_dcs;
  vector<string> table;
  string line,tempstring,subline,subline1,subline2,subline3,subline4;
  int tindex1,tindex2,tindex3,tindex4;
  double dcs;

  while(getline(input,line)) {
    if (line.empty() || line.find("W")!=string::npos) continue;
    table.push_back(line);

    tindex1 = line.find_first_not_of(" ");
    subline=line.substr(tindex1);
    tindex2 = subline.find_first_of(" ");
    subline1 = subline.substr(0, tindex2);
     
    subline2 = subline.substr(tindex2);
    tindex3 = subline2.find_first_not_of(" ");
    subline3 = subline2.substr(tindex3);
    tindex4 = subline3.find_first_of(" ");
    subline4 = subline3.substr(0, tindex4);
    istringstream(subline4)>>dcs;
    vect_dcs.push_back(dcs);
  }
  return vect_dcs;
}
//// Function "readSdcs" to read differential cross section of table of the second step ////


//// Function "calcdcsKL_adjust" to calculate cross section for piminus + p ->K + Lambda ////
double calcdcsKL_adjust(vector<double> vect, double w, double ct, int &status){
  double w_min=1626, w_max=2405;
  int w_order,ct_order;
  double wmin,ctmin;
  int index_wmin_ctmin,index_wmin_ctmax,index_wmax_ctmin,index_wmax_ctmax;
  const int w_binnum=49;
  const int ct_binnum=50;
  double w_binwidth=(w_max-w_min)/(double)(w_binnum-1);
  double ct_binwidth=2./(double)(ct_binnum-1);
  double wmin_dcs,wmax_dcs,calcdcs;

  int index_ctmin,index_ctmax;
    
  if (w<w_min){
    ct_order=(int)((ct+1)/ct_binwidth); //order from 0
    ctmin=-1+ct_order*ct_binwidth;
    index_ctmin=ct_order;
    index_ctmax=ct_order+1;
        
    calcdcs=(vect[index_ctmax]-vect[index_ctmin])/ct_binwidth*(ct-ctmin)+vect[index_ctmin];
        
    status=0;
    return calcdcs;
  }

  else if(w>w_max){            
    ct_order=(int)((ct+1)/ct_binwidth); //order from 0
    ctmin=-1+ct_order*ct_binwidth;
    index_ctmin=vect.size()-(w_binnum-1)*ct_binnum+ct_order;
    index_ctmax=vect.size()-(w_binnum-1)*ct_binnum+ct_order+1;
            
    calcdcs=(vect[index_ctmax]-vect[index_ctmin])/ct_binwidth*(ct-ctmin)+vect[index_ctmin];
    status=0;
    return calcdcs;
  }

  else{
    w_order=(int)((w-w_min)/w_binwidth); //order from 0
    ct_order=(int)((ct+1)/ct_binwidth); //order from 0
    wmin=w_min+w_order*w_binwidth;
    ctmin=-1+ct_order*ct_binwidth;
    index_wmin_ctmin=w_order*ct_binnum+ct_order;
    index_wmin_ctmax=w_order*ct_binnum+ct_order+1;
    index_wmax_ctmin=(w_order+1)*ct_binnum+ct_order;
    index_wmax_ctmax=(w_order+1)*ct_binnum+ct_order+1;

    wmin_dcs=(vect[index_wmin_ctmax]-vect[index_wmin_ctmin])/ct_binwidth*(ct-ctmin)+vect[index_wmin_ctmin];
    wmax_dcs=(vect[index_wmax_ctmax]-vect[index_wmax_ctmin])/ct_binwidth*(ct-ctmin)+vect[index_wmax_ctmin];
    calcdcs = (wmax_dcs-wmin_dcs)/w_binwidth*(w-wmin)+wmin_dcs;
    status=1;

    return calcdcs;
  }
}
//// Function "calcdcsKL_adjust" to calculate cross section for piminus + p ->K + Lambda ////


//// Function "readSTOxOzPforGPKL" to read Sigma, T, Ox, Oz and P from a linear polarization observable table for gamma+p -> kaon+Lambda ////
vector<vector<double> > readSTOxOzPforGPKL(string filename){
  ifstream infile(filename.c_str());
  double w,costheta,sigma,t,ox,oz,p;
  vector<double> vectobs;
  vector<vector<double> > vectreturn;
  vectobs.clear();
  vectreturn.clear();
  for(int i=0;i<192;i++){
    infile>>w>>costheta>>sigma>>t>>ox>>oz>>p;
    vectobs.push_back(sigma);
    vectobs.push_back(t);
    vectobs.push_back(ox);
    vectobs.push_back(oz);
    vectobs.push_back(p);
    vectreturn.push_back(vectobs);
    vectobs.clear();
  }
  infile.close();

  return vectreturn;
}
//// Function "readSTOxOzPforGPKL" to read Sigma, T, Ox, Oz and P from a linear polarization observable table for gamma+p -> kaon+Lambda ////

//// Function "calcSTOxOzPforGPKL" to calculate Sigma, T, Ox, Oz and P for gamma+p -> kaon+Lambda
vector<double> calcSTOxOzPforGPKL(double w, double ct, int &obs_status, vector<vector<double> > vectObs){

  vector<double> vectSigma, vectT, vectOx, vectOz, vectP;
  for(int i=0;i<vectObs.size();i++){
    vectSigma.push_back(vectObs[i][0]);
    vectT.push_back(vectObs[i][1]);
    vectOx.push_back(vectObs[i][2]);
    vectOz.push_back(vectObs[i][3]);
    vectP.push_back(vectObs[i][4]);
  }

  double w_min=1720, w_max=2180;
  double ct_min=-0.65, ct_max=0.75;
  int w_order,ct_order;
  double wmin,ctmin;
  int index_wmin_ctmin,index_wmin_ctmax,index_wmax_ctmin,index_wmax_ctmax;
  const int w_binnum=24;
  const int ct_binnum=8;
  double w_binwidth=(w_max-w_min)/(double)(w_binnum-1);
  double ct_binwidth=(ct_max-ct_min)/(double)(ct_binnum-1);
  double wminSigma,wmaxSigma;
  double wminT,wmaxT;
  double wminOx,wmaxOx;
  double wminOz,wmaxOz;
  double wminP,wmaxP;
  double calcSigma,calcT,calcOx,calcOz,calcP;

  int index_ctmin,index_ctmax;
    
  if (w<w_min){
    if(ct<ct_min){
      calcSigma=vectSigma[0];
      calcT=vectT[0];
      calcOx=vectOx[0];
      calcOz=vectOz[0];
      calcP=vectP[0];
    }
    else if(ct>ct_max){
      calcSigma=vectSigma[ct_binnum-1];
      calcT=vectT[ct_binnum-1];
      calcOx=vectOx[ct_binnum-1];
      calcOz=vectOz[ct_binnum-1];
      calcP=vectP[ct_binnum-1];
    }

    else{
      ct_order=(int)((ct-ct_min)/ct_binwidth); //order from 0
      ctmin=ct_min+ct_order*ct_binwidth;
      index_ctmin=ct_order;
      index_ctmax=ct_order+1;   
      calcSigma=(vectSigma[index_ctmax]-vectSigma[index_ctmin])/ct_binwidth*(ct-ctmin)+vectSigma[index_ctmin];
      calcT=(vectT[index_ctmax]-vectT[index_ctmin])/ct_binwidth*(ct-ctmin)+vectT[index_ctmin];
      calcOx=(vectOx[index_ctmax]-vectOx[index_ctmin])/ct_binwidth*(ct-ctmin)+vectOx[index_ctmin];
      calcOz=(vectOz[index_ctmax]-vectOz[index_ctmin])/ct_binwidth*(ct-ctmin)+vectOz[index_ctmin];
      calcP=(vectP[index_ctmax]-vectP[index_ctmin])/ct_binwidth*(ct-ctmin)+vectP[index_ctmin];
    }
        
    obs_status=0;
    vector<double> vect;
    vect.push_back(calcSigma);
    vect.push_back(calcT);
    vect.push_back(calcOx);
    vect.push_back(calcOz);
    vect.push_back(calcP);
    return vect;
  }

  else if(w>w_max){            
    if(ct<ct_min){
      calcSigma=vectSigma[(w_binnum-1)*ct_binnum];
      calcT=vectT[(w_binnum-1)*ct_binnum];
      calcOx=vectOx[(w_binnum-1)*ct_binnum];
      calcOz=vectOz[(w_binnum-1)*ct_binnum];
      calcP=vectP[(w_binnum-1)*ct_binnum];
    }
    else if(ct>ct_max){
      calcSigma=vectSigma[w_binnum*ct_binnum-1];
      calcT=vectT[w_binnum*ct_binnum-1];
      calcOx=vectOx[w_binnum*ct_binnum-1];
      calcOz=vectOz[w_binnum*ct_binnum-1];
      calcP=vectP[w_binnum*ct_binnum-1];
    }

    else{
      ct_order=(int)((ct-ct_min)/ct_binwidth); //order from 0
      ctmin=ct_min+ct_order*ct_binwidth;
      index_ctmin=(w_binnum-1)*ct_binnum+ct_order;
      index_ctmax=(w_binnum-1)*ct_binnum+ct_order+1;   
      calcSigma=(vectSigma[index_ctmax]-vectSigma[index_ctmin])/ct_binwidth*(ct-ctmin)+vectSigma[index_ctmin];
      calcT=(vectT[index_ctmax]-vectT[index_ctmin])/ct_binwidth*(ct-ctmin)+vectT[index_ctmin];
      calcOx=(vectOx[index_ctmax]-vectOx[index_ctmin])/ct_binwidth*(ct-ctmin)+vectOx[index_ctmin];
      calcOz=(vectOz[index_ctmax]-vectOz[index_ctmin])/ct_binwidth*(ct-ctmin)+vectOz[index_ctmin];
      calcP=(vectP[index_ctmax]-vectP[index_ctmin])/ct_binwidth*(ct-ctmin)+vectP[index_ctmin];
    }
        
    obs_status=0;
    vector<double> vect;
    vect.push_back(calcSigma);
    vect.push_back(calcT);
    vect.push_back(calcOx);
    vect.push_back(calcOz);
    vect.push_back(calcP);
    return vect;
  }

  else{
    if(ct<ct_min){
      w_order=(int)((w-w_min)/w_binwidth); //order from 0
      wmin=w_min+w_binwidth*w_order;
      index_wmin_ctmin=w_order*ct_binnum;
      index_wmax_ctmin=(w_order+1)*ct_binnum;
      calcSigma=(vectSigma[index_wmax_ctmin]-vectSigma[index_wmin_ctmin])/w_binwidth*(w-wmin)+vectSigma[index_wmin_ctmin];
      calcT=(vectT[index_wmax_ctmin]-vectT[index_wmin_ctmin])/w_binwidth*(w-wmin)+vectT[index_wmin_ctmin];
      calcOx=(vectOx[index_wmax_ctmin]-vectOx[index_wmin_ctmin])/w_binwidth*(w-wmin)+vectOx[index_wmin_ctmin];
      calcOz=(vectOz[index_wmax_ctmin]-vectOz[index_wmin_ctmin])/w_binwidth*(w-wmin)+vectOz[index_wmin_ctmin];
      calcP=(vectP[index_wmax_ctmin]-vectP[index_wmin_ctmin])/w_binwidth*(w-wmin)+vectP[index_wmin_ctmin];
    }

    else if(ct>ct_max){
      w_order=(int)((w-w_min)/w_binwidth); //order from 0
      wmin=w_min+w_binwidth*w_order;
      index_wmin_ctmax=w_order*ct_binnum+ct_binnum-1;
      index_wmax_ctmax=(w_order+1)*ct_binnum+ct_binnum-1;
      calcSigma=(vectSigma[index_wmax_ctmax]-vectSigma[index_wmin_ctmax])/w_binwidth*(w-wmin)+vectSigma[index_wmin_ctmax];
      calcT=(vectT[index_wmax_ctmax]-vectT[index_wmin_ctmax])/w_binwidth*(w-wmin)+vectT[index_wmin_ctmax];
      calcOx=(vectOx[index_wmax_ctmax]-vectOx[index_wmin_ctmax])/w_binwidth*(w-wmin)+vectOx[index_wmin_ctmax];
      calcOz=(vectOz[index_wmax_ctmax]-vectOz[index_wmin_ctmax])/w_binwidth*(w-wmin)+vectOz[index_wmin_ctmax];
      calcP=(vectP[index_wmax_ctmax]-vectP[index_wmin_ctmax])/w_binwidth*(w-wmin)+vectP[index_wmin_ctmax];
    }

    else{
      w_order=(int)((w-w_min)/w_binwidth); //order from 0
      ct_order=(int)((ct-ct_min)/ct_binwidth); //order from 0
      wmin=w_min+w_order*w_binwidth;
      ctmin=ct_min+ct_order*ct_binwidth;
      index_wmin_ctmin=w_order*ct_binnum+ct_order;
      index_wmin_ctmax=w_order*ct_binnum+ct_order+1;
      index_wmax_ctmin=(w_order+1)*ct_binnum+ct_order;
      index_wmax_ctmax=(w_order+1)*ct_binnum+ct_order+1;

      wminSigma=(vectSigma[index_wmin_ctmax]-vectSigma[index_wmin_ctmin])/ct_binwidth*(ct-ctmin)+vectSigma[index_wmin_ctmin];
      wmaxSigma=(vectSigma[index_wmax_ctmax]-vectSigma[index_wmax_ctmin])/ct_binwidth*(ct-ctmin)+vectSigma[index_wmax_ctmin];
      calcSigma = (wmaxSigma-wminSigma)/w_binwidth*(w-wmin)+wminSigma;

      wminT=(vectT[index_wmin_ctmax]-vectT[index_wmin_ctmin])/ct_binwidth*(ct-ctmin)+vectT[index_wmin_ctmin];
      wmaxT=(vectT[index_wmax_ctmax]-vectT[index_wmax_ctmin])/ct_binwidth*(ct-ctmin)+vectT[index_wmax_ctmin];
      calcT = (wmaxT-wminT)/w_binwidth*(w-wmin)+wminT;

      wminOx=(vectOx[index_wmin_ctmax]-vectOx[index_wmin_ctmin])/ct_binwidth*(ct-ctmin)+vectOx[index_wmin_ctmin];
      wmaxOx=(vectOx[index_wmax_ctmax]-vectOx[index_wmax_ctmin])/ct_binwidth*(ct-ctmin)+vectOx[index_wmax_ctmin];
      calcOx = (wmaxOx-wminOx)/w_binwidth*(w-wmin)+wminOx;

      wminOz=(vectOz[index_wmin_ctmax]-vectOz[index_wmin_ctmin])/ct_binwidth*(ct-ctmin)+vectOz[index_wmin_ctmin];
      wmaxOz=(vectOz[index_wmax_ctmax]-vectOz[index_wmax_ctmin])/ct_binwidth*(ct-ctmin)+vectOz[index_wmax_ctmin];
      calcOz = (wmaxOz-wminOz)/w_binwidth*(w-wmin)+wminOz;

      wminP=(vectP[index_wmin_ctmax]-vectP[index_wmin_ctmin])/ct_binwidth*(ct-ctmin)+vectP[index_wmin_ctmin];
      wmaxP=(vectP[index_wmax_ctmax]-vectP[index_wmax_ctmin])/ct_binwidth*(ct-ctmin)+vectP[index_wmax_ctmin];
      calcP = (wmaxP-wminP)/w_binwidth*(w-wmin)+wminP;
    }

    obs_status=1;
    vector<double> vect;
    vect.push_back(calcSigma);
    vect.push_back(calcT);
    vect.push_back(calcOx);
    vect.push_back(calcOz);
    vect.push_back(calcP);
    return vect;
  }
}
//// Function "calcSTOxOzPforGPKL" to calculate Sigma, T, Ox, Oz and P for gamma+p -> kaon+Lambda



/* ----------- Calculate differential cross section ---------------- */

void PrintUsage (char *processName){
  cout << "\nUsage: " << processName << " [options] inputFile\n\n";
  cout << "\toptions:\n";
  cout << "\t-R <filename>\tdirect ROOT output to <filename>\n";

  cout << "\t-M[#]\t\tProcess  only (#) number of events\n";
  cout << "\t-h\t\tThis information\n\n";
  cout << "\t-n\t\tHow many processed channels? Default: 9\n";
  cout << "\t-c\t\tProcessed channel number. Must be followed after option 'n'. Default: 123456789\n";
  cout << "\tFor example: -n -c 2457 means 4 channels 2, 4, 5 and 7 will be processed.\n\n";


  //JLA  explain usage of -E
  cout << "\t-E <expression>\tPhoton energy\n";
  cout << "\texpression (no whitespace): \n";
  cout << "\t[x]           \tMonoenergetic E = [x] GeV\n";
  cout << "\t[x]:[y]       \tPlain distribution from [x] to [y] GeV\n";
  cout << "\tmono:[x]      \tMonoenergetic as above\n";
  cout << "\tplain:[x]:[y] \tPlain distribution as above\n";
  cout << "\tbrems:[x]:[y] \tBremsstrahlung distribution 1/E\n";
  cout << "\thisto:<filename>:<hstname>\tDistribution according\n";
  cout << "\t              \to one-dim. ROOT-histogram <hstname> in <filename>\n" << endl;

  cout << "\t-p\tDetermine if implementing unpolarized(u) or circularly-polarized(c) or linearly-polarizied(l)  differential-cross-section(cs) for the first step. For all channels, ncs means not to implement differential cross section. For channels 2 and 8, ucs means implementing unpolarized DCS, ccs means implementing circularly polarized DCS and lcs means implementing linearly polarized DCS. For channels 1 and 5, ucs means implementing unpolarized DCS when argument is ucs or ccs since no circular ploarization observable tables can be applied for these two channels so that circularly polarized DCS can not be implemented, and lcs means implementing linearly polarized DCS. For channels 3, 4 and 9, unpolarized DCS is implemented when argument is ucs, ccs or lcs since no ploarization observable tables are applied for these three channels so that polarized DCS are not implemented. For channes 6 and 7, no corss section is implemented when argument is ncs, ucs, ccs or lcs since no unpolarized DCS tables can be applied for these two channels so that unpolarized and polarized DCS can not be implemented.\n";
  cout << "\t ncs (default)\tDon't implement differential cross section\n";
  cout << "\t ucs \tImplement unpolarized differential cross section\n";
  cout << "\t ccs\tImplement circularly polarized differential cross section\n"<<endl;
  cout << "\t lcs\tImplement linearly polarized differential cross section\n"<<endl;;

  cout << "\t-q\tDetermine if implement unpolarized cross section for the second step of channels 1 and 5\n";
  cout << "\t ncs (default)\tDon't implement differential cross section\n";
  cout << "\t ocs \tImplement unpolarized differential cross section using oringinal differential cross section table (longer running time)\n";
  cout << "\t acs\tImplement unpolarized differential cross section using adjusted differential cross section table (shorter running time)\n"<<endl;

  cout << "\t-s\tDetermine if the maximum of random values for unpolarized differential cross section of the first step of channels 1, 2, 3, 4, 5, 8 and 9 is different or not, and also determine if the other maximum of random values for linearly polarized differential cross section of the first step of channels 1, 2, 5 and 8 is different or not.\n";
  cout << "\t diff (default)\tDifferent channels have different maximum values dependent of their own different cross section table and observable table (shorter running time)\n";
  cout << "\t same\tAll channels have the same maximum value, which is the maximum value of all tables for the case of unpolarized differential cross section and is the other maximum value of linearly polairzed differential cross section. (much longer running time)\n"<<endl;

  exit (0);
}
