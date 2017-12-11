using namespace std;
#include <unistd.h>
#include <string>

extern "C"
{

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <signal.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <ntypes.h>
#include <bostypes.h>
#include <clas_cern.h>
#include <scalers.h>
#include <utility.h>
#include <printBOS.h>
#include <time.h>
  //Paul//
#include <bosfun.h>
#include <bosio.h>
#include <biofun.h>
  //    //
#define speed_light 29979245800 // cm/s
#define lifetime_Lambda 0.00000000026 // s
#define lifetime_SigmaStar 0.000000000000000000000018 // s
#define lifetime_SigmaStarMinus 0.0000000000000000000000167 // s

#define mass_protPDG 0.93827201
#define mass_neutPDG 0.939565378
#define mass_piplusPDG 0.13957018
#define mass_piminusPDG 0.1394557018
#define mass_pizeroPDG 0.1349766
#define mass_kaonPDG 0.493677
#define mass_lambdaPDG 1.115683
#define mass_photonPDG 0.0
#define mass_sigmaPDG 1.192642
#define mass_sigmastarPDG 1.3837
#define mass_sigmastarminusPDG 1.3872

    
    
#define prot_numPDG 2212
#define neut_numPDG 2112
#define piplus_numPDG 211
#define piminus_numPDG -211
#define pizero_numPDG 111
#define kaon_numPDG 321
#define lambda_numPDG 3122
#define photon_numPDG 22
#define sigma_numPDG 3212
#define sigmastar_numPDG 3214
#define sigmastarminus_numPDG 3114
    
  BOSbank bcs_;
}

#include <iostream>
#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TNtuple.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <TF1.h>
#include <TSystem.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TStyle.h>
#include <sstream>

void Display_Help();

main (int argc, char *argv[]){
  ////////////////////////////////////////////     ////////////////////////////////////////////
  cout << "Author: Tongtong Cao (caot@jlab.org)" << endl;
  cout << "Author of origninal version: Nicholas Zachariou" << endl;
  cout << "Run with '-H' or '-h' switch for help." << endl;
  int loc_i, locMaxNumEventsPerFile = -1, locNumEvents, loc_numbroot, numbroot=0;
  vector<string> locStringInputs; //vectors (can change their size and allocate by push_back command)
  string locTempString;
  istringstream locIStream;
  for(loc_i = 1; loc_i < argc; loc_i++){
    if(argv[loc_i][0] != '-'){
      locStringInputs.push_back(argv[loc_i]);
      numbroot++;
    }
    else{	//it's a switch:
      switch(argv[loc_i][1]){
      case 'h'://display help
	Display_Help();
	return 0;
      case 'H'://display help
	Display_Help();
	return 0;
      case 'M'://set MaxNumEventsPerFile
	locTempString = argv[loc_i];
	locTempString = locTempString.substr(2, locTempString.length() - 2);	//strip "-M"
	locIStream.str(locTempString);	//stores locTempString to locIStreat
	if(!(locIStream >> locMaxNumEventsPerFile))	//converts locIStream into integer and assigns it to LocMaxNumEventsPerFile
	  cout << "ERROR: COMMAND LINE INPUT NOT RECOGNIZED.  MaxNumEventsPerFile NOT SET." << endl;
	break;
      default:
	break;
      }//kills switch
    }
  }//kills loop over arguments
  if(locStringInputs.size() == 0){	//Checks if there is an input bos file
    cout << "ERROR: NEED INPUT ROOT FILE." << endl;
    return 2;
  }
  /////////////////////////////////////////////////////////////////////////////

    
  ////// Load root files in a TChain for running through all events//////
  TChain *tree = new TChain("mytree");
  for (loc_numbroot=0; loc_numbroot<numbroot; loc_numbroot++){
    cout <<" Additing to TChain, Root file: "<<locStringInputs[loc_numbroot].c_str()<<endl;
    tree->Add(locStringInputs[loc_numbroot].c_str());
  }
  ///////////////////////////////////////////////////////////////////////

    
  ////////////// Set up BOS file -- Open and initialize/////////////////
  BOSbank bThreadBOSCommonBlock_Output;
  int bThreadBOSIOptr_Output;
  //initialize BOS common block
  bosInit(bThreadBOSCommonBlock_Output.iw, 700000); //if more than one thread accessing the same file, should have one per thread (as well as locks around BOS access!)
    
  //open the file
  //string locOutputFileName = "rename.bos";
  string locOutputFileName = locStringInputs[0].substr(0, locStringInputs[0].length()-4);
  locOutputFileName.append("evt");
    
  char* locCommand_Char = new char[locOutputFileName.size() + 1];
  strncpy(locCommand_Char, locOutputFileName.c_str(), locOutputFileName.size() + 1);
  bool locOutputFileOpenFlag = !bosOpen(locCommand_Char, "w", &bThreadBOSIOptr_Output);
  delete locCommand_Char;

  /////////////////// Set up Root Branches to read //////////////////////
  int rlabel,rcirpol,rlinpol;
  TLorentzVector *rtarget=0;
  TLorentzVector *rbeam=0;
  TLorentzVector *rQp1=0;
  TLorentzVector *rQp2=0;
  TLorentzVector *rStarget=0;
  TLorentzVector *rSbeam=0;
  TLorentzVector *rRp1=0;
  TLorentzVector *rRp2=0;
  TLorentzVector *rDp1=0;
  TLorentzVector *rDp2=0;
    
  tree->SetBranchAddress("target",&rtarget);
  tree->SetBranchAddress("beam",&rbeam);
  tree->SetBranchAddress("Qp1",&rQp1);
  tree->SetBranchAddress("Qp2",&rQp2);
  tree->SetBranchAddress("Starget",&rStarget);
  tree->SetBranchAddress("Sbeam",&rSbeam);
  tree->SetBranchAddress("Rp1",&rRp1);
  tree->SetBranchAddress("Rp2",&rRp2);
  tree->SetBranchAddress("Dp1",&rDp1);
  tree->SetBranchAddress("Dp2",&rDp2);
  tree->SetBranchAddress("label",&rlabel);
  tree->SetBranchAddress("cirpol",&rcirpol);
  tree->SetBranchAddress("linpol",&rlinpol);
    
  float vertex_x, vertex_y, vertex_z;
  TF1 *myfuncx = new TF1("myfuncx","TMath::BreitWigner(x, 0.1398, 0.366507 )",-2,2); //Breit-Wigner to simulate x-vertex from
  TF1 *myfuncy = new TF1("myfuncy","TMath::BreitWigner(x, -0.0989511, 0.408273 )",-2,2); //Breit-Wigner to simulate y-vertex from
  TRandom *r3 = new TRandom3(0); // In this case the seed is guaranteed to be unique in space and time
  gRandom=r3;
    
  locNumEvents=0;
  int loc_NumEntries = tree->GetEntries();
  if (locMaxNumEventsPerFile>0)
    cout << "Main file has " << loc_NumEntries<<" entries. Will process "<<locMaxNumEventsPerFile<<" events" << endl;
  else 
    cout << "Main file has " << loc_NumEntries<<" entries. "<< endl;
    
  for(loc_i = 0; loc_i < loc_NumEntries; loc_i++){
    if((locMaxNumEventsPerFile <= locNumEvents) && (locMaxNumEventsPerFile > 0))//past max
      break;
    tree->GetEvent(loc_i);
    if (loc_i%10000==0)
      cout<<loc_i<<"  events processed"<<endl;
    vertex_x=myfuncx->GetRandom();
    vertex_y=myfuncy->GetRandom();
    vertex_z=r3->Uniform(-40, 0);
        
    int locNumRowsHead=1, locNumRowsMCVX=2, locNumRowsMCTK=10;
    //for each event: create a BOS bank
    string locBankList;
    locBankList.clear();
        
    clasHEAD_t* locHEADBank = (clasHEAD_t*)makeBank(&bThreadBOSCommonBlock_Output, "HEAD", 0, sizeof(head_t)/sizeof(int), locNumRowsHead);
    locBankList += "HEAD";
    clasMCTK_t* locMCTKBank = (clasMCTK_t*)makeBank(&bThreadBOSCommonBlock_Output, "MCTK", 0, sizeof(mctk_t)/sizeof(int), locNumRowsMCTK);
    locBankList += "MCTK";
    clasMCVX_t* locMCVXBank = (clasMCVX_t*)makeBank(&bThreadBOSCommonBlock_Output, "MCVX", 0, sizeof(mcvx_t)/sizeof(int), locNumRowsMCVX);
    locBankList += "MCVX";
        
    bosNformat(bThreadBOSCommonBlock_Output.iw, (char*)"HEAD", "8I");
    bosNformat(bThreadBOSCommonBlock_Output.iw, (char*)"MCTK","6F,5I");
    bosNformat(bThreadBOSCommonBlock_Output.iw, (char*)"MCVX","4F,I");
        
    locHEADBank->head[0].version = 1;
    locHEADBank->head[0].nrun = 53600;
    locHEADBank->head[0].nevent = locNumEvents++;
    locHEADBank->head[0].type = -4;  //Negative for simulations
    locHEADBank->head[0].evtclass = 15; //0-15 for physics events
    locHEADBank->head[0].trigbits = 0;
    locHEADBank->head[0].time =0;

    //// vertex is smeared with ffread card for beam position and sigma
    // Reaction vertex
    locMCVXBank->mcvx[0].x = vertex_x;
    locMCVXBank->mcvx[0].y = vertex_y;
    locMCVXBank->mcvx[0].z = vertex_z;
    locMCVXBank->mcvx[0].tof = 0.0;
    locMCVXBank->mcvx[0].flag = 0;
      
    //-----------FIRST TRACK is Photon -------------//
    locMCTKBank->mctk[0].cx = rbeam->Px()/rbeam->Rho();
    locMCTKBank->mctk[0].cy = rbeam->Py()/rbeam->Rho();;
    locMCTKBank->mctk[0].cz = rbeam->Pz()/rbeam->Rho();;
    locMCTKBank->mctk[0].pmom = rbeam->Rho();
    locMCTKBank->mctk[0].mass = mass_photonPDG;
    locMCTKBank->mctk[0].charge = 0;
    locMCTKBank->mctk[0].id = photon_numPDG ;
    locMCTKBank->mctk[0].beg_vtx = 0;
    locMCTKBank->mctk[0].end_vtx = 1;
    locMCTKBank->mctk[0].flag = 1;
    locMCTKBank->mctk[0].parent = 0;
    //--------------------------------------------//
        
    if (rlabel == 1){
      //Decay vertex
      locMCVXBank->mcvx[1].x = vertex_x+r3->Exp(lifetime_Lambda)*rRp1->Px()/mass_lambdaPDG*speed_light;
      locMCVXBank->mcvx[1].y = vertex_y+r3->Exp(lifetime_Lambda)*rRp1->Py()/mass_lambdaPDG*speed_light;
      locMCVXBank->mcvx[1].z = vertex_z+r3->Exp(lifetime_Lambda)*rRp1->Pz()/mass_lambdaPDG*speed_light;
      locMCVXBank->mcvx[1].tof = 0.0;
      locMCVXBank->mcvx[1].flag = 0;

      locMCTKBank->mctk[1].cx = rtarget->Px()/rtarget->Rho();
      locMCTKBank->mctk[1].cy = rtarget->Py()/rtarget->Rho();
      locMCTKBank->mctk[1].cz = rtarget->Pz()/rtarget->Rho();
      locMCTKBank->mctk[1].pmom = rtarget->Rho();
      locMCTKBank->mctk[1].mass = rtarget->M();
      locMCTKBank->mctk[1].charge = 0;
      locMCTKBank->mctk[1].id = neut_numPDG;
      locMCTKBank->mctk[1].beg_vtx = 0;
      locMCTKBank->mctk[1].end_vtx = 1;
      locMCTKBank->mctk[1].flag = rlabel;
      locMCTKBank->mctk[1].parent = 0;
        
      locMCTKBank->mctk[2].cx = rQp1->Px()/rQp1->Rho();
      locMCTKBank->mctk[2].cy = rQp1->Py()/rQp1->Rho();
      locMCTKBank->mctk[2].cz = rQp1->Pz()/rQp1->Rho();
      locMCTKBank->mctk[2].pmom = rQp1->Rho();
      locMCTKBank->mctk[2].mass = mass_pizeroPDG;
      locMCTKBank->mctk[2].charge = 0;
      locMCTKBank->mctk[2].id = pizero_numPDG;
      locMCTKBank->mctk[2].beg_vtx = 0;
      locMCTKBank->mctk[2].end_vtx = 1;
      locMCTKBank->mctk[2].flag = rlabel;
      locMCTKBank->mctk[2].parent = 0;
        
      locMCTKBank->mctk[3].cx = rQp2->Px()/rQp2->Rho();
      locMCTKBank->mctk[3].cy = rQp2->Py()/rQp2->Rho();
      locMCTKBank->mctk[3].cz = rQp2->Pz()/rQp2->Rho();
      locMCTKBank->mctk[3].pmom = rQp2->Rho();
      locMCTKBank->mctk[3].mass = mass_neutPDG;
      locMCTKBank->mctk[3].charge = 0;
      locMCTKBank->mctk[3].id = neut_numPDG;
      locMCTKBank->mctk[3].beg_vtx = 1;
      locMCTKBank->mctk[3].end_vtx = 0;
      locMCTKBank->mctk[3].flag = rlabel;
      locMCTKBank->mctk[3].parent = 0;
        
      locMCTKBank->mctk[4].cx = rSbeam->Px()/rSbeam->Rho();
      locMCTKBank->mctk[4].cy = rSbeam->Py()/rSbeam->Rho();
      locMCTKBank->mctk[4].cz = rSbeam->Pz()/rSbeam->Rho();
      locMCTKBank->mctk[4].pmom = rSbeam->Rho();
      locMCTKBank->mctk[4].mass = mass_pizeroPDG;
      locMCTKBank->mctk[4].charge = 0;
      locMCTKBank->mctk[4].id = pizero_numPDG;
      locMCTKBank->mctk[4].beg_vtx = 0;
      locMCTKBank->mctk[4].end_vtx = 1;
      locMCTKBank->mctk[4].flag = rlinpol*10+rcirpol;
      locMCTKBank->mctk[4].parent = 0;
        
      locMCTKBank->mctk[5].cx = rStarget->Px()/rStarget->Rho();
      locMCTKBank->mctk[5].cy = rStarget->Py()/rStarget->Rho();
      locMCTKBank->mctk[5].cz = rStarget->Pz()/rStarget->Rho();
      locMCTKBank->mctk[5].pmom = rStarget->Rho();
      locMCTKBank->mctk[5].mass = rStarget->M();
      locMCTKBank->mctk[5].charge = 1;
      locMCTKBank->mctk[5].id = prot_numPDG;
      locMCTKBank->mctk[5].beg_vtx = 0;
      locMCTKBank->mctk[5].end_vtx = 1;
      locMCTKBank->mctk[5].flag = 1;
      locMCTKBank->mctk[5].parent = 0;
        
      locMCTKBank->mctk[6].cx = rRp1->Px()/rRp1->Rho();
      locMCTKBank->mctk[6].cy = rRp1->Py()/rRp1->Rho();
      locMCTKBank->mctk[6].cz = rRp1->Pz()/rRp1->Rho();
      locMCTKBank->mctk[6].pmom = rRp1->Rho();
      locMCTKBank->mctk[6].mass = mass_lambdaPDG;
      locMCTKBank->mctk[6].charge = 0;
      locMCTKBank->mctk[6].id = lambda_numPDG;
      locMCTKBank->mctk[6].beg_vtx = 1;
      locMCTKBank->mctk[6].end_vtx = 0;
      locMCTKBank->mctk[6].flag = 1;
      locMCTKBank->mctk[6].parent = 0;
            
      locMCTKBank->mctk[7].cx = rRp2->Px()/rRp2->Rho();
      locMCTKBank->mctk[7].cy = rRp2->Py()/rRp2->Rho();
      locMCTKBank->mctk[7].cz = rRp2->Pz()/rRp2->Rho();
      locMCTKBank->mctk[7].pmom = rRp2->Rho();
      locMCTKBank->mctk[7].mass = mass_kaonPDG;
      locMCTKBank->mctk[7].charge = 1;
      locMCTKBank->mctk[7].id = kaon_numPDG;
      locMCTKBank->mctk[7].beg_vtx = 1;
      locMCTKBank->mctk[7].end_vtx = 0;
      locMCTKBank->mctk[7].flag = 1;
      locMCTKBank->mctk[7].parent = 0;

      locMCTKBank->mctk[8].cx = rDp1->Px()/rDp1->Rho();
      locMCTKBank->mctk[8].cy = rDp1->Py()/rDp1->Rho();
      locMCTKBank->mctk[8].cz = rDp1->Pz()/rDp1->Rho();
      locMCTKBank->mctk[8].pmom = rDp1->Rho();
      locMCTKBank->mctk[8].mass = mass_protPDG;
      locMCTKBank->mctk[8].charge = 1;
      locMCTKBank->mctk[8].id = prot_numPDG;
      locMCTKBank->mctk[8].beg_vtx = 2;
      locMCTKBank->mctk[8].end_vtx = 0;
      locMCTKBank->mctk[8].flag = 1;
      locMCTKBank->mctk[8].parent = 0;

      locMCTKBank->mctk[9].cx = rDp2->Px()/rDp2->Rho();
      locMCTKBank->mctk[9].cy = rDp2->Py()/rDp2->Rho();
      locMCTKBank->mctk[9].cz = rDp2->Pz()/rDp2->Rho();
      locMCTKBank->mctk[9].pmom = rDp2->Rho();
      locMCTKBank->mctk[9].mass = mass_piminusPDG;
      locMCTKBank->mctk[9].charge = -1;
      locMCTKBank->mctk[9].id = piminus_numPDG ;
      locMCTKBank->mctk[9].beg_vtx = 2;
      locMCTKBank->mctk[9].end_vtx = 0;
      locMCTKBank->mctk[9].flag = 1;
      locMCTKBank->mctk[9].parent = 0;
    }
    else if (rlabel == 2 ){
      //Decay vertex
      locMCVXBank->mcvx[1].x = vertex_x+r3->Exp(lifetime_Lambda)*rQp1->Px()/mass_lambdaPDG*speed_light;
      locMCVXBank->mcvx[1].y = vertex_y+r3->Exp(lifetime_Lambda)*rQp1->Py()/mass_lambdaPDG*speed_light;
      locMCVXBank->mcvx[1].z = vertex_z+r3->Exp(lifetime_Lambda)*rQp1->Pz()/mass_lambdaPDG*speed_light;
      locMCVXBank->mcvx[1].tof = 0.0;
      locMCVXBank->mcvx[1].flag = 0;

      locMCTKBank->mctk[1].cx = rtarget->Px()/rtarget->Rho();
      locMCTKBank->mctk[1].cy = rtarget->Py()/rtarget->Rho();
      locMCTKBank->mctk[1].cz = rtarget->Pz()/rtarget->Rho();
      locMCTKBank->mctk[1].pmom = rtarget->Rho();
      locMCTKBank->mctk[1].mass = rtarget->M();
      locMCTKBank->mctk[1].charge = 1;
      locMCTKBank->mctk[1].id = prot_numPDG;
      locMCTKBank->mctk[1].beg_vtx = 0;
      locMCTKBank->mctk[1].end_vtx = 1;
      locMCTKBank->mctk[1].flag = rlabel;
      locMCTKBank->mctk[1].parent = 0;
            
      locMCTKBank->mctk[2].cx = rQp1->Px()/rQp1->Rho();
      locMCTKBank->mctk[2].cy = rQp1->Py()/rQp1->Rho();
      locMCTKBank->mctk[2].cz = rQp1->Pz()/rQp1->Rho();
      locMCTKBank->mctk[2].pmom = rQp1->Rho();
      locMCTKBank->mctk[2].mass = mass_lambdaPDG;
      locMCTKBank->mctk[2].charge = 0;
      locMCTKBank->mctk[2].id = lambda_numPDG;
      locMCTKBank->mctk[2].beg_vtx = 0;
      locMCTKBank->mctk[2].end_vtx = 2;
      locMCTKBank->mctk[2].flag = rlabel;
      locMCTKBank->mctk[2].parent = 0;
            
      locMCTKBank->mctk[3].cx = rQp2->Px()/rQp2->Rho();
      locMCTKBank->mctk[3].cy = rQp2->Py()/rQp2->Rho();
      locMCTKBank->mctk[3].cz = rQp2->Pz()/rQp2->Rho();
      locMCTKBank->mctk[3].pmom = rQp2->Rho();
      locMCTKBank->mctk[3].mass = mass_kaonPDG;
      locMCTKBank->mctk[3].charge = 1;
      locMCTKBank->mctk[3].id = kaon_numPDG;
      locMCTKBank->mctk[3].beg_vtx = 0;
      locMCTKBank->mctk[3].end_vtx = 1;
      locMCTKBank->mctk[3].flag = rlabel;
      locMCTKBank->mctk[3].parent = 0;
            
      locMCTKBank->mctk[4].cx = rSbeam->Px()/rSbeam->Rho();
      locMCTKBank->mctk[4].cy = rSbeam->Py()/rSbeam->Rho();
      locMCTKBank->mctk[4].cz = rSbeam->Pz()/rSbeam->Rho();
      locMCTKBank->mctk[4].pmom = rSbeam->Rho();
      locMCTKBank->mctk[4].mass = mass_kaonPDG;
      locMCTKBank->mctk[4].charge = 1;
      locMCTKBank->mctk[4].id = kaon_numPDG;
      locMCTKBank->mctk[4].beg_vtx = 0;
      locMCTKBank->mctk[4].end_vtx = 1;
      locMCTKBank->mctk[4].flag = rlinpol*10+rcirpol;
      locMCTKBank->mctk[4].parent = 0;
            
      locMCTKBank->mctk[5].cx = rStarget->Px()/rStarget->Rho();
      locMCTKBank->mctk[5].cy = rStarget->Py()/rStarget->Rho();
      locMCTKBank->mctk[5].cz = rStarget->Pz()/rStarget->Rho();
      locMCTKBank->mctk[5].pmom = rStarget->Rho();
      locMCTKBank->mctk[5].mass = rStarget->M();
      locMCTKBank->mctk[5].charge = 0;
      locMCTKBank->mctk[5].id = neut_numPDG;
      locMCTKBank->mctk[5].beg_vtx = 0;
      locMCTKBank->mctk[5].end_vtx = 1;
      locMCTKBank->mctk[5].flag = 1;
      locMCTKBank->mctk[5].parent = 0;
        
      locMCTKBank->mctk[6].cx = rRp1->Px()/rRp1->Rho();
      locMCTKBank->mctk[6].cy = rRp1->Py()/rRp1->Rho();
      locMCTKBank->mctk[6].cz = rRp1->Pz()/rRp1->Rho();
      locMCTKBank->mctk[6].pmom = rRp1->Rho();
      locMCTKBank->mctk[6].mass = mass_kaonPDG;
      locMCTKBank->mctk[6].charge = 1;
      locMCTKBank->mctk[6].id = kaon_numPDG;
      locMCTKBank->mctk[6].beg_vtx = 1;
      locMCTKBank->mctk[6].end_vtx = 0;
      locMCTKBank->mctk[6].flag = 1;
      locMCTKBank->mctk[6].parent = 0;
        
      locMCTKBank->mctk[7].cx = rRp2->Px()/rRp2->Rho();
      locMCTKBank->mctk[7].cy = rRp2->Py()/rRp2->Rho();
      locMCTKBank->mctk[7].cz = rRp2->Pz()/rRp2->Rho();
      locMCTKBank->mctk[7].pmom = rRp2->Rho();
      locMCTKBank->mctk[7].mass = mass_neutPDG;
      locMCTKBank->mctk[7].charge = 0;
      locMCTKBank->mctk[7].id = neut_numPDG;
      locMCTKBank->mctk[7].beg_vtx = 1;
      locMCTKBank->mctk[7].end_vtx = 0;
      locMCTKBank->mctk[7].flag = 1;
      locMCTKBank->mctk[7].parent = 0;

      locMCTKBank->mctk[8].cx = rDp1->Px()/rDp1->Rho();
      locMCTKBank->mctk[8].cy = rDp1->Py()/rDp1->Rho();
      locMCTKBank->mctk[8].cz = rDp1->Pz()/rDp1->Rho();
      locMCTKBank->mctk[8].pmom = rDp1->Rho();
      locMCTKBank->mctk[8].mass = mass_protPDG;
      locMCTKBank->mctk[8].charge = 1;
      locMCTKBank->mctk[8].id = prot_numPDG;
      locMCTKBank->mctk[8].beg_vtx = 2;
      locMCTKBank->mctk[8].end_vtx = 0;
      locMCTKBank->mctk[8].flag = 1;
      locMCTKBank->mctk[8].parent = 0;

      locMCTKBank->mctk[9].cx = rDp2->Px()/rDp2->Rho();
      locMCTKBank->mctk[9].cy = rDp2->Py()/rDp2->Rho();
      locMCTKBank->mctk[9].cz = rDp2->Pz()/rDp2->Rho();
      locMCTKBank->mctk[9].pmom = rDp2->Rho();
      locMCTKBank->mctk[9].mass = mass_piminusPDG;
      locMCTKBank->mctk[9].charge = -1;
      locMCTKBank->mctk[9].id = piminus_numPDG ;
      locMCTKBank->mctk[9].beg_vtx = 2;
      locMCTKBank->mctk[9].end_vtx = 0;
      locMCTKBank->mctk[9].flag = 1;
      locMCTKBank->mctk[9].parent = 0;
    }
    else if (rlabel == 3 ){
      //Decay vertex
      locMCVXBank->mcvx[1].x = vertex_x+r3->Exp(lifetime_Lambda)*rRp1->Px()/mass_lambdaPDG*speed_light;
      locMCVXBank->mcvx[1].y = vertex_y+r3->Exp(lifetime_Lambda)*rRp1->Py()/mass_lambdaPDG*speed_light;
      locMCVXBank->mcvx[1].z = vertex_z+r3->Exp(lifetime_Lambda)*rRp1->Pz()/mass_lambdaPDG*speed_light;
      locMCVXBank->mcvx[1].tof = 0.0;
      locMCVXBank->mcvx[1].flag = 0;

      locMCTKBank->mctk[1].cx = rtarget->Px()/rtarget->Rho();
      locMCTKBank->mctk[1].cy = rtarget->Py()/rtarget->Rho();
      locMCTKBank->mctk[1].cz = rtarget->Pz()/rtarget->Rho();
      locMCTKBank->mctk[1].pmom = rtarget->Rho();
      locMCTKBank->mctk[1].mass = rtarget->M();
      locMCTKBank->mctk[1].charge = 1;
      locMCTKBank->mctk[1].id = prot_numPDG;
      locMCTKBank->mctk[1].beg_vtx = 0;
      locMCTKBank->mctk[1].end_vtx = 1;
      locMCTKBank->mctk[1].flag = rlabel;
      locMCTKBank->mctk[1].parent = 0;
            
      locMCTKBank->mctk[2].cx = rQp1->Px()/rQp1->Rho();
      locMCTKBank->mctk[2].cy = rQp1->Py()/rQp1->Rho();
      locMCTKBank->mctk[2].cz = rQp1->Pz()/rQp1->Rho();
      locMCTKBank->mctk[2].pmom = rQp1->Rho();
      locMCTKBank->mctk[2].mass = mass_lambdaPDG;
      locMCTKBank->mctk[2].charge = 0;
      locMCTKBank->mctk[2].id = lambda_numPDG;
      locMCTKBank->mctk[2].beg_vtx = 0;
      locMCTKBank->mctk[2].end_vtx = 1;
      locMCTKBank->mctk[2].flag = rlabel;
      locMCTKBank->mctk[2].parent = 0;
            
      locMCTKBank->mctk[3].cx = rQp2->Px()/rQp2->Rho();
      locMCTKBank->mctk[3].cy = rQp2->Py()/rQp2->Rho();
      locMCTKBank->mctk[3].cz = rQp2->Pz()/rQp2->Rho();
      locMCTKBank->mctk[3].pmom = rQp2->Rho();
      locMCTKBank->mctk[3].mass = mass_kaonPDG;
      locMCTKBank->mctk[3].charge = 1;
      locMCTKBank->mctk[3].id = kaon_numPDG;
      locMCTKBank->mctk[3].beg_vtx = 1;
      locMCTKBank->mctk[3].end_vtx = 0;
      locMCTKBank->mctk[3].flag = rlabel;
      locMCTKBank->mctk[3].parent = 0;
            
      locMCTKBank->mctk[4].cx = rSbeam->Px()/rSbeam->Rho();
      locMCTKBank->mctk[4].cy = rSbeam->Py()/rSbeam->Rho();
      locMCTKBank->mctk[4].cz = rSbeam->Pz()/rSbeam->Rho();
      locMCTKBank->mctk[4].pmom = rSbeam->Rho();
      locMCTKBank->mctk[4].mass = mass_lambdaPDG;
      locMCTKBank->mctk[4].charge = 0;
      locMCTKBank->mctk[4].id = lambda_numPDG;
      locMCTKBank->mctk[4].beg_vtx = 0;
      locMCTKBank->mctk[4].end_vtx = 1;
      locMCTKBank->mctk[4].flag = rlinpol*10+rcirpol;
      locMCTKBank->mctk[4].parent = 0;
            
      locMCTKBank->mctk[5].cx = rStarget->Px()/rStarget->Rho();
      locMCTKBank->mctk[5].cy = rStarget->Py()/rStarget->Rho();
      locMCTKBank->mctk[5].cz = rStarget->Pz()/rStarget->Rho();
      locMCTKBank->mctk[5].pmom = rStarget->Rho();
      locMCTKBank->mctk[5].mass = rStarget->M();
      locMCTKBank->mctk[5].charge = 0;
      locMCTKBank->mctk[5].id = neut_numPDG;
      locMCTKBank->mctk[5].beg_vtx = 0;
      locMCTKBank->mctk[5].end_vtx = 1;
      locMCTKBank->mctk[5].flag = 1;
      locMCTKBank->mctk[5].parent = 0;
            
      locMCTKBank->mctk[6].cx = rRp1->Px()/rRp1->Rho();
      locMCTKBank->mctk[6].cy = rRp1->Py()/rRp1->Rho();
      locMCTKBank->mctk[6].cz = rRp1->Pz()/rRp1->Rho();
      locMCTKBank->mctk[6].pmom = rRp1->Rho();
      locMCTKBank->mctk[6].mass = mass_lambdaPDG;
      locMCTKBank->mctk[6].charge = 0;
      locMCTKBank->mctk[6].id = lambda_numPDG;
      locMCTKBank->mctk[6].beg_vtx = 1;
      locMCTKBank->mctk[6].end_vtx = 0;
      locMCTKBank->mctk[6].flag = 1;
      locMCTKBank->mctk[6].parent = 0;
            
      locMCTKBank->mctk[7].cx = rRp2->Px()/rRp2->Rho();
      locMCTKBank->mctk[7].cy = rRp2->Py()/rRp2->Rho();
      locMCTKBank->mctk[7].cz = rRp2->Pz()/rRp2->Rho();
      locMCTKBank->mctk[7].pmom = rRp2->Rho();
      locMCTKBank->mctk[7].mass = mass_neutPDG;
      locMCTKBank->mctk[7].charge = 0;
      locMCTKBank->mctk[7].id = neut_numPDG;
      locMCTKBank->mctk[7].beg_vtx = 1;
      locMCTKBank->mctk[7].end_vtx = 0;
      locMCTKBank->mctk[7].flag = 1;
      locMCTKBank->mctk[7].parent = 0;

      locMCTKBank->mctk[8].cx = rDp1->Px()/rDp1->Rho();
      locMCTKBank->mctk[8].cy = rDp1->Py()/rDp1->Rho();
      locMCTKBank->mctk[8].cz = rDp1->Pz()/rDp1->Rho();
      locMCTKBank->mctk[8].pmom = rDp1->Rho();
      locMCTKBank->mctk[8].mass = mass_protPDG;
      locMCTKBank->mctk[8].charge = 1;
      locMCTKBank->mctk[8].id = prot_numPDG;
      locMCTKBank->mctk[8].beg_vtx = 2;
      locMCTKBank->mctk[8].end_vtx = 0;
      locMCTKBank->mctk[8].flag = 1;
      locMCTKBank->mctk[8].parent = 0;

      locMCTKBank->mctk[9].cx = rDp2->Px()/rDp2->Rho();
      locMCTKBank->mctk[9].cy = rDp2->Py()/rDp2->Rho();
      locMCTKBank->mctk[9].cz = rDp2->Pz()/rDp2->Rho();
      locMCTKBank->mctk[9].pmom = rDp2->Rho();
      locMCTKBank->mctk[9].mass = mass_piminusPDG;
      locMCTKBank->mctk[9].charge = -1;
      locMCTKBank->mctk[9].id = piminus_numPDG ;
      locMCTKBank->mctk[9].beg_vtx = 2;
      locMCTKBank->mctk[9].end_vtx = 0;
      locMCTKBank->mctk[9].flag = 1;
      locMCTKBank->mctk[9].parent = 0;
    }
    else if (rlabel == 4 ){
      //Only one vertex

      locMCTKBank->mctk[1].cx = rtarget->Px()/rtarget->Rho();
      locMCTKBank->mctk[1].cy = rtarget->Py()/rtarget->Rho();
      locMCTKBank->mctk[1].cz = rtarget->Pz()/rtarget->Rho();
      locMCTKBank->mctk[1].pmom = rtarget->Rho();
      locMCTKBank->mctk[1].mass = rtarget->M();
      locMCTKBank->mctk[1].charge = 1;
      locMCTKBank->mctk[1].id = prot_numPDG;
      locMCTKBank->mctk[1].beg_vtx = 0;
      locMCTKBank->mctk[1].end_vtx = 1;
      locMCTKBank->mctk[1].flag = rlabel;
      locMCTKBank->mctk[1].parent = 0;
            
      locMCTKBank->mctk[2].cx = rQp1->Px()/rQp1->Rho();
      locMCTKBank->mctk[2].cy = rQp1->Py()/rQp1->Rho();
      locMCTKBank->mctk[2].cz = rQp1->Pz()/rQp1->Rho();
      locMCTKBank->mctk[2].pmom = rQp1->Rho();
      locMCTKBank->mctk[2].mass = mass_sigmaPDG;
      locMCTKBank->mctk[2].charge = 0;
      locMCTKBank->mctk[2].id = sigma_numPDG;
      locMCTKBank->mctk[2].beg_vtx = 0;
      locMCTKBank->mctk[2].end_vtx = 1;
      locMCTKBank->mctk[2].flag = rlabel;
      locMCTKBank->mctk[2].parent = 0;
            
      locMCTKBank->mctk[3].cx = rQp2->Px()/rQp2->Rho();
      locMCTKBank->mctk[3].cy = rQp2->Py()/rQp2->Rho();
      locMCTKBank->mctk[3].cz = rQp2->Pz()/rQp2->Rho();
      locMCTKBank->mctk[3].pmom = rQp2->Rho();
      locMCTKBank->mctk[3].mass = mass_kaonPDG;
      locMCTKBank->mctk[3].charge = 1;
      locMCTKBank->mctk[3].id = kaon_numPDG;
      locMCTKBank->mctk[3].beg_vtx = 1;
      locMCTKBank->mctk[3].end_vtx = 0;
      locMCTKBank->mctk[3].flag = rlabel;
      locMCTKBank->mctk[3].parent = 0;
            
      locMCTKBank->mctk[4].cx = rSbeam->Px()/rSbeam->Rho();
      locMCTKBank->mctk[4].cy = rSbeam->Py()/rSbeam->Rho();
      locMCTKBank->mctk[4].cz = rSbeam->Pz()/rSbeam->Rho();
      locMCTKBank->mctk[4].pmom = rSbeam->Rho();
      locMCTKBank->mctk[4].mass = mass_sigmaPDG;
      locMCTKBank->mctk[4].charge = 0;
      locMCTKBank->mctk[4].id = sigma_numPDG;
      locMCTKBank->mctk[4].beg_vtx = 0;
      locMCTKBank->mctk[4].end_vtx = 1;
      locMCTKBank->mctk[4].flag = rlinpol*10+rcirpol;
      locMCTKBank->mctk[4].parent = 0;
            
      locMCTKBank->mctk[5].cx = rStarget->Px()/rStarget->Rho();
      locMCTKBank->mctk[5].cy = rStarget->Py()/rStarget->Rho();
      locMCTKBank->mctk[5].cz = rStarget->Pz()/rStarget->Rho();
      locMCTKBank->mctk[5].pmom = rStarget->Rho();
      locMCTKBank->mctk[5].mass = rStarget->M();
      locMCTKBank->mctk[5].charge = 0;
      locMCTKBank->mctk[5].id = neut_numPDG;
      locMCTKBank->mctk[5].beg_vtx = 0;
      locMCTKBank->mctk[5].end_vtx = 1;
      locMCTKBank->mctk[5].flag = 1;
      locMCTKBank->mctk[5].parent = 0;
            
      locMCTKBank->mctk[6].cx = rRp1->Px()/rRp1->Rho();
      locMCTKBank->mctk[6].cy = rRp1->Py()/rRp1->Rho();
      locMCTKBank->mctk[6].cz = rRp1->Pz()/rRp1->Rho();
      locMCTKBank->mctk[6].pmom = rRp1->Rho();
      locMCTKBank->mctk[6].mass = mass_sigmaPDG;
      locMCTKBank->mctk[6].charge = 0;
      locMCTKBank->mctk[6].id = sigma_numPDG;
      locMCTKBank->mctk[6].beg_vtx = 1;
      locMCTKBank->mctk[6].end_vtx = 0;
      locMCTKBank->mctk[6].flag = 1;
      locMCTKBank->mctk[6].parent = 0;
            
      locMCTKBank->mctk[7].cx = rRp2->Px()/rRp2->Rho();
      locMCTKBank->mctk[7].cy = rRp2->Py()/rRp2->Rho();
      locMCTKBank->mctk[7].cz = rRp2->Pz()/rRp2->Rho();
      locMCTKBank->mctk[7].pmom = rRp2->Rho();
      locMCTKBank->mctk[7].mass = mass_neutPDG;
      locMCTKBank->mctk[7].charge = 0;
      locMCTKBank->mctk[7].id = neut_numPDG;
      locMCTKBank->mctk[7].beg_vtx = 1;
      locMCTKBank->mctk[7].end_vtx = 0;
      locMCTKBank->mctk[7].flag = 1;
      locMCTKBank->mctk[7].parent = 0;

      locMCTKBank->mctk[8].cx = rbeam->Px()/rbeam->Rho();
      locMCTKBank->mctk[8].cy = rbeam->Py()/rbeam->Rho();;
      locMCTKBank->mctk[8].cz = rbeam->Pz()/rbeam->Rho();;
      locMCTKBank->mctk[8].pmom = rbeam->Rho();
      locMCTKBank->mctk[8].mass = mass_photonPDG;
      locMCTKBank->mctk[8].charge = 0;
      locMCTKBank->mctk[8].id = photon_numPDG ;
      locMCTKBank->mctk[8].beg_vtx = 0;
      locMCTKBank->mctk[8].end_vtx = 1;
      locMCTKBank->mctk[8].flag = 1;
      locMCTKBank->mctk[8].parent = 0;

      locMCTKBank->mctk[9].cx = rbeam->Px()/rbeam->Rho();
      locMCTKBank->mctk[9].cy = rbeam->Py()/rbeam->Rho();;
      locMCTKBank->mctk[9].cz = rbeam->Pz()/rbeam->Rho();;
      locMCTKBank->mctk[9].pmom = rbeam->Rho();
      locMCTKBank->mctk[9].mass = mass_photonPDG;
      locMCTKBank->mctk[9].charge = 0;
      locMCTKBank->mctk[9].id = photon_numPDG ;
      locMCTKBank->mctk[9].beg_vtx = 0;
      locMCTKBank->mctk[9].end_vtx = 1;
      locMCTKBank->mctk[9].flag = 1;
      locMCTKBank->mctk[9].parent = 0;
    }
    else if (rlabel == 5 ){
      //Decay vertex
      locMCVXBank->mcvx[1].x = vertex_x+r3->Exp(lifetime_Lambda)*rRp1->Px()/mass_lambdaPDG*speed_light;
      locMCVXBank->mcvx[1].y = vertex_y+r3->Exp(lifetime_Lambda)*rRp1->Py()/mass_lambdaPDG*speed_light;
      locMCVXBank->mcvx[1].z = vertex_z+r3->Exp(lifetime_Lambda)*rRp1->Pz()/mass_lambdaPDG*speed_light;
      locMCVXBank->mcvx[1].tof = 0.0;
      locMCVXBank->mcvx[1].flag = 0;

      locMCTKBank->mctk[1].cx = rtarget->Px()/rtarget->Rho();
      locMCTKBank->mctk[1].cy = rtarget->Py()/rtarget->Rho();
      locMCTKBank->mctk[1].cz = rtarget->Pz()/rtarget->Rho();
      locMCTKBank->mctk[1].pmom = rtarget->Rho();
      locMCTKBank->mctk[1].mass = rtarget->M();
      locMCTKBank->mctk[1].charge = 1;
      locMCTKBank->mctk[1].id = prot_numPDG;
      locMCTKBank->mctk[1].beg_vtx = 0;
      locMCTKBank->mctk[1].end_vtx = 1;
      locMCTKBank->mctk[1].flag = rlabel;
      locMCTKBank->mctk[1].parent = 0;
            
      locMCTKBank->mctk[2].cx = rQp1->Px()/rQp1->Rho();
      locMCTKBank->mctk[2].cy = rQp1->Py()/rQp1->Rho();
      locMCTKBank->mctk[2].cz = rQp1->Pz()/rQp1->Rho();
      locMCTKBank->mctk[2].pmom = rQp1->Rho();
      locMCTKBank->mctk[2].mass = mass_piplusPDG;
      locMCTKBank->mctk[2].charge = 1;
      locMCTKBank->mctk[2].id = piplus_numPDG;
      locMCTKBank->mctk[2].beg_vtx = 0;
      locMCTKBank->mctk[2].end_vtx = 1;
      locMCTKBank->mctk[2].flag = rlabel;
      locMCTKBank->mctk[2].parent = 0;
            
      locMCTKBank->mctk[3].cx = rQp2->Px()/rQp2->Rho();
      locMCTKBank->mctk[3].cy = rQp2->Py()/rQp2->Rho();
      locMCTKBank->mctk[3].cz = rQp2->Pz()/rQp2->Rho();
      locMCTKBank->mctk[3].pmom = rQp2->Rho();
      locMCTKBank->mctk[3].mass = mass_neutPDG;
      locMCTKBank->mctk[3].charge = 0;
      locMCTKBank->mctk[3].id = neut_numPDG;
      locMCTKBank->mctk[3].beg_vtx = 1;
      locMCTKBank->mctk[3].end_vtx = 0;
      locMCTKBank->mctk[3].flag = rlabel;
      locMCTKBank->mctk[3].parent = 0;
            
      locMCTKBank->mctk[4].cx = rSbeam->Px()/rSbeam->Rho();
      locMCTKBank->mctk[4].cy = rSbeam->Py()/rSbeam->Rho();
      locMCTKBank->mctk[4].cz = rSbeam->Pz()/rSbeam->Rho();
      locMCTKBank->mctk[4].pmom = rSbeam->Rho();
      locMCTKBank->mctk[4].mass = mass_piplusPDG;
      locMCTKBank->mctk[4].charge = 1;
      locMCTKBank->mctk[4].id = piplus_numPDG;
      locMCTKBank->mctk[4].beg_vtx = 0;
      locMCTKBank->mctk[4].end_vtx = 1;
      locMCTKBank->mctk[4].flag = rlinpol*10+rcirpol;
      locMCTKBank->mctk[4].parent = 0;
            
      locMCTKBank->mctk[5].cx = rStarget->Px()/rStarget->Rho();
      locMCTKBank->mctk[5].cy = rStarget->Py()/rStarget->Rho();
      locMCTKBank->mctk[5].cz = rStarget->Pz()/rStarget->Rho();
      locMCTKBank->mctk[5].pmom = rStarget->Rho();
      locMCTKBank->mctk[5].mass = rStarget->M();
      locMCTKBank->mctk[5].charge = 0;
      locMCTKBank->mctk[5].id = neut_numPDG;
      locMCTKBank->mctk[5].beg_vtx = 0;
      locMCTKBank->mctk[5].end_vtx = 1;
      locMCTKBank->mctk[5].flag = 1;
      locMCTKBank->mctk[5].parent = 0;
        
      locMCTKBank->mctk[6].cx = rRp1->Px()/rRp1->Rho();
      locMCTKBank->mctk[6].cy = rRp1->Py()/rRp1->Rho();
      locMCTKBank->mctk[6].cz = rRp1->Pz()/rRp1->Rho();
      locMCTKBank->mctk[6].pmom = rRp1->Rho();
      locMCTKBank->mctk[6].mass = mass_lambdaPDG;
      locMCTKBank->mctk[6].charge = 0;
      locMCTKBank->mctk[6].id = lambda_numPDG;
      locMCTKBank->mctk[6].beg_vtx = 1;
      locMCTKBank->mctk[6].end_vtx = 0;
      locMCTKBank->mctk[6].flag = 1;
      locMCTKBank->mctk[6].parent = 0;
            
      locMCTKBank->mctk[7].cx = rRp2->Px()/rRp2->Rho();
      locMCTKBank->mctk[7].cy = rRp2->Py()/rRp2->Rho();
      locMCTKBank->mctk[7].cz = rRp2->Pz()/rRp2->Rho();
      locMCTKBank->mctk[7].pmom = rRp2->Rho();
      locMCTKBank->mctk[7].mass = mass_kaonPDG;
      locMCTKBank->mctk[7].charge = 1;
      locMCTKBank->mctk[7].id = kaon_numPDG;
      locMCTKBank->mctk[7].beg_vtx = 1;
      locMCTKBank->mctk[7].end_vtx = 0;
      locMCTKBank->mctk[7].flag = 1;
      locMCTKBank->mctk[7].parent = 0;

      locMCTKBank->mctk[8].cx = rDp1->Px()/rDp1->Rho();
      locMCTKBank->mctk[8].cy = rDp1->Py()/rDp1->Rho();
      locMCTKBank->mctk[8].cz = rDp1->Pz()/rDp1->Rho();
      locMCTKBank->mctk[8].pmom = rDp1->Rho();
      locMCTKBank->mctk[8].mass = mass_protPDG;
      locMCTKBank->mctk[8].charge = 1;
      locMCTKBank->mctk[8].id = prot_numPDG;
      locMCTKBank->mctk[8].beg_vtx = 2;
      locMCTKBank->mctk[8].end_vtx = 0;
      locMCTKBank->mctk[8].flag = 1;
      locMCTKBank->mctk[8].parent = 0;

      locMCTKBank->mctk[9].cx = rDp2->Px()/rDp2->Rho();
      locMCTKBank->mctk[9].cy = rDp2->Py()/rDp2->Rho();
      locMCTKBank->mctk[9].cz = rDp2->Pz()/rDp2->Rho();
      locMCTKBank->mctk[9].pmom = rDp2->Rho();
      locMCTKBank->mctk[9].mass = mass_piminusPDG;
      locMCTKBank->mctk[9].charge = -1;
      locMCTKBank->mctk[9].id = piminus_numPDG ;
      locMCTKBank->mctk[9].beg_vtx = 2;
      locMCTKBank->mctk[9].end_vtx = 0;
      locMCTKBank->mctk[9].flag = 1;
      locMCTKBank->mctk[9].parent = 0;
    }
    else if (rlabel == 6 ){
      //Decay vertex
      locMCVXBank->mcvx[1].x = vertex_x+r3->Exp(lifetime_SigmaStar)*rQp1->Px()/mass_sigmastarPDG*speed_light;
      locMCVXBank->mcvx[1].y = vertex_y+r3->Exp(lifetime_SigmaStar)*rQp1->Py()/mass_sigmastarPDG*speed_light;
      locMCVXBank->mcvx[1].z = vertex_z+r3->Exp(lifetime_SigmaStar)*rQp1->Pz()/mass_sigmastarPDG*speed_light;
      locMCVXBank->mcvx[1].tof = 0.0;
      locMCVXBank->mcvx[1].flag = 0;

      locMCTKBank->mctk[1].cx = rtarget->Px()/rtarget->Rho();
      locMCTKBank->mctk[1].cy = rtarget->Py()/rtarget->Rho();
      locMCTKBank->mctk[1].cz = rtarget->Pz()/rtarget->Rho();
      locMCTKBank->mctk[1].pmom = rtarget->Rho();
      locMCTKBank->mctk[1].mass = rtarget->M();
      locMCTKBank->mctk[1].charge = 1;
      locMCTKBank->mctk[1].id = prot_numPDG;
      locMCTKBank->mctk[1].beg_vtx = 0;
      locMCTKBank->mctk[1].end_vtx = 1;
      locMCTKBank->mctk[1].flag = rlabel;
      locMCTKBank->mctk[1].parent = 0;
            
      locMCTKBank->mctk[2].cx = rQp1->Px()/rQp1->Rho();
      locMCTKBank->mctk[2].cy = rQp1->Py()/rQp1->Rho();
      locMCTKBank->mctk[2].cz = rQp1->Pz()/rQp1->Rho();
      locMCTKBank->mctk[2].pmom = rQp1->Rho();
      locMCTKBank->mctk[2].mass = mass_sigmastarPDG;
      locMCTKBank->mctk[2].charge = 0;
      locMCTKBank->mctk[2].id = sigmastar_numPDG;
      locMCTKBank->mctk[2].beg_vtx = 0;
      locMCTKBank->mctk[2].end_vtx = 2;
      locMCTKBank->mctk[2].flag = rlabel;
      locMCTKBank->mctk[2].parent = 0;
            
      locMCTKBank->mctk[3].cx = rQp2->Px()/rQp2->Rho();
      locMCTKBank->mctk[3].cy = rQp2->Py()/rQp2->Rho();
      locMCTKBank->mctk[3].cz = rQp2->Pz()/rQp2->Rho();
      locMCTKBank->mctk[3].pmom = rQp2->Rho();
      locMCTKBank->mctk[3].mass = mass_kaonPDG;
      locMCTKBank->mctk[3].charge = 1;
      locMCTKBank->mctk[3].id = kaon_numPDG;
      locMCTKBank->mctk[3].beg_vtx = 1;
      locMCTKBank->mctk[3].end_vtx = 0;
      locMCTKBank->mctk[3].flag = rlabel;
      locMCTKBank->mctk[3].parent = 0;
            
      locMCTKBank->mctk[4].cx = rbeam->Px()/rbeam->Rho();
      locMCTKBank->mctk[4].cy = rbeam->Py()/rbeam->Rho();;
      locMCTKBank->mctk[4].cz = rbeam->Pz()/rbeam->Rho();;
      locMCTKBank->mctk[4].pmom = rbeam->Rho();
      locMCTKBank->mctk[4].mass = mass_photonPDG;
      locMCTKBank->mctk[4].charge = 0;
      locMCTKBank->mctk[4].id = photon_numPDG ;
      locMCTKBank->mctk[4].beg_vtx = 0;
      locMCTKBank->mctk[4].end_vtx = 1;
      locMCTKBank->mctk[4].flag = rlinpol*10+rcirpol;
      locMCTKBank->mctk[4].parent = 0;
            
      locMCTKBank->mctk[5].cx = rStarget->Px()/rStarget->Rho();
      locMCTKBank->mctk[5].cy = rStarget->Py()/rStarget->Rho();
      locMCTKBank->mctk[5].cz = rStarget->Pz()/rStarget->Rho();
      locMCTKBank->mctk[5].pmom = rStarget->Rho();
      locMCTKBank->mctk[5].mass = rStarget->M();
      locMCTKBank->mctk[5].charge = 1;
      locMCTKBank->mctk[5].id = neut_numPDG;
      locMCTKBank->mctk[5].beg_vtx = 1;
      locMCTKBank->mctk[5].end_vtx = 0;
      locMCTKBank->mctk[5].flag = 1;
      locMCTKBank->mctk[5].parent = 0;
            
      locMCTKBank->mctk[6].cx = rbeam->Px()/rbeam->Rho();
      locMCTKBank->mctk[6].cy = rbeam->Py()/rbeam->Rho();;
      locMCTKBank->mctk[6].cz = rbeam->Pz()/rbeam->Rho();;
      locMCTKBank->mctk[6].pmom = rbeam->Rho();
      locMCTKBank->mctk[6].mass = mass_photonPDG;
      locMCTKBank->mctk[6].charge = 0;
      locMCTKBank->mctk[6].id = photon_numPDG ;
      locMCTKBank->mctk[6].beg_vtx = 0;
      locMCTKBank->mctk[6].end_vtx = 1;
      locMCTKBank->mctk[6].flag = 1;
      locMCTKBank->mctk[6].parent = 0;
            
      locMCTKBank->mctk[7].cx = rbeam->Px()/rbeam->Rho();
      locMCTKBank->mctk[7].cy = rbeam->Py()/rbeam->Rho();;
      locMCTKBank->mctk[7].cz = rbeam->Pz()/rbeam->Rho();;
      locMCTKBank->mctk[7].pmom = rbeam->Rho();
      locMCTKBank->mctk[7].mass = mass_photonPDG;
      locMCTKBank->mctk[7].charge = 0;
      locMCTKBank->mctk[7].id = photon_numPDG ;
      locMCTKBank->mctk[7].beg_vtx = 0;
      locMCTKBank->mctk[7].end_vtx = 1;
      locMCTKBank->mctk[7].flag = 1;
      locMCTKBank->mctk[7].parent = 0;

      locMCTKBank->mctk[8].cx = rDp1->Px()/rRp1->Rho();
      locMCTKBank->mctk[8].cy = rDp1->Py()/rRp1->Rho();
      locMCTKBank->mctk[8].cz = rDp1->Pz()/rRp1->Rho();
      locMCTKBank->mctk[8].pmom = rDp1->Rho();
      locMCTKBank->mctk[8].mass = mass_lambdaPDG;
      locMCTKBank->mctk[8].charge = 0;
      locMCTKBank->mctk[8].id = lambda_numPDG;
      locMCTKBank->mctk[8].beg_vtx = 2;
      locMCTKBank->mctk[8].end_vtx = 0;
      locMCTKBank->mctk[8].flag = 1;
      locMCTKBank->mctk[8].parent = 0;
            
      locMCTKBank->mctk[9].cx = rDp2->Px()/rRp2->Rho();
      locMCTKBank->mctk[9].cy = rDp2->Py()/rRp2->Rho();
      locMCTKBank->mctk[9].cz = rDp2->Pz()/rRp2->Rho();
      locMCTKBank->mctk[9].pmom = rDp2->Rho();
      locMCTKBank->mctk[9].mass = mass_pizeroPDG;
      locMCTKBank->mctk[9].charge = 0;
      locMCTKBank->mctk[9].id = pizero_numPDG;
      locMCTKBank->mctk[9].beg_vtx = 2;
      locMCTKBank->mctk[9].end_vtx = 0;
      locMCTKBank->mctk[9].flag = 1;
      locMCTKBank->mctk[9].parent = 0;
    }
    else if (rlabel == 7 ){
      //Decay vertex
      locMCVXBank->mcvx[1].x = vertex_x+r3->Exp(lifetime_SigmaStarMinus)*rQp1->Px()/mass_sigmastarminusPDG*speed_light;
      locMCVXBank->mcvx[1].y = vertex_y+r3->Exp(lifetime_SigmaStarMinus)*rQp1->Py()/mass_sigmastarminusPDG*speed_light;
      locMCVXBank->mcvx[1].z = vertex_z+r3->Exp(lifetime_SigmaStarMinus)*rQp1->Pz()/mass_sigmastarminusPDG*speed_light;
      locMCVXBank->mcvx[1].tof = 0.0;
      locMCVXBank->mcvx[1].flag = 0;

      locMCTKBank->mctk[1].cx = rtarget->Px()/rtarget->Rho();
      locMCTKBank->mctk[1].cy = rtarget->Py()/rtarget->Rho();
      locMCTKBank->mctk[1].cz = rtarget->Pz()/rtarget->Rho();
      locMCTKBank->mctk[1].pmom = rtarget->Rho();
      locMCTKBank->mctk[1].mass = rtarget->M();
      locMCTKBank->mctk[1].charge = 0;
      locMCTKBank->mctk[1].id = neut_numPDG;
      locMCTKBank->mctk[1].beg_vtx = 0;
      locMCTKBank->mctk[1].end_vtx = 1;
      locMCTKBank->mctk[1].flag = rlabel;
      locMCTKBank->mctk[1].parent = 0;
            
      locMCTKBank->mctk[2].cx = rQp1->Px()/rQp1->Rho();
      locMCTKBank->mctk[2].cy = rQp1->Py()/rQp1->Rho();
      locMCTKBank->mctk[2].cz = rQp1->Pz()/rQp1->Rho();
      locMCTKBank->mctk[2].pmom = rQp1->Rho();
      locMCTKBank->mctk[2].mass = mass_sigmastarminusPDG;
      locMCTKBank->mctk[2].charge = -1;
      locMCTKBank->mctk[2].id = sigmastarminus_numPDG;
      locMCTKBank->mctk[2].beg_vtx = 0;
      locMCTKBank->mctk[2].end_vtx = 2;
      locMCTKBank->mctk[2].flag = rlabel;
      locMCTKBank->mctk[2].parent = 0;
            
      locMCTKBank->mctk[3].cx = rQp2->Px()/rQp2->Rho();
      locMCTKBank->mctk[3].cy = rQp2->Py()/rQp2->Rho();
      locMCTKBank->mctk[3].cz = rQp2->Pz()/rQp2->Rho();
      locMCTKBank->mctk[3].pmom = rQp2->Rho();
      locMCTKBank->mctk[3].mass = mass_kaonPDG;
      locMCTKBank->mctk[3].charge = 1;
      locMCTKBank->mctk[3].id = kaon_numPDG;
      locMCTKBank->mctk[3].beg_vtx = 1;
      locMCTKBank->mctk[3].end_vtx = 0;
      locMCTKBank->mctk[3].flag = rlabel;
      locMCTKBank->mctk[3].parent = 0;
            
      locMCTKBank->mctk[4].cx = rbeam->Px()/rbeam->Rho();
      locMCTKBank->mctk[4].cy = rbeam->Py()/rbeam->Rho();;
      locMCTKBank->mctk[4].cz = rbeam->Pz()/rbeam->Rho();;
      locMCTKBank->mctk[4].pmom = rbeam->Rho();
      locMCTKBank->mctk[4].mass = mass_photonPDG;
      locMCTKBank->mctk[4].charge = 0;
      locMCTKBank->mctk[4].id = photon_numPDG ;
      locMCTKBank->mctk[4].beg_vtx = 0;
      locMCTKBank->mctk[4].end_vtx = 1;
      locMCTKBank->mctk[4].flag = rlinpol*10+rcirpol;
      locMCTKBank->mctk[4].parent = 0;
            
      locMCTKBank->mctk[5].cx = rStarget->Px()/rStarget->Rho();
      locMCTKBank->mctk[5].cy = rStarget->Py()/rStarget->Rho();
      locMCTKBank->mctk[5].cz = rStarget->Pz()/rStarget->Rho();
      locMCTKBank->mctk[5].pmom = rStarget->Rho();
      locMCTKBank->mctk[5].mass = rStarget->M();
      locMCTKBank->mctk[5].charge = 1;
      locMCTKBank->mctk[5].id = prot_numPDG;
      locMCTKBank->mctk[5].beg_vtx = 1;
      locMCTKBank->mctk[5].end_vtx = 0;
      locMCTKBank->mctk[5].flag = 1;
      locMCTKBank->mctk[5].parent = 0;
            
      locMCTKBank->mctk[6].cx = rbeam->Px()/rbeam->Rho();
      locMCTKBank->mctk[6].cy = rbeam->Py()/rbeam->Rho();;
      locMCTKBank->mctk[6].cz = rbeam->Pz()/rbeam->Rho();;
      locMCTKBank->mctk[6].pmom = rbeam->Rho();
      locMCTKBank->mctk[6].mass = mass_photonPDG;
      locMCTKBank->mctk[6].charge = 0;
      locMCTKBank->mctk[6].id = photon_numPDG ;
      locMCTKBank->mctk[6].beg_vtx = 0;
      locMCTKBank->mctk[6].end_vtx = 1;
      locMCTKBank->mctk[6].flag = 1;
      locMCTKBank->mctk[6].parent = 0;
            
      locMCTKBank->mctk[7].cx = rbeam->Px()/rbeam->Rho();
      locMCTKBank->mctk[7].cy = rbeam->Py()/rbeam->Rho();;
      locMCTKBank->mctk[7].cz = rbeam->Pz()/rbeam->Rho();;
      locMCTKBank->mctk[7].pmom = rbeam->Rho();
      locMCTKBank->mctk[7].mass = mass_photonPDG;
      locMCTKBank->mctk[7].charge = 0;
      locMCTKBank->mctk[7].id = photon_numPDG ;
      locMCTKBank->mctk[7].beg_vtx = 0;
      locMCTKBank->mctk[7].end_vtx = 1;
      locMCTKBank->mctk[7].flag = 1;
      locMCTKBank->mctk[7].parent = 0;

      locMCTKBank->mctk[8].cx = rDp1->Px()/rRp1->Rho();
      locMCTKBank->mctk[8].cy = rDp1->Py()/rRp1->Rho();
      locMCTKBank->mctk[8].cz = rDp1->Pz()/rRp1->Rho();
      locMCTKBank->mctk[8].pmom = rDp1->Rho();
      locMCTKBank->mctk[8].mass = mass_lambdaPDG;
      locMCTKBank->mctk[8].charge = 0;
      locMCTKBank->mctk[8].id = lambda_numPDG;
      locMCTKBank->mctk[8].beg_vtx = 2;
      locMCTKBank->mctk[8].end_vtx = 0;
      locMCTKBank->mctk[8].flag = 1;
      locMCTKBank->mctk[8].parent = 0;
            
      locMCTKBank->mctk[9].cx = rDp2->Px()/rRp2->Rho();
      locMCTKBank->mctk[9].cy = rDp2->Py()/rRp2->Rho();
      locMCTKBank->mctk[9].cz = rDp2->Pz()/rRp2->Rho();
      locMCTKBank->mctk[9].pmom = rDp2->Rho();
      locMCTKBank->mctk[9].mass = mass_piminusPDG;
      locMCTKBank->mctk[9].charge = -1;
      locMCTKBank->mctk[9].id = piminus_numPDG;
      locMCTKBank->mctk[9].beg_vtx = 2;
      locMCTKBank->mctk[9].end_vtx = 0;
      locMCTKBank->mctk[9].flag = 1;
      locMCTKBank->mctk[9].parent = 0;
    }
    else if (rlabel == 8 ){
      //Decay vertex
      locMCVXBank->mcvx[1].x = vertex_x+r3->Exp(lifetime_Lambda)*rQp1->Px()/mass_lambdaPDG*speed_light;
      locMCVXBank->mcvx[1].y = vertex_y+r3->Exp(lifetime_Lambda)*rQp1->Py()/mass_lambdaPDG*speed_light;
      locMCVXBank->mcvx[1].z = vertex_z+r3->Exp(lifetime_Lambda)*rQp1->Pz()/mass_lambdaPDG*speed_light;
      locMCVXBank->mcvx[1].tof = 0.0;
      locMCVXBank->mcvx[1].flag = 0;

      locMCTKBank->mctk[1].cx = rtarget->Px()/rtarget->Rho();
      locMCTKBank->mctk[1].cy = rtarget->Py()/rtarget->Rho();
      locMCTKBank->mctk[1].cz = rtarget->Pz()/rtarget->Rho();
      locMCTKBank->mctk[1].pmom = rtarget->Rho();
      locMCTKBank->mctk[1].mass = rtarget->M();
      locMCTKBank->mctk[1].charge = 1;
      locMCTKBank->mctk[1].id = prot_numPDG;
      locMCTKBank->mctk[1].beg_vtx = 0;
      locMCTKBank->mctk[1].end_vtx = 1;
      locMCTKBank->mctk[1].flag = rlabel;
      locMCTKBank->mctk[1].parent = 0;
            
      locMCTKBank->mctk[2].cx = rQp1->Px()/rQp1->Rho();
      locMCTKBank->mctk[2].cy = rQp1->Py()/rQp1->Rho();
      locMCTKBank->mctk[2].cz = rQp1->Pz()/rQp1->Rho();
      locMCTKBank->mctk[2].pmom = rQp1->Rho();
      locMCTKBank->mctk[2].mass = mass_lambdaPDG;
      locMCTKBank->mctk[2].charge = 0;
      locMCTKBank->mctk[2].id = lambda_numPDG;
      locMCTKBank->mctk[2].beg_vtx = 0;
      locMCTKBank->mctk[2].end_vtx = 2;
      locMCTKBank->mctk[2].flag = rlabel;
      locMCTKBank->mctk[2].parent = 0;
            
      locMCTKBank->mctk[3].cx = rQp2->Px()/rQp2->Rho();
      locMCTKBank->mctk[3].cy = rQp2->Py()/rQp2->Rho();
      locMCTKBank->mctk[3].cz = rQp2->Pz()/rQp2->Rho();
      locMCTKBank->mctk[3].pmom = rQp2->Rho();
      locMCTKBank->mctk[3].mass = mass_kaonPDG;
      locMCTKBank->mctk[3].charge = 1;
      locMCTKBank->mctk[3].id = kaon_numPDG;
      locMCTKBank->mctk[3].beg_vtx = 1;
      locMCTKBank->mctk[3].end_vtx = 0;
      locMCTKBank->mctk[3].flag = rlabel;
      locMCTKBank->mctk[3].parent = 0;
            
      locMCTKBank->mctk[4].cx = rbeam->Px()/rbeam->Rho();
      locMCTKBank->mctk[4].cy = rbeam->Py()/rbeam->Rho();;
      locMCTKBank->mctk[4].cz = rbeam->Pz()/rbeam->Rho();;
      locMCTKBank->mctk[4].pmom = rbeam->Rho();
      locMCTKBank->mctk[4].mass = mass_photonPDG;
      locMCTKBank->mctk[4].charge = 0;
      locMCTKBank->mctk[4].id = photon_numPDG ;
      locMCTKBank->mctk[4].beg_vtx = 0;
      locMCTKBank->mctk[4].end_vtx = 1;
      locMCTKBank->mctk[4].flag = rlinpol*10+rcirpol;
      locMCTKBank->mctk[4].parent = 0;
            
      locMCTKBank->mctk[5].cx = rStarget->Px()/rStarget->Rho();
      locMCTKBank->mctk[5].cy = rStarget->Py()/rStarget->Rho();
      locMCTKBank->mctk[5].cz = rStarget->Pz()/rStarget->Rho();
      locMCTKBank->mctk[5].pmom = rStarget->Rho();
      locMCTKBank->mctk[5].mass = rStarget->M();
      locMCTKBank->mctk[5].charge = 1;
      locMCTKBank->mctk[5].id = neut_numPDG;
      locMCTKBank->mctk[5].beg_vtx = 1;
      locMCTKBank->mctk[5].end_vtx = 0;
      locMCTKBank->mctk[5].flag = 1;
      locMCTKBank->mctk[5].parent = 0;
            
      locMCTKBank->mctk[6].cx = rbeam->Px()/rbeam->Rho();
      locMCTKBank->mctk[6].cy = rbeam->Py()/rbeam->Rho();;
      locMCTKBank->mctk[6].cz = rbeam->Pz()/rbeam->Rho();;
      locMCTKBank->mctk[6].pmom = rbeam->Rho();
      locMCTKBank->mctk[6].mass = mass_photonPDG;
      locMCTKBank->mctk[6].charge = 0;
      locMCTKBank->mctk[6].id = photon_numPDG ;
      locMCTKBank->mctk[6].beg_vtx = 0;
      locMCTKBank->mctk[6].end_vtx = 1;
      locMCTKBank->mctk[6].flag = 1;
      locMCTKBank->mctk[6].parent = 0;
            
      locMCTKBank->mctk[7].cx = rbeam->Px()/rbeam->Rho();
      locMCTKBank->mctk[7].cy = rbeam->Py()/rbeam->Rho();;
      locMCTKBank->mctk[7].cz = rbeam->Pz()/rbeam->Rho();;
      locMCTKBank->mctk[7].pmom = rbeam->Rho();
      locMCTKBank->mctk[7].mass = mass_photonPDG;
      locMCTKBank->mctk[7].charge = 0;
      locMCTKBank->mctk[7].id = photon_numPDG ;
      locMCTKBank->mctk[7].beg_vtx = 0;
      locMCTKBank->mctk[7].end_vtx = 1;
      locMCTKBank->mctk[7].flag = 1;
      locMCTKBank->mctk[7].parent = 0;
            
      locMCTKBank->mctk[8].cx = rDp1->Px()/rDp1->Rho();
      locMCTKBank->mctk[8].cy = rDp1->Py()/rDp1->Rho();
      locMCTKBank->mctk[8].cz = rDp1->Pz()/rDp1->Rho();
      locMCTKBank->mctk[8].pmom = rDp1->Rho();
      locMCTKBank->mctk[8].mass = mass_protPDG;
      locMCTKBank->mctk[8].charge = 1;
      locMCTKBank->mctk[8].id = prot_numPDG;
      locMCTKBank->mctk[8].beg_vtx = 2;
      locMCTKBank->mctk[8].end_vtx = 0;
      locMCTKBank->mctk[8].flag = 1;
      locMCTKBank->mctk[8].parent = 0;

      locMCTKBank->mctk[9].cx = rDp2->Px()/rDp2->Rho();
      locMCTKBank->mctk[9].cy = rDp2->Py()/rDp2->Rho();
      locMCTKBank->mctk[9].cz = rDp2->Pz()/rDp2->Rho();
      locMCTKBank->mctk[9].pmom = rDp2->Rho();
      locMCTKBank->mctk[9].mass = mass_piminusPDG;
      locMCTKBank->mctk[9].charge = -1;
      locMCTKBank->mctk[9].id = piminus_numPDG ;
      locMCTKBank->mctk[9].beg_vtx = 2;
      locMCTKBank->mctk[9].end_vtx = 0;
      locMCTKBank->mctk[9].flag = 1;
      locMCTKBank->mctk[9].parent = 0;
    }
    else if (rlabel == 9 ){
      locMCTKBank->mctk[1].cx = rtarget->Px()/rtarget->Rho();
      locMCTKBank->mctk[1].cy = rtarget->Py()/rtarget->Rho();
      locMCTKBank->mctk[1].cz = rtarget->Pz()/rtarget->Rho();
      locMCTKBank->mctk[1].pmom = rtarget->Rho();
      locMCTKBank->mctk[1].mass = rtarget->M();
      locMCTKBank->mctk[1].charge = 1;
      locMCTKBank->mctk[1].id = prot_numPDG;
      locMCTKBank->mctk[1].beg_vtx = 0;
      locMCTKBank->mctk[1].end_vtx = 1;
      locMCTKBank->mctk[1].flag = rlabel;
      locMCTKBank->mctk[1].parent = 0;
            
      locMCTKBank->mctk[2].cx = rQp1->Px()/rQp1->Rho();
      locMCTKBank->mctk[2].cy = rQp1->Py()/rQp1->Rho();
      locMCTKBank->mctk[2].cz = rQp1->Pz()/rQp1->Rho();
      locMCTKBank->mctk[2].pmom = rQp1->Rho();
      locMCTKBank->mctk[2].mass = mass_sigmaPDG;
      locMCTKBank->mctk[2].charge = 0;
      locMCTKBank->mctk[2].id = sigma_numPDG;
      locMCTKBank->mctk[2].beg_vtx = 1;
      locMCTKBank->mctk[2].end_vtx = 0;
      locMCTKBank->mctk[2].flag = rlabel;
      locMCTKBank->mctk[2].parent = 0;
            
      locMCTKBank->mctk[3].cx = rQp2->Px()/rQp2->Rho();
      locMCTKBank->mctk[3].cy = rQp2->Py()/rQp2->Rho();
      locMCTKBank->mctk[3].cz = rQp2->Pz()/rQp2->Rho();
      locMCTKBank->mctk[3].pmom = rQp2->Rho();
      locMCTKBank->mctk[3].mass = mass_kaonPDG;
      locMCTKBank->mctk[3].charge = 1;
      locMCTKBank->mctk[3].id = kaon_numPDG;
      locMCTKBank->mctk[3].beg_vtx = 1;
      locMCTKBank->mctk[3].end_vtx = 0;
      locMCTKBank->mctk[3].flag = rlabel;
      locMCTKBank->mctk[3].parent = 0;
            
      locMCTKBank->mctk[4].cx = rbeam->Px()/rbeam->Rho();
      locMCTKBank->mctk[4].cy = rbeam->Py()/rbeam->Rho();;
      locMCTKBank->mctk[4].cz = rbeam->Pz()/rbeam->Rho();;
      locMCTKBank->mctk[4].pmom = rbeam->Rho();
      locMCTKBank->mctk[4].mass = mass_photonPDG;
      locMCTKBank->mctk[4].charge = 0;
      locMCTKBank->mctk[4].id = photon_numPDG ;
      locMCTKBank->mctk[4].beg_vtx = 0;
      locMCTKBank->mctk[4].end_vtx = 1;
      locMCTKBank->mctk[4].flag = rlinpol*10+rcirpol;
      locMCTKBank->mctk[4].parent = 0;
            
      locMCTKBank->mctk[5].cx = rStarget->Px()/rStarget->Rho();
      locMCTKBank->mctk[5].cy = rStarget->Py()/rStarget->Rho();
      locMCTKBank->mctk[5].cz = rStarget->Pz()/rStarget->Rho();
      locMCTKBank->mctk[5].pmom = rStarget->Rho();
      locMCTKBank->mctk[5].mass = rStarget->M();
      locMCTKBank->mctk[5].charge = 0;
      locMCTKBank->mctk[5].id = neut_numPDG;
      locMCTKBank->mctk[5].beg_vtx = 1;
      locMCTKBank->mctk[5].end_vtx = 0;
      locMCTKBank->mctk[5].flag = 1;
      locMCTKBank->mctk[5].parent = 0;
            
      locMCTKBank->mctk[6].cx = rbeam->Px()/rbeam->Rho();
      locMCTKBank->mctk[6].cy = rbeam->Py()/rbeam->Rho();;
      locMCTKBank->mctk[6].cz = rbeam->Pz()/rbeam->Rho();;
      locMCTKBank->mctk[6].pmom = rbeam->Rho();
      locMCTKBank->mctk[6].mass = mass_photonPDG;
      locMCTKBank->mctk[6].charge = 0;
      locMCTKBank->mctk[6].id = photon_numPDG ;
      locMCTKBank->mctk[6].beg_vtx = 0;
      locMCTKBank->mctk[6].end_vtx = 1;
      locMCTKBank->mctk[6].flag = 1;
      locMCTKBank->mctk[6].parent = 0;
            
      locMCTKBank->mctk[7].cx = rbeam->Px()/rbeam->Rho();
      locMCTKBank->mctk[7].cy = rbeam->Py()/rbeam->Rho();;
      locMCTKBank->mctk[7].cz = rbeam->Pz()/rbeam->Rho();;
      locMCTKBank->mctk[7].pmom = rbeam->Rho();
      locMCTKBank->mctk[7].mass = mass_photonPDG;
      locMCTKBank->mctk[7].charge = 0;
      locMCTKBank->mctk[7].id = photon_numPDG ;
      locMCTKBank->mctk[7].beg_vtx = 0;
      locMCTKBank->mctk[7].end_vtx = 1;
      locMCTKBank->mctk[7].flag = 1;
      locMCTKBank->mctk[7].parent = 0;
            
      locMCTKBank->mctk[8].cx = rbeam->Px()/rbeam->Rho();
      locMCTKBank->mctk[8].cy = rbeam->Py()/rbeam->Rho();;
      locMCTKBank->mctk[8].cz = rbeam->Pz()/rbeam->Rho();;
      locMCTKBank->mctk[8].pmom = rbeam->Rho();
      locMCTKBank->mctk[8].mass = mass_photonPDG;
      locMCTKBank->mctk[8].charge = 0;
      locMCTKBank->mctk[8].id = photon_numPDG ;
      locMCTKBank->mctk[8].beg_vtx = 0;
      locMCTKBank->mctk[8].end_vtx = 1;
      locMCTKBank->mctk[8].flag = 1;
      locMCTKBank->mctk[8].parent = 0;

      locMCTKBank->mctk[9].cx = rbeam->Px()/rbeam->Rho();
      locMCTKBank->mctk[9].cy = rbeam->Py()/rbeam->Rho();;
      locMCTKBank->mctk[9].cz = rbeam->Pz()/rbeam->Rho();;
      locMCTKBank->mctk[9].pmom = rbeam->Rho();
      locMCTKBank->mctk[9].mass = mass_photonPDG;
      locMCTKBank->mctk[9].charge = 0;
      locMCTKBank->mctk[9].id = photon_numPDG ;
      locMCTKBank->mctk[9].beg_vtx = 0;
      locMCTKBank->mctk[9].end_vtx = 1;
      locMCTKBank->mctk[9].flag = 1;
      locMCTKBank->mctk[9].parent = 0;
    }
        
    //for each event: write, drop, clean
    char* locCommand_Char = new char[locBankList.size() + 1];
    strncpy(locCommand_Char, locBankList.c_str(), locBankList.size() + 1);
    bosWrite(bThreadBOSIOptr_Output, bThreadBOSCommonBlock_Output.iw, locCommand_Char);
    bosLdrop(bThreadBOSCommonBlock_Output.iw, locCommand_Char);
    bosNgarbage(bThreadBOSCommonBlock_Output.iw);
        
        
    rtarget->Clear();
    rbeam->Clear();
    rQp1->Clear();
    rQp2->Clear();
    rStarget->Clear();
    rSbeam->Clear();
    rRp1->Clear();
    rRp2->Clear();
    rDp1->Clear();
    rDp2->Clear();
  }

  //cleanup at the end
  bosWrite(bThreadBOSIOptr_Output, bThreadBOSCommonBlock_Output.iw, "0");
    
  //close the file
  bosClose(bThreadBOSIOptr_Output);


}

void Display_Help(){
  cout << "DISPLAY HELP:" << endl;
  cout << "Command Line: ExecutableFile InputRootFile" << endl;
  cout << "The optional '-M' flag is used for specifying the maximum number of events to be evaluated from each file." << endl;
}


