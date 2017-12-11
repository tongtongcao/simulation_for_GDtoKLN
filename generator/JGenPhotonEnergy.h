#include <unistd.h>
#include <stdlib.h>
#include <fstream>
#include <stddef.h>
#include <iostream>
#include <signal.h>
#include <string.h>
#include <TROOT.h>
#include <TRandom3.h>
#include <TH2.h>
#include <TFile.h>
#include <TMath.h>

using namespace std;
using std::cout;
using std::endl;

enum Edistrib_t { mono, plain, brems, histo, file };
enum Jexcept { XnumberArgument, XmaxLowerMin, XfileNotOpen, XhistNotFound }; 

class JGenPhotonEnergy {
  Edistrib_t type;
  Double_t Evalue;
  Double_t Emin;
  Double_t Emax;
  Double_t Ebeam;
  Double_t* fSum; 
  Int_t nbins;
  Double_t wcont[800];
  Double_t eaver[800];
  Double_t x1, x2;
  TH1F*       hgen;
  TDirectory* genDir;
  void constructor (Edistrib_t type_, Double_t Emin_, Double_t Emax_);
  void constructor (char* filename, char* hstname);
  void constructor (char* filenameW, char* filenameE, double Eb);

public:
  Double_t ComputeSum();
  Double_t GetRndm();
  Double_t ConvertEcounter(Int_t ec);
  JGenPhotonEnergy (char * copt); 
  JGenPhotonEnergy (Double_t Efixed) {
    constructor (mono, Efixed, Efixed);}
  JGenPhotonEnergy (Double_t Emin_, Double_t Emax_) {
    constructor (plain, Emin_, Emax_); }
  JGenPhotonEnergy (Edistrib_t type_, Double_t Emin_, Double_t Emax_) {
    constructor (type_, Emin_, Emax_); }
  Double_t Generate ();
  Double_t GetE () { return Evalue; }
};
