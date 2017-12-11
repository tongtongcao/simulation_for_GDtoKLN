#include "JGenPhotonEnergy.h"

// constructor parses -E option
JGenPhotonEnergy::JGenPhotonEnergy (char * coption_) {
  int   narg = 0;
  char  coption [256];
  char* copt [20];
  char* s = coption;
  char* keyword [] = { "mono", "plain", "brems", "histo", "file", NULL };
  int   jKey = 0;

  strncpy (coption, coption_, 255);   // copy line before manipulating
  int   len = strlen(coption);

  while ((s=strtok(s, ":")) != NULL) {  // delimiter is  ':'
    copt [narg++] = s;
    if ((len -= (strlen (s) + 1)) <= 0) break;
    s += strlen (s) + 1;
  }

  while ( keyword [jKey] && strcasecmp (keyword [jKey], copt [0])) {
    jKey++;
  }

  if (! keyword [jKey]) {   // no keyword found, assume float argument
    switch (narg) {
    case 1: 
      constructor (mono, atof(copt[0]), atof(copt[0]));
      break;
    case 2:
      constructor (plain, atof(copt[0]), atof(copt[1]));
      break;
    default:
      throw XnumberArgument;
    }
  }
  else {                    // keyword found, check number arguments
    switch ((Edistrib_t) jKey) {
    case mono:
      if (narg != 2) throw XnumberArgument;
      constructor (mono, atof(copt[1]), atof(copt[1]));
      break;
    case histo:
      if (narg != 3) throw XnumberArgument;
      constructor (copt[1], copt[2]);
      break;
    case file:
      if (narg != 4) throw XnumberArgument;
      constructor (copt[1], copt[2], atof(copt[3]));
      break;
    default:
      if (narg != 4) throw XnumberArgument;
      constructor ((Edistrib_t) jKey, atof(copt[1]), atof(copt[2]));
      break;
    }
  }
}

void JGenPhotonEnergy::constructor (Edistrib_t type_, 
				    Double_t Emin_, Double_t Emax_) {
  type = type_;
  Emin = Emin_;
  Emax = Emax_;

  if (Emax < Emin) throw XmaxLowerMin;

  switch (type) {
  case plain:
    x1 = Emax - Emin;
    break;
  case brems:
    x1 = log (Emin);
    x2 = log (Emax) - x1;
    break;
  default:
    break;
  }
}

void JGenPhotonEnergy::constructor (char* filename, char* hstname) {
  type = histo;
  genDir = new TDirectory ("gendir", "gendir");

  TFile f(filename);                    // open a file switches directory  
  if (!f.IsOpen ()) throw XfileNotOpen;
  TH1F* dummy = (TH1F*) f.Get(hstname); // access to histogram
  if (dummy == NULL) throw XhistNotFound;

  genDir -> cd ();                     // change dir back to root-memory
  hgen = new TH1F();                   // allocate memory
  *hgen = * ((TH1F*) dummy->Clone ()); // copy histogram
  hgen -> ComputeIntegral ();          // integral needed for GetRandom()
}

void JGenPhotonEnergy::constructor (char* filenameW, char* filenameE, double Eb) {
  type = file;
  Ebeam = Eb;
  
  Int_t Ec,Ecmin,Ecmax,num,numec;
  Double_t weight, eweight;
  Double_t egamma;
  Double_t random, Eaverage;
  Double_t egmin, egmax, eg, edev;

  ifstream fweight;
  fweight.open(filenameW,ios::in);

  ifstream fecounter;
  fecounter.open(filenameE,ios::in);

  int i=0;
  num = 0;
  while (fweight >> weight) {
   wcont[i]=weight;
   num++;
   i++;
             }
  nbins = num;
  std::cout << "Nbins=" << nbins << endl;
  ComputeSum();

  i=0;
  while (fecounter >> numec >> egmin >> egmax >> eg >> edev) {
   eaver[i] = eg;
   i++;
             }

}

Double_t JGenPhotonEnergy::Generate () {
  switch (type) {
  case mono:
    Evalue = Emin;
    break;
  case plain:
    Evalue = Emin + x1 * gRandom->Rndm();
    break;
  case brems:
    Evalue = TMath::Exp (x1 + x2 * gRandom->Rndm(0));
    break;
  case histo:
    Evalue = hgen->GetRandom();
    break;
  case file:
    Int_t Ecounter = int (GetRndm());
    Evalue = ConvertEcounter(Ecounter);
    break;
  }
  return Evalue;
}

Double_t JGenPhotonEnergy::ComputeSum() {
 Int_t bin, ibin;
  
 if (fSum) delete [] fSum;

 fSum = new Double_t[nbins+2]; 
 ibin=0;
 fSum[ibin] = 0;
 for (bin=1;bin<=nbins;bin++) {
   ibin++;
   fSum[ibin] = fSum[ibin-1] + wcont[bin-1];
 }

 fSum[nbins+1] = fSum[nbins];

//   - Normalize sum to 1
 if (fSum[nbins] ==0) {
     printf("ComputeSum", "Sum = zero"); return 0;
 }

// for (bin=0;bin<=nbins+1;bin++) cout << fSum[bin] << endl; 

 for (bin=1;bin<=nbins;bin++)  fSum[bin] /= fSum[nbins];
 return fSum[nbins];
}

Double_t JGenPhotonEnergy::GetRndm() {
 Double_t integral;
 Double_t LowEdge, BinWidth;
 Double_t tgen;

 if (fSum) {
  integral=ComputeSum();
 }
 if (integral == 0 || fSum == 0) return 0;

 Double_t r1 = gRandom->Rndm();
 Int_t ibin = TMath::BinarySearch(nbins,fSum,r1);

 BinWidth = (ibin+1)-(ibin);
 LowEdge = (ibin+1)-(BinWidth/2.);

 tgen = LowEdge + BinWidth*(r1-fSum[ibin])/(fSum[ibin+1] - fSum[ibin]);
 //cout << ibin << " " << BinWidth << " " << LowEdge << " " << endl;
// cout << r1 << endl;
 return(tgen);
}

Double_t JGenPhotonEnergy::ConvertEcounter(Int_t ec) {
 Double_t egamma;
 
 egamma = eaver[ec-1]*Ebeam;
 return(egamma);
}
