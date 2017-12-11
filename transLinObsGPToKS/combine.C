void combine(){

  ifstream input1("SigmaKS");
  ifstream input2("TKS");
  ifstream input3("OxKS");
  ifstream input4("OzKS");
  ifstream input5("PKS");

  ofstream output("STOxOzPforGPKS");

  double w,costheta,sigma,t,ox,oz,p;
  for(int i=0;i<88;i++){
    input1>>w>>costheta>>sigma;
    input2>>w>>costheta>>t;
    input3>>w>>costheta>>ox;
    input4>>w>>costheta>>oz;
    input5>>w>>costheta>>p;

    output<<setw(15)<<setiosflags(ios::left)<<w<<setw(15)<<setiosflags(ios::left)<<costheta<<setw(15)<<setiosflags(ios::left)<<sigma<<setw(15)<<setiosflags(ios::left)<<t<<setw(15)<<setiosflags(ios::left)<<ox<<setw(15)<<setiosflags(ios::left)<<oz<<setw(15)<<setiosflags(ios::left)<<p<<endl;
  }

  input1.close();
  input2.close();
  input3.close();
  input4.close();
  input5.close();

  output.close();


}
