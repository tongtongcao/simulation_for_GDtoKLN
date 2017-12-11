void combine(){

  ifstream input1("SigmaKL");
  ifstream input2("TKL");
  ifstream input3("OxKL");
  ifstream input4("OzKL");
  ifstream input5("PKL");

  ofstream output("STOxOzPforGPKL");

  double w,costheta,sigma,t,ox,oz,p;
  for(int i=0;i<192;i++){
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
