void separate(){

  string str;
  double w,wLow,wHigh,cosTheta,cosThetaLow,cosThetaHigh, obs,obsLow,obsHigh;

  ifstream input("gammaPToKLambda_G8.dat");

  ofstream outputSigma("Sigma");
  ofstream outputOx("Ox");
  ofstream outputP("P");
  ofstream outputT("T");
  ofstream outputOz("Oz");
  while(!input.eof()){
    input>>str>>w>>wLow>>wHigh>>cosTheta>>cosThetaLow>>cosThetaHigh>>obs>>obsLow>>obsHigh;
    if(str=="B") outputSigma<<setw(15)<<setiosflags(ios::left)<<w+(-wLow+wHigh)/2<<setw(15)<<setiosflags(ios::left)<<cosTheta+(-cosThetaLow+cosThetaHigh)/2<<setw(15)<<setiosflags(ios::left)<<obs<<endl;
    else if(str=="Ox") outputOx<<setw(15)<<setiosflags(ios::left)<<w+(-wLow+wHigh)/2<<setw(15)<<setiosflags(ios::left)<<cosTheta+(-cosThetaLow+cosThetaHigh)/2<<setw(15)<<setiosflags(ios::left)<<obs<<endl;
    else if(str=="R") outputP<<setw(15)<<setiosflags(ios::left)<<w+(-wLow+wHigh)/2<<setw(15)<<setiosflags(ios::left)<<cosTheta+(-cosThetaLow+cosThetaHigh)/2<<setw(15)<<setiosflags(ios::left)<<obs<<endl;
    else if(str=="T") outputT<<setw(15)<<setiosflags(ios::left)<<w+(-wLow+wHigh)/2<<setw(15)<<setiosflags(ios::left)<<cosTheta+(-cosThetaLow+cosThetaHigh)/2<<setw(15)<<setiosflags(ios::left)<<obs<<endl;
    else if(str=="Oz") outputOz<<setw(15)<<setiosflags(ios::left)<<w+(-wLow+wHigh)/2<<setw(15)<<setiosflags(ios::left)<<cosTheta+(-cosThetaLow+cosThetaHigh)/2<<setw(15)<<setiosflags(ios::left)<<obs<<endl;

  }
  outputSigma.close();
  outputOx.close();
  outputP.close();
  outputT.close();
  outputOz.close();

}
