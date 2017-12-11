void transferT(){

  ifstream input("T");

  ofstream output("TKL");

  double w, costheta,obs;
  double costhetacalc,obscalc;
  int wIncreaseIndex=0;
  int jindex;
  double min;
  int flag;

  vector<double> vectw, vectcostheta,vectobs;
  while(!input.eof()){
    input>>w>>costheta>>obs;
    if(w>1.72+wIncreaseIndex*0.02-0.01 && w<1.72+wIncreaseIndex*0.02+0.01){
      vectw.push_back(w);
      vectcostheta.push_back(costheta);
      vectobs.push_back(obs);
    }
    else{
      for(int i=0;i<8;i++){
	costhetacalc=-0.65+0.2*i;
	flag=0;
	for(int j=0;j<vectcostheta.size();j++){
	  if(vectcostheta[j]>costhetacalc-0.05 && vectcostheta[j]<costhetacalc+0.05){
	    output<<setw(15)<<setiosflags(ios::left)<<w-0.02<<setw(15)<<setiosflags(ios::left)<<costhetacalc<<setw(15)<<setiosflags(ios::left)<<vectobs[j]<<endl;
	    flag=1;
	    break;
	  }
	}
	if(flag==0){
	  if(costhetacalc<vectcostheta[0]){
	    output<<setw(15)<<setiosflags(ios::left)<<w-0.02<<setw(15)<<setiosflags(ios::left)<<costhetacalc<<setw(15)<<setiosflags(ios::left)<<vectobs[0]<<endl;
	  }
	  else if(costhetacalc>vectcostheta[vectcostheta.size()-1]){
	    output<<setw(15)<<setiosflags(ios::left)<<w-0.02<<setw(15)<<setiosflags(ios::left)<<costhetacalc<<setw(15)<<setiosflags(ios::left)<<vectobs[vectcostheta.size()-1]<<endl;
	  }
	  else{
	    min=2;
	    for(int j=0;j<vectcostheta.size();j++){
	      if(min>fabs(vectcostheta[j]-costhetacalc)){
		min=fabs(vectcostheta[j]-costhetacalc);
		jindex=j;
	      }
	    }
	    if(fabs(vectcostheta[jindex]-costhetacalc)<=0.100000000001){
	      obscalc=(vectobs[jindex+1]-vectobs[jindex])/(vectcostheta[jindex+1]-vectcostheta[jindex])*(costhetacalc-vectcostheta[jindex])+vectobs[jindex];
	      output<<setw(15)<<setiosflags(ios::left)<<w-0.02<<setw(15)<<setiosflags(ios::left)<<costhetacalc<<setw(15)<<setiosflags(ios::left)<<obscalc<<endl;
	    }
	    else{
	      obscalc=(vectobs[jindex-1]-vectobs[jindex])/(vectcostheta[jindex-1]-vectcostheta[jindex])*(costhetacalc-vectcostheta[jindex])+vectobs[jindex];
	      output<<setw(15)<<setiosflags(ios::left)<<w-0.02<<setw(15)<<setiosflags(ios::left)<<costhetacalc<<setw(15)<<setiosflags(ios::left)<<obscalc<<endl;
	    }
	  }
	}
      }
	  
      wIncreaseIndex++;
      vectw.clear();
      vectcostheta.clear();
      vectobs.clear();
      vectw.push_back(w);
      vectcostheta.push_back(costheta);
      vectobs.push_back(obs);
    }
  }

  for(int i=0;i<8;i++){
    costhetacalc=-0.65+0.2*i;
    flag=0;
    for(int j=0;j<vectcostheta.size();j++){
      if(vectcostheta[j]>costhetacalc-0.05 && vectcostheta[j]<costhetacalc+0.05){
	output<<setw(15)<<setiosflags(ios::left)<<w<<setw(15)<<setiosflags(ios::left)<<costhetacalc<<setw(15)<<setiosflags(ios::left)<<vectobs[j]<<endl;
	flag=1;
	break;
      }
    }
    if(flag==0){
      if(costhetacalc<vectcostheta[0]){
	output<<setw(15)<<setiosflags(ios::left)<<w<<setw(15)<<setiosflags(ios::left)<<costhetacalc<<setw(15)<<setiosflags(ios::left)<<vectobs[0]<<endl;
      }
      else if(costhetacalc>vectcostheta[vectcostheta.size()-1]){
	output<<setw(15)<<setiosflags(ios::left)<<w<<setw(15)<<setiosflags(ios::left)<<costhetacalc<<setw(15)<<setiosflags(ios::left)<<vectobs[vectcostheta.size()-1]<<endl;
      }
      else{
	min=2;
	for(int j=0;j<vectcostheta.size();j++){
	  if(min>fabs(vectcostheta[j]-costhetacalc)){
	    min=fabs(vectcostheta[j]-costhetacalc);
	    jindex=j;
	  }
	}
	if(fabs(vectcostheta[jindex]-costhetacalc)<=0.100000000001){
	  obscalc=(vectobs[jindex+1]-vectobs[jindex])/(vectcostheta[jindex+1]-vectcostheta[jindex])*(costhetacalc-vectcostheta[jindex])+vectobs[jindex];
	  output<<setw(15)<<setiosflags(ios::left)<<w<<setw(15)<<setiosflags(ios::left)<<costhetacalc<<setw(15)<<setiosflags(ios::left)<<obscalc<<endl;
	}
	else{
	  obscalc=(vectobs[jindex-1]-vectobs[jindex])/(vectcostheta[jindex-1]-vectcostheta[jindex])*(costhetacalc-vectcostheta[jindex])+vectobs[jindex];
	  output<<setw(15)<<setiosflags(ios::left)<<w<<setw(15)<<setiosflags(ios::left)<<costhetacalc<<setw(15)<<setiosflags(ios::left)<<obscalc<<endl;
	}
      }
    }
  }
	 


output.close();
}
    

