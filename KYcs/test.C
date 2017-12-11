double test(){
  double Wv;
  double ctv;
  //cin>>Wv>>ctv;
  int Sstatus;


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
   
  double temp=0;
  for(int i=0;i<W_binnum;i++){
    for(int j=0;j<ct_binnum;j++){
      if(dcs[i][j]>temp) temp=dcs[i][j];
    }
  }
  cout<<temp<<endl; 

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
