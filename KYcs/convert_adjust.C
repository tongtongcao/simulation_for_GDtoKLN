#include <iostream>  
#include <fstream>
#include <sstream>
#include <string> 

using namespace std;
void convert_adjust(){

  {
    ifstream input("filename_KL");
    ofstream output("csKL_adjust");

    string str,strsub;
    int ch_first,ch_last;
    const int W_binnum=49;
    const int ct_binnum=50;
    double W[W_binnum],W_min=1626,W_max=2405,W_binwidth=(W_max-W_min)/(double)(W_binnum-1);

    double ct[W_binnum][ct_binnum],dcs[W_binnum][ct_binnum],ct_binwidth=2/(double)(ct_binnum-1);

    ostringstream sstr;
    for(int i=0;i<W_binnum;i++){
      getline(input,str);
      if (str.empty()) continue;
      ifstream infile(str.c_str());
      ch_first=str.find_last_of("L")+1;
      ch_last=str.find_last_of(".")-1;
      W[i]=atof(str.substr(ch_first,ch_last).c_str());
      for(int j=0;j<ct_binnum;j++){
	infile>>ct[i][j]>>dcs[i][j];
      }
    }

    double W_adjust[W_binnum],ct_adjust[W_binnum][ct_binnum],dcs_adjust[W_binnum][ct_binnum];
    for(int i=0;i<W_binnum;i++){
      W_adjust[i]=W_min+W_binwidth*i;
      for(int j=0;j<ct_binnum;j++){
	ct_adjust[i][j]=-1+ct_binwidth*j;
      }
    }

    double W1=0,W2=0,ct11=0,ct12=0,ct21=0,ct22=0,dcs11=0,dcs12=0,dcs21=0,dcs22=0,dcs1=0,dcs2=0;
    int ith;
    for(int u=0;u<W_binnum;u++){
      for(int i=0;i<W_binnum;i++){
	if(W_adjust[u]>=W[i] && W_adjust[u]<=W[i+1]){
	  W1=W[i];
	  W2=W[i+1];
	  ith=i;
	  break;
	}
      }
      for(int v=0;v<ct_binnum;v++){
	if(v==0){
	  dcs1=dcs[ith][0];
	  dcs2=dcs[ith+1][0];
	  dcs_adjust[u][v]=(W_adjust[u]-W1)/(W2-W1)*(dcs2-dcs1)+dcs1;
	  break;
	}
	else{
	  for(int j=0;j<ct_binnum;j++){
	    if(ct_adjust[u][v]>=ct[ith][j] && ct_adjust[u][v]<=ct[ith][j+1]){
	      ct11=ct[ith][j];
	      ct12=ct[ith][j+1];
	      dcs11=dcs[ith][j];
	      dcs12=dcs[ith][j+1];
	      break;
	    }
	  }
	  for(int j=0;j<ct_binnum;j++){
	    if(ct_adjust[u][v]>=ct[ith+1][j] && ct_adjust[u][v]<=ct[ith+1][j+1]){
	      ct21=ct[ith+1][j];
	      ct22=ct[ith+1][j+1];
	      dcs21=dcs[ith+1][j];
	      dcs22=dcs[ith+1][j+1];
	      break;
	    }
	  }
	 
	  dcs1=(ct_adjust[u][v]-ct11)/(ct12-ct11)*(dcs12-dcs11)+dcs11;
	  dcs2=(ct_adjust[u][v]-ct21)/(ct22-ct21)*(dcs22-dcs21)+dcs21;
	  dcs_adjust[u][v]=(W_adjust[u]-W1)/(W2-W1)*(dcs2-dcs1)+dcs1;
	  break;
	}
      }
    }

    for(int i=0;i<W_binnum;i++){
      sstr<<"W="<<W_adjust[i]<<endl;
      output<<sstr.str();
      sstr.str("");
      for(int j=0;j<ct_binnum;j++){
	output<<ct_adjust[i][j]<<"        "<<dcs_adjust[i][j]<<endl;
      }
      output<<endl;
    }

    input.close();
    output.close();
  }

{
    ifstream input("filename_KpSp");
    ofstream output("csKpSp_adjust");

    string str,strsub;
    int ch_first,ch_last;
    const int W_binnum=35;
    const int ct_binnum=50;
    double W[W_binnum],W_min=1700,W_max=2355,W_binwidth=(W_max-W_min)/(double)(W_binnum-1);

    double ct[W_binnum][ct_binnum],dcs[W_binnum][ct_binnum],ct_binwidth=2/(double)(ct_binnum-1);

    ostringstream sstr;
    for(int i=0;i<W_binnum;i++){
      getline(input,str);
      if (str.empty()) continue;
      ifstream infile(str.c_str());
      ch_first=str.find_last_of("p")+1;
      ch_last=str.find_last_of(".")-1;
      W[i]=atof(str.substr(ch_first,ch_last).c_str());
      for(int j=0;j<ct_binnum;j++){
	infile>>ct[i][j]>>dcs[i][j];
      }
    }

    double W_adjust[W_binnum],ct_adjust[W_binnum][ct_binnum],dcs_adjust[W_binnum][ct_binnum];
    for(int i=0;i<W_binnum;i++){
      W_adjust[i]=W_min+W_binwidth*i;
      for(int j=0;j<ct_binnum;j++){
	ct_adjust[i][j]=-1+ct_binwidth*j;
      }
    }

    double W1=0,W2=0,ct11=0,ct12=0,ct21=0,ct22=0,dcs11=0,dcs12=0,dcs21=0,dcs22=0,dcs1=0,dcs2=0;
    int ith;
    for(int u=0;u<W_binnum;u++){
      for(int i=0;i<W_binnum;i++){
	if(W_adjust[u]>=W[i] && W_adjust[u]<=W[i+1]){
	  W1=W[i];
	  W2=W[i+1];
	  ith=i;
	  break;
	}
      }
      for(int v=0;v<ct_binnum;v++){
	if(v==0){
	  dcs1=dcs[ith][0];
	  dcs2=dcs[ith+1][0];
	  dcs_adjust[u][v]=(W_adjust[u]-W1)/(W2-W1)*(dcs2-dcs1)+dcs1;
	}
	else{
	  for(int j=0;j<ct_binnum;j++){
	    if(ct_adjust[u][v]>=ct[ith][j] && ct_adjust[u][v]<=ct[ith][j+1]){
	      ct11=ct[ith][j];
	      ct12=ct[ith][j+1];
	      dcs11=dcs[ith][j];
	      dcs12=dcs[ith][j+1];
	      break;
	    }
	  }
	  for(int j=0;j<ct_binnum;j++){
	    if(ct_adjust[u][v]>=ct[ith+1][j] && ct_adjust[u][v]<=ct[ith+1][j+1]){
	      ct21=ct[ith+1][j];
	      ct22=ct[ith+1][j+1];
	      dcs21=dcs[ith+1][j];
	      dcs22=dcs[ith+1][j+1];
	      break;
	    }
	  }
	 
	  dcs1=(ct_adjust[u][v]-ct11)/(ct12-ct11)*(dcs12-dcs11)+dcs11;
	  dcs2=(ct_adjust[u][v]-ct21)/(ct22-ct21)*(dcs22-dcs21)+dcs21;
	  dcs_adjust[u][v]=(W_adjust[u]-W1)/(W2-W1)*(dcs2-dcs1)+dcs1;
	}
      }
    }

    for(int i=0;i<W_binnum;i++){
      sstr<<"W="<<W_adjust[i]<<endl;
      output<<sstr.str();
      sstr.str("");
      for(int j=0;j<ct_binnum;j++){
	output<<ct_adjust[i][j]<<"        "<<dcs_adjust[i][j]<<endl;
      }
      output<<endl;
    }

    input.close();
    output.close();
  }

{
    ifstream input("filename_KpSm");
    ofstream output("csKpSm_adjust");

    string str,strsub;
    int ch_first,ch_last;
    const int W_binnum=15;
    const int ct_binnum=50;
    double W[W_binnum],W_min=1739,W_max=2405,W_binwidth=(W_max-W_min)/(double)(W_binnum-1);

    double ct[W_binnum][ct_binnum],dcs[W_binnum][ct_binnum],ct_binwidth=2/(double)(ct_binnum-1);

    ostringstream sstr;
    for(int i=0;i<W_binnum;i++){
      getline(input,str);
      if (str.empty()) continue;
      ifstream infile(str.c_str());
      ch_first=str.find_last_of("m")+1;
      ch_last=str.find_last_of(".")-1;
      W[i]=atof(str.substr(ch_first,ch_last).c_str());
      for(int j=0;j<ct_binnum;j++){
	infile>>ct[i][j]>>dcs[i][j];
      }
    }

    double W_adjust[W_binnum],ct_adjust[W_binnum][ct_binnum],dcs_adjust[W_binnum][ct_binnum];
    for(int i=0;i<W_binnum;i++){
      W_adjust[i]=W_min+W_binwidth*i;
      for(int j=0;j<ct_binnum;j++){
	ct_adjust[i][j]=-1+ct_binwidth*j;
      }
    }

    double W1=0,W2=0,ct11=0,ct12=0,ct21=0,ct22=0,dcs11=0,dcs12=0,dcs21=0,dcs22=0,dcs1=0,dcs2=0;
    int ith;
    for(int u=0;u<W_binnum;u++){
      for(int i=0;i<W_binnum;i++){
	if(W_adjust[u]>=W[i] && W_adjust[u]<=W[i+1]){
	  W1=W[i];
	  W2=W[i+1];
	  ith=i;
	  break;
	}
      }
      for(int v=0;v<ct_binnum;v++){
	if(v==0){
	  dcs1=dcs[ith][0];
	  dcs2=dcs[ith+1][0];
	  dcs_adjust[u][v]=(W_adjust[u]-W1)/(W2-W1)*(dcs2-dcs1)+dcs1;
	}
	else{
	  for(int j=0;j<ct_binnum;j++){
	    if(ct_adjust[u][v]>=ct[ith][j] && ct_adjust[u][v]<=ct[ith][j+1]){
	      ct11=ct[ith][j];
	      ct12=ct[ith][j+1];
	      dcs11=dcs[ith][j];
	      dcs12=dcs[ith][j+1];
	      break;
	    }
	  }
	  for(int j=0;j<ct_binnum;j++){
	    if(ct_adjust[u][v]>=ct[ith+1][j] && ct_adjust[u][v]<=ct[ith+1][j+1]){
	      ct21=ct[ith+1][j];
	      ct22=ct[ith+1][j+1];
	      dcs21=dcs[ith+1][j];
	      dcs22=dcs[ith+1][j+1];
	      break;
	    }
	  }
	 
	  dcs1=(ct_adjust[u][v]-ct11)/(ct12-ct11)*(dcs12-dcs11)+dcs11;
	  dcs2=(ct_adjust[u][v]-ct21)/(ct22-ct21)*(dcs22-dcs21)+dcs21;
	  dcs_adjust[u][v]=(W_adjust[u]-W1)/(W2-W1)*(dcs2-dcs1)+dcs1;
	}
      }
    }

    for(int i=0;i<W_binnum;i++){
      sstr<<"W="<<W_adjust[i]<<endl;
      output<<sstr.str();
      sstr.str("");
      for(int j=0;j<ct_binnum;j++){
	output<<ct_adjust[i][j]<<"        "<<dcs_adjust[i][j]<<endl;
      }
      output<<endl;
    }

    input.close();
    output.close();
  }

{
    ifstream input("filename_KzSz");
    ofstream output("csKzSz_adjust");

    string str,strsub;
    int ch_first,ch_last;
    const int W_binnum=30;
    const int ct_binnum=50;
    double W[W_binnum],W_min=1694,W_max=2405,W_binwidth=(W_max-W_min)/(double)(W_binnum-1);

    double ct[W_binnum][ct_binnum],dcs[W_binnum][ct_binnum],ct_binwidth=2/(double)(ct_binnum-1);

    ostringstream sstr;
    for(int i=0;i<W_binnum;i++){
      getline(input,str);
      if (str.empty()) continue;
      ifstream infile(str.c_str());
      ch_first=str.find_last_of("z")+1;
      ch_last=str.find_last_of(".")-1;
      W[i]=atof(str.substr(ch_first,ch_last).c_str());
      for(int j=0;j<ct_binnum;j++){
	infile>>ct[i][j]>>dcs[i][j];
      }
    }

    double W_adjust[W_binnum],ct_adjust[W_binnum][ct_binnum],dcs_adjust[W_binnum][ct_binnum];
    for(int i=0;i<W_binnum;i++){
      W_adjust[i]=W_min+W_binwidth*i;
      for(int j=0;j<ct_binnum;j++){
	ct_adjust[i][j]=-1+ct_binwidth*j;
      }
    }

    double W1=0,W2=0,ct11=0,ct12=0,ct21=0,ct22=0,dcs11=0,dcs12=0,dcs21=0,dcs22=0,dcs1=0,dcs2=0;
    int ith;
    for(int u=0;u<W_binnum;u++){
      for(int i=0;i<W_binnum;i++){
	if(W_adjust[u]>=W[i] && W_adjust[u]<=W[i+1]){
	  W1=W[i];
	  W2=W[i+1];
	  ith=i;
	  break;
	}
      }
      for(int v=0;v<ct_binnum;v++){
	if(v==0){
	  dcs1=dcs[ith][0];
	  dcs2=dcs[ith+1][0];
	  dcs_adjust[u][v]=(W_adjust[u]-W1)/(W2-W1)*(dcs2-dcs1)+dcs1;
	}
	else{
	  for(int j=0;j<ct_binnum;j++){
	    if(ct_adjust[u][v]>=ct[ith][j] && ct_adjust[u][v]<=ct[ith][j+1]){
	      ct11=ct[ith][j];
	      ct12=ct[ith][j+1];
	      dcs11=dcs[ith][j];
	      dcs12=dcs[ith][j+1];
	      break;
	    }
	  }
	  for(int j=0;j<ct_binnum;j++){
	    if(ct_adjust[u][v]>=ct[ith+1][j] && ct_adjust[u][v]<=ct[ith+1][j+1]){
	      ct21=ct[ith+1][j];
	      ct22=ct[ith+1][j+1];
	      dcs21=dcs[ith+1][j];
	      dcs22=dcs[ith+1][j+1];
	      break;
	    }
	  }
	 
	  dcs1=(ct_adjust[u][v]-ct11)/(ct12-ct11)*(dcs12-dcs11)+dcs11;
	  dcs2=(ct_adjust[u][v]-ct21)/(ct22-ct21)*(dcs22-dcs21)+dcs21;
	  dcs_adjust[u][v]=(W_adjust[u]-W1)/(W2-W1)*(dcs2-dcs1)+dcs1;
	}
      }
    }

    for(int i=0;i<W_binnum;i++){
      sstr<<"W="<<W_adjust[i]<<endl;
      output<<sstr.str();
      sstr.str("");
      for(int j=0;j<ct_binnum;j++){
	output<<ct_adjust[i][j]<<"        "<<dcs_adjust[i][j]<<endl;
      }
      output<<endl;
    }

    input.close();
    output.close();
  }
}
