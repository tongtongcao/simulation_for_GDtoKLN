#include <iostream>  
#include <fstream>
#include <sstream>
#include <string> 

using namespace std;
void convert_original(){

  {
    ifstream input("filename_KL");
    ofstream output("csKL_orignial");

    string str,strsub;
    int ch_first,ch_last;
    double W;

    ostringstream sstr;

    double ct[50],dcs[50];

    while(getline(input,str)) {
      if (str.empty()) continue;
      ifstream infile(str.c_str());
      cout<<str.c_str()<<endl;
      ch_first=str.find_last_of("L")+1;
      ch_last=str.find_last_of(".")-1;
      W=atof(str.substr(ch_first,ch_last).c_str());
      sstr<<"W="<<W<<endl;
      output<<sstr.str();
      sstr.str("");
      for(int i=0;i<50;i++){
	infile>>ct[i]>>dcs[i];
	output<<ct[i]<<"        "<<dcs[i]<<endl;
      }
      output<<endl;
      infile.close();
    }

    input.close();
    output.close();
  }

  {
    ifstream input("filename_KpSp");
    ofstream output("csKpSp_orignial");

    string str,strsub;
    int ch_first,ch_last;
    double W;

    ostringstream sstr;

    double ct[50],dcs[50];

    while(getline(input,str)) {
      if (str.empty()) continue;
      ifstream infile(str.c_str());
      cout<<str.c_str()<<endl;
      ch_first=str.find_last_of("p")+1;
      ch_last=str.find_last_of(".")-1;
      W=atof(str.substr(ch_first,ch_last).c_str());
      sstr<<"W="<<W<<endl;
      output<<sstr.str();
      sstr.str("");
      for(int i=0;i<50;i++){
	infile>>ct[i]>>dcs[i];
	output<<ct[i]<<"        "<<dcs[i]<<endl;
      }
      output<<endl;
      infile.close();
    }

    input.close();
    output.close();
  }

  {
    ifstream input("filename_KpSm");
    ofstream output("csKpSm_orignial");

    string str,strsub;
    int ch_first,ch_last;
    double W;

    ostringstream sstr;

    double ct[50],dcs[50];

    while(getline(input,str)) {
      if (str.empty()) continue;
      ifstream infile(str.c_str());
      cout<<str.c_str()<<endl;
      ch_first=str.find_last_of("m")+1;
      ch_last=str.find_last_of(".")-1;
      W=atof(str.substr(ch_first,ch_last).c_str());
      sstr<<"W="<<W<<endl;
      output<<sstr.str();
      sstr.str("");
      for(int i=0;i<50;i++){
	infile>>ct[i]>>dcs[i];
	output<<ct[i]<<"        "<<dcs[i]<<endl;
      }
      output<<endl;
      infile.close();
    }

    input.close();
    output.close();
  }

  {
    ifstream input("filename_KzSz");
    ofstream output("csKzSz_orignial");

    string str,strsub;
    int ch_first,ch_last;
    double W;

    ostringstream sstr;

    double ct[50],dcs[50];

    while(getline(input,str)) {
      if (str.empty()) continue;
      ifstream infile(str.c_str());
      cout<<str.c_str()<<endl;
      ch_first=str.find_last_of("z")+1;
      ch_last=str.find_last_of(".")-1;
      W=atof(str.substr(ch_first,ch_last).c_str());
      sstr<<"W="<<W<<endl;
      output<<sstr.str();
      sstr.str("");
      for(int i=0;i<50;i++){
	infile>>ct[i]>>dcs[i];
	output<<ct[i]<<"        "<<dcs[i]<<endl;
      }
      output<<endl;
      infile.close();
    }

    input.close();
    output.close();
  }
}
