//Takes file list from RootList and runs the BOSWrite code to create a bos file from the root file

#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <iostream>
#include <iomanip>

using namespace std;
void automate() {
    string rootfilename;
    cout <<"Opening  file RootList"<<endl;
    
    ifstream filename;
    filename.open("RootList");
    if (!filename) {
        cerr << "Can't open file run_files"<< endl;
        exit(1);
    }
    ostringstream command;
   	while( filename >> rootfilename){
        command <<"~/Hyperon/jgen_v3/PutInBOS/build/bin/BOSwrite_v2.0 "<<rootfilename;
        system(command.str().c_str());
        command.clear();
        command.str("");
    }
    filename.close();
}