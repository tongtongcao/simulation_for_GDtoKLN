#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
using namespace std;

void Spin(){
  ostringstream sstr;
  sstr<<"php Spin.php";
  system(sstr.str().c_str());
  sstr.str("");

  ofstream output("Spin.txt");
  double E;
  ifstream input;
  string str;
  for (int i = 0; i < 86; i++) {
    E = 900+20*i;
    sstr<<E<<"S.txt";
    input.open(sstr.str().c_str());
    sstr.str("");

    sstr<<"E="<<E;
    output<<sstr.str().c_str()<<endl;
    sstr.str("");

    while(getline(input,str)){
      if(str.size()==54)
	output<<str<<endl;
    }
    output<<endl;
    input.close();

    sstr<<"rm "<<E<<"S.txt";
    system(sstr.str().c_str());
    sstr.str("");
  }
}
