void test1(){
  ifstream input("csKL_adjust");

  vector<double> vect_dcs;
  vector<string> table;
  string line,tempstring,subline,subline1,subline2,subline3,subline4;
  int tindex1,tindex2,tindex3,tindex4;
  double dcs;
  int i=0;

  while(getline(input,line)) {
    if (line.empty() || line.find("W")!=string::npos) continue;
    table.push_back(line);

    tindex1 = line.find_first_not_of(" ");
    subline=line.substr(tindex1);
    tindex2 = subline.find_first_of(" ");
    subline1 = subline.substr(0, tindex2);
     
    subline2 = subline.substr(tindex2);
    tindex3 = subline2.find_first_not_of(" ");
    subline3 = subline2.substr(tindex3);
    tindex4 = subline3.find_first_of(" ");
    subline4 = subline3.substr(0, tindex4);
    istringstream(subline4)>>dcs;
    vect_dcs.push_back(dcs);
    cout<<i+1<<"  "<<dcs<<endl;
    i++;
    if(i%50==0) i=0;
  }

  double temp=0;
  for(int i=0;i<vect_dcs.size();i++){
    if(temp<vect_dcs[i]) temp=vect_dcs[i];
  }
  cout<<temp<<endl;
}
