#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include "TGraph.h"

using namespace std;

int Am241_PureC13RawData()
{
  string indir = "/darkside/users/hqian/AmC13/";
  string infile = indir + "Am241_Pure13CAlphaNSpec.root";
  TFile *indata = new TFile(infile.c_str());
  vector<string> isotopes = Isotopes();
  vector<TGraph*> graph;
  for(size_t i=0; i<isotopes.size(); i++)
    {
      graph.push_back((TGraph*) indata->Get(isotopes.at(i).c_str()));
      if(GetPoints(isotopes.at(i),graph.at(i),indir))
	cout<<isotopes.at(i)<<" is saved!"<<endl;
    }
  return 1;
  //  TGraph *Am241_graph = (TGraph*) indata->Get("Am241");

}

int GetPoints(string name,TGraph* Am241_graph,string indir)
  {  
  Int_t N = Am241_graph->GetN(); 
  Double_t* Energy = Am241_graph->GetX();
  Double_t* Probability = Am241_graph->GetY();

  string outdir = indir;
  string outdataname = outdir+name+"C13.txt";
  ofstream outdata (outdataname.c_str());
  if(outdata.is_open())
    {
      outdata<<"Energy\t"<<"Neutron Yield"<<"\n";
      for(int i=0; i<N; i++)
	outdata<<*(Energy+i)<<"\t"<<*(Probability+i)<<"\n";
    }
  else {
    cout<<"Cannot open the file"<<endl;
    return 0;
  }
  outdata.close();
  
  return 1;
}

vector<string> Isotopes()
{
  vector<string> isotopes;
  isotopes.push_back("Total");
  isotopes.push_back("Am241");
  
  return isotopes;
}
