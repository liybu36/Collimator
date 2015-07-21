#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include "TApplication.h"
#include "TRint.h"
#include "TProof.h"
#include "TDSet.h"
#include "TROOT.h"
#include "reconneutronClass.h"

using namespace std;
string Time = "June10";
TRint* theApp;

void reconneutronMain(int start, int end)
{
  reconneutronClass *r = new reconneutronClass;
  string indir="/darkside/users/hqian/AmC13/MCNeutron/";
  string outdir = indir;
  string inmiddle = "outneutron";
  string intail = "clustered.root";
  TString outfile;
  outfile.Form("%s_%d%d_%s.root",inmiddle.c_str(),start,end,Time.c_str());
  cout<<indir<<endl;
  cout<<outdir<<outfile<<endl;
  r->SetInPath(indir);
  r->SetOutPath(outdir);
  r->SetInMiddle(inmiddle);  
  r->SetOutFile(outfile);
  TChain *t=new TChain("Recon");
  r->ReadDataFile(t,start,end,intail);
  cout<<"events= "<<t->GetEntries()<<endl;  
  r->SetTChain(t);
  r->Init();
  r->BookHistograms();
  r->FillHistograms();
  r->SaveHistograms();
  std::cout << "==> Saved Plots." << std::endl;
  
}

#ifndef __CINT__
int main(int argc, char** argv){
  theApp = new TRint("App",&argc,argv,NULL,0);
  //  gROOT->ProcessLine("#pragma link C++ struct Event+;");
  //  gROOT->ProcessLine("#pragma link C++ struct Plots+;");
  //  gROOT->ProcessLine("#pragma link C++ class reconneutronClass+;");
  //  gROOT->ProcessLine("#pragma link C++ class std::vector<Plots>+;");
  int start, end;
  if(theApp->Argc() == 2)
    {
      start = atoi(theApp->Argv(1));
      end = start;
    }
  else if(theApp->Argc() == 3)
    {
      start = atoi(theApp->Argv(1));
      end = atoi(theApp->Argv(2));
    }
  else{
    cout<<"Usage: ./SladOD startfile endfile"<<endl;
    cout<<"Usage: ./SladOD startfile"<<endl;
  }  
  reconneutronMain(start,end);
  
  std::cout << "==> Application finished." << std::endl;
  return 0; 
}
#endif /*__CINT__*/
