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

using namespace std;
#include "./dstreeneutronSelector/dstreeneutronSelector.h"

string Time = "_June1";
TRint* theApp;

void Readdatafile(TChain *t,string dirname,int start,int end)
{ 
  string middle="outnvetoambe";
  string last=".root";
  for(int i=start; i<=end; i++)
    {
      TString filename;
      filename.Form("%s%06d%s",middle.c_str(),i,last.c_str());
      filename.Prepend(dirname.c_str());
      ifstream NameCheck;
      NameCheck.open(filename.Data());
      if(!NameCheck.good())
	continue;
      else{
	TFile *f = new TFile(filename);
	if(f->IsZombie())
	  continue;
	else{
	  t->Add(filename);
	  cout<<"Processing Data file: "<<filename<<endl;
	}
      }
    }

}

void Process(TChain* chain, TString label){
  TString option = label;
#define use_TProof
#ifdef use_TProof
  chain->SetProof();
  TProof* pr = TProof::Open("workers=9");
  pr->SetParameter("PROOF_Packetizer","TPacketizer");
  pr->SetParameter("PROOF_MaxSlavesPerNode",8);
  chain->Process("./dstreeneutronSelector/dstreeneutronSelector.C+",option.Data(),-1,0);  
#else
  dstreeneutronSelector *selector = new dstreeneutronSelector();
  selector->SetTRint(theApp);
  chain->Process(selector,option.Data());
  gSystem->Exit(0);
#endif

}

#ifndef __CINT__
int main(int argc, char** argv){
  theApp = new TRint("App",&argc,argv,NULL,0);
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
  
  string dirname="/darkside/users/hqian/neutron0G163tiltlow/cluster_data/";
  // string dirname="/ds50/data/user/hqian36/Collimator/neutron0G163tilt/cluster_data/";
  
  string label;
  label=dirname+Time;

  TChain *t=new TChain("Recon");
  Readdatafile(t,dirname,start,end);
  cout<<"Start Using the Selector..."<<endl;
  Process(t, label);
  
  std::cout << "==> Application finished." << std::endl;
  return 0; 
}
#endif /*__CINT__*/
