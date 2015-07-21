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
#include "./reconneutronSelector/reconneutronSelector.h"

string Time = "June1";
TRint* theApp;

void Readdatafile(TChain *t,string dirname,int start,int end)
{ 
  string middle="outneutron";
  string last="_clustered.root";
  for(int i=start; i<=end; i++)
    {
      TString filename;
      if(i==0)
	filename.Form("%s%s",middle.c_str(),last.c_str());
      else 
	filename.Form("%s_v%d%s",middle.c_str(),i,last.c_str());	
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
  //#define use_TProof
#ifdef use_TProof
  chain->SetProof();
  TProof* pr = TProof::Open("workers=2");
  pr->SetParameter("PROOF_Packetizer","TPacketizer");
  pr->SetParameter("PROOF_MaxSlavesPerNode",4);
  //  chain->Process("./dstreeneutronSelector/dstreeneutronSelector.C+",option.Data(),-1,0);  
  chain->Process("./reconneutronSelector/reconneutronSelector.C+",option.Data(),-1,0);  
#else
  reconneutronSelector *selector = new reconneutronSelector();
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
  
  string dirname="/darkside/users/hqian/AmC13/MCNeutron/";
  
  TString label;
  label.Form("MCNeutron_%d_%s.root",end,Time.c_str());
  label.Prepend(dirname);

  TChain *t=new TChain("Recon");
  Readdatafile(t,dirname,start,end);
  cout<<"events= "<<t->GetEntries()<<endl;
  cout<<"Start Using the Selector..."<<endl;
  Process(t, label);
  
  std::cout << "==> Application finished." << std::endl;
  return 0; 
}
#endif /*__CINT__*/
