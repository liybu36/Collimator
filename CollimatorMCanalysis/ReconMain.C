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
#include "./ReconSelector/ReconSelector.h"

//const string Time = "_Oct13AM";
const int datafiles=7;
int Volume = 9;

TRint* theApp;

void Readdatafile(TChain *t, int volume)
{   string dirname="/darkside/users/hqian/neutron0G163tiltlow/cluster_data/";
  //   string dirname="/ds50/data/user/hqian36/Collimator/neutron0G163tilt/cluster_data/";
  string filename;
  string middle="outnvetoambe";
  string last="_cylinder_clustered.root";
  stringstream oss;
  stringstream oss2;
  for(int i=0; i<datafiles; i++)
    {
      if(i==0)
        {
          oss<<volume;
          filename=dirname+middle+oss.str()+"L"+last;
          //      filename<<dirname<<middle<<volume<<"L"<<last;
          //sprintf(filename,"%s%s%iL%s",dirname.c_str(),middle.c_str(),volume,last.c_str());
          t->Add(filename.c_str());
          cout<<"Processing Data file: "<<filename<<endl;
          oss.str("");
        }
      else
	{
	  oss<<volume;
	  oss2<<i;
	  filename=dirname+middle+oss.str()+"L_v"+oss2.str()+last;
	  t->Add(filename.c_str());
	  cout<<"Processing Data file: "<<filename<<endl;
	  oss.str("");
	  oss2.str("");
	}
    }
}

void Process(TChain* chain, TString label){
  TString option = label;
  chain->SetProof();
  TProof* pr = TProof::Open("workers=6");
  //Print more information
  //  pr->SetLogLevel(2,TProofDebug::kPacketizer);
  // pr->SetParameter("PROOF_Packetizer","TPacketizer");
  // pr->SetParameter("PROOF_MaxSlavesPerNode",8);

  //  ReconSelector *selector = new ReconSelector();
  //  Bool_t withfriends = kFALSE;
  //  TDSet *dataset = new TDSet(*chain, withfriends);
  // dataset->Add("/darkside/users/hqian/neutron0G163tilt/cluster_data/outnvetoambe9L_cylinder_clustered.root");
  // pr->ShowDataSets();
  // dataset->Process("./ReconSelector/ReconSelector.C++",option.Data(),100,100);  
  chain->Process("./ReconSelector/ReconSelector.C+",option.Data(),-1,0);  
  //  gSystem->Exit(0);
  //  chain->Process(pr,option.Data(),100,0);    
}

void ReconMain(){ 
  string label;
  if(Volume == 9) label="collimator";
  else label="no_collimator";

  TChain *t=new TChain("Recon");
  t->Add("/darkside/users/hqian/neutron0G163tilt/cluster_data/outnvetoambe9L_test_cylinder_clustered.root");
  //data for BG
  //  t->Add("/home/hqian/montecarlo/g4ds10/Linux-g++/neutron0G163tilt/cluster_data/outnvetoambe0L_v5_cylinder_clustered.root");
  cout<<"Start Using the Selector..."<<endl;
  Process(t, label);
  cout<<"Successfully finished the app."<<endl;
   
}

#ifndef __CINT__
//main function
int main(){
  theApp = new TRint("App",NULL,NULL);
  //  theApp->Connect("keypressed(Int_t)","TSystem",gSystem,"ExitLoop()");
  string label;
  if(Volume == 9) label="collimator";
  else label="no_collimator";

  TChain *t=new TChain("Recon");
  //  Readdatafile(t,Volume);
  //data for collimator
  //  t->Add("/darkside/users/hqian/neutron0G163tiltlow/cluster_data/outnvetoambe9L_test_cylinder_clustered.root");
  t->Add("/darkside/users/hqian/neutron0G163tiltlow/cluster_data/outneutronsiso9_cylinder_clustered.root");
  //data for BG
  //  t->Add("/home/hqian/montecarlo/g4ds10/Linux-g++/neutron0G163tilt/cluster_data/outnvetoambe0L_v5_cylinder_clustered.root");
  cout<<"Start Using the Selector..."<<endl;
  Process(t, label);
  cout<<"Successfully finished the app."<<endl;
  return 1;

}
#endif /*__CINT__*/
