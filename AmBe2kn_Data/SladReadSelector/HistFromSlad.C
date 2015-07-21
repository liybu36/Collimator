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
#include "TProofMgr.h"
#include "TProofLog.h"
//#include "TROOT.h"
#include "TDSet.h"

using namespace std;

#include "./SladReadSelector/SladReadSelector.h"

TRint* theApp;
void Process(TChain* chain, TString fHistOutName);

void HistFromSlad(TString fInNameEvent) {
  gROOT->SetBatch(kTRUE);

  TChain* chain = new TChain("events");
  chain->Add(fInNameEvent.Data());

  TChain *chainFriend(NULL);
  TString fInNameXY = fInNameEvent;
  fInNameXY.Replace(fInNameXY.Length()-5, 5, "_masas_xy.root"); //replace .root w/ _masas_xy.root
  chainFriend = new TChain("masas_xy"); chainFriend->Add(fInNameXY.Data());
  chain->AddFriend(chainFriend);
  TString fInNameJasonXY = fInNameEvent;
  fInNameJasonXY.Replace(fInNameJasonXY.Length()-5, 5, "_xylocator_xy.root"); //replace .root w/ _xylocator_xy.root
  chainFriend = new TChain("xylocator_xy"); chainFriend->Add(fInNameJasonXY.Data());
  chain->AddFriend(chainFriend);
  TString fInNameS2 = fInNameEvent;
  fInNameS2.Replace(fInNameS2.Length()-5, 5, "_s2.root"); //replace .root w/ _s2.root
  chainFriend = new TChain("s2_fraction"); chainFriend->Add(fInNameS2.Data());
  chain->AddFriend(chainFriend);
  TString fInNameAllPulses = fInNameEvent;
  fInNameAllPulses.Replace(fInNameAllPulses.Length()-5, 5, "_allpulses.root"); //replace .root w/ _s2.root
  chainFriend = new TChain("pulse_info"); chainFriend->Add(fInNameAllPulses.Data());
  chain->AddFriend(chainFriend);

  TString fHistOutName(fInNameEvent);
  fHistOutName.Replace(fHistOutName.Length()-5, 5, "_Hist.root"); //replace .root w/ _s2.root

  Process(chain, fHistOutName);

  delete chain;

}

//#define PROOF

void Process(TChain* chain, TString fHistOutName){
  std::cout << "outhistfile: " << fHistOutName.Data() << std::endl;

  TString home(gSystem->Getenv("PWD"));
  TString option = fHistOutName;

#ifdef PROOF
//  TProof::Open("");
  TProof* pr = TProof::Open("workers=4");
  pr->Exec(Form("gSystem->Load(\"%s/SladDict.so\")", home.Data()));

  // tell the chain that we want to use PROOF
  chain->SetProof();
#endif

#if 0 // Compile
//  chain->Process("./SladReadSelector/SladReadSelector.C+", option.Data(), -1, 0);
  chain->Process("./SladReadSelector/SladReadSelector.C+", option.Data(), 400, 0);
#else
  SladReadSelector *selector = new SladReadSelector();
//  selector->SetDriftFieldConf(200);

#if 1 // not Debug
  chain->Process(selector, option.Data());
//  chain->Process(selector, option.Data(), 50000, 0);//nEvents,0);
#else
#ifndef PROOF
  chain->Process(selector, option.Data(), 400,0);
#else

  Bool_t withfriends = kTRUE;
  TDSet *dataset = new TDSet(*chain, withfriends); //TDSet(const TChain& chain, Bool_t withfriends = kTRUE)

  dataset->Process(selector, option.Data());//nEvents,0);
  TProofMgr* mgr = pr->GetManager() ;
  // save log
  mgr->GetSessionLogs()->Save("*", "log_all-workers.txt") ;

#endif // PROOF

#endif // Debug
#endif // Compile
//  delete selector; // somehow this line causes seg fal.

  std::cout << "outhistfile: " << fHistOutName.Data() << std::endl;

}

#ifndef __CINT__
int main(int argc, char **argv) {
	theApp = new TRint("App", &argc, argv, NULL, 0);
	theApp->Connect("KeyPressed(Int_t)","TSystem",gSystem,"ExitLoop()");
	if ( theApp->Argc() == 2 ) {
		std::cout << "\n==========> HistFromSlad <=============" << std::endl;
	        HistFromSlad(theApp->Argv(1));

	} else {
		std::cout << "Usage:" <<theApp->Argc()<< std::endl;
                std::cout << "./HistFromSlad fInNameEvent.root" << std::endl;
		return 0;
	}

	std::cout << "==> Application finished." << std::endl;
	return 0;
}
#endif /* __CINT __ */
