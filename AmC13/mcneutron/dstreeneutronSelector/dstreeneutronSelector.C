#define dstreeneutronSelector_cxx
// The class definition in dstreeneutronSelector.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("dstreeneutronSelector.C")
// Root > T->Process("dstreeneutronSelector.C","some options")
// Root > T->Process("dstreeneutronSelector.C+")
//

#include "dstreeneutronSelector.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>

#include <TH2.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TBox.h>
#include <TStyle.h>

using namespace std;

void dstreeneutronSelector::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

void dstreeneutronSelector::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

  Info("SlaveBegin()","beginning .....");
  TString option = GetOption();
  BookHistograms();

}

Bool_t dstreeneutronSelector::Process(Long64_t entry)
{
  fChain->GetEntry(entry);
  FillHistograms();

  return kTRUE;
}

void dstreeneutronSelector::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void dstreeneutronSelector::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

  string label = GetOption();
  TList* list = GetOutputList();
  id = dynamic_cast<TH1F*>(list->FindObject(Form("id")));


  string output = label;
  TFile outfile(output.c_str(), "RECREATE");
  id->Write();

  outfile.Write();
  outfile.Close();

  Info("Terminate()","terminating ...");
  
}

void dstreeneutronSelector::FillHistograms()
{



}

void dstreeneutronSelector::BookHistograms()
{
  string label = GetOption();
  TList* list = GetOutputList();

  
}

