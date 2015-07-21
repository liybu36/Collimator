#define ReconSelector_cxx
// The class definition in ReconSelector.h has been generated automatically
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
// Root > T->Process("ReconSelector.C")
// Root > T->Process("ReconSelector.C","some options")
// Root > T->Process("ReconSelector.C+")
//
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "ReconSelector.h"
#include <TH2.h>
#include <TStyle.h>
#include <TH2F.h>
#include <TCanvas.h>

using namespace std;

const string Time = "_Oct22AM";
//const int datafiles=9;
//int Volume = 9;

void ReconSelector::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

void ReconSelector::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).
  Info("ReconSelector::SlaveBegin()","beginning .....");
  TString option = GetOption();
  std::cout<<"GetOption()=  "<<option.Data()<<std::endl;
  BookHistograms();
  
}

Bool_t ReconSelector::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either ReconSelector::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.
  Int_t chainentry = fChain->GetChainEntryNumber(entry);
  if(chainentry%(fChain->GetEntries()/10)==0) printf("Processing Entry number %ld [%ld %% of %lld]\n", (long int)chainentry, (long int)(chainentry/(fChain->GetEntries()/100)), fChain->GetEntries());
  fChain->GetEntry(entry);
  ++fNumberOfEvents;
  FillHistograms();
   return kTRUE;
}

void ReconSelector::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void ReconSelector::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
  /*  string label;
  if(Volume == 9) label="collimator";
  else label="no_collimator";  
  */
  string label = GetOption();
  cout<<"fNumberOfEvents= "<<fNumberOfEvents<<endl;
  Info("ReconSelector::Terminate()","Closing ....");
    string outdirname="/darkside/users/hqian/neutron0G163tiltlow/cluster_data/";
  //  string outdirname="/ds50/data/user/hqian36/Collimator/neutron0G163tilt/cluster_data/";
  string outputname=outdirname+"nuclear0G163tilt_recoil_"+label+Time+"_clustered"+".root";

  TList *list = GetOutputList();

  TCanvas *c6=new TCanvas("c6",label.c_str(),1000,600);
  c6->SetLogy();
  QuenchingFactor=dynamic_cast<TH1F*>(list->FindObject(Form("QuenchingFactor")));
  QuenchingFactor->Draw();
  c6->SaveAs(Form("ambe_QuenchingFactor_%s%s.png",label.c_str(),Time.c_str()));

  Hist1DTPCedep=dynamic_cast<TH1F*>(list->FindObject(Form("Hist1DTPCedep")));
  Hist1DTPCedep->Draw();
  c6->SaveAs(Form("ambe_Hist1DTPCedep_%s%s.png",label.c_str(),Time.c_str()));

  Hist1DTPCedep_nuclear=dynamic_cast<TH1F*>(list->FindObject(Form("Hist1DTPCedep_nuclear")));
  Hist1DTPCedep_nuclear->Draw();
  c6->SaveAs(Form("ambe_Hist1DTPCedep_nuclear_%s%s.png",label.c_str(),Time.c_str()));

  Hist1DTPCedep_electron=dynamic_cast<TH1F*>(list->FindObject(Form("Hist1DTPCedep_electron")));
  Hist1DTPCedep_electron->Draw();
  c6->SaveAs(Form("ambe_Hist1DTPCedep_electron_%s%s.png",label.c_str(),Time.c_str()));

  Hist1DNVedep=dynamic_cast<TH1F*>(list->FindObject(Form("Hist1DNVedep")));
  Hist1DNVedep->Draw();
  c6->SaveAs(Form("ambe_Hist1DNVedep_%s%s.png",label.c_str(),Time.c_str()));

  Hist1DNVedep_nuclear=dynamic_cast<TH1F*>(list->FindObject(Form("Hist1DNVedep_nuclear")));
  Hist1DNVedep_nuclear->Draw();
  c6->SaveAs(Form("ambe_Hist1DNVedep_nuclear_%s%s.png",label.c_str(),Time.c_str()));

  Hist1DNVedep_electron=dynamic_cast<TH1F*>(list->FindObject(Form("Hist1DNVedep_electron")));
  Hist1DNVedep_electron->Draw();
  c6->SaveAs(Form("ambe_Hist1DNVedep_electron_%s%s.png",label.c_str(),Time.c_str()));

  TCanvas *c1=new TCanvas("c1",label.c_str(),1000,600);
  Hist2D=dynamic_cast<TH2F*>(list->FindObject(Form("Hist2D")));
  Hist2D->Draw("colz");
  c1->SaveAs(Form("ambe_Hist2D_%s%s.png",label.c_str(),Time.c_str()));
 
  Hist2DYZ=dynamic_cast<TH2F*>(list->FindObject(Form("Hist2DYZ")));
  Hist2DYZ->Draw("colz");
  c1->SaveAs(Form("ambe_Hist2DYZ_%s%s.png",label.c_str(),Time.c_str()));

  Hist2DXZ=dynamic_cast<TH2F*>(list->FindObject(Form("Hist2DXZ")));
  Hist2DXZ->Draw("colz");
  c1->SaveAs(Form("ambe_Hist2DXZ_%s%s.png",label.c_str(),Time.c_str()));

  Hist2DXYslice1=dynamic_cast<TH2F*>(list->FindObject(Form("Hist2DXYslice1")));
  Hist2DXYslice1->Draw("colz");
  c1->SaveAs(Form("ambe_Hist2DXYslice1_%s%s.png",label.c_str(),Time.c_str()));

  Hist2DXYslice2=dynamic_cast<TH2F*>(list->FindObject(Form("Hist2DXYslice2")));
  Hist2DXYslice2->Draw("colz");
  c1->SaveAs(Form("ambe_Hist2DXYslice2_%s%s.png",label.c_str(),Time.c_str()));

  Hist2DXYslice3=dynamic_cast<TH2F*>(list->FindObject(Form("Hist2DXYslice3")));
  Hist2DXYslice3->Draw("colz");
  c1->SaveAs(Form("ambe_Hist2DXYslice3_%s%s.png",label.c_str(),Time.c_str()));

  Hist2DGap=dynamic_cast<TH2F*>(list->FindObject(Form("Hist2DGap")));
  c1->SetLogz();
  Hist2DGap->Draw("colz");
  c1->SaveAs(Form("ambe_Hist2DGap_%s%s.png",label.c_str(),Time.c_str()));
  
  TCanvas *c2=new TCanvas("c2",label.c_str(),1000,600);
  Hist2Dnuclear=dynamic_cast<TH2F*>(list->FindObject(Form("Hist2Dnuclear")));
  Hist2Dnuclear->Draw("colz");
  c2->SaveAs(Form("ambe_Hist2Dnuclear_%s%s.png",label.c_str(),Time.c_str()));

  Hist2DYZnuclear=dynamic_cast<TH2F*>(list->FindObject(Form("Hist2DYZnuclear")));
  Hist2DYZnuclear->Draw("colz");
  c2->SaveAs(Form("ambe_Hist2DYZnuclear_%s%s.png",label.c_str(),Time.c_str()));

  Hist2DXZnuclear=dynamic_cast<TH2F*>(list->FindObject(Form("Hist2DXZnuclear")));
  Hist2DXZnuclear->Draw("colz");
  c2->SaveAs(Form("ambe_Hist2DXZnuclear_%s%s.png",label.c_str(),Time.c_str()));

  Hist2DXYnuclearslice1=dynamic_cast<TH2F*>(list->FindObject(Form("Hist2DXYnuclearslice1")));
  Hist2DXYnuclearslice1->Draw("colz");
  c2->SaveAs(Form("ambe_Hist2DXYnuclearslice1_%s%s.png",label.c_str(),Time.c_str()));

  Hist2DXYnuclearslice2=dynamic_cast<TH2F*>(list->FindObject(Form("Hist2DXYnuclearslice2")));
  Hist2DXYnuclearslice2->Draw("colz");
  c2->SaveAs(Form("ambe_Hist2DXYnuclearslice2_%s%s.png",label.c_str(),Time.c_str()));

  Hist2DXYnuclearslice3=dynamic_cast<TH2F*>(list->FindObject(Form("Hist2DXYnuclearslice3")));
  Hist2DXYnuclearslice3->Draw("colz");
  c2->SaveAs(Form("ambe_Hist2DXYnuclearslice3_%s%s.png",label.c_str(),Time.c_str()));

  TCanvas *c4=new TCanvas("c4",label.c_str(),1000,600);
  Hist2Dsinglenuclear=dynamic_cast<TH2F*>(list->FindObject(Form("Hist2Dsinglenuclear")));
  Hist2Dsinglenuclear->Draw("colz");
  c4->SaveAs(Form("ambe_Hist2Dsinglenuclear_%s%s.png",label.c_str(),Time.c_str()));

  Hist2DXYsinglenuclearslice1=dynamic_cast<TH2F*>(list->FindObject(Form("Hist2DXYsinglenuclearslice1")));
  Hist2DXYsinglenuclearslice1->Draw("colz");
  c4->SaveAs(Form("ambe_Hist2DXYsinglenuclearslice1_%s%s.png",label.c_str(),Time.c_str()));

  Hist2DXYsinglenuclearslice2=dynamic_cast<TH2F*>(list->FindObject(Form("Hist2DXYsinglenuclearslice2")));
  Hist2DXYsinglenuclearslice2->Draw("colz");
  c4->SaveAs(Form("ambe_Hist2DXYsinglenuclearslice2_%s%s.png",label.c_str(),Time.c_str()));

  Hist2DXYsinglenuclearslice3=dynamic_cast<TH2F*>(list->FindObject(Form("Hist2DXYsinglenuclearslice3")));
  Hist2DXYsinglenuclearslice3->Draw("colz");
  c4->SaveAs(Form("ambe_Hist2DXYsinglenuclearslice3_%s%s.png",label.c_str(),Time.c_str()));

  Hist2DYZsinglenuclear=dynamic_cast<TH2F*>(list->FindObject(Form("Hist2DYZsinglenuclear")));
  Hist2DYZsinglenuclear->Draw("colz");
  c4->SaveAs(Form("ambe_Hist2DYZsinglenuclear_%s%s.png",label.c_str(),Time.c_str()));

  TCanvas *c5=new TCanvas("c5",label.c_str(),1000,600);
  Hist2Dsingle=dynamic_cast<TH2F*>(list->FindObject(Form("Hist2Dsingle")));
  Hist2Dsingle->Draw("colz");
  c5->SaveAs(Form("ambe_Hist2Dsingle_%s%s.png",label.c_str(),Time.c_str()));

  Hist2DXYsingleslice1=dynamic_cast<TH2F*>(list->FindObject(Form("Hist2DXYsingleslice1")));
  Hist2DXYsingleslice1->Draw("colz");
  c5->SaveAs(Form("ambe_Hist2DXYsingleslice1_%s%s.png",label.c_str(),Time.c_str()));

  Hist2DXYsingleslice2=dynamic_cast<TH2F*>(list->FindObject(Form("Hist2DXYsingleslice2")));
  Hist2DXYsingleslice2->Draw("colz");
  c5->SaveAs(Form("ambe_Hist2DXYsingleslice2_%s%s.png",label.c_str(),Time.c_str()));

  Hist2DXYsingleslice3=dynamic_cast<TH2F*>(list->FindObject(Form("Hist2DXYsingleslice3")));
  Hist2DXYsingleslice3->Draw("colz");
  c5->SaveAs(Form("ambe_Hist2DXYsingleslice3_%s%s.png",label.c_str(),Time.c_str()));

  Hist2DYZsingle=dynamic_cast<TH2F*>(list->FindObject(Form("Hist2DYZsingle")));
  Hist2DYZsingle->Draw("colz");
  c5->SaveAs(Form("ambe_Hist2DYZsingle_%s%s.png",label.c_str(),Time.c_str()));
  
  //generate the projection plots
  float leftend = -32.5;
  int leftendbin = Hist2DGap->GetXaxis()->FindBin(leftend);
  TH1D* py1 = Hist2DGap->ProjectionY("py1",leftendbin,leftendbin);

  int slicebin = Hist2DXYslice1->GetNbinsX();
  TH1D* Hist2DXYslice1ProjectionY = Hist2DXYslice1->ProjectionY("Hist2DXYslice1ProjectionY",0,slicebin);
  TH1D* Hist2DXYnuclearslice1ProjectionY = Hist2DXYnuclearslice1->ProjectionY("Hist2DXYnuclearslice1ProjectionY",0,slicebin);
  TH1D* Hist2DXYsingleslice1ProjectionY = Hist2DXYsingleslice1->ProjectionY("Hist2DXYsingleslice1ProjectionY",0,slicebin);
  TH1D* Hist2DXYsinglenuclearslice1ProjectionY = Hist2DXYsinglenuclearslice1->ProjectionY("Hist2DXYsinglenuclearslice1ProjectionY",0,slicebin);

  TH1D* Hist2DXYslice2ProjectionY = Hist2DXYslice2->ProjectionY("Hist2DXYslice2ProjectionY",0,slicebin);
  TH1D* Hist2DXYnuclearslice2ProjectionY = Hist2DXYnuclearslice2->ProjectionY("Hist2DXYnuclearslice2ProjectionY",0,slicebin);
  TH1D* Hist2DXYsingleslice2ProjectionY = Hist2DXYsingleslice2->ProjectionY("Hist2DXYsingleslice2ProjectionY",0,slicebin);
  TH1D* Hist2DXYsinglenuclearslice2ProjectionY = Hist2DXYsinglenuclearslice2->ProjectionY("Hist2DXYsinglenuclearslice2ProjectionY",0,slicebin);

  TH1D* Hist2DXYslice3ProjectionY = Hist2DXYslice3->ProjectionY("Hist2DXYslice3ProjectionY",0,slicebin);
  TH1D* Hist2DXYnuclearslice3ProjectionY = Hist2DXYnuclearslice3->ProjectionY("Hist2DXYnuclearslice3ProjectionY",0,slicebin);
  TH1D* Hist2DXYsingleslice3ProjectionY = Hist2DXYsingleslice3->ProjectionY("Hist2DXYsingleslice3ProjectionY",0,slicebin);
  TH1D* Hist2DXYsinglenuclearslice3ProjectionY = Hist2DXYsinglenuclearslice3->ProjectionY("Hist2DXYsinglenuclearslice3ProjectionY",0,slicebin);

  int slicebinX1 = Hist2DXYslice1->GetYaxis()->FindBin(-2);
  int slicebinX2 = Hist2DXYslice1->GetYaxis()->FindBin(2);
  TH1D* Hist2DXYslice1ProjectionX = Hist2DXYslice1->ProjectionX("Hist2DXYslice1ProjectionX",slicebinX1,slicebinX2);
  TH1D* Hist2DXYnuclearslice1ProjectionX = Hist2DXYnuclearslice1->ProjectionX("Hist2DXYnuclearslice1ProjectionX",slicebinX1,slicebinX2);
  TH1D* Hist2DXYsingleslice1ProjectionX = Hist2DXYsingleslice1->ProjectionX("Hist2DXYsingleslice1ProjectionX",slicebinX1,slicebinX2);
  TH1D* Hist2DXYsinglenuclearslice1ProjectionX = Hist2DXYsinglenuclearslice1->ProjectionX("Hist2DXYsinglenuclearslice1ProjectionX",slicebinX1,slicebinX2);

  TH1D* Hist2DXYslice2ProjectionX = Hist2DXYslice2->ProjectionX("Hist2DXYslice2ProjectionX",slicebinX1,slicebinX2);
  TH1D* Hist2DXYnuclearslice2ProjectionX = Hist2DXYnuclearslice2->ProjectionX("Hist2DXYnuclearslice2ProjectionX",slicebinX1,slicebinX2);
  TH1D* Hist2DXYsingleslice2ProjectionX = Hist2DXYsingleslice2->ProjectionX("Hist2DXYsingleslice2ProjectionX",slicebinX1,slicebinX2);
  TH1D* Hist2DXYsinglenuclearslice2ProjectionX = Hist2DXYsinglenuclearslice2->ProjectionX("Hist2DXYsinglenuclearslice2ProjectionX",slicebinX1,slicebinX2);

  TH1D* Hist2DXYslice3ProjectionX = Hist2DXYslice3->ProjectionX("Hist2DXYslice3ProjectionX",slicebinX1,slicebinX2);
  TH1D* Hist2DXYnuclearslice3ProjectionX = Hist2DXYnuclearslice3->ProjectionX("Hist2DXYnuclearslice3ProjectionX",slicebinX1,slicebinX2);
  TH1D* Hist2DXYsingleslice3ProjectionX = Hist2DXYsingleslice3->ProjectionX("Hist2DXYsingleslice3ProjectionX",slicebinX1,slicebinX2);
  TH1D* Hist2DXYsinglenuclearslice3ProjectionX = Hist2DXYsinglenuclearslice3->ProjectionX("Hist2DXYsinglenuclearslice3ProjectionX",slicebinX1,slicebinX2);
  
  //Draw the Projection histograms
  TCanvas *c11=new TCanvas("c11",label.c_str(),1000,600);
  Hist2DXYslice1ProjectionY->Draw();
  c11->SaveAs(Form("ambe_Hist2DXYslice1ProjectionY_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYslice2ProjectionY->Draw();
  c11->SaveAs(Form("ambe_Hist2DXYslice2ProjectionY_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYslice3ProjectionY->Draw();
  c11->SaveAs(Form("ambe_Hist2DXYslice3ProjectionY_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYslice1ProjectionX->Draw();
  c11->SaveAs(Form("ambe_Hist2DXYslice1ProjectionX_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYslice2ProjectionX->Draw();
  c11->SaveAs(Form("ambe_Hist2DXYslice2ProjectionX_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYslice3ProjectionX->Draw();
  c11->SaveAs(Form("ambe_Hist2DXYslice3ProjectionX_%s%s.png",label.c_str(),Time.c_str()));
  py1->Draw();
  c11->SaveAs(Form("ambe_Hist2DGapProjectionY_%s%s.png",label.c_str(),Time.c_str()));

  TCanvas *c12=new TCanvas("c12",label.c_str(),1000,600);
  Hist2DXYnuclearslice1ProjectionY->Draw();
  c12->SaveAs(Form("ambe_Hist2DXYnuclearslice1ProjectionY_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYnuclearslice2ProjectionY->Draw();
  c12->SaveAs(Form("ambe_Hist2DXYnuclearslice2ProjectionY_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYnuclearslice3ProjectionY->Draw();
  c12->SaveAs(Form("ambe_Hist2DXYnuclearslice3ProjectionY_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYnuclearslice1ProjectionX->Draw();
  c12->SaveAs(Form("ambe_Hist2DXYnuclearslice1ProjectionX_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYnuclearslice2ProjectionX->Draw();
  c12->SaveAs(Form("ambe_Hist2DXYnuclearslice2ProjectionX_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYnuclearslice3ProjectionX->Draw();
  c12->SaveAs(Form("ambe_Hist2DXYnuclearslice3ProjectionX_%s%s.png",label.c_str(),Time.c_str()));

  TCanvas *c13=new TCanvas("c13",label.c_str(),1000,600);
  Hist2DXYsinglenuclearslice1ProjectionY->Draw();
  c13->SaveAs(Form("ambe_Hist2DXYsinglenuclearslice1ProjectionY_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYsinglenuclearslice2ProjectionY->Draw();
  c13->SaveAs(Form("ambe_Hist2DXYsinglenuclearslice2ProjectionY_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYsinglenuclearslice3ProjectionY->Draw();
  c13->SaveAs(Form("ambe_Hist2DXYsinglenuclearslice3ProjectionY_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYsinglenuclearslice1ProjectionX->Draw();
  c13->SaveAs(Form("ambe_Hist2DXYsinglenuclearslice1ProjectionX_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYsinglenuclearslice2ProjectionX->Draw();
  c13->SaveAs(Form("ambe_Hist2DXYsinglenuclearslice2ProjectionX_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYsinglenuclearslice3ProjectionX->Draw();
  c13->SaveAs(Form("ambe_Hist2DXYsinglenuclearslice3ProjectionX_%s%s.png",label.c_str(),Time.c_str()));

  TCanvas *c14=new TCanvas("c14",label.c_str(),1000,600);
  Hist2DXYsingleslice1ProjectionY->Draw();
  c14->SaveAs(Form("ambe_Hist2DXYsingleslice1ProjectionY_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYsingleslice2ProjectionY->Draw();
  c14->SaveAs(Form("ambe_Hist2DXYsingleslice2ProjectionY_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYsingleslice3ProjectionY->Draw();
  c14->SaveAs(Form("ambe_Hist2DXYsingleslice3ProjectionY_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYsingleslice1ProjectionX->Draw();
  c14->SaveAs(Form("ambe_Hist2DXYsingleslice1ProjectionX_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYsingleslice2ProjectionX->Draw();
  c14->SaveAs(Form("ambe_Hist2DXYsingleslice2ProjectionX_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYsingleslice3ProjectionX->Draw();
  c14->SaveAs(Form("ambe_Hist2DXYsingleslice3ProjectionX_%s%s.png",label.c_str(),Time.c_str()));

  //write the histograms into new root file
  TFile f2D(outputname.c_str(), "RECREATE");
  Hist2D->Write();
  Hist2DYZ->Write();
  Hist2DXZ->Write();
  Hist2DXYslice1->Write();
  Hist2DXYslice2->Write();
  Hist2DXYslice3->Write();
  Hist2DXYslice1ProjectionY->Write();
  Hist2DXYslice2ProjectionY->Write();
  Hist2DXYslice3ProjectionY->Write();
  Hist2DXYslice1ProjectionX->Write();
  Hist2DXYslice2ProjectionX->Write();
  Hist2DXYslice3ProjectionX->Write();
  Hist2DGap->Write();
  py1->Write();
  QuenchingFactor->Write();
  Hist2Dnuclear->Write();
  Hist2DYZnuclear->Write();
  Hist2DXZnuclear->Write();
  Hist2DXYnuclearslice1->Write();
  Hist2DXYnuclearslice2->Write();
  Hist2DXYnuclearslice3->Write();
  Hist2DXYnuclearslice1ProjectionY->Write();
  Hist2DXYnuclearslice2ProjectionY->Write();
  Hist2DXYnuclearslice3ProjectionY->Write();
  Hist2DXYnuclearslice1ProjectionX->Write();
  Hist2DXYnuclearslice2ProjectionX->Write();
  Hist2DXYnuclearslice3ProjectionX->Write();
  Hist1DNRcountsslice1->Write();
  Hist1DNRcountsslice2->Write();
  Hist1DNRcountsslice3->Write();
  Hist1DNRcountsslice4->Write();
  Hist2Dsinglenuclear->Write();
  Hist2DYZsinglenuclear->Write();
  Hist2DXYsinglenuclearslice1->Write();
  Hist2DXYsinglenuclearslice2->Write();
  Hist2DXYsinglenuclearslice3->Write();
  Hist2DXYsinglenuclearslice1ProjectionY->Write();
  Hist2DXYsinglenuclearslice2ProjectionY->Write();
  Hist2DXYsinglenuclearslice3ProjectionY->Write();
  Hist2DXYsinglenuclearslice1ProjectionX->Write();
  Hist2DXYsinglenuclearslice2ProjectionX->Write();
  Hist2DXYsinglenuclearslice3ProjectionX->Write();
  Hist2Dsingle->Write();
  Hist2DYZsingle->Write();
  Hist2DXYsingleslice1->Write();
  Hist2DXYsingleslice2->Write();
  Hist2DXYsingleslice3->Write();
  Hist2DXYsingleslice1ProjectionY->Write();
  Hist2DXYsingleslice2ProjectionY->Write();
  Hist2DXYsingleslice3ProjectionY->Write();
  Hist2DXYsingleslice1ProjectionX->Write();
  Hist2DXYsingleslice2ProjectionX->Write();
  Hist2DXYsingleslice3ProjectionX->Write();
  Hist1DTPCedep->Write();
  Hist1DTPCedep_nuclear->Write();
  Hist1DTPCedep_electron->Write();
  Hist1DNVedep->Write();
  Hist1DNVedep_nuclear->Write();
  Hist1DNVedep_electron->Write();
  f2D.Write();
  f2D.Close();

  Info("ReconSelector::Terminate()","terminating ...");
  
}

void ReconSelector::FillHistograms(){
  
  double quenchingcut = 0.43;
  double dep_enecut= 23.0; //keV
  double dep_enecut2=269.0; //keV

  vector<double> nrx;
  vector<double> nry;
  vector<double> nrz;
  vector<double> rx;
  vector<double> ry;
  vector<double> rz;

  int eff;
  bool IsNR,IsER;
  int NNR, NER, ndeps, nsing, NRcounts;
  
  double topzslice=-3.8;
  double middlezslice1=-5.8;
  double middlezslice2=-7.8;
  double bottomzslice=-9.8;

  NNR=NER=0;
  nsing=0;
  eff=0;
  cout<<"eff=  "<<eff<<endl;
  
  ndeps=0;
  NRcounts=0;

  IsNR = false;
  IsER = false;
  
  for(size_t i=0; i< et->size(); i++)
    {
      if(volume->at(i) =="p_active")
	{ ndeps++;
	  IsER=true;
	  Hist2D->Fill(ex->at(i),ey->at(i));
	  Hist2DYZ->Fill(ey->at(i),ez->at(i));
	  Hist2DXZ->Fill(ex->at(i),ez->at(i));
	  Hist1DTPCedep->Fill(edep->at(i)/1000.);
	  Hist1DTPCedep_nuclear->Fill(edep_nuclear->at(i)/1000.);
	  Hist1DTPCedep_electron->Fill(edep_electron->at(i)/1000.);
	  QuenchingFactor->Fill(quenchingfactor->at(i));
	  if(ndeps==1)
	    { eff++;
	      rx.push_back(ex->at(i));
	      ry.push_back(ey->at(i));
	      rz.push_back(ez->at(i));
	    }
	  if(ez->at(i)<topzslice && ez->at(i)>middlezslice1)
	    Hist2DXYslice1->Fill(ex->at(i),ey->at(i));
	  if(ez->at(i)<middlezslice1 && ez->at(i)>middlezslice2)
	    Hist2DXYslice2->Fill(ex->at(i),ey->at(i));
	  if(ez->at(i)>bottomzslice && ez->at(i)<middlezslice2)
	    Hist2DXYslice3->Fill(ex->at(i),ey->at(i));
	  if(quenchingfactor->at(i) < quenchingcut && edep->at(i)<dep_enecut2 && edep->at(i)>dep_enecut )
	    {
	      //  NRcounts++;
	      IsNR = true;
	      Hist2Dnuclear->Fill(ex->at(i),ey->at(i));
	      Hist2DYZnuclear->Fill(ey->at(i),ez->at(i));
	      Hist2DXZnuclear->Fill(ex->at(i),ez->at(i));
	      nrx.push_back(ex->at(i));
	      nry.push_back(ey->at(i));
	      nrz.push_back(ez->at(i));
	      if(ez->at(i)>topzslice)
		{
		  NRcounts++;
		  Hist1DNRcountsslice4->Fill(NRcounts);
		}
	      else if(ez->at(i)<topzslice && ez->at(i)>middlezslice1)
		{  Hist2DXYnuclearslice1->Fill(ex->at(i),ey->at(i));
		  NRcounts++;
		  Hist1DNRcountsslice1->Fill(NRcounts);
		}
	      else if(ez->at(i)<middlezslice1 && ez->at(i)>middlezslice2)
		{ Hist2DXYnuclearslice2->Fill(ex->at(i),ey->at(i));
		  NRcounts++;
		  Hist1DNRcountsslice2->Fill(NRcounts);
		}
	      else if(ez->at(i)>bottomzslice && ez->at(i)<middlezslice2)
		{  Hist2DXYnuclearslice3->Fill(ex->at(i),ey->at(i));
		  NRcounts++;
		  Hist1DNRcountsslice3->Fill(NRcounts);
		}
	    }
	}

      if(volume->at(i) =="p_scint")
	{
	  Hist1DNVedep->Fill(edep->at(i)/1000.);
	  Hist1DNVedep_nuclear->Fill(edep_nuclear->at(i)/1000.);
	  Hist1DNVedep_electron->Fill(edep_electron->at(i)/1000.);	  
	}

      if(ex->at(i)>-100 && ex->at(i)<0 && ey->at(i)<0.5 && ey->at(i)>-0.5)
	{
	  Hist2DGap->Fill(ex->at(i),ez->at(i));
	}
    }

  if(IsNR)
    NNR++;
  if(IsER)
    NER++;
  if(ndeps==1)
    { nsing++;
      for(int a=0; a<ndeps; a++)
	{
	  Hist2Dsingle->Fill(rx.at(a),ry.at(a));
	  Hist2DYZsingle->Fill(ry.at(a),rz.at(a));
	  if(rz.at(a)<topzslice && rz.at(a)>middlezslice1)
	    Hist2DXYsingleslice1->Fill(rx.at(a),ry.at(a));
	  if(rz.at(a)<middlezslice1 && rz.at(a)>middlezslice2)
	    Hist2DXYsingleslice2->Fill(rx.at(a),ry.at(a));
	  if(rz.at(a)>bottomzslice && rz.at(a)<middlezslice2)
	    Hist2DXYsingleslice3->Fill(rx.at(a),ry.at(a));
	}
    }

  if(ndeps==1 && IsNR )
    {
      for(int h=0; h<ndeps; h++)
	{
	  Hist2Dsinglenuclear->Fill(nrx.at(h),nry.at(h));
	  Hist2DYZsinglenuclear->Fill(nry.at(h),nrz.at(h));
	  if(nrz.at(h)<topzslice && nrz.at(h)>middlezslice1)
	    Hist2DXYsinglenuclearslice1->Fill(nrx.at(h),nry.at(h));
	  if(nrz.at(h)<middlezslice1 && nrz.at(h)>middlezslice2)
	    Hist2DXYsinglenuclearslice2->Fill(nrx.at(h),nry.at(h));
	  if(nrz.at(h)>bottomzslice && nrz.at(h)<middlezslice2)
	    Hist2DXYsinglenuclearslice3->Fill(nrx.at(h),nry.at(h));
	}
    }
	      
}

void ReconSelector::BookHistograms(){
  Info("ReconSelector::BookHistograms()","creating histograms...");
  std::cout<<"Starting ..... ReconSelector::BookHistograms()"<<std::endl;
  
  string label = GetOption();
   /*  string label;
  if(Volume == 9) label="collimator";
  else label="no_collimator";
   */
  TList* list = GetOutputList();

  string quchf = " quenchingfactor<0.43 ";
  string quche = " quenchingfactor>0.43 ";

  string zslice1 = " -5.8<dep_z<-3.8 ";
  string zslice2 = " -7.8<dep_z<-5.8 ";
  string zslice3 = " -9.8<dep_z<-7.8 ";
  
  string Hist2DTitle = "dep_y vs dep_x in TPC active volume "+ label;
  Hist2D = new TH2F("Hist2D",Hist2DTitle.c_str(), 100,-20,20,100,-20,20);
  Hist2D->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2D->GetYaxis()->SetTitle("dep_y [cm]");
  list->Add(Hist2D);
 
  string Hist2DYZTitle = "dep_y vs dep_z in TPC active volume "+ label;
  Hist2DYZ = new TH2F("Hist2DYZ",Hist2DYZTitle.c_str(), 100,-20,20,100,-25,20);
  Hist2DYZ->GetXaxis()->SetTitle("dep_y [cm]");
  Hist2DYZ->GetYaxis()->SetTitle("dep_z [cm]");
  list->Add(Hist2DYZ);

  string Hist2DXZTitle = "dep_x vs dep_z in TPC active volume "+ label;
  Hist2DXZ = new TH2F("Hist2DXZ",Hist2DXZTitle.c_str(), 100,-20,20,100,-25,20);
  Hist2DXZ->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2DXZ->GetYaxis()->SetTitle("dep_z [cm]");
  list->Add(Hist2DXZ);

  string Hist2DXYslice1Title = "dep_y vs dep_x in TPC active volume"+zslice1+ label;
  Hist2DXYslice1 = new TH2F("Hist2DXYslice1",Hist2DXYslice1Title.c_str(),100,-20,20,100,-20,20);
  Hist2DXYslice1->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2DXYslice1->GetYaxis()->SetTitle("dep_y [cm]");
  list->Add(Hist2DXYslice1);

  string Hist2DXYslice2Title = "dep_y vs dep_x in TPC active volume"+zslice2+ label;
  Hist2DXYslice2 = new TH2F("Hist2DXYslice2",Hist2DXYslice2Title.c_str(),100,-20,20,100,-20,20);
  Hist2DXYslice2->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2DXYslice2->GetYaxis()->SetTitle("dep_y [cm]");
  list->Add(Hist2DXYslice2);

  string Hist2DXYslice3Title = "dep_y vs dep_x in TPC active volume"+zslice3+ label;
  Hist2DXYslice3 = new TH2F("Hist2DXYslice3",Hist2DXYslice3Title.c_str(),100,-20,20,100,-20,20);
  Hist2DXYslice3->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2DXYslice3->GetYaxis()->SetTitle("dep_y [cm]");
  list->Add(Hist2DXYslice3);

  string Hist2DGapTitle = "dep_z vs dep_x with y slice "+ label;
  Hist2DGap = new TH2F("Hist2DGap",Hist2DGapTitle.c_str(), 150,-100,0,150,-50,100);
  Hist2DGap->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2DGap->GetYaxis()->SetTitle("dep_z [cm]");
  list->Add(Hist2DGap);
  
  string QuenchingFactorTitle = "quenchingfactor in TPC active volume "+ label;
  QuenchingFactor = new TH1F("QuenchingFactor",QuenchingFactorTitle.c_str(),100,0,1);
  QuenchingFactor->GetXaxis()->SetTitle("quenching factor");
  QuenchingFactor->GetYaxis()->SetTitle("counts");
  list->Add(QuenchingFactor);

  //quecnhing cuts for nuclear-like recoil
  string Hist2DnuclearTitle = "dep_y vs dep_x in TPC active volume"+quchf+ label;
  Hist2Dnuclear = new TH2F("Hist2Dnuclear",Hist2DnuclearTitle.c_str(), 100,-20,20,100,-20,20);
  Hist2Dnuclear->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2Dnuclear->GetYaxis()->SetTitle("dep_y [cm]");
  list->Add(Hist2Dnuclear);

  string Hist2DYZnuclearTitle = "dep_y vs dep_z in TPC active volume"+quchf+ label;
  Hist2DYZnuclear = new TH2F("Hist2DYZnuclear",Hist2DYZnuclearTitle.c_str(), 100,-20,20,100,-25,20);
  Hist2DYZnuclear->GetXaxis()->SetTitle("dep_y [cm]");
  Hist2DYZnuclear->GetYaxis()->SetTitle("dep_z [cm]");
  list->Add(Hist2DYZnuclear);

  string Hist2DXZnuclearTitle = "dep_x vs dep_z in TPC active volume"+quchf+ label;
  Hist2DXZnuclear = new TH2F("Hist2DXZnuclear",Hist2DXZnuclearTitle.c_str(), 100,-20,20,100,-25,20);
  Hist2DXZnuclear->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2DXZnuclear->GetYaxis()->SetTitle("dep_z [cm]");
  list->Add(Hist2DXZnuclear);

  string Hist2DXYnuclearslice1Title = "dep_y vs dep_x in TPC active volume"+zslice1+quchf+ label;
  Hist2DXYnuclearslice1 = new TH2F("Hist2DXYnuclearslice1",Hist2DXYnuclearslice1Title.c_str(),100,-20,20,100,-20,20);
  Hist2DXYnuclearslice1->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2DXYnuclearslice1->GetYaxis()->SetTitle("dep_y [cm]");
  list->Add(Hist2DXYnuclearslice1);

  string Hist2DXYnuclearslice2Title = "dep_y vs dep_x in TPC active volume"+zslice2+quchf+ label;
  Hist2DXYnuclearslice2 = new TH2F("Hist2DXYnuclearslice2",Hist2DXYnuclearslice2Title.c_str(),100,-20,20,100,-20,20);
  Hist2DXYnuclearslice2->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2DXYnuclearslice2->GetYaxis()->SetTitle("dep_y [cm]");
  list->Add(Hist2DXYnuclearslice2);

  string Hist2DXYnuclearslice3Title = "dep_y vs dep_x in TPC active volume"+zslice3+quchf+ label;
  Hist2DXYnuclearslice3 = new TH2F("Hist2DXYnuclearslice3",Hist2DXYnuclearslice3Title.c_str(),100,-20,20,100,-20,20);
  Hist2DXYnuclearslice3->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2DXYnuclearslice3->GetYaxis()->SetTitle("dep_y [cm]");
  list->Add(Hist2DXYnuclearslice3);

  string Hist1DNRcountsslice1Title = "NR counts in slice1 in TPC active volume "+ label;
  Hist1DNRcountsslice1 = new TH1F("Hist1DNRcountsslice1",Hist1DNRcountsslice1Title.c_str(),100,0,20);
  Hist1DNRcountsslice1->GetXaxis()->SetTitle("NR counts");
  list->Add(Hist1DNRcountsslice1);
  
  string Hist1DNRcountsslice2Title = "NR counts in slice2 in TPC active volume "+ label;
  Hist1DNRcountsslice2 = new TH1F("Hist1DNRcountsslice2",Hist1DNRcountsslice2Title.c_str(),100,0,20);
  Hist1DNRcountsslice2->GetXaxis()->SetTitle("NR counts");
  list->Add(Hist1DNRcountsslice2);

  string Hist1DNRcountsslice3Title = "NR counts in slice3 in TPC active volume "+ label;
  Hist1DNRcountsslice3 = new TH1F("Hist1DNRcountsslice3",Hist1DNRcountsslice3Title.c_str(),100,0,20);
  Hist1DNRcountsslice3->GetXaxis()->SetTitle("NR counts");
  list->Add(Hist1DNRcountsslice3);

  string Hist1DNRcountsslice4Title = "NR counts in slice4 in TPC active volume "+ label;
  Hist1DNRcountsslice4 = new TH1F("Hist1DNRcountsslice4",Hist1DNRcountsslice4Title.c_str(),100,0,20);
  Hist1DNRcountsslice4->GetXaxis()->SetTitle("NR counts");
  list->Add(Hist1DNRcountsslice4);

  //single nuclear sacttering in TPC active volume
  string Hist2DsinglenuclearTitle = "dep_y vs dep_x for single nuclear scattering in TPC active volume"+quchf+ label;
  Hist2Dsinglenuclear = new TH2F("Hist2Dsinglenuclear",Hist2DsinglenuclearTitle.c_str(), 100,-20,20,100,-20,20);
  Hist2Dsinglenuclear->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2Dsinglenuclear->GetYaxis()->SetTitle("dep_y [cm]");
  list->Add(Hist2Dsinglenuclear);

  string Hist2DYZsinglenuclearTitle = "dep_z vs dep_y for single nuclear scattering in TPC active volume"+quchf+ label;
  Hist2DYZsinglenuclear = new TH2F("Hist2DYZsinglenuclear",Hist2DYZsinglenuclearTitle.c_str(),100,-20,20,100,-25,20);
  Hist2DYZsinglenuclear->GetXaxis()->SetTitle("dep_y [cm]");
  Hist2DYZsinglenuclear->GetYaxis()->SetTitle("dep_z [cm]");
  list->Add(Hist2DYZsinglenuclear);

  string Hist2DXYsinglenuclearslice1Title = "dep_y vs dep_x for single nuclear scattering in TPC active volume"+zslice1+quchf+ label;
  Hist2DXYsinglenuclearslice1 = new TH2F("Hist2DXYsinglenuclearslice1",Hist2DXYsinglenuclearslice1Title.c_str(),100,-20,20,100,-20,20);
  Hist2DXYsinglenuclearslice1->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2DXYsinglenuclearslice1->GetYaxis()->SetTitle("dep_y [cm]");
  list->Add(Hist2DXYsinglenuclearslice1);

  string Hist2DXYsinglenuclearslice2Title = "dep_y vs dep_x for single nuclear scattering in TPC active volume"+zslice2+quchf+ label;
  Hist2DXYsinglenuclearslice2 = new TH2F("Hist2DXYsinglenuclearslice2",Hist2DXYsinglenuclearslice2Title.c_str(),100,-20,20,100,-20,20);
  Hist2DXYsinglenuclearslice2->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2DXYsinglenuclearslice2->GetYaxis()->SetTitle("dep_y [cm]");
  list->Add(Hist2DXYsinglenuclearslice2);

  string Hist2DXYsinglenuclearslice3Title = "dep_y vs dep_x for single nuclear scattering in TPC active volume"+zslice3+quchf+ label;
  Hist2DXYsinglenuclearslice3 = new TH2F("Hist2DXYsinglenuclearslice3",Hist2DXYsinglenuclearslice3Title.c_str(),100,-20,20,100,-20,20);
  Hist2DXYsinglenuclearslice3->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2DXYsinglenuclearslice3->GetYaxis()->SetTitle("dep_y [cm]");
  list->Add(Hist2DXYsinglenuclearslice3);

  //single sacttering in TPC active volume
  string Hist2DsingleTitle = "dep_y vs dep_x for single scattering in TPC active volume "+ label;
  Hist2Dsingle = new TH2F("Hist2Dsingle",Hist2DsingleTitle.c_str(), 100,-20,20,100,-20,20);
  Hist2Dsingle->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2Dsingle->GetYaxis()->SetTitle("dep_y [cm]");
  list->Add(Hist2Dsingle);

  string Hist2DYZsingleTitle = "dep_z vs dep_y for single scattering in TPC active volume "+ label;
  Hist2DYZsingle = new TH2F("Hist2DYZsingle",Hist2DYZsingleTitle.c_str(),100,-20,20,100,-25,20);
  Hist2DYZsingle->GetXaxis()->SetTitle("dep_y [cm]");
  Hist2DYZsingle->GetYaxis()->SetTitle("dep_z [cm]");
  list->Add(Hist2DYZsingle);

  string Hist2DXYsingleslice1Title = "dep_y vs dep_x for single scattering in TPC active volume"+zslice1+ label;
  Hist2DXYsingleslice1 = new TH2F("Hist2DXYsingleslice1",Hist2DXYsingleslice1Title.c_str(),100,-20,20,100,-20,20);
  Hist2DXYsingleslice1->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2DXYsingleslice1->GetYaxis()->SetTitle("dep_y [cm]");
  list->Add(Hist2DXYsingleslice1);

  string Hist2DXYsingleslice2Title = "dep_y vs dep_x for single scattering in TPC active volume"+zslice2+ label;
  Hist2DXYsingleslice2 = new TH2F("Hist2DXYsingleslice2",Hist2DXYsingleslice2Title.c_str(),100,-20,20,100,-20,20);
  Hist2DXYsingleslice2->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2DXYsingleslice2->GetYaxis()->SetTitle("dep_y [cm]");
  list->Add(Hist2DXYsingleslice2);

  string Hist2DXYsingleslice3Title = "dep_y vs dep_x for single scattering in TPC active volume"+zslice3+ label;
  Hist2DXYsingleslice3 = new TH2F("Hist2DXYsingleslice3",Hist2DXYsingleslice3Title.c_str(),100,-20,20,100,-20,20);
  Hist2DXYsingleslice3->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2DXYsingleslice3->GetYaxis()->SetTitle("dep_y [cm]");
  list->Add(Hist2DXYsingleslice3);

  //energy spectrum inside TPC
  string Hist1DTPCedepTitle = "Energy Spectrum in TPC active volume "+ label;
  Hist1DTPCedep = new TH1F("Hist1DTPCedep",Hist1DTPCedepTitle.c_str(), 1000,0,20);
  Hist1DTPCedep->GetXaxis()->SetTitle("Energy [MeV]");
  Hist1DTPCedep->GetYaxis()->SetTitle("Counts");
  list->Add(Hist1DTPCedep);
  
  string Hist1DTPCedep_nuclearTitle = "Energy Spectrum of nuclear recoils in TPC active volume "+ label;
  Hist1DTPCedep_nuclear = new TH1F("Hist1DTPCedep_nuclear",Hist1DTPCedep_nuclearTitle.c_str(), 1000,0,20);
  Hist1DTPCedep_nuclear->GetXaxis()->SetTitle("Energy [MeV]");
  Hist1DTPCedep_nuclear->GetYaxis()->SetTitle("Counts");
  list->Add(Hist1DTPCedep_nuclear);
  
  string Hist1DTPCedep_electronTitle = "Energy Spectrum of electron recoils in TPC active volume "+ label;
  Hist1DTPCedep_electron = new TH1F("Hist1DTPCedep_electron",Hist1DTPCedep_electronTitle.c_str(), 1000,0,20);
  Hist1DTPCedep_electron->GetXaxis()->SetTitle("Energy [MeV]");
  Hist1DTPCedep_electron->GetYaxis()->SetTitle("Counts");
  list->Add(Hist1DTPCedep_electron);
  
  //energy spectrum inside NV
  string Hist1DNVedepTitle = "Energy Spectrum in NV active volume "+ label;
  Hist1DNVedep = new TH1F("Hist1DNVedep",Hist1DNVedepTitle.c_str(), 1000,0,20);
  Hist1DNVedep->GetXaxis()->SetTitle("Energy [MeV]");
  Hist1DNVedep->GetYaxis()->SetTitle("Counts");
  list->Add(Hist1DNVedep);
  
  string Hist1DNVedep_nuclearTitle = "Energy Spectrum of nuclear recoils in NV "+ label;
  Hist1DNVedep_nuclear = new TH1F("Hist1DNVedep_nuclear",Hist1DNVedep_nuclearTitle.c_str(), 1000,0,20);
  Hist1DNVedep_nuclear->GetXaxis()->SetTitle("Energy [MeV]");
  Hist1DNVedep_nuclear->GetYaxis()->SetTitle("Counts");
  list->Add(Hist1DNVedep_nuclear);
  
  string Hist1DNVedep_electronTitle = "Energy Spectrum of electron recoils in NV "+ label;
  Hist1DNVedep_electron = new TH1F("Hist1DNVedep_electron",Hist1DNVedep_electronTitle.c_str(), 1000,0,20);
  Hist1DNVedep_electron->GetXaxis()->SetTitle("Energy [MeV]");
  Hist1DNVedep_electron->GetYaxis()->SetTitle("Counts");
  list->Add(Hist1DNVedep_electron);
  

}
