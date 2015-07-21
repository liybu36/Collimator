#define SladReadSelector_cxx
// The class definition in SladReadSelector.h has been generated automatically
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
// root> T->Process("SladReadSelector.C")
// root> T->Process("SladReadSelector.C","some options")
// root> T->Process("SladReadSelector.C+")
//

#include "SladReadSelector.h"
#include <TH2.h>
#include <TStyle.h>
#include <TMath.h>

using namespace std;

void SladReadSelector::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
//   Info("SladReadSelector::Begin()",Form("GetOption() = %s",option.Data()));
//   std::cout << "GetOption() = " << option.Data() << std::endl;
   TObjArray* optarr = option.Tokenize(":");
   TString outdir = ((TObjString*)optarr->At(0))->GetString();
   SetOutputFile(outdir);

}

void SladReadSelector::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).
  Info("SladReadSelector::SlaveBegin()","begining ...");

//  TString option = GetOption();

   BookHistograms();
   SetCutG4Outliers();

   TString profilename(fInput->FindObject("profilename")->GetTitle());
   this->SetProfileName(profilename);

   fMyFcn = new XYReconstructor();
   fMyFcn->LoadProfile(fProfileName.Data());
   fMyFcn->SetMinimizer();

#define OUTTREE

#ifdef OUTTREE
//   fOut_Tree = new TTree("TAlphas","Tree for high energy alpha candidates");
   fOut_Tree = new TTree("TOutliers","Tree for outliers (DM candidates)");
   fOut_Tree->Branch("run_id",&run_id,"run_id/I");
   fOut_Tree->Branch("event_id",&event_id,"event_id/I");
   fOut_Tree->Branch("total_s1",&total_s1,"total_s1/F");
   fOut_Tree->Branch("total_s1_corr",&total_s1_corr,"total_s1_corr/F");
   fOut_Tree->Branch("total_f90",&total_f90,"total_f90/F");
   fOut_Tree->Branch("total_s2",&total_s2,"total_s2/F");
   fOut_Tree->Branch("total_s2_corr",&total_s2_corr,"total_s2_corr/F");
   fOut_Tree->Branch("total_s2_true",&total_s2_true,"total_s2_true/F");
   fOut_Tree->Branch("t_drift",&t_drift,"t_drift/D");
   fOut_Tree->Branch("radius",&radius,"radius/F");
   fOut_Tree->Branch("s1_max_frac",&s1_max_frac,"s1_max_frac/F");
   fOut_Tree->Branch("mf99",&mf99,"mf99/F");
//   fOut_Tree->Branch("s2",&s2,"s2[38]/D");

   fOutput->Add(fOut_Tree);
#endif
}

Bool_t SladReadSelector::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either SladReadSelector::GetEntry() or TBranch::GetEntry()
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
  fChain->GetTree()->GetEntry(entry);

  hRun_dT_Trigger->Fill(livetime+inhibittime, run_id);
  hRun_livetime->Fill(livetime, run_id);
  hRun_inhibittime->Fill(inhibittime, run_id);

//  Double_t radius=sqrt(pow(masas_x,2)+pow(masas_y,2));
//  Double_t total_s2_true = (total_s2_corr*max(1,33.48/(33.48 - 1.099*radius - 0.0288*pow(radius,2) + 0.01662*pow(radius,3) - 0.002384*pow(radius,4) + 8.481E-5*pow(radius,5))));


//  hTotal_livetime->Fill(0.5,livetime-std::max(1.35e-3 - inhibittime, 0.));
  hTotal_livetime->Fill(0.5,livetime);
//  if(IsValidEvent()) {
//      hTotal_livetime->Fill(0.5,livetime-std::max(1.35e-3 - inhibittime, 0.));
//  }
  if(IsEventOK()) FillHistograms();

   return kTRUE;
}

void SladReadSelector::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void SladReadSelector::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
  Info("SladReadSelector::Terminate()","closining ...");

  Int_t nhists = 0;
  TList *list = GetOutputList();
//  list->Print();
  TIter next(list);
  TFile *hOutF = new TFile(fOutName.Data(), "RECREATE");
  TH1 *h;
  TObject* obj;
  while ((obj = (TObject*) next()))
  {
      if(!obj) continue;
      if (obj->InheritsFrom(TH1::Class()))
      {
          h = (TH1*) obj;
//          if(h->InheritsFrom(TH2::Class())){
//          }
          h->Write();

//          if (Debug)cout << "Write " << h->GetName() << " into the file" << endl;
//          delete h;
          nhists++;
      } else if (obj->InheritsFrom(TTree::Class())) {
          fOut_Tree = (TTree*) obj;
          fOut_Tree->Write();
      }
  }
  fCutg->Write();

  hOutF->Close();
}

Bool_t SladReadSelector::IsValidEvent(){
  /*  // Define cuts
    string CX1     = "(nchannel==38)";                                              // Number of channels
    string CX2     = "(SumChannelHasNoBaseline==0)";                                // Baseline
    string CX3     = "(1.35E-3<(lifetime+inhibittime))";                            // Livetime and inhibit time
    string CX4     = "(lifetime<1.)";                                               // Stalled DAQ
    string CX5     = "(veto_present)";                                              // Veto present
  */

  Bool_t isEventValid = true;
  for (Int_t i=1; i<6; i++) {
      if(selection_results[i]==0) continue; // ==0 means event is accepted

      isEventValid = false;
      break;
  }
  return isEventValid;
}

Bool_t SladReadSelector::IsEventOK_Default(){
/*  // Define cuts
  string CX1     = "(nchannel==38)";                                              // Number of channels
  string CX2     = "(SumChannelHasNoBaseline==0)";                                // Baseline
  string CX3     = "(1.35E-3<(lifetime+inhibittime))";                            // Livetime and inhibit time
  string CX4     = "(lifetime<1.)";                                               // Stalled DAQ
  string CX5     = "(veto_present)";                                              // Veto present
  string CX6     = "(50.<abs(closest_time_wrt_prompt_time_25pe_thr))";            // Coincident veto signal
  //string CX6     = "((50.<abs(closest_time_wrt_prompt_time)) || ((abs(closest_time_wrt_prompt_time)<=50.) && (multiplicity_closest_time_wrt_prompt_time<25)))";
  // Coincident veto signal
  string CX7     = "(max_charge<150)";                                            // High energy veto signal
  string CX8     = "(n_phys_pulses==2||(has_s3==1&&n_phys_pulses==3))";           // Number of pulses
  string CX9     = "((run_id<7344 &&-0.25<s1_start_time&&s1_start_time<-0.15) || (run_id>7640&&-6.10<s1_start_time&&s1_start_time<-6.00))";
  // Trigger time
  string CX10    = "(is_saturated_pulse0==0)";                                    // No S1 saturation
  string CX11    = "(s1_max_frac<s1_max_frac_threshold)";                         // S1 maximum fraction
  string CX12    = "(total_s2_f90_fixed<0.20)";                                   // Valid S2
  string CX13    = "(30.<total_s2_corr)";                                         // Substantial S2
  string CX14    = "((60.<total_s1_corr)&&(total_s1_corr<200.))";                 // S1 range
  string CX15    = "((25.0<tdrift)&&(tdrift<353.3))";                             // Drift time */

  Bool_t isEventOK = true;
  for (Int_t i=1; i<16; i++) {
      if(i==11) continue; // without s1 max fraction cut
      if(i==10) continue; // without s1 saturation cut
      if(i==14) continue; // without s1 cut
      if(i==15) continue; // without t_drift cut
      if(i==6) continue; // without coincidental veto signal cut
      if(i==7) continue; // without high energy veto signal cut
      if(selection_results[i]==0) continue; // ==0 means event is accepted

      isEventOK = false;
      break;
  }
  return isEventOK;
}

Bool_t SladReadSelector::IsEventOK(){
  // Define cuts // from cutdef.h on 09/18/2014
/*#define X1      (nchannel==38)                                                                                                          // Number of channels
#define X2      (SumChannelHasNoBaseline==0)                                                                                            // Baseline
#define X3      (1.35E-3<(lifetime+inhibittime))                                                                                        // Livetime and inhibit time
#define X4      (lifetime<1.)                                                                                                           // Stalled DAQ
#define X5      (veto_present)  */                                                                                                        // Veto present
#define X6      (roi_total_prompt_charge<10.)                                                                                            // Coincident veto signal
#define X7      (lsv_max_window_charge<80.&&lsv_late_window_charge<110.&&wt_charge<200.)                                                   // High energy veto signal
//#define X8      (n_phys_pulses==2||(has_s3==1&&n_phys_pulses==3))                                                                       // Number of pulses
//#define X9      ((run_id<7344 &&-0.25<s1_start_time&&s1_start_time<-0.15) || (run_id>7640&&-6.10<s1_start_time&&s1_start_time<-6.00))   // Trigger time
//#define X10     (is_saturated_pulse0==0)                                                                                                // S1 saturation
#define X11     (s1_max_frac<mf99)                                                                                                      // S1 maximum fraction
//#define X11     (s1_max_frac<0.4)                                                                                                      // S1 maximum fraction
//#define X12     (total_s2_f90_fixed<0.20)                                                                                               // Valid S2
#define X13     (100.>total_s2_true)                                                                                                    // Minimal S2
#define X14     ((60.<total_s1_corr)&&(total_s1_corr<200.))                                                                             // S1 range
#define X15     ((40.0<t_drift)&&(t_drift<334.5))                                                                                         // Drift time

  Bool_t isEventOK = true;
#if 0
  for (Int_t i=1; i<16; i++) {
      if (i==6) {
          if(X6) continue;  // coincidental veto signal cut
      } else if (i==7) {
          if(X7) continue;  // high energy veto signal cut
      } else if (i==10) {
          continue; // without s1 saturation cut
//          if(X10) continue; // s1 saturation cut
      } else if (i==11) {
          continue; // s1 max fraction cut
//          if(X11) continue; // s1 max fraction cut
      } else if (i==13) {
          continue;
//          if(X13) continue; // valid s2 cut
      } else if (i==14) {
//          continue;         // without s1 cut
          if(X14) continue;   // s1 cut
      } else if (i==15) {
          continue;         // without drift time cut
//          if(X15) continue; // drift time cut
      } else if(selection_results[i]==0) continue; // ==0 means event is accepted
#else //for single extracted electron study
      for (Int_t i=1; i<16; i++) {
          if (i==3) {
              continue;  // without short livetime cut
          } else if (i==6) {
              continue;  // coincidental veto signal cut
          } else if (i==7) {
              continue;  // high energy veto signal cut
//          } else if (i==8) {
//              if (n_phys_pulses==2) continue; // only 2 pulses
          } else if (i==9) {
              continue;  // without start time cut
          }else if (i==10) {
              continue; // without s1 saturation cut
          } else if (i==11) {
              continue; // without s1 max fraction cut
          } else if (i==13) {
              continue; // without valid s2 cut
          } else if (i==14) {
              continue; // without s1 cut
          } else if (i==15) {
              if((374.<t_drift)&&(t_drift<390.)) continue; // drift time cut
          } else if(selection_results[i]==0) continue; // ==0 means event is accepted
      isEventOK = false;
      break;
  }

  if(total_f90>0.05) isEventOK = false; //small f90 only ->S2 tirrigered event
#endif

  return isEventOK;
}


void SladReadSelector::FillHistograms(){
  pS2OverS1corr_S1corr_XY->Fill(masas_x, masas_y, total_s1_corr, total_s2/total_s1_corr);
  pS2_S1corr_XY->Fill(masas_x, masas_y, total_s1_corr, total_s2);

  Double_t r = TMath::Sqrt(masas_x*masas_x+masas_y*masas_y);
  hS2overS1corr_R_S1corr->Fill(total_s1_corr, r, total_s2/total_s1_corr); // s2_over_s1_corr is value before S2 charge normalization

  hTDrift->Fill(t_drift);
  hRun_TDrift->Fill(t_drift, run_id);
  hTotalS1corr_TDrift->Fill(t_drift, total_s1_corr);
  hTotalF90_TDrift->Fill(t_drift, total_f90);

  hS1Total_Ext_vs_F90->Fill(total_s1, total_f90);
  hTotalS1corr_TotalF90->Fill(total_s1_corr, total_f90);

  hR_TDrift_Max->Fill(r, t_drift);
  hR_TDrift_Min->Fill(r, t_drift);

  if(total_s1>6.e+3 && total_s1<=50.e+3 && total_f90>0.4 && total_f90<=0.8){
      hR_TDrift->Fill(r, t_drift);
#ifdef OUTTREE
//      fOut_Tree->Fill();
#endif
  }

//  hR_TotalS1_TotalS2->Fill(r, total_s1, total_s2);
  hR_TotalS1_TotalS2->Fill(r, TMath::Log10(total_s1), total_s2);

  Double_t S2corr_fac = fMyFcn->XYCorrection4S2(masas_x, masas_y);
  hR_TotalS1_TotalS2_corr->Fill(r, TMath::Log10(total_s1), S2corr_fac*total_s2);

  if(fCutg->IsInside(total_s1, total_f90)){
#ifdef OUTTREE
      fOut_Tree->Fill();
#endif

  }
}

void SladReadSelector::BookHistograms()
{
//  if(Debug)
    Info("SladReadSelector::BookHistograms()","booking histogram...");

//  Double_t fDriftTimeMax = 373.3;
  Int_t nBinX(180), nBinY(180);//, nBinR(180);
  Double_t xmin(-20.), xmax(20.), ymin(-20.), ymax(20.);//, rmin(-36.), rmax(36.);
  Int_t nBinS1(80);
  Double_t S1min(0.), S1max(800);
//  Double_t minDriftT(0.), maxDriftT((int)(fDriftTimeMax*1.2));
  TList  *list = GetOutputList();

  hTotal_livetime = new TH1F("hTotal_livetime", "total_livetime", 1, 0, 1);
  list->Add(hTotal_livetime);

  hRun_livetime = new TH2D("hRun_livetime", "Run ID vs livetime;livetime [s];run_id", 250, 0., 2.5e-3, 3200, 5299.5, 8499.5);
  list->Add(hRun_livetime);

  hRun_inhibittime = new TH2D("hRun_inhibittime", "Run ID vs inhibittime;inhibittime [s];run_id", 250, 0., 2.5e-3, 3200, 5299.5, 8499.5);
  list->Add(hRun_inhibittime);

  hRun_dT_Trigger = new TH2D("hRun_dT_Trigger", "Run ID vs time from previous trigger (inhibit + live time);inhibittime + livetime [s];run_id", 250, 0., 2.5e-3, 3200, 5299.5, 8499.5);
  list->Add(hRun_dT_Trigger);

  pS2OverS1corr_S1corr_XY = new TProfile3D("pS2OverS1corr_S1corr_XY", "total_s2/total_s1_corr vs total_s1_corr vs xy.;x [cm];y [cm];total_s1_corr [PE]", nBinX, xmin, xmax, nBinY, ymin, ymax, nBinS1, S1min, S1max);
  list->Add(pS2OverS1corr_S1corr_XY);

  pS2_S1corr_XY = new TProfile3D("pS2_S1corr_XY", "total_s2 vs total_s1_corr vs xy.;x [cm];y [cm];total_s1_corr [PE]", nBinX, xmin, xmax, nBinY, ymin, ymax, nBinS1, S1min, S1max);
  list->Add(pS2_S1corr_XY);

  hS2overS1corr_R_S1corr = new TH3D("hS2overS1corr_R_S1corr", "total_s2/total_s1_corr vs Distance from TPC center vs total_s1_corr;total_s1_corr [PE];r [cm];total_s2/total_s1_cor", nBinS1, S1min, S1max, 200, 0., 20, 600, 0, 150);
  list->Add(hS2overS1corr_R_S1corr);

  hTDrift = new TH1D("hTDrift", "Drift time; t_drift [#mus]", 400, 0., 400.);
  list->Add(hTDrift);

  hRun_TDrift = new TH2D("hRun_TDrift", "Run ID vs Drift time; t_drift [#mus];run_id", 800, 0., 400., 3200, 5299.5, 8499.5);
  list->Add(hRun_TDrift);

  hTotalS1corr_TDrift = new TH2D("hTotalS1corr_TDrift", "Total S1(corr.) vs Drift time; t_drift [#mus];total_s1_corr [PE]", 400, 0., 400., 400, 0, 400);
  list->Add(hTotalS1corr_TDrift);

  hTotalF90_TDrift = new TH2D("hTotalF90_TDrift", "Total F90 vs Drift time; t_drift [#mus];total_f90", 400, 0., 400., 110, 0, 1.1);
  list->Add(hTotalF90_TDrift);

  hS1Total_Ext_vs_F90 = new TH2D("hS1Total_Ext_vs_F90", "S1 Spectrum vs f90;total_s1 [PE];total_f90", 35000, 0., 7e+4, 110, 0, 1.1);
  list->Add(hS1Total_Ext_vs_F90);

  hR_TDrift = new TH2D("hR_TDrift", "t_drift vs R;R [cm];t_drift [#mus]", 200, 0., 20, 400, 0., 400.);
  list->Add(hR_TDrift);

  hR_TDrift_Max = new TH2D("hR_TDrift_Max", "t_drift vs R around t_drift_max;R [cm];t_drift [#mus]", 200, 0., 20, 100, 370, 380.);
  list->Add(hR_TDrift_Max);

  hR_TDrift_Min = new TH2D("hR_TDrift_Min", "t_drift vs R at small t_drift;R [cm];t_drift [#mus]", 200, 0., 20, 100, 0, 10.);
  list->Add(hR_TDrift_Min);

  hTotalS1corr_TotalF90 = new TH2D("hTotalS1corr_TotalF90", "F90 vs S1 (corrected for z-dependence); total_s1_corr [PE]; total_f90", 1000, 0, 1000, 110, 0, 1.1);
  list->Add(hTotalS1corr_TotalF90);

//  hR_TotalS1_TotalS2 = new TH3D("hR_TotalS1_TotalS2", "total_s2 (S3) vs total_s1 (S2) vs R;r [cm];total_s1 (S2) [PE];total_s2 (S3) [PE]", 6, 0., 18., 10, 0., 5000., 160, -10., 150);
  hR_TotalS1_TotalS2 = new TH3D("hR_TotalS1_TotalS2", "total_s2 (S3) vs total_s1 (S2) vs R;r [cm];Log10(total_s1) (S2) [PE];total_s2 (S3) [PE]", 7, 0., 21., 10,-1, +4., 160, -10., 150);
  list->Add(hR_TotalS1_TotalS2);

  hR_TotalS1_TotalS2_corr = new TH3D("hR_TotalS1_TotalS2_corr", "total_s2_corr (S3, XY corrected) vs total_s1 (S2) vs R;r [cm];Log10(total_s1) (S2) [PE];total_s2_xycorr (S3) [PE]", 7, 0., 21., 10,-1, +4., 160, -10., 150);
  list->Add(hR_TotalS1_TotalS2_corr);


}

void SladReadSelector::SetCutG4Outliers(){
  fCutg = new TCutG("outlierscut",9);
  fCutg->SetVarX("total_s1_corr");
  fCutg->SetVarY("total_f90");
  fCutg->SetPoint(0,40,0.75);
  fCutg->SetPoint(1,60,0.66);
  fCutg->SetPoint(2,80,0.62);
  fCutg->SetPoint(3,100,0.58);
  fCutg->SetPoint(4,120,0.56);
  fCutg->SetPoint(5,160,0.52);
  fCutg->SetPoint(6,1000,0.52);
  fCutg->SetPoint(7,1000,1.1);
  fCutg->SetPoint(8,40,1.1);
  fCutg->SetPoint(9,40,0.7);

  fOutput->Add(fCutg);

}
