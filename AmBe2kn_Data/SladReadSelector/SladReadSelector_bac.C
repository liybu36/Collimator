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
#include <iostream>

#include "../XYReconstructor/PMTGeom.h"

using namespace std;

enum TimeCategory_t {RUNT, LT, IT, LTS1, LTS3};
const Int_t nTimeBins = 5;
const TString timelabels[nTimeBins] = { "Run Time", "Live Time", "Inhibit Time", "LT_{S1}", "LT_{S3}"};

enum EvntCut_t {TOTAL, ACCEPTED, NCHANNEL, BASELINE, BASIC, LIVETIME, NPULSE, TRIGERWINDOW, SMALLF90, VALIDS2};
const Int_t nEventBins = 20;
const TString eventlabels[nEventBins] = {"Total", "Accepted", "Number of Channels", "Baseline Cut", "After Basic Cuts", "Short livetime", "# of pulse Cut", "Trigger Time Cut", "Small F90", "Valid S2"};

#define JASONS_XY

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

#define CorrXY

void SladReadSelector::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).
  Info("SladReadSelector::SlaveBegin()","begining ...");

//  TString option = GetOption();

   SetCutG4Outliers();

#ifdef CorrXY
   TString fname = "S2CorrectionFactor_Kr_pS2VsXY_Kr_Jason_Masa.root";
#else
   TString fname = "S2CorrectionFactor_Kr_pR2_theta_Kr.root";
#endif
   TFile *g = (TFile*)gROOT->GetListOfFiles()->FindObject(fname.Data());
   if (!g || !g->IsOpen()) {
       g = new TFile(fname.Data());
       TString hname = "hS2KrVsXY_norm_center";

       fhS2Correction_factor = (TH2D*) g->Get(hname.Data());
       if(!fhS2Correction_factor) std::cout<<"profile: "<<hname.Data()<<" is not found."<<std::endl;
       fhS2Correction_factor->SetDirectory(0);

       g->Close();
   }
//   TString profilename(fInput->FindObject("profilename")->GetTitle());
//   this->SetProfileName(profilename);

//   fMyFcn = new XYReconstructor();
//   fMyFcn->LoadProfile(fProfileName.Data());
//   fMyFcn->SetMinimizer();

//#define OUTTREE

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
//   fOut_Tree->Branch("total_s2_true",&total_s2_true,"total_s2_true/F");
   fOut_Tree->Branch("t_drift",&t_drift,"t_drift/D");
   fOut_Tree->Branch("radius",&radius,"radius/F");
   fOut_Tree->Branch("s1_max_frac",&s1_max_frac,"s1_max_frac/F");
   fOut_Tree->Branch("max_s1_frac_cut_threshold99",&max_s1_frac_cut_threshold99,"max_s1_frac_cut_threshold99/F");
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

  if(!fIsFieldSet || (chainentry==0 && (!CheckFieldSetting(acqui_window) || GetExpectedAcquisitionWindow(fEdrift)+50.<acqui_window))) {
      Info("SladReadSelector::Process()","Field setting is wrong... Please check it.");
      fIsFieldSet = false;
      Int_t Edrift = GetFieldFromAcquisitionWindow(acqui_window);
      //Drift field setting
      SetDriftFieldConf(Edrift, run_id);
      BookHistograms();
//          return kFALSE;
  }


  if (total_s1!=0 && fPrevious_total_s1 == total_s1) {
      cout<<"total_s1 is repeated!!"<<endl;
      cout<<"run: "<<run_id<<" event: "<<event_id<<" total_s1: "<<total_s1<<endl;
      return kFALSE;
  }

  fPrevious_total_s1 = total_s1;

  hRun_dT_Trigger->Fill(livetime+inhibittime, run_id);
  hRun_livetime->Fill(livetime, run_id);
  hRun_inhibittime->Fill(inhibittime, run_id);


    hRunTime->Fill(LT, livetime);
    hRunTime->Fill(IT, inhibittime);
    if(livetime>=200.e-6) // longer than 200 us
       hRunTime->Fill(LTS1, livetime);
    else
       hRunTime->Fill(LTS3, livetime);

    hEventConter->Fill(TOTAL);



//  Double_t radius=sqrt(pow(masas_x,2)+pow(masas_y,2));
//  Double_t total_s2_true = (total_s2_corr*max(1,33.48/(33.48 - 1.099*radius - 0.0288*pow(radius,2) + 0.01662*pow(radius,3) - 0.002384*pow(radius,4) + 8.481E-5*pow(radius,5))));


//  hTotal_livetime->Fill(0.5,livetime-std::max(1.35e-3 - inhibittime, 0.));
//  hTotal_livetime->Fill(0.5,livetime);
  if(IsValidEvent()) {
      hTotal_livetime->Fill(0.5,livetime-std::max(1.35e-3 - inhibittime, 0.));
  }

  if (npulses==3 && (fEdrift!=200 && fEdrift!=0)) has_s3 = HasS3(s2_start_time);
//  if (npulses==3){
//      Bool_t Has_s3 = HasS3(s2_start_time);
//      if (Has_s3!=has_s3) cout<<"diferent S3 identification"<<endl;
//  }
  if (npulses>=2 && (fEdrift!=200 && fEdrift!=0)) total_s1_corr = total_s1*s1_corr_factor(fDriftTimeMax, t_drift);

  const Double_t electron_lifetime = 4733;//update from Guangyong on 02/03/2014//3338.;//updated to the value from Guangyong with adjusted laser window on 01/27/2014//3330; // � 62.4 //3076; //[us]
  total_s2_corr = (run_id>=6553)? total_s2 : total_s2/TMath::Exp(-t_drift/electron_lifetime); // only apply S2 z_correction after run6553.

//  if(total_s1>296. && total_s1<297) cout<<"run: "<<run_id<<" event: "<<event_id<<" total_s1: "<<total_s1<<endl;


  if(IsEventOK()){
    FillHistograms();
  }
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

  hOutF->Save();
  hOutF->Close();
}

Bool_t SladReadSelector::IsValidEvent(){
#define X1      (nchannel==38)                                      // Number of channels                                                                    // Number of channels
#define X2      (baseline_not_found==0)                        // Baseline                                                                    // Baseline
#define X3      (1.35E-3<(livetime+inhibittime))                    // Livetime and inhibit time                                                                    // Livetime and inhibit time
#define X4      (lifetime<1.)                                       // Stalled DAQ                                                                    // Stalled DAQ
#define X5      (veto_present)                                      // Veto present
//  if(X1 && X2 && X3 ) return true;
  if (!X1) {
          hEventConter->Fill(1);
          return false;
  }
  else if (!X2) {
          hEventConter->Fill(2);
          return false;
  }
  else if (!X3) {
          hEventConter->Fill(3);
          return false;
  }
  else return true;
}


Bool_t SladReadSelector::IsEventOK(){
  // Define cuts // from cutdef.h on 09/18/2014
                                                                                                    // Veto present
#define X6      (roi_total_prompt_charge<10.)                                                                                            // Coincident veto signal
#define X7      (lsv_max_window_charge<80.&&lsv_late_window_charge<110.&&wt_charge<200.)                                                   // High energy veto signal
#define X8      (npulses==2||(has_s3==1&&npulses==3))                                                                       // Number of pulses
#define X9      ((run_id<7344 &&-0.25<s1_start_time&&s1_start_time<=-0.15) || (run_id>=7344&&run_id<7641&&-4.10<s1_start_time&&s1_start_time<=-4.00) || (run_id>7640&&-6.10<s1_start_time&&s1_start_time<=-6.00))   // Trigger time
#define X10     (is_saturated_pulse0==0)                                                                                                // S1 saturation
#define X11     (s1_max_frac < max_s1_frac_cut_threshold99)                                                                                                      // S1 maximum fraction
//#define X11     (s1_max_frac<0.4)                                                                                                      // S1 maximum fraction
#define X12     (total_s2_f90 < 0.20 && total_s2 > 30.) //(total_s2_f90<0.20)                                                                                               // Valid S2
#define X13     (100. > total_s2_true)                                                                                                    // Minimal S2
#define X14     ((60. < total_s1_corr)&&(total_s1_corr < 200.))                                                                             // S1 range
//#define X15     ((40.0<t_drift)&&(t_drift<334.5))                                                                                         // Drift time
#define X15     ((fT_drift_cut_min < t_drift)&&( t_drift <=fT_drift_cut_max))                                                                                         // Drift time
//cout<<"fT_drift_cut_min: "<<fT_drift_cut_min<<" fT_drift_cut_max: "<<fT_drift_cut_max<<" t_drift: "<<t_drift<<endl;
//cout<<X15<<" "<<(fT_drift_cut_min < t_drift)<<" "<<( t_drift < fT_drift_cut_max)<<endl;
//	if (!X12 && npulses<2) {
//		  hEventConter->Fill(15);
////		  return false;
//	  }
//	if (!X12) {
//		  hEventConter->Fill(12);
//		  return false;
//	  } else
	if (!X8 && fNRequiredPulses==2) {
	  hEventConter->Fill(8);
	  return false;
  }
  else if (!X9) {
	  hEventConter->Fill(9);
	  return false;
  }
  else if (!X10) {
	  hEventConter->Fill(10);
	  return false;
  }
//  else if (!X11) {
//	  hEventConter->Fill(11);
//	  return false;
//  }
  else if (!X12 && fNRequiredPulses==2) {
	  hEventConter->Fill(12);
	  return false;
  }
  else {
#ifdef JASONS_XY
          Double_t PosiX              = xyl_best_x;//masas_x;//xyl_best_x;
          Double_t PosiY              = xyl_best_y;//masas_y;//xyl_best_y;
#else
          Double_t PosiX              = masas_x;
          Double_t PosiY              = masas_y;
#endif
	  Double_t r2                 = PosiX*PosiX+PosiY*PosiY;

	  hTDrift_All->Fill(t_drift);
	  hTDrift_Max->Fill(t_drift);
	  hR_TDrift_Max->Fill(r2, t_drift);
	  hR_TDrift_Min->Fill(r2, t_drift);

	  hRun_TDrift->Fill(t_drift, run_id);
	  hTotalS1corr_TDrift->Fill(t_drift, total_s1_corr);
	  hTotalF90_TDrift->Fill(t_drift, total_f90);
	  hR_TDrift->Fill(r2, t_drift);

	  if (!X15 && fNRequiredPulses==2){
		  hEventConter->Fill(15);
		  return false;
	  }
	  return true;
  }

}

void SladReadSelector::FillHistograms(){

	  Double_t s2_noncorr_s1_corr = total_s2/total_s1_corr;
#ifdef JASONS_XY
	  Double_t PosiX              = xyl_best_x;//masas_x;//xyl_best_x;
	  Double_t PosiY              = xyl_best_y;//masas_y;//xyl_best_y;
#else
	  Double_t PosiX              = masas_x;
          Double_t PosiY              = masas_y;
#endif
	  Double_t r2                 = PosiX*PosiX+PosiY*PosiY;
	  Float_t  varianceXY         = GetXYVariance();
	  Double_t reco_r             = TMath::Sqrt(r2);
	  Double_t reco_theta         = TMath::ATan2(PosiY, PosiX);

	  Double_t  S2XYcorrfac         = this->XYCorrection4S2(PosiX, PosiY);//1.;

  hTDrift->Fill(t_drift);

  pS2OverS1corr_S1corr_XY->Fill(PosiX, PosiY, total_s1_corr, total_s2/total_s1_corr);
  pS2_S1corr_XY->Fill(PosiX, PosiY, total_s1_corr, total_s2);

  hS2overS1corr_R_S1corr->Fill(total_s1_corr, reco_r, total_s2/total_s1_corr); // s2_over_s1_corr is value before S2 charge normalization

  hS1Total_Ext_vs_F90->Fill(total_s1, total_f90);
  hTotalS1corr_TotalF90->Fill(total_s1_corr, total_f90);


/*  if(total_s1>6.e+3 && total_s1<=50.e+3 && total_f90>0.4 && total_f90<=0.8){ // alpha particle selection
#ifdef OUTTREE
//      fOut_Tree->Fill();
#endif
  }
*/
//  hR_TotalS1_TotalS2->Fill(r, total_s1, total_s2);
  hR_TotalS1_TotalS2->Fill(reco_r, TMath::Log10(total_s1), total_s2);

  hR_TotalS1_TotalS2_corr->Fill(reco_r, TMath::Log10(total_s1), S2XYcorrfac*total_s2);

/*  if(fCutg->IsInside(total_s1, total_f90)){
#ifdef OUTTREE
      fOut_Tree->Fill();
#endif

  }
*/
//  pS2OverS1VsXY->Fill( PosiX,  PosiY, s2_noncorr_s1_corr);// s2_over_s1_corr is value before S2 charge normalization
  hRecPosiXY->Fill(PosiX, PosiY);
  hR_Jason_Masa->Fill(TMath::Sqrt(masas_x*masas_x+masas_y*masas_y), TMath::Sqrt(xyl_best_x*xyl_best_x+xyl_best_y*xyl_best_y));
  hR2_Jason_Masa->Fill(masas_x*masas_x+masas_y*masas_y, xyl_best_x*xyl_best_x+xyl_best_y*xyl_best_y);

//  if(reco_theta>1.5 && reco_theta<=2.4 &&  t_drift>100.)
//    hS2overS1VsR->Fill(reco_r, s2_noncorr_s1_corr);// s2_over_s1_corr); // s2_over_s1_corr is value before S2 charge normalization
  r_hist_wo_norm->Fill(reco_r);
  theta_hist_wo_norm->Fill(reco_theta);
  r_hist->Fill(reco_r);
  theta_hist->Fill(reco_theta);
//  if(reco_theta>1.5 && reco_theta<=2.4)
//  if(reco_theta>1.5 && reco_theta<=2.4)
  r_t_drift_hist              ->Fill( t_drift, reco_r);
  theta_t_drift_hist          ->Fill( t_drift, reco_theta);
  r_total_s1_corr_hist        ->Fill( total_s1_corr, reco_r);
  theta_total_s1_corr_hist    ->Fill( total_s1_corr, reco_theta);
  r2_t_drift_hist             ->Fill( t_drift, r2);

  Double_t LogS2overS1Corr = TMath::Log10(S2XYcorrfac*s2_noncorr_s1_corr);//TMath::Log10(S2XYcorrfac* s2_over_s1_corr);
  hTotalS1vsLogS2overS1->Fill( total_s1_corr, LogS2overS1Corr);
  hTotalS1vsS2overS1->Fill( total_s1_corr, S2XYcorrfac*s2_noncorr_s1_corr);

  hTotalS1->Fill( total_s1);
  hTotalS1vsF90->Fill( total_s1_corr,  total_f90);
  hTotalS1CorrvsDriftTime->Fill( total_s1_corr,  t_drift);
  hTotalS2CorrvsS1->Fill( total_s1,  total_s2_corr);
//  if (reco_r<3)hTotalS2CorrvsS1Corr->Fill( total_s1_corr,  total_s2_corr);
  Float_t x_cut(3.);
  if (PosiX>=0 &&PosiX<x_cut &&PosiY>=0 &&PosiY<x_cut)
	  hTotalS2CorrvsS1Corr->Fill( total_s1_corr,  total_s2_corr);
  hTotalS2XYCorrvsS1Corr->Fill( total_s1_corr, S2XYcorrfac* total_s2_corr);

  if(reco_r<10. &&  t_drift>135. && t_drift<=235.) hTotalS1Corr_Fiducial->Fill( total_s1_corr);

  ht_drift_S2overS1->Fill(t_drift, s2_noncorr_s1_corr);
  ht_drift_S2XYCorroverS1->Fill(t_drift, S2XYcorrfac*s2_noncorr_s1_corr);

  Double_t TBAsym = (total_s1_top-total_s1_bottom)/ total_s1;
  ht_drift_TBAsym->Fill(TBAsym,  t_drift);
  hTotalS1_TBAsym->Fill(TBAsym,  total_s1);
  hDiffTotalS1_TBAsym->Fill(TBAsym, ( total_s1- total_s1_corr)/ total_s1);

  hVarianceXYVsR->Fill(reco_r,  varianceXY);
  hRvsS2BotoverS2Top->Fill(total_s2_bottom/total_s2_top, reco_r);
  hVarianceXYvsS2BotoverS2Top->Fill(total_s2_bottom/total_s2_top,  varianceXY);


  Double_t par[] = {16.6, 0.07669, -2.14e-5, 4.234, -0.02978};
  Double_t x = 200.; // normalized at 200 PE
  Double_t scale_point = par[0]+ par[1]*x+ par[2]*x*x+ TMath::Exp(par[3]+par[4]*x);
  x =  total_s1_corr;
  Double_t scale = par[0]+ par[1]*x+ par[2]*x*x+ TMath::Exp(par[3]+par[4]*x); //mean value of S2/S1 at center at total_s1_corr
  hScaledS2overS1vsS1->Fill( total_s1_corr, s2_noncorr_s1_corr*scale_point/scale);
  hScaledS2XYCorroverS1vsS1->Fill( total_s1_corr, S2XYcorrfac* total_s2_corr/total_s1_corr*scale_point/scale);
  hS2overS1VsR->Fill(reco_r, s2_noncorr_s1_corr*scale_point/scale);// s2_over_s1_corr); // s2_over_s1_corr is value before S2 charge normalization
  pS2OverS1VsXY->Fill(PosiX, PosiY, s2_noncorr_s1_corr*scale_point/scale);// s2_over_s1_corr is value before S2 charge normalization


//  if(S2XYcorrfac*s2_noncorr_s1_corr>scale) hVarianceXYVsR->Fill(reco_r,  varianceXY); // If S2/S1 is fluctuate upward compare to mean value, fill this.
//  if(S2XYcorrfac*s2_noncorr_s1_corr<scale) hVarianceXYVsR->Fill(reco_r,  varianceXY); // If S2/S1 is fluctuate downward compare to mean value, fill this.

 //TBAsymmetry correction
//  Double_t par2[] = {-0.0397956, -0.27216, 0.794036, 1.70427, -3.98323, -8.50783, -2.66051};
//
//  x = TBAsym;
//  Double_t diff_total_s1 = par2[0]+(par2[1]+(par2[2]+(par2[3]+(par2[4]+par2[5]*x)*x)*x)*x)*x; //(total_s1-total_s1_corr)/total_s1
//  Double_t total_s1_TBAcorr =  total_s1*(1.-diff_total_s1);
  Double_t total_s1_TBAcorr =  total_s1*s1_TBAcorr_factor(total_s1_top, total_s1_bottom, total_s1);
  hTotalS1_TBAcorr_TBAsym->Fill(TBAsym, total_s1_TBAcorr);

  Float_t S1corr_mean_Kr_peak = 290.6;
  Float_t S1corr_sigma_Kr_peak = 18.45;
  Float_t nsigma = 3.;
  Float_t s1_min = fS1corr_mean_Kr_peak-nsigma*fS1corr_sigma_Kr_peak;//232.;//S1corr_mean_Kr_peak-nsigma*S1corr_sigma_Kr_peak;
  Float_t s1_max = fS1corr_mean_Kr_peak+nsigma*fS1corr_sigma_Kr_peak;//352.;//S1corr_mean_Kr_peak+nsigma*S1corr_sigma_Kr_peak;

  if(total_s1_corr >= s1_min && total_s1_corr < s1_max) {
      pS2VsXY_Kr->Fill(PosiX, PosiY, total_s2_corr);
      pR2_theta_Kr->Fill(r2, reco_theta, total_s2_corr);
      ht_drift_S2XYCorroverS1_Kr->Fill(t_drift, S2XYcorrfac*s2_noncorr_s1_corr);
      ht_drift_S2XYCorr_Kr->Fill(t_drift, S2XYcorrfac*total_s2);
  }

}

void SladReadSelector::BookHistograms()
{
//  if(Debug)
    Info("SladReadSelector::BookHistograms()","booking histogram...");

    TList  *list = GetOutputList();


    hRunTime = new TH1D("hRunTime", "Time summary [ns];Time Categories", nTimeBins, -0.5, nTimeBins-0.5);
    for (Int_t i=0; i<nTimeBins; i++) hRunTime->GetXaxis()->SetBinLabel(i+1, timelabels[i].Data());
    hRunTime->SetFillColor(kCyan+1);
    list->Add(hRunTime);

    hEventConter = new TH1D("hEventConter", "Event Counter", nEventBins, -0.5, nEventBins-0.5);
//    for (Int_t i=0; i<nEventBins; i++) hEventConter->GetXaxis()->SetBinLabel(i+1, eventlabels[i].Data());
    for (Int_t i=0; i<nEventBins; i++) hEventConter->GetXaxis()->SetBinLabel(i+1, Form("CX%d",i));
    hEventConter->SetFillColor(kCyan+1);
    list->Add(hEventConter);


    Int_t nBinX(100), nBinY(100), nVarianceXY(280);
    Double_t xmin(-20.), xmax(20.), ymin(-20.), ymax(20.), varmin(80.), varmax(220.);
//    Int_t nBinS2(250);
//    Double_t S2min(0.), S2max(0.5);
    Double_t minDriftT(0.), maxDriftT((int)(fDriftTimeMax*1.2));
    Double_t minRunId(TMath::Max(run_id-1000, 0)+0.5), maxRunId(run_id+1000+0.5);
    Int_t nRunId = maxRunId - minRunId;

  Int_t nBinS1(80);
  Double_t S1min(0.), S1max(800);
  Int_t nS1Bins(1000), nS2Bins(1000), nF90Bins(1200), nDriftT(400);//nDriftT((int)(fDriftTimeMax*1.2));
  Float_t minS1(0.), maxS1(1000.), minS2(0.), maxS2(40.e+3), minF90(0.), maxF90(1.2);
//  Float_t minS2overS1(0.), maxS2overS1(1.e+2);



  hTotal_livetime = new TH1F("hTotal_livetime", "total_livetime", 1, 0, 1);
  list->Add(hTotal_livetime);

  hRun_livetime = new TH2D("hRun_livetime", "Run ID vs livetime;livetime [s];run_id", 250, 0., 2.5e-3, nRunId, minRunId, maxRunId);
  list->Add(hRun_livetime);

  hRun_inhibittime = new TH2D("hRun_inhibittime", "Run ID vs inhibittime;inhibittime [s];run_id", 250, 0., 2.5e-3, nRunId, minRunId, maxRunId);
  list->Add(hRun_inhibittime);

  hRun_dT_Trigger = new TH2D("hRun_dT_Trigger", "Run ID vs time from previous trigger (inhibit + live time);inhibittime + livetime [s];run_id", 250, 0., 2.5e-3, nRunId, minRunId, maxRunId);
  list->Add(hRun_dT_Trigger);

  pS2OverS1corr_S1corr_XY = new TProfile3D("pS2OverS1corr_S1corr_XY", "total_s2/total_s1_corr vs total_s1_corr vs xy.;x [cm];y [cm];total_s1_corr [PE]", nBinX, xmin, xmax, nBinY, ymin, ymax, nBinS1, S1min, S1max);
  list->Add(pS2OverS1corr_S1corr_XY);

  pS2_S1corr_XY = new TProfile3D("pS2_S1corr_XY", "total_s2 vs total_s1_corr vs xy.;x [cm];y [cm];total_s1_corr [PE]", nBinX, xmin, xmax, nBinY, ymin, ymax, nBinS1, S1min, S1max);
  list->Add(pS2_S1corr_XY);

  hS2overS1corr_R_S1corr = new TH3D("hS2overS1corr_R_S1corr", "total_s2/total_s1_corr vs Distance from TPC center vs total_s1_corr;total_s1_corr [PE];r [cm];total_s2/total_s1_cor", nBinS1, S1min, S1max, 200, 0., 20, 600, 0, 150);
  list->Add(hS2overS1corr_R_S1corr);

  hTDrift = new TH1D("hTDrift", "Drift time; t_drift [#mus]", 400, minDriftT, maxDriftT);
  list->Add(hTDrift);

  hTDrift_All = new TH1D("hTDrift_All", "Drift time without t_drift cut; t_drift [#mus]", 400, minDriftT, maxDriftT);
  list->Add(hTDrift_All);

  hTDrift_Max = new TH1D("hTDrift_Max", "Drift time without t_drift cut around t_drift_max; t_drift [#mus]",  100, fDriftTimeMax-fT_drift_cut_min, fDriftTimeMax+fT_drift_cut_min);
  list->Add(hTDrift_Max);

  hRun_TDrift = new TH2D("hRun_TDrift", "Run ID vs Drift time; t_drift [#mus];run_id", 400, minDriftT, maxDriftT, nRunId, minRunId, maxRunId);
  list->Add(hRun_TDrift);

  hTotalS1corr_TDrift = new TH2D("hTotalS1corr_TDrift", "Total S1(corr.) vs Drift time; t_drift [#mus];total_s1_corr [PE]", 400, minDriftT, maxDriftT, 400, 0, 400);
  list->Add(hTotalS1corr_TDrift);

  hTotalF90_TDrift = new TH2D("hTotalF90_TDrift", "Total F90 vs Drift time; t_drift [#mus];total_f90", 400, minDriftT, maxDriftT, 110, 0, 1.1);
  list->Add(hTotalF90_TDrift);

  hS1Total_Ext_vs_F90 = new TH2D("hS1Total_Ext_vs_F90", "S1 Spectrum vs f90;total_s1 [PE];total_f90", 35000, 0., 7e+4, 110, 0, 1.1);
  list->Add(hS1Total_Ext_vs_F90);

  hR_TDrift = new TH2D("hR_TDrift", "t_drift vs R^{2};R^{2} [cm];t_drift [#mus]", 100, 0, 400, 400, minDriftT, maxDriftT);
  list->Add(hR_TDrift);

  hR_TDrift_Max = new TH2D("hR_TDrift_Max", "t_drift vs R^{2} around t_drift_max;R^{2} [cm];t_drift [#mus]", 100, 0, 400,  100, fDriftTimeMax-fT_drift_cut_min, fDriftTimeMax+fT_drift_cut_min);
  list->Add(hR_TDrift_Max);

  hR_TDrift_Min = new TH2D("hR_TDrift_Min", "t_drift vs R^{2} at small t_drift;R^{2} [cm];t_drift [#mus]", 100, 0, 400, 100, 0, 2.*fT_drift_cut_min);
  list->Add(hR_TDrift_Min);

  hTotalS1corr_TotalF90 = new TH2D("hTotalS1corr_TotalF90", "F90 vs S1 (corrected for z-dependence); total_s1_corr [PE]; total_f90", 1000, 0, 1000, 110, 0, 1.1);
  list->Add(hTotalS1corr_TotalF90);

//  hR_TotalS1_TotalS2 = new TH3D("hR_TotalS1_TotalS2", "total_s2 (S3) vs total_s1 (S2) vs R;r [cm];total_s1 (S2) [PE];total_s2 (S3) [PE]", 6, 0., 18., 10, 0., 5000., 160, -10., 150);
  hR_TotalS1_TotalS2 = new TH3D("hR_TotalS1_TotalS2", "total_s2 (S3) vs total_s1 (S2) vs R;r [cm];Log10(total_s1) (S2) [PE];total_s2 (S3) [PE]", 7, 0., 21., 10,-1, +4., 160, -10., 150);
  list->Add(hR_TotalS1_TotalS2);

  hR_TotalS1_TotalS2_corr = new TH3D("hR_TotalS1_TotalS2_corr", "total_s2_corr (S3, XY corrected) vs total_s1 (S2) vs R;r [cm];Log10(total_s1) (S2) [PE];total_s2_xycorr (S3) [PE]", 7, 0., 21., 10,-1, +4., 160, -10., 150);
  list->Add(hR_TotalS1_TotalS2_corr);

  hRecPosiXY = new TH2D("hRecPosiXY", "Reconstructed Position in x-y ; x [cm]; y [cm]", nBinX, xmin, xmax, nBinY, ymin, ymax);
  list->Add(hRecPosiXY);
  hR_Jason_Masa = new TH2D("hR_Jason_Masa", "Reconstructed Position in R Jason's vs Masa ; R_{Masa} [cm]; R_{Jason} [cm]", nBinX, 0, xmax, nBinY, 0, ymax);
  list->Add(hR_Jason_Masa);
  hR2_Jason_Masa = new TH2D("hR2_Jason_Masa", "Reconstructed Position in R^{2} Jason's vs Masa ; R^{2}_{Masa} [cm]; R^{2}_{Jason} [cm]", 100, 0, 400, 100, 0, 400);
  list->Add(hR2_Jason_Masa);

  hS2overS1VsR = new TH2D("hS2overS1VsR", "Total S2/S1 vs Distance from TPC center  (corrected for z-dependence); r [cm]; total_s2_corr/total_s1_cor", 200, 0., xmax, 500, 0, 100);
   list->Add(hS2overS1VsR);

  r_hist_wo_norm              = new TH1F("r_hist_wo_norm", "radial distribution; r [cm]", 40, 0, 20);
  list->Add(r_hist_wo_norm);
  theta_hist_wo_norm          = new TH1F("theta_hist_wo_norm", "azimuthal distribution; #theta [rad]", 40,-3.5, 3.5);// -TMath::Pi(), TMath::Pi());
  list->Add(theta_hist_wo_norm);
  r_hist                      = new TH1F("r_hist", "r; r [cm]", 20, 0, 20);
  list->Add(r_hist);
  theta_hist                  = new TH1F("theta_hist", "theta ; #theta [rad]", 20, -3.5, 3.5);//-TMath::Pi(), TMath::Pi());
  list->Add(theta_hist);
  r_t_drift_hist              = new TH2F("r_t_drift_hist", "r vs t_drift;t_drift [#mus]; r [cm]", 100, minDriftT, maxDriftT, 100, 0, 20);
  list->Add(r_t_drift_hist);
  theta_t_drift_hist          = new TH2F("theta_t_drift_hist", "theta vs t_drift;t_drift [#mus]; #theta [rad]", 40, minDriftT, maxDriftT, 100,-3.5, 3.5);// -TMath::Pi(), TMath::Pi());
  list->Add(theta_t_drift_hist);
  r_total_s1_corr_hist        = new TH2F("r_total_s1_corr_hist", "r vs total_s1_corr;total_s1_corr [PE]; r [cm]", 200, 0., 1000., 20, 0, 20);
  list->Add(r_total_s1_corr_hist);
  theta_total_s1_corr_hist    = new TH2F("theta_total_s1_corr_hist", "theta vs total_s1_corr;total_s1_corr [PE];#theta [rad]", 200, 0., 1000., 20, -3.5, 3.5);//-TMath::Pi(), TMath::Pi());
  list->Add(theta_total_s1_corr_hist);

  r2_t_drift_hist             = new TH2F("r2_t_drift_hist", "r^{2} vs t_drift;t_drift [#mus]; r^{2} [cm]", 100, minDriftT, maxDriftT, 100, 0, 400);
  list->Add(r2_t_drift_hist);


  pS2OverS1VsXY = new TProfile2D("pS2OverS1VsXY", "S2/S1 vs xy.; x [cm]; y [cm]; total_s2/total_s1_corr", nBinX, xmin, xmax, nBinY, ymin, ymax);//, 500, 0., 1.);
  list->Add(pS2OverS1VsXY);

  pS2VsXY_Kr = new TProfile2D("pS2VsXY_Kr", "Averaged S2 vs xy. around Kr peak; x [cm]; y [cm]; total_s2_corr [PE]", 14, -21, 21, 14, -21, 21);//, 20, xmin, xmax, 20, ymin, ymax);//, 500, 0., 1.);
//  pS2VsXY_Kr = new TProfile2D("pS2VsXY_Kr", "Averaged S2 vs xy. around Kr peak; x [cm]; y [cm]; total_s2_corr [PE]", 40, -20, 20, 40, -20, 20);//, 20, xmin, xmax, 20, ymin, ymax);//, 500, 0., 1.);
  list->Add(pS2VsXY_Kr);

  pR2_theta_Kr             = new TProfile2D("pR2_theta_Kr", "theta vs r^{2};r^{2} [cm];#theta [rad]", 20, 0, 400, 10,-TMath::Pi(), TMath::Pi());
  list->Add(pR2_theta_Kr);

  hTotalS1vsLogS2overS1 = new TH2D("hTotalS1vsLogS2overS1", "log_{10}(S2_{xy corr.}/S1_{corr.}) vs S1_{corr.}; total_s1_corr [PE]; log_{10}(total_s2_xy_corr/total_s1_corr)", 200, 0., 1000., 200, -1, 3);
  list->Add(hTotalS1vsLogS2overS1);

  hTotalS1vsS2overS1 = new TH2D("hTotalS1vsS2overS1", "S2_{xy corr.}/S1_{corr.} vs S1_{corr.}; total_s1_corr [PE]; total_s2_xy_corr/total_s1_corr", 600, 0., 3000., 200, 0, 200);
  list->Add(hTotalS1vsS2overS1);


  hTotalS1 = new TH1D("hTotalS1", "S1; total_s1 [PE]", nS1Bins, minS1, maxS1);
  list->Add(hTotalS1);
  hTotalS1vsF90 = new TH2D("hTotalS1vsF90", "F90 vs S1; total_s1_corr [PE]; total_f90", nS1Bins, minS1, maxS1, nF90Bins, minF90, maxF90);
  list->Add(hTotalS1vsF90);
  hTotalS1CorrvsDriftTime = new TH2D("hTotalS1CorrvsDriftTime", "Drift time vs S1_{Corr.}; total_s1_corr [PE]; t_drift [#mu s]", nS1Bins, minS1, maxS1, nDriftT, minDriftT, maxDriftT);
  list->Add(hTotalS1CorrvsDriftTime);
  hTotalS2CorrvsS1 = new TH2D("hTotalS2CorrvsS1", "S2_{Corr.} vs S1; total_s1 [PE]; total_s2_corr [PE]", nS1Bins, minS1, maxS1, nS2Bins, minS2, maxS2);
  list->Add(hTotalS2CorrvsS1);
  hTotalS2CorrvsS1Corr = new TH2D("hTotalS2CorrvsS1Corr", "S2_{Corr.} vs S1_{Corr.}; total_s1_corr [PE]; total_s2_corr [PE]", nS1Bins, minS1, maxS1, nS2Bins, minS2, maxS2);
  list->Add(hTotalS2CorrvsS1Corr);
  hTotalS2XYCorrvsS1Corr = new TH2D("hTotalS2XYCorrvsS1Corr", "S2_{XY Corr.} vs S1_{Corr.}; total_s1_corr [PE]; total_s2_xy_corr [PE]", nS1Bins, minS1, maxS1, nS2Bins, minS2, maxS2);
  list->Add(hTotalS2XYCorrvsS1Corr);


  Int_t nEnergyBins(700);
  Double_t minEne(0.), maxEne(3500.);
  hTotalS1CorrvsEnergy = new TH2D("hTotalS1CorrvsEnergy", "S1_{Corr.} vs Energy; energy [keV]; total_s1_corr [PE]", nEnergyBins, minEne, maxEne, nS1Bins, minS1, maxS1);
  list->Add(hTotalS1CorrvsEnergy);
  hS1LYvsEnergy = new TH2D("hS1LYvsEnergy", "S1_{Corr.} Light Yield vs Energy; energy [keV]; S1_{Corr.} Light Yield [PE/keV]", nEnergyBins, minEne, maxEne, 500, 0, 10);
  list->Add(hS1LYvsEnergy);
  hTotalS2XYCorrvsEnergy = new TH2D("hTotalS2XYCorrvsEnergy", "S2_{XY Corr.} vs Energy; energy [keV]; total_s2_xy_corr [PE]", nEnergyBins, minEne, maxEne, nS2Bins, minS2, maxS2);
  list->Add(hTotalS2XYCorrvsEnergy);
  hS2XYCorrLYvsEnergy = new TH2D("hS2XYCorrLYvsEnergy", "S2_{XY Corr.} Light Yield vs Energy; energy [keV]; S2_{XY Corr.} Light Yield [PE/keV]", nEnergyBins, minEne, maxEne, 500, 0, 1000);
  list->Add(hS2XYCorrLYvsEnergy);

  hS2overS1vsEnergy = new TH2D("hS2overS1vsEnergy", "S2_{XY corr.}/S1_{corr.} vs Energy; energy [keV]; total_s2_xy_corr/total_s1_corr", nEnergyBins, minEne, maxEne, 200, 0, 200);
  list->Add(hS2overS1vsEnergy);

  hS2overS1vsS1_Center = new TH2D("hS2overS1vsS1_Center", "S2_{corr.}/S1_{corr.} vs S1_{Corr.} w/ max_S2_ch=30; total_s1_corr [PE]; total_s2_corr/total_s1_corr", nS1Bins, minS1, maxS1, 200, 0, 200);
  list->Add(hS2overS1vsS1_Center);

  hScaledS2overS1vsS1 = new TH2D("hScaledS2overS1vsS1", "Scaled S2/S1_{corr.} vs S1_{Corr.}; total_s1_corr [PE]; Scaled total_s2/total_s1_corr", nS1Bins, minS1, maxS1, 200, 0, 200);
  list->Add(hScaledS2overS1vsS1);

  hScaledS2XYCorroverS1vsS1 = new TH2D("hScaledS2XYCorroverS1vsS1", "Scaled S2_{xy corr.}/S1_{corr.} vs S1_{Corr.}; total_s1_corr [PE]; Scaled total_s2_xy_corr/total_s1_corr", nS1Bins, minS1, maxS1, 200, 0, 200);
  list->Add(hScaledS2XYCorroverS1vsS1);

  ht_drift_S2overS1 = new TH2D("ht_drift_S2overS1", "S2/S1_{corr.} vs t_drift; t_drift [#mus]; total_s2/total_s1_corr", 100, minDriftT, maxDriftT, 500, 0, 100);
  list->Add(ht_drift_S2overS1);

  ht_drift_S2XYCorroverS1 = new TH2D("ht_drift_S2XYCorroverS1", "S2_{XYcorr. not z}/S1_{corr.} vs t_drift; t_drift [#mus]; total_s2_XYcorr_not_z/total_s1_corr", 100, minDriftT, maxDriftT, 500, 0, 100);
  list->Add(ht_drift_S2XYCorroverS1);

  ht_drift_S2XYCorroverS1_Kr = new TH2D("ht_drift_S2XYCorroverS1_Kr", "S2_{XYcorr. not z}/S1_{corr.} vs t_drift around Kr peak; t_drift [#mus]; total_s2_XYcorr_not_z/total_s1_corr", 100, minDriftT, maxDriftT, 500, 0, 100);
  list->Add(ht_drift_S2XYCorroverS1_Kr);

  ht_drift_S2XYCorr_Kr = new TH2D("ht_drift_S2XYCorr_Kr", "S2_{XYcorr. not z} vs t_drift around Kr peak; t_drift [#mus]; total_s2_XYcorr_not_z [PE]", 100, minDriftT, maxDriftT, 300, 0, 30.e3);
  list->Add(ht_drift_S2XYCorr_Kr);

  ht_drift_TBAsym = new TH2D("ht_drift_TBAsym", "t_drift vs Top Bottom Asymmetry (total_s1_{top}-total_s1_{bot.})/total_s1; TBAsym; t_drift [#mus]", 200, -0.5, 0.3, 100, minDriftT, maxDriftT);
  list->Add(ht_drift_TBAsym);

  hTotalS1_TBAsym = new TH2D("hTotalS1_TBAsym", "total_s1 vs Top Bottom Asymmetry (total_s1_{top}-total_s1_{bot.})/total_s1; TBAsym; total_s1 [PE]", 200, -0.5, 0.3, nS1Bins, minS1, maxS1);
  list->Add(hTotalS1_TBAsym);

  hDiffTotalS1_TBAsym = new TH2D("hDiffTotalS1_TBAsym", "(total_s1-total_s1_corr)/total_s1 vs Top Bottom Asymmetry (total_s1_{top}-total_s1_{bot.})/total_s1; TBAsym; (total_s1-total_s1_corr)/total_s1", 200, -0.5, 0.3, 200, -0.1, 0.1);
  list->Add(hDiffTotalS1_TBAsym);

  hTotalS1_TBAcorr_TBAsym = new TH2D("hTotalS1_TBAcorr_TBAsym", "total_s1_TBAcorr vs Top Bottom Asymmetry (total_s1_{top}-total_s1_{bot.})/total_s1; TBAsym; total_s1_TBAcorr [PE]", 200, -0.5, 0.3, nS1Bins, minS1, maxS1);
  list->Add(hTotalS1_TBAcorr_TBAsym);

  hVarianceXYVsR = new TH2D("hVarianceXYVsR", "x-y variance vs Distance from TPC center; R [cm]; x-y variance", 300, 0., xmax, nVarianceXY, varmin, varmax);
  list->Add(hVarianceXYVsR);


  hRvsS2BotoverS2Top = new TH2D("hRvsS2BotoverS2Top", "Distance from TPC center vs S2_{Bottom PMT}/ S2_{TOP PMT} z-corrected; S2_{BOT. PMT}/ S2{TOP PMT} (z-corrected); R [cm]", 200, 0, 1., 300, 0., xmax);
  list->Add(hRvsS2BotoverS2Top);

  hVarianceXYvsS2BotoverS2Top = new TH2D("hVarianceXYvsS2BotoverS2Top", "x-y variance vs S2_{Bottom PMT}/ S2_{TOP PMT} z-corrected; S2_{BOT. PMT}/ S2{TOP PMT} (z-corrected); x-y variance", 200, 0, 1., nVarianceXY, varmin, varmax);
  list->Add(hVarianceXYvsS2BotoverS2Top);


  pdrVsR = new TProfile("pdrVsR", "Distance between rec. and prev. rec. positions; R [cm];dr", 36, 0, 18);
  list->Add(pdrVsR);

  hTotalS1Corr_Fiducial = new TH1D("hTotalS1Corr_Fiducial", "S1 with fiducial cut; total_s1_corr [PE]", nS1Bins, minS1, maxS1);
  list->Add(hTotalS1Corr_Fiducial);
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

void SladReadSelector::SetDriftFieldConf(Int_t Edrift, Int_t RunId)
{
  fIsFieldSet = true;
  fEdrift = Edrift;
  fNRequiredPulses=2;
  if(Edrift == 0){
      fDriftTimeMax = 10.; // usually acquisition window is 25 us
      fNRequiredPulses=1; //default is 2
      fS1corr_mean_Kr_peak = (RunId<7909)? 326.7 : 321.6;
      fS1corr_sigma_Kr_peak = 20.;

  } else if (Edrift == 50){
      fDriftTimeMax = 1250 ;
      fS1corr_mean_Kr_peak = 308.1;
      fS1corr_sigma_Kr_peak = 20.;

  } else if (Edrift == 100){
      fDriftTimeMax = (RunId<7909)? 650.4 : 665.;
      fS1corr_mean_Kr_peak = (RunId<7909)? 303.3 : 300.6;
      fS1corr_sigma_Kr_peak = 20.;

  } else if (Edrift == 150){
      fDriftTimeMax = (RunId<7909)? 464.3 : 469.;
      fS1corr_mean_Kr_peak = (RunId<7909)? 299.3 : 294.3;
      fS1corr_sigma_Kr_peak = 20.;

  } else if (Edrift == 200){
      fDriftTimeMax = (RunId<7909)? 373.3 : 376.;
      fS1corr_mean_Kr_peak = (RunId<7909)? 292.9 : 288.8;
      fS1corr_sigma_Kr_peak = 20.;

  } else {
      cout<<"Unknown Drift field value!! Please check input value. E_drift: "<<Edrift<<endl;
  }
  if(Edrift !=0){
      const Float_t t_drift_default = 373.3;
      const Double_t t_s3_sep_min_ratio = 372./t_drift_default; //[us] // scaled by 200 V/cm value
      const Double_t t_s3_sep_max_ratio = 400./t_drift_default; //[us]
      fS3SepTimeMin = t_s3_sep_min_ratio*fDriftTimeMax;
      fS3SepTimeMax = t_s3_sep_max_ratio*fDriftTimeMax;
      cout<<"fS3SepTimeMin: "<<fS3SepTimeMin<<" us, fS3SepTimeMax: "<<fS3SepTimeMax<<" us"<<endl;

      Float_t delta_t = 10./t_drift_default*fDriftTimeMax;
      fT_drift_cut_min = delta_t;
      fT_drift_cut_max = fDriftTimeMax - delta_t;
      cout<<"t_drift_cut_min: "<<fT_drift_cut_min<<" us, T_drift_cut_max: "<<fT_drift_cut_max<<" us"<<endl;
  }
  cout<<"======== XYRecSelector::SetDriftFieldConf() =========="<<endl;
  cout<<"Drift field is "<<Edrift<< " and max drift time is "<< fDriftTimeMax <<endl;
  cout<<"Number of required pulses is "<<fNRequiredPulses<<endl;
}

Bool_t SladReadSelector::CheckFieldSetting(Float_t acqui_wind){
  if (acqui_wind<fDriftTimeMax) return false;
  else return true;
}

Int_t SladReadSelector::GetFieldFromAcquisitionWindow(Float_t acqui_wind){
  if(acqui_wind<100.){ // null field
      return 0;
  } else if (acqui_wind<480.){ // 200 V/cm
      return 200;
  } else if (acqui_wind<600.){ // 150 V/cm
      return 150;
  } else if (acqui_wind<800.){ // 100 V/cm
      return 100;
  } else if (acqui_wind<2000.){ // 50 V/cm
      return 50;
  } else {
      cout<<"Too long acquisition window: "<<acqui_wind<<" us. Please check the data."<<endl;
      return -1;
  }
}

Float_t SladReadSelector::GetExpectedAcquisitionWindow(Int_t Edrift){
  switch (Edrift){ // null field
  case 0:
    return 25.;
  case 200:
    return 440.;
  case 150:
    return 530.;
  case 100:
    return 720.;
  case 50:
    return 1500.;
  default:
    cout<<"unknown drift field: "<<Edrift<<" V/cm."<<endl;
    return -1.;
  }
}


Bool_t SladReadSelector::HasS3(Float_t s2_starttime){
  if(s2_starttime!=pulse_info_start_time[1]) cout<<"s2_starttime doesn't match. Please check it."<<endl;

  Double_t t_drift2to3 = pulse_info_start_time[2] - s2_starttime;
  if (t_drift2to3 > fS3SepTimeMin && t_drift2to3 <= fS3SepTimeMax) return true;
  else return false;

}

Double_t SladReadSelector::s1_corr_factor(Double_t t_drift_max, Double_t t_drift)
{
  Double_t z = t_drift/(0.5*t_drift_max); // note normalization is to 0.5*t_drift_max
  // looked at Kr peaks in 15us t_drift windows (Run5330+5340), and fit these to [0]*z^5 + [1]*z^4 + [2]*z^3+[3]*z^2+[4]*z+[5].
  Double_t fit_par0 = 0.0407;
  Double_t fit_par1 = -0.206;
  Double_t fit_par2 = 0.407;
  Double_t fit_par3 = -0.389;
  Double_t fit_par4 = 0.247;
  Double_t fit_par5 = 0.898;
  // normalizing all points on fitted curve to expected Kr peak at t_drift_max/2
  Double_t exp_Kr_peak_at_half_t_drift_max = fit_par0 + fit_par1 + fit_par2 + fit_par3 + fit_par4 + fit_par5;
  Double_t exp_Kr_peak_at_t_drift = fit_par0*z*z*z*z*z + fit_par1*z*z*z*z + fit_par2*z*z*z + fit_par3*z*z + fit_par4*z + fit_par5;
  return exp_Kr_peak_at_half_t_drift_max/exp_Kr_peak_at_t_drift; // s1 correction factor
}

Double_t SladReadSelector::s1_TBAcorr_factor(Double_t s1_top, Double_t s1_bottom, Double_t s1)
{
  //TBAsymmetry correction
   Double_t par2[] = {-0.0397956, -0.27216, 0.794036, 1.70427, -3.98323, -8.50783, -2.66051};
   Double_t x = (s1_top-s1_bottom)/ s1;// TBAsym;
   Double_t diff_total_s1 = par2[0]+(par2[1]+(par2[2]+(par2[3]+(par2[4]+par2[5]*x)*x)*x)*x)*x; //(total_s1-total_s1_corr)/total_s1
//   Double_t total_s1_TBAcorr =  s1*(1.-diff_total_s1);
   return (1.-diff_total_s1);
}

Float_t SladReadSelector::GetXYVariance(){
  Double_t bary_x(0.), bary_y(0.), bary_x_sq(0.), bary_y_sq(0.);
  Double_t top_s2(0.), top_s2_sq(0.);
  for (Int_t ch = PMTGeom::N_CHANNELS/2; ch < PMTGeom::N_CHANNELS; ch++) {
    Double_t s2 =total_s2*s2_ch_frac[ch];

    top_s2 += s2;
    top_s2_sq += s2*s2;
    bary_x += PMTGeom::pmtUnit*PMTGeom::pmt_x[ch] * s2;
    bary_y += PMTGeom::pmtUnit*PMTGeom::pmt_y[ch] * s2;
    bary_x_sq += PMTGeom::pmtUnit*PMTGeom::pmtUnit*PMTGeom::pmt_x[ch]*PMTGeom::pmt_x[ch] * s2;
    bary_y_sq += PMTGeom::pmtUnit*PMTGeom::pmtUnit*PMTGeom::pmt_y[ch]*PMTGeom::pmt_y[ch] * s2;
  }// end loop over channels

  //calculate unbiased weighted variance of xy
  Double_t denom = top_s2*top_s2 - top_s2_sq;
  Double_t variance_x = (top_s2*bary_x_sq - bary_x*bary_x)/denom;
  Double_t variance_y = (top_s2*bary_y_sq - bary_y*bary_y)/denom;

//  Double_t biased_variance_x = (top_s2*bary_x_sq - bary_x*bary_x)/top_s2/top_s2;

//  cout<<"unbiased variance_x: "<<variance_x<<endl;//TMath::Sqrt(variance_x)<<endl;
//  cout<<"unbiased variance_y: "<<variance_y<<endl;//<<TMath::Sqrt(variance_y)<<endl;
//  cout<<"biased variance_x: "<<biased_variance_x<<endl;//TMath::Sqrt(biased_variance_x)<<endl;

  return variance_x + variance_y;


}

Double_t SladReadSelector::XYCorrection4S2(Double_t x, Double_t y) {
#ifdef CorrXY
	return 1./fhS2Correction_factor->GetBinContent(fhS2Correction_factor->FindBin(x, y));
#else
	return 1./fhS2Correction_factor->GetBinContent(fhS2Correction_factor->FindBin(x*x+y*y, TMath::ATan2(y, x)));
#endif
}

