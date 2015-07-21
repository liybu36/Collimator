//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Aug  2 11:48:34 2014 by ROOT version 5.99/06
// from TTree SladReadSelector/SladReadSelector
// found on file: period_I+III_SLAD-1.4_run-ADAM_veto-v6.root
//////////////////////////////////////////////////////////

#ifndef SladReadSelector_h
#define SladReadSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TString.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include "TProfile2D.h"
#include "TProfile3D.h"
#include "TCutG.h"

//#include "../XYReconstructor/XYReconstructor.h"

// Header file for the classes stored in the TTree if any.

class SladReadSelector : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
//   XYReconstructor *fMyFcn;
   Bool_t          fIsFieldSet;         //flag for field setting
   Float_t         fEdrift;             // drift field in V/cm
   Float_t         fDriftTimeMax;       // maximum drift time
   Float_t         fT_drift_cut_min;       // drift time minimum cut value
   Float_t         fT_drift_cut_max;       // drift time maximum cut value
   Float_t         fS3SepTimeMin;       // drift time minimum S3 window
   Float_t         fS3SepTimeMax;       // drift time maximum S3 window

   Int_t           fNRequiredPulses;    // number of required pulses
//   static const TString timelabels[];
//   static const TString eventlabels[];

   Float_t         fS1corr_mean_Kr_peak;
   Float_t         fS1corr_sigma_Kr_peak;

   Float_t         fPrevious_total_s1;  // this is check for bug


// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run_id;
   Int_t           subrun_id;
   Int_t           event_id;
   Double_t        livetime;
   Double_t        inhibittime;
   Float_t         acqui_window;
   Int_t           nchannel;
   short           baseline_not_found;
   Int_t           npulses;
   Int_t           has_s3;

   Float_t         total_s1;
   Float_t         total_s1_corr;
   Float_t         total_s1_top;
   Float_t         total_s1_bottom;
   Float_t         total_f90;
   Float_t         s1_start_time;
   Float_t         s1_end_time;
   Int_t           is_saturated_pulse0;

   Float_t         total_s2;
   Float_t         total_s2_corr;
   Float_t         total_s2_top;
   Float_t         total_s2_bottom;
   Float_t         total_s2_f90;
   Float_t         s2_start_time;
   Float_t         s2_end_time;
   Int_t           is_saturated_pulse1;

   Double_t        t_drift;
   Float_t         s1_max_frac;

   Short_t         selection_results[16];

   Float_t         masas_x;
   Float_t         masas_y;
   Float_t         masas_chi2;

   Float_t         xyl_SCM;
   Float_t         xyl_best_x;
   Float_t         xyl_best_y;
   Float_t         xyl_best_chi2;

//   Float_t         roi_total_prompt_charge;
//   Float_t         lsv_max_window_charge;
//   Float_t         lsv_late_window_charge;
//   Float_t         wt_charge;

   Float_t         max_s1_frac_cut_threshold99;
   Int_t           max_s1_frac_cut_exceeds99;

   Int_t           s2_max_chan;
   Float_t         s2_max_frac;
   Float_t         s2_ch_frac[38];

   Float_t         radius;

   Int_t           pulse_info_npulses;
   Int_t           pulse_info_pulse_id[22];   //[pulse_info_npulses]
   Float_t         pulse_info_start_time[22];   //[pulse_info_npulses]



   // List of branches
   TBranch        *b_run_id;   //!
   TBranch        *b_subrun_id;   //!
   TBranch        *b_event_id;   //!
   TBranch        *b_livetime;   //!
   TBranch        *b_inhibittime;   //!
   TBranch        *b_acqui_window;   //!
   TBranch        *b_nchannel;   //!
   TBranch        *b_SumChannelHasNoBaseline;   //!
   TBranch        *b_n_phys_pulses;   //!
   TBranch        *b_has_s3;   //!

   TBranch        *b_total_s1;   //!
   TBranch        *b_total_s1_corr;   //!
   TBranch        *b_total_s1_top;   //!
   TBranch        *b_total_s1_bottom;   //!
   TBranch        *b_total_f90;   //!
   TBranch        *b_is_saturated_pulse0;   //!
   TBranch        *b_total_s2;   //!
   TBranch        *b_total_s2_corr;   //!
   TBranch        *b_total_s2_top;   //!
   TBranch        *b_total_s2_bottom;   //!
   TBranch        *b_total_s2_f90;   //!
   TBranch        *b_is_saturated_pulse1;   //!
   TBranch        *b_s1_start_time;   //!
   TBranch        *b_s1_end_time;   //!
   TBranch        *b_s2_start_time;   //!
   TBranch        *b_s2_end_time;   //!

   TBranch        *b_tdrift;   //!
//   TBranch        *b_selection_results;   //!

   TBranch        *b_max_s1_frac_cut_threshold;   //!
   TBranch        *b_max_s1_frac_cut_exceeds;   //!

   TBranch        *b_selection_results;   //!

   TBranch        *b_s2_max_chan;   //!
   TBranch        *b_s2_max_frac;   //!
   TBranch        *b_s2_chan;   //!

   TBranch        *b_masas_x;   //!
   TBranch        *b_masas_y;   //!
   TBranch        *b_masas_chi2;   //!

   TBranch        *b_xyl_SCM;   //!
   TBranch        *b_xyl_best_x;   //!
   TBranch        *b_xyl_best_y;   //!
   TBranch        *b_xyl_best_chi2;   //!

   TBranch        *b_pulse_info_npulses;   //!
   TBranch        *b_pulse_info_pulse_id;   //!
   TBranch        *b_pulse_info_start_time;   //!


   SladReadSelector(TTree * /*tree*/ =0) :
     fChain(0),
     fIsFieldSet(false), fEdrift(-1), fDriftTimeMax(-1), fT_drift_cut_min(-1), fT_drift_cut_max(-1), fS3SepTimeMin(-1), fS3SepTimeMax(-1), fNRequiredPulses(-1),
     fPrevious_total_s1(-1),
     fOutName(""),
     fOut_Tree(0),fCutg(0),fProfileName(""), fhS2Correction_factor(0),
     hRunTime(0), hEventConter(0),
     pS2OverS1corr_S1corr_XY(0),
     pS2_S1corr_XY(0),
     hS2overS1corr_R_S1corr(0),
     hRun_TDrift(0), hTotalS1corr_TDrift(0), hTotalF90_TDrift(0), hRun_dT_Trigger(0),
     hRun_livetime(0), hRun_inhibittime(0),
     hTDrift(0), hTDrift_All(0), hTDrift_Max(0), hTotalS1corr_TotalF90(0), hS1Total_Ext_vs_F90(0),
     hR_TDrift(0), hR_TDrift_Max(0), hR_TDrift_Min(0),
     hTotal_livetime(0), hR_TotalS1_TotalS2(0), hR_TotalS1_TotalS2_corr(0),
     hRecPosiXY(0), hR_Jason_Masa(0), hR2_Jason_Masa(0), hS2overS1VsR(0), hTotalS1vsS2overS1(0),
     hTotalS1vsLogS2overS1(0),
     r_hist(0), theta_hist(0), r_hist_wo_norm(0), theta_hist_wo_norm(0),
     r_t_drift_hist(0), theta_t_drift_hist(0), r_total_s1_corr_hist(0), theta_total_s1_corr_hist(0), r2_t_drift_hist(0),
     pS2OverS1VsXY(0), pS2VsXY_Kr(0), pR2_theta_Kr(0),
     hTotalS1(0),
     hTotalS1vsF90(0), hTotalS1CorrvsDriftTime(0), hTotalS2CorrvsS1(0), hTotalS2CorrvsS1Corr(0), hTotalS2XYCorrvsS1Corr(0),
     hTotalS1CorrvsEnergy(0), hS1LYvsEnergy(0), hTotalS2XYCorrvsEnergy(0), hS2XYCorrLYvsEnergy(0), hS2overS1vsEnergy(0), ht_drift_S2overS1(0), ht_drift_S2XYCorroverS1(0),
	 ht_drift_ScaledS2XYCorroverS1_Kr(0), ht_drift_S2XYCorroverS1_Kr(0), ht_drift_S2XYCorr_Kr(0), ht_drift_S2AllCorr_Kr(0), ht_drift_S2overS1_center(0),
     ht_drift_TBAsym(0),hTotalS1_TBAsym(0), hDiffTotalS1_TBAsym(0), hTotalS1_TBAcorr_TBAsym(0),
     hVarianceXYVsR(0), hRvsS2BotoverS2Top(0), hVarianceXYvsS2BotoverS2Top(0), hS2overS1vsS1_Center(0), hScaledS2overS1vsS1(0), hScaledS2XYCorroverS1vsS1(0),
     pdrVsR(0),
     hTotalS1Corr_Fiducial(0)

   { }
   virtual ~SladReadSelector() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();

   void SetOutputFile(TString name){fOutName = name;}
   void SetProfileName(TString fProfile) { fProfileName = fProfile; }

private:
   TString fOutName;

   TTree      *fOut_Tree; //!

   TCutG      *fCutg;      //!

   TString fProfileName; //! for XY reconstructor

   TH2D *fhS2Correction_factor; // S2 xy correction factor

   Bool_t          IsValidEvent();
   Bool_t          IsEventOK_Default();
   Bool_t          IsEventOK();
   void            FillHistograms();
   void            BookHistograms();
   void            SetCutG4Outliers();
   void            SetDriftFieldConf(Int_t Edrift, Int_t RunId);
   Double_t        SetElectronLifetime(Int_t RunId);
   Bool_t          CheckFieldSetting(Float_t acqui_window);
   Int_t           GetFieldFromAcquisitionWindow(Float_t acqui_window);
   Float_t         GetExpectedAcquisitionWindow(Int_t Edrift);
   Bool_t          HasS3(Float_t s2_starttime);
   Double_t        s1_corr_factor(Double_t t_drift_max, Double_t t_drift);
   Double_t        s1_TBAcorr_factor(Double_t s1_top, Double_t s1_bottom, Double_t s1);
   Float_t         GetXYVariance();
   Double_t        XYCorrection4S2(Double_t x, Double_t y);

   TH1D* hRunTime, *hEventConter;

   TProfile3D *pS2OverS1corr_S1corr_XY, *pS2_S1corr_XY;
   TH3D *hS2overS1corr_R_S1corr;
   TH2D *hRun_TDrift, *hTotalS1corr_TDrift, *hTotalF90_TDrift, *hRun_dT_Trigger;
   TH2D *hRun_livetime, *hRun_inhibittime;
   TH1D *hTDrift, *hTDrift_All, *hTDrift_Max;
   TH2D *hTotalS1corr_TotalF90, *hS1Total_Ext_vs_F90;
   TH2D *hR_TDrift, *hR_TDrift_Max, *hR_TDrift_Min;
   TH1F* hTotal_livetime;

   TH3D *hR_TotalS1_TotalS2, *hR_TotalS1_TotalS2_corr;

   TH2D *hRecPosiXY, *hR_Jason_Masa, *hR2_Jason_Masa;
   TH2D *hS2overS1VsR; //!
   TH2D *hTotalS1vsS2overS1; //!
   TH2D *hTotalS1vsLogS2overS1;


   TH1F *r_hist, *theta_hist, *r_hist_wo_norm, *theta_hist_wo_norm; //!
   TH2F *r_t_drift_hist, *theta_t_drift_hist, *r_total_s1_corr_hist, *theta_total_s1_corr_hist, *r2_t_drift_hist; //!

   TProfile2D *pS2OverS1VsXY, *pS2VsXY_Kr, *pR2_theta_Kr;

   TH1D *hTotalS1;
   TH2D *hTotalS1vsF90, *hTotalS1CorrvsDriftTime, *hTotalS2CorrvsS1, *hTotalS2CorrvsS1Corr, *hTotalS2XYCorrvsS1Corr;
   TH2D *hTotalS1CorrvsEnergy, *hS1LYvsEnergy, *hTotalS2XYCorrvsEnergy, *hS2XYCorrLYvsEnergy, *hS2overS1vsEnergy;
   TH2D *ht_drift_S2overS1, *ht_drift_S2XYCorroverS1, *ht_drift_ScaledS2XYCorroverS1_Kr, *ht_drift_S2XYCorroverS1_Kr, *ht_drift_S2XYCorr_Kr;
   TH2D *ht_drift_S2AllCorr_Kr, *ht_drift_S2overS1_center;
   TH2D *ht_drift_TBAsym, *hTotalS1_TBAsym, *hDiffTotalS1_TBAsym, *hTotalS1_TBAcorr_TBAsym, *hVarianceXYVsR;
   TH2D *hRvsS2BotoverS2Top, *hVarianceXYvsS2BotoverS2Top;

   TH2D *hS2overS1vsS1_Center, *hScaledS2overS1vsS1, *hScaledS2XYCorroverS1vsS1;

   TProfile *pdrVsR;

   TH1D *hTotalS1Corr_Fiducial;


   ClassDef(SladReadSelector,0);
};

#endif

#ifdef SladReadSelector_cxx
void SladReadSelector::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

//   fChain->SetBranchStatus("*",0); //disable all

//   b_run_id->SetStatus(1);
   fChain->SetBranchAddress("run_id", &run_id, &b_run_id);
   fChain->SetBranchAddress("subrun_id", &subrun_id, &b_subrun_id);
   fChain->SetBranchAddress("event_id", &event_id, &b_event_id);
   fChain->SetBranchAddress("lifetime", &livetime, &b_livetime);
   fChain->SetBranchAddress("inhibittime", &inhibittime, &b_inhibittime);
   fChain->SetBranchAddress("acqui_window", &acqui_window, &b_acqui_window);
   fChain->SetBranchAddress("nchannel", &nchannel, &b_nchannel);
   fChain->SetBranchAddress("SumChannelHasNoBaseline", &baseline_not_found, &b_SumChannelHasNoBaseline);
   fChain->SetBranchAddress("n_phys_pulses", &npulses, &b_n_phys_pulses);
   fChain->SetBranchAddress("has_s3", &has_s3, &b_has_s3);

   fChain->SetBranchAddress("total_s1", &total_s1, &b_total_s1);
   fChain->SetBranchAddress("total_s1_corr", &total_s1_corr, &b_total_s1_corr);
   fChain->SetBranchAddress("total_s1_top", &total_s1_top, &b_total_s1_top);
   fChain->SetBranchAddress("total_s1_bottom", &total_s1_bottom, &b_total_s1_bottom);
   fChain->SetBranchAddress("total_f90", &total_f90, &b_total_f90); // total_f90_fixed should be same as total_f90;
   fChain->SetBranchAddress("is_saturated_pulse0", &is_saturated_pulse0, &b_is_saturated_pulse0);
   fChain->SetBranchAddress("total_s2", &total_s2, &b_total_s2);
   fChain->SetBranchAddress("total_s2_corr", &total_s2_corr, &b_total_s2_corr);
   fChain->SetBranchAddress("total_s2_top", &total_s2_top, &b_total_s2_top);
   fChain->SetBranchAddress("total_s2_bottom", &total_s2_bottom, &b_total_s2_bottom);
   fChain->SetBranchAddress("total_s2_f90", &total_s2_f90, &b_total_s2_f90);
   fChain->SetBranchAddress("is_saturated_pulse1", &is_saturated_pulse1, &b_is_saturated_pulse1);
   fChain->SetBranchAddress("s1_start_time", &s1_start_time, &b_s1_start_time);
   fChain->SetBranchAddress("s1_end_time", &s1_end_time, &b_s1_end_time);
   fChain->SetBranchAddress("s2_start_time", &s2_start_time, &b_s2_start_time);
   fChain->SetBranchAddress("s2_end_time", &s2_end_time, &b_s2_end_time);

   fChain->SetBranchAddress("tdrift", &t_drift, &b_tdrift);
   fChain->SetBranchAddress("s1_max_frac", &s1_max_frac);
   fChain->SetBranchAddress("selection_results", selection_results, &b_selection_results);

   fChain->SetBranchAddress("masas_x", &masas_x, &b_masas_x);
   fChain->SetBranchAddress("masas_y", &masas_y, &b_masas_y);
   fChain->SetBranchAddress("masas_chi2", &masas_chi2, &b_masas_chi2);

   fChain->SetBranchAddress("xyl_SCM", &xyl_SCM, &b_xyl_SCM);
   fChain->SetBranchAddress("xyl_best_x", &xyl_best_x, &b_xyl_best_x);
   fChain->SetBranchAddress("xyl_best_y", &xyl_best_y, &b_xyl_best_y);
   fChain->SetBranchAddress("xyl_best_chi2", &xyl_best_chi2, &b_xyl_best_chi2);

//   fChain->SetBranchAddress("roi_total_prompt_charge", &roi_total_prompt_charge);
//   fChain->SetBranchAddress("lsv_max_window_charge", &lsv_max_window_charge);
//   fChain->SetBranchAddress("lsv_late_window_charge", &lsv_late_window_charge);
//   fChain->SetBranchAddress("wt_charge", &wt_charge);

   fChain->SetBranchAddress("max_s1_frac_cut_threshold99", &max_s1_frac_cut_threshold99, &b_max_s1_frac_cut_threshold);
   fChain->SetBranchAddress("max_s1_frac_cut_exceeds99", &max_s1_frac_cut_exceeds99, &b_max_s1_frac_cut_exceeds);

   fChain->SetBranchAddress("s2_max_chan", &s2_max_chan, &b_s2_max_chan);
   fChain->SetBranchAddress("s2_max_frac", &s2_max_frac, &b_s2_max_frac);
   fChain->SetBranchAddress("s2_chan", s2_ch_frac, &b_s2_chan);

   fChain->SetBranchAddress("pulse_info_npulses", &pulse_info_npulses, &b_pulse_info_npulses);
   fChain->SetBranchAddress("pulse_info_pulse_id", pulse_info_pulse_id, &b_pulse_info_pulse_id);
   fChain->SetBranchAddress("pulse_info_start_time", pulse_info_start_time, &b_pulse_info_start_time);

//   fChain->SetBranchAddress("radius", &radius);
}

Bool_t SladReadSelector::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef SladReadSelector_cxx
