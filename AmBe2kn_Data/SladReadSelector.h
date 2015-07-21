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
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include "TProfile3D.h"
#include "TCutG.h"

#include "../XYReconstructor/XYReconstructor.h"

// Header file for the classes stored in the TTree if any.

class SladReadSelector : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   XYReconstructor *fMyFcn;

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           run_id;
   Int_t           event_id;
   Double_t        livetime;
   Double_t        inhibittime;
   Float_t         total_s1;
   Float_t         total_s1_corr;
   Float_t         total_f90;
   Float_t         total_s2;
   Float_t         total_s2_corr;
   Float_t         total_s2_true;
   Double_t        t_drift;
   Float_t         s1_max_frac;
   Short_t         selection_results[16];

   Float_t         masas_x;
   Float_t         masas_y;
   Float_t         masas_chi2;

   Float_t         roi_total_prompt_charge;
   Float_t         lsv_max_window_charge;
   Float_t         lsv_late_window_charge;
   Float_t         wt_charge;

   Float_t         mf99;
   Float_t         radius;



   // List of branches
   TBranch        *b_run_id;   //!
   TBranch        *b_event_id;   //!
   TBranch        *b_livetime;   //!
   TBranch        *b_inhibittime;   //!
   TBranch        *b_total_s1;   //!
   TBranch        *b_total_s1_corr;   //!
   TBranch        *b_total_f90;   //!
   TBranch        *b_total_s2;   //!
   TBranch        *b_total_s2_corr;   //!
   TBranch        *b_tdrift;   //!
   TBranch        *b_selection_results;   //!

   TBranch        *b_masas_x;   //!
   TBranch        *b_masas_y;   //!
   TBranch        *b_masas_chi2;   //!


   SladReadSelector(TTree * /*tree*/ =0) :
     fChain(0), fMyFcn(0),
     fOutName(""),
     fOut_Tree(0),fCutg(0),fProfileName(""),
     pS2OverS1corr_S1corr_XY(0),
     pS2_S1corr_XY(0),
     hS2overS1corr_R_S1corr(0),
     hRun_TDrift(0), hTotalS1corr_TDrift(0), hTotalF90_TDrift(0), hRun_dT_Trigger(0),
     hRun_livetime(0), hRun_inhibittime(0),
     hTDrift(0), hTotalS1corr_TotalF90(0), hS1Total_Ext_vs_F90(0),
     hR_TDrift(0), hR_TDrift_Max(0), hR_TDrift_Min(0),
     hTotal_livetime(0), hR_TotalS1_TotalS2(0), hR_TotalS1_TotalS2_corr(0)
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

   Bool_t          IsValidEvent();
   Bool_t          IsEventOK_Default();
   Bool_t          IsEventOK();
   void            FillHistograms();
   void            BookHistograms();
   void            SetCutG4Outliers();

   TProfile3D *pS2OverS1corr_S1corr_XY, *pS2_S1corr_XY;
   TH3D *hS2overS1corr_R_S1corr;
   TH2D *hRun_TDrift, *hTotalS1corr_TDrift, *hTotalF90_TDrift, *hRun_dT_Trigger;
   TH2D *hRun_livetime, *hRun_inhibittime;
   TH1D *hTDrift;
   TH2D *hTotalS1corr_TotalF90, *hS1Total_Ext_vs_F90;
   TH2D *hR_TDrift, *hR_TDrift_Max, *hR_TDrift_Min;
   TH1F* hTotal_livetime;

   TH3D *hR_TotalS1_TotalS2, *hR_TotalS1_TotalS2_corr;

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

   fChain->SetBranchAddress("run_id", &run_id, &b_run_id);
   fChain->SetBranchAddress("event_id", &event_id, &b_event_id);
   fChain->SetBranchAddress("lifetime", &livetime, &b_livetime);
   fChain->SetBranchAddress("inhibittime", &inhibittime, &b_inhibittime);
   fChain->SetBranchAddress("total_s1", &total_s1, &b_total_s1);
   fChain->SetBranchAddress("total_s1_corr", &total_s1_corr, &b_total_s1_corr);
   fChain->SetBranchAddress("total_f90", &total_f90, &b_total_f90); // total_f90_fixed should be same as total_f90;
   fChain->SetBranchAddress("total_s2", &total_s2, &b_total_s2);
   fChain->SetBranchAddress("total_s2_corr", &total_s2_corr, &b_total_s2_corr);
   fChain->SetBranchAddress("total_s2_true", &total_s2_true);
   fChain->SetBranchAddress("tdrift", &t_drift, &b_tdrift);
   fChain->SetBranchAddress("s1_max_frac", &s1_max_frac);
   fChain->SetBranchAddress("selection_results", selection_results, &b_selection_results);

   fChain->SetBranchAddress("masas_x", &masas_x, &b_masas_x);
   fChain->SetBranchAddress("masas_y", &masas_y, &b_masas_y);
   fChain->SetBranchAddress("masas_chi2", &masas_chi2, &b_masas_chi2);

   fChain->SetBranchAddress("roi_total_prompt_charge", &roi_total_prompt_charge);
   fChain->SetBranchAddress("lsv_max_window_charge", &lsv_max_window_charge);
   fChain->SetBranchAddress("lsv_late_window_charge", &lsv_late_window_charge);
   fChain->SetBranchAddress("wt_charge", &wt_charge);

   fChain->SetBranchAddress("mf99", &mf99);
   fChain->SetBranchAddress("radius", &radius);
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
