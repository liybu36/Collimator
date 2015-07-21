//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Jun  5 17:31:00 2015 by ROOT version 5.34/21
// from TTree Recon/Reconstructed Events
// found on file: outneutron_clustered.root
//////////////////////////////////////////////////////////

#ifndef reconneutronSelector_h
#define reconneutronSelector_h

#include <algorithm>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TNtuple.h>
#include <TString.h>
#include <TRint.h>
// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.
using namespace std;

class reconneutronSelector : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain

   // Declaration of leaf types
   vector<double>  *parent_energy;
   Int_t           event_id;
   Int_t           event_pdg;
   Float_t         source_ene;
   Float_t         source_x;
   Float_t         source_y;
   Float_t         source_z;
   Float_t         source_r;
   Float_t         source_px;
   Float_t         source_py;
   Float_t         source_pz;
   vector<TString> *volume;
   vector<TString> *particle;
   Int_t           n_active;
   vector<double>  *et;
   vector<double>  *ex;
   vector<double>  *ey;
   vector<double>  *ez;
   vector<double>  *dt;
   vector<double>  *dz;
   vector<double>  *dr;
   vector<double>  *edep;
   vector<double>  *edep_nuclear;
   vector<double>  *edep_electron;
   vector<double>  *eqch;
   vector<double>  *quenchingfactor;
   vector<bool>    *contain_alpha;

   // List of branches
   TBranch        *b_parent_energy;   //!
   TBranch        *b_event_id;   //!
   TBranch        *b_event_pdg;   //!
   TBranch        *b_source_ene;   //!
   TBranch        *b_source_x;   //!
   TBranch        *b_source_y;   //!
   TBranch        *b_source_z;   //!
   TBranch        *b_source_r;   //!
   TBranch        *b_source_px;   //!
   TBranch        *b_source_py;   //!
   TBranch        *b_source_pz;   //!
   TBranch        *b_volume;   //!
   TBranch        *b_particle;   //!
   TBranch        *b_n_active;   //!
   TBranch        *b_et;   //!
   TBranch        *b_ex;   //!
   TBranch        *b_ey;   //!
   TBranch        *b_ez;   //!
   TBranch        *b_dt;   //!
   TBranch        *b_dz;   //!
   TBranch        *b_dr;   //!
   TBranch        *b_edep;   //!
   TBranch        *b_edep_nuclear;   //!
   TBranch        *b_edep_electron;   //!
   TBranch        *b_eqch;   //!
   TBranch        *b_quenchingfactor;   //!
   TBranch        *b_contain_alpha;   //!

 reconneutronSelector(TTree * /*tree*/ =0) : fChain(0),
     source_ntuple(0),source_ene_theta(0),tpc_events(8,0),veto_events(5,0)
     { }
   virtual ~reconneutronSelector() { }
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
   void SetTRint(TRint *rint){fRint = rint;}

 private:
   TRint* fRint;
   void BookHistograms();
   void FillHistograms();

   TNtuple* source_ntuple;
   TH2F* source_ene_theta;
   vector<TString> ntuple_name;
   vector<TNtuple*> ntuple_plots;
   
   vector<TString> eqch_name;
   vector<TH1F*> eqch_hist;
   vector<TString> eqch_time_name;
   vector<TH2F*> eqch_time_hist;
   vector<TString> eqch_sourceene_name;
   vector<TH2F*> eqch_sourceene_hist;
   vector<TString> eqch_sourcetheta_name;
   vector<TH2F*> eqch_sourcetheta_hist;
   vector<TString> eqch_theta_name;
   vector<TH2F*> eqch_theta_hist;
   vector<TString> eqch_radius_name;
   vector<TH2F*> eqch_radius_hist;
   
   vector<TString> tpc_total_name;
   vector<TH1F*> tpc_total;
   vector<TString> tpc_total_nv_name;
   vector<TH2F*> tpc_total_nv;

   vector<int> tpc_events;
   vector<int> veto_events;

   int SetTail();
   vector<string> tail;
   vector<double> timelow;
   vector<double> timehigh;

   ClassDef(reconneutronSelector,0);
};

#endif

#ifdef reconneutronSelector_cxx
void reconneutronSelector::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   parent_energy = 0;
   volume = 0;
   particle = 0;
   et = 0;
   ex = 0;
   ey = 0;
   ez = 0;
   dt = 0;
   dz = 0;
   dr = 0;
   edep = 0;
   edep_nuclear = 0;
   edep_electron = 0;
   eqch = 0;
   quenchingfactor = 0;
   contain_alpha = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("parent_energy", &parent_energy, &b_parent_energy);
   fChain->SetBranchAddress("event_id", &event_id, &b_event_id);
   fChain->SetBranchAddress("event_pdg", &event_pdg, &b_event_pdg);
   fChain->SetBranchAddress("source_ene", &source_ene, &b_source_ene);
   fChain->SetBranchAddress("source_x", &source_x, &b_source_x);
   fChain->SetBranchAddress("source_y", &source_y, &b_source_y);
   fChain->SetBranchAddress("source_z", &source_z, &b_source_z);
   fChain->SetBranchAddress("source_r", &source_r, &b_source_r);
   fChain->SetBranchAddress("source_px", &source_px, &b_source_px);
   fChain->SetBranchAddress("source_py", &source_py, &b_source_py);
   fChain->SetBranchAddress("source_pz", &source_pz, &b_source_pz);
   fChain->SetBranchAddress("volume", &volume, &b_volume);
   fChain->SetBranchAddress("particle", &particle, &b_particle);
   fChain->SetBranchAddress("n_active", &n_active, &b_n_active);
   fChain->SetBranchAddress("et", &et, &b_et);
   fChain->SetBranchAddress("ex", &ex, &b_ex);
   fChain->SetBranchAddress("ey", &ey, &b_ey);
   fChain->SetBranchAddress("ez", &ez, &b_ez);
   fChain->SetBranchAddress("dt", &dt, &b_dt);
   fChain->SetBranchAddress("dz", &dz, &b_dz);
   fChain->SetBranchAddress("dr", &dr, &b_dr);
   fChain->SetBranchAddress("edep", &edep, &b_edep);
   fChain->SetBranchAddress("edep_nuclear", &edep_nuclear, &b_edep_nuclear);
   fChain->SetBranchAddress("edep_electron", &edep_electron, &b_edep_electron);
   fChain->SetBranchAddress("eqch", &eqch, &b_eqch);
   fChain->SetBranchAddress("quenchingfactor", &quenchingfactor, &b_quenchingfactor);
   fChain->SetBranchAddress("contain_alpha", &contain_alpha, &b_contain_alpha);
}

Bool_t reconneutronSelector::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef reconneutronSelector_cxx



/*
   vector<TString> tpc_gamma1d_name;
   vector<TH1F*> tpc_gamma1d;
   vector<TString> tpc_gamma2d_name;
   vector<TH2F*> tpc_gamma2d;

   vector<TString> tpc_neutron1d_name; 
   vector<TH1F*> tpc_neutron1d;
   vector<TString> tpc_neutron2d_name; 
   vector<TH2F*> tpc_neutron2d;
   
   vector<TString> nv_gamma1d_name;
   vector<TH1F*> nv_gamma1d;
   vector<TString> nv_gamma2d_name;
   vector<TH2F*> nv_gamma2d;

   vector<TString> nv_neutron1d_name;
   vector<TH1F*> nv_neutron1d;
   vector<TString> nv_neutron2d_name;
   vector<TH2F*> nv_neutron2d;
*/
