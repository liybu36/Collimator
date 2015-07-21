//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Sep 28 14:25:06 2014 by ROOT version 5.34/21
// from TTree Recon/Reconstructed Events
// found on file: outnvetoambe9L_v8_cylinder_clustered.root
//////////////////////////////////////////////////////////

#ifndef ReconSelector_h
#define ReconSelector_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH2F.h>
#include <TH1F.h>
#include <TH1D.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>
#include <vector>

// Fixed size dimensions of array or collections stored in the TTree if any.
using namespace std;

class ReconSelector : public TSelector {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t fNumberOfEvents;

   // Declaration of leaf types
   vector<double>  *parent_energy;
   vector<int>     *parent_pdg;
   Int_t           event_id;
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
   
   // List of branches
   TBranch        *b_parent_energy;   //!
   TBranch        *b_parent_pdg;   //!
   TBranch        *b_event_id;   //!
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

 ReconSelector(TTree * /*tree*/ =0) :
   fChain(0), fNumberOfEvents(0),
     Hist2DGap(0), QuenchingFactor(0), Hist2D(0), Hist2DXYslice1(0), Hist2DXYslice2(0), Hist2DXYslice3(0),
     Hist2Dnuclear(0), Hist2DXYnuclearslice1(0), Hist2DXYnuclearslice2(0), Hist2DXYnuclearslice3(0),
     Hist2Dsingle(0), Hist2DXYsingleslice1(0), Hist2DXYsingleslice2(0), Hist2DXYsingleslice3(0),
     Hist2Dsinglenuclear(0), Hist2DXYsinglenuclearslice1(0), Hist2DXYsinglenuclearslice2(0), Hist2DXYsinglenuclearslice3(0),
     Hist2DYZsingle(0), Hist2DYZsinglenuclear(0), Hist2DYZ(0), Hist2DYZnuclear(0), Hist2DXZ(0), Hist2DXZnuclear(0),
     Hist1DNRcountsslice1(0), Hist1DNRcountsslice2(0), Hist1DNRcountsslice3(0), Hist1DNRcountsslice4(0),
     Hist1DTPCedep(0), Hist1DTPCedep_nuclear(0), Hist1DTPCedep_electron(0), Hist1DNVedep(0), Hist1DNVedep_nuclear(0), Hist1DNVedep_electron(0)
     { }
   
   virtual ~ReconSelector() { }
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

 private:
   /* vector<double> nrx;
   vector<double> nry;
   vector<double> nrz;
   vector<double> rx;
   vector<double> ry;
   vector<double> rz;
   
   int eff;
   bool IsNR,IsER;
   int NNR, NER, ndeps, nsing, NRcounts;
   */
   //Create Histograms and Fill them
   void BookHistograms();
   void FillHistograms();

   // Initialize Histograms
   TH2F *Hist2DGap;
   TH1F *QuenchingFactor;
   TH2F *Hist2D, *Hist2DXYslice1, *Hist2DXYslice2, *Hist2DXYslice3;
   TH2F *Hist2Dnuclear, *Hist2DXYnuclearslice1, *Hist2DXYnuclearslice2, *Hist2DXYnuclearslice3;
   TH2F *Hist2Dsingle, *Hist2DXYsingleslice1, *Hist2DXYsingleslice2, *Hist2DXYsingleslice3;
   TH2F *Hist2Dsinglenuclear, *Hist2DXYsinglenuclearslice1, *Hist2DXYsinglenuclearslice2, *Hist2DXYsinglenuclearslice3;
   TH2F *Hist2DYZsingle, *Hist2DYZsinglenuclear, *Hist2DYZ, *Hist2DYZnuclear, *Hist2DXZ, *Hist2DXZnuclear;
   TH1F *Hist1DNRcountsslice1, *Hist1DNRcountsslice2, *Hist1DNRcountsslice3, *Hist1DNRcountsslice4;
   TH1F *Hist1DTPCedep, *Hist1DTPCedep_nuclear, *Hist1DTPCedep_electron, *Hist1DNVedep, *Hist1DNVedep_nuclear, *Hist1DNVedep_electron;
   ClassDef(ReconSelector,0);
};

#endif

#ifdef ReconSelector_cxx
void ReconSelector::Init(TTree *tree)
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
   parent_pdg = 0;
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
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("parent_energy", &parent_energy, &b_parent_energy);
   fChain->SetBranchAddress("parent_pdg", &parent_pdg, &b_parent_pdg);
   fChain->SetBranchAddress("event_id", &event_id, &b_event_id);
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
}

Bool_t ReconSelector::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // #ifdef ReconSelector_cxx
