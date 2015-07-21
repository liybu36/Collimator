/*

  2015-02-09 AFan

  Example macro to apply cuts and produce histograms using TTree::Draw().
  
  Usage:
  $ root -b -q example1.C(+)

  To run in compiled mode, use the "+".
  
 */

#include "TROOT.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TCut.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"

#include <iostream>


void make_histos(TChain* events)
{


  // Define cuts
  TCut cx1 = "nchannel.nchannel==38";
  TCut cx2 = "baseline.SumChannelHasNoBaseline==0";
  TCut cx3 = "(long_lifetime.lifetime+long_lifetime.inhibittime)>=1.35e-3";
  TCut cx4 = "long_lifetime.lifetime < 1";
  TCut basic_cuts = cx1 && cx2 && cx3 && cx4;
  TCut cx8 = "npulses.n_phys_pulses == 2 || (npulses.n_phys_pulses == 3 && npulses.has_s3)";
  TCut cx9 = "(events.run_id >= -999 && events.run_id < 7344   && s1_time.s1_start_time >= -0.25 && s1_time.s1_start_time <= -0.15) || \
              (events.run_id >= 7344 && events.run_id < 7641   && s1_time.s1_start_time >= -4.10 && s1_time.s1_start_time <= -4.00) || \
              (events.run_id >= 7641 && events.run_id < 999999 && s1_time.s1_start_time >= -6.10 && s1_time.s1_start_time <= -6.00)";


  // Open output file
  TFile* outfile = new TFile("example1.root", "recreate");
  std::cout << "Saving output to "<<outfile->GetName()<<std::endl;

  new TCanvas();
  
  // Fill histograms
  TH2F* h_s1_vs_s2 = new TH2F("h_s1_vs_s2", "; S1_{corr} [PE]; S2 [PE]", 1000, 0, 10000, 1000, 0, 1.e5);
  std::cout << "Drawing histogram "<<h_s1_vs_s2->GetName()<<"..."<<std::endl;
  events->Draw("total_s2:s1.total_s1_corr>>h_s1_vs_s2", basic_cuts && cx8 && cx9);

  
  outfile->Write();
  outfile->Close();
}

void example1()
{
    // Prevent canvases from being drawn.
  gROOT->SetBatch(kTRUE);

  TStopwatch* clock = new TStopwatch();
  clock->Start();
  
  // The main SLAD file containing the data we want.
  TString mainfile = "data/example.root";

  // Construct the name of the s2 file
  TString s2file = mainfile;
  s2file.Remove(s2file.Length()-5);
  s2file+="_s2.root";
  
  // Load the TTrees 
  TChain* events = new TChain("events");
  events->Add(mainfile);

  TChain* s2_fraction = new TChain("s2_fraction");
  s2_fraction->Add(s2file);
  events->AddFriend(s2_fraction);

  make_histos(events);

  std::cout << "Done!" << " "<<clock->RealTime()<<" s."<<std::endl;
}
