/*

  2015-02-09 AFan

  Example macro to apply cuts and produce histograms using an event loop.

  Usage:
  $ root -b -q example2.C(+)

  To run in compiled mode, use the "+".
  
 */

#include "TROOT.h"
#include "TStopwatch.h"
#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"

#include <iostream>


//------------------------------------------------------------------------------
// Simple struct to hold the variables that we extract from SLAD files.
// To add a variable:
// 1. Add it to this struct
// 2. Set the appropriate branch address in load_tree()
struct MyEvent
{
  // variables from main file
  int    run_id;
  int    event_id;
  int    nchannels;
  short  baseline_not_found;
  double live_time;
  double inhibit_time;
  int    npulses;
  int    has_s3;
  float  s1_start_time;
  float  s1;
  float  s2;

  // variables from s2 file
  float  s2_ch_frac[38];

  // variables from xy file
  float  x;
  float  y;

  // variables from allpulses file
  float pulse_start_time[100];
  float pulse_total_npe[100];

};


//------------------------------------------------------------------------------
// Set the branch address for the variables in MyEvent.
// To increase speed, we disable all branches by default, and turn on only the
// ones we want to use.
void load_tree(TChain* events, MyEvent & e)
{
  events->SetBranchStatus("*",0); //disable all
  
  events->SetBranchStatus("events.run_id", 1);
  events->SetBranchAddress("run_id", &e.run_id);

  events->SetBranchStatus("events.event_id", 1);
  events->SetBranchAddress("event_id", &e.event_id);
  
  events->SetBranchStatus("nchannel.nchannel", 1);
  events->SetBranchAddress("nchannel.nchannel", &e.nchannels);

  events->SetBranchStatus("baseline.SumChannelHasNoBaseline",1);
  events->SetBranchAddress("baseline.SumChannelHasNoBaseline", &e.baseline_not_found);

  events->SetBranchStatus("long_lifetime.lifetime",1);
  events->SetBranchAddress("long_lifetime.lifetime", &e.live_time);

  events->SetBranchStatus("long_lifetime.inhibittime",1);
  events->SetBranchAddress("long_lifetime.inhibittime", &e.inhibit_time);

  events->SetBranchStatus("npulses.n_phys_pulses",1);
  events->SetBranchAddress("npulses.n_phys_pulses", &e.npulses);

  events->SetBranchStatus("npulses.has_s3",1);
  events->SetBranchAddress("npulses.has_s3", &e.has_s3);

  events->SetBranchStatus("s1_time.s1_start_time", 1);
  events->SetBranchAddress("s1_time.s1_start_time", &e.s1_start_time);

  events->SetBranchStatus("s1.total_s1_corr", 1);
  events->SetBranchAddress("s1.total_s1_corr", &e.s1);

  events->SetBranchStatus("s2.total_s2_corr", 1);
  events->SetBranchAddress("s2.total_s2_corr", &e.s2);

  events->SetBranchStatus("s2_fraction.s2_chan", 1);
  events->SetBranchAddress("s2_fraction.s2_chan", e.s2_ch_frac);

  events->SetBranchStatus("masas_xy.masas_x", 1);
  events->SetBranchAddress("masas_xy.masas_x", &e.x);

  events->SetBranchStatus("masas_xy.masas_y", 1);
  events->SetBranchAddress("masas_xy.masas_y", &e.y);

  events->SetBranchStatus("pulse_info.pulse_info_start_time", 1);
  events->SetBranchAddress("pulse_info.pulse_info_start_time", e.pulse_start_time);

  events->SetBranchStatus("pulse_info.pulse_info_total_npe", 1);
  events->SetBranchAddress("pulse_info.pulse_info_total_npe", e.pulse_total_npe);
}


//------------------------------------------------------------------------------
// Define cuts. These definitions are taken from the 2014 analysis memo v11.
bool CX1(MyEvent const& e) { return e.nchannels==38 ; }
bool CX2(MyEvent const& e) { return e.baseline_not_found == false; }
bool CX3(MyEvent const& e) { return (e.live_time+e.inhibit_time)>=1.35e-3; }
bool CX4(MyEvent const& e) { return e.live_time < 1.; }
bool CX8(MyEvent const& e) { return e.npulses==2 || (e.npulses==3 && e.has_s3); }
bool CX9(MyEvent const& e) { return
    ((e.run_id >= -999 && e.run_id < 7344   && e.s1_start_time >= -0.25 && e.s1_start_time <= -0.15) ||
     (e.run_id >= 7344 && e.run_id < 7641   && e.s1_start_time >= -4.10 && e.s1_start_time <= -4.00) ||
     (e.run_id >= 7641 && e.run_id < 999999 && e.s1_start_time >= -6.10 && e.s1_start_time <= -6.00)); }


//------------------------------------------------------------------------------
// Main event loop is contained here.
void event_loop(TChain* events)
{
 
  // Set the branch addresses.
  MyEvent e;
  load_tree(events, e);

  // Open output file
  TFile* outfile = new TFile("example2.root", "recreate");
  std::cout << "Writing output to "<<outfile->GetName()<<std::endl;

  // Define histograms
  outfile->cd();
  TH2F* h_s1_vs_s2  = new TH2F("h_s1_vs_s2", "; S1_{corr} [PE]; S2 [PE]", 1000, 0, 10000, 1000, 0, 1.e5);
  TH2F* h_s2_ch_occ = new TH2F("h_s2_ch_occ", "; ch ID; S2_{i}/S2_{tot}", 40, 0, 40, 1000, 0, 1);
  TH2F* h_xy        = new TH2F("h_xy", "; x [cm]; y [cm]", 80, -20, 20, 80, -20, 20);
  TH1F* h_p0_start  = new TH1F("h_p0_start", "; pulse start time [#mus]", 1000, -10, 10);

  //-------------------------//
  //     MAIN EVENT LOOP     //
  //-------------------------//
  int nevents = events->GetEntries();
  for (int n=0; n<nevents; ++n) {
    if (n%100000==0) std::cout << "Processing event "<<n<<"/"<<nevents<<std::endl;
    
    events->GetEntry(n);

    // Generate cuts.
    bool cx1 = CX1(e);
    bool cx2 = CX2(e);
    bool cx3 = CX3(e);
    bool cx4 = CX4(e);
    bool cx8 = CX8(e);
    bool cx9 = CX9(e);
    bool basic_cuts = cx1 && cx2 && cx3 && cx4;

    // Apply cuts and fill histograms
    if (basic_cuts && cx8 && cx9) h_s1_vs_s2->Fill(e.s1, e.s2);
    if (basic_cuts && cx8 && cx9) h_xy->Fill(e.x, e.y);
    if (basic_cuts && cx8 && cx9) { for (int ch=0; ch<38; ++ch) h_s2_ch_occ->Fill(ch, e.s2_ch_frac[ch]); }
    if (basic_cuts && e.npulses>0) h_p0_start->Fill(e.pulse_start_time[0]);
    
  }//loop over events

  outfile->Write();
  outfile->Close();

}


//------------------------------------------------------------------------------
// Main method. Load DST files and invoke event_loop().
void example2()
{
   // Prevent canvases from being drawn.
  gROOT->SetBatch(kTRUE);

  TStopwatch* clock = new TStopwatch();
  clock->Start();

  // The main SLAD file containing the data we want.
  TString mainfile = "data/example.root";
  

  // Construct the name of the friend tree files
  TString s2file = mainfile;
  s2file.Remove(s2file.Length()-5);
  s2file+="_s2.root";

  TString xyfile = mainfile;
  xyfile.Remove(xyfile.Length()-5);
  xyfile+="_masas_xy.root";

  TString pulsefile = mainfile;
  pulsefile.Remove(pulsefile.Length()-5);
  pulsefile+="_allpulses.root";
  
  // Load and friend the TTrees 
  TChain* events = new TChain("events");
  events->Add(mainfile);

  TChain* s2_fraction = new TChain("s2_fraction");
  s2_fraction->Add(s2file);
  events->AddFriend(s2_fraction);

  TChain* xy = new TChain("masas_xy");
  xy->Add(xyfile);
  events->AddFriend(xy);

  TChain* pulse_info = new TChain("pulse_info");
  pulse_info->Add(pulsefile);
  events->AddFriend(pulse_info);

  event_loop(events);
  
  std::cout << "Done!"<<" "<<clock->RealTime()<<" s."<<std::endl;
}
