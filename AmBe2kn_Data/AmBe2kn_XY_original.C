#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <math.h>
#include "TString.h"
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TRint.h"
#include "TMath.h"
#include "TH1D.h"
#include "TCanvas.h"

TRint *theApp;
using namespace std;

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
  double tdrift;
  float  s1;
  float  s2;
  float  s1_total_f90;
  float  s2_total_f90;

  // variables from s2 file
  float  s2_ch_frac[38];

  // variables from xy file
  float  masa_x;
  float  masa_y;
  float  jason_x;
  float  jason_y;

};

bool load_tree(TChain* events, MyEvent & e)
{
  events->SetBranchStatus("*",0); //disable all
  events->SetBranchStatus("events.run_id", 1);
  events->SetBranchStatus("events.event_id", 1);
  events->SetBranchStatus("nchannel.nchannel", 1);
  events->SetBranchStatus("baseline.SumChannelHasNoBaseline",1);
  events->SetBranchStatus("long_lifetime.lifetime",1);
  events->SetBranchStatus("long_lifetime.inhibittime",1);
  events->SetBranchStatus("npulses.n_phys_pulses",1);
  events->SetBranchStatus("npulses.has_s3",1);
  events->SetBranchStatus("tdrift.tdrift",1);
  events->SetBranchStatus("s1_time.s1_start_time", 1);
  events->SetBranchStatus("s1.total_s1_corr", 1);
  events->SetBranchStatus("s1_f90.total_f90", 1);
  events->SetBranchStatus("s2.total_s2_corr", 1);
  events->SetBranchStatus("s2_fraction.s2_chan", 1);
  events->SetBranchStatus("s2_f90.total_s2_f90", 1);
  events->SetBranchStatus("masas_xy.masas_x", 1);
  events->SetBranchStatus("masas_xy.masas_y", 1);
  events->SetBranchStatus("xylocator_xy.xyl_best_x", 1);
  events->SetBranchStatus("xylocator_xy.xyl_best_y", 1);
 
  events->SetBranchAddress("run_id", &e.run_id);
  events->SetBranchAddress("event_id", &e.event_id);
  events->SetBranchAddress("nchannel.nchannel", &e.nchannels);
  events->SetBranchAddress("baseline.SumChannelHasNoBaseline", &e.baseline_not_found);
  events->SetBranchAddress("long_lifetime.lifetime", &e.live_time);
  events->SetBranchAddress("long_lifetime.inhibittime", &e.inhibit_time);
  events->SetBranchAddress("npulses.n_phys_pulses", &e.npulses);
  events->SetBranchAddress("npulses.has_s3", &e.has_s3);
  events->SetBranchAddress("tdrift.tdrift", &e.tdrift);
  events->SetBranchAddress("s1_time.s1_start_time", &e.s1_start_time);
  events->SetBranchAddress("s1.total_s1_corr", &e.s1);
  events->SetBranchAddress("s1_f90.total_f90", &e.s1_total_f90);
  events->SetBranchAddress("s2.total_s2_corr", &e.s2);
  events->SetBranchAddress("s2_fraction.s2_chan", e.s2_ch_frac);
  events->SetBranchAddress("s2_f90.total_s2_f90", &e.s2_total_f90);
  events->SetBranchAddress("masas_xy.masas_x", &e.masa_x);
  events->SetBranchAddress("masas_xy.masas_y", &e.masa_y);  
  events->SetBranchAddress("xylocator_xy.xyl_best_x", &e.jason_x);
  events->SetBranchAddress("xylocator_xy.xyl_best_y", &e.jason_y);

  return true;
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

//---------------------------------------------------------------------------------

// Main event loop is contained here.
void event_loop(TChain* events, TString outfilename)
{

  // Set the branch addresses.
  MyEvent e;
  if(!load_tree(events, e))
    cout<<"Cannot load tree!"<<endl;

  // Define histograms
  TH2F* h_s1_vs_s2  = new TH2F("h_s1_vs_s2", "; S1_{corr} [PE]; S2 [PE]", 1000, 0, 10000, 1000, 0, 1.e5);
  TH2F* h_s2_ch_occ = new TH2F("h_s2_ch_occ", "; ch ID; S2_{i}/S2_{tot}", 40, 0, 40, 1000, 0, 1);
  TH2F* h_xy        = new TH2F("h_xy", "; x [cm]; y [cm]", 80, -20, 20, 80, -20, 20);

  TH2F* s1_f90 = new TH2F("s1_f90","; total_s1_corr [PE]; s1_total_f90",1000,0,10000,100,0,1);
  TH2F* xy   = new TH2F("xy", "; x [cm]; y [cm]", 80, -20, 20, 80, -20, 20);
  TH2F* r_tdrift   = new TH2F("r_tdrift", "; tdrift [us]; r [cm]", 100, 0, 400, 80, 0, 20);
  TH2F* theta_tdrift   = new TH2F("theta_tdrift", "; tdrift [us]; theta [rad]", 100, 0, 400, 100, -3, 3);
  
  string syntax[] = {"neutron","single","singleNR"};

  vector<TH2F*> xy_slice;
  vector<TH2F*> r_tdrift_slice;
  vector<TH2F*> theta_tdrift_slice;
  vector<TH1D*> xy_py;
  vector<TH1D*> xy_px;

  for(int i=0; i<(int)(sizeof(syntax)/sizeof(string)); i++)
    {
      xy_slice.push_back(new TH2F(Form("%s_xy_slice",syntax[i].c_str()),"; x [cm]; y [cm]",80, -20, 20, 80, -20, 20));
      r_tdrift_slice.push_back(new TH2F(Form("%s_r_tdrift_slice",syntax[i].c_str()),"; tdrift [us]; r [cm]", 100, 0, 400, 80, 0, 20)); 
      theta_tdrift_slice.push_back(new TH2F(Form("%s_theta_tdrift_slice",syntax[i].c_str()),"; tdrift [us]; theta [rad]", 100, 0, 400, 100, -3, 3)); 
    }
  vector<int> counts (2,0);
  vector<float> probs (2,0);

  string tag[3] = {"all","NR","ER"};
  int start_s1 = 0;
  int end_s1 = 200;
  int step = 10;
  int NUM = (end_s1-start_s1)/step;
  vector<double> low_s1, high_s1;
  double f90_cut = 0.55;

  vector<TH2F*> f90_s1;
  vector<TH2F*> s2overs1_vs_s1;
  vector<TH2F*> s2overs1_vs_f90;
  vector<TH2F*> s2overs1_vs_f90_slice;
  vector<TH1D*> s2overs1_vs_f90_slice_py;
  vector<TH1D*> s2overs1_vs_f90_slice_px; 
  for(int i=0; i<3; i++)
    {
      f90_s1.push_back(new TH2F(Form("%s_f90_s1",tag[i].c_str()),"; total_s1_corr [PE]; s1_total_f90",1000,0,10000,100,0,1));
      s2overs1_vs_s1.push_back(new TH2F(Form("%s_s2overs1_vs_s1",tag[i].c_str()),"; total_s1_corr [PE]; log(s2/s1)",1000,0,10000,200,0,3));
      s2overs1_vs_f90.push_back(new TH2F(Form("%s_s2overs1_vs_f90",tag[i].c_str()),"; s1_total_f90 [PE]; log(s2/s1)",100,0,1,200,0,3));      
      //s2overs1_vs_f90_slice.push_back(new TH2F(Form("%s_s2overs1_vs_f90_slice",tag[i].c_str()),"; s1_total_f90 [PE]; log(s2/s1)",100,0,1,200,0,3));   
      for(int k=0; k<NUM; k++)
	{
	  low_s1.push_back(start_s1 + k*step +0.);
	  high_s1.push_back(start_s1 + (k+1)*step +0.);     
	  s2overs1_vs_f90_slice.push_back(new TH2F(Form("%s_s2overs1_vs_f90_slice%d",tag[i].c_str(),k),"; s1_total_f90 [PE]; log(s2/s1)",100,0,1,200,0,3));
	  // s2overs1_vs_f90_ER_slice.push_back(new TH2F(Form("s2overs1_vs_f90_ER_slice%d",i),"; s1_total_f90 [PE]; log(s2/s1)",100,0,f90_cut,200,0,3));
	  //s2overs1_vs_f90_NR_slice.push_back(new TH2F(Form("s2overs1_vs_f90_NR_slice%d",i),"; s1_total_f90 [PE]; log(s2/s1)",100,f90_cut,1,200,0,3)); 
	}
    }

#define CorrXY
#ifdef CorrXY
  TString fname = "S2CorrectionFactor_Kr_pS2VsXY_Kr_Jason_Masa.root";
#else
  TString fname = "S2CorrectionFactor_Kr_pR2_theta_Kr.root";
#endif
  TFile *g = (TFile*)gROOT->GetListOfFiles()->FindObject(fname.Data());
  g = new TFile(fname.Data());
  TString hname = "hS2KrVsXY_norm_center";
  TH2D *fhS2Correction_factor = (TH2D*) g->Get(hname.Data());  

  //-------------------------//
  //     MAIN EVENT LOOP     //
  //-------------------------//
  int nevents = events->GetEntries();
  cout<<"Events: "<<nevents<<endl;
  for (int n=0; n<nevents; ++n) {
    if (n%500000==0) std::cout << "Processing event "<<n<<"/"<<nevents<<std::endl;
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
    if (basic_cuts && cx8 && cx9)
      { 
	for(int ch=0; ch<38; ++ch)
	  h_s2_ch_occ->Fill(ch, e.s2_ch_frac[ch]);
	h_s1_vs_s2->Fill(e.s1, e.s2);
	h_xy->Fill(e.masa_x, e.masa_y);
	s1_f90->Fill(e.s1,e.s1_total_f90);
	float r = TMath::Sqrt(TMath::Power(e.jason_x,2)+TMath::Power(e.jason_y,2));
	float theta = TMath::ATan2(e.jason_y,e.jason_x);
	r_tdrift->Fill(e.tdrift,r);
	theta_tdrift->Fill(e.tdrift,theta);
	xy->Fill(e.jason_x,e.jason_y);
	++counts.at(0);

	f90_s1.at(0)->Fill(e.s1,e.s1_total_f90);
	double s2_corr_frac;
#ifdef CorrXY
	s2_corr_frac = 1./fhS2Correction_factor->GetBinContent(fhS2Correction_factor->FindBin(e.jason_x, e.jason_y));
#else
	s2_corr_frac = 1./fhS2Correction_factor->GetBinContent(fhS2Correction_factor->FindBin(e.jason_x*e.jason_x+e.jason_y*e.jason_y, TMath::ATan2(e.jason_y, e.jason_x)));
#endif
	double s2_xy_corr = e.s2*s2_corr_frac;
	s2overs1_vs_s1.at(0)->Fill(e.s1,TMath::Log10(s2_xy_corr*1.0/e.s1));
	s2overs1_vs_f90.at(0)->Fill(e.s1_total_f90,TMath::Log10(s2_xy_corr*1.0/e.s1));
	//	if(e.s1>0 && e.s1<0.)
	//	  s2overs1_vs_f90_slice.at(0)->Fill(e.s1_total_f90,TMath::Log10(s2_xy_corr*1.0/e.s1));

	for(int i=0; i<NUM; i++)
	  {	    
	    if(e.s1>low_s1.at(i) && e.s1<high_s1.at(i))
	      {
		s2overs1_vs_f90_slice.at(i)->Fill(e.s1_total_f90,TMath::Log10(s2_xy_corr*1.0/e.s1));
		if(e.s1_total_f90>0 && e.s1_total_f90<f90_cut)
		  s2overs1_vs_f90_slice.at(i+NUM*2)->Fill(e.s1_total_f90,TMath::Log10(s2_xy_corr*1.0/e.s1));
		if(e.s1_total_f90>f90_cut && e.s1_total_f90<1.)
		  s2overs1_vs_f90_slice.at(i+NUM*1)->Fill(e.s1_total_f90,TMath::Log10(s2_xy_corr*1.0/e.s1));		
	      }
	  }

	//Select s1 f90 to get the XY for collimator data
	if(e.s1_total_f90>0.6 && e.s1_total_f90<0.9 && e.s1>100.)
	  {
	    xy_slice.at(0)->Fill(e.jason_x,e.jason_y);
	    theta_tdrift_slice.at(0)->Fill(e.tdrift,theta);
	    ++counts.at(1);	   
	    if(theta>1. && theta<3.)
	      {
		r_tdrift_slice.at(0)->Fill(e.tdrift,r);
		xy_slice.at(1)->Fill(e.jason_x,e.jason_y);
		if(e.tdrift>60. && e.tdrift<350)	   
		  xy_slice.at(2)->Fill(e.jason_x,e.jason_y);
	      }
	  }
	
      }
  }//loop over events

  for(size_t j=0; j<xy_slice.size(); j++)
    {
      int binx11 = xy_slice.at(j)->GetXaxis()->FindBin(-10.);
      int binx12 = xy_slice.at(j)->GetXaxis()->FindBin(0.);
      int biny11 = xy_slice.at(j)->GetYaxis()->FindBin(10.);
      int biny12 = xy_slice.at(j)->GetYaxis()->FindBin(15.);
      xy_px.push_back(xy_slice.at(j)->ProjectionX(Form("%s_xy_px",syntax[j].c_str()),biny11,biny12));
      xy_py.push_back(xy_slice.at(j)->ProjectionY(Form("%s_xy_py",syntax[j].c_str()),binx11,binx12));
    }

  //  for(size_t j=0; j<s2overs1_vs_f90_slice.size(); j++)
  for(int i=0; i<3; i++)
    {
      for(int k=0; k<NUM; k++)
	{
	  int binx11 = s2overs1_vs_f90_slice.at(k+i*NUM)->GetXaxis()->FindBin(0.);
	  int binx12 = s2overs1_vs_f90_slice.at(k+i*NUM)->GetXaxis()->FindBin(1.);
	  int biny11 = s2overs1_vs_f90_slice.at(k+i*NUM)->GetYaxis()->FindBin(0.);
	  int biny12 = s2overs1_vs_f90_slice.at(k+i*NUM)->GetYaxis()->FindBin(3.);
	  s2overs1_vs_f90_slice_px.push_back(s2overs1_vs_f90_slice.at(k+i*NUM)->ProjectionX(Form("%s_s2overs1_vs_f90_slice%d_px",tag[i].c_str(),k),biny11,biny12));
	  s2overs1_vs_f90_slice_py.push_back(s2overs1_vs_f90_slice.at(k+i*NUM)->ProjectionY(Form("%s_s2overs1_vs_f90_slice%d_py",tag[i].c_str(),k),binx11,binx12));      
	}
    }
  
  vector<TCanvas*> canv;
  for(int k=0; k<NUM; k++)
    {
      s2overs1_vs_f90_slice_py.at(k+NUM)->SetLineColor(2);
      s2overs1_vs_f90_slice_py.at(k+NUM*2)->SetLineColor(4);
      TH1D* nr = (TH1D*)s2overs1_vs_f90_slice_py.at(k+NUM)->Clone("nr");
      TH1D* er = (TH1D*)s2overs1_vs_f90_slice_py.at(k+NUM*2)->Clone("er");

      canv.push_back(new TCanvas(Form("canv%d",k),"",600,400));
      canv.back()->cd();
      s2overs1_vs_f90_slice_py.at(k+NUM*2)->Draw();
      s2overs1_vs_f90_slice_py.at(k+NUM*1)->Draw("same");

      canv.push_back(new TCanvas(Form("cdf%d",k),"",600,400));
      canv.back()->cd();     
      /*
      s2overs1_vs_f90_slice_py.at(k+NUM*2)->ComputeIntegral();
      Double_t *integral2 = s2overs1_vs_f90_slice_py.at(k+NUM*2)->GetIntegral();
      s2overs1_vs_f90_slice_py.at(k+NUM*2)->SetContent(integral2);            
      s2overs1_vs_f90_slice_py.at(k+NUM*2)->Draw();

      s2overs1_vs_f90_slice_py.at(k+NUM*1)->ComputeIntegral();
      Double_t *integral1 = s2overs1_vs_f90_slice_py.at(k+NUM*1)->GetIntegral();
      s2overs1_vs_f90_slice_py.at(k+NUM*1)->SetContent(integral1);            
      s2overs1_vs_f90_slice_py.at(k+NUM*1)->Draw("same");      
      */
      nr->ComputeIntegral();
      Double_t *integral1 = nr->GetIntegral();
      nr->SetContent(integral1);
      nr->Draw();
      
      er->ComputeIntegral();
      Double_t *integral2 = er->GetIntegral();
      er->SetContent(integral2);
      er->Draw("same");
      }

  // Open output file
  TFile* outfile = new TFile(outfilename, "recreate");
  std::cout << "====> output file: "<<outfile->GetName()<<std::endl;
  h_s1_vs_s2->Write();
  h_s2_ch_occ->Write();
  h_xy->Write();
  s1_f90->Write();
  xy->Write();
  r_tdrift->Write();
  theta_tdrift->Write();
  for(size_t j=0; j<xy_slice.size(); j++)
    {
      xy_slice.at(j)->Write();
      r_tdrift_slice.at(j)->Write();
      theta_tdrift_slice.at(j)->Write();
      xy_py.at(j)->Write();
      xy_px.at(j)->Write();
    }
  for(size_t j=0; j<s2overs1_vs_f90.size(); j++)
    {
      f90_s1.at(j)->Write();
      s2overs1_vs_s1.at(j)->Write();
      s2overs1_vs_f90.at(j)->Write();
    }
  for(size_t j=0; j<s2overs1_vs_f90_slice.size(); j++)
    {
      s2overs1_vs_f90_slice.at(j)->Write();
      s2overs1_vs_f90_slice_px.at(j)->Write();
      s2overs1_vs_f90_slice_py.at(j)->Write();
    }
  for(size_t j=0; j<canv.size(); j++)
    canv.at(j)->Write();
  
  outfile->Write();
  outfile->Close();

  for(size_t i=0; i<counts.size(); i++)
    {
      probs.at(i) = 1.0*counts.at(i)/nevents;
      cout<<i<<"  "<<counts.at(i)<<"   "<<probs.at(i)<<endl;	
    }

  TString outdataname = outfilename;
  outdataname.Remove(outdataname.Length()-5);
  outdataname+=".txt";
  ofstream outdata (outdataname.Data());
  if(outdata.is_open())
    {
      for(size_t i=0; i<counts.size(); i++)
	outdata<<i<<"  "<<counts.at(i)<<"   "<<probs.at(i)<<"\n";
    }
  outdata.close();    
}

void AmBe2kn_XY(int start, int end)  
{
  TString Time ="May14";
  TString indir =  "/darkside/users/hqian/AmBe2kn_XYData/data/";

  TChain* events = new TChain("events");
  TChain* s2_fraction = new TChain("s2_fraction");
  TChain* masaxy = new TChain("masas_xy");
  TChain* jasonxy = new TChain("xylocator_xy");
   
  for(int i=start; i<=end; i++)
    { 
      TString mainfile;
      mainfile.Form("Run%06d.root",i);
      mainfile.Prepend(indir);

      ifstream NameCheck;
      NameCheck.open(mainfile.Data());
      if(!NameCheck.good())
	continue;
      else{
	TFile *f = new TFile(mainfile);
	if(f->IsZombie())
	  continue;
	else{	  
	  cout<<"Processing the data file: "<<mainfile<<endl;
	  TString s2file = mainfile;
	  s2file.Remove(s2file.Length()-5);
	  s2file+="_s2.root";
	  cout<<"Processing the data file: "<<s2file<<endl;
	  
	  TString masaxyfile = mainfile;
	  masaxyfile.Remove(masaxyfile.Length()-5);
	  masaxyfile+="_masas_xy.root";
	  cout<<"Processing the data file: "<<masaxyfile<<endl;

	  TString jasonxyfile = mainfile;
	  jasonxyfile.Remove(jasonxyfile.Length()-5);
	  jasonxyfile+="_xylocator_xy.root";	  
	  cout<<"Processing the data file: "<<jasonxyfile<<endl;

	  events->Add(mainfile);	  
	  s2_fraction->Add(s2file);
	  masaxy->Add(masaxyfile);
	  jasonxy->Add(jasonxyfile);

	}
      }
    }
  events->AddFriend(s2_fraction);	  
  events->AddFriend(masaxy);
  events->AddFriend(jasonxy);
  
  TString outdir = "/darkside/users/hqian/AmBe2kn_XYData/";
  TString outfile;
  outfile.Form("%sAmBe2knXY_%d%d_%s.root",outdir.Data(),start,end,Time.Data());
  event_loop(events,outfile);

}
  
#ifndef __CINT__
int main(int argc, char **argv)
{
  theApp = new TRint("theApp",&argc,argv,NULL,0);
  int start, end;
  if(theApp->Argc() == 2)
    {
      start = atoi(theApp->Argv(1));
      end = start;
      AmBe2kn_XY(start,end);  
     }
  else if(theApp->Argc() == 3)
    {
      start = atoi(theApp->Argv(1));
      end = atoi(theApp->Argv(2));
      AmBe2kn_XY(start,end);  
     }
  else{
    cout<<"Usage: ./DSTAway startfile endfile "<<endl;
    cout<<"Usage: ./DSTAway startfile "<<endl;
  }
  std::cout << "==> Application finished." << std::endl;
  return 0;
}
#endif /* __CINT __ */



