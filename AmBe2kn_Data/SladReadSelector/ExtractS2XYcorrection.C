#include <iostream>
#include <fstream>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TApplication.h"
#include "TRint.h"
#include "TH1.h"
#include "TH2D.h"
#include "TF1.h"
#include "TMath.h"
#include "TColor.h"
#include "TLegend.h"
#include "TObjString.h"
#include "TVirtualFitter.h"
#include <TGaxis.h>
#include <TLatex.h>
#include "TPaveStats.h"
#include "TProfile2D.h"
#include "TEllipse.h"

using namespace std;
//#include "./HistogramSelector/HistogramSelector.h"

#ifndef __CINT__
  int main(int argc, char **argv);
#endif

  const Int_t N_CHANNELS=38;
  const Double_t pmtUnit = 4.1; //in cm   center to center distance 8.2cm /2
  const Double_t tpcRadi = 17.945;//35.89/2.; //in cm   center to center distance 8.2cm /2, 35.17 is 0.98% of 35.89
  const Double_t pmt_x[N_CHANNELS] = {  0.,-1.73,-3.46,-3.46,-3.46,-1.73,-1.73, 0.,1.73,3.46,1.73, 0., 0.,-1.73,  0,1.73,1.73,3.46,3.46,-3.46,-1.73, 0.,1.73,3.46,1.73, 0.,-1.73,-3.46,-3.46,-1.73,0.,1.73,3.46,3.46,1.73,0.,-1.73, 0. };
  const Double_t pmt_y[N_CHANNELS] = { -4.,  -3.,  -2.,   0.,   2.,   1.,  -1.,-2., -3., -2., -1., 0., 2.,   3., 4.,  3.,  1.,  0.,  2.,  -2.,  -3.,-4., -3., -2., -1.,-2.,  -1.,   0.,   2.,   1.,0.,  1.,  0.,  2.,  3.,2.,   3., 4. };
  const Double_t pmtRadi = 3.81;//position of pmts // 3.81 cm is pmt radius


TRint* theApp;
vector<int> colors;
void SetColors();

Int_t GetIntegrate(TH1* h, Double_t xmin, Double_t xmax);
TH1* GetProjection(TH2 * h, Double_t xmin, Double_t xmax);
TH1* GetProjRegion(TH1* h, Double_t xmin, Double_t xmax);
void NormalizeAtCenter(TH2 &h2D);
TH2D* ExtendHisto(TH2 *h2D);
void GetNNeighborBinsAndAverage(TH2 *h2D, Int_t xbin, Int_t ybin, Int_t &nbin, Double_t &avrg, Double_t &dev);
//TH2D * SubtractBGProf(TProfile2D *p1, TProfile2D *p2, Double_t c2);
TH2D * SubtractBGProf(TProfile2D *p1, TProfile2D *p2, Double_t c2, TH2 &hErr, TH2 &hEntry);
void CompareS2XYCorrFactors(TString fJason, TString fMasa);

void CopyEventNumberAndLivetime(TString fOrig, TString fDest){

  TFile *_file0 = TFile::Open(fDest.Data(), "UPDATE");
  TFile *_file1 = TFile::Open(fOrig.Data());
  TH1 *hRunTime = (TH1*) _file1->Get("hRunTime");
  TH1 *hEventConter = (TH1*) _file1->Get("hEventConter");
  hRunTime->SetDirectory(_file0);
  hEventConter->SetDirectory(_file0);
  _file0->cd();
  hRunTime->Write();
  hEventConter->Write();
  _file0->Save();
  _file0->ls();
}

//____________________________________________________________________________________________________
Double_t getNormalization(TH1* h, Double_t xmin, Double_t xmax)
{
  /// Get normalization of histogram
  /// Normalization = 1/(Nevts * binWidth)

  const Double_t area   = h->Integral(h->FindBin(xmin), h->FindBin(xmax)) ;
//  const Double_t area   = h->Integral(h->FindBin(xmin), h->FindBin(xmax), "width") ;
  if( area == 0 ) return 1.0 ;

  const Double_t binWidth = h->GetBinWidth(1) ;
  if( binWidth == 0.0 ) return 1.0 ;

  return 1.0/(area*binWidth);
}


void ExtractS2XYcorrection(TString fInNameSource, TString fInNameBGRun) {
  SetColors();
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(1);
//  gStyle->SetOptStat(1100);
  gStyle->SetOptTitle(1);
  TGaxis::SetMaxDigits(3);

  Bool_t drawInOnePad(false);//true);
  Bool_t useBG = true;
  Bool_t SubtractBG=true; //if you want to draw Source & Bg together, set this to false
  Bool_t fit = true;

    //-----Setting for s1 vs s2 for Kr-------------------
//  TString histoName("pS2VsXY_Kr");// histogram name you want to get;
  TString histoName("pR2_theta_Kr");// histogram name you want to get;

        //Normalize by integration
//    #define NORMINTEGRAL
//  Int_t rebinX(10), rebinY(10);
  Int_t rebinX(4), rebinY(4);
        Double_t nomMin(400.), normMax(550.);
        Double_t xmin(0.), xmax(800.), ymin(0.), ymax(20.e+3);
        Double_t x_int(2.), x_mean(290.), x_sigma(19), x_fitmin(150), x_fitmax(450); // fit initial values
        Double_t y_int(160.),y_mean(8.5e+3), y_sigma(1.4e+3), y_fitmin(2.e+3), y_fitmax(16.e+3);
        Bool_t autoaxisscale_x(true), autoaxisscale_y(true);
    //------Setting for s1 vs s2 for Kr------------------

        TCanvas *c1 = new TCanvas("c1","c1",1000, 750);
        c1->cd();
      //  c1->Draw();
        if(drawInOnePad)c1->Divide(2,3);
        (drawInOnePad)? c1->cd(1): c1->cd();

  cout<<"Input Source file: "<<fInNameSource.Data()<<endl;


  TFile *hfileSource = new TFile(fInNameSource.Data());
  if(!hfileSource) cout<<"file: "<<fInNameSource.Data()<<" is not found."<<endl;

  TProfile2D* p2D_Source = (TProfile2D*) hfileSource->Get(histoName.Data());
  if(!p2D_Source) cout<<"hist: "<<histoName.Data()<<" is not found."<<endl;
  p2D_Source->SetDirectory(0);
  p2D_Source->Sumw2();
  p2D_Source->SetName(Form("%s_Source", p2D_Source->GetName()));
//  p2D_Source = (TProfile2D*)p2D_Source->Rebin2D(2,2);

//  TString s2ovrs1name = "hTotalS2CorrvsS1Corr";
//  TString s2ovrs1name = "hTotalS2XYCorrvsS1Corr";
  TString s2ovrs1name = "hTotalS2CorrvsS1";
  TH2D *hs2vss1_Source = (TH2D*) hfileSource->Get(s2ovrs1name.Data());
  if(!hs2vss1_Source) cout<<"hist: "<<histoName.Data()<<" is not found."<<endl;
  hs2vss1_Source->SetDirectory(0);
  hs2vss1_Source->Sumw2();
  hs2vss1_Source->SetName(Form("%s_Source", hs2vss1_Source->GetName()));
  hs2vss1_Source = (TH2D*)hs2vss1_Source->Rebin2D(rebinX,rebinY);

  TString hname = "hRunTime";
  TH1D *hRunTime_Source = (TH1D*) hfileSource->Get(hname.Data());
  if(!hRunTime_Source) cout<<"hist: "<<hname.Data()<<" is not found."<<endl;
  hRunTime_Source->SetDirectory(0);
  hRunTime_Source->SetName(Form("%s_Source", hRunTime_Source->GetName()));


  hname = "hEventConter";
  TH1D *hEventConter_Source = (TH1D*) hfileSource->Get(hname.Data());
  if(!hEventConter_Source) cout<<"hist: "<<hname.Data()<<" is not found."<<endl;
  hEventConter_Source->SetDirectory(0);
  hEventConter_Source->SetName(Form("%s_Source", hEventConter_Source->GetName()));
  Int_t Nevents_Source = hEventConter_Source->GetBinContent(2);

  Double_t livetime_Source(1.);
  if(((TString)hRunTime_Source->GetXaxis()->GetBinLabel(2)).EqualTo("Live Time"))
    livetime_Source = hRunTime_Source->GetBinContent(2); //bincontent is s
//  livetime_Source = hRunTime_Source->GetBinContent(2)*1.e-9; //bincontent is ns
  else
    cout<<"Bin is not Live time bin!! Please check input histogram."<<endl;

  cout<<"LiveTime is "<<livetime_Source<<" ns."<<endl;

  if(fInNameBGRun.EqualTo("")) useBG = false;

  cout<<"Input Bg file: "<<fInNameBGRun.Data()<<endl;

  TProfile2D* p2D_Bg = NULL;
  TH2D *hs2ovrs1_Bg = NULL;
  Double_t livetime_Bg(1.);
  if(useBG){
      TFile *hfileBg = new TFile(fInNameBGRun.Data());
      if(!hfileBg) cout<<"file: "<<fInNameBGRun.Data()<<" is not found."<<endl;
      p2D_Bg = (TProfile2D*) hfileBg->Get(histoName.Data());
      if(!p2D_Bg) cout<<"hist: "<<histoName.Data()<<" is not found."<<endl;
      p2D_Bg->SetDirectory(0);
      p2D_Bg->Sumw2();
      p2D_Bg->SetName(Form("%s_Bg", p2D_Bg->GetName()));
//      p2D_Bg = (TProfile2D*)p2D_Bg->Rebin2D(2,2);

      hs2ovrs1_Bg = (TH2D*) hfileBg->Get(s2ovrs1name.Data());
      if(!hs2ovrs1_Bg) cout<<"hist: "<<histoName.Data()<<" is not found."<<endl;
      hs2ovrs1_Bg->SetDirectory(0);
      hs2ovrs1_Bg->Sumw2();
      hs2ovrs1_Bg->SetName(Form("%s_Bg", hs2ovrs1_Bg->GetName()));
      hs2ovrs1_Bg = (TH2D*)hs2ovrs1_Bg->Rebin2D(rebinX,rebinY);


      hname = "hRunTime";
      TH1D *hRunTime_Bg = (TH1D*) hfileBg->Get(hname.Data());
      if(!hRunTime_Bg) cout<<"hist: "<<hname.Data()<<" is not found."<<endl;
      hRunTime_Bg->SetDirectory(0);
      hRunTime_Bg->SetName(Form("%s_Bg", hRunTime_Bg->GetName()));

      if(((TString)hRunTime_Bg->GetXaxis()->GetBinLabel(2)).EqualTo("Live Time"))
        livetime_Bg = hRunTime_Bg->GetBinContent(2); //bincontent is s
//      livetime_Bg = hRunTime_Bg->GetBinContent(2)*1.e-9; //bincontent is ns
      else
        cout<<"Bin is not Live time bin!! Please check input histogram."<<endl;
  }



  (drawInOnePad)? c1->cd(1)->SetGrid(): c1->SetGrid();
  TH2D* htmp = (TH2D*)p2D_Source->Clone("colz");
  htmp->SetStats(0);
//  htmp->SetTitleOffset(1.2, "Y");
  htmp->Draw("colz");
  cout<<"Kr+Bg mean: "<< htmp->GetBinContent(htmp->FindBin(0.5, 0.5))<<endl;
  TEllipse el(0., 0., tpcRadi);
  el.SetFillStyle(0);
  el.Draw("same");
  TLatex L;
//  L.SetNDC();
  L.SetTextSize(0.03);
  TEllipse* pmt[N_CHANNELS/2];
  for (Int_t ch = N_CHANNELS/2; ch < N_CHANNELS; ch++)
  {
      pmt[ch%(N_CHANNELS/2)] = new TEllipse(pmtUnit*pmt_x[ch], pmtUnit*pmt_y[ch], pmtRadi);//position of pmts // 3.81 cm is pmt radius
      pmt[ch%(N_CHANNELS/2)]->SetFillStyle(0);
      pmt[ch%(N_CHANNELS/2)]->SetLineStyle(2);
      pmt[ch%(N_CHANNELS/2)]->Draw("same");
      L.DrawLatex(pmtUnit*pmt_x[ch], pmtUnit*pmt_y[ch], Form("%d",ch));
  }
  if(!drawInOnePad) {
      c1->Update();
      theApp->Run();
      c1->Print(Form("Source_BeforeBGsubt_%s.pdf", histoName.Data()));
  }

//  p2D_Bg->Scale(livetime_Source/livetime_Bg);
  (drawInOnePad)? c1->cd(2)->SetGrid(): c1->SetGrid();
  htmp = (TH2D*)p2D_Bg->Clone("colz");
  htmp->SetStats(0);
//  htmp->SetTitleOffset(1.2, "Y");
  htmp->Draw("colz");
  cout<<"Bg mean: "<< htmp->GetBinContent(htmp->FindBin(0.5, 0.5))<<endl;
  el.Draw("same");
  for (Int_t ch = N_CHANNELS/2; ch < N_CHANNELS; ch++) pmt[ch%(N_CHANNELS/2)]->Draw("same");
  for (Int_t ch = N_CHANNELS/2; ch < N_CHANNELS; ch++) L.DrawLatex(pmtUnit*pmt_x[ch], pmtUnit*pmt_y[ch], Form("%d",ch));
  if(!drawInOnePad) {
      c1->Update();
      theApp->Run();
      c1->Print(Form("BG_%s.pdf", histoName.Data()));
  }

  TH2D* h2D_Source_Bgsub;
  TH2* hErr(NULL), *hEntry(NULL);
  if(SubtractBG && useBG){
    //Subtract normalized BG.
      hErr = (TH2D *)p2D_Source->ProjectionXY();
      hEntry = (TH2D *)p2D_Source->ProjectionXY();

      h2D_Source_Bgsub = SubtractBGProf(p2D_Source, p2D_Bg, -livetime_Source/livetime_Bg, *hErr, *hEntry);
      h2D_Source_Bgsub->SetTitle("Kr S2 mean after BG subtraction; x [cm]; y [cm]; total_s2_corr [PE]");
      cout<<"Scaled: "<<livetime_Source/livetime_Bg<<endl;
  }


  (drawInOnePad)? c1->cd(3)->SetGrid(): c1->SetGrid();
  htmp = (TH2D*)h2D_Source_Bgsub->Clone("colz");
  htmp->SetStats(0);
//  htmp->SetTitleOffset(1.2, "Y");
  htmp->Draw("colz");
  cout<<"Kr mean: "<< htmp->GetBinContent(htmp->FindBin(0.5, 0.5))<<endl;
  el.Draw("same");
  for (Int_t ch = N_CHANNELS/2; ch < N_CHANNELS; ch++) pmt[ch%(N_CHANNELS/2)]->Draw("same");
  for (Int_t ch = N_CHANNELS/2; ch < N_CHANNELS; ch++) L.DrawLatex(pmtUnit*pmt_x[ch], pmtUnit*pmt_y[ch], Form("%d",ch));
  if(!drawInOnePad) {
      c1->Update();
      theApp->Run();
      c1->Print(Form("Source_AfterBGsubt_%s.pdf", histoName.Data()));
  }

  if(!drawInOnePad) {
      if(!hErr) cout<<"no hist"<<endl;
      hErr->Draw("colz");
      c1->Update();
      theApp->Run();
      c1->Print(Form("hErr_%s.pdf", histoName.Data()));
      hEntry->Draw("colz");
      c1->Update();
      theApp->Run();
      c1->Print(Form("hEntry_%s.pdf", histoName.Data()));
  }


  (drawInOnePad)? c1->cd(4)->SetGrid(): c1->SetGrid();

  TH2D *hS2KrVsXY = (TH2D*)h2D_Source_Bgsub->Clone("hS2KrVsXY_norm_center");
  hS2KrVsXY = ExtendHisto((TH2 *)hS2KrVsXY);
  NormalizeAtCenter(*hS2KrVsXY);
//  hS2KrVsXY = ExtendHisto((TH2 *)hS2KrVsXY);
  hS2KrVsXY->SetName("hS2KrVsXY_norm_center");
  hS2KrVsXY->SetTitle("total_s2_corr from Kr normalized at center; x [cm]; y [cm];");
  hS2KrVsXY->SetStats(0);
  hS2KrVsXY->Draw("colz");
  el.Draw("same");
  for (Int_t ch = N_CHANNELS/2; ch < N_CHANNELS; ch++) pmt[ch%(N_CHANNELS/2)]->Draw("same");
  for (Int_t ch = N_CHANNELS/2; ch < N_CHANNELS; ch++) L.DrawLatex(pmtUnit*pmt_x[ch], pmtUnit*pmt_y[ch], Form("%d",ch));
  if(!drawInOnePad) {
      c1->Update();
      theApp->Run();
      c1->Print(Form("S2KrVsXY_norm_center_%s.pdf", histoName.Data()));
  }

 TH2D *h2D_Bg = ExtendHisto((TH2 *)p2D_Bg->ProjectionXY("h2D_Bg"));
 NormalizeAtCenter(*h2D_Bg);


  TLegend* leg = new TLegend(0.65, 0.65, 0.9, 0.85);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->SetFillColor(10);
//  leg->SetHeader("");


  Float_t S1corr_mean_Kr_peak = 290.6;
  Float_t S1corr_sigma_Kr_peak = 18.45;
  Float_t nsigma = 3.;
  Float_t s1_min = S1corr_mean_Kr_peak-nsigma*S1corr_sigma_Kr_peak;
  Float_t s1_max = S1corr_mean_Kr_peak+nsigma*S1corr_sigma_Kr_peak;


  (drawInOnePad)? c1->cd(5)->SetGrid(): c1->SetGrid();
  TH1* h2D_Source_px = hs2vss1_Source->ProjectionX();
  h2D_Source_px->SetTitle("total_s1_corr");
  h2D_Source_px->SetMarkerStyle(20);
  TH1* h2D_Bg_px = hs2ovrs1_Bg->ProjectionX();
  h2D_Bg_px->Scale(livetime_Source/livetime_Bg);
  h2D_Bg_px->SetLineColor(kRed);
  h2D_Bg_px->SetMarkerColor(kRed);
  h2D_Bg_px->SetMarkerStyle(21);
  TH1* htemp = (TH1*)h2D_Source_px->Clone("h_s1_Bgsubt");
  h2D_Source_px->SetYTitle(Form("Events / ( %3.1f PE)", h2D_Source_px->GetXaxis()->GetBinWidth(1)) );
  h2D_Source_px->SetTitleOffset(1.2, "Y");
  h2D_Source_px->Draw();
  h2D_Bg_px->Draw("same");
  htemp->Add(h2D_Bg_px, -1.);
  htemp->SetLineColor(kBlue);
  htemp->SetMarkerColor(kBlue);
  htemp->SetMarkerStyle(22);
  htemp->Draw("sameC");
  leg->AddEntry(h2D_Source_px, "Kr+^{39}Ar", "LP");
  leg->AddEntry(h2D_Bg_px, "^{39}Ar", "LP");
  leg->AddEntry(htemp, "Kr", "LP");
  leg->Draw();

  if(!drawInOnePad) {
      c1->Update();
      theApp->Run();
      c1->Print(Form("Source_S1_atCenter_%s.pdf", histoName.Data()));
  }


  (drawInOnePad)? c1->cd(6)->SetGrid(): c1->SetGrid();
//  TH1* h2D_Source_py = hs2vss1_Source->ProjectionY();
//  TH1* h2D_Bg_py = hs2ovrs1_Bg->ProjectionY();
  cout<<"s1_min: "<<s1_min<<" s1_max: "<<s1_max<<endl;
  TH1* h2D_Source_py = GetProjection(hs2vss1_Source, s1_min, s1_max);
  TH1* h2D_Bg_py = GetProjection(hs2ovrs1_Bg, s1_min, s1_max);
//  h2D_Source_px->SetTitle("total_s2_corr");

  h2D_Bg_py->Scale(livetime_Source/livetime_Bg);
  h2D_Source_py->SetMarkerStyle(20);
  h2D_Source_py->SetYTitle(Form("Events / ( %3.1f PE)", h2D_Source_py->GetXaxis()->GetBinWidth(1)) );
  h2D_Source_py->SetTitleOffset(1.2, "Y");
  h2D_Source_py->Clone()->Draw();
  cout<<"Kr+Bg mean: "<< h2D_Source_py->GetMean()<<endl;

//  h2D_Source_py->Draw();
  h2D_Bg_py->SetLineColor(kRed);
  h2D_Bg_py->SetMarkerColor(kRed);
  h2D_Bg_py->SetMarkerStyle(21);
  h2D_Bg_py->Draw("same");
  h2D_Source_py->Add(h2D_Bg_py, -1.);
  h2D_Source_py->SetLineColor(kBlue);
  h2D_Source_py->SetMarkerColor(kBlue);
  h2D_Source_py->SetMarkerStyle(22);
  h2D_Source_py->Draw("same");
  h2D_Source_py->SetName("h_s2_Bgsubt");
//  leg->AddEntry(h2D_Source_py, "Kr", "LP");
  leg->Draw();
  cout<<"Bg mean: "<< h2D_Bg_py->GetMean()<<endl;
  cout<<"Kr mean: "<< h2D_Source_py->GetMean()<<endl;

  if(!drawInOnePad) {
      c1->Update();
      theApp->Run();
      c1->Print(Form("Source_S2_atCenter_%s.pdf", histoName.Data()));
  }



  theApp->Run();

  if(drawInOnePad) c1->Print(Form("All_%s.pdf", histoName.Data()));
//  c1->cd(1)->Print(Form("SourcePeak_Edrift%sVOvercm_2d.pdf",Edrift.Data()));

  //Normalize
  htemp->Scale(1./htemp->Integral());
  h2D_Source_py->Scale(1./h2D_Source_py->Integral());


  TString outFileName;
  outFileName=Form("S2CorrectionFactor_Kr_%s.root", histoName.Data());
  TFile* f = new TFile(outFileName.Data(), "RECREATE");
  f->cd();
  p2D_Bg->Write();
  p2D_Source->Write();
  hS2KrVsXY->Write();
  h2D_Bg->Write();
  htemp->Write();
  h2D_Source_py->Write();
  c1->SetName("S2CorrMap");
  c1->SetTitle("S2 correction map");
  if(drawInOnePad) c1->Write();

  f->Write();
  f->Close();

  delete c1;
}

Int_t GetIntegrate(TH1* h, Double_t xmin, Double_t xmax){
  TAxis *xaxis = h->GetXaxis();
  Int_t lowbin = xaxis->FindBin(xmin+0.00001);
  Int_t hibin = xaxis->FindBin(xmax+0.00001);
  return h->Integral(lowbin, hibin);
}

TH1* GetProjection(TH2* h, Double_t xmin, Double_t xmax){
  TAxis *xaxis = h->GetXaxis();
  Int_t lowbin = xaxis->FindBin(xmin+0.00001);
  Int_t hibin = xaxis->FindBin(xmax+0.00001);
  cout<<"xlow: "<<xaxis->GetBinLowEdge(lowbin)<<" xup: "<<xaxis->GetBinUpEdge(hibin)<<endl;
  TH1* htmp = h->ProjectionY("_py", lowbin, hibin);
  htmp->SetTitle(Form("%s %s=[%4f, %4f]", htmp->GetXaxis()->GetTitle(), xaxis->GetTitle(), xmin, xmax));
  return htmp;
}

TH1* GetProjRegion(TH1* h, Double_t xmin, Double_t xmax){
//  TAxis *xaxis = h->GetXaxis();
//  Int_t lowbin = xaxis->FindBin(xmin+0.00001);
//  Int_t hibin = xaxis->FindBin(xmax+0.00001);
  TH1* htmp = (TH1*) h->Clone(Form("%s_fill", h->GetName()));
  htmp->SetAxisRange(xmin+0.00001, xmax+0.00001);
  htmp->SetFillColor(kGray);
  htmp->SetLineWidth(0);
  htmp->SetLineColor(kGray);
  return htmp;
}

void GetNNeighborBinsAndAverage(TH2 *h2D, Int_t xbin, Int_t ybin, Int_t &nbin, Double_t &avrg, Double_t &dev){
  nbin =0; avrg=0.; dev=0.;
  for(Int_t ibin=xbin-1; ibin<=xbin+1; ibin++){
    for(Int_t jbin=ybin-1; jbin<=ybin+1; jbin++){
        if(ibin==xbin && jbin==ybin) continue; // skip current bin
        Double_t val =h2D->GetBinContent(ibin, jbin);
//        if(val==0. || val>1.) continue; //skip empty bin
        if(val==0. || val!=val) continue; //skip empty bin // val!=val is true if val is nan
        nbin++;
        avrg += val;//h2D->GetBinContent(ibin,jbin);
        dev  += val*val;
    }
  }
  avrg/=1.*nbin;
  dev = TMath::Sqrt(dev/nbin-avrg*avrg);
}

TH2D* ExtendHisto(TH2 *h2D){
  Int_t nXbins = h2D->GetNbinsX();
  Int_t nYbins = h2D->GetNbinsY();
  Double_t xfirst = h2D->GetXaxis()->GetXmin();
  Double_t xlast = h2D->GetXaxis()->GetXmax();
  Double_t yfirst = h2D->GetYaxis()->GetXmin();
  Double_t ylast = h2D->GetYaxis()->GetXmax();
//  TH2D *h2D_Ext = (TH2D*) h2D->Clone(Form("%s_Extended", h2D->GetName()));
  TH2D *h2D_Ext = new TH2D(Form("%s_Extended", h2D->GetName()), Form("%s Extended", h2D->GetTitle()),nXbins, xfirst, xlast,  nYbins, yfirst, ylast);
  h2D->Copy(*h2D_Ext);
  Bool_t done(false);
  Int_t nloop(0);
  while(!done && nloop<100){
      done = true; nloop++;
      for(Int_t xbin=1; xbin<=nXbins; xbin++){
        for(Int_t ybin=1; ybin<=nYbins; ybin++){
//            if(h2D_Ext->GetBinContent(xbin, ybin)!=0.) continue;
//            if(h2D->GetBinContent(xbin, ybin)!=0.) {/*h2D_Ext->SetBinContent(xbin,ybin,h2D->GetBinContent(xbin, ybin)); */continue;}
            Double_t cont = h2D->GetBinContent(xbin, ybin);
//            if(cont!=0. && cont==cont) {/*h2D_Ext->SetBinContent(xbin,ybin,h2D->GetBinContent(xbin, ybin)); */continue;}
            Int_t nbin(0); Double_t avrg(0.); Double_t dev(0.);
            GetNNeighborBinsAndAverage(h2D_Ext, xbin, ybin, nbin, avrg, dev);
            if(nbin<3) {done = false; continue;}
            if(cont!=0. && cont==cont) {
                if(TMath::Abs(cont-avrg) > 3.*dev) h2D_Ext->SetBinContent(xbin,ybin,avrg);
                continue;}
//            cout<<"xbin: "<<xbin<<" ybin: "<<ybin<<" nbin: "<<nbin<<" average: "<<avrg<<endl;
            h2D_Ext->SetBinContent(xbin,ybin,avrg);
        }
      }
  }
  cout<<"Finish Extending. nloop: "<<nloop<<endl;
 return h2D_Ext;
}


void NormalizeAtCenter(TH2 &h2D){
  Int_t nXbins = h2D.GetNbinsX();
  Int_t nYbins = h2D.GetNbinsY();
  Double_t xfirst = h2D.GetXaxis()->GetXmin();
  Double_t xlast = h2D.GetXaxis()->GetXmax();
  Double_t yfirst = h2D.GetYaxis()->GetXmin();
  Double_t ylast = h2D.GetYaxis()->GetXmax();

  Int_t xbin_center = h2D.GetXaxis()->FindBin(0.);
  Int_t ybin_center = h2D.GetYaxis()->FindBin(0.);

  Double_t avg_value = h2D.GetBinContent(xbin_center,ybin_center);
  avg_value += h2D.GetBinContent(xbin_center-1,ybin_center);
  avg_value += h2D.GetBinContent(xbin_center-1,ybin_center-1);
  avg_value += h2D.GetBinContent(xbin_center,ybin_center-1);
  avg_value /= 4.;

/*  for(Int_t xbin=1; xbin<=nXbins; xbin++){
    for(Int_t ybin=1; ybin<=nYbins; ybin++){
        h2D.SetBinContent(xbin,ybin,h2D.GetBinContent(xbin,ybin)/avg_value);
    }
  }*/
  h2D.Scale(1./avg_value);

  return;
}

TH2D * SubtractBGProf(TProfile2D *p1, TProfile2D *p2, Double_t c2, TH2 &hErr, TH2 &hEntry){
  TH2D *p = (TH2D *)p1->ProjectionXY();
  Double_t ac2 = TMath::Abs(c2);
  // Make the loop over the bins to calculate the Addition
     Int_t bin;
     // create sumw2 per bin if not set
     p->Sumw2();
     // if p1 has not the sum of weight squared/bin stored use just the sum of weights
     Int_t ncells = p->GetSize();
     for (bin =0;bin< ncells;bin++) {
         if(p1->GetBinEntries(bin)==0 || p2->GetBinEntries(bin)==0) continue;
        Double_t sum   = p1->GetBinContent(bin)*p1->GetBinEntries(bin) + c2*p2->GetBinContent(bin)*p2->GetBinEntries(bin);
        Double_t entry = p1->GetBinEntries(bin) + c2*p2->GetBinEntries(bin); // this part is changed from original // http://www.slac.stanford.edu/grp/eg/minos/ROOTSYS/pro/hist/hist/src/TProfileHelper.h
        p->SetBinContent(bin, sum/entry);
        p->SetBinError(bin, sum/entry/TMath::Sqrt(entry));

        cout<<"Entry: "<<entry<<" Mean: "<<sum/entry<<" error: "<<sum/entry/TMath::Sqrt(entry)<<endl;
        hErr.SetBinContent(bin, 100./TMath::Sqrt(entry));
        hEntry.SetBinContent(bin, entry);
     }
 return p;
}

void CompareS2XYCorrFactors(TString fJason, TString fMasa){

  gStyle->SetOptStat(0);

  TCanvas *c1 = new TCanvas("c1","c1",1000, 750);
  c1->cd();
  c1->Divide(2,2);

//  TFile *_file0 = TFile::Open(fJason.Data(), "UPDATE");
  TFile *_file0 = TFile::Open(fJason.Data());
  TFile *_file1 = TFile::Open(fMasa.Data());
  TH2 *h_Jason  = (TH2*) _file0->Get("hS2KrVsXY_norm_center");
  TH2 *h_Masa   = (TH2*) _file1->Get("hS2KrVsXY_norm_center");
  TH2 *h_Masa_Bg   = (TH2*) _file1->Get("h2D_Bg");
  if(!h_Masa_Bg) cout<<"no hist"<<endl;
  h_Masa->SetName("hS2KrVsXY_norm_center_Masa");
//  c1->cd(1)->SetGrid();
//  h_Jason->Draw("colz");
  c1->cd(2)->SetGrid();
  h_Masa->Draw("colz");

  TH1* htemp = (TH1*)h_Jason->Clone("htemp");
  c1->cd(1)->SetGrid();
  htemp->Draw("colz");

//  htemp->Add(h_Masa, -1.);
//  htemp->Draw("colz");
  c1->cd(3)->SetGrid();
//  h_Jason->Add(h_Masa, -1.);
//  h_Jason->Draw("colz");

  h_Masa_Bg->Add(h_Masa, -1.);
  h_Masa_Bg->Draw("colz");

  TEllipse el(0., 0., tpcRadi);
  el.SetFillStyle(0);
  el.Draw("same");
  TLatex L;
//  L.SetNDC();
  L.SetTextSize(0.03);
  TEllipse* pmt[N_CHANNELS/2];
  for (Int_t ch = N_CHANNELS/2; ch < N_CHANNELS; ch++)
  {
      pmt[ch%(N_CHANNELS/2)] = new TEllipse(pmtUnit*pmt_x[ch], pmtUnit*pmt_y[ch], pmtRadi);//position of pmts // 3.81 cm is pmt radius
      pmt[ch%(N_CHANNELS/2)]->SetFillStyle(0);
      pmt[ch%(N_CHANNELS/2)]->SetLineStyle(2);
      pmt[ch%(N_CHANNELS/2)]->Draw("same");
      L.DrawLatex(pmtUnit*pmt_x[ch], pmtUnit*pmt_y[ch], Form("%d",ch));
  }


  c1->Update();
  theApp->Run();

  c1->Print("Diff_Jason_Masa.pdf");

  delete c1;

 //move to Jason's file
//  h_Masa->SetDirectory(_file0);
//  _file0->cd();
//  h_Masa->Write();
//  _file0->Save();
//  _file0->ls();
}

void SetColors(){


    colors.push_back(TColor::GetColor("#E6855E"));
    colors.push_back(TColor::GetColor("#44A5CB"));
    colors.push_back(TColor::GetColor("#D45D87"));
    colors.push_back(TColor::GetColor("#40BFB0"));
    colors.push_back(TColor::GetColor("#F9DB57"));
    colors.push_back(TColor::GetColor("#9D73BB"));
    colors.push_back(TColor::GetColor("#009F8C"));
    colors.push_back(TColor::GetColor("#8B90BE"));
    colors.push_back(TColor::GetColor("#F3C0AB"));

}

#ifndef __CINT__
int main(int argc, char **argv) {
	theApp = new TRint("App", &argc, argv, NULL, 0);
        theApp->Connect("KeyPressed(Int_t)","TSystem",gSystem,"ExitLoop()");
	if (theApp->Argc() == 3) {
            std::cout << "\n==========> ExtractS2XYcorrection without Bg Runs!!! <=============" << std::endl;
            std::cout << "==> Application start." << std::endl;
            ExtractS2XYcorrection(theApp->Argv(1), theApp->Argv(2));

	} else if(theApp->Argc() == 2) {
            std::cout << "\n==========> ExtractS2XYcorrection without Bg Runs!!! <=============" << std::endl;
            std::cout << "==> Application start." << std::endl;
//            CompareS2XYCorrFactors("S2CorrectionFactor_Kr_pS2VsXY_Kr_Jasons_xy.root", "S2CorrectionFactor_Kr_pS2VsXY_Kr_Masa.root");
//            CompareS2XYCorrFactors("S2CorrectionFactor_Kr_pS2VsXY_Kr.root", "S2CorrectionFactor_Kr_pS2VsXY_Kr_Masa.root");
            CompareS2XYCorrFactors("S2CorrectionFactor_Kr_pS2VsXY_Kr.root", "S2CorrectionFactor_Kr_pS2VsXY_Kr.root");

        } else {
		std::cout << "Usage:" <<theApp->Argc()<< std::endl;
		std::cout << "./ExtractS2XYcorrection fInName fInName_Bg" << std::endl;
		return 0;
	}


//        std::cout << "==> Application start." << std::endl;
//	ExtractS2XYcorrection(theApp->Argv(1), theApp->Argv(2));
	std::cout << "==> Application finished." << std::endl;
//	theApp->Run();

	return 0;
}
#endif /* __CINT __ */
