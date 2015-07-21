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
using namespace std;
//#include "./HistogramSelector/HistogramSelector.h"

#ifndef __CINT__
  int main(int argc, char **argv);
#endif

TF1 *funcS2Source;
TRint* theApp;
vector<int> colors;
void SetColors();
Int_t fitLineColor;

Double_t GetFractionOfBG(TF1 *f, Double_t xmin, Double_t xmax);
Int_t GetIntegrate(TH1* h, Double_t xmin, Double_t xmax);
TH1* GetProjection(TH2 * h, Double_t xmin, Double_t xmax);
TH1* GetProjRegion(TH1* h, Double_t xmin, Double_t xmax);
Double_t g2(Double_t *x, Double_t *par);
Double_t g2Pol1(Double_t *x, Double_t *par);
Double_t lognormal(Double_t *x, Double_t *par);
Double_t myfunc(Double_t *x, Double_t *par);
Double_t myfuncX(Double_t *x, Double_t *par);
void RemoveWhitespace(TString &s);

Double_t GetPeakError2(Double_t peak);
Double_t dMeandPar(int par, TF1* f, Double_t peak);

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

void DrawDistributionsFromAllFields(TString hname, TString f100Vovercm, TString f150Vovercm, TString f200Vovercm){
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  SetColors();
  TFile *_file0 = TFile::Open(f100Vovercm.Data());
  TFile *_file1 = TFile::Open(f150Vovercm.Data());
  TFile *_file2 = TFile::Open(f200Vovercm.Data());
  TH1 *h100Vcm = (TH1*) _file0->Get(hname.Data());
  h100Vcm->SetName("100Vcm");
  TH1 *h150Vcm = (TH1*) _file1->Get(hname.Data());
  h150Vcm->SetName("150Vcm");
  TH1 *h200Vcm = (TH1*) _file2->Get(hname.Data());
  h200Vcm->SetName("200Vcm");
  _file0->ls();
//  h100Vcm->SetDirectory(0);
//  h150Vcm->SetDirectory(0);
//  h200Vcm->SetDirectory(0);
  Int_t nrebin(2);
  h100Vcm->Rebin(nrebin);
  h150Vcm->Rebin(nrebin);
  h200Vcm->Rebin(nrebin);
//  Double_t xmin = h100Vcm->GetXaxis()->GetXmin();
//  Double_t xmax = h100Vcm->GetXaxis()->GetXmax();
//  Double_t xmin = 35.;
//  Double_t xmax = 48.;
  Double_t xmin = 250.;
  Double_t xmax = 360.;
  cout<<"xmin: "<<xmin<<" xmax: "<<xmax<<endl;
//  Double_t norm = h100Vcm->Integral();//getNormalization(h100Vcm, xmin, xmax);
  Double_t norm = getNormalization(h100Vcm, xmin, xmax);
  h100Vcm->Scale(norm);
  cout<<"norm: "<<norm<<endl;
//  norm = h150Vcm->Integral();//getNormalization(h150Vcm, xmin, xmax);
  norm = getNormalization(h150Vcm, xmin, xmax);
  h150Vcm->Scale(norm);
  cout<<"norm: "<<norm<<endl;
//  norm = h200Vcm->Integral();//getNormalization(h200Vcm, xmin, xmax);
  norm = getNormalization(h200Vcm, xmin, xmax);
  h200Vcm->Scale(norm);
  cout<<"norm: "<<norm<<endl;
  gStyle->SetOptStat(0);

  h100Vcm->SetTitle("100 V/cm");
  h150Vcm->SetTitle("150 V/cm");
  h200Vcm->SetTitle("200 V/cm");
  h100Vcm->GetYaxis()->SetTitle("A.U.");

  h100Vcm->SetLineColor(colors.at(1));
  h150Vcm->SetLineColor(colors.at(2));
  h200Vcm->SetLineColor(colors.at(3));
  h100Vcm->SetMarkerColor(colors.at(1));
  h150Vcm->SetMarkerColor(colors.at(2));
  h200Vcm->SetMarkerColor(colors.at(3));
  h100Vcm->SetMarkerStyle(20);
  h150Vcm->SetMarkerStyle(21);
  h200Vcm->SetMarkerStyle(22);

  TCanvas *c1 = new TCanvas("c1","c1",1000, 750);
  c1->SetGrid();
  c1->cd();

  h100Vcm->Draw("C");
//  h150Vcm->Draw("Csame");
  h200Vcm->Draw("Csame");
  TLegend* lg= c1->BuildLegend();
  lg->SetHeader("^{39}Ar");
  c1->Update();
  c1->Print(Form("%s.pdf",hname.Data()));

}


void ExtractBaSpectum(TString Edrift, TString fInNameSource, TString fInNameBGRun) {
  SetColors();
  gStyle->SetOptFit(1);
//  gStyle->SetOptStat(0);
  gStyle->SetOptStat(1100);
  gStyle->SetOptTitle(0);
  TGaxis::SetMaxDigits(3);
  fitLineColor= colors.at(1);

  Bool_t drawInOnePad(true);
  Bool_t useBG = true;
  Bool_t SubtractBG=true; //if you want to draw Source & Bg together, set this to false
  Bool_t fit = true;

  //  Double_t nomMin(600), normMax(3500);
  //  Double_t nomMin(400.), normMax(500.);

  //  Int_t rebinX(4), rebinY(2);
/*
//-----Setting for s1 vs s2 for Ba-------------------
    TString histoName("hTotalS2XYCorrvsS1Corr");// histogram name you want to get;
    //Normalize by integration
#define NORMINTEGRAL
    Int_t rebinX(5), rebinY(4);
//    Double_t nomMin(2000.), normMax(3000.);
    Double_t nomMin(2500.), normMax(3000.); // 150 Ba
    Double_t xmin(0.), xmax(4.e+3), ymin(0.), ymax(40.e+3);
    Double_t x_int(1.e+6), x_mean(0.7e+3), x_sigma(0.4e+3), x_fitmin(0.4e+3), x_fitmax(1.e+3); // fit initial values
    Double_t y_int(1.e+6),y_mean(24.e+3), y_sigma(5.e+3), y_fitmin(15.e+3), y_fitmax(30.e+3);
    Bool_t autoaxisscale_x(false), autoaxisscale_y(false);
//------Setting for s1 vs s2 for Ba------------------
*/
/*    //-----Setting for Energy vs s1 for Ba-------------------
    TString histoName("hTotalS1CorrvsEnergy");// histogram name you want to get;
    //Normalize by integration
#define NORMINTEGRAL
  Int_t rebinX(5), rebinY(4);
        Double_t nomMin(400.), normMax(500.);
        Double_t xmin(0.), xmax(5.e+2), ymin(0.), ymax(4.e+3);
        Double_t x_int(1.e+6), x_mean(0.7e+2), x_sigma(0.4e+2), x_fitmin(0.4e+2), x_fitmax(1.e+2); // fit initial values
        Double_t y_int(1.e+6), y_mean(0.7e+3), y_sigma(0.4e+3), y_fitmin(0.4e+3), y_fitmax(1.e+3);               // fit initial values
        Bool_t autoaxisscale_x(false), autoaxisscale_y(false);
    //------Setting for Energy vs s1 for Ba------------------
*/

    //-----Setting for s1 vs s2 for Kr-------------------
        TString histoName("hTotalS2XYCorrvsS1Corr");// histogram name you want to get;
        //Normalize by integration
//    #define NORMINTEGRAL
        Int_t rebinX(2), rebinY(4);
        Double_t nomMin(400.), normMax(550.);
//        Double_t nomMin(2500.), normMax(3000.); // 150 Ba
        Double_t xmin(0.), xmax(800.), ymin(0.), ymax(20.e+3);
        Double_t x_int(2.), x_mean(290.), x_sigma(19), x_fitmin(150), x_fitmax(450); // fit initial values
//        Double_t x_int(2.), x_mean(290.), x_sigma(30), x_fitmin(150), x_fitmax(450); // fit initial values
        Double_t y_int(160.),y_mean(4.1e+3), y_sigma(0.83e+3), y_fitmin(2.e+3), y_fitmax(6.e+3); // 50 V/cm
//        Double_t y_int(160.),y_mean(8.5e+3), y_sigma(1.4e+3), y_fitmin(5.e+3), y_fitmax(10.e+3); // 200 V/cm
//        Double_t y_int(160.),y_mean(8.e+3), y_sigma(1.4e+3), y_fitmin(5.e+3), y_fitmax(12.e+3); // 200 V/cm
//        Double_t y_int(160.),y_mean(6.e+3), y_sigma(1.4e+3), y_fitmin(4.e+3), y_fitmax(8.e+3); // 100 V/cm
//        Double_t y_int(160.),y_mean(6.8e+3), y_sigma(1.35e+3), y_fitmin(2.e+3), y_fitmax(16.e+3); // 150 V/cm
        Bool_t autoaxisscale_x(true), autoaxisscale_y(true);
    //------Setting for s1 vs s2 for Kr------------------

/*          //-----Setting for Energy vs s1 for Ba-------------------
            TString histoName("hTotalS1CorrvsEnergy");// histogram name you want to get;
            //Normalize by integration
//        #define NORMINTEGRAL
          Int_t rebinX(1), rebinY(1);
                Double_t nomMin(400.), normMax(500.);
                Double_t xmin(0.), xmax(1.e+2), ymin(0.), ymax(1.e+3);
                Double_t x_int(0.5), x_mean(41.5), x_sigma(5), x_fitmin(30.), x_fitmax(50.); // fit initial values
                Double_t y_int(2), y_mean(290), y_sigma(20), y_fitmin(150), y_fitmax(450);               // fit initial values
                Bool_t autoaxisscale_x(true), autoaxisscale_y(true);
            //------Setting for Energy vs s1 for Kr------------------
*/
  //  TString histoName("hTotalS2vsS1");// histogram name you want to get;
  //  TString histoName("hTotalS2CorrvsS1");// histogram name you want to get;
  //  TString histoName("hTotalS1CorrvsS2CorroverS1Corr");// histogram name you want to get;
  //  TString histoName("hTotalS2CorrvsS1Corr");// histogram name you want to get;
//    TString histoName("hTotalS2XYCorrvsS1Corr");// histogram name you want to get;
  //  TString histoName("hTotalS1CorrvsEnergy");// histogram name you want to get;
  //  TString histoName("hTotalS2XYCorrvsEnergy");// histogram name you want to get;
  //  TString histoName("hTotalS1vsF90");// histogram name you want to get;
  //  TString histoName("total_s1_f90_hist");// histogram name you want to get;
  //  TString histoName("total_s1_corr_f90_hist");// histogram name you want to get;
  //  TString histoName("t_drift_s1_corr_hist");// histogram name you want to get;

//  gStyle->SetStatX(0.52);
//  gStyle->SetStatY(0.6);
//  gStyle->SetStatW(0.4);
//  gStyle->SetStatH(0.3);

//  TString fHistOutName(fInName);
//  if (fHistOutName.Contains('/'))
//    fHistOutName.Replace(0, fHistOutName.Last('/') + 1, fOutDir + "/Hist");
//  else
//    fHistOutName.Prepend(fOutDir + "/Hist");
//  std::cout << "outhistfile: " << fHistOutName.Data() << std::endl;
  cout<<"Input Source file: "<<fInNameSource.Data()<<endl;

  if(!fInNameSource.Contains("Drift"+ Edrift +"Vcm"))
    cout<<"Drift field is not matched with file name!!"<<endl;
  if(!fInNameBGRun.EqualTo("") && !fInNameBGRun.Contains("Drift"+ Edrift +"Vcm"))
    cout<<"Drift field is not matched with file name!!"<<endl;


  TFile *hfileSource = new TFile(fInNameSource.Data());
  if(!hfileSource) cout<<"file: "<<fInNameSource.Data()<<" is not found."<<endl;

  TH2D *h2D_Source = (TH2D*) hfileSource->Get(histoName.Data());
  if(!h2D_Source) cout<<"hist: "<<histoName.Data()<<" is not found."<<endl;
  h2D_Source->SetDirectory(0);
  h2D_Source->Sumw2();
  h2D_Source->SetName(Form("%s_Source", h2D_Source->GetName()));

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
//    livetime_Source = hRunTime_Source->GetBinContent(2); //bincontent is s
  livetime_Source = hRunTime_Source->GetBinContent(2)*1.e-9; //bincontent is ns
  else
    cout<<"Bin is not Live time bin!! Please check input histogram."<<endl;

  cout<<"LiveTime is "<<livetime_Source<<" ns."<<endl;

  if(fInNameBGRun.EqualTo("")) useBG = false;

  cout<<"Input Bg file: "<<fInNameBGRun.Data()<<endl;

  TH2D *h2D_Bg = NULL;
  Double_t livetime_Bg(1.);
  if(useBG){
      TFile *hfileBg = new TFile(fInNameBGRun.Data());
      if(!hfileBg) cout<<"file: "<<fInNameBGRun.Data()<<" is not found."<<endl;
      h2D_Bg = (TH2D*) hfileBg->Get(histoName.Data());
      if(!h2D_Bg) cout<<"hist: "<<histoName.Data()<<" is not found."<<endl;
      h2D_Bg->SetDirectory(0);
      h2D_Bg->Sumw2();
      h2D_Bg->SetName(Form("%s_Bg", h2D_Bg->GetName()));
/*      h2D_Bg->SetXTitle("total_s1 [PE]");
      h2D_Bg->SetYTitle(ytitle);
      if(histoName.Contains("_corr")) h2D_Bg->SetXTitle("total_s1_corr [PE]");
*/
      hname = "hRunTime";
      TH1D *hRunTime_Bg = (TH1D*) hfileBg->Get(hname.Data());
      if(!hRunTime_Bg) cout<<"hist: "<<hname.Data()<<" is not found."<<endl;
      hRunTime_Bg->SetDirectory(0);
      hRunTime_Bg->SetName(Form("%s_Bg", hRunTime_Bg->GetName()));

      if(((TString)hRunTime_Bg->GetXaxis()->GetBinLabel(2)).EqualTo("Live Time"))
//        livetime_Bg = hRunTime_Bg->GetBinContent(2); //bincontent is s
      livetime_Bg = hRunTime_Bg->GetBinContent(2)*1.e-9; //bincontent is ns
      else
        cout<<"Bin is not Live time bin!! Please check input histogram."<<endl;
  }



#ifdef NORMINTEGRAL
  // normalized by integration
  TH1* htmp_Source_S1 = h2D_Source->ProjectionX();
  Double_t int_Source = GetIntegrate(htmp_Source_S1, nomMin, normMax);

  TH1* htmp_Bg_S1 = h2D_Bg->ProjectionX();
  Double_t int_Bg = GetIntegrate(htmp_Bg_S1, nomMin, normMax);

  Double_t normfactor = int_Source/int_Bg;

  if(normfactor!=0.)
    h2D_Bg->Scale(normfactor);
  cout<<"normfactor: "<<normfactor<<endl;

#endif
  //normalized by live time
  if(livetime_Source!=0.)
    h2D_Source->Scale(1./livetime_Source);

#ifdef NORMINTEGRAL
  if(livetime_Bg!=0. && useBG)
    h2D_Bg->Scale(1./livetime_Source);
#else
  if(livetime_Bg!=0. && useBG)
    h2D_Bg->Scale(1./livetime_Bg);
#endif

  TH2F* h2D_Source_scaled = (TH2F*) h2D_Source->Clone(histoName+"_Source_scaled");


  if(SubtractBG && useBG){
    //Subtract normalized BG.
      h2D_Source->Add(h2D_Bg,-1.);
  }
  h2D_Source = (TH2D*)h2D_Source->Rebin2D(rebinX,rebinY);
  if(!SubtractBG && useBG) h2D_Bg = (TH2D*)h2D_Bg->Rebin2D(rebinX,rebinY);

  TCanvas *c1 = new TCanvas("c1","c1",1000, 750);
  c1->cd();
//  c1->Draw();
  if(drawInOnePad)c1->Divide(2,2);
  (drawInOnePad)? c1->cd(1): c1->cd();
//  h2D_Source->SetAxisRange(0, 1500, "X");
//  h2D_Source->Draw("colz");
//  TH2D* htmp = (TH2D*)h2D_Source->DrawClone("colz");
  TH2D* htmp = (TH2D*)h2D_Source->Clone("colz");
  htmp->SetStats(0);
  htmp->SetTitleOffset(1.2, "Y");
//  htmp->Rebin2D(2,2);
//  htmp->SetAxisRange(0, 800);
//  htmp->SetAxisRange(0, 12000, "Y");
//  htmp->SetAxisRange(0, 100);
//  htmp->SetAxisRange(0, 800, "Y");
  htmp->SetAxisRange(xmin, xmax); // energy variable
  htmp->SetAxisRange(ymin, ymax, "Y"); // s1
  htmp->Draw("colz");
  if(!drawInOnePad) {
      c1->Update();
      theApp->Run();
      c1->Print(Form("SourcePeak_Edrift%sVOvercm_%s.pdf",Edrift.Data(), histoName.Data()));
  }

  (drawInOnePad)? c1->cd(2)->SetGrid(): c1->SetGrid();
  TH1* h2D_Source_px = h2D_Source->ProjectionX();
  h2D_Source_px->SetTitle("");
  h2D_Source_px->SetTitleOffset(1.2, "Y");
  h2D_Source_px->SetLineWidth(2);
//  h2D_Source_px->Scale(1./h2D_Source_px->GetXaxis()->GetBinWidth(1));
  //  h2D_Source_px->SetYTitle(Form("Events/livetime/binwidth (%2.1f PE)", h2D_Source_px->GetXaxis()->GetBinWidth(1)));
    h2D_Source_px->SetYTitle(Form("Events/livetime (%2.1f PE)", h2D_Source_px->GetXaxis()->GetBinWidth(1)));
  h2D_Source_px->Draw();

  if(!SubtractBG && useBG){ //draw Source and BG together
      TH1* h2D_Bg_px = h2D_Bg->ProjectionX();
      h2D_Bg_px->SetLineColor(kRed);
      if(h2D_Bg_px) h2D_Bg_px->Draw("same");
      theApp->Run();
  }


  Double_t s1_mean = h2D_Source_px->GetBinCenter(h2D_Source_px->GetMaximumBin());
  Double_t s1_sigma = 16.;
  Double_t s1_meanErr(0.), s1_sigmaErr(0.);
//  TF1 *fS1Source = (useBG && SubtractBG)? new TF1("fS1Source",g2,x_fitmin,x_fitmax, 3): new TF1("fS1Source",g2Pol1,x_fitmin,x_fitmax, 5);
//  (useBG && SubtractBG)? fS1Source->SetParameters(2.2e+8, x_mean, x_sigma): fS1Source->SetParameters(1., x_mean, x_sigma, 0.04, 3.e-5);
//  (useBG && SubtractBG)? fS1Source->SetParNames("Counts", "Mean", "Sigma"): fS1Source->SetParNames("Counts", "Mean", "Sigma", "Bg const.", "Bg slope");
  TF1 *fS1Source = new TF1("fS1Source",g2,x_fitmin,x_fitmax, 3);
  fS1Source->SetParameters(x_int, x_mean, x_sigma);
  fS1Source->SetParNames("Counts", "Mean", "Sigma");
//  Double_t fitMin = (useBG)? s1_mean-nsigma*s1_sigma: s1_mean-nsigma*s1_sigma
  fS1Source->SetLineColor(kRed);
  fS1Source->SetLineColor(fitLineColor);
//  if (useBG && SubtractBG) h2D_Source_px->Fit(fS1Source, "EM", "", s1_mean-nsigmaFit*s1_sigma, s1_mean+nsigmaFit*s1_sigma);
//  else
  if(fit) {
      h2D_Source_px->Fit(fS1Source, "EMR");

  Double_t nsigmafit = 3.;
  s1_mean = fS1Source->GetParameter(1);
  s1_sigma = fS1Source->GetParameter(2);
  h2D_Source_px->Fit(fS1Source, "EM", "", s1_mean-nsigmafit*s1_sigma, s1_mean+nsigmafit*s1_sigma);


  s1_mean = fS1Source->GetParameter(1);
  s1_sigma = fS1Source->GetParameter(2);
  s1_meanErr = fS1Source->GetParError(1);
  s1_sigmaErr = fS1Source->GetParError(2);
  if(autoaxisscale_x) h2D_Source_px->SetAxisRange(s1_mean-5.*s1_sigma, s1_mean+8.5*s1_sigma, "X");
  }

  Double_t nsigmaCut = 1.5;
  Double_t xProjmin = s1_mean-nsigmaCut*s1_sigma;
  Double_t xProjmax = s1_mean+nsigmaCut*s1_sigma;
#if 1
#ifndef NORMINTEGRAL
  TH1* hfill = GetProjRegion(h2D_Source_px, xProjmin, xProjmax);
#else
  TH1* hfill = GetProjRegion(h2D_Source_px, nomMin, normMax);
#endif
  hfill->Draw("same hist");
  h2D_Source_px->Draw("same");
#endif
  TLegend* leg = new TLegend(0.65, 0.15, 0.9, 0.3);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetFillColor(10);
//  leg->SetHeader("");

  if(!SubtractBG && useBG){
      TH1* h2D_Bg_px = h2D_Bg->ProjectionX();
      h2D_Bg_px->SetLineColor(kRed);
      h2D_Bg_px->Draw("same");
//      leg->AddEntry(h2D_Source_px, "Source+Ar39", "L");
      leg->AddEntry(h2D_Source_px, "Ba+Ar39", "L");
      leg->AddEntry(h2D_Bg_px, "Ar39", "L");
      if(fit) leg->AddEntry(fS1Source, "Fit func.", "L");
#ifdef NORMINTEGRAL
      leg->AddEntry(hfill, "Norm. range", "F");
#endif

      leg->Draw();
//      theApp->Run();
  }
  gPad->Update();
  TPaveStats *st = (TPaveStats*)h2D_Source_px->FindObject("stats");
  st->SetX1NDC(0.52);st->SetX2NDC(0.9);
  st->SetY1NDC(0.6);st->SetY2NDC(0.9);
  st->Draw();

  h2D_Source_px->GetListOfFunctions()->Print();
  (drawInOnePad)? c1->cd(2)->Update(): c1->Update();
  if(!drawInOnePad) {
      c1->Update();
      theApp->Run();
      c1->Print(Form("SourcePeak_Edrift%sVOvercm_%s_px.pdf",Edrift.Data(), histoName.Data()));
  }
  (drawInOnePad)? c1->cd(3)->SetGrid(): c1->SetGrid();
  TH1* h2D_Source_py = GetProjection( h2D_Source, xProjmin, xProjmax);
//  TH1* h2D_Source_py = h2D_Source->ProjectionY();
  if(h2D_Source_py) {
      h2D_Source_py->SetYTitle(Form("Events/livetime/ (%1.2f)", h2D_Source_py->GetXaxis()->GetBinWidth(1)));
      h2D_Source_py->SetTitleOffset(1.3, "Y");
      h2D_Source_py->SetLineWidth(2);
      h2D_Source_py->SetFillColor(kGray);
//      h2D_Source_py->Draw("HIST");
      h2D_Source_py->Draw();
  }
  if(!drawInOnePad) {
      c1->Update();
      theApp->Run();
      c1->Print(Form("SourcePeak_Edrift%sVOvercm_%s_py.pdf",Edrift.Data(), histoName.Data()));
  }


  Double_t bgFrac(0.);
  if(!useBG || !SubtractBG) bgFrac = GetFractionOfBG(fS1Source, xProjmin, xProjmax);
  TF1 *fS2Source = new TF1("fS2Source",g2,y_fitmin,y_fitmax, 3); // 150V/cm for xy-corrected S2 from Source
//    fS2Source->SetNpx(100);
    fS2Source->SetParNames("int", "mean", "sigma");
    fS2Source->SetParameters(58., y_mean, y_sigma);
//    fS2Source->SetParameters(1.6e-3, 0.25, 0.03);// for f90 fit
//    fS2Source->SetParLimits(0, 1, 1.e+3);
    fS2Source->SetLineColor(fitLineColor);

    h2D_Source_py->Fit(fS2Source, "EMR");//, "", s1_mean-nsigma*s1_sigma, s1_mean+nsigma*s1_sigma);
    Double_t s2_mean = fS2Source->GetParameter(1);
    Double_t s2_sigma = fS2Source->GetParameter(2);
    Double_t nsigmafit = 3.;
    h2D_Source_py->Fit(fS2Source, "EM", "", s2_mean-nsigmafit*s2_sigma, s2_mean+nsigmafit*s2_sigma);
    fS2Source->SetRange(s2_mean-nsigmafit*s2_sigma, s2_mean+nsigmafit*s2_sigma);
    fS2Source->Draw("same");
    s2_mean = fS2Source->GetParameter(1);
    s2_sigma = fS2Source->GetParameter(2);

    if(autoaxisscale_y) h2D_Source_py->SetAxisRange(s2_mean-5.*s2_sigma, s2_mean+8.5*s2_sigma, "X");

    gPad->Update();
    TPaveStats* st2 = (TPaveStats*)h2D_Source_py->FindObject("stats");
    st2->SetX1NDC(0.52);st2->SetX2NDC(0.9);
    st2->SetY1NDC(0.6);st2->SetY2NDC(0.9);
    st2->Draw();


  TLegend* leg2 = new TLegend(0.65, 0.25, 0.9, 0.4);
  leg2->SetBorderSize(0);
  leg2->SetTextSize(0.03);
  leg2->SetFillColor(10);
  if(!SubtractBG && useBG){
      TH1* h2D_Bg_py = GetProjection( h2D_Bg, xProjmin, xProjmax);
      h2D_Bg_py->SetLineColor(kRed);
      if(h2D_Bg_py) h2D_Bg_py->Draw("same");
      leg2->AddEntry(h2D_Source_py, "Source+Ar39", "L");
      leg2->AddEntry(h2D_Bg_py, "Ar39", "L");
      leg2->AddEntry(fS2Source, "Fit func.", "L");
      leg2->Draw();
      if(!drawInOnePad) theApp->Run();
  }
//  gPad->Update();
//  TPaveStats *st2 = (TPaveStats*)h2D_Source_py->FindObject("stats");
//  st2->SetX1NDC(0.52);st2->SetX2NDC(0.9);
//  st2->SetY1NDC(0.6);st2->SetY2NDC(0.9);
//  st2->Draw();



  (drawInOnePad)? c1->cd(4): c1->cd();

  TLatex L;
  L.SetNDC();
  L.SetTextSize(0.04);
  Double_t xposi = 0.01;
//  L.DrawLatex(xposi, 0.75, Form("Live Time(Source events) [s]: %4f", livetime_Source));
  L.DrawLatex(xposi, 0.7, Form("# of Event (Source events): %d", Nevents_Source));
  L.DrawLatex(xposi, 0.6, Form("S1 mean: %f #pm %f", s1_mean, s1_meanErr));
  L.DrawLatex(xposi, 0.55, Form("S2 mean: %f #pm %f", fS2Source->GetParameter(1), fS2Source->GetParError(1)));

  if(!useBG || !SubtractBG) L.DrawLatex(xposi, 0.2 , Form("Bg fraction under S1 Source peak: %f", bgFrac));
  L.DrawLatex(xposi, 0.15, Form("LY(S1_{mean}/41.5keV): %f #pm %f", s1_mean/41.5, s1_meanErr/41.5));




  theApp->Run();

/*  ofstream out;
  TString outputFilename = Form("S2vsS1_Edrift%sVOvercm.txt", Edrift.Data());
  while (true) {
    cout << "Storing output text to " << outputFilename.Data() << endl;
    out.open(outputFilename.Data());
    if (out.is_open()) break;
    cout << "Unable to open file " << outputFilename.Data() << endl;
    break;
  }
  out <<"#E_drift[V/cm]  S1_Mean[npe]  S1_MeanErr  S1_Sigma  S1_SigmaErr  S2_Mean_Fit  S2_Mean_FitErr  S2_Sigma_Fit  S2_Sigma_FitErr  S2_Mean_Hist  S2_Mean_HistErr  S2_Peak_Fit  S2_Peak_FitErr"<<endl;
  out << Edrift.Data() << " "
  << s1_mean << " "
  << s1_meanErr << " "
  << s1_sigma << " "
  << s1_sigmaErr << " "
  << MeanX << " "
  << 0.    << " "
  << fS2SourceX->Variance(fitS2Min, fitS2Max) << " "
  << 0.    << " "
  << meanFromHist << " "
  << meanFromHistErr << " "
  << peak << " "
  << peakErr << " "
  << endl;

  TString xtitle = h2D_Source_px->GetXaxis()->GetTitle();
  RemoveWhitespace(xtitle);
  TString ytitle = h2D_Source_py->GetXaxis()->GetTitle();
  RemoveWhitespace(ytitle);
  out <<"X_axis:  "<<xtitle.Data()<<"  Y_axis:  "<<ytitle.Data()<<endl;
*/
  if(drawInOnePad) c1->Print(Form("SourcePeak_Edrift%sVOvercm%s.pdf",Edrift.Data(), histoName.Data()));
//  c1->cd(1)->Print(Form("SourcePeak_Edrift%sVOvercm_2d.pdf",Edrift.Data()));

  TString outFileName;
  outFileName=Form("DS_Edrift%sVOvercm_%s.root",Edrift.Data(), histoName.Data());
  TFile* f = new TFile(outFileName.Data(), "RECREATE");
  f->cd();
//  h2D_Bg->SetName("h2D_Bg");
//  h2D_Bg->SetTitle("total_f90 vs total_s1 ^{39}Ar");
  h2D_Bg->Write();
//  h2D_Source_scaled->SetName("hF90VsTotalS1_SourceAr");
//  h2D_Source_scaled->SetTitle("total_f90 vs total_s1 ^{39}Ar+^{83m}Source");
  h2D_Source_scaled->Write();
//  h2D_Source->SetName("h2D_Source");
//  h2D_Source->SetTitle("total_f90 vs total_s1 ^{83m}Source");
  h2D_Source->Write();
  c1->SetName("FitResult");
  c1->SetTitle("Fit Result");
  if(drawInOnePad) c1->Write();

  h2D_Source_px->Write();
  h2D_Source_py->Write();

  h2D_Bg->ProjectionX()->Write();

  f->Write();
  f->Close();

  delete c1;
}

void RemoveWhitespace(TString &s){
  TObjArray *tokens = s.Tokenize(" ");
  Int_t n = tokens->GetEntries();
  TString tmp=((TObjString*) tokens->At(0))->GetString();
  for (Int_t i = 1; i < n; i++) {
      tmp+="_";
      tmp += ((TObjString*) tokens->At(i))->GetString();
  }
  cout<<"output string: "<<tmp.Data()<<endl;
  s=tmp;
  return;

}

Double_t GetFractionOfBG(TF1 *f, Double_t xmin, Double_t xmax){
  Double_t par[10];
  f->GetParameters(par);

  TF1 *fS1Source = new TF1("fS1Source",g2,f->GetXmin(),f->GetXmax(), 3);
  fS1Source->SetParameters(par);
  TF1 *fS1Bg = new TF1("fS1Bg","pol1",f->GetXmin(),f->GetXmax());
  fS1Bg->SetParameters(&par[3]);

  Double_t SourceInt = fS1Source->Integral(xmin, xmax);
  Double_t BgInt = fS1Bg->Integral(xmin, xmax);

  cout<<"Total entries in the signal region: "<<SourceInt+BgInt<<" Source: "<<SourceInt<<" Bg: "<<BgInt<<endl;
  cout<<"Source fraction in the signal region: "<<SourceInt/(SourceInt+BgInt)<<endl;
  return BgInt/(SourceInt+BgInt);
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


//-----------------------------------------------------------------------------------------------------
Double_t g2(Double_t *x, Double_t *par) {
        if(par[2]<=0.) return -99999.;
   Double_t r1 = Double_t((x[0]-par[1])/par[2]);
   return par[0]/TMath::Sqrt(TMath::TwoPi())/par[2]*TMath::Exp(-0.5*r1*r1);
}

//-----------------------------------------------------------------------------------------------------
Double_t g2Pol1(Double_t *x, Double_t *par) {
   Double_t pol1 = par[3] + x[0]*par[4];
   return g2(x,par) + pol1;
}

//-----------------------------------------------------------------------------------------------------
Double_t lognormal(Double_t *x, Double_t *par) {
  //LogNormal(Double_t x, Double_t sigma, Double_t theta = 0, Double_t m = 1)
  Double_t sigma = par[2]; //sigma is the shape parameter
  Double_t theta = par[1]; //theta is the location parameter
  Double_t m = par[0]; //m is the scale parameter
//  Double_t m = par[3]; //m is the scale parameter


  if(x[0]<theta) return 0.;
  if(m<0. || sigma<0.) return 0.;
//  return TMath::LogNormal(x[0], sigma, theta, m);

  return par[3]*TMath::Exp(TMath::Power(TMath::Log((x[0]-theta)/m), 2.) / (2.*sigma*sigma))/(x[0]-theta)/sigma/TMath::Sqrt(2.*TMath::Pi());
}

//-----------------------------------------------------------------------------------------------------
Double_t myfunc(Double_t *x, Double_t *par) {
  //LogNormal(Double_t x, Double_t sigma, Double_t theta = 0, Double_t m = 1)
//  Double_t amp = par[0]; //m is the scale parameter
//  Double_t shift = par[1]; //sigma is the shape parameter
//  Double_t m = par[3]; //m is the scale parameter


//  if(x[0]<theta) return 0.;
//  if(m<0. || sigma<0.) return 0.;
////  return TMath::LogNormal(x[0], sigma, theta, m);

//  return amp*TMath::Exp(par[2]*TMath::Power((x[0]-shift), par[3]));
//  return par[0]*TMath::Exp(x[0]/par[2])-par[3]*TMath::Exp(-(x[0]-par[1])/par[4]);

  if(x[0] > -(par[0]+par[3])/par[4]) return 0.;
  Double_t val = par[0]+TMath::TanH(par[2]*(x[0]-par[1]))*(par[3]+par[4]*x[0]);

  return (val>0.)? val : 0.;

}

Double_t myfuncX(Double_t *x, Double_t *par) {

  return x[0]*myfunc(x, par);
}

Double_t GetPeakError2(Double_t peak){
        Int_t npar =funcS2Source->GetNpar();
        Double_t *matr=(TVirtualFitter::GetFitter())->GetCovarianceMatrix();
        if(!matr) {
            cout<<"Covariant matrix is not avilable!! Please ckeck fit results."<<endl;
            return 0.;
        }
        // loop on the parameter and calculate the errors
        Double_t *ig = new Double_t[npar];
        Double_t *sum_vector = new Double_t[npar];

        for (int i=0; i < npar; ++i) {
                // check that parameter error is not zero - otherwise skip it
                // should check the limits
                double deriv = 0;
                deriv = dMeandPar(i, funcS2Source, peak);
                ig[i] = deriv;
        }

        Double_t err2(0.);
        for (Int_t irow=0; irow<npar; irow++) {
                sum_vector[irow]=0;
                for (Int_t icol=0; icol<npar; icol++){
                      cout<<matr[irow*npar+icol]<<" "<<ig[icol]<<" ";
                        sum_vector[irow]+=matr[irow*npar+icol]*ig[icol];
                }
              cout<<endl;
                err2 += ig[irow] * sum_vector[irow];
        }

//      double err2 = covMatrix.Similarity(ig);
        delete ig; delete sum_vector;

        return err2;
}

//-----------------------------------------------------------------------------------------------------
Double_t dMeandPar(int par, TF1* f, Double_t peak){
//  cout<<"here!!"<<endl;

  Double_t Par = f->GetParameter(par);
  Double_t dPar = f->GetParError(par);
  f->SetParameter(par, Par+ dPar);
  Double_t dPeak = f->GetMaximumX()-peak;
  f->SetParameter(par, Par);
  return dPeak/dPar;

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
	if ( theApp->Argc() == 4 ) {
		std::cout << "\n==========> ExtractBaSpectum!!! <=============" << std::endl;
		std::cout << "==> Application start." << std::endl;
	        ExtractBaSpectum(theApp->Argv(1), theApp->Argv(2), theApp->Argv(3));

	} else if(theApp->Argc() == 3) {
            std::cout << "\n==========> ExtractBaSpectum without Bg Runs!!! <=============" << std::endl;
            std::cout << "==> Application start." << std::endl;
            ExtractBaSpectum(theApp->Argv(1), theApp->Argv(2), "");

	} else {
		std::cout << "Usage:" <<theApp->Argc()<< std::endl;
		std::cout << "./ExtractBaSpectum Edrift fInName fInName_Bg" << std::endl;
		return 0;
	}


//        std::cout << "==> Application start." << std::endl;
//	ExtractBaSpectum(theApp->Argv(1), theApp->Argv(2));
	std::cout << "==> Application finished." << std::endl;
//	theApp->Run();

	return 0;
}
#endif /* __CINT __ */
