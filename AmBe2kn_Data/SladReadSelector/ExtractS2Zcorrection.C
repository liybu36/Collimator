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

Int_t fEdrift(0);
Double_t fDriftTimeMax(0.), ft_drift_fitmin(0.), ft_drift_fitmax(0.);

void SetDriftFieldConf(Int_t Edrift)
{
  fEdrift = Edrift;
  Double_t DriftTimeMax_default = 373.3;
  Double_t t_drift_fitmin_default = 80.;
  Double_t t_drift_fitmax_default = 330.;
  if(Edrift == 0){
      fDriftTimeMax = 10.; // usually acquisition window is 25 us
  } else if (Edrift == 50){
      fDriftTimeMax = 1250 ;
      ft_drift_fitmin = t_drift_fitmin_default*fDriftTimeMax/DriftTimeMax_default;
      ft_drift_fitmax = t_drift_fitmax_default*fDriftTimeMax/DriftTimeMax_default;
  } else if (Edrift == 100){
      fDriftTimeMax =  665.;
      ft_drift_fitmin = t_drift_fitmin_default*fDriftTimeMax/DriftTimeMax_default;
      ft_drift_fitmax = t_drift_fitmax_default*fDriftTimeMax/DriftTimeMax_default;
  } else if (Edrift == 150){
      fDriftTimeMax = 469.;
      ft_drift_fitmin = t_drift_fitmin_default*fDriftTimeMax/DriftTimeMax_default;
      ft_drift_fitmax = t_drift_fitmax_default*fDriftTimeMax/DriftTimeMax_default;
  } else if (Edrift == 200){
      fDriftTimeMax = 373.3;
      ft_drift_fitmin = t_drift_fitmin_default;
      ft_drift_fitmax = t_drift_fitmax_default;
  } else {
      cout<<"Unknown Drift field value!! Please check input value. E_drift: "<<Edrift<<endl;
  }

}

void electron_lifetime(){
  TCanvas *c1 = new TCanvas("c1","c1",1000, 750);
  c1->cd();
  c1->Draw();
  c1->SetGrid();


  const Int_t nfields =4;
  Int_t fields[] = {50, 100, 150, 200}; // fields in V/cm
  Double_t elifetime[] = {5226, 4344, 5863, 5574}; // fields in V/cm
  Double_t elifetime_err[] = {452.6, 383.8, 1882.1, 978.7}; // fields in V/cm


  TH1D *helectronlife = new TH1D("helectronlife", "Electron lifetime vs Drift field; field [V/cm]; electron lifetime [#mus]", nfields, -0.5, nfields-0.5);
  for (Int_t i=0; i<nfields; i++) helectronlife->GetXaxis()->SetBinLabel(i+1, Form("%d", fields[i]));

  for (Int_t i=0; i<nfields; i++){
      helectronlife->SetBinContent(i+1, elifetime[i]);
      helectronlife->SetBinError(i+1, elifetime_err[i]);
  }

  helectronlife->Draw();

  TF1 *f1 = new TF1("f1", "pol0", -0.5, nfields-0.5);
  helectronlife->Fit("f1", "REM");
  c1->Update();
  theApp->Run();
  c1->Print("Electron_lifetime_fit.pdf");



}


void ExtractS2Xcorrection(TString fInNameSource, TString fInNameBGRun) {
  SetColors();
  gStyle->SetOptFit(1);
//  gStyle->SetOptStat(1);
  gStyle->SetOptStat(1100);
  gStyle->SetOptTitle(1);
  TGaxis::SetMaxDigits(3);

  Bool_t drawInOnePad(false);//true);
  Bool_t useBG = true;//false;
  Bool_t SubtractBG=true; //if you want to draw Source & Bg together, set this to false
  Bool_t fit = true;

    //-----Setting for s2 vs tdrift for Kr-------------------
  TString histoName("ht_drift_S2XYCorr_Kr");// histogram name you want to get;
//  TString histoName("ht_drift_ScaledS2XYCorroverS1_Kr");// histogram name you want to get;
        //Normalize by integration
//    #define NORMINTEGRAL
//  Int_t rebinX(10), rebinY(10);
  Int_t rebinX(1), rebinY(1);
        Double_t nomMin(400.), normMax(550.);
        Double_t xmin(0.), xmax(800.), ymin(0.), ymax(20.e+3);
        Double_t x_int(2.), x_mean(290.), x_sigma(19), x_fitmin(150), x_fitmax(450); // fit initial values
        Double_t y_int(160.),y_mean(8.5e+3), y_sigma(1.4e+3), y_fitmin(2.e+3), y_fitmax(16.e+3);
        Bool_t autoaxisscale_x(true), autoaxisscale_y(true);

#if 1
      std::cout << "Which field??? Type 'q' to quit or 'number' to set field." << std::endl;
      std::string input = "";
      std::cin>>input;
      if(input=="q") return;
      TString field = input;
      SetDriftFieldConf(field.Atoi());
#endif


        //------Setting for s1 vs s2 for Kr------------------

        TCanvas *c1 = new TCanvas("c1","c1",1000, 750);
        c1->cd();
      //  c1->Draw();
        if(drawInOnePad)c1->Divide(2,3);
        (drawInOnePad)? c1->cd(1): c1->cd();

        TLegend* leg = new TLegend(0.65, 0.65, 0.9, 0.85);
          leg->SetBorderSize(0);
          leg->SetTextSize(0.03);
          leg->SetFillColor(0);
        //  leg->SetHeader("");

  cout<<"Input Source file: "<<fInNameSource.Data()<<endl;


  TFile *hfileSource = new TFile(fInNameSource.Data());
  if(!hfileSource) cout<<"file: "<<fInNameSource.Data()<<" is not found."<<endl;

  TH2D* h2D_Source = (TH2D*) hfileSource->Get(histoName.Data());
  if(!h2D_Source) cout<<"hist: "<<histoName.Data()<<" is not found."<<endl;
  h2D_Source->SetDirectory(0);
  h2D_Source->Sumw2();
  h2D_Source->SetName(Form("%s_Source", h2D_Source->GetName()));
  h2D_Source = (TProfile2D*)h2D_Source->Rebin2D(rebinX,rebinY);

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

  TH2D* h2D_Bg = NULL;
  Double_t livetime_Bg(1.);
  if(useBG){
      TFile *hfileBg = new TFile(fInNameBGRun.Data());
      if(!hfileBg) cout<<"file: "<<fInNameBGRun.Data()<<" is not found."<<endl;
      h2D_Bg = (TProfile2D*) hfileBg->Get(histoName.Data());
      if(!h2D_Bg) cout<<"hist: "<<histoName.Data()<<" is not found."<<endl;
      h2D_Bg->SetDirectory(0);
      h2D_Bg->Sumw2();
      h2D_Bg->SetName(Form("%s_Bg", h2D_Bg->GetName()));
      h2D_Bg = (TProfile2D*)h2D_Bg->Rebin2D(rebinX,rebinY);

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
  TH2D* htmp = (TH2D*)h2D_Source->Clone("colz");
  htmp->SetStats(0);
//  htmp->SetTitleOffset(1.2, "Y");
  htmp->SetTitle(Form("Kr S2_{XY corr not z} vs t_drift before BG subtraction @ %s V/cm", field.Data()));

  htmp->Draw("colz");

  if(!drawInOnePad) {
      c1->Update();
      theApp->Run();
      c1->Print(Form("Source_BeforeBGsubt_%s_%sVcm.pdf", histoName.Data(), field.Data()));
  }

  TString hname_fit;
  //  Int_t col_mean(colors.at(4)), col_sig(colors.at(2)), col_prof(colors.at(0)), col_fit(colors.at(6));
    Int_t col_mean(kBlack), col_sig(kBlack), col_prof(kGreen), col_fit(kRed);
    TF1 *fexpo = new TF1("fexpo", "[0]*exp(-x/[1])", ft_drift_fitmin, ft_drift_fitmax);
    fexpo->SetParNames("constant","e_lifetime");
    fexpo->SetParameters(1.,2.e+04);
    fexpo->SetLineColor(col_fit);

  if(useBG){
//  h2D_Bg->Scale(livetime_Source/livetime_Bg);
  (drawInOnePad)? c1->cd(2)->SetGrid(): c1->SetGrid();
  htmp = (TH2D*)h2D_Bg->Clone("colz");
  htmp->SetStats(0);
  htmp->SetTitle(Form("^{39}Ar S2_{XY corr not z} vs t_drift around Kr s1 peak @ %s V/cm", field.Data()));
//  htmp->SetTitleOffset(1.2, "Y");
  htmp->Draw("colz");

  h2D_Bg->FitSlicesY();
  hname_fit = Form("%s_1", h2D_Bg->GetName());
  TH1D *h2D_Bg_mean = (TH1D*)gDirectory->Get(hname_fit.Data());
  if(!h2D_Bg_mean) cout<<"hist: "<<hname_fit.Data()<<" is not found."<<endl;
  hname_fit = Form("%s_2", h2D_Bg->GetName());
  TH1D *h2D_Bg_sigma = (TH1D*)gDirectory->Get(hname_fit.Data());
  if(!h2D_Bg_sigma) cout<<"hist: "<<hname_fit.Data()<<" is not found."<<endl;

  TProfile *h2D_Bg_pfx = h2D_Bg->ProfileX();//"s");

  h2D_Bg_mean->Draw("sames");
  h2D_Bg_mean->SetLineColor(col_mean);
  h2D_Bg_mean->SetMarkerStyle(20);
  h2D_Bg_mean->SetMarkerColor(col_mean);

  h2D_Bg_sigma->Draw("same");
  h2D_Bg_sigma->SetLineColor(col_sig);
  h2D_Bg_sigma->SetMarkerStyle(21);
  h2D_Bg_sigma->SetMarkerColor(col_sig);

  h2D_Bg_pfx->Draw("same");
  h2D_Bg_pfx->SetLineColor(col_prof);
  h2D_Bg_pfx->SetMarkerStyle(22);
  h2D_Bg_pfx->SetMarkerColor(col_prof);

  //  TF1 *fexpo = new TF1("fexpo", "expo", ft_drift_fitmin, ft_drift_fitmax);
  h2D_Bg_mean->Fit("fexpo", "EMR+");
  fexpo->Draw("same");

  leg->AddEntry(h2D_Bg_pfx, "Mean from profile", "LP");
  leg->AddEntry(h2D_Bg_mean, "Gauss. fit mean", "LP");
  leg->AddEntry(h2D_Bg_sigma, "Gauss. fit sigma", "LP");
  leg->AddEntry(fexpo, "Fit to fit mean", "LP");
  leg->Draw();

  if(!drawInOnePad) {
      c1->Update();
      theApp->Run();
      c1->Print(Form("BG_%s_%sVcm.pdf", histoName.Data(), field.Data()));
  }

  }

  (drawInOnePad)? c1->cd(3)->SetGrid(): c1->SetGrid();
  leg->Clear();

  TH2* hErr(NULL), *hEntry(NULL);
  if(SubtractBG && useBG){
    //Subtract normalized BG.
      h2D_Source->Add(h2D_Bg, -livetime_Source/livetime_Bg);
      h2D_Source->SetTitle(Form("Kr S2_{XY corr not z} vs t_drift after BG subtraction around Kr s1 peak @ %s V/cm", field.Data()));
      cout<<"Scaled: "<<livetime_Source/livetime_Bg<<endl;
  }

  htmp = (TH2D*)h2D_Source->Clone("colz");
  htmp->SetStats(0);
//  htmp->SetTitleOffset(1.2, "Y");
  htmp->Draw("colz");

  h2D_Source->FitSlicesY();
  hname_fit = Form("%s_1", h2D_Source->GetName());
  TH1D *h2D_Source_mean = (TH1D*)gDirectory->Get(hname_fit.Data());
  if(!h2D_Source_mean) cout<<"hist: "<<hname_fit.Data()<<" is not found."<<endl;
  hname_fit = Form("%s_2", h2D_Source->GetName());
  TH1D *h2D_Source_sigma = (TH1D*)gDirectory->Get(hname_fit.Data());
  if(!h2D_Source_sigma) cout<<"hist: "<<hname_fit.Data()<<" is not found."<<endl;

  TProfile *h2D_Source_pfx = h2D_Source->ProfileX();//"s");

  h2D_Source_mean->Draw("sames");
  h2D_Source_mean->SetLineColor(col_mean);
  h2D_Source_mean->SetMarkerStyle(20);
  h2D_Source_mean->SetMarkerColor(col_mean);

  h2D_Source_sigma->Draw("same");
  h2D_Source_sigma->SetLineColor(col_sig);
  h2D_Source_sigma->SetMarkerStyle(21);
  h2D_Source_sigma->SetMarkerColor(col_sig);

  h2D_Source_pfx->Draw("same");
  h2D_Source_pfx->SetLineColor(col_prof);
  h2D_Source_pfx->SetMarkerStyle(22);
  h2D_Source_pfx->SetMarkerColor(col_prof);

////  TF1 *fexpo = new TF1("fexpo", "expo", ft_drift_fitmin, ft_drift_fitmax);
//  TF1 *fexpo = new TF1("fexpo", "[0]*exp(-x/[1])", ft_drift_fitmin, ft_drift_fitmax);
//  fexpo->SetParNames("constant","e_lifetime");
//  fexpo->SetParameters(1.,2.e+04);
  h2D_Source_mean->Fit("fexpo", "EMR+");
  fexpo->Draw("same");

  leg->AddEntry(h2D_Source_pfx, "Mean from profile", "LP");
  leg->AddEntry(h2D_Source_mean, "Gauss. fit mean", "LP");
  leg->AddEntry(h2D_Source_sigma, "Gauss. fit sigma", "LP");
  leg->AddEntry(fexpo, "Fit to fit mean", "LP");
  leg->Draw();

  if(!drawInOnePad) {
      c1->Update();
      theApp->Run();
      c1->Print(Form("Source_AfterBGsubt_%s_%sVcm.pdf", histoName.Data(), field.Data()));
  }


/*



  (drawInOnePad)? c1->cd(5)->SetGrid(): c1->SetGrid();


  Float_t S1corr_mean_Kr_peak = 290.6;
  Float_t S1corr_sigma_Kr_peak = 18.45;
  Float_t nsigma = 3.;
  Float_t s1_min = S1corr_mean_Kr_peak-nsigma*S1corr_sigma_Kr_peak;
  Float_t s1_max = S1corr_mean_Kr_peak+nsigma*S1corr_sigma_Kr_peak;

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

*/

/*  TString outFileName;
  outFileName=Form("S2vsTDrift_Kr_%s.root", histoName.Data());
  TFile* f = new TFile(outFileName.Data(), "RECREATE");
  f->cd();
  h2D_Bg->Write();
  h2D_Source->Write();
  h2D_Bg->Write();
  htemp->Write();
  h2D_Source_py->Write();
  c1->SetName("S2CorrMap");
  c1->SetTitle("S2 correction map");
  if(drawInOnePad) c1->Write();

  f->Write();
  f->Close();
*/
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


void SetColors(){
  colors.push_back(TColor::GetColor("#7F7F7E")); // Gray
  colors.push_back(TColor::GetColor("#44A5CB")); // Light Blue
  colors.push_back(TColor::GetColor("#D45D87")); // Red
  colors.push_back(TColor::GetColor("#F9DB57")); // Yellow
  colors.push_back(TColor::GetColor("#9D73BB")); // Purple
  colors.push_back(TColor::GetColor("#009F8C")); // Green
  colors.push_back(TColor::GetColor("#FF8600")); // Orange
  colors.push_back(TColor::GetColor("#E8300F")); // Orange Red
  colors.push_back(TColor::GetColor("#8B90BE")); // Light Purple
  colors.push_back(TColor::GetColor("#F3C0AB"));
}

#ifndef __CINT__
int main(int argc, char **argv) {
	theApp = new TRint("App", &argc, argv, NULL, 0);
        theApp->Connect("KeyPressed(Int_t)","TSystem",gSystem,"ExitLoop()");
	if (theApp->Argc() == 3) {
            std::cout << "\n==========> ExtractS2Xcorrection without Bg Runs!!! <=============" << std::endl;
            std::cout << "==> Application start." << std::endl;
            ExtractS2Xcorrection(theApp->Argv(1), theApp->Argv(2));

	} else if (theApp->Argc() == 1) {
	    electron_lifetime();
	} else {
		std::cout << "Usage:" <<theApp->Argc()<< std::endl;
		std::cout << "./ExtractS2Xcorrection fInName fInName_Bg" << std::endl;
		return 0;
	}


//        std::cout << "==> Application start." << std::endl;
//	ExtractS2Xcorrection(theApp->Argv(1), theApp->Argv(2));
	std::cout << "==> Application finished." << std::endl;
//	theApp->Run();

	return 0;
}
#endif /* __CINT __ */
