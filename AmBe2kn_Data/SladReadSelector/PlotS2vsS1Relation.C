#include <iostream>
#include <fstream>
#include "TCanvas.h"
#include "TFile.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TApplication.h"
#include "TRint.h"
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include "TColor.h"
#include <TGaxis.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <Rtypes.h>
#include "TPaveStats.h"
#include "TObjString.h"

using namespace std;

bool readValuesFromFile(TString inputFName, vector<vector<TString> > &allInfos);
Double_t Pol1(Double_t *x, Double_t *par);

#ifndef __CINT__
  int main(int argc, char **argv);
#endif

TRint* theApp;
vector<int> colors;
void SetColors();

void PlotS2vsS1Relation() {
  gStyle->SetOptFit(1);
  TGaxis::SetMaxDigits(3);
  SetColors();

  TCanvas *c1 = new TCanvas("c1","c1",1000, 750);
  c1->cd();
  c1->Draw();
  c1->SetGrid();

/*  const Int_t nFields = 4;
  TString EDrift[nFields] = {"200", "150", "100", "0"}; // in V/cm

  //Defalt Values
  // E_drift = 200, 150, 100, 0
//  Double_t S1[nFields]  = {398.3, 406.3, 414.6}; // in npe
//  Double_t S2FromHist[nFields]  = {7024.8, 6512.0, 6315.4};
  Double_t S1[nFields]            = {293.1, 306.9, 306.8, 335.8}; // in npe
  Double_t S1Err[nFields]         = {0.5, 0.4, 0.6, 0.5};
  Double_t S2FromHist[nFields]    = {5113.3, 5311.1, 4642.3, 0};
  Double_t S2FromHistErr[nFields] = {48.8, 29.4, 23.8, 0};
  Double_t S2Func[nFields]  = {5120.2, 5260.4, 4622.1, 0};
  Double_t S2Peak[nFields]  = {2855.3, 3041.6, 1793.8, 0};
  Double_t S2PeakErr[nFields]  = {0., 0., 0., 0};
*/
/*
  const Int_t nFields = 2;
  TString EDrift[nFields] = {"200", "100"}; // in V/cm
  Double_t S1[nFields]            = {290.57, 301.66}; // in npe
  Double_t S1Err[nFields]         = {0.13, 0.18}; //error on s1
//  Double_t S1Err[nFields]         = {18.45, 18.61}; //sigma
  Double_t S2FromHist[nFields]  = {8524.36, 5964.25};
  Double_t S2FromHistErr[nFields]  = {10.08, 9.91}; //error on s2
//  Double_t S2FromHistErr[nFields]  = {1472, 932.8}; //sigma of s2

  const Int_t nFields_old = 3;
  Double_t S1_old[]            = {289.78, 296.89, 302.33}; // in npe
  Double_t S1_oldErr[]         = {0.56, 0.50, 0.71}; //error on s1
  Double_t S2Peak[]            = {8428.67, 7227.37, 5770.15};
  Double_t S2PeakErr[]         = {40.09, 40.69, 33.04};


  const Int_t npoints_all = nFields_old+nFields;
  Double_t S1_all[]            = {290.57, 301.66, 289.78, 296.89, 302.33}; // in npe
  Double_t S1_allErr[]         = {0.13, 0.18, 0.56, 0.50, 0.71}; //error on s1
  Double_t S2all[]             = {8524.36, 5964.25, 8428.67, 7227.37, 5770.15};
  Double_t S2allErr[]          = {10.08, 9.91, 40.09, 40.69, 33.04};
*/
//  const Int_t nFields = 3;
//  TString EDrift[nFields] = {"200", "100", "0"}; // in V/cm
//  Double_t S1[nFields]            = {292.9, 303.9, 326.7}; // in npe
//  Double_t S1Err[nFields]         = {0.1, 0.2, 0.2}; //error on s1
////  Double_t S1Err[nFields]         = {18.45, 18.61}; //sigma
//  Double_t S2FromHist[nFields]  = {8531, 5957, 0};
//  Double_t S2FromHistErr[nFields]  = {10.3, 20.4, 0}; //error on s2
////  Double_t S2FromHistErr[nFields]  = {1472, 932.8}; //sigma of s2

  const Int_t nFields = 3;
  TString EDrift[nFields] = {"200", "100", "0"}; // in V/cm
  Double_t S1[nFields]            = {290.7, 301.7, 326.7}; // in npe
  Double_t S1Err[nFields]         = {0.1, 0.2, 0.2}; //error on s1
//  Double_t S1Err[nFields]         = {18.45, 18.61}; //sigma
  Double_t S2FromHist[nFields]  = {8171., 5840., 0};
  Double_t S2FromHistErr[nFields]  = {9.0, 9.6, 0}; //error on s2
//  Double_t S2FromHistErr[nFields]  = {1472, 932.8}; //sigma of s2

  const Int_t nFields_old = 4;
  Double_t S1_old[]            = {291.3, 299.3, 304.2, 329.5}; // in npe
  Double_t S1_oldErr[]         = {0.5, 0.5, 0.7, 0.9}; //error on s1
  Double_t S2Peak[]            = {8407, 7260, 5777, 0};
  Double_t S2PeakErr[]         = {48.1, 43.3, 48.6, 0};


  const Int_t npoints_all = nFields_old+nFields;
  Double_t S1_all[]            = {290.57, 301.66, 289.78, 296.89, 302.33}; // in npe
  Double_t S1_allErr[]         = {0.13, 0.18, 0.56, 0.50, 0.71}; //error on s1
  Double_t S2all[]             = {8524.36, 5964.25, 8428.67, 7227.37, 5770.15};
  Double_t S2allErr[]          = {10.08, 9.91, 40.09, 40.69, 33.04};

  const Int_t nFields_Feb2015 = 4;
  TString EDrift_Feb2015[nFields_Feb2015] = {"200", "100", "50", "0"}; // in V/cm
  Double_t S1_Feb2015[nFields_Feb2015]    = {288.8, 300.6, 308.1, 321.6}; // in npe // TBA corrected not usual z correction
  Double_t S1Err_Feb2015[nFields_Feb2015] = {0.6, 0.8, 0.6, 0.7}; //error on s1
  Double_t S2_Feb2015[nFields_Feb2015]    = {8586, 5724, 4191, 0};
  Double_t S2Err_Feb2015[nFields_Feb2015] = {89.6, 78.0, 35.1, 0}; //error on s2
//  Double_t S2FromHistErr[nFields]  = {1472, 932.8}; //sigma of s2

  const Int_t nFields_Mar2015 = 4;
  TString EDrift_Mar2015[nFields_Mar2015] = {"200", "150", "100", "50"}; // in V/cm
  Double_t S1_Mar2015[nFields_Mar2015]    = {289.2, 296.9, 302.4, 311.4}; // in npe // TBA corrected not usual z correction
  Double_t S1Err_Mar2015[nFields_Mar2015] = {0.18, 0.4, 0.2, 0.3}; //error on s1
  Double_t S2_Mar2015[nFields_Mar2015]    = {7916., 6820., 5520., 4128.};
  Double_t S2Err_Mar2015[nFields_Mar2015] = {14.5, 30.7, 14.1, 13.3}; //error on s2


//
//  for (Int_t i=0; i<npoints_all; i++){
//      S1_allErr[i]=TMath::Sqrt(S1_allErr[i]*S1_allErr[i]+ 0.5*0.5*1.e-4*S1_all[i]*S1_all[i]);
//      S2allErr[i]=TMath::Sqrt(S2allErr[i]*S2allErr[i]+ 0.5*0.5*1.e-4*S2all[i]*S2all[i]);
//  }

  TString xtitle("total_s1_corr [PE]"), ytitle("total_s2_xyz_corr [PE]");
/*
  //Read from File created with S2vsS1Fitter.C
  for(Int_t iField=0; iField<nFields; iField++){
      TString outputFilename = Form("S2vsS1_Edrift%sVOvercm.txt", EDrift[iField].Data());
      cout<<"\nOpen file: "<<outputFilename.Data()<<endl;
      vector<vector<TString> > allInfos;
      bool isFileOK = readValuesFromFile(outputFilename, allInfos);
      const int nInfos = (const int)allInfos.size(); // This has to be 2. One is for labels and the other is for values.
      if (nInfos!=3)cout<<"+++++++ Something Wrong!! Please check input file: "<<outputFilename.Data()<<endl;
      //loop over entries and create graphs      Int_t iInfo = 0;
      vector<TString> const& labels = allInfos.at(0);
      vector<TString> const& values = allInfos.at(1);
      const int nEntries = (const int)values.size();
      if(nEntries!=labels.size()) cout<<"number of values does not match with number of labels!!"<<endl;
      for(int i = 0; i < nEntries; i++) {
          cout<<labels.at(i).Data()<<"   "<<values.at(i).Data()<<endl;
          if(labels.at(i).Contains("E_drift") && !values.at(i).EqualTo(EDrift[iField]))
            cout<<"E_drift field does not match!! Please check input file: "<<outputFilename.Data()<<endl;

          if(labels.at(i).EqualTo("S1_Mean[npe]")) S1[iField] = values.at(i).Atof();
          else if(labels.at(i).EqualTo("S1_MeanErr")) S1Err[iField] = values.at(i).Atof();
          else if(labels.at(i).EqualTo("S2_Mean_Fit")&&!EDrift[iField].EqualTo("0")) S2Func[iField] = values.at(i).Atof();
//          else if(labels.at(i).EqualTo("S2_Mean_FitErr")) S2FuncErr[iField] = values.at(i).Atof();
          else if(labels.at(i).EqualTo("S2_Mean_Hist")&&!EDrift[iField].EqualTo("0")) S2FromHist[iField] = values.at(i).Atof();
          else if(labels.at(i).EqualTo("S2_Mean_HistErr")&&!EDrift[iField].EqualTo("0")) S2FromHistErr[iField] = values.at(i).Atof();
          else if(labels.at(i).EqualTo("S2_Peak_Fit")&&!EDrift[iField].EqualTo("0")) S2Peak[iField] = values.at(i).Atof();
          else if(labels.at(i).EqualTo("S2_Peak_FitErr")&&!EDrift[iField].EqualTo("0")) S2PeakErr[iField] = values.at(i).Atof();
          else continue;
      }

      vector<TString> const& axistitle = allInfos.at(2);
      if(axistitle.at(0).EqualTo("X_axis:")) xtitle = axistitle.at(1);
      if(axistitle.at(2).EqualTo("Y_axis:")) ytitle = axistitle.at(3);
  }
*/
//  double_t min(0.), max(500.), ymin(0.), ymax(30.e+3);
  Double_t min(270.), max(360.), ymin(0.), ymax(9.e+3);
  TH1F *frame = (TH1F*)c1->DrawFrame(min, ymin, max, ymax);
  frame->SetXTitle(xtitle.Data());//"Fit Width [GeV/c^{2}]");
  frame->SetYTitle(ytitle.Data());//"Fit Mass [GeV/c^{2}]");
//  frame->SetTitleOffset(titleOffset, "X");
//  frame->SetTitleOffset(titleOffset+0.1, "Y");

  TLegend* leg = new TLegend(0.6, 0.5, 0.9, 0.7);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetFillColor(10);
//  leg->SetHeader("TOP PMTs");
//  leg-> SetNColumns(2);

  TGraphErrors *gr = new TGraphErrors(nFields,S1,S2FromHist, S1Err, S2FromHistErr);
  gr->SetName("grS1vsS2");
  gr->SetTitle(Form("%s vs %s", ytitle.Data(), xtitle.Data()));
  gr->GetXaxis()->SetTitle(xtitle.Data());
  gr->GetYaxis()->SetTitle(ytitle.Data());
  gr->SetLineColor(colors[1]);//kRed);
  gr->SetLineWidth(2);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1.3);
  gr->SetMarkerColor(gr->GetLineColor());//kRed);

  gr->Draw("p same");//ap");

  leg->AddEntry(gr, "Kr in Feb. 2014", "Lp");
  TF1 *fpol1 = new TF1("fpol1",Pol1, 280, 320, 2);
  fpol1->SetParameters(0.16, 28.);
  fpol1->SetParNames("#epsilon_{1}", "#epsilon_{2}");
  fpol1->SetLineColor(gr->GetLineColor());
  gr->Fit(fpol1, "EMR");
  fpol1->SetRange(min, max);
  fpol1->Draw("same");

  gPad->Update();
  TPaveStats *st = (TPaveStats*)gr->FindObject("stats");
  st->SetX1NDC(0.52);st->SetX2NDC(0.9);
  st->SetY1NDC(0.6);st->SetY2NDC(0.8);
  st->SetTextColor(gr->GetLineColor());
  st->Draw();


/*  TGraphErrors *grFunc = new TGraphErrors(npoints_all,S1_all,S2all, S1_allErr, S2allErr);
  grFunc->SetName("grFuncS1vsS2");
//  grFunc->SetLineColor(colors[1]);//kOrange);
//  grFunc->SetMarkerStyle(21);
  grFunc->SetMarkerSize(0);
//  grFunc->SetMarkerColor(colors[1]);//kOrange);
  grFunc->Draw("p same");
//  leg->AddEntry(grFunc, "S2 mean from Func.", "Lp");
*/
  TGraphErrors *grPeak = new TGraphErrors(nFields_old,S1_old,S2Peak, S1_oldErr, S2PeakErr);
  grPeak->SetName("grPeakS1vsS2");
  grPeak->SetLineColor(colors[2]);//kCyan);
  grPeak->SetMarkerStyle(22);
  grPeak->SetMarkerSize(1.3);
  grPeak->SetMarkerColor(grPeak->GetLineColor());//kCyan);
//  leg->AddEntry(grPeak, "Kr in Oct. 2013", "Lp");
//  grPeak->Draw("p same");

  TGraphErrors *grFeb2015 = new TGraphErrors(nFields_Feb2015,S1_Feb2015,S2_Feb2015, S1Err_Feb2015, S2Err_Feb2015);
  grFeb2015->SetName("grFeb2015");
  grFeb2015->SetLineColor(colors[3]);//kCyan);
  grFeb2015->SetMarkerStyle(22);
  grFeb2015->SetMarkerSize(1.3);
  grFeb2015->SetMarkerColor(grFeb2015->GetLineColor());//kCyan);
//  leg->AddEntry(grFeb2015, "Kr in Feb. 2015", "Lp");
//  grFeb2015->Draw("p same");

  TGraphErrors *grMar2015 = new TGraphErrors(nFields_Mar2015,S1_Mar2015,S2_Mar2015, S1Err_Mar2015, S2Err_Mar2015);
  grMar2015->SetName("grMar2015");
  grMar2015->SetLineColor(colors[2]);//kCyan);
  grMar2015->SetMarkerStyle(22);
  grMar2015->SetMarkerSize(1.3);
  grMar2015->SetMarkerColor(grMar2015->GetLineColor());//kCyan);
  leg->AddEntry(grMar2015, "Kr in Mar. 2015", "Lp");
  grMar2015->Draw("p sames");

/*
//  TF1 *fpol12 = new TF1("fpol12","pol1",S1[0],S1[nFields-1]);
//  TF1 *fpol12 = new TF1("fpol12",Pol1,S1[0],S1[nFields-1], 2);
  TF1 *fpol12 = new TF1("fpol12",Pol1, 280, 320, 2);
//  fpol12->SetParameters(S2Peak[0], (S2Peak[0]-S2Peak[nFields-1])/(S1[0]-S1[nFields-1]));
  fpol12->SetParameters(0.16, 28.);
//  fpol12->SetParameters(6.e-4, 2.e+2);
//  fpol12->SetParameters(6.1e+2, 0.122);
//  fpol12->SetParNames("W/g_{2}", "g_{2}/g_{1}");
  fpol12->SetParNames("#epsilon_{1}", "#epsilon_{2}");
  fpol12->SetLineColor(grPeak->GetLineColor());
  grPeak->Fit(fpol12, "EMR");
  fpol12->SetRange(min, max);
  fpol12->Draw("same");

  gPad->Update();
  TPaveStats *st2 = (TPaveStats*)grPeak->FindObject("stats");
  st2->SetX1NDC(0.52);st2->SetX2NDC(0.9);
  st2->SetY1NDC(0.6);st2->SetY2NDC(0.8);
  st2->SetTextColor(grPeak->GetLineColor());
  st2->Draw();


  TF1 *fpol14 = new TF1("fpol14",Pol1, 280, 320, 2);
  fpol14->SetParameters(0.16, 28.);
//  fpol14->SetParameters(6.e-4, 2.e+2);
//  fpol14->SetParNames("W/g_{2}", "g_{2}/g_{1}");
  fpol14->SetParNames("#epsilon_{1}", "#epsilon_{2}");
  fpol14->SetLineColor(grFeb2015->GetLineColor());
  fpol14->SetLineStyle(2);
  grFeb2015->Fit(fpol14, "EMR");
  fpol14->SetRange(min, max);
  fpol14->Draw("same");

  gPad->Update();
  TPaveStats *st4 = (TPaveStats*)grFeb2015->FindObject("stats");
  st4->SetX1NDC(0.52);st4->SetX2NDC(0.9);
  st4->SetY1NDC(0.6);st4->SetY2NDC(0.8);
  st4->SetTextColor(grFeb2015->GetLineColor());
  st4->Draw();
*/
  TF1 *fpol15 = new TF1("fpol15",Pol1, 280, 320, 2);
  fpol15->SetParameters(0.16, 28.);
//  fpol15->SetParameters(6.e-4, 2.e+2);
//  fpol15->SetParNames("W/g_{2}", "g_{2}/g_{1}");
  fpol15->SetParNames("#epsilon_{1}", "#epsilon_{2}");
  fpol15->SetLineColor(grMar2015->GetLineColor());
  fpol15->SetLineStyle(2);
  grMar2015->Fit(fpol15, "EMR0");
  fpol15->SetRange(min, max);
  fpol15->Draw("same");
  gPad->Update();
  TPaveStats *st5 = (TPaveStats*)grMar2015->FindObject("stats");
  st5->SetX1NDC(0.52);st5->SetX2NDC(0.9);
  st5->SetY1NDC(0.6);st5->SetY2NDC(0.8);
  st5->SetTextColor(grMar2015->GetLineColor());
  st5->Draw();


/*
  TF1 *fpol13 = new TF1("fpol13",Pol1, min, max, 2); // previous result
  fpol13->SetParameters(0.16, 28.);
//  fpol13->SetParameters(6.e-4, 2.e+2);
//  fpol13->SetParNames("W/g_{2}", "g_{2}/g_{1}");
  fpol13->SetParNames("#epsilon_{1}", "#epsilon_{2}");
  fpol13->SetLineColor(kBlack);
  fpol13->SetLineStyle(2);
  grFunc->Fit(fpol13, "EMR");

  fpol13->Draw("same");
  leg->AddEntry(fpol13, "Combine fit", "L");

  gPad->Update();
  TPaveStats *st3 = (TPaveStats*)grFunc->FindObject("stats");
  st3->SetX1NDC(0.52);st3->SetX2NDC(0.9);
  st3->SetY1NDC(0.6);st3->SetY2NDC(0.8);
  st3->SetTextColor(fpol13->GetLineColor());
  st3->Draw();
*/
//  leg->Draw();

  TLatex L;
  L.SetNDC();
  L.SetTextSize(0.03);
  Double_t xposi = 0.65;
//  L.DrawLatex(xposi, 0.25, "S2_{Corr.}=E_{Kr}g_{2}/W-g_{2}/g_{1}S1");
  L.DrawLatex(xposi, 0.25, "S2_{Corr.}=#epsilon_{2}(E_{Kr}(1+f)/W-S1/#epsilon_{1})");
  L.DrawLatex(xposi, 0.3, "W=23.6 eV");
  L.DrawLatex(xposi, 0.35, "f=0.21");
  L.SetNDC(0);
 L.DrawLatex(S1_old[0]+3.5, S2Peak[0], "E_{drift} = 200V/cm");
  L.DrawLatex(S1_old[1]+3, S2Peak[1], "E_{drift} = 150V/cm");
  L.DrawLatex(S1_old[2]+3, S2Peak[2], "E_{drift} = 100V/cm");
  L.DrawLatex(S1_Feb2015[2]+6, S2_Feb2015[2], "E_{drift} = 50V/cm");
//  L.DrawLatex(xposi, 0.7, Form("%s: %4.3e", fpol13->GetParName(0), fpol13->GetParameter(0)));
//  L.DrawLatex(xposi, 0.65, Form("%s: %4.3e", fpol13->GetParName(1), fpol13->GetParameter(1)));

  leg->Draw();

  c1->Update();
#ifndef __CINT__
  theApp->Run();
#endif /* __CINT __ */

  TString printName = Form("S1vsS2.pdf");
  c1->Print(printName.Data());

  delete c1;
}

bool readValuesFromFile(TString inputFName, vector<vector<TString> > &allInfos){
  ifstream inFile(inputFName.Data());
  TString line;
  vector<TString> values;

  while(inFile.is_open() && line.ReadLine(inFile)) {
      TObjArray *tokens = line.Tokenize(" ");
      Int_t nEntriesPerLine = tokens->GetEntries();

      for (Int_t i = 0; i < nEntriesPerLine; i++) {
          TString value = ((TObjString*) tokens->At(i))->GetString();
          values.push_back(value);
      }
      allInfos.push_back(values);
      values.clear();
      tokens->Delete(); //Remove all objects from the array AND delete all heap based objects.
      delete tokens;
  }

  inFile.close();
  return true;

}

/*
//-----------------------------------------------------------------------------------------------------
Double_t Pol1(Double_t *x, Double_t *par) {
  Double_t Woverg2 = par[0];
  Double_t g2overg1 = par[1];
  Double_t KrEnergy = 41.5; //keV Kr peak
  Double_t s1 = x[0];
  Double_t s2 = KrEnergy/Woverg2 - g2overg1*s1;
   return s2;
}
*/
//-----------------------------------------------------------------------------------------------------
Double_t Pol1(Double_t *x, Double_t *par) {
  Double_t W = 23.6e-3; //keV Energy required to make an ion-electron pair
  Double_t f = 0.21; // ratio of excited atoms to ions
  Double_t epsi1 = par[0]; // PE(S1) per excited atoms
  Double_t epsi2 = par[1]; // PE(S2) per extracted electron
  Double_t KrEnergy = 41.5; //keV Kr peak
  Double_t s1 = x[0];
  Double_t s2 = epsi2*(KrEnergy*(1.+f)/W - s1/epsi1);
   return s2;
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
	if ( theApp->Argc() == 1 ) {
		std::cout << "\n==========> PlotS2vsS1Relation <=============" << std::endl;
		std::cout << "==> Application start." << std::endl;
		PlotS2vsS1Relation();

	} else {
		std::cout << "Usage:" <<theApp->Argc()<< std::endl;
		std::cout << "./PlotS2vsS1Relation" << std::endl;
		return 0;
	}


//        std::cout << "==> Application start." << std::endl;
//	S2vsS1Fitter(theApp->Argv(1), theApp->Argv(2));
	std::cout << "==> Application finished." << std::endl;
//	theApp->Run();

	return 0;
}
#endif /* __CINT __ */
