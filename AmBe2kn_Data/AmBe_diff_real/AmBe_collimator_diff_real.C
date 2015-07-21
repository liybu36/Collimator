#include <vector>
#include <iostream>

#include "TH2F.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TROOT.h"
#include "TRint.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TAxis.h"

#ifdef __MAKECINT__
#pragma link C++ class std::vector < std::vector<double> >+;
#pragma link C++ class std::vector < std::vector<TH2F*> >+;
#pragma link C++ class std::vector <TFile*>+;
#endif

TRint *theApp;
using namespace std;

void AmBe_collimator_diff_real()
{
  string indir = "/darkside/users/hqian/AmBe2kn_XYData/";

 #define Low
#ifdef Low
  string collimator = indir +  "AmBe2knXY_1046510480_May12.root";
  string bg = indir + "AmBe2knXY_1048310493_May12.root";
#else
  string collimator = indir +  "AmBe2knXY_1044410456_May12.root";
  string bg = indir + "AmBe2knXY_1043410441_May12.root";
#endif

  const int N =2;
  string name[N] = {collimator, bg};
  std::vector<TFile*> f;
  std::vector< std::vector<TH2F*> > Hist2DsingleNR;
  std::vector< std::vector<TH1D*> > HistsingleNRpy;
  std::vector< std::vector<TH1D*> > HistsingleNRpx;
  std::vector< std::vector<double> > sum2D;
  std::vector< std::vector<double> > sum1Dx;
  std::vector< std::vector<double> > sum1Dy;

  for(int i=0; i<N; i++)
    {
      f.push_back(new TFile(name[i].c_str()));
      std::vector<TH2F*> histsinglenr;
      histsinglenr.push_back((TH2F*)f.at(i)->Get("neutron_theta_tdrift_slice"));         
      histsinglenr.push_back((TH2F*)f.at(i)->Get("neutron_xy_slice"));         
      histsinglenr.push_back((TH2F*)f.at(i)->Get("single_xy_slice"));         
      histsinglenr.push_back((TH2F*)f.at(i)->Get("singleNR_xy_slice"));         
      Hist2DsingleNR.push_back(histsinglenr);
 
      std::vector<TH1D*> histpy;
      histpy.push_back((TH1D*)f.at(i)->Get("neutron_xy_py"));
      histpy.push_back((TH1D*)f.at(i)->Get("single_xy_py"));
      histpy.push_back((TH1D*)f.at(i)->Get("singleNR_xy_py"));
      HistsingleNRpy.push_back(histpy);

      std::vector<TH1D*> histpx;
      histpx.push_back((TH1D*)f.at(i)->Get("neutron_xy_px"));
      histpx.push_back((TH1D*)f.at(i)->Get("single_xy_px"));
      histpx.push_back((TH1D*)f.at(i)->Get("singleNR_xy_px"));
      HistsingleNRpx.push_back(histpx);
      
      std::vector<double> sum_temp2d;
      std::vector<double> sum_temp1dx;
      std::vector<double> sum_temp1dy;

      for(size_t j=0; j<Hist2DsingleNR.at(i).size(); j++)
	{
	  int rebin2D = 1;
	  Hist2DsingleNR[i][j]->Rebin2D(rebin2D,rebin2D);
	  int binx11(0),binx12(0),biny11(0),biny12(0);
	  if(j<1)
	    {
	      binx11 = Hist2DsingleNR[i][j]->GetXaxis()->FindBin(100.);
	      binx12 = Hist2DsingleNR[i][j]->GetXaxis()->FindBin(150.);
	      biny11 = Hist2DsingleNR[i][j]->GetYaxis()->FindBin(-1.);
	      biny12 = Hist2DsingleNR[i][j]->GetYaxis()->FindBin(0.);
	    }
	  else
	    {
	      binx11 = Hist2DsingleNR[i][j]->GetXaxis()->FindBin(-10.);
	      binx12 = Hist2DsingleNR[i][j]->GetXaxis()->FindBin(-5.);
	      biny11 = Hist2DsingleNR[i][j]->GetYaxis()->FindBin(5.);
	      biny12 = Hist2DsingleNR[i][j]->GetYaxis()->FindBin(10.);
	    }

	  Hist2DsingleNR[i][j]->Sumw2();
	  sum_temp2d.push_back(Hist2DsingleNR[i][j]->Integral(binx11,binx12,biny11,biny12));
	}	  
      sum2D.push_back(sum_temp2d);

      for(size_t j=0; j<HistsingleNRpy.at(i).size(); j++)
	{
	  int rebin1D = 2;
	  HistsingleNRpy[i][j]->Rebin(rebin1D);
	  int binx21 = HistsingleNRpy[i][j]->GetXaxis()->FindBin(0.);
	  int binx22 = HistsingleNRpy[i][j]->GetXaxis()->FindBin(5.);

	  HistsingleNRpy[i][j]->Sumw2();
	  sum_temp1dy.push_back(HistsingleNRpy[i][j]->Integral(binx21,binx22));
	}
      sum1Dy.push_back(sum_temp1dy);

      for(size_t j=0; j<HistsingleNRpx.at(i).size(); j++)
	{
	  int rebin1D = 2;
	  HistsingleNRpx[i][j]->Rebin(rebin1D);
	  int binx31 = HistsingleNRpx[i][j]->GetXaxis()->FindBin(0.);
	  int binx32 = HistsingleNRpx[i][j]->GetXaxis()->FindBin(5.);

	  HistsingleNRpx[i][j]->Sumw2();
	  sum_temp1dx.push_back(HistsingleNRpx[i][j]->Integral(binx31,binx32));
	}
      sum1Dx.push_back(sum_temp1dx);
      
    }
  
  std::vector<TCanvas*> canv;  
  for(size_t j=0; j<Hist2DsingleNR.at(0).size(); j++)
    {
      Hist2DsingleNR[1][j]->Scale(sum2D[0][j]/sum2D[1][j]);      
      Hist2DsingleNR[0][j]->Add(Hist2DsingleNR[1][j],-1);
      canv.push_back(new TCanvas(Form("c1%d",static_cast<int>(j)),"",600,400));
      canv.at(j)->cd();
      Hist2DsingleNR[1][j]->Draw("colz");
    }
  for(size_t j=0; j<HistsingleNRpy.at(0).size(); j++)
    {
      HistsingleNRpy[1][j]->Scale(sum1Dy[0][j]/sum1Dy[1][j]);
      std::cout<<sum1Dy[0][j]/sum1Dy[1][j]<<std::endl;
      HistsingleNRpy[0][j]->SetLineColor(2); //kRed
      HistsingleNRpy[1][j]->SetLineColor(4); //kBlue
      canv.push_back(new TCanvas(Form("c2%d",static_cast<int>(j)),"",600,400));
      canv.back()->cd();
      HistsingleNRpy[0][j]->Draw();
      HistsingleNRpy[1][j]->Draw("same");
    }
  for(size_t j=0; j<HistsingleNRpx.at(0).size(); j++)
    {
      HistsingleNRpx[1][j]->Scale(sum1Dx[0][j]/sum1Dx[1][j]);
      std::cout<<sum1Dx[0][j]/sum1Dx[1][j]<<std::endl;
      HistsingleNRpx[0][j]->SetLineColor(2); //kRed
      HistsingleNRpx[1][j]->SetLineColor(4); //kBlue
      canv.push_back(new TCanvas(Form("c3%d",static_cast<int>(j)),"",600,400));
      canv.back()->cd();
      HistsingleNRpx[0][j]->Draw();
      HistsingleNRpx[1][j]->Draw("same");
    }

  
#ifdef Low
  string output1=indir+"diffAmBe2knXY_1046510493_May12.root";
#else
  string output1=indir+"diffAmBe2knXY_1043410456_May12.root";
#endif
  TFile out1(output1.c_str(),"RECREATE");
  for(size_t j=0; j<canv.size(); j++)
    canv.at(j)->Write();
  out1.Write();
  out1.Close();


}


#ifndef __CINT__
int main() {
  theApp = new TRint("theApp",NULL,NULL,NULL,0);
  //  gROOT->ProcessLine("#pragma link C++ class std::vector < std::vector<double> >");
  AmBe_collimator_diff_real();

  return 1;
}
#endif
