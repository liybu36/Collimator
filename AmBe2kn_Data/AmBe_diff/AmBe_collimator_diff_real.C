#include <vector>

#include "TH2F.h"
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

#define High
#ifdef High
  string collimator = indir +  "AmBe2knXY_1046510480_May12.root";
  string nocollimator = indir + "AmBe2knXY_1048310493_May12.root";
#else
  string collimator = indir +  "AmBe2knXY_1044410456_May12.root";
  string nocollimator = indir + "AmBe2knXY_1043410441_May12.root";
#endif

  const int N =2;
  string name[N] = {collimator, nocollimator};
  std::vector<TFile*> f;
  std::vector< std::vector<TH2F*> > Hist2DsingleNR;
  std::vector< std::vector<double> > sum;

  for(int i=0; i<N; i++)
    {
      f.push_back(new TFile(name[i].c_str()));
      std::vector<TH2F*> histsinglenr;
      histsinglenr.push_back((TH2F*)f.at(i)->Get("neutron_theta_tdrift_slice"));   
      
      Hist2DsingleNR.push_back(histsinglenr);
      std::vector<double> sum_temp;
  
      for(size_t j=0; j<Hist2DsingleNR.at(i).size(); j++)
	{
	  int rebin2D = 1;
	  int binx11 = Hist2DsingleNR[i][j]->GetXaxis()->FindBin(100.);
	  int binx12 = Hist2DsingleNR[i][j]->GetXaxis()->FindBin(150.);
	  int biny11 = Hist2DsingleNR[i][j]->GetYaxis()->FindBin(-1.);
	  int biny12 = Hist2DsingleNR[i][j]->GetYaxis()->FindBin(0.);
      
	  Hist2DsingleNR[i][j]->Rebin2D(rebin2D,rebin2D);
	  Hist2DsingleNR[i][j]->Sumw2();
	  sum_temp.push_back(Hist2DsingleNR[i][j]->Integral(binx11,binx12,biny11,biny12));
	}	  
      sum.push_back(sum_temp);
    }
  
  std::vector<TCanvas*> canv;  
  for(size_t j=0; j<Hist2DsingleNR.at(0).size(); j++)
    {
      Hist2DsingleNR[1][j]->Scale(sum[1][j]/sum[0][j]);      
      Hist2DsingleNR[1][j]->Add(Hist2DsingleNR[0][j],-1);
      canv.push_back(new TCanvas(Form("c%d",static_cast<int>(j)),"",600,400));
      canv.at(j)->cd();
      Hist2DsingleNR[1][j]->Draw("colz");
    }

#ifdef High
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
