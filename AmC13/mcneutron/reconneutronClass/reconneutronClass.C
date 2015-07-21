#include "reconneutronClass.h"

using namespace std;
ClassImp(reconneutronClass);

bool reconneutronClass::VerifyDataFile(TString filename)
{
  ifstream NameCheck;
  NameCheck.open(filename.Data());
  if(!NameCheck.good())
    return false; 
  else{
    TFile *f = new TFile(filename);
    if(f->IsZombie())
      return false;
    else{
      cout<<"Processing Data file: "<<filename<<endl;	
      return true;
    }
  }
}

void reconneutronClass::ReadDataFile(TChain *t,int start,int end,string last)
{
  string dirname = GetInPath();
  string middle=GetInMiddle();
  for(int i=start; i<=end; i++)
    {
      TString filename;
      if(i==0)
	filename.Form("%s_%s",middle.c_str(),last.c_str());
      else
	filename.Form("%s_v%d_%s",middle.c_str(),i,last.c_str());
      filename.Prepend(dirname.c_str());
      if(!VerifyDataFile(filename))
	continue;
      else t->Add(filename);     
    }
}

void reconneutronClass::Init()
{
  int time_bins = 1.e+6;
  head.push_back("tpc");
  head.push_back("tpc_NR");
  head.push_back("tpc_ER");
  head.push_back("nv");
  head.push_back("tpc_single");
  head.push_back("tpc_singleNR");
  head.push_back("tpc_singleER");
  head.push_back("nv_NR_coin");
  head.push_back("nv_ER_coin");
  head.push_back("nv_singleNR_coin");
  head.push_back("nv_singleER_coin");

  eqch_bins.push_back(35000);
  eqch_bins.push_back(3500);
  eqch_bins.push_back(8000);
  eqch_bins.push_back(35000);
  eqch_bins.push_back(6000);
  eqch_bins.push_back(3500);
  eqch_bins.push_back(6000);
  eqch_bins.push_back(35000);
  eqch_bins.push_back(35000);
  eqch_bins.push_back(35000);
  eqch_bins.push_back(35000);

  total_bins.push_back(50000);
  total_bins.push_back(35000);
  total_bins.push_back(10000);
  total_bins.push_back(5000);
  total_bins.push_back(15000);
  total_bins.push_back(10000);

  LoadTree(fChain, e);
};

void reconneutronClass::BookHistograms()
{
  source_ntuple = new TNtuple("source_ntuple","","x:y:z:r:theta:px:py:pz:ene");
  for(size_t i=0; i<head.size(); i++)
    {
      Plots t;

      t.ntuple_plots = new TNtuple(head[i].c_str(),"","quenchingfactor:x:y:z:r:theta:eqch:ene:time");
      t.eqch_sourceene_hist = new TH2F(Form("%s_eqch_sourceene",head[i].c_str()),"; Energy [keVee]; Energy [keV]",eqch_bins[i]/10,0,eqch_bins[i]*1.0,500,0,8000);
      t.eqch_sourcetheta_hist = new TH2F(Form("%s_eqch_sourcetheta",head[i].c_str()),"; Energy [keVee]; Theta [rad]",eqch_bins[i]/10,0,eqch_bins[i]*1.0,200,-4,4);
      t.eqch_hist = new TH1F(Form("%s_eqch",head[i].c_str()),"; Energy [keVee]",eqch_bins[i]/10,0,eqch_bins[i]*1.0);
      t.eqch_time_hist = new TH2F(Form("%s_eqch_time",head[i].c_str()),"; Energy [keVee]; Time [ns]",eqch_bins[i]/10,0,eqch_bins[i]*1.0,200,0,time_bins);
      t.eqch_theta_hist = new TH2F(Form("%s_eqch_theta",head[i].c_str()),"; Energy [keVee]; Theta [rad]",eqch_bins[i]/10,0,eqch_bins[i]*1.0,200,-4,4);
      t.eqch_radius_hist = new TH2F(Form("%s_eqch_radius",head[i].c_str()),"; Energy [keVee]; Radius [cm]",eqch_bins[i]/10,0,eqch_bins[i]*1.0,200,0,200);

      p.push_back(t);
    }
}

bool CompareTime(double a, double b)
{
  return (a<b);
}

void reconneutronClass::FillHistograms()
{
  int nentries = fChain->GetEntries();
  for(int i=0; i<nentries; i++)
    {
      fChain->GetEntry(i);
      if(!(100*i%nentries))
	std::cout<<"Processing event "<<i<<std::endl;
      
      double quenchingcut = 0.43;
      std::vector<double> tpcNR_trigger;
      std::vector<double> tpcER_trigger;
      std::vector<double> tpcsingleNR_trigger;
      std::vector<double> tpcsingleER_trigger;
      bool veto = false;
      int single = 0;
      bool IsNR = false;
      bool IsER = false;
      
      e.sourcer = TMath::Sqrt(TMath::Power(e.source_x,2)+TMath::Power(e.source_y,2));
      e.sourcetheta = TMath::ATan2(e.source_x,e.source_y);
      source_ntuple->Fill(e.source_x,e.source_y,e.source_z,e.sourcer,e.sourcetheta,e.source_px,e.source_py,e.source_pz,e.source_ene);

      for(size_t j=0; j<e.et->size(); j++)
	{
	  double er = TMath::Sqrt(TMath::Power(e.ex->at(j),2)+TMath::Power(e.ey->at(j),2));
	  double theta = TMath::ATan2(e.ex->at(j),e.ey->at(j));
	  p.at(0).ntuple_plots->Fill(e.quenchingfactor->at(j),e.ex->at(j),e.ey->at(j),e.ez->at(j),er,theta,e.eqch->at(j),e.edep->at(j),e.et->at(j)*1.e+9);
	  p.at(0).eqch_hist->Fill(e.eqch->at(j));
	  p.at(0).eqch_time_hist->Fill(e.eqch->at(j),e.et->at(j)*1.e+9);
	  p.at(0).eqch_sourceene_hist->Fill(e.eqch->at(j),e.source_ene*1000.);
	  p.at(0).eqch_sourcetheta_hist->Fill(e.eqch->at(j),e.sourcetheta);
	  p.at(0).eqch_theta_hist->Fill(e.eqch->at(j),theta);
	  p.at(0).eqch_radius_hist->Fill(e.eqch->at(j),er);

	  if(e.volume->at(j) == "p_active")
	    {
	      ++single;
 	      if(e.quenchingfactor->at(j)<quenchingcut)
		{		
		  IsNR = true;
		}
	      if(e.quenchingfactor->at(j)>quenchingcut)
		{		
		  IsER = true;
		}
	    }
	  if(e.volume->at(j) == "p_scint")
	    { 
	      veto = true;	    
	    }
	}      
    }
}

void reconneutronClass::LoopOverEvent(Long64_t entry)
{
  

}

void reconneutronClass::SaveHistograms()
{
  string output = GetOutPath() + GetOutFile().Data();
  TFile outfile(output.c_str(), "RECREATE");
  source_ntuple->Write();

  for(size_t j=0; j<p.size(); j++)
    {
      p.at(j).ntuple_plots->Write();
      p.at(j).eqch_sourceene_hist->Write();
      p.at(j).eqch_sourcetheta_hist->Write();
      p.at(j).eqch_hist->Write();
      p.at(j).eqch_time_hist->Write();
      p.at(j).eqch_theta_hist->Write();
      p.at(j).eqch_radius_hist->Write();

    }

  outfile.Write();
  outfile.Close();
}
