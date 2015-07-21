#define reconneutronSelector_cxx
// The class definition in reconneutronSelector.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("reconneutronSelector.C")
// Root > T->Process("reconneutronSelector.C","some options")
// Root > T->Process("reconneutronSelector.C+")
//

#include "reconneutronSelector.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
#include <vector>

#include <TH2.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TBox.h>
#include <TStyle.h>
#include <TMath.h>
#include <TFile.h>

using namespace std;

void reconneutronSelector::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   SetTail();
}

void reconneutronSelector::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

  Info("SlaveBegin()","beginning .....");
  TString option = GetOption();
  BookHistograms();
  
}

Bool_t reconneutronSelector::Process(Long64_t entry)
{
  
  fChain->GetEntry(entry);
  //  if(entry%10000) std::cout<<"entry "<<entry<<std::endl;

  FillHistograms();

  return kTRUE;
}

void reconneutronSelector::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void reconneutronSelector::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
  string label = GetOption();
  TList* list = GetOutputList();
  source_ntuple = dynamic_cast<TNtuple*>(list->FindObject(Form("source_ntuple")));  
  source_ene_theta = dynamic_cast<TH2F*>(list->FindObject(Form("source_ene_theta")));  

  for(size_t j=0; j<ntuple_name.size(); j++)
    {
      ntuple_plots.at(j) = dynamic_cast<TNtuple*>(list->FindObject(Form(ntuple_name.at(j))));
      eqch_hist.at(j) = dynamic_cast<TH1F*>(list->FindObject(Form(eqch_name.at(j))));
      eqch_time_hist.at(j) = dynamic_cast<TH2F*>(list->FindObject(Form(eqch_time_name.at(j))));
      eqch_sourceene_hist.at(j) = dynamic_cast<TH2F*>(list->FindObject(Form(eqch_sourceene_name.at(j))));
      eqch_sourcetheta_hist.at(j) = dynamic_cast<TH2F*>(list->FindObject(Form(eqch_sourcetheta_name.at(j))));
      eqch_theta_hist.at(j) = dynamic_cast<TH2F*>(list->FindObject(Form(eqch_theta_name.at(j))));
      eqch_radius_hist.at(j) = dynamic_cast<TH2F*>(list->FindObject(Form(eqch_radius_name.at(j))));

    }
  for(size_t j=0; j<tpc_total_name.size(); j++)
    tpc_total.at(j) = dynamic_cast<TH1F*>(list->FindObject(Form(tpc_total_name.at(j))));
  for(size_t j=0; j<tpc_total_nv_name.size(); j++)
    tpc_total_nv.at(j) = dynamic_cast<TH2F*>(list->FindObject(Form(tpc_total_nv_name.at(j))));
  int entries = fChain->GetEntries();
  std::cout<<"Total Events "<<entries<<std::endl;
  for(size_t j=0; j<tpc_events.size(); j++)
    std::cout<<"tpc "<<j<<"  "<<tpc_events.at(j)<<"  events  "<<1.0*tpc_events.at(j)/entries<<std::endl;
  for(size_t j=0; j<veto_events.size(); j++)
    std::cout<<"veto "<<j<<"  "<<veto_events.at(j)<<"  events  "<<1.0*veto_events.at(j)/entries<<std::endl;

  string output = label;
  TFile outfile(output.c_str(), "RECREATE");
  source_ntuple->Write();
  source_ene_theta->Write();
  for(size_t j=0; j<tpc_total_name.size(); j++)
    tpc_total.at(j)->Write();
  for(size_t j=0; j<tpc_total_nv_name.size(); j++)
    tpc_total_nv.at(j)->Write();
  for(size_t j=0; j<ntuple_name.size(); j++)
    {
      ntuple_plots.at(j)->Write();
      eqch_hist.at(j)->Write();
      eqch_time_hist.at(j)->Write();
      eqch_sourceene_hist.at(j)->Write();
      eqch_sourcetheta_hist.at(j)->Write();
      eqch_theta_hist.at(j)->Write();
      eqch_radius_hist.at(j)->Write();
    }
  outfile.Write();
  outfile.Close();

  Info("Terminate()","terminating ...");
}

int reconneutronSelector::SetTail()
{  
  tail.push_back("prompt");
  tail.push_back("signal");
  tail.push_back("first");
  tail.push_back("late");

  timelow.push_back(0.);
  timelow.push_back(20.e+3);
  timelow.push_back(20.e+3);
  timelow.push_back(120.e+3);

  timehigh.push_back(20.e+3);
  timehigh.push_back(120.e+3);
  timehigh.push_back(120.e+3);
  timehigh.push_back(220.e+3);

  return (int)tail.size();
}

bool CompareTime(double a, double b)
{
  return (a<b);
}

void reconneutronSelector::FillHistograms()
{
  double quenchingcut = 0.43;
  std::vector<double> tpcNR_trigger;
  std::vector<double> tpcER_trigger; 
  std::vector<double> tpcsingleNR_trigger;
  std::vector<double> tpcsingleER_trigger; 
  bool veto = false;
  int single = 0;
  bool IsNR = false;
  bool IsER = false;
  std::vector<double> total (6,0.);
  std::vector<bool> tpc_tag (8,false);
  std::vector<bool> veto_tag (5,false);

  double sourcer = TMath::Sqrt(TMath::Power(source_x,2)+TMath::Power(source_y,2));
  double sourcetheta = TMath::ATan2(source_y,source_x);
  source_ntuple->Fill(source_x,source_y,source_z,sourcer,sourcetheta,source_px,source_py,source_pz,source_ene);
  source_ene_theta->Fill(source_ene*1000.,sourcetheta);

  for(size_t j=0; j<et->size(); j++)
    {
      double er = TMath::Sqrt(TMath::Power(ex->at(j),2)+TMath::Power(ey->at(j),2));
      double theta = TMath::ATan2(ey->at(j),ex->at(j));
      ntuple_plots.at(0)->Fill(quenchingfactor->at(j),ex->at(j),ey->at(j),ez->at(j),er,theta,eqch->at(j),edep->at(j),et->at(j)*1.e+9);
      eqch_hist.at(0)->Fill(eqch->at(j));
      eqch_time_hist.at(0)->Fill(eqch->at(j),et->at(j)*1.e+9);
      eqch_sourceene_hist.at(0)->Fill(eqch->at(j),source_ene*1000.);
      eqch_sourcetheta_hist.at(0)->Fill(eqch->at(j),sourcetheta);
      eqch_theta_hist.at(0)->Fill(eqch->at(j),theta);
      eqch_radius_hist.at(0)->Fill(eqch->at(j),er);
      total.at(0) += edep->at(j);
      total.at(1) += eqch->at(j);

      if(volume->at(j) == "p_active")
	{
	  tpc_tag.at(0) = true;
	  ++single;
	  if(quenchingfactor->at(j)<quenchingcut)
	    {
	      tpcNR_trigger.push_back(et->at(j)*1.e+9);
	      ntuple_plots.at(1)->Fill(quenchingfactor->at(j),ex->at(j),ey->at(j),ez->at(j),er,theta,eqch->at(j),edep->at(j),et->at(j)*1.e+9);
	      eqch_hist.at(1)->Fill(eqch->at(j));
	      eqch_time_hist.at(1)->Fill(eqch->at(j),et->at(j)*1.e+9);
	      eqch_sourceene_hist.at(1)->Fill(eqch->at(j),source_ene*1000.);
	      eqch_sourcetheta_hist.at(1)->Fill(eqch->at(j),sourcetheta);
	      eqch_theta_hist.at(1)->Fill(eqch->at(j),theta);
	      eqch_radius_hist.at(1)->Fill(eqch->at(j),er);
	      total.at(2) += edep->at(j);
	      total.at(3) += eqch->at(j);
	      tpc_tag.at(1) = true;
	      IsNR = true;
	    }
	  if(quenchingfactor->at(j)>quenchingcut)
	    {
	      tpcER_trigger.push_back(et->at(j)*1.e+9);
	      ntuple_plots.at(2)->Fill(quenchingfactor->at(j),ex->at(j),ey->at(j),ez->at(j),er,theta,eqch->at(j),edep->at(j),et->at(j)*1.e+9);
	      eqch_hist.at(2)->Fill(eqch->at(j));
	      eqch_time_hist.at(2)->Fill(eqch->at(j),et->at(j)*1.e+9);
	      eqch_sourceene_hist.at(2)->Fill(eqch->at(j),source_ene*1000.);
	      eqch_sourcetheta_hist.at(2)->Fill(eqch->at(j),sourcetheta);
	      eqch_theta_hist.at(2)->Fill(eqch->at(j),theta);
	      eqch_radius_hist.at(2)->Fill(eqch->at(j),er);
	      total.at(4) += edep->at(j);
	      total.at(5) += eqch->at(j);
	      tpc_tag.at(2) = true;
	      IsER = true;
	    }	  	  
	}
      if(volume->at(j) == "p_scint")
	{ veto = true;
	  ntuple_plots.at(3)->Fill(quenchingfactor->at(j),ex->at(j),ey->at(j),ez->at(j),er,theta,eqch->at(j),edep->at(j),et->at(j)*1.e+9);
	  eqch_hist.at(3)->Fill(eqch->at(j));
	  eqch_time_hist.at(3)->Fill(eqch->at(j),et->at(j)*1.e+9);
	  eqch_sourceene_hist.at(3)->Fill(eqch->at(j),source_ene*1000.);
	  eqch_sourcetheta_hist.at(3)->Fill(eqch->at(j),sourcetheta);
	  eqch_theta_hist.at(3)->Fill(eqch->at(j),theta);
	  eqch_radius_hist.at(3)->Fill(eqch->at(j),er);
	  veto_tag.at(0) = true;
	}      
    }
  for(int i=0; i<6; i++)
    tpc_total.at(i)->Fill(total.at(i));

  for(size_t j=0;j<et->size(); j++)
    {
      double er = TMath::Sqrt(TMath::Power(ex->at(j),2)+TMath::Power(ey->at(j),2));
      double theta = TMath::ATan2(ey->at(j),ex->at(j));
      
      if(volume->at(j) == "p_scint")
	{
	  tpc_total_nv.at(0)->Fill(total.at(1),eqch->at(j));
	  tpc_total_nv.at(1)->Fill(total.at(3),eqch->at(j));
	  tpc_total_nv.at(2)->Fill(total.at(5),eqch->at(j));
	}     
      if(volume->at(j) == "p_active")
	{
	  if(single==1)
	    {
	      ntuple_plots.at(4)->Fill(quenchingfactor->at(j),ex->at(j),ey->at(j),ez->at(j),er,theta,eqch->at(j),edep->at(j),et->at(j)*1.e+9);
	      eqch_hist.at(4)->Fill(eqch->at(j));
	      eqch_time_hist.at(4)->Fill(eqch->at(j),et->at(j)*1.e+9);
	      eqch_sourceene_hist.at(4)->Fill(eqch->at(j),source_ene*1000.);
	      eqch_sourcetheta_hist.at(4)->Fill(eqch->at(j),sourcetheta);
	      eqch_theta_hist.at(4)->Fill(eqch->at(j),theta);
	      eqch_radius_hist.at(4)->Fill(eqch->at(j),er);		  
	      tpc_tag.at(3) = true;
	    }
	  if(IsNR && !IsER && quenchingfactor->at(j)<quenchingcut)
	    { tpc_tag.at(4) = true;
	      if(single==1)
		{
		  tpcsingleNR_trigger.push_back(et->at(j)*1.e+9);
		  ntuple_plots.at(5)->Fill(quenchingfactor->at(j),ex->at(j),ey->at(j),ez->at(j),er,theta,eqch->at(j),edep->at(j),et->at(j)*1.e+9);
		  eqch_hist.at(5)->Fill(eqch->at(j));
		  eqch_time_hist.at(5)->Fill(eqch->at(j),et->at(j)*1.e+9);
		  eqch_sourceene_hist.at(5)->Fill(eqch->at(j),source_ene*1000.);
		  eqch_sourcetheta_hist.at(5)->Fill(eqch->at(j),sourcetheta);
		  eqch_theta_hist.at(5)->Fill(eqch->at(j),theta);
		  eqch_radius_hist.at(5)->Fill(eqch->at(j),er);		  
		  tpc_tag.at(5) = true;
		}
	    }
	  if(IsER && !IsNR && quenchingfactor->at(j)>quenchingcut)
	    { tpc_tag.at(6) = true;
	      if(single==1)
		{		  
		  tpcsingleER_trigger.push_back(et->at(j)*1.e+9);
		  ntuple_plots.at(6)->Fill(quenchingfactor->at(j),ex->at(j),ey->at(j),ez->at(j),er,theta,eqch->at(j),edep->at(j),et->at(j)*1.e+9);
		  eqch_hist.at(6)->Fill(eqch->at(j));
		  eqch_time_hist.at(6)->Fill(eqch->at(j),et->at(j)*1.e+9);
		  eqch_sourceene_hist.at(6)->Fill(eqch->at(j),source_ene*1000.);
		  eqch_sourcetheta_hist.at(6)->Fill(eqch->at(j),sourcetheta);
		  eqch_theta_hist.at(6)->Fill(eqch->at(j),theta);
		  eqch_radius_hist.at(6)->Fill(eqch->at(j),er);		  
		  tpc_tag.at(7) = true;
		}
	    }       
	}
    }
  
  if(tpcNR_trigger.size()>1)
    std::sort (tpcNR_trigger.begin(), tpcNR_trigger.end(), CompareTime);
  if(tpcER_trigger.size()>1)
    std::sort (tpcER_trigger.begin(), tpcER_trigger.end(), CompareTime);
  if(tpcsingleNR_trigger.size()>1)
    std::sort (tpcsingleNR_trigger.begin(), tpcsingleNR_trigger.end(), CompareTime);
  if(tpcsingleER_trigger.size()>1)
    std::sort (tpcsingleER_trigger.begin(), tpcsingleER_trigger.end(), CompareTime);
  
  int NUM = 4;
  int offset = 7;
  for(int k=0;veto&&k<NUM; k++)
    {
      double trigger = 0;
      if(k==0 && tpcNR_trigger.size()>0)
	trigger = tpcNR_trigger.front();
      else if(k==1 && tpcER_trigger.size()>0)
	trigger = tpcER_trigger.front();             
      else if(k==2 && tpcsingleNR_trigger.size()>0)
	trigger = tpcsingleNR_trigger.front();
      else if(k==3 && tpcsingleER_trigger.size()>0)
	trigger = tpcsingleER_trigger.front();             
      else continue;

      for(size_t j=0; j<et->size(); j++)
	{
	  double er = TMath::Sqrt(TMath::Power(ex->at(j),2)+TMath::Power(ey->at(j),2));
	  double theta = TMath::ATan2(ey->at(j),ex->at(j));
	  double gps = et->at(j)*1.e+9 - trigger;
	  if(volume->at(j) == "p_scint")
	    {
	      size_t max=tail.size()+1;
	      ntuple_plots.at(offset+k*max)->Fill(quenchingfactor->at(j),ex->at(j),ey->at(j),ez->at(j),er,theta,eqch->at(j),edep->at(j),gps);
	      eqch_hist.at(offset+k*max)->Fill(eqch->at(j));
	      eqch_time_hist.at(offset+k*max)->Fill(eqch->at(j),gps);
	      eqch_sourceene_hist.at(offset+k*max)->Fill(eqch->at(j),source_ene*1000.);
	      eqch_sourcetheta_hist.at(offset+k*max)->Fill(eqch->at(j),sourcetheta);
	      eqch_theta_hist.at(offset+k*max)->Fill(eqch->at(j),theta);
	      eqch_radius_hist.at(offset+k*max)->Fill(eqch->at(j),er);
	      veto_tag.at(k+1) = true;
	      
	      for(size_t h=0; h<tail.size(); h++)
		{
		  if(gps>timelow[h] && gps<timehigh[h])
		    {
		      ntuple_plots.at(offset+k*max+h+1)->Fill(quenchingfactor->at(j),ex->at(j),ey->at(j),ez->at(j),er,theta,eqch->at(j),edep->at(j),gps);
		      eqch_hist.at(offset+k*max+h+1)->Fill(eqch->at(j));
		      eqch_time_hist.at(offset+k*max+h+1)->Fill(eqch->at(j),gps);
		      eqch_sourceene_hist.at(offset+k*max+h+1)->Fill(eqch->at(j),source_ene*1000.);
		      eqch_sourcetheta_hist.at(offset+k*max+h+1)->Fill(eqch->at(j),sourcetheta);
		      eqch_theta_hist.at(offset+k*max+h+1)->Fill(eqch->at(j),theta);
		      eqch_radius_hist.at(offset+k*max+h+1)->Fill(eqch->at(j),er);
		    }
		}
	    }      
	}
    }
  
  for(size_t j=0; j<tpc_tag.size(); j++)
    {
      if(tpc_tag.at(j))
	++tpc_events.at(j);
    }
  for(size_t j=0; j<veto_tag.size(); j++)
    {
      if(veto_tag.at(j))
	++veto_events.at(j);
    }
}

void reconneutronSelector::BookHistograms()
{
  string label = GetOption();
  TList* list = GetOutputList();
  const int M = 11;
  string head[M]   = {"tpc","tpc_NR","tpc_ER","nv","tpc_single","tpc_singleNR","tpc_singleER"
		      ,"nv_NR_coin","nv_ER_coin","nv_singleNR_coin","nv_singleER_coin"};
  int eqch_bins[M] = {35000,3500,8000,35000,6000,3500,6000,35000,35000,35000,35000};
  int time_bins = 1.e+6;
  int total_bins[6] = {50000,35000,10000,5000,15000,10000};
  //  const int N = 4;
  //  string tail[N] = {"prompt","signal","first","late"};

  source_ntuple = new TNtuple("source_ntuple","","x:y:z:r:theta:px:py:pz:ene");
  list->Add(source_ntuple);  
  source_ene_theta = new TH2F("source_ene_theta","; Energy [keV]; theta [rad]",1000,0,8000,100,-4,4);
  list->Add(source_ene_theta);

  for(int i=0; i<3; ++i)
    {
      tpc_total_name.push_back(Form("%s_totalene",head[i].c_str()));
      tpc_total.push_back(new TH1F(tpc_total_name.back(),"; Energy [keV]",total_bins[i*2]/10,0,total_bins[i*2]));
      list->Add(tpc_total.back());
      
      tpc_total_name.push_back(Form("%s_totaleqch",head[i].c_str()));
      tpc_total.push_back(new TH1F(tpc_total_name.back(),"; Energy [keVee]",total_bins[i*2+1]/10,0,total_bins[i*2+1]));
      list->Add(tpc_total.back());

      tpc_total_nv_name.push_back(Form("%s_totaleqch_nveqch",head[i].c_str()));
      tpc_total_nv.push_back(new TH2F(tpc_total_nv_name.back(),";TPC Energy [keVee]; Veto Energy [keVee]",total_bins[i*2+1]/10,0,total_bins[i*2+1],3500,0,35000));
      list->Add(tpc_total_nv.back());
    }
  for(int i=0; i<M; i++)
    {
      ntuple_name.push_back(head[i].c_str());
      ntuple_plots.push_back(new TNtuple(ntuple_name.back(),"","quenchingfactor:x:y:z:r:theta:eqch:ene:time"));
      list->Add(ntuple_plots.back());

      eqch_sourceene_name.push_back(Form("%s_eqch_sourceene",head[i].c_str()));
      eqch_sourceene_hist.push_back(new TH2F(eqch_sourceene_name.back(),"; Energy [keVee]; Energy [keV]",eqch_bins[i]/10,0,eqch_bins[i]*1.0,500,0,8000));
      list->Add(eqch_sourceene_hist.back());

      eqch_sourcetheta_name.push_back(Form("%s_eqch_sourcetheta",head[i].c_str()));
      eqch_sourcetheta_hist.push_back(new TH2F(eqch_sourcetheta_name.back(),"; Energy [keVee]; Theta [rad]",eqch_bins[i]/10,0,eqch_bins[i]*1.0,200,-4,4));
      list->Add(eqch_sourcetheta_hist.back());
      
      eqch_name.push_back(Form("%s_eqch",head[i].c_str()));
      eqch_hist.push_back(new TH1F(eqch_name.back(),"; Energy [keVee]",eqch_bins[i]/10,0,eqch_bins[i]*1.0));
      list->Add(eqch_hist.back());

      eqch_time_name.push_back(Form("%s_eqch_time",head[i].c_str()));
      eqch_time_hist.push_back(new TH2F(eqch_time_name.back(),"; Energy [keVee]; Time [ns]",eqch_bins[i]/10,0,eqch_bins[i]*1.0,200,0,time_bins));
      list->Add(eqch_time_hist.back());
      
      eqch_theta_name.push_back(Form("%s_eqch_theta",head[i].c_str()));
      eqch_theta_hist.push_back(new TH2F(eqch_theta_name.back(),"; Energy [keVee]; Theta [rad]",eqch_bins[i]/10,0,eqch_bins[i]*1.0,200,-4,4));
      list->Add(eqch_theta_hist.back());

      eqch_radius_name.push_back(Form("%s_eqch_radius",head[i].c_str()));
      eqch_radius_hist.push_back(new TH2F(eqch_radius_name.back(),"; Energy [keVee]; Radius [cm]",eqch_bins[i]/10,0,eqch_bins[i]*1.0,200,0,200));
      list->Add(eqch_radius_hist.back());

      for(size_t j=0; i>6 && j<tail.size(); j++)
	{	  
	  ntuple_name.push_back(Form("%s_%s",head[i].c_str(),tail[j].c_str()));
	  ntuple_plots.push_back(new TNtuple(ntuple_name.back(),"","quenchingfactor:x:y:z:r:theta:eqch:ene:time"));
	  list->Add(ntuple_plots.back());
	  
	  eqch_sourceene_name.push_back(Form("%s_eqch_sourceene_%s",head[i].c_str(),tail[j].c_str()));
	  eqch_sourceene_hist.push_back(new TH2F(eqch_sourceene_name.back(),"; Energy [keVee]; Energy [keV]",eqch_bins[i]/10,0,eqch_bins[i]*1.0,500,0,8000));
	  list->Add(eqch_sourceene_hist.back());

	  eqch_sourcetheta_name.push_back(Form("%s_eqch_sourcetheta_%s",head[i].c_str(),tail[j].c_str()));
	  eqch_sourcetheta_hist.push_back(new TH2F(eqch_sourcetheta_name.back(),"; Energy [keVee]; Theta [rad]",eqch_bins[i]/10,0,eqch_bins[i]*1.0,200,-4,4));
	  list->Add(eqch_sourcetheta_hist.back());
	  
	  eqch_name.push_back(Form("%s_eqch_%s",head[i].c_str(),tail[j].c_str()));
	  eqch_hist.push_back(new TH1F(eqch_name.back(),"; Energy [keVee]",eqch_bins[i]/10,0,eqch_bins[i]*1.0));
	  list->Add(eqch_hist.back());

	  eqch_time_name.push_back(Form("%s_eqch_time_%s",head[i].c_str(),tail[j].c_str()));
	  eqch_time_hist.push_back(new TH2F(eqch_time_name.back(),"; Energy [keVee]; Time [ns]",eqch_bins[i]/10,0,eqch_bins[i]*1.0,200,timelow[j],timehigh[j]));
	  list->Add(eqch_time_hist.back());
      
	  eqch_theta_name.push_back(Form("%s_eqch_theta_%s",head[i].c_str(),tail[j].c_str()));
	  eqch_theta_hist.push_back(new TH2F(eqch_theta_name.back(),"; Energy [keVee]; Theta [rad]",eqch_bins[i]/10,0,eqch_bins[i]*1.0,200,-4,4));
	  list->Add(eqch_theta_hist.back());
	  
	  eqch_radius_name.push_back(Form("%s_eqch_radius_%s",head[i].c_str(),tail[j].c_str()));
	  eqch_radius_hist.push_back(new TH2F(eqch_radius_name.back(),"; Energy [keVee]; Radius [cm]",eqch_bins[i]/10,0,eqch_bins[i]*1.0,200,0,200));
	  list->Add(eqch_radius_hist.back());

	}

    }
}







  /*
  const int N=9;
  string tail[N]  = {"quenchingfactor","x","y","z","r","theta","eqch","ene","time"};
  string unit[N]  = {"Quenching Factor","X [cm]","Y [cm]","Z [cm]","Radius [cm]","theta [rad]","Energy [keVee]","Energy [keV]","Time [ns]"};
  for(int i=0; i<N; i++)
    {
      tpc_gamma1d_name.push_back(Form("tpc_gamma_%s",tail[i].c_str()));
      tpc_gamma1d.push_back(new TH1F(tpc_gamma1d_name.back(),"",bins[i],low[i],high[i]));
      tpc_gamma1d.at(i)->GetXaxis()->SetTitle(unit[i].c_str());
      list->Add(tpc_gamma1d.at(i));
      
      tpc_neutron1d_name.push_back(Form("tpc_neutron_%s",tail[i].c_str()));
      tpc_neutron1d.push_back(new TH1F(tpc_neutron1d_name.back(),"",bins[i],low[i],high[i]));
      tpc_neutron1d.at(i)->GetXaxis()->SetTitle(unit[i].c_str());
      list->Add(tpc_neutron1d.at(i));
      
      nv_gamma1d_name.push_back(Form("nv_gamma_%s",tail[i].c_str()));
      nv_gamma1d.push_back(new TH1F(nv_gamma1d_name.back(),"",bins[i],low[i],high[i]));
      tpc_gamma1d.at(i)->GetXaxis()->SetTitle(unit[i].c_str());
      list->Add(nv_gamma1d.at(i));
      
      nv_neutron1d_name.push_back(Form("nv_neutron_%s",tail[i].c_str()));
      nv_neutron1d.push_back(new TH1F(nv_neutron1d_name.back(),"",bins[i],low[i],high[i]));
      tpc_neutron1d.at(i)->GetXaxis()->SetTitle(unit[i].c_str());
      list->Add(nv_neutron1d.at(i));            
    }

      tpc_neutron1d.at(0)->Fill(quenchingfactor->at(j));
      tpc_neutron1d.at(1)->Fill(ex->at(j));
      tpc_neutron1d.at(2)->Fill(ey->at(j));
      tpc_neutron1d.at(3)->Fill(ez->at(j));
      tpc_neutron1d.at(4)->Fill(er);
      tpc_neutron1d.at(5)->Fill(theta);
      tpc_neutron1d.at(6)->Fill(eqch->at(j));
      tpc_neutron1d.at(7)->Fill(edep->at(j));
      tpc_neutron1d.at(8)->Fill(et->at(j)*1.e+9);
	

      tpc_gamma1d.at(0)->Fill(quenchingfactor->at(j));
      tpc_gamma1d.at(1)->Fill(ex->at(j));
      tpc_gamma1d.at(2)->Fill(ey->at(j));
      tpc_gamma1d.at(3)->Fill(ez->at(j));
      tpc_gamma1d.at(4)->Fill(er);
      tpc_gamma1d.at(5)->Fill(theta);
      tpc_gamma1d.at(6)->Fill(eqch->at(j));
      tpc_gamma1d.at(7)->Fill(edep->at(j));
      tpc_gamma1d.at(8)->Fill(et->at(j)*1.e+9);

      nv_neutron1d.at(0)->Fill(quenchingfactor->at(j));
      nv_neutron1d.at(1)->Fill(ex->at(j));
      nv_neutron1d.at(2)->Fill(ey->at(j));
      nv_neutron1d.at(3)->Fill(ez->at(j));
      nv_neutron1d.at(4)->Fill(er);
      nv_neutron1d.at(5)->Fill(theta);
      nv_neutron1d.at(6)->Fill(eqch->at(j));
      nv_neutron1d.at(7)->Fill(edep->at(j));
      nv_neutron1d.at(8)->Fill(et->at(j)*1.e+9);
      
      nv_gamma1d.at(0)->Fill(quenchingfactor->at(j));
      nv_gamma1d.at(1)->Fill(ex->at(j));
      nv_gamma1d.at(2)->Fill(ey->at(j));
      nv_gamma1d.at(3)->Fill(ez->at(j));
      nv_gamma1d.at(4)->Fill(er);
      nv_gamma1d.at(5)->Fill(theta);
      nv_gamma1d.at(6)->Fill(eqch->at(j));
      nv_gamma1d.at(7)->Fill(edep->at(j));
      nv_gamma1d.at(8)->Fill(et->at(j)*1.e+9);
	  if(quenchingfactor->at(j)<quenchingcut)
	    {
	      ntuple_plots.at(4)->Fill(quenchingfactor->at(j),ex->at(j),ey->at(j),ez->at(j),er,theta,eqch->at(j),edep->at(j),et->at(j)*1.e+9);
	      eqch_hist.at(4)->Fill(eqch->at(j));
	      eqch_time_hist.at(4)->Fill(eqch->at(j),et->at(j)*1.e+9);
	      eqch_sourceene_hist.at(4)->Fill(eqch->at(j),source_ene*1000.);
	      eqch_sourcetheta_hist.at(4)->Fill(eqch->at(j),sourcetheta);
	      eqch_theta_hist.at(4)->Fill(eqch->at(j),theta);
	      eqch_radius_hist.at(4)->Fill(eqch->at(j),er);

	    }
	  if(quenchingfactor->at(j)>quenchingcut)
	    {
	      ntuple_plots.at(5)->Fill(quenchingfactor->at(j),ex->at(j),ey->at(j),ez->at(j),er,theta,eqch->at(j),edep->at(j),et->at(j)*1.e+9);
	      eqch_hist.at(5)->Fill(eqch->at(j));
	      eqch_time_hist.at(5)->Fill(eqch->at(j),et->at(j)*1.e+9);
	      eqch_sourceene_hist.at(5)->Fill(eqch->at(j),source_ene*1000.);
	      eqch_sourcetheta_hist.at(5)->Fill(eqch->at(j),sourcetheta);
	      eqch_theta_hist.at(5)->Fill(eqch->at(j),theta);
	      eqch_radius_hist.at(5)->Fill(eqch->at(j),er);

	    }
      
  for(size_t j=0; j<tpc_gamma1d_name.size(); j++)
    {
      tpc_gamma1d.at(j) = dynamic_cast<TH1F*>(list->FindObject(Form(tpc_gamma1d_name.at(j))));
    }
  for(size_t j=0; j<tpc_neutron1d_name.size(); j++)
    {
      tpc_neutron1d.at(j) = dynamic_cast<TH1F*>(list->FindObject(Form(tpc_neutron1d_name.at(j))));
    }
  for(size_t j=0; j<nv_gamma1d_name.size(); j++)
    {
      nv_gamma1d.at(j) = dynamic_cast<TH1F*>(list->FindObject(Form(nv_gamma1d_name.at(j))));
    }
  for(size_t j=0; j<nv_neutron1d_name.size(); j++)
    {
      nv_neutron1d.at(j) = dynamic_cast<TH1F*>(list->FindObject(Form(nv_neutron1d_name.at(j))));
    }
  for(size_t j=0; j<tpc_gamma1d.size(); j++)
    tpc_gamma1d.at(j)->Write();
  for(size_t j=0; j<tpc_neutron1d.size(); j++)
    tpc_neutron1d.at(j)->Write();
  for(size_t j=0; j<nv_gamma1d.size(); j++)
    nv_gamma1d.at(j)->Write();
  for(size_t j=0; j<nv_neutron1d.size(); j++)
    nv_neutron1d.at(j)->Write();
	 
	  if(quenchingfactor->at(j)<quenchingcut)
	    {
	      ntuple_plots.at(7+k*NUM)->Fill(quenchingfactor->at(j),ex->at(j),ey->at(j),ez->at(j),er,theta,eqch->at(j),edep->at(j),et->at(j)*1.e+9);
	      eqch_hist.at(7+k*NUM)->Fill(eqch->at(j));
	      eqch_time_hist.at(7+k*NUM)->Fill(eqch->at(j),et->at(j)*1.e+9);
	      eqch_sourceene_hist.at(7+k*NUM)->Fill(eqch->at(j),source_ene*1000.);
	      eqch_sourcetheta_hist.at(7+k*NUM)->Fill(eqch->at(j),sourcetheta);
	      eqch_theta_hist.at(7+k*NUM)->Fill(eqch->at(j),theta);
	      eqch_radius_hist.at(7+k*NUM)->Fill(eqch->at(j),er);

	    }
	  if(quenchingfactor->at(j)>quenchingcut)
	    {
	      ntuple_plots.at(8+k*NUM)->Fill(quenchingfactor->at(j),ex->at(j),ey->at(j),ez->at(j),er,theta,eqch->at(j),edep->at(j),et->at(j)*1.e+9);
	      eqch_hist.at(8+k*NUM)->Fill(eqch->at(j));
	      eqch_time_hist.at(8+k*NUM)->Fill(eqch->at(j),et->at(j)*1.e+9);
	      eqch_sourceene_hist.at(8+k*NUM)->Fill(eqch->at(j),source_ene*1000.);
	      eqch_sourcetheta_hist.at(8+k*NUM)->Fill(eqch->at(j),sourcetheta);
	      eqch_theta_hist.at(8+k*NUM)->Fill(eqch->at(j),theta);
	      eqch_radius_hist.at(8+k*NUM)->Fill(eqch->at(j),er);

	    }

  */
