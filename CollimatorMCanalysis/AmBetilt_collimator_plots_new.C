#include "TMath.h"
#include <string>
#include <sstream>
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
using namespace std;

const string Time = "_Sep27AM";
const int datafiles=7;
int Volume = 0;

void AmBetilt_collimator_plots_small(){
  string label;
  if(Volume == 9) label="collimator";
  else label="no_collimator";

  TChain *t=new TChain("Recon");
  Readdatafile(t,Volume);
  // t->Add("/darkside/users/hqian/neutron0G163tilt/cluster_data/outnvetoambe9L_test_cylinder_clustered.root");
  //    t->Add("/darkside/users/hqian/neutron0G163tilt/cluster_data/outnvetoambe9L_v7_cylinder_clustered.root");
  //data for BG
  //  t->Add("/home/hqian/montecarlo/g4ds10/Linux-g++/neutron0G163tilt/cluster_data/outnvetoambe0L_v5_cylinder_clustered.root");
  MakePlots(t,label);
  MakeCounts(t,label);
}

void Readdatafile(TChain *t, int volume)
{   string dirname="/darkside/users/hqian/neutron0G163tilt/cluster_data/";
  //   string dirname="/ds50/data/user/hqian36/Collimator/neutron0G163tilt/cluster_data/";
    string filename;
    string middle="outnvetoambe";
    string last="_cylinder_clustered.root";
    stringstream oss;
    stringstream oss2;
    for(int i=0; i<datafiles; i++)
    {
        if(i==0)
        {
	  oss<<volume;
	  filename=dirname+middle+oss.str()+"L"+last;
	  //	  filename<<dirname<<middle<<volume<<"L"<<last;
	  //sprintf(filename,"%s%s%iL%s",dirname.c_str(),middle.c_str(),volume,last.c_str());
	  t->Add(filename.c_str());
	  cout<<"Processing Data file: "<<filename<<endl;
	  oss.str("");
        }
        else
	  { 
	    oss<<volume;
	    oss2<<i;
	    filename=dirname+middle+oss.str()+"L_v"+oss2.str()+last;
            t->Add(filename.c_str());	
	    cout<<"Processing Data file: "<<filename<<endl;
	    oss.str("");
	    oss2.str("");
	  }
    }
}

void MakePlots(TChain *t, string label){

  TCut nuclearenergy="edep_nuclear>0.025";
  //  TCut tpcactive="ez<14.6 && ez>-22.2 && sqrt(ex**2+ey**2)<17.77";
  TCut tpcactive="volume==\"p_active\"";
  TCut yslice="ey>-0.5 && ey<0.5";
  TCut gap="ex>-100 && ex<0";
  TCut quchcutn="quenchingfactor<0.43";
  TCut quchcute="quenchingfactor>0.43";
  
  string quchf = " quenchingfactor<0.43 ";
  string quche = " quenchingfactor>0.43 ";
  
  string zslice1 = " 2.2<dep_z<4.6 ";
  string zslice2 = " 0<dep_z<2.2 ";
  string zslice3 = " -2.2<dep_z<0 ";

  //without quenching cuts
  string Hist2DTitle = "dep_y vs dep_x in TPC active volume "+ label;
  TH2F* Hist2D = new TH2F("Hist2D",Hist2DTitle.c_str(), 100,-20,20,100,-20,20);
  Hist2D->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2D->GetYaxis()->SetTitle("dep_y [cm]");

  string Hist2DXYslice1Title = "dep_y vs dep_x in TPC active volume"+zslice1+ label;
  TH2F* Hist2DXYslice1 = new TH2F("Hist2DXYslice1",Hist2DXYslice1Title.c_str(),100,-20,20,100,-20,20);
  Hist2DXYslice1->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2DXYslice1->GetYaxis()->SetTitle("dep_y [cm]");

  string Hist2DXYslice2Title = "dep_y vs dep_x in TPC active volume"+zslice2+ label;
  TH2F* Hist2DXYslice2 = new TH2F("Hist2DXYslice2",Hist2DXYslice2Title.c_str(),100,-20,20,100,-20,20);
  Hist2DXYslice2->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2DXYslice2->GetYaxis()->SetTitle("dep_y [cm]");

  string Hist2DXYslice3Title = "dep_y vs dep_x in TPC active volume"+zslice3+ label;
  TH2F* Hist2DXYslice3 = new TH2F("Hist2DXYslice3",Hist2DXYslice3Title.c_str(),100,-20,20,100,-20,20);
  Hist2DXYslice3->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2DXYslice3->GetYaxis()->SetTitle("dep_y [cm]");

  string Hist2DGapTitle = "dep_z vs dep_x with y slice "+ label;
  TH2F* Hist2DGap = new TH2F("Hist2DGap",Hist2DGapTitle.c_str(), 150,-100,0,150,-50,100);
  Hist2DGap->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2DGap->GetYaxis()->SetTitle("dep_z [cm]");

  string QuenchingFactorTitle = "quenchingfactor in TPC active volume "+ label;
  TH1F* QuenchingFactor = new TH1F("QuenchingFactor",QuenchingFactorTitle.c_str(),100,0,1);
  QuenchingFactor->GetXaxis()->SetTitle("quenching factor");
  QuenchingFactor->GetYaxis()->SetTitle("counts");

  //quecnhing cuts for nuclear-like recoil
  string Hist2DnuclearTitle = "dep_y vs dep_x in TPC active volume"+quchf+ label;
  TH2F* Hist2Dnuclear = new TH2F("Hist2Dnuclear",Hist2DnuclearTitle.c_str(), 100,-20,20,100,-20,20);
  Hist2Dnuclear->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2Dnuclear->GetYaxis()->SetTitle("dep_y [cm]");

  string Hist2DXYnuclearslice1Title = "dep_y vs dep_x in TPC active volume"+zslice1+quchf+ label;
  TH2F* Hist2DXYnuclearslice1 = new TH2F("Hist2DXYnuclearslice1",Hist2DXYnuclearslice1Title.c_str(),100,-20,20,100,-20,20);
  Hist2DXYnuclearslice1->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2DXYnuclearslice1->GetYaxis()->SetTitle("dep_y [cm]");

  string Hist2DXYnuclearslice2Title = "dep_y vs dep_x in TPC active volume"+zslice2+quchf+ label;
  TH2F* Hist2DXYnuclearslice2 = new TH2F("Hist2DXYnuclearslice2",Hist2DXYnuclearslice2Title.c_str(),100,-20,20,100,-20,20);
  Hist2DXYnuclearslice2->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2DXYnuclearslice2->GetYaxis()->SetTitle("dep_y [cm]");

  string Hist2DXYnuclearslice3Title = "dep_y vs dep_x in TPC active volume"+zslice3+quchf+ label;
  TH2F* Hist2DXYnuclearslice3 = new TH2F("Hist2DXYnuclearslice3",Hist2DXYnuclearslice3Title.c_str(),100,-20,20,100,-20,20);
  Hist2DXYnuclearslice3->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2DXYnuclearslice3->GetYaxis()->SetTitle("dep_y [cm]");

  string Hist1DNRcountsslice1Title = "NR counts in slice1 in TPC active volume "+ label;
  TH1F* Hist1DNRcountsslice1 = new TH1F("Hist1DNRcountsslice1",Hist1DNRcountsslice1Title.c_str(),100,0,20);
  Hist1DNRcountsslice1->GetXaxis()->SetTitle("NR counts");

  string Hist1DNRcountsslice2Title = "NR counts in slice2 in TPC active volume "+ label;
  TH1F* Hist1DNRcountsslice2 = new TH1F("Hist1DNRcountsslice2",Hist1DNRcountsslice2Title.c_str(),100,0,20);
  Hist1DNRcountsslice2->GetXaxis()->SetTitle("NR counts");

  string Hist1DNRcountsslice3Title = "NR counts in slice3 in TPC active volume "+ label;
  TH1F* Hist1DNRcountsslice3 = new TH1F("Hist1DNRcountsslice3",Hist1DNRcountsslice3Title.c_str(),100,0,20);
  Hist1DNRcountsslice3->GetXaxis()->SetTitle("NR counts");

  string Hist1DNRcountsslice4Title = "NR counts in slice4 in TPC active volume "+ label;
  TH1F* Hist1DNRcountsslice4 = new TH1F("Hist1DNRcountsslice4",Hist1DNRcountsslice4Title.c_str(),100,0,20);
  Hist1DNRcountsslice4->GetXaxis()->SetTitle("NR counts");

  //single nuclear sacttering in TPC active volume
  string Hist2DsinglenuclearTitle = "dep_y vs dep_x for single nuclear scattering in TPC active volume"+quchf+ label;
  TH2F* Hist2Dsinglenuclear = new TH2F("Hist2Dsinglenuclear",Hist2DsinglenuclearTitle.c_str(), 100,-20,20,100,-20,20);
  Hist2Dsinglenuclear->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2Dsinglenuclear->GetYaxis()->SetTitle("dep_y [cm]");

  string Hist2DYZsinglenuclearTitle = "dep_z vs dep_y for single nuclear scattering in TPC active volume"+quchf+ label;
  TH2F* Hist2DYZsinglenuclear = new TH2F("Hist2DYZsinglenuclear",Hist2DYZsinglenuclearTitle.c_str(),100,-20,20,100,-25,20);
  Hist2DYZsinglenuclear->GetXaxis()->SetTitle("dep_y [cm]");
  Hist2DYZsinglenuclear->GetYaxis()->SetTitle("dep_z [cm]");
  
  string Hist2DXYsinglenuclearslice1Title = "dep_y vs dep_x for single nuclear scattering in TPC active volume"+zslice1+quchf+ label;
  TH2F* Hist2DXYsinglenuclearslice1 = new TH2F("Hist2DXYsinglenuclearslice1",Hist2DXYsinglenuclearslice1Title.c_str(),100,-20,20,100,-20,20);
  Hist2DXYsinglenuclearslice1->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2DXYsinglenuclearslice1->GetYaxis()->SetTitle("dep_y [cm]");

  string Hist2DXYsinglenuclearslice2Title = "dep_y vs dep_x for single nuclear scattering in TPC active volume"+zslice2+quchf+ label;
  TH2F* Hist2DXYsinglenuclearslice2 = new TH2F("Hist2DXYsinglenuclearslice2",Hist2DXYsinglenuclearslice2Title.c_str(),100,-20,20,100,-20,20);
  Hist2DXYsinglenuclearslice2->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2DXYsinglenuclearslice2->GetYaxis()->SetTitle("dep_y [cm]");
    
  string Hist2DXYsinglenuclearslice3Title = "dep_y vs dep_x for single nuclear scattering in TPC active volume"+zslice3+quchf+ label;
  TH2F* Hist2DXYsinglenuclearslice3 = new TH2F("Hist2DXYsinglenuclearslice3",Hist2DXYsinglenuclearslice3Title.c_str(),100,-20,20,100,-20,20);
  Hist2DXYsinglenuclearslice3->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2DXYsinglenuclearslice3->GetYaxis()->SetTitle("dep_y [cm]");

  //single sacttering in TPC active volume
  string Hist2DsingleTitle = "dep_y vs dep_x for single scattering in TPC active volume "+ label;
  TH2F* Hist2Dsingle = new TH2F("Hist2Dsingle",Hist2DsingleTitle.c_str(), 100,-20,20,100,-20,20);
  Hist2Dsingle->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2Dsingle->GetYaxis()->SetTitle("dep_y [cm]");

  string Hist2DYZsingleTitle = "dep_z vs dep_y for single scattering in TPC active volume "+ label;
  TH2F* Hist2DYZsingle = new TH2F("Hist2DYZsingle",Hist2DYZsingleTitle.c_str(),100,-20,20,100,-25,20);
  Hist2DYZsingle->GetXaxis()->SetTitle("dep_y [cm]");
  Hist2DYZsingle->GetYaxis()->SetTitle("dep_z [cm]");

  string Hist2DXYsingleslice1Title = "dep_y vs dep_x for single scattering in TPC active volume"+zslice1+ label;
  TH2F* Hist2DXYsingleslice1 = new TH2F("Hist2DXYsingleslice1",Hist2DXYsingleslice1Title.c_str(),100,-20,20,100,-20,20);
  Hist2DXYsingleslice1->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2DXYsingleslice1->GetYaxis()->SetTitle("dep_y [cm]");

  string Hist2DXYsingleslice2Title = "dep_y vs dep_x for single scattering in TPC active volume"+zslice2+ label;
  TH2F* Hist2DXYsingleslice2 = new TH2F("Hist2DXYsingleslice2",Hist2DXYsingleslice2Title.c_str(),100,-20,20,100,-20,20);
  Hist2DXYsingleslice2->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2DXYsingleslice2->GetYaxis()->SetTitle("dep_y [cm]");
    
  string Hist2DXYsingleslice3Title = "dep_y vs dep_x for single scattering in TPC active volume"+zslice3+ label;
  TH2F* Hist2DXYsingleslice3 = new TH2F("Hist2DXYsingleslice3",Hist2DXYsingleslice3Title.c_str(),100,-20,20,100,-20,20);
  Hist2DXYsingleslice3->GetXaxis()->SetTitle("dep_x [cm]");
  Hist2DXYsingleslice3->GetYaxis()->SetTitle("dep_y [cm]");

  int nEntries = t->GetEntries();
  cout << "Entries= " << nEntries << endl;
  vector<double> *ex=0;
  vector<double> *ey=0;
  vector<double> *ez=0;
  vector<double> *et=0;
  vector<double> *edep=0;
  vector<double> *edep_nuclear=0;
  vector<double> *edep_electron=0;
  vector<double> *eqch=0;
  vector<double> *quenchingfactor=0;
  vector<TString>* volume=0;
  vector<double> nrx=0;
  vector<double> nry=0;
  vector<double> nrz=0;
  vector<double> rx=0;
  vector<double> ry=0;
  vector<double> rz=0;
  
  t->SetBranchAddress("volume", &volume);
  t->SetBranchAddress("et", &et);
  t->SetBranchAddress("ex", &ex);
  t->SetBranchAddress("ey", &ey);
  t->SetBranchAddress("ez", &ez);
  t->SetBranchAddress("edep", &edep);
  t->SetBranchAddress("edep_nuclear", &edep_nuclear);
  t->SetBranchAddress("edep_electron", &edep_electron);
  t->SetBranchAddress("eqch", &eqch);
  t->SetBranchAddress("quenchingfactor", &quenchingfactor);
  
  double quenchingcut = 0.43;
  double dep_enecut= 23.0; //keV
  double dep_enecut2=269.0; //keV

  int eff=0;   
  bool IsNR,IsER;
  int NNR, NER, ndeps, nsing, NRcounts;;
  NNR=NER=0;
  nsing=0;
  for(int j=0; j<nEntries; j++)
    { 
      if(!(j%100000)) std::cout<<"Processing Event "<<j<<", " <<Int_t(100.*j/nEntries)<<"% Completed"<<std::endl;
      ndeps=0;
      NRcounts=0;
      t->GetEntry(j);
      IsNR = false;
      IsER = false;
      nrx.clear();
      nry.clear();
      nrz.clear();
      rx.clear();
      ry.clear();
      rz.clear();
      for(int i=0; i< et->size(); i++)
        {
	  if(volume->at(i) =="p_active")
            { ndeps++;
	      IsER=true;
	      Hist2D->Fill(ex->at(i),ey->at(i));
	      QuenchingFactor->Fill(quenchingfactor->at(i));
	      if(ndeps==1)
		{ eff++;
		  rx.push_back(ex->at(i));
		  ry.push_back(ey->at(i));
		  rz.push_back(ez->at(i));
		}           
	      if(ez->at(i)<4.6 && ez->at(i)>2.2)
		Hist2DXYslice1->Fill(ex->at(i),ey->at(i));
	      if(ez->at(i)<2.2 && ez->at(i)>0)
		Hist2DXYslice2->Fill(ex->at(i),ey->at(i));
	      if(ez->at(i)>-2.2 && ez->at(i)<0)
		Hist2DXYslice3->Fill(ex->at(i),ey->at(i));
	      if(quenchingfactor->at(i) < quenchingcut && edep->at(i)<dep_enecut2 && edep->at(i)>dep_enecut )
                {
		  //  NRcounts++;
		  IsNR = true;
		  Hist2Dnuclear->Fill(ex->at(i),ey->at(i));
		  nrx.push_back(ex->at(i));
		  nry.push_back(ey->at(i));
		  nrz.push_back(ez->at(i));
		  if(ez->at(i)>4.6)
		  {
		    NRcounts++;
		    Hist1DNRcountsslice4->Fill(NRcounts);		    		   
		  } 
		 else if(ez->at(i)<4.6 && ez->at(i)>2.2)
		    {  Hist2DXYnuclearslice1->Fill(ex->at(i),ey->at(i));
		      NRcounts++;
		      Hist1DNRcountsslice1->Fill(NRcounts);
		    }
		 else if(ez->at(i)<2.2 && ez->at(i)>0)
		    { Hist2DXYnuclearslice2->Fill(ex->at(i),ey->at(i));	
		      NRcounts++;
		      Hist1DNRcountsslice2->Fill(NRcounts);
		    }
		 else if(ez->at(i)>-2.2 && ez->at(i)<0)
		    {  Hist2DXYnuclearslice3->Fill(ex->at(i),ey->at(i));	
		      NRcounts++;
		      Hist1DNRcountsslice3->Fill(NRcounts);		    
		    }
		}
	    }

	  if(ex->at(i)>-100 && ex->at(i)<0 && ey->at(i)<0.5 && ey->at(i)>-0.5)
            {
	      Hist2DGap->Fill(ex->at(i),ez->at(i));
	    }
        }

      if(IsNR)
	NNR++;
      if(IsER)
	NER++;
      if(ndeps==1)
	{ nsing++;
	  for(int a=0; a<ndeps; a++)
	    {
	      Hist2Dsingle->Fill(rx.at(a),ry.at(a));
	      Hist2DYZsingle->Fill(ry.at(a),rz.at(a));
	      if(rz.at(a)<4.6 && rz.at(a)>2.2)
		Hist2DXYsingleslice1->Fill(rx.at(a),ry.at(a));
	      if(rz.at(a)<2.2 && rz.at(a)>0)
		Hist2DXYsingleslice2->Fill(rx.at(a),ry.at(a));	  
	      if(rz.at(a)>-2.2 && rz.at(a)<0)
		Hist2DXYsingleslice3->Fill(rx.at(a),ry.at(a));	  
	    }
	}

      if(ndeps==1 && IsNR )
	{ 
	  for(int h=0; h<ndeps; h++)
	    {
	      Hist2Dsinglenuclear->Fill(nrx.at(h),nry.at(h));
	      Hist2DYZsinglenuclear->Fill(nry.at(h),nrz.at(h));
	      if(nrz.at(h)<4.6 && nrz.at(h)>2.2)
		Hist2DXYsinglenuclearslice1->Fill(nrx.at(h),nry.at(h));
	      if(nrz.at(h)<2.2 && nrz.at(h)>0)
		Hist2DXYsinglenuclearslice2->Fill(nrx.at(h),nry.at(h));	  
	      if(nrz.at(h)>-2.2 && nrz.at(h)<0)
		Hist2DXYsinglenuclearslice3->Fill(nrx.at(h),nry.at(h));	  
	    }
	}
    }
    
  cout<<"NNR= " <<NNR <<"     "<<"NER= "<<NER<<"  "<<"nsing=  "<<nsing<<endl;
  cout<<"eff= "<<eff<<endl;
  float leftend = -33.5;
  int leftendbin = Hist2DGap->GetXaxis()->FindBin(leftend);
  //    int leftendbin = 102;
  TH1D* py1 = Hist2DGap->ProjectionY("py1",leftendbin,leftendbin);
  
  int slicebin = Hist2DXYslice1->GetNbinsX();
  TH1D* Hist2DXYslice1ProjectionY = Hist2DXYslice1->ProjectionY("Hist2DXYslice1ProjectionY",0,slicebin);
  TH1D* Hist2DXYnuclearslice1ProjectionY = Hist2DXYnuclearslice1->ProjectionY("Hist2DXYnuclearslice1ProjectionY",0,slicebin);
  TH1D* Hist2DXYsingleslice1ProjectionY = Hist2DXYsingleslice1->ProjectionY("Hist2DXYsingleslice1ProjectionY",0,slicebin);
  TH1D* Hist2DXYsinglenuclearslice1ProjectionY = Hist2DXYsinglenuclearslice1->ProjectionY("Hist2DXYsinglenuclearslice1ProjectionY",0,slicebin);
    
  TH1D* Hist2DXYslice2ProjectionY = Hist2DXYslice2->ProjectionY("Hist2DXYslice2ProjectionY",0,slicebin);
  TH1D* Hist2DXYnuclearslice2ProjectionY = Hist2DXYnuclearslice2->ProjectionY("Hist2DXYnuclearslice2ProjectionY",0,slicebin);
  TH1D* Hist2DXYsingleslice2ProjectionY = Hist2DXYsingleslice2->ProjectionY("Hist2DXYsingleslice2ProjectionY",0,slicebin);
  TH1D* Hist2DXYsinglenuclearslice2ProjectionY = Hist2DXYsinglenuclearslice2->ProjectionY("Hist2DXYsinglenuclearslice2ProjectionY",0,slicebin);

  TH1D* Hist2DXYslice3ProjectionY = Hist2DXYslice3->ProjectionY("Hist2DXYslice3ProjectionY",0,slicebin);
  TH1D* Hist2DXYnuclearslice3ProjectionY = Hist2DXYnuclearslice3->ProjectionY("Hist2DXYnuclearslice3ProjectionY",0,slicebin);
  TH1D* Hist2DXYsingleslice3ProjectionY = Hist2DXYsingleslice3->ProjectionY("Hist2DXYsingleslice3ProjectionY",0,slicebin);
  TH1D* Hist2DXYsinglenuclearslice3ProjectionY = Hist2DXYsinglenuclearslice3->ProjectionY("Hist2DXYsinglenuclearslice3ProjectionY",0,slicebin);

  int slicebinX1 = Hist2DXYslice1->GetYaxis()->FindBin(-2);
  int slicebinX2 = Hist2DXYslice1->GetYaxis()->FindBin(2);
  TH1D* Hist2DXYslice1ProjectionX = Hist2DXYslice1->ProjectionX("Hist2DXYslice1ProjectionX",slicebinX1,slicebinX2);
  TH1D* Hist2DXYnuclearslice1ProjectionX = Hist2DXYnuclearslice1->ProjectionX("Hist2DXYnuclearslice1ProjectionX",slicebinX1,slicebinX2);
  TH1D* Hist2DXYsingleslice1ProjectionX = Hist2DXYsingleslice1->ProjectionX("Hist2DXYsingleslice1ProjectionX",slicebinX1,slicebinX2);
  TH1D* Hist2DXYsinglenuclearslice1ProjectionX = Hist2DXYsinglenuclearslice1->ProjectionX("Hist2DXYsinglenuclearslice1ProjectionX",slicebinX1,slicebinX2);
    
  TH1D* Hist2DXYslice2ProjectionX = Hist2DXYslice2->ProjectionX("Hist2DXYslice2ProjectionX",slicebinX1,slicebinX2);
  TH1D* Hist2DXYnuclearslice2ProjectionX = Hist2DXYnuclearslice2->ProjectionX("Hist2DXYnuclearslice2ProjectionX",slicebinX1,slicebinX2);
  TH1D* Hist2DXYsingleslice2ProjectionX = Hist2DXYsingleslice2->ProjectionX("Hist2DXYsingleslice2ProjectionX",slicebinX1,slicebinX2);
  TH1D* Hist2DXYsinglenuclearslice2ProjectionX = Hist2DXYsinglenuclearslice2->ProjectionX("Hist2DXYsinglenuclearslice2ProjectionX",slicebinX1,slicebinX2);

  TH1D* Hist2DXYslice3ProjectionX = Hist2DXYslice3->ProjectionX("Hist2DXYslice3ProjectionX",slicebinX1,slicebinX2);
  TH1D* Hist2DXYnuclearslice3ProjectionX = Hist2DXYnuclearslice3->ProjectionX("Hist2DXYnuclearslice3ProjectionX",slicebinX1,slicebinX2);
  TH1D* Hist2DXYsingleslice3ProjectionX = Hist2DXYsingleslice3->ProjectionX("Hist2DXYsingleslice3ProjectionX",slicebinX1,slicebinX2);
  TH1D* Hist2DXYsinglenuclearslice3ProjectionX = Hist2DXYsinglenuclearslice3->ProjectionX("Hist2DXYsinglenuclearslice3ProjectionX",slicebinX1,slicebinX2);
    
  TCanvas *c1=new TCanvas("c1",label.c_str(),1000,600);
  Hist2D->Draw("colz");
  c1->SaveAs(Form("ambe_Hist2D_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYslice1->Draw("colz");
  c1->SaveAs(Form("ambe_Hist2DXYslice1_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYslice2->Draw("colz");
  c1->SaveAs(Form("ambe_Hist2DXYslice2_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYslice3->Draw("colz");
  c1->SaveAs(Form("ambe_Hist2DXYslice3_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYslice1ProjectionY->Draw();
  c1->SaveAs(Form("ambe_Hist2DXYslice1ProjectionY_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYslice2ProjectionY->Draw();
  c1->SaveAs(Form("ambe_Hist2DXYslice2ProjectionY_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYslice3ProjectionY->Draw();
  c1->SaveAs(Form("ambe_Hist2DXYslice3ProjectionY_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYslice1ProjectionX->Draw();
  c1->SaveAs(Form("ambe_Hist2DXYslice1ProjectionX_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYslice2ProjectionX->Draw();
  c1->SaveAs(Form("ambe_Hist2DXYslice2ProjectionX_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYslice3ProjectionX->Draw();
  c1->SaveAs(Form("ambe_Hist2DXYslice3ProjectionX_%s%s.png",label.c_str(),Time.c_str()));
  py1->Draw();
  c1->SaveAs(Form("ambe_Hist2DGapProjectionY_%s%s.png",label.c_str(),Time.c_str()));
  c1->SetLogz();
  Hist2DGap->Draw("colz");
  c1->SaveAs(Form("ambe_Hist2DGap_%s%s.png",label.c_str(),Time.c_str()));
    
  TCanvas *c2=new TCanvas("c2",label.c_str(),1000,600);
  Hist2Dnuclear->Draw("colz");
  c2->SaveAs(Form("ambe_Hist2Dnuclear_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYnuclearslice1->Draw("colz");
  c2->SaveAs(Form("ambe_Hist2DXYnuclearslice1_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYnuclearslice2->Draw("colz");
  c2->SaveAs(Form("ambe_Hist2DXYnuclearslice2_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYnuclearslice3->Draw("colz");
  c2->SaveAs(Form("ambe_Hist2DXYnuclearslice3_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYnuclearslice1ProjectionY->Draw();
  c2->SaveAs(Form("ambe_Hist2DXYnuclearslice1ProjectionY_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYnuclearslice2ProjectionY->Draw();
  c2->SaveAs(Form("ambe_Hist2DXYnuclearslice2ProjectionY_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYnuclearslice3ProjectionY->Draw();
  c2->SaveAs(Form("ambe_Hist2DXYnuclearslice3ProjectionY_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYnuclearslice1ProjectionX->Draw();
  c2->SaveAs(Form("ambe_Hist2DXYnuclearslice1ProjectionX_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYnuclearslice2ProjectionX->Draw();
  c2->SaveAs(Form("ambe_Hist2DXYnuclearslice2ProjectionX_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYnuclearslice3ProjectionX->Draw();
  c2->SaveAs(Form("ambe_Hist2DXYnuclearslice3ProjectionX_%s%s.png",label.c_str(),Time.c_str()));
    
  TCanvas *c4=new TCanvas("c4",label.c_str(),1000,600);
  Hist2Dsinglenuclear->Draw("colz");
  c4->SaveAs(Form("ambe_Hist2Dsinglenuclear_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYsinglenuclearslice1->Draw("colz");
  c4->SaveAs(Form("ambe_Hist2DXYsinglenuclearslice1_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYsinglenuclearslice2->Draw("colz");
  c4->SaveAs(Form("ambe_Hist2DXYsinglenuclearslice2_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYsinglenuclearslice3->Draw("colz");
  c4->SaveAs(Form("ambe_Hist2DXYsinglenuclearslice3_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYsinglenuclearslice1ProjectionY->Draw();
  c4->SaveAs(Form("ambe_Hist2DXYsinglenuclearslice1ProjectionY_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYsinglenuclearslice2ProjectionY->Draw();
  c4->SaveAs(Form("ambe_Hist2DXYsinglenuclearslice2ProjectionY_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYsinglenuclearslice3ProjectionY->Draw();
  c4->SaveAs(Form("ambe_Hist2DXYsinglenuclearslice3ProjectionY_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYsinglenuclearslice1ProjectionX->Draw();
  c4->SaveAs(Form("ambe_Hist2DXYsinglenuclearslice1ProjectionX_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYsinglenuclearslice2ProjectionX->Draw();
  c4->SaveAs(Form("ambe_Hist2DXYsinglenuclearslice2ProjectionX_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYsinglenuclearslice3ProjectionX->Draw();
  c4->SaveAs(Form("ambe_Hist2DXYsinglenuclearslice3ProjectionX_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DYZsinglenuclear->Draw("colz");
  c4->SaveAs(Form("ambe_Hist2DYZsinglenuclear_%s%s.png",label.c_str(),Time.c_str()));

  TCanvas *c5=new TCanvas("c5",label.c_str(),1000,600);
  Hist2Dsingle->Draw("colz");
  c5->SaveAs(Form("ambe_Hist2Dsingle_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYsingleslice1->Draw("colz");
  c5->SaveAs(Form("ambe_Hist2DXYsingleslice1_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYsingleslice2->Draw("colz");
  c5->SaveAs(Form("ambe_Hist2DXYsingleslice2_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYsingleslice3->Draw("colz");
  c5->SaveAs(Form("ambe_Hist2DXYsingleslice3_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYsingleslice1ProjectionY->Draw();
  c5->SaveAs(Form("ambe_Hist2DXYsingleslice1ProjectionY_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYsingleslice2ProjectionY->Draw();
  c5->SaveAs(Form("ambe_Hist2DXYsingleslice2ProjectionY_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYsingleslice3ProjectionY->Draw();
  c5->SaveAs(Form("ambe_Hist2DXYsingleslice3ProjectionY_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYsingleslice1ProjectionX->Draw();
  c5->SaveAs(Form("ambe_Hist2DXYsingleslice1ProjectionX_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYsingleslice2ProjectionX->Draw();
  c5->SaveAs(Form("ambe_Hist2DXYsingleslice2ProjectionX_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DXYsingleslice3ProjectionX->Draw();
  c5->SaveAs(Form("ambe_Hist2DXYsingleslice3ProjectionX_%s%s.png",label.c_str(),Time.c_str()));
  Hist2DYZsingle->Draw("colz");
  c5->SaveAs(Form("ambe_Hist2DYZsingle_%s%s.png",label.c_str(),Time.c_str()));

  /*  TCanvas *c1=new TCanvas("c1",label.c_str(),1000,600);
      c1->SetLogz();
  //  c1->SaveAs(Form("AmBe_yz_%s.png",label.c_str()));
  //without quenching cuts
  t->Draw("ey:ex>>Hist2D",tpcactive,"colz");
  //  c1->SaveAs(Form("ambe_Hist2D_%s%s.png",label.c_str(),Time.c_str()));
  t->Draw("ez:ey>>Hist2DYZ",tpcactive, "colz");
  //  c1->SaveAs(Form("ambe_Hist2DYZ_%s%s.png",label.c_str(),Time.c_str()));
  t->Draw("ez:ex>>Hist2DXZ",tpcactive, "colz");
  t->Draw("ey:ex>>Hist2DLAr",tpcactive && zslice, "colz");
  t->Draw("ey:ex>>Hist2DGap",gap && zslice, "colz");
  t->Draw("quenchingfactor>>QuenchingFactor",tpcactive);
  //quenching cuts for nuclear-like recoil
  t->Draw("ey:ex>>Hist2Dnuclear",tpcactive && quchcutn,"colz");
  t->Draw("ez:ey>>Hist2DYZnuclear",tpcactive && quchcutn, "colz");
  t->Draw("ez:ex>>Hist2DXZnuclear",tpcactive && quchcutn, "colz");
  t->Draw("ey:ex>>Hist2DLArnuclear",tpcactive && zslice && quchcutn, "colz");
  t->Draw("ey:ex>>Hist2DGapnuclear",quchcutn && gap && zslice, "colz");
  //quenching cuts for electron-like recoil
  t->Draw("ey:ex>>Hist2Delectron",tpcactive && quchcute,"colz");
  t->Draw("ez:ey>>Hist2DYZelectron",tpcactive && quchcute, "colz");
  t->Draw("ez:ex>>Hist2DXZelectron",tpcactive && quchcute, "colz");
  t->Draw("ey:ex>>Hist2DLArelectron",tpcactive && zslice && quchcute, "colz");
  t->Draw("ey:ex>>Hist2DGapelectron",quchcute && gap && zslice, "colz");
 */

  string outdirname="/darkside/users/hqian/neutron0G163tilt/cluster_data/";
  //  string outdirname="/ds50/data/user/hqian36/Collimator/neutron0G163tilt/cluster_data/";
       
  string outputname=outdirname+"nuclear0G163tilt_recoil_"+label+Time+"_clustered"+".root";  
  TFile f2D(outputname.c_str(), "RECREATE");
  Hist2D->Write(); 
  Hist2DXYslice1->Write();
  Hist2DXYslice2->Write();
  Hist2DXYslice3->Write();
  Hist2DXYslice1ProjectionY->Write();
  Hist2DXYslice2ProjectionY->Write();
  Hist2DXYslice3ProjectionY->Write();
  Hist2DXYslice1ProjectionX->Write();
  Hist2DXYslice2ProjectionX->Write();
  Hist2DXYslice3ProjectionX->Write();
  Hist2DGap->Write();
  py1->Write();
  QuenchingFactor->Write();
  Hist2Dnuclear->Write(); 
  Hist2DXYnuclearslice1->Write();
  Hist2DXYnuclearslice2->Write();
  Hist2DXYnuclearslice3->Write();
  Hist2DXYnuclearslice1ProjectionY->Write();
  Hist2DXYnuclearslice2ProjectionY->Write();
  Hist2DXYnuclearslice3ProjectionY->Write();
  Hist2DXYnuclearslice1ProjectionX->Write();
  Hist2DXYnuclearslice2ProjectionX->Write();
  Hist2DXYnuclearslice3ProjectionX->Write();
  Hist1DNRcountsslice1->Write();
  Hist1DNRcountsslice2->Write();
  Hist1DNRcountsslice3->Write();
  Hist1DNRcountsslice4->Write();
  Hist2Dsinglenuclear->Write();
  Hist2DYZsinglenuclear->Write();
  Hist2DXYsinglenuclearslice1->Write();
  Hist2DXYsinglenuclearslice2->Write();
  Hist2DXYsinglenuclearslice3->Write();
  Hist2DXYsinglenuclearslice1ProjectionY->Write();
  Hist2DXYsinglenuclearslice2ProjectionY->Write();
  Hist2DXYsinglenuclearslice3ProjectionY->Write();
  Hist2DXYsinglenuclearslice1ProjectionX->Write();
  Hist2DXYsinglenuclearslice2ProjectionX->Write();
  Hist2DXYsinglenuclearslice3ProjectionX->Write();
  Hist2Dsingle->Write();
  Hist2DYZsingle->Write();
  Hist2DXYsingleslice1->Write();
  Hist2DXYsingleslice2->Write();
  Hist2DXYsingleslice3->Write();
  Hist2DXYsingleslice1ProjectionY->Write();
  Hist2DXYsingleslice2ProjectionY->Write();
  Hist2DXYsingleslice3ProjectionY->Write();
  Hist2DXYsingleslice1ProjectionX->Write();
  Hist2DXYsingleslice2ProjectionX->Write();
  Hist2DXYsingleslice3ProjectionX->Write();
  f2D.Write();
  f2D.Close();
}

void MakeCounts(TChain *t, string label){
    
    int nEntries = t->GetEntries();
    cout << "Entries= " << nEntries << endl;
    vector<double> *ex=0;
    vector<double> *ey=0;
    vector<double> *ez=0;
    vector<double> *et=0;
    vector<double> *edep=0;
    vector<double> *edep_nuclear=0;
    vector<double> *edep_electron=0;
    vector<double> *eqch=0;
    vector<double> *quenchingfactor=0;
    vector<TString>* volume=0;
    
    t->SetBranchAddress("volume", &volume);
    t->SetBranchAddress("et", &et);
    t->SetBranchAddress("ex", &ex);
    t->SetBranchAddress("ey", &ey);
    t->SetBranchAddress("ez", &ez);
    t->SetBranchAddress("edep", &edep);
    t->SetBranchAddress("edep_nuclear", &edep_nuclear);
    t->SetBranchAddress("edep_electron", &edep_electron);
    t->SetBranchAddress("eqch", &eqch);
    t->SetBranchAddress("quenchingfactor", &quenchingfactor);
    
    double quenchingcut = 0.43;
    double dep_enecut= 23.0; //keV
    double dep_enecut2=269.0; //keV
    
    int eff=0;   
    int ndeps;
    bool InTPC,IsNR,IsER,IsCannotUse,IsLarge40PE,IsG2; 
    int insideTPCevents, nuclearevents, purenuclearevents, electronevents, cannotuseevents, Large40PEevents,G2events, singleevents, singlenuclearevents;
    int  purecannotuseevents, pureLarge40PEevents, pureG2events;
    insideTPCevents=0;
    nuclearevents=0;
    purenuclearevents=0;
    electronevents=0;
    cannotuseevents=0;
    Large40PEevents=0;
    G2events=0;
    purecannotuseevents=0;
    pureLarge40PEevents=0;
    pureG2events=0;
    singleevents=0;
    singlenuclearevents=0;

    for(int j=0; j<nEntries; j++)
      { 
	if(!(j%100000)) std::cout<<"Processing Event "<<j<<", " <<Int_t(100.*j/nEntries)<<"% Completed"<<std::endl;
	ndeps=0;
	t->GetEntry(j);
	InTPC = false;
	IsNR = false;
	IsER=false;
	IsCannotUse=false;
	IsLarge40PE=false;
	IsG2=false;
	for(int i=0; i< et->size(); i++)
	  {
	    if(volume->at(i) =="p_active")
	      { ndeps++;
		InTPC=true;
		if(ndeps==1)
		  { eff++;
		  }           
		if(quenchingfactor->at(i) < quenchingcut )
		  {
		    IsNR = true;
		    if( edep->at(i)<dep_enecut )
		      {
			IsCannotUse = true;
		      }
		    if( edep->at(i)>dep_enecut )
		      {
			IsLarge40PE = true;
		      }
		    if(edep->at(i)>dep_enecut && edep->at(i)<dep_enecut2 )
		      {
			IsG2 = true;
		      }
		  }
	
		if(quenchingfactor->at(i) > quenchingcut )
		  {
		    IsER = true;
		  }
	      }
	  }
	if(InTPC)               insideTPCevents++;
	if(IsNR)                nuclearevents++;
	if(IsNR && !IsER)       purenuclearevents++;
	if(IsER)	        electronevents++;
	if(IsCannotUse)         cannotuseevents++;
	if(IsLarge40PE)         Large40PEevents++;
	if(IsG2)                G2events++;
	if(IsCannotUse && !IsLarge40PE && !IsG2)         purecannotuseevents++;
	if(!IsCannotUse && IsLarge40PE && !IsG2)         pureLarge40PEevents++;
	if(!IsCannotUse && IsG2)         pureG2events++;
	if(ndeps==1)            singleevents++;	
	if(ndeps==1 && IsNR && !IsER )   singlenuclearevents++;
      }
    cout<<label<<Time<<endl;
    cout<<"InTPC= " <<insideTPCevents <<"     "<<"IsNR= "<<nuclearevents<<"  "<<"IsPureNR=  "<<purenuclearevents<<"  "<<"IsER= "<<electronevents<<endl;
    cout<<"IsCannotUse= " <<cannotuseevents <<"     "<<"IsLarge40PE= "<<Large40PEevents<<"  "<<"IsG2=  "<<G2events<<endl;
    cout<<"IspureCannotUse= " <<purecannotuseevents <<"     "<<"IspureLarge40PE= "<<pureLarge40PEevents<<"  "<<"IspureG2=  "<<pureG2events<<endl;
    cout<<"IsSingle= " <<singleevents <<"     "<<"IsSinglenuclear= "<<singlenuclearevents<<endl;
    cout<<"eff= "<<eff<<endl;
    
}


