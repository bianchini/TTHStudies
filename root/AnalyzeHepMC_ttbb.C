

#include <cstdlib>
#include <iostream> 
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "TMath.h"
#include "TMatrixT.h"
#include "TMatrixTBase.h"
#include "TFile.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TF2.h"


#include "TPaveText.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TList.h"
#include "TCollection.h"
#include "TObject.h"
#include "TLegend.h"
#include "THStack.h"
#include "TCut.h"
#include "TGraphAsymmErrors.h"
#include "TGraphPainter.h"
#include "TMultiGraph.h"
#include "TArrayF.h"
#include "TLine.h"
#include "TLorentzVector.h"


#define SAVE 1

typedef struct 
{
  float et; 
  float sumet;
  float sig;
  float phi;
} metInfo;

typedef struct 
{
  int run;
  int lumi;
  int event;
  int json;
} EventInfo;

typedef struct 
{
  int signalID;
  float eventscale;
  float alphaQCD;
  float alphaQED;
  float xSec;
  float xSecErr;
  int id1;
  int id2;
  float x1;
  float x2;
  float scalePDF;
  float pdf1;
  float pdf2;
  float weights[8];
} GenEventInfo;


typedef struct
{
  float bmass;
  float bpt;
  float beta;
  float bphi;
  float bstatus;
  float wdau1mass;
  float wdau1pt;
  float wdau1eta;
  float wdau1phi;
  float wdau1id;
  float wdau2mass;
  float wdau2pt;
  float wdau2eta;
  float wdau2phi;
  float wdau2id;
} genTopInfo;

typedef struct 
{
  float mass; 
  float pt;
  float eta;
  float phi;
  float status;
  float charge;
  float momid;
} genParticleInfo;


float dR( TLorentzVector reco, TLorentzVector gen){
  return TMath::Sqrt( (reco.Eta()-gen.Eta())*(reco.Eta()-gen.Eta()) +  
		      TMath::ACos( TMath::Cos(reco.Phi()-gen.Phi()) )*
		      TMath::ACos( TMath::Cos(reco.Phi()-gen.Phi()) ) )  ;
}


using namespace std;

void analyzeHepMC( TString sample = "test", TCut cut = "type==6", float fact = 0.01 ){ 

  gStyle->SetOptStat(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleH(0.07);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleOffset(1.3,"y");

  TFile* fout = 0;
  if(SAVE) fout = TFile::Open("plots_HepMC.root", "UPDATE");

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,1);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TLegend* leg = new TLegend(0.52,0.60,0.77,0.88,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04); 

  TFile* f = TFile::Open( "./files/Sherpa/MEM/hepMC_Sherpa_gen_std_"+sample+".root" );
  if(f==0) return;
  if(f->IsZombie()) return;
  
  TH1D* h  = new TH1D( "h_"+sample,  sample+" ; w; d#sigma/dw", 12, 0, 1);
  h->Sumw2();
  h->SetLineWidth(2);

  TTree* t = (TTree*)f->Get("tree");
  if(t==0) return;

  TTreeFormula* treeformula = new TTreeFormula("cat_selection", cut , t );

  float massH = 125. ;

  // ME probability  
  float p_vsMH_s[999];
  float p_vsMT_b[999];

  // b-tagging probability
  float p_tt_bb [999];
  float p_tt_jj [999];

  // # of permutations
  int nPermut_s;
  int nPermut_b;

  // num of mass points scanned
  int nMassPoints;

  // value of the mass scan grid
  float mH_scan[999];
  float mT_scan[999];


  float SCALEsyst[3];
  
  t->SetBranchAddress("p_vsMH_s",     p_vsMH_s);
  t->SetBranchAddress("p_vsMT_b",     p_vsMT_b);
  t->SetBranchAddress("p_tt_bb",      p_tt_bb);
  t->SetBranchAddress("p_tt_jj",      p_tt_jj);
  t->SetBranchAddress("nPermut_s",    &nPermut_s);
  t->SetBranchAddress("nPermut_b",    &nPermut_b);
  t->SetBranchAddress("nMassPoints",  &nMassPoints);
  t->SetBranchAddress("mH_scan",      mH_scan);
  t->SetBranchAddress("mT_scan",      mT_scan);
  t->SetBranchAddress("SCALEsyst",    SCALEsyst   );      

  Long64_t nentries = t->GetEntries();
  int num_of_trial  = 0;

  Long64_t quantile20 = nentries*1/5;
  Long64_t quantile40 = nentries*2/5;
  Long64_t quantile60 = nentries*3/5;
  Long64_t quantile80 = nentries*4/5;

  for (Long64_t i = 0; i < nentries ; i++){

    if( i==0 ){
      cout << "\e[1;32m   [                         ] (0%)\r\e[0m";
      cout.flush();
    }
    if( i==quantile20 ){
      cout << "\e[1;32m   [#####                    ] (20%)\r\e[0m";
      cout.flush();
    }
    if( i==quantile40 ){
      cout << "\e[1;32m   [##########               ] (40%)\r\e[0m";
      cout.flush();
    }
    if( i==quantile60 ){
      cout << "\e[1;32m   [###############          ] (60%)\r\e[0m";
      cout.flush();
    }
    if( i==quantile80 ){
      cout << "\e[1;32m   [####################     ] (80%)\r\e[0m";
      cout.flush();
    }
    if( i==(nentries-1) ){
      cout << "\e[1;32m   [#########################] (100%)\r\e[0m";
      cout.flush();
      cout << endl;
    }

    t->LoadTree(i);
    if( treeformula-> EvalInstance() == 0){
      //t->GetEntry(i);
      //num_of_trial += SCALEsyst[1];
      continue;
    }

    num_of_trial += SCALEsyst[1];

    // first of all, get the entry
    t->GetEntry(i);

    // this number are the total probability under the ttH, ttbb, and ttjj hypotheses
    double p_125_all_s = 0.;
    double p_125_all_b = 0.;

    // consider only the mass scan at the nominal Higgs mass
    for(int mH_it = 0; mH_it<nMassPoints; mH_it++){	
      
      // skip other mass values...
      float mH =  mH_scan[mH_it];
      if(mH>(massH+1) || mH<(massH-1)) continue;
      
      // once we got the good Higgs mass, loop over permutations...
      for( int perm_it = 0 ; perm_it<nPermut_s ; perm_it++){		
	float ME_prob     = p_vsMH_s[mH_it*nPermut_s + perm_it];
	float bb_prob     = p_tt_bb[perm_it] ;
	p_125_all_s      += ME_prob;
      }
    }
    
      // consider only the mass scan at the nominal top mass
    for(int mT_it = 0; mT_it<nMassPoints; mT_it++){
      
      // skip other mass values...
      float mT =  mT_scan[mT_it];
      if(mT>175 || mT<173) continue;
      
      // once we got the good Top mass, loop over permutations...
      for( int perm_it = 0 ; perm_it<nPermut_b ; perm_it++){		
	float ME_prob     = p_vsMT_b[mT_it*nPermut_b + perm_it];
	float bb_prob     = p_tt_bb[perm_it] ;
	float jj_prob     = p_tt_jj[perm_it] ;
	p_125_all_b      += ME_prob;
      }
    }
    
    // prepare the sb likelihood ratio ingredients
    double wDen = p_125_all_s + fact*p_125_all_b;
    double wNum = p_125_all_s;
    
    // the sb likelihood ratio
    double w = wDen>0 ? wNum/wDen : 0.;
    
    h->Fill( w, SCALEsyst[0] );
    
  }
 
  cout << "Total entries = "       << h->GetEntries() << endl;
  cout << "Total num of trials = " << num_of_trial << endl;
  h->Scale(1./num_of_trial/h->GetBinWidth(1));
  cout << "Total x-section = "       << h->Integral("width") << endl;

  h->Draw("HISTE");

  if(SAVE && fout){
    fout->cd();    
    h->Write("", TObject::kOverwrite);
    fout->Close();
    f->Close();
  }

  return;
 
}



void analyzeHepMC_all(){

  vector<TString> samples;
  samples.push_back("default");
  samples.push_back("defaultX05");
  //samples.push_back("defaultX2");
  samples.push_back("glo_soft");
  samples.push_back("R_Mbb");
  samples.push_back("Q_CMMPS");
  samples.push_back("NNPDF");
  samples.push_back("MSTW");
  samples.push_back("CSS_KIN");
  samples.push_back("MPI_up");
  samples.push_back("MPI_down");


  for(unsigned int s = 0 ; s < samples.size() ; s++ ){
    analyzeHepMC( samples[s], "type==6", 0.015);
  }

}


void make_comparison(){

  gStyle->SetOptStat(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleH(0.07);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleOffset(1.3,"y");

  TFile* f = TFile::Open("plots_HepMC.root");

  vector<TString> samples;
  samples.push_back("default");
  samples.push_back("defaultX05");
  //samples.push_back("defaultX2");
  samples.push_back("glo_soft");
  samples.push_back("R_Mbb");
  samples.push_back("Q_CMMPS");
  samples.push_back("NNPDF");
  samples.push_back("MSTW");
  samples.push_back("CSS_KIN");
  samples.push_back("MPI_up");
  samples.push_back("MPI_down");


  TCanvas *c1 = new TCanvas("c1","",5,30,650,1000);
  c1->SetGrid(0,1);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  float histsize   = 0.5;
  float histmargin = 0.07;
  int variations = int(samples.size())-1;
  float step     = (histsize-histmargin)/variations;

  map<TString, TPad*> amap;

  for( unsigned int s = 0 ; s < samples.size() ; s++ ){

    TPad *p0 = 0;
    if(s==0) 
      p0 = new TPad("p_"+samples[s],"",0, histsize ,1,1);
    else
      p0 = new TPad("p_"+samples[s],"",0, step*s   ,1, step*(s+1));

    p0->SetGrid(0,0);
    p0->SetFillStyle(4000);
    p0->SetFillColor(10);
    p0->SetTicky();
    p0->SetObjectStat(0);
    p0->Draw(); 
    amap.insert( make_pair(samples[s], p0) );
  }

 
  
  for( unsigned int s = 0 ; s < samples.size() ; s++ ){
    TH1F* h = (TH1F*)f->Get("h_"+samples[0]);
    (amap.find(samples[s])->second)->cd();
    TH1F* h0 = (TH1F*)f->Get("h_default");

    h->Divide(h, h0, 1., 1.);
    h->GetYaxis()->SetRangeUser(-1,1);
    h->GetYaxis()->SetTitle( samples[s] );

    h->Draw("HISTE");
  }


  /*
  TLegend* leg = new TLegend(0.52,0.60,0.77,0.88,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04); 
  */
  

}
