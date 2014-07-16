

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


#define NJETSMAX 99
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
  return TMath::Sqrt( (reco.Eta()-gen.Eta())*(reco.Eta()-gen.Eta()) +  TMath::ACos( TMath::Cos(reco.Phi()-gen.Phi()) )*TMath::ACos( TMath::Cos(reco.Phi()-gen.Phi()) ) )  ;
}

void analyzeHepMC( TString sample = "test", float fact = 1000 ){ 

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

  TFile* f = TFile::Open( "/hdfs/local/bianchi/Sherpa/test/ROOT/ME_hepMC_"+sample+".root" );
  if(f==0) return;
  if(f->IsZombie()) return;
  
  TH1D* h  = new TH1D( "h_"+sample,  sample+" ; w; d#sigma/dw", 20, 0, 1);
  h->Sumw2();

  TTree* t = (TTree*)f->Get("tree");
  if(t==0) return;


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

  for (Long64_t i = 0; i < nentries ; i++){

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
	p_125_all_s      += bb_prob;//ME_prob;
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
	p_125_all_b      += jj_prob;//ME_prob;
      }
    }
    
    // prepare the sb likelihood ratio ingredients
    double wDen = p_125_all_s + fact*p_125_all_b;
    double wNum = p_125_all_s;
    
    // the sb likelihood ratio
    double w = wDen>0 ? wNum/wDen : 0.;
    
    h->Fill( w, SCALEsyst[0] );
    
  }
 
  h->Draw("HISTE");

  fout->cd();    
  h->Write("", TObject::kOverwrite);
  fout->Close();

  //f->Close();
 
}
