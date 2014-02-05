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

void analyzeHepMC( TString sample = "default-antikt04", int isWeighted = 0){

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
  if(SAVE) fout = new TFile("plots-"+sample+".root", "RECREATE");

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

  TH1D* h_njet  = new TH1D( "h_njet",  "; N_{jet} ; d#sigma/dN_{jet}",     8,0,8);
  TH1D* h_nbjet = new TH1D( "h_nbjet", "; N_{jet} ; d#sigma/dN_{jet}",     8,0,8);
  TH1D* h_ncjet = new TH1D( "h_ncjet", "; N_{jet} ; d#sigma/dN_{jet}",     8,0,8);

  TH1D* h_pttop  = new TH1D( "h_pttop",  "top pT; p_{T} (GeV) ;  d#sigma/dp_{T}",     20,0,400);
  TH1D* h_ptatop = new TH1D( "h_apttop", "tbar pT; p_{T} (GeV) ; d#sigma/dp_{T}",     20,0,400);

  TH1D* h_ptj1_bb = new TH1D( "h_ptj1_bb", "1st light jet pT in #geq2b events; p_{T} (GeV) ; d#sigma/dp_{T}",     20,0,400);

  TH1D* h_ptb1_bb = new TH1D( "h_ptb1_bb", "1st b-jet pT in #geq2b events; p_{T} (GeV) ; d#sigma/dp_{T}",     20,0,400);
  TH1D* h_ptb2_bb = new TH1D( "h_ptb2_bb", "2nd b-jet pT in #geq2b events; p_{T} (GeV) ; d#sigma/dp_{T}",     20,0,400);

  TH1D* h_mbb_bb = new TH1D( "h_mbb_bb", "M(bb) in #geq2b events; m_{bb} (GeV) ; d#sigma/dm_{bb}",     20,0,400);
  TH1D* h_dR_bb  = new TH1D( "h_dR_bb",  "#DeltaR(bb) in #geq2b events; #DeltaR (GeV) ; d#sigma/d#DeltaR",    20,0,5);

  TH1D* h_mbb_bb_cut = new TH1D( "h_mbb_bb_cut", "M(bb) in #geq2b events w/ M(bb)>100 GeV; m_{bb} (GeV) ; d#sigma/dm_{bb}",     20,0,400);
  TH1D* h_dR_bb_cut  = new TH1D( "h_dR_bb_cut", "#DeltaR(bb) in #geq2b events w/ M(bb)>100 GeV; #DeltaR (GeV) ; d#sigma/d#DeltaR",    20,0,5);


  h_njet ->Sumw2();
  h_nbjet->Sumw2();
  h_ncjet->Sumw2();
  h_pttop->Sumw2();
  h_ptatop->Sumw2();
  h_ptj1_bb->Sumw2();
  h_ptb1_bb->Sumw2();
  h_ptb2_bb->Sumw2();
  h_mbb_bb->Sumw2();
  h_dR_bb->Sumw2();
  h_mbb_bb_cut->Sumw2();
  h_dR_bb_cut->Sumw2();
  

  //TFile* f = TFile::Open("/scratch/bianchi/HBB_EDMNtuple/SherpaOpenLoops/hepMC-"+sample+".root");
  TFile* f = TFile::Open( "../bin/root/"+sample+".root" );

  if(f==0) return;
  if(f->IsZombie()) return;

  TTree* t = (TTree*)f->Get("tree");
  if(t==0) return;

  Long64_t nentries = t->GetEntries(); 

  //genParticleInfo genB, genBbar;
  genTopInfo      genTop, genTbar;
  //metInfo         METtype1p2corr;
  EventInfo       EVENT;
  GenEventInfo    GENEVENT;
  
  int naJets;
  Float_t aJet_pt           [NJETSMAX];
  Float_t aJet_eta          [NJETSMAX];
  Float_t aJet_phi          [NJETSMAX];
  Float_t aJet_e            [NJETSMAX];
  Float_t aJet_flavour      [NJETSMAX];

  Float_t aJet_nBs          [NJETSMAX];
  Float_t aJet_nCs          [NJETSMAX];
  Float_t aJet_nLs          [NJETSMAX];

  t->SetBranchAddress("EVENT",     &EVENT         );
  t->SetBranchAddress("GENEVENT",  &GENEVENT      );
  t->SetBranchAddress("aJet_pt",     aJet_pt      );
  t->SetBranchAddress("aJet_eta",    aJet_eta     );
  t->SetBranchAddress("aJet_phi",    aJet_phi     );
  t->SetBranchAddress("aJet_e",      aJet_e       );
  t->SetBranchAddress("aJet_flavour",aJet_flavour );
  t->SetBranchAddress("aJet_nBs",    aJet_nBs     );
  t->SetBranchAddress("aJet_nCs",    aJet_nCs     );
  t->SetBranchAddress("aJet_nLs",    aJet_nLs     );
  t->SetBranchAddress("naJets",      &naJets      );
  t->SetBranchAddress("genTop",      &genTop      );
  t->SetBranchAddress("genTbar",     &genTbar     );

  double ev_weight       = 0.;
  int    num_of_trial    = 0;
  int    num_of_trial_per_run    = 0;
  double xSec            = 0;
  double xSecErr         = 0;

  double xSecUniform            = 0;
  double xSecErrUniform         = 0;

  double xSecErrInv2     = 0;
  int num_of_runs = 1;

  double xSecAddPrevUniform    = 0.;
  double xSecErrAddPrevUniform = 0.;
  
  double xSecAddPrev        = 0.;
  double xSecErrInv2AddPrev = 0.;

  double xSecAddLast        = 0.;
  double xSecErrInv2AddLast = 0.;

  double xSecAddLastUniform        = 0.;
  double xSecErrAddLastUniform     = 0.;

  for (Long64_t i = 0; i < nentries ; i++){

    if( i%20000==0 ) cout << "Processing " << i << " event...(" << float(i)/nentries*100 << "%)" << endl;
    t->GetEntry(i);

    ev_weight    += GENEVENT.weights[0];
    num_of_trial += int(GENEVENT.weights[3]);

    double xSecAdd        = GENEVENT.xSec / GENEVENT.xSecErr / GENEVENT.xSecErr;
    double xSecErrInv2Add = 1. / GENEVENT.xSecErr / GENEVENT.xSecErr;

    if( EVENT.event==0 || i==nentries-1){
      cout << "Run " << num_of_runs << ": " <<  num_of_trial_per_run << " trials" << endl;
      num_of_trial_per_run = 0;
    }
    num_of_trial_per_run += int(GENEVENT.weights[3]);

    if( i!=0 && EVENT.event==0 ){      

      xSec        += xSecAddPrev;
      xSecErrInv2 += xSecErrInv2AddPrev;

      xSecUniform    += xSecAddPrevUniform;
      xSecErrUniform += xSecErrAddPrevUniform;

      num_of_runs++;
    }

    else if(i==nentries-1){
      xSecAddLast        += xSecAdd;
      xSecErrInv2AddLast += xSecErrInv2Add;

      xSecAddLastUniform        += GENEVENT.xSec;
      xSecErrAddLastUniform     += GENEVENT.xSecErr*GENEVENT.xSecErr;

      xSec        += xSecAdd;
      xSecErrInv2 += xSecErrInv2Add;

      xSecUniform    += GENEVENT.xSec;
      xSecErrUniform += GENEVENT.xSecErr*GENEVENT.xSecErr;
    }


    TLorentzVector jet_b1(0,0,0,0);
    TLorentzVector jet_b2(0,0,0,0);

    int countJets  = 0;
    int countbJets = 0;
    int countcJets = 0;
    int countnobJets = 0;
    for(int j = 0 ; j < naJets ; j++){

      float pt   = aJet_pt[j];
      float eta  = aJet_eta[j];
      float phi  = aJet_phi[j];
      float e    = aJet_e[j];
      float m2  = e*e - pt*pt*TMath::CosH(eta)*TMath::CosH(eta);
      if(m2<0) m2 = 0.; 
      float m      = TMath::Sqrt( m2 ); 

      TLorentzVector jet;
      jet.SetPtEtaPhiM(pt,eta, phi, m );

      if( jet.Pt()<25 ) continue;
      if( TMath::Abs(jet.Eta()) > 2.5 ) continue;

      countJets++;

      if( aJet_nBs[j] > 0 ){
	countbJets++;
	if( countbJets==1 ){
	  h_ptb1_bb->Fill( jet.Pt(), GENEVENT.weights[0] );
	  jet_b1 = jet;
	}
	if( countbJets==2 ){
	  h_ptb2_bb->Fill( jet.Pt(), GENEVENT.weights[0] );
	  jet_b2 = jet;
	}
      }
      else{
	countnobJets++;
	if( countnobJets==1 ) 
	  h_ptj1_bb->Fill( jet.Pt(), GENEVENT.weights[0]);
      }

      if( aJet_nCs[j] > 0 ) countcJets++;

    }

    h_njet ->Fill( countJets,  GENEVENT.weights[0] );
    h_nbjet->Fill( countbJets, GENEVENT.weights[0] );
    h_ncjet->Fill( countcJets, GENEVENT.weights[0] );
    
    if( jet_b1.Pt()>0 && jet_b2.Pt()>0 ){
      h_mbb_bb->Fill( (jet_b1 + jet_b2).M(), GENEVENT.weights[0] );
      h_dR_bb->Fill( dR(jet_b1, jet_b2), GENEVENT.weights[0]);

      if(  (jet_b1 + jet_b2).M()>100 ){
	h_mbb_bb_cut->Fill( (jet_b1 + jet_b2).M(), GENEVENT.weights[0] );
	h_dR_bb_cut->Fill( dR(jet_b1, jet_b2), GENEVENT.weights[0]);
      }

    }


    TLorentzVector top;
    top.SetPtEtaPhiM( genTop.wdau1pt,genTop.wdau1eta,genTop.wdau1phi,genTop.wdau1mass);
    TLorentzVector atop;
    atop.SetPtEtaPhiM( genTbar.wdau1pt,genTbar.wdau1eta,genTbar.wdau1phi,genTbar.wdau1mass);

    h_pttop ->Fill(  top.Pt() , GENEVENT.weights[0]);
    h_ptatop->Fill( atop.Pt() , GENEVENT.weights[0]);


    xSecAddPrev           = xSecAdd;
    xSecErrInv2AddPrev    = xSecErrInv2Add;
    xSecAddPrevUniform    = GENEVENT.xSec;
    xSecErrAddPrevUniform = GENEVENT.xSecErr*GENEVENT.xSecErr;

  }

  if( num_of_runs==1 ){
    xSec        = xSecAddLast;
    xSecErrInv2 = xSecErrInv2AddLast;

    xSecUniform        = xSecAddLastUniform;
    xSecErrUniform     = xSecErrAddLastUniform;
  }
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  xSec /=  xSecErrInv2;
  xSecErr = TMath::Sqrt( 1./xSecErrInv2 );

  xSecUniform /= num_of_runs;
  xSecErrUniform = TMath::Sqrt( xSecErrUniform ) / num_of_runs;

  cout << "Num of runs = " << num_of_runs << ". Num of trials = " << num_of_trial << endl;  
  cout << "cross-section from sum-of-weights          = ( " << ev_weight/num_of_trial << " ) pb" << endl;
  cout << "cross-section from HepMC record (ML)       = ( " << xSec << " +/- " << xSecErr << " ) pb" << endl;
  cout << "cross-section from HepMC record (uniform)  = ( " << xSecUniform << " +/- " << xSecErrUniform << " ) pb" << endl;
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


  double BR = xSec/(ev_weight/num_of_trial);
  cout << "BR = " << BR << endl;

  double scaleFactor = 1.;
  if( isWeighted )
    scaleFactor = 1./num_of_trial*BR;
  else
    scaleFactor = 1./ev_weight*xSec;
    

  h_njet  ->Scale(scaleFactor/h_njet  ->GetBinWidth(1) );
  h_nbjet ->Scale(scaleFactor/h_nbjet ->GetBinWidth(1) );
  h_ncjet ->Scale(scaleFactor/h_ncjet ->GetBinWidth(1) );

  h_ptj1_bb ->Scale(scaleFactor/h_ptj1_bb ->GetBinWidth(1) );
  h_ptb1_bb ->Scale(scaleFactor/h_ptb1_bb ->GetBinWidth(1) );
  h_ptb2_bb ->Scale(scaleFactor/h_ptb2_bb ->GetBinWidth(1) );
  
  h_mbb_bb ->Scale(scaleFactor/h_mbb_bb ->GetBinWidth(1) );
  h_dR_bb  ->Scale(scaleFactor/h_dR_bb ->GetBinWidth(1) );

  h_mbb_bb_cut ->Scale(scaleFactor/h_mbb_bb_cut ->GetBinWidth(1) );
  h_dR_bb_cut  ->Scale(scaleFactor/h_dR_bb_cut ->GetBinWidth(1) );

  h_pttop ->Scale(scaleFactor/h_mbb_bb_cut ->GetBinWidth(1) );
  h_ptatop ->Scale(scaleFactor/h_mbb_bb_cut ->GetBinWidth(1) );



  h_njet  ->SetLineColor(kBlack);
  //h_njet  ->SetFillColor(kBlack);
  h_njet  ->SetLineWidth(2);
  h_nbjet ->SetLineColor(kRed);
  //h_nbjet ->SetFillColor(kRed);
  h_nbjet ->SetLineWidth(2);
  h_ncjet ->SetLineColor(kBlue);
  h_ncjet ->SetLineWidth(2);

  h_ptj1_bb ->SetLineColor(kBlue);
  h_ptj1_bb ->SetLineWidth(2);
  h_ptb1_bb ->SetLineColor(kBlue);
  h_ptb1_bb ->SetLineWidth(2);
  h_ptb2_bb ->SetLineColor(kBlue);
  h_ptb2_bb ->SetLineWidth(2);

  h_mbb_bb ->SetLineColor(kBlue);
  h_mbb_bb ->SetLineWidth(2);
  h_dR_bb ->SetLineColor(kBlue);
  h_dR_bb ->SetLineWidth(2);

  h_mbb_bb_cut ->SetLineColor(kBlue);
  h_mbb_bb_cut ->SetLineWidth(2);
  h_dR_bb_cut ->SetLineColor(kBlue);
  h_dR_bb_cut ->SetLineWidth(2);

  h_pttop->SetLineColor(kBlue);
  h_pttop->SetLineWidth(2);
  h_ptatop->SetLineColor(kBlue);
  h_ptatop->SetLineWidth(2);

  double err;
  double integ = h_nbjet ->IntegralAndError(0, h_nbjet ->GetNbinsX()+1, err, "width");

  leg->SetHeader( Form("#sigma=(%.4f #pm %.4f) pb",integ, err) );
  leg->AddEntry( h_njet,  "any flavour", "L");
  leg->AddEntry( h_nbjet, "b-jets", "L");
  //leg->AddEntry( h_ncjet, "c-jets", "L");
  
  h_nbjet ->SetMinimum(0.0);

  //h_ncjet ->Draw( "HISTE" );
  h_nbjet ->Draw( "HIST" );
  h_njet  ->Draw( "HISTSAME" );
  TH1D* h_nbjet_clone = (TH1D*)h_nbjet->Clone("h_nbjet_clone");
  h_nbjet_clone  ->SetFillColor(kRed);
  h_nbjet_clone  ->SetFillStyle(3002);
  h_nbjet_clone  ->SetLineWidth(1);
  h_nbjet_clone  ->Draw("E2SAME");
  TH1D* h_njet_clone = (TH1D*)h_njet->Clone("h_njet_clone");
  h_njet_clone  ->SetFillColor(kBlack);
  h_njet_clone  ->SetFillStyle(3002);
  h_njet_clone  ->SetLineWidth(1);
  h_njet_clone  ->Draw("E2SAME");

  leg->AddEntry( h_njet_clone, "Stat. error", "F");

  leg->Draw();

  if(SAVE){
    fout->cd();
    h_njet->Write("", TObject::kOverwrite);
    h_nbjet->Write("", TObject::kOverwrite);
    h_ncjet->Write("", TObject::kOverwrite);
    h_ptj1_bb->Write("", TObject::kOverwrite);
    h_ptb1_bb->Write("", TObject::kOverwrite);
    h_ptb2_bb->Write("", TObject::kOverwrite);
    h_mbb_bb ->Write("", TObject::kOverwrite);
    h_dR_bb->Write("", TObject::kOverwrite);
    h_mbb_bb_cut ->Write("", TObject::kOverwrite);
    h_dR_bb_cut->Write("", TObject::kOverwrite);
    h_pttop ->Write("", TObject::kOverwrite);
    h_ptatop->Write("", TObject::kOverwrite);
    fout->Close();
  }
}


