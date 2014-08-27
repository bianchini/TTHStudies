

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
#define JETPT 40.
#define ADDTTH 1

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

void analyzeHepMC( TString sample = "test", Long64_t  tot_trials = -1, float BR = 1.,
		   string obs="Wsb",
		   TCut cut = "type==6", float fact = 0.01 , 
		   int nBins=12, float xMin = 0., float xMax = 1.){ 

  TString production = "nosmearAll_";

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
  if(SAVE) fout = TFile::Open( Form("./files/Sherpa/plots_HepMC_%s.root", obs.c_str()), "UPDATE");

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

  TFile* f = TFile::Open( "./files/Sherpa/MEM/hepMC_Sherpa_"+production+"gen_std_"+sample+".root" );
  if(f==0) return;
  if(f->IsZombie()) return;
  
  TH1D* h  = new TH1D( "h_"+sample,  sample+" ; w; d#sigma/dw", nBins, xMin, xMax);
  h->Sumw2();
  h->SetLineWidth(2);

  TTree* t = (TTree*)f->Get("tree");
  if(t==0) return;

  int doObs     = !(obs.find("Psb")!=string::npos || obs.find("Wsb")!=string::npos);
  int doObsComp = doObs && obs.find("best")!=string::npos;
  int doCMS     = (obs=="Psb") ? 1 : 0;

  TTreeFormula* treeformula    = new TTreeFormula("cat_selection",  cut , t );
  TTreeFormula* treeobservable = 0;
  if(doObs==1 && doObsComp==0){
    string formulaStr = obs;
    if( obs.find("[")!=string::npos )
      formulaStr.erase( obs.find("[") , 3 );
    treeobservable = new TTreeFormula("cat_observable", formulaStr.c_str() , t );
  }

 
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

  // jet kinematics
  float jet_pt [8];
  float jet_eta[8];
  float jet_phi[8];
  float jet_m  [8];
  int   perm_to_jet_s[99];
  int   perm_to_jet_b[99];


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
  t->SetBranchAddress("jet_pt",       jet_pt);
  t->SetBranchAddress("jet_eta",      jet_eta);
  t->SetBranchAddress("jet_phi",      jet_phi);
  t->SetBranchAddress("jet_m",        jet_m);
  t->SetBranchAddress("perm_to_jet_s",perm_to_jet_s);
  t->SetBranchAddress("perm_to_jet_b",perm_to_jet_b);
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

    // check jet pts
    //if( !(jet_pt[2]>=JETPT && jet_pt[5]>=JETPT && jet_pt[6]>=JETPT && jet_pt[7]>=JETPT && jet_pt[1]>=30) )
    //continue;


    if( doObs==0 ){

      // this number are the total probability under the ttH, ttbb, and ttjj hypotheses
      double p_125_all_s  = 0.;
      double p_125_all_b1 = 0.;
      double p_125_all_b2 = 0.;
    
      // consider only the mass scan at the nominal Higgs mass
      for(int mH_it = 0; mH_it<nMassPoints; mH_it++){	
	
	// skip other mass values...
	float mH =  mH_scan[mH_it];
	if(mH>(massH+1) || mH<(massH-1)) continue;
	
	// once we got the good Higgs mass, loop over permutations...
	for( int perm_it = 0 ; perm_it<nPermut_s ; perm_it++){		
	  float ME_prob     = p_vsMH_s[mH_it*nPermut_s + perm_it];
	  float bb_prob     = p_tt_bb[perm_it] ;
	  p_125_all_s      += ME_prob*(doCMS ? bb_prob : 1.0);
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
	  p_125_all_b1     += ME_prob*(doCMS ? bb_prob : 1.0);
	  p_125_all_b2     += ME_prob*(doCMS ? jj_prob : 1.0);
	}
      }
      
      // prepare the sb likelihood ratio ingredients
      //double wDen = p_125_all_b1*fact;
      double wDen = doCMS ? p_125_all_s + fact*(0.014056*p_125_all_b1 + 6.107426*p_125_all_b2) : (p_125_all_s + fact*p_125_all_b1);
      double wNum = p_125_all_s;
      
      // the sb likelihood ratio
      double w = wDen>0 ? wNum/wDen : 0.;
      
      h->Fill( w, SCALEsyst[0] );
      continue;
    }
    else if(doObs==1 && doObsComp==0){
      double eval  = -99.;
      if(obs.find("[")!=string::npos ){
	Int_t arraySize = treeobservable->GetNdata();
	if( arraySize<2 ){
	  cout << "Incosistent branch" << endl;
	  continue;
	}
	string pos_str = obs.substr(obs.find("[")+1, 1);
	int pos_int    = atoi(pos_str.c_str());
	eval = treeobservable->EvalInstance(pos_int);
      }
      else
	eval = treeobservable->EvalInstance();
      h->Fill( eval , SCALEsyst[0] );
      continue;
    }
    else if(doObs==1 && doObsComp==1){

      int ttH    = obs.find("best_0")!=string::npos ;

      int nPermut      =  ttH ? nPermut_s     : nPermut_b;
      float* p_vsM     =  ttH ? p_vsMH_s      : p_vsMT_b;
      
      double pmax          =   0.;
      unsigned int itermax = 999;
      
      // start the scan over the Higgs mass
      for(int mH_it = 0; mH_it<nMassPoints; mH_it++){
	
	// given mH, scan over the permutations
	for( int perm_it = 0 ; perm_it<nPermut ; perm_it++){	
	  
	  // ME probability (undet ttH hypothesis)
	  float ME_prob = p_vsM[mH_it*nPermut + perm_it];
	  
	  // b-tagging prob. (under ttbb hypothesis)
	  float bb_prob =  p_tt_bb[perm_it] ;
	  
	  // total probability
	  double p =  ME_prob*bb_prob;
	  
	  if( p>pmax ){
	    pmax    = p;
	    itermax = perm_it;
	  }
	}
      }
      if(itermax==999) continue;
      
      int permutList = ttH ? perm_to_jet_s[itermax] : perm_to_jet_b[itermax];

      int bLep_pos = permutList%1000000/100000;
      int w1_pos   = permutList%100000/10000;
      int w2_pos   = permutList%10000/1000;
      int bHad_pos = permutList%1000/100;
      int b1_pos   = permutList%100/10;
      int b2_pos   = permutList%10/1;  

      TLorentzVector TOPLEPW1(1,0,0,1);
      TLorentzVector TOPLEPW2(1,0,0,1);
      TLorentzVector MET(1,0,0,1);
      TLorentzVector TOPHADW1(1,0,0,1);
      TLorentzVector TOPHADW2(1,0,0,1);
      TLorentzVector TOPHADB (1,0,0,1);
      TLorentzVector TOPLEPB (1,0,0,1);
      TLorentzVector HIGGSB1 (1,0,0,1);
      TLorentzVector HIGGSB2 (1,0,0,1);
      
      TOPLEPW1.SetPtEtaPhiM( jet_pt[0],   jet_eta[0],   jet_phi[0],        jet_m[0] );
      // this works only for DL only
      TOPLEPW2.SetPtEtaPhiM( jet_pt[3],   jet_eta[3],   jet_phi[3],        jet_m[3]      );
      MET.SetPtEtaPhiM     ( jet_pt[1],   0.,           jet_phi[1],                   0. );
      
      TOPHADW1.SetPtEtaPhiM( jet_pt[w1_pos],   jet_eta[w1_pos],   jet_phi[w1_pos],   jet_m[w1_pos] );
      TOPHADW2.SetPtEtaPhiM( jet_pt[w2_pos],   jet_eta[w2_pos],   jet_phi[w2_pos],   jet_m[w2_pos] );
      TOPHADB. SetPtEtaPhiM( jet_pt[bHad_pos], jet_eta[bHad_pos], jet_phi[bHad_pos], jet_m[bHad_pos] );
      TOPLEPB. SetPtEtaPhiM( jet_pt[bLep_pos], jet_eta[bLep_pos], jet_phi[bLep_pos], jet_m[bLep_pos] );
      HIGGSB1. SetPtEtaPhiM( jet_pt[b1_pos],   jet_eta[b1_pos],   jet_phi[b1_pos],   jet_m[b1_pos] );
      HIGGSB2. SetPtEtaPhiM( jet_pt[b2_pos],   jet_eta[b2_pos],   jet_phi[b2_pos],   jet_m[b2_pos] );
      
      double eval  = 0.;
      if( obs=="best_0_mass_WHad"    || obs=="best_1_mass_WHad" )   eval = ( TOPHADW1+TOPHADW2 ).M();
      if( obs=="best_0_mass_TopHad"  || obs=="best_1_mass_TopHad" ) eval = ( TOPHADW1+TOPHADW2+TOPHADB ).M();
      if( obs=="best_0_mass_H"       || obs=="best_1_mass_H" )      eval = ( HIGGSB1+HIGGSB2 ).M();
      if( obs=="best_0_mass_TopLep1" || obs=="best_1_mass_TopLep1" )eval = ( MET+TOPLEPB+TOPLEPW1).M();
      if( obs=="best_0_mass_TopLep2" || obs=="best_1_mass_TopLep2" )eval = ( MET+TOPHADB+TOPLEPW2).M();
      if( obs=="best_0_e_MET"        || obs=="best_1_e_MET" )       eval = ( MET ).E();
      if( obs=="best_0_dphi_MET_WLep1"||obs=="best_1_dphi_MET_WLep1" )  eval = TMath::ACos( TMath::Cos(( MET ).DeltaPhi( TOPLEPW1 )) );
      if( obs=="best_0_dphi_MET_WLep2"||obs=="best_1_dphi_MET_WLep2" )  eval = TMath::ACos( TMath::Cos(( MET ).DeltaPhi( TOPLEPW2 )) );
      if( obs=="best_0_dphi_b1_b2" || obs=="best_1_dphi_b1_b2" )       eval = TMath::ACos( TMath::Cos(( HIGGSB2 ).DeltaPhi( HIGGSB1 )) );
      
      if( obs=="best_0_pt_WLep1" || obs=="best_1_pt_WLep1" ) eval = ( TOPLEPW1).Pt();
      if( obs=="best_0_pt_WLep2" || obs=="best_1_pt_WLep2")  eval = ( TOPLEPW2).Pt();
      if( obs=="best_0_pt_bLep"  || obs=="best_1_pt_bLep"  ) eval = ( TOPLEPB ).Pt();
      if( obs=="best_0_pt_WHad1" || obs=="best_1_pt_WHad1" ) eval = ( TOPHADW1).Pt();
      if( obs=="best_0_pt_WHad2" || obs=="best_1_pt_WHad2" ) eval = ( TOPHADW2).Pt();
      if( obs=="best_0_pt_bHad"  || obs=="best_1_pt_bHad" )  eval = ( TOPHADB ).Pt();
      if( obs=="best_0_pt_b1"    || obs=="best_1_pt_b1" )  eval = ( HIGGSB1 ).Pt();
      if( obs=="best_0_pt_b2"    || obs=="best_1_pt_b2" )  eval = ( HIGGSB2 ).Pt();
      
      if( obs=="best_0_eta_WLep1" || obs=="best_1_eta_WLep1" ) eval = ( TOPLEPW1).Eta();
      if( obs=="best_0_eta_WLep2" || obs=="best_1_eta_WLep2")  eval = ( TOPLEPW2).Eta();
      if( obs=="best_0_eta_bLep"  || obs=="best_1_eta_bLep"  ) eval = ( TOPLEPB ).Eta();
      if( obs=="best_0_eta_WHad1" || obs=="best_1_eta_WHad1" ) eval = ( TOPHADW1).Eta();
      if( obs=="best_0_eta_WHad2" || obs=="best_1_eta_WHad2" ) eval = ( TOPHADW2).Eta();
      if( obs=="best_0_eta_bHad"  || obs=="best_1_eta_bHad" )  eval = ( TOPHADB ).Eta();
      if( obs=="best_0_eta_b1"    || obs=="best_1_eta_b1" )    eval = ( HIGGSB1 ).Eta();
      if( obs=="best_0_eta_b2"    || obs=="best_1_eta_b2" )    eval = ( HIGGSB2 ).Eta();

      h->Fill( eval, SCALEsyst[0]);
      continue;
    }
    else{ /*...*/ }
    
  }
 
  cout << "Total entries = "       << h->GetEntries() << endl;

  if(tot_trials>0)
    num_of_trial = tot_trials;
  cout << "Total num of trials = " << num_of_trial << endl;
  h->Scale(1./num_of_trial/h->GetBinWidth(1)*BR);
  cout << "Total x-section = "       << h->Integral("width") << endl;

  h->Draw("HISTE");

  if(SAVE && fout){
    fout->cd();    
    h->Write("", TObject::kOverwrite);
    fout->Close();
    f->Close();
    delete c1;
  }

  return;
 
}



void analyzeHepMC_all( string obs="Wsb", float fact = 0.015, 
		       TCut cut = "type==6",
		       int nBins = 10, float xMin = 0., float xMax = 1.){

  vector<TString>  samples;
  vector<Long64_t> trials;
  vector<float>    BRs;

  samples.push_back("default");      trials.push_back(353687884);  BRs.push_back(0.0121571);
  samples.push_back("LO_full");      trials.push_back(9601859);    BRs.push_back(0.0124564);
  samples.push_back("defaultX05");   trials.push_back(232133379);  BRs.push_back(0.012151);
  samples.push_back("defaultX2");    trials.push_back(518698802);  BRs.push_back(0.0121725);
  samples.push_back("glo_soft");     trials.push_back(308303420);  BRs.push_back(0.0121661);
  samples.push_back("R_Mbb");        trials.push_back(327282417);  BRs.push_back(0.0121811);
  samples.push_back("Q_CMMPS");      trials.push_back(356945498);  BRs.push_back(0.0121679);
  samples.push_back("NNPDF");        trials.push_back(325380693);  BRs.push_back(0.0121311);
  samples.push_back("MSTW");         trials.push_back(317993869);  BRs.push_back(0.0121764);
  samples.push_back("CSS_KIN");      trials.push_back(353681469);  BRs.push_back(0.0121576);
  samples.push_back("MPI_up");       trials.push_back(353680503);  BRs.push_back(0.012157);
  samples.push_back("MPI_down");     trials.push_back(353682112);  BRs.push_back(0.012157);
  samples.push_back("ttH_LO");;      trials.push_back(3103857);    BRs.push_back(0.00712719); 

  for(unsigned int s = 0 ; s < samples.size() ; s++ ){
    analyzeHepMC( samples[s], trials[s], BRs[s], obs , cut , fact, nBins, xMin, xMax);
  }

}


void make_comparison(string obs="Wsb", string units = ""){

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

  TFile* f = TFile::Open(Form("./files/Sherpa/plots_HepMC_%s.root", obs.c_str()));
  if(f==0) return;

  vector<TString> samples;
  samples.push_back("default");
  samples.push_back("LO_full");
  samples.push_back("defaultX05");
  samples.push_back("defaultX2");
  samples.push_back("glo_soft");
  samples.push_back("R_Mbb");
  samples.push_back("Q_CMMPS");
  samples.push_back("NNPDF");
  samples.push_back("MSTW");
  samples.push_back("CSS_KIN");
  samples.push_back("MPI_up");
  samples.push_back("MPI_down");
  samples.push_back("ttH_LO");

  TH1F* h_ref = (TH1F*)((TH1F*)f->Get("h_"+samples[0]))->Clone("h_ref");


  TCanvas *c1 = new TCanvas("c1","",5,30,650,1000);
  c1->SetGrid(0,1);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  float histsize   = 0.65;
  float histmargin = 0.035;
  int variations = int(samples.size())-1-(ADDTTH ? 1 : 0);
  float step     = (histsize-histmargin)/variations;

  map<TString, TPad*> amap;
  map<TString, TH1F*> amap2;

  for( unsigned int s = 0 ; s < samples.size() ; s++ ){

    if(ADDTTH && s==(samples.size()-1)) continue;

    TPad *p0 = 0;
    if(s==0) 
      p0 = new TPad("p_"+samples[s],"",0, histsize ,1,1);
    else
      p0 = new TPad("p_"+samples[s],"",0, step*(s-1)   ,1, step*s);

    p0->SetGrid(0,0);
    p0->SetFillStyle(4000);
    p0->SetFillColor(10);
    p0->SetTicky();
    p0->SetObjectStat(0);
    p0->Draw(); 
    if(s==0) 
      amap.insert( make_pair(samples[s], p0) );
    else
      amap.insert( make_pair(samples[(samples.size()-(ADDTTH ? 1 : 0))-s], p0) );
  }

 
  
  for( unsigned int s = 0 ; s < samples.size() ; s++ ){

    //cout << "Drawing histo " << string(samples[s].Data()) << endl;

    TH1F* h = (TH1F*)f->Get("h_"+samples[s]);
    //amap2.insert( make_pair(samples[s], h) );

    if(h==0){
      cout << "Missing histo for " << string(samples[s].Data()) << endl;
      continue;
    }

    //cout << "Cd into directory... " ;
    if(ADDTTH && s==(samples.size()-1)) 
      (amap.find(samples[0])->second)->cd();
    else
      (amap.find(samples[s])->second)->cd();
    //cout << "Done" << endl;

    if(s==0){
      h->Scale(1000.);
      h->GetYaxis()->SetTitle("#frac{d#sigma}{dw} "+TString(("[fb"+( units!="" ?  "/"+units+"]" : "]")).c_str()));
      h->GetXaxis()->SetTitle("w = "+TString((obs+( units!="" ? " ["+units+"]" : "")).c_str()));
      h->GetYaxis()->SetTitleSize( 0.07  );
      h->GetYaxis()->SetTitleOffset( 0.6  );
      h->GetYaxis()->CenterTitle();
      h->GetXaxis()->SetTitleSize( 0.07  );
      h->GetXaxis()->SetTitleOffset( 0.6  );
      h->SetMinimum(0.);
      h->SetLineColor( kBlack );
      h->SetFillColor(kBlack);
      h->SetFillStyle(3004);
      h->SetTitle("t#bar{t}+b#bar{b} #rightarrow #mu^{+}#mu^{-} + jets    Sherpa+OpenLoops (MC@NLO) #sqrt{s}=8TeV");
      h->Draw("HISTE");
      amap2.insert( make_pair(samples[s], h) );
    }
    else if(ADDTTH && s==(samples.size()-1)){
      h->Scale(1000.*1.30*10);
      h->SetLineColor( 46 );
      //h->SetFillColor(46);
      //h->SetFillStyle(3005);
      h->SetLineStyle( kDashed );
      h->Draw("HISTESAME");
      amap2.insert( make_pair(samples[s], h) );
    }
    else{
      if( string(samples[s].Data()).find("LO_full")!=string::npos ) h->Scale( 2.011 );
      h->Divide(h, h_ref, 1., 1.);

      if( string(samples[s].Data()).find("defaultX")!=string::npos ) {
	h->GetYaxis()->SetRangeUser(0.40,1.60);
	h->GetYaxis()->SetNdivisions(602,kFALSE);
      }
      else{
	h->GetYaxis()->SetRangeUser(0.70,1.30);
	h->GetYaxis()->SetNdivisions(602,kFALSE);
      }

      h->GetYaxis()->SetTitle( samples[s] );
      h->GetYaxis()->SetLabelSize( 0.20  );
      h->GetYaxis()->SetTitle( Form("(%d)", s) );
      h->GetYaxis()->SetTitleSize( 0.20 );
      h->GetYaxis()->SetTitleOffset( 0.15 );
      h->SetTitle("");
      h->GetYaxis()->CenterTitle();
      h->GetXaxis()->SetTitle( "" );
      h->SetLineColor( s<=8 ? (s!=4 ? 1+s : kYellow+2) : 21+s );

      TH1F* h_clone = (TH1F*)h->Clone("h_clone");
      h_clone->SetFillColor( s<=8 ? (s!=4 ? 1+s : kYellow+2) : 21+s  );
      h_clone->SetFillStyle( 3001 );
      h_clone->SetLineWidth(1);
      amap2.insert( make_pair(samples[s], h_clone) );

      h->Draw("HIST");
      h_clone->Draw("E2SAME");
    }
   

    if(s>0 && !(ADDTTH && s==(samples.size()-1))){
      TF1* line = new TF1("line", "1", h_ref->GetXaxis()->GetXmin(),h_ref->GetXaxis()->GetXmax());
      line->SetLineStyle(kDashed);
      line->SetLineWidth(2);
      line->SetLineColor(kBlack);
      line->Draw("SAME");
    }

  }



  (amap.find(samples[0])->second)->cd(); 
  TLegend* leg = new TLegend(0.67,0.33,0.88,0.88,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04); 
  
  leg->AddEntry( (amap2.find(samples[0])->second), Form("%s",samples[0].Data()),  "L" );
  if(ADDTTH)
    leg->AddEntry( (amap2.find(samples[samples.size()-1])->second), "t#bar{t}H(#rightarrow b#bar{b}) #times 10", "L" );  

  for( unsigned int s = 1 ; s < samples.size() ; s++ ){

    if(ADDTTH && s==(samples.size()-1)) continue;
    
    if( string(samples[s].Data()).find("LO_full")!=string::npos )
      leg->AddEntry( (amap2.find(samples[s])->second), Form("(%d): LO #times k_{NLO}",s ), "F" );	
    else
      leg->AddEntry( (amap2.find(samples[s])->second), Form("(%d): %s",s, samples[s].Data()), "F" );
    
  }
 

  leg->Draw();

  if(SAVE){
    c1->cd();
    //c1->SaveAs(Form("./files/Sherpa/plot_HepMC_%s.png", obs.c_str()));
    c1->SaveAs(Form("./files/Sherpa/plot_HepMC_%s.pdf", obs.c_str()));
    f->Close();
  }
  

}



void doAll(){

  analyzeHepMC_all("Wsb",  0.01, "type==6 && numBTagM>=0", 6, 0, 1);
  make_comparison( "Wsb" );

  //return;

  //analyzeHepMC_all("Psb",    0.6, "(type==6)", 6, 0, 1);
  //make_comparison( "Psb" );

  analyzeHepMC_all("jet_pt[2]", 1, "type==6", 5, 30, 200);
  make_comparison( "jet_pt[2]" , "GeV");
  
  analyzeHepMC_all("jet_pt[5]", 1, "type==6", 5, 30, 200);
  make_comparison( "jet_pt[5]" , "GeV");
  
  analyzeHepMC_all("jet_pt[6]", 1, "type==6", 5, 30, 200);
  make_comparison( "jet_pt[6]" , "GeV");
  
  analyzeHepMC_all("jet_pt[7]", 1, "type==6", 5, 30, 200);
  make_comparison( "jet_pt[7]" , "GeV");
  
  analyzeHepMC_all("best_0_mass_H",  1, "type==6", 5, 50, 200);
  make_comparison( "best_0_mass_H" , "GeV");

  analyzeHepMC_all("best_1_mass_H",  1, "type==6", 6, 20, 200);
  make_comparison( "best_1_mass_H" , "GeV");

  //analyzeHepMC_all("best_0_pt_bLep",  1, "type==6", 6, 30, 200);
  //make_comparison( "best_0_pt_bLep" );

  //analyzeHepMC_all("best_1_pt_bLep",  1, "type==6", 6, 30, 200);
  //make_comparison( "best_1_pt_bLep" );

  //analyzeHepMC_all("best_0_pt_bHad",  1, "type==6", 6, 30, 200);
  //make_comparison( "best_0_pt_bHad" );

  //analyzeHepMC_all("best_1_pt_bHad",  1, "type==6", 6, 30, 200);
  //make_comparison( "best_1_pt_bHad" );

  //analyzeHepMC_all("best_0_pt_b1",  1, "type==6", 6, 30, 200);
  //make_comparison( "best_0_pt_b1" );

  //analyzeHepMC_all("best_1_pt_b1",  1, "type==6", 6, 30, 200);
  //make_comparison( "best_1_pt_b1" );

  analyzeHepMC_all("best_0_dphi_b1_b2",  1, "type==6", 4, 0, TMath::Pi());
  make_comparison( "best_0_dphi_b1_b2" );

  analyzeHepMC_all("best_1_dphi_b1_b2",  1, "type==6", 4, 0, TMath::Pi());
  make_comparison( "best_1_dphi_b1_b2" );

  analyzeHepMC_all("numJets",  1, "type==6", 4, 4, 8);
  make_comparison( "numJets" );

  analyzeHepMC_all("numBTagM", 1, "type==6", 3, 2, 5);
  make_comparison( "numBTagM" );

  analyzeHepMC_all("lepton_pt[0]", 1, "type==6", 5, 20, 200);
  make_comparison( "lepton_pt[0]" , "GeV");

  analyzeHepMC_all("lepton_pt[1]", 1, "type==6", 5, 20, 200);
  make_comparison( "lepton_pt[1]" , "GeV");

  analyzeHepMC_all("MET_pt", 1, "type==6", 5, 0, 200);
  make_comparison( "MET_pt", "GeV");

 

}
