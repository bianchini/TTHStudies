#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include "DataFormats/FWLite/interface/LuminosityBlock.h"
#include "DataFormats/FWLite/interface/Run.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"

#include <cstdlib>
#include <iostream> 
#include <fstream>
#include <map>
#include <string>
#include <algorithm>

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
#include "TLeaf.h"
#include "TTreeFormula.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TF2.h"

#include "TLorentzVector.h"

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


using namespace std;


typedef TMatrixT<double> TMatrixD;

// to produce txt files, enable it
#define CREATEDATACARDS   1

#define USESHIFTEDSAMPLES 1
#define USEALLSAMPLES     1

#define VERBOSE           1

// use csv calibration from BDT
#define USECSVCALIBRATION 1
#define NCSVSYS          16

// number of extra systematics
#define NTHSYS           14 // used to be 4

// test effect of matching MEt phi distribution to data
#define RESHAPEMETPHI     0

// if doing ttbb/ttjj measurement, fit for ttbb as signal
#define FITTTBB           0

// split tt+jj into tt+cc
#define SPLITCC            1

// use combine convention
#define USECOMBINENAMES    1 


string DUMMY;

TF1* corrector = 0;

// JECup & JECdown are treated as 100% correlated with JECUp/JECDown (sys==1 & sys==2) 
const string csv_sys_names[16] = {
   "CSVLFUp",
   "CSVLFDown",
   "CSVHFStats1Up",
   "CSVHFStats1Down",
   "CSVHFStats2Up",
   "CSVHFStats2Down",
   "CSVCErr1Up",
   "CSVCErr1Down",
   "CSVCErr2Up",
   "CSVCErr2Down",
   "CSVHFUp",
   "CSVHFDown",
   "CSVLFStats1Up",  
   "CSVLFStats1Down",
   "CSVLFStats2Up",  
   "CSVLFStats2Down"
};

const string th_sys_names[14] = {
  //"Q2ScaleUp",     
  //"Q2ScaleDown",   
  "CMS_ttH_topPtcorrUp",
  "CMS_ttH_topPtcorrDown",
  "Q2scale_ttH_ttbar1pUp",
  "Q2scale_ttH_ttbar1pDown",
  "Q2scale_ttH_ttbar2pUp",
  "Q2scale_ttH_ttbar2pDown",
  "Q2scale_ttH_ttbar3pUp",
  "Q2scale_ttH_ttbar3pDown",
  "Q2scale_ttH_ttbar_bbUp",
  "Q2scale_ttH_ttbar_bbDown",
  "Q2scale_ttH_ttbar_bUp",
  "Q2scale_ttH_ttbar_bDown",
  "Q2scale_ttH_ttbar_ccUp",
  "Q2scale_ttH_ttbar_ccDown",
  //"Q2ScaleLFUp",
  //"Q2ScaleLFDown"
};

typedef struct 
{
  int run;
  int lumi;
  int event;
  int json;
} EventInfo;


// function that interpolates quadratically around a global maximum

pair<double,double> getMaxValue( TH1F* hMassProb){
 
  // the mass estimator
  float est  = -99;

  // the probability at the maximum
  float prob = -99;

  // parameters of the interpolating parabola y = a*x^2 + b*x + c
  float a    = -99;
  float b    = -99;
  float c    = -99;

  // if the histogram is empty, return (0.,0.)
  if( hMassProb->GetEntries()==0 ){
    est  =  0.;
    prob =  0.;
    return make_pair(est, prob);
  }

  // find the maximum bin
  int maxBin = hMassProb->GetMaximumBin();
  int bD = -999;
  int bC = -999;
  int bU = -999;

  // is the maximum bin the last filled ?
  int isRightMostBin = 1;
  for(int k=1; k<=hMassProb->GetNbinsX(); k++){
    if( hMassProb->GetBinContent(maxBin+k)>0 ) isRightMostBin=0;
  }

  // is the maximum bin the first filled ?
  int isLeftMostBin = 1;
  for(int k=1; k<=hMassProb->GetNbinsX(); k++){
    if( hMassProb->GetBinContent(maxBin-k)>0 ) isLeftMostBin=0;
  }

  // if the maximum bin is in the middle...
  if( (maxBin < hMassProb->GetNbinsX() && maxBin > 1) && !isRightMostBin && !isLeftMostBin ){

    // bC is the maximum bin
    bC = maxBin;

    // find the first filled bin on the left of bC (bD)...
    for(int k=1; k<=hMassProb->GetNbinsX() && bD==-999; k++){
      if( hMassProb->GetBinContent(bC-k)>0 ){
	bD = bC - k;	
      }
    }

    // find the first filled bin on the right of bC (bU)...
    for(int k=1; k<=hMassProb->GetNbinsX() && bU==-999; k++){
      if( hMassProb->GetBinContent(bC+k)>0 ){
	bU = bC + k;	
      }
    }   
  }

  // if the maximum bin is the rightmost filled
  else if( (maxBin == hMassProb->GetNbinsX()) || isRightMostBin){

    // bU is the maximum bin
    bU = maxBin; 

    // find the first filled bin on the left of bU (bC)...
    for(int k=1; k<=hMassProb->GetNbinsX() && bC==-999; k++){
      if( hMassProb->GetBinContent(bU-k)>0 ) bC = bU - k;	
    }
    // find the first filled bin on the left of bC (bD)...
    for(int k=1; k<=hMassProb->GetNbinsX() && bD==-999; k++){
      if( hMassProb->GetBinContent(bC-k)>0 ) bD = bC-k;	
    }
  }

  // if the maximum bin is the leftmost filled
  else if( (maxBin == 1) || isLeftMostBin){

    // bD is the maximum bin
    bD = maxBin; 

    // find the first filled bin on the right of bD (bC)...
    for(int k=1; k<=hMassProb->GetNbinsX() && bC==-999; k++){
      if( hMassProb->GetBinContent(bD+k)>0 ) bC = bD+k;	
    }

    // find the first filled bin on the right of bC (bU)...
    for(int k=1; k<=hMassProb->GetNbinsX() && bU==-999; k++){
      if( hMassProb->GetBinContent(bC+k)>0 ) bU = bC+k;	
    }
  }

  // else, if problems, return the maximum bin
  else{
    est  =  hMassProb->GetBinCenter  (hMassProb->GetMaximumBin());
    prob =  hMassProb->GetBinContent (hMassProb->GetMaximumBin());
    return make_pair(est, prob);
  }

  // for sanity
  if( bD==-999 || bU==-999 || bC==-999){
    est  =  hMassProb->GetBinCenter  (hMassProb->GetMaximumBin());
    prob =  hMassProb->GetBinContent (hMassProb->GetMaximumBin());
    return make_pair(est, prob);
  }

  // content (y) and position (x) of the three selected bins
  double xD =  hMassProb->GetBinCenter (bD);
  double xC =  hMassProb->GetBinCenter (bC);
  double xU =  hMassProb->GetBinCenter (bU);
  double yD =  hMassProb->GetBinContent(bD);
  double yC =  hMassProb->GetBinContent(bC);
  double yU =  hMassProb->GetBinContent(bU);
  
  // solve the linear system

  TMatrixD A(3,3);
  const double elements[9] = 
    { 1, xD,  xD*xD,
      1, xC,  xC*xC,
      1, xU,  xU*xU 
    };
  A.SetMatrixArray(elements, "C");
  TMatrixD AInv(3,3);
  double det;
  AInv = A.Invert(&det);
  
  TMatrixD Y(3,1);
  const double yPos[3] = 
    { yD, yC, yU
    };
  Y.SetMatrixArray(yPos, "C");
  
  TMatrixD C(3,1);
  const double dummy[3] = 
    { 1., 1., 1.};
  C.SetMatrixArray(dummy,"C");
  C.Mult(AInv,Y);
  
  a = C(2,0);
  b = C(1,0);
  c = C(0,0);
  
  // the coordinate of the parabola maximum
  est  = -b/2/a ;

  // the parabola evaluated at the maximum
  prob = a*est*est + b*est + c;

  // if a negative mass is obtained, just use the maximum bin
  if(est<0){
    est  =  hMassProb->GetBinCenter  (hMassProb->GetMaximumBin());
    prob =  hMassProb->GetBinContent (hMassProb->GetMaximumBin());
    return make_pair(est, prob);  
  }
  
  // return the values in a pair
  return make_pair(est, prob);
  
}


void bbb( TH1* hin, TH1* hout_Down, TH1* hout_Up, int bin){

  // get the bin content shifted down and up by 1sigma
  float bin_down = hin->GetBinContent( bin ) - hin->GetBinError( bin );
  float bin_up   = hin->GetBinContent( bin ) + hin->GetBinError( bin );

  // just copy the nominal histogram...
  hout_Down->Reset();
  hout_Down->Add(hin,1.0);

  // ... and then shift down the bin content by -1sigma
  hout_Down->SetBinContent(bin, TMath::Max(bin_down,float(0.01)));

  // just copy the nominal histogram...
  hout_Up->Reset();
  hout_Up->Add(hin,1.0);

  // ... and then shift up the bin content by +1sigma
  hout_Up->SetBinContent  (bin, TMath::Max(bin_up,float(0.01)));

  return;

}


float poorman_significance( TH1F* h_sgn=0, TH1F* h_bkg=0){  

  float out = -99;

  if( h_sgn==0 || h_bkg==0 ) return out;

  int nbins = h_sgn->GetNbinsX();

  float Ns = h_sgn->Integral(1,nbins);
  float Nb = h_bkg->Integral(1,nbins);

  if( Ns<=0 || Nb<=0 ) return out;

  float qA = 0.;
  for( int b = 1; b<=nbins ; b++){
    qA -= 2*(h_bkg->GetBinContent(b))*TMath::Log( 1 + h_sgn->GetBinContent(b)/h_bkg->GetBinContent(b) );
  }
  qA += 2*Ns;

  if(qA<=0){
    cout << "Negative value of qA!!!" << endl;
    return out;
  }

  float sigmaA = 1./TMath::Sqrt( qA );

  out = 0.5*(TMath::Erf( -1./sigmaA ) + 1 );

  if( TMath::IsNaN(out) ) out = -99;

  return out;

}



void draw(vector<float> param, TTree* t = 0, TString var = "", TH1F* h = 0, TCut cut = "", bool isMC=1){

  TCut weight = isMC ? "(weight*PUweight*trigger)" : "1";

  h->Reset();
  if(param.size()==3){
    t->Draw(var+">>"+TString(h->GetName()),   "weight"*cut );
    return;
  }
  else{

    if( param.size()!=5 ){
      cout << "Error in draw()" << endl;
      return;
    }
    
    float cutval = h->GetBinLowEdge( 3 );
    TCut    cut1( Form("%s<%f", var.Data(), cutval) );    
    TString var2(Form("p_125_all_b_ttbb/((p_125_all_b_ttbb+%f*p_125_all_b_ttjj))",    param[3]));

    TCut    cut2(Form("p_125_all_b_ttbb/((p_125_all_b_ttbb+%f*p_125_all_b_ttjj))<%f", param[3],param[4]));    
    t->Draw(var2+">>"+TString(h->GetName()),   weight*(cut&&cut1&&cut2) );
    float binL    = h->Integral();
    float binLErr = h->GetEntries()>0 ? sqrt(h->GetEntries())*h->Integral()/h->GetEntries() : 0.;
    h->Reset();
    t->Draw(var2+">>"+TString(h->GetName()),   weight*(cut&&cut1&&(!cut2)) );
    float binH    = h->Integral();
    float binHErr = h->GetEntries()>0 ? sqrt(h->GetEntries())*h->Integral()/h->GetEntries() : 0.;
    h->Reset();

    t->Draw(var+">>"+TString(h->GetName()),   weight*cut );
    h->SetBinContent(1, binL);  h->SetBinError(1, binLErr);
    h->SetBinContent(2, binH);  h->SetBinError(2, binHErr);
  }

  return;
}



TF1* getMETphiCorrection(TTree* tDATA=0, TTree* tMC=0, TCut cut = ""){

  TString var("MET_phi");

  TF1* sin_MC   = new TF1("sin_MC",   "[0] + [1]*TMath::Sin(x + [2])", -TMath::Pi(), TMath::Pi() );  
  TF1* sin_DATA = new TF1("sin_DATA", "[0] + [1]*TMath::Sin(x + [2])", -TMath::Pi(), TMath::Pi() );

  TH1F* h = new TH1F("h", "", 20, -TMath::Pi(), TMath::Pi());

  tDATA->Draw(var+">>h", cut);
  h->Fit( sin_DATA , "Q");
  sin_DATA->SetParameter(1, sin_DATA->GetParameter(1)/sin_DATA->GetParameter(0));
  sin_DATA->SetParameter(0, 1.);

  h->Reset();

  tMC->Draw(var+">>h", cut);
  h->Fit( sin_MC , "Q");
  sin_MC->SetParameter(1, sin_MC->GetParameter(1)/sin_MC->GetParameter(0));
  sin_MC->SetParameter(0, 1.);

  TF1* scale = new TF1("scale", "sin_DATA/sin_MC" , -TMath::Pi(), TMath::Pi());
  
  if(VERBOSE) cout << "Normalization of scale: " <<  scale->Integral(-TMath::Pi(), TMath::Pi())/(TMath::Pi()-(-TMath::Pi())) << endl;

  return scale;

  delete sin_MC;
  delete sin_DATA;
  delete h;
}



void fill(  TTree* tFull = 0, int nparts=1, int part=0,  TH1* h = 0, TCut cut = "" , int analysis = 0, vector<float>* param = 0, TF1* xsec = 0, int isMC=1, int pos_weight1=0, int pos_weight2=0, TString category="", string sample="", float massH=125 ){

  // needed as usual to copy a tree
  //TFile* dummy = new TFile("/scratch/bianchi/dummy_"+TString(sample.c_str())+".root","RECREATE");

  Long64_t totalentries = tFull->GetEntries();
  Long64_t nentries     = (totalentries/nparts+1);
  Long64_t firstentry   = part*nentries;

  if(VERBOSE) cout << " > run over : [" <<  firstentry << "," << firstentry+nentries-1 << "]" << endl;;  

  // tree of events in the category
  TTree* t = tFull ;//(TTree*)tFull->CopyTree( cut , "", nentries, firstentry );

  TTreeFormula* treeformula = new TTreeFormula("cat_selection", cut , t );

  // the variable name
  string cat_to_string = string(category.Data());
  
  // what kind of analysis
  bool isStandardObservable = cat_to_string.find("best_")==string::npos &&
    !( cat_to_string.find(":")!=string::npos && cat_to_string.find("::")==string::npos);
  bool isBestObservable = cat_to_string.find("best_")!=string::npos;
  bool is2DObservable   = cat_to_string.find(":")!=string::npos && cat_to_string.find("::")==string::npos;
  
  TTreeFormula* treeobservable1 = 0;
  TTreeFormula* treeobservable2 = 0;
  if(analysis==4){

    // a standard 1D variable (not a best_ analysis)
    if( isStandardObservable )
      treeobservable1 = new TTreeFormula("cat_observable1",   category , t );

    // a 2D analysis
    else if( is2DObservable ){

      size_t size = cat_to_string.size();
      size_t div  = cat_to_string.find(":");

      string var1 = cat_to_string.substr(0,    div   );
      string var2 = cat_to_string.substr(div+1,size  );

      cout << var1 << ", " << var2 << endl;

      treeobservable1 = new TTreeFormula("cat_observable1",   var1.c_str() , t );
      treeobservable2 = new TTreeFormula("cat_observable2",   var2.c_str() , t );
    }
    else{ /* ... */ }

  }


  if(VERBOSE) cout << " > cut: " << string(cut.GetTitle()) << endl;

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

  // DL or SL
  int Vtype;

  // the event-dependent weight
  float weight;
  float PUweight;
  float trigger;
  float weightTopPt;
  float weightEle;
  float weightCSV[19];
  float SCALEsyst[3];

  // the MET_phi (for reweghting)
  float MET_phi;

  // event information
  EventInfo EVENT;

  // the observable (if doMEM==4)
  //float observableF;
  //int   observableI;
  
  // number of extra partons in ME  
  int n_b,n_c,n_l,n_g;
  int nSimBs, nMatchSimBs, nMatchSimCs;

  if( analysis!=4 ){
    t->SetBranchAddress("p_vsMH_s",     p_vsMH_s);
    t->SetBranchAddress("p_vsMT_b",     p_vsMT_b);
    t->SetBranchAddress("p_tt_bb",      p_tt_bb);
    t->SetBranchAddress("p_tt_jj",      p_tt_jj);
    t->SetBranchAddress("nPermut_s",    &nPermut_s);
    t->SetBranchAddress("nPermut_b",    &nPermut_b);
    t->SetBranchAddress("nMassPoints",  &nMassPoints);
    t->SetBranchAddress("mH_scan",      mH_scan);
    t->SetBranchAddress("mT_scan",      mT_scan);
  }
  if( analysis==4 ){
    t->SetBranchAddress("Vtype",        &Vtype);
    t->SetBranchAddress("jet_pt",       jet_pt);
    t->SetBranchAddress("jet_eta",      jet_eta);
    
    if( string(category.Data()).find("best_")!=string::npos ){
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
    }
  }
  t->SetBranchAddress("weight",       &weight);
  t->SetBranchAddress("PUweight",     &PUweight);
  t->SetBranchAddress("trigger",      &trigger);
  t->SetBranchAddress("EVENT",        &EVENT);

  if( corrector ){
    t->SetBranchAddress("MET_phi",      &MET_phi);
  }
  else
    MET_phi = 0.;

  if( t->GetBranch("weightTopPt") )    
    t->SetBranchAddress("weightTopPt",&weightTopPt);
  else 
    weightTopPt = 1.0;

  if( USECSVCALIBRATION && t->GetBranch("weightCSV"))
    t->SetBranchAddress("weightCSV",     weightCSV  );      
  else{
    for( int k = 0 ; k < 19 ; k++ ){
      weightCSV[k] = 1.0;
    }
  }

  if( t->GetBranch("SCALEsyst"))
    t->SetBranchAddress("SCALEsyst",    SCALEsyst   );      
  else{
    for( int k = 0 ; k < 3 ; k++ ){
      SCALEsyst[k] = 1.0;
    }
  }

  /*
  if( analysis==4 ){
    if( t->GetBranch( category ) ){
      TLeaf* leaf = t->GetBranch( category )->GetLeaf( category );
      string leanType = string( leaf->GetTypeName());
      if( leanType.find("Float")!=string::npos ){
	t->SetBranchAddress( category ,   &observableF   );
	observableI = -998;
      }
      else if( leanType.find("Int")!=string::npos ){
	t->SetBranchAddress( category ,   &observableI   );
	observableF = -998;
      }
      else {
	observableF = -998.;	
	observableI = -998 ;	
      }
    }
    else{
      observableF = -998.;
      observableI = -998.;
    }
  }
  */

  if( t->GetBranch("n_b") ){
    t->SetBranchAddress("n_b",    &n_b);   
    t->SetBranchAddress("n_c",    &n_c);   
    t->SetBranchAddress("n_l",    &n_l);   
    t->SetBranchAddress("n_g",    &n_g);   
    t->SetBranchAddress("nSimBs", &nSimBs);  
    t->SetBranchAddress("nMatchSimBs", &nMatchSimBs);  
    if(SPLITCC && t->GetBranch("nMatchSimCs")) 
      t->SetBranchAddress("nMatchSimCs", &nMatchSimCs); 
    else
      nMatchSimCs = -99;
  }
  else{
    n_b = -99;
    n_c = -99;
    n_l = -99;
    n_g = -99;
    nSimBs      = -99;
    nMatchSimBs = -99;
    nMatchSimCs = -99;
  }

  if( t->GetBranch("weightEle") )    
    t->SetBranchAddress("weightEle",&weightEle);
  else 
    weightEle = 1.0;

  int equivalent_pos_weight2 = pos_weight2;


  if(VERBOSE) cout << " > sample weight: " << (isMC>0 ? string(Form("weight*weightEle*PUweight*trigger*weightCSV[ %d ]", pos_weight1)) : "1.0" ) ;
  if(isMC==2){

    if( pos_weight2<0 ){
      if( pos_weight2==-1 ) cout << "*(1+2*(weightTopPt-1))" ;
      if( pos_weight2==-2 ) cout << "*(1+0*(weightTopPt-1))" ;
      equivalent_pos_weight2 = 0;
    }
    else cout << "*(1+1*(weightTopPt-1))";
    cout << string(Form("*SCALEsyst[ %d ]", equivalent_pos_weight2  ));

  }
  else{
    // make sure that for not ttjets samples we don't evaluate scale unc.
    equivalent_pos_weight2 = 0;
  }
  if(isMC>0 && corrector) cout << "*corrector->Eval( MET_phi )";
  cout << endl;

  Long64_t pass    = 0;
  Long64_t entries = t->GetEntries(); 
  if(VERBOSE) cout << " > total entries: " << entries << endl;


  // a temporary histogram that contains the prob. vs mass
  TH1F* hTmp      = new TH1F("hTmp", "", 500,0,500);

  Long64_t quantile20 = firstentry + nentries*1/5;
  Long64_t quantile40 = firstentry + nentries*2/5;
  Long64_t quantile60 = firstentry + nentries*3/5;
  Long64_t quantile80 = firstentry + nentries*4/5;

  // loop over tree entries
  for (Long64_t i = firstentry; i < firstentry+nentries ; i++){
    
    if( i>= entries ) continue;

    if( i==firstentry ){
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
    if( i==(firstentry+nentries-1) ){
      cout << "\e[1;32m   [#########################] (100%)\r\e[0m";
      cout.flush();
      cout << endl;
    }


    t->LoadTree(i);
    if( treeformula-> EvalInstance() == 0){
      continue;
    }
    pass++;

    // first of all, get the entry
    t->GetEntry(i);

    // determine the event weight       
    float fill_weight = isMC>0 ? weight*weightEle*PUweight*trigger*weightCSV[ pos_weight1 ] : 1.0;

    // this is used to flag the top pt systematics
    if(isMC==2 && pos_weight2<0){
      if     ( pos_weight2 == -1 ) weightTopPt = 1 + 2.*(weightTopPt-1); // the +1s band      
      else if( pos_weight2 == -2 ) weightTopPt = 1 + 0.*(weightTopPt-1); // the -1s band
      else cout << "This position for top pt reweighting is not valid..." << endl;    
    }

    if(isMC==2 && pos_weight2>2 ){
      switch( pos_weight2 ){
      case 3:
	equivalent_pos_weight2 = (n_b+n_c+n_l+n_g)==1 ?  2 : 0; //Q2 Up
	break;
      case 4:
	equivalent_pos_weight2 = (n_b+n_c+n_l+n_g)==1 ?  1 : 0; //Q2 Down
	break;
      case 5:
	equivalent_pos_weight2 = (n_b+n_c+n_l+n_g)==2 ?  2 : 0; //Q2 Up
	break;
      case 6:
	equivalent_pos_weight2 = (n_b+n_c+n_l+n_g)==2 ?  1 : 0; //Q2 Down
	break;
      case 7:
	equivalent_pos_weight2 = nSimBs>2 && nMatchSimBs>=2 ?  2 : 0; //Q2 Up
	break;
      case 8:
	equivalent_pos_weight2 = nSimBs>2 && nMatchSimBs>=2 ?  1 : 0; //Q2 Down
	break;
      case 9:
	equivalent_pos_weight2 = nSimBs>2 && nMatchSimBs<2 ?  2 : 0; //Q2 Up
	break;
      case 10:
	equivalent_pos_weight2 = nSimBs>2 && nMatchSimBs<2 ?  1 : 0; //Q2 Down
	break;
      case 11:
	equivalent_pos_weight2 = ((!SPLITCC && nSimBs<=2) || (SPLITCC && nSimBs<=2 && nMatchSimCs<1)) ?  2 : 0; //Q2 Up
	break;
      case 12:
	equivalent_pos_weight2 = ((!SPLITCC && nSimBs<=2) || (SPLITCC && nSimBs<=2 && nMatchSimCs<1)) ?  1 : 0; //Q2 Down
	break;
      case 13:
	equivalent_pos_weight2 = (n_b+n_c+n_l+n_g)==3 ?  2 : 0; //Q2 Up
	break;
      case 14:
	equivalent_pos_weight2 = (n_b+n_c+n_l+n_g)==3 ?  1 : 0; //Q2 Down
	break;
      case 15:
	equivalent_pos_weight2 = SPLITCC && nSimBs<=2 && nMatchSimCs>0 ?  2 : 0; //Q2 Up
	break;
      case 16:
	equivalent_pos_weight2 = SPLITCC && nSimBs<=2 && nMatchSimCs>0 ?  1 : 0; //Q2 Down
	break;

      default:
	equivalent_pos_weight2 = 0;
	break;
      }
    }

    // if tt+jets, apply extra top pt reweighting
    if(isMC==2) fill_weight *= (weightTopPt*SCALEsyst[ equivalent_pos_weight2 ]);

    // if MC and we want to correct for phi modeling, correct for phi
    if(isMC>0 && corrector) fill_weight *= corrector->Eval( MET_phi );


    if( analysis==4 ){

      if( isBestObservable ){

	int ttH    = string(category.Data()).find("best_0")!=string::npos ;

	int nPermut      =  ttH ? nPermut_s     : nPermut_b;
	float* m_scan    =  ttH ? mH_scan       : mT_scan;
	float* p_vsM     =  ttH ? p_vsMH_s      : p_vsMT_b;

	double pmax          =   0.;
	unsigned int itermax = 999;

	// start the scan over the Higgs mass
	for(int mH_it = 0; mH_it<nMassPoints; mH_it++){

	  float mH =  m_scan[mH_it];
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
	if( category=="best_0_mass_WHad"    || category=="best_1_mass_WHad" )   eval = ( TOPHADW1+TOPHADW2 ).M();
	if( category=="best_0_mass_TopHad"  || category=="best_1_mass_TopHad" ) eval = ( TOPHADW1+TOPHADW2+TOPHADB ).M();
	if( category=="best_0_mass_H"       || category=="best_1_mass_H" )      eval = ( HIGGSB1+HIGGSB2 ).M();
	if( category=="best_0_mass_TopLep1" || category=="best_1_mass_TopLep1" )eval = ( MET+TOPLEPB+TOPLEPW1).M();
	if( category=="best_0_mass_TopLep2" || category=="best_1_mass_TopLep2" )eval = ( MET+TOPHADB+TOPLEPW2).M();
	if( category=="best_0_e_MET"        || category=="best_1_e_MET" )       eval = ( MET ).E();
	if( category=="best_0_dphi_MET_WLep1"||category=="best_1_dphi_MET_WLep1" )  eval = TMath::ACos( TMath::Cos(( MET ).DeltaPhi( TOPLEPW1 )) );
	if( category=="best_0_dphi_MET_WLep2"||category=="best_1_dphi_MET_WLep2" )  eval = TMath::ACos( TMath::Cos(( MET ).DeltaPhi( TOPLEPW2 )) );
	if( category=="best_0_dphi_b1_b2" || category=="best_1_dphi_b1_b2" )       eval = TMath::ACos( TMath::Cos(( HIGGSB2 ).DeltaPhi( HIGGSB1 )) );

	if( category=="best_0_pt_WLep1" || category=="best_1_pt_WLep1" ) eval = ( TOPLEPW1).Pt();
	if( category=="best_0_pt_WLep2" || category=="best_1_pt_WLep2")  eval = ( TOPLEPW2).Pt();
	if( category=="best_0_pt_bLep"  || category=="best_1_pt_bLep"  ) eval = ( TOPLEPB ).Pt();
	if( category=="best_0_pt_WHad1" || category=="best_1_pt_WHad1" ) eval = ( TOPHADW1).Pt();
	if( category=="best_0_pt_WHad2" || category=="best_1_pt_WHad2" ) eval = ( TOPHADW2).Pt();
	if( category=="best_0_pt_bHad"  || category=="best_1_pt_bHad" )  eval = ( TOPHADB ).Pt();
	if( category=="best_0_pt_b1"    || category=="best_1_pt_b1" )  eval = ( HIGGSB1 ).Pt();
	if( category=="best_0_pt_b2"    || category=="best_1_pt_b2" )  eval = ( HIGGSB2 ).Pt();

	if( category=="best_0_eta_WLep1" || category=="best_1_eta_WLep1" ) eval = ( TOPLEPW1).Eta();
	if( category=="best_0_eta_WLep2" || category=="best_1_eta_WLep2")  eval = ( TOPLEPW2).Eta();
	if( category=="best_0_eta_bLep"  || category=="best_1_eta_bLep"  ) eval = ( TOPLEPB ).Eta();
	if( category=="best_0_eta_WHad1" || category=="best_1_eta_WHad1" ) eval = ( TOPHADW1).Eta();
	if( category=="best_0_eta_WHad2" || category=="best_1_eta_WHad2" ) eval = ( TOPHADW2).Eta();
	if( category=="best_0_eta_bHad"  || category=="best_1_eta_bHad" )  eval = ( TOPHADB ).Eta();
	if( category=="best_0_eta_b1"    || category=="best_1_eta_b1" )    eval = ( HIGGSB1 ).Eta();
	if( category=="best_0_eta_b2"    || category=="best_1_eta_b2" )    eval = ( HIGGSB2 ).Eta();



	h->Fill( eval, fill_weight);
	continue;
      }
      else if( string(category.Data()).find("bjet_")!=string::npos || string(category.Data()).find("leadjet_")!=string::npos){
        float eval = 0.;
	
	float lead_jet_pt = 0;
        float lead_jet_eta = -99;
        for( unsigned int i = 2; i != 8; i++){
          //      std::cout<<"i = " <<i<<", Vtype = "<<Vtype<<", jet_pt = " << jet_pt[i]<<std::endl;
	  if( ( (Vtype==2 || Vtype==3) || ((Vtype<2 || Vtype==4) && (i==2 || i>4)) ) && ( jet_pt[i] > lead_jet_pt ) )
	    {
	      lead_jet_pt = jet_pt[i];
	      lead_jet_eta = jet_eta[i];
	    }
          //std::cout<<"lead jet pt = " << lead_jet_pt << std::endl;
        }


        if ( category=="bjet_pt" )     eval = jet_pt[2];
        if ( category=="bjet_eta" )    eval = jet_eta[2];
        if ( category=="leadjet_pt")   eval = lead_jet_pt;
        if ( category=="leadjet_eta")  eval = lead_jet_eta;

        h->Fill( eval, fill_weight);
	//	std::cout<<"eval = " <<eval<<std::endl;
        continue;
      }

      else if( is2DObservable ){
	double eval1  = treeobservable1->EvalInstance();
	double eval2  = treeobservable2->EvalInstance();
	TH2F* h_cast  = dynamic_cast<TH2F*>(h); 
	//cout << eval1 << ", " << eval2 << endl;
	h_cast->Fill( eval1, eval2, fill_weight);
	continue;
      }

      double eval  = treeobservable1->EvalInstance();
      h->Fill( eval, fill_weight);
      continue;
    }

    // reset (needed because hTmp is filled with Fill()
    hTmp->Reset();

    // if doing a mass analysis... ( 0 = scan over mH ; -1 = scan over mT )
    if( abs(analysis)==1 ){


      // choose what scan to perform
      int nPermut      = analysis>0 ? nPermut_s     : nPermut_b;
      float* m_scan    = analysis>0 ? mH_scan       : mT_scan;
      float* p_vsM     = analysis>0 ? p_vsMH_s      : p_vsMT_b;

      // reset the per-permutation mass-integrated probability
      float perm_prob [nPermut];
      for( int perm_it = 0 ; perm_it<nPermut ; perm_it++){
	perm_prob [perm_it] = 0.;
      }
      
      // start the scan over the Higgs mass
      for(int mH_it = 0; mH_it<nMassPoints; mH_it++){

	float mH =  m_scan[mH_it];

	// given mH, scan over the permutations
	for( int perm_it = 0 ; perm_it<nPermut ; perm_it++){	

	  // ME probability (undet ttH hypothesis)
	  float ME_prob = p_vsM[mH_it*nPermut + perm_it];

	  // b-tagging prob. (under ttbb hypothesis)
	  float bb_prob =  p_tt_bb[perm_it] ;

	  // nornalize weights (only if xsec is provided)
	  float norm    = xsec!=0 ? 1./xsec->Eval( mH ) : 1.0;	

	  // total probability
	  double p =  ME_prob*bb_prob*norm;	

	  // add...
	  perm_prob [perm_it] += p;	  
	}
      }
      
      // now we search for the permutation with the largest integrated probability
      float maxP   = 0.;
      int   maxM    = 0;
      int   maxPerm = 0;

      // loop over the permutations
      for( int perm_it = 0 ; perm_it<nPermut ; perm_it++){	
	if(perm_prob [perm_it] > maxP ){
	  maxP    = perm_prob [perm_it];
	  maxPerm = perm_it;
	}
      }
      
      // now that we found the best permutation (maxPerm),
      // we get its mass scan...
      for(int mH_it = 0; mH_it<nMassPoints; mH_it++){
	
	float mH      =  m_scan[mH_it];
	float ME_prob = p_vsM[mH_it*nPermut + maxPerm];
	float bb_prob =  p_tt_bb[maxPerm] ;
	float norm    = xsec!=0 ? 1./xsec->Eval( mH ) : 1.0;	
	double p      =  ME_prob*bb_prob*norm;       

	// fill the histogram
	hTmp->Fill( mH, p );
      }
      
      // now that the histogram is filled, find the global maximum
      pair<double,double> bestMass = getMaxValue(hTmp);

      // this comes from a quadratic interpolation around the maximum
      double mass = bestMass.first;

      // fill the output histo with the interpolated mass value
      h->Fill(mass, fill_weight);

    }

    // if doing a ME analysis...
    else if( abs(analysis)>=2 ){

      // this number are the total probability under the ttH, ttbb, and ttjj hypotheses
      double p_125_all_s_ttbb = 0.;
      double p_125_all_b_ttbb = 0.;
      double p_125_all_b_ttjj = 0.;

      // consider only the mass scan at the nominal Higgs mass
      for(int mH_it = 0; mH_it<nMassPoints; mH_it++){	

	// skip other mass values...
	float mH =  mH_scan[mH_it];
	if(mH>(massH+1) || mH<(massH-1)) continue;
	
	// once we got the good Higgs mass, loop over permutations...
	for( int perm_it = 0 ; perm_it<nPermut_s ; perm_it++){	

	  // if running unbiased analysis...
	  if( analysis<0 ) if( perm_it!=(EVENT.event%nPermut_s) ) continue;

	  float ME_prob     = p_vsMH_s[mH_it*nPermut_s + perm_it];
	  float bb_prob     = p_tt_bb[perm_it] ;
	  p_125_all_s_ttbb += ME_prob*bb_prob;
	}
      }

      // consider only the mass scan at the nominal top mass
      for(int mT_it = 0; mT_it<nMassPoints; mT_it++){

	// skip other mass values...
	float mT =  mT_scan[mT_it];
	if(mT>175 || mT<173) continue;

	// once we got the good Top mass, loop over permutations...
	for( int perm_it = 0 ; perm_it<nPermut_b ; perm_it++){	

	  if( analysis<0 ) if( perm_it!=(EVENT.event%nPermut_s) ) continue;

	  float ME_prob     = p_vsMT_b[mT_it*nPermut_b + perm_it];
	  float bb_prob     = p_tt_bb[perm_it] ;
	  float jj_prob     = p_tt_jj[perm_it] ;
	  p_125_all_b_ttbb += ME_prob*bb_prob;
	  p_125_all_b_ttjj += ME_prob*jj_prob;
	}
      }

      // prepare the sb likelihood ratio ingredients
      double wDen = p_125_all_s_ttbb + (*param)[0]*( (*param)[1]*p_125_all_b_ttbb + (*param)[2]*p_125_all_b_ttjj );
      double wNum = p_125_all_s_ttbb;
      
      // the sb likelihood ratio
      double w = wDen>0 ? wNum/wDen : 0.;

      if( abs(analysis)==3 ){
	wDen = (*param)[1]*p_125_all_b_ttbb + (*param)[2]*p_125_all_b_ttjj ;
	wNum = (*param)[1]*p_125_all_b_ttbb;
	w    = wDen>0 ? wNum/wDen : 0.;
      }
      
      // if not splitting the first bin, just fill with the weight
      if( param->size()==3 || (abs(analysis)!=2) ) 
	h->Fill( w, fill_weight);
      
      else{
	
	// sanity check
	if( param->size()<5 ){
	  cout << "Problem in fill..." << endl;
	  return;
	}
	
	// if w<cutval, do a 2D analyses LR_sb vs LR_bj
	float cutval  = h->GetBinLowEdge( 3 );
	float cutval2 = 0.;
	if( param->size()==5 ) 
	  cutval2 = cutval;
	else if( param->size()==6 )
	  cutval2 = h->GetBinLowEdge( 5 );

	// prepare the bj likelihood ratio ingredients
	double wbjDen = p_125_all_b_ttbb + (*param)[3]*p_125_all_b_ttjj;
	double wbjNum = p_125_all_b_ttbb ;
	
	// the bj likelihood ratio
	double wbj = wbjDen>0 ? wbjNum/wbjDen : 0.;
	
	// is above the cut for 2D analysis, just do 1D analysis
	if( w >= cutval2 )  h->Fill( w, fill_weight);

	// else fill out two bins...
	else{

	  if( w <= cutval ){
	    // if LR_bj is below the cut value param[4], fill the first bin...
	    if( wbj <= (*param)[4] )  h->Fill(   cutval/4. , fill_weight);
	    // if LR_bj is above the cut value param[4], fill the second first...
	    if( wbj >  (*param)[4] )  h->Fill( 3*cutval/4. , fill_weight);
	  }
	  else if( w > cutval &&  w < cutval2 ){
	    // if LR_bj is below the cut value param[4], fill the first bin...
	    if( wbj <= (*param)[4] )  h->Fill(  cutval +   (cutval2-cutval)/4. , fill_weight);
	    // if LR_bj is above the cut value param[4], fill the second first...
	    if( wbj >  (*param)[4] )  h->Fill(  cutval + 3*(cutval2-cutval)/4. , fill_weight);
	  }
	  else{ /* ... */ }

	}  
     		
      }
      
    }
    else{ /* NOT IMPLEMENTED */ }

  } // entries


  if(VERBOSE) cout << " > passing: " << pass << string(Form(" ( <=> %.1f weighted events)", h->Integral())) << endl;

  //dummy->Close();
  //delete dummy;
  //gSystem->Exec("rm /scratch/bianchi/dummy_"+TString(sample.c_str())+".root");

  delete treeformula;
  if(treeobservable1) 
    delete treeobservable1;
  if(treeobservable2) 
    delete treeobservable2;
  delete hTmp;

  return;

}


void appendLogNSyst(TH1* h=0, TDirectory* dir=0, string name="",  float err = 0., string& line = DUMMY ){

  if(line.find("lnN")==string::npos) 
    line = name+"                  lnN  ";

  if( h->Integral()>0 && TMath::Abs(err)>0){
    line += Form("%.3f      ", 1+err);
  }  
  else if( h->Integral()<0 )
    line += "";
  else
    line += "-          ";
  
  dir->cd();
  
  string newname = string(h->GetName());//.erase(0,2);
  
  TH1F* h_Up = (TH1F*)h->Clone((newname+"_"+name+"Up").c_str());
  h_Up->Scale( 1+err );
  h_Up->Write( (newname+"_"+name+"Up").c_str(), TObject::kOverwrite);
  
  TH1F* h_Down = (TH1F*)h->Clone( (newname+"_"+name+"Down").c_str());
  h_Down->Scale( 1./(1+err) );
  h_Down->Write( (newname+"_"+name+"Down").c_str(), TObject::kOverwrite);
  
  return;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, const char* argv[])
{

  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
  /* @@@@@@@@@@@@@@@@@@@@@@@@ FWLITE @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */

 
  std::cout << "DatacardMakerFWlite" << std::endl;
  gROOT->SetBatch(true);
 
  gSystem->Load("libFWCoreFWLite");
  gSystem->Load("libDataFormatsFWLite");

  AutoLibraryLoader::enable();


  gSystem->Exec("mkdir /scratch/bianchi/");


  PythonProcessDesc builder(argv[1]);
  const edm::ParameterSet& in = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteInput");
  string name             =  ( in.getParameter<string>  ("name" ) );
  string version          =  ( in.getParameter<string>  ("version" ) );
  string extraname        =  ( in.getParameter<string>  ("extraname" ) );
  TString fname           =  TString( (in.getParameter<string>  ("fname")).c_str() );
  TString  inputpath      =  TString( (in.getParameter<string>  ("inputpath")).c_str() );
  TString  directory      =  TString( (in.getParameter<string>  ("directory")).c_str() );
  string cut              =  ( in.getParameter<string>  ("cut" ) );
  TString category        =  TString( ( in.getParameter<string>  ("category" ) ).c_str() );
  TString varname         =  TString( ( in.getParameter<string>  ("varname" ) ).c_str() );
  int doMEM               =  ( in.getParameter<int>     ("doMEM" ) );
  float fact1             =  ( in.getParameter<double>  ("fact1" ) );
  float fact2             =  ( in.getParameter<double>  ("fact2" ) );
  float factbb            =  ( in.getParameter<double>  ("factbb" ) );
  float lumiScale         =  ( in.getParameter<double>  ("lumiScale" ) );
  int nBins               =  ( in.getParameter<int>     ("nBins" ) );
  int nBinsY              =  ( in.getUntrackedParameter<int> ("nBinsY", 1 ) );
  int splitFirstBin       =  ( in.getParameter<int>     ("splitFirstBin" ) );
  vector<double> binvec   =  ( in.getParameter<vector<double> >  ("binvec" ) );
  vector<double> binvecY  =  ( in.getParameter<vector<double> >  ("binvecY" ) );
  vector<string> samples  =  ( in.getParameter<vector<string> >  ("samples" ) );
  int nparts              =  ( in.getParameter<int>    ("nparts" ) );
  int part                =  ( in.getParameter<int>    ("part" ) );
  int analysis            =  ( in.getUntrackedParameter<int>    ("analysis", -1 ) );
  int doSystematics       =  ( in.getUntrackedParameter<int>    ("doSystematics", 1 ) );
  int runOnData           =  ( in.getUntrackedParameter<int>    ("runOnData", 1 ) );
  double massH            =  ( in.getUntrackedParameter<double> ("massH", 125 ) );
  int do2Dplots           =  ( in.getUntrackedParameter<int>    ("do2Dplots", 0 ) );

  string nJob = "";
  if( argc==3 )
    nJob = "_"+string( argv[2] );


  TF1* xsec = 0;

  // doMEM ==> 
  //  1 = mass
  //  2 = MEM
  //  3 = ttbb/ttjj
  //  (if unbiased, MEM -> -MEM)

  cout << endl;
  cout << "\e[1;32m>>>>>> produceNew will be exectuted <<<<<<\e[0m" << endl << endl;

  // cout the version
  if(VERBOSE) {
    cout << "Use trees " << name << " (version " << version << "): category " << category << endl;
    cout << "FOM index = " << doMEM << endl;
  }

  bool useHistos = false;

  // decide if it is a SL or DL analysis (important because of data)
  int isDL       = false;
  if(analysis>=0){
    if(analysis==0) isDL = false;
    else            isDL = true;
  }
  else{
    isDL = (string(category.Data()).find("cat6")!=string::npos || string(category.Data()).find("cat7")!=string::npos);
  }

  // selection cut
  string basecut = cut;
  if( USEALLSAMPLES ){
    basecut = basecut+string(Form(" && syst==%d",0));
  }

  // name of the main background (needed because one could use TTJetsSemiLept, TTJetsFullLept, ...)
  TString ttjets = "TTJets_nC";

  // name tag for the nominal samples
  TString nominal = USEALLSAMPLES ? "_all" : "_"+string(fname.Data())+"_nominal";

  // strength modifiers 
  float scaleTTH      = 1.0;//4.7;
  float scaleTTJetsLF = 1.0;//3.9;
  float scaleTTJetsHF = 1.0;//3.9;

  if(VERBOSE){
    cout << "\e[1;31m*****************************************************\e[0m" << endl;
    cout << "\e[1;31mAttention!!! The following rescaling will be applied:\e[0m" << endl;
    cout << "\e[1;31mTTH:       \e[0m" << scaleTTH << endl;
    cout << "\e[1;31mTTJets LF: \e[0m" << scaleTTJetsLF << endl;
    cout << "\e[1;31mTTJets HF: \e[0m" << scaleTTJetsHF << endl;
    cout << "\e[1;31m*****************************************************\e[0m" << endl;
  }

  // if doing the ME analysis, we first derive normalization factors to
  // properly nornalize the ME weights. They are defined such that:
  //  - S_bb is defined so that 'p_125_all_s_ttbb' for ttH  is normalized to unity in the category
  //  - B_bb is defined so that 'p_125_all_b_ttbb' for ttbb is normalized to unity in the category
  //  - B_jj is defined so that 'p_125_all_b_ttjj' for ttjj is normalized to unity in the category

  float S_bb = 1.;  float Int_S_bb = 1.;
  float B_bb = 1.;  float Int_B_bb = 1.;
  float B_jj = 1.;  float Int_B_jj = 1.;

  // param contains:
  //  param[0] =  ttbb/ttjj cross-sections
  //  param[1] =  cut value on the likelihood ratio ttbb/(ttbb+ttjj)

  vector<float> param;
  param.clear();
  
  if(abs(doMEM)==2 || abs(doMEM)==3){

    TFile* f_B = TFile::Open(inputpath+"MEAnalysis"+name+nominal+version+"_"+ttjets+".root", "OPEN");
    if( f_B==0 || f_B->IsZombie() ){
      cout << "Could not find f_B" << endl;
      return 1;
    }
    TTree* tB = (TTree*)f_B->Get("tree");

    TFile* f_S = TFile::Open(inputpath+"MEAnalysis"+name+nominal+version+"_TTH125.root");
    if( f_S==0 || f_S->IsZombie() ){
      cout << "Could not find f_S" << endl;
      return 1;
    }
    TTree* tS = (TTree*)f_S->Get("tree");
    
    if( useHistos ){

      TH1F* hInt = new TH1F("hInt","",200,-100,100);  
      
      tS->Draw("log(p_125_all_s_ttbb)>>hInt", TCut((basecut+" && nSimBs>=2 && p_125_all_s_ttbb>0").c_str())); 
      S_bb     = TMath::Exp(hInt->GetMean());
      Int_S_bb = hInt->Integral();
      if(VERBOSE) cout << "S_bb = " << S_bb <<endl;
      hInt->Reset();
      
      tB->Draw("log(p_125_all_b_ttbb)>>hInt", TCut((basecut+" && nMatchSimBs>=2 && nSimBs>2 && p_125_all_b_ttbb>0").c_str())); 
      B_bb     = TMath::Exp(hInt->GetMean());
      Int_B_bb = hInt->Integral();
      if(VERBOSE) cout << "B_bb = " << B_bb <<endl;
      hInt->Reset();
      
      tB->Draw("log(p_125_all_b_ttjj)>>hInt", TCut((basecut+" && nMatchSimBs<2  && nSimBs==2 && p_125_all_b_ttjj>0").c_str())); 
      B_jj     = TMath::Exp(hInt->GetMean());
      Int_B_jj = hInt->Integral();
      if(VERBOSE) cout << "B_jj = " << B_jj <<endl;
      hInt->Reset();
      
      delete hInt;
    }
    else{
      
      // needed as usual to copy a tree
      TFile* dummy = new TFile("/scratch/bianchi/dummy_forNorm.root","RECREATE");

      float xsec_S    = 0.1302*0.569;
      float xsec_B    = 107.66+25.81+112.33;

      // weights
      float weight;
      float PUweight;
      float trigger;

      float weightCSV[19];
      for( int k = 0 ; k < 19 ; k++ ){
	weightCSV[k] = 1.0;
      }
       
      // tree of events in the category
      S_bb = 0;
      Int_S_bb = 0;
      TTree* tS_cut = (TTree*)tS->CopyTree( TCut((basecut+" && nSimBs>=2").c_str()) );
      float p_125_all_s_ttbb;
      tS_cut->SetBranchAddress(Form("p_%d_all_s_ttbb",int(massH)),   &p_125_all_s_ttbb );
      tS_cut->SetBranchAddress("weight",             &weight);
      tS_cut->SetBranchAddress("PUweight",           &PUweight);
      tS_cut->SetBranchAddress("trigger",            &trigger);
      if( USECSVCALIBRATION && tS_cut->GetBranch("weightCSV"))
	tS_cut->SetBranchAddress("weightCSV",     weightCSV  ); 
      for (Long64_t i = 0; i < tS_cut->GetEntries() ; i++){
	tS_cut->GetEntry(i);
	S_bb += p_125_all_s_ttbb;
	Int_S_bb += weight*PUweight*trigger*weightCSV[0];
      }
      if(tS_cut->GetEntries()>0) S_bb /= tS_cut->GetEntries();
      if(VERBOSE) cout << "S_bb = " << S_bb <<endl;

      B_bb     = 0;
      Int_B_bb = 0;
      TTree* tBbb_cut = (TTree*)tB->CopyTree( TCut((basecut+" && nMatchSimBs>=2 && nSimBs>2").c_str()) );
      float p_125_all_b_ttbb;
      tBbb_cut->SetBranchAddress(Form("p_%d_all_b_ttbb",int(massH)),   &p_125_all_b_ttbb );
      tBbb_cut->SetBranchAddress("weight",             &weight);
      tBbb_cut->SetBranchAddress("PUweight",           &PUweight);
      tBbb_cut->SetBranchAddress("trigger",            &trigger);
      if( USECSVCALIBRATION && tBbb_cut->GetBranch("weightCSV"))
	tBbb_cut->SetBranchAddress("weightCSV",     weightCSV  ); 
      for (Long64_t i = 0; i < tBbb_cut->GetEntries() ; i++){
	tBbb_cut->GetEntry(i);
	B_bb     += p_125_all_b_ttbb;
	Int_B_bb += weight*PUweight*trigger*weightCSV[0];
      }
      if(tBbb_cut->GetEntries()>0) B_bb /= tBbb_cut->GetEntries();
      if(VERBOSE) cout << "B_bb = " << B_bb <<endl;
	    
      B_jj     = 0;
      Int_B_jj = 0;
      TTree* tBjj_cut = (TTree*)tB->CopyTree( TCut((basecut+" && nSimBs==2").c_str()) );
      float p_125_all_b_ttjj;
      tBjj_cut->SetBranchAddress(Form("p_%d_all_b_ttjj",int(massH)),   &p_125_all_b_ttjj );
      tBjj_cut->SetBranchAddress("weight",             &weight);
      tBjj_cut->SetBranchAddress("PUweight",           &PUweight);
      tBjj_cut->SetBranchAddress("trigger",            &trigger); 
      if( USECSVCALIBRATION && tBjj_cut->GetBranch("weightCSV"))
	tBjj_cut->SetBranchAddress("weightCSV",        weightCSV); 
      for (Long64_t i = 0; i < tBjj_cut->GetEntries() ; i++){
	tBjj_cut->GetEntry(i);
	B_jj     += p_125_all_b_ttjj;
	Int_B_jj += weight*PUweight*trigger*weightCSV[0];
      }
      if(tBjj_cut->GetEntries()>0) B_jj /= tBjj_cut->GetEntries();
      if(VERBOSE) cout << "B_jj = " << B_jj <<endl;

      dummy->Close();
      delete dummy;
      gSystem->Exec("rm /scratch/bianchi/dummy_forNorm.root");
    }   

    f_S->Close();
    f_B->Close();

    if(VERBOSE){
      cout << "Ratio S_bb/B_bb="  << S_bb/B_bb << endl;
      cout << "Ratio S_bb/B_jj="  << S_bb/B_jj << endl;
      cout << "Ratio bb/(bb+jj)=" << Int_B_bb << "/(" << Int_B_bb << "+" << Int_B_jj << ")=" << Int_B_bb/(Int_B_bb+Int_B_jj) << endl;       
      cout << "Ratio ttH/ttbb="   << Int_S_bb << "/"  << Int_B_bb << "=" << Int_S_bb/Int_B_bb << endl;             
      cout << "Ratio ttbb/ttjj="  << Int_B_bb << "/"  << Int_B_jj << "=" << Int_B_bb/Int_B_jj << endl;             
    }

    param.push_back( fact1 );
    if(abs(doMEM)!=3){
      param.push_back( TMath::Abs(fact2)<=1 ? (1-fact2)*S_bb/B_bb*(Int_B_bb/(Int_B_bb+Int_B_jj)) : 1.0);
      param.push_back( TMath::Abs(fact2)<=1 ? (1+fact2)*S_bb/B_jj*(Int_B_jj/(Int_B_bb+Int_B_jj)) : 0.0);
    }
    else{
      param.push_back( TMath::Abs(fact2)<=1 ? (1-fact2)*S_bb/B_bb : 1.0);
      param.push_back( TMath::Abs(fact2)<=1 ? (1+fact2)*S_bb/B_jj : 0.0);
    }

    if(abs(doMEM)==2 && splitFirstBin>=1){
      param.push_back( B_bb/B_jj );
      param.push_back( factbb );
      if(splitFirstBin>1)
	param.push_back( factbb );
    }

  }
  else{
    param.clear();
  }


  if(RESHAPEMETPHI){

    if(VERBOSE) cout << "\e[1;31m!!! MC will be reweighted to match MET_phi distribution observed in data !!!\e[0m" << endl;

    TFile* f_MC = TFile::Open(inputpath+"MEAnalysis"+name+nominal+version+"_"+ttjets+".root", "OPEN");
    if( f_MC==0 || f_MC->IsZombie() ){
      cout << "Could not find f_MC" << endl;
      return 1;
    }
    TTree* tMC = (TTree*)f_MC->Get("tree");

    TFile* f_DATA = TFile::Open(inputpath+"MEAnalysis"+name+nominal+version+"_Run2012_SingleMu"+".root", "OPEN");
    if( f_DATA==0 || f_DATA->IsZombie() ){
      cout << "Could not find f_DATA" << endl;
      return 1;
    }
    TTree* tDATA = (TTree*)f_DATA->Get("tree");
    
    corrector = getMETphiCorrection( tDATA, tMC, TCut(basecut.c_str()) );

    f_MC->Close();
    f_DATA->Close();

  }


  // if true, read binning from input, or just take unformly distributed bins [0,1]
  bool normalbinning = ( param.size()==3 || abs(doMEM)!=2 );
  int n_split_bins   = param.size()>=5 ? param.size()-4 : 0;

  // clean output file (if any)
  //gSystem->Exec("rm ../root/datacards/"+directory+"/"+fname+"_"+name+version+extraname+".root");

  // output file with shapes
  TFile* fout = TFile::Open("../root/datacards/"+directory+"/"+fname+"_"+name+version+extraname+nJob+".root","UPDATE");

  // directory for this particular category
  TString dirName = fname+"_"+category ;
  TDirectory* dir =  fout->GetDirectory( dirName ); 
  if( !dir) dir = fout->mkdir( dirName ) ;

  

  // fix the binning
  TArrayF bins  (nBins+1);
  TArrayF binsY (nBinsY+1);
  TArrayF bins2 (nBins+1+n_split_bins);
  if(VERBOSE){
    if(!do2Dplots){
      if( normalbinning )
	cout << "Making histograms with " << nBins << " bins:" << endl;
      else
	cout << "Making histograms with " << nBins+n_split_bins << " bins:" << endl;
    }
    else 
      cout << "Making histograms with " << nBins << " x " << nBinsY << " bins:" << endl;
  }

  // if do not split, then...
  if( normalbinning ){

    // if ME analysis, just take nBins fix-size bins in [0,1]
    if(abs(doMEM)==2 || abs(doMEM)==3){
      for(int b = 0; b < nBins+1; b++)
	bins[b] = b*1.0/nBins;
    }
    // else, read from binvec what are the axis bins
    else if(binvec.size()==(nBins+1)){
      for(int b = 0; b < nBins+1; b++) bins[b] = (binvec)[b];
      if( do2Dplots ){
	if( binvecY.size()==(nBinsY+1) ) 
	  for(int b = 0; b < nBinsY+1; b++) binsY[b] = (binvecY)[b];
	else{
	  cout << "WATCH OUT!! Mismatch nBinsY/binvecY" << endl;
	  return 1;
	}
      }

    }
    // sanity check
    else{
      cout << "WATCH OUT!! Mismatch nBins/binvec" << endl;
      return 1;
    }
  }
  // else, split the first bin into two...
  else{
    for(int b = 0; b < nBins+1+n_split_bins; b++){
      if(b<=(1+n_split_bins)) 
	bins2[b] = b*0.5/(nBins);
      else
	bins2[b] = (b-n_split_bins)*1.0/(nBins);
      if(VERBOSE) cout <<  bins2[b] << ", ";
    }
    if(VERBOSE) cout << endl;
  }

  // master histogram (all other histograms are a clone of this one)
  TH1* h = 0;

  if(do2Dplots==0)
    h = new TH1F("h","Simulation #sqrt{s}=8 TeV, "+fname+"; S/(S+B); units", 
		 normalbinning ? nBins : nBins+n_split_bins ,  
		 normalbinning ? bins.GetArray() : bins2.GetArray());
  else
    h = new TH2F("h","Simulation #sqrt{s}=8 TeV, "+fname+"; S/(S+B); units", 
		 nBins, bins.GetArray(),nBinsY, binsY.GetArray() );

  // the observable when doing the ME analysis for cat2 (workaround)
  TString var("");

  // the observable when doing the ME analysis for ttbb/ttjj separation
  TString var2(""); 

  // the observable for control plots
  TString var3("");  
  
  if(abs(doMEM)==1){
    if(VERBOSE) cout << "Variable = Higgs mass estimator" << endl;
  }
  else if(abs(doMEM)==2){
    var  = TString(Form("p_%d_all_s_ttbb/(p_125_all_s_ttbb+%f*(%f*p_125_all_b_ttbb+%f*p_125_all_b_ttjj))", 
			int(massH), param[0], param[1], param[2] ));
    if(VERBOSE) cout << "Variable = " << string(var.Data()) << endl;
    
  }
  if(abs(doMEM)==3){
    var2 = TString(Form("%f*p_%d_all_b_ttbb/(%f*p_125_all_b_ttbb+%f*p_125_all_b_ttjj)", 
			param[1], int(massH), param[1], param[2] ));
    if(VERBOSE) cout << "Variable = " << string(var2.Data()) << endl;
  }
  else if(doMEM==4){
    var3 = varname;
    if(VERBOSE) cout << "Variable = " << string(var3.Data()) << endl;
  }
  else{}

  // list of MC processes that will be considered 
  // (N.B. the name must match the name of the input trees)
 
  int countRun2012samples = 0;

  // these are al the samples needed in the datacard
  vector<TString> all_datacard_samples;
  all_datacard_samples.push_back("TTH125");
  all_datacard_samples.push_back("TTV");
  all_datacard_samples.push_back("DiBoson");
  all_datacard_samples.push_back("SingleT");
  all_datacard_samples.push_back("EWK");
  all_datacard_samples.push_back("TTJetsHFbb");
  all_datacard_samples.push_back("TTJetsHFb");
  if(SPLITCC)
    all_datacard_samples.push_back("TTJetsLFcc");
  all_datacard_samples.push_back("TTJetsLF");
  all_datacard_samples.push_back("data_obs");

  vector<TString> all_datacard_samples_COMBINE;
  all_datacard_samples_COMBINE.push_back("ttH_hbb");
  all_datacard_samples_COMBINE.push_back("ttbarV");
  all_datacard_samples_COMBINE.push_back("diboson");
  all_datacard_samples_COMBINE.push_back("singlet");
  all_datacard_samples_COMBINE.push_back("ewk");
  all_datacard_samples_COMBINE.push_back("ttbarPlusBBbar");
  all_datacard_samples_COMBINE.push_back("ttbarPlusB");
  if(SPLITCC)
    all_datacard_samples_COMBINE.push_back("ttbarPlusCCbar");
  all_datacard_samples_COMBINE.push_back("ttbar");
  all_datacard_samples_COMBINE.push_back("data_obs");


  // list of names that will appear in datacards 
  // (N.B. same order as in samples)
  vector<TString> datacard_samples;

  // loop over the processes
  for( unsigned int s = 0 ; s < samples.size(); s++){
    
    string sample = samples[s]; 

    // this cut selects the category   
    TCut sample_cut(cut.c_str());

    // attach to its file a name tag (same as appearing in datacards)
    string datacard_name         = "";
    string datacard_name_COMBINE = "";

    if( sample.find("TTH")!=string::npos ){
      datacard_name = "TTH125";
      datacard_name_COMBINE = "ttH_hbb";
    }
    if( sample.find("TTV")!=string::npos ){
      datacard_name = "TTV";
      datacard_name_COMBINE = "ttbarV";
    }
    if( sample.find("DiBoson")!=string::npos ){
      datacard_name = "DiBoson";
      datacard_name_COMBINE = "diboson";
    }
    if( sample.find("SingleT")!=string::npos ){
      datacard_name = "SingleT";
      datacard_name_COMBINE = "singlet";
    }
    if( sample.find("EWK")!=string::npos ){
      datacard_name = "EWK";
      datacard_name_COMBINE = "ewk";
    }    
    if( sample.find("TTJetsBB")!=string::npos ){
      sample        = string(ttjets.Data());
      sample_cut    = sample_cut && TCut("nSimBs>2 && nMatchSimBs>=2");
      datacard_name = "TTJetsHFbb";
      datacard_name_COMBINE = "ttbarPlusBBbar";
    }
    if( sample.find("TTJetsBJ")!=string::npos ){
      sample        = string(ttjets.Data());
      sample_cut    = sample_cut && TCut("nSimBs>2 && nMatchSimBs<2");
      datacard_name = "TTJetsHFb";
      datacard_name_COMBINE = "ttbarPlusB";
    }
    if( SPLITCC && sample.find("TTJetsCC")!=string::npos ){
      sample        = string(ttjets.Data());
      sample_cut    = sample_cut && TCut("nSimBs==2 && nMatchSimCs>0");
      datacard_name = "TTJetsLFcc";
      datacard_name_COMBINE = "ttbarPlusCCbar";
    }
    if( sample.find("TTJetsJJ")!=string::npos ){
      sample        = string(ttjets.Data());
      if(SPLITCC)
	sample_cut    = sample_cut && TCut("nSimBs==2 && nMatchSimCs<1");
      else
	sample_cut    = sample_cut && TCut("nSimBs==2");
      datacard_name = "TTJetsLF";
      datacard_name_COMBINE = "ttbar";
    }
    if( sample.find("Run2012")!=string::npos ){

      datacard_name         = "data_obs";
      datacard_name_COMBINE = "data_obs";

      TCut triggerVtype0 = TCut("(Vtype==0 && ( triggerFlags[22]>0 || triggerFlags[23]>0  ))");                        // mm
      TCut triggerVtype1 = TCut("(Vtype==1 && ( triggerFlags[6]>0 ) )");                                               // ee
      TCut triggerVtype2 = TCut("(Vtype==2 && ( triggerFlags[22]>0 || triggerFlags[23]>0 || triggerFlags[47]>0 ))");   // m
      TCut triggerVtype3 = TCut("(Vtype==3 && ( triggerFlags[44]>0 ) )");                                              // e
      TCut triggerVtype4 = TCut("(Vtype==4 && ( triggerFlags[22]>0 || triggerFlags[23]>0  ))");                        // em

      if( sample.find("SingleMu")!=string::npos       &&  isDL)
	sample_cut = sample_cut && (triggerVtype0 || triggerVtype4);
      if( sample.find("SingleMu")!=string::npos       && !isDL)
	sample_cut = sample_cut && triggerVtype2;
      if( sample.find("DoubleElectron")!=string::npos &&  isDL)
	sample_cut = sample_cut && triggerVtype1;
      if( sample.find("SingleElectron")!=string::npos && !isDL)
	sample_cut = sample_cut && triggerVtype3;
     
      //sample_cut = sample_cut && (triggerVtype0 || triggerVtype1 || triggerVtype2 || triggerVtype3);

      countRun2012samples++;
    }

    if(sample.find("Run2012")==string::npos || 
       (sample.find("Run2012")!=string::npos && countRun2012samples==1)) 
      datacard_samples.push_back( TString(datacard_name.c_str()) );


    // list of samples for systematics
    vector<string> systematics;

    // provided only for tt+jets and ttH at the moment
    if(USEALLSAMPLES){

      systematics.push_back("all");  // nominal
      
      if(sample.find("Run2012")==string::npos && doSystematics){
	if( !USECSVCALIBRATION ){
	  systematics.push_back("all");// cvsUp
	  systematics.push_back("all");// csvDown
	}
	else{
	  for(int sy = 0 ; sy<NCSVSYS ; sy++)
	    systematics.push_back("all");// csv calibration sys
	}
	for(int sy = 0 ; sy<NTHSYS ; sy++)
	  systematics.push_back("all");// theory sys

	systematics.push_back("all");// JESUp
	systematics.push_back("all");// JESDown
	systematics.push_back("all");// JERUp
	systematics.push_back("all");// JERDown
      }

    }
    
    else if( !USEALLSAMPLES && USESHIFTEDSAMPLES ){
      systematics.push_back("nominal");
      if( (sample.find("TTH")!=string::npos || sample.find("TTJets")!=string::npos) ){
	systematics.push_back("csvUp");
	systematics.push_back("csvDown");
	systematics.push_back("JECUp");
	systematics.push_back("JECDown");
      }
    }
    
    else{
      systematics.push_back("nominal");
    }

    // the input file
    TFile* f = 0;

    // loop over the systematics
    for( unsigned int sy = 0 ; sy < systematics.size(); sy++){
      
      string sys        = systematics[sy];

      string syst_name = sys;
      if( USEALLSAMPLES && !USECSVCALIBRATION){
	if(sy==0)
	  syst_name = "nominal";
	else if(sy==1)
	  syst_name = "csvUp";
	else if(sy==2)
	  syst_name = "csvDown";
	else if(sy>=3 && sy <= 2+NTHSYS )
	  syst_name = string(Form("%s", th_sys_names[sy-3].c_str()));
	else if(sy==2+NTHSYS+1)
	  syst_name = "JECUp";
	else if(sy==2+NTHSYS+2)
	  syst_name = "JECDown";
	else if(sy==2+NTHSYS+3)
	  syst_name = "JERUp";
	else if(sy==2+NTHSYS+4)
	  syst_name = "JERDown";
      }
      else if( USEALLSAMPLES && USECSVCALIBRATION ){
	if(sy==0)
	  syst_name = "nominal";
	else if(sy>=1 && sy <= NCSVSYS)
	  syst_name = string(Form("CMS_ttH_%s",   csv_sys_names[sy-1].c_str()));
	else if(sy>=(NCSVSYS+1) && sy <= (NCSVSYS+NTHSYS))
	  syst_name = string(Form("%s",           th_sys_names[sy-(NCSVSYS+1)].c_str()));
	else if(sy==(NCSVSYS+NTHSYS)+1)
	  syst_name = "CMS_scale_jUp";
	else if(sy==(NCSVSYS+NTHSYS)+2)
	  syst_name = "CMS_scale_jDown";
	else if(sy==(NCSVSYS+NTHSYS)+3)
	  syst_name = "CMS_res_jUp";
	else if(sy==(NCSVSYS+NTHSYS)+4)
	  syst_name = "CMS_res_jDown";
	else{ cout << "... too many systematics ..." << endl;}
      }

      string event_type = !USEALLSAMPLES ? "_"+string(fname.Data()) : "";


      TCut syst_sample_cut = sample_cut;

      if( USEALLSAMPLES && !USECSVCALIBRATION){

	int equivalent_sy = sy;

	  // this is for the nominal and all the theory systematics
	if( sy==0 || (sy>=3 && sy <= 2+NTHSYS) )
	    equivalent_sy = 0;

	// this is for the csv systematics
	else if( sy>=1 && sy<=2 )
	  equivalent_sy = sy;

	// this is for the others
	else
	  equivalent_sy = (sy - NTHSYS ) ;
	
	// apply the cut
	syst_sample_cut = sample_cut && TCut(Form("syst==%d",equivalent_sy));

      }
      else if( USEALLSAMPLES && USECSVCALIBRATION ){

	  int equivalent_sy = sy;

	  // this is for the nominal and all the csv systematics
	  if( sy <= (NCSVSYS+NTHSYS))
	    equivalent_sy = 0;

	  // else, use the shifted events
	  else
	    equivalent_sy = (sy - (NCSVSYS+NTHSYS) + 2) ;
	    
	  // apply the cut
	  syst_sample_cut = sample_cut && TCut(Form("syst==%d",equivalent_sy));
      }
      

      // input file
      string inputfilename = string(inputpath.Data())+"MEAnalysis"+name+event_type+"_"+sys+""+version+"_"+sample+".root";
      if( USEALLSAMPLES && f==0 ) 
	f = TFile::Open(inputfilename.c_str());
      else if( !USEALLSAMPLES && f==0 )
	f = TFile::Open(inputfilename.c_str());
      
      if(f==0 || f->IsZombie()){
	cout << "Missing " << sample << " file" << endl;
	continue;
      }
      else{
	if(VERBOSE) cout << "\e[1;34mDoing......... sample=" << sample << "  process=" << datacard_name << "  analysis=" << syst_name << "\e[0m" << endl;
      }
      
      // the input tree
      TTree* tree = (TTree*)f->Get("tree");
      
      // make sure that all the sample has been processed... this may be dangerous!
      float missing_job = 1.;
      TH1F* hcounter = (TH1F*)f->Get("hcounter");
      if( hcounter ){
	float total_job = hcounter->GetBinContent(1);       
	  if(sample=="TTV")     missing_job = 2./total_job;
	  if(sample=="SingleT") missing_job = 6./total_job;
	  if(sample=="DiBoson") missing_job = 3./total_job;
	  if(sample.find("TTJets")!=string::npos) missing_job = 3./total_job;
	  if(sample.find("TTH125")!=string::npos) missing_job = 1./total_job;
	  if(sample.find("EWK")!=string::npos) missing_job    = 3./total_job;
	  if(VERBOSE) cout << "\e[1;31m scale by " << missing_job << " (" << total_job <<  " partial files)\e[0m" << endl;
      }

      // histogram for the particluar process/systematics
      TH1* h_tmp = 0;
      if(do2Dplots==0) 
	h_tmp = (TH1F*)h->Clone(("h_"+datacard_name).c_str());
      else 
	h_tmp = (TH2F*)h->Clone(("h_"+datacard_name).c_str());
      h_tmp->Reset();
      h_tmp->Sumw2();
     
      
      // decide if this is data or MC (or tt+jets, in which case apply top pt reweighting)
      int isMC = 0;
      if( sample.find("Run2012")==string::npos ) isMC++;
      if( sample.find("TTJets")!=string::npos )  isMC++;

      // decide which weight to use for filling
      int pos_weight1 = 0;
      if(USECSVCALIBRATION){

	// if doing JECUp, read element 1 of weightCSV
	if( syst_name.find("CMS_scale_jUp")  !=string::npos )
	  pos_weight1 = 1;

	// if doing JECDown, read element 2 of weightCSV
	else if( syst_name.find("CMS_scale_jDown")!=string::npos )
	  pos_weight1 = 2;

	// if doing any of the other csv systematics, shift by two units
	else if(sy>=1 && sy <= NCSVSYS)
	  pos_weight1 = sy+2;

	// else, use the nominal one (central value)
	else
	  pos_weight1 = 0;
      }

      int pos_weight2 = 0;
      //if( syst_name.find( "Q2ScaleUp"   )  !=string::npos )
      //pos_weight2 = 2;
      //if( syst_name.find( "Q2ScaleDown" )  !=string::npos )
      //pos_weight2 = 1;
      if( syst_name.find( "CMS_ttH_topPtcorrUp"     )  !=string::npos )
	pos_weight2 = -1;
      if( syst_name.find( "CMS_ttH_topPtcorrDown"   )  !=string::npos )
	pos_weight2 = -2;
      if( syst_name.find( "Q2scale_ttH_ttbar1pUp"   )  !=string::npos )
	pos_weight2 = 3;
      if( syst_name.find( "Q2scale_ttH_ttbar1pDown" )  !=string::npos )
	pos_weight2 = 4;
      if( syst_name.find( "Q2scale_ttH_ttbar2pUp"   )  !=string::npos )
	pos_weight2 = 5;
      if( syst_name.find( "Q2scale_ttH_ttbar2pDown" )  !=string::npos )
	pos_weight2 = 6;
      if( syst_name.find( "Q2scale_ttH_ttbar_bbUp"   )  !=string::npos )
	pos_weight2 = 7;
      if( syst_name.find( "Q2scale_ttH_ttbar_bbDown" )  !=string::npos )
	pos_weight2 = 8;
      if( syst_name.find( "Q2scale_ttH_ttbar_bUp"   )   !=string::npos )
	pos_weight2 = 9;
      if( syst_name.find( "Q2scale_ttH_ttbar_bDown" )   !=string::npos )
	pos_weight2 = 10;
      //if( syst_name.find( "Q2ScaleLFUp"   )    !=string::npos )
      //pos_weight2 = 11;
      //if( syst_name.find( "Q2ScaleLFDown" )    !=string::npos )
      //pos_weight2 = 12;
      if( syst_name.find( "Q2scale_ttH_ttbar3pUp"   )  !=string::npos )
	pos_weight2 = 13;
      if( syst_name.find( "Q2scale_ttH_ttbar3pDown" )  !=string::npos )
	pos_weight2 = 14;
      if( syst_name.find( "Q2scale_ttH_ttbar_ccUp"   )    !=string::npos )
	pos_weight2 = 15;
      if( syst_name.find( "Q2scale_ttH_ttbar_ccDown" )    !=string::npos )
	pos_weight2 = 16;


      // fill the histogram
      fill( tree, nparts, part, h_tmp, syst_sample_cut, doMEM, &param, xsec, isMC, pos_weight1 , pos_weight2, varname, sample+"_"+extraname , massH);
      
      if(do2Dplots==0){
	// add underflow bin to the first bin...
	int firstBin              =  1;
	float firstBinContent     =  h_tmp->GetBinContent( firstBin ); 
	float underflowBinContent =  h_tmp->GetBinContent( firstBin-1 ); 
	h_tmp->SetBinContent( firstBin, firstBinContent+underflowBinContent); 
	
	// add overflow bin into the last bin...
	int lastBin              =  h_tmp->GetNbinsX();
	float lastBinContent     =  h_tmp->GetBinContent( lastBin ); 
	float overflowBinContent =  h_tmp->GetBinContent( lastBin+1 ); 
	h_tmp->SetBinContent( lastBin, lastBinContent+overflowBinContent);    
      }
      
      // scale to the target luminosity
      h_tmp->Scale( sample.find("Run2012")==string::npos          ? lumiScale*missing_job  : 1.0 );
      h_tmp->Scale( datacard_name.find("TTJetsHFb")!=string::npos ? scaleTTJetsHF          : 1.0 );
      h_tmp->Scale( datacard_name.find("TTJetsLF") !=string::npos ? scaleTTJetsLF          : 1.0 );
      h_tmp->Scale( datacard_name.find("TTH")!=string::npos       ? scaleTTH               : 1.0 );

      
      // cd to the directory and save there
      dir->cd();

      // decide how to call your histograms      
      string save_name = USECOMBINENAMES ? datacard_name_COMBINE : datacard_name;

      // if doing nominal analysis, save histogram...
      if( sy==0 ){       

	// if doing a MC sample, just overwrite the directory
	if(sample.find("Run2012")==string::npos ){
	  h_tmp->SetName(save_name.c_str());
	  h_tmp->Write(save_name.c_str(), TObject::kOverwrite);
	}
	
	// else, need to either overwrite (if first data sample), or add (if second or more)
	else{
	  if(countRun2012samples==0)  
	    cout << "Inconsistency (1) found when filling Run2012" << endl;
	  if(countRun2012samples==1){
	    h_tmp->SetName(save_name.c_str());
	    h_tmp->Write(save_name.c_str(), TObject::kOverwrite);
	  }
	  if(countRun2012samples>=2){
	    TH1* h_tmp_tmp = 0;
	    if(do2Dplots==0) 
	      h_tmp_tmp = (TH1F*)dir->FindObjectAny(save_name.c_str());
	    else 
	      h_tmp_tmp = (TH2F*)dir->FindObjectAny(save_name.c_str());

	    if( h_tmp_tmp==0 ) 
	      cout << "Inconsistency (2) found when filling Run2012" << endl;
	    else{
	      h_tmp_tmp->Add( h_tmp , 1.0);
	      dir->cd();
	      h_tmp_tmp->SetName(save_name.c_str());
	      h_tmp_tmp->Write(save_name.c_str(), TObject::kOverwrite);
	    }
	  }

	  // close the input file
	  if( !USEALLSAMPLES )
	    f->Close();
	  
	  // don't need bbb for data
	  continue; 

	}  
      }

      // if doing systematic analysis, save histogram with correct label and continue
      // no need to code for Run2012 occurrence...
      else{
	h_tmp->SetName((save_name+"_"+syst_name).c_str());
	h_tmp->Write((save_name+"_"+syst_name).c_str(), TObject::kOverwrite);
	
	// close the input file
	if( !USEALLSAMPLES )
	  f->Close();
	
	// don't need bbb for shifted samples
	continue;      
      }
      
      
      // add bin-by-bin systematics
      for(int bin = 1; bin<=h_tmp->GetNbinsX(); bin++ ){

	if(do2Dplots) continue;

	TH1F* h_tmp_b_up   = (TH1F*)h->Clone(Form("h_tmp_%d_Up",  bin));
	TH1F* h_tmp_b_down = (TH1F*)h->Clone(Form("h_tmp_%d_Down",bin));

	// the actual function that produces the shifted templates
	bbb( h_tmp, h_tmp_b_down, h_tmp_b_up, bin);      
	dir->cd();

	// save only if the bin is filled
	if(h_tmp->GetBinContent(bin)>0){
	  h_tmp_b_up  ->SetName(Form("%s_%s%sbin%dUp",   save_name.c_str(), save_name.c_str(), category.Data(), bin));
	  h_tmp_b_up  ->Write(Form("%s_%s%sbin%dUp",     save_name.c_str(), save_name.c_str(), category.Data(), bin), TObject::kOverwrite);
	}
	if(h_tmp->GetBinContent(bin)>0){
	  h_tmp_b_down->SetName(Form("%s_%s%sbin%dDown", save_name.c_str(), save_name.c_str(), category.Data(), bin));
	  h_tmp_b_down->Write(Form("%s_%s%sbin%dDown",   save_name.c_str(), save_name.c_str(), category.Data(), bin), TObject::kOverwrite);      	  
	}
      }
      
      // close the input file
      if( !USEALLSAMPLES )
	f->Close();
      
    }   
   
    if( USEALLSAMPLES )
      f->Close();
 
  }


  // return if you don't want to create the datadacards
  if( do2Dplots==1 ){
    fout->Close();
    return 1;
  }


  // for later use, use a map between process name and histograms
  std::map<TString, TH1*> aMap;

  // for the moment, use th sum of backgrounds as 'data' (asymov dataset)
  TH1F* h_data = (TH1F*)h->Clone("h_data");
  h_data->Reset();

  // observation is the sum of all processes yields
  float observation = 0.;

  // keep track of missing files
  int errFlag = 0 ;

  // sum-up the MCs to store the total SoB
  float totBkg   = 0.;
  TH1F* h_totBkg = 0;
  float totSgn   = 0.;
  TH1F* h_totSgn = 0;

  cout << "Yields:" << endl;
  for( unsigned int s = 0 ; s < all_datacard_samples.size(); s++){

    // decide how to call your histograms      
    TString saved_name = USECOMBINENAMES ? all_datacard_samples_COMBINE[s] : all_datacard_samples[s];

    // get the histogram
    TH1F* h_tmp = (TH1F*)fout->Get(fname+"_"+category+"/"+saved_name);
    if(  h_tmp==0  ){
      errFlag++;
      continue;
    }

    // if TTH, add to signal
    if( string(all_datacard_samples[s].Data()).find("TTH")!=string::npos ){
      totSgn +=  h_tmp->Integral();
      if( !h_totSgn )
	h_totSgn = (TH1F*)h_tmp->Clone("h_totSgn");
      else
	h_totSgn->Add( h_tmp, 1.0);
    }

    // if not signal and not data, add to bkg
    else if( string(all_datacard_samples[s].Data()).find("data_obs")==string::npos ){
      totBkg +=  h_tmp->Integral();
      if( !h_totBkg )
	h_totBkg = (TH1F*)h_tmp->Clone("h_totBkg");
      else
	h_totBkg->Add( h_tmp, 1.0);
    }

    cout << " > " << string(all_datacard_samples[s].Data()) << ": \e[1;34m" << h_tmp->Integral() << "\e[0m" << endl;

    // add the histogram to the map
    aMap[all_datacard_samples[s]] = h_tmp;

    if(runOnData<=0 && 
       string(all_datacard_samples[s].Data()).find("data_obs")==string::npos){

      // add the histogram to the asymov dataset
      if( runOnData==0 && string(all_datacard_samples[s].Data()).find("TTH")==string::npos)
	h_data->Add( h_tmp );
      else if( runOnData<0 )
	h_data->Add( h_tmp );
      else{ /* ... */ }

      // add the histogram yield to the asymov dataset
      observation +=  h_tmp->Integral();

      cout << "   ( h_data++ : " << h_data->Integral()  << " )" << endl; 
    }

    else if(runOnData>0){
      if(string(all_datacard_samples[s].Data()).find("data_obs")!=string::npos){
	h_data->Add( h_tmp );
	observation +=  h_tmp->Integral();
      }
    }

    else{ /* ... */ }

  }

  if(errFlag>0){
    cout << "Missing histogram for a process. Return." << endl;
    //fout->Close();
    //return 1;
  }

  cout << "Purity:" << endl;
  cout << " > S/B = " << ": \e[1;34m" << (totBkg>0 ? totSgn/totBkg : -99) << "\e[0m" << endl;
  cout << " > p0  = " << poorman_significance( h_totSgn, h_totBkg) << endl;

  //cout << "************************" << endl;

  if(runOnData<=0){

    // approximate to integer precision
    observation = int(observation);
    h_data->Scale(observation/h_data->Integral());

    // save asymov data
    dir->cd();
    h_data->SetName("data_obs");
    h_data->Write("data_obs", TObject::kOverwrite);
  }

  // return if you don't want to create the datadacards
  if( CREATEDATACARDS==0 ){
    fout->Close();
    return 1;
  }



  // the text file for the datacard
  ofstream out(Form("../root/datacards/%s/%s_%s.txt",directory.Data(), (fname+"_"+category).Data(), (name+version+extraname).c_str()));
  out.precision(8);


  string longspace  = "              ";
  string shortspace = "          ";

  // if   null == "", theory unc are correlated among categories
  // else each category has its own uncorrelated unc.
  string null("");
  //string null(category.Data());

  // internal convention
  string line1("imax 1"); out << line1 << endl;
  string line2("jmax *"); out << line2 << endl;
  string line3("kmax *"); out << line3 << endl;
  string line4(Form("shapes *  *    %s_%s.root  $CHANNEL/$PROCESS $CHANNEL/$PROCESS_$SYSTEMATIC", fname.Data(), (name+version+extraname).c_str())); out << line4 << endl;
  out<<endl;
  string line5(Form("observation %.0f", observation )); out << line5 << endl;
  out<<endl;

  // check if all the histograms are there (and filled)
  bool isTTH125there     = aMap.find("TTH125")!=aMap.end()     ? aMap["TTH125"]->Integral()>0     : false;
  bool isTTJetsHFbbthere = aMap.find("TTJetsHFbb")!=aMap.end() ? aMap["TTJetsHFbb"]->Integral()>0 : false;
  bool isTTJetsHFbthere  = aMap.find("TTJetsHFb")!=aMap.end()  ? aMap["TTJetsHFb"]->Integral()>0  : false;
  bool isTTJetsLFccthere = aMap.find("TTJetsLFcc")!=aMap.end() ? aMap["TTJetsLFcc"]->Integral()>0 : false;
  bool isTTJetsLFthere   = aMap.find("TTJetsLF")!=aMap.end()   ? aMap["TTJetsLF"]->Integral()>0   : false;
  bool isTTVthere        = aMap.find("TTV")!=aMap.end()        ? aMap["TTV"]->Integral()>0        : false;
  bool isSingleTthere    = aMap.find("SingleT")!=aMap.end()    ? aMap["SingleT"]->Integral()>0    : false;
  bool isEWKthere        = aMap.find("EWK")!=aMap.end()        ? aMap["EWK"]->Integral()>0        : false;

  //////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////
  // add processes only if yield is non null

  int proc_num_offset = 0;

  // if do ttbb/ttjj analysis, signal is tt+bb
  if(isTTH125there && abs(doMEM)==3 && FITTTBB){
    isTTH125there   = false;
    proc_num_offset = 2;
  }

  string line("bin                         ");
  if( isTTH125there     ) line += string(Form("%s    ", (fname+"_"+category).Data()));
  if( isTTJetsHFbbthere ) line += string(Form("%s    ", (fname+"_"+category).Data()));
  if( isTTJetsHFbthere  ) line += string(Form("%s    ", (fname+"_"+category).Data()));
  if( isTTJetsLFccthere ) line += string(Form("%s    ", (fname+"_"+category).Data()));
  if( isTTJetsLFthere   ) line += string(Form("%s    ", (fname+"_"+category).Data()));
  if( isTTVthere        ) line += string(Form("%s    ", (fname+"_"+category).Data()));
  if( isSingleTthere    ) line += string(Form("%s    ", (fname+"_"+category).Data()));
  if( isEWKthere        ) line += string(Form("%s    ", (fname+"_"+category).Data()));
  out<<line;
  out<<endl; 
  line = "process                     ";

  if(!USECOMBINENAMES){
    if( isTTH125there )     line += "TTH125      ";
    if( isTTJetsHFbbthere ) line += "TTJetsHFbb    ";
    if( isTTJetsHFbthere )  line += "TTJetsHFb    ";
    if( isTTJetsLFccthere ) line += "TTJetsLFcc  ";
    if( isTTJetsLFthere )   line += "TTJetsLF    ";
    if( isTTVthere )        line += "TTV      ";
    if( isSingleTthere )    line += "SingleT     ";
    if( isEWKthere )        line += "EWK      ";
  }
  else{
    if( isTTH125there )     line += "ttH_hbb      ";
    if( isTTJetsHFbbthere ) line += "ttbarPlusBBbar    ";
    if( isTTJetsHFbthere )  line += "ttbarPlusB    ";
    if( isTTJetsLFccthere ) line += "ttbarPlusCCbar  ";
    if( isTTJetsLFthere )   line += "ttbar    ";
    if( isTTVthere )        line += "ttbarV      ";
    if( isSingleTthere )    line += "singlet     ";
    if( isEWKthere )        line += "ewk      ";
  }

  out<<line;
  out<<endl;

  line = "process                       ";
  int count_proc_num = 0;
  if( isTTH125there ){
    line += string(Form("%d          ", count_proc_num));
    count_proc_num++;
  }
  if( isTTJetsHFbbthere ){
    line += string(Form("%d          ", count_proc_num - proc_num_offset));
    count_proc_num++;
  }
  if( isTTJetsHFbthere ){
    line += string(Form("%d          ", count_proc_num - proc_num_offset)); 
    count_proc_num++; 
  }
  if( isTTJetsLFthere ) {
    line += string(Form("%d          ", count_proc_num  - proc_num_offset));
    count_proc_num++;
  }
  if( isTTJetsLFccthere ){
    line += string(Form("%d          ", count_proc_num  - proc_num_offset)); 
    count_proc_num++;
  }
  if( isTTVthere ) {
    line += string(Form("%d          ", count_proc_num  - proc_num_offset)); 
    count_proc_num++;
  }
  if( isSingleTthere ){
    line += string(Form("%d          ", count_proc_num  - proc_num_offset)); 
    count_proc_num++;
  }
  if( isEWKthere ){
    line += string(Form("%d          ", count_proc_num  - proc_num_offset));
  }
  out<<line;
  out<<endl;

  line = "rate                        ";
  if( isTTH125there )     line += string(Form("%.4f     ", aMap["TTH125"]->Integral() ));
  if( isTTJetsHFbbthere ) line += string(Form("%.4f     ", aMap["TTJetsHFbb"]->Integral() ));
  if( isTTJetsHFbthere )  line += string(Form("%.4f     ", aMap["TTJetsHFb"]->Integral() ));
  if( isTTJetsLFccthere ) line += string(Form("%.4f     ", aMap["TTJetsLFcc"]->Integral())); 
  if( isTTJetsLFthere )   line += string(Form("%.4f     ", aMap["TTJetsLF"]->Integral()));
  if( isTTVthere )        line += string(Form("%.4f     ", aMap["TTV"]->Integral()));
  if( isSingleTthere )    line += string(Form("%.4f     ", aMap["SingleT"]->Integral()));
  if( isEWKthere )        line += string(Form("%.4f     ", aMap["EWK"]->Integral()));
  out<<line;
  out<<endl;
  out<<endl;
  out<< "-----------------------------------------------------------------" << endl;

  //////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////
  // systematics

  // LUMINOSITY
  if( isTTH125there )     appendLogNSyst( aMap["TTH125"],     dir, "lumi_8TeV", 0.026, line);
  if( isTTJetsHFbbthere ) appendLogNSyst( aMap["TTJetsHFbb"], dir, "lumi_8TeV", 0.026, line);
  if( isTTJetsHFbthere )  appendLogNSyst( aMap["TTJetsHFb"],  dir, "lumi_8TeV", 0.026, line);
  if( isTTJetsLFccthere ) appendLogNSyst( aMap["TTJetsLFcc"], dir, "lumi_8TeV", 0.026, line);
  if( isTTJetsLFthere )   appendLogNSyst( aMap["TTJetsLF"],   dir, "lumi_8TeV", 0.026, line);
  if( isTTVthere )        appendLogNSyst( aMap["TTV"],        dir, "lumi_8TeV", 0.026, line);
  if( isSingleTthere )    appendLogNSyst( aMap["SingleT"],    dir, "lumi_8TeV", 0.026, line);
  if( isEWKthere )        appendLogNSyst( aMap["EWK"],        dir, "lumi_8TeV", 0.026, line);
  out<<line;
  out<<endl;

  // TRIGGER/RECO for leptons
  line="";
  if( isTTH125there )     appendLogNSyst( aMap["TTH125"],     dir, "CMS_ttH_eff_lep", 0.02*(1+isDL), line);
  if( isTTJetsHFbbthere ) appendLogNSyst( aMap["TTJetsHFbb"], dir, "CMS_ttH_eff_lep", 0.02*(1+isDL), line);
  if( isTTJetsHFbthere )  appendLogNSyst( aMap["TTJetsHFb"],  dir, "CMS_ttH_eff_lep", 0.02*(1+isDL), line);
  if( isTTJetsLFccthere ) appendLogNSyst( aMap["TTJetsLFcc"], dir, "CMS_ttH_eff_lep", 0.02*(1+isDL), line);
  if( isTTJetsLFthere )   appendLogNSyst( aMap["TTJetsLF"],   dir, "CMS_ttH_eff_lep", 0.02*(1+isDL), line);
  if( isTTVthere )        appendLogNSyst( aMap["TTV"],        dir, "CMS_ttH_eff_lep", 0.02*(1+isDL), line);
  if( isSingleTthere )    appendLogNSyst( aMap["SingleT"],    dir, "CMS_ttH_eff_lep", 0.02*(1+isDL), line);
  if( isEWKthere )        appendLogNSyst( aMap["EWK"],        dir, "CMS_ttH_eff_lep", 0.02*(1+isDL), line);
  out<<line;
  out<<endl;

  if( USESHIFTEDSAMPLES ){

    // BTAGGING
    if(!USECSVCALIBRATION){
      line = "csv                   shape  ";                              
      if( isTTH125there )     line += "1.0        ";       
      if( isTTJetsHFbbthere ) line += "1.0        ";       
      if( isTTJetsHFbthere )  line += "1.0        ";       
      if( isTTJetsLFccthere ) line += "1.0        ";  
      if( isTTJetsLFthere )   line += "1.0        ";            
      if( isTTVthere )        line += "1.0        ";
      if( isSingleTthere )    line += " -         ";
      if( isEWKthere )        line += " -         ";
      out<<line;
      out<<endl;
    }
    else{
      for( int sy = 0 ; sy < NCSVSYS ; sy++){

	string sys_name = csv_sys_names[sy];

	// do it only once per systematic
	if( sy%2==0 ){

	  if( sys_name.find("Up")!=string::npos )
	    sys_name.erase( sys_name.find("Up"),   2 );
	  else if ( sys_name.find("Down")!=string::npos )
	    sys_name.erase( sys_name.find("Down"), 4 );
	  else{
	    cout << "Error found for sys name " << sys_name << endl;
	    continue;
	  }
	  
	  line = "CMS_ttH_"+sys_name+"             shape  ";                              
	  if( isTTH125there )     line += "1.0        ";       
	  if( isTTJetsHFbbthere ) line += "1.0        ";       
	  if( isTTJetsHFbthere )  line += "1.0        ";       
	  if( isTTJetsLFccthere ) line += "1.0        ";       
	  if( isTTJetsLFthere )   line += "1.0        ";       
	  if( isTTVthere )        line += "1.0        ";
	  if( isSingleTthere )    line += "1.0        ";
	  if( isEWKthere )        line += "1.0        ";
	  out<<line;
	  out<<endl;
	  
	}
	
      }
    }

    // SCALE SYSTEMATICS ON TT+JETS
    for( int sy = 0 ; sy < NTHSYS ; sy++){

      string sys_name = th_sys_names[sy];

      // do it only once per systematic
      if( sy%2==0 ){

	bool isNotForHFbb = sys_name.find("ttbar_cc") != string::npos || (sys_name.find("ttbar_b")!=string::npos && sys_name.find("ttbar_bb")==string::npos );
	bool isNotForHFb  = sys_name.find("ttbar_cc") != string::npos || sys_name.find("ttbar_bb")!=string::npos;
	bool isNotForLF   = sys_name.find("ttbar_cc") != string::npos || sys_name.find("ttbar_b")!= string::npos;
	bool isNotForLFcc = sys_name.find("ttbar_b")  != string::npos ; //|| (sys_name.find("LF") != string::npos && sys_name.find("LFcc")== string::npos);
	

	if( sys_name.find("Up")!=string::npos )
	  sys_name.erase( sys_name.find("Up"),   2 );
	else if ( sys_name.find("Down")!=string::npos )
	  sys_name.erase( sys_name.find("Down"), 4 );
	else{
	  cout << "Error found for sys name " << sys_name << endl;
	  continue;
	}
	
	line = sys_name+"                 shape  ";                              
	if( isTTH125there )     line += " -         ";       
	if( isTTJetsHFbbthere ) {
	  if( !isNotForHFbb ) line += "1.0        ";  
	  else                line += " -         ";       
	}     	
	if( isTTJetsHFbthere){
	  if( !isNotForHFb ) line += "1.0        "; 
	  else               line += " -         ";         
	}
	if( isTTJetsLFccthere ) {
	  if( !isNotForLFcc ) line += "1.0        ";       
	  else                line += " -         "; 
	}
	if( isTTJetsLFthere ) {
	  if( !isNotForLF ) line += "1.0        ";       
	  else              line += " -         "; 
	}
	if( isTTVthere )        line += " -         ";
	if( isSingleTthere )    line += " -         ";
	if( isEWKthere )        line += " -         ";
	out<<line;
	out<<endl;
	
      }
      
    }


    // JET ENERGY SCALE
    line = "CMS_scale_j                   shape  ";                             
    if( isTTH125there )     line += "1.0        ";
    if( isTTJetsHFbbthere ) line += "1.0        ";
    if( isTTJetsHFbthere )  line += "1.0        ";
    if( isTTJetsLFccthere ) line += "1.0        ";
    if( isTTJetsLFthere )   line += "1.0        ";
    if( isTTVthere )        line += "1.0        ";
    if( isSingleTthere )    line += "1.0        ";
    if( isEWKthere )        line += " -         ";
    out<<line;
    out<<endl;

    // JET ENERGY RESOLUTION
    line = "CMS_res_j                   shape  ";                             
    if( isTTH125there )     line += "1.0        ";
    if( isTTJetsHFbbthere ) line += "1.0        ";
    if( isTTJetsHFbthere )  line += "1.0        ";
    if( isTTJetsLFccthere ) line += "1.0        ";
    if( isTTJetsLFthere )   line += "1.0        ";
    if( isTTVthere )        line += "1.0        ";
    if( isSingleTthere )    line += "1.0        ";
    if( isEWKthere )        line += " -         ";
    out<<line;
    out<<endl;

  }
  else{

    // BTAGGING
    line="";
    if( isTTH125there )     appendLogNSyst( aMap["TTH125"],     dir, "csv", 0.2, line);
    if( isTTJetsHFbbthere ) appendLogNSyst( aMap["TTJetsHFbb"], dir, "csv", 0.2, line);
    if( isTTJetsHFbthere )  appendLogNSyst( aMap["TTJetsHFb"],  dir, "csv", 0.2, line);
    if( isTTJetsLFccthere ) appendLogNSyst( aMap["TTJetsLFcc"], dir, "csv", 0.2, line);
    if( isTTJetsLFthere )   appendLogNSyst( aMap["TTJetsLF"],   dir, "csv", 0.2, line);
    if( isTTVthere )        appendLogNSyst( aMap["TTV"],        dir, "csv", 0.2, line);
    if( isSingleTthere )    appendLogNSyst( aMap["SingleT"],    dir, "csv", 0.2, line);
    if( isEWKthere )        appendLogNSyst( aMap["EWK"],        dir, "csv", 0.2, line);
    out<<line;
    out<<endl;

    // JET ENERGY SCALE 
    line="";
    if( isTTH125there )     appendLogNSyst( aMap["TTH125"],     dir, "JEC", 0.05, line);
    if( isTTJetsHFbbthere ) appendLogNSyst( aMap["TTJetsHFbb"], dir, "JEC", 0.05, line);
    if( isTTJetsHFbthere )  appendLogNSyst( aMap["TTJetsHFb"],  dir, "JEC", 0.05, line);
    if( isTTJetsLFccthere ) appendLogNSyst( aMap["TTJetsLFcc"], dir, "JEC", 0.05, line);
    if( isTTJetsLFthere )   appendLogNSyst( aMap["TTJetsLF"],   dir, "JEC", 0.05, line);
    if( isTTVthere )        appendLogNSyst( aMap["TTV"],        dir, "JEC", 0.05, line);
    if( isSingleTthere )    appendLogNSyst( aMap["SingleT"],    dir, "JEC", 0.05, line);
    if( isEWKthere )        appendLogNSyst( aMap["EWK"],        dir, "JEC", 0.05, line);
    out<<line;
    out<<endl;
   
  }

  // TTBB CROSS-SECTION
  line="";
  if( isTTH125there )     appendLogNSyst( aMap["TTH125"],     dir, "CMS_ttH_QCDscale_ttbb", 0.0,  line);
  if( isTTJetsHFbbthere ) appendLogNSyst( aMap["TTJetsHFbb"], dir, "CMS_ttH_QCDscale_ttbb", 0.618, line);
  if( isTTJetsHFbthere )  appendLogNSyst( aMap["TTJetsHFb"],  dir, "CMS_ttH_QCDscale_ttbb", 0.0,  line);
  if( isTTJetsLFccthere ) appendLogNSyst( aMap["TTJetsLFcc"], dir, "CMS_ttH_QCDscale_ttbb", 0.0,  line);
  if( isTTJetsLFthere )   appendLogNSyst( aMap["TTJetsLF"],   dir, "CMS_ttH_QCDscale_ttbb", 0.0,  line);
  if( isTTVthere )        appendLogNSyst( aMap["TTV"],        dir, "CMS_ttH_QCDscale_ttbb", 0.0,  line);
  if( isSingleTthere )    appendLogNSyst( aMap["SingleT"],    dir, "CMS_ttH_QCDscale_ttbb", 0.0,  line);
  if( isEWKthere )        appendLogNSyst( aMap["EWK"],        dir, "CMS_ttH_QCDscale_ttbb", 0.0,  line);
  out<<line;
  out<<endl;

  // TTB CROSS-SECTION
  line="";
  if( isTTH125there )     appendLogNSyst( aMap["TTH125"],     dir, "CMS_ttH_QCDscale_ttb", 0.0,  line);
  if( isTTJetsHFbbthere ) appendLogNSyst( aMap["TTJetsHFbb"], dir, "CMS_ttH_QCDscale_ttb", 0.0,  line);
  if( isTTJetsHFbthere )  appendLogNSyst( aMap["TTJetsHFb"],  dir, "CMS_ttH_QCDscale_ttb", 0.618,line);
  if( isTTJetsLFccthere ) appendLogNSyst( aMap["TTJetsLFcc"], dir, "CMS_ttH_QCDscale_ttb", 0.0,  line);
  if( isTTJetsLFthere )   appendLogNSyst( aMap["TTJetsLF"],   dir, "CMS_ttH_QCDscale_ttb", 0.0,  line);
  if( isTTVthere )        appendLogNSyst( aMap["TTV"],        dir, "CMS_ttH_QCDscale_ttb", 0.0,  line);
  if( isSingleTthere )    appendLogNSyst( aMap["SingleT"],    dir, "CMS_ttH_QCDscale_ttb", 0.0,  line);
  if( isEWKthere )        appendLogNSyst( aMap["EWK"],        dir, "CMS_ttH_QCDscale_ttb", 0.0,  line);
  out<<line;
  out<<endl;

  if( SPLITCC ){

    // TTCC CROSS-SECTION
    line="";
    if( isTTH125there )     appendLogNSyst( aMap["TTH125"],     dir, "CMS_ttH_QCDscale_ttcc", 0.0,  line);
    if( isTTJetsHFbbthere ) appendLogNSyst( aMap["TTJetsHFbb"], dir, "CMS_ttH_QCDscale_ttcc", 0.0,  line);
    if( isTTJetsHFbthere )  appendLogNSyst( aMap["TTJetsHFb"],  dir, "CMS_ttH_QCDscale_ttcc", 0.0,  line);
    if( isTTJetsLFccthere ) appendLogNSyst( aMap["TTJetsLFcc"], dir, "CMS_ttH_QCDscale_ttcc", 0.618,line);
    if( isTTJetsLFthere )   appendLogNSyst( aMap["TTJetsLF"],   dir, "CMS_ttH_QCDscale_ttcc", 0.0,  line);
    if( isTTVthere )        appendLogNSyst( aMap["TTV"],        dir, "CMS_ttH_QCDscale_ttcc", 0.0,  line);
    if( isSingleTthere )    appendLogNSyst( aMap["SingleT"],    dir, "CMS_ttH_QCDscale_ttcc", 0.0,  line);
    if( isEWKthere )        appendLogNSyst( aMap["EWK"],        dir, "CMS_ttH_QCDscale_ttcc", 0.0,  line);
    out<<line;
    out<<endl;
  }

  // TTV CROSS-SECTION
  line="";
  if( isTTH125there )     appendLogNSyst( aMap["TTH125"],     dir, "QCDscale_ttbarV", 0.0,  line);
  if( isTTJetsHFbbthere ) appendLogNSyst( aMap["TTJetsHFbb"], dir, "QCDscale_ttbarV", 0.0,  line);
  if( isTTJetsHFbthere )  appendLogNSyst( aMap["TTJetsHFb"],  dir, "QCDscale_ttbarV", 0.0,  line);
  if( isTTJetsLFccthere ) appendLogNSyst( aMap["TTJetsLFcc"], dir, "QCDscale_ttbarV", 0.0,  line);
  if( isTTJetsLFthere )   appendLogNSyst( aMap["TTJetsLF"],   dir, "QCDscale_ttbarV", 0.0,  line);
  if( isTTVthere )        appendLogNSyst( aMap["TTV"],        dir, "QCDscale_ttbarV", 0.15, line);
  if( isSingleTthere )    appendLogNSyst( aMap["SingleT"],    dir, "QCDscale_ttbarV", 0.0,  line);
  if( isEWKthere )        appendLogNSyst( aMap["EWK"],        dir, "QCDscale_ttbarV", 0.0,  line);
  out<<line;
  out<<endl;

  //  SINGLE-TOP CROSS-SECTION
  line="";
  if( isTTH125there )     appendLogNSyst( aMap["TTH125"],     dir, "QCDscale_singlet", 0.0,  line);
  if( isTTJetsHFbbthere ) appendLogNSyst( aMap["TTJetsHFbb"], dir, "QCDscale_singlet", 0.0,  line);
  if( isTTJetsHFbthere )  appendLogNSyst( aMap["TTJetsHFb"],  dir, "QCDscale_singlet", 0.0,  line);
  if( isTTJetsLFccthere ) appendLogNSyst( aMap["TTJetsLFcc"], dir, "QCDscale_singlet", 0.0,  line);
  if( isTTJetsLFthere )   appendLogNSyst( aMap["TTJetsLF"],   dir, "QCDscale_singlet", 0.0,  line);
  if( isTTVthere )        appendLogNSyst( aMap["TTV"],        dir, "QCDscale_singlet", 0.0,  line);
  if( isSingleTthere )    appendLogNSyst( aMap["SingleT"],    dir, "QCDscale_singlet", 0.02, line);
  if( isEWKthere )        appendLogNSyst( aMap["EWK"],        dir, "QCDscale_singlet", 0.0,  line);
  out<<line;
  out<<endl;

  // EWK CROSS-SECTION
  line="";
  if( isTTH125there )     appendLogNSyst( aMap["TTH125"],     dir, "QCDscale_ewk", 0.0,  line);
  if( isTTJetsHFbbthere ) appendLogNSyst( aMap["TTJetsHFbb"], dir, "QCDscale_ewk", 0.0,  line);
  if( isTTJetsHFbthere )  appendLogNSyst( aMap["TTJetsHFb"],  dir, "QCDscale_ewk", 0.0,  line);
  if( isTTJetsLFccthere ) appendLogNSyst( aMap["TTJetsLFcc"], dir, "QCDscale_ewk", 0.0,  line); 
  if( isTTJetsLFthere )   appendLogNSyst( aMap["TTJetsLF"],   dir, "QCDscale_ewk", 0.0,  line);
  if( isTTVthere )        appendLogNSyst( aMap["TTV"],        dir, "QCDscale_ewk", 0.0,  line);
  if( isSingleTthere )    appendLogNSyst( aMap["SingleT"],    dir, "QCDscale_ewk", 0.0,  line);
  if( isEWKthere )        appendLogNSyst( aMap["EWK"],        dir, "QCDscale_ewk", 0.013,line);
  out<<line;
  out<<endl;

  //  QCD SCALE UNC. ON TTH
  line="";
  if( isTTH125there )     appendLogNSyst( aMap["TTH125"],     dir, "QCDscale_ttH", 0.067,line);
  if( isTTJetsHFbbthere ) appendLogNSyst( aMap["TTJetsHFbb"], dir, "QCDscale_ttH", 0.0,  line);
  if( isTTJetsHFbthere )  appendLogNSyst( aMap["TTJetsHFb"],  dir, "QCDscale_ttH", 0.0,  line);
  if( isTTJetsLFccthere ) appendLogNSyst( aMap["TTJetsLFcc"], dir, "QCDscale_ttH", 0.0,  line);
  if( isTTJetsLFthere )   appendLogNSyst( aMap["TTJetsLF"],   dir, "QCDscale_ttH", 0.0,  line);
  if( isTTVthere )        appendLogNSyst( aMap["TTV"],        dir, "QCDscale_ttH", 0.0,  line);
  if( isSingleTthere )    appendLogNSyst( aMap["SingleT"],    dir, "QCDscale_ttH", 0.0,  line);
  if( isEWKthere )        appendLogNSyst( aMap["EWK"],        dir, "QCDscale_ttH", 0.0,  line);
  out<<line;
  out<<endl;

  /*
  //  QCD SCALE UNC. ON TT+HFbb
  line="";
  if( isTTH125there )     appendLogNSyst( aMap["TTH125"],     dir,  string(Form("QCDscale%s_TTJetsHFbb", null.c_str())), 0.0,  line);
  if( isTTJetsHFbbthere ) appendLogNSyst( aMap["TTJetsHFbb"], dir,  string(Form("QCDscale%s_TTJetsHFbb", null.c_str())), 0.17, line);
  if( isTTJetsHFbthere )  appendLogNSyst( aMap["TTJetsHFb"],  dir,  string(Form("QCDscale%s_TTJetsHFbb", null.c_str())), 0.0,  line);
  if( isTTJetsLFccthere ) appendLogNSyst( aMap["TTJetsLFcc"], dir,  string(Form("QCDscale%s_TTJetsHFbb", null.c_str())), 0.0,  line);
  if( isTTJetsLFthere )   appendLogNSyst( aMap["TTJetsLF"],   dir,  string(Form("QCDscale%s_TTJetsHFbb", null.c_str())), 0.0,  line);
  if( isTTVthere )        appendLogNSyst( aMap["TTV"],        dir,  string(Form("QCDscale%s_TTJetsHFbb", null.c_str())), 0.0,  line);
  if( isSingleTthere )    appendLogNSyst( aMap["SingleT"],    dir,  string(Form("QCDscale%s_TTJetsHFbb", null.c_str())), 0.0,  line);
  if( isEWKthere )        appendLogNSyst( aMap["EWK"],        dir,  string(Form("QCDscale%s_TTJetsHFbb", null.c_str())), 0.0,  line);
  out<<line;
  out<<endl;

  //  QCD SCALE UNC. ON TT+HFb
  line="";
  if( isTTH125there )     appendLogNSyst( aMap["TTH125"],     dir,  string(Form("QCDscale%s_TTJetsHFb", null.c_str())), 0.0,  line);
  if( isTTJetsHFbbthere ) appendLogNSyst( aMap["TTJetsHFbb"], dir,  string(Form("QCDscale%s_TTJetsHFb", null.c_str())), 0.0,  line);
  if( isTTJetsHFbthere )  appendLogNSyst( aMap["TTJetsHFb"],  dir,  string(Form("QCDscale%s_TTJetsHFb", null.c_str())), 0.07, line);
  if( isTTJetsLFccthere ) appendLogNSyst( aMap["TTJetsLFcc"], dir,  string(Form("QCDscale%s_TTJetsHFb", null.c_str())), 0.0,  line); 
  if( isTTJetsLFthere )   appendLogNSyst( aMap["TTJetsLF"],   dir,  string(Form("QCDscale%s_TTJetsHFb", null.c_str())), 0.0,  line);
  if( isTTVthere )        appendLogNSyst( aMap["TTV"],        dir,  string(Form("QCDscale%s_TTJetsHFb", null.c_str())), 0.0,  line);
  if( isSingleTthere )    appendLogNSyst( aMap["SingleT"],    dir,  string(Form("QCDscale%s_TTJetsHFb", null.c_str())), 0.0,  line);
  if( isEWKthere )        appendLogNSyst( aMap["EWK"],        dir,  string(Form("QCDscale%s_TTJetsHFb", null.c_str())), 0.0,  line);
  out<<line;
  out<<endl;

  if( SPLITCC ){
    //  QCD SCALE UNC. ON TT+LFcc
    line="";
    if( isTTH125there )     appendLogNSyst( aMap["TTH125"],     dir,  string(Form("QCDscale%s_TTJetsLFcc", null.c_str())), 0.0,  line);
    if( isTTJetsHFbbthere ) appendLogNSyst( aMap["TTJetsHFbb"], dir,  string(Form("QCDscale%s_TTJetsLFcc", null.c_str())), 0.0,  line);
    if( isTTJetsHFbthere )  appendLogNSyst( aMap["TTJetsHFb"],  dir,  string(Form("QCDscale%s_TTJetsLFcc", null.c_str())), 0.0,  line);
    if( isTTJetsLFccthere ) appendLogNSyst( aMap["TTJetsLFcc"], dir,  string(Form("QCDscale%s_TTJetsLFcc", null.c_str())), 0.11, line);
    if( isTTJetsLFthere )   appendLogNSyst( aMap["TTJetsLF"],   dir,  string(Form("QCDscale%s_TTJetsLFcc", null.c_str())), 0.0,  line);
    if( isTTVthere )        appendLogNSyst( aMap["TTV"],        dir,  string(Form("QCDscale%s_TTJetsLFcc", null.c_str())), 0.0,  line);
    if( isSingleTthere )    appendLogNSyst( aMap["SingleT"],    dir,  string(Form("QCDscale%s_TTJetsLFcc", null.c_str())), 0.0,  line);
    if( isEWKthere )        appendLogNSyst( aMap["EWK"],        dir,  string(Form("QCDscale%s_TTJetsLFcc", null.c_str())), 0.0,  line);
    out<<line;
    out<<endl;
  }
  */

  //  QCD SCALE UNC. ON TT+LF
  line="";
  if( isTTH125there )     appendLogNSyst( aMap["TTH125"],     dir,  string(Form("QCDscale%s_ttbar", null.c_str())), 0.0,  line);
  if( isTTJetsHFbbthere ) appendLogNSyst( aMap["TTJetsHFbb"], dir,  string(Form("QCDscale%s_ttbar", null.c_str())), 0.03, line);
  if( isTTJetsHFbthere )  appendLogNSyst( aMap["TTJetsHFb"],  dir,  string(Form("QCDscale%s_ttbar", null.c_str())), 0.03, line);
  if( isTTJetsLFccthere ) appendLogNSyst( aMap["TTJetsLFcc"], dir,  string(Form("QCDscale%s_ttbar", null.c_str())), 0.03, line);
  if( isTTJetsLFthere )   appendLogNSyst( aMap["TTJetsLF"],   dir,  string(Form("QCDscale%s_ttbar", null.c_str())), 0.03, line);
  if( isTTVthere )        appendLogNSyst( aMap["TTV"],        dir,  string(Form("QCDscale%s_ttbar", null.c_str())), 0.0,  line);
  if( isSingleTthere )    appendLogNSyst( aMap["SingleT"],    dir,  string(Form("QCDscale%s_ttbar", null.c_str())), 0.0,  line);
  if( isEWKthere )        appendLogNSyst( aMap["EWK"],        dir,  string(Form("QCDscale%s_ttbar", null.c_str())), 0.0,  line);
  out<<line;
  out<<endl;
  
  // PDFs gg 
  line="";
  if( isTTH125there )     appendLogNSyst( aMap["TTH125"],     dir, "pdf_gg",-0.077,line);
  if( isTTJetsHFbbthere ) appendLogNSyst( aMap["TTJetsHFbb"], dir, "pdf_gg", 0.00, line);
  if( isTTJetsHFbthere )  appendLogNSyst( aMap["TTJetsHFb"],  dir, "pdf_gg", 0.00 ,line);
  if( isTTJetsLFccthere ) appendLogNSyst( aMap["TTJetsLFcc"], dir, "pdf_gg", 0.00 ,line);
  if( isTTJetsLFthere )   appendLogNSyst( aMap["TTJetsLF"],   dir, "pdf_gg", 0.00 ,line);
  if( isTTVthere )        appendLogNSyst( aMap["TTV"],        dir, "pdf_gg",-0.077,line);
  if( isSingleTthere )    appendLogNSyst( aMap["SingleT"],    dir, "pdf_gg", 0.00, line);
  if( isEWKthere )        appendLogNSyst( aMap["EWK"],        dir, "pdf_gg", 0.00, line);
  out<<line;
  out<<endl;

  // PDFs gg for tt+jets
  line="";
  if( isTTH125there )     appendLogNSyst( aMap["TTH125"],     dir, "CMS_ttHbb_pdf_ttbar", 0.00,  line);
  if( isTTJetsHFbbthere ) appendLogNSyst( aMap["TTJetsHFbb"], dir, "CMS_ttHbb_pdf_ttbar",-0.026, line);
  if( isTTJetsHFbthere )  appendLogNSyst( aMap["TTJetsHFb"],  dir, "CMS_ttHbb_pdf_ttbar",-0.026 ,line);
  if( isTTJetsLFccthere ) appendLogNSyst( aMap["TTJetsLFcc"], dir, "CMS_ttHbb_pdf_ttbar",-0.026 ,line);
  if( isTTJetsLFthere )   appendLogNSyst( aMap["TTJetsLF"],   dir, "CMS_ttHbb_pdf_ttbar",-0.026 ,line);
  if( isTTVthere )        appendLogNSyst( aMap["TTV"],        dir, "CMS_ttHbb_pdf_ttbar", 0.0,   line);
  if( isSingleTthere )    appendLogNSyst( aMap["SingleT"],    dir, "CMS_ttHbb_pdf_ttbar", 0.00,  line);
  if( isEWKthere )        appendLogNSyst( aMap["EWK"],        dir, "CMS_ttHbb_pdf_ttbar", 0.00,  line);
  out<<line;
  out<<endl;
  
  // PDFs qq
  line="";
  if( isTTH125there )     appendLogNSyst( aMap["TTH125"],     dir, "pdf_qqbar", 0.00, line);
  if( isTTJetsHFbbthere ) appendLogNSyst( aMap["TTJetsHFbb"], dir, "pdf_qqbar", 0.00, line);
  if( isTTJetsHFbthere )  appendLogNSyst( aMap["TTJetsHFb"],  dir, "pdf_qqbar", 0.00, line);
  if( isTTJetsLFccthere ) appendLogNSyst( aMap["TTJetsLFcc"], dir, "pdf_qqbar", 0.00, line);
  if( isTTJetsLFthere )   appendLogNSyst( aMap["TTJetsLF"],   dir, "pdf_qqbar", 0.00, line);
  if( isTTVthere )        appendLogNSyst( aMap["TTV"],        dir, "pdf_qqbar", 0.00, line);
  if( isSingleTthere )    appendLogNSyst( aMap["SingleT"],    dir, "pdf_qqbar", 0.00, line);
  if( isEWKthere )        appendLogNSyst( aMap["EWK"],        dir, "pdf_qqbar", 0.048,line);
  out<<line;
  out<<endl;
  
  // PDFs gq
  line="";
  if( isTTH125there )     appendLogNSyst( aMap["TTH125"],     dir, "pdf_qg", 0.00, line);
  if( isTTJetsHFbbthere ) appendLogNSyst( aMap["TTJetsHFbb"], dir, "pdf_qg", 0.00, line);
  if( isTTJetsHFbthere )  appendLogNSyst( aMap["TTJetsHFb"],  dir, "pdf_qg", 0.00, line);
  if( isTTJetsLFccthere ) appendLogNSyst( aMap["TTJetsLFcc"], dir, "pdf_qg", 0.00, line);
  if( isTTJetsLFthere )   appendLogNSyst( aMap["TTJetsLF"],   dir, "pdf_qg", 0.00, line);
  if( isTTVthere )        appendLogNSyst( aMap["TTV"],        dir, "pdf_qg", 0.00, line);
  if( isSingleTthere )    appendLogNSyst( aMap["SingleT"],    dir, "pdf_qg", 0.047,line);
  if( isEWKthere )        appendLogNSyst( aMap["EWK"],        dir, "pdf_qg", 0.00, line);
  out<<line;
  out<<endl;

  //////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////
  // bin-by-bin
  if( isTTH125there ){
    string shapename = USECOMBINENAMES ? "ttH_hbb" : "TTH125";
    for(int bin = 1; bin<=aMap["TTH125"]->GetNbinsX(); bin++ ){
      if(aMap["TTH125"]->GetBinContent(bin)>0){
	line = string(Form("%s%sbin%d        shape  ", shapename.c_str(), category.Data(), bin));
	if( isTTH125there )     line += "1.0        ";
	if( isTTJetsHFbbthere ) line += " -         ";
	if( isTTJetsHFbthere )  line += " -         ";
	if( isTTJetsLFccthere ) line += " -         ";
	if( isTTJetsLFthere )   line += " -         ";
	if( isTTVthere )        line += " -         ";
	if( isSingleTthere )    line += " -         ";
	if( isEWKthere )        line += " -         ";
	out<<line;
	out<<endl;
      }
    }  
  }
  if( isTTJetsHFbbthere ){
    string shapename = USECOMBINENAMES ? "ttbarPlusBBbar" : "TTJetsHFbb";	
    for(int bin = 1; bin<=aMap["TTJetsHFbb"]->GetNbinsX(); bin++ ){
      if(aMap["TTJetsHFbb"]->GetBinContent(bin)>0){
	line = string(Form("%s%sbin%d      shape  ",shapename.c_str(), category.Data(), bin));
	if( isTTH125there )     line += " -         ";
	if( isTTJetsHFbbthere ) line += "1.0        ";
	if( isTTJetsHFbthere )  line += " -        ";
	if( isTTJetsLFccthere ) line += " -         ";
	if( isTTJetsLFthere )   line += " -         ";
	if( isTTVthere )        line += " -         ";
	if( isSingleTthere )    line += " -         ";
	if( isEWKthere )        line += " -         ";
	out<<line;
	out<<endl;
      }
    }  
  }
  if( isTTJetsHFbthere ){
    string shapename = USECOMBINENAMES ? "ttbarPlusB" : "TTJetsHFb";	
    for(int bin = 1; bin<=aMap["TTJetsHFb"]->GetNbinsX(); bin++ ){
      if(aMap["TTJetsHFb"]->GetBinContent(bin)>0){
	line = string(Form("%s%sbin%d      shape  ",shapename.c_str(), category.Data(), bin));
	if( isTTH125there )     line += " -         ";
	if( isTTJetsHFbbthere ) line += " -        ";
	if( isTTJetsHFbthere )  line += "1.0        ";
	if( isTTJetsLFccthere ) line += " -         ";
	if( isTTJetsLFthere )   line += " -         ";
	if( isTTVthere )        line += " -         ";
	if( isSingleTthere )    line += " -         ";
	if( isEWKthere )        line += " -         ";
	out<<line;
	out<<endl;
      }
    }  
  }
  if( isTTJetsLFccthere ){
    string shapename = USECOMBINENAMES ? "ttbarPlusCCbar" : "TTJetsLFcc";	
    for(int bin = 1; bin<=aMap["TTJetsLFcc"]->GetNbinsX(); bin++ ){
      if(aMap["TTJetsLFcc"]->GetBinContent(bin)>0){
	line = string(Form("%s%sbin%d      shape  ",shapename.c_str(), category.Data(), bin));
	if( isTTH125there )     line += " -         ";
	if( isTTJetsHFbbthere ) line += " -         ";
	if( isTTJetsHFbthere )  line += " -         ";
	if( isTTJetsLFccthere ) line += "1.0        ";
	if( isTTJetsLFthere )   line += " -         ";
	if( isTTVthere )        line += " -         ";
	if( isSingleTthere )    line += " -         ";
	if( isEWKthere )        line += " -         ";
	out<<line;
	out<<endl;
      }
    } 
  } 
  if( isTTJetsLFthere ){
    string shapename = USECOMBINENAMES ? "ttbar" : "TTJetsLF";	
    for(int bin = 1; bin<=aMap["TTJetsLF"]->GetNbinsX(); bin++ ){
      if(aMap["TTJetsLF"]->GetBinContent(bin)>0){
	line = string(Form("%s%sbin%d      shape  ",shapename.c_str(), category.Data(), bin));
	if( isTTH125there )     line += " -         ";
	if( isTTJetsHFbbthere ) line += " -         ";
	if( isTTJetsHFbthere )  line += " -         ";
	if( isTTJetsLFccthere ) line += " -         ";
	if( isTTJetsLFthere )   line += "1.0        ";
	if( isTTVthere )        line += " -         ";
	if( isSingleTthere )    line += " -         ";
	if( isEWKthere )        line += " -         ";
	out<<line;
	out<<endl;
      }
    } 
  } 
  if( isTTVthere ){
    string shapename = USECOMBINENAMES ? "ttbarV" : "TTV";	
    for(int bin = 1; bin<=aMap["TTV"]->GetNbinsX(); bin++ ){
      if(aMap["TTV"]->GetBinContent(bin)>0){
	line = string(Form("%s%sbin%d           shape  ",shapename.c_str(), category.Data(), bin));
	if( isTTH125there )      line += " -         ";
	if( isTTJetsHFbbthere )  line += " -         ";
	if( isTTJetsHFbthere )   line += " -         ";
	if( isTTJetsLFccthere )  line += " -         ";
	if( isTTJetsLFthere )    line += " -         ";
	if( isTTVthere )         line += "1.0        ";
	if( isSingleTthere )     line += " -         ";
	if( isEWKthere )         line += " -         ";
	out<<line;
	out<<endl;
      }
    } 
  }
  if( isSingleTthere ){
    string shapename = USECOMBINENAMES ? "singlet" : "SingleT";	
    for(int bin = 1; bin<=aMap["SingleT"]->GetNbinsX(); bin++ ){
      if(aMap["SingleT"]->GetBinContent(bin)>0){
	line = string(Form("%s%sbin%d       shape  ",shapename.c_str(), category.Data(), bin));
	if( isTTH125there )     line += " -         ";
	if( isTTJetsHFbbthere ) line += " -         ";
	if( isTTJetsHFbthere )  line += " -         ";
	if( isTTJetsLFccthere ) line += " -         ";
	if( isTTJetsLFthere )   line += " -         ";
	if( isTTVthere )        line += " -         ";
	if( isSingleTthere )    line += "1.0        ";
	if( isEWKthere )        line += " -         ";
	out<<line;
	out<<endl;
      }
    } 
  }
  if( isEWKthere ){
    string shapename = USECOMBINENAMES ? "ewk" : "EWK";	
    for(int bin = 1; bin<=aMap["EWK"]->GetNbinsX(); bin++ ){
      if(aMap["EWK"]->GetBinContent(bin)>0){
	line = string(Form("%s%sbin%d           shape  ",shapename.c_str(), category.Data(), bin));
	if( isTTH125there )      line += " -         ";
	if( isTTJetsHFbbthere )  line += " -         ";
	if( isTTJetsHFbthere )   line += " -         ";
	if( isTTJetsLFccthere )  line += " -         ";
	if( isTTJetsLFthere )    line += " -         ";
	if( isTTVthere )         line += " -         ";
	if( isSingleTthere )     line += " -         ";
	if( isEWKthere )         line += "1.0        ";
	out<<line;
	out<<endl;
      }
    } 
  }

  out<<endl;
  out<< "-----------------------------------------------------------------" << endl;

  //////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////// 

  // a reminder...
  cout << "Recap: " << endl;
  cout << " > Tag = " << name << endl;
  cout << " > Version = " <<  version << endl;
  cout << " > Extra Name = " << extraname << endl;
  cout << " > Inputpath = " << string(inputpath.Data()) << endl;
  cout << " > FOM " << doMEM << endl;
  cout << " > Category = " << category << endl;
  cout << " > Cut = (" << cut  << ")" << endl;
  cout << " > Histo = h(" << h->GetNbinsX() << "," << h->GetXaxis()->GetXmin() << "," <<  h->GetXaxis()->GetXmax() << ")" << endl;
  cout << " > Lumi = " << lumiScale*12.1 << " fb-1" << endl ;
  cout << " > S_bb/B_bb = " << S_bb/B_bb << ", S_bb/B_jj = " << S_bb/B_jj << ", B_bb/B_jj = " << B_bb/B_jj << endl << endl;

  // close the output file
  delete h;
  if(corrector) delete corrector;
  fout->Close();
  //delete fout;

  return 0;
}
