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

typedef TMatrixT<double> TMatrixD;

#define CREATEDATACARDS 1

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


void bbb( TH1F* hin, TH1F* hout_Down, TH1F* hout_Up, int bin){

  // get the bin content shifted down and up by 1sigma
  float bin_down = hin->GetBinContent( bin ) - hin->GetBinError( bin );
  float bin_up   = hin->GetBinContent( bin ) + hin->GetBinError( bin );

  // just copy the nominal histogram...
  hout_Down->Reset();
  hout_Down->Add(hin,1.0);

  // ... and then shift down the bin content by -1sigma
  hout_Down->SetBinContent(bin, TMath::Max(bin_down,0.01));

  // just copy the nominal histogram...
  hout_Up->Reset();
  hout_Up->Add(hin,1.0);

  // ... and then shift up the bin content by +1sigma
  hout_Up->SetBinContent  (bin, TMath::Max(bin_up,0.01));

  return;

}


void draw(vector<float> param, TTree* t = 0, TString var = "", TH1F* h = 0, TCut cut = ""){

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
    t->Draw(var2+">>"+TString(h->GetName()),   "weight"*(cut&&cut1&&cut2) );
    float binL    = h->Integral();
    float binLErr = h->GetEntries()>0 ? sqrt(h->GetEntries())*h->Integral()/h->GetEntries() : 0.;
    h->Reset();
    t->Draw(var2+">>"+TString(h->GetName()),   "weight"*(cut&&cut1&&(!cut2)) );
    float binH    = h->Integral();
    float binHErr = h->GetEntries()>0 ? sqrt(h->GetEntries())*h->Integral()/h->GetEntries() : 0.;
    h->Reset();

    t->Draw(var+">>"+TString(h->GetName()),   "weight"*cut );
    h->SetBinContent(1, binL);  h->SetBinError(1, binLErr);
    h->SetBinContent(2, binH);  h->SetBinError(2, binHErr);
  }

  return;
}


void fill(  TTree* tFull = 0, TH1F* h = 0, TCut cut = "" , int analysis = 0, vector<float>* param = 0, TF1* xsec = 0){

  // needed as usual to copy a tree
  TFile* dummy = new TFile("dummy.root","RECREATE");

  // tree of events in the category
  TTree* t = (TTree*)tFull->CopyTree( cut );

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

  // the event-dependent weight
  float weight;
  
  t->SetBranchAddress("p_vsMH_s",     p_vsMH_s);
  t->SetBranchAddress("p_vsMT_b",     p_vsMT_b);
  t->SetBranchAddress("p_tt_bb",      p_tt_bb);
  t->SetBranchAddress("p_tt_jj",      p_tt_jj);
  t->SetBranchAddress("nPermut_s",    &nPermut_s);
  t->SetBranchAddress("nPermut_b",    &nPermut_b);
  t->SetBranchAddress("nMassPoints",  &nMassPoints);
  t->SetBranchAddress("mH_scan",      mH_scan);
  t->SetBranchAddress("mT_scan",      mT_scan);
  t->SetBranchAddress("weight",       &weight);
  
  Long64_t nentries = t->GetEntries(); 
  cout << "Total entries: " << nentries << endl;
  
  // a temporary histogram that contains the prob. vs mass
  TH1F* hTmp      = new TH1F("hTmp", "", 500,0,500);

  // loop over tree entries
  for (Long64_t i = 0; i < nentries ; i++){

    t->GetEntry(i);
    
    // reset (needed because hTmp is filled with Fill()
    hTmp->Reset();

    // if doing a mass analysis... ( 0 = scan over mH ; -1 = scan over mT )
    if( analysis<=0 ){


      // choose what scan to perform
      int nPermut      = analysis==0 ? nPermut_s     : nPermut_b;
      float* m_scan    = analysis==0 ? mH_scan       : mT_scan;
      float* p_vsM     = analysis==0 ? p_vsMH_s      : p_vsMT_b;

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
	double p =  ME_prob*bb_prob*norm;       

	// fill the histogram
	hTmp->Fill( mH, p );
      }
      
      // now that the histogram is filled, find the global maximum
      pair<double,double> bestMass = getMaxValue(hTmp);

      // this comes from a quadratic interpolation around the maximum
      double mass = bestMass.first;

      // fill the output histo with the interpolated mass value
      h->Fill(mass, weight);

    }

    // if doing a ME analysis...
    else if( analysis==1 ){

      // this number are the total probability under the ttH, ttbb, and ttjj hypotheses
      double p_125_all_s_ttbb = 0.;
      double p_125_all_b_ttbb = 0.;
      double p_125_all_b_ttjj = 0.;

      // consider only the mass scan at the nominal Higgs mass
      for(int mH_it = 0; mH_it<nMassPoints; mH_it++){	

	// skip other mass values...
	float mH =  mH_scan[mH_it];
	if(mH>126 || mH<124) continue;
	
	// once we got the good Higgs mass, loop over permutations...
	for( int perm_it = 0 ; perm_it<nPermut_s ; perm_it++){	
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
      
      // if not splitting the first bin, just fill with the weight
      if( param->size()==3 ) h->Fill( w, weight);
      
      else{
	
	// sanity check
	if( param->size()!= 5 ){
	  cout << "Problem in fill..." << endl;
	  return;
	}
	
	// if w<cutval, do a 2D analyses LR_sb vs LR_bj
	float cutval = h->GetBinLowEdge( 3 );
	
	// prepare the bj likelihood ratio ingredients
	double wbjDen = p_125_all_b_ttbb + (*param)[3]*p_125_all_b_ttjj;
	double wbjNum = p_125_all_b_ttbb ;
	
	// the bj likelihood ratio
	double wbj = wbjDen>0 ? wbjNum/wbjDen : 0.;
	
	// is above the cut for 2D analysis, just do 1D analysis
	if( w >= cutval )  h->Fill( w, weight);

	// else fill out two bins...
	else{
	  // if LR_bj is below the cut value param[4], fill the first bin...
	  if( wbj <= (*param)[4] )  h->Fill(   cutval/4. , weight);
	  // if LR_bj is above the cut value param[4], fill the second first...
	  if( wbj >  (*param)[4] )  h->Fill( 3*cutval/4. , weight);
	}       		
      }
      
    }
    else{ /* NOT IMPLEMENTED */ }

  } // entries

  return;

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////


// code that produces the datacards  for limit computation 

void produceNew(// main name of the trees
		string name = "New",
		
		// version of the input trees
		string version = "_v2_reg",

		// extra name (appended at the end of the .root/.txt out files
		string extraname = "",

		// secondary name of the input trees
		TString fname = "SL", 

		// selection cut
		string cut = "type==0", 

		// a header for the output file name
		TString category = "cat0", 

		// if 1, run the ME analysis, if 0 run the mass analysis
		int doMEM = 1,

		// factors that modify the likelihood ratio 
		float fact1 = 0.02, float fact2 = 0., float factbb = 0.5,

		// luminosity normalization (if 1, samples normalized to 12.1 fb-1)
		float lumiScale = 20./12.1,

		// number of bins in the plot
		int nBins = 6 , 

		// wheteher to split the first bin into ttbb/ttjj dominated sidebands
		int splitFirstBin=0, 

		// vector of bin edges
		vector<float>* binvec = 0,

		// TF1 to normalize to xsec
		TF1* xsec = 0
		){


  // cout the version
  cout << "Use trees " << name << " (version " << version << "): category " << category << endl;

  // input files path
  TString inputpath = "gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/Trees/MEM/Nov21_2013/BTagLoose/";
  //TString inputpath = "./";

  // selection cut
  string basecut = cut;
  
  // name of the main background (needed because one could use TTJetsSemiLept, TTJetsFullLept, ...)
  TString ttjets = "_TTJets";

  // strength modifiers 
  float scaleTTH      = 1.0;//4.7;
  float scaleTTJetsLF = 1.0;//3.9;
  float scaleTTJetsHF = 1.0;//3.9;

  cout << "*****************************************************" << endl;
  cout << "Attention!!! The following rescaling will be applied:" << endl;
  cout << "TTH:       " << scaleTTH << endl;
  cout << "TTJets LF: " << scaleTTJetsLF << endl;
  cout << "TTJets HF: " << scaleTTJetsHF << endl;
  cout << "*****************************************************" << endl;


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
  
  if(doMEM){

    TFile* f_B = TFile::Open(inputpath+"MEAnalysis"+name+"_"+fname+"_nominal"+version+ttjets+".root", "OPEN");
    if( f_B==0 || f_B->IsZombie() ){
      cout << "Could not find f_B" << endl;
      return;
    }
    TTree* tB = (TTree*)f_B->Get("tree");

    TFile* f_S = TFile::Open(inputpath+"MEAnalysis"+name+"_"+fname+"_nominal"+version+"_TTH125.root");
    if( f_S==0 || f_S->IsZombie() ){
      cout << "Could not find f_S" << endl;
      return;
    }
    TTree* tS = (TTree*)f_S->Get("tree");
    
    TH1F* hInt = new TH1F("hInt","",200,-100,100);  
    
    tS->Draw("log(p_125_all_s_ttbb)>>hInt", TCut((basecut+" && nSimBs>=2 && p_125_all_s_ttbb>0").c_str())); 
    S_bb     = TMath::Exp(hInt->GetMean());
    Int_S_bb = hInt->Integral();
    cout << "S_bb = " << S_bb <<endl;
    hInt->Reset();
    
    tB->Draw("log(p_125_all_b_ttbb)>>hInt", TCut((basecut+" && nMatchSimBs>=2 && nSimBs>2 && p_125_all_b_ttbb>0").c_str())); 
    B_bb     = TMath::Exp(hInt->GetMean());
    Int_B_bb = hInt->Integral();
    cout << "B_bb = " << B_bb <<endl;
    hInt->Reset();
    
    tB->Draw("log(p_125_all_b_ttjj)>>hInt", TCut((basecut+" && nMatchSimBs<2  && nSimBs==2 && p_125_all_b_ttjj>0").c_str())); 
    B_jj     = TMath::Exp(hInt->GetMean());
    Int_B_jj = hInt->Integral();
    cout << "B_jj = " << B_jj <<endl;
    hInt->Reset();
    
    cout << "Ratio S_bb/B_bb=" << S_bb/B_bb << endl;
    cout << "Ratio S_bb/B_jj=" << S_bb/B_jj << endl;
    cout << "Ratio bb/(bb+jj)=" << Int_B_bb << "/(" << Int_B_bb << "+" << Int_B_jj << ")=" << Int_B_bb/(Int_B_bb+Int_B_jj) << endl;       

    param.push_back( fact1 );
    param.push_back( (1-fact2)*S_bb/B_bb*(Int_B_bb/(Int_B_bb+Int_B_jj)) );
    param.push_back( (1+fact2)*S_bb/B_jj*(Int_B_jj/(Int_B_bb+Int_B_jj)) );

    if(splitFirstBin){
      param.push_back( B_bb/B_jj*factbb );
      param.push_back( 0.50 );
    }

    delete hInt;
  }
  else{
    param.clear();
  }

  // if true, read binning from input, or just take unformly distributed bins [0,1]
  bool normalbinning = ( param.size()==3 || doMEM==0 );

  // output file with shapes
  TFile* fout = new TFile("datacards/"+fname+"_"+name+version+extraname+".root","UPDATE");

  // directory for this particular category
  TDirectory* dir =  fout->GetDirectory( fname+"_"+category); 
  if( !dir) dir = fout->mkdir( fname+"_"+category ) ;

  // fix the binning
  TArrayF bins (nBins+1);
  TArrayF bins2(nBins+2);
  cout << "Making histograms with " << nBins << " bins:" << endl;

  // if do not split, then...
  if( normalbinning ){

    // if ME analysis, just take nBins fix-size bins in [0,1]
    if(doMEM){
      for(int b = 0; b < nBins+1; b++)
	bins[b] = b*1.0/nBins;
    }
    // else, read from binvec what are the axis bins
    else if(binvec->size()==(nBins+1)){
      for(int b = 0; b < nBins+1; b++)
	bins[b] = (*binvec)[b];
    }
    // sanity check
    else{
      cout << "WATCH OUT!! Mismatch nBins/binvec" << endl;
      return;
    }
  }
  // else, split the first bin into two...
  else{
    for(int b = 0; b < nBins+2; b++){
      if(b<=2) bins2[b] = b*0.5/(nBins);
      else     bins2[b] = (b-1)*1.0/(nBins);
      cout <<  bins2[b] << ", ";
    }
    cout << endl;
  }

  // master histogram (all other histograms are a clone of this one)
  TH1F* h = new TH1F("h","Simulation #sqrt{s}=8 TeV, "+fname+"; S/(S+B); units", 
		     normalbinning ? nBins : nBins+1 ,  
		     normalbinning ? bins.GetArray() : bins2.GetArray());
 
  // the observable when doing the ME analysis for cat2 (workaround)
  TString var("");
  var = TString(Form("p_125_all_s_ttbb/(p_125_all_s_ttbb+%f*(%f*p_125_all_b_ttbb+%f*p_125_all_b_ttjj))", 
		     fact1, (1-fact2)*S_bb/B_bb*(Int_B_bb/(Int_B_bb+Int_B_jj)), (1+fact2)*S_bb/B_jj*(Int_B_jj/(Int_B_bb+Int_B_jj)) ));

  if(doMEM)
    cout << "Variable = " << string(var.Data()) << endl;
  else{
    cout << "Variable = Higgs mass estimator" << endl;
  }

  // list of MC processes that will be considered 
  // (N.B. the name must match the name of the input trees)
  vector<string> samples;
  samples.push_back("TTV");
  samples.push_back("SingleT");
  samples.push_back("DiBoson");
  samples.push_back("TTJetsBB");
  samples.push_back("TTJetsBJ");
  samples.push_back("TTJetsJJ");
  samples.push_back("TTH125");
  samples.push_back("EWK");

  // list of names that will appear in datacards 
  // (N.B. same order as in samples)
  vector<TString> datacard_samples;

  // loop over the processes
  for( unsigned int s = 0 ; s < samples.size(); s++){
    
    string sample = samples[s]; 

    // this cut selects the category   
    TCut sample_cut(cut.c_str());

    string datacard_name = "";
    if( sample.find("TTH")!=string::npos ){
      datacard_name = "TTH125";
    }
    if( sample.find("TTV")!=string::npos ){
      datacard_name = "TTV";
    }
    if( sample.find("DiBoson")!=string::npos ){
      datacard_name = "DiBoson";
    }
    if( sample.find("SingleT")!=string::npos ){
      datacard_name = "SingleT";
    }
    if( sample.find("EWK")!=string::npos ){
      datacard_name = "EWK";
    }    
    if( sample.find("TTJetsBB")!=string::npos ){
      sample        = "TTJets";
      sample_cut    = sample_cut && TCut("nSimBs>2 && nMatchSimBs>=2");
      datacard_name = "TTJetsHFbb";
    }
    if( sample.find("TTJetsBJ")!=string::npos ){
      sample        = "TTJets";
      sample_cut    = sample_cut && TCut("nSimBs>2 && nMatchSimBs<2");
      datacard_name = "TTJetsHFb";
    }
    if( sample.find("TTJetsJJ")!=string::npos ){
      sample        = "TTJets";
      sample_cut    = sample_cut && TCut("nSimBs==2");
      datacard_name = "TTJetsLF";
    }
    datacard_samples.push_back( TString(datacard_name.c_str()) );

    // list of samples for systematics
    vector<string> systematics;

    // provided only for tt+jets and ttH at the moment
    if(  sample.find("TTH")!=string::npos ||  sample.find("TTJets")!=string::npos ){
      systematics.push_back("nominal");
      systematics.push_back("csvUp");
      systematics.push_back("csvDown");
      systematics.push_back("JECUp");
      systematics.push_back("JECDown");
    }
    else{
      systematics.push_back("nominal");
    }

    // loop over the systematics
    for( unsigned int sy = 0 ; sy < systematics.size(); sy++){
      
      string sys = systematics[sy];

      // input file
      string inputfilename = string(inputpath.Data())+"MEAnalysis"+name+"_"+string(fname.Data())+"_"+sys+""+version+"_"+sample+".root";
      TFile* f = TFile::Open(inputfilename.c_str());
      
      if(f==0 || f->IsZombie()){
	cout << "Missing " << sample << " file" << endl;
	continue;
      }
      else{
	cout << "Doing sample " << sample << " and analysis " << sys <<  endl;
      }
      
      // the input tree
      TTree* tree = (TTree*)f->Get("tree");
      
      // histogram for the particluar process/systematics
      TH1F* h_tmp = (TH1F*)h->Clone(("h_"+datacard_name).c_str());
      h_tmp->Reset();
      h_tmp->Sumw2();


      if(doMEM){
	// workaround due to bug in the trees on cat. 2
	//if( (string(category.Data())).find("cat2") != string::npos )
	//draw( param, tree, var, h_tmp, sample_cut );
	//else
	fill( tree, h_tmp, sample_cut, 1, &param, xsec);
      }
      else{
	fill( tree, h_tmp, sample_cut, 0, &param, xsec);

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
      h_tmp->Scale(lumiScale);

      // cd to the directory and save there
      dir->cd();

      // if doing nominal analysis, save histogram...
      if( sys.find("nominal")!=string::npos ){
	h_tmp->Write(datacard_name.c_str(), TObject::kOverwrite);
      }

      // if doing systematic analysis, save histogram with correct label and continue
      else{
	h_tmp->Write((datacard_name+"_"+sys).c_str(), TObject::kOverwrite);
	continue;      
      }
      

      // add bin-by-bin systematics
      for(int bin = 1; bin<=h_tmp->GetNbinsX(); bin++ ){

	TH1F* h_tmp_b_up   = (TH1F*)h->Clone(Form("h_tmp_%d_Up",bin));
	TH1F* h_tmp_b_down = (TH1F*)h->Clone(Form("h_tmp_%d_Down",bin));

	// the actual function that produces the shifted templates
	bbb( h_tmp, h_tmp_b_up, h_tmp_b_down, bin);      
	dir->cd();

	// save only if the bin is filled
	if(h_tmp->GetBinContent(bin)>0) h_tmp_b_up  ->Write(Form("%s_%s%sbin%dUp",   datacard_name.c_str(), datacard_name.c_str(), category.Data(), bin), TObject::kOverwrite);
	if(h_tmp->GetBinContent(bin)>0) h_tmp_b_down->Write(Form("%s_%s%sbin%dDown", datacard_name.c_str(), datacard_name.c_str(), category.Data(), bin), TObject::kOverwrite);      
      }
      
      // close the input file
      f->Close();
      
    }    
    
  }

  // for later use, use a map between process name and histograms
  std::map<TString, TH1F*> aMap;

  // for the moment, use th sum of backgrounds as 'data' (asymov dataset)
  TH1F* h_data = (TH1F*)h->Clone("h_data");
  h_data->Reset();

  // observation is the sum of all processes yields
  float observation = 0.;

  // keep track of missing files
  int errFlag = 0 ;

  cout << "************************" << endl;
  for( unsigned int s = 0 ; s < datacard_samples.size(); s++){

    // get the histogram
    TH1F* h_tmp = (TH1F*)fout->Get(fname+"_"+category+"/"+datacard_samples[s]);
    if(  h_tmp==0 ){
      errFlag++;
      continue;
    }
    cout << string(datacard_samples[s].Data()) << ": " << h_tmp->Integral() << endl;

    // add the histogram to the map
    aMap[datacard_samples[s]] = h_tmp;

    // add the histogram to the asymov dataset
    h_data->Add( h_tmp );

    // add the histogram yield to the asymov dataset
    observation +=  h_tmp->Integral();
  }
  cout << "************************" << endl;

  if(errFlag>0){
    cout << "Missing histogram for a process. Return." << endl;
    return;
  }

  // approximate to integer precision
  observation = int(observation);
  h_data->Scale(observation/h_data->Integral());

  // save data
  dir->cd();
  h_data->Write("data_obs", TObject::kOverwrite);

  // return if you don't want to create the datadacards
  if( CREATEDATACARDS==0 ){
    fout->Close();
    return;
  }

 
  // the text file for the datacard
  ofstream out(Form("datacards/%s_%s.txt",(fname+"_"+category).Data(), (name+version+extraname).c_str()));
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


  //////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////
  // add processes only if yield is non null

  string line("bin                         ");
  if( aMap["TTH125"]->Integral()>0     ) line += string(Form("%s    ", (fname+"_"+category).Data()));
  if( aMap["TTJetsHFbb"]->Integral()>0 ) line += string(Form("%s    ", (fname+"_"+category).Data()));
  if( aMap["TTJetsHFb"]->Integral()>0  ) line += string(Form("%s    ", (fname+"_"+category).Data()));
  if( aMap["TTJetsLF"]->Integral()>0   ) line += string(Form("%s    ", (fname+"_"+category).Data()));
  if( aMap["TTV"]->Integral()>0        ) line += string(Form("%s    ", (fname+"_"+category).Data()));
  if( aMap["SingleT"]->Integral()>0    ) line += string(Form("%s    ", (fname+"_"+category).Data()));
  out<<line;
  out<<endl; 
  line = "process                     ";
  if( aMap["TTH125"]->Integral()>0 )     line += "TTH125      ";
  if( aMap["TTJetsHFbb"]->Integral()>0 ) line += "TTJetsHFbb    ";
  if( aMap["TTJetsHFb"]->Integral()>0 )  line += "TTJetsHFb    ";
  if( aMap["TTJetsLF"]->Integral()>0 )   line += "TTJetsLF    ";
  if( aMap["TTV"]->Integral()>0 )        line += "TTV      ";
  if( aMap["SingleT"]->Integral()>0 )    line += "SingleT     ";
  out<<line;
  out<<endl;
  line = "process                       ";
  if( aMap["TTH125"]->Integral()>0 )     line += "0          ";
  if( aMap["TTJetsHFbb"]->Integral()>0 ) line += "1          ";
  if( aMap["TTJetsHFb"]->Integral()>0 )  line += "2          ";
  if( aMap["TTJetsLF"]->Integral()>0 )   line += "3          ";
  if( aMap["TTV"]->Integral()>0 )        line += "4          ";
  if( aMap["SingleT"]->Integral()>0 )    line += "5          ";
  out<<line;
  out<<endl;
  line = "rate                        ";
  if( aMap["TTH125"]->Integral()>0 )     line += string(Form("%.4f     ", aMap["TTH125"]->Integral() ));
  if( aMap["TTJetsHFbb"]->Integral()>0 ) line += string(Form("%.4f     ", aMap["TTJetsHFbb"]->Integral() ));
  if( aMap["TTJetsHFb"]->Integral()>0 )  line += string(Form("%.4f     ", aMap["TTJetsHFb"]->Integral() ));
  if( aMap["TTJetsLF"]->Integral()>0 )   line += string(Form("%.4f     ", aMap["TTJetsLF"]->Integral()));
  if( aMap["TTV"]->Integral()>0 )        line += string(Form("%.4f     ", aMap["TTV"]->Integral()));
  if( aMap["SingleT"]->Integral()>0 )    line += string(Form("%.4f     ", aMap["SingleT"]->Integral()));
  out<<line;
  out<<endl;
  out<<endl;
  out<< "-----------------------------------------------------------------" << endl;

  //////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////
  // systematics

  // LUMINOSITY
  line = "lumi                  lnN  ";
  if( aMap["TTH125"]->Integral()>0 )     line += "1.022      ";
  if( aMap["TTJetsHFbb"]->Integral()>0 ) line += "1.022      ";
  if( aMap["TTJetsHFb"]->Integral()>0 )  line += "1.022      ";
  if( aMap["TTJetsLF"]->Integral()>0 )   line += "1.022      ";
  if( aMap["TTV"]->Integral()>0 )        line += "1.022      ";
  if( aMap["SingleT"]->Integral()>0 )    line += "1.022      ";
  out<<line;
  out<<endl;

  // BTAGGING
  line = "csv                   shape  ";                              
  if( aMap["TTH125"]->Integral()>0 )     line += "1.0        ";       
  if( aMap["TTJetsHFbb"]->Integral()>0 ) line += "1.0        ";       
  if( aMap["TTJetsHFb"]->Integral()>0 )  line += "1.0        ";       
  if( aMap["TTJetsLF"]->Integral()>0 )   line += "1.0        ";       
  if( aMap["TTV"]->Integral()>0 )        line += " -         ";
  if( aMap["SingleT"]->Integral()>0 )    line += " -         ";
  out<<line;
  out<<endl;

  // JET ENERGY SCALE
  line = "JEC                   shape  ";                             
  if( aMap["TTH125"]->Integral()>0 )     line += "1.0        ";
  if( aMap["TTJetsHFbb"]->Integral()>0 ) line += "1.0        ";
  if( aMap["TTJetsHFb"]->Integral()>0 )  line += "1.0        ";
  if( aMap["TTJetsLF"]->Integral()>0 )   line += "1.0        ";
  if( aMap["TTV"]->Integral()>0 )        line += " -         ";
  if( aMap["SingleT"]->Integral()>0 )    line += " -         ";
  out<<line;
  out<<endl;

  // TTBB CROSS-SECTION
  line = "Norm_TTbb             lnN    ";
  if( aMap["TTH125"]->Integral()>0 )     line += " -         ";
  if( aMap["TTJetsHFbb"]->Integral()>0 ) line += "1.50       ";
  if( aMap["TTJetsHFb"]->Integral()>0 )  line += " -       ";
  if( aMap["TTJetsLF"]->Integral()>0 )   line += " -         ";
  if( aMap["TTV"]->Integral()>0 )        line += " -         ";
  if( aMap["SingleT"]->Integral()>0 )    line += " -         ";
  out<<line;
  out<<endl;

  // TTB CROSS-SECTION
  line = "Norm_TTb              lnN    ";
  if( aMap["TTH125"]->Integral()>0 )     line += " -         ";
  if( aMap["TTJetsHFbb"]->Integral()>0 ) line += " -       ";
  if( aMap["TTJetsHFb"]->Integral()>0 )  line += "1.50       ";
  if( aMap["TTJetsLF"]->Integral()>0 )   line += " -         ";
  if( aMap["TTV"]->Integral()>0 )        line += " -         ";
  if( aMap["SingleT"]->Integral()>0 )    line += " -         ";
  out<<line;
  out<<endl;

  // TTV CROSS-SECTION
  line = "Norm_TTV              lnN    ";
  if( aMap["TTH125"]->Integral()>0 )     line += " -         ";
  if( aMap["TTJetsHFbb"]->Integral()>0 ) line += " -         ";
  if( aMap["TTJetsHFb"]->Integral()>0 )  line += " -         ";
  if( aMap["TTJetsLF"]->Integral()>0 )   line += " -         ";
  if( aMap["TTV"]->Integral()>0 )        line += "1.20       ";
  if( aMap["SingleT"]->Integral()>0 )    line += " -         ";
  out<<line;
  out<<endl;

  // SINGLE-TOP CROSS-SECTION
  line = "Norm_SingleT          lnN    ";
  if( aMap["TTH125"]->Integral()>0 )     line += " -         ";
  if( aMap["TTJetsHFbb"]->Integral()>0 ) line += " -         ";
  if( aMap["TTJetsHFb"]->Integral()>0 )  line += " -         ";
  if( aMap["TTJetsLF"]->Integral()>0 )   line += " -         ";
  if( aMap["TTV"]->Integral()>0 )        line += " -         ";
  if( aMap["SingleT"]->Integral()>0 )    line += "1.20       ";
  out<<line;
  out<<endl;

  // QCD SCALE UNC. ON TTH
  line = "QCDscale_TTH          lnN    ";
  if( aMap["TTH125"]->Integral()>0 )     line += "1.12       ";
  if( aMap["TTJetsHFbb"]->Integral()>0 ) line += " -         ";
  if( aMap["TTJetsHFb"]->Integral()>0 )  line += " -         ";
  if( aMap["TTJetsLF"]->Integral()>0 )   line += " -         ";
  if( aMap["TTV"]->Integral()>0 )        line += " -         ";
  if( aMap["SingleT"]->Integral()>0 )    line += " -         ";
  out<<line;
  out<<endl;

  // QCD SCALE UNC. ON TT+HF
  line = string(Form("QCDscale%s_TTJetsHF     lnN    ", null.c_str()));
  if( aMap["TTH125"]->Integral()>0 )     line += " -         ";
  if( aMap["TTJetsHFbb"]->Integral()>0 ) line += "1.35       ";
  if( aMap["TTJetsHFb"]->Integral()>0 )  line += "1.35       ";
  if( aMap["TTJetsLF"]->Integral()>0 )   line += " -         ";
  if( aMap["TTV"]->Integral()>0 )        line += " -         ";
  if( aMap["SingleT"]->Integral()>0 )    line += " -         ";
  out<<line;
  out<<endl;

 // QCD SCALE UNC. ON TT+LF
  line = string(Form("QCDscale%s_TTJetsLF     lnN    ", null.c_str()));
  if( aMap["TTH125"]->Integral()>0 )     line += " -         ";
  if( aMap["TTJetsHFbb"]->Integral()>0 ) line += " -         ";
  if( aMap["TTJetsHFb"]->Integral()>0 )  line += " -         ";
  if( aMap["TTJetsLF"]->Integral()>0 )   line += "1.35       ";
  if( aMap["TTV"]->Integral()>0 )        line += " -         ";
  if( aMap["SingleT"]->Integral()>0 )    line += " -         ";
  out<<line;
  out<<endl;

  // PDFs
  line = "pdf_gg                lnN    ";
  if( aMap["TTH125"]->Integral()>0 )     line += "1.03       ";
  if( aMap["TTJetsHFbb"]->Integral()>0 ) line += "1.03       ";
  if( aMap["TTJetsHFb"]->Integral()>0 )  line += "1.03       ";
  if( aMap["TTJetsLF"]->Integral()>0 )   line += "1.03       ";
  if( aMap["TTV"]->Integral()>0 )        line += "1.03       ";
  if( aMap["SingleT"]->Integral()>0 )    line += "1.03       ";
  out<<line;
  out<<endl;

  //////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////
  // bin-by-bin

  for(int bin = 1; bin<=aMap["TTH125"]->GetNbinsX(); bin++ ){
    if(aMap["TTH125"]->GetBinContent(bin)>0){
      line = string(Form("TTH125%sbin%d        shape  ", category.Data(), bin));
      if( aMap["TTH125"]->Integral()>0 )     line += "1.0        ";
      if( aMap["TTJetsHFbb"]->Integral()>0 ) line += " -         ";
      if( aMap["TTJetsHFb"]->Integral()>0 )  line += " -         ";
      if( aMap["TTJetsLF"]->Integral()>0 )   line += " -         ";
      if( aMap["TTV"]->Integral()>0 )        line += " -         ";
      if( aMap["SingleT"]->Integral()>0 )    line += " -         ";
      out<<line;
      out<<endl;
    }
  }  
  for(int bin = 1; bin<=aMap["TTJetsHFbb"]->GetNbinsX(); bin++ ){
    if(aMap["TTJetsHFbb"]->GetBinContent(bin)>0){
      line = string(Form("TTJetsHFbb%sbin%d      shape  ", category.Data(), bin));
      if( aMap["TTH125"]->Integral()>0 )     line += " -         ";
      if( aMap["TTJetsHFbb"]->Integral()>0 ) line += "1.0        ";
      if( aMap["TTJetsHFb"]->Integral()>0 )  line += " -        ";
      if( aMap["TTJetsLF"]->Integral()>0 )   line += " -         ";
      if( aMap["TTV"]->Integral()>0 )        line += " -         ";
      if( aMap["SingleT"]->Integral()>0 )    line += " -         ";
      out<<line;
      out<<endl;
    }
  }  
  for(int bin = 1; bin<=aMap["TTJetsHFb"]->GetNbinsX(); bin++ ){
    if(aMap["TTJetsHFb"]->GetBinContent(bin)>0){
      line = string(Form("TTJetsHFb%sbin%d      shape  ", category.Data(), bin));
      if( aMap["TTH125"]->Integral()>0 )     line += " -         ";
      if( aMap["TTJetsHFbb"]->Integral()>0 ) line += " -        ";
      if( aMap["TTJetsHFb"]->Integral()>0 )  line += "1.0        ";
      if( aMap["TTJetsLF"]->Integral()>0 )   line += " -         ";
      if( aMap["TTV"]->Integral()>0 )        line += " -         ";
      if( aMap["SingleT"]->Integral()>0 )    line += " -         ";
      out<<line;
      out<<endl;
    }
  }  
  for(int bin = 1; bin<=aMap["TTJetsLF"]->GetNbinsX(); bin++ ){
    if(aMap["TTJetsLF"]->GetBinContent(bin)>0){
      line = string(Form("TTJetsLF%sbin%d      shape  ", category.Data(), bin));
      if( aMap["TTH125"]->Integral()>0 )     line += " -         ";
      if( aMap["TTJetsHFbb"]->Integral()>0 ) line += " -         ";
      if( aMap["TTJetsHFb"]->Integral()>0 )  line += " -         ";
      if( aMap["TTJetsLF"]->Integral()>0 )   line += "1.0        ";
      if( aMap["TTV"]->Integral()>0 )        line += " -         ";
      if( aMap["SingleT"]->Integral()>0 )    line += " -         ";
      out<<line;
      out<<endl;
    }
  }  
  for(int bin = 1; bin<=aMap["TTV"]->GetNbinsX(); bin++ ){
    if(aMap["TTV"]->GetBinContent(bin)>0){
      line = string(Form("TTV%sbin%d           shape  ", category.Data(), bin));
      if( aMap["TTH125"]->Integral()>0 )      line += " -         ";
      if( aMap["TTJetsHFbb"]->Integral()>0 )  line += " -         ";
      if( aMap["TTJetsHFb"]->Integral()>0 )   line += " -         ";
      if( aMap["TTJetsLF"]->Integral()>0 )    line += " -         ";
      if( aMap["TTV"]->Integral()>0 )         line += "1.0        ";
      if( aMap["SingleT"]->Integral()>0 )     line += " -         ";
      out<<line;
      out<<endl;
    }
  } 
  for(int bin = 1; bin<=aMap["SingleT"]->GetNbinsX(); bin++ ){
    if(aMap["SingleT"]->GetBinContent(bin)>0){
      line = string(Form("SingleT%sbin%d       shape  ", category.Data(), bin));
      if( aMap["TTH125"]->Integral()>0 )     line += " -         ";
      if( aMap["TTJetsHFbb"]->Integral()>0 ) line += " -         ";
      if( aMap["TTJetsHFb"]->Integral()>0 )  line += " -         ";
      if( aMap["TTJetsLF"]->Integral()>0 )   line += " -         ";
      if( aMap["TTV"]->Integral()>0 )        line += " -         ";
      if( aMap["SingleT"]->Integral()>0 )    line += "1.0        ";
      out<<line;
      out<<endl;
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
  cout << " > Category = " << category << endl;
  cout << " > Cut = (" << cut  << ")" << endl;
  cout << " > Histo = h(" << h->GetNbinsX() << "," << h->GetXaxis()->GetXmin() << "," <<  h->GetXaxis()->GetXmax() << ")" << endl;
  cout << " > Lumi = " << lumiScale*12.1 << " fb-1" << endl << endl;

  // close the output file
  fout->Close();

  return;
}


void produceAllNew(string name = "New", string version = "_v2_reg",  string extraname = "", 
		   float LumiScale = 19.5/12.1, int doMEM = 1 ){

  vector<float> binvec;
  binvec.push_back( 55.);
  binvec.push_back( 75.);
  binvec.push_back( 95.);
  binvec.push_back(115.);
  binvec.push_back(135.);
  binvec.push_back(155.);
  binvec.push_back(175.);
  binvec.push_back(195.);

  TF1* xsec = new TF1("xsec",Form("x^(-%f)", 1. ), 20, 500);
  xsec = 0;

  if(doMEM){
    produceNew( name, version, extraname,  "SL", Form("type==0 && btag_LR>=%f",                  0.975), "cat1", doMEM, 2.2, 0.00, 0.2 , LumiScale   , 6,  1);
    produceNew( name, version, extraname,  "SL", Form("type==1 && btag_LR>=%f",                  0.975), "cat2", doMEM, 1.7, 0.15, 0.5 , LumiScale   , 6,  1);
    produceNew( name, version, extraname,  "SL", Form("type==2 && flag_type2>0  && btag_LR>=%f", 0.986), "cat3", doMEM, 2.2, 0.00, 0.5 , LumiScale   , 6,  1);
    produceNew( name, version, extraname,  "SL", Form("type==2 && flag_type2<=0 && btag_LR>=%f", 0.990), "cat4", doMEM, 2.0, 0.30, 0.5 , LumiScale   , 6,  1);
    produceNew( name, version, extraname,  "SL", Form("type==3 && btag_LR>=%f",                  0.975), "cat5", doMEM, 5.5, 0.00, 0.5 , LumiScale   , 7,  1);
    produceNew( name, version, extraname,  "DL", Form("type==6 && btag_LR>=%f",                  0.953), "cat6", doMEM, 2.2, 0.00, 0.1 , LumiScale*2 , 5,  1);
  }
  else{
    produceNew( name, version, extraname,  "SL", "type==0",                                  "cat1", doMEM, 0.0, 0.0, 0.0 , LumiScale   , binvec.size()-1,  0, &binvec);
    produceNew( name, version, extraname,  "SL", "type==2",                                  "cat2", doMEM, 0.0, 0.0, 0.0 , LumiScale   , binvec.size()-1,  0, &binvec);
    produceNew( name, version, extraname,  "SL", "type==3 && flag_type3>0 && p_125_all_s>0", "cat3", doMEM, 0.0, 0.0, 0.0 , LumiScale   , binvec.size()-1,  0, &binvec);
    produceNew( name, version, extraname,  "DL", "type==6",                                  "cat4", doMEM, 0.0, 0.0, 0.0 , LumiScale*2 , binvec.size()-1,  0, &binvec);
  }

}


void optimize_mem_category_btag(string name = "New", string version = "_v2_reg",  string extraname = "", 
				float LumiScale = 19.5/12.1){
  
  int doMEM = 1;

  vector<float> btag_LRs;
  btag_LRs.push_back(0.946);
  btag_LRs.push_back(0.948);
  btag_LRs.push_back(0.950);
  btag_LRs.push_back(0.952);
  btag_LRs.push_back(0.954);
  btag_LRs.push_back(0.956);
  btag_LRs.push_back(0.958);

  for(int it = 0; it<btag_LRs.size() ; it++){    
    //if(it<6) continue;
    //produceNew( name, version, extraname+string(Form("_%d", it)),  "SL", Form("type==0 && btag_LR>=%f", btag_LRs[it]),                                  "cat1", doMEM, 1.0, 0.0, 0.2 , LumiScale   , 6,  1);
    //produceNew( name, version, extraname+string(Form("_%d", it)),  "SL", Form("type==1 && btag_LR>=%f", btag_LRs[it]),                                  "cat2", doMEM, 1.7, 0.0, 0.5 , LumiScale   , 6,  1);
    //produceNew( name, version, extraname+string(Form("_%d", it)),  "SL", Form("type==2 && flag_type2>0 && btag_LR>=%f", btag_LRs[it]),                  "cat3", doMEM, 2.2, 0.0, 0.5 , LumiScale   , 6,  1);
    //produceNew( name, version, extraname+string(Form("_%d", it)),  "SL", Form("type==2 && flag_type2<=0 && btag_LR>=%f", btag_LRs[it]),                 "cat4", doMEM, 2.0, 0.0, 0.5 , LumiScale   , 6,  1);
    //produceNew( name, version, extraname+string(Form("_%d", it)),  "SL", Form("type==3 && flag_type3>0 && p_125_all_s>0 && btag_LR>=%f", btag_LRs[it]), "cat5", doMEM, 5.5, 0.0, 0.5 , LumiScale   , 7,  1);
    produceNew( name, "_v3_reg" , extraname+string(Form("_%d", it)),  "DL", Form("type==6 && btag_LR>=%f", btag_LRs[it]),                                 "cat6", doMEM, 2.2, 0.0, 0.1 , LumiScale*2 , 5,  1);
  }

  // cat1 0.975
  // cat2 0.975
  // cat3 0.980
  // cat4 0.990
  // cat5 0.975
  // cat6 0.950

}



void optimize_mem_category_rel(string name = "New", string version = "_v2_reg",  string extraname = "", 
			       float LumiScale = 19.5/12.1){
  
  int doMEM = 1;

  vector<float> rels;
  rels.push_back(1.0);
  rels.push_back(1.40);
  rels.push_back(1.80);
  rels.push_back(2.20);
  rels.push_back(2.60);
  rels.push_back(3.00);
  rels.push_back(3.40);

  for(int it = 0; it<rels.size() ; it++){    
    produceNew( name, version, extraname+string(Form("_%d", it)),  "SL", Form("type==0 && btag_LR>=%f", 0.975 ), "cat1", doMEM, rels[it] , 0.0, 0.2 , LumiScale   , 6,  1);
  }

  rels.clear();
  rels.push_back(0.95);
  rels.push_back(1.20);
  rels.push_back(1.45);
  rels.push_back(1.70);
  rels.push_back(1.95);
  rels.push_back(2.20);
  rels.push_back(2.45);

  for(int it = 0; it<rels.size() ; it++){    
    //produceNew( name, version, extraname+string(Form("_%d", it)),  "SL", Form("type==1 && btag_LR>=%f", 0.975 ), "cat2", doMEM, rels[it] , 0.0, 0.5 , LumiScale   , 6,  1);
  }


  rels.clear();
  rels.push_back(1.45);
  rels.push_back(1.70);
  rels.push_back(1.95);
  rels.push_back(2.20);
  rels.push_back(2.45);
  rels.push_back(2.70);
  rels.push_back(2.95);

  for(int it = 0; it<rels.size() ; it++){    
    //produceNew( name, version, extraname+string(Form("_%d", it)),  "SL", Form("type==2 && flag_type2>0 && btag_LR>=%f", 0.980 ), "cat3", doMEM, rels[it] , 0.0, 0.5 , LumiScale   , 6,  1);
  }


  rels.clear();
  rels.push_back(0.80);
  rels.push_back(1.20);
  rels.push_back(1.60);
  rels.push_back(2.00);
  rels.push_back(2.40);
  rels.push_back(2.80);
  rels.push_back(3.20);

  for(int it = 0; it<rels.size() ; it++){    
    //produceNew( name, version, extraname+string(Form("_%d", it)),  "SL", Form("type==2 && flag_type2<=0 && btag_LR>=%f", 0.990 ), "cat4", doMEM, rels[it] , 0.0, 0.5 , LumiScale   , 6,  1);
  }


  rels.clear();
  rels.push_back(4.00);
  rels.push_back(4.50);
  rels.push_back(5.00);
  rels.push_back(5.50);
  rels.push_back(6.00);
  rels.push_back(6.50);
  rels.push_back(7.00);

  for(int it = 0; it<rels.size() ; it++){    
    //produceNew( name, version, extraname+string(Form("_%d", it)),  "SL", Form("type==3 && btag_LR>=%f", 0.975 ), "cat5", doMEM, rels[it] , 0.0, 0.5 , LumiScale   , 7,  1);
  }

  rels.clear();
  rels.push_back(0.30);
  rels.push_back(0.70);
  rels.push_back(1.10);
  rels.push_back(1.50);
  rels.push_back(1.90);
  rels.push_back(2.30);
  rels.push_back(2.70);

  for(int it = 0; it<rels.size() ; it++){    
    //produceNew( name, version, extraname+string(Form("_%d", it)),  "DL", Form("type==6 && btag_LR>=%f", 0.950 ), "cat6", doMEM, rels[it] , 0.0, 0.1 , LumiScale*2   , 5,  1);
  }


  // cat1 0.975
  // cat2 0.975
  // cat3 0.980
  // cat4 0.990
  // cat5 0.975
  // cat6 0.950

}



void optimize_mem_category_rel_bbjj(string name = "New", string version = "_v2_reg",  string extraname = "", 
				    float LumiScale = 19.5/12.1){
  
  int doMEM = 1;

  vector<float> rels;
  rels.push_back(-0.45);
  rels.push_back(-0.30);
  rels.push_back(-0.15);
  rels.push_back(0.0);
  rels.push_back(0.15);
  rels.push_back(0.30);
  rels.push_back(0.45);

  for(int it = 0; it<rels.size() ; it++){    
    //produceNew( name, version, extraname+string(Form("_%d", it)),  "SL", Form("type==0 && btag_LR>=%f", 0.975 ), "cat1", doMEM, 2.2 , rels[it], 0.2 , LumiScale   , 6,  1);
  }

  for(int it = 0; it<rels.size() ; it++){    
    //produceNew( name, version, extraname+string(Form("_%d", it)),  "SL", Form("type==1 && btag_LR>=%f", 0.975 ), "cat2", doMEM, 1.7, rels[it] , 0.5 , LumiScale   , 6,  1);
  }

  for(int it = 0; it<rels.size() ; it++){    
    //produceNew( name, version, extraname+string(Form("_%d", it)),  "SL", Form("type==2 && flag_type2>0 && btag_LR>=%f", 0.986 ), "cat3", doMEM, 2.2, rels[it] , 0.5 , LumiScale   , 6,  1);
  }

  for(int it = 0; it<rels.size() ; it++){    
    produceNew( name, version, extraname+string(Form("_%d", it)),  "SL", Form("type==2 && flag_type2<=0 && btag_LR>=%f", 0.990 ), "cat4", doMEM, 2.0, rels[it], 0.5 , LumiScale   , 6,  1);
  }

  for(int it = 0; it<rels.size() ; it++){    
    produceNew( name, version, extraname+string(Form("_%d", it)),  "SL", Form("type==3 && btag_LR>=%f", 0.975 ), "cat5", doMEM, 5.5 , rels[it], 0.5 , LumiScale   , 7,  1);
  }

  for(int it = 0; it<rels.size() ; it++){    
    produceNew( name, version, extraname+string(Form("_%d", it)),  "DL", Form("type==6 && btag_LR>=%f", 0.950 ), "cat6", doMEM, 2.2 , rels[it], 0.1 , LumiScale*2   , 5,  1);
  }

}
