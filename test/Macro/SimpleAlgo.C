#include <cstdlib>
#include <iostream> 
#include <map>
#include <string>
#include <vector>

#include "TChain.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "THStack.h"
#include "TCut.h"
#include "Math/GenVector/LorentzVector.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#if not defined(__CINT__) || defined(__MAKECINT__)

#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#endif

using namespace std;

typedef  ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LV;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LV2;


struct sorterByPt {
  bool operator() (float i,float j) const { return (i>j);}
};

struct sorterByLess {
  bool operator() (float i,float j) const { return (i<j);}
};


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


typedef struct
{
  int index;
  float pt;
  float eta;
  float phi;
  float mass;
  float csv;
  float topB;
  float topW;
  float atopB;
  float atopW;
  float higgsB;
  float flavor;
  float unc;
  void reset(){
    index = -99; pt = -99; eta = -99; phi = -99; mass = -99; csv = -99; topB = -99; topW = -99;  atopW = -99;  atopB = -99;
    higgsB = -99; flavor = -99; unc = -99;
  }
  void set(int index_,  float pt_, float eta_, float phi_, float mass_,float csv_, float topB_, float topW_, float atopB_,
	   float atopW_, float higgsB_, float flavor_, float unc_){
    index = index_; 
    pt = pt_; 
    eta = eta_; 
    phi = phi_; 
    mass = mass_; 
    csv = csv_; 
    topB = topB_; 
    topW = topW_;  
    atopW = atopW_;  
    atopB = atopB_;
    higgsB = higgsB_; 
    flavor = flavor_; 
    unc = unc_;
  }
} JetByPt;


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
   float wmass;
   float wpt;
   float weta;
   float wphi;
   float wid;
   float topmass;
   float toppt;
   float topeta;
   float topphi;
   float topid;
   float bbarmass;
   float bbarpt;
   float bbareta;
   float bbarphi;
 } recoTopInfo;



float deltaR(LV reco, LV gen){
  return TMath::Sqrt( (reco.Eta()-gen.Eta())*(reco.Eta()-gen.Eta()) +  TMath::ACos( TMath::Cos(reco.Phi()-gen.Phi()) )* TMath::ACos( TMath::Cos(reco.Phi()-gen.Phi()) ) )  ;
}


void findGenMatch(int& genMatch, LV genJet, LV topBLV, LV topW1LV, LV topW2LV, LV atopBLV, LV atopW1LV, LV atopW2LV, LV genBLV, LV genBbarLV){

  if(genJet.Pt()<=0.001){
    genMatch = -99;
    return;
  }

  int numMatches = 0;

  if( topBLV.Pt()>0 && deltaR(genJet, topBLV)       < 0.3 ){
    genMatch = 1;
    numMatches++;
  }
  if( topW1LV.Pt()>0 && deltaR(genJet, topW1LV)     < 0.3 ){
    genMatch = 2;
    numMatches++;
  }
  if( topW2LV.Pt()>0 && deltaR(genJet, topW2LV)     < 0.3 ){
    genMatch = 2;
    numMatches++;
  }
  if( atopBLV.Pt()>0 && deltaR(genJet, atopBLV)     < 0.3 ){
    genMatch = -1;
    numMatches++;
  }
  if( atopW1LV.Pt()>0 && deltaR(genJet, atopW1LV)   < 0.3 ){
    genMatch = -2;
    numMatches++;
  }
  if( atopW2LV.Pt()>0 && deltaR(genJet, atopW2LV)   < 0.3 ){
    genMatch = -2;
    numMatches++;
  }
  if( genBLV.Pt()>0 && deltaR(genJet, genBLV)       < 0.3 ){
    genMatch = 0;
    numMatches++;
  }
  if( genBbarLV.Pt()>0 && deltaR(genJet, genBbarLV) < 0.3 ){
    genMatch = 0;
    numMatches++;
  }

  if(numMatches==0)
    genMatch = 3; // NO MATCH 
  if(numMatches>1)
    genMatch = 4; // AMBIGUOUS MATCH



  return;

}

void findGenMatch2(int& genMatch, LV genJet, LV topBLV, LV topW1LV, LV topW2LV, LV atopBLV, LV atopW1LV, LV atopW2LV, LV genBLV, LV genBbarLV){

  if(genJet.Pt()<=0.001){
    genMatch = -99;
    return;
  }

  int numMatches = 0;

  if( topBLV.Pt()>0 && deltaR(genJet, topBLV)       < 0.3 ){
    genMatch = 1;
    numMatches++;
  }
  if( topW1LV.Pt()>0 && deltaR(genJet, topW1LV)     < 0.3 ){
    genMatch = 2;
    numMatches++;
  }
  if( topW2LV.Pt()>0 && deltaR(genJet, topW2LV)     < 0.3 ){
    genMatch = 2;
    numMatches++;
  }
  if( atopBLV.Pt()>0 && deltaR(genJet, atopBLV)     < 0.3 ){
    genMatch = -1;
    numMatches++;
  }
  if( atopW1LV.Pt()>0 && deltaR(genJet, atopW1LV)   < 0.3 ){
    genMatch = -2;
    numMatches++;
  }
  if( atopW2LV.Pt()>0 && deltaR(genJet, atopW2LV)   < 0.3 ){
    genMatch = -2;
    numMatches++;
  }
  if( genBLV.Pt()>0 && deltaR(genJet, genBLV)       < 0.3 ){
    genMatch = 5;
    numMatches++;
  }
  if( genBbarLV.Pt()>0 && deltaR(genJet, genBbarLV) < 0.3 ){
    genMatch = -5;
    numMatches++;
  }

  if(numMatches==0)
    genMatch = 3; // NO MATCH 
  if(numMatches>1)
    genMatch = 4; // AMBIGUOUS MATCH



  return;

}


bool bookKeeper(unsigned int index1, unsigned int index2, 
		unsigned int index3, unsigned int index4, 
		vector<unsigned int>& visited){

  unsigned int combination1 = 1000*index1 + 100*index2 + 10*index3 + 1*index4;
  unsigned int combination2 = 1000*index2 + 100*index1 + 10*index3 + 1*index4;
  bool isThere = false;

  for(unsigned int k = 0; k < visited.size() ; k++ )
    if( combination1 == visited[k] || combination2 == visited[k] ) isThere = true;

  if(!isThere){
    visited.push_back( combination1 );
    visited.push_back( combination2 );
    return true;
  }
  else{
    return false;
  }

}


bool bookKeeper2(unsigned int index1, unsigned int index2, 
		 unsigned int index3, unsigned int index4, 
		 vector<unsigned int>& visited){

  unsigned int combination1 = 1000*index1 + 100*index2 + 10*index3 + 1*index4;
  unsigned int combination2 = 1000*index2 + 100*index1 + 10*index3 + 1*index4;
  unsigned int combination3 = 1000*index1 + 100*index2 + 10*index4 + 1*index3;
  unsigned int combination4 = 1000*index2 + 100*index1 + 10*index4 + 1*index3;

  bool isThere = false;

  for(unsigned int k = 0; k < visited.size() ; k++ )
    if( combination1 == visited[k] || combination2 == visited[k] || combination3 == visited[k] || combination4 == visited[k]) isThere = true;

  if(!isThere){
    visited.push_back( combination1 );
    visited.push_back( combination2 );
    return true;
  }
  else{
    return false;
  }

}



void makeTree(bool sgn = true, TString weightFile = "TMVAClassification_SelM_BDT.weights.xml"){


  TMVA::Tools::Instance();
  TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );   


  TRandom3* ran = new TRandom3();

  float LOWVALUE=0.0001;

  TString fName = "";
  if(sgn) 
    fName = "SgnTree.root";
  else
    fName = "BkgTree.root";

  TFile *fOutput           = new TFile(fName,"RECREATE");
  TTree* outTree           = new TTree("outTree","");
  float numJets20, numJets30, numJets60, numJets100, numJetsBtagL, numJetsBtagM, numJetsBtagT;
  float matchesL, matchesM, matchesT;
  float weightL, weightM, weightT;
  int type;
  float nSvs;

  float m12L,m13L,m14L,m23L,m24L,m34L;
  float m12M,m13M,m14M,m23M,m24M,m34M;
  float m12T,m13T,m14T,m23T,m24T,m34T;

  float mRanL, mPtL, mSumPtL, mPtoML, mClosL, mClos2L;
  float mRanM, mPtM, mSumPtM, mPtoMM, mClosM, mClos2M;
  float mRanT, mPtT, mSumPtT, mPtoMT, mClosT, mClos2T;

  int nClosL, nClosM, nClosT;

  float BDT;

  TBranch *typeBR = outTree->Branch("type",&type,"type/I"); 
  TBranch *nSvsBR = outTree->Branch("nSvs",&nSvs,"nSvs/F"); 
  TBranch *numJets20BR  = outTree->Branch("numJets20", &numJets20,"numJets20/F");
  TBranch *numJets30BR  = outTree->Branch("numJets30", &numJets30,"numJets30/F"); 
  TBranch *numJets60BR  = outTree->Branch("numJets60", &numJets60,"numJets60/F"); 
  TBranch *numJets100BR = outTree->Branch("numJets100",&numJets100,"numJets100/F"); 

  TBranch *numJetsBtagLBR = outTree->Branch("numJetsBtagL",&numJetsBtagL,"numJetsBtagL/F"); 
  TBranch *numJetsBtagMBR = outTree->Branch("numJetsBtagM",&numJetsBtagM,"numJetsBtagM/F"); 
  TBranch *numJetsBtagTBR = outTree->Branch("numJetsBtagT",&numJetsBtagT,"numJetsBtagT/F"); 

  TBranch *weightLBR = outTree->Branch("weightL",&weightL,"weightL/F"); 
  TBranch *weightMBR = outTree->Branch("weightM",&weightM,"weightM/F"); 
  TBranch *weightTBR = outTree->Branch("weightT",&weightT,"weightT/F"); 

  TBranch *matchesLBR = outTree->Branch("matchesL",&matchesL,"matchesL/F"); 
  TBranch *matchesMBR = outTree->Branch("matchesM",&matchesM,"matchesM/F"); 
  TBranch *matchesTBR = outTree->Branch("matchesT",&matchesT,"matchesT/F");

  TBranch *m12LBR = outTree->Branch("m12L",&m12L,"m12L/F"); 
  TBranch *m13LBR = outTree->Branch("m13L",&m13L,"m13L/F"); 
  TBranch *m14LBR = outTree->Branch("m14L",&m14L,"m14L/F"); 
  TBranch *m23LBR = outTree->Branch("m23L",&m23L,"m23L/F"); 
  TBranch *m24LBR = outTree->Branch("m24L",&m24L,"m24L/F"); 
  TBranch *m34LBR = outTree->Branch("m34L",&m34L,"m34L/F"); 

  TBranch *m12MBR = outTree->Branch("m12M",&m12M,"m12M/F"); 
  TBranch *m13MBR = outTree->Branch("m13M",&m13M,"m13M/F"); 
  TBranch *m14MBR = outTree->Branch("m14M",&m14M,"m14M/F"); 
  TBranch *m23MBR = outTree->Branch("m23M",&m23M,"m23M/F"); 
  TBranch *m24MBR = outTree->Branch("m24M",&m24M,"m24M/F"); 
  TBranch *m34MBR = outTree->Branch("m34M",&m34M,"m34M/F"); 

  TBranch *m12TBR = outTree->Branch("m12T",&m12T,"m12T/F"); 
  TBranch *m13TBR = outTree->Branch("m13T",&m13T,"m13T/F"); 
  TBranch *m14TBR = outTree->Branch("m14T",&m14T,"m14T/F"); 
  TBranch *m23TBR = outTree->Branch("m23T",&m23T,"m23T/F"); 
  TBranch *m24TBR = outTree->Branch("m24T",&m24T,"m24T/F"); 
  TBranch *m34TBR = outTree->Branch("m34T",&m34T,"m34T/F"); 

  TBranch *nClosLBR = outTree->Branch("nClosL",&nClosL,"nClosL/I"); 
  TBranch *nClosMBR = outTree->Branch("nClosM",&nClosM,"nClosM/I"); 
  TBranch *nClosTBR = outTree->Branch("nClosT",&nClosT,"nClosT/I"); 
  


  TBranch *mRanLBR    = outTree->Branch("mRanL",&mRanL,"mRanL/F"); 
  TBranch *mPtLBR     = outTree->Branch("mPtL" ,&mPtL,  "mPtL/F");  
  TBranch *mSumPtLBR  = outTree->Branch("mSumPtL" ,&mSumPtL,  "mSumPtL/F"); 
  TBranch *mPtoMLBR   = outTree->Branch("mPtoML" ,&mPtoML,  "mPtoML/F");  
  TBranch *mClosLBR   = outTree->Branch("mClosL" ,&mClosL,  "mClosL/F");  
  TBranch *mClos2LBR   = outTree->Branch("mClos2L" ,&mClos2L,  "mClos2L/F");  


  TBranch *mRanMBR    = outTree->Branch("mRanM",&mRanM,"mRanM/F"); 
  TBranch *mPtMBR     = outTree->Branch("mPtM" ,&mPtM,  "mPtM/F");  
  TBranch *mSumPtMBR  = outTree->Branch("mSumPtM" ,&mSumPtM,  "mSumPtM/F"); 
  TBranch *mPtoMMBR   = outTree->Branch("mPtoMM" ,&mPtoMM,  "mPtoMM/F"); 
  TBranch *mClosMBR   = outTree->Branch("mClosM" ,&mClosM,  "mClosM/F");
  TBranch *mClos2MBR   = outTree->Branch("mClos2M" ,&mClos2M,  "mClos2M/F");    


  TBranch *mRanTBR    = outTree->Branch("mRanT",&mRanT,"mRanT/F"); 
  TBranch *mPtTBR     = outTree->Branch("mPtT" ,&mPtT,  "mPtT/F");  
  TBranch *mSumPtTBR  = outTree->Branch("mSumPtT" ,&mSumPtT,  "mSumPtT/F"); 
  TBranch *mPtoMTBR   = outTree->Branch("mPtoMT" ,&mPtoMT,  "mPtoMT/F");  
  TBranch *mClosTBR   = outTree->Branch("mClosT" ,&mClosT,  "mClosT/F");
  TBranch *mClos2TBR   = outTree->Branch("mClos2T" ,&mClos2T,  "mClos2T/F");   


  TBranch *BDTBR      = outTree->Branch("BDT" ,&BDT,  "BDT/F"); 

  reader->AddVariable( "nSvs",         &nSvs);
  reader->AddVariable( "numJets20",    &numJets20);
  reader->AddVariable( "numJets30",    &numJets30);
  reader->AddVariable( "numJets60",    &numJets60);
  reader->AddVariable( "numJets100",   &numJets100);
  reader->AddVariable( "numJetsBtagL", &numJetsBtagL);
  reader->AddVariable( "numJetsBtagM", &numJetsBtagM);
  reader->AddVariable( "numJetsBtagT", &numJetsBtagT);

  reader->BookMVA( "BDT", TString("weights/")+weightFile ); 





  TFile* fInput = 0;
  if(sgn) 
    fInput = TFile::Open("gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_VType2/v2/DiJetPt_TTH_HToBB_M-125_8TeV-pythia6_VType2.root","READ");
  else 
    fInput = TFile::Open("gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_VType2/v2/DiJetPt_TTJets_SemiLeptMGDecays_8TeV-madgraph-part_VType2.root","READ");
  TTree* tree = (TTree*)fInput->Get("tree");

  if(!tree){
    cout << "No tree" << endl; return ;
  }

  JetByPt jet1, jet2, jet3, jet4, jet5, jet6, jet7, jet8, jet9, jet10;
  genParticleInfo genB, genBbar;
  int nSV;

  tree->SetBranchAddress("jet1",   &jet1);
  tree->SetBranchAddress("jet2",   &jet2);
  tree->SetBranchAddress("jet3",   &jet3);
  tree->SetBranchAddress("jet4",   &jet4);
  tree->SetBranchAddress("jet5",   &jet5);
  tree->SetBranchAddress("jet6",   &jet6);
  tree->SetBranchAddress("jet7",   &jet7);
  tree->SetBranchAddress("jet8",   &jet8);
  tree->SetBranchAddress("jet9",   &jet9);
  tree->SetBranchAddress("jet10",  &jet10);
  tree->SetBranchAddress("genB",   &genB);
  tree->SetBranchAddress("genBbar",&genBbar);
  tree->SetBranchAddress("nSvs",   &nSV);
  

  Long64_t nentries = tree->GetEntries();
  //if(!sgn)   nentries = 500000;

  for (Long64_t i = 0; i < nentries ; i++){
    
    if(i%10000==0) cout << i << endl;
    tree->GetEntry(i);

    LV genBLV(0.,0.,0.,0.);
    LV genBbarLV(0.,0.,0.,0.);
    if(genB.mass>0 && genB.momid==25){
      genBLV.SetPt(  genB.pt );
      genBLV.SetEta( genB.eta );
      genBLV.SetPhi( genB.phi );
      genBLV.SetM(   genB.mass );
    }
    if(genBbar.mass>0 && genBbar.momid==25){
      genBbarLV.SetPt(  genBbar.pt );
      genBbarLV.SetEta( genBbar.eta );
      genBbarLV.SetPhi( genBbar.phi );
      genBbarLV.SetM(   genBbar.mass );
    }
    
    vector<LV> myJets;
    std::map<unsigned int, JetByPt> map;
  
    map[0] = jet1;  map[1] = jet2; map[2] = jet3; map[3] = jet4; 
    map[4] = jet5;  map[5] = jet6; map[6] = jet7; map[7] = jet8; 
    map[8] = jet9;  map[9] = jet10; 
    
    myJets.push_back(LV(jet1.pt, jet1.eta, jet1.phi, jet1.mass));
    myJets.push_back(LV(jet2.pt, jet2.eta, jet2.phi, jet2.mass));
    myJets.push_back(LV(jet3.pt, jet3.eta, jet3.phi, jet3.mass));
    myJets.push_back(LV(jet4.pt, jet4.eta, jet4.phi, jet4.mass));
    myJets.push_back(LV(jet5.pt, jet5.eta, jet5.phi, jet5.mass));
    myJets.push_back(LV(jet6.pt, jet6.eta, jet6.phi, jet6.mass));
    myJets.push_back(LV(jet7.pt, jet7.eta, jet7.phi, jet7.mass));
    myJets.push_back(LV(jet8.pt, jet8.eta, jet8.phi, jet8.mass));
    myJets.push_back(LV(jet9.pt, jet9.eta, jet9.phi, jet9.mass));
    myJets.push_back(LV(jet10.pt,jet10.eta,jet10.phi,jet10.mass));
  
  
    vector<LV> myJetsFilt;
    std::map<unsigned int, JetByPt> mapFilt;
    unsigned int iter = 0;
    
    for(unsigned int k = 0; k<myJets.size() ; k++){
      
      if(myJets[k].Pt()>20 && TMath::Abs(myJets[k].Eta()) < 2.5 ){
	
	myJetsFilt.push_back( myJets[k] );
	
	switch(k){
	case 0:
	  mapFilt[iter] = jet1;
	  break;
	case 1:
	  mapFilt[iter] = jet2;
	  break;
	case 2:
	  mapFilt[iter] = jet3;
	  break;
	case 3:
	  mapFilt[iter] = jet4;
	  break;
	case 4:
	  mapFilt[iter] = jet5;
	  break;
	case 5:
	  mapFilt[iter] = jet6;
	  break;
	case 6:
	  mapFilt[iter] = jet7;
	  break;
	case 7:
	  mapFilt[iter] = jet8;
	  break;
	case 8:
	  mapFilt[iter] = jet9;
	  break;
	case 9:
	  mapFilt[iter] = jet10;
	  break;
	default:
	  break;
	}
	
	iter++;
      }
    }

    numJets20 =0; 
    numJets30 =0; 
    numJets60 =0; 
    numJets100=0; 
    numJetsBtagL=0; 
    numJetsBtagM=0; 
    numJetsBtagT=0;

    int matchBL=0;   int matchBM=0;   int matchBT=0;
    int matchBbarL=0;int matchBbarM=0;int matchBbarT=0;
    

    for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
      float csv_k = mapFilt[k].csv>0 ? mapFilt[k].csv : 0.0 ;
      float pt_k  = mapFilt[k].pt>0  ? mapFilt[k].pt  : 0.0 ;

      if( pt_k>30 && csv_k>0.244 ) numJetsBtagL++;
      if( pt_k>30 && csv_k>0.679 ) numJetsBtagM++;
      if( pt_k>30 && csv_k>0.898 ) numJetsBtagT++;
      if( pt_k>20     ) numJets20++;
      if( pt_k>30     ) numJets30++;
      if( pt_k>60     ) numJets60++;
      if( pt_k>100    ) numJets100++;


      bool isB    = genBLV.Pt()>0 ?    deltaR( myJetsFilt[k], genBLV)    < 0.3 && TMath::Abs(myJetsFilt[k].Pt()-genBLV.Pt())/genBLV.Pt()<0.30        : false;
      bool isBbar = genBbarLV.Pt()>0 ? deltaR( myJetsFilt[k], genBbarLV) < 0.3 && TMath::Abs(myJetsFilt[k].Pt()-genBbarLV.Pt())/genBbarLV.Pt()<0.30  : false;
      if( isB){ 
	if(pt_k>30 && csv_k>0.244) matchBL++;
	if(pt_k>30 && csv_k>0.679) matchBM++;
	if(pt_k>30 && csv_k>0.898) matchBT++;
      }
      if( isBbar) {
	if(pt_k>30 && csv_k>0.244) matchBbarL++;
	if(pt_k>30 && csv_k>0.679) matchBbarM++;
	if(pt_k>30 && csv_k>0.898) matchBbarT++;
      }
      
    }

    float combL      = numJetsBtagL*(numJetsBtagL-1)/2.;
    matchesL = (matchBL+matchBbarL)/2.;
    weightL = 1.0;
    if(matchBL+matchBbarL>0){
      if( combL <1 ) 
	weightL = LOWVALUE;
      else
	weightL = (combL - matchesL)>0 ? matchesL*1./(combL - matchesL) : 1.0;
    }
    else
      weightL = LOWVALUE;

    float combM      = numJetsBtagM*(numJetsBtagM-1)/2.;
    matchesM = (matchBM+matchBbarM)/2.;
    weightM = 1.0;
    if(matchBM+matchBbarM>0){
      if( combM <1 ) 
	weightM = LOWVALUE;
      else
	weightM = (combM - matchesM)>0 ? matchesM*1./(combM - matchesM): 1.0;
    }
    else
      weightM = LOWVALUE;

    float combT      = numJetsBtagT*(numJetsBtagT-1)/2.;
    matchesT = (matchBT+matchBbarT)/2.;
    weightT = 1.0;
    if(matchBT+matchBbarT>0){
      if( combT <1 ) 
	weightT = LOWVALUE;
      else
	weightT = (combT - matchesT)>0 ? matchesT*1./(combT - matchesT): 1.0;
    }
    else
      weightT = LOWVALUE;

    nSvs= nSV;


    m12L=-99;m13L=-99;m14L=-99;m23L=-99;m24L=-99;m34L=-99;
    m12M=-99;m13M=-99;m14M=-99;m23M=-99;m24M=-99;m34M=-99;
    m12T=-99;m13T=-99;m14T=-99;m23T=-99;m24T=-99;m34T=-99;

    mRanL   =-99;  mRanM    =-99;mRanT    =-99;
    mPtL    =-99;  mPtM     =-99;mPtT     =-99;
    mSumPtL =-99;  mSumPtM  =-99;mSumPtT  =-99;
    mPtoML  =-99;  mPtoMM   =-99;mPtoMT   =-99;
    mClosL  =-99;  mClosM   =-99;mClosT  =-99;
    mClos2L  =-99;  mClos2M   =-99;mClos2T  =-99;

    nClosL=0; nClosM=0; nClosT=0;

    int counterL = 0;
    int counterM = 0;
    int counterT = 0;

    if(myJetsFilt.size()>1){


      float PtMaxL=-999;
      float PtMaxM=-999;
      float PtMaxT=-999;
      float SumPtMaxL=-999;
      float SumPtMaxM=-999;
      float SumPtMaxT=-999;
      float PtoML=-999;
      float PtoMM=-999;
      float PtoMT=-999;


      for(unsigned int k = 0; k < myJetsFilt.size()-1; k++){  
	float csv_k = mapFilt[k].csv>0 ? mapFilt[k].csv : 0.0 ;
	float pt_k  = mapFilt[k].pt>0  ? mapFilt[k].pt  : 0.0 ;
	
	for(unsigned int l = k+1; l < myJetsFilt.size(); l++){  
	  float csv_l = mapFilt[l].csv>0 ? mapFilt[l].csv : 0.0 ;
	  float pt_l  = mapFilt[l].pt>0  ? mapFilt[l].pt  : 0.0 ;
	  
	  LV diJet   = myJetsFilt[k]+myJetsFilt[l];
	  float mass = diJet.M();
	  float pt   = diJet.Pt();
	  float sumpt= myJetsFilt[k].Pt()+myJetsFilt[l].Pt();
	  float ptoM = diJet.Pt()/diJet.M();
	  

	  if( csv_k>0.244 && csv_l>0.244 && pt_k>30 && pt_l>30 ){

	    if( mClosL==-99  &&  mass>105 && mass<135){
	      mClosL = mass;
	      //nClosL++;
	    }
	    if( mass>105 && mass<135) nClosL++;
	    if( mClos2L==-99 &&  mass>110 && mass<130) mClos2L = mass;

	    switch(counterL){
	    case 0:
	      m12L = mass;
	      break;
	    case 1:
	      m13L = mass;
	      break;
	    case 2:
	      m14L = mass;
	      break;
	    case 3:
	      m23L = mass;
	      break;
	    case 4:
	      m24L = mass;
	      break;
	    case 5:
	      m34L = mass;
	      break;
	    default:
	      break;
	    }

	    if( pt>PtMaxL){
	      PtMaxL = pt;
	      mPtL   = mass;
	    }
	    if( sumpt>SumPtMaxL){
	      SumPtMaxL = sumpt;
	      mSumPtL   = mass;
	    }
	    if( ptoM>PtoML){
	      PtoML     = ptoM;
	      mPtoML    = mass;
	    }

	    counterL++;
	  }
	  

	  if( csv_k>0.679 && csv_l>0.679 && pt_k>30 && pt_l>30 ){

	    if( mClosM==-99 &&  mass>105 && mass<135){
	      mClosM = mass;
	      //nClosM++;
	    }
	    if(mass>105 && mass<135) nClosM++;
	    if( mClos2M==-99 &&  mass>110 && mass<130) mClos2M = mass;

	    switch(counterM){
	    case 0:
	      m12M = mass;
	      break;
	    case 1:
	      m13M = mass;
	      break;
	    case 2:
	      m14M = mass;
	      break;
	    case 3:
	      m23M = mass;
	      break;
	    case 4:
	      m24M = mass;
	      break;
	    case 5:
	      m34M = mass;
	      break;
	    default:
	      break;
	    }

	    if( pt>PtMaxM){
	      PtMaxM = pt;
	      mPtM   = mass;
	    }
	    if( sumpt>SumPtMaxM){
	      SumPtMaxM = sumpt;
	      mSumPtM   = mass;
	    }
	    if( ptoM>PtoMM){
	      PtoMM     = ptoM;
	      mPtoMM    = mass;
	    }

	    counterM++;
	  }
	  
	  if( csv_k>0.898 && csv_l>0.898 && pt_k>30 && pt_l>30 ){

	    if( mClosT==-99 &&  mass>105 && mass<135){
	      mClosT = mass;
	      //nClosT++;
	    }
	    if(mass>105 && mass<135) nClosT++;
	    if( mClos2T==-99 &&  mass>110 && mass<130) mClos2T = mass;

	    switch(counterT){
	    case 0:
	      m12T = mass;
	      break;
	    case 1:
	      m13T = mass;
	      break;
	    case 2:
	      m14T = mass;
	      break;
	    case 3:
	      m23T = mass;
	      break;
	    case 4:
	      m24T = mass;
	      break;
	    case 5:
	      m34T = mass;
	      break;
	    default:
	      break;
	    }

	    if( pt>PtMaxT){
	      PtMaxT = pt;
	      mPtT   = mass;
	    }
	    if( sumpt>SumPtMaxT){
	      SumPtMaxT = sumpt;
	      mSumPtT   = mass;
	    }
	    if( ptoM>PtoMT){
	      PtoMT     = ptoM;
	      mPtoMT    = mass;
	    }

	    counterT++;
	  }



	}
      }
    }

    //typeBR->Fill();
    //m12BR->Fill();
    //m13BR->Fill();
    //m14BR->Fill();
    //m23BR->Fill();
    //m24BR->Fill();
    //m34BR->Fill();
    int choiceL = counterL>0 ? int(ran->Uniform(0,counterL)) : -99;
    switch(choiceL){
    case 0:
      mRanL = m12L;
      break;
    case 1:
      mRanL = m13L;
      break;
    case 2:
      mRanL = m14L;
      break;
    case 3:
      mRanL = m23L;
      break;
    case 4:
      mRanL = m24L;
      break;
    case 5:
      mRanL = m34L;
      break;
    default:
      break;
    }
    if(mClosL==-99)  mClosL  = mPtoML;
    if(mClos2L==-99) mClos2L = mPtoML;

    int choiceM = counterM>0 ? int(ran->Uniform(0,counterM)): -99;
    switch(choiceM){
    case 0:
      mRanM = m12M;
      break;
    case 1:
      mRanM = m13M;
      break;
    case 2:
      mRanM = m14M;
      break;
    case 3:
      mRanM = m23M;
      break;
    case 4:
      mRanM = m24M;
      break;
    case 5:
      mRanM = m34M;
      break;
    default:
      break;
    }
    if(mClosM==-99)  mClosM  = mPtoMM;
    if(mClos2M==-99) mClos2M = mPtoMM;

    int choiceT = counterT>0 ? int(ran->Uniform(0,counterT)): -99;
    switch(choiceT){
    case 0:
      mRanT = m12T;
      break;
    case 1:
      mRanT = m13T;
      break;
    case 2:
      mRanT = m14T;
      break;
    case 3:
      mRanT = m23T;
      break;
    case 4:
      mRanT = m24T;
      break;
    case 5:
      mRanT = m34T;
      break;
    default:
      break;
    }
    if(mClosT==-99)  mClosT  = mPtoMT;
    if(mClos2T==-99) mClos2T = mPtoMT;

    type = 0;

    BDT = reader->EvaluateMVA( "BDT" );
 
    outTree->Fill();
  }
  
  fOutput->cd();
  outTree->Write("",TObject::kOverwrite );
  fOutput->Close();

  delete ran; delete reader;

  return;

}



void analyzeTree(){


  TFile* fInputSgn = 0;
  fInputSgn = TFile::Open("gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_VType2/v2/DiJetPt_TTH_HToBB_M-125_8TeV-pythia6_VType2.root","READ");
  TTree* tree = (TTree*)fInputSgn->Get("tree");
  if(!tree){
    cout << "No tree" << endl; return ;
  }
  TH1F* hCounter = (TH1F*)fInputSgn->Get("Count");
  float wSgn = (0.1302*0.569)*12.1*1000/hCounter->GetBinContent(1);
  wSgn*=1;

  TFile* fInputBkg = 0;
  fInputBkg = TFile::Open("gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_VType2/v2/DiJetPt_TTJets_SemiLeptMGDecays_8TeV-madgraph-part_VType2.root","READ");
  tree = (TTree*)fInputBkg->Get("tree");
  if(!tree){
    cout << "No tree" << endl; return ;
  }
  hCounter = (TH1F*)fInputBkg->Get("Count");
  float wBkg = (103.0)*12.1*1000/hCounter->GetBinContent(1);

  cout << "Weight signal     = " << wSgn << endl;
  cout << "Weight bsckground = " << wBkg << endl;

  TFile* fSgn = TFile::Open("SgnTree.root","READ");
  TFile* fBkg = TFile::Open("BkgTree.root","READ");

  TTree* tSgn = (TTree*)fSgn->Get("outTree");
  TTree* tBkg = (TTree*)fBkg->Get("outTree");

  /////////////////////////////////////////////////////

  int mult=4;
  TCut cat(Form("numJetsBtagM==%d && numJetsBtagT>=0", mult));

  TH1F* hRan = new TH1F("hRan","",1,0,1000);
  TH1F* hTmp = (TH1F*)hRan->Clone("hTmp");
  hTmp->Reset();

  tBkg->Draw("mRanM>>hTmp",cat&&TCut("mRanM>105 && mRanM<135"));
  hTmp->Scale(wBkg);
  hRan->Add(hTmp,1.0);  
  hTmp->Reset();
  tSgn->Draw("mRanM>>hTmp",cat&&TCut("mRanM>105 && mRanM<135"));
  hTmp->Scale(wSgn);
  //hRan->Add(hTmp,1.0);  
  hTmp->Reset();

  float sideband = hRan->Integral();
  cout << "Num events in sideband = " << sideband << endl;
  hRan->Reset();

  tBkg->Draw("mRanM>>hTmp",cat);
  hTmp->Scale(wBkg);
  hRan->Add(hTmp,1.0);  
  hTmp->Reset();
  tSgn->Draw("mRanM>>hTmp",cat);
  hTmp->Scale(wSgn);
  //hRan->Add(hTmp,1.0);  
  hTmp->Reset();

  float total = hRan->Integral();
  cout << "Num events in total = " << total << endl;
  hRan->Reset();

  float ratio    = sideband/total;
  float ratioErr = TMath::Sqrt(ratio*(1-ratio)/total); 
  cout << "Ratio = " << ratio << " +/- " << ratioErr << endl;

  float bkgInSgn = tBkg->GetEntries(cat&&TCut("nClosM>0"));
  bkgInSgn*=wBkg;
  float sgnInSgn = tSgn->GetEntries(cat&&TCut("nClosM>0"));
  sgnInSgn*=wSgn;

  cout << "Total number of events with >=1 comb in [105,135] = " << bkgInSgn+sgnInSgn << endl;
  cout << "(Sgn = " << sgnInSgn << ", bkg = " << bkgInSgn << ")" << endl;

  float expBkg     = total*( 1 - TMath::Power(1-ratio,          mult*(mult-1)/2 ));
  float expBkgUp   = total*( 1 - TMath::Power(1-ratio+ratioErr, mult*(mult-1)/2 ));
  float expBkgDown = total*( 1 - TMath::Power(1-ratio-ratioErr, mult*(mult-1)/2 ));
  float expErr = TMath::Max(TMath::Abs(expBkg-expBkgUp),TMath::Abs(expBkg-expBkgDown));

  cout << "Expected background is " << expBkg << " +/-" << expErr << endl;
  cout << "To be compared with " << bkgInSgn << " +/-" << TMath::Sqrt(bkgInSgn/wBkg)*wBkg << endl;
  
}






void secondaryVertex(){

  TFile *fOutput           = new TFile("secondaryVertex.root","RECREATE");
  TH1F* hVtxMassR = new TH1F("hVtxMassR","",100,0,200);
  TH1F* hVtxMassW = new TH1F("hVtxMassW","",100,0,200);

  TH1F* hMass = new TH1F("hMass","",100,0,400);
  TH1F* hCount         = new TH1F("hCount", "visited"   ,10, 0 ,10 );
  TH1F* hMult         = new TH1F("hMult", "visited"   ,10, 0 ,10 );


  TFile* fInput = TFile::Open("gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_VType2/v2/DiJetPt_TTH_HToBB_M-125_8TeV-pythia6_VType2.root","READ");
  //TFile* fInput = TFile::Open("gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_VType2/v2/DiJetPt_TTJets_SemiLeptMGDecays_8TeV-madgraph-part_VType2.root","READ");
  TTree* tree = (TTree*)fInput->Get("tree");

  if(!tree){
    cout << "No tree" << endl; return ;
  }

  TH1F* hCounter = (TH1F*)fInput->Get("Count");
  float weight = (0.1302*0.569)*12.1*1000/hCounter->GetBinContent(1);
  //float weight = (103.0)*12.1*1000/hCounter->GetBinContent(1);
  cout << weight << endl;
  //weight /= (500000./tree->GetEntries());
  cout << weight << endl;

  JetByPt jet1, jet2, jet3, jet4, jet5, jet6, jet7, jet8, jet9, jet10;
  genTopInfo genTop, genTbar;
  genParticleInfo genB, genBbar;
  int nvlep;
  float vLeptonpt[999]; float vLeptoneta[999]; float vLeptonphi[999];
  float met[999];

  int nSvs;
  float Sv_massBCand[999];
  float Sv_massSv[999];
  float Sv_pt[999];
  float Sv_eta[999];
  float Sv_phi[999];
  float Sv_dist3D[999];
  float Sv_dist2D[999];
  float Sv_distSim2D[999];
  float Sv_distSig3D[999];
  float Sv_dist3D_norm[999];

  
  tree->SetBranchAddress("jet1",   &jet1);
  tree->SetBranchAddress("jet2",   &jet2);
  tree->SetBranchAddress("jet3",   &jet3);
  tree->SetBranchAddress("jet4",   &jet4);
  tree->SetBranchAddress("jet5",   &jet5);
  tree->SetBranchAddress("jet6",   &jet6);
  tree->SetBranchAddress("jet7",   &jet7);
  tree->SetBranchAddress("jet8",   &jet8);
  tree->SetBranchAddress("jet9",   &jet9);
  tree->SetBranchAddress("jet10",  &jet10);
  tree->SetBranchAddress("genTop", &genTop);
  tree->SetBranchAddress("genTbar",&genTbar);
  tree->SetBranchAddress("genB",   &genB);
  tree->SetBranchAddress("genBbar",&genBbar);
  tree->SetBranchAddress("vLepton_pt", vLeptonpt);
  tree->SetBranchAddress("vLepton_phi",vLeptonphi);
  tree->SetBranchAddress("vLepton_eta",vLeptoneta);
  tree->SetBranchAddress("nvlep",  &nvlep);
  tree->SetBranchAddress("METtype1p2corr",  met);

  tree->SetBranchAddress("nSvs",&nSvs);
  tree->SetBranchAddress("Sv_massBCand", Sv_massBCand);
  tree->SetBranchAddress("Sv_massSv", Sv_massSv);
  tree->SetBranchAddress("Sv_pt", Sv_pt);
  tree->SetBranchAddress("Sv_eta", Sv_eta);
  tree->SetBranchAddress("Sv_phi", Sv_phi);
  tree->SetBranchAddress("Sv_dist3D", Sv_dist3D);
  tree->SetBranchAddress("Sv_dist2D",Sv_dist2D);
  tree->SetBranchAddress("Sv_distSim2D", Sv_distSim2D);
  tree->SetBranchAddress("Sv_distSig3D", Sv_distSig3D);
  tree->SetBranchAddress("Sv_dist3D_norm", Sv_dist3D_norm);


  Long64_t nentries = tree->GetEntries();
  int totCounter = 0;
  
  for (Long64_t i = 0; i < /*500000*/ nentries ; i++){
    
    if(i%10000==0) cout << i << endl;
    tree->GetEntry(i);
    
    //continue;

    if(nvlep<=0) continue;

    LV topBLV(  0.,0.,0.,0.);
    LV topW1LV( 0.,0.,0.,0.); 
    LV topW2LV( 0.,0.,0.,0.);
    LV atopBLV( 0.,0.,0.,0.); 
    LV atopW1LV(0.,0.,0.,0.); 
    LV atopW2LV(0.,0.,0.,0.);
    LV genBLV(0.,0.,0.,0.);
    LV genBbarLV(0.,0.,0.,0.);
    
    
     ///////////////////////////////////////////////////////////////////
     /////////////////////////////////////////////////////////////////// 
     if(genTop.bmass>0){
       topBLV.SetPt(   genTop.bpt );
       topBLV.SetEta(  genTop.beta );
       topBLV.SetPhi(  genTop.bphi );
       topBLV.SetM(    genTop.bmass );
       topW1LV.SetPt(  genTop.wdau1pt );
       topW1LV.SetEta( genTop.wdau1eta );
       topW1LV.SetPhi( genTop.wdau1phi );
       topW1LV.SetM(   genTop.wdau1mass );
       topW2LV.SetPt(  genTop.wdau2pt );
       topW2LV.SetEta( genTop.wdau2eta );
       topW2LV.SetPhi( genTop.wdau2phi );
       topW2LV.SetM(   genTop.wdau2mass );
     }
     if(genTbar.bmass>0){
	atopBLV.SetPt(   genTbar.bpt );
	atopBLV.SetEta(  genTbar.beta );
	atopBLV.SetPhi(  genTbar.bphi );
	atopBLV.SetM(    genTbar.bmass );
	atopW1LV.SetPt(  genTbar.wdau1pt );
	atopW1LV.SetEta( genTbar.wdau1eta );
	atopW1LV.SetPhi( genTbar.wdau1phi );
	atopW1LV.SetM(   genTbar.wdau1mass );
	atopW2LV.SetPt(  genTbar.wdau2pt );
	atopW2LV.SetEta( genTbar.wdau2eta );
	atopW2LV.SetPhi( genTbar.wdau2phi );
	atopW2LV.SetM(   genTbar.wdau2mass );
     }

    if(genB.mass>0 && genB.momid==25){
       genBLV.SetPt(  genB.pt );
       genBLV.SetEta( genB.eta );
       genBLV.SetPhi( genB.phi );
       genBLV.SetM(   genB.mass );
     }
     if(genBbar.mass>0 && genBbar.momid==25){
       genBbarLV.SetPt(  genBbar.pt );
       genBbarLV.SetEta( genBbar.eta );
       genBbarLV.SetPhi( genBbar.phi );
       genBbarLV.SetM(   genBbar.mass );
     }

     vector<LV> myJets;
     std::map<unsigned int, JetByPt> map;

     map[0] = jet1;  map[1] = jet2; map[2] = jet3; map[3] = jet4; 
     map[4] = jet5;  map[5] = jet6; map[6] = jet7; map[7] = jet8; 
     map[8] = jet9;  map[9] = jet10; 

     myJets.push_back(LV(jet1.pt, jet1.eta, jet1.phi, jet1.mass));
     myJets.push_back(LV(jet2.pt, jet2.eta, jet2.phi, jet2.mass));
     myJets.push_back(LV(jet3.pt, jet3.eta, jet3.phi, jet3.mass));
     myJets.push_back(LV(jet4.pt, jet4.eta, jet4.phi, jet4.mass));
     myJets.push_back(LV(jet5.pt, jet5.eta, jet5.phi, jet5.mass));
     myJets.push_back(LV(jet6.pt, jet6.eta, jet6.phi, jet6.mass));
     myJets.push_back(LV(jet7.pt, jet7.eta, jet7.phi, jet7.mass));
     myJets.push_back(LV(jet8.pt, jet8.eta, jet8.phi, jet8.mass));
     myJets.push_back(LV(jet9.pt, jet9.eta, jet9.phi, jet9.mass));
     myJets.push_back(LV(jet10.pt,jet10.eta,jet10.phi,jet10.mass));

     
     vector<LV> myJetsFilt;
     std::map<unsigned int, JetByPt> mapFilt;
     unsigned int iter = 0;

     for(unsigned int k = 0; k<myJets.size() ; k++){

       if(myJets[k].Pt()>30 && TMath::Abs(myJets[k].Eta()) < 2.5 ){

	 myJetsFilt.push_back( myJets[k] );

	 switch(k){
	 case 0:
	   mapFilt[iter] = jet1;
	   break;
	 case 1:
	   mapFilt[iter] = jet2;
	   break;
	 case 2:
	   mapFilt[iter] = jet3;
	   break;
	 case 3:
	   mapFilt[iter] = jet4;
	   break;
	 case 4:
	   mapFilt[iter] = jet5;
	   break;
	 case 5:
	   mapFilt[iter] = jet6;
	   break;
	 case 6:
	   mapFilt[iter] = jet7;
	   break;
	 case 7:
	   mapFilt[iter] = jet8;
	   break;
	 case 8:
	   mapFilt[iter] = jet9;
	   break;
	 case 9:
	   mapFilt[iter] = jet10;
	   break;
	 default:
	   break;
	 }
	 
	 iter++;
       }
     }

     LV lepLV(vLeptonpt[0], vLeptoneta[0], vLeptonphi[0], 0.);
     LV MEtLV(met[0], 0., met[3] , 0.);

     int match1=0;
     int match2=0;
     int match3=0;
     int match4=0;

     int vertexB1=0;
     int vertexB2=0;
     int vertexB3=0;
     int vertexB4=0;

     for(int l = 0; l < nSvs; l++){
       LV vertex( Sv_pt[l], Sv_eta[l], Sv_phi[l], Sv_massSv[l]);
       if(deltaR(vertex, genBLV)<0.3){
	 match1++;
	 //phi1   = vertex.Phi();
	 //eta1   = vertex.Eta();
	 //dist13D = Sv_dist3D[l];
	 vertexB1 = l;
       }
       if(deltaR(vertex, genBbarLV)<0.3){
	 match2++;
	 //phi2   = vertex.Phi();
	 //eta2   = vertex.Eta();
	 //dist23D = Sv_dist3D[l];
	 vertexB2 = l;
       }
       if(deltaR(vertex, topBLV)<0.3){
	 match3++;
	 vertexB3 = l;
       }
       if(deltaR(vertex, atopBLV)<0.3){
	 match4++;
	 vertexB4 = l;
       }
     }
     hCount->Fill( match1+match2+match3+match4, weight);



     if( match1>0 && match2>0){
       LV vertex1( Sv_pt[vertexB1], Sv_eta[vertexB1], Sv_phi[vertexB1], Sv_massSv[vertexB1]);
       LV vertex2( Sv_pt[vertexB2], Sv_eta[vertexB2], Sv_phi[vertexB2], Sv_massSv[vertexB2]);

       float mass = (vertex1+vertex2).M();
       hVtxMassR->Fill(mass,weight);

     }
     if( match3>0 && match4>0){
       LV vertex1( Sv_pt[vertexB3], Sv_eta[vertexB3], Sv_phi[vertexB3], Sv_massSv[vertexB3]);
       LV vertex2( Sv_pt[vertexB4], Sv_eta[vertexB4], Sv_phi[vertexB4], Sv_massSv[vertexB4]);

       float mass = (vertex1+vertex2).M();
       hVtxMassW->Fill(mass,weight);

     }


     int numBTag=0;
     for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
       float csv_k = mapFilt[k].csv>0 ? mapFilt[k].csv : 0.0 ;
       if(csv_k>0.679) numBTag++;
     }
     
     std::map<float,  float , sorterByLess > matchResult;

     int mult=0;
     if(nSvs>1 && myJetsFilt.size()>=5 && numBTag>3){
       
       for(int k = 0; k < nSvs-1; k++){
	 int kJet;
	 bool kFound = false;
	 LV vertex1( Sv_pt[k], Sv_eta[k], Sv_phi[k], Sv_massSv[k]);
	 for(unsigned int m = 0; m < myJetsFilt.size(); m++){  
	   if( deltaR(vertex1, myJetsFilt[m] )<0.3 ){
	     kFound=true;
	     kJet = m;
	   }
	 }

	 for(int l = k+1; l < nSvs; l++){
	   int lJet;
	   bool lFound = false;
	   LV vertex2( Sv_pt[l], Sv_eta[l], Sv_phi[l], Sv_massSv[l]);

	   for(unsigned int m = 0; m < myJetsFilt.size(); m++){  
	     if( deltaR(vertex2, myJetsFilt[m] )<0.3 ){
	       lFound=true;
	       lJet = m;
	     }
	   }

	   if((vertex1+vertex2).M()<20 || (vertex1+vertex2).M()>60) continue;

	   if(kFound && lFound && lJet!= kJet){
	     //matchResult[ deltaR(myJetsFilt[lJet], myJetsFilt[kJet])] =
	     matchResult[ deltaR( vertex1,vertex2 )] =
	       (myJetsFilt[lJet]+ myJetsFilt[kJet]).M();

	     mult++;
	   }

	   //float mass = (vertex1+vertex2).M();
	   //hVtxMass->Fill(mass,weight);

	 }
       }
       
     }

     if(matchResult.size()>0) hMass->Fill( (matchResult.begin())->second, weight );
     hMult->Fill(mult,weight);

  }



   /////////////////////////////////////////////////////////////////////////// 
   
   
   fOutput->Write();
   fOutput->Close();
   fInput->Close();
   
  delete fOutput;
  
  return  ;

}


void algoComb( float ptCut_ = 30 , float PtoM_ = 1.5, float csvH1_ = 0.679, float csvH2_ = 0.679, float ptH1_ = 30, float ptH2_ = 30){

  TRandom3* ran = new TRandom3();

  TFile *fOutput           = new TFile("algo.root","RECREATE");
  TH1F* hGoodRank          = new TH1F("hGoodRank", "good pos"   ,100, 0 ,100 );
  TH1F* hMassJJW           = new TH1F("hMassJJW", "Higgs mass"   ,40, 0 ,400 );
  TH1F* hMassJJR           = new TH1F("hMassJJR", "Higgs mass"   ,40, 0 ,400 );

  TH1F* hPtJJW           = new TH1F("hPtJJW", "Higgs mass"   ,30, 0 ,600 );
  TH1F* hPtJJR           = new TH1F("hPtJJR", "Higgs mass"   ,30, 0 ,600 );

  TH1F* hPt1JJR           = new TH1F("hPt1JJR", "Higgs mass"   ,40, 0 ,400 );
  TH1F* hPt2JJR           = new TH1F("hPt2JJR", "Higgs mass"   ,40, 0 ,400 );
  TH1F* hPt1JJW           = new TH1F("hPt1JJW", "Higgs mass"   ,40, 0 ,400 );
  TH1F* hPt2JJW           = new TH1F("hPt2JJW", "Higgs mass"   ,40, 0 ,400 );

  TH1F* hPtoMJJW           = new TH1F("hPtoMJJW", "Higgs mass"   ,50, 0 ,5 );
  TH1F* hPtoMJJR           = new TH1F("hPtoMJJR", "Higgs mass"   ,50, 0 ,5 );

  TH2F* h2PtoMJJR           = new TH2F("h2PtoMJJR", "Higgs mass"   ,40, 0 ,400, 40, 0 ,400);
  TH2F* h2PtoMJJW           = new TH2F("h2PtoMJJW", "Higgs mass"   ,40, 0 ,400, 40, 0 ,400);

  TH2F* h2PtPtR           = new TH2F("h2PtPtR", "Higgs mass"   ,40, 0 ,200, 40, 0 ,200);
  TH2F* h2PtPtW           = new TH2F("h2PtPtW", "Higgs mass"   ,40, 0 ,200, 40, 0 ,200);

  TH2F* h2CsvW              = new TH2F("h2csvW", "Higgs mass"   ,40, 0 ,1.1, 40, 0 ,1.1);
  TH2F* h2CsvR              = new TH2F("h2csvR", "Higgs mass"   ,40, 0 ,1.1, 40, 0 ,1.1);

  TH2F* h2Mass12R            = new TH2F("h2Mass12R", "Higgs mass"   ,40, 0 ,400, 40, 0 ,400);
  TH2F* h2Mass12W            = new TH2F("h2Mass12W", "Higgs mass"   ,40, 0 ,400, 40, 0 ,400);

  TH1F* hLikeGood          = new TH1F("hLikeGood", "likelihood good comb"   ,100, 0 ,100 );
  TH1F* hLikeWrong         = new TH1F("hLikeWrong", "likelihood good comb"   ,100, 0 ,100 );
  TH1F* hNumPermut         = new TH1F("hNumPermut", "visited"   ,1000, 0 ,1000 );

  TH1F* hComb         = new TH1F("hComb", "visited"   ,10, 0 ,10 );
  TH1F* hCombR         = new TH1F("hCombR", "visited"   ,1000, 0 ,1000 );

  TH1F* hCount         = new TH1F("hCount", "visited"   ,10, 0 ,10 );

  TH1F* hMassWHadW          = new TH1F("hMassWHadW",   "Higgs mass"   ,40, 0 ,200);
  TH1F* hMassTopHadW        = new TH1F("hMassTopHadW", "Higgs mass"   ,40, 0 ,400);
  TH1F* hBCsvW              = new TH1F("hBcsvW",    "Higgs mass"   ,40, 0 ,1.1);
  TH1F* hW1CsvW             = new TH1F("hW1csvW",    "Higgs mass"   ,40, 0 ,1.1);
  TH1F* hW2CsvW             = new TH1F("hW2csvW",    "Higgs mass"   ,40, 0 ,1.1);
  TH1F* hBbarCsvW           = new TH1F("hBbarcsvW", "Higgs mass"   ,40, 0 ,1.1);
  TH1F* hMassWHadR          = new TH1F("hMassWHadR",   "Higgs mass"   ,40, 0 ,200);
  TH1F* hMassTopHadR        = new TH1F("hMassTopHadR", "Higgs mass"   ,40, 0 ,400);
  TH1F* hBCsvR              = new TH1F("hBcsvR",    "Higgs mass"   ,40, 0 ,1.1);
  TH1F* hBbarCsvR           = new TH1F("hBbarcsvR", "Higgs mass"   ,40, 0 ,1.1);
  TH1F* hW1CsvR             = new TH1F("hW1csvR",    "Higgs mass"   ,40, 0 ,1.1);
  TH1F* hW2CsvR             = new TH1F("hW2csvR",    "Higgs mass"   ,40, 0 ,1.1);

  TH1F* hPt1oPt2W           = new TH1F("hPt1oPt2W", "Higgs mass"   ,40, 0 ,2 );
  TH1F* hPt1oPt2R           = new TH1F("hPt1oPt2R", "Higgs mass"   ,40, 0 ,2 );


  TFile* fLikelihoods          = new TFile("likelihoods_bkgd.root","READ");
  TF2* fLikelihoodWTopHadrMass = (TF2*)fLikelihoods->Get("fLikelihoodWTopHadrMass");
  TF1* fLikelihoodCsvB         = (TF1*)fLikelihoods->Get("fLikelihoodCsvB");
  TF1* fLikelihoodCsvNonB      = (TF1*)fLikelihoods->Get("fLikelihoodCsvNonB");
  TF1* fLikelihoodTopLeptMass  = (TF1*)fLikelihoods->Get("fLikelihoodTopLeptMass");


  TFile* fInput = TFile::Open("gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_VType2/v2/DiJetPt_TTH_HToBB_M-125_8TeV-pythia6_VType2.root","READ");
  //TFile* fInput = TFile::Open("gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_VType2/v2/DiJetPt_TTJets_SemiLeptMGDecays_8TeV-madgraph-part_VType2.root","READ");
  TTree* tree = (TTree*)fInput->Get("tree");

  if(!tree){
    cout << "No tree" << endl; return ;
  }

  TH1F* hCounter = (TH1F*)fInput->Get("Count");
  float weight = (0.1302*0.569)*12.1*1000/hCounter->GetBinContent(1);
  //float weight = (103.0)*12.1*1000/hCounter->GetBinContent(1);
  cout << weight << endl;
  //weight /= (500000./tree->GetEntries());
  cout << weight << endl;

  JetByPt jet1, jet2, jet3, jet4, jet5, jet6, jet7, jet8, jet9, jet10;
  genTopInfo genTop, genTbar;
  genParticleInfo genB, genBbar;
  int nvlep;
  float vLeptonpt[999]; float vLeptoneta[999]; float vLeptonphi[999];
  float met[999];

  tree->SetBranchAddress("jet1",   &jet1);
  tree->SetBranchAddress("jet2",   &jet2);
  tree->SetBranchAddress("jet3",   &jet3);
  tree->SetBranchAddress("jet4",   &jet4);
  tree->SetBranchAddress("jet5",   &jet5);
  tree->SetBranchAddress("jet6",   &jet6);
  tree->SetBranchAddress("jet7",   &jet7);
  tree->SetBranchAddress("jet8",   &jet8);
  tree->SetBranchAddress("jet9",   &jet9);
  tree->SetBranchAddress("jet10",  &jet10);
  tree->SetBranchAddress("genTop", &genTop);
  tree->SetBranchAddress("genTbar",&genTbar);
  tree->SetBranchAddress("genB",   &genB);
  tree->SetBranchAddress("genBbar",&genBbar);
  tree->SetBranchAddress("vLepton_pt", vLeptonpt);
  tree->SetBranchAddress("vLepton_phi",vLeptonphi);
  tree->SetBranchAddress("vLepton_eta",vLeptoneta);
  tree->SetBranchAddress("nvlep",  &nvlep);
  tree->SetBranchAddress("METtype1p2corr",  met);

  Long64_t nentries = tree->GetEntries();
  int totCounter = 0;
  
  for (Long64_t i = 0; i < /*500000*/ nentries; i++){
    
    if(i%10000==0) cout << i << endl;
    tree->GetEntry(i);
    
    //continue;

    if(nvlep<=0) continue;


     LV topBLV(  0.,0.,0.,0.);
     LV topW1LV( 0.,0.,0.,0.); 
     LV topW2LV( 0.,0.,0.,0.);
     LV atopBLV( 0.,0.,0.,0.); 
     LV atopW1LV(0.,0.,0.,0.); 
     LV atopW2LV(0.,0.,0.,0.);
     LV genBLV(0.,0.,0.,0.);
     LV genBbarLV(0.,0.,0.,0.);
    
    
     ///////////////////////////////////////////////////////////////////
     /////////////////////////////////////////////////////////////////// 
     if(genTop.bmass>0){
       topBLV.SetPt(   genTop.bpt );
       topBLV.SetEta(  genTop.beta );
       topBLV.SetPhi(  genTop.bphi );
       topBLV.SetM(    genTop.bmass );
       topW1LV.SetPt(  genTop.wdau1pt );
       topW1LV.SetEta( genTop.wdau1eta );
       topW1LV.SetPhi( genTop.wdau1phi );
       topW1LV.SetM(   genTop.wdau1mass );
       topW2LV.SetPt(  genTop.wdau2pt );
       topW2LV.SetEta( genTop.wdau2eta );
       topW2LV.SetPhi( genTop.wdau2phi );
       topW2LV.SetM(   genTop.wdau2mass );
     }
     if(genTbar.bmass>0){
	atopBLV.SetPt(   genTbar.bpt );
	atopBLV.SetEta(  genTbar.beta );
	atopBLV.SetPhi(  genTbar.bphi );
	atopBLV.SetM(    genTbar.bmass );
	atopW1LV.SetPt(  genTbar.wdau1pt );
	atopW1LV.SetEta( genTbar.wdau1eta );
	atopW1LV.SetPhi( genTbar.wdau1phi );
	atopW1LV.SetM(   genTbar.wdau1mass );
	atopW2LV.SetPt(  genTbar.wdau2pt );
	atopW2LV.SetEta( genTbar.wdau2eta );
	atopW2LV.SetPhi( genTbar.wdau2phi );
	atopW2LV.SetM(   genTbar.wdau2mass );
     }
     if(genB.mass>0 && genB.momid==25){
       genBLV.SetPt(  genB.pt );
       genBLV.SetEta( genB.eta );
       genBLV.SetPhi( genB.phi );
       genBLV.SetM(   genB.mass );
     }
     if(genBbar.mass>0 && genBbar.momid==25){
       genBbarLV.SetPt(  genBbar.pt );
       genBbarLV.SetEta( genBbar.eta );
       genBbarLV.SetPhi( genBbar.phi );
       genBbarLV.SetM(   genBbar.mass );
     }
     ///////////////////////////////////////////////////////////////////
     ///////////////////////////////////////////////////////////////////

     vector<LV> myJets;
     std::map<unsigned int, JetByPt> map;

     map[0] = jet1;  map[1] = jet2; map[2] = jet3; map[3] = jet4; 
     map[4] = jet5;  map[5] = jet6; map[6] = jet7; map[7] = jet8; 
     map[8] = jet9;  map[9] = jet10; 

     myJets.push_back(LV(jet1.pt, jet1.eta, jet1.phi, jet1.mass));
     myJets.push_back(LV(jet2.pt, jet2.eta, jet2.phi, jet2.mass));
     myJets.push_back(LV(jet3.pt, jet3.eta, jet3.phi, jet3.mass));
     myJets.push_back(LV(jet4.pt, jet4.eta, jet4.phi, jet4.mass));
     myJets.push_back(LV(jet5.pt, jet5.eta, jet5.phi, jet5.mass));
     myJets.push_back(LV(jet6.pt, jet6.eta, jet6.phi, jet6.mass));
     myJets.push_back(LV(jet7.pt, jet7.eta, jet7.phi, jet7.mass));
     myJets.push_back(LV(jet8.pt, jet8.eta, jet8.phi, jet8.mass));
     myJets.push_back(LV(jet9.pt, jet9.eta, jet9.phi, jet9.mass));
     myJets.push_back(LV(jet10.pt,jet10.eta,jet10.phi,jet10.mass));

     
     vector<LV> myJetsFilt;
     std::map<unsigned int, JetByPt> mapFilt;

     unsigned int iter = 0;
     for(unsigned int k = 0; k<myJets.size() ; k++){

       if(myJets[k].Pt()>ptCut_ && TMath::Abs(myJets[k].Eta()) < 2.5 ){

	 myJetsFilt.push_back( myJets[k] );

	 switch(k){
	 case 0:
	   mapFilt[iter] = jet1;
	   break;
	 case 1:
	   mapFilt[iter] = jet2;
	   break;
	 case 2:
	   mapFilt[iter] = jet3;
	   break;
	 case 3:
	   mapFilt[iter] = jet4;
	   break;
	 case 4:
	   mapFilt[iter] = jet5;
	   break;
	 case 5:
	   mapFilt[iter] = jet6;
	   break;
	 case 6:
	   mapFilt[iter] = jet7;
	   break;
	 case 7:
	   mapFilt[iter] = jet8;
	   break;
	 case 8:
	   mapFilt[iter] = jet9;
	   break;
	 case 9:
	   mapFilt[iter] = jet10;
	   break;
	 default:
	   break;
	 }
	 
	 iter++;
       }
     }

     LV lepLV(vLeptonpt[0], vLeptoneta[0], vLeptonphi[0], 0.);
     LV MEtLV(met[0], 0., met[3] , 0.);

     ///////////////////////////////////////////////////////////////////////////

     //if(  !((abs(genTop.wdau1id)<6  && genTop.wdau1pt>40  && TMath::Abs(genTop.wdau1eta) <2.5 && genTop.wdau2pt>40  && TMath::Abs(genTop.wdau2eta) <2.5)  || 
     //   (abs(genTbar.wdau1id)<6 && genTbar.wdau1pt>40 && TMath::Abs(genTbar.wdau1eta)<2.5 && genTbar.wdau2pt>40 && TMath::Abs(genTbar.wdau2eta)<2.5))
     //  || !(genB.pt>30 && genBbar.pt>30)
     // ) continue;

     totCounter++;
     hCount->Fill(0.5);

     if( myJetsFilt.size()<6) continue;
     hCount->Fill(1.5);

     float a = lepLV.Px();
     float b = lepLV.Py();
     float c = lepLV.Pz();
     float f = lepLV.E();
     float x = MEtLV.Pt()*TMath::Cos(MEtLV.Phi());
     float y = MEtLV.Pt()*TMath::Sin(MEtLV.Phi());
     float m = 80.385;
     
     float delta = c*c*(-4*a*a - 4*b*b - 4*c*c + 4*f*f - 4*m*m - 8*a*x - 8*b*y)*(-4*a*a - 4*b*b - 4*c*c + 4*f*f - 4*m*m - 8*a*x - 8*b*y) - 4*(-4*c*c + 4*f*f)*(-a*a*a*a - 2*a*a*b*b - b*b*b*b - 2*a*a*c*c - 2*b*b*c*c - c*c*c*c + 2*a*a*f*f + 2*b*b*f*f + 2*c*c*f*f - f*f*f*f - 2*a*a*m*m - 2*b*b*m*m - 2*c*c*m*m + 2*f*f*m*m - m*m*m*m - 4*a*a*a*x - 4*a*b*b*x - 4*a*c*c*x + 4*a*f*f*x - 4*a*m*m*x - 4*a*a*x*x + 4*f*f*x*x - 4*a*a*b*y - 4*b*b*b*y - 4*b*c*c*y+ 4*b*f*f*y - 4*b*m*m*y - 8*a*b*x*y - 4*b*b*y*y + 4*f*f*y*y);
    
     //if(delta<0) continue;
     hCount->Fill(2.5);

     float z1 = delta>0 ? (-c*(-4*a*a - 4*b*b - 4*c*c + 4*f*f - 4*m*m - 8*a*x - 8*b*y) - TMath::Sqrt(delta)) / (2*(-4*c*c + 4*f*f)) : 0.0;
     float z2 = delta>0 ? (-c*(-4*a*a - 4*b*b - 4*c*c + 4*f*f - 4*m*m - 8*a*x - 8*b*y) + TMath::Sqrt(delta)) / (2*(-4*c*c + 4*f*f)) : 0.0;

     LV2 neutrinoLV2a(x,y,z1,TMath::Sqrt(x*x+y*y+z1*z1));
     LV neutrLV1(neutrinoLV2a);
     LV2 neutrinoLV2b(x,y,z2,TMath::Sqrt(x*x+y*y+z2*z2));
     LV neutrLV2(neutrinoLV2b);

     /////////////////////////////////////////////////////////////////////////// 

     std::map<float,  vector<unsigned int>, sorterByLess > posTagged;
     std::map<float,  int, sorterByLess > matchResult;
     vector<unsigned int> visited;

     for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
       float csv_k = mapFilt[k].csv>0 ? mapFilt[k].csv : 0.0 ;
       if( csv_k>0.80 ) continue;
       int genMatchk;
       findGenMatch(genMatchk, myJetsFilt[k] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);

       for(unsigned int l = 0; l < myJetsFilt.size(); l++){
	 if( l == k) continue;
	 float csv_l = mapFilt[l].csv>0 ? mapFilt[l].csv : 0.0 ;
	 if( csv_l>0.80 ) continue;

	 int genMatchl;
	 findGenMatch(genMatchl, myJetsFilt[l] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);
	 LV WCand     = myJetsFilt[k]+myJetsFilt[l];

	 for(unsigned int m = 0; m < myJetsFilt.size(); m++){
	   if( m==l || m==k) continue;
	   float csv_m = mapFilt[m].csv>0 ? mapFilt[m].csv : 0.0 ;
	   if( csv_m<0.30 ) continue;

	   LV bCand     = myJetsFilt[m];
	   int genMatchm;
	   findGenMatch(genMatchm, myJetsFilt[m] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);

	   for(unsigned int n = 0; n < myJetsFilt.size(); n++){
	     if( n==l || n==k || n==m) continue;
	     float csv_n = mapFilt[n].csv>0 ? mapFilt[n].csv : 0.0 ;
	     if( csv_n<0.30 ) continue;

	     LV bbarCand  = myJetsFilt[n];
	     int genMatchn;
	     findGenMatch(genMatchn, myJetsFilt[n] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);

	     if( !bookKeeper(k,l,m,n, visited) ) continue;

	     if((WCand).M()<40  || (WCand).M()>140) continue;
	     if((WCand+bCand).M()<100 || (WCand+bCand).M()>280) continue;
 
	     float likeB1     = fLikelihoodCsvB->Eval( csv_m );
	     float likeB2     = fLikelihoodCsvB->Eval( csv_n );
	     float likeW1     = fLikelihoodCsvNonB->Eval( csv_k );
	     float likeW2     = fLikelihoodCsvNonB->Eval( csv_l );
	     float likeTopHad = fLikelihoodWTopHadrMass->Eval( WCand.M(), (WCand+bCand).M());
	     float likeTopLep = delta>0 ? TMath::Max(fLikelihoodTopLeptMass->Eval( (lepLV+neutrLV1+bbarCand ).M()),
						    fLikelihoodTopLeptMass->Eval( (lepLV+neutrLV2+bbarCand ).M())  ) : 1.0;

	     float logLike  = likeW1*likeW2*likeB1*likeB2*likeTopHad*likeTopLep>0 ? -TMath::Log( likeW1*likeW2*likeB1*likeB2*likeTopHad*likeTopLep ) : 100 + ran->Uniform(0.,1);

	     matchResult[logLike] = int( (genMatchk==2 && genMatchl==2 && genMatchm==1 && genMatchn==-1) || (genMatchk==-2 && genMatchl==-2 && genMatchm==-1 && genMatchn==1) );
	     vector<unsigned int> res;
	     res.push_back(k);
	     res.push_back(l);
	     res.push_back(m);
	     res.push_back(n);
	     posTagged[logLike]   = res;

	     if( (genMatchk==2 && genMatchl==2 && genMatchm==1 && genMatchn==-1) || (genMatchk==-2 && genMatchl==-2 && genMatchm==-1 && genMatchn==1))
	       hLikeGood->Fill( logLike );
	     else
	       hLikeWrong->Fill( logLike );
	       
	   }
	 }
       }
     }

     hNumPermut->Fill(visited.size());
     
     int pos = 0;
     int firstValid1 = 0;
     int firstValid2 = 0;
     int firstValid3 = 0;
     unsigned int pos1 = 999; unsigned int pos2 = 999;
     unsigned int pos3 = 999; unsigned int pos4 = 999;
     int goodCands = 0; int goodCandsR=0;
     bool higgsFound = false;

     for(std::map<float, int>::iterator it = matchResult.begin(); it!= matchResult.end() ; it++){

       if( it->first > 15 ) continue;
       if(firstValid1==0){
	 hCount->Fill(3.5);
	 firstValid1 = 1;
       }

       vector<unsigned int> res = posTagged[ it->first ]; 
       unsigned int h1 = 999;
       unsigned int h2 = 999;

       for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
	 bool mask = false;
	 for( unsigned int l = 0; l < res.size() ; l++){
	   if( k == res[l]) mask = true;
	 }
	 float csv_k = mapFilt[k].csv>0 ? mapFilt[k].csv : 0.0 ;
	 if( h1 == 999 && csv_k< csvH1_ ) mask = true;
	 if( h1 != 999 && h2 == 999 && csv_k< csvH2_ ) mask = true;
	 float pt_k = mapFilt[k].pt  ;
	 if( h1 == 999 && pt_k < ptH1_ ) mask = true;
	 if( h1 != 999 && h2 == 999 && pt_k < ptH2_ ) mask = true;

	 if( !mask && h1 == 999)
	   h1 = k;
	 else if(!mask && h1 != 999 && h2==999)
	   h2 = k;
	 else{}
       }

       if(h1 != 999 && h2!=999 && pos<999){

	 if(firstValid2==0){
	   hCount->Fill(4.5);
	   firstValid2 = 1;
	 }


	 int genMatchh1;
	 findGenMatch(genMatchh1, myJetsFilt[h1] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);
	 int genMatchh2;
	 findGenMatch(genMatchh2, myJetsFilt[h2] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);

	 float PtoM =  (myJetsFilt[h1]+ myJetsFilt[h2]).Pt()/(myJetsFilt[h1]+ myJetsFilt[h2]).M();

//       if(pos1 == 999 && pos2 == 999){
// 	   pos1 = h1;  pos2 = h2;
// 	 }
// 	 else if((pos1 != h1 || pos2 != h2) && pos3 == 999 && pos4 == 999){
// 	   pos3 = h1;  pos4 = h2;
// 	 }
// 	 else{}

	 //if(pos1!=999 && pos2!=999 && (pos1 != h1 || pos2 != h2) && pos3 == 999 && pos4 == 999){
	 //pos3 = h1;  pos4 = h2;
	 //}
	 
	 if( PtoM>PtoM_ && (myJetsFilt[h1]+ myJetsFilt[h2]).Pt()>0 ){

	   if(firstValid3==0){
	     hCount->Fill(5.5);
	     firstValid3 = 1;
	   }
	 
	   goodCands++;
	   if( !higgsFound && genMatchh1==0 && genMatchh2==0 /*&& matchResult[it->first]*/){
	     hCombR->Fill(goodCands);
	     higgsFound=true;
	   }
	   if(pos1 == 999 && pos2 == 999){
	     pos1 = h1;  pos2 = h2;
	   }
	   else if(pos1!=999 && pos2!=999 && (pos1 != h1 || pos2 != h2) && pos3 == 999 && pos4 == 999){
	     pos3 = h1;  pos4 = h2;
	   }
	   else{}

	   //cout << "Event " << i << endl;
	   ///cout << "Permut #" << pos << ": likelihood = " << it->first << ", comb = [" << res[0] << res[1] << res[2] << res[3] << "]" <<endl;
	   //cout << "Higgs candidates at [" << h1 << h2 << "]" << " => pairing = [" << genMatchh1 << genMatchh2 << "]" << endl;

	   if( genMatchh1==0 && genMatchh2==0 && goodCands==1){
	     hMassJJR->Fill( (myJetsFilt[h1]+ myJetsFilt[h2]).M(), weight  );
	     hPtoMJJR->Fill( (myJetsFilt[h1]+ myJetsFilt[h2]).Pt()/(myJetsFilt[h1]+ myJetsFilt[h2]).M()  );
	     hPtJJR->Fill( (myJetsFilt[h1]+ myJetsFilt[h2]).Pt() );
	     hPt1JJR->Fill( (myJetsFilt[h1]).Pt() );
	     hPt2JJR->Fill( (myJetsFilt[h2]).Pt() );
	     h2PtoMJJR->Fill( (myJetsFilt[h1]+ myJetsFilt[h2]).Pt() , (myJetsFilt[h1]+ myJetsFilt[h2]).M());
	     h2PtPtR->Fill( (myJetsFilt[h1]).Pt() , (myJetsFilt[h2]).Pt() );
	     h2CsvR->Fill(  mapFilt[h1].csv,  mapFilt[h2].csv );

	     hMassWHadR->Fill( (myJetsFilt[res[0]] + myJetsFilt[res[1]]).M() );
	     hMassTopHadR->Fill( (myJetsFilt[res[0]] + myJetsFilt[res[1]] + myJetsFilt[res[2]]).M() );
	     hBCsvR->Fill(mapFilt[res[2]].csv);
	     hBbarCsvR->Fill(mapFilt[res[3]].csv);
	     hW1CsvR->Fill(mapFilt[res[0]].csv);
	     hW1CsvR->Fill(mapFilt[res[1]].csv);
	     hPt1oPt2R->Fill( (myJetsFilt[h2].Pt() / myJetsFilt[h1].Pt()) );
	   }
	   else if(!(genMatchh1==0 && genMatchh2==0 ) && goodCands==1){
	     hMassJJW->Fill( (myJetsFilt[h1]+ myJetsFilt[h2]).M() , weight );
	     hPtoMJJW->Fill( (myJetsFilt[h1]+ myJetsFilt[h2]).Pt()/(myJetsFilt[h1]+ myJetsFilt[h2]).M()  );
	     hPtJJW->Fill( (myJetsFilt[h1]+ myJetsFilt[h2]).Pt() );
	     hPt1JJW->Fill( (myJetsFilt[h1]).Pt() );
	     hPt2JJW->Fill( (myJetsFilt[h2]).Pt() );
	     h2PtoMJJW->Fill( (myJetsFilt[h1]+ myJetsFilt[h2]).Pt() , (myJetsFilt[h1]+ myJetsFilt[h2]).M());
	     h2PtPtW->Fill( (myJetsFilt[h1]).Pt() , (myJetsFilt[h2]).Pt() );
	     h2CsvW->Fill(  mapFilt[h1].csv,  mapFilt[h2].csv );

	     hMassWHadW->Fill( (myJetsFilt[res[0]] + myJetsFilt[res[1]]).M() );
	     hMassTopHadW->Fill( (myJetsFilt[res[0]] + myJetsFilt[res[1]] + myJetsFilt[res[2]]).M() );
	     hBCsvW->Fill(mapFilt[res[2]].csv);
	     hBbarCsvW->Fill(mapFilt[res[3]].csv);
	     hW1CsvW->Fill(mapFilt[res[0]].csv);
	     hW1CsvW->Fill(mapFilt[res[1]].csv);
	     hPt1oPt2W->Fill( (myJetsFilt[h2].Pt() / myJetsFilt[h1].Pt()) );
	     
	   }
	   else{}

	   //break;

	 }

       }
	 

       if( it->second == 1){
	 hGoodRank->Fill( pos );
       }
       pos++;
     }

     if(pos1!=999 && pos2 != 999 && pos3 != 999 && pos4 != 999 && goodCands>0){
       int genMatchp1;
       findGenMatch(genMatchp1, myJetsFilt[pos1] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);
       int genMatchp2;
       findGenMatch(genMatchp2, myJetsFilt[pos2] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);
       int genMatchp3;
       findGenMatch(genMatchp3, myJetsFilt[pos3] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);
       int genMatchp4;
       findGenMatch(genMatchp4, myJetsFilt[pos4] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);
  
       if( ((genMatchp1==0 && genMatchp2==0) && (genMatchp3!=0 || genMatchp4!=0)) ||
	   ((genMatchp1!=0 || genMatchp2!=0) && (genMatchp3==0 && genMatchp4==0))){
	 h2Mass12R->Fill( (myJetsFilt[pos1]+ myJetsFilt[pos2]).M(), (myJetsFilt[pos3]+ myJetsFilt[pos4]).M() );
       }
       else{
	 h2Mass12W->Fill( (myJetsFilt[pos1]+ myJetsFilt[pos2]).M(), (myJetsFilt[pos3]+ myJetsFilt[pos4]).M() );
       }

     }

     hComb->Fill(goodCands);


   }

   /////////////////////////////////////////////////////////////////////////// 
   
   
   fOutput->Write();
   fOutput->Close();
   
   fLikelihoods->Close();
   fInput->Close();
   
   delete ran; delete fOutput; delete fLikelihoods;
   
   return  ;
   
}






void algoBoost( float ptCut_ = 30 , float PtoM_ = 1.5, float csvH1_ = 0.679, float csvH2_ = 0.679, float ptH1_ = 30, float ptH2_ = 30){

  TRandom3* ran = new TRandom3();

  TFile *fOutput           = new TFile("algoBoost.root","RECREATE");
  TH1F* hCount         = new TH1F("hCount", "visited"   ,10, 0 ,10 );
  TH1F* hMassJJW           = new TH1F("hMassJJW", "Higgs mass"   ,40, 0 ,400 );
  TH1F* hMassJJR           = new TH1F("hMassJJR", "Higgs mass"   ,40, 0 ,400 );
  
  TH1F* hComb2Jet         = new TH1F("hComb2Jet", "visited"   ,10, 0 ,10 );
  TH1F* hComb2JetPt       = new TH1F("hComb2JetPt", "visited"   ,10, 0 ,10 );
  TH1F* hComb         = new TH1F("hComb", "visited"   ,10, 0 ,10 );

  TH1F* hMatch         = new TH1F("hMatch", "visited"   ,4, 0 ,4 );
  

  TFile* fLikelihoods          = new TFile("likelihoods_bkgd.root","READ");
  TF2* fLikelihoodWTopHadrMass = (TF2*)fLikelihoods->Get("fLikelihoodWTopHadrMass");
  TF1* fLikelihoodCsvB         = (TF1*)fLikelihoods->Get("fLikelihoodCsvB");
  TF1* fLikelihoodCsvNonB      = (TF1*)fLikelihoods->Get("fLikelihoodCsvNonB");
  TF1* fLikelihoodTopLeptMass  = (TF1*)fLikelihoods->Get("fLikelihoodTopLeptMass");


  TFile* fInput = TFile::Open("gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_VType2/v2/DiJetPt_TTH_HToBB_M-125_8TeV-pythia6_VType2.root","READ");
  //TFile* fInput = TFile::Open("gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_VType2/v2/DiJetPt_TTJets_SemiLeptMGDecays_8TeV-madgraph-part_VType2.root","READ");
  TTree* tree = (TTree*)fInput->Get("tree");

  if(!tree){
    cout << "No tree" << endl; return ;
  }

  TH1F* hCounter = (TH1F*)fInput->Get("Count");
  float weight = (0.1302*0.569)*12.1*1000/hCounter->GetBinContent(1);
  //float weight = (103.0)*12.1*1000/hCounter->GetBinContent(1);
  cout << weight << endl;
  //weight /= (500000./tree->GetEntries());
  cout << weight << endl;

  JetByPt jet1, jet2, jet3, jet4, jet5, jet6, jet7, jet8, jet9, jet10;
  genTopInfo genTop, genTbar;
  genParticleInfo genB, genBbar;
  int nvlep;
  float vLeptonpt[999]; float vLeptoneta[999]; float vLeptonphi[999];
  float met[999];

  tree->SetBranchAddress("jet1",   &jet1);
  tree->SetBranchAddress("jet2",   &jet2);
  tree->SetBranchAddress("jet3",   &jet3);
  tree->SetBranchAddress("jet4",   &jet4);
  tree->SetBranchAddress("jet5",   &jet5);
  tree->SetBranchAddress("jet6",   &jet6);
  tree->SetBranchAddress("jet7",   &jet7);
  tree->SetBranchAddress("jet8",   &jet8);
  tree->SetBranchAddress("jet9",   &jet9);
  tree->SetBranchAddress("jet10",  &jet10);
  tree->SetBranchAddress("genTop", &genTop);
  tree->SetBranchAddress("genTbar",&genTbar);
  tree->SetBranchAddress("genB",   &genB);
  tree->SetBranchAddress("genBbar",&genBbar);
  tree->SetBranchAddress("vLepton_pt", vLeptonpt);
  tree->SetBranchAddress("vLepton_phi",vLeptonphi);
  tree->SetBranchAddress("vLepton_eta",vLeptoneta);
  tree->SetBranchAddress("nvlep",  &nvlep);
  tree->SetBranchAddress("METtype1p2corr",  met);

  Long64_t nentries = tree->GetEntries();
  int totCounter = 0;
  
  for (Long64_t i = 0; i < /*500000*/nentries; i++){
    
    if(i%10000==0) cout << i << endl;
    tree->GetEntry(i);
    
    //continue;

    if(nvlep<=0) continue;


     LV topBLV(  0.,0.,0.,0.);
     LV topW1LV( 0.,0.,0.,0.); 
     LV topW2LV( 0.,0.,0.,0.);
     LV atopBLV( 0.,0.,0.,0.); 
     LV atopW1LV(0.,0.,0.,0.); 
     LV atopW2LV(0.,0.,0.,0.);
     LV genBLV(0.,0.,0.,0.);
     LV genBbarLV(0.,0.,0.,0.);
    
    
     ///////////////////////////////////////////////////////////////////
     /////////////////////////////////////////////////////////////////// 
     if(genTop.bmass>0){
       topBLV.SetPt(   genTop.bpt );
       topBLV.SetEta(  genTop.beta );
       topBLV.SetPhi(  genTop.bphi );
       topBLV.SetM(    genTop.bmass );
       topW1LV.SetPt(  genTop.wdau1pt );
       topW1LV.SetEta( genTop.wdau1eta );
       topW1LV.SetPhi( genTop.wdau1phi );
       topW1LV.SetM(   genTop.wdau1mass );
       topW2LV.SetPt(  genTop.wdau2pt );
       topW2LV.SetEta( genTop.wdau2eta );
       topW2LV.SetPhi( genTop.wdau2phi );
       topW2LV.SetM(   genTop.wdau2mass );
     }
     if(genTbar.bmass>0){
	atopBLV.SetPt(   genTbar.bpt );
	atopBLV.SetEta(  genTbar.beta );
	atopBLV.SetPhi(  genTbar.bphi );
	atopBLV.SetM(    genTbar.bmass );
	atopW1LV.SetPt(  genTbar.wdau1pt );
	atopW1LV.SetEta( genTbar.wdau1eta );
	atopW1LV.SetPhi( genTbar.wdau1phi );
	atopW1LV.SetM(   genTbar.wdau1mass );
	atopW2LV.SetPt(  genTbar.wdau2pt );
	atopW2LV.SetEta( genTbar.wdau2eta );
	atopW2LV.SetPhi( genTbar.wdau2phi );
	atopW2LV.SetM(   genTbar.wdau2mass );
     }
     if(genB.mass>0 && genB.momid==25){
       genBLV.SetPt(  genB.pt );
       genBLV.SetEta( genB.eta );
       genBLV.SetPhi( genB.phi );
       genBLV.SetM(   genB.mass );
     }
     if(genBbar.mass>0 && genBbar.momid==25){
       genBbarLV.SetPt(  genBbar.pt );
       genBbarLV.SetEta( genBbar.eta );
       genBbarLV.SetPhi( genBbar.phi );
       genBbarLV.SetM(   genBbar.mass );
     }
     ///////////////////////////////////////////////////////////////////
     ///////////////////////////////////////////////////////////////////

     vector<LV> myJets;
     std::map<unsigned int, JetByPt> map;

     map[0] = jet1;  map[1] = jet2; map[2] = jet3; map[3] = jet4; 
     map[4] = jet5;  map[5] = jet6; map[6] = jet7; map[7] = jet8; 
     map[8] = jet9;  map[9] = jet10; 

     myJets.push_back(LV(jet1.pt, jet1.eta, jet1.phi, jet1.mass));
     myJets.push_back(LV(jet2.pt, jet2.eta, jet2.phi, jet2.mass));
     myJets.push_back(LV(jet3.pt, jet3.eta, jet3.phi, jet3.mass));
     myJets.push_back(LV(jet4.pt, jet4.eta, jet4.phi, jet4.mass));
     myJets.push_back(LV(jet5.pt, jet5.eta, jet5.phi, jet5.mass));
     myJets.push_back(LV(jet6.pt, jet6.eta, jet6.phi, jet6.mass));
     myJets.push_back(LV(jet7.pt, jet7.eta, jet7.phi, jet7.mass));
     myJets.push_back(LV(jet8.pt, jet8.eta, jet8.phi, jet8.mass));
     myJets.push_back(LV(jet9.pt, jet9.eta, jet9.phi, jet9.mass));
     myJets.push_back(LV(jet10.pt,jet10.eta,jet10.phi,jet10.mass));

     
     vector<LV> myJetsFilt;
     std::map<unsigned int, JetByPt> mapFilt;

     unsigned int iter = 0;
     for(unsigned int k = 0; k<myJets.size() ; k++){

       if(myJets[k].Pt()>ptCut_ && TMath::Abs(myJets[k].Eta()) < 2.5 ){

	 myJetsFilt.push_back( myJets[k] );

	 switch(k){
	 case 0:
	   mapFilt[iter] = jet1;
	   break;
	 case 1:
	   mapFilt[iter] = jet2;
	   break;
	 case 2:
	   mapFilt[iter] = jet3;
	   break;
	 case 3:
	   mapFilt[iter] = jet4;
	   break;
	 case 4:
	   mapFilt[iter] = jet5;
	   break;
	 case 5:
	   mapFilt[iter] = jet6;
	   break;
	 case 6:
	   mapFilt[iter] = jet7;
	   break;
	 case 7:
	   mapFilt[iter] = jet8;
	   break;
	 case 8:
	   mapFilt[iter] = jet9;
	   break;
	 case 9:
	   mapFilt[iter] = jet10;
	   break;
	 default:
	   break;
	 }
	 
	 iter++;
       }
     }

     LV lepLV(vLeptonpt[0], vLeptoneta[0], vLeptonphi[0], 0.);
     LV MEtLV(met[0], 0., met[3] , 0.);

     ///////////////////////////////////////////////////////////////////////////

     //if(  !((abs(genTop.wdau1id)<6  && genTop.wdau1pt>40  && TMath::Abs(genTop.wdau1eta) <2.5 && genTop.wdau2pt>40  && TMath::Abs(genTop.wdau2eta) <2.5)  || 
     //   (abs(genTbar.wdau1id)<6 && genTbar.wdau1pt>40 && TMath::Abs(genTbar.wdau1eta)<2.5 && genTbar.wdau2pt>40 && TMath::Abs(genTbar.wdau2eta)<2.5))
     //  || !(genB.pt>30 && genBbar.pt>30)
     // ) continue;

     totCounter++;
     hCount->Fill(0.5);

     if( myJetsFilt.size()<4) continue;
     hCount->Fill(1.5);

     float a = lepLV.Px();
     float b = lepLV.Py();
     float c = lepLV.Pz();
     float f = lepLV.E();
     float x = MEtLV.Pt()*TMath::Cos(MEtLV.Phi());
     float y = MEtLV.Pt()*TMath::Sin(MEtLV.Phi());
     float m = 80.385;
     
     float delta = c*c*(-4*a*a - 4*b*b - 4*c*c + 4*f*f - 4*m*m - 8*a*x - 8*b*y)*(-4*a*a - 4*b*b - 4*c*c + 4*f*f - 4*m*m - 8*a*x - 8*b*y) - 4*(-4*c*c + 4*f*f)*(-a*a*a*a - 2*a*a*b*b - b*b*b*b - 2*a*a*c*c - 2*b*b*c*c - c*c*c*c + 2*a*a*f*f + 2*b*b*f*f + 2*c*c*f*f - f*f*f*f - 2*a*a*m*m - 2*b*b*m*m - 2*c*c*m*m + 2*f*f*m*m - m*m*m*m - 4*a*a*a*x - 4*a*b*b*x - 4*a*c*c*x + 4*a*f*f*x - 4*a*m*m*x - 4*a*a*x*x + 4*f*f*x*x - 4*a*a*b*y - 4*b*b*b*y - 4*b*c*c*y+ 4*b*f*f*y - 4*b*m*m*y - 8*a*b*x*y - 4*b*b*y*y + 4*f*f*y*y);
    
     //if(delta<0) continue;
     hCount->Fill(2.5);

     float z1 = delta>0 ? (-c*(-4*a*a - 4*b*b - 4*c*c + 4*f*f - 4*m*m - 8*a*x - 8*b*y) - TMath::Sqrt(delta)) / (2*(-4*c*c + 4*f*f)) : 0.0;
     float z2 = delta>0 ? (-c*(-4*a*a - 4*b*b - 4*c*c + 4*f*f - 4*m*m - 8*a*x - 8*b*y) + TMath::Sqrt(delta)) / (2*(-4*c*c + 4*f*f)) : 0.0;

     LV2 neutrinoLV2a(x,y,z1,TMath::Sqrt(x*x+y*y+z1*z1));
     LV neutrLV1(neutrinoLV2a);
     LV2 neutrinoLV2b(x,y,z2,TMath::Sqrt(x*x+y*y+z2*z2));
     LV neutrLV2(neutrinoLV2b);

     /////////////////////////////////////////////////////////////////////////// 

     std::map<float,  vector<unsigned int>, sorterByLess > posTagged;
     int firstValid1 = 0;
     int firstValid2 = 0;
     int firstValid3 = 0;
     int firstValid4 = 0;
     int num2Jets = 0;
     int num2JetsPt = 0;

     int match1=0; int match2=0;

     for(unsigned int k = 0; k < myJetsFilt.size()-1; k++){  
       float csv_k = mapFilt[k].csv>0 ? mapFilt[k].csv : 0.0 ;
       float pt_k  = mapFilt[k].pt>0 ? mapFilt[k].pt : 0.0 ;
       if( csv_k<0.5 || pt_k<20) continue;
       int genMatchk;
       findGenMatch(genMatchk, myJetsFilt[k] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);

       if(genMatchk==0) match1++;

       for(unsigned int l = k+1; l < myJetsFilt.size(); l++){
	 if( l == k) continue;
	 float csv_l = mapFilt[l].csv>0 ? mapFilt[l].csv : 0.0 ;
	 float pt_l  = mapFilt[l].pt>0 ? mapFilt[l].pt : 0.0 ;
	 if( csv_l<0.5 || pt_l<20) continue;
	 int genMatchl;
	 findGenMatch(genMatchl, myJetsFilt[l] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);

	 if(firstValid1==0){
	   hCount->Fill(3.5);
	   firstValid1 = 1;
	 }

	 if(genMatchl==0) match2++;

	 num2Jets++;



	 LV HCand     = myJetsFilt[k]+myJetsFilt[l];
	 if( HCand.Pt() < 150 ) continue;

	 num2JetsPt++;

	 if(firstValid2==0){
	   hCount->Fill(4.5);
	   firstValid2 = 1;
	 }

	 for(unsigned int m = 0; m < myJetsFilt.size(); m++){
	   if( m==l || m==k) continue;
	   float csv_m = mapFilt[m].csv>0 ? mapFilt[m].csv : 0.0 ;
	   if( csv_m>999 ) continue;
	   for(unsigned int n = m+1; n < myJetsFilt.size(); n++){
	     if( n==l || n==k || n==m) continue;
	     float csv_n = mapFilt[n].csv>0 ? mapFilt[n].csv : 0.0 ;
	     if( csv_n>999 ) continue;
	     for(unsigned int o = 0; o < myJetsFilt.size(); o++){
	       if( o==l || o==k || o==m || o==n) continue;
	       float csv_o = mapFilt[o].csv>0 ? mapFilt[o].csv : 0.0 ;
	       if( csv_o<0.5 ) continue;

	       LV WHadCand     = myJetsFilt[m]+myJetsFilt[n];
	       LV THadCand      = WHadCand+myJetsFilt[o];

	       if((WHadCand).M()<65  || (WHadCand).M()>95) continue;
	       if((THadCand).M()<150 || (THadCand).M()>200) continue;

	       if(firstValid3==0){
		 hCount->Fill(5.5);
		 firstValid3 = 1;
	       }
	       
	       float THadPt = THadCand.Pt();
	       if( THadPt>200 && TMath::Cos( THadCand.Phi() - HCand.Phi())<999){

		 bool thirdBTag = false;

		 for(unsigned int p = 0; p < myJetsFilt.size(); p++){
		   if( p==l || p==k || p==m || p==n || p==o) continue;
		   float csv_p = mapFilt[p].csv>0 ? mapFilt[p].csv : 0.0 ;
		   if( csv_p>0.5 ) thirdBTag = true;
		 }


		 if(firstValid4==0){
		   hCount->Fill(6.5);
		   firstValid4 = 1;
		 }

		 vector<unsigned int> res;
		 res.push_back(k);
		 res.push_back(l);
		 res.push_back(m);
		 res.push_back(n);
		 res.push_back(o);
		 if(thirdBTag)
		   posTagged[ TMath::Abs(WHadCand.M()-81)+TMath::Abs(THadCand.M()-175) ] =  res;

	       }
	     }	     
	   }	   
	 }
       }
     }
     hComb->Fill(posTagged.size());
     hComb2Jet->Fill(num2Jets);
     hComb2JetPt->Fill(num2JetsPt);
     
     int pos = 0;
     int goodCands = 0;

     //if(posTagged.size()>0) cout << "Event " << i << endl;
     for(std::map<float, vector<unsigned int> >::iterator it = posTagged.begin(); it!= posTagged.end() ; it++){
       vector<unsigned int> res = it->second; 
       unsigned int h1 = res[0];
       unsigned int h2 = res[1];
    

       int genMatchh1;
       findGenMatch(genMatchh1, myJetsFilt[h1] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);
       int genMatchh2;
       findGenMatch(genMatchh2, myJetsFilt[h2] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);
	 
       //cout << "Permut #" << pos << ": Pt/M = " << it->first << ", comb = [" << res[0] << res[1] << res[2] << res[3] << "]" <<endl;
       ///cout << "Higgs candidates at [" << h1 << h2 << "]" << " => pairing = [" << genMatchh1 << genMatchh2 << "]" << endl;
   

       if((myJetsFilt[h1]+ myJetsFilt[h2]).Pt() / (myJetsFilt[h1]+ myJetsFilt[h2]).M() > -999 && goodCands==0 ){

	 if( genMatchh1==0 && genMatchh2==0 )
	   hMassJJR->Fill( (myJetsFilt[h1]+ myJetsFilt[h2]).M() , weight);
	 else
	   hMassJJW->Fill( (myJetsFilt[h1]+ myJetsFilt[h2]).M() , weight);
	 goodCands++;
       }
       pos++;
     }

     hMatch->Fill(match1+match2);
    

  }


   /////////////////////////////////////////////////////////////////////////// 
   
   
   fOutput->Write();
   fOutput->Close();
   
   fLikelihoods->Close();
   fInput->Close();
   
   delete ran; delete fOutput; delete fLikelihoods;
   
   return  ;
   
}



void algoIncl( float ptCut_ = 30 , float PtoM_ = 1.5, float csvH1_ = 0.679, float csvH2_ = 0.679, float ptH1_ = 30, float ptH2_ = 30){

  TRandom3* ran = new TRandom3();

  TFile *fOutput           = new TFile("algoIncl.root","RECREATE");
  TH1F* hCount             = new TH1F("hCount", "visited"   ,10, 0 ,10 );
  TH1F* hMassJJW           = new TH1F("hMassJJW", "Higgs mass"   ,40, 0 ,400 );
  TH1F* hMassJJR           = new TH1F("hMassJJR", "Higgs mass"   ,40, 0 ,400 );

  TH1F* hMassJJOneEventR     = new TH1F("hMassJJOneEventR", "Higgs mass"   ,400, 0 ,400 );

  TH1F* hMassCat1R           = new TH1F("hMassCat1R", "Higgs mass"   ,40, 0 ,400 );
  TH1F* hMassCat2R           = new TH1F("hMassCat2R", "Higgs mass"   ,40, 0 ,400 );
  TH1F* hMassCat3R           = new TH1F("hMassCat3R", "Higgs mass"   ,40, 0 ,400 );
  TH1F* hMassCat4R           = new TH1F("hMassCat4R", "Higgs mass"   ,40, 0 ,400 );
  TH1F* hMassCat5R           = new TH1F("hMassCat5R", "Higgs mass"   ,40, 0 ,400 );
  TH1F* hMassCat6R           = new TH1F("hMassCat6R", "Higgs mass"   ,40, 0 ,400 );

  TH2F* hMassCat3R2d           = new TH2F("hMassCat3R2d", "Higgs mass"   ,40, 0 ,400 , 10, 0 ,10);
  TH2F* hMassCat3W2d           = new TH2F("hMassCat3W2d", "Higgs mass"   ,40, 0 ,400 , 10, 0 ,10);

  TH1F* hMassCat1W           = new TH1F("hMassCat1W", "Higgs mass"   ,40, 0 ,400 );
  TH1F* hMassCat2W           = new TH1F("hMassCat2W", "Higgs mass"   ,40, 0 ,400 );
  TH1F* hMassCat3W           = new TH1F("hMassCat3W", "Higgs mass"   ,40, 0 ,400 );
  TH1F* hMassCat4W           = new TH1F("hMassCat4W", "Higgs mass"   ,40, 0 ,400 );
  TH1F* hMassCat5W           = new TH1F("hMassCat5W", "Higgs mass"   ,40, 0 ,400 );
  TH1F* hMassCat6W           = new TH1F("hMassCat6W", "Higgs mass"   ,40, 0 ,400 );
  
  TH1F* hComb2Jet         = new TH1F("hComb2Jet", "visited"   ,10, 0 ,10 );
  TH1F* hComb2JetPt       = new TH1F("hComb2JetPt", "visited"   ,10, 0 ,10 );
  TH1F* hComb             = new TH1F("hComb", "visited"   ,30, 0 ,30 );

  TH1F* hMatch         = new TH1F("hMatch", "visited"   ,10, 0 ,10 );
  TH2F* hMatchDebug         = new TH2F("hMatchDebug", "visited"   ,12, -6 ,6  ,12, -6 ,6 );
  TH1F* hMatchAny         = new TH1F("hMatchAny", "visited"   ,4, 0 ,4 );
  TH1F* hMatchAll         = new TH1F("hMatchAll", "visited"   ,4, 0 ,4 );
  TH1F* hMatchAllPt         = new TH1F("hMatchAllPt", "visited"   ,4, 0 ,4 );
  TH1F* hMatchAllPtBTag         = new TH1F("hMatchAllPtBTag", "visited"   ,4, 0 ,4 );

  TH1F* hPtResHiggs         = new TH1F("hPtResHiggs", "visited"   ,100, -10 ,10 );
  TH1F* hPtResHiggsAmb         = new TH1F("hPtResHiggsAmb", "visited"   ,100, -10 ,10 );
  
  

  TH1F* hDPhiR         = new TH1F("hDPhiR", "visited"   ,20, -4 ,4 );
  TH1F* hDPhiW         = new TH1F("hDPhiW", "visited"   ,20, -4 ,4 );

  TH1F* hDRR         = new TH1F("hDRR", "visited"   ,40, 0 ,8 );
  TH1F* hDRW         = new TH1F("hDRW", "visited"   ,40, 0 ,8 );

  TH2F* hMatchComb         = new TH2F("hMatchComb", "visited"   ,10, 0 ,10, 10, 0 ,10 );
  TH1F* hCombPass6         = new TH1F("hCombPass6", "visited"   ,10, 0 ,10 );
  TH1F* hCombPass3         = new TH1F("hCombPass3", "visited"   ,10, 0 ,10 );
  TH1F* hCombPass1         = new TH1F("hCombPass1", "visited"   ,10, 0 ,10 );
  TH1F* hPtMerg         = new TH1F("hPtMerg", "visited"   ,40, 0. ,400 );
  TH1F* hMassMerg         = new TH1F("hMassMerg", "visited"   ,100, 0. ,1 );
  TH1F* hMassUnMerg         = new TH1F("hMassUnMerg", "visited"   ,100, 0. ,1 );

  TH1F* hFrac6         = new TH1F("hFrac6", "visited"   ,6, 0. ,6 );
  TH1F* hFrac3         = new TH1F("hFrac3", "visited"   ,6, 0. ,6 );
  TH1F* hFrac1         = new TH1F("hFrac1", "visited"   ,6, 0. ,6 );

  TH1F* hPt1R         = new TH1F("hPt1R", "visited"   ,40, 0. ,200 );
  TH1F* hPt2R         = new TH1F("hPt2R", "visited"   ,40, 0. ,200 );
  TH1F* hPtHR         = new TH1F("hPtHR", "visited"   ,40, 0. ,200 );
  TH1F* hPtSumR       = new TH1F("hPtSumR", "visited"   ,40, 0. ,200 );
  TH1F* hBTag1R       = new TH1F("hBTag1R", "visited"   ,40, 0. ,1.1 );
  TH1F* hBTag2R       = new TH1F("hBTag2R", "visited"   ,40, 0. ,1.1 );
  TH1F* hMass1R       = new TH1F("hMass1R", "visited"   ,40, 0. ,80 );
  TH1F* hMass2R       = new TH1F("hMass2R", "visited"   ,40, 0. ,80 );

  TH1F* hPt1W         = new TH1F("hPt1W", "visited"   ,40, 0. ,200 );
  TH1F* hPt2W         = new TH1F("hPt2W", "visited"   ,40, 0. ,200 );
  TH1F* hPtHW         = new TH1F("hPtHW", "visited"   ,40, 0. ,200 );
  TH1F* hPtSumW       = new TH1F("hPtSumW", "visited"   ,40, 0. ,200 );
  TH1F* hBTag1W       = new TH1F("hBTag1W", "visited"   ,40, 0. ,1.1 );
  TH1F* hBTag2W       = new TH1F("hBTag2W", "visited"   ,40, 0. ,1.1 );
  TH1F* hMass1W       = new TH1F("hMass1W", "visited"   ,40, 0. ,80 );
  TH1F* hMass2W       = new TH1F("hMass2W", "visited"   ,40, 0. ,80 );

  TH1F* hMassCat3It1R           = new TH1F("hMassCat3It1R", "Higgs mass"   ,40, 0 ,400 );
  TH1F* hMassCat3It2R           = new TH1F("hMassCat3It2R", "Higgs mass"   ,40, 0 ,400 );
  TH1F* hMassCat3It3R           = new TH1F("hMassCat3It3R", "Higgs mass"   ,40, 0 ,400 );
  TH1F* hMassCat3It4R           = new TH1F("hMassCat3It4R", "Higgs mass"   ,40, 0 ,400 );
  TH1F* hMassCat3It5R           = new TH1F("hMassCat3It5R", "Higgs mass"   ,40, 0 ,400 );
  TH1F* hMassCat3It6R           = new TH1F("hMassCat3It6R", "Higgs mass"   ,40, 0 ,400 );

  TH1F* hMassCat3It1W           = new TH1F("hMassCat3It1W", "Higgs mass"   ,40, 0 ,400 );
  TH1F* hMassCat3It2W           = new TH1F("hMassCat3It2W", "Higgs mass"   ,40, 0 ,400 );
  TH1F* hMassCat3It3W           = new TH1F("hMassCat3It3W", "Higgs mass"   ,40, 0 ,400 );
  TH1F* hMassCat3It4W           = new TH1F("hMassCat3It4W", "Higgs mass"   ,40, 0 ,400 );
  TH1F* hMassCat3It5W           = new TH1F("hMassCat3It5W", "Higgs mass"   ,40, 0 ,400 );
  TH1F* hMassCat3It6W           = new TH1F("hMassCat3It6W", "Higgs mass"   ,40, 0 ,400 );


  TFile* fLikelihoods          = new TFile("likelihoods_bkgd.root","READ");
  TF2* fLikelihoodWTopHadrMass = (TF2*)fLikelihoods->Get("fLikelihoodWTopHadrMass");
  TF1* fLikelihoodCsvB         = (TF1*)fLikelihoods->Get("fLikelihoodCsvB");
  TF1* fLikelihoodCsvNonB      = (TF1*)fLikelihoods->Get("fLikelihoodCsvNonB");
  TF1* fLikelihoodTopLeptMass  = (TF1*)fLikelihoods->Get("fLikelihoodTopLeptMass");


  TFile* fInput = TFile::Open("gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_VType2/v2/DiJetPt_TTH_HToBB_M-125_8TeV-pythia6_VType2.root","READ");
  //TFile* fInput = TFile::Open("gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_VType2/v2/DiJetPt_TTJets_SemiLeptMGDecays_8TeV-madgraph-part_VType2.root","READ");
  TTree* tree = (TTree*)fInput->Get("tree");

  if(!tree){
    cout << "No tree" << endl; return ;
  }

  TH1F* hCounter = (TH1F*)fInput->Get("Count");
  float weight = (0.1302*0.569)*12.1*1000/hCounter->GetBinContent(1);
  //float weight = (103.0)*12.1*1000/hCounter->GetBinContent(1);
  cout << weight << endl;
  //weight /= (500000./tree->GetEntries());
  cout << weight << endl;

  JetByPt jet1, jet2, jet3, jet4, jet5, jet6, jet7, jet8, jet9, jet10;
  genTopInfo genTop, genTbar;
  genParticleInfo genB, genBbar;
  int nvlep;
  float vLeptonpt[999]; float vLeptoneta[999]; float vLeptonphi[999];
  float met[999];

  tree->SetBranchAddress("jet1",   &jet1);
  tree->SetBranchAddress("jet2",   &jet2);
  tree->SetBranchAddress("jet3",   &jet3);
  tree->SetBranchAddress("jet4",   &jet4);
  tree->SetBranchAddress("jet5",   &jet5);
  tree->SetBranchAddress("jet6",   &jet6);
  tree->SetBranchAddress("jet7",   &jet7);
  tree->SetBranchAddress("jet8",   &jet8);
  tree->SetBranchAddress("jet9",   &jet9);
  tree->SetBranchAddress("jet10",  &jet10);
  tree->SetBranchAddress("genTop", &genTop);
  tree->SetBranchAddress("genTbar",&genTbar);
  tree->SetBranchAddress("genB",   &genB);
  tree->SetBranchAddress("genBbar",&genBbar);
  tree->SetBranchAddress("vLepton_pt", vLeptonpt);
  tree->SetBranchAddress("vLepton_phi",vLeptonphi);
  tree->SetBranchAddress("vLepton_eta",vLeptoneta);
  tree->SetBranchAddress("nvlep",  &nvlep);
  tree->SetBranchAddress("METtype1p2corr",  met);

  Long64_t nentries = tree->GetEntries();
  int totCounter = 0;
  int firstOcc = 0;
  
  for (Long64_t i = 0; i < /*500000*/ nentries; i++){
    
    if(i%10000==0) cout << i << endl;
    tree->GetEntry(i);
    
    //continue;

    if(nvlep<=0) continue;


     LV topBLV(  0.,0.,0.,0.);
     LV topW1LV( 0.,0.,0.,0.); 
     LV topW2LV( 0.,0.,0.,0.);
     LV atopBLV( 0.,0.,0.,0.); 
     LV atopW1LV(0.,0.,0.,0.); 
     LV atopW2LV(0.,0.,0.,0.);
     LV genBLV(0.,0.,0.,0.);
     LV genBbarLV(0.,0.,0.,0.);
    
    
     ///////////////////////////////////////////////////////////////////
     /////////////////////////////////////////////////////////////////// 
     if(genTop.bmass>0){
       topBLV.SetPt(   genTop.bpt );
       topBLV.SetEta(  genTop.beta );
       topBLV.SetPhi(  genTop.bphi );
       topBLV.SetM(    genTop.bmass );
       topW1LV.SetPt(  genTop.wdau1pt );
       topW1LV.SetEta( genTop.wdau1eta );
       topW1LV.SetPhi( genTop.wdau1phi );
       topW1LV.SetM(   genTop.wdau1mass );
       topW2LV.SetPt(  genTop.wdau2pt );
       topW2LV.SetEta( genTop.wdau2eta );
       topW2LV.SetPhi( genTop.wdau2phi );
       topW2LV.SetM(   genTop.wdau2mass );
     }
     if(genTbar.bmass>0){
	atopBLV.SetPt(   genTbar.bpt );
	atopBLV.SetEta(  genTbar.beta );
	atopBLV.SetPhi(  genTbar.bphi );
	atopBLV.SetM(    genTbar.bmass );
	atopW1LV.SetPt(  genTbar.wdau1pt );
	atopW1LV.SetEta( genTbar.wdau1eta );
	atopW1LV.SetPhi( genTbar.wdau1phi );
	atopW1LV.SetM(   genTbar.wdau1mass );
	atopW2LV.SetPt(  genTbar.wdau2pt );
	atopW2LV.SetEta( genTbar.wdau2eta );
	atopW2LV.SetPhi( genTbar.wdau2phi );
	atopW2LV.SetM(   genTbar.wdau2mass );
     }
     if(genB.mass>0 && genB.momid==25){
       genBLV.SetPt(  genB.pt );
       genBLV.SetEta( genB.eta );
       genBLV.SetPhi( genB.phi );
       genBLV.SetM(   genB.mass );
     }
     if(genBbar.mass>0 && genBbar.momid==25){
       genBbarLV.SetPt(  genBbar.pt );
       genBbarLV.SetEta( genBbar.eta );
       genBbarLV.SetPhi( genBbar.phi );
       genBbarLV.SetM(   genBbar.mass );
     }
     ///////////////////////////////////////////////////////////////////
     ///////////////////////////////////////////////////////////////////

     vector<LV> myJets;
     std::map<unsigned int, JetByPt> map;

     map[0] = jet1;  map[1] = jet2; map[2] = jet3; map[3] = jet4; 
     map[4] = jet5;  map[5] = jet6; map[6] = jet7; map[7] = jet8; 
     map[8] = jet9;  map[9] = jet10; 

     myJets.push_back(LV(jet1.pt, jet1.eta, jet1.phi, jet1.mass));
     myJets.push_back(LV(jet2.pt, jet2.eta, jet2.phi, jet2.mass));
     myJets.push_back(LV(jet3.pt, jet3.eta, jet3.phi, jet3.mass));
     myJets.push_back(LV(jet4.pt, jet4.eta, jet4.phi, jet4.mass));
     myJets.push_back(LV(jet5.pt, jet5.eta, jet5.phi, jet5.mass));
     myJets.push_back(LV(jet6.pt, jet6.eta, jet6.phi, jet6.mass));
     myJets.push_back(LV(jet7.pt, jet7.eta, jet7.phi, jet7.mass));
     myJets.push_back(LV(jet8.pt, jet8.eta, jet8.phi, jet8.mass));
     myJets.push_back(LV(jet9.pt, jet9.eta, jet9.phi, jet9.mass));
     myJets.push_back(LV(jet10.pt,jet10.eta,jet10.phi,jet10.mass));

     
     vector<LV> myJetsFilt;
     std::map<unsigned int, JetByPt> mapFilt;

     unsigned int iter = 0;
     for(unsigned int k = 0; k<myJets.size() ; k++){

       if(myJets[k].Pt()>ptCut_ && TMath::Abs(myJets[k].Eta()) < 2.5 ){

	 myJetsFilt.push_back( myJets[k] );

	 switch(k){
	 case 0:
	   mapFilt[iter] = jet1;
	   break;
	 case 1:
	   mapFilt[iter] = jet2;
	   break;
	 case 2:
	   mapFilt[iter] = jet3;
	   break;
	 case 3:
	   mapFilt[iter] = jet4;
	   break;
	 case 4:
	   mapFilt[iter] = jet5;
	   break;
	 case 5:
	   mapFilt[iter] = jet6;
	   break;
	 case 6:
	   mapFilt[iter] = jet7;
	   break;
	 case 7:
	   mapFilt[iter] = jet8;
	   break;
	 case 8:
	   mapFilt[iter] = jet9;
	   break;
	 case 9:
	   mapFilt[iter] = jet10;
	   break;
	 default:
	   break;
	 }
	 
	 iter++;
       }
     }

     LV lepLV(vLeptonpt[0], vLeptoneta[0], vLeptonphi[0], 0.);
     LV MEtLV(met[0], 0., met[3] , 0.);

     ///////////////////////////////////////////////////////////////////////////

     //if(  !((abs(genTop.wdau1id)<6  && genTop.wdau1pt>40  && TMath::Abs(genTop.wdau1eta) <2.5 && genTop.wdau2pt>40  && TMath::Abs(genTop.wdau2eta) <2.5)  || 
     //   (abs(genTbar.wdau1id)<6 && genTbar.wdau1pt>40 && TMath::Abs(genTbar.wdau1eta)<2.5 && genTbar.wdau2pt>40 && TMath::Abs(genTbar.wdau2eta)<2.5))
     //  || !(genB.pt>30 && genBbar.pt>30)
     // ) continue;

     totCounter++;
     hCount->Fill(0.5,  weight);

     if( myJetsFilt.size()<4) continue;
     hCount->Fill(1.5,  weight);

     float a = lepLV.Px();
     float b = lepLV.Py();
     float c = lepLV.Pz();
     float f = lepLV.E();
     float x = MEtLV.Pt()*TMath::Cos(MEtLV.Phi());
     float y = MEtLV.Pt()*TMath::Sin(MEtLV.Phi());
     float m = 80.385;
     
     float delta = c*c*(-4*a*a - 4*b*b - 4*c*c + 4*f*f - 4*m*m - 8*a*x - 8*b*y)*(-4*a*a - 4*b*b - 4*c*c + 4*f*f - 4*m*m - 8*a*x - 8*b*y) - 4*(-4*c*c + 4*f*f)*(-a*a*a*a - 2*a*a*b*b - b*b*b*b - 2*a*a*c*c - 2*b*b*c*c - c*c*c*c + 2*a*a*f*f + 2*b*b*f*f + 2*c*c*f*f - f*f*f*f - 2*a*a*m*m - 2*b*b*m*m - 2*c*c*m*m + 2*f*f*m*m - m*m*m*m - 4*a*a*a*x - 4*a*b*b*x - 4*a*c*c*x + 4*a*f*f*x - 4*a*m*m*x - 4*a*a*x*x + 4*f*f*x*x - 4*a*a*b*y - 4*b*b*b*y - 4*b*c*c*y+ 4*b*f*f*y - 4*b*m*m*y - 8*a*b*x*y - 4*b*b*y*y + 4*f*f*y*y);
    
     //if(delta<0) continue;
     hCount->Fill(2.5,  weight);

     float z1 = delta>0 ? (-c*(-4*a*a - 4*b*b - 4*c*c + 4*f*f - 4*m*m - 8*a*x - 8*b*y) - TMath::Sqrt(delta)) / (2*(-4*c*c + 4*f*f)) : 0.0;
     float z2 = delta>0 ? (-c*(-4*a*a - 4*b*b - 4*c*c + 4*f*f - 4*m*m - 8*a*x - 8*b*y) + TMath::Sqrt(delta)) / (2*(-4*c*c + 4*f*f)) : 0.0;

     LV2 neutrinoLV2a(x,y,z1,TMath::Sqrt(x*x+y*y+z1*z1));
     LV neutrLV1(neutrinoLV2a);
     LV2 neutrinoLV2b(x,y,z2,TMath::Sqrt(x*x+y*y+z2*z2));
     LV neutrLV2(neutrinoLV2b);

     /////////////////////////////////////////////////////////////////////////// 

     std::map<float,  vector<unsigned int>, sorterByLess > posTagged;
     int firstValid1 = 0;
     int firstValid2 = 0;
     int firstValid3 = 0;
     int firstValid4 = 0;
     int num2Jets = 0;
     int num2JetsPt = 0;

     int match1=0; int match2=0;

     int matchB    = 0;
     int matchBbar = 0;
     int otherMatches = 0;
     float resAmb= -99;
     float resB = 0.;
     float resBbar = 0.;
     LV newJet1(0,0,0,0);
     LV newJet2(0,0,0,0);
     LV newJet3(0,0,0,0);
     LV newJet4(0,0,0,0);

     for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
       int genMatchk;
       findGenMatch2(genMatchk, myJetsFilt[k] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);
       float csv_k = mapFilt[k].csv>0 ? mapFilt[k].csv : 0.0 ;
       float pt_k  = mapFilt[k].pt>0 ? mapFilt[k].pt : 0.0 ;
       if( csv_k<0.679 || pt_k<30) continue;
       if( deltaR( myJetsFilt[k],genBLV)<0.3 & genBLV.Pt()>0 ){
	 matchB++;
	 resB = (myJetsFilt[k].Pt()-genBLV.Pt())/genBLV.Pt();
	 if(genMatchk==4){
	   resAmb = resB;
	   hPtMerg->Fill(myJetsFilt[k].Pt(),weight);
	   hMassMerg->Fill(myJetsFilt[k].M()/myJetsFilt[k].Pt(),weight);
	   //newJet1 = LV( );
	 }
	 else
	   hMassUnMerg->Fill(myJetsFilt[k].M()/myJetsFilt[k].Pt(),weight);
       }
       if( deltaR( myJetsFilt[k],genBbarLV)<0.3 & genBbarLV.Pt()>0){
	 matchBbar++;
	 resBbar= (myJetsFilt[k].Pt()-genBbarLV.Pt())/genBbarLV.Pt();
	 if(genMatchk==4){
	   resAmb = resBbar;
	   hPtMerg->Fill(myJetsFilt[k].Pt(),weight);
	   hMassMerg->Fill(myJetsFilt[k].M()/myJetsFilt[k].Pt(),weight);
	 }
	 else
	   hMassUnMerg->Fill(myJetsFilt[k].M()/myJetsFilt[k].Pt(),weight);
       }
     }



     hPtResHiggs->Fill( resB    ,weight);
     hPtResHiggs->Fill( resBbar ,weight);
     hPtResHiggsAmb->Fill( resAmb ,weight);

     hMatchAny->Fill(matchB+matchBbar, weight);

     matchB=0;
     matchBbar=0;
     for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
       int genMatchk;
       findGenMatch2(genMatchk, myJetsFilt[k] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);
       float csv_k = mapFilt[k].csv>0 ? mapFilt[k].csv : 0.0 ;
       float pt_k  = mapFilt[k].pt>0 ? mapFilt[k].pt : 0.0 ;
       //if( csv_k<0.679 || pt_k<30) continue;
       if( deltaR( myJetsFilt[k],genBLV)<0.3 & genBLV.Pt()>0 ){
	 matchB++;
	 resB = (myJetsFilt[k].Pt()-genBLV.Pt())/genBLV.Pt();
	 if(genMatchk==4) resAmb = resB;
       }
       if( deltaR( myJetsFilt[k],genBbarLV)<0.3 & genBbarLV.Pt()>0){
	 matchBbar++;
	 resBbar= (myJetsFilt[k].Pt()-genBbarLV.Pt())/genBbarLV.Pt();
	 if(genMatchk==4) resAmb = resBbar;
       }
     }
     hMatchAll->Fill(matchB+matchBbar, weight);

     matchB=0;
     matchBbar=0;
     for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
       int genMatchk;
       findGenMatch2(genMatchk, myJetsFilt[k] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);
       float csv_k = mapFilt[k].csv>0 ? mapFilt[k].csv : 0.0 ;
       float pt_k  = mapFilt[k].pt>0 ? mapFilt[k].pt : 0.0 ;
       if( csv_k<-999 || pt_k<30) continue;
       if( deltaR( myJetsFilt[k],genBLV)<0.3 & genBLV.Pt()>0 ){
	 matchB++;
	 resB = (myJetsFilt[k].Pt()-genBLV.Pt())/genBLV.Pt();
	 if(genMatchk==4) resAmb = resB;
       }
       if( deltaR( myJetsFilt[k],genBbarLV)<0.3 & genBbarLV.Pt()>0){
	 matchBbar++;
	 resBbar= (myJetsFilt[k].Pt()-genBbarLV.Pt())/genBbarLV.Pt();
	 if(genMatchk==4) resAmb = resBbar;
       }
     }
     hMatchAllPt->Fill(matchB+matchBbar, weight);

     matchB=0;
     matchBbar=0;
     for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
       int genMatchk;
       findGenMatch2(genMatchk, myJetsFilt[k] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);
       float csv_k = mapFilt[k].csv>0 ? mapFilt[k].csv : 0.0 ;
       float pt_k  = mapFilt[k].pt>0 ? mapFilt[k].pt : 0.0 ;
       if( csv_k<0.679 || pt_k<30) continue;
       if( deltaR( myJetsFilt[k],genBLV)<0.3 & genBLV.Pt()>0 ){
	 matchB++;
	 resB = (myJetsFilt[k].Pt()-genBLV.Pt())/genBLV.Pt();
	 if(genMatchk==4) resAmb = resB;
       }
       if( deltaR( myJetsFilt[k],genBbarLV)<0.3 & genBbarLV.Pt()>0){
	 matchBbar++;
	 resBbar= (myJetsFilt[k].Pt()-genBbarLV.Pt())/genBbarLV.Pt();
	 if(genMatchk==4) resAmb = resBbar;
       }
     }
     hMatchAllPtBTag->Fill(matchB+matchBbar, weight);


     //cout << "i = " << i << ", " << myJetsFilt.size() << endl;
     for(unsigned int k = 0; k < myJetsFilt.size()-1; k++){  
       float csv_k = mapFilt[k].csv>0 ? mapFilt[k].csv : 0.0 ;
       float pt_k  = mapFilt[k].pt>0 ? mapFilt[k].pt : 0.0 ;

       //cout << k << endl;

       //if( csv_k<0.679 || pt_k<30) continue;
       if( csv_k<0.0 || pt_k<30) continue;
       int genMatchk;
       findGenMatch2(genMatchk, myJetsFilt[k] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);

       if(genMatchk==5){
	 match1=1;
       }

       for(unsigned int l = k+1; l < myJetsFilt.size(); l++){
	 if( l == k) continue;

	 float csv_l = mapFilt[l].csv>0 ? mapFilt[l].csv : 0.0 ;
	 float pt_l  = mapFilt[l].pt>0 ? mapFilt[l].pt : 0.0 ;
	 //if( csv_l<0.679 || pt_l<30) continue;
	 if( csv_l<0.0 || pt_l<30) continue;

	 //cout << l << endl;

	 int genMatchl;
	 findGenMatch2(genMatchl, myJetsFilt[l] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);

	 if(firstValid1==0){
	   hCount->Fill(3.5,  weight);
	   firstValid1 = 1;
	 }

	 if(genMatchl==-5){
	   match2=1;
	 }

	 num2Jets++;


	 LV HCand     = myJetsFilt[k]+myJetsFilt[l];
	 if( HCand.Pt() < 0 ) continue;

	 num2JetsPt++;

	 if(firstValid2==0){
	   hCount->Fill(4.5,  weight);
	   firstValid2 = 1;
	 }

	 vector<unsigned int> res;
	 res.push_back(k);
	 res.push_back(l);
	 res.push_back(0);
	 res.push_back(0);
	 res.push_back(0);
	 posTagged[ deltaR( myJetsFilt[k], myJetsFilt[l]) ] =  res;
	 //posTagged[ 1./( myJetsFilt[k].Pt()+myJetsFilt[l].Pt()) ] =  res;
	 
       }
     }


     //cout << "num2Jets=" <<num2Jets << endl;
     //cout << "match1,2" << match1 << "," << match2 << endl;
     //cout << "i1=" << i1 << endl;
     //cout << "i2=" << i2 << endl;

     if( match1==1 && match2==1 && num2Jets==1){
       //vector<unsigned int> res = (posTagged.begin())->second; 
       //unsigned int h1 = res[0];
       //unsigned int h2 = res[1];
       //cout << posTagged.size() << "=> " << h1 << ":" << h2 << "  " << (myJetsFilt[h1]+myJetsFilt[h2]).M() << endl;
     }

     hComb->Fill(posTagged.size(),  weight);
     hComb2Jet->Fill(num2Jets,  weight);
     hComb2JetPt->Fill(num2JetsPt,  weight);
     
     int pos = 0;
     int goodCands = 0;
     int goodCands6 = 0;
     int goodCands3 = 0;
     int goodCands1 = 0;

     hMassJJOneEventR->Reset();

     //if(posTagged.size()>0) cout << "Event " << i << endl;
     for(std::map<float, vector<unsigned int> >::iterator it = posTagged.begin(); it!= posTagged.end() ; it++){
       vector<unsigned int> res = it->second; 
       unsigned int h1 = res[0];
       unsigned int h2 = res[1];
    

       int genMatchh1;
       findGenMatch2(genMatchh1, myJetsFilt[h1] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);
       int genMatchh2;
       findGenMatch2(genMatchh2, myJetsFilt[h2] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);
	 
       //cout << "Permut #" << pos << ": Pt/M = " << it->first << ", comb = [" << res[0] << res[1] << res[2] << res[3] << "]" <<endl;
       ///cout << "Higgs candidates at [" << h1 << h2 << "]" << " => pairing = [" << genMatchh1 << genMatchh2 << "]" << endl;

       bool isMatched =
	 (genMatchh1==5 && genMatchh2==-5) || (genMatchh2==5 && genMatchh1==-5);

       bool cut = num2Jets>0 || 
	 (num2Jets==1 && deltaR( myJetsFilt[h1], myJetsFilt[h2]) < 999 ) ||
	 (num2Jets==3 && deltaR( myJetsFilt[h1], myJetsFilt[h2]) < 999 ) ||
	 (num2Jets==6 && deltaR( myJetsFilt[h1], myJetsFilt[h2]) < 999 ) ;

       if(num2Jets==6 && deltaR( myJetsFilt[h1], myJetsFilt[h2]) < 999 && firstOcc==0){
	 hMassJJOneEventR->Fill((myJetsFilt[h1]+ myJetsFilt[h2]).M());
       }

       if(cut) goodCands++;
       if(num2Jets==6 && deltaR( myJetsFilt[h1], myJetsFilt[h2]) < 1.6 && (myJetsFilt[h1]+ myJetsFilt[h2]).M()>0 && (myJetsFilt[h1]+ myJetsFilt[h2]).M()<1400){
	 goodCands6++;
	 hMassJJOneEventR->Fill((myJetsFilt[h1]+ myJetsFilt[h2]).M());
       }
       if(num2Jets==3 && deltaR( myJetsFilt[h1], myJetsFilt[h2]) < 1.6 && (myJetsFilt[h1]+ myJetsFilt[h2]).M()>0 && (myJetsFilt[h1]+ myJetsFilt[h2]).M()<1400){
	 goodCands3++;
	 hMassJJOneEventR->Fill((myJetsFilt[h1]+ myJetsFilt[h2]).M());
       }
       if(num2Jets==1 && deltaR( myJetsFilt[h1], myJetsFilt[h2]) < 1.6 && (myJetsFilt[h1]+ myJetsFilt[h2]).M()>0 && (myJetsFilt[h1]+ myJetsFilt[h2]).M()<1400){
	 goodCands1++;
	 hMassJJOneEventR->Fill((myJetsFilt[h1]+ myJetsFilt[h2]).M());
       }
       

       if( cut && num2Jets==6 && num2Jets<=999 && deltaR( myJetsFilt[h1], myJetsFilt[h2]) < 999 ){

	 if(isMatched){
	   if( goodCands==1 ) hMassCat3It1R->Fill((myJetsFilt[h1]+ myJetsFilt[h2]).M() , weight);
	   if( goodCands==2 ) hMassCat3It2R->Fill((myJetsFilt[h1]+ myJetsFilt[h2]).M() , weight);
	   if( goodCands==3 ) hMassCat3It3R->Fill((myJetsFilt[h1]+ myJetsFilt[h2]).M() , weight);
	   if( goodCands==4 ) hMassCat3It4R->Fill((myJetsFilt[h1]+ myJetsFilt[h2]).M() , weight);
	   if( goodCands==5 ) hMassCat3It5R->Fill((myJetsFilt[h1]+ myJetsFilt[h2]).M() , weight);
	   if( goodCands==6 ) hMassCat3It6R->Fill((myJetsFilt[h1]+ myJetsFilt[h2]).M() , weight);
	 }else{
	   if( goodCands==1 ) hMassCat3It1W->Fill((myJetsFilt[h1]+ myJetsFilt[h2]).M() , weight);
	   if( goodCands==2 ) hMassCat3It2W->Fill((myJetsFilt[h1]+ myJetsFilt[h2]).M() , weight);
	   if( goodCands==3 ) hMassCat3It3W->Fill((myJetsFilt[h1]+ myJetsFilt[h2]).M() , weight);
	   if( goodCands==4 ) hMassCat3It4W->Fill((myJetsFilt[h1]+ myJetsFilt[h2]).M() , weight);
	   if( goodCands==5 ) hMassCat3It5W->Fill((myJetsFilt[h1]+ myJetsFilt[h2]).M() , weight);
	   if( goodCands==6 ) hMassCat3It6W->Fill((myJetsFilt[h1]+ myJetsFilt[h2]).M() , weight);
	 }
	 
       }


       if( cut && goodCands==1){

	 if(num2Jets<=2 && deltaR( myJetsFilt[h1], myJetsFilt[h2]) < 1.6 && isMatched)
	   hMassCat1R->Fill((myJetsFilt[h1]+ myJetsFilt[h2]).M() , weight);
	 if(num2Jets>=3 && num2Jets<=5  && deltaR( myJetsFilt[h1], myJetsFilt[h2]) < 1.6 && isMatched)
	   hMassCat2R->Fill((myJetsFilt[h1]+ myJetsFilt[h2]).M() , weight);
	 if(num2Jets>=6 && num2Jets<=999 && deltaR( myJetsFilt[h1], myJetsFilt[h2]) < 999 && isMatched){
	   hMassCat3R->Fill((myJetsFilt[h1]+ myJetsFilt[h2]).M() , weight);

	   if((myJetsFilt[h1]+ myJetsFilt[h2]).M()>100 && (myJetsFilt[h1]+ myJetsFilt[h2]).M()<140){
	     hPt1R->Fill(myJetsFilt[h1].Pt(),weight);
	     hPt2R->Fill(myJetsFilt[h2].Pt(),weight);
	     hPtHR->Fill( (myJetsFilt[h2]+myJetsFilt[h1]).Pt(),weight);
	     hPtSumR->Fill( myJetsFilt[h1].Pt()+myJetsFilt[h2].Pt(), weight);
	     hBTag1R->Fill( mapFilt[h1].csv, weight);
	     hBTag2R->Fill( mapFilt[h2].csv, weight);
	     hMass1R->Fill( myJetsFilt[h1].M(), weight);
	     hMass2R->Fill( myJetsFilt[h2].M(), weight);
	   }



	 }
	 
	 if(num2Jets<=2 && deltaR( myJetsFilt[h1], myJetsFilt[h2]) < 1.6 && !isMatched)
	   hMassCat1W->Fill((myJetsFilt[h1]+ myJetsFilt[h2]).M() , weight);
	 if(num2Jets>=3 && num2Jets<=5  && deltaR( myJetsFilt[h1], myJetsFilt[h2]) < 1.6 && !isMatched)
	   hMassCat2W->Fill((myJetsFilt[h1]+ myJetsFilt[h2]).M() , weight);
	 if(num2Jets>=6 && num2Jets<=999 && deltaR( myJetsFilt[h1], myJetsFilt[h2]) < 999 && !isMatched){
	   hMassCat3W->Fill((myJetsFilt[h1]+ myJetsFilt[h2]).M() , weight);

	   if((myJetsFilt[h1]+ myJetsFilt[h2]).M()>100 && (myJetsFilt[h1]+ myJetsFilt[h2]).M()<140){
	     hPt1W->Fill(myJetsFilt[h1].Pt(),weight);
	     hPt2W->Fill(myJetsFilt[h2].Pt(),weight);
	     hPtHW->Fill( (myJetsFilt[h2]+myJetsFilt[h1]).Pt(),weight);
	     hPtSumW->Fill( myJetsFilt[h1].Pt()+myJetsFilt[h2].Pt(), weight);
	     hBTag1W->Fill( mapFilt[h1].csv, weight);
	     hBTag2W->Fill( mapFilt[h2].csv, weight);
	     hMass1W->Fill( myJetsFilt[h1].M(), weight);
	     hMass2W->Fill( myJetsFilt[h2].M(), weight);
	   }
	 }
	 
	


	 if( isMatched ){
	   hMassJJR->Fill( (myJetsFilt[h1]+ myJetsFilt[h2]).M() , weight);
	   hDPhiR->Fill( myJetsFilt[h1].Phi()-myJetsFilt[h2].Phi(), weight);
	   hDRR->Fill( deltaR( myJetsFilt[h1], myJetsFilt[h2]), weight );
	 }
	 else{
	   hMassJJW->Fill( (myJetsFilt[h1]+ myJetsFilt[h2]).M() , weight);
	   hDPhiW->Fill( myJetsFilt[h1].Phi()-myJetsFilt[h2].Phi(), weight);
	   hDRW->Fill( deltaR( myJetsFilt[h1], myJetsFilt[h2]), weight );
	 }
	 
       }

       pos++;
     }

     float frac =  hMassJJOneEventR->Integral()>0 ? hMassJJOneEventR->Integral(hMassJJOneEventR->FindBin(105.),hMassJJOneEventR->FindBin(130.))/ hMassJJOneEventR->Integral() : 0.0;
     if(num2Jets==1) hFrac1->Fill( frac*goodCands1, weight);
     if(num2Jets==3) hFrac3->Fill( frac*goodCands3, weight);
     if(num2Jets==6) hFrac6->Fill( frac*goodCands6, weight);

     hCombPass6->Fill(goodCands6, weight);
     hCombPass3->Fill(goodCands3, weight);
     hCombPass1->Fill(goodCands1, weight);

     if(hMassJJOneEventR->Integral()>0)
       firstOcc = 1;

     hMatch->Fill(match1+match2, weight);

     if((match1+match2)==1) hMatchDebug->Fill(match1,match2,weight);

     hMatchComb->Fill(match1+match2, num2Jets, weight);
    

  }


   /////////////////////////////////////////////////////////////////////////// 
   
   
   fOutput->Write();
   fOutput->Close();
   
   fLikelihoods->Close();
   fInput->Close();
   
   delete ran; delete fOutput; delete fLikelihoods;
   
   return  ;
   
}





std::pair<float,float> algo( float cutW_ = 4.0, float cutW2_ = 4.0, float csvCut_ = 0.60, float likeTopHadCut_ = 0.00004, float likeTopLepCut_ = 0.001, float ptCut_ = 30, float likeTopHadCut2_ = 0.00004, float likeTopLepCut2_ = 0.001){

  TRandom3* ran = new TRandom3();

  int nSteps = 5; 
  TFile *fOutput           = new TFile("algo.root","RECREATE");
  TH1F* hCounterR          = new TH1F("hCounterR", "hCounterR", nSteps, 0., nSteps); 
  TH1F* hCounterW          = new TH1F("hCounterW", "hCounterW", nSteps, 0., nSteps); 
  TH1F* hNonBLogLikeS      = new TH1F("hNonBLogLikeS", "hNonBLogLike signal",100, 0 ,100 ); 
  TH1F* hNonBLogLikeB      = new TH1F("hNonBLogLikeB", "hNonBLogLike bkg"   ,100, 0 ,100 );
  TH1F* hWHadLikeS         = new TH1F("hWHadLikeS", "hWHadLike sgn"   ,300, 0 ,30 );
  TH1F* hWHadLikeB         = new TH1F("hWHadLikeB", "hWHadLike bkg"   ,300, 0 ,30 );
  TH1F* hWHadRank          = new TH1F("hWHadRank", "good pos"   ,100, 0 ,100 );
  TH1F* hWHadLikeSRank0    = new TH1F("hWHadLikeSRank0", "hWHadLike sgn"   ,300, 0 ,30 );
  TH1F* hWHadLikeBRank0    = new TH1F("hWHadLikeBRank0", "hWHadLike bkg"   ,300, 0 ,30 );
  TH1F* hWHadLikeSRank1    = new TH1F("hWHadLikeSRank1", "hWHadLike sgn"   ,300, 0 ,30 );
  TH1F* hWHadLikeBRank1    = new TH1F("hWHadLikeBRank1", "hWHadLike bkg"   ,300, 0 ,30 );
  TH1F* hMassJJ            = new TH1F("hMassJJ", "Higgs mass"   ,80, 0 ,300 );

  TH1F* hWHadrMass = 0;
 

  TFile* fLikelihoods          = new TFile("likelihoods_bkgd.root","READ");
  TF2* fLikelihoodWTopHadrMass = (TF2*)fLikelihoods->Get("fLikelihoodWTopHadrMass");
  TF1* fLikelihoodCsvB         = (TF1*)fLikelihoods->Get("fLikelihoodCsvB");
  TF1* fLikelihoodCsvNonB      = (TF1*)fLikelihoods->Get("fLikelihoodCsvNonB");
  TF1* fLikelihoodTopLeptMass  = (TF1*)fLikelihoods->Get("fLikelihoodTopLeptMass");

  TFile* fHistos               = new TFile("histos_bkgd.root","READ");
  TH2F* hWTopHadrMass          = (TH2F*)fHistos->Get("hWTopHadrMass");
  hWHadrMass = (TH1F*)hWTopHadrMass->ProjectionX("_px");

  for(int k = 0; k < 10000 ; k++){
    double csv = ran->Uniform(-1E-03,1+1E-03);
    if(fLikelihoodCsvB->Eval(csv)>0)    hNonBLogLikeS->Fill( -2*TMath::Log( fLikelihoodCsvNonB->Eval(csv)/fLikelihoodCsvB->Eval(csv) ), fLikelihoodCsvNonB->Eval(csv));
    if(fLikelihoodCsvNonB->Eval(csv)>0) hNonBLogLikeB->Fill( -2*TMath::Log( fLikelihoodCsvNonB->Eval(csv)/fLikelihoodCsvB->Eval(csv) ), fLikelihoodCsvB->Eval(csv));
  }
  hNonBLogLikeS->Scale(1./hNonBLogLikeS->Integral());
  hNonBLogLikeB->Scale(1./hNonBLogLikeB->Integral());

  //fOutput->cd();
  //hWHadrMass->Write("hWTopHadrMass_px");
  //fOutput->Write();
  //fOutput->Close();
  //return;


  TFile* fInput = TFile::Open("gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_VType2/v2/DiJetPt_TTH_HToBB_M-125_8TeV-pythia6_VType2.root","READ");
  //TFile* fInput = TFile::Open("gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_VType2/v2/","READ");
  TTree* tree = (TTree*)fInput->Get("tree");

  if(!tree){
    cout << "No tree" << endl; return make_pair(-99.,-99);
  }

  JetByPt jet1, jet2, jet3, jet4, jet5, jet6, jet7, jet8, jet9, jet10;
  genTopInfo genTop, genTbar;
  genParticleInfo genB, genBbar;
  int nvlep;
  float vLeptonpt[999]; float vLeptoneta[999]; float vLeptonphi[999];
  float met[999];

  tree->SetBranchAddress("jet1",   &jet1);
  tree->SetBranchAddress("jet2",   &jet2);
  tree->SetBranchAddress("jet3",   &jet3);
  tree->SetBranchAddress("jet4",   &jet4);
  tree->SetBranchAddress("jet5",   &jet5);
  tree->SetBranchAddress("jet6",   &jet6);
  tree->SetBranchAddress("jet7",   &jet7);
  tree->SetBranchAddress("jet8",   &jet8);
  tree->SetBranchAddress("jet9",   &jet9);
  tree->SetBranchAddress("jet10",  &jet10);
  tree->SetBranchAddress("genTop", &genTop);
  tree->SetBranchAddress("genTbar",&genTbar);
  tree->SetBranchAddress("genB",   &genB);
  tree->SetBranchAddress("genBbar",&genBbar);
  tree->SetBranchAddress("vLepton_pt", vLeptonpt);
  tree->SetBranchAddress("vLepton_phi",vLeptonphi);
  tree->SetBranchAddress("vLepton_eta",vLeptoneta);
  tree->SetBranchAddress("nvlep",  &nvlep);
  tree->SetBranchAddress("METtype1p2corr",  met);

  Long64_t nentries = tree->GetEntries();
  int totCounter = 0;

   for (Long64_t i = 0; i < nentries; i++){
     
     //if(i%100000==0) cout << i << endl;
     tree->GetEntry(i);

     //continue;

     if(nvlep<=0) continue;


     LV topBLV(  0.,0.,0.,0.);
     LV topW1LV( 0.,0.,0.,0.); 
     LV topW2LV( 0.,0.,0.,0.);
     LV atopBLV( 0.,0.,0.,0.); 
     LV atopW1LV(0.,0.,0.,0.); 
     LV atopW2LV(0.,0.,0.,0.);
     LV genBLV(0.,0.,0.,0.);
     LV genBbarLV(0.,0.,0.,0.);
    
    
     ///////////////////////////////////////////////////////////////////
     /////////////////////////////////////////////////////////////////// 
     if(genTop.bmass>0){
       topBLV.SetPt(   genTop.bpt );
       topBLV.SetEta(  genTop.beta );
       topBLV.SetPhi(  genTop.bphi );
       topBLV.SetM(    genTop.bmass );
       topW1LV.SetPt(  genTop.wdau1pt );
       topW1LV.SetEta( genTop.wdau1eta );
       topW1LV.SetPhi( genTop.wdau1phi );
       topW1LV.SetM(   genTop.wdau1mass );
       topW2LV.SetPt(  genTop.wdau2pt );
       topW2LV.SetEta( genTop.wdau2eta );
       topW2LV.SetPhi( genTop.wdau2phi );
       topW2LV.SetM(   genTop.wdau2mass );
     }
     if(genTbar.bmass>0){
	atopBLV.SetPt(   genTbar.bpt );
	atopBLV.SetEta(  genTbar.beta );
	atopBLV.SetPhi(  genTbar.bphi );
	atopBLV.SetM(    genTbar.bmass );
	atopW1LV.SetPt(  genTbar.wdau1pt );
	atopW1LV.SetEta( genTbar.wdau1eta );
	atopW1LV.SetPhi( genTbar.wdau1phi );
	atopW1LV.SetM(   genTbar.wdau1mass );
	atopW2LV.SetPt(  genTbar.wdau2pt );
	atopW2LV.SetEta( genTbar.wdau2eta );
	atopW2LV.SetPhi( genTbar.wdau2phi );
	atopW2LV.SetM(   genTbar.wdau2mass );
     }
     if(genB.mass>0 && genB.momid==25){
       genBLV.SetPt(  genB.pt );
       genBLV.SetEta( genB.eta );
       genBLV.SetPhi( genB.phi );
       genBLV.SetM(   genB.mass );
     }
     if(genBbar.mass>0 && genBbar.momid==25){
       genBbarLV.SetPt(  genBbar.pt );
       genBbarLV.SetEta( genBbar.eta );
       genBbarLV.SetPhi( genBbar.phi );
       genBbarLV.SetM(   genBbar.mass );
     }
     ///////////////////////////////////////////////////////////////////
     ///////////////////////////////////////////////////////////////////

     vector<LV> myJets;
     std::map<unsigned int, JetByPt> map;

     map[0] = jet1;  map[1] = jet2; map[2] = jet3; map[3] = jet4; 
     map[4] = jet5;  map[5] = jet6; map[6] = jet7; map[7] = jet8; 
     map[8] = jet9;  map[9] = jet10; 

     myJets.push_back(LV(jet1.pt, jet1.eta, jet1.phi, jet1.mass));
     myJets.push_back(LV(jet2.pt, jet2.eta, jet2.phi, jet2.mass));
     myJets.push_back(LV(jet3.pt, jet3.eta, jet3.phi, jet3.mass));
     myJets.push_back(LV(jet4.pt, jet4.eta, jet4.phi, jet4.mass));
     myJets.push_back(LV(jet5.pt, jet5.eta, jet5.phi, jet5.mass));
     myJets.push_back(LV(jet6.pt, jet6.eta, jet6.phi, jet6.mass));
     myJets.push_back(LV(jet7.pt, jet7.eta, jet7.phi, jet7.mass));
     myJets.push_back(LV(jet8.pt, jet8.eta, jet8.phi, jet8.mass));
     myJets.push_back(LV(jet9.pt, jet9.eta, jet9.phi, jet9.mass));
     myJets.push_back(LV(jet10.pt,jet10.eta,jet10.phi,jet10.mass));

     
     vector<LV> myJetsFilt;
     std::map<unsigned int, JetByPt> mapFilt;

     unsigned int iter = 0;
     for(unsigned int k = 0; k<myJets.size() ; k++){

       if(myJets[k].Pt()>ptCut_ && TMath::Abs(myJets[k].Eta()) < 2.5 ){

	 myJetsFilt.push_back( myJets[k] );

	 switch(k){
	 case 0:
	   mapFilt[iter] = jet1;
	   break;
	 case 1:
	   mapFilt[iter] = jet2;
	   break;
	 case 2:
	   mapFilt[iter] = jet3;
	   break;
	 case 3:
	   mapFilt[iter] = jet4;
	   break;
	 case 4:
	   mapFilt[iter] = jet5;
	   break;
	 case 5:
	   mapFilt[iter] = jet6;
	   break;
	 case 6:
	   mapFilt[iter] = jet7;
	   break;
	 case 7:
	   mapFilt[iter] = jet8;
	   break;
	 case 8:
	   mapFilt[iter] = jet9;
	   break;
	 case 9:
	   mapFilt[iter] = jet10;
	   break;
	 default:
	   break;
	 }
	 
	 iter++;
       }
     }

     LV lepLV(vLeptonpt[0], vLeptoneta[0], vLeptonphi[0], 0.);
     LV MEtLV(met[0], 0., met[3] , 0.);

     ///////////////////////////////////////////////////////////////////////////

     //if(  !((abs(genTop.wdau1id)<6  && genTop.wdau1pt>40  && TMath::Abs(genTop.wdau1eta) <2.5 && genTop.wdau2pt>40  && TMath::Abs(genTop.wdau2eta) <2.5)  || 
     //   (abs(genTbar.wdau1id)<6 && genTbar.wdau1pt>40 && TMath::Abs(genTbar.wdau1eta)<2.5 && genTbar.wdau2pt>40 && TMath::Abs(genTbar.wdau2eta)<2.5))
     //  || !(genB.pt>30 && genBbar.pt>30)
     // ) continue;

     totCounter++;

     if( myJetsFilt.size()<2) continue;

     /////////////////////////////////////////////////////////////////////////// STEP2

     std::map<float,  int, sorterByLess> rankWHad;
     std::map<float,  pair<unsigned int,unsigned int>, sorterByLess > posWHad;

     for(unsigned int k = 0; k < myJetsFilt.size()-1; k++){  
       float csv_k = mapFilt[k].csv>0 ? mapFilt[k].csv : 0.0 ;

       for(unsigned int l = k+1; l < myJetsFilt.size(); l++){
	 float csv_l = mapFilt[l].csv>0 ? mapFilt[l].csv : 0.0 ;

	 LV diJet  = myJetsFilt[k]+myJetsFilt[l];
	 //LV diJetU = myJetsFilt[k]*(1+mapFilt[k].unc) + myJetsFilt[l]*(1+mapFilt[l].unc);
	 //LV diJetD = myJetsFilt[k]*(1-mapFilt[k].unc) + myJetsFilt[l]*(1-mapFilt[l].unc);
	 //float dM = TMath::Max( TMath::Abs(diJetU.M()-diJet.M()), TMath::Abs(diJetD.M()-diJet.M()) );
	 //dM = TMath::Sqrt(dM*dM + 2.08*2.08);
	 //TF1 *fGaus = new TF1("fGaus",Form("1/sqrt(2*TMath::Pi())/%f*exp(-0.5*((x-%f)/%f)**2)", dM, diJet.M(), dM),0,400);
	 //bool pass = (fGaus->Eval( 80.325 ) > fGaus->Eval( diJet.M()+1.96*dM ));
	 //bool pass = csv_l < 0.50 &&  csv_k < 0.50;
	 float likeB1   = fLikelihoodCsvNonB->Eval( csv_k );
	 float likeB2   = fLikelihoodCsvNonB->Eval( csv_l );
	 float likeWHad = hWHadrMass->GetBinContent( hWHadrMass->FindBin( diJet.M() ) );
	 float logLike  = likeB1*likeB2*likeWHad>0 ? -TMath::Log( likeB1*likeB2*likeWHad ) : 25 + ran->Uniform(0.,1);
	 //float logLike  = fGaus->Eval( 80.325 );

	 int genMatchk;
	 findGenMatch(genMatchk, myJetsFilt[k] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);
	 int genMatchl;
	 findGenMatch(genMatchl, myJetsFilt[l] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);

	 if( abs(genMatchk)==2 && abs(genMatchl)==2 )
	   hWHadLikeS->Fill(logLike);
	 else 
	   hWHadLikeB->Fill(logLike);

	 rankWHad[logLike] = int( abs(genMatchk)==2 && abs(genMatchl)==2 );
	 posWHad[logLike]  = make_pair(k,l);

	 //delete fGaus;	 
       }
     }

     int pos=0;
     for(std::map<float, int>::iterator it = rankWHad.begin(); it!= rankWHad.end() ; it++){
       if(it->second>0.5)
	 hWHadRank->Fill( pos );
       if(pos==0 && it->second>0.5)
	 hWHadLikeSRank0->Fill( it->first );
       if(pos==0 && it->second<0.5)
	 hWHadLikeBRank0->Fill( it->first );
       if(pos==1 && it->second>0.5)
	 hWHadLikeSRank1->Fill( it->first );
       if(pos==1 && it->second<0.5)
	 hWHadLikeBRank1->Fill( it->first );
       pos++;
     }

     if( rankWHad.size()<1 ) continue;
     if(  (rankWHad.begin())->first > cutW_ || (rankWHad.size()>1 && (++rankWHad.begin())->first < cutW2_) ) continue;

     unsigned int W1 = posWHad.size()>0 ? ((posWHad.begin())->second).first : 999;
     unsigned int W2 = posWHad.size()>0 ? ((posWHad.begin())->second).second: 999;

     if( W1==999 || W2==999 ) continue;

     LV w1LV = W1<999 ? myJetsFilt[W1] : LV(0,0,0,0);
     LV w2LV = W2<999 ? myJetsFilt[W2] : LV(0,0,0,0);

     int genMatchW1;
     findGenMatch(genMatchW1, w1LV , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);
     int genMatchW2;
     findGenMatch(genMatchW2, w2LV , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);
     if(abs(genMatchW1)==2 && abs(genMatchW2)==2 ) hCounterR->Fill(0.5);
     else  hCounterW->Fill(0.5);

     //continue;

     /////////////////////////////////////////////////////////////////////////// STEP2

     std::vector<unsigned int> unMasked;
     for(unsigned int k = 0; k < myJetsFilt.size(); k++)
       if( k!= W1 && k!=W2 && mapFilt[k].csv> csvCut_) unMasked.push_back( k );

     if( unMasked.size()<2 ) continue;
     if( abs(genMatchW1)==2 && abs(genMatchW2)==2 ) hCounterR->Fill(1.5);
     else  hCounterW->Fill(1.5);

     /////////////////////////////////////////////////////////////////////////// STEP3

     float a = lepLV.Px();
     float b = lepLV.Py();
     float c = lepLV.Pz();
     float f = lepLV.E();
     float x = MEtLV.Pt()*TMath::Cos(MEtLV.Phi());
     float y = MEtLV.Pt()*TMath::Sin(MEtLV.Phi());
     float m = 80.385;
     
     float delta = c*c*(-4*a*a - 4*b*b - 4*c*c + 4*f*f - 4*m*m - 8*a*x - 8*b*y)*(-4*a*a - 4*b*b - 4*c*c + 4*f*f - 4*m*m - 8*a*x - 8*b*y) - 4*(-4*c*c + 4*f*f)*(-a*a*a*a - 2*a*a*b*b - b*b*b*b - 2*a*a*c*c - 2*b*b*c*c - c*c*c*c + 2*a*a*f*f + 2*b*b*f*f + 2*c*c*f*f - f*f*f*f - 2*a*a*m*m - 2*b*b*m*m - 2*c*c*m*m + 2*f*f*m*m - m*m*m*m - 4*a*a*a*x - 4*a*b*b*x - 4*a*c*c*x + 4*a*f*f*x - 4*a*m*m*x - 4*a*a*x*x + 4*f*f*x*x - 4*a*a*b*y - 4*b*b*b*y - 4*b*c*c*y+ 4*b*f*f*y - 4*b*m*m*y - 8*a*b*x*y - 4*b*b*y*y + 4*f*f*y*y);
    
     if(delta<0) continue;
     if( abs(genMatchW1)==2 && abs(genMatchW2)==2 ) hCounterR->Fill(2.5);
     else  hCounterW->Fill(2.5);

     float z1 = (-c*(-4*a*a - 4*b*b - 4*c*c + 4*f*f - 4*m*m - 8*a*x - 8*b*y) - TMath::Sqrt(delta)) / (2*(-4*c*c + 4*f*f));
     float z2 = (-c*(-4*a*a - 4*b*b - 4*c*c + 4*f*f - 4*m*m - 8*a*x - 8*b*y) + TMath::Sqrt(delta)) / (2*(-4*c*c + 4*f*f));

     LV2 neutrinoLV2a(x,y,z1,TMath::Sqrt(x*x+y*y+z1*z1));
     LV neutrLV1(neutrinoLV2a);
     LV2 neutrinoLV2b(x,y,z2,TMath::Sqrt(x*x+y*y+z2*z2));
     LV neutrLV2(neutrinoLV2b);


     /////////////////////////////////////////////////////////////////////////// STEP4

     std::map<float,  pair<unsigned int,unsigned int> , sorterByLess> posHiggs;
     vector<unsigned int> unMasked2;

     int counterTopHad = 0;
     int counterTopLep = 0;

     for(unsigned int k = 0; k < unMasked.size(); k++){

       bool mask = false;

       unsigned int j1 = unMasked[k];
       LV j1LV = myJetsFilt[j1];
       
       float likeTopHad = fLikelihoodWTopHadrMass->Eval( (w1LV+w2LV).M(), (w1LV+w2LV+j1LV).M());
       float likeTopLep = TMath::Max(fLikelihoodTopLeptMass->Eval( (lepLV+neutrLV1+j1LV).M()),
				     fLikelihoodTopLeptMass->Eval( (lepLV+neutrLV2+j1LV).M())  );

       if( likeTopHad>= likeTopHadCut2_ ) counterTopHad++;
       if( likeTopLep>= likeTopLepCut2_ ) counterTopLep++;

       ///cout << likeTopHad << " - " << likeTopLep << endl;
       if( likeTopHad > likeTopHadCut_ || likeTopLep > likeTopLepCut_ ){
	 mask = true;
       }
       
       if(!mask) unMasked2.push_back( j1 );
	 
     }

     if(unMasked2.size()<2) continue;
     if(counterTopHad<1 ||  counterTopLep<1 ) continue;

     if( abs(genMatchW1)==2 && abs(genMatchW2)==2 ) hCounterR->Fill(3.5);
     else  hCounterW->Fill(3.5);

     for(unsigned int k = 0; k < unMasked2.size()-1; k++){
       for(unsigned int l = k+1; l < unMasked2.size(); l++){
	 unsigned int j1 = unMasked2[k];
	 unsigned int j2 = unMasked2[l];

	 LV h1LV = myJetsFilt[j1];
	 LV h2LV = myJetsFilt[j2];

	 //posHiggs[ (h1LV+h2LV).M() / (h1LV+h2LV).Pt()  ] = make_pair(j1,j2);
	 posHiggs[ 1./( mapFilt[j1].csv + mapFilt[j2].csv)   ] = make_pair(j1,j2);
	 //posHiggs[ TMath::Abs((h1LV+h2LV).M() - 120.)   ] = make_pair(j1,j2);
       }
     }
     
  
     unsigned int H1 = ((posHiggs.begin())->second).first;
     unsigned int H2 = ((posHiggs.begin())->second).second;

     int genMatchH1;
     findGenMatch(genMatchH1, myJetsFilt[H1] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);
     int genMatchH2;
     findGenMatch(genMatchH2, myJetsFilt[H2] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);

     if(  genMatchH1==0 && genMatchH2==0) hCounterR->Fill(4.5);
     else  hCounterW->Fill(4.5);

     //if(genMatchH1==0 && genMatchH2==0)
     hMassJJ->Fill( (myJetsFilt[H1] + myJetsFilt[H2]).M() );


   } // end LOOP


   hCounterW->Scale(1./totCounter);
   hCounterR->Scale(1./totCounter);

//    cout << "Step 1 " << endl;
//    cout << "Wrong = " << hCounterW->GetBinContent(1) << endl;
//    cout << "Right = " << hCounterR->GetBinContent(1) << endl;

//    cout << "Step 2 " << endl;
//    cout << "Wrong = " << hCounterW->GetBinContent(2) << endl;
//    cout << "Right = " << hCounterR->GetBinContent(2) << endl;

//    cout << "Step 3 " << endl;
//    cout << "Wrong = " << hCounterW->GetBinContent(3) << endl;
//    cout << "Right = " << hCounterR->GetBinContent(3) << endl;

//    cout << "Step 4 " << endl;
//    cout << "Wrong = " << hCounterW->GetBinContent(4) << endl;
//    cout << "Right = " << hCounterR->GetBinContent(4) << endl;

//    cout << "Step 5 " << endl;
//    cout << "Wrong = " << hCounterW->GetBinContent(5) << endl;
//    cout << "Right = " << hCounterR->GetBinContent(5) << " (entries = " << hCounterR->GetBinContent(5)*totCounter  << ")" << endl;


   hWHadLikeS->Scale(1./hWHadLikeS->Integral());
   hWHadLikeB->Scale(1./hWHadLikeB->Integral());

   hWHadLikeSRank0->Scale(1./hWHadLikeSRank0->Integral());
   hWHadLikeBRank0->Scale(1./hWHadLikeBRank0->Integral());
   hWHadLikeSRank1->Scale(1./hWHadLikeSRank1->Integral());
   hWHadLikeBRank1->Scale(1./hWHadLikeBRank1->Integral());

   float res1 = (hCounterR->GetBinContent(5)/hCounterW->GetBinContent(5));
   float res2 =  hCounterR->GetBinContent(5);

   fOutput->Write();
   fOutput->Close();

   fLikelihoods->Close();
   fInput->Close();

   delete ran; delete fOutput; delete fLikelihoods;

   return make_pair(res1,res2) ;

}



void algoTest(){

  float cutW[]   = {2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5};
  float cutCsv[] = {0.4,0.6,0.8};
  float likeTopHadCut[] = {0.00002, 0.00004, 0.00006, 0.00008, 0.00010, 0.0002};
  float likeTopLepCut[] = {0.0005, 0.001, 0.0015, 0.020, 0.030, 0.040};
  float ptCut[] = {20,30,40,50};

  float maxResult = -99;

  for(int i = 0; i < 7; i++){
    for(int j = 0; j < 3; j++){
      for(int k = 0; k < 6; k++){
	for(int l = 0; l < 6; l++){
	  for(int m = 0; m < 4; m++){
	    
	    std::pair<float,float> result = algo( cutW[i], 0.0, cutCsv[j], likeTopHadCut[k], likeTopLepCut[l], ptCut[m], -99, -99);
	    
	    if(result.first>maxResult && result.second>0.001 && result.second<0.0015){
	      maxResult = result.first;
	      cout << "Eff: " << result.second << ", Pur: " << result.first << " ==> " << cutW[i] << ", " << cutCsv[j] << ", " << likeTopHadCut[k] << ", " << likeTopLepCut[l] << ", " << ptCut[m] << endl;
	    }
	  }
	}
      }
    }
  }

  cout << maxResult << "" << endl;
  
}


void algoPre( float cutW_ = 4.0, float cutW2_ = 4.0, float csvCut_ = 0.60, float likeTopHadCut_ = 0.00004, float likeTopLepCut_ = 0.001, float ptCut_ = 30, float likeTopHadCut2_ = 0.00004, float likeTopLepCut2_ = 0.001){

  TRandom3* ran = new TRandom3();

  int nSteps = 5; 
  TFile *fOutput           = new TFile("algoPre.root","RECREATE");
  TH1F* hCounterR          = new TH1F("hCounterR", "hCounterR", nSteps, 0., nSteps); 
  TH1F* hCounterW          = new TH1F("hCounterW", "hCounterW", nSteps, 0., nSteps); 
  TH1F* hNonBLogLikeS      = new TH1F("hNonBLogLikeS", "hNonBLogLike signal",100, 0 ,100 ); 
  TH1F* hNonBLogLikeB      = new TH1F("hNonBLogLikeB", "hNonBLogLike bkg"   ,100, 0 ,100 );
  TH1F* hWHadLikeS         = new TH1F("hWHadLikeS", "hWHadLike sgn"   ,300, 0 ,30 );
  TH1F* hWHadLikeB         = new TH1F("hWHadLikeB", "hWHadLike bkg"   ,300, 0 ,30 );
  TH1F* hWHadRank          = new TH1F("hWHadRank", "good pos"   ,100, 0 ,100 );
  TH1F* hWHadLikeSRank0    = new TH1F("hWHadLikeSRank0", "hWHadLike sgn"   ,300, 0 ,30 );
  TH1F* hWHadLikeBRank0    = new TH1F("hWHadLikeBRank0", "hWHadLike bkg"   ,300, 0 ,30 );
  TH1F* hWHadLikeSRank1    = new TH1F("hWHadLikeSRank1", "hWHadLike sgn"   ,300, 0 ,30 );
  TH1F* hWHadLikeBRank1    = new TH1F("hWHadLikeBRank1", "hWHadLike bkg"   ,300, 0 ,30 );
  TH1F* hMassJJ            = new TH1F("hMassJJ", "Higgs mass"   ,80, 0 ,300 );

  TH1F* hWHadrMass = 0;
 

  TFile* fLikelihoods          = new TFile("likelihoods_bkgd.root","READ");
  TF2* fLikelihoodWTopHadrMass = (TF2*)fLikelihoods->Get("fLikelihoodWTopHadrMass");
  TF1* fLikelihoodCsvB         = (TF1*)fLikelihoods->Get("fLikelihoodCsvB");
  TF1* fLikelihoodCsvNonB      = (TF1*)fLikelihoods->Get("fLikelihoodCsvNonB");
  TF1* fLikelihoodTopLeptMass  = (TF1*)fLikelihoods->Get("fLikelihoodTopLeptMass");

  TFile* fHistos               = new TFile("histos_bkgd.root","READ");
  TH2F* hWTopHadrMass          = (TH2F*)fHistos->Get("hWTopHadrMass");
  hWHadrMass = (TH1F*)hWTopHadrMass->ProjectionX("_px");

  for(int k = 0; k < 10000 ; k++){
    double csv = ran->Uniform(-1E-03,1+1E-03);
    if(fLikelihoodCsvB->Eval(csv)>0)    hNonBLogLikeS->Fill( -2*TMath::Log( fLikelihoodCsvNonB->Eval(csv)/fLikelihoodCsvB->Eval(csv) ), fLikelihoodCsvNonB->Eval(csv));
    if(fLikelihoodCsvNonB->Eval(csv)>0) hNonBLogLikeB->Fill( -2*TMath::Log( fLikelihoodCsvNonB->Eval(csv)/fLikelihoodCsvB->Eval(csv) ), fLikelihoodCsvB->Eval(csv));
  }
  hNonBLogLikeS->Scale(1./hNonBLogLikeS->Integral());
  hNonBLogLikeB->Scale(1./hNonBLogLikeB->Integral());

  //fOutput->cd();
  //hWHadrMass->Write("hWTopHadrMass_px");
  //fOutput->Write();
  //fOutput->Close();
  //return;


  //TFile* fInput = TFile::Open("gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_VType2/v2/DiJetPt_TTH_HToBB_M-125_8TeV-pythia6_VType2.root","READ");
  TFile* fInput = TFile::Open("gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_VType2/v2/DiJetPt_TTJets_SemiLeptMGDecays_8TeV-madgraph-part_VType2.root","READ");
  TTree* tree = (TTree*)fInput->Get("tree");

  if(!tree){
    cout << "No tree" << endl; return;
  }

  JetByPt jet1, jet2, jet3, jet4, jet5, jet6, jet7, jet8, jet9, jet10;
  genTopInfo genTop, genTbar;
  genParticleInfo genB, genBbar;
  int nvlep;
  float vLeptonpt[999]; float vLeptoneta[999]; float vLeptonphi[999];
  float met[999];

  tree->SetBranchAddress("jet1",   &jet1);
  tree->SetBranchAddress("jet2",   &jet2);
  tree->SetBranchAddress("jet3",   &jet3);
  tree->SetBranchAddress("jet4",   &jet4);
  tree->SetBranchAddress("jet5",   &jet5);
  tree->SetBranchAddress("jet6",   &jet6);
  tree->SetBranchAddress("jet7",   &jet7);
  tree->SetBranchAddress("jet8",   &jet8);
  tree->SetBranchAddress("jet9",   &jet9);
  tree->SetBranchAddress("jet10",  &jet10);
  tree->SetBranchAddress("genTop", &genTop);
  tree->SetBranchAddress("genTbar",&genTbar);
  tree->SetBranchAddress("genB",   &genB);
  tree->SetBranchAddress("genBbar",&genBbar);
  tree->SetBranchAddress("vLepton_pt", vLeptonpt);
  tree->SetBranchAddress("vLepton_phi",vLeptonphi);
  tree->SetBranchAddress("vLepton_eta",vLeptoneta);
  tree->SetBranchAddress("nvlep",  &nvlep);
  tree->SetBranchAddress("METtype1p2corr",  met);

  Long64_t nentries = tree->GetEntries();
  int totCounter = 0;

   for (Long64_t i = 0; i < nentries; i++){
     
     if(i%10000==0) cout << i << endl;
     tree->GetEntry(i);

     if(i>1000000) break;

     if(nvlep<=0) continue;


     LV topBLV(  0.,0.,0.,0.);
     LV topW1LV( 0.,0.,0.,0.); 
     LV topW2LV( 0.,0.,0.,0.);
     LV atopBLV( 0.,0.,0.,0.); 
     LV atopW1LV(0.,0.,0.,0.); 
     LV atopW2LV(0.,0.,0.,0.);
     LV genBLV(0.,0.,0.,0.);
     LV genBbarLV(0.,0.,0.,0.);
    
    
     ///////////////////////////////////////////////////////////////////
     /////////////////////////////////////////////////////////////////// 
     if(genTop.bmass>0){
       topBLV.SetPt(   genTop.bpt );
       topBLV.SetEta(  genTop.beta );
       topBLV.SetPhi(  genTop.bphi );
       topBLV.SetM(    genTop.bmass );
       topW1LV.SetPt(  genTop.wdau1pt );
       topW1LV.SetEta( genTop.wdau1eta );
       topW1LV.SetPhi( genTop.wdau1phi );
       topW1LV.SetM(   genTop.wdau1mass );
       topW2LV.SetPt(  genTop.wdau2pt );
       topW2LV.SetEta( genTop.wdau2eta );
       topW2LV.SetPhi( genTop.wdau2phi );
       topW2LV.SetM(   genTop.wdau2mass );
     }
     if(genTbar.bmass>0){
	atopBLV.SetPt(   genTbar.bpt );
	atopBLV.SetEta(  genTbar.beta );
	atopBLV.SetPhi(  genTbar.bphi );
	atopBLV.SetM(    genTbar.bmass );
	atopW1LV.SetPt(  genTbar.wdau1pt );
	atopW1LV.SetEta( genTbar.wdau1eta );
	atopW1LV.SetPhi( genTbar.wdau1phi );
	atopW1LV.SetM(   genTbar.wdau1mass );
	atopW2LV.SetPt(  genTbar.wdau2pt );
	atopW2LV.SetEta( genTbar.wdau2eta );
	atopW2LV.SetPhi( genTbar.wdau2phi );
	atopW2LV.SetM(   genTbar.wdau2mass );
     }
     if(genB.mass>0 && genB.momid==25){
       genBLV.SetPt(  genB.pt );
       genBLV.SetEta( genB.eta );
       genBLV.SetPhi( genB.phi );
       genBLV.SetM(   genB.mass );
     }
     if(genBbar.mass>0 && genBbar.momid==25){
       genBbarLV.SetPt(  genBbar.pt );
       genBbarLV.SetEta( genBbar.eta );
       genBbarLV.SetPhi( genBbar.phi );
       genBbarLV.SetM(   genBbar.mass );
     }
     ///////////////////////////////////////////////////////////////////
     ///////////////////////////////////////////////////////////////////

     vector<LV> myJets;
     std::map<unsigned int, JetByPt> map;

     map[0] = jet1;  map[1] = jet2; map[2] = jet3; map[3] = jet4; 
     map[4] = jet5;  map[5] = jet6; map[6] = jet7; map[7] = jet8; 
     map[8] = jet9;  map[9] = jet10; 

     myJets.push_back(LV(jet1.pt, jet1.eta, jet1.phi, jet1.mass));
     myJets.push_back(LV(jet2.pt, jet2.eta, jet2.phi, jet2.mass));
     myJets.push_back(LV(jet3.pt, jet3.eta, jet3.phi, jet3.mass));
     myJets.push_back(LV(jet4.pt, jet4.eta, jet4.phi, jet4.mass));
     myJets.push_back(LV(jet5.pt, jet5.eta, jet5.phi, jet5.mass));
     myJets.push_back(LV(jet6.pt, jet6.eta, jet6.phi, jet6.mass));
     myJets.push_back(LV(jet7.pt, jet7.eta, jet7.phi, jet7.mass));
     myJets.push_back(LV(jet8.pt, jet8.eta, jet8.phi, jet8.mass));
     myJets.push_back(LV(jet9.pt, jet9.eta, jet9.phi, jet9.mass));
     myJets.push_back(LV(jet10.pt,jet10.eta,jet10.phi,jet10.mass));

     
     vector<LV> myJetsFilt;
     std::map<unsigned int, JetByPt> mapFilt;

     unsigned int iter = 0;
     for(unsigned int k = 0; k<myJets.size() ; k++){

       if(myJets[k].Pt()>ptCut_ && TMath::Abs(myJets[k].Eta()) < 2.5 ){

	 myJetsFilt.push_back( myJets[k] );

	 switch(k){
	 case 0:
	   mapFilt[iter] = jet1;
	   break;
	 case 1:
	   mapFilt[iter] = jet2;
	   break;
	 case 2:
	   mapFilt[iter] = jet3;
	   break;
	 case 3:
	   mapFilt[iter] = jet4;
	   break;
	 case 4:
	   mapFilt[iter] = jet5;
	   break;
	 case 5:
	   mapFilt[iter] = jet6;
	   break;
	 case 6:
	   mapFilt[iter] = jet7;
	   break;
	 case 7:
	   mapFilt[iter] = jet8;
	   break;
	 case 8:
	   mapFilt[iter] = jet9;
	   break;
	 case 9:
	   mapFilt[iter] = jet10;
	   break;
	 default:
	   break;
	 }
	 
	 iter++;
       }
     }

     LV lepLV(vLeptonpt[0], vLeptoneta[0], vLeptonphi[0], 0.);
     LV MEtLV(met[0], 0., met[3] , 0.);

     ///////////////////////////////////////////////////////////////////////////

     //if(  ( !((abs(genTop.wdau1id)<6  && genTop.wdau1pt>40  && TMath::Abs(genTop.wdau1eta) <2.5 && genTop.wdau2pt>40  && TMath::Abs(genTop.wdau2eta) <2.5)  || 
     // (abs(genTbar.wdau1id)<6 && genTbar.wdau1pt>40 && TMath::Abs(genTbar.wdau1eta)<2.5 && genTbar.wdau2pt>40 && TMath::Abs(genTbar.wdau2eta)<2.5))
     //    || !(genB.pt>30 && genBbar.pt>30) )
     //) continue;

     totCounter++;

     if( myJetsFilt.size()<2) continue;

     /////////////////////////////////////////////////////////////////////////// STEP2

     std::map<float,  int, sorterByLess> rankWHad;
     std::map<float,  pair<unsigned int,unsigned int>, sorterByLess > posWHad;

     for(unsigned int k = 0; k < myJetsFilt.size()-1; k++){  
       float csv_k = mapFilt[k].csv>0 ? mapFilt[k].csv : 0.0 ;

       for(unsigned int l = k+1; l < myJetsFilt.size(); l++){
	 float csv_l = mapFilt[l].csv>0 ? mapFilt[l].csv : 0.0 ;

	 LV diJet  = myJetsFilt[k]+myJetsFilt[l];
	 //LV diJetU = myJetsFilt[k]*(1+mapFilt[k].unc) + myJetsFilt[l]*(1+mapFilt[l].unc);
	 //LV diJetD = myJetsFilt[k]*(1-mapFilt[k].unc) + myJetsFilt[l]*(1-mapFilt[l].unc);
	 //float dM = TMath::Max( TMath::Abs(diJetU.M()-diJet.M()), TMath::Abs(diJetD.M()-diJet.M()) );
	 //dM = TMath::Sqrt(dM*dM + 2.08*2.08);
	 //TF1 *fGaus = new TF1("fGaus",Form("1/sqrt(2*TMath::Pi())/%f*exp(-0.5*((x-%f)/%f)**2)", dM, diJet.M(), dM),0,400);
	 //bool pass = (fGaus->Eval( 80.325 ) > fGaus->Eval( diJet.M()+1.96*dM ));
	 //bool pass = csv_l < 0.50 &&  csv_k < 0.50;
	 float likeB1   = fLikelihoodCsvNonB->Eval( csv_k );
	 float likeB2   = fLikelihoodCsvNonB->Eval( csv_l );
	 float likeWHad = hWHadrMass->GetBinContent( hWHadrMass->FindBin( diJet.M() ) );
	 float logLike  = likeB1*likeB2*likeWHad>0 ? -TMath::Log( likeB1*likeB2*likeWHad ) : 25 + ran->Uniform(0.,1);
	 //float logLike  = fGaus->Eval( 80.325 );

	 int genMatchk;
	 findGenMatch(genMatchk, myJetsFilt[k] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);
	 int genMatchl;
	 findGenMatch(genMatchl, myJetsFilt[l] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);

	 if( abs(genMatchk)==2 && abs(genMatchl)==2 )
	   hWHadLikeS->Fill(logLike);
	 else 
	   hWHadLikeB->Fill(logLike);

	 rankWHad[logLike] = int( abs(genMatchk)==2 && abs(genMatchl)==2 );
	 posWHad[logLike]  = make_pair(k,l);

	 //delete fGaus;	 
       }
     }

     int pos=0;
     for(std::map<float, int>::iterator it = rankWHad.begin(); it!= rankWHad.end() ; it++){
       if(it->second>0.5)
	 hWHadRank->Fill( pos );
       if(pos==0 && it->second>0.5)
	 hWHadLikeSRank0->Fill( it->first );
       if(pos==0 && it->second<0.5)
	 hWHadLikeBRank0->Fill( it->first );
       if(pos==1 && it->second>0.5)
	 hWHadLikeSRank1->Fill( it->first );
       if(pos==1 && it->second<0.5)
	 hWHadLikeBRank1->Fill( it->first );
       pos++;
     }

     if( rankWHad.size()<1 ) continue;
     if(  (rankWHad.begin())->first > cutW_ || (rankWHad.size()>1 && (++rankWHad.begin())->first < cutW2_) ) continue;

     unsigned int W1 = posWHad.size()>0 ? ((posWHad.begin())->second).first : 999;
     unsigned int W2 = posWHad.size()>0 ? ((posWHad.begin())->second).second: 999;

     if( W1==999 || W2==999 ) continue;

     LV w1LV = W1<999 ? myJetsFilt[W1] : LV(0,0,0,0);
     LV w2LV = W2<999 ? myJetsFilt[W2] : LV(0,0,0,0);

     int genMatchW1;
     findGenMatch(genMatchW1, w1LV , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);
     int genMatchW2;
     findGenMatch(genMatchW2, w2LV , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);
     if(abs(genMatchW1)==2 && abs(genMatchW2)==2 ) hCounterR->Fill(0.5);
     else  hCounterW->Fill(0.5);

     //continue;

     /////////////////////////////////////////////////////////////////////////// STEP2

     std::vector<unsigned int> unMasked;
     for(unsigned int k = 0; k < myJetsFilt.size(); k++)
       if( k!= W1 && k!=W2 && mapFilt[k].csv> csvCut_) unMasked.push_back( k );

     if( unMasked.size()<2 ) continue;
     if( abs(genMatchW1)==2 && abs(genMatchW2)==2 ) hCounterR->Fill(1.5);
     else  hCounterW->Fill(1.5);

     /////////////////////////////////////////////////////////////////////////// STEP3

     float a = lepLV.Px();
     float b = lepLV.Py();
     float c = lepLV.Pz();
     float f = lepLV.E();
     float x = MEtLV.Pt()*TMath::Cos(MEtLV.Phi());
     float y = MEtLV.Pt()*TMath::Sin(MEtLV.Phi());
     float m = 80.385;
     
     float delta = c*c*(-4*a*a - 4*b*b - 4*c*c + 4*f*f - 4*m*m - 8*a*x - 8*b*y)*(-4*a*a - 4*b*b - 4*c*c + 4*f*f - 4*m*m - 8*a*x - 8*b*y) - 4*(-4*c*c + 4*f*f)*(-a*a*a*a - 2*a*a*b*b - b*b*b*b - 2*a*a*c*c - 2*b*b*c*c - c*c*c*c + 2*a*a*f*f + 2*b*b*f*f + 2*c*c*f*f - f*f*f*f - 2*a*a*m*m - 2*b*b*m*m - 2*c*c*m*m + 2*f*f*m*m - m*m*m*m - 4*a*a*a*x - 4*a*b*b*x - 4*a*c*c*x + 4*a*f*f*x - 4*a*m*m*x - 4*a*a*x*x + 4*f*f*x*x - 4*a*a*b*y - 4*b*b*b*y - 4*b*c*c*y+ 4*b*f*f*y - 4*b*m*m*y - 8*a*b*x*y - 4*b*b*y*y + 4*f*f*y*y);
    
     if(delta<0) continue;
     if( abs(genMatchW1)==2 && abs(genMatchW2)==2 ) hCounterR->Fill(2.5);
     else  hCounterW->Fill(2.5);

     float z1 = (-c*(-4*a*a - 4*b*b - 4*c*c + 4*f*f - 4*m*m - 8*a*x - 8*b*y) - TMath::Sqrt(delta)) / (2*(-4*c*c + 4*f*f));
     float z2 = (-c*(-4*a*a - 4*b*b - 4*c*c + 4*f*f - 4*m*m - 8*a*x - 8*b*y) + TMath::Sqrt(delta)) / (2*(-4*c*c + 4*f*f));

     LV2 neutrinoLV2a(x,y,z1,TMath::Sqrt(x*x+y*y+z1*z1));
     LV neutrLV1(neutrinoLV2a);
     LV2 neutrinoLV2b(x,y,z2,TMath::Sqrt(x*x+y*y+z2*z2));
     LV neutrLV2(neutrinoLV2b);


     /////////////////////////////////////////////////////////////////////////// STEP4

     std::map<float,  pair<unsigned int,unsigned int> , sorterByLess> posHiggs;
     vector<unsigned int> unMasked2;

     int counterTopHad = 0;
     int counterTopLep = 0;

     for(unsigned int k = 0; k < unMasked.size(); k++){

       bool mask = false;

       unsigned int j1 = unMasked[k];
       LV j1LV = myJetsFilt[j1];
       
       float likeTopHad = fLikelihoodWTopHadrMass->Eval( (w1LV+w2LV).M(), (w1LV+w2LV+j1LV).M());
       float likeTopLep = TMath::Max(fLikelihoodTopLeptMass->Eval( (lepLV+neutrLV1+j1LV).M()),
				     fLikelihoodTopLeptMass->Eval( (lepLV+neutrLV2+j1LV).M())  );

       if( likeTopHad>= likeTopHadCut2_ ) counterTopHad++;
       if( likeTopLep>= likeTopLepCut2_ ) counterTopLep++;

       ///cout << likeTopHad << " - " << likeTopLep << endl;
       if( likeTopHad > likeTopHadCut_ || likeTopLep > likeTopLepCut_ ){
	 mask = true;
       }
       
       if(!mask) unMasked2.push_back( j1 );
	 
     }

     if(unMasked2.size()<2) continue;
     if(counterTopHad<1 ||  counterTopLep<1 ) continue;

     if( abs(genMatchW1)==2 && abs(genMatchW2)==2 ) hCounterR->Fill(3.5);
     else  hCounterW->Fill(3.5);

     for(unsigned int k = 0; k < unMasked2.size()-1; k++){
       for(unsigned int l = k+1; l < unMasked2.size(); l++){
	 unsigned int j1 = unMasked2[k];
	 unsigned int j2 = unMasked2[l];

	 LV h1LV = myJetsFilt[j1];
	 LV h2LV = myJetsFilt[j2];

	 //posHiggs[ (h1LV+h2LV).M() / (h1LV+h2LV).Pt()  ] = make_pair(j1,j2);
	 posHiggs[ 1./( mapFilt[j1].csv + mapFilt[j2].csv)   ] = make_pair(j1,j2);
	 //posHiggs[ TMath::Abs((h1LV+h2LV).M() - 120.)   ] = make_pair(j1,j2);
       }
     }
     
  
     unsigned int H1 = ((posHiggs.begin())->second).first;
     unsigned int H2 = ((posHiggs.begin())->second).second;

     int genMatchH1;
     findGenMatch(genMatchH1, myJetsFilt[H1] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);
     int genMatchH2;
     findGenMatch(genMatchH2, myJetsFilt[H2] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);

     if(  genMatchH1==0 && genMatchH2==0) hCounterR->Fill(4.5);
     else  hCounterW->Fill(4.5);

     //if(genMatchH1==0 && genMatchH2==0)
     hMassJJ->Fill( (myJetsFilt[H1] + myJetsFilt[H2]).M() );


   } // end LOOP


   hCounterW->Scale(1./totCounter);
   hCounterR->Scale(1./totCounter);

   cout << "Step 1 " << endl;
   cout << "Wrong = " << hCounterW->GetBinContent(1) << endl;
   cout << "Right = " << hCounterR->GetBinContent(1) << endl;
   
   cout << "Step 2 " << endl;
   cout << "Wrong = " << hCounterW->GetBinContent(2) << endl;
   cout << "Right = " << hCounterR->GetBinContent(2) << endl;
   
   cout << "Step 3 " << endl;
   cout << "Wrong = " << hCounterW->GetBinContent(3) << endl;
   cout << "Right = " << hCounterR->GetBinContent(3) << endl;
   
   cout << "Step 4 " << endl;
   cout << "Wrong = " << hCounterW->GetBinContent(4) << endl;
   cout << "Right = " << hCounterR->GetBinContent(4) << endl;

   cout << "Step 5 " << endl;
   cout << "Wrong = " << hCounterW->GetBinContent(5) << endl;
   cout << "Right = " << hCounterR->GetBinContent(5) << " (entries = " << hCounterR->GetBinContent(5)*totCounter  << ")" << endl;


   hWHadLikeS->Scale(1./hWHadLikeS->Integral());
   hWHadLikeB->Scale(1./hWHadLikeB->Integral());

   hWHadLikeSRank0->Scale(1./hWHadLikeSRank0->Integral());
   hWHadLikeBRank0->Scale(1./hWHadLikeBRank0->Integral());
   hWHadLikeSRank1->Scale(1./hWHadLikeSRank1->Integral());
   hWHadLikeBRank1->Scale(1./hWHadLikeBRank1->Integral());

   float res1 = (hCounterR->GetBinContent(5)/hCounterW->GetBinContent(5));
   float res2 =  hCounterR->GetBinContent(5);

   fOutput->Write();
   fOutput->Close();

   fLikelihoods->Close();
   fInput->Close();

   delete ran; delete fOutput; delete fLikelihoods;

   return ;

}



void Acc(int doGen=1, int jets30=6 ,  int bTagJets30 = 4, int excl=0){

  TFile *fOutput       = new TFile("Acc.root","RECREATE");
  TH1F* hCount         = new TH1F("hCount", "visited"   ,10, 0 ,10 );
  TH1F* hMult          = new TH1F("hMult", "visited"   ,10, 0 ,10 );
  TH2F* h2Mult         = new TH2F("h2Mult", "visited"   ,10, 0 ,10 ,10, 0 ,10 );
  TH1F* hMassWThere     = new TH1F("hMassWThere", "visited"   ,40,0, 400 );
  TH1F* hMassWLost     = new TH1F("hMassWLost", "visited"   ,40,0, 400 );

  TH1F* hCsvW     = new TH1F("hCsvW", "visited"   ,40,0, 1 );

  TH2F* hCsvWvsId     = new TH2F("hCsvWvsId", "visited"   ,40,0, 1 , 23, 0, 23);

  TH1F* hCsvB     = new TH1F("hCsvB", "visited"   ,40,0, 1 );
  TH1F* hCsvH     = new TH1F("hCsvH", "visited"   ,40,0, 1 );
  

  string relation = excl ? "=" : ">=" ;
  TH1F* hDistr = new TH1F("hDistr", Form("Njets%s%d, Nbjets%s%d ; ;", relation.c_str(), jets30, relation.c_str(),bTagJets30  )   ,20, 0 ,20 );

  TH1F* h46MatchBTag30   = new TH1F("h46MatchBTag30", "visited"   ,10, 0 ,10 );
  TH1F* h46MatchBTag20   = new TH1F("h46MatchBTag20", "visited"   ,10, 0 ,10 );
  TH1F* h46MatchNoTag30  = new TH1F("h46MatchNoTag30", "visited"   ,10, 0 ,10 );
  TH1F* h46MatchNoTag20  = new TH1F("h46MatchNoTag20", "visited"   ,10, 0 ,10 );

  TH1F* h45MatchBTag30   = new TH1F("h45MatchBTag30", "visited"   ,10, 0 ,10 );
  TH1F* h45MatchBTag20   = new TH1F("h45MatchBTag20", "visited"   ,10, 0 ,10 );
  TH1F* h45MatchNoTag30  = new TH1F("h45MatchNoTag30", "visited"   ,10, 0 ,10 );
  TH1F* h45MatchNoTag20  = new TH1F("h45MatchNoTag20", "visited"   ,10, 0 ,10 );

  TH1F* h35MatchBTag30   = new TH1F("h35MatchBTag30", "visited"   ,10, 0 ,10 );
  TH1F* h35MatchBTag20   = new TH1F("h35MatchBTag20", "visited"   ,10, 0 ,10 );
  TH1F* h35MatchNoTag30  = new TH1F("h35MatchNoTag30", "visited"   ,10, 0 ,10 );
  TH1F* h35MatchNoTag20  = new TH1F("h35MatchNoTag20", "visited"   ,10, 0 ,10 );

  TH1F* h36MatchBTag30   = new TH1F("h36MatchBTag30", "visited"   ,10, 0 ,10 );
  TH1F* h36MatchBTag20   = new TH1F("h36MatchBTag20", "visited"   ,10, 0 ,10 );
  TH1F* h36MatchNoTag30  = new TH1F("h36MatchNoTag30", "visited"   ,10, 0 ,10 );
  TH1F* h36MatchNoTag20  = new TH1F("h36MatchNoTag20", "visited"   ,10, 0 ,10 );



  TH1F* hcosThetaStarTopB =  new TH1F("hcosThetaStarTopB", ""     ,50, -1.1 ,1.1 );
  TH1F* hcosThetaStarTopW1 =  new TH1F("hcosThetaStarTopW1", ""   ,50, -1.1 ,1.1 );
  TH1F* hcosThetaStarTopW2 =  new TH1F("hcosThetaStarTopW2", ""   ,50, -1.1 ,1.1 );
  
  TH1F* hcosChiStarTopB     =  new TH1F("hcosChiStarTopB", ""      ,50, -1.1 ,1.1 );
  TH1F* hcosChiStarATopB    =  new TH1F("hcosChiStarATopB", ""     ,50, -1.1 ,1.1 );

  TH1F* hcosChiStarTopHadB     =  new TH1F("hcosChiStarTopHadB", ""      ,50, -1.1 ,1.1 );
  TH1F* hcosChiStarATopHadB    =  new TH1F("hcosChiStarATopHadB", ""     ,50, -1.1 ,1.1 );

  TH1F* hcosThetaStarHelTopW1 =  new TH1F("hcosThetaStarHelTopW1", ""   ,50, -1.1 ,1.1 );
  TH1F* hcosThetaStarHelTopW2 =  new TH1F("hcosThetaStarHelTopW2", ""   ,50, -1.1 ,1.1 );


  TH1F* hcosThetaStarATopB =  new TH1F("hcosThetaStarATopB", ""     ,50, -1.1 ,1.1 );
  TH1F* hcosThetaStarATopW1 =  new TH1F("hcosThetaStarATopW1", ""   ,50, -1.1 ,1.1 );
  TH1F* hcosThetaStarATopW2 =  new TH1F("hcosThetaStarATopW2", ""   ,50, -1.1 ,1.1 );

  TH1F* hcosThetaStarHelATopW1 =  new TH1F("hcosThetaStarHelATopW1", ""   ,50, -1.1 ,1.1 );
  TH1F* hcosThetaStarHelATopW2 =  new TH1F("hcosThetaStarHelATopW2", ""   ,50, -1.1 ,1.1 );

  TH1F* hcosThetaStarHiggsB    =  new TH1F("hcosThetaStarHiggsB", ""   ,50, -1.1 ,1.1 );
  TH1F* hcosThetaStarHiggsBbar =  new TH1F("hcosThetaStarHiggsBbar", ""   ,50, -1.1 ,1.1 );

  TH1F* hcosThetaStarWrong    =  new TH1F("hcosThetaStarWrong", ""   ,50, -1.1 ,1.1 );

  TH1F* hPt =  new TH1F("hPt", ""   ,60, 0 ,6 );

  TH2F* hPtMiss =  new TH2F("hPtMiss", ""   ,120, 0 ,600, 100,-5,5 );
  TH2F* hPtTrue =  new TH2F("hPtTrue", ""   ,120, 0 ,600, 100,-5,5 );
  TH1F* hcosMiss1 =  new TH1F("hcosMiss1", ""     ,50, -1.1 ,1.1 );
  TH1F* hcosMiss2 =  new TH1F("hcosMiss2", ""     ,50, -1.1 ,1.1 );
  TH1F* hcosMiss3 =  new TH1F("hcosMiss3", ""     ,50, -1.1 ,1.1 );

  //TH1F* h =  new TH1F("", ""   ,50, -1.1 ,1.1 );

  TFile* fInput = TFile::Open("gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_VType2/v2/DiJetPt_TTH_HToBB_M-125_8TeV-pythia6_VType2.root","READ");
  //TFile* fInput = TFile::Open("gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_VType2/v2/DiJetPt_TTJets_SemiLeptMGDecays_8TeV-madgraph-part_VType2.root","READ");
  TTree* tree = (TTree*)fInput->Get("tree");

  if(!tree){
    cout << "No tree" << endl; return ;
  }

  TH1F* hCounter = (TH1F*)fInput->Get("Count");
  //float weight = (0.1302*0.569)*12.1*1000/hCounter->GetBinContent(1);
  //cout << weight << endl;

  float weight = (103.0)*12.1*1000/hCounter->GetBinContent(1);
  //weight /= (500000./tree->GetEntries());
  cout << weight << endl;

  JetByPt jet1, jet2, jet3, jet4, jet5, jet6, jet7, jet8, jet9, jet10;
  genTopInfo genTop, genTbar;
  genParticleInfo genB, genBbar;
  recoTopInfo recoTop, recoTbar;

  int nvlep;
  float vLeptonpt[999]; float vLeptoneta[999]; float vLeptonphi[999];
  float met[999];

  int nSvs;
  float Sv_massBCand[999];
  float Sv_massSv[999];
  float Sv_pt[999];
  float Sv_eta[999];
  float Sv_phi[999];
  float Sv_dist3D[999];
  float Sv_dist2D[999];
  float Sv_distSim2D[999];
  float Sv_distSig3D[999];
  float Sv_dist3D_norm[999];

  
  tree->SetBranchAddress("jet1",   &jet1);
  tree->SetBranchAddress("jet2",   &jet2);
  tree->SetBranchAddress("jet3",   &jet3);
  tree->SetBranchAddress("jet4",   &jet4);
  tree->SetBranchAddress("jet5",   &jet5);
  tree->SetBranchAddress("jet6",   &jet6);
  tree->SetBranchAddress("jet7",   &jet7);
  tree->SetBranchAddress("jet8",   &jet8);
  tree->SetBranchAddress("jet9",   &jet9);
  tree->SetBranchAddress("jet10",  &jet10);
  tree->SetBranchAddress("genTop", &genTop);
  tree->SetBranchAddress("genTbar",&genTbar);
  tree->SetBranchAddress("recoTop", &recoTop);
  tree->SetBranchAddress("recoTbar",&recoTbar);


  tree->SetBranchAddress("genB",   &genB);
  tree->SetBranchAddress("genBbar",&genBbar);
  tree->SetBranchAddress("vLepton_pt", vLeptonpt);
  tree->SetBranchAddress("vLepton_phi",vLeptonphi);
  tree->SetBranchAddress("vLepton_eta",vLeptoneta);
  tree->SetBranchAddress("nvlep",  &nvlep);
  tree->SetBranchAddress("METtype1p2corr",  met);

  tree->SetBranchAddress("nSvs",&nSvs);
  tree->SetBranchAddress("Sv_massBCand", Sv_massBCand);
  tree->SetBranchAddress("Sv_massSv", Sv_massSv);
  tree->SetBranchAddress("Sv_pt", Sv_pt);
  tree->SetBranchAddress("Sv_eta", Sv_eta);
  tree->SetBranchAddress("Sv_phi", Sv_phi);
  tree->SetBranchAddress("Sv_dist3D", Sv_dist3D);
  tree->SetBranchAddress("Sv_dist2D",Sv_dist2D);
  tree->SetBranchAddress("Sv_distSim2D", Sv_distSim2D);
  tree->SetBranchAddress("Sv_distSig3D", Sv_distSig3D);
  tree->SetBranchAddress("Sv_dist3D_norm", Sv_dist3D_norm);


  Long64_t nentries = tree->GetEntries();
  //nentries = 500000;
  int totCounter = 0;
  int doneTest = 0;

  
  for (Long64_t i = 0; i < nentries ; i++){
    
    if(i%10000==0) cout << i << endl;
    tree->GetEntry(i);
    
    //continue;

    if(nvlep<=0) continue;

    LV topBLV(  0.,0.,0.,0.);
    LV topW1LV( 0.,0.,0.,0.); 
    LV topW2LV( 0.,0.,0.,0.);
    LV atopBLV( 0.,0.,0.,0.); 
    LV atopW1LV(0.,0.,0.,0.); 
    LV atopW2LV(0.,0.,0.,0.);
    LV genBLV(0.,0.,0.,0.);
    LV genBbarLV(0.,0.,0.,0.);
    
    
     ///////////////////////////////////////////////////////////////////
     /////////////////////////////////////////////////////////////////// 
    if(doGen){

      if(genTop.bmass>0){
	topBLV.SetPt(   genTop.bpt );
	topBLV.SetEta(  genTop.beta );
	topBLV.SetPhi(  genTop.bphi );
	topBLV.SetM(    genTop.bmass );
	topW1LV.SetPt(  genTop.wdau1pt );
	topW1LV.SetEta( genTop.wdau1eta );
	topW1LV.SetPhi( genTop.wdau1phi );
	topW1LV.SetM(   genTop.wdau1mass );
	topW2LV.SetPt(  genTop.wdau2pt );
	topW2LV.SetEta( genTop.wdau2eta );
	topW2LV.SetPhi( genTop.wdau2phi );
	topW2LV.SetM(   genTop.wdau2mass );
      }
      if(genTbar.bmass>0){
	atopBLV.SetPt(   genTbar.bpt );
	atopBLV.SetEta(  genTbar.beta );
	atopBLV.SetPhi(  genTbar.bphi );
	atopBLV.SetM(    genTbar.bmass );
	atopW1LV.SetPt(  genTbar.wdau1pt );
	atopW1LV.SetEta( genTbar.wdau1eta );
	atopW1LV.SetPhi( genTbar.wdau1phi );
	atopW1LV.SetM(   genTbar.wdau1mass );
	atopW2LV.SetPt(  genTbar.wdau2pt );
	atopW2LV.SetEta( genTbar.wdau2eta );
	atopW2LV.SetPhi( genTbar.wdau2phi );
	atopW2LV.SetM(   genTbar.wdau2mass );
     }

    }
    else{

      if(recoTop.bmass>0 && (  (
				(abs(genTop.wdau1id)<6 && abs(genTop.wdau2id)<6 && recoTop.topid==1 && recoTop.wid==1) || 
				(abs(genTop.wdau1id)>=11 && abs(genTop.wdau2id)>=11 && recoTop.bstatus==1)
				) 
			       ) ){
       topBLV.SetPt(   recoTop.bpt );
       topBLV.SetEta(  recoTop.beta );
       topBLV.SetPhi(  recoTop.bphi );
       topBLV.SetM(    recoTop.bmass );
       if(abs(genTop.wdau1id)<6 && abs(genTop.wdau2id)<6 && recoTop.topid==1 && recoTop.wid==1){
	 topW1LV.SetPt(  recoTop.wdau1pt );
	 topW1LV.SetEta( recoTop.wdau1eta );
	 topW1LV.SetPhi( recoTop.wdau1phi );
	 topW1LV.SetM(   recoTop.wdau1mass );
	 topW2LV.SetPt(  recoTop.wdau2pt );
	 topW2LV.SetEta( recoTop.wdau2eta );
	 topW2LV.SetPhi( recoTop.wdau2phi );
	 topW2LV.SetM(   recoTop.wdau2mass );
       }
       else{
	 topW1LV.SetPt(  genTop.wdau1pt );
	 topW1LV.SetEta( genTop.wdau1eta );
	 topW1LV.SetPhi( genTop.wdau1phi );
	 topW1LV.SetM(   genTop.wdau1mass );
	 topW2LV.SetPt(  genTop.wdau2pt );
	 topW2LV.SetEta( genTop.wdau2eta );
	 topW2LV.SetPhi( genTop.wdau2phi );
	 topW2LV.SetM(   genTop.wdau2mass );
       }
      }

      if(recoTbar.bmass>0 && (  (
				 (abs(genTbar.wdau1id)<6 && abs(genTbar.wdau2id)<6 && recoTbar.topid==1 && recoTbar.wid==1) || 
				 (abs(genTbar.wdau1id)>=11 && abs(genTbar.wdau2id)>=11 && recoTbar.bstatus==1)
				 ) 
			       ) ){
	atopBLV.SetPt(   recoTbar.bpt );
	atopBLV.SetEta(  recoTbar.beta );
	atopBLV.SetPhi(  recoTbar.bphi );
	atopBLV.SetM(    recoTbar.bmass );
	if( (abs(genTbar.wdau1id)<6 && abs(genTbar.wdau2id)<6 && recoTbar.topid==1 && recoTbar.wid==1) ){
	  atopW1LV.SetPt(  recoTbar.wdau1pt );
	  atopW1LV.SetEta( recoTbar.wdau1eta );
	  atopW1LV.SetPhi( recoTbar.wdau1phi );
	  atopW1LV.SetM(   recoTbar.wdau1mass );
	  atopW2LV.SetPt(  recoTbar.wdau2pt );
	  atopW2LV.SetEta( recoTbar.wdau2eta );
	  atopW2LV.SetPhi( recoTbar.wdau2phi );
	  atopW2LV.SetM(   recoTbar.wdau2mass );
	}
	else{
	  atopW1LV.SetPt(  genTbar.wdau1pt );
	  atopW1LV.SetEta( genTbar.wdau1eta );
	  atopW1LV.SetPhi( genTbar.wdau1phi );
	  atopW1LV.SetM(   genTbar.wdau1mass );
	  atopW2LV.SetPt(  genTbar.wdau2pt );
	  atopW2LV.SetEta( genTbar.wdau2eta );
	  atopW2LV.SetPhi( genTbar.wdau2phi );
	  atopW2LV.SetM(   genTbar.wdau2mass );
	}

     }
    }

    if(genB.mass>0 && genB.momid==25){
       genBLV.SetPt(  genB.pt );
       genBLV.SetEta( genB.eta );
       genBLV.SetPhi( genB.phi );
       genBLV.SetM(   genB.mass );
     }
     if(genBbar.mass>0 && genBbar.momid==25){
       genBbarLV.SetPt(  genBbar.pt );
       genBbarLV.SetEta( genBbar.eta );
       genBbarLV.SetPhi( genBbar.phi );
       genBbarLV.SetM(   genBbar.mass );
     }

     //cout <<  topBLV.Pt() << "," << topW1LV.Pt() << "," << topW2LV.Pt() << "," << atopBLV.Pt() <<  "," << atopW1LV.Pt() << endl;

     if(!(topBLV.Pt()>0 &&
	  topW1LV.Pt()>0 &&
	  topW2LV.Pt()>0 &&
	  atopBLV.Pt()>0 &&
	  atopW1LV.Pt()>0 &&
	  atopW2LV.Pt()>0 
	  //&& genBLV.Pt()>0 &&
	  //genBbarLV.Pt()>0 
	  && (abs(genTop.wdau1id)<6 || abs(genTbar.wdau1id)<6 )
	  )
	) continue;
     
     vector<LV> myJets;
     std::map<unsigned int, JetByPt> map;

     map[0] = jet1;  map[1] = jet2; map[2] = jet3; map[3] = jet4; 
     map[4] = jet5;  map[5] = jet6; map[6] = jet7; map[7] = jet8; 
     map[8] = jet9;  map[9] = jet10; 

     myJets.push_back(LV(jet1.pt, jet1.eta, jet1.phi, jet1.mass));
     myJets.push_back(LV(jet2.pt, jet2.eta, jet2.phi, jet2.mass));
     myJets.push_back(LV(jet3.pt, jet3.eta, jet3.phi, jet3.mass));
     myJets.push_back(LV(jet4.pt, jet4.eta, jet4.phi, jet4.mass));
     myJets.push_back(LV(jet5.pt, jet5.eta, jet5.phi, jet5.mass));
     myJets.push_back(LV(jet6.pt, jet6.eta, jet6.phi, jet6.mass));
     myJets.push_back(LV(jet7.pt, jet7.eta, jet7.phi, jet7.mass));
     myJets.push_back(LV(jet8.pt, jet8.eta, jet8.phi, jet8.mass));
     myJets.push_back(LV(jet9.pt, jet9.eta, jet9.phi, jet9.mass));
     myJets.push_back(LV(jet10.pt,jet10.eta,jet10.phi,jet10.mass));

     
     vector<LV> myJetsFilt;
     std::map<unsigned int, JetByPt> mapFilt;
     unsigned int iter = 0;

     for(unsigned int k = 0; k<myJets.size() ; k++){

       if(myJets[k].Pt()>20 && TMath::Abs(myJets[k].Eta()) < 2.5 ){

	 myJetsFilt.push_back( myJets[k] );

	 switch(k){
	 case 0:
	   mapFilt[iter] = jet1;
	   break;
	 case 1:
	   mapFilt[iter] = jet2;
	   break;
	 case 2:
	   mapFilt[iter] = jet3;
	   break;
	 case 3:
	   mapFilt[iter] = jet4;
	   break;
	 case 4:
	   mapFilt[iter] = jet5;
	   break;
	 case 5:
	   mapFilt[iter] = jet6;
	   break;
	 case 6:
	   mapFilt[iter] = jet7;
	   break;
	 case 7:
	   mapFilt[iter] = jet8;
	   break;
	 case 8:
	   mapFilt[iter] = jet9;
	   break;
	 case 9:
	   mapFilt[iter] = jet10;
	   break;
	 default:
	   break;
	 }
	 
	 iter++;
       }
     }

     LV lepLV(vLeptonpt[0], vLeptoneta[0], vLeptonphi[0], 0.);
     LV MEtLV(met[0], 0., met[3] , 0.);

     int numBTag=0;
     for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
       float csv_k = mapFilt[k].csv>0 ? mapFilt[k].csv : 0.0 ;
       float pt_k = mapFilt[k].pt>0 ? mapFilt[k].pt : 0.0 ;
       if(csv_k>0.679 && pt_k>30) numBTag++;
     }
     int numJets30=0;
     for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
       float pt_k = mapFilt[k].pt>0 ? mapFilt[k].pt : 0.0 ;
       if(pt_k>30) numJets30++;
     }
     
     hMult->Fill(0.5, weight);
     h2Mult->Fill(numJets30, numBTag, weight);

     if(numBTag==3 && numJets30==5)
       hMult->Fill(1.5, weight);
     if(numBTag==3 && numJets30>=6)
       hMult->Fill(2.5, weight);
     if(numBTag==4 && numJets30==5)
       hMult->Fill(3.5, weight);
     if(numBTag==4 && numJets30>=6)
       hMult->Fill(4.5, weight);

     if( (!excl && ( (jets30>=6 && numJets30>=0/*jets30*/) || (jets30<6 && numJets30>=jets30)) && numBTag>=bTagJets30) ||
	 (excl  && ( (jets30>=6 && numJets30>=jets30) || (jets30<6 && numJets30==jets30)) && numBTag==bTagJets30) 
	 ){


       /////////////////////////////////////////////

       if(abs(genTop.wdau1id)<6 && topBLV.Pt()>40 && topW1LV.Pt()>40){

	 doneTest++;
	 if(doneTest==2){


	 float ptS  = 1.0;
	 float etaS = 0.1;
	 float phiS = 0.1;
	 
	 hPtTrue->Fill(topW2LV.Pt(), topW2LV.Eta());
	 cout << topW2LV.Pt() << ", " << topW2LV.Eta() << ", " << topW2LV.Phi() << endl;
	 cout << topW1LV.Pt() << ", " << topW1LV.Eta() << ", " << topW1LV.Phi() << endl;
	 cout << topBLV.Pt() << ", " <<topBLV.Eta() << ", " << topBLV.Phi() << endl;

	 for(int scan1 = 0; scan1<100    ; scan1++){
	   for(int scan2 = 0; scan2<100  ; scan2++){
	     for(int scan3 = 0; scan3<30 ; scan3++){
	       
	       float pt  =  5.0 + ptS*scan1;
	       float eta = -5  + etaS*scan2; 
	       float phi = -TMath::Pi() +  (2*TMath::Pi())/30.*scan3;
	       

	       LV jetLV(pt, eta, phi, 0.);
	       float TopMass =  (topBLV+topW1LV+jetLV).M();
	       float WMass   =  (topW1LV+jetLV).M();

	       LV topLV  = topBLV+topW1LV+jetLV;
	       LV topWLV = topW1LV+jetLV;
	       
	       TVector3 boost(topLV.Px()/topLV.E(), topLV.Py()/topLV.E(), topLV.Pz()/topLV.E());
	       
	       TLorentzVector bHad;
	       TLorentzVector w1Had;
	       TLorentzVector w2Had;
	       
	       bHad.SetPxPyPzE( topBLV.Px(), topBLV.Py(), topBLV.Pz(), topBLV.E());
	       w1Had.SetPxPyPzE( topW1LV.Px(), topW1LV.Py(), topW1LV.Pz(), topW1LV.E());
	       w2Had.SetPxPyPzE( jetLV.Px(), jetLV.Py(), jetLV.Pz(), jetLV.E());
	       
	       bHad.Boost(-boost);
	       w1Had.Boost(-boost);
	       w2Had.Boost(-boost);

	       double cosThetaStarTopB  = TMath::Cos(bHad.Vect().Angle(boost));
	       double cosThetaStarTopW1 = TMath::Cos(w1Had.Vect().Angle(boost));
	       double cosThetaStarTopW2 = TMath::Cos(w2Had.Vect().Angle(boost));

	       hcosMiss1->Fill(cosThetaStarTopB, weight);
	       hcosMiss2->Fill(cosThetaStarTopW1, weight);
	       hcosMiss3->Fill(cosThetaStarTopW2, weight);

	   
	       if( pt>55 && pt < 58 && eta>-2.3 && eta<-2. && phi>2.5 && phi<2.7){
		 cout << TopMass << ", " << WMass << endl;
	       }

	       if(TMath::Abs(WMass-81.)<15. && TMath::Abs(TopMass-175.)<20){
		 hPtMiss->Fill( pt , eta);
		 //doneTest++;
	       }
	       
	     }
	   }
	 }
	 }
       }
	  
       /////////////////////////////////////////////

       int matchW1 = 0;
       int matchW2 = 0;
       int matchB1 = 0;
       int matchB2 = 0;
       int matchH1 = 0;
       int matchH2 = 0;

       int twoBsMerg = 0;
       int twoWsMerg = 0;
       int oneBWMerg = 0;

       LV w1(0.,0.,0.,0.); unsigned int w1index=999;
       LV w2(0.,0.,0.,0.); unsigned int w2index=999;

       for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
	 float pt_k = mapFilt[k].pt>0 ? mapFilt[k].pt : 0.0 ;
	 float eta_k = mapFilt[k].eta;
	 if(pt_k<30 || TMath::Abs(eta_k)>2.5) continue;
	 float csv_k = mapFilt[k].csv>0 ? mapFilt[k].csv : 0.0 ;
	 float pdgid_k = abs(mapFilt[k].flavor);

	 if(csv_k<0.679 && w1.Pt()<0.001) w1 = myJetsFilt[k]; 
	 else if(csv_k<0.679 && w2.Pt()<0.001 && w1.Pt()>0) w2 = myJetsFilt[k]; 
	 else{}

	 int matchByTopB=0;
	 int matchByTopW=0;
	 int matchByHiggsB=0;

	 //cout << myJetsFilt[k].Pt() << " - " << topBLV.Pt() << endl;

	 int genMatch;
	 findGenMatch2(genMatch, myJetsFilt[k] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);
	 if( deltaR(myJetsFilt[k], topBLV)<0.3 ){
	   if(csv_k>0.679){
	     matchB1++;
	     hCsvB->Fill(csv_k);
	   }
	     matchByTopB++;
	 }
	 if( (abs(genTop.wdau1id)<6 && deltaR(myJetsFilt[k], topW1LV)<0.3) || 
	     (abs(genTbar.wdau1id)<6 && deltaR(myJetsFilt[k], atopW1LV)<0.3) ){
	   if(csv_k<0.679){
	     matchW1++;
	     w1index = k;
	     hCsvW->Fill(csv_k);
	     hCsvWvsId->Fill(csv_k,pdgid_k);
	   }
	   matchByTopW++;
	 }
	 if( (abs(genTop.wdau1id)<6 && deltaR(myJetsFilt[k], topW2LV)<0.3) || 
	     (abs(genTbar.wdau1id)<6 && deltaR(myJetsFilt[k], atopW2LV)<0.3) ){
	   if(csv_k<0.679){
	     matchW2++;
	     w2index = k;
	     hCsvW->Fill(csv_k);
	     hCsvWvsId->Fill(csv_k,pdgid_k);
	   }
	   matchByTopW++;
	 }
	 if( deltaR(myJetsFilt[k], atopBLV)<0.3 ){
	   if(csv_k>0.679){
	     matchB2++;
	     hCsvB->Fill(csv_k);
	   }
	   matchByTopB++;
	 }
	 if( deltaR(myJetsFilt[k], genBLV)<0.3 ){
	   if(csv_k>0.679){
	     matchH1++;
	     hCsvH->Fill(csv_k);
	   }
	   matchByHiggsB++;
	 }
	 if( deltaR(myJetsFilt[k], genBbarLV)<0.3 ){
	   if(csv_k>0.679){
	     matchH2++;
	     hCsvH->Fill(csv_k);
	   }
	   matchByHiggsB++;
	 }

	 if(matchByTopW>1) twoWsMerg=1;
	 if(matchByTopW==1 && (matchByTopB>0 || matchByHiggsB>0)) oneBWMerg=1;
	 if(matchByTopB>1 || matchByHiggsB>1 || (matchByTopB>0 && matchByHiggsB>0)) twoBsMerg=1;

       }

       int nomerge =  twoWsMerg==0 && twoBsMerg==0 && oneBWMerg==0;
       int merge1 = twoWsMerg>0;
       int merge2 = twoBsMerg>0;
       int merge3 = oneBWMerg>0;

       // if bkg
       //matchH1=0; matchH2=0;

       //if( !(w1index!=999 && w2index!=999 && ((myJetsFilt[w1index]+myJetsFilt[w2index]).M()<1000 && (myJetsFilt[w1index]+myJetsFilt[w2index]).M()>0)  )) continue;
       //if( !(abs(genTop.wdau1id)==4 ||  abs(genTop.wdau2id)==4 || abs(genTbar.wdau1id)==4 ||  abs(genTbar.wdau2id)==4) ) continue;
       //if( !(w1.Pt()>0 && w2.Pt()>0 &&  (w1+w2).M() > 60 &&  (w1+w2).M()<100) ) continue;

       // cout << matchW1 << ", " << matchW2 << ", " << matchB1 << ", " <<  matchB2 << ", " << matchH1 << ", " <<  matchH2 << endl;
       //cout << nomerge << ", " << merge1 << ", " << merge2 << ", "<< merge2 << endl;

       if( matchW1==1 && matchW2==1 && matchB1==1 && matchB2==1 && matchH1==1 && matchH2==1 && nomerge )/*&& (w1index!=999 && w2index!=999 && (myJetsFilt[w1index]+myJetsFilt[w2index]).M()<120 && (myJetsFilt[w1index]+myJetsFilt[w2index]).M()>30*/ {
	 hDistr->Fill(  0.5 ); // ideal
	 hMassWThere->Fill( (w1+w2).M() );
       }
       else if( ( (matchW1<1 && matchW2==1) || (matchW1==1 && matchW2<1)) && matchB1==1 && matchB2==1 && matchH1==1 && matchH2==1 && nomerge){
	 hDistr->Fill(  1.5 ); // no W
	 hMassWLost->Fill( (w1+w2).M() );
       }
       else if( matchW1==1 && matchW2==1 && (matchB1<1 || matchB2<1) && matchH1==1 && matchH2==1 && nomerge) 
	 hDistr->Fill(  2.5 ); // no top bs
       else if( matchW1==1 && matchW2==1 && matchB1==1 && matchB2==1 && (matchH1<1 || matchH2<1) && nomerge) 
	 hDistr->Fill(  3.5 ); // no higgs bs
       else if( matchW1==1 && matchW2==1 && matchB1==1 && matchB2==1 && matchH1==1 && matchH2==1 && merge1) 
	 hDistr->Fill(  4.5 ); // ideal
       else if( (matchW1<1 || matchW2<1) && matchB1==1 && matchB2==1 && matchH1==1 && matchH2==1 && merge1) 
	 hDistr->Fill(  5.5 ); // no W
       else if( matchW1==1 && matchW2==1 && (matchB1<1 || matchB2<1) && matchH1==1 && matchH2==1 && merge1) 
	 hDistr->Fill(  6.5 ); // no top bs
       else if( matchW1==1 && matchW2==1 && matchB1==1 && matchB2==1 && (matchH1<1 || matchH2<1) && merge1) 
	 hDistr->Fill(  7.5 ); // no higgs bs
       else if( matchW1==1 && matchW2==1 && matchB1==1 && matchB2==1 && matchH1==1 && matchH2==1 && merge2) 
	 hDistr->Fill(  8.5 ); // ideal
       else if( (matchW1<1 || matchW2<1) && matchB1==1 && matchB2==1 && matchH1==1 && matchH2==1 && merge2) 
	 hDistr->Fill(  9.5 ); // no W
       else if( matchW1==1 && matchW2==1 && (matchB1<1 || matchB2<1) && matchH1==1 && matchH2==1 && merge2) 
	 hDistr->Fill(  10.5 ); // no top bs
       else if( matchW1==1 && matchW2==1 && matchB1==1 && matchB2==1 && (matchH1<1 || matchH2<1) && merge2) 
	 hDistr->Fill(  11.5 ); // no higgs bs
       else if( matchW1==1 && matchW2==1 && matchB1==1 && matchB2==1 && matchH1==1 && matchH2==1 && merge3) 
	 hDistr->Fill(  12.5 ); // ideal
       else if( (matchW1<1 || matchW2<1) && matchB1==1 && matchB2==1 && matchH1==1 && matchH2==1 && merge3) 
	 hDistr->Fill(  13.5 ); // no W
       else if( matchW1==1 && matchW2==1 && (matchB1<1 || matchB2<1) && matchH1==1 && matchH2==1 && merge3) 
	 hDistr->Fill(  14.5 ); // no top bs
       else if( matchW1==1 && matchW2==1 && matchB1==1 && matchB2==1 && (matchH1<1 || matchH2<1) && merge3) 
	 hDistr->Fill(  15.5 ); // no higgs bs
       else if( (matchW1<1 || matchW2<1) && matchB1==1 && matchB2==1 && (matchH1<1 || matchH2<1)) 
	 hDistr->Fill(  16.5 ); // no W
       else if( (matchW1<1 || matchW2<1) && (matchB1<1 || matchB2<1) && matchH1==1 && matchH2==1) 
	 hDistr->Fill(  17.5 ); // no W
       else if( (matchW1<1 && matchW2<1) && matchB1==1 && matchB2==1 && matchH1==1 && matchH2==1 && nomerge) 
	 hDistr->Fill(  18.5 ); // no W
       else{
	 hDistr->Fill(  19.5 );

	 //cout << matchW1 << ", " << matchW2 << ", " << matchB1 << ", " <<  matchB2 << ", " << matchH1 << ", " <<  matchH2 << endl;
	 //cout << nomerge << ", " << merge1 << ", " << merge2 << ", "<< merge2 << endl;

       }

     }


     if(abs(genTop.wdau1id)<6 && topBLV.Pt()>0 && topW1LV.Pt()>0 && topW2LV.Pt()>0){
       LV topLV  = topBLV+topW1LV+topW2LV;
       LV topWLV = topW1LV+topW2LV;
       
       TVector3 boost(topLV.Px()/topLV.E(), topLV.Py()/topLV.E(), topLV.Pz()/topLV.E());
       
       TLorentzVector bHad;
       TLorentzVector w1Had;
       TLorentzVector w2Had;

       bHad.SetPxPyPzE( topBLV.Px(), topBLV.Py(), topBLV.Pz(), topBLV.E());
       w1Had.SetPxPyPzE( topW1LV.Px(), topW1LV.Py(), topW1LV.Pz(), topW1LV.E());
       w2Had.SetPxPyPzE( topW2LV.Px(), topW2LV.Py(), topW2LV.Pz(), topW2LV.E());

       bHad.Boost(-boost);
       w1Had.Boost(-boost);
       w2Had.Boost(-boost);

       double cosThetaStarTopB  = TMath::Cos(bHad.Vect().Angle(boost));
       double cosThetaStarTopW1 = TMath::Cos(w1Had.Vect().Angle(boost));
       double cosThetaStarTopW2 = TMath::Cos(w2Had.Vect().Angle(boost));

       hcosThetaStarTopB->Fill(cosThetaStarTopB, weight);
       hcosThetaStarTopW1->Fill(cosThetaStarTopW1, weight);
       hcosThetaStarTopW2->Fill(cosThetaStarTopW2, weight);

       TVector3 boostW(topWLV.Px()/topWLV.E(), topWLV.Py()/topWLV.E(), topWLV.Pz()/topWLV.E());
       w1Had.SetPxPyPzE( topW1LV.Px(), topW1LV.Py(), topW1LV.Pz(), topW1LV.E());
       w2Had.SetPxPyPzE( topW2LV.Px(), topW2LV.Py(), topW2LV.Pz(), topW2LV.E());
       w1Had.Boost(-boostW);
       w2Had.Boost(-boostW);
       double cosThetaStarHelTopW1 = TMath::Cos(w1Had.Vect().Angle(boostW));
       double cosThetaStarHelTopW2 = TMath::Cos(w2Had.Vect().Angle(boostW));
       hcosThetaStarHelTopW1->Fill(cosThetaStarHelTopW1, weight);
       hcosThetaStarHelTopW2->Fill(cosThetaStarHelTopW2, weight);
     }

     if(abs(genTbar.wdau1id)<6 && atopBLV.Pt()>0 && atopW1LV.Pt()>0 && atopW2LV.Pt()>0){
       LV topLV  = atopBLV+atopW1LV+atopW2LV;
       LV topWLV = atopW1LV+atopW2LV;

       TVector3 boost(topLV.Px()/topLV.E(), topLV.Py()/topLV.E(), topLV.Pz()/topLV.E());
       
       TLorentzVector bHad;
       TLorentzVector w1Had;
       TLorentzVector w2Had;

       bHad.SetPxPyPzE(  atopBLV.Px(),  atopBLV.Py(),  atopBLV.Pz(),  atopBLV.E());
       w1Had.SetPxPyPzE( atopW1LV.Px(), atopW1LV.Py(), atopW1LV.Pz(), atopW1LV.E());
       w2Had.SetPxPyPzE( atopW2LV.Px(), atopW2LV.Py(), atopW2LV.Pz(), atopW2LV.E());

       bHad.Boost( -boost);
       w1Had.Boost(-boost);
       w2Had.Boost(-boost);

       double cosThetaStarATopB  = TMath::Cos(bHad.Vect().Angle(boost));
       double cosThetaStarATopW1 = TMath::Cos(w1Had.Vect().Angle(boost));
       double cosThetaStarATopW2 = TMath::Cos(w2Had.Vect().Angle(boost));

       hcosThetaStarATopB->Fill(cosThetaStarATopB, weight);
       hcosThetaStarATopW1->Fill(cosThetaStarATopW1, weight);
       hcosThetaStarATopW2->Fill(cosThetaStarATopW2, weight);


       TVector3 boostW(topWLV.Px()/topWLV.E(), topWLV.Py()/topWLV.E(), topWLV.Pz()/topWLV.E());
       w1Had.SetPxPyPzE( atopW1LV.Px(), atopW1LV.Py(), atopW1LV.Pz(), atopW1LV.E());
       w2Had.SetPxPyPzE( atopW2LV.Px(), atopW2LV.Py(), atopW2LV.Pz(), atopW2LV.E());
       w1Had.Boost(-boostW);
       w2Had.Boost(-boostW);
       double cosThetaStarHelATopW1 = TMath::Cos(w1Had.Vect().Angle(boostW));
       double cosThetaStarHelATopW2 = TMath::Cos(w2Had.Vect().Angle(boostW));
       hcosThetaStarHelATopW1->Fill(cosThetaStarHelATopW1, weight);
       hcosThetaStarHelATopW2->Fill(cosThetaStarHelATopW2, weight);

     }

     //cout << genTbar.wdau1id << " - " << genTbar.wdau2id << endl;
     // chi* variables pag. 56
     if( genTbar.wdau1id==13 && atopBLV.Pt()>0 && atopW1LV.Pt()>0 && atopW2LV.Pt()>0 && genBLV.Pt()>0){
       LV topLV  = atopBLV+atopW1LV+atopW2LV; //<----------
       LV topWLV = atopW1LV+atopW2LV;

       TVector3 boost(topLV.Px()/topLV.E(), topLV.Py()/topLV.E(), topLV.Pz()/topLV.E());
       
       TLorentzVector bHad;
       TLorentzVector w1Had;
       TLorentzVector w2Had;

       TVector3 boostW(topWLV.Px()/topWLV.E(), topWLV.Py()/topWLV.E(), topWLV.Pz()/topWLV.E());

       w1Had.SetPxPyPzE( atopW1LV.Px(), atopW1LV.Py(), atopW1LV.Pz(), atopW1LV.E());
       w2Had.SetPxPyPzE( atopW2LV.Px(), atopW2LV.Py(), atopW2LV.Pz(), atopW2LV.E());
       bHad.SetPxPyPzE(  atopBLV.Px(),  atopBLV.Py(),  atopBLV.Pz(),  atopBLV.E());
       //bHad.SetPxPyPzE(  genBLV.Px(),  genBLV.Py(),  genBLV.Pz(),  genBLV.E());

       w1Had.Boost(-boostW);
       w2Had.Boost(-boostW);
       bHad.Boost( -boostW);
       double cosChiStarATop = TMath::Cos(w1Had.Vect().Angle( -bHad.Vect() ));
       hcosChiStarATopB->Fill(cosChiStarATop, weight);
     }
     if( genTop.wdau1id==-13 && topBLV.Pt()>0 && topW1LV.Pt()>0 && topW2LV.Pt()>0){
       LV topLV  = topBLV+topW1LV+topW2LV;
       LV topWLV = topW1LV+topW2LV;

       TVector3 boost(topLV.Px()/topLV.E(), topLV.Py()/topLV.E(), topLV.Pz()/topLV.E());
       
       TLorentzVector bHad;
       TLorentzVector w1Had;
       TLorentzVector w2Had;

       TVector3 boostW(topWLV.Px()/topWLV.E(), topWLV.Py()/topWLV.E(), topWLV.Pz()/topWLV.E());

       w1Had.SetPxPyPzE( topW1LV.Px(), topW1LV.Py(), topW1LV.Pz(), topW1LV.E());
       w2Had.SetPxPyPzE( topW2LV.Px(), topW2LV.Py(), topW2LV.Pz(), topW2LV.E());
       bHad.SetPxPyPzE(  topBLV.Px(),  topBLV.Py(),  topBLV.Pz(),  topBLV.E());

       w1Had.Boost(-boostW);
       w2Had.Boost(-boostW);
       bHad.Boost( -boostW);
       double cosChiStarTop = TMath::Cos(w1Had.Vect().Angle( -bHad.Vect() ));
       hcosChiStarTopB->Fill(cosChiStarTop, weight);
     }

     if( abs(genTop.wdau1id)<6 && topBLV.Pt()>0 && topW1LV.Pt()>0 && topW2LV.Pt()>0){
       LV topLV  = topBLV+topW1LV+topW2LV;
       LV topWLV = topW1LV+topW2LV;

       TVector3 boost(topLV.Px()/topLV.E(), topLV.Py()/topLV.E(), topLV.Pz()/topLV.E());
       
       TLorentzVector bHad;
       TLorentzVector w1Had;
       TLorentzVector w2Had;

       TVector3 boostW(topWLV.Px()/topWLV.E(), topWLV.Py()/topWLV.E(), topWLV.Pz()/topWLV.E());

       w1Had.SetPxPyPzE( topW1LV.Px(), topW1LV.Py(), topW1LV.Pz(), topW1LV.E());
       w2Had.SetPxPyPzE( topW2LV.Px(), topW2LV.Py(), topW2LV.Pz(), topW2LV.E());
       bHad.SetPxPyPzE(  topBLV.Px(),  topBLV.Py(),  topBLV.Pz(),  topBLV.E());

       w1Had.Boost(-boostW);
       w2Had.Boost(-boostW);
       bHad.Boost( -boostW);
       double cosChiStarTop1 = TMath::Cos(w1Had.Vect().Angle( -bHad.Vect() ));
       double cosChiStarTop2 = TMath::Cos(w2Had.Vect().Angle( -bHad.Vect() ));
       hcosChiStarTopHadB->Fill(cosChiStarTop1, weight);
       hcosChiStarTopHadB->Fill(cosChiStarTop2, weight);
     }
     if( abs(genTbar.wdau1id)<6 && atopBLV.Pt()>0 && atopW1LV.Pt()>0 && atopW2LV.Pt()>0){
       LV topLV  = atopBLV+atopW1LV+atopW2LV;
       LV topWLV = atopW1LV+atopW2LV;

       TVector3 boost(topLV.Px()/topLV.E(), topLV.Py()/topLV.E(), topLV.Pz()/topLV.E());
       
       TLorentzVector bHad;
       TLorentzVector w1Had;
       TLorentzVector w2Had;

       TVector3 boostW(topWLV.Px()/topWLV.E(), topWLV.Py()/topWLV.E(), topWLV.Pz()/topWLV.E());

       w1Had.SetPxPyPzE( atopW1LV.Px(), atopW1LV.Py(), atopW1LV.Pz(), atopW1LV.E());
       w2Had.SetPxPyPzE( atopW2LV.Px(), atopW2LV.Py(), atopW2LV.Pz(), atopW2LV.E());
       bHad.SetPxPyPzE(  atopBLV.Px(),  atopBLV.Py(),  atopBLV.Pz(),  atopBLV.E());

       w1Had.Boost(-boostW);
       w2Had.Boost(-boostW);
       bHad.Boost( -boostW);
       double cosChiStarATop1 = TMath::Cos(w1Had.Vect().Angle( -bHad.Vect() ));
       double cosChiStarATop2 = TMath::Cos(w2Had.Vect().Angle( -bHad.Vect() ));
       hcosChiStarATopHadB->Fill(cosChiStarATop1, weight);
       hcosChiStarATopHadB->Fill(cosChiStarATop2, weight);
     }


     if(genBLV.Pt()>0 && genBbarLV.Pt()>0){

       unsigned int index1=999;
       for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
	 if(deltaR(myJetsFilt[k], genBLV)<0.30)
	   index1=k;
       }
       unsigned int index2=999;
       for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
	 if(deltaR(myJetsFilt[k], genBbarLV)<0.30)
	   index2=k;
       }
       LV jet1 =  index1!=999 ? myJetsFilt[index1] : genBLV;
       LV jet2 =  index2!=999 ? myJetsFilt[index2] : genBbarLV;

       LV higgsLV = jet1+jet2;//genBLV+genBbarLV;
       
       TVector3 boost(higgsLV.Px()/higgsLV.E(), higgsLV.Py()/higgsLV.E(), higgsLV.Pz()/higgsLV.E());
      
       TLorentzVector b1Had;
       TLorentzVector b2Had;

       b1Had.SetPxPyPzE(  jet1.Px(),  jet1.Py(),  jet1.Pz(),  jet1.E());
       b2Had.SetPxPyPzE(  jet2.Px(),  jet2.Py(),  jet2.Pz(),  jet2.E());

       b1Had.Boost(-boost);
       b2Had.Boost(-boost);

       double cosThetaStarHiggsB    = TMath::Cos(b1Had.Vect().Angle(boost));
       double cosThetaStarHiggsBbar = TMath::Cos(b2Had.Vect().Angle(boost));

       if(higgsLV.M()>90 && higgsLV.M()<150 && higgsLV.Pt()>100 && index1!=999 && index2!=999){
	 hcosThetaStarHiggsB->Fill(TMath::Abs(cosThetaStarHiggsB), weight);
	 if(TMath::Abs(cosThetaStarHiggsB)>0.85) hPt->Fill( deltaR(jet1,jet2));
       }
       hcosThetaStarHiggsBbar->Fill(cosThetaStarHiggsBbar, weight);
     }


     if(abs(genTbar.wdau1id)<6 && atopBLV.Pt()>0 && atopW1LV.Pt()>0 && atopW1LV.Pt()>0 && genBLV.Pt()>0 && genBbarLV.Pt()>0){
       LV topLV  = genBLV+atopW1LV+atopW2LV; // wrong pairing
       LV topWLV = atopW1LV+atopW2LV;

       TVector3 boost(topLV.Px()/topLV.E(), topLV.Py()/topLV.E(), topLV.Pz()/topLV.E());

       TLorentzVector w1Had;
       TLorentzVector w2Had;

       TVector3 boostW(topWLV.Px()/topWLV.E(), topWLV.Py()/topWLV.E(), topWLV.Pz()/topWLV.E());
       w1Had.SetPxPyPzE( atopW1LV.Px(), atopW1LV.Py(), atopW1LV.Pz(), atopW1LV.E());
       w2Had.SetPxPyPzE( atopW2LV.Px(), atopW2LV.Py(), atopW2LV.Pz(), atopW2LV.E());
       w1Had.Boost(-boostW);
       w2Had.Boost(-boostW);
       double cosThetaStarHelATopW1 = TMath::Cos(w1Had.Vect().Angle(boost));
       double cosThetaStarHelATopW2 = TMath::Cos(w2Had.Vect().Angle(boost));
       //hcosThetaStarWrong->Fill(cosThetaStarHelATopW1, weight);
       //hcosThetaStarWrong->Fill(cosThetaStarHelATopW2, weight);

     }

     if(genBLV.Pt()>0 && genBbarLV.Pt()>0 && atopBLV.Pt()>0 && topBLV.Pt()>0){

       unsigned int index1=999;
       for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
	 if(deltaR(myJetsFilt[k], atopBLV)<0.30)
	   index1=k;
       }
       unsigned int index2=999;
       for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
	 if(deltaR(myJetsFilt[k], topBLV)<0.30)
	   index2=k;
       }
       LV jet1 =  index1!=999 ? myJetsFilt[index1] : atopBLV;
       LV jet2 =  index2!=999 ? myJetsFilt[index2] : topBLV;

       LV higgsLV = jet1+jet2;//topBLV;
       
       TVector3 boost(higgsLV.Px()/higgsLV.E(), higgsLV.Py()/higgsLV.E(), higgsLV.Pz()/higgsLV.E());
      
       TLorentzVector b1Had;
       TLorentzVector b2Had;

       b1Had.SetPxPyPzE(  jet1.Px(),  jet1.Py(),  jet1.Pz(),  jet1.E());
       b2Had.SetPxPyPzE(  jet2.Px(),  jet2.Py(),  jet2.Pz(),  jet2.E());

       b1Had.Boost(-boost);
       b2Had.Boost(-boost);

       double cosThetaStarHiggsB    = TMath::Cos(b1Had.Vect().Angle(boost));
       double cosThetaStarHiggsBbar = TMath::Cos(b2Had.Vect().Angle(boost));

       if(higgsLV.M()>90 && higgsLV.M()<150  && higgsLV.Pt()>100 && index1!=999 && index2!=999)
	 hcosThetaStarWrong->Fill(TMath::Abs(cosThetaStarHiggsB), weight);

     }
     

     int matchB    = 0;
     int matchBbar = 0;
     for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
       float csv_k = mapFilt[k].csv>0 ? mapFilt[k].csv : 0.0 ;
       float pt_k  = mapFilt[k].pt>0 ? mapFilt[k].pt : 0.0 ;
       if( !( csv_k>0.679 && pt_k>30 ) ) continue;
       int genMatchH;
       findGenMatch2(genMatchH, myJetsFilt[k] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);
       bool isB    = deltaR( myJetsFilt[k], genBLV)    < 0.5;
       bool isBbar = deltaR( myJetsFilt[k], genBbarLV) < 0.5;
       if( genMatchH == 5 || isB)     matchB++;
       if( genMatchH == -5 || isBbar) matchBbar++;
     }
     if(numBTag==3 && numJets30==5) h35MatchBTag30->Fill(matchB+matchBbar,weight);
     if(numBTag==3 && numJets30>=6) h36MatchBTag30->Fill(matchB+matchBbar,weight);
     if(numBTag==4 && numJets30==5) h45MatchBTag30->Fill(matchB+matchBbar,weight);
     if(numBTag==4 && numJets30>=6) h46MatchBTag30->Fill(matchB+matchBbar,weight);

     matchB    = 0;
     matchBbar = 0;
     for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
       float csv_k = mapFilt[k].csv>0 ? mapFilt[k].csv : 0.0 ;
       float pt_k  = mapFilt[k].pt>0 ? mapFilt[k].pt : 0.0 ;
       if( !( csv_k>0.679 && pt_k>20 ) ) continue;
       int genMatchH;
       findGenMatch2(genMatchH, myJetsFilt[k] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);
       bool isB    = deltaR( myJetsFilt[k], genBLV)    < 0.5;
       bool isBbar = deltaR( myJetsFilt[k], genBbarLV) < 0.5;
       if( genMatchH == 5 || isB)     matchB++;
       if( genMatchH == -5 || isBbar) matchBbar++;   
     }
     if(numBTag==3 && numJets30==5) h35MatchBTag20->Fill(matchB+matchBbar,weight);
     if(numBTag==3 && numJets30>=6) h36MatchBTag20->Fill(matchB+matchBbar,weight);
     if(numBTag==4 && numJets30==5) h45MatchBTag20->Fill(matchB+matchBbar,weight);
     if(numBTag==4 && numJets30>=6) h46MatchBTag20->Fill(matchB+matchBbar,weight);


     matchB    = 0;
     matchBbar = 0;
     for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
       float csv_k = mapFilt[k].csv>0 ? mapFilt[k].csv : 0.0 ;
       float pt_k  = mapFilt[k].pt>0 ? mapFilt[k].pt : 0.0 ;
       if( !( csv_k>=0.00 && pt_k>30 ) ) continue;
       int genMatchH;
       findGenMatch2(genMatchH, myJetsFilt[k] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);
       bool isB    = deltaR( myJetsFilt[k], genBLV)    < 0.5;
       bool isBbar = deltaR( myJetsFilt[k], genBbarLV) < 0.5;
       if( genMatchH == 5 || isB)     matchB++;
       if( genMatchH == -5 || isBbar) matchBbar++;   
     }
     if(numBTag==3 && numJets30==5) h35MatchNoTag30->Fill(matchB+matchBbar,weight);
     if(numBTag==3 && numJets30>=6) h36MatchNoTag30->Fill(matchB+matchBbar,weight);
     if(numBTag==4 && numJets30==5) h45MatchNoTag30->Fill(matchB+matchBbar,weight);
     if(numBTag==4 && numJets30>=6) h46MatchNoTag30->Fill(matchB+matchBbar,weight);

     matchB    = 0;
     matchBbar = 0;
     for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
       float csv_k = mapFilt[k].csv>0 ? mapFilt[k].csv : 0.0 ;
       float pt_k  = mapFilt[k].pt>0 ? mapFilt[k].pt : 0.0 ;
       if( !( csv_k>=0.00 && pt_k>20 ) ) continue;
       int genMatchH;
       findGenMatch2(genMatchH, myJetsFilt[k] , topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);
       bool isB    = deltaR( myJetsFilt[k], genBLV)    < 0.5;
       bool isBbar = deltaR( myJetsFilt[k], genBbarLV) < 0.5;
       if( genMatchH == 5 || isB)     matchB++;
       if( genMatchH == -5 || isBbar) matchBbar++;   
     }
     if(numBTag==3 && numJets30==5) h35MatchNoTag20->Fill(matchB+matchBbar,weight);
     if(numBTag==3 && numJets30>=6) h36MatchNoTag20->Fill(matchB+matchBbar,weight);
     if(numBTag==4 && numJets30==5) h45MatchNoTag20->Fill(matchB+matchBbar,weight);
     if(numBTag==4 && numJets30>=6) h46MatchNoTag20->Fill(matchB+matchBbar,weight);
     

     if(false && numBTag==4 && numJets30>=6 && matchB==0){
       cout << "Event " << i << endl;
       cout << "HiggsB: Pt= " << genBLV.Pt() << ", Eta= " << genBLV.Eta() << ", Phi= " << genBLV.Phi() << endl;
       for(unsigned int k = 0; k < myJetsFilt.size(); k++){
	 cout << "Pt" << k << "= " << myJetsFilt[k].Pt();
	 cout << " Eta" << k << "= " << myJetsFilt[k].Eta();
	 cout << " Phi" << k << "= " << myJetsFilt[k].Phi() << endl;
	 bool isB    = deltaR( myJetsFilt[k], genBLV)    < 0.5;
	 if(isB) cout << "Matched to jet #" << k << endl;
       }
     }
     if(false && numBTag==4 && numJets30>=6 && matchBbar==0){
       cout << "Event " << i << endl;
       cout << "HiggsBbar: Pt= " << genBbarLV.Pt() << ", Eta= " << genBbarLV.Eta() << ", Phi= " << genBbarLV.Phi() << endl;
       for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
	 cout << "Pt" << k << "= " << myJetsFilt[k].Pt();
	 cout << " Eta" << k << "= " << myJetsFilt[k].Eta();
	 cout << " Phi" << k << "= " << myJetsFilt[k].Phi() << endl;
	 bool isB    = deltaR( myJetsFilt[k], genBbarLV)    < 0.5;
	 if(isB) cout << "Matched to jet #" << k << endl;
       }
     }

  }



   /////////////////////////////////////////////////////////////////////////// 

  float total = hDistr->Integral() ;
  cout << "*************** Distribution *****************" << endl;
  cout << "Total: " << total << endl;
  cout << "No merge, all partons: "  << hDistr->GetBinContent(1)/total << endl;
  cout << "No merge, ==1 W lost : "  << hDistr->GetBinContent(2)/total << endl;
  cout << "No merge, >1 W lost  : "  << hDistr->GetBinContent(19)/total << endl;
  cout << "No merge, one b lost : "  << hDistr->GetBinContent(3)/total << endl;
  cout << "No merge, one H lost : "  << hDistr->GetBinContent(4)/total << endl;
  cout << "Ws merge, all partons: "  << hDistr->GetBinContent(5)/total << endl;
  cout << "Ws merge, >=1 W lost : "  << hDistr->GetBinContent(6)/total << endl;
  cout << "Ws merge, one b lost : "  << hDistr->GetBinContent(7)/total << endl;
  cout << "Ws merge, one H lost : "  << hDistr->GetBinContent(8)/total << endl;
  cout << "Bs merge, all partons: "  << hDistr->GetBinContent(9)/total << endl;
  cout << "Bs merge, >=1 W lost : "  << hDistr->GetBinContent(10)/total << endl;
  cout << "Bs merge, one b lost : "  << hDistr->GetBinContent(11)/total << endl;
  cout << "Bs merge, one H lost : "  << hDistr->GetBinContent(12)/total << endl;
  cout << "WB merge, all partons: "  << hDistr->GetBinContent(13)/total << endl;
  cout << "WB merge, >=1 W lost : "  << hDistr->GetBinContent(14)/total << endl;
  cout << "WB merge, one b lost : "  << hDistr->GetBinContent(15)/total << endl;
  cout << "WB merge, one H lost : "  << hDistr->GetBinContent(16)/total << endl;
  cout << ">=1 W >=1 H lost     : "  << hDistr->GetBinContent(17)/total << endl;
  cout << ">=1 W >=1 b lost     : "  << hDistr->GetBinContent(18)/total << endl;
  cout << "Other                : "  << hDistr->GetBinContent(20)/total << endl;

  fOutput->Write();
  fOutput->Close();
  fInput->Close();
  
  delete fOutput;
  
  return  ;

}




void algoTopDown(float seedBtag     = 0.898, float seedPt=80,      float seedEta=1.5, float seedDR=1.5,
		 float seedPtnBtag  = 0.244, float seedPtnPt=30,
		 float antiSeedBtag = 0.679, float antiSeed1Pt=80, float antiSeed2Pt=30,  float antiSeedDR=2.0){


  TFile *fOutput       = new TFile("AlgoTopDown.root","RECREATE");
  TH1F* hCount         = new TH1F("hCount", "visited"   ,10, 0 ,10 );
  TH1F* hSeeds         = new TH1F("hSeeds", "visited"   ,10, 0 ,10 );
  TH1F* hNumResJets20  = new TH1F("hNumResJets20", "visited"   ,10, 0 ,10 );
  TH1F* hNumResJets30  = new TH1F("hNumResJets30", "visited"   ,10, 0 ,10 );
  TH1F* hNumResBTagJets20  = new TH1F("hNumResBTagJets20", "visited"   ,10, 0 ,10 );
  TH1F* hNumResBTagJets30  = new TH1F("hNumResBTagJets20", "visited"   ,10, 0 ,10 );
  TH1F* hMass           = new TH1F("hMass", "Higgs mass"   ,40, 0 ,400 );
  TH1F* hMatch         = new TH1F("hMatch", "visited"   ,4, 0 ,4 );
  TH1F* hMatchAll      = new TH1F("hMatchAll", "visited"   ,4, 0 ,4 );
  
  TFile* fInput = TFile::Open("gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_VType2/v2/DiJetPt_TTH_HToBB_M-125_8TeV-pythia6_VType2.root","READ");
  //TFile* fInput = TFile::Open("gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_VType2/v2/DiJetPt_TTJets_SemiLeptMGDecays_8TeV-madgraph-part_VType2.root","READ");
  TTree* tree = (TTree*)fInput->Get("tree");

  if(!tree){
    cout << "No tree" << endl; return ;
  }

  TH1F* hCounter = (TH1F*)fInput->Get("Count");
  float weight = (0.1302*0.569)*12.1*1000/hCounter->GetBinContent(1);
  cout << weight << endl;

  //float weight = (103.0)*12.1*1000/hCounter->GetBinContent(1);
  //weight /= (500000./tree->GetEntries());
  //cout << weight << endl;

  JetByPt jet1, jet2, jet3, jet4, jet5, jet6, jet7, jet8, jet9, jet10;
  genTopInfo genTop, genTbar;
  genParticleInfo genB, genBbar;
  int nvlep;
  float vLeptonpt[999]; float vLeptoneta[999]; float vLeptonphi[999];
  float met[999];

  int nSvs;
  float Sv_massBCand[999];
  float Sv_massSv[999];
  float Sv_pt[999];
  float Sv_eta[999];
  float Sv_phi[999];
  float Sv_dist3D[999];
  float Sv_dist2D[999];
  float Sv_distSim2D[999];
  float Sv_distSig3D[999];
  float Sv_dist3D_norm[999];

  
  tree->SetBranchAddress("jet1",   &jet1);
  tree->SetBranchAddress("jet2",   &jet2);
  tree->SetBranchAddress("jet3",   &jet3);
  tree->SetBranchAddress("jet4",   &jet4);
  tree->SetBranchAddress("jet5",   &jet5);
  tree->SetBranchAddress("jet6",   &jet6);
  tree->SetBranchAddress("jet7",   &jet7);
  tree->SetBranchAddress("jet8",   &jet8);
  tree->SetBranchAddress("jet9",   &jet9);
  tree->SetBranchAddress("jet10",  &jet10);
  tree->SetBranchAddress("genTop", &genTop);
  tree->SetBranchAddress("genTbar",&genTbar);
  tree->SetBranchAddress("genB",   &genB);
  tree->SetBranchAddress("genBbar",&genBbar);
  tree->SetBranchAddress("vLepton_pt", vLeptonpt);
  tree->SetBranchAddress("vLepton_phi",vLeptonphi);
  tree->SetBranchAddress("vLepton_eta",vLeptoneta);
  tree->SetBranchAddress("nvlep",  &nvlep);
  tree->SetBranchAddress("METtype1p2corr",  met);

  tree->SetBranchAddress("nSvs",&nSvs);
  tree->SetBranchAddress("Sv_massBCand", Sv_massBCand);
  tree->SetBranchAddress("Sv_massSv", Sv_massSv);
  tree->SetBranchAddress("Sv_pt", Sv_pt);
  tree->SetBranchAddress("Sv_eta", Sv_eta);
  tree->SetBranchAddress("Sv_phi", Sv_phi);
  tree->SetBranchAddress("Sv_dist3D", Sv_dist3D);
  tree->SetBranchAddress("Sv_dist2D",Sv_dist2D);
  tree->SetBranchAddress("Sv_distSim2D", Sv_distSim2D);
  tree->SetBranchAddress("Sv_distSig3D", Sv_distSig3D);
  tree->SetBranchAddress("Sv_dist3D_norm", Sv_dist3D_norm);


  Long64_t nentries = tree->GetEntries();
  int totCounter = 0;
  
  for (Long64_t i = 0; i < /*500000*/ nentries ; i++){
    
    if(i%10000==0) cout << i << endl;
    tree->GetEntry(i);
    
    //continue;

    if(nvlep<=0) continue;

    LV topBLV(  0.,0.,0.,0.);
    LV topW1LV( 0.,0.,0.,0.); 
    LV topW2LV( 0.,0.,0.,0.);
    LV atopBLV( 0.,0.,0.,0.); 
    LV atopW1LV(0.,0.,0.,0.); 
    LV atopW2LV(0.,0.,0.,0.);
    LV genBLV(0.,0.,0.,0.);
    LV genBbarLV(0.,0.,0.,0.);
    
    
     ///////////////////////////////////////////////////////////////////
     /////////////////////////////////////////////////////////////////// 
     if(genTop.bmass>0){
       topBLV.SetPt(   genTop.bpt );
       topBLV.SetEta(  genTop.beta );
       topBLV.SetPhi(  genTop.bphi );
       topBLV.SetM(    genTop.bmass );
       topW1LV.SetPt(  genTop.wdau1pt );
       topW1LV.SetEta( genTop.wdau1eta );
       topW1LV.SetPhi( genTop.wdau1phi );
       topW1LV.SetM(   genTop.wdau1mass );
       topW2LV.SetPt(  genTop.wdau2pt );
       topW2LV.SetEta( genTop.wdau2eta );
       topW2LV.SetPhi( genTop.wdau2phi );
       topW2LV.SetM(   genTop.wdau2mass );
     }
     if(genTbar.bmass>0){
	atopBLV.SetPt(   genTbar.bpt );
	atopBLV.SetEta(  genTbar.beta );
	atopBLV.SetPhi(  genTbar.bphi );
	atopBLV.SetM(    genTbar.bmass );
	atopW1LV.SetPt(  genTbar.wdau1pt );
	atopW1LV.SetEta( genTbar.wdau1eta );
	atopW1LV.SetPhi( genTbar.wdau1phi );
	atopW1LV.SetM(   genTbar.wdau1mass );
	atopW2LV.SetPt(  genTbar.wdau2pt );
	atopW2LV.SetEta( genTbar.wdau2eta );
	atopW2LV.SetPhi( genTbar.wdau2phi );
	atopW2LV.SetM(   genTbar.wdau2mass );
     }

    if(genB.mass>0 && genB.momid==25){
       genBLV.SetPt(  genB.pt );
       genBLV.SetEta( genB.eta );
       genBLV.SetPhi( genB.phi );
       genBLV.SetM(   genB.mass );
     }
     if(genBbar.mass>0 && genBbar.momid==25){
       genBbarLV.SetPt(  genBbar.pt );
       genBbarLV.SetEta( genBbar.eta );
       genBbarLV.SetPhi( genBbar.phi );
       genBbarLV.SetM(   genBbar.mass );
     }

     if(!(topBLV.Pt()>0 &&
	  topW1LV.Pt()>0 &&
	  topW2LV.Pt()>0 &&
	  atopBLV.Pt()>0 &&
	  atopW1LV.Pt()>0 &&
	  atopW2LV.Pt()>0 &&
	  genBLV.Pt()>0 &&
	  genBbarLV.Pt()>0)) continue;
     
     vector<LV> myJets;
     std::map<unsigned int, JetByPt> map;

     map[0] = jet1;  map[1] = jet2; map[2] = jet3; map[3] = jet4; 
     map[4] = jet5;  map[5] = jet6; map[6] = jet7; map[7] = jet8; 
     map[8] = jet9;  map[9] = jet10; 

     myJets.push_back(LV(jet1.pt, jet1.eta, jet1.phi, jet1.mass));
     myJets.push_back(LV(jet2.pt, jet2.eta, jet2.phi, jet2.mass));
     myJets.push_back(LV(jet3.pt, jet3.eta, jet3.phi, jet3.mass));
     myJets.push_back(LV(jet4.pt, jet4.eta, jet4.phi, jet4.mass));
     myJets.push_back(LV(jet5.pt, jet5.eta, jet5.phi, jet5.mass));
     myJets.push_back(LV(jet6.pt, jet6.eta, jet6.phi, jet6.mass));
     myJets.push_back(LV(jet7.pt, jet7.eta, jet7.phi, jet7.mass));
     myJets.push_back(LV(jet8.pt, jet8.eta, jet8.phi, jet8.mass));
     myJets.push_back(LV(jet9.pt, jet9.eta, jet9.phi, jet9.mass));
     myJets.push_back(LV(jet10.pt,jet10.eta,jet10.phi,jet10.mass));

     
     vector<LV> myJetsFilt;
     std::map<unsigned int, JetByPt> mapFilt;
     unsigned int iter = 0;

     for(unsigned int k = 0; k<myJets.size() ; k++){

       if(myJets[k].Pt()>20 && TMath::Abs(myJets[k].Eta()) < 2.5 ){

	 myJetsFilt.push_back( myJets[k] );

	 switch(k){
	 case 0:
	   mapFilt[iter] = jet1;
	   break;
	 case 1:
	   mapFilt[iter] = jet2;
	   break;
	 case 2:
	   mapFilt[iter] = jet3;
	   break;
	 case 3:
	   mapFilt[iter] = jet4;
	   break;
	 case 4:
	   mapFilt[iter] = jet5;
	   break;
	 case 5:
	   mapFilt[iter] = jet6;
	   break;
	 case 6:
	   mapFilt[iter] = jet7;
	   break;
	 case 7:
	   mapFilt[iter] = jet8;
	   break;
	 case 8:
	   mapFilt[iter] = jet9;
	   break;
	 case 9:
	   mapFilt[iter] = jet10;
	   break;
	 default:
	   break;
	 }
	 
	 iter++;
       }
     }

     LV lepLV(vLeptonpt[0], vLeptoneta[0], vLeptonphi[0], 0.);
     LV MEtLV(met[0], 0., met[3] , 0.);

     

     hCount->Fill(0.5);

     unsigned int j1Seed  = 999;
     unsigned int j2Seed  = 999;
     unsigned int j1ASeed = 999;
     unsigned int j2ASeed = 999;
     int match=0;     int matchAny=0;
     int numSeeds=0;
     int numResJets20=0;
     int numResJets30=0;
     int numResBTagJets20=0;
     int numResBTagJets30=0;

     vector<unsigned int> visited;

     for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
       float csv_k = mapFilt[k].csv>0 ? mapFilt[k].csv : 0.0 ;
       float pt_k  = mapFilt[k].pt>0 ? mapFilt[k].pt : 0.0 ;
       float eta_k = mapFilt[k].eta;
       if(csv_k>seedBtag && pt_k>seedPt && TMath::Abs(eta_k)<seedEta){ // valid seed
	 j1Seed = k;
	  for(unsigned int l = 0; l < myJetsFilt.size(); l++){
	    if(l==j1Seed) continue;
	    float csv_l = mapFilt[l].csv>0 ? mapFilt[l].csv : 0.0 ;
	    float pt_l  = mapFilt[l].pt>0 ? mapFilt[l].pt : 0.0 ;
	    float eta_l = mapFilt[l].eta;
	    if(deltaR( myJetsFilt[j1Seed], myJetsFilt[l])< seedDR && csv_l>seedPtnBtag && pt_l>seedPtnPt){
	      j2Seed = l;

	      for(unsigned int m = 0; m < myJetsFilt.size()-1; m++){
		if( m==j1Seed || m==j2Seed) continue;
		float csv_m = mapFilt[m].csv>0 ? mapFilt[m].csv : 0.0 ;
		float pt_m  = mapFilt[m].pt>0 ? mapFilt[m].pt : 0.0 ;
		float eta_m = mapFilt[m].eta;

		for(unsigned int n = m+1; n < myJetsFilt.size(); n++){
		  if( n==j1Seed || n==j2Seed) continue;
		  float csv_n = mapFilt[n].csv>0 ? mapFilt[n].csv : 0.0 ;
		  float pt_n  = mapFilt[n].pt>0 ? mapFilt[n].pt : 0.0 ;
		  float eta_n = mapFilt[n].eta;

		  if( deltaR(myJetsFilt[n],myJetsFilt[m])>antiSeedDR &&  csv_m>antiSeedBtag && 
		      csv_n>antiSeedBtag && pt_m>antiSeed1Pt && pt_n>antiSeed2Pt){
		    j1ASeed = m;
		    j2ASeed = n;
		    if( bookKeeper2(j1Seed,j2Seed,j1ASeed,j2ASeed, visited ) ) numSeeds++;
		    if( numSeeds==1){
		      int isJ1B    = genBLV.Pt()>0 ?    deltaR( myJetsFilt[j1Seed], genBLV)    < 0.3 && TMath::Abs(myJetsFilt[j1Seed].Pt()-genBLV.Pt())/genBLV.Pt()<0.25        : 0;
		      int isJ1Bbar = genBbarLV.Pt()>0 ? deltaR( myJetsFilt[j1Seed], genBbarLV) < 0.3 && TMath::Abs(myJetsFilt[j1Seed].Pt()-genBbarLV.Pt())/genBbarLV.Pt()<0.25  : 0;
		      int isJ2B    = genBLV.Pt()>0 ?    deltaR( myJetsFilt[j2Seed], genBLV)    < 0.3 && TMath::Abs(myJetsFilt[j2Seed].Pt()-genBLV.Pt())/genBLV.Pt()<0.25        : 0;
		      int isJ2Bbar = genBbarLV.Pt()>0 ? deltaR( myJetsFilt[j2Seed], genBbarLV) < 0.3 && TMath::Abs(myJetsFilt[j2Seed].Pt()-genBbarLV.Pt())/genBbarLV.Pt()<0.25  : 0;
		    
		      //hMatch->Fill(isJ1B+isJ1Bbar+isJ2B+isJ2Bbar);
		      match = isJ1B+isJ1Bbar+isJ2B+isJ2Bbar;
		      hMass->Fill( (myJetsFilt[j1Seed] +  myJetsFilt[j2Seed]).M() );
		    }
		    int isAnyB=0; int isAnyBbar=0;
		    for(unsigned int o = 0; o < myJetsFilt.size(); o++){

		      float csv_o = mapFilt[o].csv>0 ? mapFilt[o].csv : 0.0 ;
		      float pt_o  = mapFilt[o].pt>0 ? mapFilt[o].pt : 0.0 ;

		      if(pt_o>seedPtnPt && csv_o>seedBtag &&  genBLV.Pt()>0     &&   deltaR( myJetsFilt[o], genBLV)    < 0.3 && TMath::Abs(myJetsFilt[o].Pt()-genBLV.Pt())/genBLV.Pt()<0.25) isAnyB++;
		      if(pt_o>seedPtnPt && csv_o>seedBtag &&  genBbarLV.Pt()>0  &&   deltaR( myJetsFilt[o], genBbarLV) < 0.3 && TMath::Abs(myJetsFilt[o].Pt()-genBbarLV.Pt())/genBbarLV.Pt()<0.25) isAnyBbar++;

		      if(o==j1Seed || o==j2Seed || o==j1ASeed || o==j2ASeed) continue;
		      if(numSeeds!=1) continue;
		      
		      if(csv_o>-999  && pt_o>20)  numResJets20++;
		      if(csv_o>-999  && pt_o>30)  numResJets30++;
		      if(csv_o>0.679 && pt_o>20)  numResBTagJets20++;
		      if(csv_o>0.679 && pt_o>30)  numResBTagJets30++;
		      
		    }

		    //if( deltaR( genBLV, genBbarLV)> seedDR || TMath::Max(genBLV.Pt(),genBbarLV.Pt())<seedPt ) isAnyB = -isAnyBbar;
		    //if( numSeeds==1) hMatchAll->Fill(isAnyB+isAnyBbar);
		    matchAny=isAnyB+isAnyBbar;
		  }

		}

	      }

	    }

	  }
	 
       }
     }

     if(numSeeds==1){
       hMatch->Fill(match);
       hMatchAll->Fill(matchAny);
     }

     hSeeds->Fill(numSeeds , weight);
     hNumResJets20->Fill( numResJets20, weight);
     hNumResJets30->Fill( numResJets30, weight);
     hNumResBTagJets20->Fill( numResBTagJets20, weight);
     hNumResBTagJets30->Fill( numResBTagJets30, weight);
     
  }

  fOutput->Write();
  fOutput->Close();
  fInput->Close();
  
  delete fOutput;
  
  return  ;

}
