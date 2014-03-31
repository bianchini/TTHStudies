#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include <cstdlib>
#include <iostream> 
#include <fstream>
#include <map>
#include <string>

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
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TLegend.h"
#include "TList.h"
#include "THStack.h"
#include "TCut.h"
#include "TArrayF.h"
#include "TObjArray.h"
#include "TVector3.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TKey.h"
#include "TMultiGraph.h"
#include "Bianchi/TTHStudies/interface/Samples.h"
#include "Bianchi/TTHStudies/interface/Test.h"
#include "Math/GenVector/LorentzVector.h"

#include <algorithm>

#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"

#include "DataFormats/FWLite/interface/LuminosityBlock.h"
#include "DataFormats/FWLite/interface/Run.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "PhysicsTools/FWLite/interface/TFileService.h"

// event shapes
#include "PhysicsTools/CandUtils/interface/EventShapeVariables.h"
#include "TopQuarkAnalysis/TopTools/interface/TopologyWorker.h"

#define BTAGTHR     0.679
#define GENJETDR    0.30

#define LEADLEPISOZEE   0.10;
#define TRAILLEPISOZEE  0.20;
#define LEADLEPPTZEE    30;
#define TRAILLEPPTZEE   20;
#define LEADLEPETAZEE   2.5;
#define TRAILLEPETAZEE  2.5;

#define LEADLEPISOZMM   0.10;
#define TRAILLEPISOZMM  0.20;
#define LEADLEPPTZMM    30;
#define TRAILLEPPTZMM   20;
#define LEADLEPETAZMM   2.1;
#define TRAILLEPETAZMM  2.4;

using namespace std;

typedef  ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LV;

struct sorterByPt {
  bool operator() (float i,float j) const { return (i>j);}
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

typedef struct
 {
   float hdau1mass;
   float hdau1pt;
   float hdau1eta;
   float hdau1phi;
   float hdau1id;
   float hdau2mass;
   float hdau2pt;
   float hdau2eta;
   float hdau2phi;
   float hdau2id;
   float hmass;
   float hpt;
   float heta;
   float hphi;
   float hid;
 } recoHiggsInfo;



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


  
void resetHiggsInfo(recoHiggsInfo &info){
   info.hdau1mass = -99;
   info.hdau1pt = -99;
   info.hdau1eta = -99;
   info.hdau1phi = -99;
   info.hdau1id = -99;
   info.hdau2mass = -99;
   info.hdau2pt = -99;
   info.hdau2eta = -99;
   info.hdau2phi = -99;
   info.hdau2id = -99;
   info.hmass = -99;
   info.hpt = -99;
   info.heta = -99;
   info.hphi = -99;
   info.hid = -99;
}



void resetTopInfo(recoTopInfo &info){
  info.bmass = -99;
  info.bpt = -99;
  info.beta = -99;
  info.bphi = -99;
  info.bstatus = -99;
  info.wdau1mass = -99;
  info.wdau1pt = -99;
  info.wdau1eta = -99;
  info.wdau1phi = -99;
  info.wdau1id = -99;
  info.wdau2mass = -99;
  info.wdau2pt = -99;
  info.wdau2eta = -99;
  info.wdau2phi = -99;
  info.wdau2id = -99;
  info.topmass = -99;
  info.toppt  = -99;
  info.topeta = -99;
  info.topphi = -99;
  info.topid  = -99;
  info.bbarmass = -99;
  info.bbarpt  = -99;
  info.bbareta = -99;
  info.bbarphi = -99;
}


float bestHiggsMass( std::vector<LV> leptons, std::vector<LV> mets, 
		     std::vector<LV> jetsUntag, std::vector<LV> jetsBtag, 
		     std::vector<float> jetUncUntag, 
		     std::vector<float> jetUncBtag, 
		     bool verbose){

  if( leptons.size() != 1){
    if(verbose) cout << "Exactly one lepton is required... return" << endl;
    return -99.;
  }
  if( mets.size() < 1){
    if(verbose) cout << "MET is required... return" << endl;
    return -99.;
  }
  unsigned int untagJets   = jetsUntag.size();
  unsigned int btagJets    = jetsBtag.size();
  unsigned int jecuncBtag  = jetUncBtag.size();
  unsigned int jecuncUntag = jetUncUntag.size();

  if( untagJets < 2 || btagJets<4){
    if(verbose) cout << "2 untag jets and 4-b tag jets are required... return" << endl;
    return -99.;
  }
  if( jecuncBtag != btagJets || jecuncUntag!=untagJets ){
     if(verbose) cout << "Need JEC for all input jets... return" << endl;
    return -99.;
  }

  float WMass   = 80.385;
  float WMass2  = WMass*WMass;
  float WMassG  = 2.085;
  float WMassG2 = WMassG*WMassG;
  float TMass   = 173.5;
  float TMassG  = 2.0;
  float TMassG2 = TMassG*TMassG;

  LV lepton = leptons[0];
  LV met    = mets[0];

  // solve for neutrino pz
  float El    = lepton.E();  
  float El_z  = lepton.Pz();
  float Delta = WMass2 - El*El + (lepton+met).Pt()*(lepton+met).Pt() + El_z*El_z;
  float a     = 4*(El*El - El_z*El_z);
  float b     = -4*Delta*El_z;
  float c     = 4*met.Pt()*met.Pt()*El*El - Delta*Delta;

  float Enu_z1 = 0.0;
  float Enu_z2 = 0.0;

  if( (b*b - 4*a*c) <0){
    if(verbose) cout << "Negative determinant!!!" << endl;
    Enu_z1 = 0.0;
    Enu_z2 = 0.0;
  }
  else{
    Enu_z1 = (-b + TMath::Sqrt(b*b - 4*a*c))/(2.*a);
    Enu_z2 = (-b - TMath::Sqrt(b*b - 4*a*c))/(2.*a);
  }

  LV neut1 = met;
  neut1.SetEta( TMath::ACosH( TMath::Sqrt( met.Pt()*met.Pt() + Enu_z1*Enu_z1 )/met.Pt() ) );
  LV neut2 = met;
  neut2.SetEta( TMath::ACosH( TMath::Sqrt( met.Pt()*met.Pt() + Enu_z2*Enu_z2 )/met.Pt() ) );

  vector<LV> neutrinos;
  vector<unsigned int> mask;
  neutrinos.push_back( neut1 );
  neutrinos.push_back( neut2 );

  float minChi2 = 99999.;

  for(unsigned int nu = 0 ; nu<neutrinos.size() ; nu++){

    LV neut = neutrinos[nu];
    LV W1p4 = neut + lepton;

    for(unsigned int i = 0; i<untagJets-1; i++ ){
      for(unsigned int j = i+1; j<untagJets; j++ ){

	LV jet1 = jetsUntag[i];
	LV jet2 = jetsUntag[j];

	LV W2p4 = jet1+jet2;

	for(unsigned int k = 0; k<btagJets-1; k++ ){
	  LV bjet1 = jetsBtag[k];
	  for(unsigned int l = k+1; l<btagJets; l++ ){
	    LV bjet2 = jetsBtag[l];

	    LV top11 = W1p4+bjet1;
	    LV top22 = W2p4+bjet2;
	    LV top12 = W1p4+bjet2;
	    LV top21 = W2p4+bjet1;

	    float dMW2_2 =  W2p4.M()*W2p4.M()*( jetUncUntag[i] + jetUncUntag[j] );

	    float dMtop11_2 = W1p4.M()*W1p4.M()*( 0.0                              ) + bjet1.M()*bjet1.M()*( jetUncBtag[k] ) + 2*( 0.0                             + jetUncBtag[k])*( bjet1.Dot(W1p4) );
	    float dMtop22_2 = W2p4.M()*W2p4.M()*( jetUncUntag[i] + jetUncUntag[j]  ) + bjet2.M()*bjet2.M()*( jetUncBtag[l] ) + 2*( jetUncUntag[i] + jetUncUntag[j] + jetUncBtag[l])*( bjet2.Dot(W2p4) );
	    float dMtop12_2 = W1p4.M()*W1p4.M()*( 0.0                              ) + bjet2.M()*bjet2.M()*( jetUncBtag[l] ) + 2*( 0.0                             + jetUncBtag[l])*( bjet2.Dot(W1p4) );
	    float dMtop21_2 = W2p4.M()*W2p4.M()*( jetUncUntag[i] + jetUncUntag[j]  ) + bjet1.M()*bjet1.M()*( jetUncBtag[k] ) + 2*( jetUncUntag[i] + jetUncUntag[j] + jetUncBtag[k])*( bjet1.Dot(W2p4) );

	    float chi2_11 = 
	      TMath::Power(TMass - top11.M() , 2)/TMath::Sqrt( TMath::Power(dMtop11_2/2./top11.M(),2) + TMassG2 ) + 
	      TMath::Power(WMass - W2p4.M(),   2)/TMath::Sqrt( TMath::Power(dMW2_2/2./W2p4.M()    ,2) + WMassG2 ); 

	    float chi2_12 = 
	      TMath::Power(TMass - top12.M() , 2)/TMath::Sqrt( TMath::Power(dMtop12_2/2./top12.M(),2) + TMassG2 ) + 
	      TMath::Power(WMass - W2p4.M(),   2)/TMath::Sqrt( TMath::Power(dMW2_2/2./W2p4.M()    ,2) + WMassG2 ); 

	    float chi2_21 = 
	      TMath::Power(TMass - top21.M() , 2)/TMath::Sqrt( TMath::Power(dMtop21_2/2./top21.M(),2) + TMassG2 ) + 
	      TMath::Power(WMass - W2p4.M(),   2)/TMath::Sqrt( TMath::Power(dMW2_2/2./W2p4.M()    ,2) + WMassG2 ); 

	    float chi2_22 = 
	      TMath::Power(TMass - top22.M() , 2)/TMath::Sqrt( TMath::Power(dMtop22_2/2./top22.M(),2) + TMassG2 ) + 
	      TMath::Power(WMass - W2p4.M(),   2)/TMath::Sqrt( TMath::Power(dMW2_2/2./W2p4.M()    ,2) + WMassG2 );

	    float minComb = TMath::Min( TMath::Min( chi2_11, chi2_12), TMath::Min( chi2_21, chi2_22)); 

	    if( minComb  < minChi2){
	      minChi2 = minComb;
	      mask.clear();
	      mask.push_back(k);
	      mask.push_back(l);
	    }

	  }
	}

      }
    }
  }

  if( mask.size() == 2){
    LV HiggsB1, HiggsB2; int counter = -1;
    for(unsigned int k = 0; k<btagJets; k++){
      if( k == mask[0] ||  k == mask[1] ) continue;
      counter++;
      if( counter  == 0) HiggsB1 = jetsBtag[k];
      if( counter  == 1) HiggsB2 = jetsBtag[k];
    }

    return (HiggsB1+HiggsB2).M();

  }
  else{
    if(verbose) cout << "Could not find a pair of btag jets..." << endl;
    return -99.;
  }

  return -99.;

}


//  -99 => jet not matched to partons
//    3  => jet not matched to top or higgs decay products
//    0  => jet matched to higgs b quarks
// +/-1  => jet matched to top/antitop b quarks
// +/-2  => jet matched to top/antitop W quarks
//    4  => ambiguous match

void findGenMatch(int& genMatch, LV genJet, LV topBLV, LV topW1LV, LV topW2LV, LV atopBLV, LV atopW1LV, LV atopW2LV, LV genBLV, LV genBbarLV){

  if(genJet.Pt()<=0.001){
    genMatch = -99;
    return;
  }

  int numMatches = 0;

  if( topBLV.Pt()>0 && Geom::deltaR(genJet, topBLV)       < GENJETDR ){
    genMatch = 1;
    numMatches++;
  }
  if( topW1LV.Pt()>0 && Geom::deltaR(genJet, topW1LV)     < GENJETDR ){
    genMatch = 2;
    numMatches++;
  }
  if( topW2LV.Pt()>0 && Geom::deltaR(genJet, topW2LV)     < GENJETDR ){
    genMatch = 2;
    numMatches++;
  }
  if( atopBLV.Pt()>0 && Geom::deltaR(genJet, atopBLV)     < GENJETDR ){
    genMatch = -1;
    numMatches++;
  }
  if( atopW1LV.Pt()>0 && Geom::deltaR(genJet, atopW1LV)   < GENJETDR ){
    genMatch = -2;
    numMatches++;
  }
  if( atopW2LV.Pt()>0 && Geom::deltaR(genJet, atopW2LV)   < GENJETDR ){
    genMatch = -2;
    numMatches++;
  }
  if( genBLV.Pt()>0 && Geom::deltaR(genJet, genBLV)       < GENJETDR ){
    genMatch = 0;
    numMatches++;
  }
  if( genBbarLV.Pt()>0 && Geom::deltaR(genJet, genBbarLV) < GENJETDR ){
    genMatch = 0;
    numMatches++;
  }

  if(numMatches==0)
    genMatch = 3; // NO MATCH 
  if(numMatches>1)
    genMatch = 4; // AMBIGUOUS MATCH



  return;

}

void setTopDecay(float& topB_,  float& topW_, 
		 float& atopB_, float& atopW_, float& higgsB_, 
		 int genMatch){

  topB_  = int(genMatch==1);
  topW_  = int(genMatch==2);
  atopB_ = int(genMatch==-1);
  atopW_ = int(genMatch==-2);
  higgsB_= int(genMatch==0);

  if(genMatch==3 || genMatch==4){
    topB_  = int(genMatch);
    topW_  = int(genMatch);
    atopB_ = int(genMatch);
    atopW_ = int(genMatch);
    higgsB_= int(genMatch);
  }

  return;
}


bool jsonContainsEvent (const std::vector< edm::LuminosityBlockRange > &jsonVec,
                        unsigned long run, unsigned long lumi
			)
{
  // if the jsonVec is empty, then no JSON file was provided so all
  // events should pass
  if (jsonVec.empty())
    {
      return true;
    }
  bool (* funcPtr) (edm::LuminosityBlockRange const &,
		    edm::LuminosityBlockID const &) = &edm::contains;
  edm::LuminosityBlockID lumiID (run,lumi); 
  std::vector< edm::LuminosityBlockRange >::const_iterator iter = 
    std::find_if (jsonVec.begin(), jsonVec.end(),
		  boost::bind(funcPtr, _1, lumiID) );
  return jsonVec.end() != iter;

}


int main(int argc, const char* argv[])
{

  std::cout << "Step3Clone" << std::endl;
  gROOT->SetBatch(true);
 
  gSystem->Load("libFWCoreFWLite");
  gSystem->Load("libDataFormatsFWLite");

  AutoLibraryLoader::enable();

  PythonProcessDesc builder(argv[1]);
 
  const edm::ParameterSet& in = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteInput");

  std::string pathToFile( in.getParameter<std::string>("pathToFile" ) );
  std::string outPath(    in.getParameter<std::string>("outPath" ) );

  std::string ordering(   in.getParameter<std::string>("ordering" ) );
  double lumi(            in.getParameter<double>("lumi" ) );
  bool verbose(           in.getParameter<bool>("verbose" ) );

  int eventsToDebug = -1;
  if(in.exists("eventsToDebug")) eventsToDebug = in.getParameter<int>("eventsToDebug") ;
  int eventsToDebugCounter = 0;


  std::vector<edm::LuminosityBlockRange> jsonVector;
  if ( in.exists("lumisToProcess") ) 
    {
      std::vector<edm::LuminosityBlockRange> const & lumisTemp =
	in.getUntrackedParameter<std::vector<edm::LuminosityBlockRange> > ("lumisToProcess");
      jsonVector.resize( lumisTemp.size() );
      copy( lumisTemp.begin(), lumisTemp.end(), jsonVector.begin() );
    }
  

  edm::VParameterSet samples;
  if( in.exists("samples") && argc==2 ) samples = in.getParameter<edm::VParameterSet>("samples") ;
  else{
    if(argc==4){
      string fName = string((argv[2])); 
      string fNick = string((argv[3]));
      string cut   = "DUMMY";
      edm::ParameterSet pset;
      pset.addParameter("skip",false);
      pset.addParameter("name",    fName);
      pset.addParameter("nickName",fNick);
      pset.addParameter("color",1);
      pset.addParameter("xSec",-1.);
      pset.addParameter("update",false);
      pset.addParameter("cut",cut);
      samples.push_back( pset );
    }
    else{
      cout << "You need to specify name and nickname for your sample" << endl;
      return 0;
    }
  }

 

  bool openAllFiles = false;
  Samples* mySamples = new Samples(openAllFiles, pathToFile, ordering, samples, lumi, verbose);
  map<string, TH1F*> mapHist;
  vector<string> mySampleFiles;

  if(mySamples->IsOk()){

    cout << "Ok!" << endl;
    mySampleFiles = mySamples->Files();

    for( unsigned int i = 0 ; i < mySampleFiles.size(); i++){
      string sampleName       = mySampleFiles[i];

      if(verbose){
	cout << mySampleFiles[i] << " ==> " << mySamples->GetXSec(sampleName) 
	     << " pb,"
	     << " ==> weight = "            << mySamples->GetWeight(sampleName) << endl;
      }
    }
  }
  else{
    cout << "Problems... leaving" << endl;
    return 0;
  }


  for(unsigned int i = 0 ; i < mySampleFiles.size(); i++){
    
    string currentName       = mySampleFiles[i];
    string fileName          = mySamples->GetFileName(currentName);

    string outputName = outPath+fileName+".root"; 

    //string cleanSE  = "srmrm "+outputName;
    //cout << cleanSE << endl;
    //gSystem->Exec(cleanSE.c_str());
    
    TFile* fs      = TFile::Open(outputName.c_str(), "RECREATE");

    mySamples->OpenFile( currentName );
    TIter nextkey((mySamples->GetFile( currentName ))->GetListOfKeys());
    TH1F *key;
    while ( (key = (TH1F*)nextkey()) ) {
      string name(key->GetName());
      if( name.find("tree")!=string::npos) continue;
      TH1F* h = (TH1F*)(mySamples->GetFile( currentName ))->FindObjectAny( key->GetName());
      fs->cd();
      h->Write(key->GetName());
    }
    cout << "Start copying tree " <<  endl;
    TTree* outTree = (mySamples->GetTree( currentName, "tree"))->CloneTree();
    if(!outTree){
      cout << "Null tree... continue" << endl;
      continue;
    }
    cout << "Done..." << endl;
    

    //////////////////////////////////NEW VARIABLES///////////////////////////////
    //////////////////////////////////////////////////////////////////////////////

    int myJson_;
    int hJetRank_[100];
    int aJetRank_[100];

    JetByPt jet1_, jet2_, jet3_, jet4_, jet5_, jet6_, jet7_, jet8_, jet9_, jet10_;

    int nLF_,    nC_,    nB_;
    int nLFTop_, nCTop_, nBTop_;
    int numOfBs_, numOfBsAcc_;
    int numOfBsFlav_, numOfBsFlavAcc_, numOfCsFlav_, numOfCsFlavAcc_;

    int numJets40_; int numJets40bTag_;
    int numJets30_; int numJets30bTag_;
    int numJets20_; int numJets20bTag_;
    float isotropy_, circularity_, sphericity_, aplanarity_, Cparam_, Dparam_;
    float thrust0_, thrust1_, thrust2_;
    float sphericity2_, aplanarity2_; 
    float h10_,h20_,h30_,h40_,h50_,h60_;
    float aveCsv_;
    float aveDeltaRbTag_;
    float aveMbTag_, aveMunTag_;
    float closestJJbTagMass_;
    float varCsv_;
    float maxCsv_;  float minCsv_;
    float sumAllPt_;
    float sumAllJet20Pt_; float sumAllJet30Pt_; float sumAllJet40Pt_;
    float massAll_;
    float massLJ_;
    float M3_;
    float MHT30_, MHT20_, MHT40_;
    float minDeltaRLJ_, minDeltaRBtag_;
    float firstBtag_, secondBtag_, thirdBtag_, fourthBtag_;
    float bestHiggsMass_;

    recoTopInfo recoTop_, recoTbar_, genMatchTop_, genMatchTbar_;
    recoHiggsInfo recoHiggs_, genMatchHiggs_;


    TBranch *hJetRankBR = outTree->Branch("hJetRank",hJetRank_,"hJetRank[nhJets]/I");
    TBranch *aJetRankBR = outTree->Branch("aJetRank",aJetRank_,"aJetRank[naJets]/I");
    
    TBranch *jet1BR = outTree->Branch("jet1",&jet1_,"index/I:pt/F:eta/F:phi/F:mass/F:csv/F:topB/F:topW/F:atopB/F:atopW/F:higgsB/F:flavor/F:unc/F");
    TBranch *jet2BR = outTree->Branch("jet2",&jet2_,"index/I:pt/F:eta/F:phi/F:mass/F:csv/F:topB/F:topW/F:atopB/F:atopW/F:higgsB/F:flavor/F:unc/F");
    TBranch *jet3BR = outTree->Branch("jet3",&jet3_,"index/I:pt/F:eta/F:phi/F:mass/F:csv/F:topB/F:topW/F:atopB/F:atopW/F:higgsB/F:flavor/F:unc/F");
    TBranch *jet4BR = outTree->Branch("jet4",&jet4_,"index/I:pt/F:eta/F:phi/F:mass/F:csv/F:topB/F:topW/F:atopB/F:atopW/F:higgsB/F:flavor/F:unc/F");
    TBranch *jet5BR = outTree->Branch("jet5",&jet5_,"index/I:pt/F:eta/F:phi/F:mass/F:csv/F:topB/F:topW/F:atopB/F:atopW/F:higgsB/F:flavor/F:unc/F");
    TBranch *jet6BR = outTree->Branch("jet6",&jet6_,"index/I:pt/F:eta/F:phi/F:mass/F:csv/F:topB/F:topW/F:atopB/F:atopW/F:higgsB/F:flavor/F:unc/F");
    TBranch *jet7BR = outTree->Branch("jet7",&jet7_,"index/I:pt/F:eta/F:phi/F:mass/F:csv/F:topB/F:topW/F:atopB/F:atopW/F:higgsB/F:flavor/F:unc/F");
    TBranch *jet8BR = outTree->Branch("jet8",&jet8_,"index/I:pt/F:eta/F:phi/F:mass/F:csv/F:topB/F:topW/F:atopB/F:atopW/F:higgsB/F:flavor/F:unc/F");
    TBranch *jet9BR = outTree->Branch("jet9",&jet9_,"index/I:pt/F:eta/F:phi/F:mass/F:csv/F:topB/F:topW/F:atopB/F:atopW/F:higgsB/F:flavor/F:unc/F");
    TBranch *jet10BR= outTree->Branch("jet10",&jet10_,"index/I:pt/F:eta/F:phi/F:mass/F:csv/F:topB/F:topW/F:atopB/F:atopW/F:higgsB/F:flavor/F:unc/F");


    TBranch *nLFBR     = outTree->Branch("nLF",   &nLF_,   "nLF/I");
    TBranch *nCBR      = outTree->Branch("nC",    &nC_,    "nC/I");
    TBranch *nBBR      = outTree->Branch("nB",    &nB_,    "nB/I");
    TBranch *nLFTopBR  = outTree->Branch("nLFTop",&nLFTop_,"nLFTop/I");
    TBranch *nCTopBR   = outTree->Branch("nCTop", &nCTop_, "nCTop/I");
    TBranch *nBTopBR   = outTree->Branch("nBTop", &nBTop_, "nBTop/I");

    TBranch *numJets40BR     = outTree->Branch("numJets40",&numJets40_,"numJets40/I");
    TBranch *numJets30BR     = outTree->Branch("numJets30",&numJets30_,"numJets30/I");
    TBranch *numJets20BR     = outTree->Branch("numJets20",&numJets20_,"numJets20/I");
    TBranch *numJets40bTagBR = outTree->Branch("numJets40bTag",&numJets40bTag_,"numJets40bTag/I");
    TBranch *numJets30bTagBR = outTree->Branch("numJets30bTag",&numJets30bTag_,"numJets30bTag/I");
    TBranch *numJets20bTagBR = outTree->Branch("numJets20bTag",&numJets20bTag_,"numJets20bTag/I");

    TBranch *isotropyBR      = outTree->Branch("isotropy",&isotropy_,"isotropy/F");
    TBranch *circularityBR   = outTree->Branch("circularity",&circularity_,"circularity/F");
    TBranch *sphericityBR    = outTree->Branch("sphericity",&sphericity_,"sphericity/F");
    TBranch *aplanarityBR    = outTree->Branch("aplanarity",&aplanarity_,"aplanarity/F");
    TBranch *CparamBR        = outTree->Branch("Cparam",&Cparam_,"Cparam/F");
    TBranch *DparamBR        = outTree->Branch("Dparam",&Dparam_,"Dparam/F");

    TBranch *thrust0BR = outTree->Branch("thrust0",&thrust0_,"thrust0/F");
    TBranch *thrust1BR = outTree->Branch("thrust1",&thrust1_,"thrust1/F");
    TBranch *thrust2BR = outTree->Branch("thrust2",&thrust2_,"thrust2/F");
    TBranch *sphericity2BR = outTree->Branch("sphericity2",&sphericity2_,"sphericity2/F");
    TBranch *aplanarity2BR = outTree->Branch("aplanarity2",&aplanarity2_,"aplanarity2/F");
    TBranch *h10BR = outTree->Branch("h10",&h10_,"h10/F");
    TBranch *h20BR = outTree->Branch("h20",&h20_,"h20/F");
    TBranch *h30BR = outTree->Branch("h30",&h30_,"h30/F");
    TBranch *h40BR = outTree->Branch("h40",&h40_,"h40/F");
    TBranch *h50BR = outTree->Branch("h50",&h50_,"h50/F");
    TBranch *h60BR = outTree->Branch("h60",&h60_,"h60/F");

    TBranch *aveCsvBR            = outTree->Branch("aveCsv",&aveCsv_,"aveCsv/F");
    TBranch *aveDeltaRbTagBR     = outTree->Branch("aveDeltaRbTag",&aveDeltaRbTag_,"aveDeltaRbTag/F");
    TBranch *aveMbTagBR          = outTree->Branch("aveMbTag",&aveMbTag_,"aveMbTag/F");
    TBranch *aveMunTagBR         = outTree->Branch("aveMunTag",&aveMunTag_,"aveMunTag/F");
    TBranch *closestJJbTagMassBR = outTree->Branch("closestJJbTagMass",&closestJJbTagMass_,"closestJJbTagMass/F");

    TBranch *varCsvBR            = outTree->Branch("varCsv",&varCsv_,"varCsv/F");
    TBranch *maxCsvBR            = outTree->Branch("maxCsv",&maxCsv_,"maxCsv/F");
    TBranch *minCsvBR            = outTree->Branch("minCsv",&minCsv_,"minCsv/F");

    TBranch *sumAllPtBR          = outTree->Branch("sumAllPt",&sumAllPt_,"sumAllPt/F");
    TBranch *sumAllJet20PtBR     = outTree->Branch("sumAllJet20Pt",&sumAllJet20Pt_,"sumAllJet20Pt/F");
    TBranch *sumAllJet30PtBR     = outTree->Branch("sumAllJet30Pt",&sumAllJet30Pt_,"sumAllJet30Pt/F");
    TBranch *sumAllJet40PtBR     = outTree->Branch("sumAllJet40Pt",&sumAllJet40Pt_,"sumAllJet40Pt/F");

    TBranch *massAllBR           = outTree->Branch("massAll",&massAll_,"massAll/F");
    TBranch *massLJBR            = outTree->Branch("massLJ",&massLJ_,"massLJ/F");
    TBranch *M3BR                = outTree->Branch("M3",&M3_,"M3/F");

    TBranch *MHT20BR = outTree->Branch("MHT20",&MHT20_,"MHT20/F");
    TBranch *MHT30BR = outTree->Branch("MHT30",&MHT30_,"MHT30/F");
    TBranch *MHT40BR = outTree->Branch("MHT40",&MHT40_,"MHT40/F");

    TBranch *minDeltaRLJBR   = outTree->Branch("minDeltaRLJ",&minDeltaRLJ_,"minDeltaRLJ/F");
    TBranch *minDeltaRBtagBR = outTree->Branch("minDeltaRBtag",&minDeltaRBtag_,"minDeltaRBtag/F");
    
    TBranch *firstBtagBR  = outTree->Branch("firstBtag", &firstBtag_,"firstBtag/F");
    TBranch *secondBtagBR = outTree->Branch("secondBtag",&secondBtag_,"secondBtag/F");
    TBranch *thirdBtagBR  = outTree->Branch("thirdBtag", &thirdBtag_,"thirdBtag/F");
    TBranch *fourthBtagBR = outTree->Branch("fourthBtag",&fourthBtag_,"fourthBtag/F");

    TBranch *bestHiggsMassBR  = outTree->Branch("bestHiggsMass",&bestHiggsMass_,"bestHiggsMass/F");

    TBranch *myJsonBR         = outTree->Branch("myJson",&myJson_,"myJson/F");
    
    TBranch *numOfBsBR        = outTree->Branch("numOfBs",   &numOfBs_,   "numOfBs/I");
    TBranch *numOfBsAccBR     = outTree->Branch("numOfBsAcc",&numOfBsAcc_,"numOfBsAcc/I");
    TBranch *numOfBsFlavBR    = outTree->Branch("numOfBsFlav",   &numOfBsFlav_,   "numOfBsFlav/I");
    TBranch *numOfBsFlavAccBR = outTree->Branch("numOfBsFlavAcc",&numOfBsFlavAcc_,"numOfBsFlavAcc/I");
    TBranch *numOfCsFlavBR    = outTree->Branch("numOfCsFlav",   &numOfCsFlav_,   "numOfCsFlav/I");
    TBranch *numOfCsFlavAccBR = outTree->Branch("numOfCsFlavAcc",&numOfCsFlavAcc_,"numOfCsFlavAcc/I");
    
    TBranch *recoTopBR       = outTree->Branch("recoTop",      &recoTop_,"bmass/F:bpt/F:beta:bphi/F:bstatus/F:wdau1mass/F:wdau1pt/F:wdau1eta:wdau1phi/F:wdau1id/F:wdau2mass/F:wdau2pt/F:wdau2eta:wdau2phi/F:wdau2id/F:wmass/F:wpt/F:weta/F:wphi/F:wid/F:topmass/F:toppt/F:topeta/F:topphi/F:topid/F:bbarmass/F:bbarpt/F:bbareta/F:bbarphi/F");
    TBranch *genMatchTopBR   = outTree->Branch("genMatchTop",  &genMatchTop_,"bmass/F:bpt/F:beta:bphi/F:bstatus/F:wdau1mass/F:wdau1pt/F:wdau1eta:wdau1phi/F:wdau1id/F:wdau2mass/F:wdau2pt/F:wdau2eta:wdau2phi/F:wdau2id/F:wmass/F:wpt/F:weta/F:wphi/F:wid/F:topmass/F:toppt/F:topeta/F:topphi/F:topid/F:bbarmass/F:bbarpt/F:bbareta/F:bbarphi/F");
    TBranch *recoTbarBR      = outTree->Branch("recoTbar",     &recoTbar_,"bmass/F:bpt/F:beta:bphi/F:bstatus/F:wdau1mass/F:wdau1pt/F:wdau1eta:wdau1phi/F:wdau1id/F:wdau2mass/F:wdau2pt/F:wdau2eta:wdau2phi/F:wdau2id/F:wmass/F:wpt/F:weta/F:wphi/F:wid/F:topmass/F:toppt/F:topeta/F:topphi/F:topid/F:bbarmass/F:bbarpt/F:bbareta/F:bbarphi/F");
    TBranch *genMatchTbarBR  = outTree->Branch("genMatchTbar", &genMatchTbar_,"bmass/F:bpt/F:beta:bphi/F:bstatus/F:wdau1mass/F:wdau1pt/F:wdau1eta:wdau1phi/F:wdau1id/F:wdau2mass/F:wdau2pt/F:wdau2eta:wdau2phi/F:wdau2id/F:wmass/F:wpt/F:weta/F:wphi/F:wid/F:topmass/F:toppt/F:topeta/F:topphi/F:topid/F:bbarmass/F:bbarpt/F:bbareta/F:bbarphi/F");

    TBranch *recoHiggsBR     = outTree->Branch("recoHiggs",    &recoHiggs_,"hdau1mass/F:hdau1pt/F:hdau1eta:hdau1phi/F:hdau1id/F:hdau2mass/F:hdau2pt/F:hdau2eta:hdau2phi/F:hdau2id/F:hmass/F:hpt/F:heta/F:hphi/F:hid/F");
    TBranch *genMatchHiggsBR = outTree->Branch("genMatchHiggs",    &genMatchHiggs_,"hdau1mass/F:hdau1pt/F:hdau1eta:hdau1phi/F:hdau1id/F:hdau2mass/F:hdau2pt/F:hdau2eta:hdau2phi/F:hdau2id/F:hmass/F:hpt/F:heta/F:hphi/F:hid/F");



    //TBranch *BR = outTree->Branch("",&_,"/F");
    //////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////





    //////////////////////////////////INPUT VARIABLES/////////////////////////////
    //////////////////////////////////////////////////////////////////////////////

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
      int run;
      int lumi;
      int event;
      int json;
    } EventInfo;

    int nhJets;
    int naJets;
    int nvlep;
    int Vtype;
    int nSimBs;
    //int run, event, lumi;
    EventInfo ev;

    float met[999];
    float hJetspt[999]; float hJetseta[999]; float hJetsphi[999]; float hJetse[999]; float hJetscsv[999]; float hJetsunc[999]; float hJetsflavor[999]; float hJetsgenpt[999];float hJetsgeneta[999];float hJetsgenphi[999]; float hJetspuJetIdL[999]; bool hJetsid[999];
    float aJetspt[999]; float aJetseta[999]; float aJetsphi[999]; float aJetse[999]; float aJetscsv[999]; float aJetsunc[999]; float aJetsflavor[999]; float aJetsgenpt[999];float aJetsgeneta[999];float aJetsgenphi[999]; float aJetspuJetIdL[999]; bool aJetsid[999];
    float vLeptonpt[999]; float vLeptoneta[999]; float vLeptonphi[999];
    float vLeptonpfCombRelIso[999];
    float SimBspt[999]; float SimBseta[999]; float SimBsphi[999]; float SimBsmass[999];

    genTopInfo genTop, genTbar;
    genParticleInfo genB, genBbar;

    outTree->SetBranchAddress("nhJets",   &nhJets);
    outTree->SetBranchAddress("naJets",   &naJets);
    outTree->SetBranchAddress("Vtype",    &Vtype);
    outTree->SetBranchAddress("nSimBs",   &nSimBs);

    outTree->SetBranchAddress("hJet_pt",   hJetspt);
    outTree->SetBranchAddress("hJet_eta",  hJetseta);
    outTree->SetBranchAddress("hJet_phi",  hJetsphi);
    outTree->SetBranchAddress("hJet_e",    hJetse);
    outTree->SetBranchAddress("hJet_csv",  hJetscsv);
    outTree->SetBranchAddress("hJet_JECUnc", hJetsunc);
    outTree->SetBranchAddress("hJet_flavour",hJetsflavor);
    outTree->SetBranchAddress("hJet_genPt",  hJetsgenpt);
    outTree->SetBranchAddress("hJet_genEta", hJetsgeneta);
    outTree->SetBranchAddress("hJet_genPhi", hJetsgenphi);
    outTree->SetBranchAddress("hJet_puJetIdL", hJetspuJetIdL);
    outTree->SetBranchAddress("hJet_id", hJetsid);


    outTree->SetBranchAddress("aJet_pt",   aJetspt);
    outTree->SetBranchAddress("aJet_eta",  aJetseta);
    outTree->SetBranchAddress("aJet_phi",  aJetsphi);
    outTree->SetBranchAddress("aJet_e",    aJetse);
    outTree->SetBranchAddress("aJet_csv",  aJetscsv);
    outTree->SetBranchAddress("aJet_JECUnc", aJetsunc);
    outTree->SetBranchAddress("aJet_flavour",aJetsflavor);
    outTree->SetBranchAddress("aJet_genPt",  aJetsgenpt);
    outTree->SetBranchAddress("aJet_genEta", aJetsgeneta);
    outTree->SetBranchAddress("aJet_genPhi", aJetsgenphi);
    outTree->SetBranchAddress("aJet_puJetIdL", aJetspuJetIdL);
    outTree->SetBranchAddress("aJet_id", aJetsid);

    outTree->SetBranchAddress("vLepton_pt", vLeptonpt);
    outTree->SetBranchAddress("vLepton_phi",vLeptonphi);
    outTree->SetBranchAddress("vLepton_eta",vLeptoneta);
    outTree->SetBranchAddress("vLepton_pfCombRelIso",vLeptonpfCombRelIso);
    outTree->SetBranchAddress("nvlep",      &nvlep);

    outTree->SetBranchAddress("METtype1p2corr",  met);

    outTree->SetBranchAddress("genTop", &genTop);
    outTree->SetBranchAddress("genTbar",&genTbar);
    outTree->SetBranchAddress("genB",   &genB);
    outTree->SetBranchAddress("genBbar",&genBbar);

    outTree->SetBranchAddress("SimBs_pt",  SimBspt);
    outTree->SetBranchAddress("SimBs_eta" ,SimBseta );
    outTree->SetBranchAddress("SimBs_phi" ,SimBsphi );
    outTree->SetBranchAddress("SimBs_mass",SimBsmass );

    //outTree->SetBranchAddress("EVENT.run",  &run);
    //outTree->SetBranchAddress("EVENT.event",&event);
    //outTree->SetBranchAddress("EVENT.lumi", &lumi);

    outTree->SetBranchAddress("EVENT",  &ev);
  
    //////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////



    Long64_t nentries = outTree->GetEntries();
 
    for (Long64_t i = 0; i < nentries; i++){

      if(i%10000==0) cout << i << endl;

      jet1_.reset();
      jet2_.reset();
      jet3_.reset();
      jet4_.reset();
      jet5_.reset();
      jet6_.reset();
      jet7_.reset();
      jet8_.reset();
      jet9_.reset();
      jet10_.reset();
    
      nLF_    = 0; nC_    = 0; nB_    = 0;
      nLFTop_ = 0; nCTop_ = 0; nBTop_ = 0;

      numOfBs_=0; numOfBsAcc_=0;
      numOfBsFlav_=0;  numOfBsFlavAcc_=0;
      numOfCsFlav_=0;  numOfCsFlavAcc_=0;

      thrust0_       = -99;  thrust1_ = -99; thrust2_ = -99; 
      sphericity2_   = -99;  aplanarity2_ = -99;
      h10_           = -99;  h20_ = -99; h30_ = -99; h40_ = -99; h50_ = -99; h60_ = -99;
      massLJ_ = -99; 
      varCsv_ = 0.;
      firstBtag_     = -99;  secondBtag_ = -99; thirdBtag_ = -99; fourthBtag_ = -99;
      sumAllPt_      = 0.;   sumAllJet20Pt_  = 0.; sumAllJet30Pt_  = 0.; sumAllJet40Pt_  = 0.;
      numJets20_     = 0;    numJets30_     = 0; numJets40_     = 0;
      numJets30bTag_ = 0;    numJets20bTag_ = 0; numJets40bTag_ = 0;
      minDeltaRLJ_   = -99;  M3_ = -99;
      minDeltaRBtag_ = -99.; bestHiggsMass_ = -99.;

      resetTopInfo(recoTop_);
      resetTopInfo(genMatchTop_);
      resetTopInfo(recoTbar_);
      resetTopInfo(genMatchTbar_);
      resetHiggsInfo(recoHiggs_);
      resetHiggsInfo(genMatchHiggs_);


      ////////////////////////////////////////////////

      outTree->GetEntry(i);

      ////////////////////////////////////////////////


      std::vector<math::RhoEtaPhiVector> allJetsForShapeVar;
      TObjArray jetArrayForShapeVar;
      TObjArray leptArrayForShapeVar;
      jetArrayForShapeVar.SetOwner();
      leptArrayForShapeVar.SetOwner();

      std::vector<LV> allBs;
      for(int k = 0; k < nSimBs; k++){
	LV bLV(SimBspt[k], SimBseta[k], SimBsphi[k], SimBsmass[k]);
	allBs.push_back(bLV);
      }
     


      std::vector<LV> allLeptons;
      std::vector<LV> allMets;
      std::vector<LV> allBtagJetsForShapeVar;
      std::vector<LV> allUntagJetsForShapeVar;
      std::vector<float> csvBtag;
      vector<float> jetUncUntag;
      vector<float> jetUncBtag;

      std::map<float, int, sorterByPt> hMapPt;
      std::map<float, int, sorterByPt> aMapPt;
      std::map<float, int, sorterByPt> allMapPt20;
      std::map<float, int, sorterByPt> allMapPt30;
      std::map<float, int, sorterByPt> bTagMap30;

      std::vector<LV> topHiggsDecay;

      LV topBLV(  0.,0.,0.,0.);
      LV topW1LV( 0.,0.,0.,0.); 
      LV topW2LV( 0.,0.,0.,0.);
      LV atopBLV( 0.,0.,0.,0.); 
      LV atopW1LV(0.,0.,0.,0.); 
      LV atopW2LV(0.,0.,0.,0.);

      LV recotopBLV(  0.,0.,0.,0.);
      LV recotopW1LV( 0.,0.,0.,0.); 
      LV recotopW2LV( 0.,0.,0.,0.);
      LV recoatopBLV( 0.,0.,0.,0.); 
      LV recoatopW1LV(0.,0.,0.,0.); 
      LV recoatopW2LV(0.,0.,0.,0.);

      LV recohiggsB1LV(0.,0.,0.,0.); 
      LV recohiggsB2LV(0.,0.,0.,0.);

      int topB_flag   = 0;
      int topW1_flag  = 0;
      int topW2_flag  = 0;
      int atopB_flag  = 0;
      int atopW1_flag = 0;
      int atopW2_flag = 0;
      int higgsB1_flag= 0;
      int higgsB2_flag= 0;

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

	topHiggsDecay.push_back(topBLV);
	if(abs(genTop.wdau1id)<6) topHiggsDecay.push_back(topW1LV);
	if(abs(genTop.wdau2id)<6) topHiggsDecay.push_back(topW2LV);
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

	topHiggsDecay.push_back(atopBLV);
	if(abs(genTbar.wdau1id)<6) topHiggsDecay.push_back(atopW1LV);
	if(abs(genTbar.wdau2id)<6) topHiggsDecay.push_back(atopW2LV);
      }
      LV genBLV(0.,0.,0.,0.);
      LV genBbarLV(0.,0.,0.,0.);
      if(genB.mass>0 && genB.momid==25){
	genBLV.SetPt(  genB.pt );
	genBLV.SetEta( genB.eta );
	genBLV.SetPhi( genB.phi );
	genBLV.SetM(   genB.mass );
	topHiggsDecay.push_back(genBLV);
      }
      if(genBbar.mass>0 && genBbar.momid==25){
	genBbarLV.SetPt(  genBbar.pt );
	genBbarLV.SetEta( genBbar.eta );
	genBbarLV.SetPhi( genBbar.phi );
	genBbarLV.SetM(   genBbar.mass );
	topHiggsDecay.push_back(genBbarLV);
      }


      float sumBtag  = 0.;
      LV MHT20LV(0.,0.,0.,0.);
      LV MHT30LV(0.,0.,0.,0.);
      LV MHT40LV(0.,0.,0.,0.);

      for(int i = 0; i < nhJets; i++){

	hMapPt[ hJetspt[i]]   = i;

	LV genJetLV(hJetsgenpt[i], hJetsgeneta[i], hJetsgenphi[i], 0.0);
	for(unsigned int b = 0; b < allBs.size(); b++){
	  LV bLV = allBs[b];
	  if( genJetLV.Pt()>20. && Geom::deltaR(bLV,genJetLV )< 0.5 ){
	    numOfBs_++;
	    if( TMath::Abs(genJetLV.Eta())<2.5 )
	      numOfBsAcc_++;
	  }
	}
	if( genJetLV.Pt()>20. && TMath::Abs( hJetsflavor[i])==5 ){
	  numOfBsFlav_++;	
	  if( TMath::Abs(genJetLV.Eta())<2.5 ) numOfBsFlavAcc_++;
	}
	if( genJetLV.Pt()>20. && TMath::Abs( hJetsflavor[i])==4 ){
	  numOfCsFlav_++;	
	  if( TMath::Abs(genJetLV.Eta())<2.5 ) numOfCsFlavAcc_++;
	}


	bool jetID = hJetspuJetIdL[i]>0.5 && TMath::Abs(hJetseta[i])<5.0 && hJetsid[i];
	if(!jetID){
	  if(verbose) cout << "hJet " << i << "does not pass" << endl;
	  continue;
	}


	float jetMass2 = hJetse[i]*hJetse[i] -  TMath::Power(hJetspt[i]*TMath::CosH(hJetseta[i]) ,2);
	LV jetLV(hJetspt[i], hJetseta[i], hJetsphi[i], TMath::Sqrt(jetMass2));

	if( topBLV.Pt()>0   && Geom::deltaR(topBLV, jetLV) <GENJETDR ){
	  recotopBLV   = jetLV;
	  topB_flag++;
	}
	else if( abs(genTop.wdau1id)<6 && topW1LV.Pt()>0  && Geom::deltaR(topW1LV,jetLV) <GENJETDR ){
	  recotopW1LV  = jetLV;
	  topW1_flag++;
	}
	else if( abs(genTop.wdau2id)<6 && topW2LV.Pt()>0  && Geom::deltaR(topW2LV,jetLV) <GENJETDR ){
	  recotopW2LV  = jetLV;
	  topW2_flag++;
	}
	else if( atopBLV.Pt()>0  && Geom::deltaR(atopBLV, jetLV)<GENJETDR ){
	  recoatopBLV  = jetLV;
	  atopB_flag++;
	}
	else if( abs(genTbar.wdau1id)<6 && atopW1LV.Pt()>0 && Geom::deltaR(atopW1LV,jetLV)<GENJETDR ){
	  recoatopW1LV = jetLV;
	  atopW1_flag++;
	}
	else if( abs(genTbar.wdau2id)<6 && atopW2LV.Pt()>0 && Geom::deltaR(atopW2LV,jetLV)<GENJETDR ){
	  recoatopW2LV = jetLV;
	  atopW2_flag++;
	}
	else if( genBLV.Pt()>0 && Geom::deltaR(genBLV,jetLV)<GENJETDR ){
	  recohiggsB1LV = jetLV;
	  higgsB1_flag++;
	}
	else if( genBbarLV.Pt()>0 && Geom::deltaR(genBbarLV,jetLV)<GENJETDR ){
	  recohiggsB2LV = jetLV;
	  higgsB2_flag++;
	}
	else{}


	if( hJetspt[i] > 20.){
	  numJets20_++;
	  if( hJetscsv[i] > BTAGTHR) numJets20bTag_++;
	  MHT20LV += jetLV;
	  sumAllJet20Pt_ +=  hJetspt[i];
	  allMapPt20[ hJetspt[i]] = i;
	}
	if( hJetspt[i]>30.){
	  numJets30_++;
	  MHT30LV += jetLV;
	  if( hJetscsv[i] > BTAGTHR){
	    numJets30bTag_++;
	    sumBtag += hJetscsv[i];
	    allBtagJetsForShapeVar.push_back(  jetLV );
	    csvBtag.push_back( hJetscsv[i] );
	    jetUncBtag.push_back( hJetsunc[i] );
	  }
	  else{
	    allUntagJetsForShapeVar.push_back( jetLV );
	    jetUncUntag.push_back( hJetsunc[i] );
	  }
	  allMapPt30[ hJetspt[i]] = i;
	  sumAllJet30Pt_ +=  hJetspt[i];

	  allJetsForShapeVar.push_back( math::RhoEtaPhiVector(hJetspt[i], hJetseta[i], hJetsphi[i]) );
	  TVector3 *tv3 = new TVector3(jetLV.Px(),jetLV.Py(),jetLV.Pz());
	  jetArrayForShapeVar.Add( tv3 );

	}
	if( hJetspt[i]>40.){
	  numJets40_++;
	  if( hJetscsv[i] > BTAGTHR) numJets40bTag_++;
	  MHT40LV += jetLV;
	  sumAllJet40Pt_ +=  hJetspt[i];
	}
      }


      for(int i = 0; i < naJets; i++){

	aMapPt[ aJetspt[i]] = i;

	LV genJetLV(aJetsgenpt[i], aJetsgeneta[i], aJetsgenphi[i], 0.0);
	for(unsigned int b = 0; b < allBs.size(); b++){
	  LV bLV = allBs[b];
	  if( genJetLV.Pt()>20. && Geom::deltaR(bLV,genJetLV )< 0.5 ){
	    numOfBs_++;
	    if( TMath::Abs(genJetLV.Eta())<2.5 )
	      numOfBsAcc_++;
	  }
	}
	if( genJetLV.Pt()>20. && TMath::Abs( aJetsflavor[i])==5 ){
	  numOfBsFlav_++;	
	  if( TMath::Abs(genJetLV.Eta())<2.5 ) numOfBsFlavAcc_++;
	}
	if( genJetLV.Pt()>20. && TMath::Abs( aJetsflavor[i])==4 ){
	  numOfCsFlav_++;	
	  if( TMath::Abs(genJetLV.Eta())<2.5 ) numOfCsFlavAcc_++;
	}


	bool jetID = aJetspuJetIdL[i]>0.5 && TMath::Abs(aJetseta[i])<5.0 && aJetsid[i];
	if(!jetID){
	  if(verbose) cout << "aJet " << i << "does not pass" << endl;
	  continue;
	}

	float jetMass2 = aJetse[i]*aJetse[i] -  TMath::Power(aJetspt[i]*TMath::CosH(aJetseta[i]),2);
	LV jetLV(aJetspt[i], aJetseta[i], aJetsphi[i],  TMath::Sqrt(jetMass2));

	if( topBLV.Pt()>0   && Geom::deltaR(topBLV, jetLV) <GENJETDR ){
          recotopBLV   = jetLV;
          topB_flag++;
        }
        else if( topW1LV.Pt()>0  && Geom::deltaR(topW1LV,jetLV) <GENJETDR ){
          recotopW1LV  = jetLV;
          topW1_flag++;
        }
        else if( topW2LV.Pt()>0  && Geom::deltaR(topW2LV,jetLV) <GENJETDR ){
          recotopW2LV  = jetLV;
          topW2_flag++;
        }
        else if( atopBLV.Pt()>0  && Geom::deltaR(atopBLV, jetLV)<GENJETDR ){
          recoatopBLV  = jetLV;
          atopB_flag++;
        }
        else if( atopW1LV.Pt()>0 && Geom::deltaR(atopW1LV,jetLV)<GENJETDR ){
          recoatopW1LV = jetLV;
          atopW1_flag++;
        }
        else if( atopW2LV.Pt()>0 && Geom::deltaR(atopW2LV,jetLV)<GENJETDR ){
          recoatopW2LV = jetLV;
          atopW2_flag++;
        }
	else if( genBLV.Pt()>0 && Geom::deltaR(genBLV,jetLV)<GENJETDR ){
	  recohiggsB1LV = jetLV;
	  higgsB1_flag++;
	}
	else if( genBbarLV.Pt()>0 && Geom::deltaR(genBbarLV,jetLV)<GENJETDR ){
	  recohiggsB2LV = jetLV;
	  higgsB2_flag++;
	}
	else{}


	if( aJetspt[i] > 20.){
	  numJets20_++;
	  if( aJetscsv[i] > BTAGTHR) numJets20bTag_++;
	  MHT20LV += jetLV;
	  sumAllJet20Pt_ +=  aJetspt[i];
	  allMapPt20[ aJetspt[i]] = -i-1;// distinguish jet collection
	}
	if( aJetspt[i]>30.){
	  numJets30_++;
	  MHT30LV += jetLV;
	  if( aJetscsv[i] > BTAGTHR){
	    numJets30bTag_++;
	    sumBtag += aJetscsv[i];
	    allBtagJetsForShapeVar.push_back( jetLV );
	    csvBtag.push_back( aJetscsv[i] );
	    jetUncBtag.push_back( aJetsunc[i] );
	  }
	  else{
	    allUntagJetsForShapeVar.push_back( jetLV );
	    jetUncUntag.push_back( aJetsunc[i] );
	  }
	  allMapPt30[ aJetspt[i]] = -i-1; // distinguish jet collection
	  sumAllJet30Pt_ +=  aJetspt[i];

	  allJetsForShapeVar.push_back( math::RhoEtaPhiVector(aJetspt[i], aJetseta[i], aJetsphi[i]) );
	  TVector3 *tv3 = new TVector3(jetLV.Px(),jetLV.Py(),jetLV.Pz());
	  jetArrayForShapeVar.Add( tv3 );
	  
	}
	if( aJetspt[i]>40.){
	  numJets40_++;
	  if( aJetscsv[i] > BTAGTHR) numJets40bTag_++;
	  MHT40LV += jetLV;
	  sumAllJet40Pt_ +=  aJetspt[i];
	}
      }





      MHT20_ = MHT20LV.Pt()>0 ? MHT20LV.Pt() : -99;
      MHT30_ = MHT30LV.Pt()>0 ? MHT30LV.Pt() : -99;
      MHT40_ = MHT40LV.Pt()>0 ? MHT40LV.Pt() : -99;
      aveCsv_ = numJets30bTag_>0 ? sumBtag/numJets30bTag_ : -99;

      float sumDeltaRbTag = 0.;
      float sumMassBtag   = 0.;
      float sumMassUntag  = 0.;
      int numOfComb       = 0;

      if(allBtagJetsForShapeVar.size()>1){
	float minDeltaR = 999.;
	for(unsigned k = 0; k< allBtagJetsForShapeVar.size()-1 ; k++){
	  for(unsigned l = k+1; l<allBtagJetsForShapeVar.size() ; l++){

	    LV first  = allBtagJetsForShapeVar[k];
	    LV second = allBtagJetsForShapeVar[l];
	    float deltaR = TMath::Sqrt( TMath::Power(first.Eta()-second.Eta(),2) + TMath::Power(first.Phi()-second.Phi(),2)  );

	    sumDeltaRbTag += deltaR;
	    sumMassBtag   += (first+second).M();
	    if(deltaR < minDeltaR){
	      minDeltaRBtag_ = deltaR;
	      closestJJbTagMass_ = (first+second).M();
	      minDeltaR = deltaR;
	    }
	    numOfComb++;
	  }
	}
      }

      aveDeltaRbTag_ = numOfComb>0 ? sumDeltaRbTag/numOfComb : -99;
      aveMbTag_      = numOfComb>0 ? sumMassBtag/numOfComb   : -99;


      float maxCsv = -999.; float minCsv = 999.;
      for(unsigned k = 0; k< csvBtag.size() ; k++){
	varCsv_ +=  TMath::Power(csvBtag[k]-aveCsv_,2);
	if( csvBtag[k]>maxCsv ) maxCsv = csvBtag[k];
	if( csvBtag[k]<minCsv ) minCsv = csvBtag[k];
      }
      maxCsv_ = maxCsv;
      minCsv_ = minCsv;
      varCsv_ = numJets30bTag_>1 ? varCsv_/(numJets30bTag_-1) : -99;
  
      numOfComb =  0;
      if(allUntagJetsForShapeVar.size()>1){
	for(unsigned k = 0; k< allUntagJetsForShapeVar.size()-1 ; k++){
	  for(unsigned l = k+1; l<allUntagJetsForShapeVar.size() ; l++){
	    LV first  = allUntagJetsForShapeVar[k];
	    LV second = allUntagJetsForShapeVar[l];
	    sumMassUntag   += (first+second).M();
	    numOfComb++;
	  }
	}
      }
      aveMunTag_     = numOfComb>0 ? sumMassUntag/numOfComb : -99;
      

      LV allP4(0.,0.,0.,0.);
      LV MEtP4(met[0], 0., met[3] , 0.);
      allMets.push_back(MEtP4);    
    
      float leadLepIso,  leadLepPt,  leadLeptEta;
      float trailLepIso, trailLepPt, trailLeptEta;

      switch(Vtype){
      case 0:
	leadLepIso   = LEADLEPISOZMM;
	leadLepPt    = LEADLEPPTZMM;
	leadLeptEta  = LEADLEPETAZMM;
	trailLepIso  = TRAILLEPISOZMM;
	trailLepPt   = TRAILLEPPTZMM;
	trailLeptEta = TRAILLEPETAZMM;
	break;
      case 1:
	leadLepIso   = LEADLEPISOZEE;
	leadLepPt    = LEADLEPPTZEE;
	leadLeptEta  = LEADLEPETAZEE;
	trailLepIso  = TRAILLEPISOZEE;
	trailLepPt   = TRAILLEPPTZEE;
	trailLeptEta = TRAILLEPETAZEE;
	break;
      case 2:
	leadLepIso   = LEADLEPISOZMM;
	leadLepPt    = LEADLEPPTZMM;
	leadLeptEta  = LEADLEPETAZMM;
	trailLepIso  = TRAILLEPISOZMM;
	trailLepPt   = TRAILLEPPTZMM;
	trailLeptEta = TRAILLEPETAZMM;
	break;
      case 3:
	leadLepIso   = LEADLEPISOZEE;
	leadLepPt    = LEADLEPPTZEE;
	leadLeptEta  = LEADLEPETAZEE;
	trailLepIso  = TRAILLEPISOZEE;
	trailLepPt   = TRAILLEPPTZEE;
	trailLeptEta = TRAILLEPETAZEE;
	break;
      default:
	leadLepIso   = LEADLEPISOZMM;
	leadLepPt    = LEADLEPPTZMM;
	leadLeptEta  = LEADLEPETAZMM;
	trailLepIso  = TRAILLEPISOZMM;
	trailLepPt   = TRAILLEPPTZMM;
	trailLeptEta = TRAILLEPETAZMM;
	break;
      }


      for(int k = 0; k<nvlep ; k++){
	bool tight = vLeptonpfCombRelIso[k]<leadLepIso  && vLeptonpt[k]>leadLepPt  && TMath::Abs(vLeptoneta[k])<leadLeptEta; 
	bool loose = vLeptonpfCombRelIso[k]<trailLepIso && vLeptonpt[k]>trailLepPt && TMath::Abs(vLeptoneta[k])<trailLeptEta;
	if( tight || loose){
	  sumAllPt_ += vLeptonpt[k];
	  LV lepLV(vLeptonpt[k], vLeptoneta[k], vLeptonphi[k], 0.);
	  allLeptons.push_back( lepLV );
	  TVector3 *tv3 = new TVector3(lepLV.Px(),lepLV.Py(),lepLV.Pz());
	  leptArrayForShapeVar.Add( tv3 );
	}
      }



      for(unsigned int l = 0; l< allLeptons.size() ; l++){

	LV lep = allLeptons[l];
	if( abs(genTop.wdau1id)>=11 && topW1LV.Pt()>0  && Geom::deltaR(topW1LV,lep) <GENJETDR ){
	  recotopW1LV  = lep;
	  topW1_flag++;
	}
	if( abs(genTop.wdau2id)>=11 && topW2LV.Pt()>0  && Geom::deltaR(topW2LV,lep) <GENJETDR ){
	  recotopW2LV  = lep;
	  topW2_flag++;
	}
	if( abs(genTbar.wdau1id)>=11 && atopW1LV.Pt()>0 && Geom::deltaR(atopW1LV,lep)<GENJETDR ){
	  recoatopW1LV = lep;
	  atopW1_flag++;
	}
	if( abs(genTbar.wdau2id)>=11 && atopW2LV.Pt()>0 && Geom::deltaR(atopW2LV,lep)<GENJETDR ){
	  recoatopW2LV = lep;
	  atopW2_flag++;
	}
      }
      





      float minDeltaR = 999.;
      for(unsigned k = 0; k< allBtagJetsForShapeVar.size() ; k++){
	for(unsigned l = 0; l<allLeptons.size() ; l++){
	  LV first  = allBtagJetsForShapeVar[k];
	  LV second = allLeptons[l];
	  float deltaR = TMath::Sqrt( TMath::Power(first.Eta()-second.Eta(),2) + TMath::Power(first.Phi()-second.Phi(),2)  );
	  if(deltaR < minDeltaR){
	    massLJ_ = (first+second).M();
	    minDeltaR = deltaR;
	  }
	}
      }


      if(allBtagJetsForShapeVar.size()>0 && allUntagJetsForShapeVar.size()>1){
	float highestPt = -999;
	for(unsigned k = 0; k< allBtagJetsForShapeVar.size() ; k++){
	  LV zero = allBtagJetsForShapeVar[k];
	  for(unsigned l = 0; l<allUntagJetsForShapeVar.size()-1 ; l++){
	    for(unsigned m = l+1; m<allUntagJetsForShapeVar.size() ; m++){
	      LV first  = allUntagJetsForShapeVar[l];
	      LV second = allUntagJetsForShapeVar[m];
	      if( (zero+first+second).Pt()> highestPt){
		highestPt = (zero+first+second).Pt();
		M3_ = (zero+first+second).M();
	      }
	    }
	  }

	}
      }



      if(allLeptons.size()>0 && (allBtagJetsForShapeVar.size()>0 || allUntagJetsForShapeVar.size()>0)){

	LV first = allLeptons[0];
	float minDeltaR = 999.;

	for(unsigned k = 0; k< allBtagJetsForShapeVar.size() ; k++){
	  LV second = allBtagJetsForShapeVar[k];
	  float deltaR = TMath::Sqrt( TMath::Power(first.Eta()-second.Eta(),2) + TMath::Power(first.Phi()-second.Phi(),2)  );
	  if( deltaR < minDeltaR){
	    minDeltaRLJ_ = deltaR;
	    minDeltaR    = deltaR;
	  }
	}

	for(unsigned k = 0; k< allUntagJetsForShapeVar.size() ; k++){
	  LV second = allUntagJetsForShapeVar[k];
	  float deltaR = TMath::Sqrt( TMath::Power(first.Eta()-second.Eta(),2) + TMath::Power(first.Phi()-second.Phi(),2)  );
	  if( deltaR < minDeltaR){
	    minDeltaRLJ_ = deltaR;
	    minDeltaR    = deltaR;
	  }
	}
      }
      

      for(unsigned k = 0; k< allLeptons.size() ; k++){
	allP4 += allLeptons[k];
      }
      for(unsigned k = 0; k< allUntagJetsForShapeVar.size() ; k++){
	sumAllPt_ += allUntagJetsForShapeVar[k].Pt();
	allP4     += allUntagJetsForShapeVar[k];
      }
      for(unsigned k = 0; k< allBtagJetsForShapeVar.size() ; k++){
	sumAllPt_ += allBtagJetsForShapeVar[k].Pt();
	allP4     += allBtagJetsForShapeVar[k];
      }

      sumAllPt_ += MEtP4.Pt();
      allP4     += MEtP4;

      massAll_ = allP4.M();


      if( allLeptons.size() == 1 && allMets.size()>0 && allUntagJetsForShapeVar.size()>1 && allBtagJetsForShapeVar.size()>=4)
	bestHiggsMass_ = bestHiggsMass( allLeptons, allMets, allUntagJetsForShapeVar, allBtagJetsForShapeVar, jetUncUntag, jetUncBtag, verbose);


      int h = 0; int a = 0;
      for(std::map<float, int>::iterator it = hMapPt.begin(); it!=hMapPt.end(); it++, h++) hJetRank_[h] = it->second ;
      for(std::map<float, int>::iterator it = aMapPt.begin(); it!=aMapPt.end(); it++, a++) aJetRank_[a] = it->second ;


      //////////////////////////////////////////////////////////////////////////////////////////////
      int all = 0;
      for(std::map<float, int>::iterator it = allMapPt20.begin(); it!=allMapPt20.end(); it++, all++){

	int index = it->second;
	float pt  = -99; float eta = -99; float phi  = -99; float csv  = -99; float mass   = -99; float unc    = -99;
	float topB= -99; float topW= -99; float atopB= -99; float atopW= -99; float higgsB = -99; float flavor = -99;


	// N.B. 
	// ((index_i)<0)*aJet_flavour[(TMath::Max(-1-index_i,0))]+((index_i)>=0)*hJet_flavour[(TMath::Max(index_i,0))]
	// will plot the flavour of the ith jet ordered by pt

	LV genJet;
	if( index >= 0 ){
	  genJet.SetPt(  hJetsgenpt[ index] );
	  genJet.SetEta( hJetsgeneta[index] );
	  genJet.SetPhi( hJetsgenphi[index] );
	  genJet.SetM( 0.0 );
	}else{
	  genJet.SetPt(  aJetsgenpt[ -index-1] );
	  genJet.SetEta( aJetsgeneta[-index-1] );
	  genJet.SetPhi( aJetsgenphi[-index-1] );
	  genJet.SetM( 0.0 );
	}

	int genMatch;
	findGenMatch(genMatch, genJet, topBLV, topW1LV, topW2LV, atopBLV, atopW1LV, atopW2LV, genBLV, genBbarLV);

	float jetMass2 = (index >= 0 ) ?  hJetse[index]*hJetse[index] -  TMath::Power(hJetspt[index]*TMath::CosH(hJetseta[index]) ,2) :
	  aJetse[-index-1]*aJetse[-index-1] -  TMath::Power(aJetspt[-index-1]*TMath::CosH(aJetseta[-index-1]) ,2);

	pt    = (index >= 0 ) ?  hJetspt[index]  : aJetspt[ -index-1];
	eta   = (index >= 0 ) ?  hJetseta[index] : aJetseta[-index-1];
	phi   = (index >= 0 ) ?  hJetsphi[index] : aJetsphi[-index-1];
	csv   = (index >= 0 ) ?  hJetscsv[index] : aJetscsv[-index-1];
	flavor= (index >= 0 ) ?  hJetsflavor[index] : aJetsflavor[-index-1];
	mass  = TMath::Sqrt(jetMass2);
	unc   = (index >= 0 ) ?  hJetsunc[index]    : hJetsunc[-index-1];

	setTopDecay(topB, topW, atopB, atopW, higgsB, genMatch);
       
	if(TMath::Abs(flavor) == 21 || (TMath::Abs(flavor) > 0 && TMath::Abs(flavor)<4) ){
	  nLF_++;
	  if( topW ) nLFTop_++;
	}
	if(TMath::Abs(flavor) == 4 ){
	  nC_++;
	  if( topW ) nCTop_++;
	}
	if(TMath::Abs(flavor) == 5 ){
	  nB_++;
	  if( topB ) nBTop_++;
	}

	if(all==0)      jet1_.set( index, pt, eta, phi, mass, csv, topB, topW, atopB, atopW, higgsB, flavor, unc );
	else if(all==1) jet2_.set( index, pt, eta, phi, mass, csv, topB, topW, atopB, atopW, higgsB, flavor, unc );
	else if(all==2) jet3_.set( index, pt, eta, phi, mass, csv, topB, topW, atopB, atopW, higgsB, flavor, unc );
	else if(all==3) jet4_.set( index, pt, eta, phi, mass, csv, topB, topW, atopB, atopW, higgsB, flavor, unc );
	else if(all==4) jet5_.set( index, pt, eta, phi, mass, csv, topB, topW, atopB, atopW, higgsB, flavor, unc );
	else if(all==5) jet6_.set( index, pt, eta, phi, mass, csv, topB, topW, atopB, atopW, higgsB, flavor, unc );
	else if(all==6) jet7_.set( index, pt, eta, phi, mass, csv, topB, topW, atopB, atopW, higgsB, flavor, unc );
	else if(all==7) jet8_.set( index, pt, eta, phi, mass, csv, topB, topW, atopB, atopW, higgsB, flavor, unc );
	else if(all==8) jet9_.set( index, pt, eta, phi, mass, csv, topB, topW, atopB, atopW, higgsB, flavor, unc );
	else if(all==9) jet10_.set(index, pt, eta, phi, mass, csv, topB, topW, atopB, atopW, higgsB, flavor, unc );
	else{}

      }

      ////////////////////////////////////////////////////////////
      // TTH and TTJets specific
      ////////////////////////////////////////////////////////////

      bool overlap = false;
      if(topHiggsDecay.size()>1){
	for(unsigned int k = 0; k < topHiggsDecay.size()-1; k++){
	  for(unsigned int l = k+1; l < topHiggsDecay.size(); l++){
	    LV first  = topHiggsDecay[k];
	    LV second = topHiggsDecay[l];
	    float deltaR = Geom::deltaR(first,second);
	    if( deltaR < 0.5) overlap = true;
	  }
	}
      }
      bool overlapHiggs = false;
      if(genBLV.Pt()>0 && genBbarLV.Pt()>0){
	for(unsigned int k = 0; k < topHiggsDecay.size(); k++){
	  float deltaR1 = Geom::deltaR(topHiggsDecay[k],genBLV);
	  float deltaR2 = Geom::deltaR(topHiggsDecay[k],genBbarLV);
	  if( deltaR1 < 0.5 && deltaR1 > 0.01) overlapHiggs = true; // match Higgs b with other top decay
	  if( deltaR2 < 0.5 && deltaR2 > 0.01) overlapHiggs = true; // match Higgs bbar with other top decay
	}
      }
      

      recoHiggs_.hdau1mass = higgsB1_flag>0 ? recohiggsB1LV.M()   : -99;
      recoHiggs_.hdau1pt   = higgsB1_flag>0 ? recohiggsB1LV.Pt()  : -99;
      recoHiggs_.hdau1eta  = higgsB1_flag>0 ? recohiggsB1LV.Eta() : -99;
      recoHiggs_.hdau1phi  = higgsB1_flag>0 ? recohiggsB1LV.Phi() : -99;
      recoHiggs_.hdau1id   = higgsB1_flag;
      recoHiggs_.hdau2mass = higgsB2_flag>0 ? recohiggsB2LV.M()   : -99;
      recoHiggs_.hdau2pt   = higgsB2_flag>0 ? recohiggsB2LV.Pt()  : -99;
      recoHiggs_.hdau2eta  = higgsB2_flag>0 ? recohiggsB2LV.Eta() : -99;
      recoHiggs_.hdau2phi  = higgsB2_flag>0 ? recohiggsB2LV.Phi() : -99;
      recoHiggs_.hdau2id   = higgsB2_flag;
      recoHiggs_.hmass     = (higgsB1_flag>0 && higgsB2_flag>0) ? (recohiggsB1LV+recohiggsB2LV).M()   : -99;
      recoHiggs_.hpt       = (higgsB1_flag>0 && higgsB2_flag>0) ? (recohiggsB1LV+recohiggsB2LV).Pt()  : -99;
      recoHiggs_.heta      = (higgsB1_flag>0 && higgsB2_flag>0) ? (recohiggsB1LV+recohiggsB2LV).Eta() : -99;
      recoHiggs_.hphi      = (higgsB1_flag>0 && higgsB2_flag>0) ? (recohiggsB1LV+recohiggsB2LV).Phi() : -99;
      recoHiggs_.hid       = (higgsB1_flag + higgsB2_flag - 2);

      genMatchHiggs_.hdau1mass = higgsB1_flag>0 ? genBLV.M()   : -99;
      genMatchHiggs_.hdau1pt   = higgsB1_flag>0 ? genBLV.Pt()  : -99;
      genMatchHiggs_.hdau1eta  = higgsB1_flag>0 ? genBLV.Eta() : -99;
      genMatchHiggs_.hdau1phi  = higgsB1_flag>0 ? genBLV.Phi() : -99;
      genMatchHiggs_.hdau1id   = higgsB1_flag;
      genMatchHiggs_.hdau2mass = higgsB2_flag>0 ? genBbarLV.M()   : -99;
      genMatchHiggs_.hdau2pt   = higgsB2_flag>0 ? genBbarLV.Pt()  : -99;
      genMatchHiggs_.hdau2eta  = higgsB2_flag>0 ? genBbarLV.Eta() : -99;
      genMatchHiggs_.hdau2phi  = higgsB2_flag>0 ? genBbarLV.Phi() : -99;
      genMatchHiggs_.hdau2id   = higgsB2_flag;
      genMatchHiggs_.hmass     = (higgsB1_flag>0 && higgsB2_flag>0) ? (genBLV+genBbarLV).M()   : -99;
      genMatchHiggs_.hpt       = (higgsB1_flag>0 && higgsB2_flag>0) ? (genBLV+genBbarLV).Pt()  : -99;
      genMatchHiggs_.heta      = (higgsB1_flag>0 && higgsB2_flag>0) ? (genBLV+genBbarLV).Eta() : -99;
      genMatchHiggs_.hphi      = (higgsB1_flag>0 && higgsB2_flag>0) ? (genBLV+genBbarLV).Phi() : -99;
      genMatchHiggs_.hid       = !overlapHiggs ? (higgsB1_flag + higgsB2_flag - 2) : -(higgsB1_flag + higgsB2_flag - 2);



      if( (abs(genTop.wdau1id)<6 && abs(genTop.wdau2id)<6) || (abs(genTop.wdau1id)>=11 && abs(genTop.wdau2id)>=11)
	  ){ // top->Wb->jjb || top->Wb->lnub
	
	recoTop_.bmass     = topB_flag>0  ? recotopBLV.M()   : -99;
	recoTop_.bpt       = topB_flag>0  ? recotopBLV.Pt()  : -99;
	recoTop_.beta      = topB_flag>0  ? recotopBLV.Eta() : -99;
	recoTop_.bphi      = topB_flag>0  ? recotopBLV.Phi() : -99;
	recoTop_.bstatus   = topB_flag;
	recoTop_.wdau1mass = topW1_flag>0 ? recotopW1LV.M()  : -99;
	recoTop_.wdau1pt   = topW1_flag>0 ? recotopW1LV.Pt() : -99;
	recoTop_.wdau1eta  = topW1_flag>0 ? recotopW1LV.Eta(): -99;
	recoTop_.wdau1phi  = topW1_flag>0 ? recotopW1LV.Phi(): -99;
	recoTop_.wdau1id   = topW1_flag;
	recoTop_.wdau2mass = topW2_flag>0 ? recotopW2LV.M()  : -99;
        recoTop_.wdau2pt   = topW2_flag>0 ? recotopW2LV.Pt() : -99;
        recoTop_.wdau2eta  = topW2_flag>0 ? recotopW2LV.Eta(): -99;
        recoTop_.wdau2phi  = topW2_flag>0 ? recotopW2LV.Phi(): -99;
        recoTop_.wdau2id   = topW2_flag;
	recoTop_.wmass     = (topW1_flag>0 && topW2_flag>0) ? (recotopW1LV+recotopW2LV).M()  : -99;
	recoTop_.wpt       = (topW1_flag>0 && topW2_flag>0) ? (recotopW1LV+recotopW2LV).Pt() : -99;
	recoTop_.weta      = (topW1_flag>0 && topW2_flag>0) ? (recotopW1LV+recotopW2LV).Eta(): -99;
	recoTop_.wphi      = (topW1_flag>0 && topW2_flag>0) ? (recotopW1LV+recotopW2LV).Phi(): -99;
	recoTop_.wid       = (topW1_flag + topW2_flag - 1);
	recoTop_.topmass   = (topW1_flag>0 && topW2_flag>0 && topB_flag>0) ? (recotopW1LV+recotopW2LV+recotopBLV).M()  : -99;
        recoTop_.toppt     = (topW1_flag>0 && topW2_flag>0 && topB_flag>0) ? (recotopW1LV+recotopW2LV+recotopBLV).Pt() : -99;
        recoTop_.topeta    = (topW1_flag>0 && topW2_flag>0 && topB_flag>0) ? (recotopW1LV+recotopW2LV+recotopBLV).Eta(): -99;
        recoTop_.topphi    = (topW1_flag>0 && topW2_flag>0 && topB_flag>0) ? (recotopW1LV+recotopW2LV+recotopBLV).Phi(): -99;
        recoTop_.topid     = (topW1_flag + topW2_flag + topB_flag - 2 );
	recoTop_.bbarmass  = topB_flag>0 && atopB_flag>0 ? (recoatopBLV+recotopBLV).M() : -99;
	recoTop_.bbarpt    = topB_flag>0 && atopB_flag>0 ? (recoatopBLV+recotopBLV).Pt() : -99;
	recoTop_.bbareta   = topB_flag>0 && atopB_flag>0 ? (recoatopBLV+recotopBLV).Eta() : -99;
	recoTop_.bbarphi   = topB_flag>0 && atopB_flag>0 ? (recoatopBLV+recotopBLV).Phi() : -99;

	genMatchTop_.bmass     = topB_flag>0  ? topBLV.M()   : -99;
        genMatchTop_.bpt       = topB_flag>0  ? topBLV.Pt()  : -99;
        genMatchTop_.beta      = topB_flag>0  ? topBLV.Eta() : -99;
        genMatchTop_.bphi      = topB_flag>0  ? topBLV.Phi() : -99;
        genMatchTop_.bstatus   = topB_flag;
        genMatchTop_.wdau1mass = topW1_flag>0 ? topW1LV.M()  : -99;
        genMatchTop_.wdau1pt   = topW1_flag>0 ? topW1LV.Pt() : -99;
        genMatchTop_.wdau1eta  = topW1_flag>0 ? topW1LV.Eta(): -99;
        genMatchTop_.wdau1phi  = topW1_flag>0 ? topW1LV.Phi(): -99;
        genMatchTop_.wdau1id   = topW1_flag;
        genMatchTop_.wdau2mass = topW2_flag>0 ? topW2LV.M()  : -99;
        genMatchTop_.wdau2pt   = topW2_flag>0 ? topW2LV.Pt() : -99;
        genMatchTop_.wdau2eta  = topW2_flag>0 ? topW2LV.Eta(): -99;
        genMatchTop_.wdau2phi  = topW2_flag>0 ? topW2LV.Phi(): -99;
        genMatchTop_.wdau2id   = topW2_flag;
        genMatchTop_.wmass     = (topW1_flag>0 && topW2_flag>0) ? (topW1LV+topW2LV).M()  : -99;
        genMatchTop_.wpt       = (topW1_flag>0 && topW2_flag>0) ? (topW1LV+topW2LV).Pt() : -99;
        genMatchTop_.weta      = (topW1_flag>0 && topW2_flag>0) ? (topW1LV+topW2LV).Eta(): -99;
        genMatchTop_.wphi      = (topW1_flag>0 && topW2_flag>0) ? (topW1LV+topW2LV).Phi(): -99;
        genMatchTop_.wid       = (topW1_flag + topW2_flag -1 );
        genMatchTop_.topmass   = (topW1_flag>0 && topW2_flag>0 && topB_flag>0) ? (topW1LV+topW2LV+topBLV).M()  : -99;
	genMatchTop_.toppt     = (topW1_flag>0 && topW2_flag>0 && topB_flag>0) ? (topW1LV+topW2LV+topBLV).Pt() : -99;
        genMatchTop_.topeta    = (topW1_flag>0 && topW2_flag>0 && topB_flag>0) ? (topW1LV+topW2LV+topBLV).Eta(): -99;
        genMatchTop_.topphi    = (topW1_flag>0 && topW2_flag>0 && topB_flag>0) ? (topW1LV+topW2LV+topBLV).Phi(): -99;
        genMatchTop_.topid     = !overlap ? (topW1_flag + topW2_flag + topB_flag - 2 ) : -(topW1_flag + topW2_flag + topB_flag - 2 );
	genMatchTop_.bbarmass      = topB_flag>0 && atopB_flag>0 ? (topBLV+atopBLV).M()   : -99;
	genMatchTop_.bbarpt        = topB_flag>0 && atopB_flag>0 ? (topBLV+atopBLV).Pt()  : -99;
	genMatchTop_.bbareta       = topB_flag>0 && atopB_flag>0 ? (topBLV+atopBLV).Eta() : -99;
	genMatchTop_.bbarphi       = topB_flag>0 && atopB_flag>0 ? (topBLV+atopBLV).Phi() : -99;
      }
      
      if( (abs(genTbar.wdau1id)<6 && abs(genTbar.wdau2id)<6) || (abs(genTbar.wdau1id)>=11 && abs(genTbar.wdau2id)>=11)
          ){ // tbar->Wb->jjb || tbar->Wb->lnub

        recoTbar_.bmass     = atopB_flag>0  ? recoatopBLV.M()   : -99;
        recoTbar_.bpt       = atopB_flag>0  ? recoatopBLV.Pt()  : -99;
        recoTbar_.beta      = atopB_flag>0  ? recoatopBLV.Eta() : -99;
        recoTbar_.bphi      = atopB_flag>0  ? recoatopBLV.Phi() : -99;
        recoTbar_.bstatus   = atopB_flag;
        recoTbar_.wdau1mass = atopW1_flag>0 ? recoatopW1LV.M()  : -99;
        recoTbar_.wdau1pt   = atopW1_flag>0 ? recoatopW1LV.Pt() : -99;
        recoTbar_.wdau1eta  = atopW1_flag>0 ? recoatopW1LV.Eta(): -99;
        recoTbar_.wdau1phi  = atopW1_flag>0 ? recoatopW1LV.Phi(): -99;
        recoTbar_.wdau1id   = atopW1_flag;
        recoTbar_.wdau2mass = atopW2_flag>0 ? recoatopW2LV.M()  : -99;
        recoTbar_.wdau2pt   = atopW2_flag>0 ? recoatopW2LV.Pt() : -99;
        recoTbar_.wdau2eta  = atopW2_flag>0 ? recoatopW2LV.Eta(): -99;
        recoTbar_.wdau2phi  = atopW2_flag>0 ? recoatopW2LV.Phi(): -99;
        recoTbar_.wdau2id   = atopW2_flag;
        recoTbar_.wmass     = (atopW1_flag>0 && atopW2_flag>0) ? (recoatopW1LV+recoatopW2LV).M()  : -99;
        recoTbar_.wpt       = (atopW1_flag>0 && atopW2_flag>0) ? (recoatopW1LV+recoatopW2LV).Pt() : -99;
        recoTbar_.weta      = (atopW1_flag>0 && atopW2_flag>0) ? (recoatopW1LV+recoatopW2LV).Eta(): -99;
        recoTbar_.wphi      = (atopW1_flag>0 && atopW2_flag>0) ? (recoatopW1LV+recoatopW2LV).Phi(): -99;
        recoTbar_.wid       = (atopW1_flag + atopW2_flag - 1);
        recoTbar_.topmass   = (atopW1_flag>0 && atopW2_flag>0 && atopB_flag>0) ? (recoatopW1LV+recoatopW2LV+recoatopBLV).M()  : -99;
        recoTbar_.toppt     = (atopW1_flag>0 && atopW2_flag>0 && atopB_flag>0) ? (recoatopW1LV+recoatopW2LV+recoatopBLV).Pt() : -99;
        recoTbar_.topeta    = (atopW1_flag>0 && atopW2_flag>0 && atopB_flag>0) ? (recoatopW1LV+recoatopW2LV+recoatopBLV).Eta(): -99;
        recoTbar_.topphi    = (atopW1_flag>0 && atopW2_flag>0 && atopB_flag>0) ? (recoatopW1LV+recoatopW2LV+recoatopBLV).Phi(): -99;
	recoTbar_.bbarmass  = topB_flag>0 && atopB_flag>0 ? (recoatopBLV+recotopBLV).M()   : -99;
	recoTbar_.bbarpt    = topB_flag>0 && atopB_flag>0 ? (recoatopBLV+recotopBLV).Pt()  : -99;
	recoTbar_.bbareta   = topB_flag>0 && atopB_flag>0 ? (recoatopBLV+recotopBLV).Eta() : -99;
	recoTbar_.bbarphi   = topB_flag>0 && atopB_flag>0 ? (recoatopBLV+recotopBLV).Phi() : -99;

	genMatchTbar_.bmass     = atopB_flag>0  ? atopBLV.M()   : -99;
        genMatchTbar_.bpt       = atopB_flag>0  ? atopBLV.Pt()  : -99;
        genMatchTbar_.beta      = atopB_flag>0  ? atopBLV.Eta() : -99;
        genMatchTbar_.bphi      = atopB_flag>0  ? atopBLV.Phi() : -99;
        genMatchTbar_.bstatus   = atopB_flag;
        genMatchTbar_.wdau1mass = atopW1_flag>0 ? atopW1LV.M()  : -99;
        genMatchTbar_.wdau1pt   = atopW1_flag>0 ? atopW1LV.Pt() : -99;
        genMatchTbar_.wdau1eta  = atopW1_flag>0 ? atopW1LV.Eta(): -99;
        genMatchTbar_.wdau1phi  = atopW1_flag>0 ? atopW1LV.Phi(): -99;
        genMatchTbar_.wdau1id   = atopW1_flag;
        genMatchTbar_.wdau2mass = atopW2_flag>0 ? atopW2LV.M()  : -99;
        genMatchTbar_.wdau2pt   = atopW2_flag>0 ? atopW2LV.Pt() : -99;
        genMatchTbar_.wdau2eta  = atopW2_flag>0 ? atopW2LV.Eta(): -99;
        genMatchTbar_.wdau2phi  = atopW2_flag>0 ? atopW2LV.Phi(): -99;
        genMatchTbar_.wdau2id   = atopW2_flag;
        genMatchTbar_.wmass     = (atopW1_flag>0 && atopW2_flag>0) ? (atopW1LV+atopW2LV).M()  : -99;
        genMatchTbar_.wpt       = (atopW1_flag>0 && atopW2_flag>0) ? (atopW1LV+atopW2LV).Pt() : -99;
        genMatchTbar_.weta      = (atopW1_flag>0 && atopW2_flag>0) ? (atopW1LV+atopW2LV).Eta(): -99;
        genMatchTbar_.wphi      = (atopW1_flag>0 && atopW2_flag>0) ? (atopW1LV+atopW2LV).Phi(): -99;
        genMatchTbar_.wid       = (atopW1_flag + atopW2_flag -1 );
        genMatchTbar_.topmass   = (atopW1_flag>0 && atopW2_flag>0 && atopB_flag>0) ? (atopW1LV+atopW2LV+atopBLV).M()  : -99;
        genMatchTbar_.toppt     = (atopW1_flag>0 && atopW2_flag>0 && atopB_flag>0) ? (atopW1LV+atopW2LV+atopBLV).Pt() : -99;
        genMatchTbar_.topeta    = (atopW1_flag>0 && atopW2_flag>0 && atopB_flag>0) ? (atopW1LV+atopW2LV+atopBLV).Eta(): -99;
        genMatchTbar_.topphi    = (atopW1_flag>0 && atopW2_flag>0 && atopB_flag>0) ? (atopW1LV+atopW2LV+atopBLV).Phi(): -99;
        genMatchTbar_.topid     = !overlap ? (atopW1_flag + atopW2_flag + atopB_flag - 2 ) : -(atopW1_flag + atopW2_flag + atopB_flag - 2 );
	genMatchTbar_.bbarmass  = topB_flag>0 && atopB_flag>0 ? (topBLV+atopBLV).M()   : -99;
	genMatchTbar_.bbarpt    = topB_flag>0 && atopB_flag>0 ? (topBLV+atopBLV).Pt()  : -99;
	genMatchTbar_.bbareta   = topB_flag>0 && atopB_flag>0 ? (topBLV+atopBLV).Eta() : -99;
	genMatchTbar_.bbarphi   = topB_flag>0 && atopB_flag>0 ? (topBLV+atopBLV).Phi() : -99;
      }


      if(eventsToDebug>0 && eventsToDebugCounter<eventsToDebug){

	bool isTopToHad   = abs(genTop.wdau1id)<6 && abs(genTop.wdau2id)<6 && abs(genTbar.wdau1id)>=11;
	bool isTopInAcc   = genTop.wdau1pt>30 && genTop.wdau2pt>30 && genTop.bpt>30 && genTbar.bpt>30  && TMath::Abs(genTop.wdau1eta)<5.0 && TMath::Abs(genTop.wdau2eta)<5.0 && TMath::Abs(genTop.beta)<5.0 && TMath::Abs(genTbar.beta)<5.0;
	bool isHiggsInAcc = genB.pt>30 && genBbar.pt>30 && abs(genB.eta)<5.0 && TMath::Abs(genBbar.eta)<5.0;

	if(isTopToHad && isTopInAcc && isHiggsInAcc){

	  eventsToDebugCounter++;

	  cout << "************ DEBUG" << endl;
	  cout << "In event " << ev.event << ", the hadronic top and the higgs decay in acceptance...number of reconstructed jets: " << numJets30_ << endl;
	  cout << "The partons are " << endl;
	  cout << string(Form("W(1). (%.0f,%.1f,%.1f)",genTop.wdau1pt,genTop.wdau1eta,genTop.wdau1phi)) << endl;
	  cout << string(Form("W(2). (%.0f,%.1f,%.1f)",genTop.wdau2pt,genTop.wdau2eta,genTop.wdau2phi)) << endl;
	  cout << string(Form("b.    (%.0f,%.1f,%.1f)",genTop.bpt,genTop.beta,genTop.bphi)) << endl;
	  cout << string(Form("bbar. (%.0f,%.1f,%.1f)",genTbar.bpt,genTbar.beta,genTbar.bphi)) << endl;
	  cout << string(Form("H(1). (%.0f,%.1f,%.1f)",genB.pt,genB.eta,genB.phi)) << endl;
	  cout << string(Form("H(2). (%.0f,%.1f,%.1f)",genBbar.pt,genBbar.eta,genBbar.phi)) << endl;

	  float minDeltaR = 999.;
	  if(topHiggsDecay.size()>1){
	    for(unsigned int k = 0; k < topHiggsDecay.size()-1; k++){
	      for(unsigned int l = k+1; l < topHiggsDecay.size(); l++){
		LV first  = topHiggsDecay[k];
		LV second = topHiggsDecay[l];
		float deltaR = Geom::deltaR(first,second);
		if( deltaR < minDeltaR) minDeltaR = deltaR;
	      }
	    }
	  }
	  cout << " ==> The smallest deltaR between all partons is " << string(Form("%.1f",minDeltaR)) << endl;

	  if(numJets30_>=5 && numJets30_<=7){
	    cout << numJets30_ << " jets have been reconstructed: " << endl;
	    cout << string(Form("1. (%.0f,%.1f,%.1f)",jet1_.pt,jet1_.eta,jet1_.phi)) << endl;
	    cout << string(Form("2. (%.0f,%.1f,%.1f)",jet2_.pt,jet2_.eta,jet2_.phi)) << endl;
	    cout << string(Form("3. (%.0f,%.1f,%.1f)",jet3_.pt,jet3_.eta,jet3_.phi)) << endl;
	    cout << string(Form("4. (%.0f,%.1f,%.1f)",jet4_.pt,jet4_.eta,jet4_.phi)) << endl;
	    cout << string(Form("5. (%.0f,%.1f,%.1f)",jet5_.pt,jet5_.eta,jet5_.phi)) << endl;
	    if(numJets30_>=6)  cout << string(Form("6. (%.0f,%.1f,%.1f)",jet6_.pt,jet7_.eta,jet6_.phi)) << endl;
	    if(numJets30_==7)  cout << string(Form("7. (%.0f,%.1f,%.1f)",jet7_.pt,jet6_.eta,jet7_.phi)) << endl;

	    int topB_hits   = 0;
	    int topW1_hits  = 0;
	    int topW2_hits  = 0;
	    int atopB_hits  = 0;
	    int higgsB1_hits= 0;
	    int higgsB2_hits= 0;

	    vector<LV> myJets; 
	    myJets.push_back(LV(jet1_.pt,jet1_.eta,jet1_.phi,jet1_.mass));
	    myJets.push_back(LV(jet2_.pt,jet2_.eta,jet2_.phi,jet2_.mass));
	    myJets.push_back(LV(jet3_.pt,jet3_.eta,jet3_.phi,jet3_.mass));
	    myJets.push_back(LV(jet4_.pt,jet4_.eta,jet4_.phi,jet4_.mass));
	    myJets.push_back(LV(jet5_.pt,jet5_.eta,jet5_.phi,jet5_.mass));
	    if(numJets30_>=6)  myJets.push_back(LV(jet6_.pt,jet6_.eta,jet6_.phi,jet6_.mass));
	    if(numJets30_==7)  myJets.push_back(LV(jet7_.pt,jet7_.eta,jet7_.phi,jet7_.mass));

	    map<unsigned int, int> recoJetHits;
	    for(unsigned int k = 0 ; k <myJets.size() ; k++  ) recoJetHits[k] = 0;

	    for(unsigned int k = 0 ; k <myJets.size() ; k++  ){
	      if(Geom::deltaR( myJets[k],topW1LV ) < 0.5 ){
		cout << "W(1) is dR = " << string(Form("%.1f",Geom::deltaR( myJets[k],topW1LV ) )) << " away from jet #" << k+1 << endl;
		topW1_hits++;
		recoJetHits[k]++;
	      }
	      if(Geom::deltaR( myJets[k],topW2LV ) < 0.5 ){
		cout << "W(2) is dR = " << string(Form("%.1f",Geom::deltaR( myJets[k],topW2LV ) )) << " away from jet #" << k+1 << endl;
		topW2_hits++;
		recoJetHits[k]++;
	      }
	      if(Geom::deltaR( myJets[k],topBLV ) < 0.5 ){
		cout << "b    is dR = " << string(Form("%.1f",Geom::deltaR( myJets[k],topBLV ) )) << " away from jet #" << k+1 << endl;
		topB_hits++;
		recoJetHits[k]++;
	      }
	      if(Geom::deltaR( myJets[k],atopBLV ) < 0.5 ){
		cout << "bbar is dR = " << string(Form("%.1f",Geom::deltaR( myJets[k],atopBLV ) )) << " away from jet #" << k+1 << endl;
		atopB_hits++;
		recoJetHits[k]++;
	      }
	      if(Geom::deltaR( myJets[k],genBLV ) < 0.5 ){
		cout << "H(1) is dR = " << string(Form("%.1f",Geom::deltaR( myJets[k], genBLV) )) << " away from jet #" << k+1 << endl;
		higgsB1_hits++;
		recoJetHits[k]++;
	      }
	      if(Geom::deltaR( myJets[k],genBbarLV ) < 0.5 ){
		cout << "H(2) is dR = " << string(Form("%.1f",Geom::deltaR( myJets[k], genBbarLV) )) << " away from jet #" << k+1 << endl;
		higgsB2_hits++;
		recoJetHits[k]++;
	      }
	    }
	    cout << "Gen Summary:" << endl;
	    cout << "W(1) has " << topW1_hits << "matches" << endl;
	    if( topW1_hits==0){
	      cout << "... looking at all jets>20 GeV..." << endl;
	      for(int i = 0; i < nhJets; i++){
		float jetMass2 = hJetse[i]*hJetse[i] -  TMath::Power(hJetspt[i]*TMath::CosH(hJetseta[i]) ,2);
		LV jetLV(hJetspt[i], hJetseta[i], hJetsphi[i], TMath::Sqrt(jetMass2));
		if(Geom::deltaR( jetLV, topW1LV) < 0.5 ){
		  cout << "W(1) can be matched to this jet: " << string(Form("(%.0f,%.1f,%.1f)",jetLV.Pt(),jetLV.Eta(),jetLV.Phi())) <<  endl;
		}
	      }
	      for(int i = 0; i < naJets; i++){
		float jetMass2 = aJetse[i]*aJetse[i] -  TMath::Power(aJetspt[i]*TMath::CosH(aJetseta[i]) ,2);
		LV jetLV(aJetspt[i], aJetseta[i], aJetsphi[i], TMath::Sqrt(jetMass2));
		if(Geom::deltaR( jetLV, topW1LV) < 0.5 ){
		  cout << "W(1) can be matched to this jet: " << string(Form("(%.0f,%.1f,%.1f)",jetLV.Pt(),jetLV.Eta(),jetLV.Phi())) << endl;
		}
	      }
	    }
	    cout << "W(2) has " << topW2_hits << "matches" << endl;
	    if( topW2_hits==0){
	      cout << "... looking at all jets>20 GeV..." << endl;
	      for(int i = 0; i < nhJets; i++){
		float jetMass2 = hJetse[i]*hJetse[i] -  TMath::Power(hJetspt[i]*TMath::CosH(hJetseta[i]) ,2);
		LV jetLV(hJetspt[i], hJetseta[i], hJetsphi[i], TMath::Sqrt(jetMass2));
		if(Geom::deltaR( jetLV, topW2LV) < 0.5 ){
		  cout << "W(2) can be matched to this jet: " << string(Form("(%.0f,%.1f,%.1f)",jetLV.Pt(),jetLV.Eta(),jetLV.Phi())) << endl;
		}
	      }
	      for(int i = 0; i < naJets; i++){
		float jetMass2 = aJetse[i]*aJetse[i] -  TMath::Power(aJetspt[i]*TMath::CosH(aJetseta[i]) ,2);
		LV jetLV(aJetspt[i], aJetseta[i], aJetsphi[i], TMath::Sqrt(jetMass2));
		if(Geom::deltaR( jetLV, topW2LV) < 0.5 ){
		  cout << "W(2) can be matched to this jet: " << string(Form("(%.0f,%.1f,%.1f)",jetLV.Pt(),jetLV.Eta(),jetLV.Phi())) << endl;
		}
	      }
	    }
	    cout << "b    has " << topB_hits << "matches" << endl;
	    if( topB_hits==0){
	      cout << "... looking at all jets>20 GeV..." << endl;
	      for(int i = 0; i < nhJets; i++){
		float jetMass2 = hJetse[i]*hJetse[i] -  TMath::Power(hJetspt[i]*TMath::CosH(hJetseta[i]) ,2);
		LV jetLV(hJetspt[i], hJetseta[i], hJetsphi[i], TMath::Sqrt(jetMass2));
		if(Geom::deltaR( jetLV, topBLV) < 0.5 ){
		  cout << "b can be matched to this jet: " << string(Form("(%.0f,%.1f,%.1f)",jetLV.Pt(),jetLV.Eta(),jetLV.Phi())) << endl;
		}
	      }
	      for(int i = 0; i < naJets; i++){
		float jetMass2 = aJetse[i]*aJetse[i] -  TMath::Power(aJetspt[i]*TMath::CosH(aJetseta[i]) ,2);
		LV jetLV(aJetspt[i], aJetseta[i], aJetsphi[i], TMath::Sqrt(jetMass2));
		if(Geom::deltaR( jetLV, topBLV) < 0.5 ){
		  cout << "b can be matched to this jet: " << string(Form("(%.0f,%.1f,%.1f)",jetLV.Pt(),jetLV.Eta(),jetLV.Phi())) << endl;
		}
	      }
	    }
	    cout << "bbar has " << atopB_hits << "matches" << endl;
	    if( atopB_hits==0){
	      cout << "... looking at all jets>20 GeV..." << endl;
	      for(int i = 0; i < nhJets; i++){
		float jetMass2 = hJetse[i]*hJetse[i] -  TMath::Power(hJetspt[i]*TMath::CosH(hJetseta[i]) ,2);
		LV jetLV(hJetspt[i], hJetseta[i], hJetsphi[i], TMath::Sqrt(jetMass2));
		if(Geom::deltaR( jetLV, atopBLV) < 0.5 ){
		  cout << "bbar can be matched to this jet: " << string(Form("(%.0f,%.1f,%.1f)",jetLV.Pt(),jetLV.Eta(),jetLV.Phi())) << endl;
		}
	      }
	      for(int i = 0; i < naJets; i++){
		float jetMass2 = aJetse[i]*aJetse[i] -  TMath::Power(aJetspt[i]*TMath::CosH(aJetseta[i]) ,2);
		LV jetLV(aJetspt[i], aJetseta[i], aJetsphi[i], TMath::Sqrt(jetMass2));
		if(Geom::deltaR( jetLV, atopBLV) < 0.5 ){
		  cout << "bbar can be matched to this jet: " << string(Form("(%.0f,%.1f,%.1f)",jetLV.Pt(),jetLV.Eta(),jetLV.Phi())) << endl;
		}
	      }
	    }
	    cout << "H(1) has " << higgsB1_hits << "matches" << endl;
	    if( higgsB1_hits==0){
	      cout << "... looking at all jets>20 GeV..." << endl;
	      for(int i = 0; i < nhJets; i++){
		float jetMass2 = hJetse[i]*hJetse[i] -  TMath::Power(hJetspt[i]*TMath::CosH(hJetseta[i]) ,2);
		LV jetLV(hJetspt[i], hJetseta[i], hJetsphi[i], TMath::Sqrt(jetMass2));
		if(Geom::deltaR( jetLV, genBLV) < 0.5 ){
		  cout << "H(1) can be matched to this jet: " << string(Form("(%.0f,%.1f,%.1f)",jetLV.Pt(),jetLV.Eta(),jetLV.Phi())) << endl;
		}
	      }
	      for(int i = 0; i < naJets; i++){
		float jetMass2 = aJetse[i]*aJetse[i] -  TMath::Power(aJetspt[i]*TMath::CosH(aJetseta[i]) ,2);
		LV jetLV(aJetspt[i], aJetseta[i], aJetsphi[i], TMath::Sqrt(jetMass2));
		if(Geom::deltaR( jetLV, genBLV) < 0.5 ){
		  cout << "H(1) can be matched to this jet: " << string(Form("(%.0f,%.1f,%.1f)",jetLV.Pt(),jetLV.Eta(),jetLV.Phi())) << endl;
		}
	      }
	    }
	    cout << "H(2) has " << higgsB2_hits << "matches" << endl;
	    if( higgsB2_hits==0){
	      cout << "... looking at all jets>20 GeV..." << endl;
	      for(int i = 0; i < nhJets; i++){
		float jetMass2 = hJetse[i]*hJetse[i] -  TMath::Power(hJetspt[i]*TMath::CosH(hJetseta[i]) ,2);
		LV jetLV(hJetspt[i], hJetseta[i], hJetsphi[i], TMath::Sqrt(jetMass2));
		if(Geom::deltaR( jetLV, genBbarLV) < 0.5 ){
		  cout << "H(2) can be matched to this jet: " << string(Form("(%.0f,%.1f,%.1f)",jetLV.Pt(),jetLV.Eta(),jetLV.Phi())) << endl;
		}
	      }
	      for(int i = 0; i < naJets; i++){
		float jetMass2 = aJetse[i]*aJetse[i] -  TMath::Power(aJetspt[i]*TMath::CosH(aJetseta[i]) ,2);
		LV jetLV(aJetspt[i], aJetseta[i], aJetsphi[i], TMath::Sqrt(jetMass2));
		if(Geom::deltaR( jetLV, genBbarLV) < 0.5 ){
		  cout << "H(2) can be matched to this jet: " << string(Form("(%.0f,%.1f,%.1f)",jetLV.Pt(),jetLV.Eta(),jetLV.Phi())) << endl;
		}
	      }
	    }
	    cout << "Reco Summary:" << endl;
	    for(unsigned int k = 0 ; k <myJets.size() ; k++ ){
	      cout << "Jet #" << k+1 << "has " << (recoJetHits[k]) << " matches";
	      int hitsk = recoJetHits[k];
	      if( (hitsk == 0) ){
		cout << " ==> this jet is not matched to any parton... flavor ";
		if(k==0) cout << jet1_.flavor << endl;
		if(k==1) cout << jet2_.flavor << endl;
		if(k==2) cout << jet3_.flavor << endl;
		if(k==3) cout << jet4_.flavor << endl;
		if(k==4) cout << jet5_.flavor << endl;
		if(k==5) cout << jet6_.flavor << endl;
		if(k==6) cout << jet7_.flavor << endl;
	      }
	      if(  (hitsk == 1)  ) cout << " ==> 1 <-> 1 match" << endl;
	      if(  (hitsk == 2)  ) cout << " ==> probably this jet merges two partons" << endl;
	    }
	  } // numJets30_
	}
      }



      for(unsigned k = 0; k< csvBtag.size() ; k++){
	bTagMap30[ csvBtag[k] ] = k;
      }
      int count=0;
      for(std::map<float, int>::iterator it = bTagMap30.begin(); it!=bTagMap30.end(); it++, count++){
	if(count==0) firstBtag_  = it->first;
	if(count==1) secondBtag_ = it->first;
	if(count==2) thirdBtag_  = it->first;
	if(count==3) fourthBtag_ = it->first;
      }


      EventShapeVariables* eventShapes = new EventShapeVariables(allJetsForShapeVar);
      isotropy_   = eventShapes->isotropy();
      circularity_= eventShapes->circularity();
      sphericity_ = eventShapes->sphericity();
      aplanarity_ = eventShapes->aplanarity();
      Cparam_     = eventShapes->C();
      Dparam_     = eventShapes->D();

      TopologyWorker* topologicalWorker = new  TopologyWorker(false);

      if(jetArrayForShapeVar.GetEntries()>0){
	leptArrayForShapeVar.Clear();
	//topologicalWorker->setVerbose(true);
	topologicalWorker->setPartList(&jetArrayForShapeVar, &leptArrayForShapeVar);
	thrust0_ = (topologicalWorker->thrust())[0] ;
	thrust1_ = (topologicalWorker->thrust())[1] ;
	thrust2_ = (topologicalWorker->thrust())[2] ;
	sphericity2_ = topologicalWorker->get_sphericity();
	aplanarity2_ = topologicalWorker->get_aplanarity();
	h10_ = topologicalWorker->get_h10();
	h20_ = topologicalWorker->get_h20();
	h30_ = topologicalWorker->get_h30();
	h40_ = topologicalWorker->get_h40();
	h50_ = topologicalWorker->get_h50();
	h60_ = topologicalWorker->get_h60();
      } 
      jetArrayForShapeVar.Clear(); 
      leptArrayForShapeVar.Clear();



      bool isInJson = jsonContainsEvent(jsonVector, ev.run, ev.lumi);
      myJson_ = isInJson ? 1 : 0;
      //if(!isInJson && verbose) cout << "Event " << run << ", " << lumi << " is not in json" << endl;






      //////////////////////////////////FILL////////////////////////////////////////
      //////////////////////////////////////////////////////////////////////////////

      hJetRankBR->Fill();
      aJetRankBR->Fill();
      numJets40BR->Fill();
      numJets30BR->Fill();
      numJets20BR->Fill();
      numJets40bTagBR->Fill();
      numJets30bTagBR->Fill();
      numJets20bTagBR->Fill();

      jet1BR->Fill();
      jet2BR->Fill();
      jet3BR->Fill();
      jet4BR->Fill();
      jet5BR->Fill();
      jet6BR->Fill();
      jet7BR->Fill();
      jet8BR->Fill();
      jet9BR->Fill();
      jet10BR->Fill();

      nLFBR->Fill();
      nCBR->Fill();
      nBBR->Fill();
      nLFTopBR->Fill();
      nCTopBR->Fill();
      nBTopBR->Fill();

      isotropyBR->Fill();
      circularityBR->Fill();
      sphericityBR->Fill();
      aplanarityBR->Fill();
      CparamBR->Fill();
      DparamBR->Fill();
      thrust0BR->Fill();
      thrust1BR->Fill();
      thrust2BR->Fill();
      sphericity2BR->Fill();
      aplanarity2BR->Fill();
      h10BR->Fill();
      h20BR->Fill();
      h30BR->Fill();
      h40BR->Fill();
      h50BR->Fill();
      h60BR->Fill();

      aveCsvBR->Fill();
      aveDeltaRbTagBR->Fill();
      aveMbTagBR->Fill();
      aveMunTagBR->Fill();
      closestJJbTagMassBR->Fill();
      varCsvBR->Fill();
      maxCsvBR->Fill();
      minCsvBR->Fill();
      sumAllPtBR->Fill();
      sumAllJet20PtBR->Fill();
      sumAllJet30PtBR->Fill();
      sumAllJet40PtBR->Fill();
      massAllBR->Fill();
      massLJBR->Fill();
      M3BR->Fill();
      MHT20BR->Fill();
      MHT30BR->Fill();
      MHT40BR->Fill();

      minDeltaRLJBR->Fill();
      minDeltaRBtagBR->Fill();

      firstBtagBR->Fill();
      secondBtagBR->Fill();
      thirdBtagBR->Fill();
      fourthBtagBR->Fill();

      bestHiggsMassBR->Fill();

      myJsonBR->Fill();
      numOfBsBR->Fill();
      numOfBsAccBR->Fill();
      numOfBsFlavBR->Fill();
      numOfBsFlavAccBR->Fill();
      numOfCsFlavBR->Fill();
      numOfCsFlavAccBR->Fill();

      recoTopBR->Fill();
      genMatchTopBR->Fill();
      recoTbarBR->Fill();
      genMatchTbarBR->Fill();
      recoHiggsBR->Fill();
      genMatchHiggsBR->Fill();

      delete eventShapes;
      delete topologicalWorker;

    }

    fs->cd();
    outTree->Write("",TObject::kOverwrite );
    fs->Close();

    mySamples->GetFile( currentName )->Close();

   }



  return 0;


}
