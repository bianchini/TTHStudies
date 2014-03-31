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


// -99 => jet not matched to partons
// -1  => jet not matched to top or higgs decay products
//  0  => jet matched to higgs b quarks
//  1  => jet matched to top b quarks
//  2  => jet matched to top W quarks

void findGenMatch(int& genMatch, LV genJet, LV topBLV, LV topW1LV, LV topW2LV, LV atopBLV, LV atopW1LV, LV atopW2LV, LV genBLV, LV genBbarLV){

  if(genJet.Pt()<=0.001){
    genMatch = -99;
    return;
  }

  if( topBLV.Pt()>0 && Geom::deltaR(genJet, topBLV)            < GENJETDR ){
    genMatch = 1;
  }
  else if( topW1LV.Pt()>0 && Geom::deltaR(genJet, topW1LV)     < GENJETDR ){
    genMatch = 2;
  }
  else if( topW2LV.Pt()>0 && Geom::deltaR(genJet, topW2LV)     < GENJETDR ){
    genMatch = 2;
  }
  else if( atopBLV.Pt()>0 && Geom::deltaR(genJet, atopBLV)     < GENJETDR ){
    genMatch = 1;
  }
  else if( atopW1LV.Pt()>0 && Geom::deltaR(genJet, atopW1LV)   < GENJETDR ){
    genMatch = 2;
  }
  else if( atopW2LV.Pt()>0 && Geom::deltaR(genJet, atopW2LV)   < GENJETDR ){
    genMatch = 2;
  }
  else if( genBLV.Pt()>0 && Geom::deltaR(genJet, genBLV)       < GENJETDR ){
    genMatch = 0;
  }
  else if( genBbarLV.Pt()>0 && Geom::deltaR(genJet, genBbarLV) < GENJETDR ){
    genMatch = 0;
  }
  else{
    genMatch = -1; 
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

  std::cout << "Step3" << std::endl;
  gROOT->SetBatch(true);
 
  gSystem->Load("libFWCoreFWLite");
  gSystem->Load("libDataFormatsFWLite");

  AutoLibraryLoader::enable();

  PythonProcessDesc builder(argv[1]);
 
  const edm::ParameterSet& in = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteInput");

  std::string pathToFile( in.getParameter<std::string>("pathToFile" ) );
  std::string ordering(   in.getParameter<std::string>("ordering" ) );
  double lumi(            in.getParameter<double>("lumi" ) );
  bool verbose(           in.getParameter<bool>("verbose" ) );


  std::vector<edm::LuminosityBlockRange> jsonVector;
  if ( in.exists("lumisToProcess") ) 
    {
      std::vector<edm::LuminosityBlockRange> const & lumisTemp =
	in.getUntrackedParameter<std::vector<edm::LuminosityBlockRange> > ("lumisToProcess");
      jsonVector.resize( lumisTemp.size() );
      copy( lumisTemp.begin(), lumisTemp.end(), jsonVector.begin() );
    }
  

  //const edm::VParameterSet& samples = in.getParameter<edm::VParameterSet>("samples") ;
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
      pset.addParameter("update",true);
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
    
    //fwlite::TFileService fs = fwlite::TFileService( ("TThNTuples_"+currentName+".root").c_str());
    //TTree* outTree = fs.make<TTree>("outTree","TTH Tree");
    //TIter nextkey((mySamples->GetFile( currentName ))->GetListOfKeys());
    //TH1F *key;
    //while ( (key = (TH1F*)nextkey()) ) {
    //string name(key->GetName());
    //if( name.find("tree")!=string::npos) continue;
    //TH1F* h = (TH1F*)(mySamples->GetFile( currentName ))->FindObjectAny( key->GetName());
    //h->Write(key->GetName());
    //}
    //outTree = (mySamples->GetTree( currentName, "tree"))->CopyTree("","");

    mySamples->OpenFile( currentName );
    cout << "Opening file " << currentName << endl;
    TTree* outTree = mySamples->GetTree( currentName, "tree");
    cout << "Done!!" << endl;

    //////////////////////////////////NEW VARIABLES///////////////////////////////
    //////////////////////////////////////////////////////////////////////////////

    int hJetRank_[100];
    int aJetRank_[100];
    int index1_,  index2_,  index3_,  index4_,  index5_,  index6_;
    float pt1_,  pt2_,  pt3_,  pt4_,  pt5_,  pt6_;
    float eta1_, eta2_, eta3_, eta4_, eta5_, eta6_;
    float phi1_, phi2_, phi3_, phi4_, phi5_, phi6_;
    float csv1_, csv2_, csv3_, csv4_, csv5_, csv6_;
    float topB1_,  topB2_,  topB3_,  topB4_,  topB5_,  topB6_;
    float topW1_,  topW2_,  topW3_,  topW4_,  topW5_,  topW6_;
    float higgsB1_,  higgsB2_,  higgsB3_,  higgsB4_,  higgsB5_,  higgsB6_;
    float flavor1_, flavor2_, flavor3_, flavor4_, flavor5_, flavor6_;
    int nLF_,    nC_,    nB_;
    int nLFTop_, nCTop_, nBTop_;

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

    TBranch *hJetRankBR = outTree->Branch("hJetRank",hJetRank_,"hJetRank[nhJets]/I");
    TBranch *aJetRankBR = outTree->Branch("aJetRank",aJetRank_,"aJetRank[naJets]/I");

    TBranch *index1BR = outTree->Branch("index1",&index1_,"index1/I");
    TBranch *index2BR = outTree->Branch("index2",&index2_,"index2/I");
    TBranch *index3BR = outTree->Branch("index3",&index3_,"index3/I");
    TBranch *index4BR = outTree->Branch("index4",&index4_,"index4/I");
    TBranch *index5BR = outTree->Branch("index5",&index5_,"index5/I");
    TBranch *index6BR = outTree->Branch("index6",&index6_,"index6/I");

    TBranch *pt1BR = outTree->Branch("pt1",&pt1_,"pt1/F");
    TBranch *pt2BR = outTree->Branch("pt2",&pt2_,"pt2/F");
    TBranch *pt3BR = outTree->Branch("pt3",&pt3_,"pt3/F");
    TBranch *pt4BR = outTree->Branch("pt4",&pt4_,"pt4/F");
    TBranch *pt5BR = outTree->Branch("pt5",&pt5_,"pt5/F");
    TBranch *pt6BR = outTree->Branch("pt6",&pt6_,"pt6/F");

    TBranch *eta1BR = outTree->Branch("eta1",&eta1_,"eta1/F");
    TBranch *eta2BR = outTree->Branch("eta2",&eta2_,"eta2/F");
    TBranch *eta3BR = outTree->Branch("eta3",&eta3_,"eta3/F");
    TBranch *eta4BR = outTree->Branch("eta4",&eta4_,"eta4/F");
    TBranch *eta5BR = outTree->Branch("eta5",&eta5_,"eta5/F");
    TBranch *eta6BR = outTree->Branch("eta6",&eta6_,"eta6/F");
 
    TBranch *phi1BR = outTree->Branch("phi1",&phi1_,"phi1/F");
    TBranch *phi2BR = outTree->Branch("phi2",&phi2_,"phi2/F");
    TBranch *phi3BR = outTree->Branch("phi3",&phi3_,"phi3/F");
    TBranch *phi4BR = outTree->Branch("phi4",&phi4_,"phi4/F");
    TBranch *phi5BR = outTree->Branch("phi5",&phi5_,"phi5/F");
    TBranch *phi6BR = outTree->Branch("phi6",&phi6_,"phi6/F");

    TBranch *csv1BR = outTree->Branch("csv1",&csv1_,"csv1/F");
    TBranch *csv2BR = outTree->Branch("csv2",&csv2_,"csv2/F");
    TBranch *csv3BR = outTree->Branch("csv3",&csv3_,"csv3/F");
    TBranch *csv4BR = outTree->Branch("csv4",&csv4_,"csv4/F");
    TBranch *csv5BR = outTree->Branch("csv5",&csv5_,"csv5/F");
    TBranch *csv6BR = outTree->Branch("csv6",&csv6_,"csv6/F");

    TBranch *topB1BR = outTree->Branch("topB1",&topB1_,"topB1/F");
    TBranch *topB2BR = outTree->Branch("topB2",&topB2_,"topB2/F");
    TBranch *topB3BR = outTree->Branch("topB3",&topB3_,"topB3/F");
    TBranch *topB4BR = outTree->Branch("topB4",&topB4_,"topB4/F");
    TBranch *topB5BR = outTree->Branch("topB5",&topB5_,"topB5/F");
    TBranch *topB6BR = outTree->Branch("topB6",&topB6_,"topB6/F");

    TBranch *topW1BR = outTree->Branch("topW1",&topW1_,"topW1/F");
    TBranch *topW2BR = outTree->Branch("topW2",&topW2_,"topW2/F");
    TBranch *topW3BR = outTree->Branch("topW3",&topW3_,"topW3/F");
    TBranch *topW4BR = outTree->Branch("topW4",&topW4_,"topW4/F");
    TBranch *topW5BR = outTree->Branch("topW5",&topW5_,"topW5/F");
    TBranch *topW6BR = outTree->Branch("topW6",&topW6_,"topW6/F");

    TBranch *higgsB1BR = outTree->Branch("higgsB1",&higgsB1_,"higgsB1/F");
    TBranch *higgsB2BR = outTree->Branch("higgsB2",&higgsB2_,"higgsB2/F");
    TBranch *higgsB3BR = outTree->Branch("higgsB3",&higgsB3_,"higgsB3/F");
    TBranch *higgsB4BR = outTree->Branch("higgsB4",&higgsB4_,"higgsB4/F");
    TBranch *higgsB5BR = outTree->Branch("higgsB5",&higgsB5_,"higgsB5/F");
    TBranch *higgsB6BR = outTree->Branch("higgsB6",&higgsB6_,"higgsB6/F");

    TBranch *flavor1BR = outTree->Branch("flavor1",&flavor1_,"flavor1/F");
    TBranch *flavor2BR = outTree->Branch("flavor2",&flavor2_,"flavor2/F");
    TBranch *flavor3BR = outTree->Branch("flavor3",&flavor3_,"flavor3/F");
    TBranch *flavor4BR = outTree->Branch("flavor4",&flavor4_,"flavor4/F");
    TBranch *flavor5BR = outTree->Branch("flavor5",&flavor5_,"flavor5/F");
    TBranch *flavor6BR = outTree->Branch("flavor6",&flavor6_,"flavor6/F");

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

    TBranch *bestHiggsMassBR = outTree->Branch("bestHiggsMass",&bestHiggsMass_,"bestHiggsMass/F");

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
    //int run, event, lumi;
    EventInfo ev;

    float met[999];
    float hJetspt[999]; float hJetseta[999]; float hJetsphi[999]; float hJetse[999]; float hJetscsv[999]; float hJetsunc[999]; float hJetsflavor[999]; float hJetsgenpt[999];float hJetsgeneta[999];float hJetsgenphi[999];
    float aJetspt[999]; float aJetseta[999]; float aJetsphi[999]; float aJetse[999]; float aJetscsv[999]; float aJetsunc[999]; float aJetsflavor[999]; float aJetsgenpt[999];float aJetsgeneta[999];float aJetsgenphi[999];
    float vLeptonpt[999]; float vLeptoneta[999]; float vLeptonphi[999];
    float vLeptonpfCombRelIso[999];

    genTopInfo genTop, genTbar;
    genParticleInfo genB, genBbar;

    outTree->SetBranchAddress("nhJets",   &nhJets);
    outTree->SetBranchAddress("naJets",   &naJets);
    outTree->SetBranchAddress("Vtype",    &Vtype);

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

    // outTree->SetBranchAddress("EVENT.run",  &run);
    //outTree->SetBranchAddress("EVENT.event",&event);
    //outTree->SetBranchAddress("EVENT.lumi", &lumi);

    outTree->SetBranchAddress("EVENT",  &ev);
    //////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////////////////////////

    int counter  = 0;
    int counterJ = 0;

    Long64_t nentries = outTree->GetEntries();
 
    for (Long64_t i = 0; i < nentries; i++){

      if(i%10000==0) cout << i << endl;

      index1_  = -99; index2_  = -99; index3_  = -99; index4_  = -99; index5_  = -99; index6_  = -99;
      pt1_  = -99; pt2_  = -99; pt3_  = -99; pt4_  = -99; pt5_  = -99; pt6_  = -99;
      eta1_ = -99; eta2_ = -99; eta3_ = -99; eta4_ = -99; eta5_ = -99; eta6_ = -99;
      phi1_ = -99; phi2_ = -99; phi3_ = -99; phi4_ = -99; phi5_ = -99; phi6_ = -99;
      csv1_ = -99; csv2_ = -99; csv3_ = -99; csv4_ = -99; csv5_ = -99; csv6_ = -99;

      topB1_  = -99; topB2_  = -99; topB3_  = -99; topB4_  = -99; topB5_  = -99; topB6_  = -99;
      topW1_  = -99; topW2_  = -99; topW3_  = -99; topW4_  = -99; topW5_  = -99; topW6_  = -99;
      higgsB1_= -99; higgsB2_= -99; higgsB3_= -99; higgsB4_= -99; higgsB5_= -99; higgsB6_= -99;
      flavor1_= -99; flavor2_= -99; flavor3_= -99; flavor4_= -99; flavor5_= -99; flavor6_= -99;

      nLF_    = 0; nC_    = 0; nB_    = 0;
      nLFTop_ = 0; nCTop_ = 0; nBTop_ = 0;

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

      outTree->GetEntry(i);

      std::vector<math::RhoEtaPhiVector> allJetsForShapeVar;
      TObjArray jetArrayForShapeVar;
      TObjArray leptArrayForShapeVar;
      jetArrayForShapeVar.SetOwner();
      leptArrayForShapeVar.SetOwner();

      std::vector<LV> allLeptons;
      std::vector<LV> allMets;
      std::vector<LV> allBtagJetsForShapeVar;
      std::vector<LV> allUntagJetsForShapeVar;
      std::vector<float> csvBtag;
      vector<float> jetUncUntag;
      vector<float> jetUncBtag;

      std::map<float, int, sorterByPt> hMapPt;
      std::map<float, int, sorterByPt> aMapPt;
      std::map<float, int, sorterByPt> allMapPt30;
      std::map<float, int, sorterByPt> bTagMap30;

      LV topBLV(  0.,0.,0.,0.);
      LV topW1LV( 0.,0.,0.,0.); 
      LV topW2LV( 0.,0.,0.,0.);
      LV atopBLV( 0.,0.,0.,0.); 
      LV atopW1LV(0.,0.,0.,0.); 
      LV atopW2LV(0.,0.,0.,0.);

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
      LV genBLV(0.,0.,0.,0.);
      LV genBbarLV(0.,0.,0.,0.);
      if(genB.mass>0){
	genBLV.SetPt(  genB.pt );
	genBLV.SetEta( genB.eta );
	genBLV.SetPhi( genB.phi );
	genBLV.SetM(   genB.mass );
      }
      if(genBbar.mass>0){
	genBbarLV.SetPt(  genBbar.pt );
	genBbarLV.SetEta( genBbar.eta );
	genBbarLV.SetPhi( genBbar.phi );
	genBbarLV.SetM(   genBbar.mass );
      }


      float sumBtag  = 0.;
      LV MHT20LV(0.,0.,0.,0.);
      LV MHT30LV(0.,0.,0.,0.);
      LV MHT40LV(0.,0.,0.,0.);

      for(int i = 0; i < nhJets; i++){

	hMapPt[ hJetspt[i]]   = i;

	float jetMass2 = hJetse[i]*hJetse[i] -  TMath::Power(hJetspt[i]*TMath::CosH(hJetseta[i]) ,2);
	LV jetLV(hJetspt[i], hJetseta[i], hJetsphi[i], TMath::Sqrt(jetMass2));

	if( hJetspt[i] > 20.){
	  numJets20_++;
	  if( hJetscsv[i] > BTAGTHR) numJets20bTag_++;
	  MHT20LV += jetLV;
	  sumAllJet20Pt_ +=  hJetspt[i];
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

	float jetMass2 = aJetse[i]*aJetse[i] -  TMath::Power(aJetspt[i]*TMath::CosH(aJetseta[i]),2);
	LV jetLV(aJetspt[i], aJetseta[i], aJetsphi[i],  TMath::Sqrt(jetMass2));

	if( aJetspt[i] > 20.){
	  numJets20_++;
	  if( aJetscsv[i] > BTAGTHR) numJets20bTag_++;
	  MHT20LV += jetLV;
	  sumAllJet20Pt_ +=  hJetspt[i];
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
	  sumAllJet30Pt_ +=  hJetspt[i];

	  allJetsForShapeVar.push_back( math::RhoEtaPhiVector(aJetspt[i], aJetseta[i], aJetsphi[i]) );
	  TVector3 *tv3 = new TVector3(jetLV.Px(),jetLV.Py(),jetLV.Pz());
	  jetArrayForShapeVar.Add( tv3 );
	  
	}
	if( hJetspt[i]>40.){
	  numJets40_++;
	  if( aJetscsv[i] > BTAGTHR) numJets40bTag_++;
	  MHT40LV += jetLV;
	  sumAllJet40Pt_ +=  hJetspt[i];
	}
      }





      MHT20_ = MHT20LV.Pt()>0 ? MHT20LV.Pt() : -99;
      MHT30_ = MHT30LV.Pt()>0 ? MHT30LV.Pt() : -99;
      MHT40_ = MHT40LV.Pt()>0 ? MHT40LV.Pt() : -99;
      aveCsv_ = numJets30bTag_>0 ? sumBtag/numJets30bTag_ : -99;

      float sumDeltaRbTag = 0.;
      float sumMassBtag   = 0.;
      float sumMassUntag  = 0.;

      if(allBtagJetsForShapeVar.size()>1){
	float minDeltaR = 999.;
	for(unsigned k = 0; k< allBtagJetsForShapeVar.size()-1 ; k++){
	  for(unsigned l = k+1; l<allBtagJetsForShapeVar.size() ; l++){
	    LV first  = allBtagJetsForShapeVar[k];
	    LV second = allBtagJetsForShapeVar[l];
	    float deltaR = TMath::Sqrt( TMath::Power(first.Eta()-second.Eta(),2) + TMath::Power(first.Phi()-second.Phi(),2)  );
	    sumDeltaRbTag += deltaR;
	    if(deltaR < minDeltaR){
	      minDeltaRBtag_ = deltaR;
	      closestJJbTagMass_ = (first+second).M();
	      minDeltaR = deltaR;
	    }
	    sumMassBtag   += (first+second).M();
	  }
	}
      }

      aveDeltaRbTag_ = numJets30bTag_>1 ? sumDeltaRbTag/(numJets30bTag_/2*(numJets30bTag_-1)) : -99;
      aveMbTag_      = numJets30bTag_>1 ? sumMassBtag/(numJets30bTag_/2*(numJets30bTag_-1))   : -99;


      float maxCsv = -999.; float minCsv = 999.;
      for(unsigned k = 0; k< csvBtag.size() ; k++){
	varCsv_ +=  TMath::Power(csvBtag[k]-aveCsv_,2);
	if( csvBtag[k]>maxCsv ) maxCsv = csvBtag[k];
	if( csvBtag[k]<minCsv ) minCsv = csvBtag[k];
      }
      maxCsv_ = maxCsv;
      minCsv_ = minCsv;
      varCsv_ = numJets30bTag_>1 ? varCsv_/(numJets30bTag_-1) : -99;
  
      if(allUntagJetsForShapeVar.size()>1){
	for(unsigned k = 0; k< allUntagJetsForShapeVar.size()-1 ; k++){
	  for(unsigned l = k+1; l<allUntagJetsForShapeVar.size() ; l++){
	    LV first  = allUntagJetsForShapeVar[k];
	    LV second = allUntagJetsForShapeVar[l];
	    sumMassUntag   += (first+second).M();
	  }
	}
      }
      aveMunTag_     = (numJets30_-numJets30bTag_)>1 ? sumMassUntag/((numJets30_-numJets30bTag_)/2*((numJets30_-numJets30bTag_)-1)) : -99;
      

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

      int all = 0;
      for(std::map<float, int>::iterator it = allMapPt30.begin(); it!=allMapPt30.end(); it++, all++){

	int index = it->second;

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

	if(all==0){
	  index1_ = index;
	  pt1_    = (index >= 0 ) ?  hJetspt[index]  : aJetspt[ -index-1];
	  eta1_   = (index >= 0 ) ?  hJetseta[index] : aJetseta[-index-1];
	  phi1_   = (index >= 0 ) ?  hJetsphi[index] : aJetsphi[-index-1];
	  csv1_   = (index >= 0 ) ?  hJetscsv[index] : aJetscsv[-index-1];
	  flavor1_= (index >= 0 ) ?  hJetsflavor[index] : aJetsflavor[-index-1];
	  topB1_  = int(genMatch==1);
	  topW1_  = int(genMatch==2);
	  higgsB1_= int(genMatch==0);
	  if(TMath::Abs(flavor1_) == 21 || (TMath::Abs(flavor1_) > 0 && TMath::Abs(flavor1_)<4) ){
	    nLF_++;
	    if( topW1_ ) nLFTop_++;
	  }
	  if(TMath::Abs(flavor1_) == 4 ){
	    nC_++;
	    if( topW1_ ) nCTop_++;
	  }
	  if(TMath::Abs(flavor1_) == 5 ){
	    nB_++;
	    if( topB1_ ) nBTop_++;
	  }
	}
	else if(all==1){
	  index2_ = index;
	  pt2_    = (index >= 0 ) ?  hJetspt[index]  : aJetspt[ -index-1];
	  eta2_   = (index >= 0 ) ?  hJetseta[index] : aJetseta[-index-1];
	  phi2_   = (index >= 0 ) ?  hJetsphi[index] : aJetsphi[-index-1];
	  csv2_   = (index >= 0 ) ?  hJetscsv[index] : aJetscsv[-index-1];
	  flavor2_= (index >= 0 ) ?  hJetsflavor[index] : aJetsflavor[-index-1];
	  topB2_  = int(genMatch==1);
	  topW2_  = int(genMatch==2);
	  higgsB2_= int(genMatch==0);
	  if(TMath::Abs(flavor2_) == 21 || (TMath::Abs(flavor2_) > 0 && TMath::Abs(flavor2_)<4) ){
	    nLF_++;
	    if( topW2_ ) nLFTop_++;
	  }
	  if(TMath::Abs(flavor2_) == 4 ){
	    nC_++;
	    if( topW2_ ) nCTop_++;
	  }
	  if(TMath::Abs(flavor2_) == 5 ){
	    nB_++;
	    if( topB2_ ) nBTop_++;
	  }
	}
	else if(all==2){
	  index3_ = index;
	  pt3_    = (index >= 0 ) ?  hJetspt[index]  : aJetspt[ -index-1];
	  eta3_   = (index >= 0 ) ?  hJetseta[index] : aJetseta[-index-1];
	  phi3_   = (index >= 0 ) ?  hJetsphi[index] : aJetsphi[-index-1];
	  csv3_   = (index >= 0 ) ?  hJetscsv[index] : aJetscsv[-index-1];
	  flavor3_= (index >= 0 ) ?  hJetsflavor[index] : aJetsflavor[-index-1];
	  topB3_  = int(genMatch==1);
	  topW3_  = int(genMatch==2);
	  higgsB3_= int(genMatch==0);
	  if(TMath::Abs(flavor3_) == 21 || (TMath::Abs(flavor3_) > 0 && TMath::Abs(flavor3_)<4) ){
	    nLF_++;
	    if( topW3_ ) nLFTop_++;
	  }
	  if(TMath::Abs(flavor3_) == 4 ){
	    nC_++;
	    if( topW3_ ) nCTop_++;
	  }
	  if(TMath::Abs(flavor3_) == 5 ){
	    nB_++;
	    if( topB3_ ) nBTop_++;
	  }
	}
	else if(all==3){
	  index4_ = index;
	  pt4_    = (index >= 0 ) ?  hJetspt[index]  : aJetspt[ -index-1];
	  eta4_   = (index >= 0 ) ?  hJetseta[index] : aJetseta[-index-1];
	  phi4_   = (index >= 0 ) ?  hJetsphi[index] : aJetsphi[-index-1];
	  csv4_   = (index >= 0 ) ?  hJetscsv[index] : aJetscsv[-index-1];
	  flavor4_= (index >= 0 ) ?  hJetsflavor[index] : aJetsflavor[-index-1];
	  topB4_  = int(genMatch==1);
	  topW4_  = int(genMatch==2);
	  higgsB4_= int(genMatch==0);
	  if(TMath::Abs(flavor4_) == 21 || (TMath::Abs(flavor4_) > 0 && TMath::Abs(flavor4_)<4) ){
	    nLF_++;
	    if( topW4_ ) nLFTop_++;
	  }
	  if(TMath::Abs(flavor4_) == 4 ){
	    nC_++;
	    if( topW4_ ) nCTop_++;
	  }
	  if(TMath::Abs(flavor4_) == 5 ){
	    nB_++;
	    if( topB4_ ) nBTop_++;
	  }
	}
	else if(all==4){
	  index5_ = index;
	  pt5_    = (index >= 0 ) ?  hJetspt[index]  : aJetspt[ -index-1];
	  eta5_   = (index >= 0 ) ?  hJetseta[index] : aJetseta[-index-1];
	  phi5_   = (index >= 0 ) ?  hJetsphi[index] : aJetsphi[-index-1];
	  csv5_   = (index >= 0 ) ?  hJetscsv[index] : aJetscsv[-index-1];
	  flavor5_= (index >= 0 ) ?  hJetsflavor[index] : aJetsflavor[-index-1];
	  topB5_  = int(genMatch==1);
	  topW5_  = int(genMatch==2);
	  higgsB5_= int(genMatch==0); 
	  if(TMath::Abs(flavor5_) == 21 || (TMath::Abs(flavor5_) > 0 && TMath::Abs(flavor5_)<4) ){
	    nLF_++;
	    if( topW5_ ) nLFTop_++;
	  }
	  if(TMath::Abs(flavor5_) == 4 ){
	    nC_++;
	    if( topW5_ ) nCTop_++;
	  }
	  if(TMath::Abs(flavor5_) == 5 ){
	    nB_++;
	    if( topB5_ ) nBTop_++;
	  }
	}
	else if(all==5){
	  index6_ = index;
	  pt6_    = (index >= 0 ) ?  hJetspt[index]  : aJetspt[ -index-1];
	  eta6_   = (index >= 0 ) ?  hJetseta[index] : aJetseta[-index-1];
	  phi6_   = (index >= 0 ) ?  hJetsphi[index] : aJetsphi[-index-1];
	  csv6_   = (index >= 0 ) ?  hJetscsv[index] : aJetscsv[-index-1];
	  flavor6_= (index >= 0 ) ?  hJetsflavor[index] : aJetsflavor[-index-1];
	  topB6_  = int(genMatch==1);
	  topW6_  = int(genMatch==2);
	  higgsB6_= int(genMatch==0);
	  if(TMath::Abs(flavor6_) == 21 || (TMath::Abs(flavor6_) > 0 && TMath::Abs(flavor6_)<4) ){
	    nLF_++;
	    if( topW6_ ) nLFTop_++;
	  }
	  if(TMath::Abs(flavor6_) == 4 ){
	    nC_++;
	    if( topW6_ ) nCTop_++;
	  }
	  if(TMath::Abs(flavor6_) == 5 ){
	    nB_++;
	    if( topB6_ ) nBTop_++;
	  }
	}
	else{}


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


      counter++;
      bool isInJson = jsonContainsEvent(jsonVector, ev.run, ev.lumi);
      if(!isInJson){
	counterJ++;
	//cout << "Run  " << ev.run << ", " << ev.lumi << " is not in json" << endl;
      }

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
      index1BR->Fill();
      index2BR->Fill();
      index3BR->Fill();
      index4BR->Fill();
      index5BR->Fill();
      index6BR->Fill();
      pt1BR->Fill();
      pt2BR->Fill();
      pt3BR->Fill();
      pt4BR->Fill();
      pt5BR->Fill();
      pt6BR->Fill();
      eta1BR->Fill();
      eta2BR->Fill();
      eta3BR->Fill();
      eta4BR->Fill();
      eta5BR->Fill();
      eta6BR->Fill();
      phi1BR->Fill();
      phi2BR->Fill();
      phi3BR->Fill();
      phi4BR->Fill();
      phi5BR->Fill();
      phi6BR->Fill();
      csv1BR->Fill();
      csv2BR->Fill();
      csv3BR->Fill();
      csv4BR->Fill();
      csv5BR->Fill();
      csv6BR->Fill();
      topB1BR->Fill();
      topB2BR->Fill();
      topB3BR->Fill();
      topB4BR->Fill();
      topB5BR->Fill();
      topB6BR->Fill();
      topW1BR->Fill();
      topW2BR->Fill();
      topW3BR->Fill();
      topW4BR->Fill();
      topW5BR->Fill();
      topW6BR->Fill();
      higgsB1BR->Fill();
      higgsB2BR->Fill();
      higgsB3BR->Fill();
      higgsB4BR->Fill();
      higgsB5BR->Fill();
      higgsB6BR->Fill();
      flavor1BR->Fill();
      flavor2BR->Fill();
      flavor3BR->Fill();
      flavor4BR->Fill();
      flavor5BR->Fill();
      flavor6BR->Fill();

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

      delete eventShapes;
      delete topologicalWorker;

    }

    cout << counter << " => " << counterJ << endl;

    outTree->Write("",TObject::kOverwrite );

    mySamples->GetFile( currentName )->Close();

   }


  return 0;


}
