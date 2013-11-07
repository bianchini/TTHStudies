#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include <cstdlib>
#include <iostream> 
#include <fstream>
#include <map>
#include <string>

#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraph.h"
#include "Math/GenVector/LorentzVector.h"
#include "TLorentzVector.h"
#include "TVectorD.h"

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

#include "TopQuarkAnalysis/TopKinFitter/interface/TtSemiLepKinFitter.h"
#include "TopQuarkAnalysis/TopKinFitter/plugins/TtSemiLepKinFitProducer.h"

#include "Bianchi/TTHStudies/interface/Samples.h"



#define GENJETDR 0.3
#define VERBOSE  false

using namespace std;

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LV;
typedef struct 
{
  float et; 
  float sumet;
  float sig;
  float phi;
} metInfo;


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
  float mass; 
  float pt;
  float eta;
  float phi;
  float status;
  float charge;
  float momid;
} genParticleInfo;


struct sorter {
  bool operator() (Float_t i, Float_t j) const { return (i>j);}
};


Float_t deltaR(LV j1, LV j2){
  return TMath::Sqrt(TMath::Power(j1.eta()-j2.eta(),2) + TMath::Power(j1.phi()-j2.phi(),2));
}

Int_t match(LV j1, LV j2){
  return deltaR(j1,j2) < GENJETDR;
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


int main(int argc, const char* argv[])
{

  std::cout << "TopKinFitter_SL" << std::endl;
  gROOT->SetBatch(true);
 
  gSystem->Load("libFWCoreFWLite");
  gSystem->Load("libDataFormatsFWLite");

  AutoLibraryLoader::enable();

  fwlite::TFileService fs = fwlite::TFileService("topKinFitterSL.root");
  TTree* outTreeR    = fs.make<TTree>("outTreeR","tree right");
  TTree* outTreeW    = fs.make<TTree>("outTreeW","tree wrong");
  TTree* outTree     = fs.make<TTree>("outTree","event tree");

  TTree* genTree     = fs.make<TTree>("genTree","event tree");
  TTree* genJetLightTree  = fs.make<TTree>("genJetLightTree","event tree");
  TTree* genJetHeavyTree  = fs.make<TTree>("genJetHeavyTree","event tree");

  float chi2R, probR, helHadR, helLepR, bTagR;
  float xaR;
  float xbR;
  float totalPtR  ;
  float planeCosThetaR  ;
  float higgsOrientR    ;
  float topEnFracR ;
  float atopEnFracR ;
  float higgsPtStarR   ;
  float topMassR  ;
  float atopMassR ;
  float higgsMassR ;
  float topsCosDPhiR;
  float topLepPtR;
  float topLepEtaR;
  float topLepPhiR;
  float topHadPtR;
  float topHadEtaR;
  float topHadPhiR;

  float chi2W, probW, helHadW, helLepW, bTagW;
  float xaW;
  float xbW;
  float totalPtW  ;
  float planeCosThetaW  ;
  float higgsOrientW    ;
  float topEnFracW ;
  float atopEnFracW ;
  float higgsPtStarW   ;
  float topMassW  ;
  float atopMassW ;
  float higgsMassW ;
  float topsCosDPhiW;

  int nResJets, matches, nPermut, nPermutAll, status ;
  int TTnoW, fullTT, fullH;
  int WTag;
  JetByPt hjet1, hjet2;
  JetByPt w1jet, w2jet;
  genTopInfo missingTop;

  float X1;
  float X2;
  float PtTTH  ;
  float PzTTH  ;
  float PhiTTH  ;
  float BetaTTH  ;
  float GammaTTH    ;
  float DeltaTTH    ;
  float PtHStar ;
  float DphiHStar ;
  float GammaTT;
  float MassTT;
  float BetaW  ;
  float GammaW  ;
  float DeltaW  ;
  float BetaWLep  ;
  float GammaWLep  ;
  float DeltaWLep  ;
  float PxT, PyT, PzT;
  float PxH, PyH, PzH;
  //float M12;

  float ptRecoLight;
  float phiRecoLight;
  float etaRecoLight;
  float csvRecoLight;
  float ptLight;
  float etaLight;
  float phiLight;
  float massRecoLight;
  float ptRecoHeavy;
  float phiRecoHeavy;
  float etaRecoHeavy;
  float csvRecoHeavy;
  float ptHeavy;
  float etaHeavy;
  float phiHeavy;
  float massRecoHeavy;

  outTreeR->Branch("chi2",    &chi2R,    "chi2/F");
  outTreeR->Branch("prob",    &probR,    "prob/F");
  outTreeR->Branch("helHad",  &helHadR,  "helHad/F");
  outTreeR->Branch("helLep",  &helLepR,  "helLep/F");
  outTreeR->Branch("bTag",    &bTagR,    "bTag/F");
  outTreeR->Branch("xa",           &xaR,    "xa/F");
  outTreeR->Branch("xb",           &xbR,    "xb/F");
  outTreeR->Branch("totalPt",      &totalPtR,    "totalPt/F");
  outTreeR->Branch("planeCosTheta",&planeCosThetaR,    "planeCosTheta/F");
  outTreeR->Branch("higgsOrient",  &higgsOrientR,    "higgsOrient/F");
  outTreeR->Branch("topEnFrac",    &topEnFracR,    "topEnFrac/F");
  outTreeR->Branch("atopEnFrac",   &atopEnFracR,    "atopEnFrac/F");
  outTreeR->Branch("higgsPtStar",    &higgsPtStarR,    "higgsPtStar/F");
  outTreeR->Branch("topMass",    &topMassR,    "topMass/F");
  outTreeR->Branch("atopMass",   &atopMassR,    "atopMass/F");
  outTreeR->Branch("higgsMass",  &higgsMassR,   "higgsMass/F");
  outTreeR->Branch("topsCosDPhi",  &topsCosDPhiR,    "topsCosDPhi/F");
  outTreeR->Branch("topLepPt",     &topLepPtR,    "topLepPt/F");
  outTreeR->Branch("topLepEta",    &topLepEtaR,    "topLepEta/F");
  outTreeR->Branch("topLepPhi",    &topLepPhiR,    "topLepPhi/F");
  outTreeR->Branch("topHadPt",     &topHadPtR,    "topHadPt/F");
  outTreeR->Branch("topHadEta",    &topHadEtaR,    "topHadEta/F");
  outTreeR->Branch("topHadPhi",    &topHadPhiR,    "topHadPhi/F");



  outTreeW->Branch("chi2",    &chi2W,    "chi2/F");
  outTreeW->Branch("prob",    &probW,    "prob/F");
  outTreeW->Branch("helHad",  &helHadW,  "helHad/F");
  outTreeW->Branch("helLep",  &helLepW,  "helLep/F");
  outTreeW->Branch("bTag",    &bTagW,    "bTag/F");
  outTreeW->Branch("xa",    &xaW,    "xa/F");
  outTreeW->Branch("xb",    &xbW,    "xb/F");
  outTreeW->Branch("totalPt",    &totalPtW,    "totalPt/F");
  outTreeW->Branch("planeCosTheta",    &planeCosThetaW,    "planeCosTheta/F");
  outTreeW->Branch("higgsOrient",    &higgsOrientW,    "higgsOrient/F");
  outTreeW->Branch("topEnFrac",    &topEnFracW,    "topEnFrac/F");
  outTreeW->Branch("atopEnFrac",    &atopEnFracW,    "atopEnFrac/F");
  outTreeW->Branch("higgsPtStar",    &higgsPtStarR,    "higgsPtStar/F");
  outTreeW->Branch("topMass",    &topMassW,    "topMass/F");
  outTreeW->Branch("atopMass",   &atopMassW,    "atopMass/F");
  outTreeW->Branch("higgsMass",  &higgsMassW,   "higgsMass/F");
  outTreeW->Branch("topsCosDPhi",    &topsCosDPhiW,    "topsCosDPhi/F");
  //outTreeR->Branch("",    &R,    "/F");

  outTree->Branch("hjet1",    &hjet1,       "index/I:pt/F:eta/F:phi/F:mass/F:csv/F:topB/F:topW/F:atopB/F:atopW/F:higgsB/F:flavor/F:unc/F");
  outTree->Branch("hjet2",    &hjet2,       "index/I:pt/F:eta/F:phi/F:mass/F:csv/F:topB/F:topW/F:atopB/F:atopW/F:higgsB/F:flavor/F:unc/F");
  outTree->Branch("w1jet",    &w1jet,       "index/I:pt/F:eta/F:phi/F:mass/F:csv/F:topB/F:topW/F:atopB/F:atopW/F:higgsB/F:flavor/F:unc/F");
  outTree->Branch("w2jet",    &w2jet,       "index/I:pt/F:eta/F:phi/F:mass/F:csv/F:topB/F:topW/F:atopB/F:atopW/F:higgsB/F:flavor/F:unc/F");
  outTree->Branch("nResJets",&nResJets,     "nResJets/I");
  outTree->Branch("matches", &matches,      "matches/I");
  outTree->Branch("nPermut", &nPermut,      "nPermut/I");
  outTree->Branch("nPermutAll", &nPermutAll,"nPermutAll/I");
  outTree->Branch("status",     &status,    "status/I");
  outTree->Branch("fullTT",     &fullTT,    "fullTT/I");
  outTree->Branch("fullH",      &fullH,     "fullH/I");
  outTree->Branch("TTnoW",      &TTnoW,     "TTnoW/I");
  outTree->Branch("WTag",       &WTag,      "WTag/I");
  outTree->Branch("genTop",     &missingTop,  "bmass/F:bpt/F:beta:bphi/F:bstatus/F:wdau1mass/F:wdau1pt/F:wdau1eta:wdau1phi/F:wdau1id/F:wdau2mass/F:wdau2pt/F:wdau2eta:wdau2phi/F:wdau2id/F");


  genTree->Branch("X1",       &X1,     "X1/F");
  genTree->Branch("X2",       &X2,     "X2/F");
  genTree->Branch("PtTTH",    &PtTTH,  "PtTTH/F");
  genTree->Branch("PzTTH",    &PzTTH,  "PzTTH/F");
  genTree->Branch("PhiTTH",   &PhiTTH, "PhiTTH/F");
  genTree->Branch("BetaTTH",  &BetaTTH,     "BetaTTH/F");
  genTree->Branch("GammaTTH", &GammaTTH,     "GammaTTH/F");
  genTree->Branch("GammaTT",  &GammaTT,      "GammaTT/F");
  genTree->Branch("DeltaTTH", &DeltaTTH,     "DeltaTTH/F");
  genTree->Branch("PtHStar",  &PtHStar,       "PtHStar/F");
  genTree->Branch("DphiHStar",&DphiHStar,    "DphiHStar/F");
  genTree->Branch("MassTT",   &MassTT,       "MassTT/F");
  genTree->Branch("BetaW",    &BetaW,        "BetaW/F");
  genTree->Branch("GammaW",   &GammaW,       "GammaW/F");
  genTree->Branch("DeltaW",   &DeltaW,       "DeltaW/F");
  genTree->Branch("BetaWLep", &BetaWLep,     "BetaWLep/F");
  genTree->Branch("GammaWLep",&GammaWLep,    "GammaWLep/F");
  genTree->Branch("DeltaWLep",&DeltaWLep,    "DeltaWLep/F");


  genTree->Branch("PxT",   &PxT,       "PxT/F");
  genTree->Branch("PyT",   &PyT,       "PyT/F");
  genTree->Branch("PzT",   &PzT,       "PzT/F");
  genTree->Branch("PxH",   &PxH,       "PxH/F");
  genTree->Branch("PyH",   &PyH,       "PyH/F");
  genTree->Branch("PzH",   &PzH,       "PzH/F");

  genJetLightTree->Branch("ptReco",   &ptRecoLight, "ptReco/F");
  genJetLightTree->Branch("etaReco",  &etaRecoLight, "etaReco/F");
  genJetLightTree->Branch("phiReco",  &phiRecoLight, "phiReco/F");
  genJetLightTree->Branch("csvReco",   &csvRecoLight, "csvReco/F");
  genJetLightTree->Branch("pt",       &ptLight,     "pt/F");
  genJetLightTree->Branch("eta",      &etaLight,    "eta/F");
  genJetLightTree->Branch("phi",      &phiLight,    "phi/F");
  genJetLightTree->Branch("massReco",     &massRecoLight,   "massReco/F");

  genJetHeavyTree->Branch("ptReco",   &ptRecoHeavy, "ptReco/F");
  genJetHeavyTree->Branch("etaReco",  &etaRecoHeavy, "etaReco/F");
  genJetHeavyTree->Branch("phiReco",  &phiRecoHeavy, "phiReco/F");
  genJetHeavyTree->Branch("csvReco",   &csvRecoHeavy, "csvReco/F");
  genJetHeavyTree->Branch("pt",       &ptHeavy,     "pt/F");
  genJetHeavyTree->Branch("eta",      &etaHeavy,    "eta/F");
  genJetHeavyTree->Branch("phi",      &phiHeavy,    "phi/F");
  genJetHeavyTree->Branch("massReco",     &massRecoHeavy,   "massReco/F");


  PythonProcessDesc builder(argv[1]);
 
  const edm::ParameterSet& in = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteInput");

  const edm::VParameterSet& samples = in.getParameter<edm::VParameterSet>("samples") ;

  std::string pathToFile( in.getParameter<std::string>("pathToFile" ) );
  std::string ordering(   in.getParameter<std::string>("ordering" ) );
  std::vector< edm::ParameterSet > udscResolutions(in.getParameter<edm::VParameterSet>("udscResolutions") );
  std::vector< edm::ParameterSet > bResolutions   (in.getParameter<edm::VParameterSet>("bResolutions") );
  std::vector< edm::ParameterSet > lepResolutions (in.getParameter<edm::VParameterSet>("lepResolutions") );
  std::vector<double> metResolutionsCoeff         (in.getParameter<std::vector<double> >("metResolutionsCoeff") );
  std::vector<double> helicityCoeff               (in.getParameter<std::vector<double> >("helicityCoeff") );
  std::vector<double> chi2Coeff                   (in.getParameter<std::vector<double> >("chi2Coeff") );
  std::vector<double> sgnMassCoeff                (in.getParameter<std::vector<double> >("sgnMassCoeff") );

  if(helicityCoeff.size()!=9){
    cout << "Nine parameters are required for the helicity angle param... return" << endl;
    return 0;
  }

  double lumi              (in.getParameter<double>("lumi") );
  bool verbose             (in.getParameter<bool>("verbose") );
  bool computeCMSVariables (in.getParameter<bool>("computeCMSVariables") );
  double udscPtCut(in.getParameter<double>("udscPtCut") );
  double bPtCut   (in.getParameter<double>("bPtCut") );
  std::string likelihoodFile(in.getParameter<std::string>("likelihoodFile") );
  vector<std::string> likelihoods(in.getParameter<vector<std::string> >("likelihoods") );

  int useKin   = 0;
  int useBtag  = 0;
  int useHel   = 0;
  int useHMass = 0;

  for(unsigned int l = 0; l < likelihoods.size() ; l++){
    if( likelihoods[l].find("KIN")   !=string::npos ) useKin   =1;
    if( likelihoods[l].find("BTAG")  !=string::npos ) useBtag  =1;
    if( likelihoods[l].find("HEL")   !=string::npos ) useHel   =1;
    if( likelihoods[l].find("HMASS") !=string::npos ) useHMass =1;
  }

  TFile* fLikelihoods          = TFile::Open(likelihoodFile.c_str(),"READ");
  TF2* fLikelihoodWTopHadrMass = (TF2*)fLikelihoods->Get("fLikelihoodWTopHadrMass");
  TF1* fLikelihoodCsvB         = (TF1*)fLikelihoods->Get("fLikelihoodCsvB");
  TF1* fLikelihoodCsvNonB      = (TF1*)fLikelihoods->Get("fLikelihoodCsvNonB");
  TF1* fLikelihoodTopLeptMass  = (TF1*)fLikelihoods->Get("fLikelihoodTopLeptMass");

  TF1* fHelTopLept = new TF1("fHelTopLept",Form("%f + %f*x + %f*x*x + %f*x*x*x",              helicityCoeff[0], helicityCoeff[1], helicityCoeff[2], helicityCoeff[3]),-1,1);
  TF1* fHelTopHadr = new TF1("fHelTopHadr",Form("%f + %f*x + %f*x*x + %f*x*x*x + %f*x*x*x*x", helicityCoeff[4], helicityCoeff[5], helicityCoeff[6], helicityCoeff[7], helicityCoeff[8]),-1,1);

  float ndf = chi2Coeff.size()>0 ? chi2Coeff[0] : 3.;
  TF1* fChi2       = new TF1("fChi2", Form("1./TMath::Gamma(%f/2)/2.^(%f/2)*TMath::Exp(-x/2)*(x^(%f/2-1))", chi2Coeff[0],chi2Coeff[1],chi2Coeff[2]), 0, 100);

  TF1 *fSgnMass = new TF1("fSgnMass","gaus",0,1000);
  fSgnMass->SetParameters(1./TMath::Sqrt(2*TMath::Pi())/sgnMassCoeff[1],sgnMassCoeff[0],sgnMassCoeff[1]);

  bool openAllFiles  = false;
  Samples* mySamples = new Samples(openAllFiles, pathToFile, ordering, samples, lumi, verbose);
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

    mySamples->OpenFile( currentName );
    cout << "Opening file " << currentName << endl;
    TTree* currentTree       = mySamples->GetTree( currentName, "tree");
    cout << "Done!!" << endl;
    
    JetByPt jet1, jet2, jet3, jet4, jet5, jet6, jet7, jet8, jet9, jet10;
    genParticleInfo genB, genBbar;
    genTopInfo genTop, genTbar;

    Float_t vLepton_mass[99];
    Float_t vLepton_pt  [99];
    Float_t vLepton_eta [99];
    Float_t vLepton_phi [99];
    Float_t vLepton_charge [99];
    metInfo METtype1p2corr;

    currentTree->SetBranchAddress("jet1",   &jet1);
    currentTree->SetBranchAddress("jet2",   &jet2);
    currentTree->SetBranchAddress("jet3",   &jet3);
    currentTree->SetBranchAddress("jet4",   &jet4);
    currentTree->SetBranchAddress("jet5",   &jet5);
    currentTree->SetBranchAddress("jet6",   &jet6);
    currentTree->SetBranchAddress("jet7",   &jet7);
    currentTree->SetBranchAddress("jet8",   &jet8);
    currentTree->SetBranchAddress("jet9",   &jet9);
    currentTree->SetBranchAddress("jet10",  &jet10);
    currentTree->SetBranchAddress("genB",   &genB);
    currentTree->SetBranchAddress("genBbar",&genBbar);
    currentTree->SetBranchAddress("genTop", &genTop);
    currentTree->SetBranchAddress("genTbar",&genTbar);
    currentTree->SetBranchAddress("vLepton_mass"  ,vLepton_mass);
    currentTree->SetBranchAddress("vLepton_pt"    ,vLepton_pt);
    currentTree->SetBranchAddress("vLepton_eta"   ,vLepton_eta);
    currentTree->SetBranchAddress("vLepton_phi"   ,vLepton_phi);
    currentTree->SetBranchAddress("vLepton_charge",vLepton_charge);
    
    currentTree->SetBranchAddress("METtype1p2corr",&METtype1p2corr);


    Long64_t nentries = currentTree->GetEntries();
    
    for (Long64_t i = 0; i < nentries ; i++){
      
      if(i%1000==0) cout << i << endl;
      currentTree->GetEntry(i);

      //continue;
      
      LV topBLV   (0.,0.,0.,0.);
      LV topW1LV  (0.,0.,0.,0.); 
      LV topW2LV  (0.,0.,0.,0.);
      LV atopBLV  (0.,0.,0.,0.); 
      LV atopW1LV (0.,0.,0.,0.); 
      LV atopW2LV (0.,0.,0.,0.);
      LV genBLV   (0.,0.,0.,0.);
      LV genBbarLV(0.,0.,0.,0.);
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
      vector<LV> myJetsFilt;
      std::map<unsigned int, JetByPt> mapFilt;
      
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
      
      int numJets20 =0;     
      int numJets30 =0; 
      int numJets60 =0; 
      int numJets20BtagL=0; 
      int numJets20BtagM=0; 
      int numJets20BtagT=0;
      int numJets30BtagL=0; 
      int numJets30BtagM=0; 
      int numJets30BtagT=0;
      int numJets60BtagL=0; 
      int numJets60BtagM=0; 
      int numJets60BtagT=0;

      for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
	float csv_k = mapFilt[k].csv>0 ? mapFilt[k].csv : 0.0 ;
	float pt_k  = mapFilt[k].pt>0  ? mapFilt[k].pt  : 0.0 ;
	
	if( pt_k>20 && csv_k>0.244 ) numJets20BtagL++;
	if( pt_k>20 && csv_k>0.679 ) numJets20BtagM++;
	if( pt_k>20 && csv_k>0.898 ) numJets20BtagT++;
	if( pt_k>30 && csv_k>0.244 ) numJets30BtagL++;
	if( pt_k>30 && csv_k>0.679 ) numJets30BtagM++;
	if( pt_k>30 && csv_k>0.898 ) numJets30BtagT++;
	if( pt_k>60 && csv_k>0.244 ) numJets60BtagL++;
	if( pt_k>60 && csv_k>0.679 ) numJets60BtagM++;
	if( pt_k>60 && csv_k>0.898 ) numJets60BtagT++;
	if( pt_k>20     ) numJets20++;
	if( pt_k>30     ) numJets30++;
	if( pt_k>60     ) numJets60++;
	
	//bool isB    = genBLV.Pt()>0 ?    deltaR( myJetsFilt[k], genBLV)    < 0.3 && TMath::Abs(myJetsFilt[k].Pt()-genBLV.Pt())/genBLV.Pt()<0.30        : false;
	//bool isBbar = genBbarLV.Pt()>0 ? deltaR( myJetsFilt[k], genBbarLV) < 0.3 && TMath::Abs(myJetsFilt[k].Pt()-genBbarLV.Pt())/genBbarLV.Pt()<0.30  : false;      
      }

      if( myJetsFilt.size()<4 || numJets30<6 || numJets30BtagM<4) continue;


      ///////////////////////////////////////////////////////////////////////////

      ptRecoLight=-99;
      etaRecoLight=-99;
      phiRecoLight=-99;
      csvRecoLight=-99;
      ptLight=-99;
      etaLight=-99;
      phiLight=-99;
      massRecoLight=-99;
      ptRecoHeavy=-99;
      etaRecoHeavy=-99;
      phiRecoHeavy=-99;
      csvRecoHeavy=-99;
      ptHeavy=-99;
      etaHeavy=-99;
      phiHeavy=-99;
      massRecoHeavy=-99;

      if(abs(genTop.wdau1id)<6 &&  abs(genTbar.wdau1id)>6){

	unsigned int indexB = 999;
	for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
	  if(deltaR(myJetsFilt[k],topBLV)<0.3)
	    indexB = k;
	}
	if(indexB!=999){
	  ptRecoHeavy   = myJetsFilt[indexB].Pt();
	  etaRecoHeavy   = myJetsFilt[indexB].Eta();
	  phiRecoHeavy   = myJetsFilt[indexB].Phi();
	  massRecoHeavy = myJetsFilt[indexB].M();
	  csvRecoHeavy  =  mapFilt[indexB].csv>0 ? mapFilt[indexB].csv : 0.0 ;
	}else{ 
	  ptRecoHeavy=-99;
	  etaRecoHeavy=-99;
	  phiRecoHeavy=-99;
	  csvRecoHeavy=-99;
	  massRecoHeavy=-99;
	}
	ptHeavy       = topBLV.Pt();
	etaHeavy      = topBLV.Eta();
	phiHeavy      = topBLV.Phi();
	genJetHeavyTree->Fill();
	

	unsigned int indexW1 = 999;
	for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
	  if(deltaR(myJetsFilt[k],topW1LV)<0.3)
	    indexW1 = k;
	}
	if(indexW1!=999){
	  ptRecoLight   = myJetsFilt[indexW1].Pt();
	  etaRecoLight   = myJetsFilt[indexW1].Eta();
	  phiRecoLight   = myJetsFilt[indexW1].Phi();
	  massRecoLight = myJetsFilt[indexW1].M();
	  csvRecoLight  =  mapFilt[indexW1].csv>0 ? mapFilt[indexW1].csv : 0.0 ;
	}
	else{
	  ptRecoLight=-99;
	  csvRecoLight=-99;
	  massRecoLight=-99;
	  etaRecoLight=-99;
	  phiRecoLight=-99;
	}
	ptLight       = topW1LV.Pt();
	etaLight      = topW1LV.Eta();
	phiLight      = topW1LV.Phi();
	genJetLightTree->Fill();
	

	unsigned int indexW2 = 999;
	for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
	  if(deltaR(myJetsFilt[k],topW2LV)<0.3)
	    indexW2 = k;
	}
	if(indexW2!=999){
	  ptRecoLight   = myJetsFilt[indexW2].Pt();
	  massRecoLight = myJetsFilt[indexW2].M();
	  etaRecoLight   = myJetsFilt[indexW2].Eta();
	  phiRecoLight   = myJetsFilt[indexW2].Phi();
	  csvRecoLight  =  mapFilt[indexW2].csv>0 ? mapFilt[indexW2].csv : 0.0 ;
	}else{
	  ptRecoLight=-99;
	  csvRecoLight=-99;
	  massRecoLight=-99;
	  etaRecoLight=-99;
	  phiRecoLight=-99;
	}	
	ptLight       = topW2LV.Pt();
	etaLight      = topW2LV.Eta();
	phiLight      = topW2LV.Phi();
	genJetLightTree->Fill();
	
	unsigned int indexBbar = 999;
	for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
	  if(deltaR(myJetsFilt[k],atopBLV)<0.3)
	    indexBbar = k;
	}
	if(indexBbar!=999){
	  ptRecoHeavy   = myJetsFilt[indexBbar].Pt();
	  massRecoHeavy = myJetsFilt[indexBbar].M();
	  csvRecoHeavy  =  mapFilt[indexBbar].csv>0 ? mapFilt[indexBbar].csv : 0.0 ;
	  etaRecoHeavy   = myJetsFilt[indexBbar].Eta();
	  phiRecoHeavy   = myJetsFilt[indexBbar].Phi();
	}else{
	  ptRecoHeavy=-99;
	  csvRecoHeavy=-99;
	  massRecoHeavy=-99;
	  etaRecoHeavy=-99;
	  phiRecoHeavy=-99;
	}
	ptHeavy       = atopBLV.Pt();
	etaHeavy      = atopBLV.Eta();
	phiHeavy      = atopBLV.Phi();
	genJetHeavyTree->Fill();
	
	unsigned int indexH1 = 999;
	for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
	  if(deltaR(myJetsFilt[k],genBLV)<0.3)
	    indexH1 = k;
	}
	if(indexH1!=999){
	  ptRecoHeavy   = myJetsFilt[indexH1].Pt();
	  massRecoHeavy = myJetsFilt[indexH1].M();
	  csvRecoHeavy  =  mapFilt[indexH1].csv>0 ? mapFilt[indexH1].csv : 0.0 ;
	  etaRecoHeavy   = myJetsFilt[indexH1].Eta();
	  phiRecoHeavy   = myJetsFilt[indexH1].Phi();
	}
	else{
	  ptRecoHeavy=-99;
	  csvRecoHeavy=-99;
	  massRecoHeavy=-99;
	  etaRecoHeavy=-99;
	  phiRecoHeavy=-99;
	}
	ptHeavy       = genBLV.Pt();
	etaHeavy      = genBLV.Eta();
	phiHeavy      = genBLV.Phi();
	genJetHeavyTree->Fill();
	
	unsigned int indexH2 = 999;
	for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
	  if(deltaR(myJetsFilt[k],genBbarLV)<0.3)
	    indexH2 = k;
	}
	if(indexH2!=999){
	  ptRecoHeavy   = myJetsFilt[indexH2].Pt();
	  massRecoHeavy = myJetsFilt[indexH2].M();
	  csvRecoHeavy  =  mapFilt[indexH2].csv>0 ? mapFilt[indexH2].csv : 0.0 ;
	  etaRecoHeavy   = myJetsFilt[indexH2].Eta();
	  phiRecoHeavy   = myJetsFilt[indexH2].Phi();
	}
	else{
	  ptRecoHeavy=-99;
	  csvRecoHeavy=-99;
	  massRecoHeavy=-99;
	  etaRecoHeavy=-99;
	  phiRecoHeavy=-99;
	}
	ptHeavy       = genBbarLV.Pt();
	etaHeavy      = genBbarLV.Eta();
	phiHeavy      = genBbarLV.Phi();
	genJetHeavyTree->Fill();
       
      }


      ////////////////////////////  TOP KIN FITTER    ////////////////////////////

      std::vector< TtSemiLepKinFitter::Constraint > constraints; 
      constraints.push_back(TtSemiLepKinFitter::kWHadMass);
      constraints.push_back(TtSemiLepKinFitter::kWLepMass);
      constraints.push_back(TtSemiLepKinFitter::kTopHadMass);
      constraints.push_back(TtSemiLepKinFitter::kTopLepMass);
      //constraints.push_back(TtSemiLepKinFitter::kEqualTopMasses);
      //constraints.push_back(TtSemiLepKinFitter::kNeutrinoMass); // causes crash

      std::vector< double > jetEnergyResolutionScaleFactors;
      std::vector< double > jetEnergyResolutionEtaBinning;
      jetEnergyResolutionScaleFactors.push_back(1.);
      jetEnergyResolutionEtaBinning.push_back(0.);
      jetEnergyResolutionEtaBinning.push_back(-1.);


      std::vector< edm::ParameterSet > metResolutions;
      edm::ParameterSet paramSetMet;  
      paramSetMet.addParameter("et",string(Form("sqrt(%f + %f*%.5f)",      metResolutionsCoeff[0],  metResolutionsCoeff[1],  METtype1p2corr.sumet)));
      paramSetMet.addParameter("eta",string("999"));
      paramSetMet.addParameter("phi",string(Form("sqrt(%f + %f*%.5f)/%.5f",metResolutionsCoeff[2],  metResolutionsCoeff[3],  METtype1p2corr.sumet,METtype1p2corr.et)));
      metResolutions .push_back(paramSetMet);
	
      int maxNbIter = 200;
      TtSemiLepKinFitter *fitter = 
	new TtSemiLepKinFitter(TtSemiLepKinFitter::kEtEtaPhi,TtSemiLepKinFitter::kEtEtaPhi,TtSemiLepKinFitter::kEtEtaPhi,
			       maxNbIter,5e-5,1e-4,constraints,80.4,173.,
			       &udscResolutions,&bResolutions,&lepResolutions,&metResolutions,
			       &jetEnergyResolutionScaleFactors,&jetEnergyResolutionEtaBinning);


      LV leptonLV(vLepton_pt[0],vLepton_eta[0],vLepton_phi[0],vLepton_mass[0]);
      float nuPx = METtype1p2corr.et*TMath::Cos(METtype1p2corr.phi);
      float nuPy = METtype1p2corr.et*TMath::Sin(METtype1p2corr.phi);
      float nuE  = TMath::Sqrt(nuPx*nuPx+nuPy*nuPy);

      TLorentzVector p4Lepton(leptonLV.Px(),leptonLV.Py(),leptonLV.Pz(),leptonLV.E());
      TLorentzVector p4Neutrino(nuPx,nuPy,0,nuE);
      TLorentzVector p4HadB(0.,0.,0.,0.);
      TLorentzVector p4HadP(0.,0.,0.,0.);
      TLorentzVector p4HadQ(0.,0.,0.,0.);
      TLorentzVector p4LepB(0.,0.,0.,0.);

      unsigned int indices[myJetsFilt.size()];
      for(unsigned int k = 0; k<myJetsFilt.size() ; k++)
	indices[k]=k; 
      vector<unsigned int> visited;

      
      ////////////////////////////////////////
      unsigned int mmatchW1 = 999;
      unsigned int mmatchW2 = 999;
      unsigned int mmatchB1 = 999;
      unsigned int mmatchB2 = 999;
      unsigned int mmatchH1 = 999;
      unsigned int mmatchH2 = 999;
      unsigned int indexLight1=999;
      unsigned int indexLight2=999;

      for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
	 float pt_k = mapFilt[k].pt>0 ? mapFilt[k].pt : 0.0 ;
	 float eta_k = mapFilt[k].eta;
	 if(pt_k<30 || TMath::Abs(eta_k)>2.5) continue;
	 float csv_k = mapFilt[k].csv>0 ? mapFilt[k].csv : 0.0 ;

	 if( deltaR(myJetsFilt[k], topBLV)<0.3 ){
	   if(csv_k>0.679) mmatchB1=k;
	 }
	 else if( (abs(genTop.wdau1id)<6 && deltaR(myJetsFilt[k], topW1LV)<0.3) || 
	     (abs(genTbar.wdau1id)<6 && deltaR(myJetsFilt[k], atopW1LV)<0.3) ){
	   if(csv_k<0.679){
	     mmatchW1    = k;
	     indexLight1 = k;
	   }
	 }
	 else if( (abs(genTop.wdau1id)<6 && deltaR(myJetsFilt[k], topW2LV)<0.3) || 
	     (abs(genTbar.wdau1id)<6 && deltaR(myJetsFilt[k], atopW2LV)<0.3) ){
	   if(csv_k<0.679){
	     mmatchW2    = k;
	     indexLight2 = k;
	   }
	 }
	 else if( deltaR(myJetsFilt[k], atopBLV)<0.3 ){
	   if(csv_k>0.679) mmatchB2=k;
	 }
	 else if( deltaR(myJetsFilt[k], genBLV)<0.3 ){
	   if(csv_k>0.679) mmatchH1=k;
	 }
	 else if( deltaR(myJetsFilt[k], genBbarLV)<0.3 ){
	   if(csv_k>0.679) mmatchH2=k;
	 }
	 else{}
       }
      if( mmatchW1!=999 && mmatchW2!=999 && mmatchB1!=999 && mmatchB2!=999 ){
	fullTT = 1;
      }
      else
	fullTT = 0;
      if( mmatchH1!=999 && mmatchH2!=999) 
	fullH = 1;
      else
	fullH = 0;
      if( (mmatchW1==999 || mmatchW1==999) && mmatchB1!=999 && mmatchB2!=999) 
	TTnoW = 1;
      else
	TTnoW = 0;


      WTag = 0;
      for(unsigned int ii = 0; ii < myJetsFilt.size()-1; ii++){
	for(unsigned int jj = ii+1; jj < myJetsFilt.size(); jj++){
	  for(unsigned int swapJ = 0; swapJ<2 ; swapJ++){
	    unsigned int indexWdau1Cand = swapJ==0 ? ii : jj;//indices[0];
	    unsigned int indexWdau2Cand = swapJ==0 ? jj : ii;//indices[1];
	    
	    if( !( (mapFilt[indexWdau1Cand ]).pt >  udscPtCut && (mapFilt[indexWdau2Cand ]).pt > udscPtCut ) )
	      continue;
	    
	    if( !( (mapFilt[indexWdau1Cand ]).csv < 0.679   && (mapFilt[indexWdau2Cand ]).csv <0.679) ) 
	      continue;
	    
	    float WMass   = (myJetsFilt[indexWdau1Cand]+myJetsFilt[indexWdau2Cand]).M();
	    if( WMass  >  40 && WMass  < 120 ){
	      WTag++;
	      if(indexLight1==999 && indexLight2==999){
		indexLight1 = ii;
		indexLight2 = jj;
	      }
	    }
	  } 
	}
      }
      
      w1jet.reset(); 
      w2jet.reset();
      if(indexLight1!=999) w1jet        = mapFilt[indexLight1];
      if(indexLight2!=999) w2jet        = mapFilt[indexLight2];
      missingTop   = abs(genTop.wdau1id)<6 ? genTop : genTbar;

      



      ////////////////////////////////////////

      int counterPer = 0;
      nPermutAll     = 0;
      status         = 0;
      float minLogP  = 999;
      vector<unsigned int> bestRank; 
      bestRank.clear();
      bestRank.push_back(0);	bestRank.push_back(0);
      bestRank.push_back(0);	bestRank.push_back(0);
      int foundMatch=0;
      
      //do {
      for(unsigned int ii = 0; ii < myJetsFilt.size()-1; ii++){
	for(unsigned int jj = ii+1; jj < myJetsFilt.size(); jj++){
	  for(unsigned int kk = 0; kk < myJetsFilt.size()-1; kk++){
	    if(kk==ii || kk==jj) continue;
	    for(unsigned int ll = kk+1; ll < myJetsFilt.size(); ll++){
	      if(ll==ii || ll==jj) continue;
	      for(unsigned int swapB = 0; swapB<2 ; swapB++){
		for(unsigned int swapJ = 0; swapJ<2 ; swapJ++){
		
		  float logp = 999;

		  unsigned int indexWdau1Cand = swapJ==0 ? ii : jj;//indices[0];
		  unsigned int indexWdau2Cand = swapJ==0 ? jj : ii;//indices[1];
		  unsigned int indexbCand     = swapB==0 ? kk : ll;//indices[2];
		  unsigned int indexbbarCand  = swapB==0 ? ll : kk;//indices[3];
		  
		  int matchW1 = 0; // W had
		  int matchW2 = 0; // W had
		  int matchB1 = 0; // b from top had
		  int matchB2 = 0; // b from top lept
		  //int matchH1 = 0; // higgs
		  //int matchH2 = 0; // higgs

		  if( (topBLV.Pt()>0 && topW1LV.Pt()>0 && topW2LV.Pt()>0 && atopBLV.Pt()>0 && atopW1LV.Pt()>0 && atopW2LV.Pt()>0 )){
		    
		    if(abs(genTop.wdau1id)<6 &&  abs(genTbar.wdau1id)>6){
		      if(deltaR(myJetsFilt[indexbCand],    topBLV)<0.3){
			matchB1++;
		      }
		      if( deltaR(myJetsFilt[indexWdau1Cand], topW1LV)<0.3 && deltaR(myJetsFilt[indexWdau2Cand], topW2LV) <0.3){
			matchW1++;
			matchW2++;
		      }
		      if(deltaR(myJetsFilt[indexbbarCand],   atopBLV)<0.3)
			matchB2++;
		    }
		    else if(abs(genTop.wdau1id)>6 &&  abs(genTbar.wdau1id)<6){
		      if(deltaR(myJetsFilt[indexbCand],    atopBLV)<0.3)
			matchB1++;
		      if( deltaR(myJetsFilt[indexWdau1Cand], atopW1LV)<0.3 && deltaR(myJetsFilt[indexWdau2Cand], atopW2LV) <0.3){
			matchW1++;
			matchW2++;
		      }
		      if(deltaR(myJetsFilt[indexbbarCand],   topBLV)<0.3)
			matchB2++;
		    }
		    else{}


		    //if( (deltaR(myJetsFilt[indexbCand],    topBLV)<0.3 && abs(genTop.wdau1id) <6) || 
		    //(deltaR(myJetsFilt[indexbCand],   atopBLV)<0.3 && abs(genTbar.wdau1id)<6)
		    //){
		    //matchB1++;
		    //}
		    //if( (abs(genTop.wdau1id)<6  && deltaR(myJetsFilt[indexWdau1Cand], topW1LV)<0.3 && deltaR(myJetsFilt[indexWdau2Cand], topW2LV) <0.3) || 
		    //(abs(genTbar.wdau1id)<6 && deltaR(myJetsFilt[indexWdau1Cand],atopW1LV)<0.3 && deltaR(myJetsFilt[indexWdau2Cand], atopW2LV)<0.3) ){
		    //matchW1++;
		    //matchW2++;
		    //}
		    //if( (deltaR(myJetsFilt[indexbbarCand],    topBLV)<0.3 && abs(genTop.wdau1id) >6) || 
		    //(deltaR(myJetsFilt[indexbbarCand],   atopBLV)<0.3 && abs(genTbar.wdau1id)>6) ){
		    //matchB2++;
		    //}
		    
		    //cout << matchW1 << ", " << matchW2 << ", " << matchB1 << ", " << matchB2 << endl;  
		  }
		  
		  
		  //if( !bookKeeper( indexWdau1Cand, indexWdau2Cand, indexbCand, indexbbarCand, visited) ) 
		  //continue;
		  nPermutAll++;
		  
		  //if(verbose || (i==1930 && indexWdau1Cand==5 && indexWdau2Cand==6 && indexbCand==2 && indexbbarCand==1)){
		  //cout << "Perm #" << counterPer << " [" << indexWdau1Cand << "," << indexWdau2Cand << "," << indexbCand << "," << indexbbarCand << " => STEP1" <<endl;
		  //cout << (mapFilt[indexWdau1Cand ]).pt << "," << (mapFilt[indexWdau2Cand ]).pt  << ", " << (mapFilt[indexbCand ]).pt   << "," << (mapFilt[indexbbarCand ]).pt  << endl;
		  //}
		  
		  if( !( (mapFilt[indexWdau1Cand ]).pt >  udscPtCut && (mapFilt[indexWdau2Cand ]).pt > udscPtCut  &&
			 (mapFilt[indexbCand ]).pt     >  bPtCut    && (mapFilt[indexbbarCand ]).pt  > bPtCut ) )
		    continue;
		  
		  //if(verbose || (i==1930 && indexWdau1Cand==5 && indexWdau2Cand==6 && indexbCand==2 && indexbbarCand==1))
		  //cout << "Perm #" << counterPer << " [" << indexWdau1Cand << "," << indexWdau2Cand << "," << indexbCand << "," << indexbbarCand << " => STEP2" <<endl;
		  
		  if( !((mapFilt[indexbCand ]).csv     > 0.679 && (mapFilt[indexbbarCand ]).csv  >  0.679 &&
			(mapFilt[indexWdau1Cand ]).csv < 0.679   && (mapFilt[indexWdau2Cand ]).csv <0.679) ) 
		    continue;
		  
		  //if(verbose || (i==1930 && indexWdau1Cand==5 && indexWdau2Cand==6 && indexbCand==2 && indexbbarCand==1))
		  //cout << "Perm #" << counterPer << " [" << indexWdau1Cand << "," << indexWdau2Cand << "," << indexbCand << "," << indexbbarCand << " => STEP3" <<endl;
		  
		  
		  float WMass   = (myJetsFilt[indexWdau1Cand]+myJetsFilt[indexWdau2Cand]).M();
		  float TopMass = (myJetsFilt[indexWdau1Cand]+myJetsFilt[indexWdau2Cand]+myJetsFilt[indexbCand]).M();
		  
		  if( WMass  <  40 || WMass  > 120 ) continue;
		  if( TopMass< 120 || TopMass> 230 ) continue;
		  
		  
		  //if(verbose || (i==1930 && indexWdau1Cand==5 && indexWdau2Cand==6 && indexbCand==2 && indexbbarCand==1))
		  //cout << "Perm #" << counterPer << " [" << indexWdau1Cand << "," << indexWdau2Cand << "," << indexbCand << "," << indexbbarCand << " => STEP4" <<endl;	
		  
		  counterPer++;
	
		  p4HadB.SetPx( (myJetsFilt[indexbCand]).Px() );
		  p4HadB.SetPy( (myJetsFilt[indexbCand]).Py() );
		  p4HadB.SetPz( (myJetsFilt[indexbCand]).Pz() );
		  p4HadB.SetE(  (myJetsFilt[indexbCand]).E()  );
		  
		  p4HadP.SetPx( (myJetsFilt[indexWdau1Cand]).Px() );
		  p4HadP.SetPy( (myJetsFilt[indexWdau1Cand]).Py() );
		  p4HadP.SetPz( (myJetsFilt[indexWdau1Cand]).Pz() );
		  p4HadP.SetE(  (myJetsFilt[indexWdau1Cand]).E()  );

		  p4HadQ.SetPx( (myJetsFilt[indexWdau2Cand]).Px() );
		  p4HadQ.SetPy( (myJetsFilt[indexWdau2Cand]).Py() );
		  p4HadQ.SetPz( (myJetsFilt[indexWdau2Cand]).Pz() );
		  p4HadQ.SetE(  (myJetsFilt[indexWdau2Cand]).E()  );
		  
		  p4LepB.SetPx( (myJetsFilt[indexbbarCand]).Px() );
		  p4LepB.SetPy( (myJetsFilt[indexbbarCand]).Py() );
		  p4LepB.SetPz( (myJetsFilt[indexbbarCand]).Pz() );
		  p4LepB.SetE(  (myJetsFilt[indexbbarCand]).E()  );
		  
		  int statusPerm = fitter->fit(p4HadP, p4HadQ, p4HadB, p4LepB, p4Lepton, p4Neutrino, -1, CovarianceMatrix::kMuon);
		  if(statusPerm!=0) continue;
		  status++;
		  
		  //if(verbose || (i==1930 && indexWdau1Cand==5 && indexWdau2Cand==6 && indexbCand==2 && indexbbarCand==1))
		  //cout << "Perm #" << counterPer << " [" << indexWdau1Cand << "," << indexWdau2Cand << "," << indexbCand << "," << indexbbarCand << " => STEP5" <<endl;	
		  helHadR = -99; helLepR = -99; bTagR = -99;
		  helHadW = -99; helLepW = -99; bTagW = -99;



		  float chi2 = fitter->fitS();
		  float p    = fitter->fitProb();
		  if(useKin && p>0.)
		    logp =  -TMath::Log(p);
		    //logp   = TMath::Exp(-chi2);

		  float bTag = -99;
		  if(useBtag && fLikelihoodCsvB && fLikelihoodCsvNonB){
		    float btag1   = fLikelihoodCsvB->Eval( (mapFilt[indexbCand    ]).csv);
		    float btag2   = fLikelihoodCsvB->Eval( (mapFilt[indexbbarCand ]).csv);
		    float buntag1 = fLikelihoodCsvNonB->Eval( (mapFilt[indexWdau1Cand ]).csv);
		    float buntag2 = fLikelihoodCsvNonB->Eval( (mapFilt[indexWdau2Cand ]).csv);

		    bTag = btag1*btag2*buntag1*buntag2;
		    if( useKin && btag1>0. && btag2>0. && buntag1>0. && buntag2>0.){
		      logp += (-TMath::Log(btag1));
		      logp += (-TMath::Log(btag2));
		      logp += (-TMath::Log(buntag1));
		      logp += (-TMath::Log(buntag2));
		    }
		    else if( !useKin && btag1>0. && btag2>0. && buntag1>0. && buntag2>0.){
		      logp =  (-TMath::Log(btag1));
		      logp += (-TMath::Log(btag2));
		      logp += (-TMath::Log(buntag1));
		      logp += (-TMath::Log(buntag2));
		    }
		    else
		      logp = 999;	  
		  }

	

		  pat::Particle resHadP     = fitter->fittedHadP();
		  pat::Particle resHadQ     = fitter->fittedHadQ();
		  pat::Particle resHadB     = fitter->fittedHadB();
		  pat::Particle resLepB     = fitter->fittedLepB();
		  pat::Particle resLepton   = fitter->fittedLepton();
		  pat::Particle resNeutrino = fitter->fittedNeutrino();
	
		  TLorentzVector lvHadP(resHadP.px(),resHadP.py(),resHadP.pz(),resHadP.energy());
		  TLorentzVector lvHadQ(resHadQ.px(),resHadQ.py(),resHadQ.pz(),resHadQ.energy());
		  TLorentzVector lvHadB(resHadB.px(),resHadB.py(),resHadB.pz(),resHadB.energy());
		  TLorentzVector lvLepB(resLepB.px(),resLepB.py(),resLepB.pz(),resLepB.energy());
		  TLorentzVector lvLepton(resLepton.px(),resLepton.py(),resLepton.pz(),resLepton.energy());
		  TLorentzVector lvNeutrino(resNeutrino.px(),resNeutrino.py(),resNeutrino.pz(),resNeutrino.energy());
		  TLorentzVector lvHadW   = lvHadP   + lvHadQ ;
		  TLorentzVector lvHadTop = lvHadW   + lvHadB ;
		  TLorentzVector lvLepW   = lvLepton + lvNeutrino ;
		  TLorentzVector lvLepTop = lvLepW   + lvLepB ;
		  
		  float helHad = -99;
		  float helLep = -99;
		  if(useHel && fHelTopLept && fHelTopHadr){
		    
		    TLorentzVector bHad;
		    TLorentzVector w1Had;
		    TLorentzVector w2Had;
		    w1Had.SetPxPyPzE( lvHadP.Px(), lvHadP.Py(), lvHadP.Pz(), lvHadP.E());
		    w2Had.SetPxPyPzE( lvHadQ.Px(), lvHadQ.Py(), lvHadQ.Pz(), lvHadQ.E());
		    bHad.SetPxPyPzE(  lvHadB.Px(), lvHadB.Py(), lvHadB.Pz(), lvHadB.E());
		    TVector3 boostWHad(lvHadW.Px()/lvHadW.E(), lvHadW.Py()/lvHadW.E(), lvHadW.Pz()/lvHadW.E());
		    w1Had.Boost(-boostWHad);
		    w2Had.Boost(-boostWHad);
		    bHad.Boost( -boostWHad);
		    double cosChiStarTopHad = TMath::Cos(w1Had.Vect().Angle( -bHad.Vect() ));
		    helHad = fHelTopHadr->Eval( cosChiStarTopHad );
		    
		    TLorentzVector bLep;
		    TLorentzVector w1Lep;
		    TLorentzVector w2Lep;
		    w1Lep.SetPxPyPzE( lvLepton.Px(), lvLepton.Py(), lvLepton.Pz(), lvLepton.E());
		    w2Lep.SetPxPyPzE( lvNeutrino.Px(), lvNeutrino.Py(), lvNeutrino.Pz(), lvNeutrino.E());
		    bLep.SetPxPyPzE(  lvLepB.Px(), lvLepB.Py(), lvLepB.Pz(), lvLepB.E());
		    TVector3 boostWLep(lvLepW.Px()/lvLepW.E(), lvLepW.Py()/lvLepW.E(), lvLepW.Pz()/lvLepW.E());
		    w1Lep.Boost(-boostWLep);
		    w2Lep.Boost(-boostWLep);
		    bLep.Boost( -boostWLep);
		    double cosChiStarTopLep = TMath::Cos(w1Lep.Vect().Angle( -bLep.Vect() ));
		    helLep = fHelTopLept->Eval( cosChiStarTopLep );
		    
		    if( (useKin || useBtag) && helHad>0. && helLep>0. ){
		      logp += (-TMath::Log(helHad));
		      logp += (-TMath::Log(helLep));
		    }
		    else if( !(useKin || useBtag) && helHad>0. && helLep>0.){
		      logp =  (-TMath::Log(helHad));
		      logp += (-TMath::Log(helLep));
		    }
		    else
		      logp = 999;
		  }

		  float higgsMassProb = -99;
		  int counterPairPerm=0; 
		  unsigned int h1Perm=999; 
		  unsigned int h2Perm=999;
		  for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
		    float csv_k = mapFilt[k].csv>0 ? mapFilt[k].csv : 0.0 ;
		    float pt_k  = mapFilt[k].pt>0  ? mapFilt[k].pt  : 0.0 ;
		    if( !(csv_k > 0.679 && pt_k>bPtCut) ) continue;
		    
		    if( k==indexbCand || k==indexWdau1Cand || k==indexWdau2Cand || k==indexbbarCand) continue;
		    
		    if(counterPairPerm==0){
		      h1Perm   = k;
		    }
		    else if(counterPairPerm==1){
		      h2Perm   = k;
		    }
		    else{}
		    counterPairPerm++;
		  }
		  
		  if(useHMass && fSgnMass){

		    if(h1Perm!=999 && h2Perm!=999){
		      higgsMassProb = fSgnMass->Eval( (myJetsFilt[h1Perm]+myJetsFilt[h2Perm]).M() );

		      if( (useKin || useBtag || useHel) )  
			logp += (-TMath::Log(higgsMassProb));
		      else
			logp  = (-TMath::Log(higgsMassProb));			
		    }
		    else{
		      logp = 999;
		    }
		  }

		  

		  /////////////////////////////////////////////////////
		  float xa            = -99;
		  float xb            = -99;
		  float totalPt       = -99;
		  float planeCosTheta = -99;
		  float higgsOrient   = -99;
		  float topEnFrac     = -99;
		  float atopEnFrac    = -99;
		  float higgsPtStar   = -99;
		  float topMass       = -99;
		  float atopMass      = -99;
		  float higgsMass     = -99;
		  float topsCosDPhi   = -99;
		  float topLepPt  = -99;
		  float topLepEta = -99;
		  float topLepPhi = -99;
		  float topHadPt  = -99;
		  float topHadEta = -99;
		  float topHadPhi = -99;

		  if(computeCMSVariables && h1Perm!=999 && h2Perm!=999){

		    TLorentzVector TOP   = (vLepton_charge[0]==1)  ? lvLepTop : lvHadTop; 
		    TLorentzVector ATOP  = (vLepton_charge[0]==-1) ? lvLepTop : lvHadTop; 
		    TLorentzVector DIJET;
		    DIJET.SetPxPyPzE( (myJetsFilt[h1Perm]+myJetsFilt[h2Perm]).Px(), (myJetsFilt[h1Perm]+myJetsFilt[h2Perm]).Py(), (myJetsFilt[h1Perm]+myJetsFilt[h2Perm]).Pz(), (myJetsFilt[h1Perm]+myJetsFilt[h2Perm]).E());

		    TLorentzVector TOT = TOP+ATOP+DIJET;
		    TVector3 boostToCMS(TOT.Px()/TOT.E(), TOT.Py()/TOT.E(), TOT.Pz()/TOT.E());

		    if(abs(genTop.wdau1id)>6 &&  abs(genTbar.wdau1id)<6){
		      topLepPt   = (lvLepTop.Pt() -(topBLV+topW1LV+topW2LV).Pt() ) /(topBLV+topW1LV+topW2LV).Pt();
		      topLepEta  = (lvLepTop.Eta()-(topBLV+topW1LV+topW2LV).Eta() )/(topBLV+topW1LV+topW2LV).Eta();
		      topLepPhi  = (lvLepTop.Phi() -(topBLV+topW1LV+topW2LV).Phi() ) /(topBLV+topW1LV+topW2LV).Phi();
		      topHadPt   = (lvHadTop.Pt() -(atopBLV+atopW1LV+atopW2LV).Pt() ) /(atopBLV+atopW1LV+atopW2LV).Pt();
		      topHadEta =  (lvHadTop.Eta()-(atopBLV+atopW1LV+atopW2LV).Eta() )/(atopBLV+atopW1LV+atopW2LV).Eta();
		      topHadPhi =  (lvHadTop.Phi() -(atopBLV+atopW1LV+atopW2LV).Phi() ) /(atopBLV+atopW1LV+atopW2LV).Phi();
		    }
		    else if(abs(genTop.wdau1id)<6 &&  abs(genTbar.wdau1id)>6){
		      topLepPt   = (lvLepTop.Pt() -(atopBLV+atopW1LV+atopW2LV).Pt() ) /(atopBLV+atopW1LV+atopW2LV).Pt();
		      topLepEta  = (lvLepTop.Eta()-(atopBLV+atopW1LV+atopW2LV).Eta() )/(atopBLV+atopW1LV+atopW2LV).Eta();
		      topLepPhi  = (lvLepTop.Phi() -(atopBLV+atopW1LV+atopW2LV).Phi() ) /(atopBLV+atopW1LV+atopW2LV).Phi();
		      topHadPt   = (lvHadTop.Pt() -(topBLV+topW1LV+topW2LV).Pt() ) /(topBLV+topW1LV+topW2LV).Pt();
		      topHadEta =  (lvHadTop.Eta()-(topBLV+topW1LV+topW2LV).Eta() )/(topBLV+topW1LV+topW2LV).Eta();
		      topHadPhi =  (lvHadTop.Phi() -(topBLV+topW1LV+topW2LV).Phi() ) /(topBLV+topW1LV+topW2LV).Phi();
		    }
		    else{}

		    totalPt = TOT.Pt();
		    xa = (  TOT.Pz() + TOT.E() )/TMath::Sqrt(8000*8000.);
		    xb = ( -TOT.Pz() + TOT.E() )/TMath::Sqrt(8000*8000.);

		    TOP.Boost(-boostToCMS);
		    ATOP.Boost(-boostToCMS);
		    DIJET.Boost(-boostToCMS);

		    if(TOP.Mag()<=0 || ATOP.Mag()<=0) cout << "Error: null momentum" << endl;
		    TVector3 decayPlane    = ((TOP.Vect()).Cross(ATOP.Vect())).Unit();
		    TVector3 decayPlanePhi = decayPlane.Cross( TVector3(0,0,1) );
		    planeCosTheta = TMath::Cos( decayPlane.Angle( TVector3(0,0,1)) );
		    higgsOrient   = TMath::Cos( (DIJET.Vect()).Angle( decayPlanePhi ) );
		    topEnFrac     = TOP.E()/TOT.E();
		    atopEnFrac    = ATOP.E()/TOT.E();

		    higgsPtStar = DIJET.Vect().Mag();
		    higgsMass   = DIJET.M();
		    topMass     = TOP.M();
		    atopMass    = ATOP.M();

		    topsCosDPhi = TMath::Cos( (TOP.Vect()).Angle( ATOP.Vect() ) );

		  }

		  /////////////////////////////////////////////////////

		  
		  if(logp<minLogP){
		    minLogP = logp;
		    bestRank.clear();
		    bestRank.push_back( indexWdau1Cand);
		    bestRank.push_back( indexWdau2Cand);
		    bestRank.push_back( indexbCand);
		    bestRank.push_back( indexbbarCand);
		  }
		  
		  if(verbose) cout << "Perm #" << counterPer << " [" << indexWdau1Cand << "," << indexWdau2Cand << "," << indexbCand << "," << indexbbarCand << "] => log(p) = " << logp << endl;
		  
		  
		  if( matchB1>0 && matchW1>0 && matchW2>0 && matchB2>0){
		    foundMatch = 1;
		    chi2R = chi2;
		    probR = p;
		    helHadR = helHad;
		    helLepR = helLep;
		    bTagR   = bTag;
		    xaR=xa;
		    xbR=xb;
		    totalPtR = totalPt;
		    planeCosThetaR = planeCosTheta;
		    higgsOrientR   = higgsOrient;
		    topEnFracR     = topEnFrac;
		    atopEnFracR    = atopEnFrac;
		    higgsPtStarR = higgsPtStar;
		    higgsMassR   = higgsMass;
		    topMassR     = topMass;
		    atopMassR    = atopMass;
		    topsCosDPhiR =topsCosDPhi;
		    
		    topLepPtR  = topLepPt;
		    topLepEtaR = topLepEta;
		    topLepPhiR = topLepPhi;
		    topHadPtR  = topHadPt;
		    topHadEtaR = topHadEta;
		    topHadPhiR = topHadPhi;
		    
		    outTreeR->Fill();

		    
		  }
		  else{
		    chi2W = chi2;
		    probW = p;
		    helHadW = helHad;
		    helLepW = helLep;
		    bTagW   = bTag;

		    xaW=xa;
		    xbW=xb;
		    totalPtW = totalPt;
		    planeCosThetaW = planeCosTheta;
		    higgsOrientW   = higgsOrient;
		    topEnFracW     = topEnFrac;
		    atopEnFracW    = atopEnFrac;
		    higgsPtStarW = higgsPtStar;
		    higgsMassW   = higgsMass;
		    topMassW     = topMass;
		    atopMassW    = atopMass;
		    topsCosDPhiW =topsCosDPhi;
		    outTreeW->Fill();
		  }
	
		

		  
		}}}}}}
      //while ( std::next_permutation( indices, indices+myJetsFilt.size()  ) );
      

      
      //if(foundMatch==0 && fullTT==1){
      //cout << "Event " << i << " <==> " << counterPer << " / " << nPermutAll << ", status :" << status << endl;
      //cout  << mmatchW1 << "," << mmatchW2 << "," << mmatchB1 << "," << mmatchB2 << endl;
      //cout << "Csv : " << mapFilt[mmatchW1].csv << "," << mapFilt[mmatchW2].csv << "," << mapFilt[mmatchB1].csv << "," << mapFilt[mmatchB2].csv << endl;
      //cout << "Mass W " <<  (myJetsFilt[mmatchW1]+myJetsFilt[mmatchW2]).M() << endl;
      //}

      if(verbose){
	cout << "Number of permutation:  " << counterPer << endl;
	cout << "Smallest logLikelihood: " << minLogP
	     << " [" << bestRank[0] << "," <<  bestRank[1] << "," << bestRank[2] << "," <<  bestRank[3] << "]"
	     << endl;
      }



      hjet1.reset(); hjet2.reset();
      int counter=0; unsigned int h1=999; unsigned int h2=999;
      for(unsigned int k = 0; k < myJetsFilt.size(); k++){  

	if( status<1 /*minLogP>10*/ ) continue;
	float csv_k = mapFilt[k].csv>0 ? mapFilt[k].csv : 0.0 ;
	float pt_k  = mapFilt[k].pt>0  ? mapFilt[k].pt  : 0.0 ;
	if( !(csv_k > 0.679 && pt_k>30) ) continue;

	if( k==bestRank[0] || k==bestRank[1] ||  k==bestRank[2]  ||  k==bestRank[3]   ) continue;

	if(counter==0){
	  //jet1 = mapFilt[k];
	  hjet1.index = mapFilt[k].index; 
	  hjet1.pt    = mapFilt[k].pt; 
	  hjet1.eta   = mapFilt[k].eta; 
	  hjet1.phi   = mapFilt[k].phi; 
	  hjet1.mass  = mapFilt[k].mass; 
	  hjet1.csv   = mapFilt[k].csv; 
	  hjet1.topB  = mapFilt[k].topB; 
	  hjet1.topW  = mapFilt[k].topW; 
	  hjet1.atopW = mapFilt[k].atopW; 
	  hjet1.higgsB= mapFilt[k].higgsB; 
	  hjet1.flavor= mapFilt[k].flavor; 
	  hjet1.unc   = mapFilt[k].unc; 
	  h1   = k;
	}
	else if(counter==1){
	  //jet2 = mapFilt[k];
	  hjet2.index = mapFilt[k].index; 
	  hjet2.pt    = mapFilt[k].pt; 
	  hjet2.eta   = mapFilt[k].eta; 
	  hjet2.phi   = mapFilt[k].phi; 
	  hjet2.mass  = mapFilt[k].mass; 
	  hjet2.csv   = mapFilt[k].csv; 
	  hjet2.topB  = mapFilt[k].topB; 
	  hjet2.topW  = mapFilt[k].topW; 
	  hjet2.atopW = mapFilt[k].atopW; 
	  hjet2.higgsB= mapFilt[k].higgsB; 
	  hjet2.flavor= mapFilt[k].flavor; 
	  hjet2.unc   = mapFilt[k].unc; 
	  h2 = k;
	}
	else{}

	counter++;
      }
      int matchB    = 0;
      int matchBbar = 0;
      if(h1!=999 && h2!=999 && genBLV.Pt()>0 ){
	if(deltaR(myJetsFilt[h1], genBLV)<0.3) matchB++;
	if(deltaR(myJetsFilt[h2], genBLV)<0.3) matchB++;
      }
      if(h1!=999 && h2!=999 && genBbarLV.Pt()>0 ){
	if(deltaR(myJetsFilt[h1], genBbarLV)<0.3) matchBbar++;
	if(deltaR(myJetsFilt[h2], genBbarLV)<0.3) matchBbar++;
      }

      
      nPermut  = counterPer;
      nResJets = counter;
      matches  = matchB+matchBbar;
      outTree->Fill();


      /////////////////////////


      X1 = -99;
      X2 = -99;
      PtTTH   = -99;
      PhiTTH   = -99;
      BetaTTH   = -99;
      GammaTTH     = -99;
      GammaTT      = -99;
      DeltaTTH     = -99;
      PtHStar  = -99;
      MassTT   = -99;
      DphiHStar  = -99;
      BetaW   = -99;
      GammaW   = -99;
      DeltaW   = -99;
      BetaWLep = -99;
      GammaWLep = -99;
      DeltaWLep = -99;

      if(computeCMSVariables && genBLV.Pt()>0 && genBbarLV.Pt()>0 && topBLV.Pt()>0 && topW1LV.Pt()>0 && topW2LV.Pt()>0 && atopBLV.Pt()>0 && atopW1LV.Pt()>0 && atopW2LV.Pt()>0){

	TLorentzVector TOPLEP;
	TLorentzVector TOPHAD;
	TLorentzVector TOPHADW1;
	TLorentzVector TOPHADW2;
	TLorentzVector TOPHADB;
	TLorentzVector TOPLEPW1;
	TLorentzVector TOPLEPW2;
	TLorentzVector TOPLEPB;
	TLorentzVector HIGGS;

	int properEvent = 1;

	if( abs(genTop.wdau1id)>6 && abs(genTbar.wdau1id)<6){
	  TOPLEP.SetPxPyPzE( (topBLV+topW1LV+topW2LV).Px(), (topBLV+topW1LV+topW2LV).Py(), (topBLV+topW1LV+topW2LV).Pz(), (topBLV+topW1LV+topW2LV).E() );
	  TOPHAD.SetPxPyPzE( (atopBLV+atopW1LV+atopW2LV).Px(), (atopBLV+atopW1LV+atopW2LV).Py(), (atopBLV+atopW1LV+atopW2LV).Pz(), (atopBLV+atopW1LV+atopW2LV).E() );
	  TOPHADW1.SetPxPyPzE( atopW1LV.Px(), atopW1LV.Py(), atopW1LV.Pz(), atopW1LV.E());
	  TOPHADW2.SetPxPyPzE( atopW2LV.Px(), atopW2LV.Py(), atopW2LV.Pz(), atopW2LV.E());
	  TOPHADB.SetPxPyPzE( atopBLV.Px(),  atopBLV.Py(),   atopBLV.Pz(),  atopBLV.E());
	  TOPLEPW1.SetPxPyPzE( topW1LV.Px(), topW1LV.Py(), topW1LV.Pz(), topW1LV.E());
	  TOPLEPW2.SetPxPyPzE( topW2LV.Px(), topW2LV.Py(), topW2LV.Pz(), topW2LV.E());
	  TOPLEPB.SetPxPyPzE(  topBLV.Px(),  topBLV.Py(),   topBLV.Pz(), topBLV.E());
	}
	else if(abs(genTop.wdau1id)<6 && abs(genTbar.wdau1id)>6){
	  TOPHAD.SetPxPyPzE( (topBLV+topW1LV+topW2LV).Px(), (topBLV+topW1LV+topW2LV).Py(), (topBLV+topW1LV+topW2LV).Pz(), (topBLV+topW1LV+topW2LV).E() );
	  TOPLEP.SetPxPyPzE( (atopBLV+atopW1LV+atopW2LV).Px(), (atopBLV+atopW1LV+atopW2LV).Py(), (atopBLV+atopW1LV+atopW2LV).Pz(), (atopBLV+atopW1LV+atopW2LV).E() );
	  TOPHADW1.SetPxPyPzE( topW1LV.Px(), topW1LV.Py(), topW1LV.Pz(), topW1LV.E());
	  TOPHADW2.SetPxPyPzE( topW2LV.Px(), topW2LV.Py(), topW2LV.Pz(), topW2LV.E());
	  TOPHADB.SetPxPyPzE( topBLV.Px(),  topBLV.Py(),   topBLV.Pz(),  topBLV.E());
	  TOPLEPW1.SetPxPyPzE( atopW1LV.Px(), atopW1LV.Py(), atopW1LV.Pz(), atopW1LV.E());
	  TOPLEPW2.SetPxPyPzE( atopW2LV.Px(), atopW2LV.Py(), atopW2LV.Pz(), atopW2LV.E());
	  TOPLEPB.SetPxPyPzE( atopBLV.Px(),  atopBLV.Py(),   atopBLV.Pz(),  atopBLV.E());
	}
	
	else{properEvent=0;}
	HIGGS.SetPxPyPzE( (genBLV+genBbarLV).Px(), (genBLV+genBbarLV).Py(),(genBLV+genBbarLV).Pz(),(genBLV+genBbarLV).E());

	if(properEvent){


	  TLorentzVector TOT = TOPLEP+TOPHAD+HIGGS;
	  TVector3 boostToCMS(TOT.Px()/TOT.E(), TOT.Py()/TOT.E(), TOT.Pz()/TOT.E());
	  TVector3 boostToTopHadCMS(TOPHAD.Px()/TOPHAD.E(), TOPHAD.Py()/TOPHAD.E(), TOPHAD.Pz()/TOPHAD.E());
	  TVector3 boostToTopLepCMS(TOPLEP.Px()/TOPLEP.E(), TOPLEP.Py()/TOPLEP.E(), TOPLEP.Pz()/TOPLEP.E());

	  PtTTH  = TOT.Pt();
	  PhiTTH = TMath::Cos(TOT.Phi());
	  PzTTH  = TOT.Pz();
	  X1 = (  TOT.Pz() + TOT.E() )/TMath::Sqrt(8000*8000.);
	  X2 = ( -TOT.Pz() + TOT.E() )/TMath::Sqrt(8000*8000.);
	  
	//trasform
	float tmp1 = X1;
	float tmp2 = X2;
	X1 = TMath::Sqrt(tmp1*tmp2);
	X2 = tmp1-tmp2;

	TOPLEP.Boost(-boostToCMS);
	TOPHAD.Boost(-boostToCMS);
	HIGGS .Boost(-boostToCMS);

	PxT = TOPLEP.Px();
	PyT = TOPLEP.Py();
	PzT = TOPLEP.Pz();
	PxH = HIGGS.Px();
	PyH = HIGGS.Py();
	PzH = HIGGS.Pz();

	if(TOPLEP.Mag()<=0 || TOPHAD.Mag()<=0) cout << "Error: null momentum" << endl;

	TVector3 decayPlane    = ((TOPHAD.Vect()).Cross(TOPLEP.Vect())).Unit();
	TVector3 aux1 = (boostToCMS.Cross( TVector3(0,0,1) )).Unit();
	TVector3 aux2 = (boostToCMS.Cross( aux1 )).Unit();
	TVector3 aux3 = boostToCMS.Unit();

	TVector3 aux4( aux1.Px()*TMath::Cos(decayPlane.Angle(aux1)) , aux1.Py()*TMath::Cos(decayPlane.Angle(aux1)), aux1.Pz()*TMath::Cos(decayPlane.Angle(aux1)));
	TVector3 aux5( aux2.Px()*TMath::Cos(decayPlane.Angle(aux2)) , aux2.Py()*TMath::Cos(decayPlane.Angle(aux2)), aux2.Pz()*TMath::Cos(decayPlane.Angle(aux2)));

	//TVector3 decayPlanePhi = decayPlane.Cross( TVector3(0,0,1) );

	BetaTTH  = TMath::Cos(decayPlane.Angle(aux3));
	DeltaTTH = TMath::Cos((aux4+aux5).Angle(aux1));
	GammaTTH = TMath::Cos((HIGGS.Vect()).Angle( boostToCMS ));

	PtHStar   = HIGGS.P()/TOT.M();
	DphiHStar = TMath::Cos((HIGGS.Vect()).Angle(TOPLEP.Vect()));

	TOPLEP.Boost( -(TOPLEP+TOPHAD).BoostVector() ) ;
	TOPHAD.Boost( -(TOPLEP+TOPHAD).BoostVector() ) ;

	GammaTT  = TMath::Cos( (TOPLEP.Vect()).Angle( (TOPLEP+TOPHAD).BoostVector() ) );
	MassTT   = (TOPLEP+TOPHAD).M();

	//////////////////////////
	int chLep = 0;
	if(abs(genTop.wdau1id)==11  || abs(genTop.wdau1id)==13  || abs(genTop.wdau1id)==15 ||
	   abs(genTbar.wdau1id)==11 || abs(genTbar.wdau1id)==13 || abs(genTbar.wdau1id)==15  )
	  chLep = 1;
	else if(abs(genTop.wdau2id) ==11  || abs(genTop.wdau2id)==13  || abs(genTop.wdau2id)==15 ||
		abs(genTbar.wdau2id)==11  || abs(genTbar.wdau2id)==13 || abs(genTbar.wdau2id)==15  )
	  chLep = 2;
	else{}

	if( abs(genTop.wdau1id)>6 && abs(genTbar.wdau1id)<6){
	  TOPHAD.SetPxPyPzE( (atopBLV+atopW1LV+atopW2LV).Px(), (atopBLV+atopW1LV+atopW2LV).Py(), (atopBLV+atopW1LV+atopW2LV).Pz(), (atopBLV+atopW1LV+atopW2LV).E() );
	  TOPHADW1.SetPxPyPzE( atopW1LV.Px(), atopW1LV.Py(), atopW1LV.Pz(), atopW1LV.E());
	  TOPHADW2.SetPxPyPzE( atopW2LV.Px(), atopW2LV.Py(), atopW2LV.Pz(), atopW2LV.E());
	  TOPHADB.SetPxPyPzE( atopBLV.Px(),  atopBLV.Py(),   atopBLV.Pz(),  atopBLV.E());
	  TOPLEP.SetPxPyPzE( (topBLV+topW1LV+topW2LV).Px(), (topBLV+topW1LV+topW2LV).Py(), (topBLV+topW1LV+topW2LV).Pz(), (topBLV+topW1LV+topW2LV).E() );
	  TOPLEPW1.SetPxPyPzE( topW1LV.Px(), topW1LV.Py(), topW1LV.Pz(), topW1LV.E());
	  TOPLEPW2.SetPxPyPzE( topW2LV.Px(), topW2LV.Py(), topW2LV.Pz(), topW2LV.E());
	  TOPLEPB.SetPxPyPzE( topBLV.Px(),  topBLV.Py(),   topBLV.Pz(),  topBLV.E());
	}
	else{
	  TOPHAD.SetPxPyPzE( (topBLV+topW1LV+topW2LV).Px(), (topBLV+topW1LV+topW2LV).Py(), (topBLV+topW1LV+topW2LV).Pz(), (topBLV+topW1LV+topW2LV).E() );
	  TOPHADW1.SetPxPyPzE( topW1LV.Px(), topW1LV.Py(), topW1LV.Pz(), topW1LV.E());
	  TOPHADW2.SetPxPyPzE( topW2LV.Px(), topW2LV.Py(), topW2LV.Pz(), topW2LV.E());
	  TOPHADB.SetPxPyPzE( topBLV.Px(),  topBLV.Py(),   topBLV.Pz(),  topBLV.E());
	  TOPLEP.SetPxPyPzE( (atopBLV+atopW1LV+atopW2LV).Px(), (atopBLV+atopW1LV+atopW2LV).Py(), (atopBLV+atopW1LV+atopW2LV).Pz(), (atopBLV+atopW1LV+atopW2LV).E() );
	  TOPLEPW1.SetPxPyPzE( atopW1LV.Px(), atopW1LV.Py(), atopW1LV.Pz(), atopW1LV.E());
	  TOPLEPW2.SetPxPyPzE( atopW2LV.Px(), atopW2LV.Py(), atopW2LV.Pz(), atopW2LV.E());
	  TOPLEPB.SetPxPyPzE( atopBLV.Px(),  atopBLV.Py(),   atopBLV.Pz(),  atopBLV.E());
	}

	//TOPHAD.Boost(-boostToTopHadCMS);
	TOPHADW1.Boost(-boostToTopHadCMS);
	TOPHADW2.Boost(-boostToTopHadCMS);
	TOPHADB.Boost( -boostToTopHadCMS);

	TVector3 decayTPlane = (TOPHADB.Vect()).Cross( (TOPHADW1).Vect()  ).Unit();

	TVector3 boostToWHadCMS((TOPHADW1+TOPHADW2).Px()/(TOPHADW1+TOPHADW2).E(), (TOPHADW1+TOPHADW2).Py()/(TOPHADW1+TOPHADW2).E(), (TOPHADW1+TOPHADW2).Pz()/(TOPHADW1+TOPHADW2).E());
	TOPHADW1.Boost(-boostToWHadCMS);
	TOPHADW2.Boost(-boostToWHadCMS);

	//cout << TOPHADB.P() << endl;
	//cout << "W1: " << TOPHADW1.E() << " W1: " << TOPHADW2.E() << endl;

	BetaW     = TMath::Cos((TOPHADB.Vect()).Angle( boostToTopHadCMS  ));
	DeltaW    = TMath::Cos( decayTPlane.Angle( boostToTopHadCMS )  );
	GammaW    = TMath::Cos( (TOPHADW1.Vect()).Angle(boostToWHadCMS) );


	///////////////////////////

	TOPLEPW1.Boost(-boostToTopLepCMS);
	TOPLEPW2.Boost(-boostToTopLepCMS);
	TOPLEPB.Boost( -boostToTopLepCMS);

	TVector3 boostToWLepCMS((TOPLEPW1+TOPLEPW2).Px()/(TOPLEPW1+TOPLEPW2).E(), (TOPLEPW1+TOPLEPW2).Py()/(TOPLEPW1+TOPLEPW2).E(), (TOPLEPW1+TOPLEPW2).Pz()/(TOPLEPW1+TOPLEPW2).E());
	TOPLEPW1.Boost(-boostToWLepCMS);
	TOPLEPW2.Boost(-boostToWLepCMS);

	//cout << TOPLEPB.P() << endl;
	//cout << "W1: " << TOPLEPW1.E() << " W1: " << TOPLEPW2.E() << endl;

	BetaWLep     = TMath::Cos( (TOPLEPB.Vect()).Angle( boostToTopLepCMS  ));
	DeltaWLep    = 0.0;
	if( chLep==1)
	  GammaWLep    = TMath::Cos( (TOPLEPW1.Vect()).Angle(boostToWLepCMS) );
	else if( chLep==2)
	  GammaWLep    = TMath::Cos( (TOPLEPW2.Vect()).Angle(boostToWLepCMS) );
	else
	  GammaWLep = -99;
	///////



	///////
	genTree->Fill();
	}
      }



      /////////////////////////
      

      delete fitter;

    } // event loop

    mySamples->GetFile( currentName )->Close();
    
  } // samples loop


  fLikelihoods->Close();
  delete fHelTopLept; delete fHelTopHadr; delete fChi2; delete fSgnMass;
 
  return 0;

}
