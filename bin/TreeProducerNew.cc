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

#include "Bianchi/TTHStudies/interface/Samples.h"
#include "Bianchi/TTHStudies/interface/HelperFunctions.h"


#define GENJETDR 0.3
#define VERBOSE  false

using namespace std;


Float_t deltaR(TLorentzVector j1, TLorentzVector j2){
  float res = TMath::Sqrt(TMath::Power(j1.Eta()-j2.Eta(),2) + TMath::Power( TMath::ACos(TMath::Cos(j1.Phi()-j2.Phi())) ,2));
  return res;
}


int main(int argc, const char* argv[])
{

  std::cout << "TreeProducerNew" << std::endl;
  gROOT->SetBatch(true);
 
  gSystem->Load("libFWCoreFWLite");
  gSystem->Load("libDataFormatsFWLite");

  AutoLibraryLoader::enable();

  fwlite::TFileService fs = fwlite::TFileService("./root/treeProducerTEST.root");

  Float_t readerVars[15];
  TMVA::Reader* reader = getTMVAReader("BDTG", readerVars);

  TTree* genTree          = fs.make<TTree>("genTree","event tree");
  TTree* genJetLightTree  = fs.make<TTree>("genJetLightTree","event tree");
  TTree* genJetGluonTree  = fs.make<TTree>("genJetGluonTree","event tree");
  TTree* genJetHeavyTree  = fs.make<TTree>("genJetHeavyTree","event tree");
  TTree* genEventTree     = fs.make<TTree>("genEventTree","event tree");

  float BetaW  ;
  float GammaW  ;
  float BetaWLep  ;
  float GammaWLep  ;

  float sumEt;
  float et;
  float etReco;
  float phi;
  float phiReco;
  float px, pxReco, py, pyReco;
  float puWeight;

  float ptRecoHeavy;
  float ptRecoRegHeavy;
  float phiRecoHeavy;
  float etaRecoHeavy;
  float csvRecoHeavy;
  float ptHeavy;
  float genPtHeavy;
  float etaHeavy;
  float phiHeavy;
  float massRecoHeavy;
  int   flavorHeavy;

  float ptRecoLight;
  float phiRecoLight;
  float etaRecoLight;
  float csvRecoLight;
  float ptLight;
  float genPtLight;
  float etaLight;
  float phiLight;
  float massRecoLight;
  int   flavorLight;

  float ptRecoGluon;
  float phiRecoGluon;
  float etaRecoGluon;
  float csvRecoGluon;
  float genPtGluon;
  float massRecoGluon;
  int   flavorGluon;

  genTree->Branch("BetaW",               &BetaW,       "BetaW/F");
  genTree->Branch("GammaW",              &GammaW,      "GammaW/F");
  genTree->Branch("BetaWLep",            &BetaWLep,    "BetaWLep/F");
  genTree->Branch("GammaWLep",           &GammaWLep,   "GammaWLep/F");

  genEventTree->Branch("sumEt",          &sumEt,       "sumEt/F");
  genEventTree->Branch("et",             &et,          "et/F");
  genEventTree->Branch("etReco",         &etReco,      "etReco/F");
  genEventTree->Branch("phi",            &phi,         "phi/F");
  genEventTree->Branch("phiReco",        &phiReco,     "phiReco/F");
  genEventTree->Branch("px",             &px,          "px/F");
  genEventTree->Branch("pxReco",         &pxReco,      "pxReco/F");
  genEventTree->Branch("py",             &py,          "py/F");
  genEventTree->Branch("pyReco",         &pyReco,      "pyReco/F");
  genEventTree->Branch("puWeight",       &puWeight,    "puWeight/F");

  genJetHeavyTree->Branch("ptReco",      &ptRecoHeavy,  "ptReco/F");
  genJetHeavyTree->Branch("ptRecoReg",   &ptRecoRegHeavy,"ptRecoReg/F");
  genJetHeavyTree->Branch("etaReco",     &etaRecoHeavy, "etaReco/F");
  genJetHeavyTree->Branch("phiReco",     &phiRecoHeavy, "phiReco/F");
  genJetHeavyTree->Branch("csvReco",     &csvRecoHeavy, "csvReco/F");
  genJetHeavyTree->Branch("pt",          &ptHeavy,      "pt/F");
  genJetHeavyTree->Branch("ptGen",       &genPtHeavy,   "ptGen/F");
  genJetHeavyTree->Branch("eta",         &etaHeavy,     "eta/F");
  genJetHeavyTree->Branch("phi",         &phiHeavy,     "phi/F");
  genJetHeavyTree->Branch("massReco",    &massRecoHeavy,"massReco/F");
  genJetHeavyTree->Branch("flavor",      &flavorHeavy,  "flavor/I");

  genJetLightTree->Branch("ptReco",      &ptRecoLight,  "ptReco/F");
  genJetLightTree->Branch("etaReco",     &etaRecoLight, "etaReco/F");
  genJetLightTree->Branch("phiReco",     &phiRecoLight, "phiReco/F");
  genJetLightTree->Branch("csvReco",     &csvRecoLight, "csvReco/F");
  genJetLightTree->Branch("pt",          &ptLight,      "pt/F");
  genJetLightTree->Branch("ptGen",       &genPtLight,   "ptGen/F");
  genJetLightTree->Branch("eta",         &etaLight,     "eta/F");
  genJetLightTree->Branch("phi",         &phiLight,     "phi/F");
  genJetLightTree->Branch("massReco",    &massRecoLight,"massReco/F");
  genJetLightTree->Branch("flavor",      &flavorLight,  "flavor/I");

  genJetGluonTree->Branch("ptReco",      &ptRecoGluon,  "ptReco/F");
  genJetGluonTree->Branch("etaReco",     &etaRecoGluon, "etaReco/F");
  genJetGluonTree->Branch("phiReco",     &phiRecoGluon, "phiReco/F");
  genJetGluonTree->Branch("massReco",    &massRecoGluon,"massReco/F");
  genJetGluonTree->Branch("csvReco",     &csvRecoGluon, "csvReco/F");
  genJetGluonTree->Branch("ptGen",       &genPtGluon,   "ptGen/F");
  genJetGluonTree->Branch("flavor",      &flavorGluon,  "flavor/I");


  PythonProcessDesc builder(argv[1]);
 
  const edm::ParameterSet& in = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteInput");

  const edm::VParameterSet& samples = in.getParameter<edm::VParameterSet>("samples") ;

  std::string pathToFile   (in.getParameter<std::string>("pathToFile" ) );
  std::string ordering     (in.getParameter<std::string>("ordering" ) );
  double lumi              (in.getParameter<double>("lumi") );
  bool verbose             (in.getParameter<bool>("verbose") );

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

  string currentName0       = mySampleFiles[0];
  mySamples->OpenFile( currentName0 );
  TTree* currentTree0       = mySamples->GetTree( currentName0, "tree");

  for(unsigned int i = 0 ; i < mySampleFiles.size(); i++){
    
    string currentName       = mySampleFiles[i];

    mySamples->OpenFile( currentName );
    cout << "Opening file " << currentName << endl;
    TTree* currentTree       = mySamples->GetTree( currentName, "tree");
    cout << "Done!!" << endl;
    
    int nvlep, nhJets, naJets;

    Float_t hJet_pt           [999];
    Float_t hJet_eta          [999];
    Float_t hJet_phi          [999];
    Float_t hJet_e            [999];
    Float_t hJet_puJetIdL     [999];
    Float_t hJet_csv_nominal  [999];
    Float_t hJet_csv_upBC     [999];
    Float_t hJet_csv_downBC   [999];
    Float_t hJet_csv_upL      [999];
    Float_t hJet_csv_downL    [999];
    Float_t hJet_JECUnc       [999];
    Float_t hJet_genPt        [999];
    Float_t hJet_genEta       [999];
    Float_t hJet_genPhi       [999];
    Float_t hJet_flavour      [999];
    Float_t aJet_pt           [999];
    Float_t aJet_eta          [999];
    Float_t aJet_phi          [999];
    Float_t aJet_e            [999];
    Float_t aJet_puJetIdL     [999];
    Float_t aJet_csv_nominal  [999];
    Float_t aJet_csv_upBC     [999];
    Float_t aJet_csv_downBC   [999];
    Float_t aJet_csv_upL      [999];
    Float_t aJet_csv_downL    [999];
    Float_t aJet_JECUnc       [999];
    Float_t aJet_genPt        [999];
    Float_t aJet_genEta       [999];
    Float_t aJet_genPhi       [999];
    Float_t aJet_flavour      [999];

    genParticleInfo genB, genBbar;
    genTopInfo      genTop, genTbar;
    metInfo         METtype1p2corr;
    float           PUweight;

  
    currentTree->SetBranchAddress("nhJets",      &nhJets);
    currentTree->SetBranchAddress("naJets",      &naJets);
    currentTree->SetBranchAddress("hJet_pt",          hJet_pt);    
    currentTree->SetBranchAddress("hJet_eta",         hJet_eta);    
    currentTree->SetBranchAddress("hJet_phi",         hJet_phi);    
    currentTree->SetBranchAddress("hJet_e",           hJet_e);    
    currentTree->SetBranchAddress("hJet_puJetIdL",    hJet_puJetIdL);
    currentTree->SetBranchAddress("hJet_csv_nominal", hJet_csv_nominal);
    currentTree->SetBranchAddress("hJet_csv_upBC",    hJet_csv_upBC);
    currentTree->SetBranchAddress("hJet_csv_downBC",  hJet_csv_downBC);
    currentTree->SetBranchAddress("hJet_csv_upL",     hJet_csv_upL);
    currentTree->SetBranchAddress("hJet_csv_downL",   hJet_csv_downL);
    currentTree->SetBranchAddress("hJet_JECUnc",      hJet_JECUnc);
    currentTree->SetBranchAddress("hJet_genPt",       hJet_genPt);
    currentTree->SetBranchAddress("hJet_genEta",      hJet_genEta);
    currentTree->SetBranchAddress("hJet_genPhi",      hJet_genPhi);
    currentTree->SetBranchAddress("hJet_flavour",     hJet_flavour);
    currentTree->SetBranchAddress("aJet_pt",          aJet_pt);    
    currentTree->SetBranchAddress("aJet_eta",         aJet_eta);    
    currentTree->SetBranchAddress("aJet_phi",         aJet_phi);    
    currentTree->SetBranchAddress("aJet_e",           aJet_e);    
    currentTree->SetBranchAddress("aJet_puJetIdL",    aJet_puJetIdL);
    currentTree->SetBranchAddress("aJet_csv_nominal", aJet_csv_nominal);
    currentTree->SetBranchAddress("aJet_csv_upBC",    aJet_csv_upBC);
    currentTree->SetBranchAddress("aJet_csv_downBC",  aJet_csv_downBC);
    currentTree->SetBranchAddress("aJet_csv_upL",     aJet_csv_upL);
    currentTree->SetBranchAddress("aJet_csv_downL",   aJet_csv_downL);
    currentTree->SetBranchAddress("aJet_JECUnc",      aJet_JECUnc);
    currentTree->SetBranchAddress("aJet_genPt",       aJet_genPt);
    currentTree->SetBranchAddress("aJet_genEta",      aJet_genEta);
    currentTree->SetBranchAddress("aJet_genPhi",      aJet_genPhi);
    currentTree->SetBranchAddress("aJet_flavour",     aJet_flavour);

    currentTree->SetBranchAddress("genB",            &genB);
    currentTree->SetBranchAddress("genBbar",         &genBbar);
    currentTree->SetBranchAddress("genTop",          &genTop);
    currentTree->SetBranchAddress("genTbar",         &genTbar);
    currentTree->SetBranchAddress("METtype1p2corr",  &METtype1p2corr);
    currentTree->SetBranchAddress("PUweight",        &PUweight);

    int printed = 0;
    Long64_t nentries = currentTree->GetEntries();
    
    for (Long64_t i = 0; i < nentries ; i++){
      
      if(i%5000==0) cout << i << endl;
      if(i>20000){
	i = nentries;
	continue;
      }

      currentTree->GetEntry(i);

    
      // read top decay products from input files   
      TLorentzVector topBLV   (1,0,0,1);
      TLorentzVector topW1LV  (1,0,0,1); 
      TLorentzVector topW2LV  (1,0,0,1);
      TLorentzVector atopBLV  (1,0,0,1); 
      TLorentzVector atopW1LV (1,0,0,1); 
      TLorentzVector atopW2LV (1,0,0,1);
      TLorentzVector genBLV   (1,0,0,1);
      TLorentzVector genBbarLV(1,0,0,1);

      // set-up top decay products (if available from input file...)
      if(genTop.bmass>0){
	topBLV.SetPtEtaPhiM(  genTop.bpt,genTop.beta,genTop.bphi,genTop.bmass );
	topW1LV.SetPtEtaPhiM( genTop.wdau1pt,genTop.wdau1eta, genTop.wdau1phi,genTop.wdau1mass);
	topW2LV.SetPtEtaPhiM( genTop.wdau2pt,genTop.wdau2eta, genTop.wdau2phi,genTop.wdau2mass);
      }
      if(genTbar.bmass>0){
	atopBLV.SetPtEtaPhiM(  genTbar.bpt,genTbar.beta,genTbar.bphi,genTbar.bmass );
	atopW1LV.SetPtEtaPhiM( genTbar.wdau1pt,genTbar.wdau1eta, genTbar.wdau1phi,genTbar.wdau1mass);
	atopW2LV.SetPtEtaPhiM( genTbar.wdau2pt,genTbar.wdau2eta,genTbar.wdau2phi,genTbar.wdau2mass);
      }
      if(genB.mass>0 && (genB.momid==25 || genB.momid==23)){
	genBLV.SetPtEtaPhiM(genB.pt,genB.eta ,genB.phi, genB.mass );
      }
      if(genBbar.mass>0 && (genBbar.momid==25 || genBbar.momid==23)){
	genBbarLV.SetPtEtaPhiM(genBbar.pt,genBbar.eta ,genBbar.phi, genBbar.mass );
      }
  
  
      // define LV for the 6 (8) particles in ttbb (ttH) events
      TLorentzVector TOPHADW1(1,0,0,1);
      TLorentzVector TOPHADW2(1,0,0,1);
      TLorentzVector TOPHADB (1,0,0,1);
      TLorentzVector TOPLEPW1(1,0,0,1);
      TLorentzVector TOPLEPW2(1,0,0,1);
      TLorentzVector TOPLEPB (1,0,0,1);
      TLorentzVector HIGGSB1 (1,0,0,1);
      TLorentzVector HIGGSB2 (1,0,0,1);

      HIGGSB1.SetPtEtaPhiM( genBLV.Pt(),    genBLV.Eta(),    genBLV.Phi(),    genBLV.M());
      HIGGSB2.SetPtEtaPhiM( genBbarLV.Pt(), genBbarLV.Eta(), genBbarLV.Phi(), genBbarLV.M());

      if( abs(genTop.wdau1id)>6 && abs(genTbar.wdau1id)<6){
	TOPHADW1.SetPxPyPzE( atopW1LV.Px(), atopW1LV.Py(), atopW1LV.Pz(), atopW1LV.E());
	TOPHADW2.SetPxPyPzE( atopW2LV.Px(), atopW2LV.Py(), atopW2LV.Pz(), atopW2LV.E());
	TOPHADB.SetPxPyPzE( atopBLV.Px(),  atopBLV.Py(),   atopBLV.Pz(),  atopBLV.E());
	if( abs(genTop.wdau1id)==11 || abs(genTop.wdau1id)==13 || abs(genTop.wdau1id)==15 ){
	  TOPLEPW1.SetPxPyPzE  ( topW1LV.Px(), topW1LV.Py(), topW1LV.Pz(), topW1LV.E());
	  TOPLEPW2.SetPxPyPzE  ( topW2LV.Px(), topW2LV.Py(), topW2LV.Pz(), topW2LV.E());
	}
	else{
	  TOPLEPW1.SetPxPyPzE  ( topW2LV.Px(), topW2LV.Py(), topW2LV.Pz(), topW2LV.E());
	  TOPLEPW2.SetPxPyPzE  ( topW1LV.Px(), topW1LV.Py(), topW1LV.Pz(), topW1LV.E());
	}
	TOPLEPB.SetPxPyPzE(  topBLV.Px(),  topBLV.Py(),   topBLV.Pz(), topBLV.E());
      }
      else if(abs(genTop.wdau1id)<6 && abs(genTbar.wdau1id)>6){
	TOPHADW1.SetPxPyPzE( topW1LV.Px(), topW1LV.Py(), topW1LV.Pz(), topW1LV.E());
	TOPHADW2.SetPxPyPzE( topW2LV.Px(), topW2LV.Py(), topW2LV.Pz(), topW2LV.E());
	TOPHADB.SetPxPyPzE( topBLV.Px(),  topBLV.Py(),   topBLV.Pz(),  topBLV.E());
	if( abs(genTbar.wdau1id)==11 || abs(genTbar.wdau1id)==13 || abs(genTbar.wdau1id)==15 ){
	  TOPLEPW1.SetPxPyPzE  ( atopW1LV.Px(), atopW1LV.Py(), atopW1LV.Pz(), atopW1LV.E());
	  TOPLEPW2.SetPxPyPzE  ( atopW2LV.Px(), atopW2LV.Py(), atopW2LV.Pz(), atopW2LV.E());
	}
	else{
	  TOPLEPW1.SetPxPyPzE  ( atopW2LV.Px(), atopW2LV.Py(), atopW2LV.Pz(), atopW2LV.E());
	  TOPLEPW2.SetPxPyPzE  ( atopW1LV.Px(), atopW1LV.Py(), atopW1LV.Pz(), atopW1LV.E());
	}
	TOPLEPB.SetPxPyPzE( atopBLV.Px(),  atopBLV.Py(),   atopBLV.Pz(),  atopBLV.E());
      }      
      else if(abs(genTop.wdau1id)>6 && abs(genTbar.wdau1id)>6){
	if( abs(genTop.wdau1id)==11 || abs(genTop.wdau1id)==13 || abs(genTop.wdau1id)==15 ){
	  TOPHADW1.SetPxPyPzE  ( topW1LV.Px(), topW1LV.Py(), topW1LV.Pz(), topW1LV.E());
	  TOPHADW2.SetPxPyPzE  ( topW2LV.Px(), topW2LV.Py(), topW2LV.Pz(), topW2LV.E());
	}
	else{
	  TOPHADW1.SetPxPyPzE  ( topW2LV.Px(), topW2LV.Py(), topW2LV.Pz(), topW2LV.E());
	  TOPHADW2.SetPxPyPzE  ( topW1LV.Px(), topW1LV.Py(), topW1LV.Pz(), topW1LV.E());
	}
	TOPHADB.SetPxPyPzE( topBLV.Px(),  topBLV.Py(),   topBLV.Pz(),  topBLV.E());
	if( abs(genTbar.wdau1id)==11 || abs(genTbar.wdau1id)==13 || abs(genTbar.wdau1id)==15 ){
	  TOPLEPW1.SetPxPyPzE  ( atopW1LV.Px(), atopW1LV.Py(), atopW1LV.Pz(), atopW1LV.E());
	  TOPLEPW2.SetPxPyPzE  ( atopW2LV.Px(), atopW2LV.Py(), atopW2LV.Pz(), atopW2LV.E());
	}
	else{
	  TOPLEPW1.SetPxPyPzE  ( atopW2LV.Px(), atopW2LV.Py(), atopW2LV.Pz(), atopW2LV.E());
	  TOPLEPW2.SetPxPyPzE  ( atopW1LV.Px(), atopW1LV.Py(), atopW1LV.Pz(), atopW1LV.E());
	}
	TOPLEPB.SetPxPyPzE( atopBLV.Px(),  atopBLV.Py(),   atopBLV.Pz(),  atopBLV.E());
      }      
      else{}

      // conisder only semileptonic decays ...
      if(  !((abs(genTop.wdau1id)<6 &&  abs(genTbar.wdau1id)>6) ||  
	     (abs(genTop.wdau1id)>6 &&  abs(genTbar.wdau1id)<6))  )
	continue; 

      //////////////////////////////////////////////////
      // met
      //////////////////////////////////////////////////

      sumEt   = METtype1p2corr.sumet;
      etReco  = METtype1p2corr.et;
      phiReco = METtype1p2corr.phi;

      TVector3 recoMEt( etReco*TMath::Cos(phiReco), etReco*TMath::Sin(phiReco), 0.0) ;
      pxReco   = recoMEt.Px();
      pyReco   = recoMEt.Py();
      puWeight = PUweight;

      px       = TMath::Abs(TOPLEPW2.Py())>0 ? TOPLEPW2.Px()  : -999;
      py       = TMath::Abs(TOPLEPW2.Py())>0 ? TOPLEPW2.Py()  : -999;
      et       = TMath::Abs(TOPLEPW2.Py())>0 ? TOPLEPW2.Pt()  : -999;
      phi      = TMath::Abs(TOPLEPW2.Py())>0 ? TOPLEPW2.Phi() : -999;

      genEventTree->Fill();

      //////////////////////////////////////////////////
      // jets
      //////////////////////////////////////////////////

      vector<TLorentzVector> genHeavy;
      if( TMath::Abs(HIGGSB1.Py())>0) genHeavy.push_back( HIGGSB1);
      if( TMath::Abs(HIGGSB2.Py())>0) genHeavy.push_back( HIGGSB2);
      if( TMath::Abs(TOPHADB.Py())>0) genHeavy.push_back( TOPHADB);
      if( TMath::Abs(TOPLEPB.Py())>0) genHeavy.push_back( TOPLEPB);
      vector<TLorentzVector> genLight;
      if(  TMath::Abs(TOPHADW1.Py())>0) genLight.push_back( TOPHADW1);
      if(  TMath::Abs(TOPHADW2.Py())>0) genLight.push_back( TOPHADW2);
   

      // loop over jet collections
      for(int coll = 0 ; coll < 2 ; coll++){
	
	// loop over jets
	for(int hj = 0; hj < (coll==0 ? nhJets : naJets); hj++){	 

	  float ptGen = -99.;
	  if(coll==0 && hJet_genPt[hj]>0.){
	    ptGen = hJet_genPt[hj]*TMath::CosH(hJet_genEta[hj]);
	  }
	  if(coll==1 && aJet_genPt[hj]>0.){
	    ptGen = aJet_genPt[hj]*TMath::CosH(aJet_genEta[hj]);	  
	  }

	  float pt     = (coll==0) ? hJet_pt [hj]  : aJet_pt [hj];
	  float eta    = (coll==0) ? hJet_eta[hj]  : aJet_eta[hj];
	  float phi    = (coll==0) ? hJet_phi[hj]  : aJet_phi[hj];
	  float e      = (coll==0) ? hJet_e  [hj]  : aJet_e  [hj];
	  float m2     = e*e - pt*pt*TMath::CosH(eta)*TMath::CosH(eta);
	  if(m2<0) m2 = 0.; 
	  float m      = TMath::Sqrt( m2 ); 

	  // This is filled only for hJets
	  int id       = (coll==0) ? hJet_puJetIdL[hj] : aJet_puJetIdL[hj];
	  float JECUnc = (coll==0) ? hJet_JECUnc  [hj] : aJet_JECUnc  [hj];
	  int flavor   = (coll==0) ? hJet_flavour [hj] : aJet_flavour [hj];

	  // only jets passing ID...
	  if( id < 0.5 ) continue;	

	  TLorentzVector p4;
	  p4.SetPtEtaPhiM( pt, eta, phi, m );

	  // only jets above the pt cut
	  if( p4.Pt()<20 ) continue;

	  // only jets in acceptance
	  if( TMath::Abs(p4.Eta())>2.5 ) continue;

	  // for csv systematics
	  float csv_nominal =  (coll==0) ? hJet_csv_nominal[hj] : aJet_csv_nominal[hj];
	  float csv = TMath::Max(csv_nominal, float(0.0));

	  // matches a gen b quark from t->bW or H->bb ?
	  int matchesB = 0; 
	  unsigned int posB = 999;
	  for(unsigned int k = 0; k < genHeavy.size(); k++){  
	    if( deltaR(genHeavy[k], p4 ) < 0.3 ){
	      matchesB++;
	      posB = k;
	    }
	  }
	  // matches a gen udcs quark from W->qq' ?
	  int matchesL = 0;
	  unsigned int posL = 999;
	  for(unsigned int k = 0; k < genLight.size(); k++){  
	    if( deltaR(genLight[k], p4 ) < 0.3 ){
	      matchesL++;
	      posL = k;
	    }
	  }

	  // test regression
	  float ptReg = -99.;
	  if (coll==0){
	    ptReg = getRegressionEnergy("BDTG", reader, readerVars, currentTree0, i, (coll==0 ? hj : -hj-1), 0 );
	  }

	  // 1st case: the jet is matched to flavor 5
	  if( abs(flavor)==5 ){

	    csvRecoHeavy  = csv;
	    genPtHeavy    = ptGen;
	    flavorHeavy   = abs(flavor);
	    ptRecoHeavy   = e;
	    etaRecoHeavy  = eta;
	    phiRecoHeavy  = phi;
	    massRecoHeavy = m;
	    ptRecoRegHeavy = ptReg*TMath::CosH(eta);
	 
	    // no ambiguous matches...
	    if( matchesB==1 && matchesL==0 && posB!=999){
	      ptHeavy       = genHeavy[posB].E();
	      etaHeavy      = genHeavy[posB].Eta();
	      phiHeavy      = genHeavy[posB].Phi();
	    }else{
	      ptHeavy  = -99;
	      etaHeavy = -99;
	      phiHeavy = -99;
	    }


	    if( TMath::Abs(eta+1.865071) < 1e-04 ){
	      cout << "Entry " << i << endl;
	      cout << "Coll " << coll << endl;
	      cout << " hj " <<  hj << endl;
	      cout << eta << ", " << csv << ", " << genPtHeavy << ", " << ptRecoHeavy << endl;
	    }

	    genJetHeavyTree->Fill();
	  }

	  // 2nd case: the jet is matched to flavor < 5
	  else if( abs(flavor)<5 ){

	    csvRecoLight  = csv;
	    genPtLight    = ptGen;
	    flavorLight   = abs(flavor);
	    ptRecoLight   = e;
	    etaRecoLight  = eta;
	    phiRecoLight  = phi;
	    massRecoLight = m;
	 
	    // no ambiguous matches...
	    if( matchesB==0 && matchesL==1 && posL!=999){
	      ptLight       = genLight[posL].E();
	      etaLight      = genLight[posL].Eta();
	      phiLight      = genLight[posL].Phi();
	    }else{
	      ptLight  = -99;
	      etaLight = -99;
	      phiLight = -99;
	    }

	    genJetLightTree->Fill();
	  }

	  // 3rd case: gluon
	  else if( abs(flavor)==21 ){

	    csvRecoGluon  = csv;
	    genPtGluon    = ptGen;
	    flavorGluon   = abs(flavor);
	    ptRecoGluon   = e;
	    etaRecoGluon  = eta;
	    phiRecoGluon  = phi;
	    massRecoGluon = m;

	    genJetGluonTree->Fill();
	  }

	  else{}

	  
	}
      }

      //////////////////////////////////////////////////
      // gen decay
      //////////////////////////////////////////////////

      if(  TMath::Abs(TOPHADW2.Py())>0 &&  TMath::Abs(TOPHADW1.Py())>0 && TMath::Abs(TOPHADB.Py())>0 &&
	   TMath::Abs(TOPLEPW2.Py())>0 &&  TMath::Abs(TOPLEPW1.Py())>0 && TMath::Abs(TOPLEPB.Py())>0 ){

	TLorentzVector TOPLEP = TOPLEPW1+TOPLEPW2+TOPLEPB;
	TLorentzVector TOPHAD = TOPHADW1+TOPHADW2+TOPHADB;

	TVector3 boostToTopHadCMS = TOPHAD.BoostVector();
	TVector3 boostToTopLepCMS = TOPLEP.BoostVector();

	// first deal with TOPHAD...

	TOPHADW1.Boost(-boostToTopHadCMS);
	TOPHADW2.Boost(-boostToTopHadCMS);
	TOPHADB.Boost( -boostToTopHadCMS);
      
	TVector3 boostToWHadCMS = (TOPHADW1+TOPHADW2).BoostVector();
	TOPHADW1.Boost(-boostToWHadCMS);
	TOPHADW2.Boost(-boostToWHadCMS);	  
	  
	BetaW     = TMath::Cos((TOPHADB.Vect()).Angle( boostToTopHadCMS  ));
	GammaW    = i%2==0 ? TMath::Cos( (TOPHADW1.Vect()).Angle(boostToWHadCMS) ) : TMath::Cos( (TOPHADW2.Vect()).Angle(boostToWHadCMS) );  //average over W flavor (unobserved)
	  

	// then deal with TOPLEP...

	TOPLEPW1.Boost(-boostToTopLepCMS);
	TOPLEPW2.Boost(-boostToTopLepCMS);
	TOPLEPB.Boost( -boostToTopLepCMS);
	
	TVector3 boostToWLepCMS = (TOPLEPW1+TOPLEPW2).BoostVector();
	TOPLEPW1.Boost(-boostToWLepCMS);
	TOPLEPW2.Boost(-boostToWLepCMS);
	
	BetaWLep     = TMath::Cos( (TOPLEPB.Vect()).Angle( boostToTopLepCMS  ));       
	GammaWLep    = TMath::Cos( (TOPLEPW1.Vect()).Angle(boostToWLepCMS) );

	genTree->Fill();
      }
      
    } // event loop

    mySamples->GetFile( currentName )->Close();
    
  } // samples loop

  return 0;

}
