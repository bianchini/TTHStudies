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

  string target = argc>2 ? string(argv[2]) : "";
  string extraname = "TEST"+target;
  fwlite::TFileService fs = fwlite::TFileService("./root/treeProducer"+extraname+".root");

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

  Float_t regVarsFHeavy[14];
  Int_t   regVarsIHeavy[2];

  float eRecoHeavy;
  float ptRecoHeavy;
  float phiRecoHeavy;
  float etaRecoHeavy;
  float massRecoHeavy;
  float csvRecoHeavy;
  float eRecoRegHeavy;
  float ptRecoRegHeavy;
  float ePartHeavy;
  float ptPartHeavy;
  float etaPartHeavy;
  float phiPartHeavy;
  float eGenHeavy;
  float ptGenHeavy;
  float etaGenHeavy;
  float phiGenHeavy;
  int   flavorHeavy;

  float eRecoLight;
  float ptRecoLight;
  float phiRecoLight;
  float etaRecoLight;
  float massRecoLight;
  float csvRecoLight;
  float eRecoRegLight;
  float ptRecoRegLight;
  float ePartLight;
  float ptPartLight;
  float etaPartLight;
  float phiPartLight;
  float eGenLight;
  float ptGenLight;
  float etaGenLight;
  float phiGenLight;
  int   flavorLight;

  float eRecoGluon;
  float ptRecoGluon;
  float phiRecoGluon;
  float etaRecoGluon;
  float massRecoGluon;
  float csvRecoGluon;
  float eRecoRegGluon;
  float ptRecoRegGluon;
  float ePartGluon;
  float ptPartGluon;
  float etaPartGluon;
  float phiPartGluon;
  float eGenGluon;
  float ptGenGluon;
  float etaGenGluon;
  float phiGenGluon;
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

  addRegressionBranchesPerJet(genJetHeavyTree, regVarsFHeavy, regVarsIHeavy);
  genJetHeavyTree->Branch("e_rec",       &eRecoHeavy,    "e_rec/F");
  genJetHeavyTree->Branch("pt_rec",      &ptRecoHeavy,   "pt_rec/F");
  genJetHeavyTree->Branch("eta_rec",     &etaRecoHeavy,  "eta_rec/F");
  genJetHeavyTree->Branch("phi_rec",     &phiRecoHeavy,  "phi_rec/F");
  genJetHeavyTree->Branch("csv_rec",     &csvRecoHeavy,  "csv_rec/F");
  genJetHeavyTree->Branch("mass_rec",    &massRecoHeavy, "mass_rec/F");
  genJetHeavyTree->Branch("e_rec_reg",   &eRecoRegHeavy, "e_rec_reg/F");
  genJetHeavyTree->Branch("pt_rec_reg",  &ptRecoRegHeavy,"pt_rec_reg/F");
  genJetHeavyTree->Branch("e_part",      &ePartHeavy,    "e_part/F");
  genJetHeavyTree->Branch("pt_part",     &ptPartHeavy,   "pt_part/F");
  genJetHeavyTree->Branch("eta_part",    &etaPartHeavy,  "eta_part/F");
  genJetHeavyTree->Branch("phi_part",    &phiPartHeavy,  "phi_part/F");
  genJetHeavyTree->Branch("e_gen",       &eGenHeavy,     "e_gen/F");
  genJetHeavyTree->Branch("pt_gen",      &ptGenHeavy,    "pt_gen/F");
  genJetHeavyTree->Branch("eta_gen",     &etaGenHeavy,   "eta_gen/F");
  genJetHeavyTree->Branch("phi_gen",     &phiGenHeavy,   "phi_gen/F");
  genJetHeavyTree->Branch("flavor",      &flavorHeavy,   "flavor/I");

  genJetLightTree->Branch("e_rec",       &eRecoLight,    "e_rec/F");
  genJetLightTree->Branch("pt_rec",      &ptRecoLight,   "pt_rec/F");
  genJetLightTree->Branch("eta_rec",     &etaRecoLight,  "eta_rec/F");
  genJetLightTree->Branch("phi_rec",     &phiRecoLight,  "phi_rec/F");
  genJetLightTree->Branch("csv_rec",     &csvRecoLight,  "csv_rec/F");
  genJetLightTree->Branch("mass_rec",    &massRecoLight, "mass_rec/F");
  genJetLightTree->Branch("e_rec_reg",   &eRecoRegLight, "e_rec_reg/F");
  genJetLightTree->Branch("pt_rec_reg",  &ptRecoRegLight,"pt_rec_reg/F");
  genJetLightTree->Branch("e_part",      &ePartLight,    "e_part/F");
  genJetLightTree->Branch("pt_part",     &ptPartLight,   "pt_part/F");
  genJetLightTree->Branch("eta_part",    &etaPartLight,  "eta_part/F");
  genJetLightTree->Branch("phi_part",    &phiPartLight,  "phi_part/F");
  genJetLightTree->Branch("e_gen",       &eGenLight,     "e_gen/F");
  genJetLightTree->Branch("pt_gen",      &ptGenLight,    "pt_gen/F");
  genJetLightTree->Branch("eta_gen",     &etaGenLight,   "eta_gen/F");
  genJetLightTree->Branch("phi_gen",     &phiGenLight,   "phi_gen/F");
  genJetLightTree->Branch("flavor",      &flavorLight,   "flavor/I");

  genJetGluonTree->Branch("e_rec",       &eRecoGluon,    "e_rec/F");
  genJetGluonTree->Branch("pt_rec",      &ptRecoGluon,   "pt_rec/F");
  genJetGluonTree->Branch("eta_rec",     &etaRecoGluon,  "eta_rec/F");
  genJetGluonTree->Branch("phi_rec",     &phiRecoGluon,  "phi_rec/F");
  genJetGluonTree->Branch("csv_rec",     &csvRecoGluon,  "csv_rec/F");
  genJetGluonTree->Branch("mass_rec",    &massRecoGluon, "mass_rec/F");
  genJetGluonTree->Branch("e_rec_reg",   &eRecoRegGluon, "e_rec_reg/F");
  genJetGluonTree->Branch("pt_rec_reg",  &ptRecoRegGluon,"pt_rec_reg/F");
  genJetGluonTree->Branch("e_part",      &ePartGluon,    "e_part/F");
  genJetGluonTree->Branch("pt_part",     &ptPartGluon,   "pt_part/F");
  genJetGluonTree->Branch("eta_part",    &etaPartGluon,  "eta_part/F");
  genJetGluonTree->Branch("phi_part",    &phiPartGluon,  "phi_part/F");
  genJetGluonTree->Branch("e_gen",       &eGenGluon,     "e_gen/F");
  genJetGluonTree->Branch("pt_gen",      &ptGenGluon,    "pt_gen/F");
  genJetGluonTree->Branch("eta_gen",     &etaGenGluon,   "eta_gen/F");
  genJetGluonTree->Branch("phi_gen",     &phiGenGluon,   "phi_gen/F");
  genJetGluonTree->Branch("flavor",      &flavorGluon,   "flavor/I");

  PythonProcessDesc builder(argv[1]);
 
  const edm::ParameterSet& in = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteInput");

  const edm::VParameterSet& samples = in.getParameter<edm::VParameterSet>("samples") ;

  std::string pathToFile   (in.getParameter<std::string>("pathToFile" ) );
  std::string ordering     (in.getParameter<std::string>("ordering" ) );
  double lumi              (in.getParameter<double>     ("lumi") );
  bool verbose             (in.getParameter<bool>       ("verbose") );
  bool evalReg             (in.getParameter<bool>       ("evalReg") );
  int  maxnum              (in.getParameter<int>        ("maxnum") );

  Float_t readerVars[12];
  TMVA::Reader* reader = evalReg ? getTMVAReader("./root/weights/", "target-"+target, "BDTG", readerVars) : 0;


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
    
    int nvlep, nhJets, naJets, Vtype;

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

  
    currentTree->SetBranchAddress("Vtype",            &Vtype);
    currentTree->SetBranchAddress("nhJets",           &nhJets);
    currentTree->SetBranchAddress("naJets",           &naJets);
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
    cout << "Total entries: " << nentries << endl;

    for (Long64_t i = 0; i < nentries ; i++){
      
      if(i%5000==0) cout << i << "  (" << float(i)/float(nentries)*100. << " % completed)" << endl;

      if(i>maxnum){
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

	  float eGen   = -99.;
	  float ptGen  = -99.;
	  float etaGen = -99.;
	  float phiGen = -99.;
	  if(coll==0 && hJet_genPt[hj]>0.){
	    eGen   = hJet_genPt[hj]*TMath::CosH(hJet_genEta[hj]);
	    ptGen  = hJet_genPt[hj];
	    etaGen = hJet_genEta[hj];
	    phiGen = hJet_genPhi[hj];
	  }
	  if(coll==1 && aJet_genPt[hj]>0.){
	    eGen   = aJet_genPt[hj]*TMath::CosH(aJet_genEta[hj]);
	    ptGen  = aJet_genPt[hj];
	    etaGen = aJet_genEta[hj];
	    phiGen = aJet_genPhi[hj];
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

	  // 1st case: the jet is matched to flavor 5
	  if( abs(flavor)==5 ){

	    fillRegressionBranchesPerJet( regVarsFHeavy, regVarsIHeavy, currentTree0, i, (coll==0 ? hj : -hj-1), int(verbose));

	    eRecoHeavy    = e;
	    ptRecoHeavy   = pt;
	    etaRecoHeavy  = eta;
	    phiRecoHeavy  = phi;
	    massRecoHeavy = m;
	    csvRecoHeavy  = csv;
	    flavorHeavy   = abs(flavor);

	    eGenHeavy     = eGen;
	    ptGenHeavy    = ptGen;
	    etaGenHeavy   = etaGen;
	    phiGenHeavy   = phiGen;

	    // no ambiguous matches...
	    if( matchesB==1 && matchesL==0 && posB!=999 ){
	      ePartHeavy        = genHeavy[posB].E();
	      ptPartHeavy       = genHeavy[posB].Pt();
	      etaPartHeavy      = genHeavy[posB].Eta();
	      phiPartHeavy      = genHeavy[posB].Phi();
	    }
	    else{
	      ePartHeavy   = -99;
	      ptPartHeavy  = -99;
	      etaPartHeavy = -99;
	      phiPartHeavy = -99;
	    }

	    if( evalReg ){
	      float output  = -99;
	      getRegressionEnergy(output, "BDTG", reader, readerVars, currentTree0, i, (coll==0 ? hj : -hj-1), int(verbose));
	      eRecoRegHeavy  = output*TMath::CosH( etaRecoHeavy );
	      ptRecoRegHeavy = output;
	    }
	    else{
	      eRecoRegHeavy  = -99;
	      ptRecoRegHeavy = -99;
	    }

	    genJetHeavyTree->Fill();
	  }

	  // 2nd case: the jet is matched to flavor < 5
	  else if( abs(flavor)<5 ){

	    eRecoLight    = e;
	    ptRecoLight   = pt;
	    etaRecoLight  = eta;
	    phiRecoLight  = phi;
	    massRecoLight = m;
	    csvRecoLight  = csv;
	    flavorLight   = abs(flavor);

	    eGenLight     = eGen;
	    ptGenLight    = ptGen;
	    etaGenLight   = etaGen;
	    phiGenLight   = phiGen;

	    // no ambiguous matches...
	    if( matchesB==0 && matchesL==1 && posL!=999 ){
	      ePartLight        = genLight[posL].E();
	      ptPartLight       = genLight[posL].Pt();
	      etaPartLight      = genLight[posL].Eta();
	      phiPartLight      = genLight[posL].Phi();
	    }
	    else{
	      ePartLight   = -99;
	      ptPartLight  = -99;
	      etaPartLight = -99;
	      phiPartLight = -99;
	    }

	    genJetLightTree->Fill();

	  }

	  // 3rd case: gluon
	  else if( abs(flavor)==21 ){

	    eRecoGluon    = e;
	    ptRecoGluon   = pt;
	    etaRecoGluon  = eta;
	    phiRecoGluon  = phi;
	    massRecoGluon = m;
	    csvRecoGluon  = csv;
	    flavorGluon   = abs(flavor);

	    eGenGluon     = eGen;
	    ptGenGluon    = ptGen;
	    etaGenGluon   = etaGen;
	    phiGenGluon   = phiGen;	  

	    genJetGluonTree->Fill();
	  }

	  else{}

	  
	}
      }

      //////////////////////////////////////////////////
      // gen decay
      //////////////////////////////////////////////////

      if(  TMath::Abs(TOPHADW2.Py())>0 &&  TMath::Abs(TOPHADW1.Py())>0 && TMath::Abs(TOPHADB.Py())>0 &&
	   TMath::Abs(TOPLEPW2.Py())>0 &&  TMath::Abs(TOPLEPW1.Py())>0 && TMath::Abs(TOPLEPB.Py())>0 &&
	   (Vtype==2 || Vtype==3)
	   ){

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

	/*
	cout << "*** check..." << endl;
	cout << "- TOPLEPW1 M=" << TOPLEPW1.M() << endl;
	cout << "- TOPLEPW2 M=" << TOPLEPW2.M() << endl;
	cout << "- GammaWLep = " << GammaWLep << endl; 
	cout << "************" << endl;
	*/

	genTree->Fill();
      }
      
    } // event loop

    mySamples->GetFile( currentName )->Close();
    
  } // samples loop

  return 0;

}
