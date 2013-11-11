#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include <cstdlib>
#include <iostream> 
#include <fstream>
#include <map>
#include <string>
#include <algorithm>

#include "TSystem.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "THStack.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraph.h"
#include "Math/GenVector/LorentzVector.h"
#include "TLorentzVector.h"
#include "TVectorD.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TMatrixT.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GSLMCIntegrator.h"
#include "Math/AllIntegrationTypes.h"
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

#include "Bianchi/TTHStudies/interface/MEIntegratorNew.h"
#include "Bianchi/TTHStudies/interface/Samples.h"
#include "Bianchi/TTHStudies/interface/MECombination.h"
#include "Bianchi/TTHStudies/interface/HelperFunctions.h"

// dR distance for reco-gen matching
#define GENJETDR  0.3

// maximum number of VEGAS evaluation
#define MAX_REEVAL_TRIES 3

// 3.141...
#define PI TMath::Pi()

// maximum number of permutations per event
#define NMAXPERMUT 30

// maximum number of mass points for likelihood scan
#define NMAXMASS   20

// maximum number of four-vectors per event saved in output 
#define NMAXJETS   10

// if one, MC jets are already corrected for JER
#define JERCORRECTED 1


using namespace std;


int main(int argc, const char* argv[])
{

  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
  /* @@@@@@@@@@@@@@@@@@@@@@@@ FWLITE @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */

 
  std::cout << "MEAnalysisNew" << std::endl;
  gROOT->SetBatch(true);
 
  gSystem->Load("libFWCoreFWLite");
  gSystem->Load("libDataFormatsFWLite");

  AutoLibraryLoader::enable();


  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
  /* @@@@@@@@@@@@@@@@@@@@@@@ CONFIGURATION @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */


  PythonProcessDesc builder(argv[1]);
  const edm::ParameterSet& in = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteInput");

  // SAMPLES
  const edm::VParameterSet& samples = in.getParameter<edm::VParameterSet>("samples") ;
  std::string outFileName( in.getParameter<std::string>  ("outFileName" ) );
  std::string pathToFile ( in.getParameter<std::string>  ("pathToFile" ) );
  std::string ordering   ( in.getParameter<std::string>  ("ordering" ) );
  std::string pathToTF   ( in.getParameter<std::string>  ("pathToTF") );
  std::string pathToCP   ( in.getParameter<std::string>  ("pathToCP") );
  bool   verbose         ( in.getParameter<bool>         ("verbose" ) );
  
  // PARAMETERS
  double lumi               ( in.getParameter<double>("lumi") );
  float  MH                 ( in.getUntrackedParameter<double> ("MH",     125.));
  float  MT                 ( in.getUntrackedParameter<double> ("MT",    174.3));
  float  MW                 ( in.getUntrackedParameter<double> ("MW",    80.19));
  float  MwL                ( in.getUntrackedParameter<double> ("MwL",      60));
  float  MwH                ( in.getUntrackedParameter<double> ("MwH",     100));
  float  btag_prob_cut_6jets( in.getUntrackedParameter<double> ("btag_prob_cut_6jets",    0.));
  float  btag_prob_cut_5jets( in.getUntrackedParameter<double> ("btag_prob_cut_5jets",    0.));
  float  btag_prob_cut_4jets( in.getUntrackedParameter<double> ("btag_prob_cut_4jets",    0.));
  double maxChi2_           ( in.getUntrackedParameter<double> ("maxChi2", 2.5));
  vector<double> massesH    ( in.getParameter<vector<double> > ("massesH"));
  vector<double> massesT    ( in.getParameter<vector<double> > ("massesT"));

  // FLAGS
  int   switchoffOL      ( in.getUntrackedParameter<int>    ("switchoffOL",     0));
  int   speedup          ( in.getUntrackedParameter<int>    ("speedup",         0));
  int   doubleGaussianB  ( in.getUntrackedParameter<int>    ("doubleGaussianB", 1));
  int   useBtag          ( in.getUntrackedParameter<int>    ("useBtag",         0));
  int   selectByBTagShape( in.getUntrackedParameter<int>    ("selectByBTagShape",0));
  int   doTypeBTag4      ( in.getUntrackedParameter<int>    ("doTypeBTag4", 0));
  int   doTypeBTag5      ( in.getUntrackedParameter<int>    ("doTypeBTag5", 0));
  int   doTypeBTag6      ( in.getUntrackedParameter<int>    ("doTypeBTag6", 0));
  int   doType0          ( in.getUntrackedParameter<int>    ("doType0", 0));
  int   doType1          ( in.getUntrackedParameter<int>    ("doType1", 0));
  int   doType2          ( in.getUntrackedParameter<int>    ("doType2", 0));
  int   doType3          ( in.getUntrackedParameter<int>    ("doType3", 0));
  int   doType6          ( in.getUntrackedParameter<int>    ("doType6", 0));
  int   doType7          ( in.getUntrackedParameter<int>    ("doType7", 0));
  int   doType0ByBTagShape(in.getUntrackedParameter<int>    ("doType0ByBTagShape", 0));
  int   doType1ByBTagShape(in.getUntrackedParameter<int>    ("doType1ByBTagShape", 0));
  int   doType2ByBTagShape(in.getUntrackedParameter<int>    ("doType2ByBTagShape", 0));
  int   doType3ByBTagShape(in.getUntrackedParameter<int>    ("doType3ByBTagShape", 0));
  int   doType6ByBTagShape(in.getUntrackedParameter<int>    ("doType6ByBTagShape", 0));
  int   useME            ( in.getParameter<int>             ("useME")     );
  int   useJac           ( in.getParameter<int>             ("useJac")    );
  int   useMET           ( in.getParameter<int>             ("useMET")    );
  int   useTF            ( in.getParameter<int>             ("useTF")     );
  int   usePDF           ( in.getParameter<int>             ("usePDF")    );
  int   norm             ( in.getUntrackedParameter<int>    ("norm",      0));
  int   hypo             ( in.getUntrackedParameter<int>    ("hypo",      0));
  int   SoB              ( in.getUntrackedParameter<int>    ("SoB",       1));
  int   doJERbias        ( in.getUntrackedParameter<int>    ("doJERbias", 0));

  int   doCSVup          ( in.getUntrackedParameter<int>    ("doCSVup",   0));
  int   doCSVdown        ( in.getUntrackedParameter<int>    ("doCSVdown", 0));
  int   doJECup          ( in.getUntrackedParameter<int>    ("doJECup",   0));
  int   doJECdown        ( in.getUntrackedParameter<int>    ("doJECdown", 0));
  int   doJERup          ( in.getUntrackedParameter<int>    ("doJERup",   0));
  int   doJERdown        ( in.getUntrackedParameter<int>    ("doJERdown", 0));
  int   fixNumEvJob      ( in.getUntrackedParameter<int>    ("fixNumEvJob",1));
  vector<string> functions(in.getParameter<vector<string> > ("functions"));
  vector<int>    evLimits (in.getParameter<vector<int> >    ("evLimits"));

  int   print            ( in.getParameter<int>             ("printout")   );
  int   debug            ( in.getParameter<int>             ("debug")      );


  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
  /* @@@@@@@@@@@@@@@@@@@@@@@@@ INITIALIZE @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */

  // flag to discriminate data/MC
  bool isMC = true;

  // upper and lower bounds to be processed
  int evLow  = evLimits[0];
  int evHigh = evLimits[1];

  TStopwatch* clock = new TStopwatch();

  // file with b-tag pdf
  TFile* fCP = TFile::Open(pathToCP.c_str(),"READ");

  // b-tag pdf for b-quark ('b'), c-quark ('c'), and light jets ('l')
  map<string,TH1F*> btagger; 
  if( useBtag && fCP!=0 ){
    btagger["b_Bin0"] = fCP->Get("csv_b_Bin0__csvReco")!=0 ? (TH1F*)fCP->Get("csv_b_Bin0__csvReco") : 0;
    btagger["b_Bin1"] = fCP->Get("csv_b_Bin1__csvReco")!=0 ? (TH1F*)fCP->Get("csv_b_Bin1__csvReco") : 0;
    btagger["c_Bin0"] = fCP->Get("csv_c_Bin0__csvReco")!=0 ? (TH1F*)fCP->Get("csv_c_Bin0__csvReco") : 0;
    btagger["c_Bin1"] = fCP->Get("csv_c_Bin1__csvReco")!=0 ? (TH1F*)fCP->Get("csv_c_Bin1__csvReco") : 0;
    btagger["l_Bin0"] = fCP->Get("csv_l_Bin0__csvReco")!=0 ? (TH1F*)fCP->Get("csv_l_Bin0__csvReco") : 0;
    btagger["l_Bin1"] = fCP->Get("csv_l_Bin1__csvReco")!=0 ? (TH1F*)fCP->Get("csv_l_Bin1__csvReco") : 0;
  }
  else if( useBtag && fCP!=0 ){
    cout << "Cound not find " << pathToCP << ": exit" << endl;
    return 0;
  }

  // Higgs mass values for scan
  const int nHiggsMassPoints  = massesH.size();
  double mH[nHiggsMassPoints];
  for( unsigned int m = 0; m < massesH.size() ; m++)
    mH[m] = massesH[m];

  // Top mass values for scan
  const int nTopMassPoints  = massesT.size();
  double mT[nTopMassPoints]; 
  for( unsigned int m = 0; m < massesT.size() ; m++)
    mT[m] = massesT[m];


  // not supported...
  if( nHiggsMassPoints>1 && nTopMassPoints>1){
    cout << "Cannot handle two mass scans at the same time... return." << endl;
    return 1;
  }
  if( nHiggsMassPoints>NMAXMASS || nTopMassPoints>NMAXMASS){
    cout << "Too many mass points required... return" << endl;
    return 1;
  }


  // configure MEIntegrator
  MEIntegratorNew* meIntegrator = new MEIntegratorNew( pathToTF , 4 , int(verbose));
  if( norm == 0)
    meIntegrator->setWeightNorm( MEIntegratorNew::None );
  else if( norm==1 )
    meIntegrator->setWeightNorm( MEIntegratorNew::xSec );
  else if( norm==2 )
    meIntegrator->setWeightNorm( MEIntegratorNew::Acc );
  else{
    cout << "Unsupported normalization... exit" << endl;
    delete meIntegrator;
    return 0;
  }

  // set normalization formulas ( not used if norm==0 )
  meIntegrator->setNormFormulas( TString(functions[0].c_str()),  
				 TString(functions[1].c_str()),  
				 TString(functions[2].c_str()),
				 TString(functions[3].c_str()),  
				 TString(functions[4].c_str()),  
				 TString(functions[5].c_str())
				 );
  
  // initialize top and W mass
  meIntegrator->setTopMass( MT , MW );

  // configure ME calculation
  meIntegrator->setUseME (useME);
  meIntegrator->setUseJac(useJac);
  meIntegrator->setUseMET(useMET);
  meIntegrator->setUseTF (useTF);
  meIntegrator->setUsePDF(usePDF);

  // use double-gaussian for b quark energy TF
  meIntegrator->setUseRefinedTF(doubleGaussianB);

  // use nominal TF from ROOT file
  meIntegrator->initTFparameters(1.0,1.0,1.0,1.0, 1.0);
  if(switchoffOL){
    meIntegrator->switchOffOL(); 
    cout << "*** Switching off OpenLoops to speed-up the calculation ***" << endl;
  }



  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
  /* @@@@@@@@@@@@@@@@@@@@@@@ OUTPUT @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */


  // clean output file (if any)
  gSystem->Exec(("rm "+outFileName).c_str());

  // output file
  TFile* fout_tmp = TFile::Open(outFileName.c_str(),"UPDATE");

  // total event counter for normalization
  TH1F*  hcounter = new TH1F("hcounter","",1,0,1);
  int events_     = 0;
  
  // output tree
  TTree* tree  = new TTree("tree","");

  // counts how many events have been analyzed (to have fixed-size jobs)
  int counter_;
  // number of jet-quark permutations (for signal and bkg hypotheses)
  int nPermut_, nPermut_alt_;
  // number of (Higgs/Top) mass points
  int nMassPoints_;
  // total number integration per event (nPermut*nMassPoints)
  int nTotInteg_,nTotInteg_alt_;
  // number of matches to higgs quarks among tagged jets
  int matchesH_;
  // number of matches to W quarks among un-tagged jets
  int matchesW_;
  // number of matches to top quarks among tagged jets
  int matchesT_;
  // number of matches to higgs quarks among all jets
  int matchesHAll_;
  // number of matches to W quarks among all jets
  int matchesWAll_;
  // number of matches to top quarks among all jets
  int matchesTAll_;
  // count how many quarks from W decay overlap by dR<0.5
  int overlapLight_;
  // count how many b-quarks overlap by dR<0.5
  int overlapHeavy_;
  // integration type
  int type_;
  // num. of b-hadrons and c-quarks
  int nSimBs_; //, nC_, nCTop_;
  // num of b-hadrons inside jets
  int nMatchSimBs_;
  // type-dependent flags
  int flag_type0_;
  int flag_type1_;
  int flag_type2_;
  int flag_type3_;
  int flag_type4_;
  int flag_type6_;

  // event-wise **ME** probability (summed over permutations)
  // ...for ttH...
  float probAtSgn_;
  // ...for ttbb...
  float probAtSgn_alt_;

  // event-wise **ME X btag** probability (summed over permutations)
  // ...for ttH...
  float probAtSgn_ttbb_;
  // ...for tt(bb,jj,bj,cc)...
  float probAtSgn_alt_ttbb_;
  float probAtSgn_alt_ttjj_;
  float probAtSgn_alt_ttbj_;
  float probAtSgn_alt_ttcc_;

  // per-permutation and mass value **ME** probability
  float probAtSgn_permut_       [NMAXPERMUT*NMAXMASS];
  float probAtSgn_alt_permut_   [NMAXPERMUT*NMAXMASS];
  float probAtSgnErr_permut_    [NMAXPERMUT*NMAXMASS];
  float probAtSgnErr_alt_permut_[NMAXPERMUT*NMAXMASS];

  // per-permutation **btag** probability
  float probAtSgn_bb_permut_[NMAXPERMUT];
  float probAtSgn_bj_permut_[NMAXPERMUT];
  float probAtSgn_cc_permut_[NMAXPERMUT];
  float probAtSgn_jj_permut_[NMAXPERMUT];

  // masses to be scanned
  float mH_scan_[NMAXMASS];
  float mT_scan_[NMAXMASS];

  // a flag for the event type
  int Vtype_;
  // event-dependent weight (for normalization)
  float weight_;
  // cpu time
  float time_;
  // event information
  EventInfo EVENT_;
  // num of PVs
  int nPVs_;
  // pu reweighting
  float PUweight_, PUweightP_, PUweightM_;
  // trigger turn-on
  float trigger_;
  // trigger bits (data)
  bool triggerFlags_[70];
  // number of gen jets
  float lheNj_;

  // lepton kinematic (at most two leptons)
  int   nLep_;
  float lepton_pt_    [2];
  float lepton_eta_   [2];
  float lepton_phi_   [2];
  float lepton_m_     [2];
  float lepton_charge_[2];
  float lepton_rIso_  [2];
  int   lepton_type_  [2];

  // met kinematic
  float MET_pt_;
  float MET_phi_;
  float MET_sumEt_;

  // jet kinematics (as passed via **jets** collection)
  int nJet_;
  float jet_pt_  [NMAXJETS];
  float jet_eta_ [NMAXJETS];
  float jet_phi_ [NMAXJETS];
  float jet_m_   [NMAXJETS];
  float jet_csv_ [NMAXJETS];

  // number of selected jets passing CSV L,M,T
  int numBTagL_, numBTagM_, numBTagT_;
  // nummber of selected jets
  int numJets_;
  // btag likelihood ratio
  float btag_LR_;

  // permutation -> jets association
  int   perm_to_jet_    [NMAXPERMUT];
  int   perm_to_jet_alt_[NMAXPERMUT];
  // permutation -> gen association
  int   perm_to_gen_     [NMAXPERMUT];
  int   perm_to_gen_alt_ [NMAXPERMUT];

  tree->Branch("counter",      &counter_,       "counter/I");
  tree->Branch("nPermut_s",    &nPermut_,       "nPermut_s/I");
  tree->Branch("nPermut_b",    &nPermut_alt_,   "nPermut_b/I");
  tree->Branch("nMassPoints",  &nMassPoints_,   "nMassPoints/I");  
  tree->Branch("nTotInteg_s",  &nTotInteg_,     "nTotInteg_s/I");  
  tree->Branch("nTotInteg_b",  &nTotInteg_alt_, "nTotInteg_b/I");  
  tree->Branch("matchesH",     &matchesH_,      "matchesH/I");
  tree->Branch("matchesW",     &matchesW_,      "matchesW/I");
  tree->Branch("matchesT",     &matchesT_,      "matchesT/I");
  tree->Branch("matchesHAll",  &matchesHAll_,   "matchesHAll/I");
  tree->Branch("matchesWAll",  &matchesWAll_,   "matchesWAll/I");
  tree->Branch("matchesTAll",  &matchesTAll_,   "matchesTAll/I");
  tree->Branch("overlapLight", &overlapLight_,  "overlapLight/I");
  tree->Branch("overlapHeavy", &overlapHeavy_,  "overlapHeavy/I");
  tree->Branch("type",         &type_,          "type/I");
  tree->Branch("nSimBs",       &nSimBs_,        "nSimBs/I");
  tree->Branch("nMatchSimBs",  &nMatchSimBs_,   "nMatchSimBs/I");
  tree->Branch("weight",       &weight_,        "weight/F");
  tree->Branch("time",         &time_,          "time/F");
  tree->Branch("flag_type0",   &flag_type0_,    "flag_type0/I");
  tree->Branch("flag_type1",   &flag_type1_,    "flag_type1/I");
  tree->Branch("flag_type2",   &flag_type2_,    "flag_type2/I");
  tree->Branch("flag_type3",   &flag_type3_,    "flag_type3/I");
  tree->Branch("flag_type4",   &flag_type4_,    "flag_type4/I");
  tree->Branch("flag_type6",   &flag_type6_,    "flag_type6/I");
  tree->Branch("Vtype",        &Vtype_,         "type/I");
  tree->Branch("EVENT",        &EVENT_,         "run/I:lumi/I:event/I:json/I");
  tree->Branch("nPVs",         &nPVs_,          "nPVs/I");
  tree->Branch("PUweight",     &PUweight_,      "PUweight/F");
  tree->Branch("PUweightP",    &PUweightP_,     "PUweightP/F");
  tree->Branch("PUweightM",    &PUweightM_,     "PUweightM/F");
  tree->Branch("lheNj",        &lheNj_,         "lheNj/F");
  tree->Branch("trigger",      &trigger_,       "trigger/F");
  tree->Branch("triggerFlags", triggerFlags_,   "triggerFlags[70]/b");


  // marginalized over permutations (all <=> Sum_{p=0}^{nTotPermut}( mH=MH, mT=MT ) )
  tree->Branch(Form("p_%d_all_s",     int(MH)),   &probAtSgn_,           Form("p_%d_all_s/F",              int(MH)) );
  tree->Branch(Form("p_%d_all_b",     int(MH)),   &probAtSgn_alt_,       Form("p_%d_all_b/F",              int(MH)) );
  tree->Branch(Form("p_%d_all_s_ttbb",int(MH)),   &probAtSgn_ttbb_,      Form("p_%d_all_s_ttbb/F",         int(MH)) );
  tree->Branch(Form("p_%d_all_b_ttbb",int(MH)),   &probAtSgn_alt_ttbb_,  Form("p_%d_all_b_ttbb/F",         int(MH)) );
  tree->Branch(Form("p_%d_all_b_ttjj",int(MH)),   &probAtSgn_alt_ttjj_,  Form("p_%d_all_b_ttjj/F",         int(MH)) );
  tree->Branch(Form("p_%d_all_b_ttbj",int(MH)),   &probAtSgn_alt_ttbj_,  Form("p_%d_all_b_ttbj/F",         int(MH)) );
  tree->Branch(Form("p_%d_all_b_ttcc",int(MH)),   &probAtSgn_alt_ttcc_,  Form("p_%d_all_b_ttcc/F",         int(MH)) );

  // differential in permutations and mass value. E.g.:
  //  p_vsMH_s[0],          ..., p_vsMH_s[  nPermut-1] => prob. per-permutation and for mH = masses[0]
  //  p_vsMH_s[nPermut], ...,    p_vsMH_s[2*nPermut-1] => prob. per-permutation and for mH = masses[1]
  //  ....
  tree->Branch("p_vsMH_s",     probAtSgn_permut_,       "p_vsMH_s[nTotInteg_s]/F" );
  tree->Branch("p_vsMT_b",     probAtSgn_alt_permut_,   "p_vsMT_b[nTotInteg_b]/F");
  tree->Branch("p_vsMH_Err_s", probAtSgnErr_permut_,    "p_vsMH_Err_s[nTotInteg_s]/F" );
  tree->Branch("p_msMT_Err_b", probAtSgnErr_alt_permut_,"p_vsMT_Err_b[nTotInteg_b]/F" );

  // differential in permutation
  tree->Branch("p_tt_bb",      probAtSgn_bb_permut_,    "p_tt_bb[nPermut_s]/F");
  tree->Branch("p_tt_bj",      probAtSgn_bj_permut_,    "p_tt_bj[nPermut_b]/F");
  tree->Branch("p_tt_cc",      probAtSgn_cc_permut_,    "p_tt_cc[nPermut_b]/F");
  tree->Branch("p_tt_jj",      probAtSgn_jj_permut_,    "p_tt_jj[nPermut_b]/F");
  
  // # of mass points scanned
  tree->Branch("mH_scan",       mH_scan_,          "mH_scan[nMassPoints]/F");
  tree->Branch("mT_scan",       mT_scan_,          "mT_scan[nMassPoints]/F");
 
  // lepton kinematics
  tree->Branch("nLep",                    &nLep_,        "nLep/I");
  tree->Branch("lepton_pt",               lepton_pt_,    "lepton_pt[nLep]/F");
  tree->Branch("lepton_eta",              lepton_eta_,   "lepton_eta[nLep]/F");
  tree->Branch("lepton_phi",              lepton_phi_,   "lepton_phi[nLep]/F");
  tree->Branch("lepton_m",                lepton_m_,     "lepton_m[nLep]/F");
  tree->Branch("lepton_charge",           lepton_charge_,"lepton_charge[nLep]/F");
  tree->Branch("lepton_rIso",             lepton_rIso_,  "lepton_rIso[nLep]/F");
  tree->Branch("lepton_type",             lepton_type_,  "lepton_type[nLep]/I");  

  // MET kinematics
  tree->Branch("MET_pt",                  &MET_pt_,    "MET_pt/F");
  tree->Branch("MET_phi",                 &MET_phi_,   "MET_phi/F");
  tree->Branch("MET_sumEt",               &MET_sumEt_, "MET_sumEt/F");

  // jet kinematics
  tree->Branch("nJet",                    &nJet_,      "nJet/I");
  tree->Branch("jet_pt",                  jet_pt_,     "jet_pt[nJet]/F");
  tree->Branch("jet_eta",                 jet_eta_,    "jet_eta[nJet]/F");
  tree->Branch("jet_phi",                 jet_phi_,    "jet_phi[nJet]/F");
  tree->Branch("jet_m",                   jet_m_,      "jet_m[nJet]/F");
  tree->Branch("jet_csv",                 jet_csv_,    "jet_csv[nJet]/F");

  // Jet multiplicity
  tree->Branch("numBTagL",                &numBTagL_,    "numBTagL/I");
  tree->Branch("numBTagM",                &numBTagM_,    "numBTagM/I");
  tree->Branch("numBTagT",                &numBTagT_,    "numBTagT/I");
  tree->Branch("numJets",                 &numJets_,     "numJets/I");
  tree->Branch("btag_LR",                 &btag_LR_,     "btag_LR/F");

  // a map that associates to each permutation [p=0...nTotPermut] to the corresponding jet,
  // indexed according to the order in the jet_* collection
  //  E.g.: perm_to_jets[0] = 234567 <==> the first permutation (element '0' of p_vsMH_s/p_vsMT_b )
  //        associates element '2' of jets_* to the bLep, element '3' to W1Had, '4' to 'W2Had', '5' to bHad, and '6','7'
  //        to non-top radiation
  tree->Branch("perm_to_jet_s",             perm_to_jet_,    "perm_to_jet[nPermut_s]/I");
  tree->Branch("perm_to_gen_s" ,            perm_to_gen_,    "perm_to_gen[nPermut_s]/I");
  tree->Branch("perm_to_jet_b",             perm_to_jet_alt_,"perm_to_jet[nPermut_b]/I");
  tree->Branch("perm_to_gen_b" ,            perm_to_gen_alt_,"perm_to_gen[nPermut_b]/I");


  //tree->Branch("",                        _,  "[]/F");

  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
  /* @@@@@@@@@@@@@@@@@@@@@@@@@ OPEN FILES @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
 

  // read input files  
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

  // loop over input files
  for(unsigned int sample = 0 ; sample < mySampleFiles.size(); sample++){
    
    string currentName       = mySampleFiles[sample];

    if(currentName.find("Data")!=string::npos) isMC = false;

    mySamples->OpenFile( currentName );
    cout << "Opening file " << currentName << endl;
    TTree* currentTree       = mySamples->GetTree( currentName, "tree");
    cout << "Done!!" << endl;
    float scaleFactor        = mySamples->GetWeight(currentName);

    // variables to be used from the input files
    genParticleInfo genB, genBbar;
    genTopInfo      genTop, genTbar;
    metInfo         METtype1p2corr;
    EventInfo       EVENT;
    int             nvlep, nSimBs, /*nC, nCTop,*/ nhJets, naJets, nPVs, Vtype;
    float           PUweight, PUweightP, PUweightM;
    float           lheNj;
    float           weightTrig2012;
    UChar_t         triggerFlags[70];
    Int_t   vLepton_type      [2];
    Float_t vLepton_mass      [2];
    Float_t vLepton_pt        [2];
    Float_t vLepton_eta       [2];
    Float_t vLepton_phi       [2];
    Float_t vLepton_charge    [2];
    Float_t vLepton_pfCorrIso [2];
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
    float SimBsmass           [999];
    float SimBspt             [999];
    float SimBseta            [999];
    float SimBsphi            [999];
    //int nSvs;
    //float SvmassSv [999];
    //float Svpt     [999];
    //float Sveta    [999];
    //float Svphi    [999];
  
    currentTree->SetBranchAddress("EVENT",            &EVENT);
    currentTree->SetBranchAddress("PUweight",         &PUweight);
    currentTree->SetBranchAddress("PUweightP",        &PUweightP);
    currentTree->SetBranchAddress("PUweightM",        &PUweightM);
    currentTree->SetBranchAddress("lheNj",            &lheNj); 
    currentTree->SetBranchAddress("weightTrig2012",   &weightTrig2012); 
    currentTree->SetBranchAddress("triggerFlags",     triggerFlags); 
    currentTree->SetBranchAddress("Vtype",            &Vtype);        
    currentTree->SetBranchAddress("nhJets",           &nhJets);
    currentTree->SetBranchAddress("naJets",           &naJets);
    currentTree->SetBranchAddress("nSimBs",           &nSimBs);
    currentTree->SetBranchAddress("nvlep",            &nvlep);
    currentTree->SetBranchAddress("nPVs",             &nPVs);
    currentTree->SetBranchAddress("genB",             &genB);
    currentTree->SetBranchAddress("genBbar",          &genBbar);
    currentTree->SetBranchAddress("genTop",           &genTop);
    currentTree->SetBranchAddress("genTbar",          &genTbar);
    currentTree->SetBranchAddress("METtype1p2corr",   &METtype1p2corr);
    currentTree->SetBranchAddress("vLepton_charge",   vLepton_charge);
    currentTree->SetBranchAddress("vLepton_mass"  ,   vLepton_mass);
    currentTree->SetBranchAddress("vLepton_pt"    ,   vLepton_pt);
    currentTree->SetBranchAddress("vLepton_eta"   ,   vLepton_eta);
    currentTree->SetBranchAddress("vLepton_phi"   ,   vLepton_phi);
    currentTree->SetBranchAddress("vLepton_charge",   vLepton_charge);
    currentTree->SetBranchAddress("vLepton_pfCorrIso",vLepton_pfCorrIso);
    currentTree->SetBranchAddress("vLepton_type",     vLepton_type);
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
    currentTree->SetBranchAddress("SimBs_mass",   SimBsmass);
    currentTree->SetBranchAddress("SimBs_pt",     SimBspt);
    currentTree->SetBranchAddress("SimBs_eta",    SimBseta);
    currentTree->SetBranchAddress("SimBs_phi",    SimBsphi);
    //currentTree->SetBranchAddress("nSvs",        &nSvs);
    //currentTree->SetBranchAddress("Sv_massSv",   SvmassSv);
    //currentTree->SetBranchAddress("Sv_pt",       Svpt);
    //currentTree->SetBranchAddress("Sv_eta",      Sveta);
    //currentTree->SetBranchAddress("Sv_phi",      Svphi);

 
    /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
    /* @@@@@@@@@@@@@@@@@@@@@@@@@ EVENT LOOP @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
    
    // loop over entries
    int counter = 0;
    Long64_t nentries = currentTree->GetEntries();

    if( evHigh<0 ) evHigh = nentries;
    for (Long64_t i = 0; i < nentries ; i++){
      
      // if fixed-size job and above upper bound, continue...
      if(counter>evHigh && fixNumEvJob) continue;

      // otherwise, if outside the event window, continue...
      if(!fixNumEvJob && !(i>=evLow && i<evHigh) ) continue;
      events_++;

      if(i%5000==0) cout << i << endl;

      // read event...
      currentTree->GetEntry(i);
      
      // reset variables
      probAtSgn_          =  0.;
      probAtSgn_alt_      =  0.;
      probAtSgn_ttbb_     =  0.;
      probAtSgn_alt_ttbb_ =  0.;
      probAtSgn_alt_ttbj_ =  0.;
      probAtSgn_alt_ttcc_ =  0.;
      probAtSgn_alt_ttjj_ =  0.;

      nMassPoints_        = TMath::Max(nHiggsMassPoints,nTopMassPoints);
      flag_type0_         = -99;
      flag_type1_         = -99;
      flag_type2_         = -99;
      flag_type3_         = -99;
      flag_type4_         = -99;
      flag_type6_         = -99;

      btag_LR_            = -99;

      for( int k = 0; k < 2; k++){
	lepton_pt_[k] = -99; lepton_eta_[k] = -99; lepton_phi_[k] = -99; lepton_m_[k] = -99; lepton_charge_[k] = -99; lepton_rIso_[k] = -99; lepton_type_[k] = -99;
      }
      for( int k = 0; k < NMAXJETS; k++){
	jet_pt_[k] = -99; jet_eta_[k] = -99; jet_phi_[k] = -99; jet_m_[k] = -99;  jet_csv_[k] = -99;
      }
      for( int k = 0; k < NMAXPERMUT; k++){
	perm_to_jet_    [k] = -99;
	perm_to_gen_    [k] = -99;
	perm_to_jet_alt_[k] = -99;
	perm_to_gen_alt_[k] = -99;
      }

      // save the values into the tree (save mH[0] in case no scan is required)
      for( unsigned int m = 0; m < (unsigned int)nHiggsMassPoints ; m++){
	mH_scan_[m] = mH[m];
      }
      for( unsigned int t = 0; t < (unsigned int)nTopMassPoints ; t++){
	mT_scan_[t] = mT[t];
      }
      
      
      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ GEN PARTICLES @@@@@@@@@@@@@@@@@@@@@@@@@@  */

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
  
      // dummy cut (for the moment)
      bool properEventSL = (genBLV.Pt()>0 && genBbarLV.Pt()>0 && topBLV.Pt()>0 && topW1LV.Pt()>0 && topW2LV.Pt()>0 && atopBLV.Pt()>0 && atopW1LV.Pt()>0 && atopW2LV.Pt()>0);
      bool properEventDL = (genBLV.Pt()>0 && genBbarLV.Pt()>0 && topBLV.Pt()>0 && topW1LV.Pt()>0 && topW2LV.Pt()>0 && atopBLV.Pt()>0 && atopW1LV.Pt()>0 && atopW2LV.Pt()>0);

      if(!(properEventSL || properEventDL)){
	cout << "A dummy cut has failed..." << endl;
	continue;
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


      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ GEN BHADRONS @@@@@@@@@@@@@@@@@@@@@@@@@@@  */


      // find out number of b-hadrons in the event...
      nMatchSimBs_ = 0;
      for(int l = 0; l<nSimBs; l++){
	TLorentzVector Bs(1,0,0,1);
	Bs.SetPtEtaPhiM( SimBspt[l], SimBseta[l], SimBsphi[l], SimBsmass[l]);	    
	
	if( topBLV.Pt()>10 && deltaR(topBLV,  Bs)<0.5 ) continue;
	if(atopBLV.Pt()>10 && deltaR(atopBLV, Bs)<0.5 ) continue;
	
	for(int hj = 0; hj<nhJets; hj++){
	  TLorentzVector hJLV(1,0,0,1);
	  if(hJet_genPt[hj]>10) 
	    hJLV.SetPtEtaPhiM( hJet_genPt[hj], hJet_genEta[hj], hJet_genPhi[hj], 0.0);
	  if( hJLV.Pt()>20 && TMath::Abs(hJLV.Eta())<5 && deltaR(Bs, hJLV)<0.5 ) nMatchSimBs_++;
	}
	for(int aj = 0; aj<naJets; aj++){
	  TLorentzVector aJLV(1,0,0,1);
	  if(aJet_genPt[aj]>10) 
	    aJLV.SetPtEtaPhiM( aJet_genPt[aj], aJet_genEta[aj], aJet_genPhi[aj], 0.0);
	  if( aJLV.Pt()>20 && TMath::Abs(aJLV.Eta())<5 && deltaR(Bs, aJLV)<0.5 ) nMatchSimBs_++;
	}	
      }
      if( nSimBs>=2){
	for(int l = 0; l<nSimBs-1; l++){
	  TLorentzVector Bs1(1,0,0,1);
	  Bs1.SetPtEtaPhiM( SimBspt[l], SimBseta[l], SimBsphi[l], SimBsmass[l]);	    
	  if( topBLV.Pt()>10 && deltaR(topBLV,  Bs1)<0.5 ) continue;
	  if(atopBLV.Pt()>10 && deltaR(atopBLV, Bs1)<0.5 ) continue;
	  for(int m = l+1; m<nSimBs; m++){
	    TLorentzVector Bs2(1,0,0,1);
		Bs2.SetPtEtaPhiM( SimBspt[m], SimBseta[m], SimBsphi[m], SimBsmass[m]);	    	    
		if( topBLV.Pt()>10 && deltaR(topBLV,  Bs2)<0.5 ) continue;
		if(atopBLV.Pt()>10 && deltaR(atopBLV, Bs2)<0.5 ) continue;
		if( deltaR(Bs1,Bs2)<0.50 ) nMatchSimBs_--;
	  }
	}
      }

      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@ LEPTON SELECTION @@@@@@@@@@@@@@@@@@@@@@@@@  */
      

      // charged leptons and MET
      TLorentzVector leptonLV, leptonLV2;
      TLorentzVector neutrinoLV;


      ///////////////////////////////////
      //         SL events             //
      ///////////////////////////////////

      properEventSL = false;      
      if( nvlep==1 && (Vtype==2 || Vtype==3) ){

	// first lepton...
	leptonLV.SetPtEtaPhiM(vLepton_pt[0],vLepton_eta[0],vLepton_phi[0],vLepton_mass[0]);      

	// cut on lepton (SL)
	int lepSelVtype2 =  (Vtype==2 && vLepton_type[0]==13 && leptonLV.Pt()>30 && TMath::Abs(leptonLV.Eta())<2.1 && vLepton_pfCorrIso[0]<0.10);
	int lepSelVtype3 =  (Vtype==3 && vLepton_type[0]==11 && leptonLV.Pt()>30 && TMath::Abs(leptonLV.Eta())<2.5 && !(TMath::Abs(leptonLV.Eta())>1.442 &&  TMath::Abs(leptonLV.Eta())<1.566) && vLepton_pfCorrIso[0]<0.10) ;

	// OR of four trigger paths:  "HLT_Mu40_eta2p1_v.*", "HLT_IsoMu24_eta2p1_v.*", "HLT_Mu40_v.*",  "HLT_IsoMu24_v.*"
	int trigVtype2 =  (Vtype==2 && ( triggerFlags[22]>0 || triggerFlags[23]>0 || triggerFlags[14]>0 ||triggerFlags[21]>0 ));

	// OR of one trigger paths:   "HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v.*"
	int trigVtype3 =  (Vtype==3 && ( triggerFlags[3]>0 ) );
  
	// ID && trigger
	properEventSL = (lepSelVtype2 && (isMC ? 1 : trigVtype2)) || (lepSelVtype3 && (isMC ? 1 : trigVtype3));	 


	// save lepton kinematics...
	nLep_ = 1;

	lepton_pt_     [0] = leptonLV.Pt();
	lepton_eta_    [0] = leptonLV.Eta();
	lepton_phi_    [0] = leptonLV.Phi();
	lepton_m_      [0] = leptonLV.M();
	lepton_charge_ [0] = vLepton_charge   [0];
	lepton_rIso_   [0] = vLepton_pfCorrIso[0];
	lepton_type_   [0] = vLepton_type[0];

      }


      ///////////////////////////////////
      //         DL events             //
      ///////////////////////////////////

      properEventDL = false;
      if( nvlep>=2 && (Vtype==0 || Vtype==1)){
	
	// first lepton...
	leptonLV.SetPtEtaPhiM (vLepton_pt[0],vLepton_eta[0],vLepton_phi[0],vLepton_mass[0]);
	// second lepton...
	leptonLV2.SetPtEtaPhiM(vLepton_pt[1],vLepton_eta[1],vLepton_phi[1],vLepton_mass[1]);

	// cut on leptons (DL)
	int lepSelVtype0 = ( Vtype==0 && vLepton_type[0]==13 && vLepton_type[1]==13 && 
			    ( (leptonLV.Pt() >20 && TMath::Abs(leptonLV.Eta()) <2.1 && vLepton_pfCorrIso[0]<0.10 &&
			       leptonLV2.Pt()>10 && TMath::Abs(leptonLV2.Eta())<2.4 && vLepton_pfCorrIso[1]<0.20) ||
			      (leptonLV.Pt() >10 && TMath::Abs(leptonLV.Eta()) <2.4 && vLepton_pfCorrIso[0]<0.20 &&
			       leptonLV2.Pt()>20 && TMath::Abs(leptonLV2.Eta())<2.1 && vLepton_pfCorrIso[1]<0.10) )
			    ) && vLepton_charge[0]*vLepton_charge[1]<0;

	int lepSelVtype1 = ( Vtype==1 && vLepton_type[0]==11 && vLepton_type[0]==11 && 
			    ( (leptonLV.Pt() >20 && TMath::Abs(leptonLV.Eta()) <2.5 && vLepton_pfCorrIso[0]<0.10 &&
			       leptonLV2.Pt()>10 && TMath::Abs(leptonLV2.Eta())<2.5 && vLepton_pfCorrIso[1]<0.20) ||
			      (leptonLV.Pt() >10 && TMath::Abs(leptonLV.Eta()) <2.5 && vLepton_pfCorrIso[0]<0.20 &&
			       leptonLV2.Pt()>20 && TMath::Abs(leptonLV2.Eta())<2.5 && vLepton_pfCorrIso[1]<0.10) )
			    ) && vLepton_charge[0]*vLepton_charge[1]<0;
	
	// OR of four trigger paths:  "HLT_Mu40_eta2p1_v.*", "HLT_IsoMu24_eta2p1_v.*", "HLT_Mu40_v.*",  "HLT_IsoMu24_v.*"
	int trigVtype0 =  (Vtype==0 && ( triggerFlags[22]>0 || triggerFlags[23]>0 || triggerFlags[14]>0 ||triggerFlags[21]>0 ));

	// OR of two trigger paths:    "HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v.*", "HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v.*"     
	int trigVtype1 =  (Vtype==1 && ( triggerFlags[5]>0 || triggerFlags[6]>0 ) );		 
			
	// ID && trigger
	properEventDL = (lepSelVtype0 && (isMC ? 1 : trigVtype0)) || (lepSelVtype1 && (isMC ? 1 : trigVtype1));
	  

	// save lepton(s) kinematics into the tree...
	nLep_ = 2;

	// lep 1...
	lepton_pt_     [0] = leptonLV.Pt();
	lepton_eta_    [0] = leptonLV.Eta();
	lepton_phi_    [0] = leptonLV.Phi();
	lepton_m_      [0] = leptonLV.M();
	lepton_charge_ [0] = vLepton_charge   [0];
	lepton_rIso_   [0] = vLepton_pfCorrIso[0];
	lepton_type_   [0] = vLepton_type[0];	

	// lep 2...
	lepton_pt_     [1] = leptonLV2.Pt();
	lepton_eta_    [1] = leptonLV2.Eta();
	lepton_phi_    [1] = leptonLV2.Phi();
	lepton_m_      [1] = leptonLV2.M();
	lepton_charge_ [1] = vLepton_charge   [1];
	lepton_rIso_   [1] = vLepton_pfCorrIso[1];
	lepton_type_   [1] = vLepton_type[1];

      }

      ////////////////////////////////////////////////////////////////////////

      // MET
      float nuPx = METtype1p2corr.et*TMath::Cos(METtype1p2corr.phi);
      float nuPy = METtype1p2corr.et*TMath::Sin(METtype1p2corr.phi);
      float nuE  = TMath::Sqrt(nuPx*nuPx+nuPy*nuPy);
      neutrinoLV.SetPxPyPzE(nuPx,nuPy,0. ,nuE);

      // save MET kinematics into the tree...
      MET_pt_    = neutrinoLV.Pt();
      MET_phi_   = neutrinoLV.Phi();
      MET_sumEt_ = METtype1p2corr.sumet; 

      // continue if leptons do not satisfy cuts
      if( !(properEventSL || properEventDL) ) continue;


      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
      /* @@@@@@@@@@@@@@@@@@@@@@@@ JET SELECTION @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */

      // this container will hold the jets
      std::vector<JetObservable> jet_map;

      // loop over jet collections
      for(int coll = 0 ; coll < 2 ; coll++){

	// loop over jets
	for(int hj = 0; hj < (coll==0 ? nhJets : naJets); hj++){

	  float ptGen = -99.;
	  if(coll==0 && hJet_genPt[hj]>0.) ptGen = hJet_genPt[hj];
	  if(coll==1 && aJet_genPt[hj]>0.) ptGen = aJet_genPt[hj];

	  float pt     = (coll==0) ? hJet_pt [hj]  : aJet_pt [hj];
	  float eta    = (coll==0) ? hJet_eta[hj]  : aJet_eta[hj];
	  float phi    = (coll==0) ? hJet_phi[hj]  : aJet_phi[hj];
	  float e      = (coll==0) ? hJet_e  [hj]  : aJet_e  [hj];
	  float m2     = e*e - pt*pt*TMath::CosH(eta)*TMath::CosH(eta);
	  if(m2<0) m2 = 0.; 
	  float m      = TMath::Sqrt( m2 ); 

	  int id       = (coll==0) ? hJet_puJetIdL[hj] : aJet_puJetIdL[hj];
	  float JECUnc = (coll==0) ? hJet_JECUnc  [hj] : aJet_JECUnc  [hj];

	  // for JEC/JER systematics (N.B. this assumes that the jets are already corrected)
	  float shift     = 1.0;
	  if     ( doJECup   )  shift *= (1+JECUnc);
	  else if( doJECdown )  shift *= (1-JECUnc);
	  else if( doJERup   )  shift *= (ptGen>0. ?  1+resolutionBias(TMath::Abs(eta), +1, JERCORRECTED)*(1-ptGen/pt)   : 1.0);
	  else if( doJERdown )  shift *= (ptGen>0. ?  1+resolutionBias(TMath::Abs(eta), -1, JERCORRECTED)*(1-ptGen/pt)   : 1.0);
	  else{}

	  // if correct for bias in jet resolution (for sanity, enforce isMC)
	  if( isMC && doJERbias && !doJERup && !doJERdown) 
	    shift *= (ptGen>0. ?  1+resolutionBias(TMath::Abs(eta), 0, JERCORRECTED)*(1-ptGen/pt)   : 1.0);

	  // change energy/mass by shift
	  pt *= shift;
	  m  *= shift;

	  // only jets in acceptance...
	  if( TMath::Abs(eta)> 2.5 ) continue;

	  // only jets passing ID...
	  if( id < 0.5 ) continue;	

	  // only jets above pt cut...
	  if( pt < 30  ) continue;	  

	  TLorentzVector p4;
	  p4.SetPtEtaPhiM( pt, eta, phi, m );

	  // for csv systematics
	  float csv_nominal =  (coll==0) ? hJet_csv_nominal[hj] : aJet_csv_nominal[hj];
	  float csv_upBC    =  (coll==0) ? hJet_csv_upBC   [hj] : aJet_csv_upBC   [hj];
	  float csv_downBC  =  (coll==0) ? hJet_csv_downBC [hj] : aJet_csv_downBC [hj];
	  float csv_upL     =  (coll==0) ? hJet_csv_upL    [hj] : aJet_csv_upL    [hj];
	  float csv_downL   =  (coll==0) ? hJet_csv_downL  [hj] : aJet_csv_downL  [hj];
	  float csv = csv_nominal;
	  if     ( doCSVup  ) csv =  TMath::Max(csv_upBC,   csv_upL);
	  else if( doCSVdown) csv =  TMath::Min(csv_downBC, csv_downL);
	  else{}

	  // the b-tagger output 
	  //  ==> Min needed because in csvUp, csv can exceed 1..., 
	  //      Max needed because we crunch  [-inf,0[ -> {0.}
	  csv         =  TMath::Min( TMath::Max( csv, float(0.)), float(0.999999) );

	  // the jet observables (p4 and csv)
	  JetObservable myJet;
	  myJet.p4  = p4;
	  myJet.csv = csv; 
	  
	  // use pt to order jet collection
	  //jet_map[ p4.Pt() ] = myJet;
	  jet_map.push_back( myJet );

	  //if( debug )
	  //cout << "Jet #" << coll << "-" << hj << " => (" << pt << "," << eta << "," << phi << "," << m << "), ID=" << id << endl;
	  
	}
      }


      // order jet list by Pt
      std::sort( jet_map.begin(), jet_map.end(), JetObservableListerByPt() );
      
      // if use btag shape, order by decreasing CSV 
      //       <=> when considering only a subset of the jets, this ensures that the combination
      //           obtained from CSVM only jets is among those considered
      if( selectByBTagShape ) 
	std::sort( jet_map.begin(), jet_map.end(), JetObservableListerByCSV() );

   
      // fill arrays of jets
      std::vector<TLorentzVector>  jets_p4;
      std::vector<double>          jets_csv;
      std::vector<double>          jets_csv_prob_b;
      std::vector<double>          jets_csv_prob_c;
      std::vector<double>          jets_csv_prob_j;
      int jetsAboveCut = 0;

      //for( jet_map_it = jet_map.begin() ; jet_map_it != jet_map.end(); jet_map_it++){
      for(unsigned int jj = 0; jj < jet_map.size() ; jj++ ){
	
	// the four-vector
	TLorentzVector p4 = jet_map[jj].p4; //(jet_map_it->second).p4;
	
	// the csv value
	float csv = jet_map[jj].csv;

	// count jets above 40 GeV
	if( p4.Pt()>40 ) jetsAboveCut++;
	
	// store jet p4...
	jets_p4.push_back ( p4  );
	// store csv
	jets_csv.push_back( csv );
	
	// needed to find appropriate csv PDF
	string bin = "";
	if( TMath::Abs( p4.Eta() ) <= 1.0 ) 
	  bin = "Bin0";
	if( TMath::Abs( p4.Eta() ) >  1.0 ) 
	  bin = "Bin1";

	// store PDF(csv)
	jets_csv_prob_b.push_back( btagger["b_"+bin]!=0 ? btagger["b_"+bin]->GetBinContent( btagger["b_"+bin]->FindBin( csv ) ) : 1.);
	jets_csv_prob_c.push_back( btagger["c_"+bin]!=0 ? btagger["c_"+bin]->GetBinContent( btagger["c_"+bin]->FindBin( csv ) ) : 1.);
	jets_csv_prob_j.push_back( btagger["l_"+bin]!=0 ? btagger["l_"+bin]->GetBinContent( btagger["l_"+bin]->FindBin( csv ) ) : 1.);

      }

      // continue if not enough jets
      if( jetsAboveCut<4 ){
	//cout << "Less then 4 jets.. continue" << endl;
	continue;
      }
      

      /* @@@@@@@@@@@@@@@@@@@@@@@@ JET ORDERING @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */

      // jet multiplicity
      int numJets30UntagM = 0; 
      int numJets30UntagL = 0; 
      int numJets30BtagL  = 0; 
      int numJets30BtagM  = 0; 
      int numJets30BtagT  = 0;

      // this vector contains the indices of up to 6 jets passing the CSVM(L) b-tag selection... 
      // (type 7 uses two working points)
      vector<unsigned int> btag_indices;
      vector<unsigned int> btagLoose_indices;

      // this vector contains the indices of up to 6 jets failing the CSVM(L) b-tag selection... 
      vector<unsigned int> buntag_indices;
      vector<unsigned int> buntagLoose_indices;

      // this vector contains the indices of up to 12 jets passing ANY b-tag selection... 
      vector<unsigned int> banytag_indices;  
    

      for(unsigned int k = 0; k < jets_p4.size(); k++){  

	float csv_k = jets_csv[k];
	float pt_k  = (jets_p4[k]).Pt();

	// passes CSVL...
	int btag_L = csv_k>0.244 ;	
	// passes CSVM...
	int btag_M = csv_k>0.679 ;
	// passes CSVT...
	int btag_T = csv_k>0.898 ;	

	// any btag value...
	if( pt_k>30 ) banytag_indices.push_back( k );

	// passing at least CSVL...
	if( pt_k>30 &&  btag_L ){
	  numJets30BtagL ++;
	  btagLoose_indices.push_back( k );
	}
	// passing at least CSVM...
	if( pt_k>30 &&  btag_M ){
	  numJets30BtagM ++;
	  btag_indices.push_back( k );
	}
	// passing at least CSVT...
	if( pt_k>30 &&  btag_T ){
	  numJets30BtagT ++;
	}

	// failing at most CSVL...
	if( pt_k>30 && !btag_L ){
	  numJets30UntagL++;
	  buntagLoose_indices.push_back( k );
	}
	// failing at most CSVM...
	if( pt_k>30 && !btag_M ){
	  numJets30UntagM++;
	  buntag_indices.push_back( k );
	}

      }

      numBTagL_ = numJets30BtagL;
      numBTagM_ = numJets30BtagM;
      numBTagT_ = numJets30BtagT;
      numJets_  = ( numJets30BtagM + numJets30UntagM);


      /* @@@@@@@@@@@@@@@@@@@@@@@@ JET RE-ORDERING BY CSV PROB @@@@@@@@@@@@@@@@@@@@@@@@@@  */      

      // in case the jet collections is reshuffled using the b-probability, keep a copy of old selection
      // (for debugging purposes)
      vector<unsigned int> btag_indices_backup   = btag_indices;
      vector<unsigned int> buntag_indices_backup = buntag_indices;

      // variable that flags events passing or failing the b-probability cut
      int passes_btagshape = 0;

      if( selectByBTagShape && useBtag &&                    // run this only if specified AND...
	  ((properEventSL && banytag_indices.size()>=5) ||   // [ run this only if SL and at least 5 jets OR...
	   (properEventDL && banytag_indices.size()>=4)) ){  //   run this only if DL and at least 4 jets ]
	
	// map that transforms the indices into the permutation convention
	std::map< unsigned int, unsigned int> btag_map;
	btag_map.clear();
	
	// number of inequivalent permutations
	int nS, nB;
	
	// sum of the b-tag probability over all permutations
	float p_bb     =  0.;
	float p_jj     =  0.;

	// permutation with the largest probability 
	float max_p_bb = -1;
	unsigned int selected_comb = 999;
	
	// a flag to keep track of the jet multiplicity
	int btag_flag = -999;
	
	// if numJets>6, the first six jets ordered by decreasing pt, and
	// then by decreasing csv are considered...
	if( properEventSL && banytag_indices.size()>=6 ){
	  btag_map[2] = banytag_indices[0]; 
	  btag_map[3] = banytag_indices[1]; 
	  btag_map[4] = banytag_indices[2];
	  btag_map[5] = banytag_indices[3]; 
	  btag_map[6] = banytag_indices[4]; 
	  btag_map[7] = banytag_indices[5];
	  nS = 15; 
	  nB = 15;
	  btag_flag = 0;
	}
	else if( properEventSL && banytag_indices.size()==5 ){
	  btag_map[2] = banytag_indices[0]; 
	  btag_map[3] = banytag_indices[1]; 
	  btag_map[4] = banytag_indices[1];
	  btag_map[5] = banytag_indices[2]; 
	  btag_map[6] = banytag_indices[3]; 
	  btag_map[7] = banytag_indices[4];
	  nS =  5; 
	  nB = 10;
	  btag_flag = 1;
	}
	// if numJets>4, the first four jets ordered by decreasing pt, and
	// then by decreasing csv are considered...
	else if( properEventDL && banytag_indices.size()>=4 ){
	  btag_map[2] = banytag_indices[0]; 
	  btag_map[3] = banytag_indices[0]; 
	  btag_map[4] = banytag_indices[0];
	  btag_map[5] = banytag_indices[1]; 
	  btag_map[6] = banytag_indices[2]; 
	  btag_map[7] = banytag_indices[3];
	  nS = 1; 
	  nB = 6;
	  btag_flag = 2;
	}
	else{ 
	  cout << "Inconsistency in selectByBTagShape... continue" << endl;
	  continue;
	}		

	// loop over hypothesis [ TTH, TTbb ]
	for(int hyp = 0 ; hyp<2;  hyp++){
	  
	  // list of permutations
	  int* permutList = 0;	   
	  if( btag_flag == 0 ) permutList = hyp==0 ?  permutations_6J_S : permutations_6J_B;	    
	  if( btag_flag == 1 ) permutList = hyp==0 ?  permutations_5J_S : permutations_5J_B;
	  if( btag_flag == 2 ) permutList = hyp==0 ?  permutations_4J_S : permutations_4J_B;
	  
	  // loop over permutations
	  for(unsigned int pos = 0; pos < (unsigned int)( hyp==0 ? nS : nB ) ; pos++){
	    
	    // index of the four jets associated to b-quarks or W->qq
	    int bLep_pos = (permutList[pos])%1000000/100000;
	    int w1_pos   = (permutList[pos])%100000/10000;
	    int w2_pos   = (permutList[pos])%10000/1000;
	    int bHad_pos = (permutList[pos])%1000/100;
	    int b1_pos   = (permutList[pos])%100/10;
	    int b2_pos   = (permutList[pos])%10/1; 
	    
	    double p_b_bLep =  jets_csv_prob_b[  btag_map[bLep_pos] ];
	    double p_b_bHad =  jets_csv_prob_b[  btag_map[bHad_pos] ];
	    double p_b_b1   =  jets_csv_prob_b[  btag_map[b1_pos]   ];
	    double p_j_b1   =  jets_csv_prob_j[  btag_map[b1_pos]   ];
	    double p_b_b2   =  jets_csv_prob_b[  btag_map[b2_pos]   ];
	    double p_j_b2   =  jets_csv_prob_j[  btag_map[b2_pos]   ];
	    double p_j_w1   =  jets_csv_prob_j[  btag_map[w1_pos]   ];
	    double p_j_w2   =  jets_csv_prob_j[  btag_map[w2_pos]   ];
	    
	    // the total probability
	    float p_pos = 0.;

	    // if sgn (ttbb)
	    if( hyp==0 ){
	      p_pos =  p_b_bLep * p_b_bHad * p_b_b1 * p_b_b2 * p_j_w1 * p_j_w2; 
	      p_bb += p_pos;

	      // look for a maximum
	      if(  p_pos > max_p_bb ){
		max_p_bb      = p_pos;
		selected_comb = pos;
	      }
	    }

	    // if bkg (ttjj)
	    if( hyp==1 ){
	      p_pos =  p_b_bLep * p_b_bHad * p_j_b1 * p_j_b2 * p_j_w1 * p_j_w2;
	      p_jj += p_pos;
	    }
	    
	  }
	  
	} // end loop over hypothesis

	// normalize the probabilities...
	p_bb /= nS;
	p_jj /= nB;

	// LR of ttbb vs ttjj hypotheses as variable to select events
	btag_LR_ = (p_bb+p_jj)>0 ? p_bb/(p_bb+p_jj) : 0. ;

	// depending on event type, check if the event passes the cut:
	// if it does, check which combination yields the **largest** ttbb probability
	int* permutListS = 0;	   	 
	switch( btag_flag ){
	case 0:
	  passes_btagshape = ( btag_LR_ >= btag_prob_cut_6jets && selected_comb!=999);
	  permutListS      = permutations_6J_S;
	  break;
	case 1:
	  passes_btagshape = ( btag_LR_ >= btag_prob_cut_5jets && selected_comb!=999);
	  permutListS      = permutations_5J_S;
	  break;
	case 2:
	  passes_btagshape = ( btag_LR_ >= btag_prob_cut_4jets && selected_comb!=999);
	  permutListS      = permutations_4J_S;
	  break;
	default:
	  break;
	}

	// if the event passes the cut, reshuffle jets into btag/buntag vectors
	// N.B. ==> jets will be sorted by descending btag probaility!!!
	if( passes_btagshape ){

	  // reset old vector
	  btag_indices.clear();

	  // these are the four jets that are most compatible with four b-quarks
	  btag_indices.push_back(  btag_map[(permutListS[selected_comb])%1000000/100000] );
	  btag_indices.push_back(  btag_map[(permutListS[selected_comb])%1000/100]       );
	  btag_indices.push_back(  btag_map[(permutListS[selected_comb])%100/10]         );
	  btag_indices.push_back(  btag_map[(permutListS[selected_comb])%10/1]           );
	  
	  // all other jets go into this collection
	  buntag_indices.clear();
	  for( unsigned int jj = 0 ; jj<banytag_indices.size(); jj++){

	    // check for a match among the b-tagged jets
	    int isTagged = 0;
	    for( unsigned int bb = 0 ; bb<btag_indices.size(); bb++ ) 
	      if( banytag_indices[jj]==btag_indices[bb] ) isTagged = 1;

	    // if not matches, push back
	    if( isTagged == 0) 
	      buntag_indices.push_back( banytag_indices[jj] );	    
	  }

	}
	
      } 
     


      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@ EVENT SELECTION @@@@@@@@@@@@@@@@@@@@@@@@@  */


      ////////////////////////////////////////////////////
      // SEMILEPTONIC EVENTS                            //
      //        FULLLEPTONIC EVENTS                     //
      ////////////////////////////////////////////////////

      // categories defined by jet and btagged jet multiplicity (Nj,Nb)
      bool analyze_type0       = properEventSL && numJets30BtagM==4 && numJets30UntagM==2  && doType0;
      bool analyze_type1       = properEventSL && numJets30BtagM==4 && numJets30UntagM==2  && doType1;
      bool analyze_type2       = properEventSL && numJets30BtagM==4 && numJets30UntagM==1  && doType2;
      bool analyze_type3       = properEventSL && numJets30BtagM==4 && numJets30UntagM >2  && doType3;
      bool analyze_type6       = properEventDL && numJets30BtagM==4                        && doType6;
      bool analyze_type7       = properEventDL && numJets30BtagM==3 && numJets30BtagL==4   && doType7;

      // categories defined by jet multiplicity (Nj)
      bool analyze_typeBTag6   = properEventSL && (numJets_==6)                            && doTypeBTag6;
      bool analyze_typeBTag5   = properEventSL && (numJets_==5)                            && doTypeBTag5;
      bool analyze_typeBTag4   = properEventDL && (numJets_==4)                            && doTypeBTag4;

      // categories defined by jet multiplicity (Nj), but passing a minimum cut on the btag likelihood
      bool analyze_type0_BTag  = properEventSL && (numJets_==6)                            && doType0ByBTagShape && passes_btagshape;
      bool analyze_type1_BTag  = properEventSL && (numJets_==6)                            && doType1ByBTagShape && passes_btagshape;
      bool analyze_type2_BTag  = properEventSL && (numJets_==5)                            && doType2ByBTagShape && passes_btagshape;
      bool analyze_type3_BTag  = properEventSL && (numJets_ >6)                            && doType3ByBTagShape && passes_btagshape;
      bool analyze_type6_BTag  = properEventDL && (numJets_>=4)                            && doType6ByBTagShape && passes_btagshape;


      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@ GEN MATCH @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */

      // temporary variables to count num. of matches
      int hMatches = 0;
      int tMatches = 0;
      int wMatches = 0;

      for( unsigned int w = 0; w<btag_indices.size(); w++){

	// use the Py() component to assess if the gen particle is present in the original tree
	// (it is 0.0 otherwise)
	if     (  TMath::Abs(HIGGSB1.Py()) >0  && deltaR(jets_p4[ btag_indices[w] ],HIGGSB1)  < GENJETDR ) hMatches++;
	else if(  TMath::Abs(HIGGSB2.Py()) >0  && deltaR(jets_p4[ btag_indices[w] ],HIGGSB2)  < GENJETDR ) hMatches++;    
	if     (  TMath::Abs(TOPHADB.Py()) >0  && deltaR(jets_p4[ btag_indices[w] ],TOPHADB)  < GENJETDR ) tMatches++;
	else if(  TMath::Abs(TOPLEPB.Py()) >0  && deltaR(jets_p4[ btag_indices[w] ],TOPLEPB)  < GENJETDR ) tMatches++;
	if     (  TMath::Abs(TOPHADW1.Py())>0  && deltaR(jets_p4[ btag_indices[w] ],TOPHADW1) < GENJETDR ) wMatches++;
	else if(  TMath::Abs(TOPHADW2.Py())>0  && deltaR(jets_p4[ btag_indices[w] ],TOPHADW2) < GENJETDR ) wMatches++;    
      }
      matchesH_    = hMatches;
      matchesHAll_ = hMatches;
      matchesT_    = tMatches;
      matchesTAll_ = tMatches;
      matchesWAll_ = wMatches;  
      wMatches = 0;

      for( unsigned int w = 0; w<buntag_indices.size(); w++){

	// use the Py() component to assess if the gen particle is present in the original tree
	// (it is 0.0 otherwise)	
	if     (   TMath::Abs(HIGGSB1.Py())>0    && deltaR(jets_p4[ buntag_indices[w] ],HIGGSB1)  < GENJETDR ) matchesHAll_++;
	else if(   TMath::Abs(HIGGSB2.Py())>0    && deltaR(jets_p4[ buntag_indices[w] ],HIGGSB2)  < GENJETDR ) matchesHAll_++;
	if     (   TMath::Abs(TOPHADB.Py())>0    && deltaR(jets_p4[ buntag_indices[w] ],TOPHADB)  < GENJETDR ) matchesTAll_++;
	else if(   TMath::Abs(TOPLEPB.Py())>0    && deltaR(jets_p4[ buntag_indices[w] ],TOPLEPB)  < GENJETDR ) matchesTAll_++;
	if     (   TMath::Abs(TOPHADW1.Py())>0   && deltaR(jets_p4[ buntag_indices[w] ],TOPHADW1) < GENJETDR ) wMatches++;
	else if(   TMath::Abs(TOPHADW2.Py())>0   && deltaR(jets_p4[ buntag_indices[w] ],TOPHADW2) < GENJETDR ) wMatches++;    
      }
      matchesW_    =  wMatches; 
      matchesWAll_ += wMatches;

      // gen level overlap
      vector<TLorentzVector> genHeavy;
      if( TMath::Abs(HIGGSB1.Py())>0) genHeavy.push_back( HIGGSB1);
      if( TMath::Abs(HIGGSB2.Py())>0) genHeavy.push_back( HIGGSB2);
      if( TMath::Abs(TOPHADB.Py())>0) genHeavy.push_back( TOPHADB);
      if( TMath::Abs(TOPLEPB.Py())>0) genHeavy.push_back( TOPLEPB);
      int overlapH = 0;
      if(genHeavy.size()>1){
	for(unsigned int k = 0; k < genHeavy.size()-1; k++){  
	  for(unsigned int l = k+1; l < genHeavy.size(); l++){  
	    if( deltaR(genHeavy[k], genHeavy[l])<0.5 ) overlapH++;
	  }	    
	}
      }
      overlapHeavy_ = overlapH;
      vector<TLorentzVector> genLight;
      if(  TMath::Abs(TOPHADW1.Py())>0) genLight.push_back( TOPHADW1);
      if(  TMath::Abs(TOPHADW2.Py())>0) genLight.push_back( TOPHADW2);
      int overlapL = 0;
      for(unsigned int k = 0; k < genLight.size(); k++){  
	for(unsigned int l = 0; l < genHeavy.size(); l++){  
	  if( deltaR(genLight[k], genHeavy[l])<0.5 ) overlapL++;
	}	    
      }
      if( genLight.size()>1 && deltaR(genLight[0], genLight[1])<0.5  ) overlapL++;
      overlapLight_  = overlapL;      
      


      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ANALYSIS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */


      //  input 4-vectors
      vector<TLorentzVector> jets;

      // internal map: [ position in "jets" ] -> [ position in "jets_p4" ]
      std::map< unsigned int, unsigned int> pos_to_index;

      // consider th event only if of the desired type
      if( analyze_typeBTag6  || analyze_typeBTag5  || analyze_typeBTag4  || 
	  analyze_type0      || analyze_type1      || analyze_type2      || analyze_type3      || analyze_type6       || analyze_type7 ||
	  analyze_type0_BTag || analyze_type1_BTag || analyze_type2_BTag || analyze_type3_BTag || analyze_type6_BTag){	
	
	// find out which two untagged jets come from W->qq'
	unsigned int ind1 = 999;
	unsigned int ind2 = 999;
	
	if( analyze_typeBTag4 ){
	  
	  counter++;      
	  if(fixNumEvJob && !(counter>=evLow && counter<=evHigh) ) continue;
	  if(print) cout << "Processing event # " << counter << " (TYPE DL-4 JETS)." << " Event ID " << EVENT.event << endl;	   

	  /////////////////////////////////////////////////////
	  type_       = -3;
	  nPermut_    =  1;
	  nPermut_alt_=  6;
	  meIntegrator->setIntType( MEIntegratorNew::DL );
	  /////////////////////////////////////////////////////	  

	}
	else if( analyze_typeBTag5 ){
	  
	  counter++;      
	  if(fixNumEvJob && !(counter>=evLow && counter<=evHigh) ) continue;
	  if(print) cout << "Processing event # " << counter << " (TYPE SL-5 JETS)." << " Event ID " << EVENT.event << endl;	   

	  /////////////////////////////////////////////////////
	  type_       = -1;
	  nPermut_    =  5;
	  nPermut_alt_= 10;
	  meIntegrator->setIntType( MEIntegratorNew::SL1wj );
	  /////////////////////////////////////////////////////	  

	}
	else if( analyze_typeBTag6 ){
	  
	  counter++;      
	  if(fixNumEvJob && !(counter>=evLow && counter<=evHigh) ) continue;
	  if(print) cout << "Processing event # " << counter << " (TYPE SL6 JETS)." << " Event ID " << EVENT.event << endl;	   

	  /////////////////////////////////////////////////////
	  type_       =  -2;
	  nPermut_    =  15;
	  nPermut_alt_=  15;	  
	  meIntegrator->setIntType( MEIntegratorNew::SL2wj );
	  /////////////////////////////////////////////////////	  

	}
	else if( (analyze_type0 || analyze_type1 || analyze_type0_BTag || analyze_type1_BTag) ){
	  
	  // sanity check:
	  if(buntag_indices.size() != 2){
	    cout << "Inconsistency found for (analyze_type0 || analyze_type1)... continue" << endl;
	    continue;
	  }

	  // use untagged mass to assign to type 0 OR type 1
	  float WMass = (jets_p4[ buntag_indices[0] ]+jets_p4[  buntag_indices[1] ]).M();
	  
	  // set index for untagged jets
	  ind1 = buntag_indices[0];
	  ind2 = buntag_indices[1];
	  
	  if( (WMass>MwL && WMass<MwH)  && (analyze_type0 || analyze_type0_BTag) ){
	    
	    counter++;      
	    if( fixNumEvJob && !(counter>=evLow && counter<=evHigh) ) continue;
	    if(print) cout << "Processing event # " << counter << " (TYPE 0), mW=" << WMass << " GeV." << " Event ID " << EVENT.event << endl;
	    
	    /////////////////////////////////////////////////////	      
	    type_       =  0;
	    nPermut_    = 12;
	    nPermut_alt_= 12;
	    meIntegrator->setIntType( MEIntegratorNew::SL2wj );	    
	    /////////////////////////////////////////////////////
	  }
	  else if( !( WMass>MwL && WMass<MwH) && (analyze_type1 || analyze_type1_BTag)){
	    
	    counter++;      
	    if( fixNumEvJob && !(counter>=evLow && counter<=evHigh)  ) continue;
	    if(print) cout << "Processing event # " << counter << " (TYPE 1), mW=" << WMass << " GeV." << " Event ID " << EVENT.event << endl;	   
	    
	    /////////////////////////////////////////////////////	  
	    type_       =  1;
	    nPermut_    = 24;
	    nPermut_alt_= 24;
	    meIntegrator->setIntType( MEIntegratorNew::SL1wj );
	    /////////////////////////////////////////////////////
	  }
	  else{ continue; }

	}
	else if( analyze_type2 || analyze_type2_BTag ){	   
	  
	  // sanity check:
	  if(buntag_indices.size() != 1){
	    cout << "Inconsistency found for analyze_type2... continue" << endl;
	    continue;
	  }

	  counter++;
	  if(fixNumEvJob && !(counter>=evLow && counter<=evHigh) ) continue;
	  if(print) cout << "Processing event # " << counter << " (TYPE 1)." << " Event ID " << EVENT.event << endl;
	  
	  // set index for untagged jet
	  ind1 = buntag_indices[0];
	  ind2 = buntag_indices[0];
	  
	  /////////////////////////////////////////////////////
	  type_       =  2;
	  nPermut_    = 12;
	  nPermut_alt_= 12;
	  meIntegrator->setIntType( MEIntegratorNew::SL1wj );
	  /////////////////////////////////////////////////////
	}
	else if( analyze_type3 || analyze_type3_BTag ){
	  
	  // sanity check:
	  if(buntag_indices.size() <= 2 ){
	    cout << "Inconsistency found for analyze_type3... continue" << endl;
	    continue;
	  }

	  counter++;      
	  if(fixNumEvJob && !(counter>=evLow && counter<=evHigh) ) continue;
	  if(print) cout << "Processing event # " << counter << " (TYPE 3)." << " Event ID " << EVENT.event << endl;	   
	  
	  // find out which are ind1 and ind2...
	  float minDiff     = 99999.;
	  for(unsigned int uj1 = 0; uj1<buntag_indices.size()-1; uj1++){
	    for(unsigned int uj2 = uj1+1; uj2<buntag_indices.size(); uj2++){
	      
	      float WMass12 = (jets_p4[ buntag_indices[uj1] ]+jets_p4[ buntag_indices[uj2] ]).M();
	      if( TMath::Abs(WMass12-80.1)<minDiff ){
		minDiff = TMath::Abs(WMass12-MW);
		ind1 = buntag_indices[uj1];
		ind2 = buntag_indices[uj2];
	      }
	      
	    }
	  }	  
	  /////////////////////////////////////////////////////
	  type_       =  3;
	  nPermut_    = 12;
	  nPermut_alt_= 12;
	  meIntegrator->setIntType( MEIntegratorNew::SL2wj );
	  /////////////////////////////////////////////////////	  
	}
	else if( analyze_type6 || analyze_type6_BTag ){
	  
	  counter++;      
	  if(fixNumEvJob && !(counter>=evLow && counter<=evHigh) ) continue;
	  if(print) cout << "Processing event # " << counter << " (TYPE 6)." << " Event ID " << EVENT.event << endl;	   

	  /////////////////////////////////////////////////////
	  type_       =  6;
	  nPermut_    = 12;
	  nPermut_alt_= 12;
	  meIntegrator->setIntType( MEIntegratorNew::DL );
	  /////////////////////////////////////////////////////	  

	}
	else if( analyze_type7 ){
	  
	  counter++;      
	  if(fixNumEvJob && !(counter>=evLow && counter<=evHigh) ) continue;
	  if(print) cout << "Processing event # " << counter << " (TYPE 7)." << " Event ID " << EVENT.event << endl;	   
	  
	  /////////////////////////////////////////////////////
	  type_       =  7;
	  nPermut_    = 12;
	  nPermut_alt_= 12;
	  meIntegrator->setIntType( MEIntegratorNew::DL );
	  /////////////////////////////////////////////////////	  

	}
	else{ 
	  cout << "Inconsistency in the analysis... continue." << endl; 
	  continue; 
	}


	// sanity-check
	if(type_>=0 && type_<=3 && (ind1==999 || ind2==999)){
	  cout << "Inconsistency found: ind1 or ind2 are not set...continue." << endl;
	  continue;
	}

	// DEBUG
	if(debug){
	  cout << "*** Event ID " << EVENT.event << " *** " << endl;
	  cout << " ==> SL=" << int(properEventSL) << ", DL=" << properEventDL << endl;
	  cout << "     NJets " << numJets_ << " (" << numBTagM_ << " tagged)" << endl;
	  cout << "     b-tagged: " << endl;
	  for( unsigned int jj = 0; jj<btag_indices_backup.size(); jj++)
	    cout << "     (" 
		 << jets_p4[ btag_indices_backup[jj] ].Pt() << "," 
		 << jets_p4[ btag_indices_backup[jj] ].Eta() << ","
		 << jets_p4[ btag_indices_backup[jj] ].Phi() << "," 
		 << jets_p4[ btag_indices_backup[jj] ].M() << "), CSV= " 
		 << jets_csv[btag_indices_backup[jj] ] << endl;
	  cout << "     b-untagged: " << endl;
	  for( unsigned int jj = 0; jj<buntag_indices_backup.size(); jj++)
	    cout << "     (" 
		 << jets_p4[ buntag_indices_backup[jj] ].Pt() << "," 
		 << jets_p4[ buntag_indices_backup[jj] ].Eta() << ","
		 << jets_p4[ buntag_indices_backup[jj] ].Phi() << "," 
		 << jets_p4[ buntag_indices_backup[jj] ].M() << "), CSV= " 
		 << jets_csv[buntag_indices_backup[jj] ] << endl;
	  cout << "     btag probability is " << btag_LR_ << endl;
	  if(passes_btagshape){
	    cout << "     @@@@@ the jet collection has been re-ordered according to btag probability @@@@@@" << endl;
	    cout << "     b-tagged: " << endl;
	    for( unsigned int jj = 0; jj < btag_indices.size(); jj++)
	      cout << "     (" << jets_p4[ btag_indices[jj] ].Pt() << "," << jets_p4[ btag_indices[jj] ].Eta() << ","
		   << jets_p4[ btag_indices[jj] ].Phi() << "," << jets_p4[ btag_indices[jj] ].M() << "), CSV= " << jets_csv[ btag_indices[jj] ] << endl;
	    cout << "     b-untagged: " << endl;
	    for( unsigned int jj = 0; jj<buntag_indices.size(); jj++)
	      cout << "     (" << jets_p4[ buntag_indices[jj] ].Pt() << "," << jets_p4[ buntag_indices[jj] ].Eta() << ","
		   << jets_p4[ buntag_indices[jj] ].Phi() << "," << jets_p4[ buntag_indices[jj] ].M() << "), CSV= " << jets_csv[buntag_indices[jj] ] << endl;
	  }	
	  
	}
	
	

	// total number of integrations
	nTotInteg_      = nPermut_    * nMassPoints_;
	nTotInteg_alt_  = nPermut_alt_* nMassPoints_;
	
	
	// setup jet collection
	jets.clear();
	jets.push_back( leptonLV     );  
	jets.push_back( neutrinoLV   );  

	// b1,...,w1,w2 are indices for jets_p4 collection;
	// This is a map between the internal ordering bLep=2, W1Had=3, ..., higgs2 = 7, and jets_p4
	pos_to_index.clear();
	  
	if( type_==-3){	  
	  jets.push_back( jets_p4[ banytag_indices[0] ]);
	  jets.push_back( leptonLV2    ); 
	  jets.push_back( neutrinoLV   );       // dummy 
	  jets.push_back( jets_p4[ banytag_indices[1] ]);
	  jets.push_back( jets_p4[ banytag_indices[2] ]);
	  jets.push_back( jets_p4[ banytag_indices[3] ]);
	  
	  pos_to_index[2] = banytag_indices[0];
	  pos_to_index[3] = banytag_indices[0]; // dummy
	  pos_to_index[4] = banytag_indices[0]; // dummy
	  pos_to_index[5] = banytag_indices[1];
	  pos_to_index[6] = banytag_indices[2];
	  pos_to_index[7] = banytag_indices[3];	 
	}
	else if( type_==-2){
	  jets.push_back( jets_p4[ banytag_indices[0] ]);
	  jets.push_back( jets_p4[ banytag_indices[1] ]);  
	  jets.push_back( jets_p4[ banytag_indices[2] ]);  
	  jets.push_back( jets_p4[ banytag_indices[3] ]);
	  jets.push_back( jets_p4[ banytag_indices[4] ]);
	  jets.push_back( jets_p4[ banytag_indices[5] ]);
	  
	  pos_to_index[2] = banytag_indices[0];
	  pos_to_index[3] = banytag_indices[1];
	  pos_to_index[4] = banytag_indices[2];
	  pos_to_index[5] = banytag_indices[3];
	  pos_to_index[6] = banytag_indices[4];
	  pos_to_index[7] = banytag_indices[5];	 
	}
	else if( type_==-1){
	  jets.push_back( jets_p4[ banytag_indices[0] ]);
	  jets.push_back( jets_p4[ banytag_indices[1] ]);  
	  jets.push_back( jets_p4[ banytag_indices[1] ]); // dummy 
	  jets.push_back( jets_p4[ banytag_indices[2] ]);
	  jets.push_back( jets_p4[ banytag_indices[3] ]);
	  jets.push_back( jets_p4[ banytag_indices[4] ]);
	  
	  pos_to_index[2] = banytag_indices[0];
	  pos_to_index[3] = banytag_indices[1];
	  pos_to_index[4] = banytag_indices[1];           // dummy 
	  pos_to_index[5] = banytag_indices[2];
	  pos_to_index[6] = banytag_indices[3];
	  pos_to_index[7] = banytag_indices[4];	 
	}
	else if( type_<=3 && type_>=0){
	  jets.push_back( jets_p4[ btag_indices[0] ]  );
	  jets.push_back( jets_p4[ ind1 ]);  
	  jets.push_back( jets_p4[ ind2 ]);  
	  jets.push_back( jets_p4[ btag_indices[1] ]  );
	  jets.push_back( jets_p4[ btag_indices[2] ]  );
	  jets.push_back( jets_p4[ btag_indices[3] ]  );
	  
	  pos_to_index[2] = btag_indices[0];
	  pos_to_index[3] = ind1;
	  pos_to_index[4] = ind2;
	  pos_to_index[5] = btag_indices[1];
	  pos_to_index[6] = btag_indices[2];
	  pos_to_index[7] = btag_indices[3];	 
	}
	else if( type_==6 ){
	  jets.push_back( jets_p4[ btag_indices[0] ] );
	  jets.push_back( leptonLV2   );  
	  jets.push_back( neutrinoLV  );      // dummy  
	  jets.push_back( jets_p4[ btag_indices[1] ] );
	  jets.push_back( jets_p4[ btag_indices[2] ] );
	  jets.push_back( jets_p4[ btag_indices[3] ] );

	  pos_to_index[2] = btag_indices[0];
	  pos_to_index[3] = btag_indices[0];  // dummy
	  pos_to_index[4] = btag_indices[0];  // dummy
	  pos_to_index[5] = btag_indices[1];
	  pos_to_index[6] = btag_indices[2];
	  pos_to_index[7] = btag_indices[3];
	}
	else if( type_==7 ){
	  jets.push_back( jets_p4[ btagLoose_indices[0] ] );
	  jets.push_back( leptonLV2   );  
	  jets.push_back( neutrinoLV  );           // dummy  
	  jets.push_back( jets_p4[ btagLoose_indices[1] ] );
	  jets.push_back( jets_p4[ btagLoose_indices[2] ] );
	  jets.push_back( jets_p4[ btagLoose_indices[3] ] );

	  pos_to_index[2] = btagLoose_indices[0];
	  pos_to_index[3] = btagLoose_indices[0];  // dummy
	  pos_to_index[4] = btagLoose_indices[0];  // dummy
	  pos_to_index[5] = btagLoose_indices[1];
	  pos_to_index[6] = btagLoose_indices[2];
	  pos_to_index[7] = btagLoose_indices[3];
	}
	else{ /* ... */ }
	



	// save jet kinematics into the tree...
	nJet_ = 8;
	for(int q = 0; q < nJet_ ; q++ ){
	  // kinematics
	  jet_pt_ [q] = jets[q].Pt(); 
	  jet_eta_[q] = jets[q].Eta(); 
	  jet_phi_[q] = jets[q].Phi(); 	    
	  jet_m_  [q] = jets[q].M(); 
	  jet_csv_[q] = q>1 ? jets_csv[ pos_to_index[q] ] : -99.;
	}
	
	// set all prob. to 0.0;
	for(int p = 0 ; p < nTotInteg_; p++){
	  probAtSgn_permut_       [p] = 0.;
	  probAtSgnErr_permut_    [p] = 0.;
	}
	for(int p = 0 ; p < nPermut_; p++){
	  probAtSgn_bb_permut_ [p] = 0.;
	}
	for(int p = 0 ; p < nTotInteg_alt_; p++){
	  probAtSgn_alt_permut_   [p] = 0.;
	  probAtSgnErr_alt_permut_[p] = 0.;
	}
	for(int p = 0 ; p < nPermut_alt_; p++){
	  probAtSgn_bj_permut_ [p] = 0.;
	  probAtSgn_cc_permut_ [p] = 0.;
	  probAtSgn_jj_permut_ [p] = 0.;
	}
	
	/////////////////////////////////////////////////////////////
	
	// check if there is a tag-untag pair that satisfies the "cs-tag" 
	for( unsigned int w = 0; w<btag_indices.size(); w++){
	  
	  float m1 = ( jets_p4[btag_indices[w]] + jets_p4[ind1] ).M();
	  float m2 = ( jets_p4[btag_indices[w]] + jets_p4[ind2] ).M();
	  
	  if( ((m1>(MwL+5) && m1<(MwH-5)) || (m2>(MwL+5) && m2<(MwH-5))) && type_== 0 ){
	    flag_type0_ = 0; 
	    if( jets_csv[btag_indices[w]]<0.95 ) flag_type0_ = 1; 
	    if( jets_csv[btag_indices[w]]<0.90 ) flag_type0_ = 2; 
	    if( jets_csv[btag_indices[w]]<0.85 ) flag_type0_ = 3;
	    if( jets_csv[btag_indices[w]]<0.80 ) flag_type0_ = 4;  
	  }
	  if( ((m1>(MwL+5) && m1<(MwH-5)) || (m2>(MwL+5) && m2<(MwH-5))) && type_== 1 ){
	    flag_type1_ = 0; 
	    if( jets_csv[btag_indices[w]]<0.95 ) flag_type1_ = 1; 
	    if( jets_csv[btag_indices[w]]<0.90 ) flag_type1_ = 2; 
	    if( jets_csv[btag_indices[w]]<0.85 ) flag_type1_ = 3;
	    if( jets_csv[btag_indices[w]]<0.80 ) flag_type1_ = 4;  
	  }
	  if( ((m1>(MwL+5) && m1<(MwH-5)) || (m2>(MwL+5) && m2<(MwH-5))) && type_== 2 ){
	    flag_type2_ = 0; 
	    if( jets_csv[btag_indices[w]]<0.95 ) flag_type2_ = 1; 
	    if( jets_csv[btag_indices[w]]<0.90 ) flag_type2_ = 2; 
	    if( jets_csv[btag_indices[w]]<0.85 ) flag_type2_ = 3;
	    if( jets_csv[btag_indices[w]]<0.80 ) flag_type2_ = 4;  
	  }
	}
	// for type 3, the W-tag is different...
	float WMass = type_==3 ? (jets_p4[ ind1 ]+jets_p4[ ind2 ]).M() : -999.;
	if( WMass>MwL && WMass<MwH )  flag_type3_ = 1;
	
	
	/////////////////////////////////////////////////////////////
	
	// init reco particles
	meIntegrator->setJets(&jets);	

	// init MET stuff
	meIntegrator->setSumEt( METtype1p2corr.sumet );
	meIntegrator->setMEtCov(-99,-99,0);
	
	// specify if topLep has pdgid +6 or -6
	meIntegrator->setTopFlags( vLepton_charge[0]==1 ? +1 : -1 , vLepton_charge[0]==1 ? -1 : +1 );
	  
	// start the clock...	  
	clock->Start();
	
	// loop over Higgs mass values...
	for(int m = 0; m < nHiggsMassPoints ; m++){
	  meIntegrator->setMass( mH[m] );
	  
	  // loop over Top mass values...
	  for(int t = 0; t < nTopMassPoints ; t++){
	    meIntegrator->setTopMass( mT[t] , MW );
	    
	    // these are used for bookkeeping
	    double maxP_s = 0.;
	    double maxP_b = 0.;
	    
	    // number of permutations for which p has been calculated
	    int num_s = 0;
	    int num_b = 0;

	    // loop over hypothesis [ TTH, TTbb ]
	    for(int hyp = 0 ; hyp<2;  hyp++){
		    
	      // choose which permutations to consider;
	      int* permutList = 0;
	      if     ( type_ == -3 ) permutList = hyp==0 ?  permutations_4J_S      : permutations_4J_B;
	      else if( type_ == -2 ) permutList = hyp==0 ?  permutations_6J_S      : permutations_6J_B;
	      else if( type_ == -1 ) permutList = hyp==0 ?  permutations_5J_S      : permutations_5J_B;
	      else if( type_ ==  0 ) permutList = hyp==0 ?  permutations_TYPE0_S   : permutations_TYPE0_B;
	      else if( type_ ==  1 ) permutList = hyp==0 ?  permutations_TYPE1_S   : permutations_TYPE1_B;
	      else if( type_ ==  2 ) permutList = hyp==0 ?  permutations_TYPE2_S   : permutations_TYPE2_B;
	      else if( type_ ==  3 ) permutList = hyp==0 ?  permutations_TYPE0_S   : permutations_TYPE0_B;
	      else if( type_ >=  6 ) permutList = hyp==0 ?  permutations_TYPE6_S   : permutations_TYPE6_B;
	      else{ 
		cout << "No permutations found...continue." << endl; 
		continue; 
	      }
	      

	      // loop over permutations
	      for(unsigned int pos = 0; pos < (unsigned int)( hyp==0 ? nPermut_ : nPermut_alt_ ) ; pos++){
		
		// consider permutation #pos & save permutation-to-jet mas into the tree...
		meIntegrator->initVersors( permutList[pos] );
		if( hyp==0 ) perm_to_jet_    [pos] =  permutList[pos];
		if( hyp==1 ) perm_to_jet_alt_[pos] =  permutList[pos];
		
		// index of the four jets associated to b-quarks or W->qq
		int bLep_pos = (permutList[pos])%1000000/100000;
		int w1_pos   = (permutList[pos])%100000/10000;
		int w2_pos   = (permutList[pos])%10000/1000;
		int bHad_pos = (permutList[pos])%1000/100;
		int b1_pos   = (permutList[pos])%100/10;
		int b2_pos   = (permutList[pos])%10/1;       
		
		// find if this particular permutation matches the expectation
		int bLep_match = 0;
		if(TMath::Abs(TOPLEPB.Py())>0 && deltaR( jets_p4[ pos_to_index[bLep_pos] ], TOPLEPB ) < GENJETDR){
		  bLep_match = 1;
		}
		int w1_match = 0;  
		int w2_match = 0;
		if( (TMath::Abs(TOPHADW1.Py())>0 && deltaR( jets_p4[ pos_to_index[w1_pos] ], TOPHADW1 ) < GENJETDR) || 
		    (TMath::Abs(TOPHADW2.Py())>0 && deltaR( jets_p4[ pos_to_index[w1_pos] ], TOPHADW2 ) < GENJETDR)){
		  w1_match = 1;
		}
		if( (TMath::Abs(TOPHADW1.Py())>0 && deltaR( jets_p4[ pos_to_index[w2_pos] ], TOPHADW1 ) < GENJETDR) || 
		    (TMath::Abs(TOPHADW2.Py())>0 && deltaR( jets_p4[ pos_to_index[w2_pos] ], TOPHADW2 ) < GENJETDR)){
		  w2_match = 1;
		}
		int bHad_match = 0;
		if(TMath::Abs(TOPHADB.Py())>0 && deltaR( jets_p4[ pos_to_index[bHad_pos] ], TOPHADB ) < GENJETDR){
		  bHad_match = 1;
		}
		int b1_match = 0;  
		int b2_match = 0;
		if( (TMath::Abs(HIGGSB1.Py())>0 && deltaR( jets_p4[ pos_to_index[b1_pos] ], HIGGSB1 ) < GENJETDR) || 
		    (TMath::Abs(HIGGSB2.Py())>0 && deltaR( jets_p4[ pos_to_index[b1_pos] ], HIGGSB2 ) < GENJETDR)){
		  b1_match = 1;
		}
		if( (TMath::Abs(HIGGSB1.Py())>0 && deltaR( jets_p4[ pos_to_index[b2_pos] ], HIGGSB1 ) < GENJETDR) || 
		    (TMath::Abs(HIGGSB2.Py())>0 && deltaR( jets_p4[ pos_to_index[b2_pos] ], HIGGSB2 ) < GENJETDR)){
		  b2_match = 1;
		}
		
		// save an integer with this convention:
		//  > 1st digit = 1/0 if bLep candidate is correct/wrong
		//  > 2nd digit = 1/0 if matched/not matched to W-quark
		//  > ... 
		if( hyp==0 ) perm_to_gen_    [pos] = 100000*bLep_match + 10000*w1_match + 1000*w2_match + 100*bHad_match + 10*b1_match + 1*b2_match;
		if( hyp==1 ) perm_to_gen_alt_[pos] = 100000*bLep_match + 10000*w1_match + 1000*w2_match + 100*bHad_match + 10*b1_match + 1*b2_match;
		
		
		// check invariant mass of jet system:
		double mass, massLow, massHigh;
		bool skip        = !( meIntegrator->compatibilityCheck    (0.95, /*print*/ 0, mass, massLow, massHigh ) );
		bool skip_WHad   = false;
		bool skip_TopHad = false;
		if( type_==0 || type_==3 ){
		  skip_WHad   = !( meIntegrator->compatibilityCheck_WHad  (0.98, /*print*/ 0, mass, massLow, massHigh ) );
		  skip_TopHad = !( meIntegrator->compatibilityCheck_TopHad(0.98, /*print*/ 0, mass, massLow, massHigh ) );
		}
		
		
		// if use btag, determine the b-tag probability density
		if( useBtag ){	       			
		  double p_b_bLep =  jets_csv_prob_b[ pos_to_index[bLep_pos] ];
		  double p_b_bHad =  jets_csv_prob_b[ pos_to_index[bHad_pos] ];
		  double p_b_b1   =  jets_csv_prob_b[ pos_to_index[b1_pos] ];
		  double p_c_b1   =  jets_csv_prob_c[ pos_to_index[b1_pos] ];
		  double p_j_b1   =  jets_csv_prob_j[ pos_to_index[b1_pos] ];
		  double p_b_b2   =  jets_csv_prob_b[ pos_to_index[b2_pos] ];
		  double p_c_b2   =  jets_csv_prob_c[ pos_to_index[b2_pos] ];
		  double p_j_b2   =  jets_csv_prob_j[ pos_to_index[b2_pos] ];
		  double p_j_w1   =  1.0;
		  double p_j_w2   =  1.0;
		  
		  // the b-tag probability for the two untagged jets...
		  if( type_==0 || type_==3 || type_==-2){
		    p_j_w1 = jets_csv_prob_j[ pos_to_index[w1_pos] ];
		    p_j_w2 = jets_csv_prob_j[ pos_to_index[w2_pos] ];
		  }
		  // the b-tag probability for the one untagged jet...
		  else if( type_==1 || type_==2 || type_==-1){
		    p_j_w1 = jets_csv_prob_j[ pos_to_index[w1_pos] ];
		    p_j_w2 = 1.0;
		  }
		  // there are untagged jets...
		  else{
		    p_j_w1 = 1.0;
		    p_j_w2 = 1.0;
		  }
		  
		  // DEBUG
		  if( debug && 0 ){
		    cout << "Hyp=" << hyp <<"  [BTag M="   << numBTagM_ << "]" << endl;
		    cout << " bLep: p_b("   << jets_csv[pos_to_index[bLep_pos]] << ")=" << p_b_bLep;
		    cout << " bHad: p_b("   << jets_csv[pos_to_index[bHad_pos]] << ")=" << p_b_bHad;
		    cout << " b1  : p_b("   << jets_csv[pos_to_index[b1_pos]]   << ")=" << p_b_b1;
		    cout << " b2  : p_b("   << jets_csv[pos_to_index[b2_pos]]   << ")=" << p_b_b2;
		    if(!(type_>=6 || type_==-3)) 
		      cout << " w1  : p_j(" << jets_csv[pos_to_index[w1_pos]]   << ")=" << p_j_w1;
		    if(!(type_>=6 || type_==-3 || type_==1 || type_==2)) 
		      cout << " w2  : p_j(" << jets_csv[pos_to_index[w2_pos]]   << ")=" << p_j_w2 << endl;
		    cout << " P = "         <<  p_b_bLep * p_b_bHad * p_b_b1 * p_b_b2 * p_j_w1 * p_j_w2 << endl;
		  }
		  
		  // fill arrays with per-permutation probability
		  if(hyp==0){
		    probAtSgn_bb_permut_[pos] =  p_b_bLep * p_b_bHad * p_b_b1 * p_b_b2 * p_j_w1 * p_j_w2;
		  }
		  if(hyp==1){
		    probAtSgn_bj_permut_[pos] =  p_b_bLep * p_b_bHad * (p_b_b1 * p_j_b2 + p_j_b1 * p_b_b2 )*0.5 * p_j_w1 * p_j_w2;
		    probAtSgn_cc_permut_[pos] =  p_b_bLep * p_b_bHad * p_c_b1 * p_c_b2 * p_j_w1 * p_j_w2;
		    probAtSgn_jj_permut_[pos] =  p_b_bLep * p_b_bHad * p_j_b1 * p_j_b2 * p_j_w1 * p_j_w2;		  
		  }

		}
		
		// if doing scan over b-tag only, don't need amplitude...
		if( type_<0  ) continue;
		
		
		// if type 0/3 and incompatible with MW or MT (and we are not scanning vs MT) continue
		// ( this applies to both hypotheses )
		if( nTopMassPoints==1 && (skip_WHad || skip_TopHad) ){    
		  if(print) cout << "Skip (THad check failed)...          Perm. #" << pos << endl;
		  continue;
		}	      
		
		// retrieve intergration boundaries from meIntegrator
		pair<double, double> range_x0 = (meIntegrator->getW1JetEnergyCI(0.95));
		pair<double, double> range_x1 =  make_pair(-1,1);
		pair<double, double> range_x2 =  make_pair(-PI,PI);	    
		pair<double, double> range_x3 =  make_pair(-1,1);
		pair<double, double> range_x4 =  useMET ? (meIntegrator->getNuPhiCI(0.95)) : make_pair(-PI,PI);
		pair<double, double> range_x5 = (meIntegrator->getB1EnergyCI(0.95));
		pair<double, double> range_x6 = (meIntegrator->getB2EnergyCI(0.95));
	      
		// boundaries
		double x0L = range_x0.first; double x0U = range_x0.second;
		double x1L = range_x1.first; double x1U = range_x1.second;
		double x2L = range_x2.first; double x2U = range_x2.second;
		double x3L = range_x3.first; double x3U = range_x3.second;
		double x4L = range_x4.first; double x4U = range_x4.second;
		double x5L = range_x5.first; double x5U = range_x5.second;
		double x6L = range_x6.first; double x6U = range_x6.second;
	      
		// these hold for the sgn integration and type0...
		double xLmode0_s[4] = {x0L, x3L, x4L, x5L};
		double xUmode0_s[4] = {x0U, x3U, x4U, x5U};
		// these hold for the bkg integration and type0...
		double xLmode0_b[5] = {x0L, x3L, x4L, x5L, x6L};
		double xUmode0_b[5] = {x0U, x3U, x4U, x5U, x6U};		 
		
		// these hold for the sgn integration and type1...
		double xLmode1_s[6] = {x0L, x1L, x2L, x3L, x4L, x5L};
		double xUmode1_s[6] = {x0U, x1U, x2U, x3U, x4U, x5U};
		// these hold for the bkg integration and type1...
		double xLmode1_b[7] = {x0L, x1L, x2L, x3L, x4L, x5L, x6L};
		double xUmode1_b[7] = {x0U, x1U, x2U, x3U, x4U, x5U, x6U};
		
		// these hold for the sgn integration and type2...
		double xLmode2_s[6] = {x0L, x1L, x2L, x3L, x4L, x5L};
		double xUmode2_s[6] = {x0U, x1U, x2U, x3U, x4U, x5U};
		// these hold for the bkg integration and type2...	      
		double xLmode2_b[7] = {x0L, x1L, x2L, x3L, x4L, x5L, x6L};
		double xUmode2_b[7] = {x0U, x1U, x2U, x3U, x4U, x5U, x6U};	     	     	   	       
		
		// these hold for the sgn integration and type3...
		double xLmode3_s[4] = {x0L, x3L, x4L, x5L};
		double xUmode3_s[4] = {x0U, x3U, x4U, x5U};
		// these hold for the bkg integration and type3...
		double xLmode3_b[5] = {x0L, x3L, x4L, x5L, x6L};
		double xUmode3_b[5] = {x0U, x3U, x4U, x5U, x6U};
	      
		// these hold for the sgn integration and type6...
		double xLmode6_s[5] = {x1L, x2L, x1L, x2L, x5L};
		double xUmode6_s[5] = {x1U, x2U, x1U, x2U, x5U};
		// these hold for the bkg integration and type6...
		double xLmode6_b[6] = {x1L, x2L, x1L, x2L, x5L, x6L};
		double xUmode6_b[6] = {x1U, x2U, x1U, x2U, x5U, x6U};
	      
		// these hold for the sgn integration and type7...
		double xLmode7_s[5] = {x1L, x2L, x1L, x2L, x5L};
		double xUmode7_s[5] = {x1U, x2U, x1U, x2U, x5U};
		// these hold for the bkg integration and type7...
		double xLmode7_b[6] = {x1L, x2L, x1L, x2L, x5L, x6L};
		double xUmode7_b[6] = {x1U, x2U, x1U, x2U, x5U, x6U};
	      
		// number of integration variables (TTH hypothesis)
		int nParam;
		if     ( type_==0 )  nParam = 4; // Eq, eta_nu, phi_nu, Eb
		else if( type_==1 )  nParam = 6; // Eq, eta_q', phi_q', eta_nu, phi_nu, Eb
		else if( type_==2 )  nParam = 6; // Eq, eta_q', phi_q', eta_nu, phi_nu, Eb
		else if( type_==3 )  nParam = 4; // Eq, eta_nu, phi_nu, Eb
		else if( type_==6 )  nParam = 5; // eta_nu, phi_nu, eta_nu', phi_nu', Eb
		else if( type_==7 )  nParam = 5; // eta_nu, phi_nu, eta_nu', phi_nu', Eb
		else{
		  cout << "No type match...contine." << endl;
		  continue;
		}
		
		// the per-permutation probability...
		double p    = 0.;	
		// the per-permutation probability error...
		double pErr = 0.;	       	
		
	     	  	
		// if doing higgs mass scan, don't consider bkg hypo
		if( nHiggsMassPoints>1 && hyp==1 ) continue;
		
		// if doing top mass scan, don't consider sgn hypo
		if( nTopMassPoints>1 && hyp==0 ) continue;
		
		// if consider only one hypothesis (SoB=0) 
		// and the current hypo is not the desired one, continue...
		if( SoB==0 && hyp!=hypo) continue;
		
		// if current hypo is TTH, but M(b1b2) incompatible with 125
		// (and we are not scanning vs MH) continue...
		if( hyp==0 && nHiggsMassPoints==1 && skip){
		  if(print){
		    cout << "Skip    hypo " << (hyp==0 ? "ttH " : "ttbb") 
			 << " [MH=" << mH[m] << ", MT=" << mT[t]
			 << "] Perm. #" << pos;
		    cout << " => p=" << p << endl;
		  }
		  continue;
		}
		
		// increment counters
		if(hyp==0) num_s++;
		if(hyp==1) num_b++;

		// if NOT doing top mass scan, and hyp==1
		// and no permutations for hyp==0 accepted, continue (the weight will be 0 anyway)
		if( nTopMassPoints==1 && hyp==1 && num_s==0){
		  if(print){
		    cout << "Skip    hypo " << (hyp==0 ? "ttH " : "ttbb") 
			 << " because no valid ttH permutations found" << endl;
		  }
		  continue;
		}

		// setup hypothesis
		if(print){
		  cout << "Testing hypo " << (hyp==0 ? "ttH " : "ttbb") 
		       << " [MH=" << mH[m] << ", MT=" << mT[t]
		       << "] Perm. #" << pos;	       
		}
		meIntegrator->setHypo(hyp);
		  
		// initial number of function calles
		int intPoints = 4000;
		if( type_==0 )  intPoints =  2000;
		if( type_==1 )  intPoints =  4000;
		if( type_==2 )  intPoints =  4000;
		if( type_==3 )  intPoints =  2000;
		if( type_==6 )  intPoints = 10000;
		if( type_==7 )  intPoints = 10000;
		
		  // count how many time the integration is rerun per permutation
		int ntries = 0;
		
		// skip ME calculation... for debugging
		if(speedup==0){
		  
		  // refinement: redo integration if it returned a bad chi2
		  while( ntries < MAX_REEVAL_TRIES){
		    
		    // integrand
		    ROOT::Math::Functor toIntegrate(meIntegrator, &MEIntegratorNew::Eval, nParam+hyp);
		    // VEGAS integrator
		    ROOT::Math::GSLMCIntegrator ig2( ROOT::Math::IntegrationMultiDim::kVEGAS , 1.e-12, 1.e-5, intPoints);
		    ig2.SetFunction(toIntegrate);
		    // setup # of parameters
		    meIntegrator->SetPar( nParam+hyp );	 
		    
		    // the integration ranges depend on hyp and type
		    if     ( type_==0 ){
		      p = (hyp==0 ? ig2.Integral(xLmode0_s, xUmode0_s) : ig2.Integral(xLmode0_b, xUmode0_b));
		    }
		    else if( type_==1 ){
		      p = (hyp==0 ? ig2.Integral(xLmode1_s, xUmode1_s) : ig2.Integral(xLmode1_b, xUmode1_b));
		    }
		    else if( type_==2 ){
		      p = (hyp==0 ? ig2.Integral(xLmode2_s, xUmode2_s) : ig2.Integral(xLmode2_b, xUmode2_b));
		    }
		    else if( type_==3 ){
		      p = (hyp==0 ? ig2.Integral(xLmode3_s, xUmode3_s) : ig2.Integral(xLmode3_b, xUmode3_b));
		    }
		    else if( type_==6 ){
		      p = (hyp==0 ? ig2.Integral(xLmode6_s, xUmode6_s) : ig2.Integral(xLmode6_b, xUmode6_b));
		    }
		    else if( type_==7 ){
		      p = (hyp==0 ? ig2.Integral(xLmode7_s, xUmode7_s) : ig2.Integral(xLmode7_b, xUmode7_b));
		    }
		      else{ /* ... */ }		    
		    
		    // chi2 of the integration
		    double chi2 =  ig2.ChiSqr();
		    
		    // error from the last VEGAS iteration
		    pErr =  ig2.Error();
		    
		    // check if the actual permutation returned a small or large number...
		    // if the value is less than 10% of the largest found that far, skip 
		    // the improvement
		    if( hyp==0 ){
		      if( p>maxP_s ) maxP_s = p;		  
		      else if( p<0.1*maxP_s) ntries = MAX_REEVAL_TRIES+1;
		    }
		    else{
		      if( p>maxP_b ) maxP_b = p;		  
		      else if( p<0.1*maxP_b) ntries = MAX_REEVAL_TRIES+1;	
		    }
		    
		    // if the chi2 is bad, increse # of points and repeat the integration
		    if( chi2 > maxChi2_ ){
		      ntries++;
		      intPoints *= 1.5;
		    }
		    // otherwise, just go to the next permutation...
		    else 
		      ntries = MAX_REEVAL_TRIES+1;			
		  }	
		}
		else{
		  // can still be interested in b-tagging, so set p=1...
		  p = 1.;
		}
		  
		// to avoid problems
		if( TMath::IsNaN(p) )    p    = 0.;
		if( TMath::IsNaN(pErr) ) pErr = 0.;
		if(print) cout << " => p=(" << p << " +/- " << pErr << ")" << endl;
		
		/////////////////////////////////////////////////////////////
		  
		if( useBtag ){
		  
		  // total ME*btag probability for nominal MH and MT, sgn hypo
		  if( hyp==0 && mH[m]<MH+1.5 && mH[m]>MH-1.5 && mT[t]<MT+1.5 && mT[t]>MT-1.5){
		    probAtSgn_ttbb_         += ( p * probAtSgn_bb_permut_[pos] );
		  }
		  
		  // total ME*btag probability for nominal MH and MT, bkg hypo
		  if( hyp==1 && mH[m]<MH+1.5 && mH[m]>MH-1.5 && mT[t]<MT+1.5 && mT[t]>MT-1.5){
		    probAtSgn_alt_ttbb_     += ( p * probAtSgn_bb_permut_[pos] );
		    probAtSgn_alt_ttbj_     += ( p * probAtSgn_bj_permut_[pos] );
		    probAtSgn_alt_ttcc_     += ( p * probAtSgn_cc_permut_[pos] );
		    probAtSgn_alt_ttjj_     += ( p * probAtSgn_jj_permut_[pos] );		      
		  }		      		      
		}
		
		/////////////////////////////////////////////////////////////
		
		// save TTH prob. per permutation AND mH
		if( hyp==0 && mT[t]<MT+1.5 && mT[t]>MT-1.5){
		  probAtSgn_permut_   [(unsigned int)(pos+m*nPermut_)] = p;
		  probAtSgnErr_permut_[(unsigned int)(pos+m*nPermut_)] = pErr;
		}
		// save TTbb prob. per permutation AND mT
		if( hyp==1 && mH[m]<MH+1.5 && mH[m]>MH-1.5 ){
		  probAtSgn_alt_permut_   [(unsigned int)(pos+t*nPermut_alt_)] = p;
		  probAtSgnErr_alt_permut_[(unsigned int)(pos+t*nPermut_alt_)] = pErr;
		}
		
		// total and per-permutation ME probability for nominal MH and MT
		if( mH[m]<MH+1.5 && mH[m]>MH-1.5 && mT[t]<MT+1.5 && mT[t]>MT-1.5){
		  if(hyp==0){
		    probAtSgn_     += p;
		  }
		  else{
		    probAtSgn_alt_ += p;
		  }
		}
		
		/////////////////////////////////////////////////////////////
		  
	      }  // hypothesis		  
	    }  // permutations
	  }  // nTopMassPoints
	}  // nHiggsMassPoints	    
	
	// stop clock and reset
	clock->Stop();
	time_ = clock->RealTime();
	clock->Reset();	    
	
	// this time is integrated over all hypotheses, permutations, and mass scan values
	if(print) cout << "Done in " << time_ << " sec" << endl;
      }

      ///////////////////////////////////////////////////
      // ALL THE REST...                               //
      ///////////////////////////////////////////////////

      else{
	continue;
      }
      
      counter_     = counter;
      nSimBs_      = nSimBs;
      weight_      = scaleFactor;
      nPVs_        = nPVs;
      Vtype_       = Vtype;

      EVENT_.run   = EVENT.run;
      EVENT_.lumi  = EVENT.lumi;
      EVENT_.event = EVENT.event;
      EVENT_.json  = EVENT.json;

      PUweight_    = PUweight;
      PUweightP_   = PUweightP;
      PUweightM_   = PUweightM;

      lheNj_       = lheNj;
      
      trigger_     = weightTrig2012;
      for(int k = 0; k < 70 ; k++){
	triggerFlags_[k] = triggerFlags[k];
      }
      
      // fill the tree...
      tree->Fill();      

    } // nentries

    // this histogram keeps track of the fraction of analyzed events per sample
    hcounter->SetBinContent(1,float(events_)/nentries);    

  } // samples



  // save the tree and the counting histo in the ROOT file
  fout_tmp->cd();
  hcounter->Write("",TObject::kOverwrite );
  tree->Write("",TObject::kOverwrite );
  fout_tmp->Close();

  // delete MEIntegrator and other allocated objects
  cout << "Delete meIntegrator..." << endl;
  delete meIntegrator;
  delete clock;
  cout << "Finished!!!" << endl;
  
  // Done!!!
  return 0;
}
