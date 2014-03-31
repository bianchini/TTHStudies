#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include <cstdlib>
#include <iostream> 
#include <fstream>
#include <map>
#include <string>

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


#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooLandau.h"
#include "RooUniform.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooFitResult.h"
#include "RooKeysPdf.h"
#include "RooBifurGauss.h"
#include "RooVoigtian.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooGaussian.h"
#include "RooAbsCategory.h"
#include "RooNumConvPdf.h"
#include "RooFFTConvPdf.h"
#include "RooTruthModel.h"
#include "RooDataHist.h"
#include "RooGenericPdf.h"
#include "RooPolynomial.h"
#include "RooBernstein.h"
#include "RooNDKeysPdf.h"
#include "RooChiSquarePdf.h"

#include "Bianchi/TTHStudies/interface/MEIntegratorNew.h"
#include "Bianchi/TTHStudies/interface/Samples.h"
#include "Bianchi/TTHStudies/interface/CalcME.h"

#define GENJETDR 0.3
#define N_TRIES_MAX 10
#define MAX_REEVAL_TRIES   3
#define VERBOSE  false
#define PI TMath::Pi()

using namespace std;
using namespace RooFit;

typedef struct 
{
  float et; 
  float sumet;
  float sig;
  float phi;
} metInfo;


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


void commitVariables( TLorentzVector b1_, TLorentzVector b2_, float csv1_, float csv2_,
		      float& dR,   float& dPhi, float& dEta, float& Mjj, 
		      float& Pt,   float& Eta,
		      float& pt1,  float& pt2,  float& eta1, float& eta2,
		      float& csv1, float& csv2){
  
  dR   = deltaR(b1_,b2_);
  dPhi = TMath::ACos( TMath::Cos(b1_.Phi()-b2_.Phi()) ); 
  dEta = TMath::Abs(b1_.Eta() - b2_.Eta());
  Mjj  = (b1_+b2_).M();
  Pt   = (b1_+b2_).Pt();
  Eta  = (b1_+b2_).Eta();
  pt1  = b1_.Pt();
  pt2  = b2_.Pt();
  eta1 = b1_.Eta();
  eta2 = b2_.Eta();
  csv1 = csv1_;
  csv2 = csv2_;

  return;

}





int main(int argc, const char* argv[])
{

  std::cout << "MEAnalysis" << std::endl;
  gROOT->SetBatch(true);
 
  gSystem->Load("libFWCoreFWLite");
  gSystem->Load("libDataFormatsFWLite");

  AutoLibraryLoader::enable();

  PythonProcessDesc builder(argv[1]);
  const edm::ParameterSet& in = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteInput");

  const edm::VParameterSet& samples = in.getParameter<edm::VParameterSet>("samples") ;
  std::string outFileName  ( in.getParameter<std::string>  ("outFileName" ) );
  std::string pathToFile   ( in.getParameter<std::string>  ("pathToFile" ) );
  std::string ordering     ( in.getParameter<std::string>  ("ordering" ) );
  std::string pathToTF     ( in.getParameter<std::string>  ("pathToTF"   ) );
  std::string pathToCP     ( in.getParameter<std::string>  ("pathToCP"   ) );

  bool verbose             ( in.getParameter<bool>  ("verbose") );
  double lumi              ( in.getParameter<double>("lumi") );

  float met       ( in.getParameter<double> ("met")        );

  float MwL       ( in.getUntrackedParameter<double> ("MwL", 60)  );
  float MwH       ( in.getUntrackedParameter<double> ("MwH", 100)  );

  float MhL       ( in.getUntrackedParameter<double> ("MhL", 100)  );
  float MhH       ( in.getUntrackedParameter<double> ("MhH", 140)  );

  int   switchoffOL( in.getUntrackedParameter<int>    ("switchoffOL", 0));
  int   speedup    ( in.getUntrackedParameter<int>    ("speedup", 0));

  int   doType0    ( in.getUntrackedParameter<int>    ("doType0", 0));
  int   doType1    ( in.getUntrackedParameter<int>    ("doType1", 0));
  int   doType2    ( in.getUntrackedParameter<int>    ("doType2", 0));
  int   doType3    ( in.getUntrackedParameter<int>    ("doType3", 0));
  int   doType4    ( in.getUntrackedParameter<int>    ("doType4", 0));
  int   doType6    ( in.getUntrackedParameter<int>    ("doType6", 0));
  int   doType7    ( in.getUntrackedParameter<int>    ("doType7", 0));

  int doubleGaussianB  ( in.getUntrackedParameter<int> ("doubleGaussianB", 0));
  int useBtag          ( in.getUntrackedParameter<int> ("useBtag", 0));


  int   useME      ( in.getParameter<int>    ("useME")      );
  int   useJac     ( in.getParameter<int>    ("useJac")     );
  int   useMET     ( in.getParameter<int>    ("useMET")     );
  int   useTF      ( in.getParameter<int>    ("useTF")      );
  int   usePDF     ( in.getParameter<int>    ("usePDF")     );

  int   printP4         ( in.getParameter<int>    ("printP4")          );
  int   norm            ( in.getUntrackedParameter<int>    ("norm", 0) );

  double maxChi2_       ( in.getUntrackedParameter<double>  ("maxChi2", 2.5));

  int   hypo                 ( in.getUntrackedParameter<int>    ("hypo", 0)     );
  int   SoB                  ( in.getUntrackedParameter<int>    ("SoB",  1)     );

  vector<double> masses       ( in.getParameter<vector<double> >  ("masses")        );
  vector<string> functions    ( in.getParameter<vector<string> >  ("functions")     );

  vector<int> evLimits( in.getParameter<vector<int> >  ("evLimits") );
  int evLow  = evLimits[0];
  int evHigh = evLimits[1];
  int fixNumberOfEventsPerJob  ( in.getUntrackedParameter<int>    ("fixNumberOfEventsPerJob", 1));

  int  doCSVup             ( in.getUntrackedParameter<int>    ("doCSVup",   0) );
  int  doCSVdown           ( in.getUntrackedParameter<int>    ("doCSVdown", 0) );
  int  doJECup             ( in.getUntrackedParameter<int>    ("doJECup",   0) );
  int  doJECdown           ( in.getUntrackedParameter<int>    ("doJECdown", 0) );

  gSystem->Exec(("rm "+outFileName).c_str());

  TStopwatch* clock = new TStopwatch();
  TRandom3*   ran   = new TRandom3();
  ran->SetSeed(4321);

  TFile* fCP = TFile::Open(pathToCP.c_str(),"READ");

  double svs[4] = {0.27, 0.73, 0.93, 0.07};
  map<string,TH1F*> btagger; 
  if( useBtag && fCP!=0 ){
    btagger["b_Bin0"] = fCP->Get("csv_b_Bin0__csvReco")!=0 ? (TH1F*)fCP->Get("csv_b_Bin0__csvReco") : 0;
    btagger["b_Bin1"] = fCP->Get("csv_b_Bin1__csvReco")!=0 ? (TH1F*)fCP->Get("csv_b_Bin1__csvReco") : 0;
    btagger["l_Bin0"] = fCP->Get("csv_l_Bin0__csvReco")!=0 ? (TH1F*)fCP->Get("csv_l_Bin0__csvReco") : 0;
    btagger["l_Bin1"] = fCP->Get("csv_l_Bin1__csvReco")!=0 ? (TH1F*)fCP->Get("csv_l_Bin1__csvReco") : 0;
  }
  else if( useBtag && fCP!=0 ){
    cout << "Cound not find " << pathToCP << ": exit" << endl;
    return 0;
  }

  const int nMassPoints  = masses.size();
  double mH[nMassPoints];
  for( unsigned int m = 0; m < masses.size() ; m++)
    mH[m] = masses[m];


  int permutations_SL2wj[12] =  
    {234567, 534267,     // CORRECT 
     634725, 734625,     // ALL WRONG
     534762, 734562,     // ONE WRONG
     234765, 734265,     // ONE WRONG
     634572, 534672,     // ONE WRONG
     234675, 634275      // ONE WRONG
    };

  int permutations_SL1wj[12] =  
    {234567, 534267,     // CORRECT 
     634725, 734625,     // ALL WRONG
     534762, 734562,     // ONE WRONG
     234765, 734265,     // ONE WRONG
     634572, 534672,     // ONE WRONG
     234675, 634275      // ONE WRONG
    };
  int permutations_SLNoBHad[3] =  
    {
      234567,     // CORRECT 
      734562,     // ONE WRONG
      634572      // ONE WRONG
    };
  int permutations_SLNoBLep[3] =
    {
      234567,      // CORRECT 
      234765,      // ONE WRONG
      234675       // ONE WRONG
    };
  int permutations_SLNoHiggs[6] =
    {
      234567,      // ONE CORRECT
      534267,      // ONE CORRECT
      534627,      // ALL WRONG
      634527,      // ALL WRONG
      234657,      // ALL WRONG
      634257       // ALL WRONG
    };

  int permutations_SL3b[12] =  
    {234567, 534267,     // CORRECT 
     634725, 734625,     // ALL WRONG
     534762, 734562,     // ONE WRONG
     234765, 734265,     // ONE WRONG
     634572, 534672,     // ONE WRONG
     234675, 634275      // ONE WRONG
    };

  int permutations_SL3b_equivalent[12] =  
    {234576, 534276,     // CORRECT 
     634752, 734652,     // ALL WRONG
     534726, 734526,     // ONE WRONG
     234756, 734256,     // ONE WRONG
     634527, 534627,     // ONE WRONG
     234657, 634257      // ONE WRONG
    };

  int permutations_DL[12] =  
    {234567, 534267,     // CORRECT 
     634725, 734625,     // ALL WRONG
     534762, 734562,     // ONE WRONG
     234765, 734265,     // ONE WRONG
     634572, 534672,     // ONE WRONG
     234675, 634275      // ONE WRONG
    };
  

  ///////////////////////////////////////////////////////

  TFile* fout_tmp = TFile::Open(outFileName.c_str(),"UPDATE");

  TH1F*  hcounter = new TH1F("hcounter","",1,0,1);
  int events_ = 0;

  TTree* tree  = new TTree("tree","");

  int counter_;
  int matchesH_;
  int matchesW_;
  int matchesT_;
  int matchesHAll_;
  int matchesWAll_;
  int matchesTAll_;
  int overlapLight_;
  int overlapHeavy_;
  int type_;
  int nSimBs_, nC_, nCTop_;
  int nMatchSimBs_;

  int flag_type0_;
  int flag_type1_;
  int flag_type2_;
  int flag_type3_;
  int flag_type4_;
  int flag_type6_;

  float time_;
  float mL_ [2];
  float ptL_[2];
  float csvL_[2];
  float csvH_[4];
  float probAtSgn_;
  float probAtSgn_alt_;

  float probAtSgn_ttbb_;
  float probAtSgn_alt_ttbb_;
  float probAtSgn_alt_ttjj_;

  float weight_;

  float dR_,dPhi_,dEta_,Mjj_,Pt_,Eta_,pt1_,pt2_,eta1_,eta2_,csv1_,csv2_;


  tree->Branch("counter",  &counter_, "counter/I");
  tree->Branch("matchesH", &matchesH_, "matchesH/I");
  tree->Branch("matchesW", &matchesW_, "matchesW/I");
  tree->Branch("matchesT", &matchesT_, "matchesT/I");
  tree->Branch("matchesHAll", &matchesHAll_, "matchesHAll/I");
  tree->Branch("matchesWAll", &matchesWAll_, "matchesWAll/I");
  tree->Branch("matchesTAll", &matchesTAll_, "matchesTAll/I");
  tree->Branch("overlapLight", &overlapLight_, "overlapLight/I");
  tree->Branch("overlapHeavy", &overlapHeavy_, "overlapHeavy/I");
  tree->Branch("csvL",     csvL_,     "csvL[2]/F");
  tree->Branch("mL",       mL_,       "mL[2]/F");
  tree->Branch("ptL",      ptL_,      "ptL[2]/F");
  tree->Branch("csvH",     csvH_,     "csvH[2]/F");

  tree->Branch("type",     &type_,    "type/I");
  tree->Branch("nSimBs",   &nSimBs_,  "nSimBs/I");
  tree->Branch("nMatchSimBs",   &nMatchSimBs_,  "nMatchSimBs/I");
  tree->Branch("nC",       &nC_,      "nC/I");
  tree->Branch("nCTop",    &nCTop_,   "nCTop/I");

  tree->Branch("weight",   &weight_,  "weight/F");
  tree->Branch("time",     &time_,    "time/F");


  tree->Branch("flag_type0",  &flag_type0_, "flag_type0/I");
  tree->Branch("flag_type1",  &flag_type1_, "flag_type1/I");
  tree->Branch("flag_type2",  &flag_type2_, "flag_type2/I");
  tree->Branch("flag_type3",  &flag_type3_, "flag_type3/I");
  tree->Branch("flag_type4",  &flag_type4_, "flag_type4/I");
  tree->Branch("flag_type6",  &flag_type6_, "flag_type6/I");

  tree->Branch("dR",       &dR_,    "dR/F");
  tree->Branch("dPhi",     &dPhi_,  "dPhi/F");
  tree->Branch("dEta",     &dEta_,  "dEta/F");
  tree->Branch("Mjj",      &Mjj_,   "Mjj/F");
  tree->Branch("Pt",       &Pt_,    "Pt/F");
  tree->Branch("Eta",      &Eta_,   "Eta/F");
  tree->Branch("pt1",      &pt1_,   "pt1/F");
  tree->Branch("pt2",      &pt2_,   "pt2/F");
  tree->Branch("eta1",     &eta1_,  "eta1/F");
  tree->Branch("eta2",     &eta2_,  "eta2/F");
  tree->Branch("csv1",     &csv1_,  "csv1/F");
  tree->Branch("csv2",     &csv2_,  "csv2/F");

  tree->Branch(Form("p_%d_all_s",int(met)),   &probAtSgn_,     Form("p_%d_all_s/F",int(met)) );
  tree->Branch(Form("p_%d_all_b",int(met)),   &probAtSgn_alt_, Form("p_%d_all_b/F",int(met)) );

  tree->Branch(Form("p_%d_all_s_ttbb",int(met)),   &probAtSgn_ttbb_,     Form("p_%d_all_s_ttbb/F",int(met)) );
  tree->Branch(Form("p_%d_all_b_ttbb",int(met)),   &probAtSgn_alt_ttbb_, Form("p_%d_all_b_ttbb/F",int(met)) );
  tree->Branch(Form("p_%d_all_b_ttjj",int(met)),   &probAtSgn_alt_ttjj_, Form("p_%d_all_b_ttjj/F",int(met)) );

  ///////////////////////////////////////////////////////


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
  meIntegrator->setNormFormulas( TString(functions[0].c_str()),  
				 TString(functions[1].c_str()),  
				 TString(functions[2].c_str()),
				 TString(functions[3].c_str()),  
				 TString(functions[4].c_str()),  
				 TString(functions[5].c_str())
				 );

  meIntegrator->setTopMass( 174.3 , 80.19);
  meIntegrator->setUseME (useME);
  meIntegrator->setUseJac(useJac);
  meIntegrator->setUseMET(useMET);
  meIntegrator->setUseTF (useTF);
  meIntegrator->setUsePDF(usePDF);

  meIntegrator->setUseRefinedTF(doubleGaussianB);

  meIntegrator->initTFparameters(1.0,1.0,1.0,1.0, 1.0);

  if(switchoffOL){
    meIntegrator->switchOffOL(); 
    cout << "*** Switching off OpenLoops to speed-up the calculation ***" << endl;
  }



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


  ////////////////////////////////////////////////////////



  for(unsigned int sample = 0 ; sample < mySampleFiles.size(); sample++){
    
    string currentName       = mySampleFiles[sample];

    mySamples->OpenFile( currentName );
    cout << "Opening file " << currentName << endl;
    TTree* currentTree       = mySamples->GetTree( currentName, "tree");
    cout << "Done!!" << endl;
    float scaleFactor        = mySamples->GetWeight(currentName);


    genParticleInfo genB, genBbar;
    genTopInfo genTop, genTbar;
    metInfo METtype1p2corr;
    JetByPt jet1, jet2, jet3, jet4, jet5, jet6, jet7, jet8, jet9, jet10;
    int nvlep, nSimBs, nC, nCTop, nhJets, naJets;
    Float_t vLepton_mass[999];
    Float_t vLepton_pt  [999];
    Float_t vLepton_eta [999];
    Float_t vLepton_phi [999];
    Float_t vLepton_charge    [999];
    Float_t vLepton_pfCorrIso [999];
    Int_t   vLepton_type      [999];


    Float_t hJet_csv_nominal[999];
    Float_t hJet_csv_upBC[999];
    Float_t hJet_csv_downBC[999];
    Float_t hJet_csv_upL[999];
    Float_t hJet_csv_downL[999];
    Float_t hJet_JECUnc[999];

    Float_t hJet_genPt [999];
    Float_t hJet_genEta[999];
    Float_t hJet_genPhi[999];

    Float_t aJet_csv_nominal[999];
    Float_t aJet_csv_upBC[999];
    Float_t aJet_csv_downBC[999];
    Float_t aJet_csv_upL[999];
    Float_t aJet_csv_downL[999];
    Float_t aJet_JECUnc[999];

    Float_t aJet_genPt [999];
    Float_t aJet_genEta[999];
    Float_t aJet_genPhi[999];

    int nSvs;
    float SvmassSv[999];
    float Svpt    [999];
    float Sveta   [999];
    float Svphi   [999];

    float SimBsmass[999];
    float SimBspt  [999];
    float SimBseta [999];
    float SimBsphi [999];

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
    currentTree->SetBranchAddress("METtype1p2corr",&METtype1p2corr);
    currentTree->SetBranchAddress("vLepton_charge",vLepton_charge);
    currentTree->SetBranchAddress("vLepton_mass"  ,vLepton_mass);
    currentTree->SetBranchAddress("vLepton_pt"    ,vLepton_pt);
    currentTree->SetBranchAddress("vLepton_eta"   ,vLepton_eta);
    currentTree->SetBranchAddress("vLepton_phi"   ,vLepton_phi);
    currentTree->SetBranchAddress("vLepton_charge",vLepton_charge);
    currentTree->SetBranchAddress("vLepton_pfCorrIso",vLepton_pfCorrIso);
    currentTree->SetBranchAddress("vLepton_type",  vLepton_type);

    currentTree->SetBranchAddress("hJet_csv_nominal", hJet_csv_nominal);
    currentTree->SetBranchAddress("hJet_csv_upBC",    hJet_csv_upBC);
    currentTree->SetBranchAddress("hJet_csv_downBC",  hJet_csv_downBC);
    currentTree->SetBranchAddress("hJet_csv_upL",     hJet_csv_upL);
    currentTree->SetBranchAddress("hJet_csv_downL",   hJet_csv_downL);
    currentTree->SetBranchAddress("hJet_JECUnc",      hJet_JECUnc);
    currentTree->SetBranchAddress("hJet_genPt",       hJet_genPt);
    currentTree->SetBranchAddress("hJet_genEta",      hJet_genEta);
    currentTree->SetBranchAddress("hJet_genPhi",      hJet_genPhi);

    currentTree->SetBranchAddress("aJet_csv_nominal", aJet_csv_nominal);
    currentTree->SetBranchAddress("aJet_csv_upBC",    aJet_csv_upBC);
    currentTree->SetBranchAddress("aJet_csv_downBC",  aJet_csv_downBC);
    currentTree->SetBranchAddress("aJet_csv_upL",     aJet_csv_upL);
    currentTree->SetBranchAddress("aJet_csv_downL",   aJet_csv_downL);
    currentTree->SetBranchAddress("aJet_JECUnc",      aJet_JECUnc);
    currentTree->SetBranchAddress("aJet_genPt",       aJet_genPt);
    currentTree->SetBranchAddress("aJet_genEta",      aJet_genEta);
    currentTree->SetBranchAddress("aJet_genPhi",      aJet_genPhi);

    currentTree->SetBranchAddress("nhJets",      &nhJets);
    currentTree->SetBranchAddress("naJets",      &naJets);

    currentTree->SetBranchAddress("nSimBs",  &nSimBs);
    currentTree->SetBranchAddress("nvlep",   &nvlep);
    currentTree->SetBranchAddress("nC",      &nC);
    currentTree->SetBranchAddress("nCTop",   &nCTop);

    currentTree->SetBranchAddress("nSvs",      &nSvs);
    currentTree->SetBranchAddress("Sv_massSv", SvmassSv);
    currentTree->SetBranchAddress("Sv_pt",     Svpt);
    currentTree->SetBranchAddress("Sv_eta",    Sveta);
    currentTree->SetBranchAddress("Sv_phi",    Svphi);

    currentTree->SetBranchAddress("SimBs_mass",   SimBsmass);
    currentTree->SetBranchAddress("SimBs_pt",     SimBspt);
    currentTree->SetBranchAddress("SimBs_eta",    SimBseta);
    currentTree->SetBranchAddress("SimBs_phi",    SimBsphi);

    
    int counter = 0;
    Long64_t nentries = currentTree->GetEntries();

    if( evHigh<0 ) evHigh = nentries;
    for (Long64_t i = 0; i < nentries ; i++){
      
      if(counter>evHigh && fixNumberOfEventsPerJob) continue;
      if(!fixNumberOfEventsPerJob && !(i>=evLow && i<evHigh) ) continue;
      events_++;


      if(i%5000==0) cout << i << endl;
      currentTree->GetEntry(i);
      
      //continue;
      
      TLorentzVector topBLV(1,0,0,1)   ;
      TLorentzVector topW1LV(1,0,0,1)  ; 
      TLorentzVector topW2LV(1,0,0,1)  ;
      TLorentzVector atopBLV(1,0,0,1)  ; 
      TLorentzVector atopW1LV(1,0,0,1) ; 
      TLorentzVector atopW2LV(1,0,0,1) ;
      TLorentzVector genBLV(1,0,0,1)   ;
      TLorentzVector genBbarLV(1,0,0,1);

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
  
      TLorentzVector TOPLEP(1,0,0,1);
      TLorentzVector TOPHAD(1,0,0,1);
      TLorentzVector TOPHADW1(1,0,0,1);
      TLorentzVector TOPHADW2(1,0,0,1);
      TLorentzVector TOPHADB(1,0,0,1);
      TLorentzVector TOPLEPW1(1,0,0,1);
      TLorentzVector TOPLEPW2(1,0,0,1);
      TLorentzVector TOPLEPB(1,0,0,1);
      TLorentzVector HIGGS(1,0,0,1);

      float neutCosTheta;
      int chargeTopLep = 0;
      bool isSL = false;
      bool isDL = false;
      int flavHADW1 = 0;
      int flavHADW2 = 0;

      // dummy
      bool properEventSL = (genBLV.Pt()>0 && genBbarLV.Pt()>0 && topBLV.Pt()>0 && topW1LV.Pt()>0 && topW2LV.Pt()>0 && atopBLV.Pt()>0 && atopW1LV.Pt()>0 && atopW2LV.Pt()>0);
      bool properEventDL = (genBLV.Pt()>0 && genBbarLV.Pt()>0 && topBLV.Pt()>0 && topW1LV.Pt()>0 && topW2LV.Pt()>0 && atopBLV.Pt()>0 && atopW1LV.Pt()>0 && atopW2LV.Pt()>0);

      if(!(properEventSL || properEventDL)) continue;

      HIGGS.SetPxPyPzE( (genBLV+genBbarLV).Px(), (genBLV+genBbarLV).Py(),(genBLV+genBbarLV).Pz(),(genBLV+genBbarLV).E());
      if( abs(genTop.wdau1id)>6 && abs(genTbar.wdau1id)<6){
	TOPLEP.SetPxPyPzE( (topBLV+topW1LV+topW2LV).Px(), (topBLV+topW1LV+topW2LV).Py(), (topBLV+topW1LV+topW2LV).Pz(), (topBLV+topW1LV+topW2LV).E() );
	TOPHAD.SetPxPyPzE( (atopBLV+atopW1LV+atopW2LV).Px(), (atopBLV+atopW1LV+atopW2LV).Py(), (atopBLV+atopW1LV+atopW2LV).Pz(), (atopBLV+atopW1LV+atopW2LV).E() );
	TOPHADW1.SetPxPyPzE( atopW1LV.Px(), atopW1LV.Py(), atopW1LV.Pz(), atopW1LV.E());
	TOPHADW2.SetPxPyPzE( atopW2LV.Px(), atopW2LV.Py(), atopW2LV.Pz(), atopW2LV.E());
	flavHADW1 = genTbar.wdau1id;
	flavHADW2 = genTbar.wdau2id;
	TOPHADB.SetPxPyPzE( atopBLV.Px(),  atopBLV.Py(),   atopBLV.Pz(),  atopBLV.E());
	if( abs(genTop.wdau1id)==11 || abs(genTop.wdau1id)==13 || abs(genTop.wdau1id)==15 ){
	  chargeTopLep = genTop.wdau1id>0 ? 1 : -1;
	  TOPLEPW1.SetPxPyPzE  ( topW1LV.Px(), topW1LV.Py(), topW1LV.Pz(), topW1LV.E());
	  TOPLEPW2.SetPtEtaPhiM( topW2LV.Pt(), 0.0, topW2LV.Phi(), 0.0);
	  neutCosTheta = TMath::Cos( topW2LV.Vect().Theta() );
	}
	else{
	  chargeTopLep = genTop.wdau2id>0 ? 1 : -1;
	  TOPLEPW2.SetPtEtaPhiM( topW1LV.Pt(), 0.0, topW1LV.Phi(), 0.0);
	  TOPLEPW1.SetPxPyPzE  ( topW2LV.Px(), topW2LV.Py(), topW2LV.Pz(), topW2LV.E());
	  neutCosTheta = TMath::Cos( topW1LV.Vect().Theta() );
	}
	TOPLEPB.SetPxPyPzE(  topBLV.Px(),  topBLV.Py(),   topBLV.Pz(), topBLV.E());
	isSL = true;
      }
      else if(abs(genTop.wdau1id)<6 && abs(genTbar.wdau1id)>6){
	TOPHAD.SetPxPyPzE( (topBLV+topW1LV+topW2LV).Px(), (topBLV+topW1LV+topW2LV).Py(), (topBLV+topW1LV+topW2LV).Pz(), (topBLV+topW1LV+topW2LV).E() );
	TOPLEP.SetPxPyPzE( (atopBLV+atopW1LV+atopW2LV).Px(), (atopBLV+atopW1LV+atopW2LV).Py(), (atopBLV+atopW1LV+atopW2LV).Pz(), (atopBLV+atopW1LV+atopW2LV).E() );
	TOPHADW1.SetPxPyPzE( topW1LV.Px(), topW1LV.Py(), topW1LV.Pz(), topW1LV.E());
	TOPHADW2.SetPxPyPzE( topW2LV.Px(), topW2LV.Py(), topW2LV.Pz(), topW2LV.E());
	flavHADW1 = genTop.wdau1id;
	flavHADW2 = genTop.wdau2id;
	TOPHADB.SetPxPyPzE( topBLV.Px(),  topBLV.Py(),   topBLV.Pz(),  topBLV.E());
	if( abs(genTbar.wdau1id)==11 || abs(genTbar.wdau1id)==13 || abs(genTbar.wdau1id)==15 ){
	  chargeTopLep = genTbar.wdau1id>0 ? 1 : -1;	  
	  TOPLEPW1.SetPxPyPzE  ( atopW1LV.Px(), atopW1LV.Py(), atopW1LV.Pz(), atopW1LV.E());
	  TOPLEPW2.SetPtEtaPhiM( atopW2LV.Pt(), 0.0, atopW2LV.Phi(), 0.0);
	  neutCosTheta = TMath::Cos( atopW2LV.Vect().Theta() );
	}
	else{
	  chargeTopLep = genTbar.wdau2id>0 ? 1 : -1;	  
	  TOPLEPW2.SetPtEtaPhiM( atopW1LV.Pt(), 0.0, atopW1LV.Phi(), 0.0);
	  TOPLEPW1.SetPxPyPzE  ( atopW2LV.Px(), atopW2LV.Py(), atopW2LV.Pz(), atopW2LV.E());
	  neutCosTheta = TMath::Cos( atopW1LV.Vect().Theta() );
	}
	TOPLEPB.SetPxPyPzE( atopBLV.Px(),  atopBLV.Py(),   atopBLV.Pz(),  atopBLV.E());
	isSL = true;
      }      
      else if(abs(genTop.wdau1id)>6 && abs(genTbar.wdau1id)>6){
	TOPHAD.SetPxPyPzE( (topBLV+topW1LV+topW2LV).Px(), (topBLV+topW1LV+topW2LV).Py(), (topBLV+topW1LV+topW2LV).Pz(), (topBLV+topW1LV+topW2LV).E() );
	TOPLEP.SetPxPyPzE( (atopBLV+atopW1LV+atopW2LV).Px(), (atopBLV+atopW1LV+atopW2LV).Py(), (atopBLV+atopW1LV+atopW2LV).Pz(), (atopBLV+atopW1LV+atopW2LV).E() );
	if( abs(genTop.wdau1id)==11 || abs(genTop.wdau1id)==13 || abs(genTop.wdau1id)==15 ){
	  TOPHADW1.SetPxPyPzE  ( topW1LV.Px(), topW1LV.Py(), topW1LV.Pz(), topW1LV.E());
	  TOPHADW2.SetPtEtaPhiM( topW2LV.Pt(), 0.0, topW2LV.Phi(), 0.0);
	}
	else{
	  TOPHADW2.SetPtEtaPhiM( topW1LV.Pt(), 0.0, topW1LV.Phi(), 0.0);
	  TOPHADW1.SetPxPyPzE  ( topW2LV.Px(), topW2LV.Py(), topW2LV.Pz(), topW2LV.E());
	}
	TOPHADB.SetPxPyPzE( topBLV.Px(),  topBLV.Py(),   topBLV.Pz(),  topBLV.E());
	if( abs(genTbar.wdau1id)==11 || abs(genTbar.wdau1id)==13 || abs(genTbar.wdau1id)==15 ){
	  chargeTopLep = genTbar.wdau1id>0 ? 1 : -1;	  
	  TOPLEPW1.SetPxPyPzE  ( atopW1LV.Px(), atopW1LV.Py(), atopW1LV.Pz(), atopW1LV.E());
	  TOPLEPW2.SetPtEtaPhiM( atopW2LV.Pt(), 0.0, atopW2LV.Phi(), 0.0);
	  neutCosTheta = TMath::Cos( atopW2LV.Vect().Theta() );
	}
	else{
	  chargeTopLep = genTbar.wdau2id>0 ? 1 : -1;	  
	  TOPLEPW2.SetPtEtaPhiM( atopW1LV.Pt(), 0.0, atopW1LV.Phi(), 0.0);
	  TOPLEPW1.SetPxPyPzE  ( atopW2LV.Px(), atopW2LV.Py(), atopW2LV.Pz(), atopW2LV.E());
	  neutCosTheta = TMath::Cos( atopW1LV.Vect().Theta() );
	}
	TOPLEPB.SetPxPyPzE( atopBLV.Px(),  atopBLV.Py(),   atopBLV.Pz(),  atopBLV.E());
	isDL = true;
      }      
      else{
	//properEventSL=false;
	//properEventDL=false;
      }

      if(isDL) TOPLEPW2 = TOPLEPW2+TOPHADW2;

      /////////////////////////////////////////////////////

      TLorentzVector leptonLV, leptonLV2;
      leptonLV.SetPtEtaPhiM(vLepton_pt[0],vLepton_eta[0],vLepton_phi[0],vLepton_mass[0]);
      properEventSL = 
	(vLepton_type[0]==13 && leptonLV.Pt()>30 && TMath::Abs(leptonLV.Eta())<2.1 && vLepton_pfCorrIso[0]<0.10) ||
	(vLepton_type[0]==11 && leptonLV.Pt()>30 && TMath::Abs(leptonLV.Eta())<2.5 && !(TMath::Abs(leptonLV.Eta())>1.442 &&  TMath::Abs(leptonLV.Eta())<1.566) && vLepton_pfCorrIso[0]<0.10) ;
      

      properEventDL = false;
      if( nvlep>=2 ){
	leptonLV.SetPtEtaPhiM (vLepton_pt[0],vLepton_eta[0],vLepton_phi[0],vLepton_mass[0]);
	leptonLV2.SetPtEtaPhiM(vLepton_pt[1],vLepton_eta[1],vLepton_phi[1],vLepton_mass[1]);

	properEventDL = (
			 ( vLepton_type[0]==13 && vLepton_type[1]==13 && 
			   ( (leptonLV.Pt() >20 && TMath::Abs(leptonLV.Eta()) <2.1 && vLepton_pfCorrIso[0]<0.10 &&
			      leptonLV2.Pt()>10 && TMath::Abs(leptonLV2.Eta())<2.4 && vLepton_pfCorrIso[1]<0.20) ||
			     (leptonLV.Pt() >10 && TMath::Abs(leptonLV.Eta()) <2.4 && vLepton_pfCorrIso[0]<0.20 &&
			      leptonLV2.Pt()>20 && TMath::Abs(leptonLV2.Eta())<2.1 && vLepton_pfCorrIso[1]<0.10) )
			   ) ||
			 ( vLepton_type[0]==11 && vLepton_type[0]==11 && 
			   ( (leptonLV.Pt() >20 && TMath::Abs(leptonLV.Eta()) <2.5 && vLepton_pfCorrIso[0]<0.10 &&
			      leptonLV2.Pt()>10 && TMath::Abs(leptonLV2.Eta())<2.5 && vLepton_pfCorrIso[1]<0.20) ||
			     (leptonLV.Pt() >10 && TMath::Abs(leptonLV.Eta()) <2.5 && vLepton_pfCorrIso[0]<0.20 &&
			      leptonLV2.Pt()>20 && TMath::Abs(leptonLV2.Eta())<2.5 && vLepton_pfCorrIso[1]<0.10) )
			   )
			 )
	  &&  vLepton_charge[0]*vLepton_charge[1]<0 ;
      }

      float nuPx = METtype1p2corr.et*TMath::Cos(METtype1p2corr.phi);
      float nuPy = METtype1p2corr.et*TMath::Sin(METtype1p2corr.phi);
      float nuE  = TMath::Sqrt(nuPx*nuPx+nuPy*nuPy);
      TLorentzVector neutrinoLV;
      neutrinoLV.SetPxPyPzE(nuPx,nuPy,0. ,nuE);

      /////////////////////////////////////////////////////

      vector<TLorentzVector> myJets;
      std::map<unsigned int, JetByPt> map;
      vector<TLorentzVector> myJetsFilt;
      std::map<unsigned int, JetByPt> mapFilt;
      
      map[0] = jet1;  map[1] = jet2; map[2] = jet3; map[3] = jet4; 
      map[4] = jet5;  map[5] = jet6; map[6] = jet7; map[7] = jet8; 
      map[8] = jet9;  map[9] = jet10; 
    

      for(unsigned int m = 0; m <10; m++){
	TLorentzVector jetLV; 

	int index_m     =  map[m].index;
	float JECUnc_m  =  index_m>=0 ? hJet_JECUnc[index_m] : aJet_JECUnc[-index_m-1];
	float shift     = 1.0;
	if     ( doJECup   )  shift *= (1+JECUnc_m);
	else if( doJECdown )  shift *= (1-JECUnc_m);
	else{}

	if( (map[m].pt*shift) > 20 ){
	  jetLV.SetPtEtaPhiM( map[m].pt*shift, map[m].eta, map[m].phi, map[m].mass*shift );
	  myJets.push_back(jetLV);
	}

      }



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
      


      int numJets20UntagM =0;     
      int numJets30UntagM =0; 
      int numJets60UntagM=0; 
      int numJets20UntagL =0;     
      int numJets30UntagL =0; 
      int numJets60UntagL=0; 
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

	int index_k =  mapFilt[k].index;
	
	float csv_k_nominal = index_k>=0 ? hJet_csv_nominal[index_k]  : aJet_csv_nominal[-index_k-1];
	float csv_k_upBC    = index_k>=0 ? hJet_csv_upBC[index_k]     : aJet_csv_upBC[-index_k-1];
	float csv_k_downBC  = index_k>=0 ? hJet_csv_downBC[index_k]   : aJet_csv_downBC[-index_k-1];
	float csv_k_upL     = index_k>=0 ? hJet_csv_upL[index_k]      : aJet_csv_upL[-index_k-1];
	float csv_k_downL   = index_k>=0 ? hJet_csv_downL[index_k]    : aJet_csv_downL[-index_k-1];
	float JECUnc_k      = index_k>=0 ? hJet_JECUnc[index_k]       : aJet_JECUnc[index_k];

	//float csv_k = mapFilt[k].csv>0 ? mapFilt[k].csv : 0.0 ;

	int btag_L = csv_k_nominal>0.244 ;	
	if     ( doCSVup  ) btag_L = ( TMath::Max(csv_k_upBC,   csv_k_upL)  >0.244 );
	else if( doCSVdown) btag_L = ( TMath::Min(csv_k_downBC, csv_k_downL)>0.244 );
	else{}

	int btag_M = csv_k_nominal>0.679 ;
	if     ( doCSVup  ) btag_M = ( TMath::Max(csv_k_upBC,   csv_k_upL)  >0.679 );
	else if( doCSVdown) btag_M = ( TMath::Min(csv_k_downBC, csv_k_downL)>0.679 );
	else{}

	int btag_T = csv_k_nominal>0.898 ;	
	if     ( doCSVup  ) btag_T = ( TMath::Max(csv_k_upBC,   csv_k_upL)  >0.898 );
	else if( doCSVdown) btag_T = ( TMath::Min(csv_k_downBC, csv_k_downL)>0.898 );
	else{}

	float pt_k  = myJetsFilt[k].Pt();

	if( pt_k>20 && btag_L ) numJets20BtagL++;
	if( pt_k>20 && btag_M ) numJets20BtagM++;
	if( pt_k>20 && btag_T ) numJets20BtagT++;
	if( pt_k>30 && btag_L ) numJets30BtagL++;
	if( pt_k>30 && btag_M ) numJets30BtagM++;
	if( pt_k>30 && btag_T ) numJets30BtagT++;
	if( pt_k>60 && btag_L ) numJets60BtagL++;
	if( pt_k>60 && btag_M ) numJets60BtagM++;
	if( pt_k>60 && btag_T ) numJets60BtagT++;
	if( pt_k>20 && !btag_M ) numJets20UntagM++;
	if( pt_k>30 && !btag_M ) numJets30UntagM++;
	if( pt_k>60 && !btag_M ) numJets60UntagM++; 
	if( pt_k>20 && !btag_L ) numJets20UntagL++;
	if( pt_k>30 && !btag_L ) numJets30UntagL++;
	if( pt_k>60 && !btag_L ) numJets60UntagL++;      
      }


      ////////////////////////////////////////////////////

      //cout << "####################" << endl;
      //for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
      //cout << myJetsFilt[k].Pt() << " vs " << mapFilt[k].pt << endl;
      //}

      bool atLeastFourJets = false;
      if( myJetsFilt.size()>=4 ){
	if( myJetsFilt[0].Pt()>40 && myJetsFilt[1].Pt()>40 && myJetsFilt[2].Pt()>40 && myJetsFilt[3].Pt()>40 )
	  atLeastFourJets = true;
      }

      /////////////////////////////////////////////////////

      if( !(properEventSL || properEventDL) || !atLeastFourJets ) continue;

      probAtSgn_          = 0.;
      probAtSgn_alt_      = 0.;
      probAtSgn_ttbb_     = 0.;
      probAtSgn_alt_ttbb_ = 0.;
      probAtSgn_alt_ttjj_ = 0.;


      matchesH_   = -99;
      matchesW_   = -99;
      matchesT_   = -99;
      matchesHAll_= -99;
      matchesWAll_= -99;
      matchesTAll_= -99;
      overlapHeavy_=-99;
      overlapLight_=-99;
      type_        = -99;
      csvL_[0]     = -99; 
      csvL_[1]     = -99;
      ptL_[0]      = -99; 
      ptL_[1]      = -99;
      mL_[0]       = -99; 
      mL_[1]       = -99;
      csvH_[0]     = -99; 
      csvH_[1]     = -99; 
      csvH_[2]     = -99; 
      csvH_[3]     = -99;  

      flag_type0_ = -99;
      flag_type1_ = -99;
      flag_type2_ = -99;
      flag_type3_ = -99;
      flag_type4_ = -99;
      flag_type6_ = -99;
      nSimBs_     = -99;
      nMatchSimBs_= -99;
      nC_         = -99;
      nCTop_      = -99;

      time_       = -99;

      dR_=-99;dPhi_=-99;dEta_=-99;Mjj_=-99;Pt_=-99;Eta_=-99;pt1_=-99;pt2_=-99;eta1_=-99;eta2_=-99;csv1_=-99;csv2_=-99;

      /////////////////////////////////////////////////////

      //which light jets?
      unsigned int w1=999;
      unsigned int w2=999;
      unsigned int w3=999;
      unsigned int w4=999;
      unsigned int w5=999;
      unsigned int b1=999;
      unsigned int b2=999;
      unsigned int b3=999;
      unsigned int b4=999;

      unsigned int wL1=999;
      unsigned int wL2=999;
      unsigned int wL3=999;
      unsigned int bL1=999;
      unsigned int bL2=999;
      unsigned int bL3=999;
      unsigned int bL4=999;
      
      for(unsigned int k = 0; k < myJetsFilt.size(); k++){  

	//float csv_k = mapFilt[k].csv>0 ? mapFilt[k].csv : 0.0 ;

	int index_k =  mapFilt[k].index;
	
	float csv_k_nominal = index_k>=0 ? hJet_csv_nominal[index_k]  : aJet_csv_nominal[-index_k-1];
	float csv_k_upBC    = index_k>=0 ? hJet_csv_upBC[index_k]     : aJet_csv_upBC[-index_k-1];
	float csv_k_downBC  = index_k>=0 ? hJet_csv_downBC[index_k]   : aJet_csv_downBC[-index_k-1];
	float csv_k_upL     = index_k>=0 ? hJet_csv_upL[index_k]      : aJet_csv_upL[-index_k-1];
	float csv_k_downL   = index_k>=0 ? hJet_csv_downL[index_k]    : aJet_csv_downL[-index_k-1];
	float JECUnc_k      = index_k>=0 ? hJet_JECUnc[index_k]       : aJet_JECUnc[index_k];

	//float csv_k = mapFilt[k].csv>0 ? mapFilt[k].csv : 0.0 ;

	int btag_L = csv_k_nominal>0.244 ;	
	if     ( doCSVup  ) btag_L = ( TMath::Max(csv_k_upBC,   csv_k_upL)  >0.244 );
	else if( doCSVdown) btag_L = ( TMath::Min(csv_k_downBC, csv_k_downL)>0.244 );
	else{}

	int btag_M = csv_k_nominal>0.679 ;
	if     ( doCSVup  ) btag_M = ( TMath::Max(csv_k_upBC,   csv_k_upL)  >0.679 );
	else if( doCSVdown) btag_M = ( TMath::Min(csv_k_downBC, csv_k_downL)>0.679 );
	else{}

	int btag_T = csv_k_nominal>0.898 ;	
	if     ( doCSVup  ) btag_T = ( TMath::Max(csv_k_upBC,   csv_k_upL)  >0.898 );
	else if( doCSVdown) btag_T = ( TMath::Min(csv_k_downBC, csv_k_downL)>0.898 );
	else{}

	if     (  btag_M && b1==999 && b2==999 && b3==999 && b4==999) b1 = k;	    
	else if(  btag_M && b1!=999 && b2==999 && b3==999 && b4==999) b2 = k;	    
	else if(  btag_M && b1!=999 && b2!=999 && b3==999 && b4==999) b3 = k;	    
	else if(  btag_M && b1!=999 && b2!=999 && b3!=999 && b4==999) b4 = k;	    
	else if( !btag_M && w1==999 && w2==999 && w3==999 && w4==999 && w5==999) w1 = k;
	else if( !btag_M && w1!=999 && w2==999 && w3==999 && w4==999 && w5==999) w2 = k; 
	else if( !btag_M && w1!=999 && w2!=999 && w3==999 && w4==999 && w5==999) w3 = k; 
	else if( !btag_M && w1!=999 && w2!=999 && w3!=999 && w4==999 && w5==999) w4 = k; 
	else if( !btag_M && w1!=999 && w2!=999 && w3!=999 && w4!=999 && w5==999) w5 = k; 
	else{}
      }


      for(unsigned int k = 0; k < myJetsFilt.size(); k++){  

	//float csv_k = mapFilt[k].csv>0 ? mapFilt[k].csv : 0.0 ;

	int index_k =  mapFilt[k].index;
	
	float csv_k_nominal = index_k>=0 ? hJet_csv_nominal[index_k]  : aJet_csv_nominal[-index_k-1];
	float csv_k_upBC    = index_k>=0 ? hJet_csv_upBC[index_k]     : aJet_csv_upBC[-index_k-1];
	float csv_k_downBC  = index_k>=0 ? hJet_csv_downBC[index_k]   : aJet_csv_downBC[-index_k-1];
	float csv_k_upL     = index_k>=0 ? hJet_csv_upL[index_k]      : aJet_csv_upL[-index_k-1];
	float csv_k_downL   = index_k>=0 ? hJet_csv_downL[index_k]    : aJet_csv_downL[-index_k-1];
	float JECUnc_k      = index_k>=0 ? hJet_JECUnc[index_k]       : aJet_JECUnc[index_k];
	
	int btag_L = csv_k_nominal>0.244 ;	
	if     ( doCSVup  ) btag_L = ( TMath::Max(csv_k_upBC,   csv_k_upL)  >0.244 );
	else if( doCSVdown) btag_L = ( TMath::Min(csv_k_downBC, csv_k_downL)>0.244 );
	else{}

	int btag_M = csv_k_nominal>0.679 ;
	if     ( doCSVup  ) btag_M = ( TMath::Max(csv_k_upBC,   csv_k_upL)  >0.679 );
	else if( doCSVdown) btag_M = ( TMath::Min(csv_k_downBC, csv_k_downL)>0.679 );
	else{}

	int btag_T = csv_k_nominal>0.898 ;	
	if     ( doCSVup  ) btag_T = ( TMath::Max(csv_k_upBC,   csv_k_upL)  >0.898 );
	else if( doCSVdown) btag_T = ( TMath::Min(csv_k_downBC, csv_k_downL)>0.898 );
	else{}

	if     (  btag_L && bL1==999 && bL2==999 && bL3==999 && bL4==999) bL1 = k;	    
	else if(  btag_L && bL1!=999 && bL2==999 && bL3==999 && bL4==999) bL2 = k;	    
	else if(  btag_L && bL1!=999 && bL2!=999 && bL3==999 && bL4==999) bL3 = k;	    
	else if(  btag_L && bL1!=999 && bL2!=999 && bL3!=999 && bL4==999) bL4 = k;	    
	else if( !btag_L && wL1==999 && wL2==999) wL1 = k;
	else if( !btag_L && wL1!=999 && wL2==999) wL2 = k; 
	else if( !btag_L && wL1!=999 && wL2!=999) wL3 = k; 
	else{}
      }


      vector<TLorentzVector> jets;
      std::map< unsigned int, unsigned int> pos_to_index;


      if(properEventSL && numJets30BtagM==4 && numJets30UntagM>=2 && (doType0 || doType1 || doType3)){

	if(numJets30UntagM==2 && (doType0 || doType1)){

	  counter++;      
	  if( fixNumberOfEventsPerJob && !(counter>=evLow && counter<=evHigh) ) continue;
	  cout << "Processing event # " << counter << " (SL[4,2])" << endl;
	 
	  if( w1==999 || w2==999 || b1==999 || b2==999 || b3==999 || b4==999){
	    cout << "Inconsistentcy type 0/1 found!" << endl;
	    cout << w1 << ", " << w2 << "," << b1 << "," << b2 << "," << b3 << "," << b4 <<  endl;
	    continue;
	  }

	  //cout << "####################" << endl;
	  //for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
	  //int index_k =  mapFilt[k].index;
	  //float csv_k_nominal = index_k>=0 ? hJet_csv_nominal[index_k]  : aJet_csv_nominal[-index_k-1];
	  //cout << k << ":  Pt = " << myJetsFilt[k].Pt() << ", csv =  " << csv_k_nominal << endl;	   
	  //}
	  //cout << w1 << "," << w2 << "," << b1 << "," << b2 << "," << b3 << "," << b4 <<  endl;
	  //cout << "####################" << endl;


	  int hMatches = 0;
	  if     (  deltaR(myJetsFilt[b1],genBLV)<0.3 )    hMatches++;
	  else if(  deltaR(myJetsFilt[b1],genBbarLV)<0.3)  hMatches++;
	  if     (  deltaR(myJetsFilt[b2],genBLV)<0.3 )    hMatches++;
	  else if(  deltaR(myJetsFilt[b2],genBbarLV)<0.3)  hMatches++;
	  if     (  deltaR(myJetsFilt[b3],genBLV)<0.3 )    hMatches++;
	  else if(  deltaR(myJetsFilt[b3],genBbarLV)<0.3)  hMatches++;
	  if     (  deltaR(myJetsFilt[b4],genBLV)<0.3 )    hMatches++;
	  else if(  deltaR(myJetsFilt[b4],genBbarLV)<0.3)  hMatches++;
	  matchesH_    = hMatches;
	  matchesHAll_ = hMatches;
	  if     (  deltaR(myJetsFilt[w1],genBLV)<0.3 )    matchesHAll_++;
	  else if(  deltaR(myJetsFilt[w1],genBbarLV)<0.3)  matchesHAll_++;
	  if     (  deltaR(myJetsFilt[w2],genBLV)<0.3 )    matchesHAll_++;
	  else if(  deltaR(myJetsFilt[w2],genBbarLV)<0.3)  matchesHAll_++;

	  int tMatches = 0;
	  if     (  deltaR(myJetsFilt[b1],TOPHADB)<0.3 )   tMatches++;
	  else if(  deltaR(myJetsFilt[b1],TOPLEPB)<0.3)    tMatches++;
	  if     (  deltaR(myJetsFilt[b2],TOPHADB)<0.3 )   tMatches++;
	  else if(  deltaR(myJetsFilt[b2],TOPLEPB)<0.3)    tMatches++;
	  if     (  deltaR(myJetsFilt[b3],TOPHADB)<0.3 )   tMatches++;
	  else if(  deltaR(myJetsFilt[b3],TOPLEPB)<0.3)    tMatches++;
	  if     (  deltaR(myJetsFilt[b4],TOPHADB)<0.3 )   tMatches++;
	  else if(  deltaR(myJetsFilt[b4],TOPLEPB)<0.3)    tMatches++;
	  matchesT_    = tMatches;
	  matchesTAll_ = tMatches;
	  if     (  deltaR(myJetsFilt[w1],TOPHADB)<0.3 )   matchesTAll_++;
	  else if(  deltaR(myJetsFilt[w1],TOPLEPB)<0.3)    matchesTAll_++;
	  if     (  deltaR(myJetsFilt[w2],TOPHADB)<0.3 )   matchesTAll_++;
	  else if(  deltaR(myJetsFilt[w2],TOPLEPB)<0.3)    matchesTAll_++;


	  int wMatches = 0;
	  if     (  deltaR(myJetsFilt[w1],TOPHADW1)<0.3 ){
	    wMatches++;
	    mL_[0] = 1;
	  }
	  else if(  deltaR(myJetsFilt[w1],TOPHADW2)<0.3){
	    wMatches++;
	    mL_[0] = 1;
	  }
	  else{ mL_[0] = 0; }
	  if     (  deltaR(myJetsFilt[w2],TOPHADW1)<0.3 ){
	    wMatches++;
	    mL_[1] = 1;
	  }
	  else if(  deltaR(myJetsFilt[w2],TOPHADW2)<0.3){
	    wMatches++;
	    mL_[1] = 1;
	  }
	  else{  mL_[1] = 0; }
	  matchesW_    = wMatches;
	  matchesWAll_ = wMatches;	  
	  if     (  deltaR(myJetsFilt[b1],TOPHADW1)<0.3 )   matchesWAll_++;
	  else if(  deltaR(myJetsFilt[b1],TOPHADW2)<0.3)    matchesWAll_++;
	  if     (  deltaR(myJetsFilt[b2],TOPHADW1)<0.3 )   matchesWAll_++;
	  else if(  deltaR(myJetsFilt[b2],TOPHADW2)<0.3)    matchesWAll_++;
	  if     (  deltaR(myJetsFilt[b3],TOPHADW1)<0.3 )   matchesWAll_++;
	  else if(  deltaR(myJetsFilt[b3],TOPHADW2)<0.3)    matchesWAll_++;
	  if     (  deltaR(myJetsFilt[b4],TOPHADW1)<0.3 )   matchesWAll_++;
	  else if(  deltaR(myJetsFilt[b4],TOPHADW2)<0.3)    matchesWAll_++;

	  vector<TLorentzVector> genHeavy;
	  genHeavy.push_back( genBLV);
	  genHeavy.push_back( genBbarLV);
	  genHeavy.push_back( TOPHADB);
	  genHeavy.push_back( TOPLEPB);
	  int overlapH = 0;
	  for(unsigned int k = 0; k < genHeavy.size()-1; k++){  
	    for(unsigned int l = k+1; l < genHeavy.size(); l++){  
	      if( deltaR(genHeavy[k], genHeavy[l])<0.5 ) overlapH++;
	    }	    
	  }

	  overlapHeavy_ = overlapH;
	  vector<TLorentzVector> genLight;
	  genLight.push_back( TOPHADW1);
	  genLight.push_back( TOPHADW2);
	  int overlapL = 0;
	  for(unsigned int k = 0; k < genLight.size(); k++){  
	    for(unsigned int l = 0; l < genHeavy.size(); l++){  
	      if( deltaR(genLight[k], genHeavy[l])<0.5 ) overlapL++;
	    }	    
	  }
	  if(  deltaR(genLight[0], genLight[1])<0.5  ) overlapL++;
	  overlapLight_  = overlapL;

	  csvL_[0] = (mapFilt[w1].csv>0 ? mapFilt[w1].csv : 0.0 );
	  csvL_[1] = (mapFilt[w2].csv>0 ? mapFilt[w2].csv : 0.0 );
	  ptL_[0]  = (mapFilt[w1].pt>0  ? mapFilt[w1].pt  : 0.0 );
	  ptL_[1]  = (mapFilt[w2].pt>0  ? mapFilt[w2].pt  : 0.0 );

	  csvH_[0] = (mapFilt[b1].csv>0 ? mapFilt[b1].csv : 0.0 );
	  csvH_[1] = (mapFilt[b2].csv>0 ? mapFilt[b2].csv : 0.0 );
	  csvH_[2] = (mapFilt[b3].csv>0 ? mapFilt[b3].csv : 0.0 );
	  csvH_[3] = (mapFilt[b4].csv>0 ? mapFilt[b4].csv : 0.0 );


	  float WMass = (myJetsFilt[w1]+myJetsFilt[w2]).M();

	  jets.clear();
	  jets.push_back( leptonLV  );  
	  jets.push_back( neutrinoLV  );  
	  jets.push_back( myJetsFilt[b1]   );
	  jets.push_back( myJetsFilt[w1]   );  
	  jets.push_back( myJetsFilt[w2]   );  
	  jets.push_back( myJetsFilt[b2]   );
	  jets.push_back( myJetsFilt[b3]   );
	  jets.push_back( myJetsFilt[b4]   );
	  
	  pos_to_index.clear();
	  pos_to_index[2] = b1;
	  pos_to_index[3] = w1;
	  pos_to_index[4] = w2;
	  pos_to_index[5] = b2;
	  pos_to_index[6] = b3;
	  pos_to_index[7] = b4;

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

	    //if( Bs.Pt()>20 && TMath::Abs(Bs.Eta())<2.5 ) nMatchSimBs_++;
	    //if     ( deltaR( Bs , myJetsFilt[b1] )<0.2 ) nMatchSimBs_++;
	    //else if( deltaR( Bs , myJetsFilt[b2] )<0.2 ) nMatchSimBs_++;
	    //else if( deltaR( Bs , myJetsFilt[b3] )<0.2 ) nMatchSimBs_++;
	    //else if( deltaR( Bs , myJetsFilt[b4] )<0.2 ) nMatchSimBs_++;
	    //else if( deltaR( Bs , myJetsFilt[w1] )<0.2 ) nMatchSimBs_++;
	    //else if( deltaR( Bs , myJetsFilt[w2] )<0.2 ) nMatchSimBs_++;
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

	  meIntegrator->setJets(&jets);

	  if( WMass>MwL && WMass<MwH ){
	    cout << "W mass = " << WMass << " ... OK: SL(4,2)" << endl;
	    type_ = 0;
	    meIntegrator->setIntType( MEIntegratorNew::SL2wj );	    
	  }
	  else{
	    cout << "W mass = " << WMass << " ... NO: SL(4,1)" << endl;
	    type_ = 1;
	    meIntegrator->setIntType( MEIntegratorNew::SL1wj );
	  }
	  
	  /////////////////////////////////////////////////////////////
	  int index_bLep          =  mapFilt[b1].index;
	  float csv_bLep_nominal  =  index_bLep>=0 ? hJet_csv_nominal[index_bLep]  : aJet_csv_nominal[-index_bLep-1];
	  float csv_bLep_upBC     =  index_bLep>=0 ? hJet_csv_upBC   [index_bLep]  : aJet_csv_upBC   [-index_bLep-1];
	  float csv_bLep_downBC   =  index_bLep>=0 ? hJet_csv_downBC [index_bLep]  : aJet_csv_downBC [-index_bLep-1];
	  float csv_bLep_upL      =  index_bLep>=0 ? hJet_csv_upL    [index_bLep]  : aJet_csv_upL    [-index_bLep-1];
	  float csv_bLep_downL    =  index_bLep>=0 ? hJet_csv_downL  [index_bLep]  : aJet_csv_downL  [-index_bLep-1];
	  float csv_bLep = csv_bLep_nominal;
	  if( doCSVup )   csv_bLep = TMath::Max(csv_bLep_upBC,   csv_bLep_upL);
	  if( doCSVdown ) csv_bLep = TMath::Min(csv_bLep_downBC, csv_bLep_downL);   

	  if     ( csv_bLep< 1.1 ){
	    float m1    = (meIntegrator->jetAt(2)+meIntegrator->jetAt(3)).M();
	    float m2    = (meIntegrator->jetAt(2)+meIntegrator->jetAt(4)).M();
	    if( ((m1>(MwL+5) && m1<(MwH-5)) || (m2>(MwL+5) && m2<(MwH-5))) && type_== 0 ){
	      flag_type0_ = 0; 
	      if( csv_bLep<0.95 ) flag_type0_ = 1; 
	      if( csv_bLep<0.90 ) flag_type0_ = 2; 
	      if( csv_bLep<0.85 ) flag_type0_ = 3;
	      if( csv_bLep<0.80 ) flag_type0_ = 4;  
	    }
	    if( ((m1>(MwL+5) && m1<(MwH-5)) || (m2>(MwL+5) && m2<(MwH-5))) && type_== 1 ){
	      flag_type1_ = 0; 
	      if( csv_bLep<0.95 ) flag_type1_ = 1; 
	      if( csv_bLep<0.90 ) flag_type1_ = 2; 
	      if( csv_bLep<0.85 ) flag_type1_ = 3;
	      if( csv_bLep<0.80 ) flag_type1_ = 4;  
	    }
	  }
	  int index_bHad          =  mapFilt[b2].index;
	  float csv_bHad_nominal  =  index_bHad>=0 ? hJet_csv_nominal[index_bHad]  : aJet_csv_nominal[-index_bHad-1];
	  float csv_bHad_upBC     =  index_bHad>=0 ? hJet_csv_upBC   [index_bHad]  : aJet_csv_upBC   [-index_bHad-1];
	  float csv_bHad_downBC   =  index_bHad>=0 ? hJet_csv_downBC [index_bHad]  : aJet_csv_downBC [-index_bHad-1];
	  float csv_bHad_upL      =  index_bHad>=0 ? hJet_csv_upL    [index_bHad]  : aJet_csv_upL    [-index_bHad-1];
	  float csv_bHad_downL    =  index_bHad>=0 ? hJet_csv_downL  [index_bHad]  : aJet_csv_downL  [-index_bHad-1];	  
	  float csv_bHad = csv_bHad_nominal;
	  if( doCSVup )   csv_bHad = TMath::Max(csv_bHad_upBC,   csv_bHad_upL);
	  if( doCSVdown ) csv_bHad = TMath::Min(csv_bHad_downBC, csv_bHad_downL);
	  if     ( csv_bHad< 1.1 ){
	    float m1    = (meIntegrator->jetAt(5)+meIntegrator->jetAt(3)).M();
	    float m2    = (meIntegrator->jetAt(5)+meIntegrator->jetAt(4)).M();
	    if( ((m1>(MwL+5) && m1<(MwH-5)) || (m2>(MwL+5) && m2<(MwH-5))) && type_== 0 ){
	      flag_type0_ = 0; 
	      if( csv_bHad<0.95 ) flag_type0_ = 1; 
	      if( csv_bHad<0.90 ) flag_type0_ = 2; 
	      if( csv_bHad<0.85 ) flag_type0_ = 3;
	      if( csv_bHad<0.80 ) flag_type0_ = 4;  
	    }
	    if( ((m1>(MwL+5) && m1<(MwH-5)) || (m2>(MwL+5) && m2<(MwH-5))) && type_== 1 ){
	      flag_type1_ = 0; 
	      if( csv_bHad<0.95 ) flag_type1_ = 1; 
	      if( csv_bHad<0.90 ) flag_type1_ = 2; 
	      if( csv_bHad<0.85 ) flag_type1_ = 3;
	      if( csv_bHad<0.80 ) flag_type1_ = 4;  
	    }
	  }
	  int index_b1          =  mapFilt[b3].index;
	  float csv_b1_nominal  =  index_b1>=0 ? hJet_csv_nominal[index_b1]  : aJet_csv_nominal[-index_b1-1];
	  float csv_b1_upBC     =  index_b1>=0 ? hJet_csv_upBC   [index_b1]  : aJet_csv_upBC   [-index_b1-1];
	  float csv_b1_downBC   =  index_b1>=0 ? hJet_csv_downBC [index_b1]  : aJet_csv_downBC [-index_b1-1];
	  float csv_b1_upL      =  index_b1>=0 ? hJet_csv_upL    [index_b1]  : aJet_csv_upL    [-index_b1-1];
	  float csv_b1_downL    =  index_b1>=0 ? hJet_csv_downL  [index_b1]  : aJet_csv_downL  [-index_b1-1];
	  float csv_b1 = csv_b1_nominal;
	  if( doCSVup )   csv_b1 = TMath::Max(csv_b1_upBC,   csv_b1_upL);
	  if( doCSVdown ) csv_b1 = TMath::Min(csv_b1_downBC, csv_b1_downL);	   
	  if     ( csv_b1< 1.1 ){
	    float m1    = (meIntegrator->jetAt(6)+meIntegrator->jetAt(3)).M();
	    float m2    = (meIntegrator->jetAt(6)+meIntegrator->jetAt(4)).M();
	    if( ((m1>(MwL+5) && m1<(MwH-5)) || (m2>(MwL+5) && m2<(MwH-5))) && type_== 0 ){
	      flag_type0_ = 0; 
	      if( csv_b1<0.95 ) flag_type0_ = 1; 
	      if( csv_b1<0.90 ) flag_type0_ = 2; 
	      if( csv_b1<0.85 ) flag_type0_ = 3;
	      if( csv_b1<0.80 ) flag_type0_ = 4; 
	    }
	    if( ((m1>(MwL+5) && m1<(MwH-5)) || (m2>(MwL+5) && m2<(MwH-5))) && type_== 1 ){
	      flag_type1_ = 0; 
	      if( csv_b1<0.95 ) flag_type1_ = 1; 
	      if( csv_b1<0.90 ) flag_type1_ = 2; 
	      if( csv_b1<0.85 ) flag_type1_ = 3;
	      if( csv_b1<0.80 ) flag_type1_ = 4;  	
	    }
	  }
	  int index_b2          =  mapFilt[b4].index;
	  float csv_b2_nominal  =  index_b2>=0 ? hJet_csv_nominal[index_b2]  : aJet_csv_nominal[-index_b2-1];
	  float csv_b2_upBC     =  index_b2>=0 ? hJet_csv_upBC   [index_b2]  : aJet_csv_upBC   [-index_b2-1];
	  float csv_b2_downBC   =  index_b2>=0 ? hJet_csv_downBC [index_b2]  : aJet_csv_downBC [-index_b2-1];
	  float csv_b2_upL      =  index_b2>=0 ? hJet_csv_upL    [index_b2]  : aJet_csv_upL    [-index_b2-1];
	  float csv_b2_downL    =  index_b2>=0 ? hJet_csv_downL  [index_b2]  : aJet_csv_downL  [-index_b2-1];
	  float csv_b2 = csv_b2_nominal;
	  if( doCSVup )   csv_b2 = TMath::Max(csv_b2_upBC,   csv_b2_upL);
	  if( doCSVdown ) csv_b2 = TMath::Min(csv_b2_downBC, csv_b2_downL);	  
	  if     ( csv_b2< 1.1 ){
	    float m1    = (meIntegrator->jetAt(7)+meIntegrator->jetAt(3)).M();
	    float m2    = (meIntegrator->jetAt(7)+meIntegrator->jetAt(4)).M();
	    if( ((m1>(MwL+5) && m1<(MwH-5)) || (m2>(MwL+5) && m2<(MwH-5))) && type_== 0 ){
	      flag_type0_ = 0; 
	      if( csv_b2<0.95 ) flag_type0_ = 1; 
	      if( csv_b2<0.90 ) flag_type0_ = 2; 
	      if( csv_b2<0.85 ) flag_type0_ = 3;
	      if( csv_b2<0.80 ) flag_type0_ = 4; 
	    }
	    if( ((m1>(MwL+5) && m1<(MwH-5)) || (m2>(MwL+5) && m2<(MwH-5))) && type_== 1 ){
	      flag_type1_ = 0; 
	      if( csv_b2<0.95 ) flag_type1_ = 1; 
	      if( csv_b2<0.90 ) flag_type1_ = 2; 
	      if( csv_b2<0.85 ) flag_type1_ = 3;
	      if( csv_b2<0.80 ) flag_type1_ = 4; 
	    }
	  }
	  /////////////////////////////////////////////////////////////


	  meIntegrator->setSumEt( METtype1p2corr.sumet );
	  meIntegrator->setMEtCov(-99,-99,0);
	  meIntegrator->setTopFlags( vLepton_charge[0]==1 ? +1 : -1 , vLepton_charge[0]==1 ? -1 : +1 );
	  
	  clock->Start();

	  for(int iter=0; iter<2; iter++){

	    if( type_ == 0 && iter==1) continue;
	    if( type_ == 1 && iter==1){
	      jets.clear();
	      jets.push_back( leptonLV  );  
	      jets.push_back( neutrinoLV  );  
	      jets.push_back( myJetsFilt[b1]   );
	      jets.push_back( myJetsFilt[w2]   );   // swap
	      jets.push_back( myJetsFilt[w1]   );  
	      jets.push_back( myJetsFilt[b2]   );
	      jets.push_back( myJetsFilt[b3]   );
	      jets.push_back( myJetsFilt[b4]   );

	      pos_to_index.clear();
	      pos_to_index[2] = b1;
	      pos_to_index[3] = w2;
	      pos_to_index[4] = w1;
	      pos_to_index[5] = b2;
	      pos_to_index[6] = b3;
	      pos_to_index[7] = b4;

	      meIntegrator->setJets(&jets);	
	    }
	    cout << "Doing iteration " << iter 
		 << ": W jet (pt,eta,phi) = (" << (iter==0 ? myJetsFilt[w1].Pt() :  myJetsFilt[w2].Pt() )  
		 << "," << (iter==0 ? myJetsFilt[w1].Eta() :  myJetsFilt[w2].Eta() ) 
		 << "," << (iter==0 ? myJetsFilt[w1].Phi() :  myJetsFilt[w2].Phi() ) 
		 << ")" << endl; 


	    for(int m = 0; m < nMassPoints ; m++){
	      meIntegrator->setMass( mH[m] );
	      
	      double maxP_s = 0.;
	      double maxP_b = 0.;

	      for(unsigned int pos = 0; pos < 12 ; pos++){

		meIntegrator->initVersors( permutations_SL2wj[pos] );
		double mass, massLow, massHigh;
		bool skip = !(meIntegrator->compatibilityCheck(0.95, /*printP4*/ 0, mass, massLow, massHigh )) ;
		bool skip_WHad   = !( meIntegrator->compatibilityCheck_WHad  (0.98, /*printP4*/ 0, mass, massLow, massHigh ) );
		bool skip_TopHad = !( meIntegrator->compatibilityCheck_TopHad(0.98, /*printP4*/ 0, mass, massLow, massHigh ) );
		if( type_==0 && (skip_WHad || skip_TopHad) ) continue;
		
		pair<double, double> range_x0 = (meIntegrator->getW1JetEnergyCI(0.95));
		pair<double, double> range_x1 =  make_pair(-1,1);
		pair<double, double> range_x2 =  make_pair(-PI,PI);	    
		pair<double, double> range_x3 =  make_pair(-1,1);
		pair<double, double> range_x4 = useMET ? (meIntegrator->getNuPhiCI(0.95)) : make_pair(-PI,PI);
		pair<double, double> range_x5 = (meIntegrator->getB1EnergyCI(0.95));
		pair<double, double> range_x6 = (meIntegrator->getB2EnergyCI(0.95));
		
		double x0L = range_x0.first;
		double x0U = range_x0.second;
		double x1L = range_x1.first;
		double x1U = range_x1.second;
		double x2L = range_x2.first;
		double x2U = range_x2.second;
		double x3L = range_x3.first;
		double x3U = range_x3.second;
		double x4L = range_x4.first;
		double x4U = range_x4.second;
		double x5L = range_x5.first;
		double x5U = range_x5.second;
		double x6L = range_x6.first;
		double x6U = range_x6.second;
	
		double xLmode0[4]   = {x0L, x3L, x4L, x5L};
		double xUmode0[4]   = {x0U, x3U, x4U, x5U};
		double xLmode0_b[5] = {x0L, x3L, x4L, x5L, x6L};
		double xUmode0_b[5] = {x0U, x3U, x4U, x5U, x6U};
	  	
		double xLmode1[6]   = {x0L, x1L, x2L, x3L, x4L, x5L};
		double xUmode1[6]   = {x0U, x1U, x2U, x3U, x4U, x5U};
		double xLmode1_b[7] = {x0L, x1L, x2L, x3L, x4L, x5L, x6L};
		double xUmode1_b[7] = {x0U, x1U, x2U, x3U, x4U, x5U, x6U};	       
		
		double p = 0.;

		if(speedup==0){
		  
		  for(int hyp = 0 ; hyp<2; hyp++){

		    if( SoB==0 && hyp!=hypo) continue;
		    if( SoB==1 && hyp==0 && skip) continue;

		    cout << "Testing hypothesis " << (hyp==0 ? "S" : "B") << " for permutation " << pos << endl; 
		    meIntegrator->setHypo(hyp);

		    int ntries = 0;
		    int intPoints = type_==0 ? 2000 : 4000;
		    while( ntries < MAX_REEVAL_TRIES){

		      ROOT::Math::Functor toIntegrate(meIntegrator, &MEIntegratorNew::Eval, (type_ == 0 ? 4+hyp : 6+hyp));
		      ROOT::Math::GSLMCIntegrator ig2( ROOT::Math::IntegrationMultiDim::kVEGAS , 1.e-12, 1.e-5, intPoints);
		      ig2.SetFunction(toIntegrate);
		      meIntegrator->SetPar( (type_ == 0 ? 4+hyp : 6+hyp) );	 

		      if( type_==0 ){
			p = (hyp==0 ? ig2.Integral(xLmode0, xUmode0) : ig2.Integral(xLmode0_b, xUmode0_b));
		      }
		      else{
			p = (hyp==0 ? ig2.Integral(xLmode1, xUmode1) : ig2.Integral(xLmode1_b, xUmode1_b));
		      }

		      double chi2 =  ig2.ChiSqr();
		      if( hyp==0 ){
			if( p>maxP_s ) maxP_s = p;		  
			else if( p<0.1*maxP_s) ntries = MAX_REEVAL_TRIES+1;
		      }
		      else{
			if( p>maxP_b ) maxP_b = p;		  
			else if( p<0.1*maxP_b) ntries = MAX_REEVAL_TRIES+1;	
		      }

		      if( chi2 > maxChi2_ ){
			ntries++;
			intPoints *= 1.5;
		      }
		      else ntries = MAX_REEVAL_TRIES+1;
		      
		    }		   
				    
		
		    if( TMath::IsNaN(p) ) p = 0.;


		    if( useBtag ){
		      
		      int bLep_pos    = (permutations_SL2wj[pos])/100000;
		      int index_bLep  =  mapFilt[pos_to_index[bLep_pos]].index;
		      float csv_bLep_nominal  =  index_bLep>=0 ? hJet_csv_nominal[index_bLep]  : aJet_csv_nominal[-index_bLep-1];
		      float csv_bLep_upBC     =  index_bLep>=0 ? hJet_csv_upBC   [index_bLep]  : aJet_csv_upBC   [-index_bLep-1];
		      float csv_bLep_downBC   =  index_bLep>=0 ? hJet_csv_downBC [index_bLep]  : aJet_csv_downBC [-index_bLep-1];
		      float csv_bLep_upL      =  index_bLep>=0 ? hJet_csv_upL    [index_bLep]  : aJet_csv_upL    [-index_bLep-1];
		      float csv_bLep_downL    =  index_bLep>=0 ? hJet_csv_downL  [index_bLep]  : aJet_csv_downL  [-index_bLep-1];
		      double csv_bLep = csv_bLep_nominal;
		      if( doCSVup )   csv_bLep = TMath::Max(csv_bLep_upBC,   csv_bLep_upL);
		      if( doCSVdown ) csv_bLep = TMath::Min(csv_bLep_downBC, csv_bLep_downL);

		      int bHad_pos    = (permutations_SL2wj[pos])%1000/100;
		      int index_bHad  =  mapFilt[pos_to_index[bHad_pos]].index;
		      //float csv_bHad  =  index_bHad>=0 ? hJet_csv_nominal[index_bHad]  : aJet_csv_nominal[-index_bHad-1];
		      float csv_bHad_nominal  =  index_bHad>=0 ? hJet_csv_nominal[index_bHad]  : aJet_csv_nominal[-index_bHad-1];
		      float csv_bHad_upBC     =  index_bHad>=0 ? hJet_csv_upBC   [index_bHad]  : aJet_csv_upBC   [-index_bHad-1];
		      float csv_bHad_downBC   =  index_bHad>=0 ? hJet_csv_downBC [index_bHad]  : aJet_csv_downBC [-index_bHad-1];
		      float csv_bHad_upL      =  index_bHad>=0 ? hJet_csv_upL    [index_bHad]  : aJet_csv_upL    [-index_bHad-1];
		      float csv_bHad_downL    =  index_bHad>=0 ? hJet_csv_downL  [index_bHad]  : aJet_csv_downL  [-index_bHad-1];
		      double csv_bHad = csv_bHad_nominal;
		      if( doCSVup )   csv_bHad = TMath::Max(csv_bHad_upBC,   csv_bHad_upL);
		      if( doCSVdown ) csv_bHad = TMath::Min(csv_bHad_downBC, csv_bHad_downL);

		      int b1_pos    =   (permutations_SL2wj[pos])%100/10;
		      int index_b1  =  mapFilt[pos_to_index[b1_pos]].index;
		      //float csv_b1  =  index_b1>=0 ? hJet_csv_nominal[index_b1]  : aJet_csv_nominal[-index_b1-1];
		      float csv_b1_nominal  =  index_b1>=0 ? hJet_csv_nominal[index_b1]  : aJet_csv_nominal[-index_b1-1];
		      float csv_b1_upBC     =  index_b1>=0 ? hJet_csv_upBC   [index_b1]  : aJet_csv_upBC   [-index_b1-1];
		      float csv_b1_downBC   =  index_b1>=0 ? hJet_csv_downBC [index_b1]  : aJet_csv_downBC [-index_b1-1];
		      float csv_b1_upL      =  index_b1>=0 ? hJet_csv_upL    [index_b1]  : aJet_csv_upL    [-index_b1-1];
		      float csv_b1_downL    =  index_b1>=0 ? hJet_csv_downL  [index_b1]  : aJet_csv_downL  [-index_b1-1];
		      double csv_b1 = csv_b1_nominal;
		      if( doCSVup )   csv_b1 = TMath::Max(csv_b1_upBC,   csv_b1_upL);
		      if( doCSVdown ) csv_b1 = TMath::Min(csv_b1_downBC, csv_b1_downL);

		      int b2_pos    = (permutations_SL2wj[pos])%10;
		      int index_b2  =  mapFilt[pos_to_index[b2_pos]].index;
		      //float csv_b2  =  index_b2>=0 ? hJet_csv_nominal[index_b2]  : aJet_csv_nominal[-index_b2-1];
		      float csv_b2_nominal  =  index_b2>=0 ? hJet_csv_nominal[index_b2]  : aJet_csv_nominal[-index_b2-1];
		      float csv_b2_upBC     =  index_b2>=0 ? hJet_csv_upBC   [index_b2]  : aJet_csv_upBC   [-index_b2-1];
		      float csv_b2_downBC   =  index_b2>=0 ? hJet_csv_downBC [index_b2]  : aJet_csv_downBC [-index_b2-1];
		      float csv_b2_upL      =  index_b2>=0 ? hJet_csv_upL    [index_b2]  : aJet_csv_upL    [-index_b2-1];
		      float csv_b2_downL    =  index_b2>=0 ? hJet_csv_downL  [index_b2]  : aJet_csv_downL  [-index_b2-1];
		      double csv_b2 = csv_b2_nominal;
		      if( doCSVup )   csv_b2 = TMath::Max(csv_b2_upBC,   csv_b2_upL);
		      if( doCSVdown ) csv_b2 = TMath::Min(csv_b2_downBC, csv_b2_downL);

		      
		      string bin_bLep = TMath::Abs( (meIntegrator->jetAt(2)).Eta() )<1.0 ? "Bin0" : "Bin1"; 
		      double p_b_bLep =  btagger["b_"+bin_bLep]!=0 ?  btagger["b_"+bin_bLep]->GetBinContent( btagger["b_"+bin_bLep]->FindBin( TMath::Min(csv_bLep, 0.999999) ) ) : 1.;
		      
		      string bin_bHad = TMath::Abs( (meIntegrator->jetAt(5)).Eta() )<1.0 ? "Bin0" : "Bin1"; 
		      double p_b_bHad =  btagger["b_"+bin_bHad]!=0 ?  btagger["b_"+bin_bHad]->GetBinContent( btagger["b_"+bin_bHad]->FindBin( TMath::Min(csv_bHad, 0.999999) ) ) : 1.;

		      string bin_b1   = TMath::Abs( (meIntegrator->jetAt(6)).Eta() )<1.0 ? "Bin0" : "Bin1"; 
		      double p_b_b1   =  btagger["b_"+bin_b1]!=0 ?  btagger["b_"+bin_b1]->GetBinContent( btagger["b_"+bin_b1]->FindBin( TMath::Min(csv_b1, 0.999999) ) ) : 1.;
		      double p_j_b1   =  btagger["l_"+bin_b1]!=0 ?  btagger["l_"+bin_b1]->GetBinContent( btagger["l_"+bin_b1]->FindBin( TMath::Min(csv_b1, 0.999999) ) ) : 1.;

		      string bin_b2   = TMath::Abs( (meIntegrator->jetAt(7)).Eta() )<1.0 ? "Bin0" : "Bin1"; 
		      double p_b_b2   =  btagger["b_"+bin_b2]!=0 ?  btagger["b_"+bin_b2]->GetBinContent( btagger["b_"+bin_b2]->FindBin( TMath::Min(csv_b2, 0.999999) ) ) : 1.;
		      double p_j_b2   =  btagger["l_"+bin_b2]!=0 ?  btagger["l_"+bin_b2]->GetBinContent( btagger["l_"+bin_b2]->FindBin( TMath::Min(csv_b2, 0.999999) ) ) : 1.;
		      
		      //cout << "bLep_pos " << bLep_pos << " -- correspond to " << index_bLep << " with csv="<< csv_bLep << ": prob under b hyp: " << p_b_bLep << endl;
		      //cout << "bHad_pos " << bHad_pos << " -- correspond to " << index_bHad << " with csv="<< csv_bHad << ": prob under b hyp: " << p_b_bHad <<  endl;
		      //cout << "b1_pos " << b1_pos << " -- correspond to " << index_b1 << " with csv="<< csv_b1         << ": prob under b hyp: " << p_b_b1 << ", under j hyp: " << p_j_b1 <<  endl;
		      //cout << "b2_pos " << b2_pos << " -- correspond to " << index_b2 << " with csv="<< csv_b2         << ": prob under b hyp: " << p_b_b2 << ", under j hyp: " << p_j_b2 << endl;		      

		      /*
		      int sv_bLep = 0;
		      int sv_bHad = 0;
		      int sv_b1   = 0;
		      int sv_b2   = 0;
		      for(int l = 0; l < nSvs; l++){
			TLorentzVector vertex(1,0,0,1);
			vertex.SetPtEtaPhiM( Svpt[l], Sveta[l], Svphi[l], SvmassSv[l]);
			if(deltaR(vertex, meIntegrator->jetAt(2) )<0.5)      sv_bLep = 1;
			else if(deltaR(vertex, meIntegrator->jetAt(5) )<0.5) sv_bHad = 1;
			else if(deltaR(vertex, meIntegrator->jetAt(6) )<0.5) sv_b1   = 1;
			else if(deltaR(vertex, meIntegrator->jetAt(7) )<0.5) sv_b2   = 1;
			else{}
		      }
		      double p_b_sv = svs[sv_bLep]*svs[sv_bHad]*svs[sv_b1]  *svs[sv_b2];
		      double p_j_sv = svs[sv_bLep]*svs[sv_bHad]*svs[2+sv_b1]*svs[2+sv_b2];
		      */

		      if( mH[m]<met+0.5 && mH[m]>met-0.5){
			if(hyp==0){
			  probAtSgn_ttbb_         += ( p * p_b_bLep * p_b_bHad * p_b_b1 * p_b_b2 /**  p_b_sv*/);
			}
			else{
			  probAtSgn_alt_ttbb_     += ( p * p_b_bLep * p_b_bHad * p_b_b1 * p_b_b2 /**  p_b_sv*/);
			  probAtSgn_alt_ttjj_     += ( p * p_b_bLep * p_b_bHad * p_j_b1 * p_j_b2 /**  p_j_sv*/);
			}
		      }
		      
		    }


		    if( p>= maxP_s && hyp==0 ){

		      int lastTwoDigits = (permutations_SL2wj[pos])%100;
		      int higgs1 = lastTwoDigits/(10);
		      int higgs2 = lastTwoDigits%10;

		      //cout << "pos=" << permutations_SL2wj[pos] << " => " << higgs1 << "," << higgs2 << endl;

		      if( higgs1 == 2)
			higgs1 = b1;		      
		      else if(  higgs1 == 5 )
			higgs1 = b2;
		      else if(  higgs1 == 6 )
			higgs1 = b3;
		      else if(  higgs1 == 7 )
			higgs1 = b4;
		      else{ cout << "What happened?" << endl; }
		      if( higgs2 == 2)
			higgs2 = b1;		      
		      else if(  higgs2 == 5 )
			higgs2 = b2;
		      else if(  higgs2 == 6 )
			higgs2 = b3;
		      else if(  higgs2 == 7 )
			higgs2 = b4;
		      else{ cout << "What happened?" << endl; }
			
		      //cout << "Corrspond to " << higgs1 << "," << higgs2 << endl;

		      int index1 =  mapFilt[higgs1].index;		      
		      int index2 =  mapFilt[higgs2].index;		      

		      float csv_1 = index1>=0 ? hJet_csv_nominal[index1]  : aJet_csv_nominal[-index1-1];
		      float csv_2 = index2>=0 ? hJet_csv_nominal[index2]  : aJet_csv_nominal[-index2-1];

		      //cout << "Deug 1: " << mapFilt[higgs1].pt << " --- " << (meIntegrator->jetAt(6)).Pt() << endl;
		      //cout << "Deug 2: " << mapFilt[higgs2].pt << " --- " << (meIntegrator->jetAt(7)).Pt() << endl;

		      commitVariables( meIntegrator->jetAt(6), meIntegrator->jetAt(7), csv_1, csv_2,
				       dR_,dPhi_,dEta_,Mjj_,Pt_,Eta_,pt1_,pt2_,eta1_,eta2_,csv1_,csv2_);
		    }

		  
		    if( mH[m]<met+0.5 && mH[m]>met-0.5){
		      if(hyp==0)  probAtSgn_     += p;
		      else        probAtSgn_alt_ += p;
		    }

		  }		  
		}
		
	      }
	    }
	  }

	  clock->Stop();
	  time_ = clock->RealTime();
	  clock->Reset();	    
	  
	}
	else if(numJets30UntagM>=3 && doType3){
	  
	  //continue;

	  counter++;      
	  if(fixNumberOfEventsPerJob && !(counter>=evLow && counter<=evHigh) ) continue;
	  cout << "Processing event # " << counter << " (SL[4,3])" << endl;
	 
	  
	  if( w1==999 || w2==999 || w3==999 || b1==999 || b2==999 || b3==999 || b4==999){
	    cout << "Inconsistentcy found!" << endl;
	    continue;
	  }

	  /*
	  float WMass12 = (myJetsFilt[w1]+myJetsFilt[w2]).M();
	  float WMass13 = (myJetsFilt[w1]+myJetsFilt[w3]).M();
	  float WMass23 = (myJetsFilt[w2]+myJetsFilt[w3]).M();
	  
	  unsigned int ind1 = 999;
	  unsigned int ind2 = 999;
	  if( TMath::Abs(WMass12-80.1) < TMath::Abs(WMass13-80.1) &&
	      TMath::Abs(WMass12-80.1) < TMath::Abs(WMass23-80.1)){
	    ind1 = w1;
	    ind2 = w2;
	    cout << "Best W mass: " << WMass12 << endl;
	  }
	  else if( TMath::Abs(WMass13-80.1) < TMath::Abs(WMass12-80.1) &&
		   TMath::Abs(WMass13-80.1) < TMath::Abs(WMass23-80.1) ){
	    ind1 = w1;
	    ind2 = w3;
	    cout << "Best W mass: " << WMass13 << endl;
	  }
	  else if( TMath::Abs(WMass23-80.1) < TMath::Abs(WMass12-80.1) &&
		   TMath::Abs(WMass23-80.1) < TMath::Abs(WMass13-80.1) ){
	    ind1 = w2;
	    ind2 = w3;
	    cout << "Best W mass: " << WMass23 << endl;
	  }
	  else{
	    cout << "Inconsistency found in numJets30UntagM==3... continue" << endl;
	    continue;
	  }
	  */


	  
	  vector<unsigned int> untaggedjets;
	  if( w1!=999 ) untaggedjets.push_back(w1);
	  if( w2!=999 ) untaggedjets.push_back(w2);
	  if( w3!=999 ) untaggedjets.push_back(w3);
	  if( w4!=999 ) untaggedjets.push_back(w4);
	  if( w5!=999 ) untaggedjets.push_back(w5);

	  unsigned int ind1 = 999;
	  unsigned int ind2 = 999;
	  float minDiff     = 9999.;
	  for(unsigned int uj1 = 0; uj1<untaggedjets.size()-1; uj1++){
	    for(unsigned int uj2 = uj1+1; uj2<untaggedjets.size(); uj2++){
	      float WMass12 = (myJetsFilt[ untaggedjets[uj1] ]+myJetsFilt[ untaggedjets[uj2] ]).M();
	      if( TMath::Abs(WMass12-80.1)<minDiff ){
		minDiff = TMath::Abs(WMass12-80.1);
		ind1 = untaggedjets[uj1];
		ind2 = untaggedjets[uj2];
	      }
	    }
	  }
	  bool wtag = ind1!=999 && ind2!=999;
	  if(wtag){
	    float WMass = (myJetsFilt[ ind1 ]+myJetsFilt[ ind2 ]).M();
	    cout << "Best W mass: " << WMass << endl;
	    if( WMass>MwL && WMass<MwH ) 
	      wtag = true;
	    else 
	      wtag = false;
	  }
	  if(  wtag  ) flag_type3_ = 1;
	  



	  jets.clear();
	  jets.push_back( leptonLV  );  
	  jets.push_back( neutrinoLV  );  
	  jets.push_back( myJetsFilt[b1]   );
	  jets.push_back( myJetsFilt[ind1]   );  
	  jets.push_back( myJetsFilt[ind2]   );  
	  jets.push_back( myJetsFilt[b2]   );
	  jets.push_back( myJetsFilt[b3]   );
	  jets.push_back( myJetsFilt[b4]   );

	  pos_to_index.clear();
	  pos_to_index[2] = b1;
	  pos_to_index[3] = ind1;
	  pos_to_index[4] = ind2;
	  pos_to_index[5] = b2;
	  pos_to_index[6] = b3;
	  pos_to_index[7] = b4;

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


	    //if( Bs.Pt()>20 && TMath::Abs(Bs.Eta())<2.5 ) nMatchSimBs_++;
	    //if     ( deltaR( Bs , myJetsFilt[b1] )<0.2 ) nMatchSimBs_++;
	    //else if( deltaR( Bs , myJetsFilt[b2] )<0.2 ) nMatchSimBs_++;
	    //else if( deltaR( Bs , myJetsFilt[b3] )<0.2 ) nMatchSimBs_++;
	    //else if( deltaR( Bs , myJetsFilt[b4] )<0.2 ) nMatchSimBs_++;
	    //else{
	    //for(unsigned int uj = 0; uj<untaggedjets.size(); uj++){
	    //if ( deltaR( Bs , myJetsFilt[ untaggedjets[uj] ] )<0.2 ) nMatchSimBs_++;
	    //}
	    //}
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

	  meIntegrator->setJets(&jets);

	  type_ = 3;

	  meIntegrator->setIntType( MEIntegratorNew::SL2wj );	    
	  meIntegrator->setSumEt( METtype1p2corr.sumet );
	  meIntegrator->setMEtCov(-99,-99,0);
	  meIntegrator->setTopFlags( vLepton_charge[0]==1 ? +1 : -1 , vLepton_charge[0]==1 ? -1 : +1 );

	  clock->Start();

	  for(int m = 0; m < nMassPoints ; m++){
	    meIntegrator->setMass( mH[m] );
	      
	    double maxP_s = 0.;
	    double maxP_b = 0.;

	    for(unsigned int pos = 0; pos < 12 ; pos++){

	      meIntegrator->initVersors( permutations_SL2wj[pos] );
	      double mass, massLow, massHigh;
	      bool skip = !(meIntegrator->compatibilityCheck(0.95, /*printP4*/ 0, mass, massLow, massHigh )) ;
	      bool skip_WHad   = !( meIntegrator->compatibilityCheck_WHad  (0.98, /*printP4*/ 0, mass, massLow, massHigh ) );
	      bool skip_TopHad = !( meIntegrator->compatibilityCheck_TopHad(0.98, /*printP4*/ 0, mass, massLow, massHigh ) );
	      if( (skip_WHad || skip_TopHad) ) continue;

	      pair<double, double> range_x0 = (meIntegrator->getW1JetEnergyCI(0.95));
	      pair<double, double> range_x1 =  make_pair(-1,1);
	      pair<double, double> range_x2 =  make_pair(-PI,PI);	    
	      pair<double, double> range_x3 =  make_pair(-1,1);
	      pair<double, double> range_x4 = useMET ? (meIntegrator->getNuPhiCI(0.95)) : make_pair(-PI,PI);
	      pair<double, double> range_x5 = (meIntegrator->getB1EnergyCI(0.95));
	      pair<double, double> range_x6 = (meIntegrator->getB2EnergyCI(0.95));
	      
	      double x0L = range_x0.first;
	      double x0U = range_x0.second;
	      double x1L = range_x1.first;
	      double x1U = range_x1.second;
	      double x2L = range_x2.first;
	      double x2U = range_x2.second;
	      double x3L = range_x3.first;
	      double x3U = range_x3.second;
	      double x4L = range_x4.first;
	      double x4U = range_x4.second;
	      double x5L = range_x5.first;
	      double x5U = range_x5.second;
	      double x6L = range_x6.first;
	      double x6U = range_x6.second;
	
	      double xLmode0[4]   = {x0L, x3L, x4L, x5L};
	      double xUmode0[4]   = {x0U, x3U, x4U, x5U};
	      double xLmode0_b[5] = {x0L, x3L, x4L, x5L, x6L};
	      double xUmode0_b[5] = {x0U, x3U, x4U, x5U, x6U};	     
		
	      double p = 0.;
	      
	      if(speedup==0){
		
		for(int hyp = 0 ; hyp<2; hyp++){
		  
		  if(SoB==0 && hyp!=hypo) continue;
		  if(SoB==1 && hyp==0 && skip) continue;
		  
		  cout << "Testing hypothesis " << (hyp==0 ? "S" : "B") << " for permutation " << pos << endl; 
		  meIntegrator->setHypo(hyp);
		  
		  int ntries = 0;
		  int intPoints = 2000;
		  while( ntries < MAX_REEVAL_TRIES){
		    
		    ROOT::Math::Functor toIntegrate(meIntegrator, &MEIntegratorNew::Eval,  4+hyp );
		    ROOT::Math::GSLMCIntegrator ig2( ROOT::Math::IntegrationMultiDim::kVEGAS , 1.e-12, 1.e-5, intPoints);
		    ig2.SetFunction(toIntegrate);
		    meIntegrator->SetPar( 4+hyp );	 
		    
		    p = (hyp==0 ? ig2.Integral(xLmode0, xUmode0) : ig2.Integral(xLmode0_b, xUmode0_b));
		    
		    double chi2 =  ig2.ChiSqr();
		    if( hyp==0 ){
		      if( p>maxP_s ) maxP_s = p;		  
		      else if( p<0.1*maxP_s) ntries = MAX_REEVAL_TRIES+1;
		    }
		    else{
		      if( p>maxP_b ) maxP_b = p;		  
		      else if( p<0.1*maxP_b) ntries = MAX_REEVAL_TRIES+1;	
		    }
		    
		    if( chi2 > maxChi2_ ){
		      ntries++;
		      intPoints *= 1.5;
		    }
		    else ntries = MAX_REEVAL_TRIES+1;
		    
		  }		   
				    
		  if( TMath::IsNaN(p) ) p = 0.;


		  if( useBtag ){
		      
		    int bLep_pos    = (permutations_SL2wj[pos])/100000;
		    int index_bLep  =  mapFilt[pos_to_index[bLep_pos]].index;
		    float csv_bLep_nominal  =  index_bLep>=0 ? hJet_csv_nominal[index_bLep]  : aJet_csv_nominal[-index_bLep-1];
		    float csv_bLep_upBC     =  index_bLep>=0 ? hJet_csv_upBC   [index_bLep]  : aJet_csv_upBC   [-index_bLep-1];
		    float csv_bLep_downBC   =  index_bLep>=0 ? hJet_csv_downBC [index_bLep]  : aJet_csv_downBC [-index_bLep-1];
		    float csv_bLep_upL      =  index_bLep>=0 ? hJet_csv_upL    [index_bLep]  : aJet_csv_upL    [-index_bLep-1];
		    float csv_bLep_downL    =  index_bLep>=0 ? hJet_csv_downL  [index_bLep]  : aJet_csv_downL  [-index_bLep-1];
		    double csv_bLep = csv_bLep_nominal;
		    if( doCSVup )   csv_bLep = TMath::Max(csv_bLep_upBC,   csv_bLep_upL);
		    if( doCSVdown ) csv_bLep = TMath::Min(csv_bLep_downBC, csv_bLep_downL);
		    
		    int bHad_pos    = (permutations_SL2wj[pos])%1000/100;
		    int index_bHad  =  mapFilt[pos_to_index[bHad_pos]].index;
		    //float csv_bHad  =  index_bHad>=0 ? hJet_csv_nominal[index_bHad]  : aJet_csv_nominal[-index_bHad-1];
		    float csv_bHad_nominal  =  index_bHad>=0 ? hJet_csv_nominal[index_bHad]  : aJet_csv_nominal[-index_bHad-1];
		    float csv_bHad_upBC     =  index_bHad>=0 ? hJet_csv_upBC   [index_bHad]  : aJet_csv_upBC   [-index_bHad-1];
		    float csv_bHad_downBC   =  index_bHad>=0 ? hJet_csv_downBC [index_bHad]  : aJet_csv_downBC [-index_bHad-1];
		    float csv_bHad_upL      =  index_bHad>=0 ? hJet_csv_upL    [index_bHad]  : aJet_csv_upL    [-index_bHad-1];
		    float csv_bHad_downL    =  index_bHad>=0 ? hJet_csv_downL  [index_bHad]  : aJet_csv_downL  [-index_bHad-1];
		    double csv_bHad = csv_bHad_nominal;
		    if( doCSVup )   csv_bHad = TMath::Max(csv_bHad_upBC,   csv_bHad_upL);
		    if( doCSVdown ) csv_bHad = TMath::Min(csv_bHad_downBC, csv_bHad_downL);
		    
		    int b1_pos    =   (permutations_SL2wj[pos])%100/10;
		    int index_b1  =  mapFilt[pos_to_index[b1_pos]].index;
		    //float csv_b1  =  index_b1>=0 ? hJet_csv_nominal[index_b1]  : aJet_csv_nominal[-index_b1-1];
		    float csv_b1_nominal  =  index_b1>=0 ? hJet_csv_nominal[index_b1]  : aJet_csv_nominal[-index_b1-1];
		    float csv_b1_upBC     =  index_b1>=0 ? hJet_csv_upBC   [index_b1]  : aJet_csv_upBC   [-index_b1-1];
		    float csv_b1_downBC   =  index_b1>=0 ? hJet_csv_downBC [index_b1]  : aJet_csv_downBC [-index_b1-1];
		    float csv_b1_upL      =  index_b1>=0 ? hJet_csv_upL    [index_b1]  : aJet_csv_upL    [-index_b1-1];
		    float csv_b1_downL    =  index_b1>=0 ? hJet_csv_downL  [index_b1]  : aJet_csv_downL  [-index_b1-1];
		    double csv_b1 = csv_b1_nominal;
		    if( doCSVup )   csv_b1 = TMath::Max(csv_b1_upBC,   csv_b1_upL);
		    if( doCSVdown ) csv_b1 = TMath::Min(csv_b1_downBC, csv_b1_downL);
		    
		    int b2_pos    = (permutations_SL2wj[pos])%10;
		    int index_b2  =  mapFilt[pos_to_index[b2_pos]].index;
		    //float csv_b2  =  index_b2>=0 ? hJet_csv_nominal[index_b2]  : aJet_csv_nominal[-index_b2-1];
		    float csv_b2_nominal  =  index_b2>=0 ? hJet_csv_nominal[index_b2]  : aJet_csv_nominal[-index_b2-1];
		    float csv_b2_upBC     =  index_b2>=0 ? hJet_csv_upBC   [index_b2]  : aJet_csv_upBC   [-index_b2-1];
		    float csv_b2_downBC   =  index_b2>=0 ? hJet_csv_downBC [index_b2]  : aJet_csv_downBC [-index_b2-1];
		    float csv_b2_upL      =  index_b2>=0 ? hJet_csv_upL    [index_b2]  : aJet_csv_upL    [-index_b2-1];
		    float csv_b2_downL    =  index_b2>=0 ? hJet_csv_downL  [index_b2]  : aJet_csv_downL  [-index_b2-1];
		    double csv_b2 = csv_b2_nominal;
		    if( doCSVup )   csv_b2 = TMath::Max(csv_b2_upBC,   csv_b2_upL);
		    if( doCSVdown ) csv_b2 = TMath::Min(csv_b2_downBC, csv_b2_downL);
	    

		    string bin_bLep = TMath::Abs( (meIntegrator->jetAt(2)).Eta() )<1.0 ? "Bin0" : "Bin1"; 
		    double p_b_bLep =  btagger["b_"+bin_bLep]!=0 ?  btagger["b_"+bin_bLep]->GetBinContent( btagger["b_"+bin_bLep]->FindBin( TMath::Min(csv_bLep, 0.999999) ) ) : 1.;
		    
		    string bin_bHad = TMath::Abs( (meIntegrator->jetAt(5)).Eta() )<1.0 ? "Bin0" : "Bin1"; 
		    double p_b_bHad =  btagger["b_"+bin_bHad]!=0 ?  btagger["b_"+bin_bHad]->GetBinContent( btagger["b_"+bin_bHad]->FindBin( TMath::Min(csv_bHad, 0.999999) ) ) : 1.;
		    
		    string bin_b1   = TMath::Abs( (meIntegrator->jetAt(6)).Eta() )<1.0 ? "Bin0" : "Bin1"; 
		    double p_b_b1   =  btagger["b_"+bin_b1]!=0 ?  btagger["b_"+bin_b1]->GetBinContent( btagger["b_"+bin_b1]->FindBin( TMath::Min(csv_b1, 0.999999) ) ) : 1.;
		    double p_j_b1   =  btagger["l_"+bin_b1]!=0 ?  btagger["l_"+bin_b1]->GetBinContent( btagger["l_"+bin_b1]->FindBin( TMath::Min(csv_b1, 0.999999) ) ) : 1.;
		    
		    string bin_b2   = TMath::Abs( (meIntegrator->jetAt(7)).Eta() )<1.0 ? "Bin0" : "Bin1"; 
		    double p_b_b2   =  btagger["b_"+bin_b2]!=0 ?  btagger["b_"+bin_b2]->GetBinContent( btagger["b_"+bin_b2]->FindBin( TMath::Min(csv_b2, 0.999999) ) ) : 1.;
		    double p_j_b2   =  btagger["l_"+bin_b2]!=0 ?  btagger["l_"+bin_b2]->GetBinContent( btagger["l_"+bin_b2]->FindBin( TMath::Min(csv_b2, 0.999999) ) ) : 1.;
		   
		     
		    //cout << "bLep_pos " << bLep_pos << " -- correspond to " << index_bLep << " with csv="<< csv_bLep << ": prob under b hyp: " << p_b_bLep << endl;
		    //cout << "bHad_pos " << bHad_pos << " -- correspond to " << index_bHad << " with csv="<< csv_bHad << ": prob under b hyp: " << p_b_bHad <<  endl;
		    //cout << "b1_pos " << b1_pos << " -- correspond to " << index_b1 << " with csv="<< csv_b1         << ": prob under b hyp: " << p_b_b1 << ", under j hyp: " << p_j_b1 <<  endl;
		    //cout << "b2_pos " << b2_pos << " -- correspond to " << index_b2 << " with csv="<< csv_b2         << ": prob under b hyp: " << p_b_b2 << ", under j hyp: " << p_j_b2 << endl;		      

		    /*
		    int sv_bLep = 0;
		    int sv_bHad = 0;
		    int sv_b1   = 0;
		    int sv_b2   = 0;
		    for(int l = 0; l < nSvs; l++){
		      TLorentzVector vertex(1,0,0,1);
		      vertex.SetPtEtaPhiM( Svpt[l], Sveta[l], Svphi[l], SvmassSv[l]);
		      if(deltaR(vertex, meIntegrator->jetAt(2) )<0.5)      sv_bLep = 1;
		      else if(deltaR(vertex, meIntegrator->jetAt(5) )<0.5) sv_bHad = 1;
		      else if(deltaR(vertex, meIntegrator->jetAt(6) )<0.5) sv_b1   = 1;
		      else if(deltaR(vertex, meIntegrator->jetAt(7) )<0.5) sv_b2   = 1;
		      else{}
		    }
		    double p_b_sv = svs[sv_bLep]*svs[sv_bHad]*svs[sv_b1]  *svs[sv_b2];
		    double p_j_sv = svs[sv_bLep]*svs[sv_bHad]*svs[2+sv_b1]*svs[2+sv_b2];
		    */

		    if( mH[m]<met+0.5 && mH[m]>met-0.5){
		      if(hyp==0){
			probAtSgn_ttbb_         += ( p * p_b_bLep * p_b_bHad * p_b_b1 * p_b_b2 /**  p_b_sv*/);
		      }
		      else{
			probAtSgn_alt_ttbb_     += ( p * p_b_bLep * p_b_bHad * p_b_b1 * p_b_b2 /**  p_b_sv*/);
			probAtSgn_alt_ttjj_     += ( p * p_b_bLep * p_b_bHad * p_j_b1 * p_j_b2 /**  p_j_sv*/);
		      }
		    }		    
		  }

		  
		  if( p>= maxP_s && hyp==0 ){
		    int lastTwoDigits = (permutations_SL2wj[pos])%100;
		    int higgs1 = lastTwoDigits/(10);
		    int higgs2 = lastTwoDigits%10;		    
		    //cout << "pos=" << permutations_SL2wj[pos] << " => " << higgs1 << "," << higgs2 << endl;		    
		    if( higgs1 == 2)
		      higgs1 = b1;		      
		    else if(  higgs1 == 5 )
		      higgs1 = b2;
		    else if(  higgs1 == 6 )
		      higgs1 = b3;
		    else if(  higgs1 == 7 )
		      higgs1 = b4;
		    else{ cout << "What happened?" << endl; }
		    if( higgs2 == 2)
		      higgs2 = b1;		      
		    else if(  higgs2 == 5 )
		      higgs2 = b2;
		    else if(  higgs2 == 6 )
		      higgs2 = b3;
		    else if(  higgs2 == 7 )
		      higgs2 = b4;
		    else{ cout << "What happened?" << endl; }		    
		    //cout << "Corrspond to " << higgs1 << "," << higgs2 << endl;		    
		    int index1 =  mapFilt[higgs1].index;		      
		    int index2 =  mapFilt[higgs2].index;		      		    
		    float csv_1 = index1>=0 ? hJet_csv_nominal[index1]  : aJet_csv_nominal[-index1-1];
		    float csv_2 = index2>=0 ? hJet_csv_nominal[index2]  : aJet_csv_nominal[-index2-1];		    
		    //cout << "Deug 1: " << mapFilt[higgs1].pt << " --- " << (meIntegrator->jetAt(6)).Pt() << endl;
		    //cout << "Deug 2: " << mapFilt[higgs2].pt << " --- " << (meIntegrator->jetAt(7)).Pt() << endl;		    
		    commitVariables( meIntegrator->jetAt(6), meIntegrator->jetAt(7), csv_1, csv_2,
				     dR_,dPhi_,dEta_,Mjj_,Pt_,Eta_,pt1_,pt2_,eta1_,eta2_,csv1_,csv2_);
		  }

		  if( mH[m]<met+0.5 && mH[m]>met-0.5){
		    if(hyp==0)  probAtSgn_     += p;
		    else        probAtSgn_alt_ += p;
		  }
		}		  
	      }	
		
	    }
	  }
	  
	  clock->Stop();
	  time_ = clock->RealTime();
	  clock->Reset();	 
	  
	}
	else{ 
	  continue;
	}


      }
      else if(properEventSL && numJets30BtagM==4 && numJets30UntagM==1 && doType2){

	counter++;
	if(fixNumberOfEventsPerJob && !(counter>=evLow && counter<=evHigh) ) continue;
	cout << "Processing event # " << counter << " (SL[4,1])" << endl;
	 
	if( w1==999  || b1==999 || b2==999 || b3==999 || b4==999){
	  cout << "Inconsistentcy type 2 found!" << endl;
	  cout << w1 << ", " << w2 << "," << b1 << "," << b2 << "," << b3 << "," << b4 <<  endl;
	  continue;
	}

	/*
	float WMass1 = (myJetsFilt[w1]+myJetsFilt[b1]).M();
	float WMass2 = (myJetsFilt[w1]+myJetsFilt[b2]).M();
	float WMass3 = (myJetsFilt[w1]+myJetsFilt[b3]).M();
	float WMass4 = (myJetsFilt[w1]+myJetsFilt[b4]).M();
	
	bool wtag = ( (WMass1>(MwL+5) && WMass1<(MwH-5)) || 
		      (WMass2>(MwL+5) && WMass2<(MwH-5)) || 
		      (WMass3>(MwL+5) && WMass3<(MwH-5)) || 
		      (WMass4>(MwL+5) && WMass4<(MwH-5))
		      );
	
	
	if(  wtag  ) flag_type2_ = 1;
	*/


	int hMatches = 0;
	if     (  deltaR(myJetsFilt[b1],genBLV)<0.3 )    hMatches++;
	else if(  deltaR(myJetsFilt[b1],genBbarLV)<0.3)  hMatches++;
	if     (  deltaR(myJetsFilt[b2],genBLV)<0.3 )    hMatches++;
	else if(  deltaR(myJetsFilt[b2],genBbarLV)<0.3)  hMatches++;
	if     (  deltaR(myJetsFilt[b3],genBLV)<0.3 )    hMatches++;
	else if(  deltaR(myJetsFilt[b3],genBbarLV)<0.3)  hMatches++;
	if     (  deltaR(myJetsFilt[b4],genBLV)<0.3 )    hMatches++;
	else if(  deltaR(myJetsFilt[b4],genBbarLV)<0.3)  hMatches++;
	matchesH_    = hMatches;
	matchesHAll_ = hMatches;
	if     (  deltaR(myJetsFilt[w1],genBLV)<0.3 )    matchesHAll_++;
	else if(  deltaR(myJetsFilt[w1],genBbarLV)<0.3)  matchesHAll_++;
	
	int tMatches = 0;
	if     (  deltaR(myJetsFilt[b1],TOPHADB)<0.3 )   tMatches++;
	else if(  deltaR(myJetsFilt[b1],TOPLEPB)<0.3)    tMatches++;
	if     (  deltaR(myJetsFilt[b2],TOPHADB)<0.3 )   tMatches++;
	else if(  deltaR(myJetsFilt[b2],TOPLEPB)<0.3)    tMatches++;
	if     (  deltaR(myJetsFilt[b3],TOPHADB)<0.3 )   tMatches++;
	else if(  deltaR(myJetsFilt[b3],TOPLEPB)<0.3)    tMatches++;
	if     (  deltaR(myJetsFilt[b4],TOPHADB)<0.3 )   tMatches++;
	else if(  deltaR(myJetsFilt[b4],TOPLEPB)<0.3)    tMatches++;
	matchesT_    = tMatches;
	matchesTAll_ = tMatches;
	if     (  deltaR(myJetsFilt[w1],TOPHADB)<0.3 )   matchesTAll_++;
	else if(  deltaR(myJetsFilt[w1],TOPLEPB)<0.3)    matchesTAll_++;


	int wMatches = 0;
	if     (  deltaR(myJetsFilt[w1],TOPHADW1)<0.3 ){
	  wMatches++;
	  mL_[0] = 1;
	}
	else if(  deltaR(myJetsFilt[w1],TOPHADW2)<0.3){
	  wMatches++;
	  mL_[0] = 1;
	}
	else{ mL_[0] = 0; }
	matchesW_    = wMatches;
	matchesWAll_ = wMatches;	  
	if     (  deltaR(myJetsFilt[b1],TOPHADW1)<0.3 )   matchesWAll_++;
	else if(  deltaR(myJetsFilt[b1],TOPHADW2)<0.3)    matchesWAll_++;
	if     (  deltaR(myJetsFilt[b2],TOPHADW1)<0.3 )   matchesWAll_++;
	else if(  deltaR(myJetsFilt[b2],TOPHADW2)<0.3)    matchesWAll_++;
	if     (  deltaR(myJetsFilt[b3],TOPHADW1)<0.3 )   matchesWAll_++;
	else if(  deltaR(myJetsFilt[b3],TOPHADW2)<0.3)    matchesWAll_++;
	if     (  deltaR(myJetsFilt[b4],TOPHADW1)<0.3 )   matchesWAll_++;
	else if(  deltaR(myJetsFilt[b4],TOPHADW2)<0.3)    matchesWAll_++;
	
	vector<TLorentzVector> genHeavy;
	genHeavy.push_back( genBLV);
	genHeavy.push_back( genBbarLV);
	genHeavy.push_back( TOPHADB);
	genHeavy.push_back( TOPLEPB);
	int overlapH = 0;
	for(unsigned int k = 0; k < genHeavy.size()-1; k++){  
	  for(unsigned int l = k+1; l < genHeavy.size(); l++){  
	    if( deltaR(genHeavy[k], genHeavy[l])<0.5 ) overlapH++;
	  }	    
	}
	overlapHeavy_ = overlapH;
	vector<TLorentzVector> genLight;
	genLight.push_back( TOPHADW1);
	genLight.push_back( TOPHADW2);
	int overlapL = 0;
	for(unsigned int k = 0; k < genLight.size(); k++){  
	  for(unsigned int l = 0; l < genHeavy.size(); l++){  
	    if( deltaR(genLight[k], genHeavy[l])<0.5 ) overlapL++;
	  }	    
	}
	if(  deltaR(genLight[0], genLight[1])<0.5  ) overlapL++;
	overlapLight_  = overlapL;
	
	csvL_[0] = (mapFilt[w1].csv>0 ? mapFilt[w1].csv : 0.0 );
	//csvL_[1] = (mapFilt[w2].csv>0 ? mapFilt[w2].csv : 0.0 );
	ptL_[0]  = (mapFilt[w1].pt>0  ? mapFilt[w1].pt  : 0.0 );
	//ptL_[1]  = (mapFilt[w2].pt>0  ? mapFilt[w2].pt  : 0.0 );
	
	csvH_[0] = (mapFilt[b1].csv>0 ? mapFilt[b1].csv : 0.0 );
	csvH_[1] = (mapFilt[b2].csv>0 ? mapFilt[b2].csv : 0.0 );
	csvH_[2] = (mapFilt[b3].csv>0 ? mapFilt[b3].csv : 0.0 );
	csvH_[3] = (mapFilt[b4].csv>0 ? mapFilt[b4].csv : 0.0 );

	jets.clear();
	jets.push_back( leptonLV  );  
	jets.push_back( neutrinoLV  );  
	jets.push_back( myJetsFilt[b1]   );
	jets.push_back( myJetsFilt[w1]   );  
	jets.push_back( myJetsFilt[w1]   );  
	jets.push_back( myJetsFilt[b2]   );
	jets.push_back( myJetsFilt[b3]   );
	jets.push_back( myJetsFilt[b4]   );

	pos_to_index.clear();
	pos_to_index[2] = b1;
	pos_to_index[3] = w1;
	pos_to_index[4] = w1;
	pos_to_index[5] = b2;
	pos_to_index[6] = b3;
	pos_to_index[7] = b4;
	  
	meIntegrator->setJets(&jets);

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

	  //if( Bs.Pt()>20 && TMath::Abs(Bs.Eta())<2.5 ) nMatchSimBs_++;
	  //if     ( deltaR( Bs , myJetsFilt[b1] )<0.2 ) nMatchSimBs_++;
	  //else if( deltaR( Bs , myJetsFilt[b2] )<0.2 ) nMatchSimBs_++;
	  //else if( deltaR( Bs , myJetsFilt[b3] )<0.2 ) nMatchSimBs_++;
	  //else if( deltaR( Bs , myJetsFilt[b4] )<0.2 ) nMatchSimBs_++;
	  //else if( deltaR( Bs , myJetsFilt[w1] )<0.2 ) nMatchSimBs_++;
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
	
	type_ = 2;

	meIntegrator->setIntType( MEIntegratorNew::SL1wj );


	/////////////////////////////////////////////////////////////
	int index_bLep          =  mapFilt[b1].index;
	float csv_bLep_nominal  =  index_bLep>=0 ? hJet_csv_nominal[index_bLep]  : aJet_csv_nominal[-index_bLep-1];
	float csv_bLep_upBC     =  index_bLep>=0 ? hJet_csv_upBC   [index_bLep]  : aJet_csv_upBC   [-index_bLep-1];
	float csv_bLep_downBC   =  index_bLep>=0 ? hJet_csv_downBC [index_bLep]  : aJet_csv_downBC [-index_bLep-1];
	float csv_bLep_upL      =  index_bLep>=0 ? hJet_csv_upL    [index_bLep]  : aJet_csv_upL    [-index_bLep-1];
	float csv_bLep_downL    =  index_bLep>=0 ? hJet_csv_downL  [index_bLep]  : aJet_csv_downL  [-index_bLep-1];
	double csv_bLep = csv_bLep_nominal;
	if( doCSVup )   csv_bLep = TMath::Max(csv_bLep_upBC,   csv_bLep_upL);
	if( doCSVdown ) csv_bLep = TMath::Min(csv_bLep_downBC, csv_bLep_downL);   
	
	if     ( csv_bLep< 1.1 ){
	  float m1    = (meIntegrator->jetAt(2)+meIntegrator->jetAt(3)).M();
	  if(  m1>(MwL+5) && m1<(MwH-5) ){
	    flag_type2_ = 0; 
	    if( csv_bLep<0.95 ) flag_type2_ = 1; 
	    if( csv_bLep<0.90 ) flag_type2_ = 2; 
	    if( csv_bLep<0.85 ) flag_type2_ = 3;
	    if( csv_bLep<0.80 ) flag_type2_ = 4;  
	  }
	}
	int index_bHad          =  mapFilt[b2].index;
	float csv_bHad_nominal  =  index_bHad>=0 ? hJet_csv_nominal[index_bHad]  : aJet_csv_nominal[-index_bHad-1];
	float csv_bHad_upBC     =  index_bHad>=0 ? hJet_csv_upBC   [index_bHad]  : aJet_csv_upBC   [-index_bHad-1];
	float csv_bHad_downBC   =  index_bHad>=0 ? hJet_csv_downBC [index_bHad]  : aJet_csv_downBC [-index_bHad-1];
	float csv_bHad_upL      =  index_bHad>=0 ? hJet_csv_upL    [index_bHad]  : aJet_csv_upL    [-index_bHad-1];
	float csv_bHad_downL    =  index_bHad>=0 ? hJet_csv_downL  [index_bHad]  : aJet_csv_downL  [-index_bHad-1];	  
	double csv_bHad = csv_bHad_nominal;
	if( doCSVup )   csv_bHad = TMath::Max(csv_bHad_upBC,   csv_bHad_upL);
	if( doCSVdown ) csv_bHad = TMath::Min(csv_bHad_downBC, csv_bHad_downL);
	if     ( csv_bHad< 1.1 ){
	  float m1    = (meIntegrator->jetAt(5)+meIntegrator->jetAt(3)).M();
	  if( m1>(MwL+5) && m1<(MwH-5) ){
	    flag_type2_ = 0; 
	    if( csv_bHad<0.95 ) flag_type2_ = 1; 
	    if( csv_bHad<0.90 ) flag_type2_ = 2; 
	    if( csv_bHad<0.85 ) flag_type2_ = 3;
	    if( csv_bHad<0.80 ) flag_type2_ = 4;  
	  }	
	}
	int index_b1          =  mapFilt[b3].index;
	float csv_b1_nominal  =  index_b1>=0 ? hJet_csv_nominal[index_b1]  : aJet_csv_nominal[-index_b1-1];
	float csv_b1_upBC     =  index_b1>=0 ? hJet_csv_upBC   [index_b1]  : aJet_csv_upBC   [-index_b1-1];
	float csv_b1_downBC   =  index_b1>=0 ? hJet_csv_downBC [index_b1]  : aJet_csv_downBC [-index_b1-1];
	float csv_b1_upL      =  index_b1>=0 ? hJet_csv_upL    [index_b1]  : aJet_csv_upL    [-index_b1-1];
	float csv_b1_downL    =  index_b1>=0 ? hJet_csv_downL  [index_b1]  : aJet_csv_downL  [-index_b1-1];
	double csv_b1 = csv_b1_nominal;
	if( doCSVup )   csv_b1 = TMath::Max(csv_b1_upBC,   csv_b1_upL);
	if( doCSVdown ) csv_b1 = TMath::Min(csv_b1_downBC, csv_b1_downL);	   
	if     ( csv_b1< 1.1 ){
	  float m1    = (meIntegrator->jetAt(6)+meIntegrator->jetAt(3)).M();
	  if( m1>(MwL+5) && m1<(MwH-5)){
	    flag_type2_ = 0; 
	    if( csv_b1<0.95 ) flag_type2_ = 1; 
	    if( csv_b1<0.90 ) flag_type2_ = 2; 
	    if( csv_b1<0.85 ) flag_type2_ = 3;
	    if( csv_b1<0.80 ) flag_type2_ = 4; 
	  }
	}
	int index_b2          =  mapFilt[b4].index;
	float csv_b2_nominal  =  index_b2>=0 ? hJet_csv_nominal[index_b2]  : aJet_csv_nominal[-index_b2-1];
	float csv_b2_upBC     =  index_b2>=0 ? hJet_csv_upBC   [index_b2]  : aJet_csv_upBC   [-index_b2-1];
	float csv_b2_downBC   =  index_b2>=0 ? hJet_csv_downBC [index_b2]  : aJet_csv_downBC [-index_b2-1];
	float csv_b2_upL      =  index_b2>=0 ? hJet_csv_upL    [index_b2]  : aJet_csv_upL    [-index_b2-1];
	float csv_b2_downL    =  index_b2>=0 ? hJet_csv_downL  [index_b2]  : aJet_csv_downL  [-index_b2-1];
	double csv_b2 = csv_b2_nominal;
	if( doCSVup )   csv_b2 = TMath::Max(csv_b2_upBC,   csv_b2_upL);
	if( doCSVdown ) csv_b2 = TMath::Min(csv_b2_downBC, csv_b2_downL);	  
	if     ( csv_b2< 1.1 ){
	  float m1    = (meIntegrator->jetAt(7)+meIntegrator->jetAt(3)).M();
	  if( m1>(MwL+5) && m1<(MwH-5) ){
	    flag_type2_ = 0; 
	    if( csv_b2<0.95 ) flag_type2_ = 1; 
	    if( csv_b2<0.90 ) flag_type2_ = 2; 
	    if( csv_b2<0.85 ) flag_type2_ = 3;
	    if( csv_b2<0.80 ) flag_type2_ = 4; 
	  }
	}
	/////////////////////////////////////////////////////////////


	meIntegrator->setSumEt( METtype1p2corr.sumet );
	meIntegrator->setMEtCov(-99,-99,0);
	meIntegrator->setTopFlags( vLepton_charge[0]==1 ? +1 : -1 , vLepton_charge[0]==1 ? -1 : +1 );
	  
	clock->Start();

	for(int m = 0; m < nMassPoints ; m++){
	  meIntegrator->setMass( mH[m] );
	    
	  double maxP_s = 0.;
	  double maxP_b = 0.;

	  for(unsigned int pos = 0; pos < 12 ; pos++){

	    meIntegrator->initVersors( permutations_SL1wj[pos] );
	    double mass, massLow, massHigh;
	    bool skip = !(meIntegrator->compatibilityCheck(0.95, /*printP4*/ 0, mass, massLow, massHigh )) ;


	    pair<double, double> range_x0 = (meIntegrator->getW1JetEnergyCI(0.95));
	    pair<double, double> range_x1 =  make_pair(-1,1);
	    pair<double, double> range_x2 =  make_pair(-PI,PI);	    
	    pair<double, double> range_x3 =  make_pair(-1,1);
	    pair<double, double> range_x4 = useMET ? (meIntegrator->getNuPhiCI(0.95)) : make_pair(-PI,PI);
	    pair<double, double> range_x5 = (meIntegrator->getB1EnergyCI(0.95));
	    pair<double, double> range_x6 = (meIntegrator->getB2EnergyCI(0.95));

	    double x0L = range_x0.first;
	    double x0U = range_x0.second;
	    double x1L = range_x1.first;
	    double x1U = range_x1.second;
	    double x2L = range_x2.first;
	    double x2U = range_x2.second;
	    double x3L = range_x3.first;
	    double x3U = range_x3.second;
	    double x4L = range_x4.first;
	    double x4U = range_x4.second;
	    double x5L = range_x5.first;
	    double x5U = range_x5.second;
	    double x6L = range_x6.first;
	    double x6U = range_x6.second;
	
	    double xLmode1[6]   = {x0L, x1L, x2L, x3L, x4L, x5L};
	    double xUmode1[6]   = {x0U, x1U, x2U, x3U, x4U, x5U};
	    double xLmode1_b[7] = {x0L, x1L, x2L, x3L, x4L, x5L, x6L};
	    double xUmode1_b[7] = {x0U, x1U, x2U, x3U, x4U, x5U, x6U};	     	     
	      
	    double p = 0.;
	    if(speedup==0){
	      
	      for(int hyp = 0 ; hyp<2; hyp++){
		
		if(SoB==0 && hyp!=hypo) continue;
		if(SoB==1 && hyp==0 && skip) continue;

		cout << "Testing hypothesis " << (hyp==0 ? "S" : "B") << " for permutation " << pos << endl; 
		meIntegrator->setHypo(hyp);

		int ntries = 0;
		int intPoints = 4000;		    
		while( ntries < MAX_REEVAL_TRIES){
		  
		  ROOT::Math::Functor toIntegrate(meIntegrator, &MEIntegratorNew::Eval,  6+hyp);
		  ROOT::Math::GSLMCIntegrator ig2( ROOT::Math::IntegrationMultiDim::kVEGAS , 1.e-12, 1.e-5, intPoints);
		  ig2.SetFunction(toIntegrate);
		  meIntegrator->SetPar(  6+hyp );

		  p = (hyp==0 ? ig2.Integral(xLmode1, xUmode1) : ig2.Integral(xLmode1_b, xUmode1_b));

		  double chi2 =  ig2.ChiSqr();
		  if( hyp==0 ){
		    if( p>maxP_s ) maxP_s = p;		  
		    else if( p<0.1*maxP_s) ntries = MAX_REEVAL_TRIES+1;
		  }
		  else{
		    if( p>maxP_b ) maxP_b = p;		  
		    else if( p<0.1*maxP_b) ntries = MAX_REEVAL_TRIES+1;	
		  }

		  if( chi2 > maxChi2_ ){
		    ntries++;
		    intPoints *= 1.5;
		  }
		  else ntries = MAX_REEVAL_TRIES+1;

		}
	    
		if( TMath::IsNaN(p) ) p = 0.;
		

		if( useBtag ){
		  
		  int bLep_pos    = (permutations_SL1wj[pos])/100000;
		  int index_bLep  =  mapFilt[pos_to_index[bLep_pos]].index;
		  float csv_bLep_nominal  =  index_bLep>=0 ? hJet_csv_nominal[index_bLep]  : aJet_csv_nominal[-index_bLep-1];
		  float csv_bLep_upBC     =  index_bLep>=0 ? hJet_csv_upBC   [index_bLep]  : aJet_csv_upBC   [-index_bLep-1];
		  float csv_bLep_downBC   =  index_bLep>=0 ? hJet_csv_downBC [index_bLep]  : aJet_csv_downBC [-index_bLep-1];
		  float csv_bLep_upL      =  index_bLep>=0 ? hJet_csv_upL    [index_bLep]  : aJet_csv_upL    [-index_bLep-1];
		  float csv_bLep_downL    =  index_bLep>=0 ? hJet_csv_downL  [index_bLep]  : aJet_csv_downL  [-index_bLep-1];
		  double csv_bLep = csv_bLep_nominal;
		  if( doCSVup )   csv_bLep = TMath::Max(csv_bLep_upBC,   csv_bLep_upL);
		  if( doCSVdown ) csv_bLep = TMath::Min(csv_bLep_downBC, csv_bLep_downL);
		      
		  int bHad_pos    = (permutations_SL1wj[pos])%1000/100;
		  int index_bHad  =  mapFilt[pos_to_index[bHad_pos]].index;
		  //float csv_bHad  =  index_bHad>=0 ? hJet_csv_nominal[index_bHad]  : aJet_csv_nominal[-index_bHad-1];
		  float csv_bHad_nominal  =  index_bHad>=0 ? hJet_csv_nominal[index_bHad]  : aJet_csv_nominal[-index_bHad-1];
		  float csv_bHad_upBC     =  index_bHad>=0 ? hJet_csv_upBC   [index_bHad]  : aJet_csv_upBC   [-index_bHad-1];
		  float csv_bHad_downBC   =  index_bHad>=0 ? hJet_csv_downBC [index_bHad]  : aJet_csv_downBC [-index_bHad-1];
		  float csv_bHad_upL      =  index_bHad>=0 ? hJet_csv_upL    [index_bHad]  : aJet_csv_upL    [-index_bHad-1];
		  float csv_bHad_downL    =  index_bHad>=0 ? hJet_csv_downL  [index_bHad]  : aJet_csv_downL  [-index_bHad-1];
		  double csv_bHad = csv_bHad_nominal;
		  if( doCSVup )   csv_bHad = TMath::Max(csv_bHad_upBC,   csv_bHad_upL);
		  if( doCSVdown ) csv_bHad = TMath::Min(csv_bHad_downBC, csv_bHad_downL);
		  
		  int b1_pos    =   (permutations_SL1wj[pos])%100/10;
		  int index_b1  =  mapFilt[pos_to_index[b1_pos]].index;
		  //float csv_b1  =  index_b1>=0 ? hJet_csv_nominal[index_b1]  : aJet_csv_nominal[-index_b1-1];
		  float csv_b1_nominal  =  index_b1>=0 ? hJet_csv_nominal[index_b1]  : aJet_csv_nominal[-index_b1-1];
		  float csv_b1_upBC     =  index_b1>=0 ? hJet_csv_upBC   [index_b1]  : aJet_csv_upBC   [-index_b1-1];
		  float csv_b1_downBC   =  index_b1>=0 ? hJet_csv_downBC [index_b1]  : aJet_csv_downBC [-index_b1-1];
		  float csv_b1_upL      =  index_b1>=0 ? hJet_csv_upL    [index_b1]  : aJet_csv_upL    [-index_b1-1];
		  float csv_b1_downL    =  index_b1>=0 ? hJet_csv_downL  [index_b1]  : aJet_csv_downL  [-index_b1-1];
		  double csv_b1 = csv_b1_nominal;
		  if( doCSVup )   csv_b1 = TMath::Max(csv_b1_upBC,   csv_b1_upL);
		  if( doCSVdown ) csv_b1 = TMath::Min(csv_b1_downBC, csv_b1_downL);
		  
		  int b2_pos    = (permutations_SL1wj[pos])%10;
		  int index_b2  =  mapFilt[pos_to_index[b2_pos]].index;
		  //float csv_b2  =  index_b2>=0 ? hJet_csv_nominal[index_b2]  : aJet_csv_nominal[-index_b2-1];
		  float csv_b2_nominal  =  index_b2>=0 ? hJet_csv_nominal[index_b2]  : aJet_csv_nominal[-index_b2-1];
		  float csv_b2_upBC     =  index_b2>=0 ? hJet_csv_upBC   [index_b2]  : aJet_csv_upBC   [-index_b2-1];
		  float csv_b2_downBC   =  index_b2>=0 ? hJet_csv_downBC [index_b2]  : aJet_csv_downBC [-index_b2-1];
		  float csv_b2_upL      =  index_b2>=0 ? hJet_csv_upL    [index_b2]  : aJet_csv_upL    [-index_b2-1];
		  float csv_b2_downL    =  index_b2>=0 ? hJet_csv_downL  [index_b2]  : aJet_csv_downL  [-index_b2-1];
		  double csv_b2 = csv_b2_nominal;
		  if( doCSVup )   csv_b2 = TMath::Max(csv_b2_upBC,   csv_b2_upL);
		  if( doCSVdown ) csv_b2 = TMath::Min(csv_b2_downBC, csv_b2_downL);
		  		  
		  
		  string bin_bLep = TMath::Abs( (meIntegrator->jetAt(2)).Eta() )<1.0 ? "Bin0" : "Bin1"; 
		  double p_b_bLep =  btagger["b_"+bin_bLep]!=0 ?  btagger["b_"+bin_bLep]->GetBinContent( btagger["b_"+bin_bLep]->FindBin( TMath::Min(csv_bLep, 0.999999) ) ) : 1.;
		  
		  string bin_bHad = TMath::Abs( (meIntegrator->jetAt(5)).Eta() )<1.0 ? "Bin0" : "Bin1"; 
		  double p_b_bHad =  btagger["b_"+bin_bHad]!=0 ?  btagger["b_"+bin_bHad]->GetBinContent( btagger["b_"+bin_bHad]->FindBin( TMath::Min(csv_bHad, 0.999999) ) ) : 1.;
		  
		  string bin_b1   = TMath::Abs( (meIntegrator->jetAt(6)).Eta() )<1.0 ? "Bin0" : "Bin1"; 
		  double p_b_b1   =  btagger["b_"+bin_b1]!=0 ?  btagger["b_"+bin_b1]->GetBinContent( btagger["b_"+bin_b1]->FindBin( TMath::Min(csv_b1, 0.999999) ) ) : 1.;
		  double p_j_b1   =  btagger["l_"+bin_b1]!=0 ?  btagger["l_"+bin_b1]->GetBinContent( btagger["l_"+bin_b1]->FindBin( TMath::Min(csv_b1, 0.999999) ) ) : 1.;
		  
		  string bin_b2   = TMath::Abs( (meIntegrator->jetAt(7)).Eta() )<1.0 ? "Bin0" : "Bin1"; 
		  double p_b_b2   =  btagger["b_"+bin_b2]!=0 ?  btagger["b_"+bin_b2]->GetBinContent( btagger["b_"+bin_b2]->FindBin( TMath::Min(csv_b2, 0.999999) ) ) : 1.;
		  double p_j_b2   =  btagger["l_"+bin_b2]!=0 ?  btagger["l_"+bin_b2]->GetBinContent( btagger["l_"+bin_b2]->FindBin( TMath::Min(csv_b2, 0.999999) ) ) : 1.;
		 
		  
		  //cout << "bLep_pos " << bLep_pos << " -- correspond to " << index_bLep << " with csv="<< csv_bLep << ": prob under b hyp: " << p_b_bLep << endl;
		  //cout << "bHad_pos " << bHad_pos << " -- correspond to " << index_bHad << " with csv="<< csv_bHad << ": prob under b hyp: " << p_b_bHad <<  endl;
		  //cout << "b1_pos " << b1_pos << " -- correspond to " << index_b1 << " with csv="<< csv_b1         << ": prob under b hyp: " << p_b_b1 << ", under j hyp: " << p_j_b1 <<  endl;
		  //cout << "b2_pos " << b2_pos << " -- correspond to " << index_b2 << " with csv="<< csv_b2         << ": prob under b hyp: " << p_b_b2 << ", under j hyp: " << p_j_b2 << endl;		      
		  
		  /*
		  int sv_bLep = 0;
		  int sv_bHad = 0;
		  int sv_b1   = 0;
		  int sv_b2   = 0;
		  for(int l = 0; l < nSvs; l++){
		    TLorentzVector vertex(1,0,0,1);
		    vertex.SetPtEtaPhiM( Svpt[l], Sveta[l], Svphi[l], SvmassSv[l]);
		    if(deltaR(vertex, meIntegrator->jetAt(2) )<0.5)      sv_bLep = 1;
		    else if(deltaR(vertex, meIntegrator->jetAt(5) )<0.5) sv_bHad = 1;
		    else if(deltaR(vertex, meIntegrator->jetAt(6) )<0.5) sv_b1   = 1;
		    else if(deltaR(vertex, meIntegrator->jetAt(7) )<0.5) sv_b2   = 1;
		    else{}
		  }
		  double p_b_sv = svs[sv_bLep]*svs[sv_bHad]*svs[sv_b1]  *svs[sv_b2];
		  double p_j_sv = svs[sv_bLep]*svs[sv_bHad]*svs[2+sv_b1]*svs[2+sv_b2];
		  */

		  if( mH[m]<met+0.5 && mH[m]>met-0.5){
		    if(hyp==0){
		      probAtSgn_ttbb_         += ( p * p_b_bLep * p_b_bHad * p_b_b1 * p_b_b2 /**  p_b_sv*/);
		    }
		    else{
		      probAtSgn_alt_ttbb_     += ( p * p_b_bLep * p_b_bHad * p_b_b1 * p_b_b2 /**  p_b_sv*/);
		      probAtSgn_alt_ttjj_     += ( p * p_b_bLep * p_b_bHad * p_j_b1 * p_j_b2 /**  p_j_sv*/);
		    }
		  }
		  
		}

		if( p>= maxP_s && hyp==0 ){
		  int lastTwoDigits = (permutations_SL1wj[pos])%100;
		  int higgs1 = lastTwoDigits/(10);
		  int higgs2 = lastTwoDigits%10;		    
		  //cout << "pos=" << permutations_SL2wj[pos] << " => " << higgs1 << "," << higgs2 << endl;		    
		  if( higgs1 == 2)
		    higgs1 = b1;		      
		  else if(  higgs1 == 5 )
		    higgs1 = b2;
		  else if(  higgs1 == 6 )
		    higgs1 = b3;
		  else if(  higgs1 == 7 )
		    higgs1 = b4;
		  else{ cout << "What happened?" << endl; }
		  if( higgs2 == 2)
		    higgs2 = b1;		      
		  else if(  higgs2 == 5 )
		    higgs2 = b2;
		  else if(  higgs2 == 6 )
		    higgs2 = b3;
		  else if(  higgs2 == 7 )
		    higgs2 = b4;
		  else{ cout << "What happened?" << endl; }		    
		  //cout << "Corrspond to " << higgs1 << "," << higgs2 << endl;		    
		  int index1 =  mapFilt[higgs1].index;		      
		  int index2 =  mapFilt[higgs2].index;		      		    
		  float csv_1 = index1>=0 ? hJet_csv_nominal[index1]  : aJet_csv_nominal[-index1-1];
		  float csv_2 = index2>=0 ? hJet_csv_nominal[index2]  : aJet_csv_nominal[-index2-1];		    
		  //cout << "Deug 1: " << mapFilt[higgs1].pt << " --- " << (meIntegrator->jetAt(6)).Pt() << endl;
		  //cout << "Deug 2: " << mapFilt[higgs2].pt << " --- " << (meIntegrator->jetAt(7)).Pt() << endl;		    
		  commitVariables( meIntegrator->jetAt(6), meIntegrator->jetAt(7), csv_1, csv_2,
				   dR_,dPhi_,dEta_,Mjj_,Pt_,Eta_,pt1_,pt2_,eta1_,eta2_,csv1_,csv2_);
		}

		if( mH[m]<met+0.5 && mH[m]>met-0.5){
		  if(hyp==0)  probAtSgn_     += p;
		  else        probAtSgn_alt_ += p;
		}

	      }
	    }
	  }
	}

	clock->Stop();
	time_ = clock->RealTime();
	clock->Reset();	 
	
      }
      else if(properEventSL && numJets30BtagM==3 && numJets30UntagM==2 && doType4){

	
	if( w1==999 || w2==999 || b1==999 || b2==999 || b3==999){
	  cout << "Inconsistentcy type 3 found!" << endl;
	  cout << w1 << ", " << w2 << "," << b1 << "," << b2 << "," << b3 << "," << b4 <<  endl;
	  continue;
	}

	float WMass  = (myJetsFilt[w1]+myJetsFilt[w2]).M();
	float Mass12 = (myJetsFilt[b1]+myJetsFilt[b2]).M();
	float Mass13 = (myJetsFilt[b1]+myJetsFilt[b3]).M();
	float Mass23 = (myJetsFilt[b2]+myJetsFilt[b3]).M();

	if( !( (Mass12>MhL && Mass12<MhH) || (Mass13>MhL && Mass13<MhH) || (Mass23>MhL && Mass23<MhH) ) ) continue;
	if( !(WMass>MwL && WMass<MwH) ) continue;

	counter++;
	if(fixNumberOfEventsPerJob && !(counter>=evLow && counter<=evHigh) ) continue;
	cout << "Processing event # " << counter << " (SL[3,2])" << endl;

	  

	int hMatches = 0;
	if     (  deltaR(myJetsFilt[b1],genBLV)<0.3 )    hMatches++;
	else if(  deltaR(myJetsFilt[b1],genBbarLV)<0.3)  hMatches++;
	if     (  deltaR(myJetsFilt[b2],genBLV)<0.3 )    hMatches++;
	else if(  deltaR(myJetsFilt[b2],genBbarLV)<0.3)  hMatches++;
	if     (  deltaR(myJetsFilt[b3],genBLV)<0.3 )    hMatches++;
	else if(  deltaR(myJetsFilt[b3],genBbarLV)<0.3)  hMatches++;
	matchesH_    = hMatches;
	matchesHAll_ = hMatches;
	if     (  deltaR(myJetsFilt[w1],genBLV)<0.3 )    matchesHAll_++;
	else if(  deltaR(myJetsFilt[w1],genBbarLV)<0.3)  matchesHAll_++;
	if     (  deltaR(myJetsFilt[w2],genBLV)<0.3 )    matchesHAll_++;
	else if(  deltaR(myJetsFilt[w2],genBbarLV)<0.3)  matchesHAll_++;
	
	int tMatches = 0;
	if     (  deltaR(myJetsFilt[b1],TOPHADB)<0.3 )   tMatches++;
	else if(  deltaR(myJetsFilt[b1],TOPLEPB)<0.3)    tMatches++;
	if     (  deltaR(myJetsFilt[b2],TOPHADB)<0.3 )   tMatches++;
	else if(  deltaR(myJetsFilt[b2],TOPLEPB)<0.3)    tMatches++;
	if     (  deltaR(myJetsFilt[b3],TOPHADB)<0.3 )   tMatches++;
	else if(  deltaR(myJetsFilt[b3],TOPLEPB)<0.3)    tMatches++;
	matchesT_    = tMatches;
	matchesTAll_ = tMatches;
	if     (  deltaR(myJetsFilt[w1],TOPHADB)<0.3 )   matchesTAll_++;
	else if(  deltaR(myJetsFilt[w1],TOPLEPB)<0.3)    matchesTAll_++;
	if     (  deltaR(myJetsFilt[w2],TOPHADB)<0.3 )   matchesTAll_++;
	else if(  deltaR(myJetsFilt[w2],TOPLEPB)<0.3)    matchesTAll_++;
	
	
	int wMatches = 0;
	if     (  deltaR(myJetsFilt[w1],TOPHADW1)<0.3 ){
	  wMatches++;
	  mL_[0] = 1;
	}
	else if(  deltaR(myJetsFilt[w1],TOPHADW2)<0.3){
	  wMatches++;
	  mL_[0] = 1;
	}
	else{ mL_[0] = 0; }
	if     (  deltaR(myJetsFilt[w2],TOPHADW1)<0.3 ){
	  wMatches++;
	  mL_[1] = 1;
	}
	else if(  deltaR(myJetsFilt[w2],TOPHADW2)<0.3){
	  wMatches++;
	  mL_[1] = 1;
	}
	else{  mL_[1] = 0; }
	matchesW_    = wMatches;
	matchesWAll_ = wMatches;	  
	if     (  deltaR(myJetsFilt[b1],TOPHADW1)<0.3 )   matchesWAll_++;
	else if(  deltaR(myJetsFilt[b1],TOPHADW2)<0.3)    matchesWAll_++;
	if     (  deltaR(myJetsFilt[b2],TOPHADW1)<0.3 )   matchesWAll_++;
	else if(  deltaR(myJetsFilt[b2],TOPHADW2)<0.3)    matchesWAll_++;
	if     (  deltaR(myJetsFilt[b3],TOPHADW1)<0.3 )   matchesWAll_++;
	else if(  deltaR(myJetsFilt[b3],TOPHADW2)<0.3)    matchesWAll_++;
	
	vector<TLorentzVector> genHeavy;
	genHeavy.push_back( genBLV);
	genHeavy.push_back( genBbarLV);
	genHeavy.push_back( TOPHADB);
	genHeavy.push_back( TOPLEPB);
	int overlapH = 0;
	for(unsigned int k = 0; k < genHeavy.size()-1; k++){  
	  for(unsigned int l = k+1; l < genHeavy.size(); l++){  
	    if( deltaR(genHeavy[k], genHeavy[l])<0.5 ) overlapH++;
	  }	    
	}
	overlapHeavy_ = overlapH;
	vector<TLorentzVector> genLight;
	genLight.push_back( TOPHADW1);
	genLight.push_back( TOPHADW2);
	int overlapL = 0;
	for(unsigned int k = 0; k < genLight.size(); k++){  
	  for(unsigned int l = 0; l < genHeavy.size(); l++){  
	    if( deltaR(genLight[k], genHeavy[l])<0.5 ) overlapL++;
	  }	    
	}
	if(  deltaR(genLight[0], genLight[1])<0.5  ) overlapL++;
	overlapLight_  = overlapL;

	csvL_[0] = (mapFilt[w1].csv>0 ? mapFilt[w1].csv : 0.0 );
	csvL_[1] = (mapFilt[w2].csv>0 ? mapFilt[w2].csv : 0.0 );
	ptL_[0]  = (mapFilt[w1].pt>0  ? mapFilt[w1].pt  : 0.0 );
	ptL_[1]  = (mapFilt[w2].pt>0  ? mapFilt[w2].pt  : 0.0 );
	
	csvH_[0] = (mapFilt[b1].csv>0 ? mapFilt[b1].csv : 0.0 );
	csvH_[1] = (mapFilt[b2].csv>0 ? mapFilt[b2].csv : 0.0 );
	csvH_[2] = (mapFilt[b3].csv>0 ? mapFilt[b3].csv : 0.0 );
	
	type_ = 4;

	meIntegrator->setSumEt( METtype1p2corr.sumet );
	meIntegrator->setMEtCov(-99,-99,0);
	meIntegrator->setTopFlags( vLepton_charge[0]==1 ? +1 : -1 , vLepton_charge[0]==1 ? -1 : +1 );
	

	for(int iter=0; iter<2; iter++){

	  cout << "Doing iter " << iter << endl;

	  if( iter==0 ){
	    jets.clear();
	    jets.push_back( leptonLV  );  
	    jets.push_back( neutrinoLV  );  
	    jets.push_back( myJetsFilt[b1]   );
	    jets.push_back( myJetsFilt[w1]   );
	    jets.push_back( myJetsFilt[w2]   );  
	    jets.push_back( myJetsFilt[b1]   );    // dummy - NOT CONSIDRED
	    jets.push_back( myJetsFilt[b2]   );
	    jets.push_back( myJetsFilt[b3]   );
	    meIntegrator->setJets(&jets);	
	    meIntegrator->setIntType( MEIntegratorNew::SLNoBHad );	    
	  }
	  else{
	    jets.clear();
	    jets.push_back( leptonLV  );  
	    jets.push_back( neutrinoLV  );  
	    jets.push_back( myJetsFilt[b1]   );    // dummy - NOT CONSIDRED
	    jets.push_back( myJetsFilt[w1]   );  
	    jets.push_back( myJetsFilt[w2]   );  
	    jets.push_back( myJetsFilt[b1]   );
	    jets.push_back( myJetsFilt[b2]   );
	    jets.push_back( myJetsFilt[b3]   );
	    meIntegrator->setJets(&jets);
	    meIntegrator->setIntType( MEIntegratorNew::SLNoBLep );	    
	  }


	  for(int m = 0; m < nMassPoints ; m++){
	    meIntegrator->setMass( mH[m] );
	      
	    double maxP_s = 0.;
	    double maxP_b = 0.;
	    
	    for(unsigned int pos = 0; pos < 3 ; pos++){

	      if(iter==0) 
		meIntegrator->initVersors( permutations_SLNoBHad[pos] );
	      else
		meIntegrator->initVersors( permutations_SLNoBLep[pos] );

	      double mass, massLow, massHigh;
	      bool skip = !(meIntegrator->compatibilityCheck(0.95, /*printP4*/ 0, mass, massLow, massHigh )) ;
		
	      pair<double, double> range_x0 = (meIntegrator->getW1JetEnergyCI(0.95));
	      pair<double, double> range_x1 =  make_pair(-1,1);
	      pair<double, double> range_x2 =  make_pair(-PI,PI);	    
	      pair<double, double> range_x3 =  make_pair(-1,1);
	      pair<double, double> range_x4 = useMET ? (meIntegrator->getNuPhiCI(0.95)) : make_pair(-PI,PI);
	      pair<double, double> range_x5 = (meIntegrator->getB1EnergyCI(0.95));
	      pair<double, double> range_x6 = (meIntegrator->getB2EnergyCI(0.95));
	      
	      double x0L = range_x0.first;
	      double x0U = range_x0.second;
	      double x1L = range_x1.first;
	      double x1U = range_x1.second;
	      double x2L = range_x2.first;
	      double x2U = range_x2.second;
	      double x3L = range_x3.first;
	      double x3U = range_x3.second;
	      double x4L = range_x4.first;
	      double x4U = range_x4.second;
	      double x5L = range_x5.first;
	      double x5U = range_x5.second;
	      double x6L = range_x6.first;
	      double x6U = range_x6.second;
		

	      double xLmode2[6]   = {x0L, x1L, x2L, x3L, x4L, x5L};
	      double xUmode2[6]   = {x0U, x1U, x2U, x3U, x4U, x5U};
	      double xLmode2_b[7] = {x0L, x1L, x2L, x3L, x4L, x5L, x6L};
	      double xUmode2_b[7] = {x0U, x1U, x2U, x3U, x4U, x5U, x6U};

	      clock->Start();
		
	      double p = 0.;

	      if(speedup==0){
		
		for(int hyp = 0 ; hyp<2; hyp++){
		  
		  if(SoB==0 && hyp!=hypo) continue;
		  if(SoB==1 && hyp==0 && skip) continue;
		  
		  cout << "Testing hypothesis " << (hyp==0 ? "S" : "B") << " for permutation " << pos << endl; 
		  meIntegrator->setHypo(hyp);
		  
		  int ntries = 0;
		  int intPoints = 5000;
		  while( ntries < MAX_REEVAL_TRIES){
		    
		    ROOT::Math::Functor toIntegrate(meIntegrator, &MEIntegratorNew::Eval, 6+hyp);
		    ROOT::Math::GSLMCIntegrator ig2( ROOT::Math::IntegrationMultiDim::kVEGAS , 1.e-12, 1.e-5, intPoints);
		    ig2.SetFunction(toIntegrate);
		    meIntegrator->SetPar( 6+hyp );	 
		    
		    p = (hyp==0 ? ig2.Integral(xLmode2, xUmode2) : ig2.Integral(xLmode2_b, xUmode2_b));

		    double chi2 =  ig2.ChiSqr();
		    if( hyp==0 ){
		      if( p>maxP_s ) maxP_s = p;		  
		      else if( p<0.1*maxP_s) ntries = MAX_REEVAL_TRIES+1;
		    }
		    else{
		      if( p>maxP_b ) maxP_b = p;		  
		      else if( p<0.1*maxP_b) ntries = MAX_REEVAL_TRIES+1;	
		    }
		    
		    if( chi2 > maxChi2_ ){
		      ntries++;
		      intPoints *= 1.5;
		    }
		    else ntries = MAX_REEVAL_TRIES+1;
		      
		  }		   
				    
		  if( TMath::IsNaN(p) ) p = 0.;
		  
		  if( mH[m]<met+0.5 && mH[m]>met-0.5){
		    if(hyp==0)  probAtSgn_     += p;
		    else        probAtSgn_alt_ += p;
		  }
		}		  
	      }

	      clock->Stop();
	      time_ = clock->RealTime();
	      clock->Reset();	 

		
	    }
	  }
	  
	}


	
      }
      else if(properEventSL && numJets30BtagM==3 && numJets30BtagL==4 && numJets30UntagL>=2 && false){

	w1 = wL1;
	w2 = wL2;
	w3 = wL3;
	b1 = bL1;
	b2 = bL2;
	b3 = bL3;
	b4 = bL4;

	if(numJets30UntagL==2 /*&& false*/){

	  counter++;      
	  if(fixNumberOfEventsPerJob && !(counter>=evLow && counter<=evHigh) ) continue;
	  cout << "Processing event # " << counter << endl;
	 
	  if( w1==999 || w2==999 || b1==999 || b2==999 || b3==999 || b4==999){
	    cout << "Inconsistentcy type 0/1 found!" << endl;
	    cout << w1 << ", " << w2 << "," << b1 << "," << b2 << "," << b3 << "," << b4 <<  endl;
	    continue;
	  }

	  int hMatches = 0;
	  if     (  deltaR(myJetsFilt[b1],genBLV)<0.3 )    hMatches++;
	  else if(  deltaR(myJetsFilt[b1],genBbarLV)<0.3)  hMatches++;
	  if     (  deltaR(myJetsFilt[b2],genBLV)<0.3 )    hMatches++;
	  else if(  deltaR(myJetsFilt[b2],genBbarLV)<0.3)  hMatches++;
	  if     (  deltaR(myJetsFilt[b3],genBLV)<0.3 )    hMatches++;
	  else if(  deltaR(myJetsFilt[b3],genBbarLV)<0.3)  hMatches++;
	  if     (  deltaR(myJetsFilt[b4],genBLV)<0.3 )    hMatches++;
	  else if(  deltaR(myJetsFilt[b4],genBbarLV)<0.3)  hMatches++;
	  matchesH_    = hMatches;
	  matchesHAll_ = hMatches;
	  if     (  deltaR(myJetsFilt[w1],genBLV)<0.3 )    matchesHAll_++;
	  else if(  deltaR(myJetsFilt[w1],genBbarLV)<0.3)  matchesHAll_++;
	  if     (  deltaR(myJetsFilt[w2],genBLV)<0.3 )    matchesHAll_++;
	  else if(  deltaR(myJetsFilt[w2],genBbarLV)<0.3)  matchesHAll_++;

	  int tMatches = 0;
	  if     (  deltaR(myJetsFilt[b1],TOPHADB)<0.3 )   tMatches++;
	  else if(  deltaR(myJetsFilt[b1],TOPLEPB)<0.3)    tMatches++;
	  if     (  deltaR(myJetsFilt[b2],TOPHADB)<0.3 )   tMatches++;
	  else if(  deltaR(myJetsFilt[b2],TOPLEPB)<0.3)    tMatches++;
	  if     (  deltaR(myJetsFilt[b3],TOPHADB)<0.3 )   tMatches++;
	  else if(  deltaR(myJetsFilt[b3],TOPLEPB)<0.3)    tMatches++;
	  if     (  deltaR(myJetsFilt[b4],TOPHADB)<0.3 )   tMatches++;
	  else if(  deltaR(myJetsFilt[b4],TOPLEPB)<0.3)    tMatches++;
	  matchesT_    = tMatches;
	  matchesTAll_ = tMatches;
	  if     (  deltaR(myJetsFilt[w1],TOPHADB)<0.3 )   matchesTAll_++;
	  else if(  deltaR(myJetsFilt[w1],TOPLEPB)<0.3)    matchesTAll_++;
	  if     (  deltaR(myJetsFilt[w2],TOPHADB)<0.3 )   matchesTAll_++;
	  else if(  deltaR(myJetsFilt[w2],TOPLEPB)<0.3)    matchesTAll_++;


	  int wMatches = 0;
	  if     (  deltaR(myJetsFilt[w1],TOPHADW1)<0.3 ){
	    wMatches++;
	    mL_[0] = 1;
	  }
	  else if(  deltaR(myJetsFilt[w1],TOPHADW2)<0.3){
	    wMatches++;
	    mL_[0] = 1;
	  }
	  else{ mL_[0] = 0; }
	  if     (  deltaR(myJetsFilt[w2],TOPHADW1)<0.3 ){
	    wMatches++;
	    mL_[1] = 1;
	  }
	  else if(  deltaR(myJetsFilt[w2],TOPHADW2)<0.3){
	    wMatches++;
	    mL_[1] = 1;
	  }
	  else{  mL_[1] = 0; }
	  matchesW_    = wMatches;
	  matchesWAll_ = wMatches;	  
	  if     (  deltaR(myJetsFilt[b1],TOPHADW1)<0.3 )   matchesWAll_++;
	  else if(  deltaR(myJetsFilt[b1],TOPHADW2)<0.3)    matchesWAll_++;
	  if     (  deltaR(myJetsFilt[b2],TOPHADW1)<0.3 )   matchesWAll_++;
	  else if(  deltaR(myJetsFilt[b2],TOPHADW2)<0.3)    matchesWAll_++;
	  if     (  deltaR(myJetsFilt[b3],TOPHADW1)<0.3 )   matchesWAll_++;
	  else if(  deltaR(myJetsFilt[b3],TOPHADW2)<0.3)    matchesWAll_++;
	  if     (  deltaR(myJetsFilt[b4],TOPHADW1)<0.3 )   matchesWAll_++;
	  else if(  deltaR(myJetsFilt[b4],TOPHADW2)<0.3)    matchesWAll_++;

	  vector<TLorentzVector> genHeavy;
	  genHeavy.push_back( genBLV);
	  genHeavy.push_back( genBbarLV);
	  genHeavy.push_back( TOPHADB);
	  genHeavy.push_back( TOPLEPB);
	  int overlapH = 0;
	  for(unsigned int k = 0; k < genHeavy.size()-1; k++){  
	    for(unsigned int l = k+1; l < genHeavy.size(); l++){  
	      if( deltaR(genHeavy[k], genHeavy[l])<0.5 ) overlapH++;
	    }	    
	  }
	  overlapHeavy_ = overlapH;
	  vector<TLorentzVector> genLight;
	  genLight.push_back( TOPHADW1);
	  genLight.push_back( TOPHADW2);
	  int overlapL = 0;
	  for(unsigned int k = 0; k < genLight.size(); k++){  
	    for(unsigned int l = 0; l < genHeavy.size(); l++){  
	      if( deltaR(genLight[k], genHeavy[l])<0.5 ) overlapL++;
	    }	    
	  }
	  if(  deltaR(genLight[0], genLight[1])<0.5  ) overlapL++;
	  overlapLight_  = overlapL;

	  csvL_[0] = (mapFilt[w1].csv>0 ? mapFilt[w1].csv : 0.0 );
	  csvL_[1] = (mapFilt[w2].csv>0 ? mapFilt[w2].csv : 0.0 );
	  ptL_[0]  = (mapFilt[w1].pt>0  ? mapFilt[w1].pt  : 0.0 );
	  ptL_[1]  = (mapFilt[w2].pt>0  ? mapFilt[w2].pt  : 0.0 );

	  csvH_[0] = (mapFilt[b1].csv>0 ? mapFilt[b1].csv : 0.0 );
	  csvH_[1] = (mapFilt[b2].csv>0 ? mapFilt[b2].csv : 0.0 );
	  csvH_[2] = (mapFilt[b3].csv>0 ? mapFilt[b3].csv : 0.0 );
	  csvH_[3] = (mapFilt[b4].csv>0 ? mapFilt[b4].csv : 0.0 );


	  float WMass = (myJetsFilt[w1]+myJetsFilt[w2]).M();

	  jets.clear();
	  jets.push_back( leptonLV  );  
	  jets.push_back( neutrinoLV  );  
	  jets.push_back( myJetsFilt[b1]   );
	  jets.push_back( myJetsFilt[w1]   );  
	  jets.push_back( myJetsFilt[w2]   );  
	  jets.push_back( myJetsFilt[b2]   );
	  jets.push_back( myJetsFilt[b3]   );
	  jets.push_back( myJetsFilt[b4]   );
	  
	  meIntegrator->setJets(&jets);

	  if( WMass>MwL && WMass<MwH ){
	    cout << "W mass = " << WMass << " ... OK: SL(4,2)" << endl;
	    type_ = -1;
	    meIntegrator->setIntType( MEIntegratorNew::SL2wj );	    
	  }
	  else{
	    cout << "W mass = " << WMass << " ... NO: SL(4,1)" << endl;
	    type_ = -2;
	    meIntegrator->setIntType( MEIntegratorNew::SL1wj );
	  }
	  
	  meIntegrator->setSumEt( METtype1p2corr.sumet );
	  meIntegrator->setMEtCov(-99,-99,0);
	  meIntegrator->setTopFlags( vLepton_charge[0]==1 ? +1 : -1 , vLepton_charge[0]==1 ? -1 : +1 );
	  
	  for(int iter=0; iter<2; iter++){

	    if( type_ == -1 && iter==1) continue;
	    if( type_ == -2 && iter==1){
	      jets.clear();
	      jets.push_back( leptonLV  );  
	      jets.push_back( neutrinoLV  );  
	      jets.push_back( myJetsFilt[b1]   );
	      jets.push_back( myJetsFilt[w2]   );   // swap
	      jets.push_back( myJetsFilt[w1]   );  
	      jets.push_back( myJetsFilt[b2]   );
	      jets.push_back( myJetsFilt[b3]   );
	      jets.push_back( myJetsFilt[b4]   );
	      meIntegrator->setJets(&jets);	
	    }
	    cout << "Doing iteration " << iter 
		 << ": W jet (pt,eta,phi) = (" << (iter==0 ? myJetsFilt[w1].Pt() :  myJetsFilt[w2].Pt() )  
		 << "," << (iter==0 ? myJetsFilt[w1].Eta() :  myJetsFilt[w2].Eta() ) 
		 << "," << (iter==0 ? myJetsFilt[w1].Phi() :  myJetsFilt[w2].Phi() ) 
		 << ")" << endl; 


	    for(int m = 0; m < nMassPoints ; m++){
	      meIntegrator->setMass( mH[m] );
	      
	      double maxP_s = 0.;
	      double maxP_b = 0.;

	      for(unsigned int pos = 0; pos < 12 ; pos++){

		meIntegrator->initVersors( permutations_SL2wj[pos] );
		double mass, massLow, massHigh;
		bool skip = !(meIntegrator->compatibilityCheck(0.95, /*printP4*/ 0, mass, massLow, massHigh )) ;
		bool skip_WHad   = !( meIntegrator->compatibilityCheck_WHad  (0.98, /*printP4*/ 0, mass, massLow, massHigh ) );
		bool skip_TopHad = !( meIntegrator->compatibilityCheck_TopHad(0.98, /*printP4*/ 0, mass, massLow, massHigh ) );
		if( type_==-1 && (skip_WHad || skip_TopHad) ) continue;
		
		pair<double, double> range_x0 = (meIntegrator->getW1JetEnergyCI(0.95));
		pair<double, double> range_x1 =  make_pair(-1,1);
		pair<double, double> range_x2 =  make_pair(-PI,PI);	    
		pair<double, double> range_x3 =  make_pair(-1,1);
		pair<double, double> range_x4 = useMET ? (meIntegrator->getNuPhiCI(0.95)) : make_pair(-PI,PI);
		pair<double, double> range_x5 = (meIntegrator->getB1EnergyCI(0.95));
		pair<double, double> range_x6 = (meIntegrator->getB2EnergyCI(0.95));
		
		double x0L = range_x0.first;
		double x0U = range_x0.second;
		double x1L = range_x1.first;
		double x1U = range_x1.second;
		double x2L = range_x2.first;
		double x2U = range_x2.second;
		double x3L = range_x3.first;
		double x3U = range_x3.second;
		double x4L = range_x4.first;
		double x4U = range_x4.second;
		double x5L = range_x5.first;
		double x5U = range_x5.second;
		double x6L = range_x6.first;
		double x6U = range_x6.second;
	
		double xLmode0[4]   = {x0L, x3L, x4L, x5L};
		double xUmode0[4]   = {x0U, x3U, x4U, x5U};
		double xLmode0_b[5] = {x0L, x3L, x4L, x5L, x6L};
		double xUmode0_b[5] = {x0U, x3U, x4U, x5U, x6U};
	  	
		double xLmode1[6]   = {x0L, x1L, x2L, x3L, x4L, x5L};
		double xUmode1[6]   = {x0U, x1U, x2U, x3U, x4U, x5U};
		double xLmode1_b[7] = {x0L, x1L, x2L, x3L, x4L, x5L, x6L};
		double xUmode1_b[7] = {x0U, x1U, x2U, x3U, x4U, x5U, x6U};

		
		clock->Start();
		
		double p = 0.;

		if(speedup==0){
		  
		  for(int hyp = 0 ; hyp<2; hyp++){

		    if( SoB==0 && hyp!=hypo) continue;
		    if( SoB==1 && hyp==0 && skip) continue;

		    cout << "Testing hypothesis " << (hyp==0 ? "S" : "B") << " for permutation " << pos << endl; 
		    meIntegrator->setHypo(hyp);

		    int ntries = 0;
		    int intPoints = type_==-1 ? 2000 : 4000;
		    while( ntries < MAX_REEVAL_TRIES){

		      ROOT::Math::Functor toIntegrate(meIntegrator, &MEIntegratorNew::Eval, (type_ == -1 ? 4+hyp : 6+hyp));
		      ROOT::Math::GSLMCIntegrator ig2( ROOT::Math::IntegrationMultiDim::kVEGAS , 1.e-12, 1.e-5, intPoints);
		      ig2.SetFunction(toIntegrate);
		      meIntegrator->SetPar( (type_ == -1 ? 4+hyp : 6+hyp) );	 

		      if( type_== -1 ){
			p = (hyp==0 ? ig2.Integral(xLmode0, xUmode0) : ig2.Integral(xLmode0_b, xUmode0_b));
		      }
		      else{
			p = (hyp==0 ? ig2.Integral(xLmode1, xUmode1) : ig2.Integral(xLmode1_b, xUmode1_b));
		      }

		      double chi2 =  ig2.ChiSqr();
		      if( hyp==0 ){
			if( p>maxP_s ) maxP_s = p;		  
			else if( p<0.1*maxP_s) ntries = MAX_REEVAL_TRIES+1;
		      }
		      else{
			if( p>maxP_b ) maxP_b = p;		  
			else if( p<0.1*maxP_b) ntries = MAX_REEVAL_TRIES+1;	
		      }

		      if( chi2 > maxChi2_ ){
			ntries++;
			intPoints *= 1.5;
		      }
		      else ntries = MAX_REEVAL_TRIES+1;
		      
		    }		   
				    
		    clock->Stop();
		
		    if( TMath::IsNaN(p) ) p = 0.;
		    clock->Reset();
		  
		    if( mH[m]<met+0.5 && mH[m]>met-0.5){
		      if(hyp==0)  probAtSgn_     += p;
		      else        probAtSgn_alt_ += p;
		    }
		  }		  
		}		
	      }
	    }
	  }	  
	}
	else if(numJets30UntagL==3 /*&& false*/){
	  
	  //continue;

	  counter++;      
	  if(fixNumberOfEventsPerJob && !(counter>=evLow && counter<=evHigh) ) continue;
	  cout << "Processing event # " << counter << endl;
	 
	  
	  if( w1==999 || w2==999 || w3==999 || b1==999 || b2==999 || b3==999 || b4==999){
	    cout << "Inconsistentcy found!" << endl;
	    continue;
	  }

	  float WMass12 = (myJetsFilt[w1]+myJetsFilt[w2]).M();
	  float WMass13 = (myJetsFilt[w1]+myJetsFilt[w3]).M();
	  float WMass23 = (myJetsFilt[w2]+myJetsFilt[w3]).M();

	  unsigned int ind1 = 999;
	  unsigned int ind2 = 999;
	  if( TMath::Abs(WMass12-80.1) < TMath::Abs(WMass13-80.1) &&  
	      TMath::Abs(WMass12-80.1) < TMath::Abs(WMass23-80.1)){
	    ind1 = w1;
	    ind2 = w2;
	    cout << "Best W mass: " << WMass12 << endl;
	  }
	  else if( TMath::Abs(WMass13-80.1) < TMath::Abs(WMass12-80.1) &&  
		   TMath::Abs(WMass13-80.1) < TMath::Abs(WMass23-80.1) ){
	    ind1 = w1;
	    ind2 = w3;
	    cout << "Best W mass: " << WMass13 << endl;
	  }
	  else if( TMath::Abs(WMass23-80.1) < TMath::Abs(WMass12-80.1) &&  
		   TMath::Abs(WMass23-80.1) < TMath::Abs(WMass13-80.1) ){
	    ind1 = w2;
	    ind2 = w3;
	    cout << "Best W mass: " << WMass23 << endl;
	  }
	  else{
	    cout << "Inconsistency found in numJets30UntagM==3... continue" << endl;
	    continue;
	  }
	  

	  jets.clear();
	  jets.push_back( leptonLV  );  
	  jets.push_back( neutrinoLV  );  
	  jets.push_back( myJetsFilt[b1]   );
	  jets.push_back( myJetsFilt[ind1]   );  
	  jets.push_back( myJetsFilt[ind2]   );  
	  jets.push_back( myJetsFilt[b2]   );
	  jets.push_back( myJetsFilt[b3]   );
	  jets.push_back( myJetsFilt[b4]   );
	
	  meIntegrator->setJets(&jets);

	  type_ = -3;

	  meIntegrator->setIntType( MEIntegratorNew::SL2wj );	    
	  meIntegrator->setSumEt( METtype1p2corr.sumet );
	  meIntegrator->setMEtCov(-99,-99,0);
	  meIntegrator->setTopFlags( vLepton_charge[0]==1 ? +1 : -1 , vLepton_charge[0]==1 ? -1 : +1 );

	  for(int m = 0; m < nMassPoints ; m++){
	    meIntegrator->setMass( mH[m] );
	      
	    double maxP_s = 0.;
	    double maxP_b = 0.;

	    for(unsigned int pos = 0; pos < 12 ; pos++){

	      meIntegrator->initVersors( permutations_SL2wj[pos] );
	      double mass, massLow, massHigh;
	      bool skip = !(meIntegrator->compatibilityCheck(0.95, /*printP4*/ 0, mass, massLow, massHigh )) ;
	      bool skip_WHad   = !( meIntegrator->compatibilityCheck_WHad  (0.98, /*printP4*/ 0, mass, massLow, massHigh ) );
	      bool skip_TopHad = !( meIntegrator->compatibilityCheck_TopHad(0.98, /*printP4*/ 0, mass, massLow, massHigh ) );
	      if( (skip_WHad || skip_TopHad) ) continue;

	      pair<double, double> range_x0 = (meIntegrator->getW1JetEnergyCI(0.95));
	      pair<double, double> range_x1 =  make_pair(-1,1);
	      pair<double, double> range_x2 =  make_pair(-PI,PI);	    
	      pair<double, double> range_x3 =  make_pair(-1,1);
	      pair<double, double> range_x4 = useMET ? (meIntegrator->getNuPhiCI(0.95)) : make_pair(-PI,PI);
	      pair<double, double> range_x5 = (meIntegrator->getB1EnergyCI(0.95));
	      pair<double, double> range_x6 = (meIntegrator->getB2EnergyCI(0.95));
	      
	      double x0L = range_x0.first;
	      double x0U = range_x0.second;
	      double x1L = range_x1.first;
	      double x1U = range_x1.second;
	      double x2L = range_x2.first;
	      double x2U = range_x2.second;
	      double x3L = range_x3.first;
	      double x3U = range_x3.second;
	      double x4L = range_x4.first;
	      double x4U = range_x4.second;
	      double x5L = range_x5.first;
	      double x5U = range_x5.second;
	      double x6L = range_x6.first;
	      double x6U = range_x6.second;
	
	      double xLmode0[4]   = {x0L, x3L, x4L, x5L};
	      double xUmode0[4]   = {x0U, x3U, x4U, x5U};
	      double xLmode0_b[5] = {x0L, x3L, x4L, x5L, x6L};
	      double xUmode0_b[5] = {x0U, x3U, x4U, x5U, x6U};
	      
	      clock->Start();
		
	      double p = 0.;
	      
	      if(speedup==0){
		
		for(int hyp = 0 ; hyp<2; hyp++){
		  
		  if(SoB==0 && hyp!=hypo) continue;
		  if(SoB==1 && hyp==0 && skip) continue;
		  
		  cout << "Testing hypothesis " << (hyp==0 ? "S" : "B") << " for permutation " << pos << endl; 
		  meIntegrator->setHypo(hyp);
		  
		  int ntries = 0;
		  int intPoints = 2000;
		  while( ntries < MAX_REEVAL_TRIES){
		    
		    ROOT::Math::Functor toIntegrate(meIntegrator, &MEIntegratorNew::Eval,  4+hyp );
		    ROOT::Math::GSLMCIntegrator ig2( ROOT::Math::IntegrationMultiDim::kVEGAS , 1.e-12, 1.e-5, intPoints);
		    ig2.SetFunction(toIntegrate);
		    meIntegrator->SetPar( 4+hyp );	 
		    
		    p = (hyp==0 ? ig2.Integral(xLmode0, xUmode0) : ig2.Integral(xLmode0_b, xUmode0_b));
		    
		    double chi2 =  ig2.ChiSqr();
		    if( hyp==0 ){
		      if( p>maxP_s ) maxP_s = p;		  
		      else if( p<0.1*maxP_s) ntries = MAX_REEVAL_TRIES+1;
		    }
		    else{
		      if( p>maxP_b ) maxP_b = p;		  
		      else if( p<0.1*maxP_b) ntries = MAX_REEVAL_TRIES+1;	
		    }
		    
		    if( chi2 > maxChi2_ ){
		      ntries++;
		      intPoints *= 1.5;
		    }
		    else ntries = MAX_REEVAL_TRIES+1;
		    
		  }		   
				    
		  clock->Stop();
		
		  if( TMath::IsNaN(p) ) p = 0.;
		  clock->Reset();
		  
		  if( mH[m]<met+0.5 && mH[m]>met-0.5){
		    if(hyp==0)  probAtSgn_     += p;
		    else        probAtSgn_alt_ += p;
		  }
		}		  
	      }		
	    }
	  }
	  
	  
	}
	else{ 
	  continue;
	}



      }
      else if(properEventDL && (numJets30BtagM==4 && doType6) || (numJets30BtagM==3 && numJets30BtagL==4 && doType7) ){
	
	counter++;      
	if(fixNumberOfEventsPerJob && !(counter>=evLow && counter<=evHigh) ) continue;
	cout << "Processing event # " << counter << " (DL[4,X])" << endl;
	
	if( numJets30BtagM==4 && (b1==999 || b2==999 || b3==999 || b4==999) ){
	  cout << "Inconsistentcy type 6 found!" << endl;
	  cout << w1 << ", " << w2 << "," << b1 << "," << b2 << "," << b3 << "," << b4 <<  endl;
	  continue;
	}
	if( numJets30BtagM==3 && numJets30BtagL==4 && (bL1==999 || bL2==999 || bL3==999 || bL4==999) ){
	  cout << "Inconsistentcy type 6 found!" << endl;
	  cout << w1 << ", " << w2 << "," << bL1 << "," << bL2 << "," << bL3 << "," << bL4 <<  endl;
	  continue;
	}

	if( numJets30BtagM==4 ){
	  jets.clear();
	  jets.push_back( leptonLV  );  
	  jets.push_back( neutrinoLV  );  
	  jets.push_back( myJetsFilt[b1]   );
	  jets.push_back( leptonLV2 );  
	  jets.push_back( neutrinoLV       ); // dummy  
	  jets.push_back( myJetsFilt[b2]   );
	  jets.push_back( myJetsFilt[b3]   );
	  jets.push_back( myJetsFilt[b4]   );
	  type_ = 6;

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

	    //if( Bs.Pt()>20 && TMath::Abs(Bs.Eta())<2.5 ) nMatchSimBs_++;
	    //if     ( deltaR( Bs , myJetsFilt[b1] )<0.2 ) nMatchSimBs_++;
	    //else if( deltaR( Bs , myJetsFilt[b2] )<0.2 ) nMatchSimBs_++;
	    //else if( deltaR( Bs , myJetsFilt[b3] )<0.2 ) nMatchSimBs_++;
	    //else if( deltaR( Bs , myJetsFilt[b4] )<0.2 ) nMatchSimBs_++;
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

	}
	else if( numJets30BtagM==3 && numJets30BtagL==4){
	  jets.clear();
	  jets.push_back( leptonLV  );  
	  jets.push_back( neutrinoLV  );  
	  jets.push_back( myJetsFilt[bL1]   );
	  jets.push_back( leptonLV2 );  
	  jets.push_back( neutrinoLV       ); // dummy  
	  jets.push_back( myJetsFilt[bL2]   );
	  jets.push_back( myJetsFilt[bL3]   );
	  jets.push_back( myJetsFilt[bL4]   );
	  type_ = 7;

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

	    //if( Bs.Pt()>20 && TMath::Abs(Bs.Eta())<2.5 ) nMatchSimBs_++;
	    //if     ( deltaR( Bs , myJetsFilt[bL1] )<0.2 ) nMatchSimBs_++;
	    //else if( deltaR( Bs , myJetsFilt[bL2] )<0.2 ) nMatchSimBs_++;
	    //else if( deltaR( Bs , myJetsFilt[bL3] )<0.2 ) nMatchSimBs_++;
	    //else if( deltaR( Bs , myJetsFilt[bL4] )<0.2 ) nMatchSimBs_++;
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

	}
	else if( numJets30BtagM==2 && numJets30BtagL==4){
	  jets.clear();
	  jets.push_back( leptonLV  );  
	  jets.push_back( neutrinoLV  );  
	  jets.push_back( myJetsFilt[bL1]   );
	  jets.push_back( leptonLV2 );  
	  jets.push_back( neutrinoLV       ); // dummy  
	  jets.push_back( myJetsFilt[bL2]   );
	  jets.push_back( myJetsFilt[bL3]   );
	  jets.push_back( myJetsFilt[bL4]   );
	  type_ = 8;
	}
	else if( numJets30BtagM==1 && numJets30BtagL==4 ){
	  jets.clear();
	  jets.push_back( leptonLV  );  
	  jets.push_back( neutrinoLV  );  
	  jets.push_back( myJetsFilt[bL1]   );
	  jets.push_back( leptonLV2 );  
	  jets.push_back( neutrinoLV       ); // dummy  
	  jets.push_back( myJetsFilt[bL2]   );
	  jets.push_back( myJetsFilt[bL3]   );
	  jets.push_back( myJetsFilt[bL4]   );
	  type_ = 9;
	}
	else if(numJets30BtagM==0 && numJets30BtagL==4){
	  jets.clear();
	  jets.push_back( leptonLV  );  
	  jets.push_back( neutrinoLV  );  
	  jets.push_back( myJetsFilt[bL1]   );
	  jets.push_back( leptonLV2 );  
	  jets.push_back( neutrinoLV       ); // dummy  
	  jets.push_back( myJetsFilt[bL2]   );
	  jets.push_back( myJetsFilt[bL3]   );
	  jets.push_back( myJetsFilt[bL4]   );
	  type_ = 10;
	}
	else{
	  cout << "Inconsistentcy if elseif type 6 found!" << endl;
	  continue;
	}

	 
	meIntegrator->setJets(&jets);
	meIntegrator->setIntType( MEIntegratorNew::DL );	
	meIntegrator->setSumEt( METtype1p2corr.sumet );
	meIntegrator->setMEtCov(-99,-99,0);
	meIntegrator->setTopFlags( vLepton_charge[0]==1 ? +1 : -1 , vLepton_charge[0]==1 ? -1 : +1 );
		
	clock->Start();

	for(int m = 0; m < nMassPoints ; m++){

	  meIntegrator->setMass( mH[m] );
	    
	  double maxP_s = 0.;
	  double maxP_b = 0.;
	  
	  for(unsigned int pos = 0; pos < 12 ; pos++){
	    
	    meIntegrator->initVersors( permutations_DL[pos] );
	    double mass, massLow, massHigh;
	    bool skip = !(meIntegrator->compatibilityCheck(0.95, /*printP4*/ 0, mass, massLow, massHigh )) ;
	    
	    pair<double, double> range_x0 = (meIntegrator->getW1JetEnergyCI(0.95));
	    pair<double, double> range_x1 =  make_pair(-1,1);
	    pair<double, double> range_x2 =  make_pair(-PI,PI);	    
	    pair<double, double> range_x3 =  make_pair(-1,1);
	    pair<double, double> range_x4 = useMET ? (meIntegrator->getNuPhiCI(0.95)) : make_pair(-PI,PI);
	    pair<double, double> range_x5 = (meIntegrator->getB1EnergyCI(0.95));
	    pair<double, double> range_x6 = (meIntegrator->getB2EnergyCI(0.95));

	    double x0L = range_x0.first;
	    double x0U = range_x0.second;
	    double x1L = range_x1.first;
	    double x1U = range_x1.second;
	    double x2L = range_x2.first;
	    double x2U = range_x2.second;
	    double x3L = range_x3.first;
	    double x3U = range_x3.second;
	    double x4L = range_x4.first;
	    double x4U = range_x4.second;
	    double x5L = range_x5.first;
	    double x5U = range_x5.second;
	    double x6L = range_x6.first;
	    double x6U = range_x6.second;
	
	    double xLmode6[5]   = {x1L, x2L, x1L, x2L, x5L};
	    double xUmode6[5]   = {x1U, x2U, x1U, x2U, x5U};
	    double xLmode6_b[6] = {x1L, x2L, x1L, x2L, x5L, x6L};
	    double xUmode6_b[6] = {x1U, x2U, x1U, x2U, x5U, x6U};
	      
	    double p = 0.;
	    if(speedup==0){
	      
	      for(int hyp = 0 ; hyp<2; hyp++){
		
		if(SoB==0 && hyp!=hypo) continue;
		if(SoB==1 && hyp==0 && skip) continue;

		cout << "Testing hypothesis " << (hyp==0 ? "S" : "B") << " for permutation " << pos << endl; 
		meIntegrator->setHypo(hyp);

		int ntries = 0;
		int intPoints = 10000;		    
		while( ntries < MAX_REEVAL_TRIES){
		  
		  ROOT::Math::Functor toIntegrate(meIntegrator, &MEIntegratorNew::Eval,  5+hyp);
		  ROOT::Math::GSLMCIntegrator ig2( ROOT::Math::IntegrationMultiDim::kVEGAS , 1.e-12, 1.e-5, intPoints);
		  ig2.SetFunction(toIntegrate);
		  meIntegrator->SetPar(  5+hyp );

		  p = (hyp==0 ? ig2.Integral(xLmode6, xUmode6) : ig2.Integral(xLmode6_b, xUmode6_b));

		  double chi2 =  ig2.ChiSqr();
		  if( hyp==0 ){
		    if( p>maxP_s ) maxP_s = p;		  
		    else if( p<0.1*maxP_s) ntries = MAX_REEVAL_TRIES+1;
		  }
		  else{
		    if( p>maxP_b ) maxP_b = p;		  
		    else if( p<0.1*maxP_b) ntries = MAX_REEVAL_TRIES+1;	
		  }

		  if( chi2 > maxChi2_ ){
		    ntries++;
		    intPoints *= 1.5;
		  }
		  else ntries = MAX_REEVAL_TRIES+1;

		}
		
		if( TMath::IsNaN(p) ) p = 0.;

  
		if( p>= maxP_s && hyp==0 ){
		  int lastTwoDigits = (permutations_DL[pos])%100;
		  int higgs1 = lastTwoDigits/(10);
		  int higgs2 = lastTwoDigits%10;		    
		  //cout << "pos=" << permutations_SL2wj[pos] << " => " << higgs1 << "," << higgs2 << endl;		    
		  if( higgs1 == 2)
		    higgs1 = b1;		      
		  else if(  higgs1 == 5 )
		    higgs1 = b2;
		  else if(  higgs1 == 6 )
		    higgs1 = b3;
		  else if(  higgs1 == 7 )
		    higgs1 = b4;
		  else{ cout << "What happened?" << endl; }
		  if( higgs2 == 2)
		    higgs2 = b1;		      
		  else if(  higgs2 == 5 )
		    higgs2 = b2;
		  else if(  higgs2 == 6 )
		    higgs2 = b3;
		  else if(  higgs2 == 7 )
		    higgs2 = b4;
		  else{ cout << "What happened?" << endl; }		    
		  //cout << "Corrspond to " << higgs1 << "," << higgs2 << endl;		    
		  int index1 =  mapFilt[higgs1].index;		      
		  int index2 =  mapFilt[higgs2].index;		      		    
		  float csv_1 = index1>=0 ? hJet_csv_nominal[index1]  : aJet_csv_nominal[-index1-1];
		  float csv_2 = index2>=0 ? hJet_csv_nominal[index2]  : aJet_csv_nominal[-index2-1];		    
		  //cout << "Deug 1: " << mapFilt[higgs1].pt << " --- " << (meIntegrator->jetAt(6)).Pt() << endl;
		  //cout << "Deug 2: " << mapFilt[higgs2].pt << " --- " << (meIntegrator->jetAt(7)).Pt() << endl;		    
		  commitVariables( meIntegrator->jetAt(6), meIntegrator->jetAt(7), csv_1, csv_2,
				   dR_,dPhi_,dEta_,Mjj_,Pt_,Eta_,pt1_,pt2_,eta1_,eta2_,csv1_,csv2_);
		}

		if( mH[m]<met+0.5 && mH[m]>met-0.5){
		  if(hyp==0)  probAtSgn_     += p;
		  else        probAtSgn_alt_ += p;
		}
		
	      }
	    }	
	  }
	}

	clock->Stop();
	time_ = clock->RealTime();
	clock->Reset();	 

      }
      else{
	continue;
      }
      /////////////////////////////////////////////////////
      /*
      // overload with smeared vectors
      TOPLEPW1  =  meIntegrator->jetAt(0);
      TOPLEPW2  =  meIntegrator->jetAt(1);
      TOPLEPB   =  meIntegrator->jetAt(2);
      TOPHADW1  =  meIntegrator->jetAt(3);
      TOPHADW2  =  meIntegrator->jetAt(4);
      TOPHADB   =  meIntegrator->jetAt(5);
      genBLV    =  meIntegrator->jetAt(6);
      genBbarLV =  meIntegrator->jetAt(7);


      int HADW1pass = 0;
      if(abs(flavHADW1)>0 && abs(flavHADW1)<3)
	HADW1pass = ( ran->Uniform(0,1)<=udEff ? 1 : 0);
      else if (abs(flavHADW1)==3)
	HADW1pass = ( ran->Uniform(0,1)<=sEff ? 1 : 0);
      else if (abs(flavHADW1)==4)
	HADW1pass = ( ran->Uniform(0,1)<=cEff ? 1 : 0);
      else
	HADW1pass = 0;

      int HADW2pass = 0;
      if(abs(flavHADW2)>0 && abs(flavHADW2)<3)
	HADW2pass = ( ran->Uniform(0,1)<=udEff ? 1 : 0);
      else if (abs(flavHADW2)==3)
	HADW2pass = ( ran->Uniform(0,1)<=sEff ? 1 : 0);
      else if (abs(flavHADW2)==4)
	HADW2pass = ( ran->Uniform(0,1)<=cEff ? 1 : 0);
      else
	HADW2pass = 0;

      int LEPBpass    = ( ran->Uniform(0,1)<=bEff ? 1 : 0);
      int HADBpass    = ( ran->Uniform(0,1)<=bEff ? 1 : 0);
      int genBpass    = ( ran->Uniform(0,1)<=bEff ? 1 : 0);
      int genBbarpass = ( ran->Uniform(0,1)<=bEff ? 1 : 0);

      int LEPW1in   =  TOPLEPW1.Pt() >20 && TMath::Abs(TOPLEPW1.Eta()) <2.1;
      int LEPW2in   =  TOPHADW1.Pt() >20 && TMath::Abs(TOPHADW1.Eta()) <2.1;
      int HADW1in   =  TOPHADW1.Pt() >30 && TMath::Abs(TOPHADW1.Eta()) <2.5;
      int HADW2in   =  TOPHADW2.Pt() >30 && TMath::Abs(TOPHADW2.Eta()) <2.5;     
      int LEPBin    =  TOPLEPB.Pt()  >30 && TMath::Abs(TOPLEPB.Eta())  <2.5;
      int HADBin    =  TOPHADB.Pt()  >30 && TMath::Abs(TOPHADB.Eta())  <2.5;
      int genBin    =  genBLV.Pt()   >30 && TMath::Abs(genBLV.Eta())   <2.5;
      int genBbarin =  genBbarLV.Pt()>30 && TMath::Abs(genBbarLV.Eta())<2.5;

      if(printP4){
	cout << "LEPBpass    " << LEPBpass << endl;
	cout << "HADW1pass   " << HADW1pass <<  endl;
	cout << "HADW2pass   " << HADW2pass << endl;
	cout << "HADBpass    " << HADBpass <<  endl;
	cout << "genBpass    " << genBpass << endl;
	cout << "genBbarpass " << genBbarpass <<  endl;
	cout << "LEPBin    " << LEPBin << endl;
	cout << "HADW1in   " << HADW1in <<  endl;
	cout << "HADW2in   " << HADW2in << endl;
	cout << "HADBin    " << HADBin <<  endl;
	cout << "genBin    " << genBin << endl;
	cout << "genBbarin " << genBbarin <<  endl;
	cout << "LEPW1in   " << LEPW1in <<  endl;
	cout << "LEPW2in   " << LEPW2in << endl;
	cout << "******* END INPUT ******" << endl;   
      }
      */
      /////////////////////////////////////////////////////
	counter_   = counter;
	nSimBs_    = nSimBs;
	nC_        = nC;
	nCTop_     = nCTop;
	weight_    = scaleFactor;

	tree->Fill();      
    } // nentries

    hcounter->SetBinContent(1,float(events_)/nentries);    

  } // samples

  fout_tmp->cd();

  hcounter->Write("",TObject::kOverwrite );
  tree->Write("",TObject::kOverwrite );
  fout_tmp->Close();

  cout << "Delete meIntegrator..." << endl;
  delete meIntegrator;
  delete ran; delete clock;
  cout << "Finished!!!" << endl;
  
  return 0;
}





/*
/////////////////////////////////////////////////////////////
float sv_bLep = -99;
float sv_bHad = -99;
float sv_b1   = -99;
float sv_b2   = -99;
for(int l = 0; l < nSvs; l++){
TLorentzVector vertex(1,0,0,1);
vertex.SetPtEtaPhiM( Svpt[l], Sveta[l], Svphi[l], SvmassSv[l]);
if     (sv_bLep<0 && deltaR(vertex, meIntegrator->jetAt(2) )<0.5) sv_bLep = SvmassSv[l] ;
else if(sv_bHad<0 && deltaR(vertex, meIntegrator->jetAt(5) )<0.5) sv_bHad = SvmassSv[l];
else if(sv_b1<0   && deltaR(vertex, meIntegrator->jetAt(6) )<0.5) sv_b1   = SvmassSv[l];
else if(sv_b2<0   && deltaR(vertex, meIntegrator->jetAt(7) )<0.5) sv_b2   = SvmassSv[l];
else{}
}
if( sv_bLep<0 || (sv_bLep<1.8) ){
float m1    = (meIntegrator->jetAt(2)+meIntegrator->jetAt(3)).M();
float m2    = (meIntegrator->jetAt(2)+meIntegrator->jetAt(4)).M();
if( ((m1>(MwL+5) && m1<(MwH-5)) || (m2>(MwL+5) && m2<(MwH-5))) && type_== 0 ) flag_type0_ = 1; 
if( ((m1>(MwL+5) && m1<(MwH-5)) || (m2>(MwL+5) && m2<(MwH-5))) && type_== 1 ) flag_type1_ = 1; 
}
if( sv_bHad<0 || (sv_bHad<1.8) ){
float m1 = (meIntegrator->jetAt(5)+meIntegrator->jetAt(3)).M();
float m2 = (meIntegrator->jetAt(5)+meIntegrator->jetAt(4)).M();
if( ((m1>(MwL+5) && m1<(MwH-5)) || (m2>(MwL+5) && m2<(MwH-5))) && type_== 0 ) flag_type0_ = 1; 
if( ((m1>(MwL+5) && m1<(MwH-5)) || (m2>(MwL+5) && m2<(MwH-5))) && type_== 1 ) flag_type1_ = 1; 
}
if( sv_b1<0 || (sv_b1<1.8) ){
float m1 = (meIntegrator->jetAt(6)+meIntegrator->jetAt(3)).M();
float m2 = (meIntegrator->jetAt(6)+meIntegrator->jetAt(4)).M();
if( ((m1>(MwL+5) && m1<(MwH-5)) || (m2>(MwL+5) && m2<(MwH-5))) && type_== 0 ) flag_type0_ = 1; 
if( ((m1>(MwL+5) && m1<(MwH-5)) || (m2>(MwL+5) && m2<(MwH-5))) && type_== 1 ) flag_type1_ = 1; 
}
if( sv_b2<0 || (sv_b2<1.8) ){
float m1 = (meIntegrator->jetAt(7)+meIntegrator->jetAt(3)).M();
float m2 = (meIntegrator->jetAt(7)+meIntegrator->jetAt(4)).M();
if( ((m1>(MwL+5) && m1<(MwH-5)) || (m2>(MwL+5) && m2<(MwH-5))) && type_== 0 ) flag_type0_ = 1; 
if( ((m1>(MwL+5) && m1<(MwH-5)) || (m2>(MwL+5) && m2<(MwH-5))) && type_== 1 ) flag_type1_ = 1; 
}
/////////////////////////////////////////////////////////////
*/
