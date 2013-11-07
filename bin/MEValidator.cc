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

#include "Math/IntegratorOptions.h"

#define GENJETDR         0.3
#define N_TRIES_MAX       10
#define MAX_REEVAL_TRIES   3
#define VERBOSE        false
#define PI       TMath::Pi()

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

pair<double,double> getMaxValue( TH1F* hMassProb){

  float est  = -99;
  float prob = -99;
  float a = -99;
  float b = -99;
  float c = -99;
  
  int maxBin = hMassProb->GetMaximumBin();
  
  if( maxBin<  hMassProb->GetNbinsX() && maxBin>1 ){
    
    double xD =  hMassProb->GetBinCenter (hMassProb->GetMaximumBin()-1);
    double xC =  hMassProb->GetBinCenter (hMassProb->GetMaximumBin());
    double xU =  hMassProb->GetBinCenter (hMassProb->GetMaximumBin()+1);
    double yD =  hMassProb->GetBinContent(hMassProb->GetMaximumBin()-1);
    double yC =  hMassProb->GetBinContent(hMassProb->GetMaximumBin());
    double yU =  hMassProb->GetBinContent(hMassProb->GetMaximumBin()+1);
    
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
      { 1., 1., 1.
      };
    C.SetMatrixArray(dummy,"C");
    C.Mult(AInv,Y);
    
    a = C(2,0);
    b = C(1,0);
    c = C(0,0);

    est  = -b/2/a ;
    prob = a*est*est + b*est + c;
  }
  else if(maxBin== hMassProb->GetNbinsX()){
    
    double xD =  hMassProb->GetBinCenter (hMassProb->GetMaximumBin()-2);
    double xC =  hMassProb->GetBinCenter (hMassProb->GetMaximumBin()-1);
    double xU =  hMassProb->GetBinCenter (hMassProb->GetMaximumBin());
    double yD =  hMassProb->GetBinContent(hMassProb->GetMaximumBin()-2);
    double yC =  hMassProb->GetBinContent(hMassProb->GetMaximumBin()-1);
    double yU =  hMassProb->GetBinContent(hMassProb->GetMaximumBin());
    
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
      { 1., 1., 1.
      };
    C.SetMatrixArray(dummy,"C");
    C.Mult(AInv,Y);
    
    a = C(2,0);
    b = C(1,0);
    c = C(0,0);

    est = -b/2/a ;
    prob = a*est*est + b*est + c;	  
  }
  else if(maxBin==1){
    
    double xD =  hMassProb->GetBinCenter (hMassProb->GetMaximumBin());
    double xC =  hMassProb->GetBinCenter (hMassProb->GetMaximumBin()+1);
    double xU =  hMassProb->GetBinCenter (hMassProb->GetMaximumBin()+2);
    double yD =  hMassProb->GetBinContent(hMassProb->GetMaximumBin());
    double yC =  hMassProb->GetBinContent(hMassProb->GetMaximumBin()+1);
    double yU =  hMassProb->GetBinContent(hMassProb->GetMaximumBin()+2);
    
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
      { 1., 1., 1.
      };
    C.SetMatrixArray(dummy,"C");
    C.Mult(AInv,Y);
    
    a = C(2,0);
    b = C(1,0);
    c = C(0,0);

    est = -b/2/a ;
    prob = a*est*est + b*est + c;
  }
  else{
    est  =  hMassProb->GetBinCenter  (hMassProb->GetMaximumBin());
    prob =  hMassProb->GetBinContent (hMassProb->GetMaximumBin());
  }

  return make_pair(est, prob);

}


int main(int argc, const char* argv[])
{

  std::cout << "MEValidator" << std::endl;
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
  int vegasPoints          ( in.getParameter<int>   ("vegasPoints"   ) );
  bool verbose             ( in.getParameter<bool>  ("verbose") );
  double lumi              ( in.getParameter<double>("lumi") );

  float met       ( in.getParameter<double> ("met")        );

  float scaleH  ( in.getParameter<double> ("scaleH")   );
  float scaleL  ( in.getParameter<double> ("scaleL")   );
  float scaleMET( in.getParameter<double> ("scaleMET") );

  int   switchoffOL( in.getUntrackedParameter<int>    ("switchoffOL", 0));
  int   useME      ( in.getParameter<int>    ("useME")      );
  int   useJac     ( in.getParameter<int>    ("useJac")     );
  int   useMET     ( in.getParameter<int>    ("useMET")     );
  int   useTF      ( in.getParameter<int>    ("useTF")      );
  int   usePDF     ( in.getParameter<int>    ("usePDF")     );

  int   doParton      ( in.getParameter<int>    ("doParton")     );
  int   doSmear       ( in.getParameter<int>    ("doSmear")     );
  int   doMassScan    ( in.getParameter<int>    ("doMassScan")     );
  int   doPermutations( in.getParameter<int>    ("doPermutations") );

  int   printP4         ( in.getParameter<int>    ("printP4")          );
  int   mode            ( in.getUntrackedParameter<int>    ("mode", 0) );
  int   norm            ( in.getUntrackedParameter<int>    ("norm", 0) );

  int   newVars         ( in.getUntrackedParameter<int>    ("newVars", 0) );

  double maxChi2_       ( in.getUntrackedParameter<double>  ("maxChi2", 10.));

  int   hypo                 ( in.getUntrackedParameter<int>    ("hypo", 0)     );
  int   SoB                  ( in.getUntrackedParameter<int>    ("SoB", 0)     );
  int   doBkgWithWorkaraound ( in.getUntrackedParameter<int>    ("doBkgWithWorkaraound", 0)     );
  vector<double> masses       ( in.getParameter<vector<double> >  ("masses")        );
  vector<string> functions    ( in.getParameter<vector<string> >  ("functions")     );

  vector<int> evLimits( in.getParameter<vector<int> >  ("evLimits") );
  int evLow  = evLimits[0];
  int evHigh = evLimits[1];

  gSystem->Exec(("rm "+outFileName).c_str());

  TStopwatch* clock = new TStopwatch();
  TRandom3*   ran   = new TRandom3();
  ran->SetSeed(4321);


  const int nMassPoints  = masses.size();
  double mH[nMassPoints];
  for( unsigned int m = 0; m < masses.size() ; m++)
    mH[m] = masses[m];


  int nPermutations;
  switch(mode){
  case 0:
    nPermutations = 12;          
    break;
  case 1:
    nPermutations = 12;
    break;
  case 2:
    nPermutations = 3;
    break;
  case 3:
    nPermutations = 3;
    break;
  case 4:
    nPermutations = 6;
    break;
  case 5:
    nPermutations = 12;
    break;
  case 6:
    nPermutations = 12;
    break;
  default:
    nPermutations = 12;
    break;
  }

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
    {234567,     // CORRECT 
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
  

  int permutations[nPermutations];
  for(int per = 0; per < nPermutations; per++){
    if(mode==0) 
      permutations[per] = permutations_SL2wj[per];
    else if(mode==1) 
      permutations[per] = permutations_SL1wj[per];
    else if(mode==2) 
      permutations[per] = permutations_SLNoBHad[per];
    else if(mode==3) 
      permutations[per] = permutations_SLNoBLep[per];
    else if(mode==4) 
      permutations[per] = permutations_SLNoHiggs[per];
    else if(mode==5) 
      permutations[per] = permutations_SL3b[per];
    else if(mode==6) 
      permutations[per] = permutations_DL[per];
  }


  //{ /*50,  55,  */60,  65,  70,  75,  
  // 80 , 85,  90, 95, 100, 105, 110, 
  // 115, 120, 125, 130, 135, 140, 145, 150, 155, 
  // 160, 165, 170, 175, 180 ,185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250 
  //}; 


  ///////////////////////////////////////////////////////

  TH1F*  hMass_gen         = new TH1F("m_good_gen","",        nMassPoints, mH[0]-2.5, mH[nMassPoints-1]+2.5);
  TH1F*  hMass_alt_gen     = new TH1F("m_good_alt_gen","",    nMassPoints, mH[0]-2.5, mH[nMassPoints-1]+2.5);
  TH1F*  hMass_rec         = new TH1F("m_good_rec","",        nMassPoints, mH[0]-2.5, mH[nMassPoints-1]+2.5);
  TH1F*  hMass_alt_rec     = new TH1F("m_good_alt_rec","",    nMassPoints, mH[0]-2.5, mH[nMassPoints-1]+2.5);
  TH1F*  hMassProb_gen     = new TH1F("m_all_gen","",         nMassPoints, mH[0]-2.5, mH[nMassPoints-1]+2.5);
  TH1F*  hMassProb_alt_gen = new TH1F("m_all_alt_gen","",     nMassPoints, mH[0]-2.5, mH[nMassPoints-1]+2.5);
  TH1F*  hMassProb_rec     = new TH1F("m_all_rec","",         nMassPoints, mH[0]-2.5, mH[nMassPoints-1]+2.5);
  TH1F*  hMassProb_alt_rec = new TH1F("m_all_alt_rec","",     nMassPoints, mH[0]-2.5, mH[nMassPoints-1]+2.5);

  //THStack*  sMassProbInt_gen  = new THStack("stack_m_all_gen", "Per permutation mass, GEN");
  //THStack*  sMassProbInt_rec  = new THStack("stack_m_all_rec", "Per permutation mass, REC");
  TH1F*  hInt     = new TH1F("hInt","",nMassPoints, mH[0]-2.5, mH[nMassPoints-1]+2.5);  
 

  TH1F*  hBestMass_gen        = new TH1F("m_best_good_gen","",   nMassPoints, mH[0]-2.5, mH[nMassPoints-1]+2.5);
  TH1F*  hBestMass_rec        = new TH1F("m_best_good_rec","",   nMassPoints, mH[0]-2.5, mH[nMassPoints-1]+2.5);
  TH1F*  hBestMassProb_gen    = new TH1F("m_best_all_gen","",    nMassPoints, mH[0]-2.5, mH[nMassPoints-1]+2.5);
  TH1F*  hBestMassProbInt_gen = new TH1F("m_bestInt_all_gen","", nMassPoints, mH[0]-2.5, mH[nMassPoints-1]+2.5);
  TH1F*  hBestMassProbMax_gen = new TH1F("m_bestMax_all_gen","", nMassPoints, mH[0]-2.5, mH[nMassPoints-1]+2.5);
  TH1F*  hBestMassProb_rec    = new TH1F("m_best_all_rec","",    nMassPoints, mH[0]-2.5, mH[nMassPoints-1]+2.5);
  TH1F*  hBestMassProbInt_rec = new TH1F("m_bestInt_all_rec","", nMassPoints, mH[0]-2.5, mH[nMassPoints-1]+2.5);
  TH1F*  hBestMassProbMax_rec = new TH1F("m_bestMax_all_rec","", nMassPoints, mH[0]-2.5, mH[nMassPoints-1]+2.5);

  TTree* tree  = new TTree("tree","");

  int counter_;
  int type_;

  float probGG_;
  float probAtSgnGG_,probAtSgnGG_alt_;
  float massGG_;
  float evalCpuGG_,evalReaGG_,evalCpuGG_alt_,evalReaGG_alt_;
  float chi2GG_;

  float probRG_;
  float probAtSgnRG_, probAtSgnRG_alt_;
  float massRG_;
  float evalCpuRG_,evalReaRG_,evalCpuRG_alt_,evalReaRG_alt_;
  float chi2RG_;

  float probGA_;
  float probAtSgnGA_,probAtSgnGA_alt_;
  float probAtGoodGA_;
  float massGA_, massIntGA_, massMaxGA_, massMaxModGA_, massIntModGA_;
  float probMaxGA_,probMaxModGA_,probIntGA_ , probIntModGA_ ;
  float evalCpuGA_,evalReaGA_,evalCpuGA_alt_,evalReaGA_alt_;

  int posIntGA_,posIntRA_,posMaxGA_,posMaxRA_;

  float probRA_;
  float probAtSgnRA_, probAtSgnRA_alt_;
  float probAtGoodRA_;
  float massRA_,massIntRA_, massMaxRA_,massMaxModRA_, massIntModRA_;
  float probMaxRA_,probMaxModRA_,probIntRA_, probIntModRA_ ;
  float evalCpuRA_,evalReaRA_, evalCpuRA_alt_,evalReaRA_alt_;
  float dPx_, dRecPx_;

  int matchByIntGA_, matchByIntModGA_, matchByMaxGA_, matchByMaxModGA_, matchByIntRA_, matchByIntModRA_, matchByMaxRA_, matchByMaxModRA_;

  tree->Branch("counter",                       &counter_,      "counter/I");
  if(newVars) tree->Branch("type",                          &type_,         "type/I");
  tree->Branch("p_best_good_gen",               &probGG_,       "p_best_good_gen/F");
   if(newVars)tree->Branch("chi2_good_gen",                 &chi2GG_,       "chi2_good_gen/F");
  tree->Branch(Form("p_%d_good_gen",int(met)),  &probAtSgnGG_, Form("p_%d_good_gen/F",int(met)) );
  if(newVars) tree->Branch(Form("p_%d_good_alt_gen",int(met)),  &probAtSgnGG_alt_, Form("p_%d_good_alt_gen/F",int(met)) );
  tree->Branch("m_best_good_gen",               &massGG_,       "m_best_good_gen/F");
  tree->Branch("cputime_good_gen",              &evalCpuGG_,    "cputime_good_gen/F");
  tree->Branch("realtime_good_gen",             &evalReaGG_,    "realtime_good_gen/F");
  if(newVars) tree->Branch("cputime_good_alt_gen",              &evalCpuGG_alt_,    "cputime_good_alt_gen/F");
  if(newVars) tree->Branch("realtime_good_alt_gen",             &evalReaGG_alt_,    "realtime_good_alt_gen/F");
  tree->Branch("p_best_good_rec",               &probRG_,       "p_best_good_rec/F");
  if(newVars)tree->Branch("chi2_good_rec",                 &chi2RG_,       "chi2_good_rec/F");
  tree->Branch(Form("p_%d_good_rec",int(met)),  &probAtSgnRG_, Form("p_%d_good_rec/F",int(met)) );
  if(newVars) tree->Branch(Form("p_%d_good_alt_rec",int(met)),  &probAtSgnRG_alt_, Form("p_%d_good_alt_rec/F",int(met)) );
  tree->Branch("m_best_good_rec",               &massRG_,       "m_best_good_rec/F");
  tree->Branch("cputime_good_rec",              &evalCpuRG_,    "cputime_good_rec/F");
  tree->Branch("realtime_good_rec",             &evalReaRG_,    "realtime_good_rec/F");
  if(newVars) tree->Branch("cputime_good_alt_rec",              &evalCpuRG_alt_,    "cputime_good_alt_rec/F");
  if(newVars) tree->Branch("realtime_good_alt_rec",             &evalReaRG_alt_,    "realtime_good_alt_rec/F");
  tree->Branch("p_best_all_gen",                &probGA_,       "p_best_all_gen/F");
  tree->Branch(Form("p_%d_all_gen",int(met)),   &probAtSgnGA_, Form("p_%d_all_gen/F",int(met)) );
  if(newVars) tree->Branch(Form("p_%d_all_alt_gen",int(met)),   &probAtSgnGA_alt_, Form("p_%d_all_alt_gen/F",int(met)) );
  tree->Branch("p_good_all_gen",                &probAtGoodGA_, "p_good_all_gen/F");
  tree->Branch("m_best_all_gen",                &massGA_,       "m_best_all_gen/F");
  tree->Branch("m_bestInt_all_gen",             &massIntGA_,    "m_bestInt_all_gen/F");
  if(newVars) tree->Branch("pos_bestInt_all_gen",           &posIntGA_,     "pos_bestInt_all_gen/I");
  if(newVars) tree->Branch("m_bestIntMod_all_gen",          &massIntModGA_, "m_bestIntMod_all_gen/F");
  tree->Branch("m_bestMax_all_gen",             &massMaxGA_,    "m_bestMax_all_gen/F");
  if(newVars) tree->Branch("pos_bestMax_all_gen",           &posMaxGA_,     "pos_bestMax_all_gen/I");
  if(newVars) tree->Branch("m_bestMaxMod_all_gen",          &massMaxModGA_, "m_bestMaxMod_all_gen/F");
  tree->Branch("p_bestInt_all_gen",             &probIntGA_,    "p_bestInt_all_gen/F");
  if(newVars) tree->Branch("p_bestIntMod_all_gen",          &probIntModGA_, "p_bestIntMod_all_gen/F");
  tree->Branch("p_bestMax_all_gen",             &probMaxGA_,    "p_bestMax_all_gen/F");
  if(newVars) tree->Branch("p_bestMaxMod_all_gen",          &probMaxModGA_, "p_bestMaxMod_all_gen/F");
  tree->Branch("cputime_all_gen",               &evalCpuGA_,    "cputime_all_gen/F");
  tree->Branch("realtime_all_gen",              &evalReaGA_,    "realtime_all_gen/F");
  if(newVars) tree->Branch("cputime_all_alt_gen",               &evalCpuGA_alt_,    "cputime_all_alt_gen/F");
  if(newVars) tree->Branch("realtime_all_alt_gen",              &evalReaGA_alt_,    "realtime_all_alt_gen/F");
  tree->Branch("p_best_all_rec",                &probRA_,       "p_best_all_rec/F");
  tree->Branch(Form("p_%d_all_rec",int(met)),   &probAtSgnRA_, Form("p_%d_all_rec/F",int(met)) );
  if(newVars) tree->Branch(Form("p_%d_all_alt_rec",int(met)),   &probAtSgnRA_alt_, Form("p_%d_all_alt_rec/F",int(met)) );
  tree->Branch("p_good_all_rec",                &probAtGoodRA_, "p_good_all_rec/F");
  tree->Branch("m_best_all_rec",                &massRA_,       "m_best_all_rec/F");
  tree->Branch("m_bestInt_all_rec",             &massIntRA_,    "m_bestInt_all_rec/F");
  if(newVars) tree->Branch("pos_bestInt_all_rec",           &posIntRA_,     "pos_bestInt_all_rec/I");
  if(newVars) tree->Branch("m_bestIntMod_all_rec",          &massIntModRA_, "m_bestIntMod_all_rec/F");
  tree->Branch("m_bestMax_all_rec",             &massMaxRA_,    "m_bestMax_all_rec/F");
  if(newVars) tree->Branch("m_bestMaxMod_all_rec",          &massMaxModRA_, "m_bestMaxMod_all_rec/F");
  tree->Branch("p_bestInt_all_rec",             &probIntRA_,    "p_bestInt_all_rec/F");
  if(newVars) tree->Branch("p_bestIntMod_all_rec",          &probIntModRA_, "p_bestIntMod_all_rec/F");
  tree->Branch("p_bestMax_all_rec",             &probMaxRA_,    "p_bestMax_all_rec/F");
  if(newVars) tree->Branch("pos_bestMax_all_rec",           &posMaxRA_,     "pos_bestMax_all_rec/I");  
  if(newVars) tree->Branch("p_bestMaxMod_all_rec",          &probMaxModRA_, "p_bestMaxMod_all_rec/F");
  tree->Branch("cputime_all_rec",               &evalCpuRA_,    "cputime_all_rec/F");
  tree->Branch("realtime_all_rec",              &evalReaRA_,    "realtime_all_rec/F");
  if(newVars) tree->Branch("cputime_all_alt_rec",               &evalCpuRA_alt_,    "cputime_all_alt_rec/F");
  if(newVars) tree->Branch("realtime_all_alt_rec",              &evalReaRA_alt_,    "realtime_all_alt_rec/F");
  tree->Branch("match_bestInt_all_gen",         &matchByIntGA_, "match_bestInt_all_gen/I");
  if(newVars) tree->Branch("match_bestIntMod_all_gen",      &matchByIntModGA_,"match_bestIntMod_all_gen/I");
  tree->Branch("match_bestMax_all_gen",         &matchByMaxGA_,   "match_bestMax_all_gen/I");
  if(newVars) tree->Branch("match_bestMaxMod_all_gen",      &matchByMaxModGA_,"match_bestMaxMod_all_gen/I");
  tree->Branch("match_bestInt_all_rec",         &matchByIntRA_,   "match_bestInt_all_rec/I");
  if(newVars) tree->Branch("match_bestIntMod_all_rec",      &matchByIntModRA_,"match_bestIntMod_all_rec/I");
  tree->Branch("match_bestMax_all_rec",         &matchByMaxRA_,   "match_bestMax_all_rec/I");
  if(newVars) tree->Branch("match_bestMaxMod_all_rec",      &matchByMaxModRA_,"match_bestMaxMod_all_rec/I");

  if(newVars) tree->Branch("dPx",                       &dPx_,      "dPx/F");
  if(newVars) tree->Branch("dRecPx",                    &dRecPx_,   "dRecPx/F");
  ///////////////////////////////////////////////////////



  int par;
  switch(mode){
  case 0:
    par = 4;  // SL all partons
    break;
  case 1:
    par = 6;  // SL, 1w lost
    break;
  case 2:
    par = 6;  // SL, bhad lost
    break;
  case 3:
    par = 6;  // SL, bLep lost
    break;
  case 4:
    par = 6;  // SL, higgs2 lost
    break;
  case 5:
    par = 6;  // SL, 3b
    break;
  case 6:
    par = 5;  // DL
    break;
  case 7:
    par = 19; // inclusive xsec
    break;
  default:
    par = 19;
    break;
  }

  MEIntegratorNew* meIntegrator = new MEIntegratorNew( pathToTF , par, int(verbose));
  if(mode == 0)
    meIntegrator->setIntType( MEIntegratorNew::SL2wj );
  else if(mode == 1)
    meIntegrator->setIntType( MEIntegratorNew::SL1wj );
  else if(mode == 2)
    meIntegrator->setIntType( MEIntegratorNew::SLNoBHad );
  else if(mode == 3)
    meIntegrator->setIntType( MEIntegratorNew::SLNoBLep );
  else if(mode == 4)
    meIntegrator->setIntType( MEIntegratorNew::SLNoHiggs );
  else if(mode == 5)
    meIntegrator->setIntType( MEIntegratorNew::SL3b );
  else if(mode == 6)
    meIntegrator->setIntType( MEIntegratorNew::DL );
  else if(mode == 7)
    meIntegrator->setIntType( MEIntegratorNew::SLXSec );
  else if(mode == 8)
    meIntegrator->setIntType( MEIntegratorNew::SLAcc );
  else if(mode == 9)
    meIntegrator->setIntType( MEIntegratorNew::SLAcc2wj );
  else if(mode == 10)
    meIntegrator->setIntType( MEIntegratorNew::SLAcc1wj );
  else if(mode == 11)
    meIntegrator->setIntType( MEIntegratorNew::SLAccNoBHad );
  else if(mode == 12)
    meIntegrator->setIntType( MEIntegratorNew::SLAccNoBLep );
  else if(mode == 13)
    meIntegrator->setIntType( MEIntegratorNew::SLAccNoHiggs );
  else{
    cout << "Unsupported mode... exit" << endl;
    delete meIntegrator;
    return 0;
  }
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

  meIntegrator->setHypo(hypo);

  meIntegrator->initTFparameters(1.0,scaleL,1.0,scaleH, scaleMET);

  if(switchoffOL){
    meIntegrator->switchOffOL(); 
    cout << "*** Switching off OpenLoops to speed-up the calculation ***" << endl;
  }


  if(doMassScan && mode==5){
    cout << "Mode 5 works only when testing the permutations... exit" << endl;
    delete meIntegrator;
    return 0;
  }


  ////////////////////////////////////////////////////////
  if(mode>6){

    double CTl = TMath::Cos(2*TMath::ATan(TMath::Exp(-2.1)));
    double CTj = TMath::Cos(2*TMath::ATan(TMath::Exp(-2.5)));

    cout << "Mode " << mode << ": cos theta l (j) will be restricted to +/-" << CTl << " (" << CTj << ")" << endl; 

    for(int m = 0; m < nMassPoints ; m++){

      ROOT::Math::Functor toIntegrate(meIntegrator, &MEIntegratorNew::Eval, par);
      ROOT::Math::GSLMCIntegrator ig2( ROOT::Math::IntegrationMultiDim::kVEGAS , 1.e-12, 1.e-5, vegasPoints);
      ig2.SetFunction(toIntegrate);
      meIntegrator->SetPar(par+hypo);

      double xLmode[19]    = {    -1, -PI,    5.,-1.,-PI,-1, -PI,   5.,  -1, -PI, -1, -PI, -1, -PI,    5,-1, -PI, -1, -PI};
      double xUmode[19]    = {     1., PI,  500., 1., PI, 1,  PI, 500.,   1,  PI,  1,  PI,  1,  PI,  500, 1,  PI,  1,  PI};
      double xLmode_b[20]  = {    -1, -PI,    5.,-1.,-PI,-1, -PI,   5.,  -1, -PI, -1, -PI, -1, -PI,    5,-1, -PI, -1, -PI,   5};
      double xUmode_b[20]  = {     1., PI,  500., 1., PI, 1,  PI, 500.,   1,  PI,  1,  PI,  1,  PI,  500, 1,  PI,  1,  PI, 500};

      double xLmode9[19]   = {  -CTl, -PI,   20.,-1.,-PI,-CTj, -PI,  20.,  -CTj, -PI, -CTj, -PI, -CTj, -PI,   20,-CTj, -PI, -CTj, -PI};
      double xUmode9[19]   = {   CTl,  PI,  500., 1., PI, CTj,  PI, 500.,   CTj,  PI,  CTj,  PI,  CTj,  PI,  500, CTj,  PI,  CTj,  PI};
      double xLmode9_b[20] = {  -CTl, -PI,   20.,-1.,-PI,-CTj, -PI,  20.,  -CTj, -PI, -CTj, -PI, -CTj, -PI,   20,-CTj, -PI, -CTj, -PI,  20};
      double xUmode9_b[20] = {   CTl,  PI,  500., 1., PI, CTj,  PI, 500.,   CTj,  PI,  CTj,  PI,  CTj,  PI,  500, CTj,  PI,  CTj,  PI, 500};

      double xLmode10[19]  = {  -CTl, -PI,   20.,-1.,-PI,-CTj, -PI,   5.,  -1, -PI, -1, -PI, -CTj, -PI,   20,-CTj, -PI, -CTj, -PI};
      double xUmode10[19]  = {   CTl,  PI,  500., 1., PI, CTj,  PI, 500.,   1,  PI,  1,  PI,  CTj,  PI,  500, CTj,  PI,  CTj,  PI};

      double xLmode11[19]  = {  -CTl, -PI,   20.,-1.,-PI,-CTj, -PI,  20.,  -CTj, -PI, -CTj, -PI, -1, -PI,   20,-CTj, -PI, -CTj, -PI};
      double xUmode11[19]  = {   CTl,  PI,  500., 1., PI, CTj,  PI, 500.,   CTj,  PI,  CTj,  PI,  1,  PI,  500, CTj,  PI,  CTj,  PI};

      double xLmode12[19]  = {  -CTl, -PI,   20.,-1.,-PI,-1, -PI,    20.,  -CTj, -PI, -CTj, -PI, -CTj, -PI,   20,-CTj, -PI, -CTj, -PI};
      double xUmode12[19]  = {   CTl,  PI,  500., 1., PI, 1,  PI,   500.,   CTj,  PI,  CTj,  PI,  CTj,  PI,  500, CTj,  PI,  CTj,  PI};

      double xLmode13[19]  = {  -CTl, -PI,   20.,-1.,-PI,-1, -PI,    20.,  -CTj, -PI, -CTj, -PI, -CTj, -PI,   20,-1, -PI, -1, -PI};
      double xUmode13[19]  = {   CTl,  PI,  500., 1., PI, 1,  PI,   500.,   CTj,  PI,  CTj,  PI,  CTj,  PI,  500, 1,  PI,  1,  PI};

      meIntegrator->setMass( mH[m] );

      clock->Start();
      double p = 0.;
      if(mode == 7 || mode == 8)
	p = ig2.Integral(xLmode, xUmode);
      else if(mode == 9)
	p = hypo==0 ? ig2.Integral(xLmode9, xUmode9) : ig2.Integral(xLmode9_b, xUmode9_b);
      else if(mode == 10)
	p = ig2.Integral(xLmode10, xUmode10);
      else if(mode == 11)
	p = ig2.Integral(xLmode11, xUmode11);
      else if(mode == 12)
	p = ig2.Integral(xLmode12, xUmode12);
      else if(mode == 13)
	p = ig2.Integral(xLmode13, xUmode13);
      else{}

      clock->Stop();

      cout << "Mass  " << mH[m] << " => p = " << p << endl;
      hInt->Fill( mH[m], p);
      cout << "[Done in " << clock->RealTime()/60. << " min]" << endl;
    }

    cout << "Delete meIntegrator..." << endl;
    delete meIntegrator;
    delete clock; delete ran;
    cout << "Finished!!!" << endl;
    
    TFile* fout_tmp0 = TFile::Open(outFileName.c_str(),"UPDATE");
    hInt->Write("hInt",TObject::kOverwrite);
    fout_tmp0->Close();
    return 0;
  }
  ////////////////////////////////////////////////////////


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

    genParticleInfo genB, genBbar;
    genTopInfo genTop, genTbar;
    metInfo METtype1p2corr;
    JetByPt jet1, jet2, jet3, jet4, jet5, jet6, jet7, jet8, jet9, jet10;
    Float_t vLepton_charge [99];

    currentTree->SetBranchAddress("genB",   &genB);
    currentTree->SetBranchAddress("genBbar",&genBbar);
    currentTree->SetBranchAddress("genTop", &genTop);
    currentTree->SetBranchAddress("genTbar",&genTbar);
    currentTree->SetBranchAddress("METtype1p2corr",&METtype1p2corr);
    currentTree->SetBranchAddress("vLepton_charge",vLepton_charge);

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



    int counter = 0;
    Long64_t nentries = currentTree->GetEntries();
    for (Long64_t i = 0; i < nentries ; i++){
      
      if( counter > evHigh ) continue;

      if(i%5000==0) cout << i << endl;
      currentTree->GetEntry(i);
      
      //continue;
      
      //////////////////////////////////////////////////////////////////
      vector<TLorentzVector> myJets;
      std::map<unsigned int, JetByPt> map;
      vector<TLorentzVector> myJetsFilt;
      std::map<unsigned int, JetByPt> mapFilt;
      
      map[0] = jet1;  map[1] = jet2; map[2] = jet3; map[3] = jet4; 
      map[4] = jet5;  map[5] = jet6; map[6] = jet7; map[7] = jet8; 
      map[8] = jet9;  map[9] = jet10; 
    

      for(unsigned int m = 0; m <10; m++){
	TLorentzVector jetLV; 
	if(map[m].pt>20){
	  jetLV.SetPtEtaPhiM( map[m].pt, map[m].eta, map[m].phi, map[m].mass );
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
      //////////////////////////////////////////////////////////////////////////


      TLorentzVector topBLV   ;
      TLorentzVector topW1LV  ; 
      TLorentzVector topW2LV  ;
      TLorentzVector atopBLV  ; 
      TLorentzVector atopW1LV ; 
      TLorentzVector atopW2LV ;
      TLorentzVector genBLV   ;
      TLorentzVector genBbarLV;

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
  

      //////////////////////////////////////////////////////////////////////////
      vector<TLorentzVector> genveto;
      genveto.push_back( topBLV );
      genveto.push_back( topW1LV);
      genveto.push_back( topW2LV);
      genveto.push_back( atopBLV);
      genveto.push_back( atopW1LV);
      genveto.push_back( atopW2LV);
      unsigned int untag1 = 999;
      unsigned int untag2 = 999;
      for(unsigned int k = 0; k < myJetsFilt.size(); k++){ 
	bool match = false;
	for(unsigned int l = 0; l < genveto.size(); l++){  
	  if( deltaR(myJetsFilt[k], genveto[l])<0.5 ) match = true;
	}
	if(!match && untag1==999 && untag2==999)      untag1 = k;
	else if(!match && untag1!=999 && untag2==999) untag2 = k;
      }
      if( doBkgWithWorkaraound && untag1!=999 && untag2!=999 ){
	genBLV    =  myJetsFilt[untag1];
	genBbarLV =  myJetsFilt[untag2];
      }
      //////////////////////////////////////////////////////////////////////////


      TLorentzVector TOPLEP;
      TLorentzVector TOPHAD;
      TLorentzVector TOPHADW1;
      TLorentzVector TOPHADW2;
      TLorentzVector TOPHADB;
      TLorentzVector TOPLEPW1;
      TLorentzVector TOPLEPW2;
      TLorentzVector TOPLEPB;
      TLorentzVector HIGGS;

      float neutCosTheta;
      int chargeTopLep = 0;
      bool isSL = false;
      bool isDL = false;

      bool properEvent = (genBLV.Pt()>0 && genBbarLV.Pt()>0 && topBLV.Pt()>0 && topW1LV.Pt()>0 && topW2LV.Pt()>0 && atopBLV.Pt()>0 && atopW1LV.Pt()>0 && atopW2LV.Pt()>0);
      if(!properEvent) continue;

      HIGGS.SetPxPyPzE( (genBLV+genBbarLV).Px(), (genBLV+genBbarLV).Py(),(genBLV+genBbarLV).Pz(),(genBLV+genBbarLV).E());
      if( abs(genTop.wdau1id)>6 && abs(genTbar.wdau1id)<6){
	TOPLEP.SetPxPyPzE( (topBLV+topW1LV+topW2LV).Px(), (topBLV+topW1LV+topW2LV).Py(), (topBLV+topW1LV+topW2LV).Pz(), (topBLV+topW1LV+topW2LV).E() );
	TOPHAD.SetPxPyPzE( (atopBLV+atopW1LV+atopW2LV).Px(), (atopBLV+atopW1LV+atopW2LV).Py(), (atopBLV+atopW1LV+atopW2LV).Pz(), (atopBLV+atopW1LV+atopW2LV).E() );
	TOPHADW1.SetPxPyPzE( atopW1LV.Px(), atopW1LV.Py(), atopW1LV.Pz(), atopW1LV.E());
	TOPHADW2.SetPxPyPzE( atopW2LV.Px(), atopW2LV.Py(), atopW2LV.Pz(), atopW2LV.E());
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
	properEvent=false;
      }

      if(mode==6) TOPLEPW2 = TOPLEPW2+TOPHADW2;

      if(mode==0)
	properEvent = isSL && ( TOPLEPW1.Pt() >20 && TMath::Abs(TOPLEPW1.Eta()) <2.1 &&
				TOPLEPB.Pt()  >30 && TMath::Abs(TOPLEPB.Eta())  <2.5 &&
				TOPHADW1.Pt() >30 && TMath::Abs(TOPHADW1.Eta()) <2.5 &&
				TOPHADW2.Pt() >30 && TMath::Abs(TOPHADW2.Eta()) <2.5 &&
				TOPHADB.Pt()  >30 && TMath::Abs(TOPHADB.Eta())  <2.5 &&
				genBLV.Pt()   >30 && TMath::Abs(genBLV.Eta())   <2.5 &&
				genBbarLV.Pt()>30 && TMath::Abs(genBbarLV.Eta())<2.5
				);
      if(mode==1)
	properEvent = isSL && (( TOPLEPW1.Pt() >20 && TMath::Abs(TOPLEPW1.Eta()) <2.1 &&
				 TOPLEPB.Pt()  >30 && TMath::Abs(TOPLEPB.Eta())  <2.5 &&
				 (TOPHADW1.Pt() < 30 || TMath::Abs(TOPHADW1.Eta()) >2.5) &&
				 TOPHADW2.Pt() >30 && TMath::Abs(TOPHADW2.Eta()) <2.5 &&
				 TOPHADB.Pt()  >30 && TMath::Abs(TOPHADB.Eta())  <2.5 &&
				 genBLV.Pt()   >30 && TMath::Abs(genBLV.Eta())   <2.5 &&
				 genBbarLV.Pt()>30 && TMath::Abs(genBbarLV.Eta())<2.5
				 ) ||
			       ( TOPLEPW1.Pt() >20 && TMath::Abs(TOPLEPW1.Eta()) <2.1 &&
				 TOPLEPB.Pt()  >30 && TMath::Abs(TOPLEPB.Eta())  <2.5 &&
				 TOPHADW1.Pt() >30 && TMath::Abs(TOPHADW1.Eta()) <2.5 &&
				 (TOPHADW2.Pt() < 30 || TMath::Abs(TOPHADW2.Eta()) >2.5) &&
				 TOPHADB.Pt()  >30 && TMath::Abs(TOPHADB.Eta())  <2.5 &&
				 genBLV.Pt()   >30 && TMath::Abs(genBLV.Eta())   <2.5 &&
				 genBbarLV.Pt()>30 && TMath::Abs(genBbarLV.Eta())<2.5
				 ) );
      if(mode==2){
	properEvent = isSL && ( TOPLEPW1.Pt() >20 && TMath::Abs(TOPLEPW1.Eta()) <2.1 &&
				TOPLEPB.Pt()  >30 && TMath::Abs(TOPLEPB.Eta())  <2.5 &&
				TOPHADW1.Pt() >30 && TMath::Abs(TOPHADW1.Eta()) <2.5 &&
				TOPHADW2.Pt() >30 && TMath::Abs(TOPHADW2.Eta()) <2.5 &&
				(TOPHADB.Pt() <30 || TMath::Abs(TOPHADB.Eta())  >2.5)&&
				genBLV.Pt()   >30 && TMath::Abs(genBLV.Eta())   <2.5 &&
				genBbarLV.Pt()>30 && TMath::Abs(genBbarLV.Eta())<2.5
				);
      }
      if(mode==3){
	properEvent = isSL && ( TOPLEPW1.Pt() >20 && TMath::Abs(TOPLEPW1.Eta()) <2.1 &&
				(TOPLEPB.Pt() <30 || TMath::Abs(TOPLEPB.Eta())  >2.5)&&
				TOPHADW1.Pt() >30 && TMath::Abs(TOPHADW1.Eta()) <2.5 &&
				TOPHADW2.Pt() >30 && TMath::Abs(TOPHADW2.Eta()) <2.5 &&
				TOPHADB.Pt()  >30 && TMath::Abs(TOPHADB.Eta())  <2.5 &&
				genBLV.Pt()   >30 && TMath::Abs(genBLV.Eta())   <2.5 &&
				genBbarLV.Pt()>30 && TMath::Abs(genBbarLV.Eta())<2.5
				);
      }
      if(mode==4){
	properEvent = isSL && (( TOPLEPW1.Pt() >20 && TMath::Abs(TOPLEPW1.Eta()) <2.1 &&
				 TOPLEPB.Pt()  >30 && TMath::Abs(TOPLEPB.Eta())  <2.5 &&
				 TOPHADW1.Pt() >30 && TMath::Abs(TOPHADW1.Eta()) <2.5 &&
				 TOPHADW2.Pt() >30 && TMath::Abs(TOPHADW2.Eta()) <2.5 &&				
				 TOPHADB.Pt()  >30 && TMath::Abs(TOPHADB.Eta())  <2.5 &&		       
				 (genBLV.Pt()  < 30 || TMath::Abs(genBLV.Eta())  >2.5)  &&
				 genBbarLV.Pt()>30 && TMath::Abs(genBbarLV.Eta())<2.5 
				 ) ||
			       ( TOPLEPW1.Pt() >20 && TMath::Abs(TOPLEPW1.Eta()) <2.1 &&
				 TOPLEPB.Pt()  >30 && TMath::Abs(TOPLEPB.Eta())  <2.5 &&
				 TOPHADW1.Pt() >30 && TMath::Abs(TOPHADW1.Eta()) <2.5 &&
				 TOPHADW2.Pt() >30 && TMath::Abs(TOPHADW2.Eta()) <2.5 &&
				 TOPHADB.Pt()  >30 && TMath::Abs(TOPHADB.Eta())  <2.5 &&			
				 genBLV.Pt()   >30 && TMath::Abs(genBLV.Eta())   <2.5 &&
				 (genBbarLV.Pt()< 30 || TMath::Abs(genBbarLV.Eta())>2.5)
				 ) );
      }
      if(mode==5){
	properEvent = isSL && (( TOPLEPW1.Pt() >20 && TMath::Abs(TOPLEPW1.Eta()) <2.1 &&
				TOPLEPB.Pt()  >30 && TMath::Abs(TOPLEPB.Eta())  <2.5 &&
				TOPHADW1.Pt() >30 && TMath::Abs(TOPHADW1.Eta()) <2.5 &&
				TOPHADW2.Pt() >30 && TMath::Abs(TOPHADW2.Eta()) <2.5 &&
				(TOPHADB.Pt() <30 || TMath::Abs(TOPHADB.Eta())  >2.5)&&
				genBLV.Pt()   >30 && TMath::Abs(genBLV.Eta())   <2.5 &&
				genBbarLV.Pt()>30 && TMath::Abs(genBbarLV.Eta())<2.5
				 ) ||
			       ( TOPLEPW1.Pt() >20 && TMath::Abs(TOPLEPW1.Eta()) <2.1 &&
				(TOPLEPB.Pt() <30 || TMath::Abs(TOPLEPB.Eta())  >2.5)&&
				 TOPHADW1.Pt() >30 && TMath::Abs(TOPHADW1.Eta()) <2.5 &&
				 TOPHADW2.Pt() >30 && TMath::Abs(TOPHADW2.Eta()) <2.5 &&
				 TOPHADB.Pt()  >30 && TMath::Abs(TOPHADB.Eta())  <2.5 &&
				 genBLV.Pt()   >30 && TMath::Abs(genBLV.Eta())   <2.5 &&
				 genBbarLV.Pt()>30 && TMath::Abs(genBbarLV.Eta())<2.5
				 ) ||
			       (( TOPLEPW1.Pt() >20 && TMath::Abs(TOPLEPW1.Eta()) <2.1 &&
				  TOPLEPB.Pt()  >30 && TMath::Abs(TOPLEPB.Eta())  <2.5 &&
				  TOPHADW1.Pt() >30 && TMath::Abs(TOPHADW1.Eta()) <2.5 &&
				  TOPHADW2.Pt() >30 && TMath::Abs(TOPHADW2.Eta()) <2.5 &&				
				  TOPHADB.Pt()  >30 && TMath::Abs(TOPHADB.Eta())  <2.5 &&		       
				  (genBLV.Pt()  < 30 || TMath::Abs(genBLV.Eta())  >2.5)  &&
				  genBbarLV.Pt()>30 && TMath::Abs(genBbarLV.Eta())<2.5 
				  ) ||
				( TOPLEPW1.Pt() >20 && TMath::Abs(TOPLEPW1.Eta()) <2.1 &&
				  TOPLEPB.Pt()  >30 && TMath::Abs(TOPLEPB.Eta())  <2.5 &&
				  TOPHADW1.Pt() >30 && TMath::Abs(TOPHADW1.Eta()) <2.5 &&
				  TOPHADW2.Pt() >30 && TMath::Abs(TOPHADW2.Eta()) <2.5 &&
				  TOPHADB.Pt()  >30 && TMath::Abs(TOPHADB.Eta())  <2.5 &&			
				  genBLV.Pt()   >30 && TMath::Abs(genBLV.Eta())   <2.5 &&
				  (genBbarLV.Pt()< 30 || TMath::Abs(genBbarLV.Eta())>2.5)
				  ) )
			       );
      }

      if(mode==6)
	properEvent = isDL && ( TOPLEPW1.Pt() >20 && TMath::Abs(TOPLEPW1.Eta()) <2.1 &&
				TOPLEPB.Pt()  >30 && TMath::Abs(TOPLEPB.Eta())  <2.5 &&
				TOPHADW1.Pt() >20 && TMath::Abs(TOPHADW1.Eta()) <2.1 &&
				TOPHADB.Pt()  >30 && TMath::Abs(TOPHADB.Eta())  <2.5 &&
				genBLV.Pt()   >30 && TMath::Abs(genBLV.Eta())   <2.5 &&
				genBbarLV.Pt()>30 && TMath::Abs(genBbarLV.Eta())<2.5
				);
      

      if(!properEvent) continue;
      counter++;

      if( !(counter>=evLow && counter<=evHigh) ) continue;
      cout << "Processing event # " << counter << endl;

      counter_ = counter;
      type_    = mode;

      vector<TLorentzVector> jets;
      jets.push_back( TOPLEPW1  );  
      jets.push_back( TOPLEPW2  );  
      jets.push_back( TOPLEPB   );
      jets.push_back( TOPHADW1  );  
      jets.push_back( TOPHADW2  );  
      jets.push_back( TOPHADB   );
      jets.push_back( genBLV    );
      jets.push_back( genBbarLV );

      int swapW = 0;
      int swapH = 0;
      if(mode==1){
	TLorentzVector w1 = jets[3];
	TLorentzVector w2 = jets[4];
	if( (w1.Pt()<30 || TMath::Abs(w1.Eta())>2.5) && (w2.Pt()>30 && TMath::Abs(w2.Eta())<2.5)){
	  jets[3] = w2;
	  jets[4] = w1;
	  swapW = 1;
	}
	else if( (w2.Pt()<30 || TMath::Abs(w2.Eta())>2.5) && (w1.Pt()>30 && TMath::Abs(w1.Eta())<2.5)){
	  jets[3] = w1;
	  jets[4] = w2;
	}
	else{
	  cout << "Inconsistentcy of mode and jet selections" << endl;
	  delete meIntegrator;
	  return 1;
	}
      }

      if(mode==4){
	TLorentzVector h1 = jets[6];
	TLorentzVector h2 = jets[7];
	if( (h1.Pt()<30 || TMath::Abs(h1.Eta())>2.5) && (h2.Pt()>30 && TMath::Abs(h2.Eta())<2.5)){
	  jets[6] = h2;
	  jets[7] = h1;
	  swapH = 1;
	}
	else if( (h2.Pt()<30 || TMath::Abs(h2.Eta())>2.5) && (h1.Pt()>30 && TMath::Abs(h1.Eta())<2.5)){
	  jets[6] = h1;
	  jets[7] = h2;
	}
	else{
	  cout << "Inconsistentcy of mode and jet selections" << endl;
	  delete meIntegrator;
	  return 1;
	}
      }

      if(mode==5){
	TLorentzVector h1 = jets[6];
	TLorentzVector h2 = jets[7];
	if( (h1.Pt()<30 || TMath::Abs(h1.Eta())>2.5) && (h2.Pt()>30 && TMath::Abs(h2.Eta())<2.5)){
	  jets[6] = h2;
	  jets[7] = h1;
	  swapH = 1;
	}
	else if( (h2.Pt()<30 || TMath::Abs(h2.Eta())>2.5) && (h1.Pt()>30 && TMath::Abs(h1.Eta())<2.5)){
	  jets[6] = h1;
	  jets[7] = h2;
	}
	else{}
      }


      
      if(printP4){
	cout << "******* INPUT ******" << endl;      
	cout << "SumEt = " << METtype1p2corr.sumet  << endl;
	cout << "lep:  jet1.SetPtEtaPhiM(" << TOPLEPW1.Pt() << "," <<  TOPLEPW1.Eta()  << "," << TOPLEPW1.Phi() << "," << TOPLEPW1.M() << ")" << endl;
	cout << "met:  jet2.SetPtEtaPhiM(" << TOPLEPW2.Pt() << "," <<  TOPLEPW2.Eta()  << "," << TOPLEPW2.Phi() << "," << TOPLEPW2.M() << ")" << endl;
	cout << "blep: jet3.SetPtEtaPhiM(" << TOPLEPB.Pt() << ","  <<  TOPLEPB.Eta()   << "," << TOPLEPB.Phi()  << "," <<  TOPLEPB.M() << ")" << endl;
	if(swapW){
	  cout << "w2:   jet5.SetPtEtaPhiM(" << TOPHADW2.Pt() << "," <<  TOPHADW2.Eta()  << "," << TOPHADW2.Phi() << "," << TOPHADW2.M() << ")" << endl;
	  cout << "w1:   jet4.SetPtEtaPhiM(" << TOPHADW1.Pt() << "," <<  TOPHADW1.Eta()  << "," << TOPHADW1.Phi() << "," << TOPHADW1.M() << ")" << endl;
	}
	else{
	  cout << "w1:   jet4.SetPtEtaPhiM(" << TOPHADW1.Pt() << "," <<  TOPHADW1.Eta()  << "," << TOPHADW1.Phi() << "," << TOPHADW1.M() << ")" << endl;
	  cout << "w2:   jet5.SetPtEtaPhiM(" << TOPHADW2.Pt() << "," <<  TOPHADW2.Eta()  << "," << TOPHADW2.Phi() << "," << TOPHADW2.M() << ")" << endl;
	}
	cout << "bhad: jet6.SetPtEtaPhiM(" << TOPHADB.Pt() << ","  <<  TOPHADB.Eta()   << "," << TOPHADB.Phi() << "," << TOPHADB.M() << ")" << endl;
	if(swapH){
	  cout << "h2:   jet8.SetPtEtaPhiM(" << genBbarLV.Pt() << ","<<  genBbarLV.Eta() << "," << genBbarLV.Phi() << "," << genBbarLV.M() << ")" << endl;	
	  cout << "h1:   jet7.SetPtEtaPhiM(" << genBLV.Pt() << ","   <<  genBLV.Eta()    << "," << genBLV.Phi() << "," << genBLV.M() << ")" << endl;
	}
	else{
	  cout << "h1:   jet7.SetPtEtaPhiM(" << genBLV.Pt() << ","   <<  genBLV.Eta()    << "," << genBLV.Phi() << "," << genBLV.M() << ")" << endl;
	  cout << "h2:   jet8.SetPtEtaPhiM(" << genBbarLV.Pt() << ","<<  genBbarLV.Eta() << "," << genBbarLV.Phi() << "," << genBbarLV.M() << ")" << endl;	
	}
	cout << "Top Lep mass = " << (TOPLEPW1+TOPLEPW2+TOPLEPB).M() << " <=> neutrino eta=0!!!" << endl;
	cout << "Top Had mass = " << (TOPHADW1+TOPHADW2+TOPHADB).M() << endl;
	cout << "Higgs mass = "   << (genBLV+genBbarLV).M() << endl;
	cout << "Neut cosTheta = " << neutCosTheta << endl;
	}


      meIntegrator->setJets(&jets);
      meIntegrator->setSumEt( METtype1p2corr.sumet );
      meIntegrator->setMEtCov(-99,-99,0);
      //meIntegrator->setTopFlags( vLepton_charge[0]==1 ? +1 : -1 , vLepton_charge[0]==1 ? -1 : +1 );
      meIntegrator->setTopFlags( chargeTopLep>0 ? -1 : +1 , chargeTopLep>0 ? +1 : -1);

      probGG_      = 0;
      probAtSgnGG_ = 0;
      probAtSgnGG_alt_ = 0;
      massGG_      = 0;
      evalCpuGG_   = 0;
      evalReaGG_   = 0;
      evalCpuGG_alt_   = 0;
      evalReaGG_alt_   = 0;
      probGA_      = 0;
      massGA_      = 0;
      massIntGA_   = 0;
      massIntModGA_= 0;
      massMaxGA_   = 0;
      massMaxModGA_= 0;
      probIntGA_   = 0;
      probIntModGA_= 0;
      probMaxGA_   = 0;
      probMaxModGA_= 0;
      evalCpuGA_   = 0;
      evalReaGA_   = 0;
      evalCpuGA_alt_   = 0;
      evalReaGA_alt_   = 0;
      probAtSgnGA_ = 0;
      probAtSgnGA_alt_ = 0;
      probAtGoodGA_= 0;
      probRG_      = 0;
      probAtSgnRG_ = 0;
      probAtSgnRG_alt_ = 0;
      massRG_      = 0;
      evalCpuRG_   = 0;
      evalReaRG_   = 0;	  
      evalCpuRG_alt_   = 0;
      evalReaRG_alt_   = 0;	  	
      probRA_      = 0;
      massRA_      = 0;
      massIntRA_   = 0;
      massIntModRA_= 0;
      massMaxRA_   = 0;
      massMaxModRA_= 0;
      probIntRA_   = 0;
      probIntModRA_= 0;
      probMaxRA_   = 0;
      probMaxModRA_= 0;
      evalCpuRA_   = 0;
      evalReaRA_   = 0;
      evalCpuRA_alt_   = 0;
      evalReaRA_alt_   = 0;
      probAtSgnRA_ = 0;
      probAtSgnRA_alt_ = 0;
      probAtGoodRA_= 0;
      matchByIntGA_   = 0;
      matchByIntModGA_= 0;
      matchByMaxGA_   = 0;
      matchByMaxModGA_= 0;
      matchByIntRA_   = 0;
      matchByIntModRA_= 0;
      matchByMaxRA_   = 0;
      matchByMaxModRA_= 0;
      chi2GG_ = -99;
      chi2RG_ = -99;

      dPx_    = -999;
      dRecPx_ = -999;

      if( doParton){

	if(mode == 5){
	  meIntegrator->setIntType( MEIntegratorNew::SL3b );
	  TLorentzVector bLep =  meIntegrator->jetAt(2);
	  TLorentzVector bHad =  meIntegrator->jetAt(5);
	  TLorentzVector hig1 =  meIntegrator->jetAt(6);
	  TLorentzVector hig2 =  meIntegrator->jetAt(7);
	  if( (bLep.Pt()<30 || TMath::Abs(bLep.Eta())>2.5) && (bHad.Pt()>30 && TMath::Abs(bHad.Eta())<2.5) && (hig1.Pt()>30 && TMath::Abs(hig1.Eta())<2.5)  && (hig2.Pt()>30 && TMath::Abs(hig2.Eta())<2.5)){
	    type_ = 3;
	  }
	  else if( (bLep.Pt()>30 && TMath::Abs(bLep.Eta())<2.5) && (bHad.Pt()<30 || TMath::Abs(bHad.Eta())>2.5) && (hig1.Pt()>30 && TMath::Abs(hig1.Eta())<2.5)  && (hig2.Pt()>30 && TMath::Abs(hig2.Eta())<2.5)){
	    type_ = 2;
	  }
	  else if( (bLep.Pt()>30 && TMath::Abs(bLep.Eta())<2.5) && (bHad.Pt()>30 && TMath::Abs(bHad.Eta())<2.5) && (hig1.Pt()>30 && TMath::Abs(hig1.Eta())<2.5)  && (hig2.Pt()<30 || TMath::Abs(hig2.Eta())>2.5)){
	    type_ = 4;
	  }
	  else if((bLep.Pt()>30 && TMath::Abs(bLep.Eta())<2.5) && (bHad.Pt()>30 && TMath::Abs(bHad.Eta())<2.5) && (hig1.Pt()<30 || TMath::Abs(hig1.Eta())>2.5)  && (hig2.Pt()>30 && TMath::Abs(hig2.Eta())<2.5)){
	    type_ = 4;
	  }
	  else{}
	}


	if( doMassScan ){
	
	  meIntegrator->initVersors(1);
		
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

	  double xLmode2[6]   = {x0L, x1L, x2L, x3L, x4L, x5L};
	  double xUmode2[6]   = {x0U, x1U, x2U, x3U, x4U, x5U};
	  double xLmode2_b[7] = {x0L, x1L, x2L, x3L, x4L, x5L, x6L};
	  double xUmode2_b[7] = {x0U, x1U, x2U, x3U, x4U, x5U, x6U};

	  double xLmode3[6]   = {x0L, x1L, x2L, x3L, x4L, x5L};
	  double xUmode3[6]   = {x0U, x1U, x2U, x3U, x4U, x5U};
	  double xLmode3_b[7] = {x0L, x1L, x2L, x3L, x4L, x5L, x6L};
	  double xUmode3_b[7] = {x0U, x1U, x2U, x3U, x4U, x5U, x6U};

	  double xLmode4[6]   = {x0L, x1L, x2L, x3L, x4L, x5L};
	  double xUmode4[6]   = {x0U, x1U, x2U, x3U, x4U, x5U};
	  double xLmode4_b[7] = {x0L, x1L, x2L, x3L, x4L, x5L, x6L};
	  double xUmode4_b[7] = {x0U, x1U, x2U, x3U, x4U, x5U, x6U};

	  double xLmode6[5]   = {x1L, x2L, x1L, x2L, x5L};
	  double xUmode6[5]   = {x1U, x2U, x1U, x2U, x5U};
	  double xLmode6_b[6] = {x1L, x2L, x1L, x2L, x5L, x6L};
	  double xUmode6_b[6] = {x1U, x2U, x1U, x2U, x5U, x6U};
	  
	  if( printP4 ){
	    cout << "Integration range: " << endl;
	    
	    
	    for(unsigned int k = 0; k < (unsigned int)par; k++){
	      if(mode==0) 
		cout << "Var " << k << ": [" << (hypo==0 ? xLmode0[k] : xLmode0_b[k]) << "," <<  (hypo==0 ? xUmode0[k] : xUmode0_b[k]) << "]" << endl;
	      else if(mode==1)
		cout << "Var " << k << ": [" << (hypo==0 ? xLmode1[k] : xLmode1_b[k])<< "," <<  (hypo==0 ? xUmode1[k] : xUmode1_b[k])<< "]" << endl;
	      else if(mode==2)
		cout << "Var " << k << ": [" << xLmode2[k] << "," <<  xUmode2[k] << "]" << endl;
	      else if(mode==3)
		cout << "Var " << k << ": [" << xLmode3[k] << "," <<  xUmode3[k] << "]" << endl;
	      else if(mode==4)
		cout << "Var " << k << ": [" << xLmode4[k] << "," <<  xUmode4[k] << "]" << endl;
	      else if(mode==6)
		cout << "Var " << k << ": [" << xLmode6[k] << "," <<  xUmode6[k] << "]" << endl;
	      else{}
	    }	   
	  }
	  
	  hMass_gen->Reset();
	  hMass_alt_gen->Reset();
	  for(int m = 0; m < nMassPoints ; m++){

	    meIntegrator->setMass( mH[m] );

	    clock->Start();

	    double p = 0.;

	    for(int hyp = 0 ; hyp<2; hyp++){

	      if(SoB==0 && hyp!=hypo) continue;
	      if( hyp==1 && false/*&& !(mode==0 || mode==1)*/){
		cout << "Unsupported at the moment..." << endl;
		continue;
	      }

	      if(printP4) cout << "Testing hypothesis " << (hyp==0 ? "S" : "B") << endl; 
	      meIntegrator->setHypo(hyp);

	      int ntries = 0;
	      int intPoints = vegasPoints;
	      double err = 0;
	      while( ntries < MAX_REEVAL_TRIES){

		ROOT::Math::Functor toIntegrate(meIntegrator, &MEIntegratorNew::Eval, par+hyp);
		ROOT::Math::GSLMCIntegrator ig2( ROOT::Math::IntegrationMultiDim::kVEGAS , 1.e-12, 1.e-5, intPoints );
		ig2.SetFunction(toIntegrate);
		meIntegrator->SetPar(par+hyp);	 
		
		if(mode==0)
		  p = (hyp==0 ? ig2.Integral(xLmode0, xUmode0) : ig2.Integral(xLmode0_b, xUmode0_b));
		else if(mode==1)
		  p = (hyp==0 ? ig2.Integral(xLmode1, xUmode1) : ig2.Integral(xLmode1_b, xUmode1_b));
		else if(mode==2)
		  p = (hyp==0 ? ig2.Integral(xLmode2, xUmode2) : ig2.Integral(xLmode2_b, xUmode2_b));
		else if(mode==3)
		  p = (hyp==0 ? ig2.Integral(xLmode3, xUmode3) : ig2.Integral(xLmode3_b, xUmode3_b)) ;
		else if(mode==4)
		  p = (hyp==0 ? ig2.Integral(xLmode4, xUmode4) : ig2.Integral(xLmode4_b, xUmode4_b));
		else if(mode==6)
		  p = (hyp==0 ? ig2.Integral(xLmode6, xUmode6) : ig2.Integral(xLmode6_b, xUmode6_b));
		else{ }
		
		cout << "Chi2 = " << ig2.ChiSqr() << endl;
		chi2GG_ =  ig2.ChiSqr();
		err = ig2.Error() ;

		if( chi2GG_ > maxChi2_ ){
		  ntries++;
		  intPoints *= 1.5;
		  cout << "VEGAS has chi2=" << chi2GG_ << " > " << maxChi2_ << ": trying again with " << intPoints << " points..." << endl;
		}
		else ntries = MAX_REEVAL_TRIES+1;
	      }

	      clock->Stop();
	      
	      if( TMath::IsNaN(p) ) p = 0.;
	      
	      if(hyp==0){
		evalCpuGG_ += clock->CpuTime();
		evalReaGG_ += clock->RealTime();
	      }
	      else{
		evalCpuGG_alt_ += clock->CpuTime();
		evalReaGG_alt_ += clock->RealTime();
	      }

	      if(printP4) cout << "Mass " << mH[m] << " => prob  = " << p << " +/- " << err << endl;
	      if(hyp==0){
		hMass_gen->Fill(mH[m], p);
	      }
	      else{
		hMass_alt_gen->Fill(mH[m], p);
	      }
	      
	      if( mH[m]<met+0.5 && mH[m]>met-0.5){
		if(hyp==0) 
		  probAtSgnGG_ += p;
		else 
		  probAtSgnGG_alt_ += p;
	      }
	      
	      p = 0.;
	      clock->Reset();
	    }


	  }
	  evalCpuGG_ /= nMassPoints;
	  evalReaGG_ /= nMassPoints;
	  evalCpuGG_alt_ /= nMassPoints;
	  evalReaGG_alt_ /= nMassPoints;
	  
	  pair<double,double> bestMass = getMaxValue(hMass_gen);
	  hBestMass_gen->Fill( bestMass.first );
	  massGG_ = bestMass.first;
	  probGG_ = bestMass.second;

	  TFile* fout_tmp = TFile::Open(outFileName.c_str(),"UPDATE");
	  if(gDirectory->GetKey(Form("Event_%d", counter)) == 0) {
	    TDirectory *dir = fout_tmp->mkdir(Form("Event_%d", counter));
	    dir->cd();
	    hMass_gen->Write("", TObject::kOverwrite);
	  }
	  else{
	    fout_tmp->cd(Form("Event_%d", counter));
	    hMass_gen->Write("", TObject::kOverwrite);
	  }
	  fout_tmp->Close();
	    
	}
	
	if(doPermutations){
	
	  hMassProb_gen->Reset();
	  hMassProb_alt_gen->Reset();
	
	  double pTotAtSgn     = 0.;
	  double pTotAtSgn_alt = 0.;
	  int combScanned = 0;
		 
	  double permutMassReco[nPermutations];
	  for(int k = 0; k < nPermutations; k++) permutMassReco[k] = 0.;

	  double permutProbInt [nPermutations];
	  for(int k = 0; k < nPermutations; k++) permutProbInt[k] = 0.;

	  double permutProbInt_alt [nPermutations];
	  for(int k = 0; k < nPermutations; k++) permutProbInt_alt[k] = 0.;

	  vector<TH1F*> histos;
	  for(unsigned int it = 0; it<(unsigned int)nPermutations; it++){
	    if(gDirectory->FindObject(Form("hMass_%d",it))!=0){
	      gDirectory->Remove(gDirectory->FindObject(Form("hMass_%d",it)));
	    }
	    TH1F* h = (TH1F*)hMass_gen->Clone(Form("hMass_%d",it));
	    h->Reset();
	    histos.push_back( h );
	  }
	  vector<TH1F*> histos_alt;
	  for(unsigned int it = 0; it<(unsigned int)nPermutations; it++){
	    if(gDirectory->FindObject(Form("hMass_alt_%d",it))!=0){
	      gDirectory->Remove(gDirectory->FindObject(Form("hMass_alt_%d",it)));
	    }
	    TH1F* h = (TH1F*)hMass_gen->Clone(Form("hMass_alt_%d",it));
	    h->Reset();
	    histos_alt.push_back( h );
	  }


	  if(gDirectory->FindObject("stack_m_all_gen")!=0){
	    gDirectory->Remove(gDirectory->FindObject("stack_m_all_gen"));
	  }
	  THStack*  sMassProbInt_gen  = new THStack("stack_m_all_gen", "Per permutation mass, GEN");
	  TLegend* leg = new TLegend(0.63,0.48,0.85,0.85,NULL,"brNDC");
	  leg->SetFillStyle(0);
	  leg->SetBorderSize(0);
	  leg->SetFillColor(10);
	  leg->SetTextSize(0.03);
	  TCanvas *c1 = new TCanvas("c_gen","canvas",10,30,650,600);

	  for(int m = 0; m < nMassPoints ; m++){
	    meIntegrator->setMass( mH[m] );

	    double maxP_s = 0.;
	    double maxP_b = 0.;
	    
	    for(unsigned int pos = 0; pos < (unsigned int)nPermutations ; pos++){
	      meIntegrator->initVersors( permutations[pos] );

	      int mode5HiggsLost = 0;
	      if( mode==5 ){ // figure out which b is missing, then set the correspinding int mode
		TLorentzVector bLep =  meIntegrator->jetAt(2);
		TLorentzVector bHad =  meIntegrator->jetAt(5);
		TLorentzVector hig1 =  meIntegrator->jetAt(6);
		TLorentzVector hig2 =  meIntegrator->jetAt(7);
		if( (bLep.Pt()<30 || TMath::Abs(bLep.Eta())>2.5) && (bHad.Pt()>30 && TMath::Abs(bHad.Eta())<2.5) && (hig1.Pt()>30 && TMath::Abs(hig1.Eta())<2.5)  && (hig2.Pt()>30 && TMath::Abs(hig2.Eta())<2.5)){
		  if(printP4) cout << "Permutation " << permutations[pos] << " has assigned a jet with (pt,eta)=(" << bLep.Pt() << "," << bLep.Eta() << ") to the bLep: setting mode SLNoBLep" << endl;
		  meIntegrator->setIntType( MEIntegratorNew::SLNoBLep );
		}
		else if( (bLep.Pt()>30 && TMath::Abs(bLep.Eta())<2.5) && (bHad.Pt()<30 || TMath::Abs(bHad.Eta())>2.5) && (hig1.Pt()>30 && TMath::Abs(hig1.Eta())<2.5)  && (hig2.Pt()>30 && TMath::Abs(hig2.Eta())<2.5)){
		   if(printP4)  cout << "Permutation " << permutations[pos] << " has assigned a jet with (pt,eta)=(" << bHad.Pt() << "," << bHad.Eta() << ") to the bHad: setting mode SLNoBHad" << endl;
		  meIntegrator->setIntType( MEIntegratorNew::SLNoBHad );
		}
		else if( (bLep.Pt()>30 && TMath::Abs(bLep.Eta())<2.5) && (bHad.Pt()>30 && TMath::Abs(bHad.Eta())<2.5) && (hig1.Pt()>30 && TMath::Abs(hig1.Eta())<2.5)  && (hig2.Pt()<30 || TMath::Abs(hig2.Eta())>2.5)){
		   if(printP4)  cout << "Permutation " << permutations[pos] << " has assigned a jet with (pt,eta)=(" << hig2.Pt() << "," << hig2.Eta() << ") to  Higgs2: setting mode SLNoHiggs" << endl;
		  meIntegrator->setIntType( MEIntegratorNew::SLNoHiggs );
		  mode5HiggsLost = 1;
		}
		else if((bLep.Pt()>30 && TMath::Abs(bLep.Eta())<2.5) && (bHad.Pt()>30 && TMath::Abs(bHad.Eta())<2.5) && (hig1.Pt()<30 || TMath::Abs(hig1.Eta())>2.5)  && (hig2.Pt()>30 && TMath::Abs(hig2.Eta())<2.5)){
		  meIntegrator->initVersors( permutations_SL3b_equivalent[pos] );
		  hig1 =  meIntegrator->jetAt(6);
		  hig2 =  meIntegrator->jetAt(7);
		   if(printP4)  cout << "Permutation " << permutations_SL3b_equivalent[pos] << " has assigned a jet with (pt,eta)=(" << hig2.Pt() << "," << hig2.Eta() << ") to  Higgs2: setting mode SLNoHiggs" << endl;
		  meIntegrator->setIntType( MEIntegratorNew::SLNoHiggs );
		  mode5HiggsLost = 1;
		}
		else{
		  cout << "Mode " << mode << ": inconsistency between mode and jet list" << endl;
		  continue;
		}
	      } 

	      
	      double mass, massLow, massHigh;
	      bool skip = !(meIntegrator->compatibilityCheck(0.95, /*printP4*/ 0, mass, massLow, massHigh )) ;
	      permutMassReco[pos] = mass;
	      if( SoB==0 && skip && !(mode==4 || (mode==5 && mode5HiggsLost))){
		continue;
	      }
	      
	      //meIntegrator->printJetList();

	      combScanned++;
	      
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

	      double xLmode2[6] = {x0L, x1L, x2L, x3L, x4L, x5L};
	      double xUmode2[6] = {x0U, x1U, x2U, x3U, x4U, x5U};
	      double xLmode2_b[7] = {x0L, x1L, x2L, x3L, x4L, x5L, x6L};
	      double xUmode2_b[7] = {x0U, x1U, x2U, x3U, x4U, x5U, x6U};

	      double xLmode3[6] = {x0L, x1L, x2L, x3L, x4L, x5L};
	      double xUmode3[6] = {x0U, x1U, x2U, x3U, x4U, x5U};
	      double xLmode3_b[7] = {x0L, x1L, x2L, x3L, x4L, x5L, x6L};
	      double xUmode3_b[7] = {x0U, x1U, x2U, x3U, x4U, x5U, x6U};

	      double xLmode4[6] = {x0L, x1L, x2L, x3L, x4L, x5L};
	      double xUmode4[6] = {x0U, x1U, x2U, x3U, x4U, x5U};
	      double xLmode4_b[7] = {x0L, x1L, x2L, x3L, x4L, x5L, x6L};
	      double xUmode4_b[7] = {x0U, x1U, x2U, x3U, x4U, x5U, x6U};

	      double xLmode6[5]   = {x1L, x2L, x1L, x2L, x5L};
	      double xUmode6[5]   = {x1U, x2U, x1U, x2U, x5U};
	      double xLmode6_b[6] = {x1L, x2L, x1L, x2L, x5L, x6L};
	      double xUmode6_b[6] = {x1U, x2U, x1U, x2U, x5U, x6U};

	      if( printP4 && false){
		cout << "Integration range: " << endl;
		for(unsigned int k = 0; k < (unsigned int)par; k++){
		  if(mode==0) 
		    cout << "Var " << k << ": [" << (hypo==0 ? xLmode0[k] : xLmode0_b[k]) << "," <<  (hypo==0 ? xUmode0[k] : xUmode0_b[k]) << "]" << endl;
		  else if(mode==1)
		    cout << "Var " << k << ": [" << (hypo==0 ? xLmode1[k] : xLmode1_b[k]) << "," <<  (hypo==0 ? xUmode1[k] : xUmode1_b[k]) << "]" << endl;
		  else if(mode==2)
		    cout << "Var " << k << ": [" << xLmode2[k] << "," <<  xUmode2[k] << "]" << endl;
		  else if(mode==3)
		    cout << "Var " << k << ": [" << xLmode3[k] << "," <<  xUmode3[k] << "]" << endl;
		  else if(mode==4 || mode==5)
		    cout << "Var " << k << ": [" << xLmode4[k] << "," <<  xUmode4[k] << "]" << endl;
		  else if(mode==6)
		    cout << "Var " << k << ": [" << xLmode6[k] << "," <<  xUmode6[k] << "]" << endl;
		  else{}
		}	   
	      }	

	      clock->Start();

	      double p = 0.;

	      for(int hyp = 0 ; hyp<2; hyp++){

		if(SoB==0 && hyp!=hypo) continue;
		if(SoB==1 && hyp==0 && skip) continue;
		if( hyp==1 && false/*&& !(mode==0 || mode==1)*/){
		  cout << "Unsupported at the moment..." << endl;
		  continue;
		}

		if(printP4) cout << "Testing hypothesis " << (hyp==0 ? "S" : "B") << endl; 
		meIntegrator->setHypo(hyp);

		int ntries = 0;
		int intPoints = vegasPoints;
		while( ntries < MAX_REEVAL_TRIES){

		  ROOT::Math::Functor toIntegrate(meIntegrator, &MEIntegratorNew::Eval, par+hyp);
		  ROOT::Math::GSLMCIntegrator ig2( ROOT::Math::IntegrationMultiDim::kVEGAS , 1.e-12, 1.e-5, intPoints );
		  ig2.SetFunction(toIntegrate);
		  meIntegrator->SetPar(par+hyp);
		  
		  if(mode==0) 
		    p = (hyp==0 ? ig2.Integral(xLmode0, xUmode0) : ig2.Integral(xLmode0_b, xUmode0_b));
		  else if(mode==1)
		    p = (hyp==0 ? ig2.Integral(xLmode1, xUmode1) : ig2.Integral(xLmode1_b, xUmode1_b));
		  else if(mode==2)
		    p = (hyp==0 ? ig2.Integral(xLmode2, xUmode2) : ig2.Integral(xLmode2_b, xUmode2_b));
		  else if(mode==3)
		    p = (hyp==0 ? ig2.Integral(xLmode3, xUmode3) : ig2.Integral(xLmode3_b, xUmode3_b)) ;
		  else if(mode==4 || mode==5)
		    p = (hyp==0 ? ig2.Integral(xLmode4, xUmode4) : ig2.Integral(xLmode4_b, xUmode4_b));
		  else if(mode==6)      
		    p = (hyp==0 ? ig2.Integral(xLmode6, xUmode6) : ig2.Integral(xLmode6_b, xUmode6_b));
		  else{ }
		  
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
		    cout << "p = " << p << "; VEGAS has chi2=" << chi2 << " > " << maxChi2_ << ": trying again with " << intPoints << " points..." << endl;
		  }
		  else ntries = MAX_REEVAL_TRIES+1;

		}

		clock->Stop();

		if( TMath::IsNaN(p) ) p = 0.;

		if(hyp==0){
		  permutProbInt[pos] += p;
		}
		else{
		  permutProbInt_alt[pos] += p;
		}

		if(hyp==0) histos[ pos ]->SetBinContent( histos[ pos ]->FindBin(mH[m]), p);
		else       histos_alt[ pos ]->SetBinContent( histos_alt[ pos ]->FindBin(mH[m]), p);

		if( pos==0 && p>pTotAtSgn && hyp==0){
		  pTotAtSgn = p;
		}	       
		if( pos==0 && p>pTotAtSgn_alt && hyp==1){
		  pTotAtSgn_alt = p;
		}

		if(hyp==0){
		  evalCpuGA_ += clock->CpuTime();
		  evalReaGA_ += clock->RealTime();
		}
		else{
		  evalCpuGA_alt_ += clock->CpuTime();
		  evalReaGA_alt_ += clock->RealTime();
		}
	      
		clock->Reset();

		if( mH[m]<met+0.5 && mH[m]>met-0.5){
		  if(hyp==0) 
		    probAtSgnGA_ += p;
		  else
		    probAtSgnGA_alt_ += p;
		}
	      
		if( hyp==0 ){
		  float old = hMassProb_gen->GetBinContent( hMassProb_gen->FindBin(mH[m]) );
		  hMassProb_gen->SetBinContent( hMassProb_gen->FindBin(mH[m]), old+p );
		}
		else{
		  float old = hMassProb_alt_gen->GetBinContent( hMassProb_alt_gen->FindBin(mH[m]) );
		  hMassProb_alt_gen->SetBinContent( hMassProb_alt_gen->FindBin(mH[m]), old+p );
		}

	      }

	    } // permutations

	    if(printP4){
	      if( SoB==0 ) cout << "M = " << mH[m]  << " => p = " << (hypo==0 ? hMassProb_gen->GetBinContent( hMassProb_gen->FindBin(mH[m]) ) : hMassProb_alt_gen->GetBinContent( hMassProb_alt_gen->FindBin(mH[m]) ) ) << endl; 
	      else{
		cout << "S: M = " << mH[m]  << " => p = " << hMassProb_gen->GetBinContent( hMassProb_gen->FindBin(mH[m]) ) << endl; 
		cout << "B: M = " << mH[m]  << " => p = " << hMassProb_alt_gen->GetBinContent( hMassProb_alt_gen->FindBin(mH[m]) ) << endl; 
	      }
	    }
	    
	  }
	  
	  double maxIntProb     = 0.;
	  unsigned int permMax  = 0;
	  int posMax = 0;
	  for(unsigned int it = 0; it<(unsigned int)nPermutations; it++){
	    if(printP4) cout << "Permut #" << it << " has int p = " << permutProbInt[it] << endl;
	    if(mode==0 || mode==1 || mode==5 || mode==6){
	      if( permutProbInt[it] > TMath::Max(permutProbInt[0],permutProbInt[1]) ) posMax++;
	    }
	    if(mode==2 || mode==3){
	      if( permutProbInt[it] > permutProbInt[0] ) posMax++;
	    }
	    if( permutProbInt[it]>maxIntProb ){
	      maxIntProb = permutProbInt[it];
	      permMax = it;
	    }	 
	  }

	  double maxIntProbMod  = 0.;
	  unsigned int permMax3 = 0;
	  if(mode==0 || mode==1 || mode==6 || mode==5){
	    for(unsigned int it = 0; it<(unsigned int)(nPermutations-1); it++){
	      if(it%2==0){
		if(printP4) cout << "Permut #" << it << "+#" << it+1 <<  "  has int p = " << permutProbInt[it] << " => " << (permutProbInt[it]+permutProbInt[it+1]) << endl;		
		if( (permutProbInt[it]+permutProbInt[it+1])>maxIntProbMod ){
		  maxIntProbMod = (permutProbInt[it]+permutProbInt[it+1]);
		  permMax3 = it;
		}
	      }
	    }
	  }
	  
	  pair<double,double> bestInt    = (getMaxValue( histos[permMax] ));
	  massIntGA_ = bestInt.first;
	  probIntGA_ = bestInt.second;
	  posIntGA_  = posMax;

	  if(gDirectory->FindObject("hSumIntMod")!=0){
	    gDirectory->Remove(gDirectory->FindObject("hSumIntMod"));
	  }
	  TH1F* hSumIntMod = (TH1F*)histos[permMax3]->Clone("hSumIntMod");
	  hSumIntMod->Add(histos[permMax3+1]);
	  pair<double,double> bestIntMod = (getMaxValue(  hSumIntMod ));
	  massIntModGA_ = bestIntMod.first;
	  probIntModGA_ = bestIntMod.second;


	  hBestMassProbInt_gen->Fill( massIntGA_ ); 

	  if(mode==0 || mode==1 || mode==6 || mode==4 || mode==5){
	    matchByIntGA_    = (permMax==0 || permMax==1) ? 1 : 0;
	    matchByIntModGA_ = (permMax3==0) ? 1 : 0;
	  }
	  else if(mode==2 || mode==3)
	    matchByIntGA_ = (permMax==0) ? 1 : 0;
	  else{}

	  double maxMaxProb = 0.;
	  unsigned int permMax2 = 0;
	  posMax = 0;

	  cout << "# of permutations considered: = " << histos.size() << endl;
	  sMassProbInt_gen->Modified();
	  for(unsigned int it = 0; it < histos.size(); it++){
	    if(mode==0 || mode==1 || mode==5 || mode==6){
	      if( histos[it]->GetMaximum() > TMath::Max(histos[0]->GetMaximum(),histos[1]->GetMaximum()) ) posMax++;
	    }
	    if(mode==2 || mode==3){
	      if( histos[it]->GetMaximum() > histos[0]->GetMaximum() ) posMax++;
	    }
	    if( histos[it]->GetMaximum()>maxMaxProb ){
	      maxMaxProb = histos[it]->GetMaximum();
	      permMax2 = it;
	    }
	    if(it>0){
	      histos[it]->SetFillStyle(3003);
	      histos[it]->SetFillColor(it+30);
	      histos[it]->SetLineColor(it+30);
	      sMassProbInt_gen->Add( histos[it] );
	      leg->AddEntry(histos[it], Form("%d (%.0f GeV)", permutations[it] , permutMassReco[it]), "F");
	    }
	  }
	  if(histos.size()>0){
	    histos[0]->SetFillStyle(3004);
	    histos[0]->SetFillColor(2);
	    histos[0]->SetLineColor(2);
	    sMassProbInt_gen->Add( histos[0] );
	    leg->AddEntry(histos[0], Form("%d (%.0f GeV)", permutations[0] , permutMassReco[0]), "F");	      
	  }

	  pair<double,double> bestMax = (getMaxValue( histos[permMax2] ));
	  massMaxGA_ =  bestMax.first ;
	  probMaxGA_ =  bestMax.second;
	  posMaxGA_  = posMax;

	  double maxMaxModProb = 0.;
	  unsigned int permMax4 = 0;

	  if(mode==0 || mode==1 || mode==6 || mode==5){
	    for(unsigned int it = 0; it < (unsigned int)(histos.size()-1); it++){
	      if(it%2==0){
		if(gDirectory->FindObject("hSum")!=0){
		  gDirectory->Remove(gDirectory->FindObject("hSum"));
		}
		TH1F* hSum = (TH1F*)histos[it]->Clone("hSum");
		hSum->Add(histos[it+1] );
		if( hSum->GetMaximum() > maxMaxModProb ){
		  maxMaxModProb =  hSum->GetMaximum();
		  permMax4 = it;
		}
	      }
	    }
	  }

	  if(gDirectory->FindObject("hSumMaxMod")!=0){
	    gDirectory->Remove(gDirectory->FindObject("hSumMaxMod"));
	  }
	  TH1F* hSumMaxMod = (TH1F*)histos[permMax4]->Clone("hSumMaxMod");
	  hSumMaxMod->Add(histos[permMax4+1]);
	  pair<double,double> bestMaxMod = (getMaxValue( hSumMaxMod ));
	  massMaxModGA_ =  bestMaxMod.first ;
	  probMaxModGA_ =  bestMaxMod.second;

	  if(mode==0 || mode==1 || mode==6 || mode==4 || mode==5){
	    matchByMaxGA_    = (permMax2==0 || permMax2==1) ? 1 : 0;
	    matchByMaxModGA_ = (permMax4==0) ? 1 : 0;
	  }
	  else if(mode==2 || mode==3)
	    matchByMaxGA_ = (permMax2==0) ? 1 : 0;
	  else{}

	  hBestMassProbMax_gen->Fill( massMaxGA_ ); 
	  histos.clear();
	  

	  evalCpuGA_ /= nMassPoints;//combScanned;
	  evalReaGA_ /= nMassPoints;//combScanned;
	  probAtGoodGA_ = pTotAtSgn;
	
	  pair<double,double> bestMass = getMaxValue(hMassProb_gen);
	  hBestMassProb_gen->Fill( bestMass.first );
	  massGA_ = bestMass.first;
	  probGA_ = bestMass.second;

	  TFile* fout_tmp = TFile::Open(outFileName.c_str(),"UPDATE");	  
	  if(gDirectory->GetKey(Form("Event_%d", counter)) == 0) {
	    TDirectory *dir = fout_tmp->mkdir(Form("Event_%d", counter));
	    dir->cd();
	    hMassProb_gen->Write("", TObject::kOverwrite);
	    c1->cd();
	    sMassProbInt_gen->Draw();
	    leg->Draw();
	    //sMassProbInt_gen->Write("", TObject::kOverwrite);
	    c1->Write("", TObject::kOverwrite);
	  }
	  else{
	    fout_tmp->cd(Form("Event_%d", counter));
	    hMassProb_gen->Write("", TObject::kOverwrite);
	    sMassProbInt_gen->Draw();
	    leg->Draw();
	    //sMassProbInt_gen->Write("", TObject::kOverwrite);
	    c1->Write("", TObject::kOverwrite);
	  }
	  fout_tmp->Close();
	  delete sMassProbInt_gen; delete leg; delete c1;

	}

      }


      if( doSmear ){

	if(mode == 5)
	  meIntegrator->setIntType( MEIntegratorNew::SL3b );
	meIntegrator->initVersors(1);

	bool passes = false;
	int ntries  = 0;
	while( ntries<N_TRIES_MAX && !passes){
	  meIntegrator->setJets(&jets);
	  passes = meIntegrator->smearByTF(30., printP4);
	  ntries++;
	  if(passes && printP4) 
	    cout << "Smearing succedes at the " << ntries << "th attempt" << endl;
	  else if(printP4 && !passes) 
	    cout << "Smearing fails... let's try again!" << endl;
	  else{}
	}
	
	if(mode == 5){
	  TLorentzVector bLep =  meIntegrator->jetAt(2);
	  TLorentzVector bHad =  meIntegrator->jetAt(5);
	  TLorentzVector hig1 =  meIntegrator->jetAt(6);
	  TLorentzVector hig2 =  meIntegrator->jetAt(7);
	  if( (bLep.Pt()<30 || TMath::Abs(bLep.Eta())>2.5) && (bHad.Pt()>30 && TMath::Abs(bHad.Eta())<2.5) && (hig1.Pt()>30 && TMath::Abs(hig1.Eta())<2.5)  && (hig2.Pt()>30 && TMath::Abs(hig2.Eta())<2.5)){
	    type_ = 3;
	  }
	  else if( (bLep.Pt()>30 && TMath::Abs(bLep.Eta())<2.5) && (bHad.Pt()<30 || TMath::Abs(bHad.Eta())>2.5) && (hig1.Pt()>30 && TMath::Abs(hig1.Eta())<2.5)  && (hig2.Pt()>30 && TMath::Abs(hig2.Eta())<2.5)){
	    type_ = 2;
	  }
	  else if( (bLep.Pt()>30 && TMath::Abs(bLep.Eta())<2.5) && (bHad.Pt()>30 && TMath::Abs(bHad.Eta())<2.5) && (hig1.Pt()>30 && TMath::Abs(hig1.Eta())<2.5)  && (hig2.Pt()<30 || TMath::Abs(hig2.Eta())>2.5)){
	    type_ = 4;
	  }
	  else if((bLep.Pt()>30 && TMath::Abs(bLep.Eta())<2.5) && (bHad.Pt()>30 && TMath::Abs(bHad.Eta())<2.5) && (hig1.Pt()<30 || TMath::Abs(hig1.Eta())>2.5)  && (hig2.Pt()>30 && TMath::Abs(hig2.Eta())<2.5)){
	    type_ = 4;
	  }
	  else{}
	}


	if( passes ){

	  ////////////////////////////////////////////
	  TLorentzVector fv1 =  meIntegrator->jetAt(0);
	  TLorentzVector fv2 =  meIntegrator->jetAt(1);
	  TLorentzVector fv3 =  meIntegrator->jetAt(2);
	  TLorentzVector fv4 =  meIntegrator->jetAt(3);
	  TLorentzVector fv5 =  meIntegrator->jetAt(4);
	  TLorentzVector fv6 =  meIntegrator->jetAt(5);
	  TLorentzVector fv7 =  meIntegrator->jetAt(6);
	  TLorentzVector fv8 =  meIntegrator->jetAt(7);
	  dPx_    = jets[1].Px() - fv2.Px();
	  dRecPx_ = -(fv1+fv2+fv3+fv4+fv5+fv6+fv7+fv8).Px() + (jets[0]+jets[1]+jets[2]+jets[3]+jets[4]+jets[5]+jets[6]+jets[7]).Px(); 
	  ////////////////////////////////////////////
	  
	  
	  if( doMassScan ){

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

	    double xLmode2[6]   = {x0L, x1L, x2L, x3L, x4L, x5L};
	    double xUmode2[6]   = {x0U, x1U, x2U, x3U, x4U, x5U};
	    double xLmode2_b[7] = {x0L, x1L, x2L, x3L, x4L, x5L, x6L};
	    double xUmode2_b[7] = {x0U, x1U, x2U, x3U, x4U, x5U, x6U};
	    
	    double xLmode3[6]   = {x0L, x1L, x2L, x3L, x4L, x5L};
	    double xUmode3[6]   = {x0U, x1U, x2U, x3U, x4U, x5U};
	    double xLmode3_b[7] = {x0L, x1L, x2L, x3L, x4L, x5L, x6L};
	    double xUmode3_b[7] = {x0U, x1U, x2U, x3U, x4U, x5U, x6U};
	    
	    double xLmode4[6]   = {x0L, x1L, x2L, x3L, x4L, x5L};
	    double xUmode4[6]   = {x0U, x1U, x2U, x3U, x4U, x5U};
	    double xLmode4_b[7] = {x0L, x1L, x2L, x3L, x4L, x5L, x6L};
	    double xUmode4_b[7] = {x0U, x1U, x2U, x3U, x4U, x5U, x6U};
	    	    
	    double xLmode6[5]   = {x1L, x2L, x1L, x2L, x5L};
	    double xUmode6[5]   = {x1U, x2U, x1U, x2U, x5U};
	    double xLmode6_b[6] = {x1L, x2L, x1L, x2L, x5L, x6L};
	    double xUmode6_b[6] = {x1U, x2U, x1U, x2U, x5U, x6U};


	    if( printP4 ){
	      cout << "Integration range: " << endl;
	      for(unsigned int k = 0; k < (unsigned int)(par); k++){
		if(mode==0) 
		  cout << "Var " << k << ": [" << (hypo==0 ? xLmode0[k] : xLmode0_b[k]) << "," <<  (hypo==0 ? xUmode0[k] : xUmode0_b[k]) << "]" << endl;
		else if(mode==1)
		  cout << "Var " << k << ": [" << (hypo==0 ? xLmode1[k] : xLmode1_b[k]) << "," <<  (hypo==0 ? xUmode1[k] : xUmode1_b[k]) << "]" << endl;
		else if(mode==2)
		  cout << "Var " << k << ": [" << xLmode2[k] << "," <<  xUmode2[k] << "]" << endl;
		else if(mode==3)
		  cout << "Var " << k << ": [" << xLmode3[k] << "," <<  xUmode3[k] << "]" << endl;
		else if(mode==4)
		  cout << "Var " << k << ": [" << xLmode4[k] << "," <<  xUmode4[k] << "]" << endl;
		else if(mode==6)
		  cout << "Var " << k << ": [" << xLmode6[k] << "," <<  xUmode6[k] << "]" << endl;
		else{}
	      }	   
	    }			   
	    
	    hMass_rec->Reset();
	    hMass_alt_rec->Reset();
	    
	    for(int m = 0; m < nMassPoints ; m++){

	      meIntegrator->setMass( mH[m] );

	      clock->Start();

	      double p = 0.;

	      for(int hyp = 0 ; hyp<2; hyp++){

		if(SoB==0 && hyp!=hypo) continue;
		if( hyp==1 && false/*&& !(mode==0 || mode==1)*/){
		  cout << "Unsupported at the moment..." << endl;
		  continue;
		}
		
		if(printP4) cout << "Testing hypothesis " << (hyp==0 ? "S" : "B") << endl; 
		meIntegrator->setHypo(hyp);
		
		int ntries = 0;
		double err = 0;
		int intPoints = vegasPoints;
		while( ntries < MAX_REEVAL_TRIES){

		  ROOT::Math::Functor toIntegrate(meIntegrator, &MEIntegratorNew::Eval, par+hyp);
		  ROOT::Math::GSLMCIntegrator ig2( ROOT::Math::IntegrationMultiDim::kVEGAS , 1.e-12, 1.e-5, intPoints);
		  ig2.SetFunction(toIntegrate);
		  meIntegrator->SetPar(par+hyp);
		  
		  if(mode==0)
		    p = (hyp==0 ? ig2.Integral(xLmode0, xUmode0) : ig2.Integral(xLmode0_b, xUmode0_b));
		  else if(mode==1)
		    p = (hyp==0 ? ig2.Integral(xLmode1, xUmode1) : ig2.Integral(xLmode1_b, xUmode1_b));
		  else if(mode==2)
		    p = (hyp==0 ? ig2.Integral(xLmode2, xUmode2) : ig2.Integral(xLmode2_b, xUmode2_b));
		  else if(mode==3)
		    p = (hyp==0 ? ig2.Integral(xLmode3, xUmode3) : ig2.Integral(xLmode3_b, xUmode3_b)) ;
		  else if(mode==4)
		    p = (hyp==0 ? ig2.Integral(xLmode4, xUmode4) : ig2.Integral(xLmode4_b, xUmode4_b));
		  else if(mode==6)
		    p = (hyp==0 ? ig2.Integral(xLmode6, xUmode6) : ig2.Integral(xLmode6_b, xUmode6_b));
		  else{ }	

		  cout << "Chi2 = " << ig2.ChiSqr() << endl;
		  chi2RG_ =  ig2.ChiSqr();
		  err = ig2.Error() ;

		  if( chi2RG_ > maxChi2_ ){
		    ntries++;
		    intPoints *= 1.5;
		    cout << "p = " << p << "; VEGAS has chi2=" << chi2RG_ << " > " << maxChi2_ << ": trying again with " << intPoints << " points..." << endl;
		  }
		  else ntries = MAX_REEVAL_TRIES+1;
		}
		
		clock->Stop();
		
		if( TMath::IsNaN(p) ) p = 0.;	     
		
		if(hyp==0){
		  evalCpuRG_ += clock->CpuTime();
		  evalReaRG_ += clock->RealTime();
		}
		else{
		  evalCpuRG_alt_ += clock->CpuTime();
		  evalReaRG_alt_ += clock->RealTime();
		}
		
		if(printP4 || true) cout << "Mass " << mH[m] << " => prob  = " << p << " +/- " <<  err << endl;
		if(hyp==0){
		  hMass_rec->Fill(mH[m], p);
		}
		else{
		  hMass_alt_rec->Fill(mH[m], p);
		}
		
		if( mH[m]<met+0.5 && mH[m]>met-0.5){
		  if(hyp==0) 
		    probAtSgnRG_ += p;
		  else 
		    probAtSgnRG_alt_ += p;
		}
		
		p = 0.;
		clock->Reset();
	      }

	    }
	    evalCpuRG_ /= nMassPoints;
	    evalReaRG_ /= nMassPoints;
	    evalCpuRG_alt_ /= nMassPoints;
	    evalReaRG_alt_ /= nMassPoints;
	    
	    pair<double,double> bestMass = getMaxValue(hMass_rec);
	    hBestMass_rec->Fill( bestMass.first );
	    massRG_ = bestMass.first;
	    probRG_ = bestMass.second;

	    TFile* fout_tmp = TFile::Open(outFileName.c_str(),"UPDATE");	  
	    if(gDirectory->GetKey(Form("Event_%d", counter)) == 0) {
	      TDirectory *dir = fout_tmp->mkdir(Form("Event_%d", counter));
	      dir->cd();
	      hMass_rec->Write("", TObject::kOverwrite);
	    }
	    else{
	      fout_tmp->cd(Form("Event_%d", counter));
	      hMass_rec->Write("", TObject::kOverwrite);
	    }
	    fout_tmp->Close();
	  }
	
	  
	  if(doPermutations){
	    
	    hMassProb_rec->Reset();
	    hMassProb_alt_rec->Reset();

	    double pTotAtSgn = 0.;
	    double pTotAtSgn_alt = 0.;

	    int combScanned = 0;
	    
	 
	    double permutMassReco[nPermutations];
	    for(int k = 0; k < nPermutations; k++) permutMassReco[k] = 0.;
	    
	    double permutProbInt [nPermutations];
	    for(int k = 0; k < nPermutations; k++) permutProbInt[k] = 0.;

	    double permutProbInt_alt [nPermutations];
	    for(int k = 0; k < nPermutations; k++) permutProbInt_alt[k] = 0.;

	    vector<TH1F*> histos;
	    for(unsigned int it = 0; it<(unsigned int)nPermutations; it++){
	      if(gDirectory->FindObject(Form("hMass_%d",it))!=0){
		gDirectory->Remove(gDirectory->FindObject(Form("hMass_%d",it)));
	      }
	      TH1F* h = (TH1F*)hMass_rec->Clone(Form("hMass_%d",it));
	      h->Reset();
	      histos.push_back( h );
	    }
	    vector<TH1F*> histos_alt;
	    for(unsigned int it = 0; it<(unsigned int)nPermutations; it++){
	      if(gDirectory->FindObject(Form("hMass_alt_%d",it))!=0){
		gDirectory->Remove(gDirectory->FindObject(Form("hMass_alt_%d",it)));
	      }
	      TH1F* h = (TH1F*)hMass_gen->Clone(Form("hMass_alt_%d",it));
	      h->Reset();
	      histos_alt.push_back( h );
	    }



	    if(gDirectory->FindObject("stack_m_all_rec")!=0){
	      gDirectory->Remove(gDirectory->FindObject("stack_m_all_rec"));
	    }
	    THStack*  sMassProbInt_rec  = new THStack("stack_m_all_rec", "Per permutation mass, REC");
	    TLegend* leg = new TLegend(0.63,0.48,0.85,0.85,NULL,"brNDC");
	    leg->SetFillStyle(0);
	    leg->SetBorderSize(0);
	    leg->SetFillColor(10);
	    leg->SetTextSize(0.03);
	    TCanvas *c1 = new TCanvas("c_rec","canvas",10,30,650,600);


	    for(int m = 0; m < nMassPoints ; m++){
	      meIntegrator->setMass( mH[m] );
	      
	      double maxP_s = 0;
	      double maxP_b = 0;

	      for(unsigned int pos = 0; pos < (unsigned int)nPermutations; pos++){
		meIntegrator->initVersors( permutations[pos] );
		
		int mode5HiggsLost = 0;
		if( mode==5 ){ // figure out which b is missing, then set the correspinding int mode
		  TLorentzVector bLep =  meIntegrator->jetAt(2);
		  TLorentzVector bHad =  meIntegrator->jetAt(5);
		  TLorentzVector hig1 =  meIntegrator->jetAt(6);
		  TLorentzVector hig2 =  meIntegrator->jetAt(7);
		  if( (bLep.Pt()<30 || TMath::Abs(bLep.Eta())>2.5) && (bHad.Pt()>30 && TMath::Abs(bHad.Eta())<2.5) && (hig1.Pt()>30 && TMath::Abs(hig1.Eta())<2.5)  && (hig2.Pt()>30 && TMath::Abs(hig2.Eta())<2.5)){
		    if(printP4) cout << "Permutation " << permutations[pos] << " has assigned a jet with (pt,eta)=(" << bLep.Pt() << "," << bLep.Eta() << ") to the bLep: setting mode SLNoBLep" << endl;
		    meIntegrator->setIntType( MEIntegratorNew::SLNoBLep );
		  }
		  else if( (bLep.Pt()>30 && TMath::Abs(bLep.Eta())<2.5) && (bHad.Pt()<30 || TMath::Abs(bHad.Eta())>2.5) && (hig1.Pt()>30 && TMath::Abs(hig1.Eta())<2.5)  && (hig2.Pt()>30 && TMath::Abs(hig2.Eta())<2.5)){
		    if(printP4)  cout << "Permutation " << permutations[pos] << " has assigned a jet with (pt,eta)=(" << bHad.Pt() << "," << bHad.Eta() << ") to the bHad: setting mode SLNoBHad" << endl;
		    meIntegrator->setIntType( MEIntegratorNew::SLNoBHad );
		  }
		  else if( (bLep.Pt()>30 && TMath::Abs(bLep.Eta())<2.5) && (bHad.Pt()>30 && TMath::Abs(bHad.Eta())<2.5) && (hig1.Pt()>30 && TMath::Abs(hig1.Eta())<2.5)  && (hig2.Pt()<30 || TMath::Abs(hig2.Eta())>2.5)){
		    if(printP4)  cout << "Permutation " << permutations[pos] << " has assigned a jet with (pt,eta)=(" << hig2.Pt() << "," << hig2.Eta() << ") to  Higgs2: setting mode SLNoHiggs" << endl;
		    meIntegrator->setIntType( MEIntegratorNew::SLNoHiggs );
		    mode5HiggsLost = 1;		  
		  }
		  else if((bLep.Pt()>30 && TMath::Abs(bLep.Eta())<2.5) && (bHad.Pt()>30 && TMath::Abs(bHad.Eta())<2.5) && (hig1.Pt()<30 || TMath::Abs(hig1.Eta())>2.5)  && (hig2.Pt()>30 && TMath::Abs(hig2.Eta())<2.5)){
		    meIntegrator->initVersors( permutations_SL3b_equivalent[pos] );
		    hig1 =  meIntegrator->jetAt(6);
		    hig2 =  meIntegrator->jetAt(7);
		    if(printP4)  cout << "Permutation " << permutations_SL3b_equivalent[pos] << " has assigned a jet with (pt,eta)=(" << hig2.Pt() << "," << hig2.Eta() << ") to  Higgs2: setting mode SLNoHiggs" << endl;
		    meIntegrator->setIntType( MEIntegratorNew::SLNoHiggs );
		    mode5HiggsLost = 1;
		  }
		  else{
		    cout << "Mode " << mode << ": inconsistency between mode and jet list" << endl;
		    continue;
		  }
		} 
		
		double mass, massLow, massHigh;
		bool skip = !(meIntegrator->compatibilityCheck(0.95, /*printP4*/ 0, mass, massLow, massHigh )) ;
		permutMassReco[pos] = mass;
		if( SoB==0 && skip && !(mode==4 || (mode==5 && mode5HiggsLost)) ){
		  continue;
		}

		
		combScanned++;
		
		pair<double, double> range_x0 = (meIntegrator->getW1JetEnergyCI(0.95));
		pair<double, double> range_x1 =  make_pair(-1,1);
		pair<double, double> range_x2 =  make_pair(-PI,PI);
		pair<double, double> range_x3 =  make_pair(-1,1);
		pair<double, double> range_x4 = useMET ? (meIntegrator->getNuPhiCI(0.95)) :  make_pair(-PI,PI); 
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

		double xLmode2[6]   = {x0L, x1L, x2L, x3L, x4L, x5L};
		double xUmode2[6]   = {x0U, x1U, x2U, x3U, x4U, x5U};
		double xLmode2_b[7] = {x0L, x1L, x2L, x3L, x4L, x5L, x6L};
		double xUmode2_b[7] = {x0U, x1U, x2U, x3U, x4U, x5U, x6U};
		
		double xLmode3[6]   = {x0L, x1L, x2L, x3L, x4L, x5L};
		double xUmode3[6]   = {x0U, x1U, x2U, x3U, x4U, x5U};
		double xLmode3_b[7] = {x0L, x1L, x2L, x3L, x4L, x5L, x6L};
		double xUmode3_b[7] = {x0U, x1U, x2U, x3U, x4U, x5U, x6U};
		
		double xLmode4[6]   = {x0L, x1L, x2L, x3L, x4L, x5L};
		double xUmode4[6]   = {x0U, x1U, x2U, x3U, x4U, x5U};
		double xLmode4_b[7] = {x0L, x1L, x2L, x3L, x4L, x5L, x6L};
		double xUmode4_b[7] = {x0U, x1U, x2U, x3U, x4U, x5U, x6U};
		
		double xLmode6[5]   = {x1L, x2L, x1L, x2L, x5L};
		double xUmode6[5]   = {x1U, x2U, x1U, x2U, x5U};
		double xLmode6_b[6] = {x1L, x2L, x1L, x2L, x5L, x6L};
		double xUmode6_b[6] = {x1U, x2U, x1U, x2U, x5U, x6U};


		if( printP4 && false){
		  cout << "Integration range: " << endl;
		  for(unsigned int k = 0; k < (unsigned int)(par); k++){
		    if(mode==0) 
		      cout << "Var " << k << ": [" << (hypo==0 ? xLmode0[k] : xLmode0_b[k]) << "," <<  (hypo==0 ? xUmode0[k] : xUmode0_b[k]) << "]" << endl;
		    else if(mode==1)
		      cout << "Var " << k << ": [" << (hypo==0 ? xLmode1[k] : xLmode1_b[k]) << "," <<  (hypo==0 ? xUmode1[k] : xUmode1_b[k]) << "]" << endl;
		    else if(mode==2)
		      cout << "Var " << k << ": [" << xLmode2[k] << "," <<  xUmode2[k] << "]" << endl;
		    else if(mode==3)
		      cout << "Var " << k << ": [" << xLmode3[k] << "," <<  xUmode3[k] << "]" << endl;
		    else if(mode==4 || mode==5)
		      cout << "Var " << k << ": [" << xLmode4[k] << "," <<  xUmode4[k] << "]" << endl;
		    else if(mode==6)
		      cout << "Var " << k << ": [" << xLmode6[k] << "," <<  xUmode6[k] << "]" << endl;
		    else{}
		  }	   
		}	
			
		clock->Start();
		
		double p = 0.;

		for(int hyp = 0 ; hyp<2; hyp++){

		  if(SoB==0 && hyp!=hypo) continue;
		  if(SoB==1 && hyp==0 && skip) continue;
		  if( hyp==1 && false/*&& !(mode==0 || mode==1)*/){
		    cout << "Unsupported at the moment..." << endl;
		    continue;
		  }

		  if(printP4) cout << "Testing hypothesis " << (hyp==0 ? "S" : "B") << endl; 
		  meIntegrator->setHypo(hyp);

		  int ntries = 0;
		  int intPoints = vegasPoints;
		  while( ntries < MAX_REEVAL_TRIES){

		    ROOT::Math::Functor toIntegrate(meIntegrator, &MEIntegratorNew::Eval, par+hyp);
		    ROOT::Math::GSLMCIntegrator ig2( ROOT::Math::IntegrationMultiDim::kVEGAS , 1.e-12, 1.e-5, intPoints);
		    ig2.SetFunction(toIntegrate);
		    meIntegrator->SetPar(par+hyp);		   
		    
		    if(mode==0)
		      p = (hyp==0 ? ig2.Integral(xLmode0, xUmode0) : ig2.Integral(xLmode0_b, xUmode0_b));
		    else if(mode==1)
		      p = (hyp==0 ? ig2.Integral(xLmode1, xUmode1) : ig2.Integral(xLmode1_b, xUmode1_b));
		    else if(mode==2)
		      p = (hyp==0 ? ig2.Integral(xLmode2, xUmode2) : ig2.Integral(xLmode2_b, xUmode2_b));
		    else if(mode==3)
		      p = (hyp==0 ? ig2.Integral(xLmode3, xUmode3) : ig2.Integral(xLmode3_b, xUmode3_b)) ;
		    else if(mode==4 || mode==5)
		      p = (hyp==0 ? ig2.Integral(xLmode4, xUmode4) : ig2.Integral(xLmode4_b, xUmode4_b));
		    else if(mode==6)
		      p = (hyp==0 ? ig2.Integral(xLmode6, xUmode6) : ig2.Integral(xLmode6_b, xUmode6_b));
		    else{ }
		    
		    double chi2 =  ig2.ChiSqr();
		    if( hyp==0 ){
		      if( p>maxP_s ) maxP_s = p;		  
		      else if( p<0.1*maxP_s){
			ntries = MAX_REEVAL_TRIES+1;
			continue;
		      }
		    }
		    else{
		      if( p>maxP_b ) maxP_b = p;		  
		      else if( p<0.1*maxP_b){
			ntries = MAX_REEVAL_TRIES+1;	
			continue;
		      }
		    }


		    if( chi2 > maxChi2_ ){
		      ntries++;
		      intPoints *= 1.5;
		      cout << (hyp==0 ? "S" : "B")  << ": p = " << p << "; VEGAS has chi2=" << chi2 << " > " << maxChi2_ << ": trying again with " << intPoints << " points..." << endl;
		    }
		    else ntries = MAX_REEVAL_TRIES+1;
		    
		  }

		  clock->Stop();
		  
		  if( TMath::IsNaN(p) ) p = 0.;

		  if(hyp==0){
		    permutProbInt[pos] += p;
		  }
		  else{
		    permutProbInt_alt[pos] += p;
		  }

		  if(hyp==0) histos[ pos ]->SetBinContent( histos[ pos ]->FindBin(mH[m]), p);
		  else       histos_alt[ pos ]->SetBinContent( histos_alt[ pos ]->FindBin(mH[m]), p);

		  if( pos==0 && p>pTotAtSgn && hyp==0){
		    pTotAtSgn = p;
		  }	       
		  if( pos==0 && p>pTotAtSgn_alt && hyp==1){
		    pTotAtSgn_alt = p;
		  }
		  
		  if(hyp==0){
		    evalCpuRA_ += clock->CpuTime();
		    evalReaRA_ += clock->RealTime();
		  }
		  else{
		    evalCpuRA_alt_ += clock->CpuTime();
		    evalReaRA_alt_ += clock->RealTime();
		  }
		
		
		  clock->Reset();


		  if( mH[m]<met+0.5 && mH[m]>met-0.5){
		    if(hyp==0) 
		      probAtSgnRA_ += p;
		    else
		      probAtSgnRA_alt_ += p;
		  }

		  if( hyp==0 ){
		    float old = hMassProb_rec->GetBinContent( hMassProb_rec->FindBin(mH[m]) );
		    hMassProb_rec->SetBinContent( hMassProb_rec->FindBin(mH[m]), old+p );
		  }
		  else{
		    float old = hMassProb_alt_rec->GetBinContent( hMassProb_alt_rec->FindBin(mH[m]) );
		    hMassProb_alt_rec->SetBinContent( hMassProb_alt_rec->FindBin(mH[m]), old+p );
		  }


		}
		
	      } // permutations
	    
	      if(printP4){
		if( SoB==0 ) cout << "M = " << mH[m]  << " => p = " << (hypo==0 ? hMassProb_rec->GetBinContent( hMassProb_rec->FindBin(mH[m]) ) : hMassProb_alt_rec->GetBinContent( hMassProb_alt_rec->FindBin(mH[m]) ) ) << endl; 
		else{
		  cout << "S: M = " << mH[m]  << " => p = " << hMassProb_rec->GetBinContent( hMassProb_rec->FindBin(mH[m]) ) << endl; 
		  cout << "B: M = " << mH[m]  << " => p = " << hMassProb_alt_rec->GetBinContent( hMassProb_alt_rec->FindBin(mH[m]) ) << endl; 
		}
	      }
	      
	    }

	    double maxIntProb = 0.;
	    unsigned int permMax = 0;
	    int posMax = 0;
	    for(unsigned int it = 0; it<(unsigned int)nPermutations; it++){
	      if(printP4) cout << "Permut #" << it << " has int p = " << permutProbInt[it] << endl;
	      if(mode==0 || mode==1 || mode==5 || mode==6){
		if( permutProbInt[it] > TMath::Max(permutProbInt[0],permutProbInt[1]) ) posMax++;
	      }
	      if(mode==2 || mode==3){
		if( permutProbInt[it] > permutProbInt[0] ) posMax++;
	      }
	      if( permutProbInt[it]>maxIntProb ){
		maxIntProb = permutProbInt[it];
		permMax = it;
	      }
	    }

	    double maxIntProbMod  = 0.;
	    unsigned int permMax3 = 0;
	    
	    if(mode==0 || mode==1 || mode==6 || mode==5){
	      for(unsigned int it = 0; it<(unsigned int)(nPermutations-1); it++){
		if(it%2==0){
		if(printP4) cout << "Permut #" << it << "+#" << it+1 <<  "  has int p = " << permutProbInt[it] << " => " << (permutProbInt[it]+permutProbInt[it+1]) << endl;		
		  if( (permutProbInt[it]+permutProbInt[it+1])>maxIntProbMod ){
		    maxIntProbMod = (permutProbInt[it]+permutProbInt[it+1]);
		    permMax3 = it;
		  }
		}
	      }
	    }

	    pair<double,double> bestInt = (getMaxValue( histos[permMax] ));
	    massIntRA_ = bestInt.first;
	    probIntRA_ = bestInt.second;
	    posIntRA_  = posMax;

	    if(gDirectory->FindObject("hSumIntMod")!=0){
	      gDirectory->Remove(gDirectory->FindObject("hSumIntMod"));
	    }
	    TH1F* hSumIntMod = (TH1F*)histos[permMax3]->Clone("hSumIntMod");
	    hSumIntMod->Add(histos[permMax3+1]);	  
	    pair<double,double> bestIntMod = (getMaxValue( hSumIntMod ));
	    massIntModRA_ = bestIntMod.first;
	    probIntModRA_ = bestIntMod.second;

	    hBestMassProbInt_rec->Fill( massIntRA_ ); 

	    if(mode==0 || mode==1 || mode==6 || mode==4 || mode==5){
	      matchByIntRA_ = (permMax==0 || permMax==1) ? 1 : 0;
	      matchByIntModRA_ = (permMax3==0) ? 1 : 0;
	    }
	    else if(mode==2 || mode==3)
	      matchByIntRA_ = (permMax==0) ? 1 : 0;
	    else{}

	    double maxMaxProb = 0.;
	    unsigned int permMax2 = 0;
	    posMax = 0;

	    //cout << "Histos size = " << histos.size() << endl;
	    sMassProbInt_rec->Modified();
	    for(unsigned int it = 0; it < histos.size(); it++){
	      if(mode==0 || mode==1 || mode==5 || mode==6){
		if( histos[it]->GetMaximum() > TMath::Max(histos[0]->GetMaximum(),histos[1]->GetMaximum()) ) posMax++;
	      }
	      if(mode==2 || mode==3){
		if( histos[it]->GetMaximum() > histos[0]->GetMaximum() ) posMax++;
	      }
	      if( histos[it]->GetMaximum()>maxMaxProb ){
		maxMaxProb = histos[it]->GetMaximum();
		permMax2 = it;
	      }
	      if(it>0){
		histos[it]->SetFillStyle(3003);
		histos[it]->SetFillColor(it+30);
		histos[it]->SetLineColor(it+30);
		sMassProbInt_rec->Add( histos[it] );
		leg->AddEntry(histos[it], Form("%d (%.0f GeV)", permutations[it] , permutMassReco[it]), "F");
	      }
	    }
	    if(histos.size()>0){
	      histos[0]->SetFillStyle(3004);
	      histos[0]->SetFillColor(2);
	      histos[0]->SetLineColor(2);
	      sMassProbInt_rec->Add( histos[0] );
	      leg->AddEntry(histos[0], Form("%d (%.0f GeV)", permutations[0] , permutMassReco[0]), "F");	      
	    }

	    pair<double,double> bestMax = (getMaxValue( histos[permMax2] ));
	    massMaxRA_ =  bestMax.first ;
	    probMaxRA_ =  bestMax.second;
	    posMaxRA_  = posMax;
	  
	    double maxMaxModProb = 0.;
	    unsigned int permMax4 = 0;

	    if(mode==0 || mode==2 || mode==6 || mode==5){
	      for(unsigned int it = 0; it < (unsigned int)(histos.size()-1); it++){
		if(it%2==0){
		  if(gDirectory->FindObject("hSum")!=0){
		    gDirectory->Remove(gDirectory->FindObject("hSum"));
		  }
		  TH1F* hSum = (TH1F*)histos[it]->Clone("hSum");
		  hSum->Add(histos[it+1] );
		  if( hSum->GetMaximum() > maxMaxModProb ){
		    maxMaxModProb =  hSum->GetMaximum();
		    permMax4 = it;
		  }
		}
	      }
	    }

	    if(gDirectory->FindObject("hSumMaxMod")!=0){
	      gDirectory->Remove(gDirectory->FindObject("hSumMaxMod"));
	    }
	    TH1F* hSumMaxMod = (TH1F*)histos[permMax4]->Clone("hSumMaxMod");
	    hSumMaxMod->Add(histos[permMax4+1]);
	    pair<double,double> bestMaxMod = (getMaxValue( hSumMaxMod ));
	    massMaxModRA_ =  bestMaxMod.first ;
	    probMaxModRA_ =  bestMaxMod.second;

	    
	    if(mode==0 || mode==1 || mode==6 || mode==4 || mode==5){
	      matchByMaxRA_    = (permMax2==0 || permMax2==1) ? 1 : 0;
	      matchByMaxModRA_ = (permMax4==0) ? 1 : 0;
	    }
	    else if(mode==2 || mode==3)
	      matchByMaxRA_ = (permMax2==0) ? 1 : 0;
	    else{}

	    hBestMassProbMax_rec->Fill( massMaxRA_ ); 
	    histos.clear();
	    
	    evalCpuRA_ /= nMassPoints;//combScanned;
	    evalReaRA_ /= nMassPoints;//combScanned;
	    probAtGoodRA_ = pTotAtSgn;
	    
	    pair<double,double> bestMass = getMaxValue(hMassProb_rec);
	    hBestMassProb_rec->Fill( bestMass.first );
	    massRA_ = bestMass.first;
	    probRA_ = bestMass.second;

	    TFile* fout_tmp = TFile::Open(outFileName.c_str(),"UPDATE");	  
	    if(gDirectory->GetKey(Form("Event_%d", counter)) == 0) {
	      TDirectory *dir = fout_tmp->mkdir(Form("Event_%d", counter));
	      dir->cd();
	      hMassProb_rec->Write("", TObject::kOverwrite);
	      c1->cd();
	      sMassProbInt_rec->Draw();
	      leg->Draw();
	      //sMassProbInt_rec->Write("", TObject::kOverwrite);
	      c1->Write("", TObject::kOverwrite);
	    }
	    else{
	      fout_tmp->cd(Form("Event_%d", counter));
	      hMassProb_rec->Write("", TObject::kOverwrite);
	      sMassProbInt_rec->Draw();
	      leg->Draw();
	      //sMassProbInt_rec->Write("", TObject::kOverwrite);
	      c1->Write("", TObject::kOverwrite);
	    }
	    fout_tmp->Close();
	    delete sMassProbInt_rec; delete leg; delete c1;

	  }
	}
	else{
	  cout << "Smeared event fails acceptance cuts... continue" << endl;
	}

      }

      tree->Fill();      
    } // nentries
    
  } // samples

  TFile* fout_tmp = TFile::Open(outFileName.c_str(),"UPDATE");
  fout_tmp->cd();
  hBestMass_gen->Write       ("m_best_good_gen",        TObject::kOverwrite);
  hBestMassProb_gen->Write   ("m_best_all_gen",         TObject::kOverwrite);
  hBestMassProbInt_gen->Write("m_bestInt_all_gen",      TObject::kOverwrite);
  hBestMassProbMax_gen->Write("m_bestMax_all_gen",      TObject::kOverwrite);
  hBestMass_rec->Write       ("m_best_good_rec",        TObject::kOverwrite);
  hBestMassProb_rec->Write   ("m_best_all_rec",         TObject::kOverwrite);
  hBestMassProbInt_rec->Write("m_bestInt_all_rec",      TObject::kOverwrite);
  hBestMassProbMax_rec->Write("m_bestMax_all_rec",      TObject::kOverwrite);

  tree->Write("",TObject::kOverwrite );
  fout_tmp->Close();

  cout << "Delete meIntegrator..." << endl;
  delete meIntegrator;
  delete ran; delete clock;
  cout << "Finished!!!" << endl;
  
  return 0;
}

