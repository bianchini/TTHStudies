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
#define VERBOSE  false

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


//namespace LHAPDF {
//    void   initPDFSet(int nset, const std::string& filename, int member=0);
//    int    numberPDF  (int nset);
//    void   usePDFMember(int nset, int member);
//    double xfx(int nset, double x, double Q, int fl);
//    double getXmin(int nset, int member);
//    double getXmax(int nset, int member);
//    double getQ2min(int nset, int member);
//    double getQ2max(int nset, int member);
//    void   extrapolate(bool extrapolate=true);
//}



int main(int argc, const char* argv[])
{

  std::cout << "TestMENew" << std::endl;
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
  //float pertW1    ( in.getParameter<double> ("pertW1")     );
  //float pertW2    ( in.getParameter<double> ("pertW2")     );
  //float pertBHad  ( in.getParameter<double> ("pertBHad")   );
  //float pertBLep  ( in.getParameter<double> ("pertBLep")   );
  float enlargeE1 ( in.getParameter<double> ("enlargeE1")  );
  float enlargeEh1( in.getParameter<double> ("enlargeEh1") );
  float enlargePt ( in.getParameter<double> ("enlargePt")  );

  float scaleH  ( in.getParameter<double> ("scaleH")   );
  float scaleL  ( in.getParameter<double> ("scaleL")   );
  float scaleMET( in.getParameter<double> ("scaleMET") );

  int   useME     ( in.getParameter<int>    ("useME")      );
  int   useJac    ( in.getParameter<int>    ("useJac")     );
  int   useMET    ( in.getParameter<int>    ("useMET")     );
  int   useTF     ( in.getParameter<int>    ("useTF")      );
  int   usePDF    ( in.getParameter<int>    ("usePDF")     );
  int   shiftMomenta    ( in.getParameter<int>    ("shiftMomenta")     );
  int   testMassScan    ( in.getParameter<int>    ("testMassScan")     );
  int   testPermutations( in.getParameter<int>    ("testPermutations") );
  int   printP4         ( in.getParameter<int>    ("printP4")          );

  vector<double> masses       ( in.getParameter<vector<double> >  ("masses")        );
  vector<string> functions    ( in.getParameter<vector<string> >  ("functions")     );
 

  vector<int> evLimits( in.getParameter<vector<int> >  ("evLimits") );
  int evLow = evLimits[0];
  int evHigh = evLimits[1];

  int   mode         ( in.getUntrackedParameter<int>    ("mode", 0)          );

  //TFile* fout = 0;
  //fout = new TFile(outFileName.c_str(),"RECREATE");
  //fout->Close();
  //cout << "Deleting file..." << endl;
  //delete fout
  gSystem->Exec(("rm "+outFileName).c_str());

  TStopwatch* clock = new TStopwatch();
  TRandom3*   ran   = new TRandom3();
  ran->SetSeed(4321);

  //int printP4          = 0;
  //int shiftMomenta     = 0;
  //int testMassScan     = 0;
  //int testPermutations = 1;


  // 21 bins 80 - 180
  const int nMassPoints  = 31;                                                                                                                        
  double mH[nMassPoints] = 
    { /*50,  55, */ 60,  65,  70,  75,
      80 , 85,  90, 95, 100, 105, 110, 
      115, 120, 125, 130, 135, 140, 145, 150, 155, 
      160, 165, 170, 175, 180 ,185, 190, 195, 200, 205, 210
      //215, 220, 225, 230, 235, 240, 245, 250 
    }; 
  

  int par;
  switch(mode){
  case 0:
    par = 4;
    break;
  case 1:
    par = 6;
    break;
  case 2:
    par = 19;
    break;
  case 3:
    par = 19;
    break;
  default:
    par = 4;
    break;
  }

  MEIntegratorNew* meIntegrator = new MEIntegratorNew( pathToTF , par, int(verbose));
  if(mode == 0)
    meIntegrator->setIntType( MEIntegratorNew::SL2wj );
  else if(mode == 1)
    meIntegrator->setIntType( MEIntegratorNew::SL1wj );
  else if(mode == 2)
    meIntegrator->setIntType( MEIntegratorNew::SLXSec );
  else if(mode == 3)
    meIntegrator->setIntType( MEIntegratorNew::SLAcc );
  else{
    cout << "Unsupported mode... exit" << endl;
    return 0;
  }
  meIntegrator->setWeightNorm( MEIntegratorNew::None ); // divide weights by xsec
  //meIntegrator->setNormFormulas( TString(functions[0].c_str()),  
  //			 TString(functions[1].c_str()),  
  //			 TString(functions[2].c_str())  );

  ///////////////////////////


  TFile* fout_tmp0 = TFile::Open(outFileName.c_str(),"UPDATE");
  TH1F*  hInt     = new TH1F("hInt","",nMassPoints, mH[0]-2.5, mH[nMassPoints-1]+2.5);  
 
  meIntegrator->setMass( met );
  meIntegrator->setTopMass( 174.3 , 80.19);
  meIntegrator->setUseME (useME);
  meIntegrator->setUseJac(useJac);
  meIntegrator->setUseMET(useMET);
  meIntegrator->setUseTF (useTF);
  meIntegrator->setUsePDF(usePDF);

  double pi = TMath::Pi();
  double p = 0.;
  if(mode==2 || mode==3){
    for(int m = 0; m < nMassPoints ; m++){

      ROOT::Math::Functor toIntegrate(meIntegrator, &MEIntegratorNew::Eval, par);
      ROOT::Math::GSLMCIntegrator ig2( ROOT::Math::IntegrationMultiDim::kVEGAS , 1.e-12, 1.e-5, vegasPoints);
      ig2.SetFunction(toIntegrate);
      meIntegrator->SetPar(par);

      double xLmode2[19] = {  -1, -pi,    5.,-1.,-pi,-1, -pi,   5.,  -1, -pi, -1, -pi, -1, -pi,    5,-1, -pi, -1, -pi};
      double xUmode2[19] = {   1., pi,  500., 1., pi, 1,  pi, 500.,   1,  pi,  1,  pi,  1,  pi,  500, 1,  pi,  1,  pi};

      bool isThere = false;
      for(unsigned int mm = 0; mm < masses.size() ; mm++){
	if( TMath::Abs(mH[m]-masses[mm])<0.1 ) isThere = true;
      }
      if(!isThere) continue;

      meIntegrator->setMass( mH[m] );
      clock->Start();
      p = ig2.Integral(xLmode2, xUmode2);
      clock->Stop();
      cout << "Mass  " << mH[m] << " => p = " << p << endl;
      hInt->Fill( mH[m], p);
      cout << "[Done in " << clock->RealTime()/60. << " min]" << endl;
    }
  }
  cout << "Delete meIntegrator..." << endl;
  delete meIntegrator;
  delete clock; delete ran;
  cout << "Finished!!!" << endl;

  hInt->Write("hInt",TObject::kOverwrite);
  fout_tmp0->Close();
  return 0;
  ///////////////////////////



  //////////////////////
  //TH1F*  hPdfLumi     = new TH1F("hPdfLumi","",42, 425 , 2525.);

  //ROOT::Math::GSLMCIntegrator ig2Pdf("vegas", 1.e-12, 1.e-5, 1000);
  //ROOT::Math::Functor toIntegratePdf(meIntegrator, &MEIntegratorNew::EvalPdf, 1); 
  //ig2Pdf.SetFunction(toIntegratePdf);

  //double xLPdf[1] = {0.0000001};
  //double xUPdf[1] = {1.0};

  //for(int q = 0; q < 8/*84*/; q++){
  //double qValue = 450 + q*25; 
  //meIntegrator->setQ(qValue);

  //double xLPdf[1] = {qValue*qValue/8000./8000.};
  //double xUPdf[1] = {1.0};

  //double pPdf = ig2Pdf.Integral(xLPdf, xUPdf);
  //cout << pPdf << endl;
  //hPdfLumi->SetBinContent(hPdfLumi->FindBin(qValue), pPdf);
  //}
  //meIntegrator->setPartonLuminosity( hPdfLumi );
  //////////////////////

  //double Q1    =  600;
  //double Q2    =  800;
  //double Q3    = 1000;
  //double Q4    = 1500;

  //double cos3 =   0;
  //TH2F* meSquared1 = new TH2F("meSquared1","; m12 ; cos1Star", 50, 2*175, Q1-met, 50, -1,1);
  //TH2F* meSquared2 = new TH2F("meSquared2","; m12 ; cos1Star", 50, 2*175, Q2-met, 50, -1,1);
  //TH2F* meSquared3 = new TH2F("meSquared3","; m12 ; cos1Star", 50, 2*175, Q3-met, 50, -1,1);
  //TH2F* meSquared4 = new TH2F("meSquared4","; m12 ; cos1Star", 50, 2*175, Q4-met, 50, -1,1);
  /*
  for(int x = 1; x <=meSquared1->GetNbinsX(); x++){
    double m12 = meSquared1->GetXaxis()->GetBinCenter(x);
    for(int y = 1; y <=meSquared1->GetNbinsY(); y++){
      double cos1Star = meSquared1->GetYaxis()->GetBinCenter(y);
      meSquared1->SetBinContent( x,y,  meIntegrator->meSquaredAtQ(Q1, m12, cos1Star, cos3) );
    }
  }

  for(int x = 1; x <=meSquared2->GetNbinsX(); x++){
    double m12 = meSquared2->GetXaxis()->GetBinCenter(x);
    for(int y = 1; y <=meSquared2->GetNbinsY(); y++){
      double cos1Star = meSquared2->GetYaxis()->GetBinCenter(y);
      meSquared2->SetBinContent( x,y,  meIntegrator->meSquaredAtQ(Q2, m12, cos1Star, cos3) );
    }
  }

  for(int x = 1; x <=meSquared3->GetNbinsX(); x++){
    double m12 = meSquared3->GetXaxis()->GetBinCenter(x);
    for(int y = 1; y <=meSquared3->GetNbinsY(); y++){
      double cos1Star = meSquared3->GetYaxis()->GetBinCenter(y);
      meSquared3->SetBinContent( x,y,  meIntegrator->meSquaredAtQ(Q3, m12, cos1Star, cos3) );
    }
  }

  for(int x = 1; x <=meSquared4->GetNbinsX(); x++){
    double m12 = meSquared4->GetXaxis()->GetBinCenter(x);
    for(int y = 1; y <=meSquared4->GetNbinsY(); y++){
      double cos1Star = meSquared4->GetYaxis()->GetBinCenter(y);
      meSquared4->SetBinContent( x,y,  meIntegrator->meSquaredAtQ(Q4, m12, cos1Star, cos3) );
    }
  }
  */

  //////////////////////

  //ROOT::Math::GSLMCIntegrator ig2("vegas", 1.e-12, 1.e-5, vegasPoints);
  //ROOT::Math::Functor toIntegrate(meIntegrator, &MEIntegratorNew::Eval, par); 
  //ig2.SetFunction(toIntegrate);
  
  //const int nMassPoints = 30;
  //double mH[nMassPoints] = {95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240};
  //const int nMassPoints  = 18;                                                                                                                        
  //double mH[nMassPoints] = {90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175}; 
 

  


  TH1F*  hMass         = new TH1F("hMass","",nMassPoints, mH[0]-2.5, mH[nMassPoints-1]+2.5);
  TH1F*  hBestMass     = new TH1F("hBestMass","",40, 50, 250);
  TH1F*  hMassProb     = new TH1F("hMassProb","",nMassPoints, mH[0]-2.5, mH[nMassPoints-1]+2.5);
  TH1F*  hBestMassProb = new TH1F("hBestMassProb","",41, 47.5, 252.5);
  TTree* tProb         = new TTree("tree","");
  float prob_;
  float mass_;
  float initCpu_,initRea_,instCpu_, instRea_,evalCpu_,evalRea_;
  tProb->Branch("prob",&prob_,"prob/F");
  tProb->Branch("mass",&mass_,"mass/F");
  tProb->Branch("initCpu",&initCpu_,"initCpu/F");
  tProb->Branch("initRea",&initRea_,"initRea/F");
  tProb->Branch("instCpu",&instCpu_,"instCpu/F");
  tProb->Branch("instRea",&instRea_,"instRea/F");
  tProb->Branch("evalCpu",&evalCpu_,"evalCpu/F");
  tProb->Branch("evalRea",&evalRea_,"evalRea/F");

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
    Float_t vLepton_charge [99];

    currentTree->SetBranchAddress("genB",   &genB);
    currentTree->SetBranchAddress("genBbar",&genBbar);
    currentTree->SetBranchAddress("genTop", &genTop);
    currentTree->SetBranchAddress("genTbar",&genTbar);
    currentTree->SetBranchAddress("METtype1p2corr",&METtype1p2corr);
    currentTree->SetBranchAddress("vLepton_charge",vLepton_charge);
    
    int counter = 0;
    Long64_t nentries = currentTree->GetEntries();
    for (Long64_t i = 0; i < nentries ; i++){
      
      if(i%5000==0) cout << i << endl;
      currentTree->GetEntry(i);
      
      //continue;
      
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
      if(genB.mass>0 && genB.momid==25){
	genBLV.SetPtEtaPhiM(genB.pt,genB.eta ,genB.phi, genB.mass );
      }
      if(genBbar.mass>0 && genBbar.momid==25){
	genBbarLV.SetPtEtaPhiM(genBbar.pt,genBbar.eta ,genBbar.phi, genBbar.mass );
      }
  
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

      bool properEvent = (genBLV.Pt()>0 && genBbarLV.Pt()>0 && topBLV.Pt()>0 && topW1LV.Pt()>0 && topW2LV.Pt()>0 && atopBLV.Pt()>0 && atopW1LV.Pt()>0 && atopW2LV.Pt()>0);

      HIGGS.SetPxPyPzE( (genBLV+genBbarLV).Px(), (genBLV+genBbarLV).Py(),(genBLV+genBbarLV).Pz(),(genBLV+genBbarLV).E());
      if( abs(genTop.wdau1id)>6 && abs(genTbar.wdau1id)<6){
	TOPLEP.SetPxPyPzE( (topBLV+topW1LV+topW2LV).Px(), (topBLV+topW1LV+topW2LV).Py(), (topBLV+topW1LV+topW2LV).Pz(), (topBLV+topW1LV+topW2LV).E() );
	TOPHAD.SetPxPyPzE( (atopBLV+atopW1LV+atopW2LV).Px(), (atopBLV+atopW1LV+atopW2LV).Py(), (atopBLV+atopW1LV+atopW2LV).Pz(), (atopBLV+atopW1LV+atopW2LV).E() );
	TOPHADW1.SetPxPyPzE( atopW1LV.Px(), atopW1LV.Py(), atopW1LV.Pz(), atopW1LV.E());
	TOPHADW2.SetPxPyPzE( atopW2LV.Px(), atopW2LV.Py(), atopW2LV.Pz(), atopW2LV.E());
	TOPHADB.SetPxPyPzE( atopBLV.Px(),  atopBLV.Py(),   atopBLV.Pz(),  atopBLV.E());
	if( abs(genTop.wdau1id)==11 || abs(genTop.wdau1id)==13 || abs(genTop.wdau1id)==15 ){
	  TOPLEPW1.SetPxPyPzE  ( topW1LV.Px(), topW1LV.Py(), topW1LV.Pz(), topW1LV.E());
	  TOPLEPW2.SetPtEtaPhiM( topW2LV.Pt(), 0.0, topW2LV.Phi(), 0.0);
	  neutCosTheta = TMath::Cos( topW2LV.Vect().Theta() );
	}
	else{
	  TOPLEPW2.SetPtEtaPhiM( topW1LV.Pt(), 0.0, topW1LV.Phi(), 0.0);
	  TOPLEPW1.SetPxPyPzE  ( topW2LV.Px(), topW2LV.Py(), topW2LV.Pz(), topW2LV.E());
	  neutCosTheta = TMath::Cos( topW1LV.Vect().Theta() );
	}
	TOPLEPB.SetPxPyPzE(  topBLV.Px(),  topBLV.Py(),   topBLV.Pz(), topBLV.E());
      }
      else if(abs(genTop.wdau1id)<6 && abs(genTbar.wdau1id)>6){
	TOPHAD.SetPxPyPzE( (topBLV+topW1LV+topW2LV).Px(), (topBLV+topW1LV+topW2LV).Py(), (topBLV+topW1LV+topW2LV).Pz(), (topBLV+topW1LV+topW2LV).E() );
	TOPLEP.SetPxPyPzE( (atopBLV+atopW1LV+atopW2LV).Px(), (atopBLV+atopW1LV+atopW2LV).Py(), (atopBLV+atopW1LV+atopW2LV).Pz(), (atopBLV+atopW1LV+atopW2LV).E() );
	TOPHADW1.SetPxPyPzE( topW1LV.Px(), topW1LV.Py(), topW1LV.Pz(), topW1LV.E());
	TOPHADW2.SetPxPyPzE( topW2LV.Px(), topW2LV.Py(), topW2LV.Pz(), topW2LV.E());
	TOPHADB.SetPxPyPzE( topBLV.Px(),  topBLV.Py(),   topBLV.Pz(),  topBLV.E());
	if( abs(genTbar.wdau1id)==11 || abs(genTbar.wdau1id)==13 || abs(genTbar.wdau1id)==15 ){
	  TOPLEPW1.SetPxPyPzE  ( atopW1LV.Px(), atopW1LV.Py(), atopW1LV.Pz(), atopW1LV.E());
	  TOPLEPW2.SetPtEtaPhiM( atopW2LV.Pt(), 0.0, atopW2LV.Phi(), 0.0);
	  neutCosTheta = TMath::Cos( atopW2LV.Vect().Theta() );
	}
	else{
	  TOPLEPW2.SetPtEtaPhiM( atopW1LV.Pt(), 0.0, atopW1LV.Phi(), 0.0);
	  TOPLEPW1.SetPxPyPzE  ( atopW2LV.Px(), atopW2LV.Py(), atopW2LV.Pz(), atopW2LV.E());
	  neutCosTheta = TMath::Cos( atopW1LV.Vect().Theta() );
	}
	TOPLEPB.SetPxPyPzE( atopBLV.Px(),  atopBLV.Py(),   atopBLV.Pz(),  atopBLV.E());
      }      
      else{
	properEvent=false;
      }




      double errWjet1 = ran->Gaus( TOPHADW1.E(),  scaleL*TOPHADW1.E());
      double errWjet2 = ran->Gaus( TOPHADW2.E(),  scaleL*TOPHADW2.E());
      double errbHad  = ran->Gaus( TOPHADB.E()  *0.93 ,  scaleH*TOPHADB.E());
      double errbLep  = ran->Gaus( TOPLEPB.E()  *0.93 ,  scaleH*TOPLEPB.E());
      double errHiggs1= ran->Gaus( genBLV.E()   *0.93,   scaleH*genBLV.E());
      double errHiggs2= ran->Gaus( genBbarLV.E()*0.93,   scaleH*genBbarLV.E());
      double errMetPx = TOPLEPW2.Px() + ran->Gaus(0.0, scaleMET);
      double errMetPy = TOPLEPW2.Py() + ran->Gaus(0.0, scaleMET);
      double errMetPz = TOPLEPW2.Pz() ;
      
      TLorentzVector TOPLEPW1scaled  (TOPLEPW1.Vect()   , TOPLEPW1.E());
      TLorentzVector TOPLEPW2scaled;
      TOPLEPW2scaled.SetPxPyPzE(errMetPx,errMetPy,errMetPz,sqrt(errMetPx*errMetPx + errMetPy*errMetPy + errMetPz*errMetPz));
      TLorentzVector TOPLEPBscaled  (TOPLEPB.Vect()  *(errbLep/TOPLEPB.E())  , errbLep);
      TLorentzVector TOPHADW1scaled (TOPHADW1.Vect() *(errWjet1/TOPHADW1.E()), errWjet1);
      TLorentzVector TOPHADW2scaled (TOPHADW2.Vect() *(errWjet2/TOPHADW2.E()), errWjet2);
      TLorentzVector TOPHADBscaled  (TOPHADB.Vect()  *(errbHad/TOPHADB.E())  , errbHad);
      TLorentzVector genBLVscaled   (genBLV.Vect()   *(errHiggs1/genBLV.E()) , errHiggs1);
      TLorentzVector genBbarLVscaled(genBbarLV.Vect()*(errHiggs2/genBbarLV.E()) , errHiggs2);


      if(shiftMomenta==0){
	if(mode==0)
	  properEvent = ( TOPLEPW1.Pt() >20 && TMath::Abs(TOPLEPW1.Eta()) <2.1 &&
			  TOPLEPB.Pt()  >30 && TMath::Abs(TOPLEPB.Eta())  <2.5 &&
			  TOPHADW1.Pt() >30 && TMath::Abs(TOPHADW1.Eta()) <2.5 &&
			  TOPHADW2.Pt() >30 && TMath::Abs(TOPHADW2.Eta()) <2.5 &&
			  TOPHADB.Pt()  >30 && TMath::Abs(TOPHADB.Eta())  <2.5 &&
			  genBLV.Pt()   >30 && TMath::Abs(genBLV.Eta())   <2.5 &&
			  genBbarLV.Pt()>30 && TMath::Abs(genBbarLV.Eta())<2.5
			  );
	else if(mode==1)
	  properEvent = (( TOPLEPW1.Pt() >20 && TMath::Abs(TOPLEPW1.Eta()) <2.1 &&
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
	else{ }
      }
      else{
	if(mode==0) 
	  properEvent = ( TOPLEPW1scaled.Pt() >20 && TMath::Abs(TOPLEPW1.Eta()) <2.1 &&
			  TOPLEPBscaled.Pt()  >30 && TMath::Abs(TOPLEPB.Eta())  <2.5 &&
			  TOPHADW1scaled.Pt() >30 && TMath::Abs(TOPHADW1.Eta()) <2.5 &&
			  TOPHADW2scaled.Pt() >30 && TMath::Abs(TOPHADW2.Eta()) <2.5 &&
			  TOPHADBscaled.Pt()  >30 && TMath::Abs(TOPHADB.Eta())  <2.5 &&
			  genBLVscaled.Pt()   >30 && TMath::Abs(genBLV.Eta())   <2.5 &&
			  genBbarLVscaled.Pt()>30 && TMath::Abs(genBbarLV.Eta())<2.5
			  );
	else if(mode==1)
	  properEvent = (( TOPLEPW1scaled.Pt() >20 && TMath::Abs(TOPLEPW1.Eta()) <2.1 &&
			   TOPLEPBscaled.Pt()  >30 && TMath::Abs(TOPLEPB.Eta())  <2.5 &&
			   (TOPHADW1scaled.Pt() <30 || TMath::Abs(TOPHADW1.Eta()) >2.5) &&
			   TOPHADW2scaled.Pt() >30 && TMath::Abs(TOPHADW2.Eta()) <2.5 &&
			   TOPHADBscaled.Pt()  >30 && TMath::Abs(TOPHADB.Eta())  <2.5 &&
			   genBLVscaled.Pt()   >30 && TMath::Abs(genBLV.Eta())   <2.5 &&
			   genBbarLVscaled.Pt()>30 && TMath::Abs(genBbarLV.Eta())<2.5
			   ) ||
			 ( TOPLEPW1scaled.Pt() >20 && TMath::Abs(TOPLEPW1.Eta()) <2.1 &&
			   TOPLEPBscaled.Pt()  >30 && TMath::Abs(TOPLEPB.Eta())  <2.5 &&
			   TOPHADW1scaled.Pt() >30 && TMath::Abs(TOPHADW1.Eta()) <2.5 &&
			   (TOPHADW2scaled.Pt() <30 || TMath::Abs(TOPHADW2.Eta()) >2.5) &&
			   TOPHADBscaled.Pt()  >30 && TMath::Abs(TOPHADB.Eta())  <2.5 &&
			   genBLVscaled.Pt()   >30 && TMath::Abs(genBLV.Eta())   <2.5 &&
			   genBbarLVscaled.Pt()>30 && TMath::Abs(genBbarLV.Eta())<2.5
			   )
			 );
	else{ }
      }


      if(!properEvent) continue;
      counter++;

      if( !(counter>=evLow && counter<=evHigh) ) continue;
      cout << "Processing event # " << counter << endl;

      //TLorentzVector jet1,jet2,jet3,jet4,jet5,jet6,jet7,jet8;
      //jet1.SetPtEtaPhiM(19.8454,            0.822914,-2.3529,0.10566);     // lep
      //jet2.SetPtEtaPhiM(91.934,             0.0,-0.0161266,   0.);         // MET    0.543608 <=> cosTheta = 0.495714;
      //jet3.SetPtEtaPhiM(175.197 * pertBLep, 0.0936224,-0.969521,4.8);      // b from top lep
      //jet4.SetPtEtaPhiM(122.914 * pertW1,   1.17562,2.17839,0.33);         // j1 from W
      //jet5.SetPtEtaPhiM(34.9939 * pertW2,   1.06994,0.879248,0.33);        // j2 from W
      //jet6.SetPtEtaPhiM(95.5899 * pertBHad, 0.829954,3.10669,4.8);         // b from top hadr
      //jet7.SetPtEtaPhiM(32.805,             0.678347,-0.0951163,4.8);      // b1 from H
      //jet8.SetPtEtaPhiM(75.4309,            2.08021,2.75503,4.8);          // b2 from H

      prob_    = 0.;
      mass_    = -99;
      initCpu_ = -99;
      initRea_ = -99;
      instCpu_ = -99;
      evalCpu_ = -99;
      evalCpu_ = -99;
      evalRea_ = -99;

      if(testMassScan){

	vector<TLorentzVector> jets;
	jets.push_back( TOPLEPW1  );  
	jets.push_back( TOPLEPW2  );  
	jets.push_back( TOPLEPB   );
	jets.push_back( TOPHADW1  );  
	jets.push_back( TOPHADW2  );  
	jets.push_back( TOPHADB   );
	jets.push_back( genBLV    );
	jets.push_back( genBbarLV );

	if(printP4 && shiftMomenta==0){
	  cout << "*******START******" << endl;      
	  cout << "SumEt = " << METtype1p2corr.sumet  << endl;
	  cout << "lep:  jet1.SetPtEtaPhiM(" << TOPLEPW1.Pt() << "," <<  TOPLEPW1.Eta() << "," << TOPLEPW1.Phi() << "," << TOPLEPW1.M() << ")" << endl;
	  cout << "met:  jet2.SetPtEtaPhiM(" << TOPLEPW2.Pt() << "," <<  TOPLEPW2.Eta() << "," << TOPLEPW2.Phi() << "," << TOPLEPW2.M() << ")" << endl;
	  cout << "blep: jet3.SetPtEtaPhiM(" << TOPLEPB.Pt() << "," <<  TOPLEPB.Eta() << "," << TOPLEPB.Phi() << "," << TOPLEPB.M() << ")" << endl;
	  cout << "w1:   jet4.SetPtEtaPhiM(" << TOPHADW1.Pt() << "," <<  TOPHADW1.Eta() << "," << TOPHADW1.Phi() << "," << TOPHADW1.M() << ")" << endl;
	  cout << "w2:   jet5.SetPtEtaPhiM(" << TOPHADW2.Pt() << "," <<  TOPHADW2.Eta() << "," << TOPHADW2.Phi() << "," << TOPHADW2.M() << ")" << endl;
	  cout << "bhad: jet6.SetPtEtaPhiM(" << TOPHADB.Pt() << "," <<  TOPHADB.Eta() << "," << TOPHADB.Phi() << "," << TOPHADB.M() << ")" << endl;
	  cout << "h1:   jet7.SetPtEtaPhiM(" << genBLV.Pt() << "," <<  genBLV.Eta() << "," << genBLV.Phi() << "," << genBLV.M() << ")" << endl;
	  cout << "h2:   jet8.SetPtEtaPhiM(" << genBbarLV.Pt() << "," <<  genBbarLV.Eta() << "," << genBbarLV.Phi() << "," << genBbarLV.M() << ")" << endl;	
	  cout << "Top Lep mass = " << (TOPLEPW1+TOPLEPW2+TOPLEPB).M() << " <=> neutrino eta=0!!!" << endl;
	  cout << "Top Had mass = " << (TOPHADW1+TOPHADW2+TOPHADB).M() << endl;
	  cout << "Higgs mass = "   << (genBLV+genBbarLV).M() << endl;
	  cout << "Neut cosTheta = " << neutCosTheta << endl;
	}


	if(shiftMomenta){

	  /*
	  meIntegrator->setJets(&jets);
	  meIntegrator->initVersors(1);
	  meIntegrator->initTF();
	  
	  double errWjet1  = (meIntegrator->getCachedTF("tfWjet1"))->GetRandom()  ;
	  double errWjet2  = (meIntegrator->getCachedTF("tfWjet2"))->GetRandom()  ;
	  double errbHad   = (meIntegrator->getCachedTF("tfbHad"))->GetRandom()   ;
	  double errbLep   = (meIntegrator->getCachedTF("tfbLep"))->GetRandom()   ;
	  double errHiggs1 = (meIntegrator->getCachedTF("tfHiggs1"))->GetRandom() ;
	  double errHiggs2 = (meIntegrator->getCachedTF("tfHiggs2"))->GetRandom() ;
	  double errMetPt  = (meIntegrator->getCachedTF("tfMetPt"))->GetRandom()  ;
	  double errMetPhi = (((TH2F*)(meIntegrator->getCachedTF("tfMetPhi")))->ProjectionY("_py", 
											    (meIntegrator->getCachedTF("tfMetPhi"))->GetXaxis()->FindBin( TOPLEPW2.Pt() ),
											    (meIntegrator->getCachedTF("tfMetPhi"))->GetXaxis()->FindBin( TOPLEPW2.Pt() )  ))->GetRandom() ;
	  
	  meIntegrator->deleteTF();
	 
	  double errWjet1 = ran->Gaus( TOPHADW1.E(),  scaleL*TOPHADW1.E());
	  double errWjet2 = ran->Gaus( TOPHADW2.E(),  scaleL*TOPHADW2.E());
	  double errbHad  = ran->Gaus( TOPHADB.E()  *0.93 ,  scaleH*TOPHADB.E());
	  double errbLep  = ran->Gaus( TOPLEPB.E()  *0.93 ,  scaleH*TOPLEPB.E());
	  double errHiggs1= ran->Gaus( genBLV.E()   *0.93,   scaleH*genBLV.E());
	  double errHiggs2= ran->Gaus( genBbarLV.E()*0.93,   scaleH*genBbarLV.E());
	  double errMetPx = TOPLEPW2.Px() + ran->Gaus(0.0, scaleMET);
	  double errMetPy = TOPLEPW2.Py() + ran->Gaus(0.0, scaleMET);
	  double errMetPz = TOPLEPW2.Pz() ;

	  TLorentzVector TOPLEPW1scaled  (TOPLEPW1.Vect()   , TOPLEPW1.E());
	  TLorentzVector TOPLEPW2scaled;
	  TOPLEPW2scaled.SetPxPyPzE(errMetPx,errMetPy,errMetPz,sqrt(errMetPx*errMetPx + errMetPy*errMetPy + errMetPz*errMetPz));
	  TLorentzVector TOPLEPBscaled  (TOPLEPB.Vect()  *(errbLep/TOPLEPB.E())  , errbLep);
	  TLorentzVector TOPHADW1scaled (TOPHADW1.Vect() *(errWjet1/TOPHADW1.E()), errWjet1);
	  TLorentzVector TOPHADW2scaled (TOPHADW2.Vect() *(errWjet2/TOPHADW2.E()), errWjet2);
	  TLorentzVector TOPHADBscaled  (TOPHADB.Vect()  *(errbHad/TOPHADB.E())  , errbHad);
	  TLorentzVector genBLVscaled   (genBLV.Vect()   *(errHiggs1/genBLV.E()) , errHiggs1);
	  TLorentzVector genBbarLVscaled(genBbarLV.Vect()*(errHiggs2/genBbarLV.E()) , errHiggs2);
	  */
  
	  if(printP4){
	    cout << "*******START******" << endl;      
	    cout << "SumEt = " << METtype1p2corr.sumet  << endl;
	    cout << "lep   (scaled): jet1.SetPtEtaPhiM(" << TOPLEPW1scaled.Pt() << "," <<  TOPLEPW1scaled.Eta() << "," << TOPLEPW1scaled.Phi() << "," << TOPLEPW1scaled.M() << ")" << endl;
	    cout << "met:  (scaled)  jet2.SetPtEtaPhiM(" << TOPLEPW2scaled.Pt() << "," <<  TOPLEPW2scaled.Eta() << "," << TOPLEPW2scaled.Phi() << "," << TOPLEPW2scaled.M() << ")" << endl;
	    cout << "blep: (scaled)  jet3.SetPtEtaPhiM(" << TOPLEPBscaled.Pt() << "," <<  TOPLEPBscaled.Eta() << "," << TOPLEPBscaled.Phi() << "," << TOPLEPBscaled.M() << ")" << endl;
	    cout << "w1:   (scaled)  jet4.SetPtEtaPhiM(" << TOPHADW1scaled.Pt() << "," <<  TOPHADW1scaled.Eta() << "," << TOPHADW1scaled.Phi() << "," << TOPHADW1scaled.M() << ")" << endl;
	    cout << "w2:   (scaled)  jet5.SetPtEtaPhiM(" << TOPHADW2scaled.Pt() << "," <<  TOPHADW2scaled.Eta() << "," << TOPHADW2scaled.Phi() << "," << TOPHADW2scaled.M() << ")" << endl;
	    cout << "bhad: (scaled)  jet6.SetPtEtaPhiM(" << TOPHADBscaled.Pt() << "," <<  TOPHADBscaled.Eta() << "," << TOPHADBscaled.Phi() << "," << TOPHADBscaled.M() << ")" << endl;
	    cout << "h1:   (scaled)  jet7.SetPtEtaPhiM(" << genBLVscaled.Pt() << "," <<  genBLVscaled.Eta() << "," << genBLVscaled.Phi() << "," << genBLVscaled.M() << ")" << endl;
	    cout << "h2:   (scaled)  jet8.SetPtEtaPhiM(" << genBbarLVscaled.Pt() << "," <<  genBbarLVscaled.Eta() << "," << genBbarLVscaled.Phi() << "," << genBbarLVscaled.M() << ")" << endl;	
	    cout << "Top Lep mass (scaled) = " << (TOPLEPW1scaled+TOPLEPW2scaled+TOPLEPBscaled).M() << " <=> neutrino eta=0!!!" << endl;
	    cout << "Top Had mass (scaled) = " << (TOPHADW1scaled+TOPHADW2scaled+TOPHADBscaled).M() << endl;
	    cout << "Higgs mass (scaled) = "   << (genBLVscaled+genBbarLVscaled).M() << endl;
	  }

	  jets.clear();
	  jets.push_back( TOPLEPW1scaled  );  
	  jets.push_back( TOPLEPW2scaled  );  
	  jets.push_back( TOPLEPBscaled   );
	  jets.push_back( TOPHADW1scaled  );  
	  jets.push_back( TOPHADW2scaled  );  
	  jets.push_back( TOPHADBscaled   );
	  jets.push_back( genBLVscaled    );
	  jets.push_back( genBbarLVscaled );

	  /*
	  bool properScaledEvent = false;
	  if(mode==0) 
	    properScaledEvent = ( TOPLEPW1scaled.Pt() >20 && TMath::Abs(TOPLEPW1.Eta()) <2.1 &&
				  TOPLEPBscaled.Pt()  >30 && TMath::Abs(TOPLEPB.Eta())  <2.5 &&
				  TOPHADW1scaled.Pt() >30 && TMath::Abs(TOPHADW1.Eta()) <2.5 &&
				  TOPHADW2scaled.Pt() >30 && TMath::Abs(TOPHADW2.Eta()) <2.5 &&
				  TOPHADBscaled.Pt()  >30 && TMath::Abs(TOPHADB.Eta())  <2.5 &&
				  genBLVscaled.Pt()   >30 && TMath::Abs(genBLV.Eta())   <2.5 &&
				  genBbarLVscaled.Pt()>30 && TMath::Abs(genBbarLV.Eta())<2.5
				  );
	  else if(mode==1)
	    properScaledEvent = (( TOPLEPW1scaled.Pt() >20 && TMath::Abs(TOPLEPW1.Eta()) <2.1 &&
				  TOPLEPBscaled.Pt()  >30 && TMath::Abs(TOPLEPB.Eta())  <2.5 &&
				  (TOPHADW1scaled.Pt() <30 || TMath::Abs(TOPHADW1.Eta()) >2.5) &&
				  TOPHADW2scaled.Pt() >30 && TMath::Abs(TOPHADW2.Eta()) <2.5 &&
				  TOPHADBscaled.Pt()  >30 && TMath::Abs(TOPHADB.Eta())  <2.5 &&
				  genBLVscaled.Pt()   >30 && TMath::Abs(genBLV.Eta())   <2.5 &&
				  genBbarLVscaled.Pt()>30 && TMath::Abs(genBbarLV.Eta())<2.5
				  ) ||
				 ( TOPLEPW1scaled.Pt() >20 && TMath::Abs(TOPLEPW1.Eta()) <2.1 &&
				   TOPLEPBscaled.Pt()  >30 && TMath::Abs(TOPLEPB.Eta())  <2.5 &&
				   TOPHADW1scaled.Pt() >30 && TMath::Abs(TOPHADW1.Eta()) <2.5 &&
				   (TOPHADW2scaled.Pt() <30 || TMath::Abs(TOPHADW2.Eta()) >2.5) &&
				   TOPHADBscaled.Pt()  >30 && TMath::Abs(TOPHADB.Eta())  <2.5 &&
				   genBLVscaled.Pt()   >30 && TMath::Abs(genBLV.Eta())   <2.5 &&
				   genBbarLVscaled.Pt()>30 && TMath::Abs(genBbarLV.Eta())<2.5
				  )
				 );
	  else{ }

	  if( shiftMomenta && !properScaledEvent ) continue;
	  */	  
	}

	if(mode==1){
	  TLorentzVector w1 = jets[3];
	  TLorentzVector w2 = jets[4];
	  if( (w1.Pt()<30 || TMath::Abs(w1.Eta())>2.5) && (w2.Pt()>30 && TMath::Abs(w2.Eta())<2.5)){
	    jets[3] = w2;
	    jets[4] = w1;
	  }
	  else if( (w2.Pt()<30 || TMath::Abs(w2.Eta())>2.5) && (w1.Pt()>30 && TMath::Abs(w1.Eta())<2.5)){
	    jets[3] = w1;
	    jets[4] = w2;
	  }
	  else{
	    cout << "Inconsistentcy of mode and jet selections" << endl;
	    return 1;
	  }
	}


	/////////////////////////////////////////////////////////////////////////////
	clock->Start();

	if(printP4) cout << ">>>" << endl;
	if(printP4) cout << "Setup jets..." << endl;
	meIntegrator->setJets(&jets);
	if(printP4) cout << "Setup H mass..." << endl;
	meIntegrator->setMass( met );
	meIntegrator->setTopMass( 174.3 , 80.19);
	meIntegrator->setSumEt( METtype1p2corr.sumet );
	if(printP4) cout << "Setup versors..." << endl;
	meIntegrator->initVersors(234567);
	if(printP4) cout << "Smear by TF (?)" << endl;
	if( !(meIntegrator->smearByTF(30., 1)) ) continue;
	if(printP4) cout << "Setup TFs..." << endl;
	meIntegrator->initTF();
	if(printP4) cout << "Build integrand..." << endl;	
	meIntegrator->setTopFlags( vLepton_charge[0]==1 ? +1 : -1 , vLepton_charge[0]==1 ? -1 : +1 );
	//meIntegrator->setTopFlags( 0, 0 );
	cout << "Lep charge = " << vLepton_charge[0] << " => Top flags = (" << (vLepton_charge[0]==1 ? +1 : -1) << "," << (vLepton_charge[0]==1 ? -1 : +1) << ")" << endl;
	meIntegrator->setUseME (useME);
	meIntegrator->setUseJac(useJac);
	meIntegrator->setUseMET(useMET);
	meIntegrator->setUseTF (useTF);
	meIntegrator->setUsePDF(usePDF);
	if(printP4) cout << ">>>" << endl;
      
	clock->Stop();

	initCpu_ = clock->CpuTime();
	initRea_ = clock->RealTime();

	clock->Reset();

	/////////////////////////////////////////////////////////////////////////////
    
	double E1low   = (meIntegrator->getCachedTF("tfWjet1"))->GetXaxis()->GetXmin() * (1-enlargeE1);
	double E1high  = (meIntegrator->getCachedTF("tfWjet1"))->GetXaxis()->GetXmax() * (1+enlargeE1);
	
	double Ptlow   = (meIntegrator->getCachedTF("tfMetPt"))->GetXaxis()->GetXmin() * (1-enlargePt);
	double Pthigh  = (meIntegrator->getCachedTF("tfMetPt"))->GetXaxis()->GetXmax() * (1+enlargePt);
	
	Ptlow  = -1; //0.781794 - 0.05;
	Pthigh = +1; //0.781794 + 0.05;
      
	double Philow   = -(meIntegrator->getCachedTF("tfMetPhi"))->GetYaxis()->GetXmax();
	double Phihigh  =  (meIntegrator->getCachedTF("tfMetPhi"))->GetYaxis()->GetXmax();
            
	//Philow  = -TMath::Pi();
	//Phihigh = +TMath::Pi();

	double Eh1low   = (meIntegrator->getCachedTF("tfHiggs1"))->GetXaxis()->GetXmin() * (1-enlargeEh1);
	double Eh1high  = (meIntegrator->getCachedTF("tfHiggs1"))->GetXaxis()->GetXmax() * (1+enlargeEh1);
      
	/////////////////////////////////////////////////////////////////////////////
      
	if(printP4){
	  cout << "E1 :       ["  << E1low  << "," << E1high  << "]" << endl;
	  if(mode==1){
	    cout << "cosThetaMiss : ["  << -1  << "," << 1  << "]" << endl;
	    cout << "PhiMiss :      ["  << -TMath::Pi() << "," << TMath::Pi() << "]" << endl;	
	  }
	  cout << "cosTheta : ["  << Ptlow  << "," << Pthigh  << "]" << endl;
	  cout << "Phi :      ["  << Philow << "," << Phihigh << "]" << endl;
	  cout << "Eh1 :      ["  << Eh1low << "," << Eh1high << "]" << endl;
	}

	/////////////////////////////////////////////////////////////////////////////
	
	double xLmode0[4] = {  E1low,  Ptlow,   Philow,   Eh1low};
	double xUmode0[4] = {  E1high, Pthigh , Phihigh,  Eh1high};

	double xLmode1[6] = {  E1low,  -1, -TMath::Pi(), Ptlow,   Philow,   Eh1low};
	double xUmode1[6] = {  E1high, +1, +TMath::Pi(), Pthigh , Phihigh,  Eh1high};

      
	clock->Start();

	ROOT::Math::Functor toIntegrate(meIntegrator, &MEIntegratorNew::Eval, par);
	ROOT::Math::GSLMCIntegrator ig2( ROOT::Math::IntegrationMultiDim::kVEGAS , 1.e-12, 1.e-5, vegasPoints);
	ig2.SetFunction(toIntegrate);
	meIntegrator->SetPar(par);

	clock->Stop();

	instCpu_ = clock->CpuTime();
	instRea_ = clock->RealTime();

	clock->Reset();

	evalCpu_ = 0.;
	evalRea_ = 0.;
	hMass->Reset();
	for(int m = 0; m < nMassPoints ; m++){
	  meIntegrator->setMass( mH[m] );
	  //meIntegrator->setTopMass( mH[m], 80.19);

	  clock->Start();
	  double p = 0.;
	  if(mode==0)
	    p = ig2.Integral(xLmode0, xUmode0);
	  else if(mode==1)
	    p = ig2.Integral(xLmode1, xUmode1);
	  else{ }

	  clock->Stop();
	  evalCpu_ = clock->CpuTime();
	  evalRea_ = clock->RealTime();
	  clock->Reset();

	  if(printP4) cout << "Mass " << mH[m] << " => prob  = " << p << endl;
	  hMass->Fill(mH[m], p);
	  if(mH[m]>122.5 && mH[m]<127.5){
	    prob_ = p;
	    tProb->Fill();
	  }
	}

	TF1* gaus = new TF1("gaus","gaus", 0,300);
	//gaus->SetParameter(1, hMass->GetBinCenter(hMass->GetMaximumBin()) );
	//gaus->SetParameter(2, 20. );
	//hBestMass->Fill( hMass->GetBinCenter(hMass->GetMaximumBin())  );
	hMass->Fit(gaus, "Q", "", 
		   hMass->GetBinCenter(hMass->GetMaximumBin()) - 20, 
		   hMass->GetBinCenter(hMass->GetMaximumBin()) + 20);
	hBestMass->Fill( gaus->GetParameter(1)  );
	delete gaus;

	if(printP4) cout << "Opening file to save histos" << endl;
	TFile* fout_tmp = TFile::Open(outFileName.c_str(),"UPDATE");
	TDirectory *dir = fout_tmp->mkdir(Form("Event_%d", counter));
	dir->cd();
	(meIntegrator->getCachedPdf("pdfGammaWHad"))->Write("pdfGammaWHad",TObject::kOverwrite);
	(meIntegrator->getCachedPdf("pdfBetaWHad")) ->Write("pdfBetaWHad",TObject::kOverwrite);
	(meIntegrator->getCachedPdf("pdfGammaWLep"))->Write("pdfGammaWLep",TObject::kOverwrite);
	(meIntegrator->getCachedPdf("pdfBetaWLep")) ->Write("pdfBetaWLep",TObject::kOverwrite);
	(meIntegrator->getCachedPdf("pdfGammaTTH")) ->Write("pdfGammaTTH",TObject::kOverwrite);
	(meIntegrator->getCachedPdf("pdf3D"))       ->Write("pdf3D",TObject::kOverwrite);
	(meIntegrator->getCachedPdf("pdf2D"))       ->Write("pdf2D",TObject::kOverwrite);
	(meIntegrator->getCachedTF("tfWjet1"))      ->Write("tfWjet1",TObject::kOverwrite);
	(meIntegrator->getCachedTF("tfWjet2"))      ->Write("tfWjet2",TObject::kOverwrite);  
	(meIntegrator->getCachedTF("tfbHad"))       ->Write("tfbHad",TObject::kOverwrite);
	(meIntegrator->getCachedTF("tfbLep"))       ->Write("tfbLep",TObject::kOverwrite);
	(meIntegrator->getCachedTF("tfHiggs1"))     ->Write("tfHiggs1",TObject::kOverwrite);
	(meIntegrator->getCachedTF("tfHiggs2"))     ->Write("tfHiggs2",TObject::kOverwrite);
	(meIntegrator->getCachedTF("tfMetPt"))      ->Write("tfMetPt",TObject::kOverwrite);
	(meIntegrator->getCachedTF("tfMetPhi"))     ->Write("tfMetPhi",TObject::kOverwrite);
	hMass->Write("hMass",TObject::kOverwrite);
	if(printP4) cout << "Closing file to save histos" << endl;
	fout_tmp->Close();
	

	meIntegrator->deleteTF();
	if(printP4) cout << "*******END********" << endl;      
      }

      prob_    = 0.;
      mass_    = -99;
      initCpu_ = -99;
      initRea_ = -99;
      instCpu_ = -99;
      evalCpu_ = -99;
      evalCpu_ = -99;
      evalRea_ = -99;

      //////////////////////////////////////////////////////////
      // test permutations
      if(testPermutations) {

	/*
	double errWjet1 = ran->Gaus( TOPHADW1.E(),         scaleL*TOPHADW1.E());
	double errWjet2 = ran->Gaus( TOPHADW2.E(),         scaleL*TOPHADW2.E());
	double errbHad  = ran->Gaus( TOPHADB.E()  *0.93 ,  scaleH*TOPHADB.E());
	double errbLep  = ran->Gaus( TOPLEPB.E()  *0.93 ,  scaleH*TOPLEPB.E());
	double errHiggs1= ran->Gaus( genBLV.E()   *0.93,   scaleH*genBLV.E());
	double errHiggs2= ran->Gaus( genBbarLV.E()*0.93,   scaleH*genBbarLV.E());
	double errMetPx = TOPLEPW2.Px() + ran->Gaus(0.0,   scaleMET);
	double errMetPy = TOPLEPW2.Py() + ran->Gaus(0.0,   scaleMET);
	double errMetPz = TOPLEPW2.Pz() ;
	
	TLorentzVector TOPLEPW1scaled  (TOPLEPW1.Vect()   , TOPLEPW1.E());
	TLorentzVector TOPLEPW2scaled;
	TOPLEPW2scaled.SetPxPyPzE(errMetPx,errMetPy,errMetPz,sqrt(errMetPx*errMetPx + errMetPy*errMetPy + errMetPz*errMetPz));
	TLorentzVector TOPLEPBscaled  (TOPLEPB.Vect()  *(errbLep/TOPLEPB.E())  , errbLep);
	TLorentzVector TOPHADW1scaled (TOPHADW1.Vect() *(errWjet1/TOPHADW1.E()), errWjet1);
	TLorentzVector TOPHADW2scaled (TOPHADW2.Vect() *(errWjet2/TOPHADW2.E()), errWjet2);
	TLorentzVector TOPHADBscaled  (TOPHADB.Vect()  *(errbHad/TOPHADB.E())  , errbHad);
	TLorentzVector genBLVscaled   (genBLV.Vect()   *(errHiggs1/genBLV.E()) , errHiggs1);
	TLorentzVector genBbarLVscaled(genBbarLV.Vect()*(errHiggs2/genBbarLV.E()) , errHiggs2);
		    
	bool properScaledEvent = false;
	if(mode==0) 
	  properScaledEvent = ( TOPLEPW1scaled.Pt() >20 && TMath::Abs(TOPLEPW1.Eta()) <2.1 &&
				TOPLEPBscaled.Pt()  >30 && TMath::Abs(TOPLEPB.Eta())  <2.5 &&
				TOPHADW1scaled.Pt() >30 && TMath::Abs(TOPHADW1.Eta()) <2.5 &&
				TOPHADW2scaled.Pt() >30 && TMath::Abs(TOPHADW2.Eta()) <2.5 &&
				TOPHADBscaled.Pt()  >30 && TMath::Abs(TOPHADB.Eta())  <2.5 &&
				genBLVscaled.Pt()   >30 && TMath::Abs(genBLV.Eta())   <2.5 &&
				genBbarLVscaled.Pt()>30 && TMath::Abs(genBbarLV.Eta())<2.5
				);
	else if(mode==1)
	  properScaledEvent = (( TOPLEPW1scaled.Pt() >20 && TMath::Abs(TOPLEPW1.Eta()) <2.1 &&
				 TOPLEPBscaled.Pt()  >30 && TMath::Abs(TOPLEPB.Eta())  <2.5 &&
				 (TOPHADW1scaled.Pt() <30 || TMath::Abs(TOPHADW1.Eta()) >2.5) &&
				 TOPHADW2scaled.Pt() >30 && TMath::Abs(TOPHADW2.Eta()) <2.5 &&
				 TOPHADBscaled.Pt()  >30 && TMath::Abs(TOPHADB.Eta())  <2.5 &&
				 genBLVscaled.Pt()   >30 && TMath::Abs(genBLV.Eta())   <2.5 &&
				 genBbarLVscaled.Pt()>30 && TMath::Abs(genBbarLV.Eta())<2.5
				 ) ||
			       ( TOPLEPW1scaled.Pt() >20 && TMath::Abs(TOPLEPW1.Eta()) <2.1 &&
				 TOPLEPBscaled.Pt()  >30 && TMath::Abs(TOPLEPB.Eta())  <2.5 &&
				 TOPHADW1scaled.Pt() >30 && TMath::Abs(TOPHADW1.Eta()) <2.5 &&
				 (TOPHADW2scaled.Pt() <30 || TMath::Abs(TOPHADW2.Eta()) >2.5) &&
				 TOPHADBscaled.Pt()  >30 && TMath::Abs(TOPHADB.Eta())  <2.5 &&
				 genBLVscaled.Pt()   >30 && TMath::Abs(genBLV.Eta())   <2.5 &&
				 genBbarLVscaled.Pt()>30 && TMath::Abs(genBbarLV.Eta())<2.5
				 )
			       );

	if( shiftMomenta && !properScaledEvent ) continue;
	*/	   

	hMassProb->Reset();
	vector<TLorentzVector> bjets;
	if(shiftMomenta==0){
	  bjets.push_back( TOPLEPB );
	  bjets.push_back( TOPHADB );
	  bjets.push_back( genBLV );
	  bjets.push_back( genBbarLV );
	}
	else{
	  bjets.push_back( TOPLEPBscaled );
	  bjets.push_back( TOPHADBscaled );
	  bjets.push_back( genBLVscaled );
	  bjets.push_back( genBbarLVscaled );
	}


 
	if(printP4 && shiftMomenta==0){
	  cout << "*******START******" << endl;      
	  cout << "SumEt = " << METtype1p2corr.sumet  << endl;
	  cout << "lep   (): jet1.SetPtEtaPhiM(" << TOPLEPW1.Pt() << "," <<  TOPLEPW1.Eta() << "," << TOPLEPW1.Phi() << "," << TOPLEPW1.M() << ")" << endl;
	  cout << "met:  ()  jet2.SetPtEtaPhiM(" << TOPLEPW2.Pt() << "," <<  TOPLEPW2.Eta() << "," << TOPLEPW2.Phi() << "," << TOPLEPW2.M() << ")" << endl;
	  cout << "blep: ()  jet3.SetPtEtaPhiM(" << TOPLEPB.Pt() << "," <<  TOPLEPB.Eta() << "," << TOPLEPB.Phi() << "," << TOPLEPB.M() << ")" << endl;
	  cout << "w1:   ()  jet4.SetPtEtaPhiM(" << TOPHADW1.Pt() << "," <<  TOPHADW1.Eta() << "," << TOPHADW1.Phi() << "," << TOPHADW1.M() << ")" << endl;
	  cout << "w2:   ()  jet5.SetPtEtaPhiM(" << TOPHADW2.Pt() << "," <<  TOPHADW2.Eta() << "," << TOPHADW2.Phi() << "," << TOPHADW2.M() << ")" << endl;
	  cout << "bhad: ()  jet6.SetPtEtaPhiM(" << TOPHADB.Pt() << "," <<  TOPHADB.Eta() << "," << TOPHADB.Phi() << "," << TOPHADB.M() << ")" << endl;
	  cout << "h1:   ()  jet7.SetPtEtaPhiM(" << genBLV.Pt() << "," <<  genBLV.Eta() << "," << genBLV.Phi() << "," << genBLV.M() << ")" << endl;
	  cout << "h2:   ()  jet8.SetPtEtaPhiM(" << genBbarLV.Pt() << "," <<  genBbarLV.Eta() << "," << genBbarLV.Phi() << "," << genBbarLV.M() << ")" << endl;	
	  cout << "Top Lep mass () = " << (TOPLEPW1+TOPLEPW2+TOPLEPB).M() << " <=> neutrino eta=0!!!" << endl;
	  cout << "Top Had mass () = " << (TOPHADW1+TOPHADW2+TOPHADB).M() << endl;
	  cout << "Higgs mass () = "   << (genBLV+genBbarLV).M() << endl;
	}
	if(printP4 && shiftMomenta==1){
	  cout << "*******START******" << endl;      
	  cout << "SumEt = " << METtype1p2corr.sumet  << endl;
	  cout << "lep   (scaled): jet1.SetPtEtaPhiM(" << TOPLEPW1scaled.Pt() << "," <<  TOPLEPW1scaled.Eta() << "," << TOPLEPW1scaled.Phi() << "," << TOPLEPW1scaled.M() << ")" << endl;
	  cout << "met:  (scaled)  jet2.SetPtEtaPhiM(" << TOPLEPW2scaled.Pt() << "," <<  TOPLEPW2scaled.Eta() << "," << TOPLEPW2scaled.Phi() << "," << TOPLEPW2scaled.M() << ")" << endl;
	  cout << "blep: (scaled)  jet3.SetPtEtaPhiM(" << TOPLEPBscaled.Pt() << "," <<  TOPLEPBscaled.Eta() << "," << TOPLEPBscaled.Phi() << "," << TOPLEPBscaled.M() << ")" << endl;
	  cout << "w1:   (scaled)  jet4.SetPtEtaPhiM(" << TOPHADW1scaled.Pt() << "," <<  TOPHADW1scaled.Eta() << "," << TOPHADW1scaled.Phi() << "," << TOPHADW1scaled.M() << ")" << endl;
	  cout << "w2:   (scaled)  jet5.SetPtEtaPhiM(" << TOPHADW2scaled.Pt() << "," <<  TOPHADW2scaled.Eta() << "," << TOPHADW2scaled.Phi() << "," << TOPHADW2scaled.M() << ")" << endl;
	  cout << "bhad: (scaled)  jet6.SetPtEtaPhiM(" << TOPHADBscaled.Pt() << "," <<  TOPHADBscaled.Eta() << "," << TOPHADBscaled.Phi() << "," << TOPHADBscaled.M() << ")" << endl;
	  cout << "h1:   (scaled)  jet7.SetPtEtaPhiM(" << genBLVscaled.Pt() << "," <<  genBLVscaled.Eta() << "," << genBLVscaled.Phi() << "," << genBLVscaled.M() << ")" << endl;
	  cout << "h2:   (scaled)  jet8.SetPtEtaPhiM(" << genBbarLVscaled.Pt() << "," <<  genBbarLVscaled.Eta() << "," << genBbarLVscaled.Phi() << "," << genBbarLVscaled.M() << ")" << endl;	
	  cout << "Top Lep mass (scaled) = " << (TOPLEPW1scaled+TOPLEPW2scaled+TOPLEPBscaled).M() << " <=> neutrino eta=0!!!" << endl;
	  cout << "Top Had mass (scaled) = " << (TOPHADW1scaled+TOPHADW2scaled+TOPHADBscaled).M() << endl;
	  cout << "Higgs mass (scaled) = "   << (genBLVscaled+genBbarLVscaled).M() << endl;
	}



	for(int m = 0; m < nMassPoints ; m++){
	  cout << "Processing " << mH[m] << endl;
	  meIntegrator->setMass( mH[m] );
	  for(unsigned int pp1=0; pp1!=bjets.size(); pp1++){
	    for(unsigned int pp2=0; pp2!=bjets.size(); pp2++){
	      if(pp2==pp1) continue;
	      for(unsigned int pp3=0; pp3!=bjets.size()-1; pp3++){
		if(pp3==pp1 || pp3==pp2) continue;
		for(unsigned int pp4=pp3+1; pp4!=bjets.size(); pp4++){
		  if(pp4==pp1 || pp4==pp2 || pp4==pp3) continue;
		  //cout << pp1 << pp2 << pp3 << pp4 << endl;

		  //if(/*pp1==0 && pp2==1 &&*/ pp3==2 && pp4==3) continue;
		  //if(/*pp1==0 && pp2==1 &&*/ pp3==3 && pp4==2) continue;
		  float massPair = (bjets[pp3]+bjets[pp4]).M(); 
		  if( TMath::Abs( massPair - mH[m]) > 2.0*scaleH*massPair ){
		    //cout << "Pair mass: " << massPair << " - mass under test: " << mH[m] << endl; 
		    continue; // 2sigmas away...
		  }


		  vector<TLorentzVector> jets;
		  if(shiftMomenta==0) {
		    jets.push_back( TOPLEPW1 );
		    jets.push_back( TOPLEPW2 );
		  }
		  else{
		    jets.push_back( TOPLEPW1scaled );
		    jets.push_back( TOPLEPW2scaled );
		  }
		  jets.push_back( bjets[pp1] );
		  if(shiftMomenta==0) {
		    jets.push_back( TOPHADW1 );
		    jets.push_back( TOPHADW2 );
		  }
		  else{
		    jets.push_back( TOPHADW1scaled );
		    jets.push_back( TOPHADW2scaled );
		  }
		  jets.push_back( bjets[pp2] );
		  jets.push_back( bjets[pp3] );
		  jets.push_back( bjets[pp4] );


		  if(mode==1){
		    TLorentzVector w1 = jets[3];
		    TLorentzVector w2 = jets[4];
		    if( (w1.Pt()<30 || TMath::Abs(w1.Eta())>2.5) && (w2.Pt()>30 && TMath::Abs(w2.Eta())<2.5)){
		      jets[3] = w2;
		      jets[4] = w1;
		    }
		    else if( (w2.Pt()<30 || TMath::Abs(w2.Eta())>2.5) && (w1.Pt()>30 && TMath::Abs(w1.Eta())<2.5)){
		      jets[3] = w1;
		      jets[4] = w2;
		    }
		    else{
		      cout << "Inconsistentcy of mode and jet selections" << endl;
		      continue;
		    }
		  }


		  meIntegrator->setSumEt( METtype1p2corr.sumet );
		  meIntegrator->setJets(&jets);
		  meIntegrator->initVersors(1);
		  meIntegrator->initTF();
		  meIntegrator->setTopFlags( vLepton_charge[0]==1 ? +1 : -1 , vLepton_charge[0]==1 ? -1 : +1 );
		  meIntegrator->setUseME (useME);
		  meIntegrator->setUseJac(useJac);
		  meIntegrator->setUseMET(useMET);
		  meIntegrator->setUseTF (useTF);
		  meIntegrator->setUsePDF(usePDF);

		  double E1low   = (meIntegrator->getCachedTF("tfWjet1"))->GetXaxis()->GetXmin() * (1-enlargeE1);
		  double E1high  = (meIntegrator->getCachedTF("tfWjet1"))->GetXaxis()->GetXmax() * (1+enlargeE1);
		  double Ptlow   = (meIntegrator->getCachedTF("tfMetPt"))->GetXaxis()->GetXmin() * (1-enlargePt);
		  double Pthigh  = (meIntegrator->getCachedTF("tfMetPt"))->GetXaxis()->GetXmax() * (1+enlargePt);
		  Ptlow  = -1; 
		  Pthigh = +1; 
		  double Philow   = -(meIntegrator->getCachedTF("tfMetPhi"))->GetYaxis()->GetXmax();
		  double Phihigh  =  (meIntegrator->getCachedTF("tfMetPhi"))->GetYaxis()->GetXmax();
		  double Eh1low   = (meIntegrator->getCachedTF("tfHiggs1"))->GetXaxis()->GetXmin() * (1-enlargeEh1);
		  double Eh1high  = (meIntegrator->getCachedTF("tfHiggs1"))->GetXaxis()->GetXmax() * (1+enlargeEh1);	  

		  double xLmode0[4] = {  E1low,  Ptlow,   Philow,   Eh1low};
		  double xUmode0[4] = {  E1high, Pthigh , Phihigh,  Eh1high};

		  double xLmode1[6] = {  E1low,   -1, -TMath::Pi(), Ptlow,   Philow,   Eh1low};
		  double xUmode1[6] = {  E1high,  +1, +TMath::Pi(), Pthigh , Phihigh,  Eh1high};

		  ROOT::Math::Functor toIntegrate(meIntegrator, &MEIntegratorNew::Eval, par);
		  ROOT::Math::GSLMCIntegrator ig2( ROOT::Math::IntegrationMultiDim::kVEGAS , 1.e-12, 1.e-5, vegasPoints);
		  ig2.SetFunction(toIntegrate);
		  meIntegrator->SetPar(par);

		  double p = 0.;
		  if(mode==0) 
		    p = ig2.Integral(xLmode0, xUmode0);
		  else if(mode==1)
		    p = ig2.Integral(xLmode1, xUmode1);
		  else{ }
		
		  float old = hMassProb->GetBinContent( hMassProb->FindBin(mH[m]) );
		  hMassProb->SetBinContent( hMassProb->FindBin(mH[m]), old+p );
	
		  meIntegrator->deleteTF();	
		}

	      }
	    }
	  }
	}// mass points

	/*
	TF1* gaus = new TF1("gaus","gaus", 50,200);
	//gaus->SetParameter(1, hMassProb->GetBinCenter(hMassProb->GetMaximumBin()) );
	//gaus->SetParameter(2, 10. );
	//hBestMass->Fill( hMass->GetBinCenter(hMass->GetMaximumBin())  );
	hMassProb->Fit(gaus, "Q", "", 
		       hMassProb->GetBinCenter(hMassProb->GetMaximumBin()) - 30, 
		       hMassProb->GetBinCenter(hMassProb->GetMaximumBin()) + 30);
	float newMean  = gaus->GetParameter(1);
	float newSigma = gaus->GetParameter(2);
	hMassProb->Fit(gaus, "Q", "", 
		       newMean - 2*newSigma,
		       newMean + 2*newSigma
		       );
	newMean  = gaus->GetParameter(1);
	newSigma = gaus->GetParameter(2);
	hMassProb->Fit(gaus, "Q", "", 
		       newMean - 2*newSigma,
		       newMean + 2*newSigma
		       );
	*/

	float est = -99;
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
	  //A.Print("");
	  TMatrixD AInv(3,3);
	  double det;
	  AInv = A.Invert(&det);
	  //AInv.Print("");
	  
	  TMatrixD Y(3,1);
	  const double yPos[3] = 
	    { yD, yC, yU
	    };
	  Y.SetMatrixArray(yPos, "C");
	  //Y.Print("");
	  
	  TMatrixD C(3,1);
	  const double dummy[3] = 
	    { 1., 1., 1.
	    };
	  C.SetMatrixArray(dummy,"C");
	  C.Mult(AInv,Y);
	  //C.Print("");
	  
	  a = C(2,0);
	  b = C(1,0);
	  c = C(0,0);
	  //cout << "a = " << a << endl;
	  //cout << "b = " << b << endl;
	  //cout << "c = " << c << endl;
	  //cout << "Estim. = " << -b/2/a << endl; 
	  est = -b/2/a ;
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
	  //A.Print("");
	  TMatrixD AInv(3,3);
	  double det;
	  AInv = A.Invert(&det);
	  //AInv.Print("");
	  
	  TMatrixD Y(3,1);
	  const double yPos[3] = 
	    { yD, yC, yU
	    };
	  Y.SetMatrixArray(yPos, "C");
	  //Y.Print("");
	  
	  TMatrixD C(3,1);
	  const double dummy[3] = 
	    { 1., 1., 1.
	    };
	  C.SetMatrixArray(dummy,"C");
	  C.Mult(AInv,Y);
	  //C.Print("");
	  
	  float a = C(2,0);
	  float b = C(1,0);
	  float c = C(0,0);
	  //cout << "a = " << a << endl;
	  //cout << "b = " << b << endl;
	  //cout << "c = " << c << endl;
	  //cout << "Estim. = " << -b/2/a << endl; 
	  est = -b/2/a ;	  
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
	  //A.Print("");
	  TMatrixD AInv(3,3);
	  double det;
	  AInv = A.Invert(&det);
	  //AInv.Print("");
	  
	  TMatrixD Y(3,1);
	  const double yPos[3] = 
	    { yD, yC, yU
	    };
	  Y.SetMatrixArray(yPos, "C");
	  //Y.Print("");
	  
	  TMatrixD C(3,1);
	  const double dummy[3] = 
	    { 1., 1., 1.
	    };
	  C.SetMatrixArray(dummy,"C");
	  C.Mult(AInv,Y);
	  //C.Print("");
	  
	  a = C(2,0);
	  b = C(1,0);
	  //c = C(0,0);

	  //cout << "a = " << a << endl;
	  //cout << "b = " << b << endl;
	  //cout << "c = " << c << endl;
	  //cout << "Estim. = " << -b/2/a << endl; 
	  est = -b/2/a ;
	}
	else{
	  est =  hMassProb->GetBinCenter (hMassProb->GetMaximumBin());
	}



	hBestMassProb->Fill( est );
	//hBestMassProb->Fill( gaus->GetParameter(1)  );
        //hBestMassProb->Fill( hMassProb->GetBinCenter(hMassProb->GetMaximumBin())  );	
	//mass_ = hMassProb->GetBinCenter(hMassProb->GetMaximumBin());
	//mass_ = gaus->GetParameter(1);
	//prob_ = gaus->Eval( gaus->GetParameter(1) );
	mass_   = est;
	prob_   = a*est*est + b*est + c;
	tProb->Fill();
	//delete gaus;

	TFile* fout_tmp = TFile::Open(outFileName.c_str(),"UPDATE");
	TDirectory *dir = fout_tmp->mkdir(Form("Event_%d", counter));
	if(dir)
	  dir->cd();
	else
	  fout_tmp->cd(Form("Event_%d", counter));
	hMassProb->Write("hMassProb",TObject::kOverwrite);
	fout_tmp->Close();
	
      }
      //////////////////////////////////////////////////////////
    
    }
    
  }
  

  TFile* fout_tmp = TFile::Open(outFileName.c_str(),"UPDATE");
  hBestMass->Write("hBestMass",TObject::kOverwrite);
  hBestMassProb->Write("hBestMassProb",TObject::kOverwrite);
  tProb->Write("",TObject::kOverwrite );
  //hPdfLumi->Write("hPdfLumi",TObject::kOverwrite);
  //meSquared1->Write("meSquared1",TObject::kOverwrite);
  //meSquared2->Write("meSquared2",TObject::kOverwrite);
  //meSquared3->Write("meSquared3",TObject::kOverwrite);
  //meSquared4->Write("meSquared4",TObject::kOverwrite);

  fout_tmp->Close();

  //cout << "Closing file..." << endl;
  //fout->Close();
  //cout << "Deleting file..." << endl;
  //delete fout;
  
  cout << "Delete meIntegrator..." << endl;
  delete meIntegrator;
  delete clock; delete ran;
  cout << "Finished!!!" << endl;
  return 0;
  

}
