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

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GSLMCIntegrator.h"

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

#include "Bianchi/TTHStudies/interface/MEIntegrator.h"

#define GENJETDR 0.3
#define VERBOSE  false

using namespace std;
using namespace RooFit;



int main(int argc, const char* argv[])
{

  std::cout << "TestME" << std::endl;
  gROOT->SetBatch(true);
 
  gSystem->Load("libFWCoreFWLite");
  gSystem->Load("libDataFormatsFWLite");
  //gSystem->Load("BianchiTTHStudiesPlugins");

  AutoLibraryLoader::enable();

  TFile* fout = 0;

  PythonProcessDesc builder(argv[1]);
  const edm::ParameterSet& in = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteInput");
  std::string outFileName( in.getParameter<std::string>  ("outFileName" ) );
  std::string pathToFile(  in.getParameter<std::string> ("pathToFile"   ) );
  int vegasPoints(  in.getParameter<int> ("vegasPoints"   ) );

  int par = 8;
  MEIntegrator* meIntegrator = new MEIntegrator( pathToFile , par);
  meIntegrator->debug();

  vector<TLorentzVector> jets;
  TLorentzVector jet1,jet2,jet3,jet4,jet5,jet6,jet7;
  jet1.SetPtEtaPhiM(100,0.0, 0.0,10);
  jet2.SetPtEtaPhiM(70, 1.0,-1.0,10);
  jet3.SetPtEtaPhiM(70,-1.0,-2.0,10);
  jet4.SetPtEtaPhiM(65,-0.5,-1.0, 0);

  jet5.SetPtEtaPhiM(49.2, -2.14,  0.421, 0.1);
  jet6.SetPtEtaPhiM(47.8, -3.16, -0.967, 0.1);
  jet7.SetPtEtaPhiM(54.7, -3.4   ,2.96,    5);

  jets.push_back( jet1 );  jets.push_back( jet2 );  jets.push_back( jet3 );
  jets.push_back( jet4 );  jets.push_back( jet5 );  jets.push_back( jet6 );
  jets.push_back( jet7 );

  vector<float> bTagging;
  bTagging.push_back( 0.0 ) ; bTagging.push_back( 0.0 ) ; bTagging.push_back( 0.0 ) ; bTagging.push_back( 0.0 ) ; 
  bTagging.push_back( 0.0 ) ; bTagging.push_back( 0.0 ) ; bTagging.push_back( 0.0 ) ; 


  double xL[8] = { 0.00, -0.6,   0.,  -TMath::Pi(), -1.,  -TMath::Pi(), -1,  -TMath::Pi()};
  double xU[8] = { 0.50,  0.6,  600.,  TMath::Pi(),  1.,   TMath::Pi(),  1,   TMath::Pi()};

  ROOT::Math::GSLMCIntegrator ig2("vegas", 1.e-12, 1.e-5, vegasPoints);
  ROOT::Math::Functor toIntegrate(meIntegrator, &MEIntegrator::Eval, par); 
  meIntegrator->SetPar(par);
  meIntegrator->setJets(&jets);
  meIntegrator->setBtag(&bTagging);
  meIntegrator->createMash();
  meIntegrator->setInputLV( TLorentzVector( 20., 20.,0.,TMath::Sqrt(120*120+20*20*2)), 
			    TLorentzVector(-20., 20.,0.,TMath::Sqrt(170*170+20*20*2)) );
  ig2.SetFunction(toIntegrate);
  double p = ig2.Integral(xL, xU);
  cout << "Prob  = " << p << endl;

  fout = new TFile("TestME.root","UPDATE");
  fout->cd();
  (meIntegrator->getCachedPdf("pdfPtTTH"))->Write("pdfPtTTH",TObject::kOverwrite);
  (meIntegrator->getCachedPdf("pdfV2"))->Write("pdfV2",TObject::kOverwrite);
  (meIntegrator->getCachedPdf("pdfPtHStar"))->Write("pdfPtHStar",TObject::kOverwrite);
  (meIntegrator->getCachedPdf("pdfBetaW"))->Write("pdfBetaW",TObject::kOverwrite);
  (meIntegrator->getCachedPdf("pdfGammaW"))->Write("pdfGammaW",TObject::kOverwrite);
  (meIntegrator->getCachedPdf("pdfCsvHeavy"))->Write("pdfCsvHeavy",TObject::kOverwrite);
  (meIntegrator->getCachedPdf("pdfCsvLight"))->Write("pdfCsvLight",TObject::kOverwrite);
  (meIntegrator->getCachedPdf("pdfV1"))->Write("pdfV1",TObject::kOverwrite);
  (meIntegrator->getMash())->Write("mash",TObject::kOverwrite);
  (meIntegrator->getDebugHisto())->Write("debugHisto1",TObject::kOverwrite);
  fout->Close();
  delete fout;


  delete meIntegrator;
  return 0;
}
