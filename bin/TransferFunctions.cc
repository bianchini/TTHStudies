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
#include "RooChebychev.h"
#include "RooPolynomial.h"


#define GENJETDR 0.3
#define VERBOSE  false

using namespace std;
using namespace RooFit;



int main(int argc, const char* argv[])
{

  std::cout << "TransferFunctions" << std::endl;
  gROOT->SetBatch(true);
 
  gSystem->Load("libFWCoreFWLite");
  gSystem->Load("libDataFormatsFWLite");

  AutoLibraryLoader::enable();

  TFile* fout = 0;
  fout = new TFile("./root/ControlPlotsNew_new.root","RECREATE");
  fout->cd();


  PythonProcessDesc builder(argv[1]);
  const edm::ParameterSet& in = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteInput");
  std::string outFileName( in.getParameter<std::string>("outFileName" ) );
  std::string pathToFile(  in.getParameter<std::string> ("pathToFile"  ) );

  TFile* file = TFile::Open(pathToFile.c_str(),"READ");
  if(!file || file->IsZombie()){
    cout << "No input file..." << endl;
    return 0;
  }
  TTree* tree = (TTree*)file->Get("genTree");
  if(!tree){
    cout << "The tree called 'genTree' is not there..."<< endl;
    return 0;
  }


  TTree* treeEvent = (TTree*)file->Get("genEventTree");
  if(!treeEvent){
    cout << "The tree called 'genEventTree' is not there..."<< endl;
    return 0;
  }

  TTree* treeJetsLight = (TTree*)file->Get("genJetLightTree");
  if(!treeJetsLight){
    cout << "The tree called 'genJetLightTree' is not there..."<< endl;
    return 0;
  }
  TTree* treeJetsHeavy = (TTree*)file->Get("genJetHeavyTree");
  if(!treeJetsHeavy){
    cout << "The tree called 'genJetHeavyTree' is not there..."<< endl;
    return 0;
  }


  float etaBinning[] = {0.0, 1.0, 2.5};

  // ACCEPTANCE
  TH1F* accLightBin0  = new TH1F("accLightBin0","",100,0,400);
  TH1F* accLightBin1  = new TH1F("accLightBin1","",100,0,400);
  TH1F* accHeavyBin0  = new TH1F("accHeavyBin0","",100,0,400);
  TH1F* accHeavyBin1  = new TH1F("accHeavyBin1","",100,0,400);


  // RESOLUTION (GAUSSIAN)
  TH1F* resolLightBin0  = new TH1F("resolLightBin0","",58,30,300);
  TH1F* respLightBin0   = new TH1F("respLightBin0","", 58,30,300);
  TH1F* resolLightBin1  = new TH1F("resolLightBin1","",58,30,300);
  TH1F* respLightBin1   = new TH1F("respLightBin1","", 58,30,300);

  TH1F* resolHeavyBin0  = new TH1F("resolHeavyBin0","",58,30,300);
  TH1F* respHeavyBin0   = new TH1F("respHeavyBin0","", 58,30,300);
  TH1F* resolHeavyBin1  = new TH1F("resolHeavyBin1","",58,30,300);
  TH1F* respHeavyBin1   = new TH1F("respHeavyBin1","", 58,30,300);


  // RESOLUTION (DOUBLE GAUSSIAN)
  TH1F* respG1HeavyBin0     = new TH1F("respG1HeavyBin0","", 58,30,300);
  TH1F* respG2HeavyBin0     = new TH1F("respG2HeavyBin0","",   58,30,300);  
  TH1F* resolG1HeavyBin0    = new TH1F("resolG1HeavyBin0","",  58,30,300);
  TH1F* resolG2HeavyBin0    = new TH1F("resolG2HeavyBin0","",  58,30,300);
  TH1F* respG1HeavyBin1     = new TH1F("respG1HeavyBin1","", 58,30,300);
  TH1F* respG2HeavyBin1     = new TH1F("respG2HeavyBin1","",   58,30,300);  
  TH1F* resolG1HeavyBin1    = new TH1F("resolG1HeavyBin1","",  58,30,300);
  TH1F* resolG2HeavyBin1    = new TH1F("resolG2HeavyBin1","",  58,30,300);

  // RESOLUTION MASS (GAUSSIAN)
  //TH1F* resolMassLightBin0  = new TH1F("resolMassLightBin0","",58,30,300);
  //TH1F* respMassLightBin0   = new TH1F("respMassLightBin0","", 58,30,300);
  //TH1F* corrMassLightBin0   = new TH1F("corrMassLightBin0","", 58,30,300);
  //TH1F* resolMassLightBin1  = new TH1F("resolMassLightBin1","",58,30,300);
  //TH1F* respMassLightBin1   = new TH1F("respMassLightBin1","", 58,30,300);
  //TH1F* corrMassLightBin1   = new TH1F("corrMassLightBin1","", 58,30,300);

  //TH1F* resolMassHeavyBin0  = new TH1F("resolMassHeavyBin0","",58,30,300);
  //TH1F* respMassHeavyBin0   = new TH1F("respMassHeavyBin0","", 58,30,300);
  //TH1F* corrMassHeavyBin0   = new TH1F("corrMassHeavyBin0","", 58,30,300);
  //TH1F* resolMassHeavyBin1  = new TH1F("resolMassHeavyBin1","",58,30,300);
  //TH1F* respMassHeavyBin1   = new TH1F("respMassHeavyBin1","", 58,30,300);
  //TH1F* corrMassHeavyBin1   = new TH1F("corrMassHeavyBin1","", 58,30,300);

  //TH2F* hCorrMassLightBin0 = new TH2F("hCorrMassLightBin0","", 25, 30,530 , 30, 0,0.4);
  //TH2F* hCorrMassLightBin1 = new TH2F("hCorrMassLightBin1","", 25, 30,530 , 30, 0,0.4);
  //TH2F* hCorrMassHeavyBin0 = new TH2F("hCorrMassHeavyBin0","", 25, 30,530 , 30, 0,0.4);
  //TH2F* hCorrMassHeavyBin1 = new TH2F("hCorrMassHeavyBin1","", 25, 30,530 , 30, 0,0.4);


  for(int k = 0; k < 2; k++){

    TH1F* accPass   = (TH1F*)accLightBin0->Clone(Form("accPass%d",k)); //accPass->Sumw2();
    TH1F* accTot    = (TH1F*)accLightBin0->Clone(Form("accTot%d",k));  //accTot->Sumw2();

    TH2F* hCorrLight = new TH2F("hCorrLight","", 58, 30,300 , 100, 0,500);
    treeJetsLight->Draw("ptReco:pt>>hCorrLight", Form("(ptReco>30)*(TMath::Abs(eta)>%f && TMath::Abs(eta)<%f)",etaBinning[k],etaBinning[k+1]) );

//    if( k==0 ){
//       treeJetsLight->Draw("massReco/ptReco:pt>>hCorrMassLightBin0", Form("(ptReco>30)*(TMath::Abs(eta)>%f && TMath::Abs(eta)<%f)",etaBinning[k],etaBinning[k+1]) );
//       for( int b = 1; b<=hCorrMassLightBin0->GetNbinsX(); b++ ){
// 	 TH1F* h = (TH1F*)hCorrMassLightBin0->ProjectionY("_py",b,b);
// 	 float norm = h->Integral()*h->GetBinWidth(1);
// 	 for( int bb = 1; bb<=hCorrMassLightBin0->GetNbinsY(); bb++ ){
// 	   hCorrMassLightBin0->SetBinContent( b,bb, h->Integral()>0 ?  hCorrMassLightBin0->GetBinContent(b,bb)/norm : 0. );
// 	 }
//       }
//     }
//     else{
//       treeJetsLight->Draw("massReco/ptReco:pt>>hCorrMassLightBin1", Form("(ptReco>30)*(TMath::Abs(eta)>%f && TMath::Abs(eta)<%f)",etaBinning[k],etaBinning[k+1]) );
//       for( int b = 1; b<=hCorrMassLightBin1->GetNbinsX(); b++ ){
// 	 TH1F* h = (TH1F*)hCorrMassLightBin1->ProjectionY("_py",b,b);
// 	 float norm = h->Integral()*h->GetBinWidth(1);
// 	 for( int bb = 1; bb<=hCorrMassLightBin1->GetNbinsY(); bb++ ){
// 	   hCorrMassLightBin1->SetBinContent( b,bb,   h->Integral()>0 ? hCorrMassLightBin1->GetBinContent(b,bb)/norm : 0.);
// 	 }
//       }
//     }

    treeJetsLight->Draw(Form("pt>>accTot%d", k), Form("(TMath::Abs(eta)>%f && TMath::Abs(eta)<%f)",            etaBinning[k],etaBinning[k+1]) );
    treeJetsLight->Draw(Form("pt>>accPass%d",k), Form("(ptReco>30)*(TMath::Abs(eta)>%f && TMath::Abs(eta)<%f)",etaBinning[k],etaBinning[k+1]) );
    accPass->Sumw2();
    accTot->Sumw2();

    for(int j = 1; j<=hCorrLight->GetNbinsX(); j++){

      TH1F* hSlice = (TH1F*)hCorrLight->ProjectionY("_py",j,j);
      if(hSlice->GetEntries()<10) continue;
      cout << hSlice->GetMean() << ", " << hSlice->GetRMS() << endl;

      //TH1F* hSliceMass = (TH1F*)hCorrMassLight->ProjectionY("_py",j,j);
      //if(hSliceMass->GetEntries()<10) continue;
      //cout << hSliceMass->GetMean() << ", " << hSliceMass->GetRMS() << endl;

      TF1* resol = new TF1("resol","gaus", 0,500);
      hSlice->Fit(resol);

      //TF1* resolMass = new TF1("resolMass","gaus", 0,500);
      //hSliceMass->Fit(resolMass);

      fout->mkdir(Form("BinLight%d_%d",j,k));
      fout->cd( Form("BinLight%d_%d",j,k) );
      hSlice    ->Write("",TObject::kOverwrite);
      //hSliceMass->Write("",TObject::kOverwrite);
	   
     
      //TF1* resolDG = new TF1("resolDG","[0]*( 0.9*exp(-0.5*(x-[1])*(x-[1])/[2]/[2]) + (1-0.9)*exp(-0.5*(x-[1])*(x-[1])/[3]/[3]))", 0,500);
      //resolDG->SetParameter(1, hSlice->GetMean());
      //resolDG->SetParameter(2, hSlice->GetRMS());
      //resolDG->SetParameter(3, hSlice->GetRMS());
      //resolDG->SetParLimits(4,0.,1);
      //hSlice->Fit(resolDG);
      //hSlice->Fit(resolDG,"","",  hSlice->GetMean()-1.5*hSlice->GetRMS(), hSlice->GetMean()+1.5*hSlice->GetRMS());


      //TH2F* hCorrEnergyMassLight = new TH2F("hCorrEnergyMassLight","", 
      //				    20, resol->GetParameter(1)-2*resol->GetParameter(2),resol->GetParameter(1)+2*resol->GetParameter(2), //hSlice->GetMean()-hSlice->GetRMS(),hSlice->GetMean()+hSlice->GetRMS()  ,
      //				    20, resolMass->GetParameter(1)-2*resolMass->GetParameter(2),resolMass->GetParameter(1)+2*resolMass->GetParameter(2)//hSliceMass->GetMean()-hSliceMass->GetRMS(),hSliceMass->GetMean()+hSliceMass->GetRMS()
      //			    );
      //treeJetsLight->Draw("massReco:ptReco>>hCorrEnergyMassLight", 
      //		  Form("(ptReco>30)*(TMath::Abs(eta)>%f && TMath::Abs(eta)<%f && pt>%f && pt<%f)",etaBinning[k],etaBinning[k+1], 
      //		       hCorrLight->GetBinLowEdge(j),  
      //		       hCorrLight->GetBinLowEdge(j)+hCorrLight->GetBinWidth(j) ) );      
      //hCorrEnergyMassLight->Write("",TObject::kOverwrite);

     
      if(k==0){
	resolLightBin0->SetBinContent(j, resol->GetParameter(2));
	resolLightBin0->SetBinError  (j, resol->GetParError (2));
	respLightBin0->SetBinContent (j, resol->GetParameter(1));
	respLightBin0->SetBinError   (j, resol->GetParError (1));

	//resolMassLightBin0->SetBinContent(j, resolMass->GetParameter(2));
	//resolMassLightBin0->SetBinError  (j, resolMass->GetParError (2));
	//respMassLightBin0->SetBinContent (j, resolMass->GetParameter(1));
	//respMassLightBin0->SetBinError   (j, resolMass->GetParError (1));
	//corrMassLightBin0->SetBinContent (j, hCorrEnergyMassLight->GetCorrelationFactor());

	//fracDGLightBin0->SetBinContent (j,  resolDG->GetParameter(4) );
	//fracDGLightBin0->SetBinError   (j,  resolDG->GetParError (4) );
	//respDGLightBin0->SetBinContent (j,  resolDG->GetParameter(1) );
	//respDGLightBin0->SetBinError   (j,  resolDG->GetParError (1) );
	//resolG1LightBin0->SetBinContent(j,  resolDG->GetParameter(2) );
	//resolG1LightBin0->SetBinError  (j,  resolDG->GetParError (2) );
	//resolG2LightBin0->SetBinContent(j,  resolDG->GetParameter(3) );
	//resolG2LightBin0->SetBinError  (j,  resolDG->GetParError (3) );
      }
      if(k==1){
	resolLightBin1->SetBinContent(j, resol->GetParameter(2));
	resolLightBin1->SetBinError  (j, resol->GetParError(2));
	respLightBin1->SetBinContent (j, resol->GetParameter(1));
	respLightBin1->SetBinError   (j, resol->GetParError(1));

	//resolMassLightBin1->SetBinContent(j, resolMass->GetParameter(2));
	//resolMassLightBin1->SetBinError  (j, resolMass->GetParError (2));
	//respMassLightBin1->SetBinContent (j, resolMass->GetParameter(1));
	//respMassLightBin1->SetBinError   (j, resolMass->GetParError (1));
	//corrMassLightBin1->SetBinContent (j, hCorrEnergyMassLight->GetCorrelationFactor());

	//fracDGLightBin1->SetBinContent (j,  resolDG->GetParameter(4) );
	//fracDGLightBin1->SetBinError   (j,  resolDG->GetParError (4) );
	//respDGLightBin1->SetBinContent (j,  resolDG->GetParameter(1) );
	//respDGLightBin1->SetBinError   (j,  resolDG->GetParError (1) );
	//resolG1LightBin1->SetBinContent(j,  resolDG->GetParameter(2) );
	//resolG1LightBin1->SetBinError  (j,  resolDG->GetParError (2) );
	//resolG2LightBin1->SetBinContent(j,  resolDG->GetParameter(3) );
	//resolG2LightBin1->SetBinError  (j,  resolDG->GetParError (3) );
  
      }

   
     
      delete resol; //delete resolDG; 
      //delete resolMass; 
    }


    if(k==0)
      accLightBin0->Divide( accPass,accTot, 1.0,1.0);
    if(k==1)
      accLightBin1->Divide( accPass,accTot, 1.0,1.0);

    accPass->Reset();
    accTot->Reset();

    TH2F* hCorrHeavy = new TH2F("hCorrHeavy","",  58, 30,300 , 100, 0,500);
    treeJetsHeavy->Draw("ptReco:pt>>hCorrHeavy", Form("(ptReco>30)*(TMath::Abs(eta)>%f && TMath::Abs(eta)<%f)",etaBinning[k],etaBinning[k+1]) );

    //TH2F* hCorrMassHeavy = new TH2F("hCorrMassHeavy","", 58, 30,300 , 40, 0,80);
    //treeJetsHeavy->Draw("massReco:pt>>hCorrMassHeavy", Form("(ptReco>30)*(TMath::Abs(eta)>%f && TMath::Abs(eta)<%f)",etaBinning[k],etaBinning[k+1]) );

    treeJetsHeavy->Draw(Form("pt>>accTot%d",k),  Form("(TMath::Abs(eta)>%f && TMath::Abs(eta)<%f)",            etaBinning[k],etaBinning[k+1]) );
    treeJetsHeavy->Draw(Form("pt>>accPass%d",k), Form("(ptReco>30)*(TMath::Abs(eta)>%f && TMath::Abs(eta)<%f)",etaBinning[k],etaBinning[k+1]) );
    accPass->Sumw2();
    accTot->Sumw2();


//     if( k==0 ){
//       treeJetsHeavy->Draw("massReco/ptReco:pt>>hCorrMassHeavyBin0", Form("(ptReco>30)*(TMath::Abs(eta)>%f && TMath::Abs(eta)<%f)",etaBinning[k],etaBinning[k+1]) );
//       for( int b = 1; b<=hCorrMassHeavyBin0->GetNbinsX(); b++ ){
// 	TH1F* h = (TH1F*)hCorrMassHeavyBin0->ProjectionY("_py",b,b);
// 	float norm = h->Integral()*h->GetBinWidth(1);
// 	for( int bb = 1; bb<=hCorrMassHeavyBin0->GetNbinsY(); bb++ ){
// 	  hCorrMassHeavyBin0->SetBinContent( b,bb, h->Integral()>0 ?  hCorrMassHeavyBin0->GetBinContent(b,bb)/norm : 0. );
// 	}
//       }
//     }
//     else{
//       treeJetsHeavy->Draw("massReco/ptReco:pt>>hCorrMassHeavyBin1", Form("(ptReco>30)*(TMath::Abs(eta)>%f && TMath::Abs(eta)<%f)",etaBinning[k],etaBinning[k+1]) );
//       for( int b = 1; b<=hCorrMassHeavyBin1->GetNbinsX(); b++ ){
// 	TH1F* h = (TH1F*)hCorrMassHeavyBin1->ProjectionY("_py",b,b);
// 	float norm = h->Integral()*h->GetBinWidth(1);
// 	for( int bb = 1; bb<=hCorrMassHeavyBin1->GetNbinsY(); bb++ ){
// 	  hCorrMassHeavyBin1->SetBinContent( b,bb, h->Integral()>0 ? hCorrMassHeavyBin1->GetBinContent(b,bb)/norm : 0.);
// 	}
//       }
//     }


    for(int j = 1; j<=hCorrHeavy->GetNbinsX(); j++){

      TH1F* hSlice = (TH1F*)hCorrHeavy ->ProjectionY("_py",j,j);
      if(hSlice->GetEntries()<8) continue;
      cout << hSlice->GetMean() << ", " << hSlice->GetRMS() << endl;

      //TH1F* hSliceMass = (TH1F*)hCorrMassHeavy->ProjectionY("_py",j,j);
      //if(hSliceMass->GetEntries()<10) continue;
      //cout << hSliceMass->GetMean() << ", " << hSliceMass->GetRMS() << endl;

      TF1* resol = new TF1("resol","gaus", 0,500);
      hSlice->Fit(resol);

      //TF1* resolMass = new TF1("resolMass","gaus", 0,500);
      //hSliceMass->Fit(resolMass);

      fout->mkdir(Form("BinHeavy%d_%d",j,k));
      fout->cd( Form("BinHeavy%d_%d",j, k) );
      hSlice->Write("",TObject::kOverwrite);
      //hSliceMass->Write("",TObject::kOverwrite);

      TF1* resolDG = 0;
      if(k==1){
	resolDG = new TF1("resolDG","[0]*( 0.65*exp(-0.5*(x-[1])*(x-[1])/[2]/[2]) + (1-0.65)*exp(-0.5*(x-[4])*(x-[4])/[3]/[3]))", 0,500);
	resolDG->SetParameter(1, hSlice->GetMean());
	resolDG->SetParameter(2, hSlice->GetRMS());
	resolDG->SetParameter(3, hSlice->GetRMS());
	resolDG->SetParLimits(4, hSlice->GetMean()-2*hSlice->GetRMS() ,hSlice->GetMean());  
      }
      else{	
	resolDG = new TF1("resolDG","[0]*( 0.65*exp(-0.5*(x-[1])*(x-[1])/[2]/[2]) + (1-0.65)*exp(-0.5*(x-[1]*0.92)*(x-[1]*0.92)/[3]/[3]))", 0,500);
	resolDG->SetParameter(1, hSlice->GetMean());
	resolDG->SetParameter(2, hSlice->GetRMS());
	resolDG->SetParameter(3, hSlice->GetRMS());
      }


      //resolDG->SetParLimits(4,0.,1);
      //resolDG->SetParameter(4, 0.9 );

      hSlice->Fit(resolDG);
      hSlice->Fit(resolDG);



      hSlice->Write("",TObject::kOverwrite);


      //TH2F* hCorrEnergyMassHeavy = new TH2F("hCorrEnergyMassHeavy","", 
      //				    20, resol->GetParameter(1)-2*resol->GetParameter(2),resol->GetParameter(1)+2*resol->GetParameter(2), //hSlice->GetMean()-hSlice->GetRMS(),hSlice->GetMean()+hSlice->GetRMS()  ,
      //				    20, resolMass->GetParameter(1)-2*resolMass->GetParameter(2),resolMass->GetParameter(1)+2*resolMass->GetParameter(2)//hSliceMass->GetMean()-hSliceMass->GetRMS(),hSliceMass->GetMean()+hSliceMass->GetRMS()
      //			    );
      //treeJetsHeavy->Draw("massReco:ptReco>>hCorrEnergyMassHeavy", 
      //		  Form("(ptReco>30)*(TMath::Abs(eta)>%f && TMath::Abs(eta)<%f && pt>%f && pt<%f)",etaBinning[k],etaBinning[k+1], 
      //		       hCorrHeavy->GetBinLowEdge(j),  
      //		       hCorrHeavy->GetBinLowEdge(j)+hCorrHeavy->GetBinWidth(j) ) );
      //hCorrEnergyMassHeavy->Write("",TObject::kOverwrite);
	  


      if(k==0){

	resolHeavyBin0->SetBinContent(j, resol->GetParameter(2));
	resolHeavyBin0->SetBinError  (j, resol->GetParError (2));
	respHeavyBin0->SetBinContent (j, resol->GetParameter(1));
	respHeavyBin0->SetBinError   (j, resol->GetParError (1));

	//resolMassHeavyBin0->SetBinContent(j, resolMass->GetParameter(2));
	//resolMassHeavyBin0->SetBinError  (j, resolMass->GetParError (2));
	//respMassHeavyBin0->SetBinContent (j, resolMass->GetParameter(1));
	//respMassHeavyBin0->SetBinError   (j, resolMass->GetParError (1));
	//corrMassHeavyBin0->SetBinContent (j, hCorrEnergyMassHeavy->GetCorrelationFactor());

	respG1HeavyBin0->SetBinContent   (j,  resolDG->GetParameter(1) );
	respG1HeavyBin0->SetBinError     (j,  resolDG->GetParError (1) );
	respG2HeavyBin0->SetBinContent   (j,  0.92* resolDG->GetParameter(1) );
	respG2HeavyBin0->SetBinError     (j,  0.92* resolDG->GetParError (1) );
	resolG1HeavyBin0->SetBinContent  (j,  resolDG->GetParameter(2) );
	resolG1HeavyBin0->SetBinError    (j,  resolDG->GetParError (2) );
	resolG2HeavyBin0->SetBinContent  (j,  resolDG->GetParameter(3) ); 
	resolG2HeavyBin0->SetBinError    (j,  resolDG->GetParError (3) );

      }
      if(k==1){
	resolHeavyBin1->SetBinContent(j, resol->GetParameter(2));
	resolHeavyBin1->SetBinError  (j, resol->GetParError(2));
	respHeavyBin1->SetBinContent (j, resol->GetParameter(1));
	respHeavyBin1->SetBinError   (j, resol->GetParError(1));

	//resolMassHeavyBin1->SetBinContent(j, resolMass->GetParameter(2));
	//resolMassHeavyBin1->SetBinError  (j, resolMass->GetParError (2));
	//respMassHeavyBin1->SetBinContent (j, resolMass->GetParameter(1));
	//respMassHeavyBin1->SetBinError   (j, resolMass->GetParError (1));
	//corrMassHeavyBin1->SetBinContent (j, hCorrEnergyMassHeavy->GetCorrelationFactor());

	respG1HeavyBin1->SetBinContent   (j,  resolDG->GetParameter(1) );
	respG1HeavyBin1->SetBinError     (j,  resolDG->GetParError (1) );
	respG2HeavyBin1->SetBinContent   (j,  resolDG->GetParameter(4) );
	respG2HeavyBin1->SetBinError     (j,  resolDG->GetParError (4) );
	resolG1HeavyBin1->SetBinContent  (j,  resolDG->GetParameter(2) );
	resolG1HeavyBin1->SetBinError    (j,  resolDG->GetParError (2) );
	resolG2HeavyBin1->SetBinContent  (j,  resolDG->GetParameter(3) );
	resolG2HeavyBin1->SetBinError    (j,  resolDG->GetParError (3) );

      }


    
      delete resol; //delete resolMass; 
      delete resolDG;
    }


    if(k==0)
      accHeavyBin0->Divide( accPass,accTot, 1.0,1.0);
    if(k==1)
      accHeavyBin1->Divide( accPass,accTot, 1.0,1.0);
      
    delete hCorrLight;  delete hCorrHeavy;
  }

  TF1* turnOnLightBin0 = new TF1("turnOnBin0","TMath::Erf([0]*x+[1])+[2]", 0,400);
  TF1* turnOnLightBin1 = new TF1("turnOnBin1","TMath::Erf([0]*x+[1])+[2]", 0,400);
  TF1* turnOnHeavyBin0 = new TF1("turnOnBin0","TMath::Erf([0]*x+[1])+[2]", 0,400);
  TF1* turnOnHeavyBin1 = new TF1("turnOnBin1","TMath::Erf([0]*x+[1])+[2]", 0,400);
  accLightBin0->Fit(turnOnLightBin0);  
  accLightBin1->Fit(turnOnLightBin1);
  accHeavyBin0->Fit(turnOnHeavyBin0);  
  accHeavyBin1->Fit(turnOnHeavyBin1);

  TF1* sigmaLightBin0 = new TF1("sigmaLightBin0","sqrt([0]*[0] + x*[1]*[1] + x*x*[2]*[2])",0,400);
  TF1* sigmaLightBin1 = new TF1("sigmaLightBin1","sqrt([0]*[0] + x*[1]*[1] + x*x*[2]*[2])",0,400);
  TF1* sigmaHeavyBin0 = new TF1("sigmaHeavyBin0","sqrt([0]*[0] + x*[1]*[1] + x*x*[2]*[2])",0,400);
  TF1* sigmaHeavyBin1 = new TF1("sigmaHeavyBin1","sqrt([0]*[0] + x*[1]*[1] + x*x*[2]*[2])",0,400);
 
  TF1* meanLightBin0  = new TF1("meanLightBin0", "x*[0]",0,400);
  TF1* meanLightBin1  = new TF1("meanLightBin1", "x*[0]",0,400);
  TF1* meanHeavyBin0  = new TF1("meanHeavyBin0", "x*[0]",0,400);
  TF1* meanHeavyBin1  = new TF1("meanHeavyBin1", "x*[0]",0,400);

  //TF1* sigmaMassLightBin0 = new TF1("sigmaMassLightBin0","sqrt([0]*[0] + x*[1]*[1] + x*x*[2]*[2])",0,400);
  //TF1* sigmaMassLightBin1 = new TF1("sigmaMassLightBin1","sqrt([0]*[0] + x*[1]*[1] + x*x*[2]*[2])",0,400);
  //TF1* sigmaMassHeavyBin0 = new TF1("sigmaMassHeavyBin0","sqrt([0]*[0] + x*[1]*[1] + x*x*[2]*[2])",0,400);
  //TF1* sigmaMassHeavyBin1 = new TF1("sigmaMassHeavyBin1","sqrt([0]*[0] + x*[1]*[1] + x*x*[2]*[2])",0,400);
 
  //TF1* meanMassLightBin0  = new TF1("meanMassLightBin0", "x*[0]+[1]",0,400);
  //TF1* meanMassLightBin1  = new TF1("meanMassLightBin1", "x*[0]+[1]",0,400);
  //TF1* meanMassHeavyBin0  = new TF1("meanMassHeavyBin0", "x*[0]+[1]",0,400);
  //TF1* meanMassHeavyBin1  = new TF1("meanMassHeavyBin1", "x*[0]+[1]",0,400);


  TF1* meanG1HeavyBin0  = new TF1("meanG1HeavyBin0", "x*[0]+[1]",0,400);
  TF1* meanG2HeavyBin0  = new TF1("meanG2HeavyBin0", "x*[0]+[1]",0,400);
  TF1* sigmaG1HeavyBin0 = new TF1("sigmaG1HeavyBin0","sqrt([0]*[0] + x*[1]*[1] + x*x*[2]*[2])",0,400);
  TF1* sigmaG2HeavyBin0 = new TF1("sigmaG2HeavyBin0","sqrt([0]*[0] + x*[1]*[1] + x*x*[2]*[2])",0,400);

  TF1* meanG1HeavyBin1  = new TF1("meanG1HeavyBin1", "x*[0]+[1]",0,400);
  TF1* meanG2HeavyBin1  = new TF1("meanG2HeavyBin1", "x*[0]+[1]",0,400);
  TF1* sigmaG1HeavyBin1 = new TF1("sigmaG1HeavyBin1","sqrt([0]*[0] + x*[1]*[1] + x*x*[2]*[2])",0,400);
  TF1* sigmaG2HeavyBin1 = new TF1("sigmaG2HeavyBin1","sqrt([0]*[0] + x*[1]*[1] + x*x*[2]*[2])",0,400);



  resolLightBin0->Fit(sigmaLightBin0,"","",60,160);
  resolLightBin1->Fit(sigmaLightBin1,"","",60,160);
  resolHeavyBin0->Fit(sigmaHeavyBin0,"","",60,160);
  resolHeavyBin1->Fit(sigmaHeavyBin1,"","",60,160);
  respLightBin0->Fit (meanLightBin0, "","",60,160);
  respLightBin1->Fit (meanLightBin1, "","",60,160);
  respHeavyBin0->Fit (meanHeavyBin0, "","",60,160);
  respHeavyBin1->Fit (meanHeavyBin1, "","",60,160);

  resolLightBin0->Fit(sigmaLightBin0,"","",40,200);
  resolLightBin1->Fit(sigmaLightBin1,"","",40,200);
  resolHeavyBin0->Fit(sigmaHeavyBin0,"","",40,200);
  resolHeavyBin1->Fit(sigmaHeavyBin1,"","",40,200);
  respLightBin0->Fit (meanLightBin0, "","",40,200);
  respLightBin1->Fit (meanLightBin1, "","",40,200);
  respHeavyBin0->Fit (meanHeavyBin0, "","",40,200);
  respHeavyBin1->Fit (meanHeavyBin1, "","",40,200);

  //resolMassLightBin0->Fit(sigmaMassLightBin0,"","",60,160);
  //resolMassLightBin1->Fit(sigmaMassLightBin1,"","",60,160);
  //resolMassHeavyBin0->Fit(sigmaMassHeavyBin0,"","",60,160);
  //resolMassHeavyBin1->Fit(sigmaMassHeavyBin1,"","",60,160);
  //respMassLightBin0->Fit (meanMassLightBin0, "","",60,160);
  //respMassLightBin1->Fit (meanMassLightBin1, "","",60,160);
  //respMassHeavyBin0->Fit (meanMassHeavyBin0, "","",60,160);
  //respMassHeavyBin1->Fit (meanMassHeavyBin1, "","",60,160);

  //resolMassLightBin0->Fit(sigmaMassLightBin0,"","",40,200);
  //resolMassLightBin1->Fit(sigmaMassLightBin1,"","",40,200);
  //resolMassHeavyBin0->Fit(sigmaMassHeavyBin0,"","",40,200);
  //resolMassHeavyBin1->Fit(sigmaMassHeavyBin1,"","",40,200);
  //respMassLightBin0->Fit (meanMassLightBin0, "","",40,200);
  //respMassLightBin1->Fit (meanMassLightBin1, "","",40,200);
  //respMassHeavyBin0->Fit (meanMassHeavyBin0, "","",40,200);
  //respMassHeavyBin1->Fit (meanMassHeavyBin1, "","",40,200);


  respG1HeavyBin0->Fit  (meanG1HeavyBin0, "","",60,160);
  respG2HeavyBin0->Fit  (meanG2HeavyBin0, "","",60,160);
  resolG1HeavyBin0->Fit (sigmaG1HeavyBin0,"","",60,160);
  resolG2HeavyBin0->Fit (sigmaG2HeavyBin0,"","",60,160);
  respG1HeavyBin1->Fit  (meanG1HeavyBin1, "","",80,200);
  respG2HeavyBin1->Fit  (meanG2HeavyBin1, "","",80,200);
  resolG1HeavyBin1->Fit (sigmaG1HeavyBin1,"","",80,200);
  resolG2HeavyBin1->Fit (sigmaG2HeavyBin1,"","",80,200);

  respG1HeavyBin0->Fit  (meanG1HeavyBin0, "","",40,200);
  respG2HeavyBin0->Fit  (meanG2HeavyBin0, "","",40,200);
  resolG1HeavyBin0->Fit (sigmaG1HeavyBin0,"","",40,200);
  resolG2HeavyBin0->Fit (sigmaG2HeavyBin0,"","",40,200);
  respG1HeavyBin1->Fit  (meanG1HeavyBin1, "","",60,300);
  respG2HeavyBin1->Fit  (meanG2HeavyBin1, "","",60,300);
  resolG1HeavyBin1->Fit (sigmaG1HeavyBin1,"","",60,300);
  resolG2HeavyBin1->Fit (sigmaG2HeavyBin1,"","",60,300);



  /////////////////////

  for(int k = 0; k < 2; k++){

    TH2F* hCorrHeavy = new TH2F("hCorrHeavy","",  58, 30,300 , 100, 0,500);
    for(int j = 1; j<=hCorrHeavy->GetNbinsX(); j++){

      float en = hCorrHeavy->GetBinCenter(j);

      //fout->cd( Form("BinHeavy%d_%d",j, k) );
      TH1F* hCorrHeavy_py = (TH1F*)fout->Get(Form("BinHeavy%d_%d/hCorrHeavy_py",j, k));

      if(!hCorrHeavy_py){
	cout << "No histo" << endl;
	continue;
      }

      TF1* resolDG_old =  hCorrHeavy_py->GetFunction("resolDG");
      if(!resolDG_old){
	cout << "No func" << endl;
	continue;
      }
      TF1* resolDG = 0;
      if(k==1){
	resolDG = new TF1("resolDG_fit","[0]*( 0.65*exp(-0.5*(x-[1])*(x-[1])/[2]/[2]) + (1-0.65)*exp(-0.5*(x-[4])*(x-[4])/[3]/[3]))", 0,500);
	resolDG->SetLineColor(kBlue);
	//hCorrHeavy_py->Fit(resolDG);
	float norm = resolDG_old->Integral(0,500);
	float m1 = meanG1HeavyBin1->Eval(en);
	float m2 = meanG2HeavyBin1->Eval(en);
	float s1 = sigmaG1HeavyBin1->Eval(en);
	float s2 = sigmaG2HeavyBin1->Eval(en);
	resolDG->SetParameter(0, norm/(0.65*sqrt(2*TMath::Pi())*s1 + (1-0.65)*sqrt(2*TMath::Pi())*s2 )  );
	resolDG->SetParameter(1, m1);
	resolDG->SetParameter(2, s1);
	resolDG->SetParameter(3, s2);
	resolDG->SetParameter(4, m2);
	//hCorrHeavy_py->Write("",TObject::kOverwrite);
	fout->cd( Form("BinHeavy%d_%d",j, k) );
	resolDG_old->Write(Form("resolDG%d_%d",    j,k),TObject::kOverwrite);
	resolDG    ->Write(Form("resolDG_fit%d_%d",j,k),TObject::kOverwrite);

      }
      else{

	resolDG = new TF1("resolDG_fit","[0]*( 0.65*exp(-0.5*(x-[1])*(x-[1])/[2]/[2]) + (1-0.65)*exp(-0.5*(x-0.92*[1])*(x-0.92*[1])/[3]/[3]))", 0,500);
	resolDG->SetLineColor(kBlue);
	//hCorrHeavy_py->Fit(resolDG);
	float norm = resolDG_old->Integral(0,500);
	float m1   = meanG1HeavyBin0->Eval(en);
	float m2   = meanG2HeavyBin0->Eval(en);
	float s1   = sigmaG1HeavyBin0->Eval(en);
	float s2   = sigmaG2HeavyBin0->Eval(en);
	resolDG->SetParameter(0, norm/(0.65*sqrt(2*TMath::Pi())*s1 + (1-0.65)*sqrt(2*TMath::Pi())*s2 )  );
	resolDG->SetParameter(1, m1);
	resolDG->SetParameter(2, s1);
	resolDG->SetParameter(3, s2);
	resolDG->SetParameter(4, m2);
	//hCorrHeavy_py->Write("",TObject::kOverwrite);

	fout->cd( Form("BinHeavy%d_%d",j, k) );
	resolDG_old->Write(Form("resolDG%d_%d",    j,k),TObject::kOverwrite);
	resolDG    ->Write(Form("resolDG_fit%d_%d",j,k),TObject::kOverwrite);

      }

    }

  }


  //////////////////////

 
  gSystem->cd("./root/");
  RooWorkspace *w = new RooWorkspace("transferFuntions","all functions for SL ttH");

  RooRealVar param0AccLightBin0("param0AccLightBin0","", turnOnLightBin0->GetParameter(0));
  RooRealVar param1AccLightBin0("param1AccLightBin0","", turnOnLightBin0->GetParameter(1));
  RooRealVar param2AccLightBin0("param2AccLightBin0","", turnOnLightBin0->GetParameter(2));
  RooRealVar param0AccLightBin1("param0AccLightBin1","", turnOnLightBin1->GetParameter(0));
  RooRealVar param1AccLightBin1("param1AccLightBin1","", turnOnLightBin1->GetParameter(1));
  RooRealVar param2AccLightBin1("param2AccLightBin1","", turnOnLightBin1->GetParameter(2));
  RooRealVar param0AccHeavyBin0("param0AccHeavyBin0","", turnOnHeavyBin0->GetParameter(0));
  RooRealVar param1AccHeavyBin0("param1AccHeavyBin0","", turnOnHeavyBin0->GetParameter(1));
  RooRealVar param2AccHeavyBin0("param2AccHeavyBin0","", turnOnHeavyBin0->GetParameter(2));
  RooRealVar param0AccHeavyBin1("param0AccHeavyBin1","", turnOnHeavyBin1->GetParameter(0));
  RooRealVar param1AccHeavyBin1("param1AccHeavyBin1","", turnOnHeavyBin1->GetParameter(1));
  RooRealVar param2AccHeavyBin1("param2AccHeavyBin1","", turnOnHeavyBin1->GetParameter(2));

  RooRealVar param0resolLightBin0("param0resolLightBin0","", sigmaLightBin0->GetParameter(0));
  RooRealVar param1resolLightBin0("param1resolLightBin0","", sigmaLightBin0->GetParameter(1));
  RooRealVar param2resolLightBin0("param2resolLightBin0","", sigmaLightBin0->GetParameter(2));
  RooRealVar param0resolLightBin1("param0resolLightBin1","", sigmaLightBin1->GetParameter(0));
  RooRealVar param1resolLightBin1("param1resolLightBin1","", sigmaLightBin1->GetParameter(1));
  RooRealVar param2resolLightBin1("param2resolLightBin1","", sigmaLightBin1->GetParameter(2));
  RooRealVar param0respLightBin0 ("param0respLightBin0","",  meanLightBin0->GetParameter(0));
  RooRealVar param1respLightBin0 ("param1respLightBin0","",  0.0);
  RooRealVar param0respLightBin1 ("param0respLightBin1","",  meanLightBin1->GetParameter(0));
  RooRealVar param1respLightBin1 ("param1respLightBin1","",  0.0);

  //RooRealVar param0resolMassLightBin0("param0resolMassLightBin0","", sigmaMassLightBin0->GetParameter(0));
  //RooRealVar param1resolMassLightBin0("param1resolMassLightBin0","", sigmaMassLightBin0->GetParameter(1));
  //RooRealVar param2resolMassLightBin0("param2resolMassLightBin0","", sigmaMassLightBin0->GetParameter(2));
  //RooRealVar param0resolMassLightBin1("param0resolMassLightBin1","", sigmaMassLightBin1->GetParameter(0));
  //RooRealVar param1resolMassLightBin1("param1resolMassLightBin1","", sigmaMassLightBin1->GetParameter(1));
  //RooRealVar param2resolMassLightBin1("param2resolMassLightBin1","", sigmaMassLightBin1->GetParameter(2));
  //RooRealVar param0respMassLightBin0 ("param0respMassLightBin0","",  meanMassLightBin0->GetParameter(0));
  //RooRealVar param1respMassLightBin0 ("param1respMassLightBin0","",  meanMassLightBin0->GetParameter(1));
  //RooRealVar param0respMassLightBin1 ("param0respMassLightBin1","",  meanMassLightBin1->GetParameter(0));
  //RooRealVar param1respMassLightBin1 ("param1respMassLightBin1","",  meanMassLightBin1->GetParameter(1));


  RooRealVar param0resolHeavyBin0("param0resolHeavyBin0","", sigmaHeavyBin0->GetParameter(0));
  RooRealVar param1resolHeavyBin0("param1resolHeavyBin0","", sigmaHeavyBin0->GetParameter(1));
  RooRealVar param2resolHeavyBin0("param2resolHeavyBin0","", sigmaHeavyBin0->GetParameter(2));
  RooRealVar param0resolHeavyBin1("param0resolHeavyBin1","", sigmaHeavyBin1->GetParameter(0));
  RooRealVar param1resolHeavyBin1("param1resolHeavyBin1","", sigmaHeavyBin1->GetParameter(1));
  RooRealVar param2resolHeavyBin1("param2resolHeavyBin1","", sigmaHeavyBin1->GetParameter(2));
  RooRealVar param0respHeavyBin0 ("param0respHeavyBin0","",  meanHeavyBin0->GetParameter(0));
  RooRealVar param1respHeavyBin0 ("param1respHeavyBin0","", 0.0);
  RooRealVar param0respHeavyBin1 ("param0respHeavyBin1","",  meanHeavyBin1->GetParameter(0));
  RooRealVar param1respHeavyBin1 ("param1respHeavyBin1","", 0.0);

  //RooRealVar param0resolMassHeavyBin0("param0resolMassHeavyBin0","", sigmaMassHeavyBin0->GetParameter(0));
  //RooRealVar param1resolMassHeavyBin0("param1resolMassHeavyBin0","", sigmaMassHeavyBin0->GetParameter(1));
  //RooRealVar param2resolMassHeavyBin0("param2resolMassHeavyBin0","", sigmaMassHeavyBin0->GetParameter(2));
  //RooRealVar param0resolMassHeavyBin1("param0resolMassHeavyBin1","", sigmaMassHeavyBin1->GetParameter(0));
  //RooRealVar param1resolMassHeavyBin1("param1resolMassHeavyBin1","", sigmaMassHeavyBin1->GetParameter(1));
  //RooRealVar param2resolMassHeavyBin1("param2resolMassHeavyBin1","", sigmaMassHeavyBin1->GetParameter(2));
  //RooRealVar param0respMassHeavyBin0 ("param0respMassHeavyBin0","",  meanMassHeavyBin0->GetParameter(0));
  //RooRealVar param1respMassHeavyBin0 ("param1respMassHeavyBin0","",  meanMassHeavyBin0->GetParameter(1));
  //RooRealVar param0respMassHeavyBin1 ("param0respMassHeavyBin1","",  meanMassHeavyBin1->GetParameter(0));
  //RooRealVar param1respMassHeavyBin1 ("param1respMassHeavyBin1","",  meanMassHeavyBin1->GetParameter(1));



  RooRealVar param0resolG1HeavyBin0("param0resolG1HeavyBin0","", sigmaG1HeavyBin0->GetParameter(0));
  RooRealVar param1resolG1HeavyBin0("param1resolG1HeavyBin0","", sigmaG1HeavyBin0->GetParameter(1));
  RooRealVar param2resolG1HeavyBin0("param2resolG1HeavyBin0","", sigmaG1HeavyBin0->GetParameter(2));
  RooRealVar param0resolG2HeavyBin0("param0resolG2HeavyBin0","", sigmaG2HeavyBin0->GetParameter(0));
  RooRealVar param1resolG2HeavyBin0("param1resolG2HeavyBin0","", sigmaG2HeavyBin0->GetParameter(1));
  RooRealVar param2resolG2HeavyBin0("param2resolG2HeavyBin0","", sigmaG2HeavyBin0->GetParameter(2));
  RooRealVar param0respG1HeavyBin0 ("param0respG1HeavyBin0","",  meanG1HeavyBin0->GetParameter(0));
  RooRealVar param1respG1HeavyBin0 ("param1respG1HeavyBin0","",  meanG1HeavyBin0->GetParameter(1));
  RooRealVar param0respG2HeavyBin0 ("param0respG2HeavyBin0","",  meanG2HeavyBin0->GetParameter(0));
  RooRealVar param1respG2HeavyBin0 ("param1respG2HeavyBin0","",  meanG2HeavyBin0->GetParameter(1));

  RooRealVar param0resolG1HeavyBin1("param0resolG1HeavyBin1","", sigmaG1HeavyBin1->GetParameter(0));
  RooRealVar param1resolG1HeavyBin1("param1resolG1HeavyBin1","", sigmaG1HeavyBin1->GetParameter(1));
  RooRealVar param2resolG1HeavyBin1("param2resolG1HeavyBin1","", sigmaG1HeavyBin1->GetParameter(2));
  RooRealVar param0resolG2HeavyBin1("param0resolG2HeavyBin1","", sigmaG2HeavyBin1->GetParameter(0));
  RooRealVar param1resolG2HeavyBin1("param1resolG2HeavyBin1","", sigmaG2HeavyBin1->GetParameter(1));
  RooRealVar param2resolG2HeavyBin1("param2resolG2HeavyBin1","", sigmaG2HeavyBin1->GetParameter(2));
  RooRealVar param0respG1HeavyBin1 ("param0respG1HeavyBin1","",  meanG1HeavyBin1->GetParameter(0));
  RooRealVar param1respG1HeavyBin1 ("param1respG1HeavyBin1","",  meanG1HeavyBin1->GetParameter(1));
  RooRealVar param0respG2HeavyBin1 ("param0respG2HeavyBin1","",  meanG2HeavyBin1->GetParameter(0));
  RooRealVar param1respG2HeavyBin1 ("param1respG2HeavyBin1","",  meanG2HeavyBin1->GetParameter(1));



  w->import( param0AccLightBin0 ) ;
  w->import( param1AccLightBin0 ) ;
  w->import( param2AccLightBin0 ) ;
  w->import( param0AccLightBin1 ) ;
  w->import( param1AccLightBin1 ) ;
  w->import( param2AccLightBin1 ) ;
  w->import( param0AccHeavyBin0 ) ;
  w->import( param1AccHeavyBin0 ) ;
  w->import( param2AccHeavyBin0 ) ;
  w->import( param0AccHeavyBin1 ) ;
  w->import( param1AccHeavyBin1 ) ;
  w->import( param2AccHeavyBin1 ) ;

  w->import(  param0resolLightBin0 );
  w->import(  param1resolLightBin0 );
  w->import(  param2resolLightBin0 );
  w->import(  param0resolLightBin1 );
  w->import(  param1resolLightBin1 );
  w->import(  param2resolLightBin1 );
  w->import(  param0respLightBin0 );
  w->import(  param1respLightBin0 );
  w->import(  param0respLightBin1 );
  w->import(  param1respLightBin1 );

  //w->import(  param0resolMassLightBin0 );
  //w->import(  param1resolMassLightBin0 );
  //w->import(  param2resolMassLightBin0 );
  //w->import(  param0resolMassLightBin1 );
  //w->import(  param1resolMassLightBin1 );
  //w->import(  param2resolMassLightBin1 );
  //w->import(  param0respMassLightBin0 );
  //w->import(  param1respLightBin0 );
  //w->import(  param0respMassLightBin1 );
  //w->import(  param1respLightBin1 );

  w->import(  param0resolHeavyBin0 );
  w->import(  param1resolHeavyBin0 );
  w->import(  param2resolHeavyBin0 );
  w->import(  param0resolHeavyBin1 );
  w->import(  param1resolHeavyBin1 );
  w->import(  param2resolHeavyBin1 );
  w->import(  param0respHeavyBin0 );
  w->import(  param1respHeavyBin0 );
  w->import(  param0respHeavyBin1 );
  w->import(  param1respHeavyBin1 );

  //w->import(  param0resolMassHeavyBin0 );
  //w->import(  param1resolMassHeavyBin0 );
  //w->import(  param2resolMassHeavyBin0 );
  //w->import(  param0resolMassHeavyBin1 );
  //w->import(  param1resolMassHeavyBin1 );
  //w->import(  param2resolMassHeavyBin1 );
  //w->import(  param0respMassHeavyBin0 );
  //w->import(  param1respHeavyBin0 );
  //w->import(  param0respMassHeavyBin1 );
  //w->import(  param1respHeavyBin1 );


  // DG
  w->import( param0resolG1HeavyBin0 );
  w->import( param1resolG1HeavyBin0 );
  w->import( param2resolG1HeavyBin0 );
  w->import( param0resolG2HeavyBin0 );
  w->import( param1resolG2HeavyBin0 );
  w->import( param2resolG2HeavyBin0 );
  w->import( param0respG1HeavyBin0 );
  w->import( param1respG1HeavyBin0 );
  //w->import( param2respG1HeavyBin0 );
  w->import( param0respG2HeavyBin0 );
  w->import( param1respG2HeavyBin0 );
  //w->import( param2respG2HeavyBin0 );
  w->import( param0resolG1HeavyBin1 );
  w->import( param1resolG1HeavyBin1 );
  w->import( param2resolG1HeavyBin1 );
  w->import( param0resolG2HeavyBin1 );
  w->import( param1resolG2HeavyBin1 );
  w->import( param2resolG2HeavyBin1 );
  w->import( param0respG1HeavyBin1 );
  w->import( param1respG1HeavyBin1 );
  //w->import( param2respG1HeavyBin1 );
  w->import( param0respG2HeavyBin1 );
  w->import( param1respG2HeavyBin1 );
  //w->import( param2respG2HeavyBin1 );







  RooRealVar flavor("flavor","csv", -6, 6);
  RooRealVar eta   ("eta","eta",    -3, 3);

  RooRealVar csvRecoLight("csvReco","csv", 0., 1.0);  //0.6, 1.0 10 bins
  csvRecoLight.setBins(25);
  RooDataSet datasetBtagLight_Bin0("datasetBtagLight_Bin0","dataset", RooArgSet(csvRecoLight,flavor,eta), Import( *treeJetsLight ), Cut("abs(flavor)!=4 && abs(flavor)!=5 && TMath::Abs(eta)<1.0 && csvReco>-1") ); //0.6
  RooDataSet datasetBtagLight_Bin1("datasetBtagLight_Bin1","dataset", RooArgSet(csvRecoLight,flavor,eta), Import( *treeJetsLight ), Cut("abs(flavor)!=4 && abs(flavor)!=5 && TMath::Abs(eta)>1.0 && TMath::Abs(eta)<2.5 && csvReco>-1") );

  RooRealVar csvRecoCharm("csvReco","csv", 0., 1.0);  //0.6, 1.0 10 bins
  csvRecoCharm.setBins(25);
  RooDataSet datasetBtagCharm_Bin0("datasetBtagCharm_Bin0","dataset", RooArgSet(csvRecoCharm,flavor,eta), Import( *treeJetsLight ), Cut("abs(flavor)==4 && TMath::Abs(eta)<1.0 && csvReco>-1") ); //0.6
  RooDataSet datasetBtagCharm_Bin1("datasetBtagCharm_Bin1","dataset", RooArgSet(csvRecoCharm,flavor,eta), Import( *treeJetsLight ), Cut("abs(flavor)==4 && TMath::Abs(eta)>1.0 && TMath::Abs(eta)<2.5 && csvReco>-1") );

  RooKeysPdf pdfCsvLight_Bin0("pdfCsvLight_Bin0","", csvRecoLight, datasetBtagLight_Bin0, RooKeysPdf::NoMirror, 1.); //1.5
  RooKeysPdf pdfCsvLight_Bin1("pdfCsvLight_Bin1","", csvRecoLight, datasetBtagLight_Bin1, RooKeysPdf::NoMirror, 1.);

  w->import( pdfCsvLight_Bin0 );
  w->import( pdfCsvLight_Bin1 );

  RooKeysPdf pdfCsvCharm_Bin0("pdfCsvCharm_Bin0","", csvRecoCharm, datasetBtagCharm_Bin0, RooKeysPdf::NoMirror, 1.);
  RooKeysPdf pdfCsvCharm_Bin1("pdfCsvCharm_Bin1","", csvRecoCharm, datasetBtagCharm_Bin1, RooKeysPdf::NoMirror, 1.);

  w->import( pdfCsvCharm_Bin0 );
  w->import( pdfCsvCharm_Bin1 );
  

  RooRealVar csvRecoHeavy("csvReco","csv", 0., 1.0); //0.679, 1.0  30 bins
  csvRecoHeavy.setBins(100);
  RooDataSet datasetBtagHeavy_Bin0("datasetBtagHeavy_Bin0","dataset", RooArgSet(csvRecoHeavy,flavor,eta), Import( *treeJetsHeavy ), Cut("TMath::Abs(eta)<1.0 && csvReco>-1") ); //0.679
  RooDataSet datasetBtagHeavy_Bin1("datasetBtagHeavy_Bin1","dataset", RooArgSet(csvRecoHeavy,flavor,eta), Import( *treeJetsHeavy ), Cut("TMath::Abs(eta)>1.0 && TMath::Abs(eta)<2.5 && csvReco>-1") );
  //RooKeysPdf pdfCsvHeavy_Bin0("pdfCsvHeavy_Bin0","", csvRecoHeavy, datasetBtagHeavy_Bin0);
  //RooKeysPdf pdfCsvHeavy_Bin1("pdfCsvHeavy_Bin1","", csvRecoHeavy, datasetBtagHeavy_Bin1);

  RooDataHist datahistBtagHeavy_Bin0("datahistBtagHeavy_Bin0","datahistBtagHeavy_Bin0", RooArgSet(csvRecoHeavy), datasetBtagHeavy_Bin0, 1.0);
  RooDataHist datahistBtagHeavy_Bin1("datahistBtagHeavy_Bin1","datahistBtagHeavy_Bin1", RooArgSet(csvRecoHeavy), datasetBtagHeavy_Bin1, 1.0);

  RooHistPdf  pdfCsvHeavy_Bin0  ("pdfCsvHeavy_Bin0", "pdfCsvHeavy_Bin0",  RooArgSet(csvRecoHeavy), datahistBtagHeavy_Bin0);
  RooHistPdf  pdfCsvHeavy_Bin1  ("pdfCsvHeavy_Bin1", "pdfCsvHeavy_Bin1",  RooArgSet(csvRecoHeavy), datahistBtagHeavy_Bin1);


  w->import( pdfCsvHeavy_Bin0 );
  w->import( pdfCsvHeavy_Bin1 );

  // Binning(30,0.679, 1.0 )

  TH1F* csv_b_Bin0 = (TH1F*)pdfCsvHeavy_Bin0.createHistogram("csv_b_Bin0", csvRecoHeavy, Binning(100,0., 1.0 ));
  csv_b_Bin0->Scale( 1./csv_b_Bin0->Integral()/csv_b_Bin0->GetBinWidth(1));
  TH1F* csv_b_Bin1 = (TH1F*)pdfCsvHeavy_Bin1.createHistogram("csv_b_Bin1", csvRecoHeavy, Binning(100,0., 1.0 ));
  csv_b_Bin1->Scale( 1./csv_b_Bin1->Integral()/csv_b_Bin1->GetBinWidth(1));

  TH1F* csv_l_Bin0 = (TH1F*)pdfCsvLight_Bin0.createHistogram("csv_l_Bin0", csvRecoLight, Binning(100,0., 1.0 ));
  csv_l_Bin0->Scale( 1./csv_l_Bin0->Integral()/csv_l_Bin0->GetBinWidth(1));
  TH1F* csv_l_Bin1 = (TH1F*)pdfCsvLight_Bin1.createHistogram("csv_l_Bin1", csvRecoLight, Binning(100,0., 1.0 ));
  csv_l_Bin1->Scale( 1./csv_l_Bin1->Integral()/csv_l_Bin1->GetBinWidth(1));

  TH1F* csv_c_Bin0 = (TH1F*)pdfCsvCharm_Bin0.createHistogram("csv_c_Bin0", csvRecoCharm, Binning(100,0., 1.0 ));
  csv_c_Bin0->Scale( 1./csv_c_Bin0->Integral()/csv_c_Bin0->GetBinWidth(1));
  TH1F* csv_c_Bin1 = (TH1F*)pdfCsvCharm_Bin1.createHistogram("csv_c_Bin1", csvRecoCharm, Binning(100,0., 1.0 ));
  csv_c_Bin1->Scale( 1./csv_c_Bin1->Integral()/csv_c_Bin1->GetBinWidth(1));


  ///////////////////////////////////////////////////////////////////////////////////////////

  float sumEtBins[11] = {500,750, 1000, 1250, 1500, 1750, 2000, 2250, 2500, 2750, 3000};
  TH1F* hWidthsPxResol = new TH1F("hWidthsPxResol", "", 10, 500, 3000);
  for(unsigned int it = 0; it<10 ; it++){
    TH1F* hPx = new TH1F("hPx","",100,-200,200);
    treeEvent->Draw("(pxReco-px)>>hPx",Form("(px>-998 && sumEt>=%f && sumEt<%f)", sumEtBins[it],  sumEtBins[it+1])); // check here!
    TF1* pxResol = new TF1("pxResol","gaus", -100,100);
    if(hPx->Integral()<10) continue;
    cout << "Fit " << sumEtBins[it] << "," << sumEtBins[it+1] << endl;
    hPx->Fit(pxResol);
    hWidthsPxResol->SetBinContent(it+1, pxResol->GetParameter(2));
    hWidthsPxResol->SetBinError  (it+1, pxResol->GetParError(2));
    delete hPx; delete pxResol;
  }
  TF1* hWidthsPxResolModel = new TF1("hWidthsPxResolModel", "sqrt([0]*[0]*x + [1]*[1])", 500, 3000);
  cout << "Fit model" << endl;
  hWidthsPxResol->Fit( hWidthsPxResolModel );
  cout << "Fited model" << endl;
  RooRealVar param0PxWidthModel("param0PxWidthModel","", hWidthsPxResolModel->GetParameter(0) );
  RooRealVar param1PxWidthModel("param1PxWidthModel","", hWidthsPxResolModel->GetParameter(1) );
  w->import( param0PxWidthModel );
  w->import( param1PxWidthModel );

  ///////////////////////////////////////////////////////////////////////////////////////////


  RooRealVar BetaWHad ("BetaW", "BetaWHad",    -1,1);  
  RooRealVar GammaWHad("GammaW","GammaWHad",   -1,1);  
  RooRealVar BetaWLep ("BetaWLep", "BetaWLep", -1,1);  
  RooRealVar GammaWLep("GammaWLep","GammaWLep",-1,1);  

  RooDataSet dataset("dataset","dataset", RooArgSet(BetaWHad,GammaWHad,BetaWLep,GammaWLep), Import( *tree ) );


  // v1)
  RooDataHist histGammaWHad("histGammaWHad","histGammaWHad", RooArgSet(GammaWHad), dataset, 1.0);
  RooRealVar p0v1("p0v1","",  1,  0,100);
  RooRealVar p1v1("p1v1","",  1,  0,100);
  RooRealVar p2v1("p2v1","",  1,  0,100);
  RooRealVar p3v1("p3v1","",  1,  0,100);
  RooBernstein pdfGammaWHad( "pdfGammaWHad", "pdfGammaWHad", GammaWHad, RooArgList(p0v1,p1v1,p2v1,p3v1));
  pdfGammaWHad.fitTo(histGammaWHad,  Minos(1), Save(1), Range(-1.,1.) ); 
  w->import(pdfGammaWHad);


  // v2)
  RooDataHist histBetaWHad("histBetaWHad","histBetaWHad", RooArgSet(BetaWHad), dataset, 1.0);
  RooRealVar p0v2("p0v2","",  1,  0,100);
  RooRealVar p1v2("p1v2","",  1,  0,100);
  RooBernstein pdfBetaWHad( "pdfBetaWHad", "pdfBetaWHad", BetaWHad, RooArgList(p0v2,p1v2));
  pdfBetaWHad.fitTo(histBetaWHad,  Minos(1), Save(1), Range(-1.,1.) ); 
  w->import(pdfBetaWHad);

  // v3)
  RooDataHist histGammaWLep("histGammaWLep","histGammaWLep", RooArgSet(GammaWLep), dataset, 1.0);
  RooRealVar p0v3("p0v3","",  1,  0,100);
  RooRealVar p1v3("p1v3","",  1,  0,100);
  RooRealVar p2v3("p2v3","",  1,  0,100);
  RooRealVar p2v4("p2v4","",  1,  0,100);
  RooBernstein pdfGammaWLep( "pdfGammaWLep", "pdfGammaWLep", GammaWLep, RooArgList(p0v3,p1v3,p2v3,p2v4));
  pdfGammaWLep.fitTo(histGammaWLep,  Minos(1), Save(1), Range(-1.,1.) ); 
  w->import(pdfGammaWLep);


  // v4)
  RooDataHist histBetaWLep("histBetaWLep","histBetaWLep", RooArgSet(BetaWLep), dataset, 1.0);
  RooRealVar p0v4("p0v4","",  1,  0,100);
  RooRealVar p1v4("p1v4","",  1,  0,100);
  RooBernstein pdfBetaWLep( "pdfBetaWLep", "pdfBetaWLep", BetaWLep, RooArgList(p0v4,p1v4));
  pdfBetaWLep.fitTo(histBetaWLep,  Minos(1), Save(1), Range(-1.,1.) ); 
  w->import(pdfBetaWLep);





  w->writeToFile(outFileName.c_str(),kTRUE);

  TCanvas *c1CsvLight_Bin0 = new TCanvas("c1CsvLight_Bin0","canvas",10,30,650,600);
  RooPlot* plotCsvLight_Bin0 = csvRecoLight.frame(Bins(10),Title("cvs light, Bin0"));
  datasetBtagLight_Bin0.plotOn(plotCsvLight_Bin0);
  pdfCsvLight_Bin0.plotOn(plotCsvLight_Bin0);

  TCanvas *c1CsvLight_Bin1 = new TCanvas("c1CsvLight_Bin1","canvas",10,30,650,600);
  RooPlot* plotCsvLight_Bin1 = csvRecoLight.frame(Bins(10),Title("cvs light, Bin1"));
  datasetBtagLight_Bin1.plotOn(plotCsvLight_Bin1);
  pdfCsvLight_Bin1.plotOn(plotCsvLight_Bin1);

  TCanvas *c1CsvCharm_Bin0 = new TCanvas("c1CsvCharm_Bin0","canvas",10,30,650,600);
  RooPlot* plotCsvCharm_Bin0 = csvRecoCharm.frame(Bins(10),Title("cvs light, Bin0"));
  datasetBtagCharm_Bin0.plotOn(plotCsvCharm_Bin0);
  pdfCsvCharm_Bin0.plotOn(plotCsvCharm_Bin0);

  TCanvas *c1CsvCharm_Bin1 = new TCanvas("c1CsvCharm_Bin1","canvas",10,30,650,600);
  RooPlot* plotCsvCharm_Bin1 = csvRecoCharm.frame(Bins(10),Title("cvs light, Bin1"));
  datasetBtagCharm_Bin1.plotOn(plotCsvCharm_Bin1);
  pdfCsvCharm_Bin1.plotOn(plotCsvCharm_Bin1);

 
  TCanvas *c1CsvHeavy_Bin0 = new TCanvas("c1CsvHeavy_Bin0","canvas",10,30,650,600);
  RooPlot* plotCsvHeavy_Bin0 = csvRecoHeavy.frame(Bins(30),Title("cvs light, Bin0"));
  datasetBtagHeavy_Bin0.plotOn(plotCsvHeavy_Bin0);
  pdfCsvHeavy_Bin0.plotOn(plotCsvHeavy_Bin0);

  TCanvas *c1CsvHeavy_Bin1 = new TCanvas("c1CsvHeavy_Bin1","canvas",10,30,650,600);
  RooPlot* plotCsvHeavy_Bin1 = csvRecoHeavy.frame(Bins(30),Title("cvs light, Bin1"));
  datasetBtagHeavy_Bin1.plotOn(plotCsvHeavy_Bin1);
  pdfCsvHeavy_Bin1.plotOn(plotCsvHeavy_Bin1);

  TCanvas *c1BetaWHad = new TCanvas("c1BetaWHad","canvas",10,30,650,600);
  RooPlot* plotBetaWHad = BetaWHad.frame(Bins(20),Title("BetaWHad"));
  histBetaWHad.plotOn(plotBetaWHad);
  pdfBetaWHad.plotOn(plotBetaWHad);

  TCanvas *c1GammaWHad = new TCanvas("c1GammaWHad","canvas",10,30,650,600);
  RooPlot* plotGammaWHad = GammaWHad.frame(Bins(20),Title("GammaWHad"));
  histGammaWHad.plotOn(plotGammaWHad);
  pdfGammaWHad.plotOn(plotGammaWHad);

  TCanvas *c1BetaWLep = new TCanvas("c1BetaWLep","canvas",10,30,650,600);
  RooPlot* plotBetaWLep = BetaWLep.frame(Bins(20),Title("BetaWLep"));
  histBetaWLep.plotOn(plotBetaWLep);
  pdfBetaWLep.plotOn(plotBetaWLep);

  TCanvas *c1GammaWLep = new TCanvas("c1GammaWLep","canvas",10,30,650,600);
  RooPlot* plotGammaWLep = GammaWLep.frame(Bins(20),Title("GammaWLep"));
  histGammaWLep.plotOn(plotGammaWLep);
  pdfGammaWLep.plotOn(plotGammaWLep);


  //TCanvas *c1V1 = new TCanvas("c1V1","canvas",10,30,650,600);
  //RooPlot* plotV1 = v1.frame(Bins(20),Title("v1"));
  //histV1.plotOn(plotV1);
  //pdfV1.plotOn(plotV1);

  //TCanvas *c1V2 = new TCanvas("c1V2","canvas",10,30,650,600);
  //RooPlot* plotV2 = v2.frame(Bins(20),Title("v2"));
  //histV2.plotOn(plotV2);
  //pdfV2.plotOn(plotV2);

  //TCanvas *c1PtTTH = new TCanvas("c1PtTTH","canvas",10,30,650,600);
  //RooPlot* plotPtTTH = PtTTH.frame(Bins(20),Title("PtTTH"));
  //histPtTTH.plotOn(plotPtTTH);
  //pdfPtTTH.plotOn(plotPtTTH);

  //TCanvas *c1PtHStar = new TCanvas("c1PtHStar","canvas",10,30,650,600);
  //RooPlot* plotPtHStar = PtHStar.frame(Bins(20),Title("PtHStar"));
  //histPtHStar.plotOn(plotPtHStar);
  //pdfPtHStar.plotOn(plotPtHStar);

  //TCanvas *c1DphiHStar = new TCanvas("c1DphiHStar","canvas",10,30,650,600);
  //RooPlot* plotDphiHStar = DphiHStar.frame(Bins(20),Title("DphiHStar"));
  //histDphiHStar.plotOn(plotDphiHStar);
  //pdfDphiHStar.plotOn(plotDphiHStar);

  
  //fout = new TFile("./root/ControlPlots_new.root","UPDATE");
  fout->cd();

  c1CsvLight_Bin0->cd();
  plotCsvLight_Bin0->Draw();
  c1CsvLight_Bin0->Write("",TObject::kOverwrite);

  c1CsvLight_Bin1->cd();
  plotCsvLight_Bin1->Draw();
  c1CsvLight_Bin1->Write("",TObject::kOverwrite);

  c1CsvCharm_Bin0->cd();
  plotCsvCharm_Bin0->Draw();
  c1CsvCharm_Bin0->Write("",TObject::kOverwrite);

  c1CsvCharm_Bin1->cd();
  plotCsvCharm_Bin1->Draw();
  c1CsvCharm_Bin1->Write("",TObject::kOverwrite);
  
  c1CsvHeavy_Bin0->cd();
  plotCsvHeavy_Bin0->Draw();
  c1CsvHeavy_Bin0->Write("",TObject::kOverwrite);

  c1CsvHeavy_Bin1->cd();
  plotCsvHeavy_Bin1->Draw();
  c1CsvHeavy_Bin1->Write("",TObject::kOverwrite);
  
  c1BetaWHad->cd();
  plotBetaWHad->Draw();
  c1BetaWHad->Write("",TObject::kOverwrite);

  c1GammaWHad->cd();
  plotGammaWHad->Draw();
  c1GammaWHad->Write("",TObject::kOverwrite);

  c1BetaWLep->cd();
  plotBetaWLep->Draw();
  c1BetaWLep->Write("",TObject::kOverwrite);

  c1GammaWLep->cd();
  plotGammaWLep->Draw();
  c1GammaWLep->Write("",TObject::kOverwrite);

  //c1V1->cd();
  //plotV1->Draw();
  //c1V1->Write("",TObject::kOverwrite);

  //c1V2->cd();
  //plotV2->Draw();
  //c1V2->Write("",TObject::kOverwrite);

  //c1PtTTH->cd();
  //plotPtTTH->Draw();
  //c1PtTTH->Write("",TObject::kOverwrite);

  //c1PtHStar->cd();
  //plotPtHStar->Draw();
  //c1PtHStar->Write("",TObject::kOverwrite);

  //c1DphiHStar->cd();
  //plotDphiHStar->Draw();
  //c1DphiHStar->Write("",TObject::kOverwrite);
  //param0DphiHStar->Write("",TObject::kOverwrite);
  //param1DphiHStar->Write("",TObject::kOverwrite);
  //param2DphiHStar->Write("",TObject::kOverwrite);
  //ratioDphiHStar->Write( "",TObject::kOverwrite);

  accLightBin0->Write("",TObject::kOverwrite);
  accLightBin1->Write("",TObject::kOverwrite);
  accHeavyBin0->Write("",TObject::kOverwrite);
  accHeavyBin1->Write("",TObject::kOverwrite);
  resolLightBin0->Write("",TObject::kOverwrite);
  respLightBin0->Write("",TObject::kOverwrite);
  resolLightBin1->Write("",TObject::kOverwrite);
  respLightBin1->Write("",TObject::kOverwrite);

  //resolMassLightBin0->Write("",TObject::kOverwrite);
  //respMassLightBin0->Write("",TObject::kOverwrite);
  //resolMassLightBin1->Write("",TObject::kOverwrite);
  //respMassLightBin1->Write("",TObject::kOverwrite);
  //corrMassLightBin0->Write("",TObject::kOverwrite);
  //corrMassLightBin1->Write("",TObject::kOverwrite);
  //resolMassHeavyBin0->Write("",TObject::kOverwrite);
  //respMassHeavyBin0->Write("",TObject::kOverwrite);
  //resolMassHeavyBin1->Write("",TObject::kOverwrite);
  //respMassHeavyBin1->Write("",TObject::kOverwrite);
  //corrMassHeavyBin0->Write("",TObject::kOverwrite);
  //corrMassHeavyBin1->Write("",TObject::kOverwrite);


  resolHeavyBin0->Write("",TObject::kOverwrite);
  respHeavyBin0->Write("",TObject::kOverwrite);
  resolHeavyBin1->Write("",TObject::kOverwrite);
  respHeavyBin1->Write("",TObject::kOverwrite);

  //fracDGLightBin0->Write("",TObject::kOverwrite);
  //respDGLightBin0->Write("",TObject::kOverwrite);
  //resolG1LightBin0->Write("",TObject::kOverwrite);
  //resolG2LightBin0->Write("",TObject::kOverwrite);
  //fracDGLightBin1->Write("",TObject::kOverwrite);
  //respDGLightBin1->Write("",TObject::kOverwrite);
  //resolG1LightBin1->Write("",TObject::kOverwrite);
  //resolG2LightBin1->Write("",TObject::kOverwrite);

  respG1HeavyBin0->Write("",TObject::kOverwrite);
  respG2HeavyBin0->Write("",TObject::kOverwrite);
  resolG1HeavyBin0->Write("",TObject::kOverwrite);
  resolG2HeavyBin0->Write("",TObject::kOverwrite);
  respG1HeavyBin1->Write("",TObject::kOverwrite);
  respG2HeavyBin1->Write("",TObject::kOverwrite);
  resolG1HeavyBin1->Write("",TObject::kOverwrite);
  resolG2HeavyBin1->Write("",TObject::kOverwrite);

  csv_b_Bin0->Write("",TObject::kOverwrite);
  csv_b_Bin1->Write("",TObject::kOverwrite);
  csv_l_Bin0->Write("",TObject::kOverwrite);
  csv_l_Bin1->Write("",TObject::kOverwrite);
  csv_c_Bin0->Write("",TObject::kOverwrite);
  csv_c_Bin1->Write("",TObject::kOverwrite);

  hWidthsPxResol->Write("",TObject::kOverwrite);
  
  //hCorrMassLightBin0->Write("",TObject::kOverwrite);
  //hCorrMassLightBin1->Write("",TObject::kOverwrite);
  //hCorrMassHeavyBin0->Write("",TObject::kOverwrite);
  //hCorrMassHeavyBin1->Write("",TObject::kOverwrite);


  fout->Close();
  delete fout;


  return 0;
}
