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
#include "TCut.h"
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

using namespace std;
using namespace RooFit;


int main(int argc, const char* argv[])
{

  std::cout << "TransferFunctions" << std::endl;
  gROOT->SetBatch(true);
 
  gSystem->Load("libFWCoreFWLite");
  gSystem->Load("libDataFormatsFWLite");

  AutoLibraryLoader::enable();

  TString extraname =  argc>2 ?  "TEST"+TString(argv[2]) : "TEST";
  TFile* fout = 0;
  fout = new TFile("./root/ControlPlots"+extraname+".root","RECREATE");
  fout->cd();

  TString reg =  argc>2 ? TString(argv[2]) : "" ;

  PythonProcessDesc builder(argv[1]);
  const edm::ParameterSet& in = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteInput");
  std::string outFileName( in.getParameter<std::string>    ("outFileName") );
  std::string pathToFile ( in.getParameter<std::string>    ("pathToFile" ) );
  int doOnlyJetTF        ( in.getUntrackedParameter<int>   ("doOnlyJetTF", 0 ) );
  double relWeightBin0   ( in.getUntrackedParameter<double>("relWeightBin0", 0.65 ));
  double relWeightBin1   ( in.getUntrackedParameter<double>("relWeightBin1", 0.65 ));
  double relShiftBin0    ( in.getUntrackedParameter<double>("relShiftBin0",  0.92 ));

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
  int nBinsY = 58;

  // RESOLUTION (GAUSSIAN)
  TH1F* resolLightBin0  = new TH1F("resolLightBin0","",nBinsY,30,300);
  TH1F* respLightBin0   = new TH1F("respLightBin0","", nBinsY,30,300);
  TH1F* resolLightBin1  = new TH1F("resolLightBin1","",nBinsY,30,300);
  TH1F* respLightBin1   = new TH1F("respLightBin1","", nBinsY,30,300);

  TH1F* resolHeavyBin0  = new TH1F("resolHeavyBin0","",nBinsY,30,300);
  TH1F* respHeavyBin0   = new TH1F("respHeavyBin0","", nBinsY,30,300);
  TH1F* resolHeavyBin1  = new TH1F("resolHeavyBin1","",nBinsY,30,300);
  TH1F* respHeavyBin1   = new TH1F("respHeavyBin1","", nBinsY,30,300);


  // RESOLUTION (DOUBLE GAUSSIAN)
  TH1F* respG1HeavyBin0     = new TH1F("respG1HeavyBin0","",   nBinsY,30,300);
  TH1F* respG2HeavyBin0     = new TH1F("respG2HeavyBin0","",   nBinsY,30,300);  
  TH1F* resolG1HeavyBin0    = new TH1F("resolG1HeavyBin0","",  nBinsY,30,300);
  TH1F* resolG2HeavyBin0    = new TH1F("resolG2HeavyBin0","",  nBinsY,30,300);
  TH1F* respG1HeavyBin1     = new TH1F("respG1HeavyBin1","",   nBinsY,30,300);
  TH1F* respG2HeavyBin1     = new TH1F("respG2HeavyBin1","",   nBinsY,30,300);  
  TH1F* resolG1HeavyBin1    = new TH1F("resolG1HeavyBin1","",  nBinsY,30,300);
  TH1F* resolG2HeavyBin1    = new TH1F("resolG2HeavyBin1","",  nBinsY,30,300);

  TString SG_formula("[0]*exp(-0.5*(x-[1])*(x-[1])/[2]/[2])");
  TString DG_formula_Bin0(Form("[0]*( %f*exp(-0.5*(x-[1])*(x-[1])/[2]/[2]) + (1-%f)*exp(-0.5*(x-[1]*%f)*(x-[1]*%f)/[3]/[3]))",     relWeightBin0, relWeightBin0, relShiftBin0, relShiftBin0 ));
  TString DG_formula_Bin1(Form("[0]*( %f*exp(-0.5*(x-[1])*(x-[1])/[2]/[2]) + (1-%f)*exp(-0.5*(x-[4])*(x-[4])/[3]/[3]))",           relWeightBin1, relWeightBin1));

  for(int k = 0; k < 2; k++){

    TCut  cut_Bin_Light(Form("(pt_rec>20)  *(TMath::Abs(eta_part)>%f && TMath::Abs(eta_part)<%f)",              etaBinning[k],etaBinning[k+1]));
    TCut  cut_Bin_Heavy(Form("(pt_rec%s>20)*(TMath::Abs(eta_part)>%f && TMath::Abs(eta_part)<%f)",  reg.Data(), etaBinning[k],etaBinning[k+1]));
    //TCut  cut_Bin_Heavy(Form("(pt_rec>30)*(TMath::Abs(eta_part)>%f && TMath::Abs(eta_part)<%f)",                etaBinning[k],etaBinning[k+1]));

    TH2F* hCorrLight = new TH2F("hCorrLight","", nBinsY, 30,300 , 100, 0,500);
    treeJetsLight->Draw("e_rec:e_part>>hCorrLight", cut_Bin_Light );

    for(int j = 1; j<=hCorrLight->GetNbinsX(); j++){

      TH1F* hSlice = (TH1F*)hCorrLight->ProjectionY("_py",j,j);
      if(hSlice->GetEntries()<10) continue;
      cout << hSlice->GetMean() << ", " << hSlice->GetRMS() << endl;

      TF1* resol = new TF1("resol","gaus", 0,500);
      hSlice->Fit(resol);

      fout->mkdir(Form("BinLight%d_%d",j,k));
      fout->cd(Form("BinLight%d_%d",j,k) );
      hSlice    ->Write("",TObject::kOverwrite);   	         
     
      if(k==0){
	resolLightBin0->SetBinContent(j, resol->GetParameter(2));
	resolLightBin0->SetBinError  (j, resol->GetParError (2));
	respLightBin0->SetBinContent (j, resol->GetParameter(1));
	respLightBin0->SetBinError   (j, resol->GetParError (1));
      }
      if(k==1){
	resolLightBin1->SetBinContent(j, resol->GetParameter(2));
	resolLightBin1->SetBinError  (j, resol->GetParError(2));
	respLightBin1->SetBinContent (j, resol->GetParameter(1));
	respLightBin1->SetBinError   (j, resol->GetParError(1));  
      }
       
      delete resol; 

    }


    TH2F* hCorrHeavy = new TH2F("hCorrHeavy","",  nBinsY, 30,300 , 100, 0,500);
    treeJetsHeavy->Draw("e_rec"+reg+":e_part>>hCorrHeavy", cut_Bin_Heavy );

    for(int j = 1; j<=hCorrHeavy->GetNbinsX(); j++){

      TH1F* hSlice = (TH1F*)hCorrHeavy ->ProjectionY("_py",j,j);
      if(hSlice->GetEntries()<8) continue;
      cout << hSlice->GetMean() << ", " << hSlice->GetRMS() << endl;
      TF1* resol = new TF1("resol","gaus", 0,500);
      hSlice->Fit(resol);

      fout->mkdir(Form("BinHeavy%d_%d",j,k));
      fout->cd( Form("BinHeavy%d_%d",j, k) );
      hSlice->Write("",TObject::kOverwrite);
  
      TF1* resolDG = 0;
      
      if(k==0){	
	resolDG = new TF1("resolDG", DG_formula_Bin0, 0,500);
	resolDG->SetParameter(1, hSlice->GetMean());
	resolDG->SetParameter(2, hSlice->GetRMS());
	resolDG->SetParameter(3, hSlice->GetRMS());
      } 
      if(k==1){
	resolDG = new TF1("resolDG", DG_formula_Bin1, 0,500);
	resolDG->SetParameter(1, hSlice->GetMean());
	resolDG->SetParameter(2, hSlice->GetRMS());
	resolDG->SetParameter(3, hSlice->GetRMS());
	resolDG->SetParLimits(4, hSlice->GetMean()-2*hSlice->GetRMS() ,hSlice->GetMean());  
      }

      hSlice->Fit(resolDG);
      hSlice->Fit(resolDG);

      hSlice->Write("",TObject::kOverwrite);

      if(k==0){
	resolHeavyBin0->SetBinContent(j, resol->GetParameter(2));
	resolHeavyBin0->SetBinError  (j, resol->GetParError (2));
	respHeavyBin0->SetBinContent (j, resol->GetParameter(1));
	respHeavyBin0->SetBinError   (j, resol->GetParError (1));
	respG1HeavyBin0->SetBinContent   (j,  resolDG->GetParameter(1) );
	respG1HeavyBin0->SetBinError     (j,  resolDG->GetParError (1) );
	respG2HeavyBin0->SetBinContent   (j,  relShiftBin0 * resolDG->GetParameter(1) );
	respG2HeavyBin0->SetBinError     (j,  relShiftBin0 * resolDG->GetParError (1) );
	//respG2HeavyBin0->SetBinContent   (j,  resolDG->GetParameter(4) );
	//respG2HeavyBin0->SetBinError     (j,  resolDG->GetParError (4) );
	resolG1HeavyBin0->SetBinContent  (j,  resolDG->GetParameter(2) );
	resolG1HeavyBin0->SetBinError    (j,  resolDG->GetParError (2) );
	resolG2HeavyBin0->SetBinContent  (j,  resolDG->GetParameter(3) ); 
	resolG2HeavyBin0->SetBinError    (j,  resolDG->GetParError (3) );
      }
      if(k==1){
	resolHeavyBin1->SetBinContent(j, resol->GetParameter(2));
	resolHeavyBin1->SetBinError  (j, resol->GetParError (2));
	respHeavyBin1->SetBinContent (j, resol->GetParameter(1));
	respHeavyBin1->SetBinError   (j, resol->GetParError (1));
	respG1HeavyBin1->SetBinContent   (j,  resolDG->GetParameter(1) );
	respG1HeavyBin1->SetBinError     (j,  resolDG->GetParError (1) );
	respG2HeavyBin1->SetBinContent   (j,  resolDG->GetParameter(4) );
	respG2HeavyBin1->SetBinError     (j,  resolDG->GetParError (4) );
	resolG1HeavyBin1->SetBinContent  (j,  resolDG->GetParameter(2) );
	resolG1HeavyBin1->SetBinError    (j,  resolDG->GetParError (2) );
	resolG2HeavyBin1->SetBinContent  (j,  resolDG->GetParameter(3) );
	resolG2HeavyBin1->SetBinError    (j,  resolDG->GetParError (3) );
      }


    
      delete resol; 
      delete resolDG;
    }


    delete hCorrLight;  delete hCorrHeavy;
  }

  TF1* sigmaLightBin0 = new TF1("sigmaLightBin0","sqrt([0]*[0] + x*[1]*[1] + x*x*[2]*[2])",0,400);
  TF1* sigmaLightBin1 = new TF1("sigmaLightBin1","sqrt([0]*[0] + x*[1]*[1] + x*x*[2]*[2])",0,400);
  TF1* sigmaHeavyBin0 = new TF1("sigmaHeavyBin0","sqrt([0]*[0] + x*[1]*[1] + x*x*[2]*[2])",0,400);
  TF1* sigmaHeavyBin1 = new TF1("sigmaHeavyBin1","sqrt([0]*[0] + x*[1]*[1] + x*x*[2]*[2])",0,400);
 
  TF1* meanLightBin0  = new TF1("meanLightBin0", "x*[0]",0,400);
  TF1* meanLightBin1  = new TF1("meanLightBin1", "x*[0]",0,400);
  TF1* meanHeavyBin0  = new TF1("meanHeavyBin0", "x*[0]",0,400);
  TF1* meanHeavyBin1  = new TF1("meanHeavyBin1", "x*[0]",0,400);

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
  resolHeavyBin0->Fit(sigmaHeavyBin0,"","",60,300);

  resolHeavyBin1->Fit(sigmaHeavyBin1,"","",40,200);
  resolHeavyBin1->Fit(sigmaHeavyBin1,"","",40,300);

  respLightBin0->Fit (meanLightBin0, "","",40,200);
  respLightBin1->Fit (meanLightBin1, "","",40,200);
  respHeavyBin0->Fit (meanHeavyBin0, "","",40,200);
  respHeavyBin1->Fit (meanHeavyBin1, "","",40,200);

  // responses
  respG1HeavyBin0->Fit  (meanG1HeavyBin0, "","",60,160);
  respG1HeavyBin0->Fit  (meanG1HeavyBin0, "","",40,200);

  respG2HeavyBin0->Fit  (meanG2HeavyBin0, "","",60,160);
  respG2HeavyBin0->Fit  (meanG2HeavyBin0, "","",40,200);

  respG1HeavyBin1->Fit  (meanG1HeavyBin1, "","",80,200);
  respG1HeavyBin1->Fit  (meanG1HeavyBin1, "","",60,300);
  respG1HeavyBin1->Fit  (meanG1HeavyBin1, "","",40,300);

  respG2HeavyBin1->Fit  (meanG2HeavyBin1, "","",80,200);
  respG2HeavyBin1->Fit  (meanG2HeavyBin1, "","",60,300);
  respG2HeavyBin1->Fit  (meanG2HeavyBin1, "","",40,300);
  
  // resolutions
  resolG1HeavyBin0->Fit (sigmaG1HeavyBin0,"","",100,300);
  //resolG1HeavyBin0->Fit (sigmaG1HeavyBin0,"","",40, 250);

  resolG2HeavyBin0->Fit (sigmaG2HeavyBin0,"","",100,300);
  //resolG2HeavyBin0->Fit (sigmaG2HeavyBin0,"","",40, 250);
 
  resolG1HeavyBin1->Fit (sigmaG1HeavyBin1,"","",30,  300);
  resolG1HeavyBin1->Fit (sigmaG1HeavyBin1,"","",60,  300);
  resolG1HeavyBin1->Fit (sigmaG1HeavyBin1,"","",80,  300);
  //resolG1HeavyBin1->Fit (sigmaG1HeavyBin1,"","",50, 300);

  resolG2HeavyBin1->Fit (sigmaG2HeavyBin1,"","",130,300);
  resolG2HeavyBin1->Fit (sigmaG2HeavyBin1,"","",130,300);
  //resolG2HeavyBin1->Fit (sigmaG2HeavyBin1,"","",50, 180);
  //resolG2HeavyBin1->Fit (sigmaG2HeavyBin1,"","",50, 250);
  //resolG2HeavyBin1->Fit (sigmaG2HeavyBin1,"","",100,300);
  //resolG2HeavyBin1->Fit (sigmaG2HeavyBin1,"","",120,300);
  //resolG2HeavyBin1->Fit (sigmaG2HeavyBin1,"","",80, 300);
  //resolG2HeavyBin1->Fit (sigmaG2HeavyBin1,"","",40, 250);



  /////////////////////

  for(int k = 0; k < 2; k++){

    TH2F* hCorrHeavy = new TH2F("hCorrHeavy","",  nBinsY, 30,300 , 100, 0,500);
    for(int j = 1; j<=hCorrHeavy->GetNbinsX(); j++){

      float en = hCorrHeavy->GetBinCenter(j);

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
      TF1* resolSG = 0;
      if(k==0){

	resolDG = new TF1("resolDG_fit", DG_formula_Bin0 , -800,800);
	resolDG->SetNpx(1000);
	resolDG->SetLineColor(kBlue);

	float norm = resolDG_old->Integral(-800,800);
	float m1   = meanG1HeavyBin0->Eval(en);
	float m2   = meanG2HeavyBin0->Eval(en);
	float s1   = sigmaG1HeavyBin0->Eval(en);
	float s2   = sigmaG2HeavyBin0->Eval(en);
	resolDG->SetParameter(0, norm/(relWeightBin0*sqrt(2*TMath::Pi())*s1 + (1-relWeightBin0)*sqrt(2*TMath::Pi())*s2 )  );
	resolDG->SetParameter(1, m1);
	resolDG->SetParameter(2, s1);
	resolDG->SetParameter(3, s2);
	resolDG->SetParameter(4, m2);

	resolSG = new TF1("resolSG_fit", SG_formula , -800, 800);
	resolSG->SetNpx(1000);
	resolSG->SetLineColor(kMagenta);	
	float m = meanHeavyBin0->Eval(en);
	float s = sigmaHeavyBin0->Eval(en);
	resolSG->SetParameter(0, norm/sqrt(2*TMath::Pi())/s );
	resolSG->SetParameter(1, m);
	resolSG->SetParameter(2, s);

	fout->cd( Form("BinHeavy%d_%d",j, k) );
	resolDG_old->Write(Form("resolDG%d_%d",    j,k),TObject::kOverwrite);
	resolDG    ->Write(Form("resolDG_fit%d_%d",j,k),TObject::kOverwrite);
	resolSG    ->Write(Form("resolSG_fit%d_%d",j,k),TObject::kOverwrite);
      }
      if(k==1){
	resolDG = new TF1("resolDG_fit", DG_formula_Bin1 , -800,800);
	resolDG->SetNpx(1000);
	resolDG->SetLineColor(kBlue);

	float norm = resolDG_old->Integral(-800,800);
	float m1 = meanG1HeavyBin1->Eval(en);
	float m2 = meanG2HeavyBin1->Eval(en);
	float s1 = sigmaG1HeavyBin1->Eval(en);
	float s2 = sigmaG2HeavyBin1->Eval(en);
	resolDG->SetParameter(0, norm/(relWeightBin1*sqrt(2*TMath::Pi())*s1 + (1-relWeightBin1)*sqrt(2*TMath::Pi())*s2 )  );
	resolDG->SetParameter(1, m1);
	resolDG->SetParameter(2, s1);
	resolDG->SetParameter(3, s2);
	resolDG->SetParameter(4, m2);

	resolSG = new TF1("resolSG_fit", SG_formula , -800, 800);
	resolSG->SetNpx(1000);
	resolSG->SetLineColor(kMagenta);	
	float m = meanHeavyBin1->Eval(en);
	float s = sigmaHeavyBin1->Eval(en);
	resolSG->SetParameter(0, norm/sqrt(2*TMath::Pi())/s );
	resolSG->SetParameter(1, m);
	resolSG->SetParameter(2, s);
	
	fout->cd( Form("BinHeavy%d_%d",j, k) );
	resolDG_old->Write(Form("resolDG%d_%d",    j,k),TObject::kOverwrite);
	resolDG    ->Write(Form("resolDG_fit%d_%d",j,k),TObject::kOverwrite);
	resolSG    ->Write(Form("resolSG_fit%d_%d",j,k),TObject::kOverwrite);

      }
 

    }

    TH2F* hCorrLight = new TH2F("hCorrLight","",  nBinsY, 30,300 , 100, 0,500);
    for(int j = 1; j<=hCorrLight->GetNbinsX(); j++){

      float en = hCorrLight->GetBinCenter(j);

      TH1F* hCorrLight_py = (TH1F*)fout->Get(Form("BinLight%d_%d/hCorrLight_py",j, k));

      if(!hCorrLight_py){
	cout << "No histo" << endl;
	continue;
      }

      TF1* resolDG_old =  hCorrLight_py->GetFunction("resol");
      if(!resolDG_old){
	cout << "No func" << endl;
	continue;
      }
      TF1* resolSG = 0;      
      if(k==0){

	resolSG = new TF1("resolSG_fit", SG_formula , -800, 800);
	resolSG->SetNpx(1000);
	resolSG->SetLineColor(kMagenta);
	float norm = resolDG_old->Integral(-800,800);	
	float m = meanLightBin0->Eval(en);
	float s = sigmaLightBin0->Eval(en);
	resolSG->SetParameter(0, norm/sqrt(2*TMath::Pi())/s );
	resolSG->SetParameter(1, m);
	resolSG->SetParameter(2, s);

	fout->cd( Form("BinLight%d_%d",j, k) );
	resolDG_old->Write(Form("resol%d_%d",    j,k),TObject::kOverwrite);
	resolSG    ->Write(Form("resolSG_fit%d_%d",j,k),TObject::kOverwrite);
      }
      if(k==1){
	resolSG = new TF1("resolSG_fit",SG_formula, -800, 800);
	resolSG->SetNpx(1000);
	resolSG->SetLineColor(kMagenta);	
	float norm = resolDG_old->Integral(-800,800);
	float m = meanLightBin1->Eval(en);
	float s = sigmaLightBin1->Eval(en);
	resolSG->SetParameter(0, norm/sqrt(2*TMath::Pi())/s );
	resolSG->SetParameter(1, m);
	resolSG->SetParameter(2, s);

	fout->cd( Form("BinLight%d_%d",j, k) );
	resolDG_old->Write(Form("resol%d_%d",    j,k),TObject::kOverwrite);
	resolSG    ->Write(Form("resolSG_fit%d_%d",j,k),TObject::kOverwrite);

      }

    }

  }
  
  if(doOnlyJetTF){
    fout->cd();  
    resolLightBin0->Write("",TObject::kOverwrite);
    respLightBin0->Write("",TObject::kOverwrite);
    resolLightBin1->Write("",TObject::kOverwrite);
    respLightBin1->Write("",TObject::kOverwrite);
    resolHeavyBin0->Write("",TObject::kOverwrite);
    respHeavyBin0->Write("",TObject::kOverwrite);
    resolHeavyBin1->Write("",TObject::kOverwrite);
    respHeavyBin1->Write("",TObject::kOverwrite);
    respG1HeavyBin0->Write("",TObject::kOverwrite);
    respG2HeavyBin0->Write("",TObject::kOverwrite);
    resolG1HeavyBin0->Write("",TObject::kOverwrite);
    resolG2HeavyBin0->Write("",TObject::kOverwrite);
    respG1HeavyBin1->Write("",TObject::kOverwrite);
    respG2HeavyBin1->Write("",TObject::kOverwrite);
    resolG1HeavyBin1->Write("",TObject::kOverwrite);
    resolG2HeavyBin1->Write("",TObject::kOverwrite);
    fout->Close();
    delete fout;
    return 0;
  }

  //////////////////////

 
  gSystem->cd("./root/");
  RooWorkspace *w = new RooWorkspace("transferFuntions","all functions for SL ttH");

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







  RooRealVar flavor("flavor","flavor",        -6, 6);
  RooRealVar eta   ("eta_part","eta_part",    -3, 3);

  RooRealVar csvRecoLight("csv_rec","csv", 0., 1.0);  //0.6, 1.0 10 bins
  csvRecoLight.setBins(25);
  RooDataSet datasetBtagLight_Bin0("datasetBtagLight_Bin0","dataset", RooArgSet(csvRecoLight,flavor,eta), Import( *treeJetsLight ), Cut("abs(flavor)!=4 && abs(flavor)!=5 && TMath::Abs(eta_part)<1.0 && csv_rec>-1") ); //0.6
  RooDataSet datasetBtagLight_Bin1("datasetBtagLight_Bin1","dataset", RooArgSet(csvRecoLight,flavor,eta), Import( *treeJetsLight ), Cut("abs(flavor)!=4 && abs(flavor)!=5 && TMath::Abs(eta_part)>1.0 && TMath::Abs(eta_part)<2.5 && csv_rec>-1") );

  RooRealVar csvRecoCharm("csv_rec","csv", 0., 1.0);  //0.6, 1.0 10 bins
  csvRecoCharm.setBins(25);
  RooDataSet datasetBtagCharm_Bin0("datasetBtagCharm_Bin0","dataset", RooArgSet(csvRecoCharm,flavor,eta), Import( *treeJetsLight ), Cut("abs(flavor)==4 && TMath::Abs(eta_part)<1.0 && csv_rec>-1") ); //0.6
  RooDataSet datasetBtagCharm_Bin1("datasetBtagCharm_Bin1","dataset", RooArgSet(csvRecoCharm,flavor,eta), Import( *treeJetsLight ), Cut("abs(flavor)==4 && TMath::Abs(eta_part)>1.0 && TMath::Abs(eta_part)<2.5 && csv_rec>-1") );

  RooKeysPdf pdfCsvLight_Bin0("pdfCsvLight_Bin0","", csvRecoLight, datasetBtagLight_Bin0, RooKeysPdf::NoMirror, 1.); //1.5
  RooKeysPdf pdfCsvLight_Bin1("pdfCsvLight_Bin1","", csvRecoLight, datasetBtagLight_Bin1, RooKeysPdf::NoMirror, 1.);

  w->import( pdfCsvLight_Bin0 );
  w->import( pdfCsvLight_Bin1 );

  RooKeysPdf pdfCsvCharm_Bin0("pdfCsvCharm_Bin0","", csvRecoCharm, datasetBtagCharm_Bin0, RooKeysPdf::NoMirror, 1.);
  RooKeysPdf pdfCsvCharm_Bin1("pdfCsvCharm_Bin1","", csvRecoCharm, datasetBtagCharm_Bin1, RooKeysPdf::NoMirror, 1.);

  w->import( pdfCsvCharm_Bin0 );
  w->import( pdfCsvCharm_Bin1 );
  

  RooRealVar csvRecoHeavy("csv_rec","csv", 0., 1.0); //0.679, 1.0  30 bins
  csvRecoHeavy.setBins(100);
  RooDataSet datasetBtagHeavy_Bin0("datasetBtagHeavy_Bin0","dataset", RooArgSet(csvRecoHeavy,flavor,eta), Import( *treeJetsHeavy ), Cut("TMath::Abs(eta_part)<1.0 && csv_rec>-1") ); //0.679
  RooDataSet datasetBtagHeavy_Bin1("datasetBtagHeavy_Bin1","dataset", RooArgSet(csvRecoHeavy,flavor,eta), Import( *treeJetsHeavy ), Cut("TMath::Abs(eta_part)>1.0 && TMath::Abs(eta_part)<2.5 && csv_rec>-1") );
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





  w->writeToFile((outFileName+string(extraname.Data())+".root").c_str(),kTRUE);

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

  resolLightBin0->Write("",TObject::kOverwrite);
  respLightBin0->Write("",TObject::kOverwrite);
  resolLightBin1->Write("",TObject::kOverwrite);
  respLightBin1->Write("",TObject::kOverwrite);

  resolHeavyBin0->Write("",TObject::kOverwrite);
  respHeavyBin0->Write("",TObject::kOverwrite);
  resolHeavyBin1->Write("",TObject::kOverwrite);
  respHeavyBin1->Write("",TObject::kOverwrite);

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
  
  fout->Close();
  delete fout;


  return 0;
}
