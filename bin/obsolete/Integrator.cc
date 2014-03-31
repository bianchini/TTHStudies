#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include <cstdlib>
#include <iostream> 
#include <fstream>
#include <map>
#include <string>

#include "TChain.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include "TH1F.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TLegend.h"
#include "THStack.h"
#include "TCut.h"
#include "TArrayF.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TMultiGraph.h"

#include "TF3.h"
#include "Math/WrappedMultiTF1.h"
#include "Math/AdaptiveIntegratorMultiDim.h"
#include "Math/AllIntegrationTypes.h"

#include <Math/IFunction.h>
#include <Math/IFunctionfwd.h>
#include "Math/IntegratorMultiDim.h"
#include <Math/Integrator.h>

#include "PhysicsTools/FWLite/interface/TFileService.h"


using namespace std;



class Integrand : public ROOT::Math::IMultiGenFunction {
  
public:

  Integrand(double, double);
  ~Integrand();

  ROOT::Math::IMultiGenFunction* Clone() const { return new Integrand(14,0.030*0.030); }
  double DoEval(const double *) const;
  unsigned int NDim() const { return 2;}


private:
  double sqrtS_;
  double pt2_;
 
  mutable long counter_;
};


Integrand::Integrand(double sqrtS, double pt2): sqrtS_(sqrtS), pt2_(pt2) , counter_(0) {}

Integrand::~Integrand(){
  cout << "Delete!" << endl;
}


double Integrand::DoEval(const double * x) const {

  counter_++;

  if( 1 - 4*pt2_/(x[0]*x[1]*sqrtS_*sqrtS_) < 0 ){
    return 0.0;
  }

  double c = TMath::Sqrt( 1 - 4*pt2_/(x[0]*x[1]*sqrtS_*sqrtS_));
  if(c<0){
    cout << "negative c!!!" << endl;
    return 0.0;
  }

  double matrixElement = 9./2*(3 - (1-c*c)/4 + 2*(1+c)/(1-c)/(1-c) + 2*(1-c)/(1+c)/(1+c)   );
  if(matrixElement<0){
    cout << "negative matrixElement!!!" << endl;
    return 0.0;
  }

  double flux          = 1./(x[0]*x[1]*x[0]*x[1]*sqrtS_*sqrtS_*sqrtS_*sqrtS_)/c;

  double gluon1        = 1./x[0]*(TMath::Power(x[0],-0.26))*TMath::Power(1-x[0],5.2)*(1-0.52*x[0]); 
  double gluon2        = 1./x[1]*(TMath::Power(x[1],-0.26))*TMath::Power(1-x[1],5.2)*(1-0.52*x[1]);

  //double gluon1 = 1.0;
  //double gluon2 = 1.0;

  double result        =  matrixElement*flux*gluon1*gluon2*2*TMath::Sqrt(pt2_);

  return result;
}


double myIntegrator(double sqrtS, double pt2){

  Integrand* myIntegrand = new Integrand(sqrtS, pt2);

  ROOT::Math::IntegratorMultiDim* integrator_ = new ROOT::Math::IntegratorMultiDim(*myIntegrand);

  integrator_->SetRelTolerance(1e-06);
  integrator_->SetFunction(*myIntegrand);

  double xMin[2] = {0.0001, 0.0001};
  double xMax[2] = {0.99, 0.99};
  double integral  = integrator_->Integral(xMin,xMax);

  cout << integral << endl;

  delete myIntegrand;

  return integral; 

}

void plot(){

  TCanvas *c1 = new TCanvas("c1","gg->gg cross-section @ 8TeV ; p_T (GeV) ; A.U.",200,10,700,500);
  c1->SetLogy();
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  gStyle->SetOptStat(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleH(0.07);
  gStyle->SetTitleFontSize(0.1);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleOffset(1.3,"y");

  c1->cd();


  float x[30] ; float y[30] ;

  for(unsigned int k = 0 ; k < 30 ; k++){
    x[k] = 0.020 + 0.030 * k;
    y[k] = myIntegrator( 8, x[k]*x[k] );
    x[k] *= 1000;
  }


  TH1F* h1 = new TH1F("h1","gg->gg cross-section @ 8TeV ; p_T (GeV) ; A.U.",30, x[0]-5, x[29]+5 );
  TGraph* graph =  new TGraph(30,x, y);

  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("gg->gg cross-section @ 8TeV");

  //h1->Draw();
  graph->SetMarkerStyle(kFullCircle);
  graph->SetMarkerSize(1.5);
  graph->SetMarkerColor(kRed);

  mg->Add(graph);
  mg->Draw("ALP3");

  c1->cd();
  gPad->Modified();
  mg->GetXaxis()->SetLimits(x[0]-5, x[29]+5);
  mg->GetXaxis()->SetTitle("p_{T} (GeV)");
  mg->GetYaxis()->SetTitle("A.U.");

  c1->SaveAs("graph.png");
  c1->SaveAs("graph.root");
}



int main(int argc, const char* argv[])
{

  std::cout << "integrator of gg cross-section" << std::endl;
  gROOT->SetBatch(true);
 
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

  //myIntegrator(  14,0.20*0.20 );

  plot();

  return 0;

}
