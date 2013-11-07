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
#include "TObjArray.h"
#include "Bianchi/TTHStudies/interface/Samples.h"
//#include "Bianchi/TTHStudies/interface/Test.h"

#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "PhysicsTools/FWLite/interface/TFileService.h"

#define VERBOSE 1

using namespace std;


struct sorterByIntegral {
  bool operator() (float i,float j) const { return (i<j);}
};


void ClearAllHisto(map<string, TH1F*> aMap){

  for(  std::map<string,TH1F*>::iterator it = aMap.begin(); it!= aMap.end() ; it++){
    delete (it->second);
  }

  return;

}


void CanvasAndLegend(TCanvas* c1,  TLegend* leg, int logy=0){

  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  c1->SetLogy(logy);

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

  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.03);
  
}


void DrawHistogramMC(TTree* tree = 0, 
		     TString variable = "", 
		     float& normalization      = *(new float()), 
		     float& normalizationError = *(new float()), 
		     float scaleFactor = 0., 
		     TH1F* h = 0, 
		     TCut cut = TCut(""),
		     int verbose = 0 ){
 
  if(tree!=0 && h!=0){
    h->Reset();
    tree->Draw(variable+">>"+TString(h->GetName()),"(PUweight*weightTrig2012)"*cut); //lheWeight
    h->Scale(scaleFactor);
    normalization      = h->Integral();
    normalizationError = h->GetEntries()>0 ? TMath::Sqrt(h->GetEntries())*(normalization/h->GetEntries()) : 0.0;
    if(verbose==0) h->Reset();
  }
  else{
    cout << "Function DrawHistogramMC has raised an error" << endl;
    return;
  }
}

void DrawHistogramData(TTree* tree = 0, 
		       TString variable = "", 
		       float& normalization      = *(new float()), 
		       float& normalizationError = *(new float()), 
		       float scaleFactor = 0., 
		       TH1F* h = 0, 
		       TCut cut = TCut(""),
		       int verbose = 0 ){
  
  if(tree!=0 && h!=0){
    h->Reset();
    tree->Draw(variable+">>"+TString(h->GetName()),cut);
    h->Scale(scaleFactor);
    normalization      = h->Integral();
    normalizationError = h->GetEntries()>0 ? TMath::Sqrt(h->GetEntries())*(normalization/h->GetEntries()) : 0.0;
    if(verbose==0) h->Reset();
  }
  else{
    cout << "Function DrawHistogramData has raised an error" << endl;
    return;
  }
}




int main(int argc, const char* argv[])
{

  std::cout << "PlotZmm" << std::endl;
  gROOT->SetBatch(true);
 
  gSystem->Load("libFWCoreFWLite");
  gSystem->Load("libDataFormatsFWLite");
  //gSystem->Load("../interface/Test_cxx.so");
  AutoLibraryLoader::enable();

  PythonProcessDesc builder(argv[1]);
 
  const edm::ParameterSet& in = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteInput");

  const edm::VParameterSet& samples = in.getParameter<edm::VParameterSet>("samples") ;
  const edm::VParameterSet& plots   = in.getParameter<edm::VParameterSet>("plots") ;

  std::string pathToFile( in.getParameter<std::string>("pathToFile" ) );
  std::string ordering(   in.getParameter<std::string>("ordering" ) );
  double lumi(       in.getParameter<double>("lumi") );
  bool debug(        in.getParameter<bool>("debug") );

  bool openAllFiles = true;
  Samples* mySamples = new Samples(openAllFiles, pathToFile, ordering, samples, lumi, VERBOSE);
  map<string, TH1F*> mapHist;
  vector<string> mySampleFiles;

  if(mySamples->IsOk()){

    cout << "Ok!" << endl;
    mySampleFiles = mySamples->Files();

    for( unsigned int i = 0 ; i < mySampleFiles.size(); i++){
      string sampleName       = mySampleFiles[i];

      if(VERBOSE){
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




  BOOST_FOREACH(edm::ParameterSet pset, plots){

    mapHist.clear();

    float totalMC   = 0.;
    float totalMCError2 = 0.;
    float totalData = 0.;
    float normalization      = 0;
    float normalizationError = 0;

    bool skip        = pset.getParameter<bool>("skip");
    if(skip) continue;

    float xLow       = pset.getParameter<double>("xLow");
    float xHigh      = pset.getParameter<double>("xHigh");
    int nBins        = pset.getParameter<int>("nBins");
    string var       = pset.getParameter<string>("variable");
    string xTitle    = pset.getParameter<string>("xTitle");
    string yTitle    = pset.getParameter<string>("yTitle");
    string histoName = pset.getParameter<string>("histoName");
    string cut       = pset.getParameter<string>("cut");
    string cutName   = pset.exists("cutName") ? pset.getParameter<string>("cutName") : "";
    int normalize    = pset.exists("normalize") ? pset.getParameter<int>("normalize") : 0;
    int logy         = pset.getParameter<int>("logy");

    for(unsigned int i = 0 ; i < mySampleFiles.size(); i++){

      string currentName       = mySampleFiles[i];

      string header = string( Form("CMS Preliminary #sqrt{s}=8 TeV  L=%.1f fb^{-1} ",lumi) );

      TH1F* h1 = new TH1F( ("h1_"+currentName).c_str(), (header+" "+cutName+"; "+xTitle+";"+yTitle).c_str(), nBins, xLow, xHigh );
      h1->SetLineColor( mySamples->GetColor(currentName) );
      h1->SetFillColor( mySamples->GetColor(currentName) );
      h1->GetXaxis()->SetLabelFont(63); //font in pixels
      h1->GetXaxis()->SetLabelSize(16); //in pixels
      h1->GetYaxis()->SetLabelFont(63); //font in pixels
      h1->GetYaxis()->SetLabelSize(16); //in pixels
      h1->GetYaxis()->SetTitleSize(0.05);
      h1->GetXaxis()->SetTitleSize(0.05);
      h1->GetYaxis()->SetTitleOffset(1.0);

      if( currentName.find("Data")!=string::npos) 
	h1->SetMarkerStyle(20);h1->SetMarkerSize(1.2);h1->SetMarkerColor(kBlack);h1->SetLineColor(kBlack);

      mapHist[currentName] = h1;


      TH1F*  currentHisto      = mapHist[currentName];
      TTree* currentTree       = mySamples->GetTree( currentName, "tree");
      
      if( currentHisto==0 || currentTree==0){ cout << "No histo/file..." << endl; continue; }
      
      TString variable(var.c_str());
     
      float scaleFactor        = mySamples->GetWeight(currentName);
      TCut sampleCut           = (mySamples->GetCut(currentName)).find("DUMMY")==string::npos ?
	TCut((mySamples->GetCut(currentName)).c_str()) : "H.HiggsFlag==1";

      TCut myCut( cut.c_str() );
      myCut = myCut&&sampleCut;
      
      if( currentName.find("Data")==string::npos ){
	DrawHistogramMC(    currentTree, variable, normalization, normalizationError, scaleFactor, currentHisto, myCut, 1);
	totalMC += normalization;
	totalMCError2 += (normalizationError*normalizationError);
	if(VERBOSE) cout << "Histo for " << currentName << " has integral " << currentHisto->Integral() << endl;
      }
      else{
	DrawHistogramData(  currentTree, variable, normalization, normalizationError, scaleFactor, currentHisto, myCut, 1);
	totalData += normalization;
	if(VERBOSE) cout << "Histo for " << currentName << " has integral " << currentHisto->Integral() << endl;
      }
      
    }
    
    
    if(VERBOSE) cout << "Toatl MC = "   << totalMC   << " +/- " << TMath::Sqrt(totalMCError2) << endl;
    if(VERBOSE) cout << "Toatl Data = " << totalData << " +/- " << TMath::Sqrt(totalData)     << endl;
    
    
    TCanvas *c1  = new TCanvas("c1","",5,30,650,600);

    //LEFT : 
    //TLegend* leg = new TLegend(0.119195,0.5611888,0.369969,0.8968531,NULL,"brNDC");
    //RIGHT: 
    //TLegend* leg = new TLegend(0.6393189,0.5594406,0.8900929,0.8951049,NULL,"brNDC");

    TLegend* leg = new TLegend(0.113,0.7353,0.874,0.9001,NULL,"brNDC");
    leg->SetNColumns(5);

    
    CanvasAndLegend(c1, leg, logy);
   
    THStack* aStack = new THStack("aStack","");
    vector<TH1F*> mcHistoArray;

    TH1F* Data = 0;
    TH1F* DiBoson = 0; TH1F* WJets = 0; TH1F* DYJets = 0;
    TH1F* top   = 0;   TH1F* singleTop = 0;
    TH1F* topLF = 0;   TH1F* topC      = 0; TH1F* topB = 0;
    TH1F* Higgs = 0;
    TH1F* TTV = 0;
     
    for(  std::map<string,TH1F*>::iterator it = mapHist.begin(); it!= mapHist.end() ; it++){
      if( (it->first).find("Data")==string::npos  ){
	if( (it->first).find("ZH")==string::npos ||
	    ((it->first).find("ZH")!=string::npos && (it->first).find("125")!=string::npos )  ){

	  if(!(it->second)){
	    cout << "Null pointer... skipping " << it->first <<  endl;
	    continue;
	  }

	  bool isHiggs     = ((it->first).find("WH")!=string::npos  || 
			      (it->first).find("ZH")!=string::npos || 
			      (it->first).find("TTH")!=string::npos );

	  bool isTTV       = ((it->first).find("TTZ")!=string::npos  || 
			      (it->first).find("TTW")!=string::npos );

	  bool isWJets     =  (it->first).find("WJets")!=string::npos;

	  bool isDYJets    =  (it->first).find("DYJets")!=string::npos;

	  bool isDiBoson   = ((it->first).find("ZZ")!=string::npos || 
			      (it->first).find("WW")!=string::npos || 
			      (it->first).find("WZ")!=string::npos);

	  bool isSingleTop = ((it->first).find("Tt")!=string::npos  || 
			      (it->first).find("TtW")!=string::npos || 
			      (it->first).find("Ts")!=string::npos  ||
			      (it->first).find("Tbar")!=string::npos );

	  bool isTTbar     = ((it->first).find("TTJets")   !=string::npos);
	  bool isTTbarLF   = ((it->first).find("TTJets")   !=string::npos && (it->first).find("_LF")!=string::npos);
	  bool isTTbarC    = ((it->first).find("TTJets")   !=string::npos && (it->first).find("_C") !=string::npos);
	  bool isTTbarB    = ((it->first).find("TTJets")   !=string::npos && (it->first).find("_B") !=string::npos);

	  if(isTTbar){
	    if( !top && !isTTbarLF && !isTTbarC && !isTTbarB)    top = (it->second);
	    else if(top && !isTTbarLF && !isTTbarC && !isTTbarB) top->Add(it->second);
	    else if( isTTbarLF ){
	      if(!topLF) topLF = (it->second);
	      else topLF->Add(it->second);
	    }
	    else if( isTTbarC ){
	      if(!topC) topC = (it->second);
	      else topC->Add(it->second);
	    }
	    else if( isTTbarB ){
	      if(!topB) topB = (it->second);
	      else topB->Add(it->second);
	    }
	    else{}
	  }
	  else if( isSingleTop){
	    if(!singleTop) singleTop = (it->second);
	    else singleTop->Add(it->second);
	  }
	  else if(isDiBoson){
	    if(!DiBoson) DiBoson = (it->second);
	    else DiBoson->Add(it->second);
	  }
	  else if(isWJets){
	    if(!WJets) WJets = (it->second);
	    else WJets->Add(it->second);
	  }
	  else if(isDYJets){
	    if(!DYJets) DYJets = (it->second);
	    else DYJets->Add(it->second);
	  }
	  else if(isTTV){
	    if(!TTV) TTV = (it->second);
	    else TTV->Add(it->second);
	  }
	  else if(isHiggs){
	    if(!Higgs) Higgs = (it->second);
	    else Higgs->Add(it->second);
	  }
	  else{
	    mcHistoArray.push_back((it->second));
	    leg->AddEntry(it->second, (it->first).c_str(), "F");
	  }

	  //aStack->Add((it->second));
	  //mcHistoArray.Add((it->second));

	}
      }
      else{
	if(!Data)  Data = (it->second);
	else Data->Add(it->second);
      }

    }

    if(top && (!topLF || !topC || !topB)){
      mcHistoArray.push_back(top);
      leg->AddEntry(top, "Top", "F");
    }
    else if( !top && topLF &&  topC && topB){
      //cout << "LF " << topLF->Integral() << ", C " << topC->Integral() << ", B " << topB->Integral() << endl;
      mcHistoArray.push_back(topLF);
      mcHistoArray.push_back(topC);
      mcHistoArray.push_back(topB);
      leg->AddEntry(topLF, "t#bar{t}+LF", "F");
      leg->AddEntry(topC,  "t#bar{t}+C", "F");
      leg->AddEntry(topB,  "t#bar{t}+B", "F");
    }
    else{}
    if(singleTop){
      mcHistoArray.push_back(singleTop);
      leg->AddEntry(singleTop, "Single top", "F");
    }
    if(DiBoson){
      mcHistoArray.push_back(DiBoson);
      leg->AddEntry(DiBoson, "Di-boson", "F");
    }
    if(WJets || DYJets){
      if(WJets && !DYJets){
	mcHistoArray.push_back(WJets);
	leg->AddEntry(WJets, "W+jets", "F");
      }
      else if(!WJets && DYJets){
	mcHistoArray.push_back(DYJets);
	leg->AddEntry(DYJets, "Z+jets", "F");
      }
      else{
	WJets->Add( DYJets );
	mcHistoArray.push_back(WJets);
	leg->AddEntry(WJets, "EWK", "F");
      }
    }
    if(TTV){
      mcHistoArray.push_back(TTV);
      leg->AddEntry(TTV, "t#bar{t}+V", "F");
    }
    if(Higgs){
      mcHistoArray.push_back(Higgs);
      Higgs->SetFillStyle(3004);
      leg->AddEntry(Higgs, "H(125)", "F");
    }


    std::map<float, TH1F*, sorterByIntegral> hMapIntegral;
    for( unsigned int k = 0 ; k < mcHistoArray.size() ; k++){
      hMapIntegral[ mcHistoArray[k]->Integral() ] = mcHistoArray[k];
     }
    
    TH1F* mcTotal = 0;
    for(std::map<float, TH1F*>::iterator it = hMapIntegral.begin(); it!=hMapIntegral.end(); it++){
      if( mcTotal==0 ) mcTotal = (TH1F*)(it->second)->Clone("mcTotal");
      else  mcTotal->Add(it->second);
      aStack->Add( it->second );
    }


    c1->cd();
    if(Data!=0){
      leg->AddEntry(Data, "Data", "PE");

      //if(mcTotal)
      //Data->SetMaximum( 1.3*TMath::Max(mcTotal->GetMaximum(), Data->GetMaximum()) );
      //else
      //Data->SetMaximum( 1.3*Data->GetMaximum() );

      Data->Sumw2();
      if(normalize){
	Data->DrawNormalized("P");
	for( int k = 0 ; k < (aStack->GetStack())->GetEntries() ; k++){
	  if(totalMC>0) ((TH1F*)(aStack->GetStack())->At(k))->Scale(1./totalMC);
	}
	aStack->Draw("HISTSAME");
	Data->DrawNormalized("PSAME");
      }
      else{

	TPad *pad1 = new TPad("pad1","pad1",0,0.3,1,1);
	pad1->SetBottomMargin(0);
	pad1->Draw();
	pad1->SetLogy(logy);
	pad1->cd();

	TH1F* HiggsScaled = 0;
	if(Higgs && mcTotal){
	  HiggsScaled = (TH1F*)Higgs->Clone("HiggsScaled");
	  HiggsScaled->SetLineWidth(2);
	  HiggsScaled->SetLineStyle(kDotted);
	  HiggsScaled->SetLineColor(kRed);
	  HiggsScaled->SetFillStyle(0);
	  HiggsScaled->Scale((mcTotal->Integral()-Higgs->Integral())/HiggsScaled->Integral());
	  Data->SetMaximum( 1.3*TMath::Max(mcTotal->GetMaximum(), TMath::Max(Data->GetMaximum(),HiggsScaled->GetMaximum() )) );
	}

	Data->Draw("PE");
	aStack->Draw("HISTSAME");
	Data->Draw("PSAME");
	if(HiggsScaled){
	  HiggsScaled->Draw("HISTSAME");
	  leg->AddEntry(HiggsScaled,"Signal norm.to bkg","L");
	}
	leg->Draw();

	c1->cd();
	TPad *pad2 = new TPad("pad2","pad2",0,0.05,1,0.298);
	pad2->SetTitle("");
	pad2->SetTopMargin(0);
	pad2->Draw();
	pad2->cd();

	if(mcTotal){
	  mcTotal->Sumw2();
	  mcTotal->SetTitle("");
	  TH1F* dataMC = (TH1F*)mcTotal->Clone("dataMC");
	  dataMC->SetMarkerStyle(20);dataMC->SetMarkerSize(1.2);dataMC->SetMarkerColor(kBlack);dataMC->SetLineColor(kBlack);
	  dataMC->GetXaxis()->SetLabelFont(63); //font in pixels
	  dataMC->GetXaxis()->SetLabelSize(16); //in pixels
	  dataMC->GetYaxis()->SetLabelFont(63); //font in pixels
	  dataMC->GetYaxis()->SetLabelSize(16); //in pixels
	  dataMC->GetYaxis()->SetTitleSize(0.10);
	  dataMC->GetXaxis()->SetTitleSize(0.13);
	  dataMC->GetYaxis()->SetTitleOffset(0.9);
	  dataMC->GetYaxis()->SetTitle("MC/data");
	  dataMC->Reset();
	  dataMC->Divide(Data, mcTotal, 1.0, 1.0);
	  dataMC->SetMinimum(0);  dataMC->SetMaximum(2);
	  dataMC->Draw("EP");

	  TF1* line = new TF1("line","1",dataMC->GetXaxis()->GetXmin(),dataMC->GetXaxis()->GetXmax());
	  line->SetLineStyle(3);
	  line->SetLineWidth(1.5);
	  line->SetLineColor(kBlack);
	  line->Draw("SAME");
	}
	c1->cd();

      }
    }
    else{
      if(normalize){
	for( int k = 0 ; k < (aStack->GetStack())->GetEntries() ; k++){
	  if(totalMC>0) ((TH1F*)(aStack->GetStack())->At(k))->Scale(1./totalMC);
	}
	aStack->Draw("HIST");
	leg->Draw();;
      }
      else{
	aStack->Draw("HIST");
	leg->Draw();
      }
      TH1F* hStack = (TH1F*)aStack->GetHistogram();
      hStack->SetTitle(cutName.c_str());
      hStack->SetXTitle(xTitle.c_str());
      hStack->SetYTitle(yTitle.c_str());
      hStack->SetTitleSize(  0.04,"X");
      hStack->SetTitleSize(  0.05,"Y");
      hStack->SetTitleOffset(0.95,"Y");
    }

    c1->SaveAs(("Plots/"+histoName+".png").c_str());
    c1->SaveAs(("Plots/"+histoName+".root").c_str());
   

    ///////////////////  clear /////////////////////////////////
    ClearAllHisto(mapHist); delete c1; delete leg; delete aStack;

    if(debug) return 0;

  }


  return 0;

}
