#include <cstdlib>
#include <iostream> 
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include "TMath.h"
#include "TMatrixT.h"
#include "TMatrixTBase.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TF2.h"

#include "TPaveText.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TKey.h"
#include "TList.h"
#include "TCollection.h"
#include "TObject.h"
#include "TLegend.h"
#include "THStack.h"
#include "TCut.h"
#include "TGraphAsymmErrors.h"
#include "TGraphPainter.h"
#include "TMultiGraph.h"
#include "TArrayF.h"
#include "TLine.h"
#include "TLorentzVector.h"

#include "TIterator.h"

#include "TRandom3.h"

#define RUNONDATA 1

#define SPLITCC   1

#define KSWOSYSTEMATICS 1

string version   =  "_all_rec_std" ;
string inputpath = "files/byLLR/Apr07_2014/";
string tag       = "New_all_rec_std";


void plot_param_3(int isB = 1, TString header = "Light jet TF", TString CP = "TEST"){

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


  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TString SG_formula("[0]*exp(-0.5*(x-[1])*(x-[1])/[2]/[2])");
  TString DG_formula(Form("[0]*( %f*exp(-0.5*(x-[1])*(x-[1])/[2]/[2]) + (1-%f)*exp(-0.5*(x-[4])*(x-[4])/[3]/[3]))",           0.65, 0.65));

  TF1* l_Bin0_1 = new TF1("l_Bin0_1", SG_formula , 0, 400); l_Bin0_1->SetLineWidth(2) ; l_Bin0_1->SetLineColor(kMagenta); l_Bin0_1->SetLineStyle(kSolid);
  TF1* l_Bin0_2 = new TF1("l_Bin0_2", SG_formula , 0, 400); l_Bin0_2->SetLineWidth(2) ; l_Bin0_2->SetLineColor(kRed);     l_Bin0_2->SetLineStyle(kSolid);
  TF1* l_Bin0_3 = new TF1("l_Bin0_3", SG_formula , 0, 400); l_Bin0_3->SetLineWidth(2) ; l_Bin0_3->SetLineColor(kBlue);    l_Bin0_3->SetLineStyle(kSolid);

  TF1* l_Bin1_1 = new TF1("l_Bin1_1", SG_formula , 0, 400); l_Bin1_1->SetLineWidth(2) ; l_Bin1_1->SetLineColor(kMagenta); l_Bin1_1->SetLineStyle(kDashed);
  TF1* l_Bin1_2 = new TF1("l_Bin1_2", SG_formula , 0, 400); l_Bin1_2->SetLineWidth(2) ; l_Bin1_2->SetLineColor(kRed);     l_Bin1_2->SetLineStyle(kDashed);
  TF1* l_Bin1_3 = new TF1("l_Bin1_3", SG_formula , 0, 400); l_Bin1_3->SetLineWidth(2) ; l_Bin1_3->SetLineColor(kBlue);    l_Bin1_3->SetLineStyle(kDashed);
  
  TF1* b_Bin0_1 = new TF1("b_Bin0_1", DG_formula , 0, 400); b_Bin0_1->SetLineWidth(2) ; b_Bin0_1->SetLineColor(kMagenta); b_Bin0_1->SetLineStyle(kSolid);
  TF1* b_Bin0_2 = new TF1("b_Bin0_2", DG_formula , 0, 400); b_Bin0_2->SetLineWidth(2) ; b_Bin0_2->SetLineColor(kRed);     b_Bin0_2->SetLineStyle(kSolid);
  TF1* b_Bin0_3 = new TF1("b_Bin0_3", DG_formula , 0, 400); b_Bin0_3->SetLineWidth(2) ; b_Bin0_3->SetLineColor(kBlue);    b_Bin0_3->SetLineStyle(kSolid);

  TF1* b_Bin1_1 = new TF1("b_Bin1_1", DG_formula , 0, 400); b_Bin1_1->SetLineWidth(2) ; b_Bin1_1->SetLineColor(kMagenta); b_Bin1_1->SetLineStyle(kDashed);
  TF1* b_Bin1_2 = new TF1("b_Bin1_2", DG_formula , 0, 400); b_Bin1_2->SetLineWidth(2) ; b_Bin1_2->SetLineColor(kRed);     b_Bin1_2->SetLineStyle(kDashed);
  TF1* b_Bin1_3 = new TF1("b_Bin1_3", DG_formula , 0, 400); b_Bin1_3->SetLineWidth(2) ; b_Bin1_3->SetLineColor(kBlue);    b_Bin1_3->SetLineStyle(kDashed);
 

  TH1F* h = new TH1F("h","CMS Simulation #sqrt{s}=8 TeV ; jet energy (GeV) ; transfer function (1/GeV)",1,5,225);
  h->GetYaxis()->SetRangeUser( 0, 0.07);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetTitleSize(0.05);

  h->GetYaxis()->SetTitleOffset(0.95);
  h->GetXaxis()->SetTitleOffset(0.90);

  TFile *f = TFile::Open("../bin/root/ControlPlots"+CP+".root","READ");


  TH1F* hBestMass = (TH1F*)f->Get("resolLightBin0");
  if( !hBestMass ){
    cout << "No such histogram" << endl;
    return;
  }
  TString fname;
  TIter iter( hBestMass->GetListOfFunctions() );
  TObject* obj = iter();
  while ( obj!=0 ){
    TF1* func = dynamic_cast<TF1*>(obj);
    fname = func->GetName();
    obj = iter();
  }
  TF1 *func = hBestMass->GetFunction(fname);
  l_Bin0_1->SetParameter(2, func->Eval(50.) ); l_Bin0_2->SetParameter(2, func->Eval(100.) ); l_Bin0_3->SetParameter(2, func->Eval(150.) );
  
  hBestMass = (TH1F*)f->Get("resolLightBin1");
  if( !hBestMass ){
    cout << "No such histogram" << endl;
    return;
  }
  TIter iter2( hBestMass->GetListOfFunctions() );
  obj = iter2();
  while ( obj!=0 ){
    TF1* func = dynamic_cast<TF1*>(obj);
    fname = func->GetName();
    obj = iter2();
  }
  func = hBestMass->GetFunction(fname);
  l_Bin1_1->SetParameter(2, func->Eval(50.) ); l_Bin1_2->SetParameter(2, func->Eval(100.) ); l_Bin1_3->SetParameter(2, func->Eval(150.) );

  hBestMass = (TH1F*)f->Get("respLightBin0");
  if( !hBestMass ){
    cout << "No such histogram" << endl;
    return;
  }
  TIter iter3( hBestMass->GetListOfFunctions() );
  obj = iter3();
  while ( obj!=0 ){
    TF1* func = dynamic_cast<TF1*>(obj);
    fname = func->GetName();
    obj = iter3();
  }
  func = hBestMass->GetFunction(fname);
  l_Bin0_1->SetParameter(1, func->Eval(50.) ); l_Bin0_2->SetParameter(1, func->Eval(100.) ); l_Bin0_3->SetParameter(1, func->Eval(150.) );
  
  hBestMass = (TH1F*)f->Get("respLightBin1");
  if( !hBestMass ){
    cout << "No such histogram" << endl;
    return;
  }
  TIter iter4( hBestMass->GetListOfFunctions() );
  obj = iter4();
  while ( obj!=0 ){
    TF1* func = dynamic_cast<TF1*>(obj);
    fname = func->GetName();
    obj = iter4();
  }
  func = hBestMass->GetFunction(fname);
  l_Bin1_1->SetParameter(1, func->Eval(50.) ); l_Bin1_2->SetParameter(1, func->Eval(100.) ); l_Bin1_3->SetParameter(1, func->Eval(150.) );

  l_Bin0_1->SetParameter(0,1./TMath::Sqrt( 2* TMath::Pi())/ l_Bin0_1->GetParameter(2) );
  l_Bin0_2->SetParameter(0,1./TMath::Sqrt( 2* TMath::Pi())/ l_Bin0_2->GetParameter(2) );
  l_Bin0_3->SetParameter(0,1./TMath::Sqrt( 2* TMath::Pi())/ l_Bin0_3->GetParameter(2) );

  l_Bin1_1->SetParameter(0,1./TMath::Sqrt( 2* TMath::Pi())/ l_Bin1_1->GetParameter(2) );
  l_Bin1_2->SetParameter(0,1./TMath::Sqrt( 2* TMath::Pi())/ l_Bin1_2->GetParameter(2) );
  l_Bin1_3->SetParameter(0,1./TMath::Sqrt( 2* TMath::Pi())/ l_Bin1_3->GetParameter(2) );


  hBestMass = (TH1F*)f->Get("respG1HeavyBin0");
  TIter iter5( hBestMass->GetListOfFunctions() );
  obj = iter5();
  while ( obj!=0 ){
    TF1* func = dynamic_cast<TF1*>(obj);
    fname = func->GetName();
    obj = iter5();
  }
  func = hBestMass->GetFunction(fname);
  b_Bin0_1->SetParameter(1, func->Eval(50.) ); b_Bin0_2->SetParameter(1, func->Eval(100.) ); b_Bin0_3->SetParameter(1, func->Eval(150.) );

  hBestMass = (TH1F*)f->Get("respG2HeavyBin0");
  TIter iter6( hBestMass->GetListOfFunctions() );
  obj = iter6();
  while ( obj!=0 ){
    TF1* func = dynamic_cast<TF1*>(obj);
    fname = func->GetName();
    obj = iter6();
  }
  func = hBestMass->GetFunction(fname);
  b_Bin0_1->SetParameter(4, func->Eval(50.) ); b_Bin0_2->SetParameter(4, func->Eval(100.) ); b_Bin0_3->SetParameter(4, func->Eval(150.) );

  hBestMass = (TH1F*)f->Get("resolG1HeavyBin0");
  TIter iter7( hBestMass->GetListOfFunctions() );
  obj = iter7();
  while ( obj!=0 ){
    TF1* func = dynamic_cast<TF1*>(obj);
    fname = func->GetName();
    obj = iter7();
  }
  func = hBestMass->GetFunction(fname);
  b_Bin0_1->SetParameter(2, func->Eval(50.) ); b_Bin0_2->SetParameter(2, func->Eval(100.) ); b_Bin0_3->SetParameter(2, func->Eval(150.) );

  hBestMass = (TH1F*)f->Get("resolG2HeavyBin0");
  TIter iter8( hBestMass->GetListOfFunctions() );
  obj = iter8();
  while ( obj!=0 ){
    TF1* func = dynamic_cast<TF1*>(obj);
    fname = func->GetName();
    obj = iter8();
  }
  func = hBestMass->GetFunction(fname);
  b_Bin0_1->SetParameter(3, func->Eval(50.) ); b_Bin0_2->SetParameter(3, func->Eval(100.) ); b_Bin0_3->SetParameter(3, func->Eval(150.) );

  hBestMass = (TH1F*)f->Get("respG1HeavyBin1");
  TIter iter9( hBestMass->GetListOfFunctions() );
  obj = iter9();
  while ( obj!=0 ){
    TF1* func = dynamic_cast<TF1*>(obj);
    fname = func->GetName();
    obj = iter9();
  }
  func = hBestMass->GetFunction(fname);
  b_Bin1_1->SetParameter(1, func->Eval(50.) ); b_Bin1_2->SetParameter(1, func->Eval(100.) ); b_Bin1_3->SetParameter(1, func->Eval(150.) );

  hBestMass = (TH1F*)f->Get("respG2HeavyBin1");
  TIter iter10( hBestMass->GetListOfFunctions() );
  obj = iter10();
  while ( obj!=0 ){
    TF1* func = dynamic_cast<TF1*>(obj);
    fname = func->GetName();
    obj = iter10();
  }
  func = hBestMass->GetFunction(fname);
  b_Bin1_1->SetParameter(4, func->Eval(50.) ); b_Bin1_2->SetParameter(4, func->Eval(100.) ); b_Bin1_3->SetParameter(4, func->Eval(150.) );

  hBestMass = (TH1F*)f->Get("resolG1HeavyBin1");
  TIter iter11( hBestMass->GetListOfFunctions() );
  obj = iter11();
  while ( obj!=0 ){
    TF1* func = dynamic_cast<TF1*>(obj);
    fname = func->GetName();
    obj = iter11();
  }
  func = hBestMass->GetFunction(fname);
  b_Bin1_1->SetParameter(2, func->Eval(50.) ); b_Bin1_2->SetParameter(2, func->Eval(100.) ); b_Bin1_3->SetParameter(2, func->Eval(150.) );

  hBestMass = (TH1F*)f->Get("resolG2HeavyBin1");
  TIter iter12( hBestMass->GetListOfFunctions() );
  obj = iter12();
  while ( obj!=0 ){
    TF1* func = dynamic_cast<TF1*>(obj);
    fname = func->GetName();
    obj = iter12();
  }
  func = hBestMass->GetFunction(fname);
  b_Bin1_1->SetParameter(3, func->Eval(50.) ); b_Bin1_2->SetParameter(3, func->Eval(100.) ); b_Bin1_3->SetParameter(3, func->Eval(150.) );


  b_Bin0_1->SetParameter(0,1./TMath::Sqrt( 2* TMath::Pi())/ ( 0.65*b_Bin0_1->GetParameter(2) + 0.35*b_Bin0_1->GetParameter(3) ) );
  b_Bin0_2->SetParameter(0,1./TMath::Sqrt( 2* TMath::Pi())/ ( 0.65*b_Bin0_2->GetParameter(2) + 0.35*b_Bin0_2->GetParameter(3) ) );
  b_Bin0_3->SetParameter(0,1./TMath::Sqrt( 2* TMath::Pi())/ ( 0.65*b_Bin0_3->GetParameter(2) + 0.35*b_Bin0_3->GetParameter(3) ) );

  b_Bin1_1->SetParameter(0,1./TMath::Sqrt( 2* TMath::Pi())/ ( 0.65*b_Bin1_1->GetParameter(2) + 0.35*b_Bin1_1->GetParameter(3) ) );
  b_Bin1_2->SetParameter(0,1./TMath::Sqrt( 2* TMath::Pi())/ ( 0.65*b_Bin1_2->GetParameter(2) + 0.35*b_Bin1_2->GetParameter(3) ) );
  b_Bin1_3->SetParameter(0,1./TMath::Sqrt( 2* TMath::Pi())/ ( 0.65*b_Bin1_3->GetParameter(2) + 0.35*b_Bin1_3->GetParameter(3) ) );


  TLegend* leg = new TLegend(0.15,0.65,0.55,0.85,NULL,"brNDC");
  leg->SetHeader( header );
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.05); 
  leg->AddEntry(l_Bin0_1 , "E_{q} = 50 GeV" , "L");
  leg->AddEntry(l_Bin0_2 , "E_{q} =100 GeV" , "L");
  leg->AddEntry(l_Bin0_3 , "E_{q} =150 GeV" , "L");

  h->Draw();
  if(isB==0){
  l_Bin0_1->Draw("SAME");
  l_Bin0_2->Draw("SAME");
  l_Bin0_3->Draw("SAME");
  l_Bin1_1->Draw("SAME");
  l_Bin1_2->Draw("SAME");
  l_Bin1_3->Draw("SAME");
  }
  if(isB){
  b_Bin0_1->Draw("SAME");
  b_Bin0_2->Draw("SAME");
  b_Bin0_3->Draw("SAME");
  b_Bin1_1->Draw("SAME");
  b_Bin1_2->Draw("SAME");
  b_Bin1_3->Draw("SAME");
  }

  leg->Draw();

  TPaveText *pt1 = new TPaveText(0.62, 0.45, 0.71, 0.54,"brNDC");
  pt1->SetFillStyle(1001);
  pt1->SetBorderSize(0);
  pt1->SetFillColor(kWhite);
  pt1->SetTextSize(0.04); 
  pt1->AddText("Solid:  |#eta|<1.0");
  pt1->AddText("Dashed: |#eta|>1.0");
  pt1->SetTextAlign(11);
  pt1->Draw();

  if(1){
    //c1->SaveAs("Plots/PAS/Plot_TF_"+CP+"_"+TString(Form("%d",isB))+".png");
    c1->SaveAs("Plots/PAS/Plot_TF_"+CP+"_"+TString(Form("%d",isB))+".pdf");
  }

  delete c1;
  return;
}

void plot_param_3All(){
  
  plot_param_3(1, "b-quarks",   "TEST");
  plot_param_3(1, "b-quarks", "TEST_std_gen");

  plot_param_3(0, "u-quarks",   "TEST");
  plot_param_3(0, "u-quarks", "TEST_std_gen");
    

}




void plot_time(TString sample = "TTH125",
	       TString title = ""
	       ){

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


  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  c1->SetLogx();

  TLegend* leg = new TLegend(0.65,0.65,0.90,0.88,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04); 


  TFile* fN_SL  = TFile::Open(inputpath+"MEAnalysisNew"+version+"_"+sample+".root");
  TTree* tN_SL = (TTree*)fN_SL->Get("tree");

  TFile* fN_DL  = TFile::Open(inputpath+"MEAnalysisNew"+version+"_"+sample+".root");
  TTree* tN_DL = (TTree*)fN_DL->Get("tree");

  TH1F* hN_type0 = new TH1F("hN_type0","CMS Simulation #sqrt{s}=8 TeV "+title+"; CPU time/event (s); normalized to unity",480,20,500);
  hN_type0->SetLineColor(kRed);
  hN_type0->SetLineWidth(3);
  hN_type0->SetFillColor(kRed);
  hN_type0->SetFillStyle(3004);
  hN_type0->SetTitleSize(0.05,"X");
  hN_type0->SetTitleSize(0.05,"Y");
  hN_type0->SetTitleOffset(0.90,"X");
  hN_type0->SetTitleOffset(0.95,"Y");
  tN_SL->Draw("time>>hN_type0","( type==0 || (type==3 && flag_type3>0)) && syst==0");

  TH1F* hN_type1 = new TH1F("hN_type1","CMS Simulation #sqrt{s}=8 TeV, "+title+"; time/event (s); normalized to unity",480,20,500);
  hN_type1->SetLineColor(kBlue);
  hN_type1->SetLineWidth(3);
  hN_type1->SetFillColor(kBlue);
  hN_type1->SetFillStyle(3004);
  tN_SL->Draw("time>>hN_type1","( type==1 || (type==3 && flag_type3<=0) ) && syst==0");

  TH1F* hN_type2 = new TH1F("hN_type2","Simulation #sqrt{s}=8 TeV, "+title+"; time/event (s); normalized to unity",480,20,500);
  hN_type2->SetLineColor(kGreen);
  hN_type2->SetLineWidth(3);
  hN_type2->SetFillColor(kGreen);
  hN_type2->SetFillStyle(3004);
  tN_SL->Draw("time>>hN_type2","type==2 && syst==0");

  TH1F* hN_type6 = new TH1F("hN_type6","CMS Simulation #sqrt{s}=8 TeV, "+title+"; time/event (s); normalized to unity",480,20,500);
  hN_type6->SetLineColor(kMagenta);
  hN_type6->SetLineWidth(3);
  hN_type6->SetFillColor(kMagenta);
  hN_type6->SetFillStyle(3005);
  tN_DL->Draw("time>>hN_type6","type==6 && syst==0");

  float max = hN_type0->GetMaximum()*1.2;

  hN_type0->SetMinimum( 0.0 );
  hN_type0->SetMaximum( max );
  hN_type0->GetXaxis()->SetMoreLogLabels();
  hN_type0->DrawNormalized("HIST");
  hN_type1->DrawNormalized("HISTSAME");
  hN_type2->DrawNormalized("HISTSAME");
  hN_type6->DrawNormalized("HISTSAME");

  leg->AddEntry(hN_type0,"SL-cat1" , "F");
  leg->AddEntry(hN_type1,"SL-cat2",  "F");
  leg->AddEntry(hN_type2,"SL-cat3",  "F");
  leg->AddEntry(hN_type6,"DL",       "F");
  leg->Draw();

  cout << max << endl;

  TLine* line = new TLine(hN_type0->GetBinLowEdge( hN_type0->FindBin(60) ), max/hN_type0->Integral()  , hN_type0->GetBinLowEdge( hN_type0->FindBin(60)), 0.);
  line->SetLineWidth(3);
  line->SetLineStyle(kDashed);
  line->SetLineColor(kBlack);
  line->Draw("SAME");

  TLine* line2 = new TLine(hN_type0->GetBinLowEdge( hN_type0->FindBin(120) ), max/hN_type0->Integral()  , hN_type0->GetBinLowEdge( hN_type0->FindBin(120)), 0.);
  line2->SetLineWidth(3);
  line2->SetLineStyle(kDashed);
  line2->SetLineColor(kBlack);
  line2->Draw("SAME");

  TLine* line3 = new TLine(hN_type0->GetBinLowEdge( hN_type0->FindBin(180) ), max/hN_type0->Integral()  , hN_type0->GetBinLowEdge( hN_type0->FindBin(180)), 0.);
  line3->SetLineWidth(3);
  line3->SetLineStyle(kDashed);
  line3->SetLineColor(kBlack);
  line3->Draw("SAME");

  leg->Draw();

  if(1){
    c1->SaveAs("Plots/PAS/Plot_Time_"+sample+".pdf");
  }

}



void plot_Slice(TString header = "Light jet TF",
		TString unitsX = "jet E (GeV)",
		TString unitsY = "TF(E|#hat{E})",
		int event = 1,
		string flavor = "",
		int bin = 0,
		TString extraLabel = "",
		TString param = ""
		){
  
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


  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TFile *f = TFile::Open("../bin/root/ControlPlotsTEST"+extraLabel+".root","READ");
  TH1F* hBestMass = (TH1F*)f->Get(TString(Form("Bin%s%d_%d/hCorr%s_py",flavor.c_str(), event, bin,flavor.c_str())));

  TF1* funcP  = 0;
  TF1* funcP2 = 0;
  if( flavor.find("DUMMY")==string::npos ){
    funcP  = (TF1*)f->Get(TString(Form("Bin%s%d_%d/resolDG_fit%d_%d",flavor.c_str(), event, bin, event, bin)));
    funcP2 = (TF1*)f->Get(TString(Form("Bin%s%d_%d/resolSG_fit%d_%d",flavor.c_str(), event, bin, event, bin)));
  }
  
  if(!funcP)  cout << "Could not find funcP" << endl;
  if(!funcP2) cout << "Could not find funcP2" << endl;

  if( !hBestMass ){
    cout << "No such histogram" << endl;
    return;
  }

 
  TString fname;
  TIter iter( hBestMass->GetListOfFunctions() );
  TObject* obj = iter();
  while ( obj!=0 ){
    TF1* func = dynamic_cast<TF1*>(obj);
    fname = func->GetName();
    obj = iter();
  }
  TF1 *func = hBestMass->GetFunction(fname);

  hBestMass->SetMarkerStyle(kFullCircle);
  hBestMass->SetMarkerColor(kBlack);
  hBestMass->SetTitle("CMS Simulation #sqrt{s}=8 TeV");
  hBestMass->SetXTitle( unitsX );
  hBestMass->SetYTitle( "a.u." );
  hBestMass->GetYaxis()->SetTitleSize(0.05);
  hBestMass->GetXaxis()->SetTitleSize(0.05);


  TLegend* leg = new TLegend(0.39,0.68,0.73,0.89,NULL,"brNDC");
  leg->SetHeader( header );
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.05); 
  //leg->AddEntry(hBestMass, unitsY);

  float normFunc = 0.;  
  if(func) {
    leg->AddEntry(func, param, "L");
    func->SetNpx(1000);
    normFunc = func->GetParameter(0);
    cout << "Bin: " << func->Integral(-800,800) << endl;
  }
  if(funcP){
    leg->AddEntry(funcP,  "parametrization (2G)", "L");
    funcP->SetNpx(1000);
    funcP->SetParameter(0, normFunc);
    cout << "Param 2G: " << funcP->Integral(-800,800) << endl;
    if( TMath::Abs(funcP->Integral(-800,800)/func->Integral(-800,800)-1)>1e-02  )
      funcP->SetParameter(0, funcP->GetParameter(0)*func->Integral(-800,800)/funcP->Integral(-800,800)  );
    cout << "RESCALE: Param 2G: " << funcP->Integral(-800,800) << endl;
  }
  if(funcP2 && flavor=="Light"){
    leg->AddEntry(funcP2, "parametrization (1G)", "L");
    funcP2->SetNpx(1000);
    funcP2->SetParameter(0, normFunc);
    cout << "Param 1G: " << funcP2->Integral(-800,800) << endl;
    if( TMath::Abs(funcP2->Integral(-800,800)/func->Integral(-800,800)-1)>1e-02  )
      funcP2->SetParameter(0, funcP2->GetParameter(0)*func->Integral(-800,800)/funcP2->Integral(-800,800)  );
    cout << "RESCALE: Param 1G: " << funcP2->Integral(-800,800) << endl;
  }

  hBestMass->GetXaxis()->SetRangeUser(0,350);
  hBestMass->SetMaximum( hBestMass->GetMaximum()*1.4);
  hBestMass->Draw("PE");
  if( func ) func->Draw("SAME");
  if( funcP )  funcP->Draw("SAME");
  if( funcP2 && flavor=="Light") funcP2->Draw("SAME");
  leg->Draw();

  if(1){
    string eventStr = event<10 ? string(Form("0%d",event)) : string(Form("%d",event));
    c1->SaveAs(Form("Plots/PAS/Plot_%s_%s_%d%s.pdf", flavor.c_str(), eventStr.c_str(), bin, extraLabel.Data()));
  }

  delete c1;

  return;
}


void plot_SliceAll(){

  plot_Slice( "u-jet. E_{q}=50 GeV, |#eta_{q}|<1.0", "jet energy (GeV)", "TF(E|#hat{E})",   5, "Light", 0, "", "fit (1G)");
  plot_Slice( "u-jet. E_{q}=100 GeV, |#eta_{q}|<1.0", "jet energy (GeV)", "TF(E|#hat{E})", 15, "Light", 0, "", "fit (1G)");
  plot_Slice( "u-jet. E_{q}=150 GeV, |#eta_{q}|<1.0", "jet energy (GeV)", "TF(E|#hat{E})", 26, "Light", 0, "", "fit (1G)");
  plot_Slice( "u-jet. E_{q}=200 GeV, |#eta_{q}|<1.0", "jet energy (GeV)", "TF(E|#hat{E})", 37, "Light", 0, "", "fit (1G)");

  plot_Slice( "u-jet. E_{q}=50 GeV, |#eta_{q}|>1.0", "jet energy (GeV)", "TF(E|#hat{E})",   5, "Light", 1, "", "fit (1G)");
  plot_Slice( "u-jet. E_{q}=100 GeV, |#eta_{q}|>1.0", "jet energy (GeV)", "TF(E|#hat{E})", 15, "Light", 1, "", "fit (1G)");
  plot_Slice( "u-jet. E_{q}=150 GeV, |#eta_{q}|>1.0", "jet energy (GeV)", "TF(E|#hat{E})", 26, "Light", 1, "", "fit (1G)");
  plot_Slice( "u-jet. E_{q}=200 GeV, |#eta_{q}|>1.0", "jet energy (GeV)", "TF(E|#hat{E})", 37, "Light", 1, "", "fit (1G)");

  plot_Slice( "b-jet. E_{q}=50 GeV, |#eta_{q}|<1.0", "jet energy (GeV)", "TF(E|#hat{E})",   5, "Heavy", 0, "", "fit (2G)");
  plot_Slice( "b-jet. E_{q}=100 GeV, |#eta_{q}|<1.0", "jet energy (GeV)", "TF(E|#hat{E})", 15, "Heavy", 0, "", "fit (2G)");
  plot_Slice( "b-jet. E_{q}=150 GeV, |#eta_{q}|<1.0", "jet energy (GeV)", "TF(E|#hat{E})", 26, "Heavy", 0, "", "fit (2G)");
  plot_Slice( "b-jet. E_{q}=200 GeV, |#eta_{q}|<1.0", "jet energy (GeV)", "TF(E|#hat{E})", 37, "Heavy", 0, "", "fit (2G)");

  plot_Slice( "b-jet. E_{q}=50 GeV, |#eta_{q}|>1.0", "jet energy (GeV)", "TF(E|#hat{E})",   5, "Heavy", 1, "", "fit (2G)");
  plot_Slice( "b-jet. E_{q}=100 GeV, |#eta_{q}|>1.0", "jet energy (GeV)", "TF(E|#hat{E})", 15, "Heavy", 1, "", "fit (2G)");
  plot_Slice( "b-jet. E_{q}=150 GeV, |#eta_{q}|>1.0", "jet energy (GeV)", "TF(E|#hat{E})", 26, "Heavy", 1, "", "fit (2G)");
  plot_Slice( "b-jet. E_{q}=200 GeV, |#eta_{q}|>1.0", "jet energy (GeV)", "TF(E|#hat{E})", 37, "Heavy", 1, "", "fit (2G)");



}




void plot_BTag( TString extraname = "",	       
		int flag = 0,
		int plotDistr = 0,
		int oldFiles = 1){

  /*
    flag
    0   tt+jets vs TTH, SL >=6 jets
    1   tt+jj   vs TTH, SL >=6 jets
    2   tt+b    vs TTH, SL >=6 jets
    3   tt+bb   vs TTH, SL >=6 jets
    4   tt+cc   vs TTH, SL >=6 jets
    5   tt+jets vs TTH, SL >=5 jets
    6   tt+jj   vs TTH, SL >=5 jets
    7   tt+b    vs TTH, SL >=5 jets
    8   tt+bb   vs TTH, SL >=5 jets
    9   tt+cc   vs TTH, SL >=5 jets
   10   tt+jets vs TTH, DL >=4 jets
   11   tt+jj   vs TTH, DL >=4 jets
   12   tt+b    vs TTH, DL >=4 jets
   13   tt+bb   vs TTH, DL >=4 jets
   14   tt+cc   vs TTH, DL >=4 jets
  */

  int isSL6j = flag<5;
  int isSL5j = flag>=5  && flag<10;
  int isDL4j = flag>=10 && flag<15;

  int nmax = 4000;
  if( isDL4j ) nmax = 4000;

  if(plotDistr==1){
    nmax = 100;
  }

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

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  if(plotDistr==1) c1->SetGrid(0,0);
  else  c1->SetGrid(1,1);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  if(plotDistr==0) c1->SetLogy(1);

  TLegend* leg = new TLegend(0.14,0.69,0.36,0.87,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04); 
  if( isSL6j ) 
    leg->SetHeader("SL, N_{jet}#geq6");
  else if( isSL5j )   
    leg->SetHeader("SL, N_{jet}=5");
  else if( isDL4j )   
    leg->SetHeader("DL, N_{jet}#geq4");

  float kappa1 = 1;

  float min_bb = 0; 

  TString titleDistr("Simulation #sqrt{s}=8 TeV;  b_{LR}; units");

  TH1F* hS = new TH1F("hS",titleDistr, nmax, 0, 1); hS->SetLineColor(kRed) ; hS->SetLineWidth(2); hS->SetFillStyle(3004); hS->SetFillColor(kRed); 
  TH1F* hB = new TH1F("hB",titleDistr, nmax, 0, 1); hB->SetLineColor(kBlue); hB->SetLineWidth(2); hB->SetFillStyle(3005); hB->SetFillColor(kBlue); 
  
  float xM[3], yM[3];

  float xWP[1], yWP[1];
  xWP[0] = -99;
  yWP[0] = -99;

  vector<string> files;
  //string path = oldFiles ? "~/public/BTagging/Nov21_2013/" : "./";
  string path = oldFiles ? "files/test_CMVA/" : "files/test_CMVA/";
  string name =  "New" ;
  string jets =  "J" ;
  string appendix = "_nC";

  if(oldFiles==0){
  
    if( isSL6j || isSL5j ){
      files.push_back(path+"MEAnalysisNew_CMVA_rec_std_TTH125.root");
      files.push_back(path+"MEAnalysisNew_CMVA_rec_std_TTJetsSemiLept"+appendix+".root");
    }
    else{
      files.push_back(path+"MEAnalysisNew_CMVA_rec_std_TTH125.root");
      files.push_back(path+"MEAnalysisNew_CMVA_rec_std_TTJetsFullLept"+appendix+".root");
    }
  
  }
  else{

    if( isSL6j || isSL5j ){
      files.push_back(path+"MEAnalysisNew_CSV_rec_std_TTH125.root");
      files.push_back(path+"MEAnalysisNew_CSV_rec_std_TTJetsSemiLept"+appendix+".root");
    }
    else{
      files.push_back(path+"MEAnalysisNew_CSV_rec_std_TTH125.root");
      files.push_back(path+"MEAnalysisNew_CSV_rec_std_TTJetsFullLept"+appendix+".root");
    }

  }
  
  /*
  else{
    if(flag<4){
      files.push_back(path+"MEAnalysis"+name+"_TTH_6"+jets+".root");
      files.push_back(path+"MEAnalysis"+name+"_TTJETS"+"_6"+jets+".root");
    }
    else if(flag>=4 && flag<8){
      files.push_back(path+"MEAnalysis"+name+"_TTH_5"+jets+".root");
      files.push_back(path+"MEAnalysis"+name+"_TTJETS"+"_5"+jets+".root");
    }
    else if(flag>=8 && flag<12){
      files.push_back(path+"MEAnalysis"+name+"_TTH_4"+jets+".root");
      files.push_back(path+"MEAnalysis"+name+"_TTJETS"+"_4"+jets+".root");
    }
  }
  */
  

  for(int j = 0; j<files.size(); j++){

    TFile *f = TFile::Open( files[j].c_str(), "READ");
    TTree* t = (TTree*)f->Get("tree");    

    float p_tt_bb[999];
    float p_tt_jj[999];
    float p_tt_bj[999];
    int nPermut_s;
    int nPermut_b;
    int numBTagM,numBTagL,numBTagT;
    int nSimBs;
    int nMatchSimBs, nMatchSimCs;
    int type;

    t->SetBranchAddress("p_tt_bb", p_tt_bb);
    t->SetBranchAddress("p_tt_jj", p_tt_jj);
    t->SetBranchAddress("p_tt_bj", p_tt_bj);
    t->SetBranchAddress("nPermut_s", &nPermut_s);
    t->SetBranchAddress("nPermut_b", &nPermut_b);
    t->SetBranchAddress("numBTagM",&numBTagM);
    t->SetBranchAddress("numBTagT",&numBTagT);
    t->SetBranchAddress("numBTagL",&numBTagL);
    t->SetBranchAddress("nSimBs",  &nSimBs);
    t->SetBranchAddress("nMatchSimBs",  &nMatchSimBs);
    t->SetBranchAddress("type",    &type);
    if( t->GetBranch("nMatchSimCs") )
      t->SetBranchAddress("nMatchSimCs",  &nMatchSimCs);
    else
      nMatchSimCs = 0;

    int counter2M = 0;
    int counter3M = 0;
    int counter4M = 0;

    Long64_t nentries = t->GetEntries(); 
    int counter = 0;

    for (Long64_t i = 0; i < nentries ; i++){

      t->GetEntry(i);

      //if( oldFiles==0 ){
      if( isSL6j  && type!=-1) continue;
      if( isSL5j  && type!=-2) continue;
      if( isDL4j  && type!=-3) continue;
      //}
      
      
      if(j==1 && flag%5==0     && nSimBs==0) continue; // tt+X
      if(j==1 && (flag-1)%5==0 && (nSimBs >2 || (nSimBs==2 && nMatchSimCs>0)))   continue; // tt+LF
      if(j==1 && (flag-2)%5==0 && (nSimBs==2 || (nSimBs >2 && nMatchSimBs>=2) )) continue; // tt+b
      if(j==1 && (flag-3)%5==0 && (nSimBs==2 || (nSimBs >2 && nMatchSimBs<2) ))  continue; // tt+bb
      if(j==1 && (flag-4)%5==0 && (nSimBs >2 || (nSimBs==2 && nMatchSimCs<1)))   continue; // tt+cc

      counter++;

      float max_p_bb = -1;
      float max_p_jj = -1;

      float next_to_max_p_bb = -1;

      float p_bb = 0.0;
      for(int m = 0 ; m<nPermut_s; m++){
	p_bb += p_tt_bb[m];
	if( p_tt_bb[m]>max_p_bb ) max_p_bb = p_tt_bb[m];
	if( p_tt_bb[m]>next_to_max_p_bb &&  p_tt_bb[m]< max_p_bb ) next_to_max_p_bb = p_tt_bb[m];
      }
      float p_jj = 0.0;
      for(int m = 0 ; m<nPermut_b; m++){
	p_jj += p_tt_jj[m];
	if( p_tt_jj[m]>max_p_jj ) max_p_jj = p_tt_jj[m];
      }
      float p_bj = 0.0;
      for(int m = 0 ; m<nPermut_b; m++){
	p_bj += p_tt_bj[m];
      }
      
      p_bb /= nPermut_s;
      p_jj /= nPermut_b;
      
      if(numBTagL>=4 && numBTagM>=3) counter2M++;
      if(numBTagM>=4 && numBTagT>1) counter3M++;
      if(numBTagM>=4) counter4M++;
     
      if(j==0) hS->Fill( p_bb>min_bb ? p_bb/(p_bb+kappa1*p_jj) : 0. );
      if(j==1) hB->Fill( p_bb>min_bb ? p_bb/(p_bb+kappa1*p_jj) : 0. );

    }

    //cout << "File " << j << endl;
    //cout << counter2M << endl;
    //cout << counter3M << endl;
    //cout << counter4M << endl;

    if( j==0 ){
      xM[0] =  float(counter2M)/float(counter);
      xM[1] =  float(counter3M)/float(counter);
      xM[2] =  float(counter4M)/float(counter);
    }
    if( j==1 ){
      yM[0] =  float(counter2M)/float(counter);
      yM[1] =  float(counter3M)/float(counter);
      yM[2] =  float(counter4M)/float(counter);
    }
    

  }
  
  float x[nmax+1], y[nmax+1];

  for(int b = 1; b <= hS->GetNbinsX(); b++ ){
    float intS = hS->Integral(b,hS->GetNbinsX())/hS->Integral();
    float intB = hB->Integral(b,hS->GetNbinsX())/hB->Integral();
    x[b-1] = intS;
    y[b-1] = intB;				
    if((b%2)==0) cout << "At x>" << hS->GetBinCenter(b) << " => Eff.=" << intS << " -- vs -- FkR.=" << intB << endl;

    if( isSL6j &&  hS->GetBinCenter(b)>0.959 &&   hS->GetBinCenter(b)<0.961 && xWP[0]<0 ){
      xWP[0] = intS;
      yWP[0] = intB;
      cout << "\e[1;32m" << intS << " --- " << intB << "\e[0m" << endl;
    }
    if( isSL5j &&  hS->GetBinCenter(b)>0.969 &&   hS->GetBinCenter(b)<0.971 && xWP[0]<0 ){
      xWP[0] = intS;
      yWP[0] = intB;
      cout << "\e[1;32m" << intS << " --- " << intB << "\e[0m" << endl;
    }
    if( isDL4j  &&  hS->GetBinCenter(b)>0.859 &&   hS->GetBinCenter(b)<0.861 && xWP[0]<0 ){
      xWP[0] = intS;
      yWP[0] = intB;
      cout << "\e[1;32m" << intS << " --- " << intB << "\e[0m" << endl;
    }

  }
  x[nmax]=0.00;
  y[nmax]=0.00;
  cout << "Eff.=" << x[nmax] << " -- vs -- FkR.=" << y[nmax] << endl;


  ///////////////////////////////////////

  TString title("Simulation #sqrt{s}=8 TeV; selection efficiency in t#bar{t}H; selection efficiency in t#bar{t}+jets");
  if( flag%5==0    ) title = "Simulation #sqrt{s}=8 TeV; selection efficiency in t#bar{t}H; selection efficiency in t#bar{t}+jets";
  if( (flag-1)%5==0) title = "Simulation #sqrt{s}=8 TeV; selection efficiency in t#bar{t}H; selection efficiency in t#bar{t}+jj";
  if( (flag-2)%5==0) title = "Simulation #sqrt{s}=8 TeV; selection efficiency in t#bar{t}H; selection efficiency in t#bar{t}+b ";
  if( (flag-3)%5==0) title = "Simulation #sqrt{s}=8 TeV; selection efficiency in t#bar{t}H; selection efficiency in t#bar{t}+b#bar{b}";
  if( (flag-4)%5==0) title = "Simulation #sqrt{s}=8 TeV; selection efficiency in t#bar{t}H; selection efficiency in t#bar{t}+c#bar{c}";

  TH1F* h = new TH1F("h",title, 100, 0.0, 0.6);  
  if(flag==0){
    h->SetMinimum(0.001);
    h->SetMaximum(0.20);
  }
  if(flag==1){
    h->SetMinimum(0.0005);
    h->SetMaximum(0.20);
  }
  if(flag==2){
    h->SetMinimum(0.001);
    h->SetMaximum(0.30);
  }
  if(flag==3){
    h->SetMinimum(0.01);
    h->SetMaximum(0.40);
  }
  if(flag==4){
    h->SetMinimum(0.0005);
    h->SetMaximum(0.20);
  }

  if(flag==5){
    h->SetMinimum(0.001);
    h->SetMaximum(0.20);
  }
  if(flag==6){
    h->SetMinimum(0.0002);
    h->SetMaximum(0.20);
  }
  if(flag==7){
    h->SetMinimum(0.001);
    h->SetMaximum(0.30);
  }
  if(flag==8){
    h->SetMinimum(0.01);
    h->SetMaximum(0.40);
  }
  if(flag==9){
    h->SetMinimum(0.0002);
    h->SetMaximum(0.20);
  }

  if(flag==10){
    h->SetMinimum(0.001);
    h->SetMaximum(0.05);
  }
  if(flag==11){
    h->SetMinimum(0.0001);
    h->SetMaximum(0.05);
  }
  if(flag==12){
    h->SetMinimum(0.001);
    h->SetMaximum(0.20);
  }
  if(flag==13){
   h->SetMinimum(0.04);
   h->SetMaximum(0.40);
  }
  if(flag==14){
    h->SetMinimum(0.0001);
    h->SetMaximum(0.05);
  }

  h->Draw();
  TGraph* hROC = new TGraph(nmax+1, x, y);
  hROC->SetLineWidth(2);
  hROC->SetLineColor(kRed);
  hROC->GetXaxis()->SetTitleSize(0.05);
  hROC->GetYaxis()->SetTitleSize(0.05);

  hROC->Draw();

  ///////////////////////////////////////
 
  TGraph* hROCM = new TGraph(3, xM, yM);
  hROCM->SetMarkerStyle(kFullCircle);
  hROCM->SetMarkerSize(1.3);
  hROCM->SetMarkerColor(kBlue);
  hROCM->Draw("PSAME");

  TGraph* hWP = new TGraph(1, xWP, yWP);
  hWP->SetMarkerStyle(kFullCross);
  hWP->SetMarkerSize(2);
  hWP->SetMarkerColor(kRed);
  hWP->Draw("PSAME");

  cout << "2CSVM:" << xM[0] << "-" << yM[0] << "  ";
  cout << "3CSVM:" << xM[1] << "-" << yM[1] << "  ";
  cout << "4CSVM:" << xM[2] << "-" << yM[2] << endl;

  if(plotDistr==0){
   leg->AddEntry(hROCM, "4L3M, 4M, 4M2T","P");
   leg->AddEntry(hROC, "b_{LR}","L");
   leg->AddEntry(hWP,  "b_{LR} WP","P");
   leg->Draw();
  }

  if(plotDistr==1){
    leg->AddEntry(hS, "t#bar{t}H","F");
    if(flag==0 || flag==5 || flag==10)
      leg->AddEntry(hB, "t#bar{t}+jets","F");
    if(flag==1 || flag==6 || flag==11)
      leg->AddEntry(hB, "t#bar{t}+jj","F");
    if(flag==2 || flag==7 || flag==12)
      leg->AddEntry(hB, "t#bar{t}+b","F");
    if(flag==3 || flag==8 || flag==13)
      leg->AddEntry(hB, "t#bar{t}+b#bar{b}","F");
    if(flag==4 || flag==9 || flag==14)
      leg->AddEntry(hB, "t#bar{t}+c#bar{c}","F");

    hS->SetMinimum(0.);
    hB->SetMinimum(0.);

    hS->GetXaxis()->SetTitleSize(0.05);
    hS->GetYaxis()->SetTitleSize(0.05);
  
    hS->Sumw2(); 
    hB->Sumw2();
    hS->DrawNormalized("HISTE");
    hB->DrawNormalized("HISTESAME");
    leg->Draw();
  }

  if(1){
    c1->SaveAs("Plots/PAS/Plot_BTag_ROC_"+extraname+".pdf");
  }

}

void plot_BTagAll(){
  
  for(int a = 0; a<15; a++){
    for(int b = 0; b<2; b++){
      plot_BTag(TString(Form("%d_%d",a,b)),a,b);
    }
  }
  return;
}




void plot_param(TString header = "Light jet TF",	      
		TString unitsX = "jet E (GeV)",
		TString unitsY = "Acceptance",
		TString hname  = "accLightBin0",
		float xLow  = 20,
		float xHigh = 100,
		float yLow  = -999,
		float yHigh = -999,
		TString param = "",
		TString extraLabel = ""
		){

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


  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TFile *f = TFile::Open("../bin/root/ControlPlotsTEST.root","READ");
  TH1F* hBestMass = (TH1F*)f->Get(hname);


  if( !hBestMass ){
    cout << "No such histogram" << endl;
    return;
  }


  TString fname;
  TIter iter( hBestMass->GetListOfFunctions() );
  TObject* obj = iter();
  while ( obj!=0 ){
    TF1* func = dynamic_cast<TF1*>(obj);
    fname = func->GetName();
    obj = iter();
  }
  TF1 *func = hBestMass->GetFunction(fname);
  if(func) func->SetRange(xLow,xHigh);

  hBestMass->SetMarkerStyle(kFullCircle);
  hBestMass->SetMarkerColor(kBlue);
  hBestMass->SetTitle("CMS Simulation #sqrt{s}=8 TeV");
  hBestMass->SetXTitle( unitsX );
  hBestMass->SetYTitle( unitsY );
  hBestMass->SetTitleSize(0.05,"X");
  hBestMass->SetTitleSize(0.05,"Y");

  hBestMass->SetTitleOffset(0.90,"X");
  hBestMass->SetTitleOffset(0.83,"Y");

  hBestMass->GetXaxis()->SetRangeUser(xLow,xHigh);
  if(func) func->GetXaxis()->SetRangeUser(xLow,xHigh);
  if(yLow!=-999 && yHigh!=-999){
    hBestMass->GetYaxis()->SetRangeUser(yLow,yHigh);
    if(func) func->GetYaxis()->SetRangeUser(yLow,yHigh);
  }

  TLegend* leg = new TLegend(0.15,0.65,0.55,0.85,NULL,"brNDC");
  leg->SetHeader( header );
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.06); 
  leg->AddEntry(hBestMass, "from fit", "P");
  if(func) leg->AddEntry(func, param, "L");


  if(func){
    hBestMass->Draw("PE");
    func->Draw("SAME");
  }
  else{
    hBestMass->SetFillColor(kRed);
    hBestMass->SetFillStyle(3004);
    hBestMass->Draw("HIST");
  }

  leg->Draw();

  if(func) cout << "f(50,100,150)=" << string(Form("(%.1f, %.1f, %.1f)",func->Eval(50), func->Eval(100), func->Eval(150))) << endl;

  if(1){
    c1->SaveAs("Plots/PAS/Plot_"+hname+"_tf_"+extraLabel+".pdf");
  }

  //delete c1;
  return;
}


void plot_paramAll(){

  vector<string> hnames;
  hnames.push_back("resolLightBin0");
  hnames.push_back("respLightBin0");
  hnames.push_back("resolLightBin1");
  hnames.push_back("respLightBin1");
  hnames.push_back("respG1HeavyBin0");
  hnames.push_back("respG2HeavyBin0");
  hnames.push_back("resolG1HeavyBin0");
  hnames.push_back("resolG2HeavyBin0");
  hnames.push_back("respG1HeavyBin1");
  hnames.push_back("respG2HeavyBin1");
  hnames.push_back("resolG1HeavyBin1");
  hnames.push_back("resolG2HeavyBin1");


  vector<string> titles;
  titles.push_back("u-jet, #sigma(E), |#eta|<1.0");
  titles.push_back("u-jet, #mu(E),  |#eta|<1.0");
  titles.push_back("u-jet, #sigma(E), |#eta|>1.0");
  titles.push_back("u-jet, #mu(E),  |#eta|>1.0");
  titles.push_back("b-jet, #mu(E), |#eta|<1.0");
  titles.push_back("b-jet, #mu'(E), |#eta|<1.0");
  titles.push_back("b-jet, #sigma(E), |#eta|<1.0");
  titles.push_back("b-jet, #sigma'(E), |#eta|<1.0");
  titles.push_back("b-jet, #mu(E), |#eta|>1.0");
  titles.push_back("b-jet, #mu'(E), |#eta|>1.0");
  titles.push_back("b-jet, #sigma(E), |#eta|>1.0");
  titles.push_back("b-jet, #sigma'(E), |#eta|>1.0");



  for(int k = 0 ; k < hnames.size(); k++){
    continue;
    cout << hnames[k] << endl;
    if( hnames[k].find("resp")!=string::npos )
      plot_param( TString(titles[k].c_str()) ,"quark energy (GeV)","GeV",TString(hnames[k]),  50 , 200, 50, 200, "parametrization", ""); 
    if( hnames[k].find("resol")!=string::npos )
      plot_param( TString(titles[k].c_str()) ,"quark energy (GeV)","GeV",TString(hnames[k]),  50 , 200, 0, 100, "parametrization", ""); 

    cout << endl;
  }

  plot_param("E_{x,y}^{miss}, #sigma_{E_{T}}", "#sum E_{T} (GeV)", "GeV", "hWidthsPxResol", 500, 3000, 0, 70, "parametrization", "");

}


void plot_category(string dir  = "Mar25_2014",
		   string type = "MEM", 
		   string cat  = "cat1",
		   string tag = "New_rec_std_byLLR_bj",
		   string header = "Cat 1",
		   string fname  = "",
		   TString title = "S/(S+B)",
		   int    drawseparation = 1,
		   float xMax = -1,
		   float  rescaleTTbb = 1.,
		   float  rescaleTTb  = 1.,
		   float  rescaleTTjj = 1.
		   ){
  
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

  TFile* f = TFile::Open(("datacards/"+dir+"/"+type+"_"+tag+".root").c_str());
  if(f==0 || f->IsZombie() ) return;


  TH1F* hS    = 0;
  TH1F* hData = 0;
  TH1F* hErr  = 0;
  TH1F* h1    = 0;


  TCanvas *c1 = new TCanvas("c1","",5,30,650,600/0.85);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TPad *p1 = new TPad("p1","",0,0.155,1,1);
  p1->SetGrid(0,0);
  p1->SetFillStyle(4000);
  p1->SetFillColor(10);
  p1->SetTicky();
  p1->SetObjectStat(0);
  p1->Draw();
  p1->cd();

  TLegend* leg = new TLegend(0.62,0.45,0.88,0.89,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04); 

  THStack* aStack = new THStack("aStack","CMS preliminary #sqrt{s}=8 TeV, L=19.5 fb^{-1}; P_{s/b} ; events ");

  /*
  vector<string> samples;
  samples.push_back("SingleT");
  samples.push_back("TTV");
  samples.push_back("EWK");
  samples.push_back("TTJetsHFbb");
  samples.push_back("TTJetsHFb");
  if( SPLITCC )
    samples.push_back("TTJetsLFcc");
  samples.push_back("TTJetsLF");
  samples.push_back("TTH125");
  samples.push_back("DiBoson");
  if(RUNONDATA)  samples.push_back("data_obs");
  */

  vector<string> samples;
  samples.push_back("singlet");
  samples.push_back("ttbarV");
  samples.push_back("ewk");
  samples.push_back("ttbarPlusBBbar");
  samples.push_back("ttbarPlusB");
  if( SPLITCC )
    samples.push_back("ttbarPlusCCbar");
  samples.push_back("ttbar");
  samples.push_back("ttH_hbb");
  samples.push_back("diboson");
  if(RUNONDATA)  samples.push_back("data_obs");

  ////////////////////////////////////////////////////////////////////////

  map<string,TH1F*> h_Ups;

  f->cd( (type+"_"+cat).c_str());
 
  //h1 = (TH1F*)((TH1F*)f->Get( (type+"_"+cat+"/TTH125").c_str() ))->Clone("h1");
  h1 = (TH1F*)((TH1F*)f->Get( (type+"_"+cat+"/ttH_hbb").c_str() ))->Clone("h1");
  h1->Reset();

  TIter iter( gDirectory->GetListOfKeys()  );
  TKey *key;
  
  while ( (key = (TKey*)iter.Next()) ){
    
    TH1F* hist =  (TH1F*)key->ReadObj();    
    string hname = string(hist->GetName());

    if( hname.find("TTJetsHFbb")!=string::npos){
      hist->Scale( rescaleTTbb);
      //cout << hname << " scaled by " << rescaleTTbb << endl;
    }
    else if(hname.find("TTJetsHFb")!=string::npos){
      hist->Scale( rescaleTTb);
      //cout << hname << " scaled by " << rescaleTTb << endl;
    }
    else if(hname=="TTJetsLF" || hname.find("TTJetsLF_")!=string::npos ){
      hist->Scale( rescaleTTjj );
      //cout << hname << " scaled by " << rescaleTTjj << endl;
    }
    
    if( //hname.find("TTH125")==string::npos && 
       hname.find("ttH_hbb")==string::npos && 
       hname.find("data")==string::npos && 
       hname.find("bin")==string::npos && 
       hname.find("Down")==string::npos &&
       hname.find("Up")!=string::npos ){

      string sysname = hname;
      sysname.erase(0, (sysname.find_first_of("_"))+1 );
      cout << sysname << endl;

      // WORKAROUND
      // if( sysname.find("pdf_qq")!=string::npos ) continue;

      TH1F* h = new TH1F( *hist );

      if( h_Ups.find(sysname) == h_Ups.end() ){
	h_Ups[ sysname ] = h;
      }
      else{
	h_Ups[ sysname ]->Add( h );
      }

    }

  }

  cout << endl << "SUMMARY: " << endl;
  for( map<string,TH1F*>::iterator it = h_Ups.begin() ; it != h_Ups.end() ; it++){
    if( (it->second) )   cout << it->first << " ==> " << (it->second)->Integral() << endl ;
  }
  

  
  ////////////////////////////////////////////////////////////////////////



  for(unsigned int sample = 0; sample < samples.size(); sample++){

    cout << "Doing " << samples[sample] << endl;

  
    TH1F* h = (TH1F*)f->Get((type+"_"+cat+"/"+samples[sample]).c_str());
    if( h==0 ){
      cout << (type+"_"+cat+"/"+samples[sample]) << " not found" << endl;
      continue;
    }
    cout << "Events = $" << h->Integral() << " \\pm " << (h->GetEntries()>0 ? sqrt(h->GetEntries())*h->Integral()/h->GetEntries() : 0.0 ) << " $" << endl;


    if( //samples[sample].find("TTH125") != string::npos ){
       samples[sample].find("ttH_hbb") != string::npos ){
      hS = (TH1F*)h->Clone("hS");
      hS->Reset();
      hS->Add(h,10.0);
      hS->SetLineWidth(3);
      hS->SetLineStyle(kDashed);
      hS->SetLineColor(kRed);

      h->SetLineColor( kRed );
      h->SetFillColor( kRed );
      h->SetFillStyle( 3002 );
      leg->AddEntry(h, "t#bar{t}H", "F");
    }
    if( //samples[sample].find("TTJetsHFbb") != string::npos ){
       samples[sample].find("ttbarPlusBBbar") != string::npos ){
      leg->AddEntry(h, "t#bar{t} + bb", "F");
      //h->Scale( rescaleTTbb );
      h->SetLineColor( 15 );
      h->SetFillColor( 15 );
    }
    else if( //samples[sample].find("TTJetsHFb") != string::npos ){
	    samples[sample].find("ttbarPlusB") != string::npos && samples[sample].find("BBbar") == string::npos ){
      leg->AddEntry(h, "t#bar{t} + b", "F");
      //h->Scale( rescaleTTb );
      h->SetLineColor( 16 );
      h->SetFillColor( 16 );
    }
    else if( //SPLITCC && samples[sample].find("TTJetsLFcc") != string::npos ){
	    SPLITCC && samples[sample].find("ttbarPlusCCbar") != string::npos ){
      leg->AddEntry(h, "t#bar{t} + cc", "F");
      //h->Scale( rescaleTTjj );                                                                                                        
      h->SetLineColor( 17 );
      h->SetFillColor( 17 );
    }
    else if( //samples[sample].find("TTJetsLF") != string::npos ){
	    samples[sample] == "ttbar" ){
      leg->AddEntry(h, "t#bar{t} + jj", "F");
      //h->Scale( rescaleTTjj );
      h->SetLineColor( 17+SPLITCC );
      h->SetFillColor( 17+SPLITCC );
    }    
    if( //samples[sample].find("SingleT") != string::npos ){
       samples[sample].find("singlet") != string::npos ){
      h->SetLineColor( kMagenta );
      h->SetFillColor( kMagenta );
      leg->AddEntry(h, "Single top", "F");
    }
    if( //samples[sample].find("EWK") != string::npos ){
       samples[sample].find("ewk") != string::npos ){
      h->SetLineColor( kGreen );
      h->SetFillColor( kGreen );
      leg->AddEntry(h, "V+jets", "F");
    }
    if( //samples[sample].find("DiBoson") != string::npos ){
       samples[sample].find("diboson") != string::npos ){
      h->SetLineColor( kYellow );
      h->SetFillColor( kYellow );
      leg->AddEntry(h, "VV", "F");
    }
    if( //samples[sample].find("TTV") != string::npos ){
       samples[sample].find("ttbarV") != string::npos ){
      h->SetLineColor( 30 );
      h->SetFillColor( 30 );
      leg->AddEntry(h, "t#bar{t}V", "F");
    }
    if( samples[sample].find("data_obs") != string::npos && RUNONDATA){
      hData = (TH1F*)h->Clone("hData");
      hData->Sumw2();
      hData->SetMarkerStyle(kFullCircle);
      hData->SetMarkerSize(1.5);
      hData->SetMarkerColor(kBlack);
      leg->AddEntry(hData, "data", "P");
    }

    if( !( /*samples[sample].find("TTH125")!=string::npos*/ samples[sample].find("ttH_hbb")!=string::npos || samples[sample].find("data_obs")!=string::npos) ){
      if( hErr==0 ){
	hErr = (TH1F*)h->Clone("hErr");
	hErr->Reset();
	hErr->SetFillStyle(0);
	hErr->SetLineColor(kBlack);
      }
      hErr->Add( h, 1.0);
    }

    if( samples[sample].find("data_obs") == string::npos )
      aStack->Add( h );
  }

  if(hErr==0 || hS==0) return;

  for(int b = 1; b <= hErr->GetNbinsX() ; b++){

    float syst_b2 = 0.;
    for( map<string,TH1F*>::iterator it = h_Ups.begin() ; it != h_Ups.end() ; it++){
      if( (it->second) ) syst_b2 += TMath::Power(hErr->GetBinContent(b)- (it->second)->GetBinContent(b), 2);
    }   
    hErr->SetBinError(b, sqrt(hErr->GetBinError(b)*hErr->GetBinError(b) + syst_b2) );
  }


  hErr->GetYaxis()->SetTitle("Events");
  if( tag.find("_sb") != string::npos )
    hErr->GetXaxis()->SetTitle("P_{s/b}");
  if( tag.find("_sb_nb") != string::npos )
    hErr->GetXaxis()->SetTitle("P_{s/b}^{1 perm}");
  if( tag.find("_bj") != string::npos )
    hErr->GetXaxis()->SetTitle("P_{b/j}");

  if(RUNONDATA==0)
    hErr->SetTitle("Simulation #sqrt{s}=8 TeV, L=19.5 fb^{-1}");
  else
    hErr->SetTitle("CMS Preliminary #sqrt{s}=8 TeV, L=19.5 fb^{-1}");
  hErr->SetTitleSize  (0.04,"X");
  hErr->SetTitleOffset(0.95,"X");
  float max =  TMath::Max(float(hErr->GetMaximum()*1.45),(hData!=0 ? hData->GetMaximum()*1.45 : -1.));
  hErr->GetYaxis()->SetRangeUser(0., max );
  hErr->GetXaxis()->SetRangeUser( hErr->GetXaxis()->GetXmin(), xMax>0 ? xMax : hErr->GetXaxis()->GetXmax());
  hErr->SetLineColor(kBlack);


  hErr->GetXaxis()->SetTitle( title );

  hErr->Draw("HIST");
  TH1D* hErr_clone = (TH1D*)hErr->Clone("hErr_clone");
  hErr_clone  ->SetFillColor(kBlack);
  hErr_clone  ->SetFillStyle(3004);
  hErr_clone  ->SetLineWidth(1);
  hErr_clone  ->Draw("E2SAME");
  leg->AddEntry(hErr_clone, "Bkg unc.", "F");

  aStack->Draw("HISTSAME");
  leg->SetHeader( header.c_str() );
  leg->AddEntry(hS, "signal x 10", "L");
  hS->Draw("HISTSAME");
  hErr_clone  ->Draw("E2SAME");

  double KS   = -99;
  double CHI2 = -99;
  if(hData!=0){
    hData->Draw("ESAME");

    TH1D* hErr_clone2 = (TH1D*)hErr->Clone("hErr_clone2");

    if( KSWOSYSTEMATICS ){
      for(int b = 1; b <= hErr_clone2->GetNbinsX() ; b++)
	hErr_clone2->SetBinError( b,  0. );
    }
    
    KS   = hErr_clone2->KolmogorovTest( hData );
    CHI2 = hErr_clone2->Chi2Test( hData );
    cout << "KS = " << KS << ", Chi2 = " << CHI2 << endl;
    //leg->AddEntry((TObject*)0, Form("KS %.3f, #chi^{2} %.3f", KS, CHI2), "");
  }
  leg->Draw();

  TPaveText *pt0 = new TPaveText(0.101, 0.839161, 0.4, 0.895105,"brNDC");
  pt0->SetFillStyle(1001);
  pt0->SetBorderSize(0);
  pt0->SetFillColor(kWhite);
  pt0->SetTextSize(0.04); 
  pt0->AddText(  Form("KS %.3f, #chi^{2} %.3f", KS, CHI2) );
  pt0->Draw();

  TLine* line = new TLine(hErr->GetBinLowEdge(3), max , hErr->GetBinLowEdge(3), 0.);
  line->SetLineWidth(4);
  line->SetLineStyle(kSolid);
  line->SetLineColor(kBlack);
  if(drawseparation) line->Draw("SAME");

  TPaveText *pt1 = new TPaveText(0.101, 0.839161, 0.198142, 0.895105,"brNDC");
  pt1->SetFillStyle(1001);
  pt1->SetBorderSize(0);
  pt1->SetFillColor(kWhite);
  pt1->SetTextSize(0.03); 
  pt1->AddText("P_{b/j}<0.5");

  TPaveText *pt2 = new TPaveText(0.191, 0.839161, 0.294118, 0.895105,"brNDC");
  pt2->SetFillStyle(1001);
  pt2->SetBorderSize(0);
  pt2->SetFillColor(kWhite);
  pt2->SetTextSize(0.03); 
  pt2->AddText("P_{b/j}>0.5");
 
  if(drawseparation) pt1->Draw();
  if(drawseparation) pt2->Draw();

  cout << "Signal = " << hS->Integral()/5. << endl;


  c1->cd();
  TPad *p2 = new TPad("p1","",0,0,1,0.175);
  p2->SetGrid(0,0);
  p2->SetFillStyle(4000);
  p2->SetFillColor(10);
  p2->SetTicky();
  p2->SetObjectStat(0);
  p2->Draw();
  p2->cd();
  
  THStack* rStack = new THStack("rStack","");
  
  TF1* one = new TF1("one","1", h1->GetXaxis()->GetXmin(),h1->GetXaxis()->GetXmax());
  
  TH1F* one_low = 0;
  one_low = (TH1F*)h1->Clone("one_low");
  
  TH1F* one_up = 0;
  one_up = (TH1F*)h1->Clone("one_up");
  
  TH1F* one_error = 0;
  one_error = (TH1F*)h1->Clone("one_error");
  
  TH1F* h_ratio = 0;
  h_ratio = (TH1F*)h1->Clone("h_ratio");
  
  for(int b=1; b<= h_ratio->GetNbinsX() ; b++){

    float ndata    = 0.;
    float ndataErr = 0.;

    if(RUNONDATA){
      ndata = hData->GetBinContent(b);
      ndataErr = hData->GetBinError(b);
    }      
    else{
      ndata = hErr->GetBinContent(b);  //should be hData
      ndataErr = hErr->GetBinError(b); //should be hData
    }
    
    float nbkg = hErr->GetBinContent(b);
    float nbkgErr = 0;
    if(hErr!=0){
      nbkgErr = hErr->GetBinError(b);
    }
    
    float nbkg_up  = hErr->GetBinContent(b) + hErr->GetBinError(b);
    float nbkg_low = hErr->GetBinContent(b) - hErr->GetBinError(b);
    
    cout << ndata << ", " << nbkg << endl;
  
    float r = nbkg >0 ? ndata/nbkg : 0; 

    float rErr = nbkg>0 >0 ? ndataErr/nbkg  : 0;

    float r_up  = nbkg_low >0 && ndata>0 ? nbkg/nbkg_low : 0;
    float r_low = nbkg_up  >0 && ndata>0 ? nbkg/nbkg_up  : 0;
    
    h_ratio->SetBinContent(b, r);
    h_ratio->SetBinError(b,rErr);
    
    one_low->SetBinContent(b, r_low);
    one_up->SetBinContent(b, r_up - r_low);
    
    one_error->SetBinContent(b, 1.0);
    one_error->SetBinError(b, (r_up-r_low)/2 );
  }
   
  h_ratio->SetMarkerStyle(kFullCircle);
  h_ratio->SetMarkerSize(1.5);
  h_ratio->SetMarkerColor(kBlack);
  h_ratio->SetTitle(0);
  h_ratio->GetYaxis()->SetRangeUser(0.0,2.0);
  h_ratio->GetXaxis()->SetTitle(0);
  h_ratio->GetXaxis()->SetLabelSize(0);
  h_ratio->GetXaxis()->SetNdivisions(0);
  h_ratio->GetXaxis()->SetTickLength(0.1);
  
  h_ratio->GetYaxis()->SetTitleSize(0.16);
  h_ratio->GetYaxis()->SetTitleOffset(0.28);
  h_ratio->GetYaxis()->SetTitle("Data/Bkg  ");
  h_ratio->GetYaxis()->SetNdivisions(202);
  h_ratio->GetYaxis()->SetLabelSize(0.16);
  
  h_ratio->Draw("EP");
  
  one->SetLineColor(kRed);
  one->SetLineWidth(2);
  
  one_error->SetFillStyle(3004);
  one_error->SetFillColor(kBlack);

  one_low->SetFillStyle(0);
  one_low->SetFillColor(18);
  one_low->SetLineColor(18);
  one_up->SetFillColor(18);
  one_up->SetLineColor(18);
  
  rStack->Add( one_low );
  rStack->Add( one_up );
  
  one_error->Draw("E2SAME");
  one->Draw("SAME");
  h_ratio->Draw("EPSAME");

  h_ratio->Print("all");
 
  cout << "REMOVE: " << endl;
  for( map<string,TH1F*>::iterator it = h_Ups.begin() ; it != h_Ups.end() ; it++){
    if( (it->second) )  delete (it->second);
  }

  if(0){
    c1->SaveAs(  ("datacards/"+dir+"/Plot_"+cat+"_"+fname+"_"+tag+".pdf").c_str() );
    f->Close();
  }

 
}


void plot_categoryAll(){

  string directory = "PreApproval/controlPlots";

  vector<string> hyps;
  hyps.push_back("0");
  hyps.push_back("1");

  for( unsigned int i = 0 ; i < hyps.size() ; i++){

    string hyp = hyps[i];

    plot_category(directory, "MEM", "best_"+hyp+"_mass_WHad",   "New_rec_std_cat1", "SL-Cat 1", "", "best W_{had} mass (GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_mass_TopHad", "New_rec_std_cat1", "SL-Cat 1", "", "best t_{had} mass (GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_mass_H",      "New_rec_std_cat1", "SL-Cat 1", "", "best Higgs mass (GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_mass_H",      "New_rec_std_cat2", "SL-Cat 2", "", "best Higgs mass (GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_mass_H",      "New_rec_std_cat3", "SL-Cat 3", "", "best Higgs mass (GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_mass_H",      "New_rec_std_cat6", "DL",       "", "best Higgs mass (GeV)", 0);
    
    plot_category(directory, "MEM", "best_"+hyp+"_mass_TopLep1","New_rec_std_cat1", "SL-Cat 1", "", "best t_{lep} mass (GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_mass_TopLep1","New_rec_std_cat2", "SL-Cat 2", "", "best t_{lep} mass (GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_mass_TopLep1","New_rec_std_cat3", "SL-Cat 3", "", "best t_{lep} mass (GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_mass_TopLep1","New_rec_std_cat6", "DL",       "", "best t_{lep} mass (GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_mass_TopLep2","New_rec_std_cat6", "DL",       "", "best t_{lep} mass (GeV)", 0);
    
    plot_category(directory, "MEM", "best_"+hyp+"_e_MET","New_rec_std_cat1", "SL-Cat 1", "", "E_{T}^{miss}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_e_MET","New_rec_std_cat2", "SL-Cat 2", "", "E_{T}^{miss}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_e_MET","New_rec_std_cat3", "SL-Cat 3", "", "E_{T}^{miss}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_e_MET","New_rec_std_cat6", "DL",       "", "E_{T}^{miss}(GeV)", 0);
    
    plot_category(directory, "MEM", "best_"+hyp+"_dphi_MET_WLep1","New_rec_std_cat1", "SL-Cat 1", "", "#Delta#phi l_{1}/E_{T}^{miss}", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_dphi_MET_WLep1","New_rec_std_cat2", "SL-Cat 2", "", "#Delta#phi l_{1}/E_{T}^{miss}", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_dphi_MET_WLep1","New_rec_std_cat3", "SL-Cat 3", "", "#Delta#phi l_{1}/E_{T}^{miss}", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_dphi_MET_WLep1","New_rec_std_cat6", "DL",       "", "#Delta#phi l_{1}/E_{T}^{miss}", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_dphi_MET_WLep2","New_rec_std_cat6", "DL",       "", "#Delta#phi l_{2}/E_{T}^{miss}", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_dphi_b1_b2","New_rec_std_cat1", "SL-Cat 1", "", "#Delta#phi b_{1}/b_{2}", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_dphi_b1_b2","New_rec_std_cat2", "SL-Cat 2", "", "#Delta#phi b_{1}/b_{2}", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_dphi_b1_b2","New_rec_std_cat3", "SL-Cat 3", "", "#Delta#phi b_{1}/b_{2}", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_dphi_b1_b2","New_rec_std_cat6", "DL",       "", "#Delta#phi b_{1}/b_{2}", 0);
        
    plot_category(directory, "MEM", "best_"+hyp+"_pt_WLep1","New_rec_std_cat1", "SL-Cat 1", "", "p_{T} l_{1}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_pt_WLep1","New_rec_std_cat2", "SL-Cat 2", "", "p_{T} l_{1}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_pt_WLep1","New_rec_std_cat3", "SL-Cat 3", "", "p_{T} l_{1}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_pt_WLep1","New_rec_std_cat6", "DL",       "", "p_{T} l_{1}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_pt_WLep2","New_rec_std_cat6", "DL",       "", "p_{T} l_{2}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_pt_bLep","New_rec_std_cat1", "SL-Cat 1", "", "p_{T} b_{l}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_pt_bLep","New_rec_std_cat2", "SL-Cat 2", "", "p_{T} b_{l}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_pt_bLep","New_rec_std_cat3", "SL-Cat 3", "", "p_{T} b_{l}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_pt_bLep","New_rec_std_cat6", "DL",       "", "p_{T} b_{l}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_pt_WHad1","New_rec_std_cat1", "SL-Cat 1", "", "p_{T} w_{1}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_pt_WHad1","New_rec_std_cat2", "SL-Cat 2", "", "p_{T} w_{1}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_pt_WHad1","New_rec_std_cat3", "SL-Cat 3", "", "p_{T} w_{1}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_pt_WHad1","New_rec_std_cat6", "DL",       "", "p_{T} w_{1}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_pt_WHad2","New_rec_std_cat1", "SL-Cat 1", "", "p_{T} w_{2}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_pt_WHad2","New_rec_std_cat2", "SL-Cat 2", "", "p_{T} w_{2}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_pt_WHad2","New_rec_std_cat3", "SL-Cat 3", "", "p_{T} w_{2}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_pt_WHad2","New_rec_std_cat6", "DL",       "", "p_{T} w_{2}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_pt_bHad","New_rec_std_cat1", "SL-Cat 1", "", "p_{T} b_{h}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_pt_bHad","New_rec_std_cat2", "SL-Cat 2", "", "p_{T} b_{h}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_pt_bHad","New_rec_std_cat3", "SL-Cat 3", "", "p_{T} b_{h}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_pt_bHad","New_rec_std_cat6", "DL",       "", "p_{T} b_{h}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_pt_b1","New_rec_std_cat1", "SL-Cat 1", "", "p_{T} b_{1}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_pt_b1","New_rec_std_cat2", "SL-Cat 2", "", "p_{T} b_{1}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_pt_b1","New_rec_std_cat3", "SL-Cat 3", "", "p_{T} b_{1}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_pt_b1","New_rec_std_cat6", "DL",       "", "p_{T} b_{1}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_pt_b2","New_rec_std_cat1", "SL-Cat 1", "", "p_{T} b_{2}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_pt_b2","New_rec_std_cat2", "SL-Cat 2", "", "p_{T} b_{2}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_pt_b2","New_rec_std_cat3", "SL-Cat 3", "", "p_{T} b_{2}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_pt_b2","New_rec_std_cat6", "DL",       "", "p_{T} b_{2}(GeV)", 0);
    
    plot_category(directory, "MEM", "best_"+hyp+"_eta_WLep1","New_rec_std_cat1", "SL-Cat 1", "", "#eta l_{1}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_eta_WLep1","New_rec_std_cat2", "SL-Cat 2", "", "#eta l_{1}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_eta_WLep1","New_rec_std_cat3", "SL-Cat 3", "", "#eta l_{1}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_eta_WLep1","New_rec_std_cat6", "DL",       "", "#eta l_{1}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_eta_WLep2","New_rec_std_cat6", "DL",       "", "#eta l_{2}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_eta_bLep","New_rec_std_cat1", "SL-Cat 1", "", "#eta b_{l}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_eta_bLep","New_rec_std_cat2", "SL-Cat 2", "", "#eta b_{l}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_eta_bLep","New_rec_std_cat3", "SL-Cat 3", "", "#eta b_{l}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_eta_bLep","New_rec_std_cat6", "DL",       "", "#eta b_{l}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_eta_WHad1","New_rec_std_cat1", "SL-Cat 1", "", "#eta w_{1}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_eta_WHad1","New_rec_std_cat2", "SL-Cat 2", "", "#eta w_{1}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_eta_WHad1","New_rec_std_cat3", "SL-Cat 3", "", "#eta w_{1}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_eta_WHad1","New_rec_std_cat6", "DL",       "", "#eta w_{1}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_eta_WHad2","New_rec_std_cat1", "SL-Cat 1", "", "#eta w_{2}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_eta_WHad2","New_rec_std_cat2", "SL-Cat 2", "", "#eta w_{2}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_eta_WHad2","New_rec_std_cat3", "SL-Cat 3", "", "#eta w_{2}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_eta_WHad2","New_rec_std_cat6", "DL",       "", "#eta w_{2}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_eta_bHad","New_rec_std_cat1", "SL-Cat 1", "", "#eta b_{h}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_eta_bHad","New_rec_std_cat2", "SL-Cat 2", "", "#eta b_{h}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_eta_bHad","New_rec_std_cat3", "SL-Cat 3", "", "#eta b_{h}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_eta_bHad","New_rec_std_cat6", "DL",       "", "#eta b_{h}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_eta_b1","New_rec_std_cat1", "SL-Cat 1", "", "#eta b_{1}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_eta_b1","New_rec_std_cat2", "SL-Cat 2", "", "#eta b_{1}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_eta_b1","New_rec_std_cat3", "SL-Cat 3", "", "#eta b_{1}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_eta_b1","New_rec_std_cat6", "DL",       "", "#eta b_{1}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_eta_b2","New_rec_std_cat1", "SL-Cat 1", "", "#eta b_{2}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_eta_b2","New_rec_std_cat2", "SL-Cat 2", "", "#eta b_{2}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_eta_b2","New_rec_std_cat3", "SL-Cat 3", "", "#eta b_{2}(GeV)", 0);
    plot_category(directory, "MEM", "best_"+hyp+"_eta_b2","New_rec_std_cat6", "DL",       "", "#eta b_{2}(GeV)", 0);

  }

  //return;

  /*
  plot_category("Apr23_2014/", "MEM", "cat1_H", "New_MH90_rec_std_sb", "SL-Cat 1 (H)", "", "P_{s/b}^{90}", 1, -1);
  plot_category("Apr23_2014/", "MEM", "cat1_L", "New_MH90_rec_std_sb", "SL-Cat 1 (L)", "", "P_{s/b}^{90}", 0, -1);

  plot_category("Apr23_2014/", "MEM", "cat2_H", "New_MH90_rec_std_sb", "SL-Cat 2 (H)", "", "P_{s/b}^{90}", 1, -1);
  plot_category("Apr23_2014/", "MEM", "cat2_L", "New_MH90_rec_std_sb", "SL-Cat 2 (L)", "", "P_{s/b}^{90}", 0, -1);

  plot_category("Apr23_2014/", "MEM", "cat3_H", "New_MH90_rec_std_sb", "SL-Cat 3 (H)", "", "P_{s/b}^{90}", 1, -1);
  plot_category("Apr23_2014/", "MEM", "cat3_L", "New_MH90_rec_std_sb", "SL-Cat 3 (L)", "", "P_{s/b}^{90}", 0, -1);

  plot_category("Apr23_2014/", "MEM", "cat6_H", "New_MH90_rec_std_sb", "DL (H)", "", "P_{s/b}^{90}", 1, -1);
  plot_category("Apr23_2014/", "MEM", "cat6_L", "New_MH90_rec_std_sb", "DL (L)", "", "P_{s/b}^{90}", 0, -1);

  return;
  */

  plot_category(directory, "MEM", "logPs", "New_rec_std_cat1_L", "SL-Cat 1 (L)", "", "log(w_{0})", 0, 50.);
  plot_category(directory, "MEM", "logPs", "New_rec_std_cat1_H", "SL-Cat 1 (H)", "", "log(w_{0})", 0, 50.);
  plot_category(directory, "MEM", "logPs", "New_rec_std_cat2_L", "SL-Cat 2 (L)", "", "log(w_{0})", 0, 50.);
  plot_category(directory, "MEM", "logPs", "New_rec_std_cat2_H", "SL-Cat 2 (H)", "", "log(w_{0})", 0, 50.);
  plot_category(directory, "MEM", "logPs", "New_rec_std_cat3_L", "SL-Cat 3 (L)", "", "log(w_{0})", 0, 50.);
  plot_category(directory, "MEM", "logPs", "New_rec_std_cat3_H", "SL-Cat 3 (H)", "", "log(w_{0})", 0, 50.);
  plot_category(directory, "MEM", "logPs", "New_rec_std_cat6_L", "DL (L)", "", "log(w_{0})", 0, 50.);
  plot_category(directory, "MEM", "logPs", "New_rec_std_cat6_H", "DL (H)", "", "log(w_{0})", 0, 50.);

  plot_category(directory, "MEM", "logPb", "New_rec_std_cat1_L", "SL-Cat 1 (L)", "", "log(w_{1})", 0, 50.);
  plot_category(directory, "MEM", "logPb", "New_rec_std_cat1_H", "SL-Cat 1 (H)", "", "log(w_{1})", 0, 50.);
  plot_category(directory, "MEM", "logPb", "New_rec_std_cat2_L", "SL-Cat 2 (L)", "", "log(w_{1})", 0, 50.);
  plot_category(directory, "MEM", "logPb", "New_rec_std_cat2_H", "SL-Cat 2 (H)", "", "log(w_{1})", 0, 50.);
  plot_category(directory, "MEM", "logPb", "New_rec_std_cat3_L", "SL-Cat 3 (L)", "", "log(w_{1})", 0, 50.);
  plot_category(directory, "MEM", "logPb", "New_rec_std_cat3_H", "SL-Cat 3 (H)", "", "log(w_{1})", 0, 50.);
  plot_category(directory, "MEM", "logPb", "New_rec_std_cat6_L", "DL (L)", "", "log(w_{1})", 0, 50.);
  plot_category(directory, "MEM", "logPb", "New_rec_std_cat6_H", "DL (H)", "", "log(w_{1})", 0, 50.);

  plot_category(directory, "MEM", "logPbb", "New_rec_std_cat1_L", "SL-Cat 1 (L)", "", "log(L_{S}^{b-tag})", 0);
  plot_category(directory, "MEM", "logPbb", "New_rec_std_cat1_H", "SL-Cat 1 (H)", "", "log(L_{S}^{b-tag})", 0);
  plot_category(directory, "MEM", "logPbb", "New_rec_std_cat2_L", "SL-Cat 2 (L)", "", "log(L_{S}^{b-tag})", 0);
  plot_category(directory, "MEM", "logPbb", "New_rec_std_cat2_H", "SL-Cat 2 (H)", "", "log(L_{S}^{b-tag})", 0);
  plot_category(directory, "MEM", "logPbb", "New_rec_std_cat3_L", "SL-Cat 3 (L)", "", "log(L_{S}^{b-tag})", 0);
  plot_category(directory, "MEM", "logPbb", "New_rec_std_cat3_H", "SL-Cat 3 (H)", "", "log(L_{S}^{b-tag})", 0);
  plot_category(directory, "MEM", "logPbb", "New_rec_std_cat6_L", "DL (L)", "", "log(L_{S}^{b-tag})", 0);
  plot_category(directory, "MEM", "logPbb", "New_rec_std_cat6_H", "DL (H)", "", "log(L_{S}^{b-tag})", 0);

  plot_category(directory, "MEM", "logPjj", "New_rec_std_cat1_L", "SL-Cat 1 (L)", "", "log(L_{B}^{b-tag})", 0);
  plot_category(directory, "MEM", "logPjj", "New_rec_std_cat1_H", "SL-Cat 1 (H)", "", "log(L_{B}^{b-tag})", 0);
  plot_category(directory, "MEM", "logPjj", "New_rec_std_cat2_L", "SL-Cat 2 (L)", "", "log(L_{B}^{b-tag})", 0);
  plot_category(directory, "MEM", "logPjj", "New_rec_std_cat2_H", "SL-Cat 2 (H)", "", "log(L_{B}^{b-tag})", 0);
  plot_category(directory, "MEM", "logPjj", "New_rec_std_cat3_L", "SL-Cat 3 (L)", "", "log(L_{B}^{b-tag})", 0);
  plot_category(directory, "MEM", "logPjj", "New_rec_std_cat3_H", "SL-Cat 3 (H)", "", "log(L_{B}^{b-tag})", 0);
  plot_category(directory, "MEM", "logPjj", "New_rec_std_cat6_L", "DL (L)", "", "log(L_{B}^{b-tag})", 0);
  plot_category(directory, "MEM", "logPjj", "New_rec_std_cat6_H", "DL (H)", "", "log(L_{B}^{b-tag})", 0);

}



void plot_comp(TString file_1 = "",
	       TString file_2 = "",
	       TString cat_1  = "cat1",
	       TString cat_2  = "cat1",
	       TString proc_1 = "TTH125",
	       TString proc_2 = "TTH125",
	       TString syst_1 = "nominal",
	       TString syst_2 = "nominal",
	       TString header_1 = "",
	       TString header_2 = "",
	       TString header = "",
	       TString fname  = "",
	       TString titleX = "P_{s/b}"
	       ){
  
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

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TLegend* leg = new TLegend(0.35913,0.695804,0.609907,0.884615,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04); 


  TH1F* h_1 = 0;
  TFile* f_1 = TFile::Open(file_1);
  if(f_1==0 || f_1->IsZombie() ){
    cout << "Could not open file1=" << string(file_1.Data()) << endl;
    return;
  }
  TString syst_1_nickname = string(syst_1.Data()).find("nominal")!=string::npos ? "" : "_"+syst_1;
  cout << string((cat_1+"/"+proc_1+syst_1_nickname).Data()) << endl;
  h_1 = (TH1F*)f_1->Get(cat_1+"/"+proc_1+syst_1_nickname);

  if(h_1==0){
    cout << "Missing histo" << endl;
    return;
  }
  h_1->SetLineWidth(3);
  h_1->SetLineStyle(kSolid);
  h_1->SetLineColor(kBlack);
  leg->AddEntry( h_1, header_1 , "L");


  TH1F* h_2 = 0;
  TFile* f_2 = TFile::Open(file_2);
  if(f_2==0 || f_2->IsZombie() ){
    cout << "Could not open file1=" << string(file_2.Data()) << endl;
    return;
  }
  TString syst_2_nickname = string(syst_2.Data()).find("nominal")!=string::npos ? "" : "_"+syst_2;
  cout << string((cat_2+"/"+proc_2+syst_2_nickname).Data()) << endl;
  h_2 = (TH1F*)f_2->Get(cat_2+"/"+proc_2+syst_2_nickname);
  if(h_2==0){
    cout << "Missing histo" << endl;
    return;
  }
  h_2->SetLineWidth(3);
  h_2->SetLineColor(kRed);
  leg->AddEntry( h_2, header_2 , "F");

  h_2->SetFillStyle(3004);
  h_2->SetFillColor(kRed);

  h_1->SetTitle("CMS Simulation #sqrt{s}=8 TeV");
  h_1->SetTitleSize(0.05,"X");
  h_1->SetTitleSize(0.05,"Y");
  h_1->SetXTitle( titleX );
  h_1->GetXaxis()->SetTitleOffset(0.85);;
  h_1->SetYTitle("units");


  h_1->Scale( 1./h_1->Integral() );
  h_2->Scale( 1./h_2->Integral() );

  h_1->SetMaximum( TMath::Max( float(h_1->GetMaximum()), float(h_2->GetMaximum()) )*1.4 );
  h_1->SetMinimum(0);
  h_1->Draw("HISTE");
  h_2->Draw("HISTESAME");

  cout << "*********" << endl;
  cout << "h_1 := " << h_1->Integral() << ", mu = " << h_1->GetMean() << endl;
  cout << "h_2 := " << h_2->Integral() << ", mu = " << h_2->GetMean() << endl;
  cout << "*********" << endl;

  leg->SetHeader( header );
  leg->Draw();

  if(1){
    c1->SaveAs("Plots/PAS/Plot_Comp_"+fname+".pdf");
    c1->SaveAs("Plots/PAS/Plot_Comp_"+fname+".png");
    delete c1;
  }

}



void plot_compAll(){


  TString dir = "May25_2014_195fb";


  plot_comp("datacards/PreApproval/controlPlots/MEM_New_rec_std_cat1.root", 
	    "datacards/PreApproval/controlPlots/MEM_New_rec_std_cat1.root", 
	    "MEM_best_1_pt_b1", "MEM_best_1_pt_b1", 
	    "ttbarPlusBBbar", "ttbar", 
	    "nominal", "nominal", 
	    "tt+bb","tt+jj", 
	    "SL Cat 1", "cat1_pt_b1_TTJetsHFbb_vs_TTJetsLF");

  plot_comp("datacards/PreApproval/controlPlots/MEM_New_rec_std_cat1.root", 
	    "datacards/PreApproval/controlPlots/MEM_New_rec_std_cat1.root", 
	    "MEM_best_1_pt_b2", "MEM_best_1_pt_b2", 
	    "ttbarPlusBBbar", "ttbar", 
	    "nominal", "nominal", 
	    "tt+bb","tt+jj", 
	    "SL Cat 1", "cat1_pt_b2_TTJetsHFbb_vs_TTJetsLF");

  plot_comp("datacards/PreApproval/controlPlots/MEM_New_rec_std_cat1.root", 
	    "datacards/PreApproval/controlPlots/MEM_New_rec_std_cat1.root", 
	    "MEM_best_1_eta_b1", "MEM_best_1_eta_b1", 
	    "ttbarPlusBBbar", "ttbar", 
	    "nominal", "nominal", 
	    "tt+bb","tt+jj", 
	    "SL Cat 1", "cat1_eta_b1_TTJetsHFbb_vs_TTJetsLF");

  plot_comp("datacards/PreApproval/controlPlots/MEM_New_rec_std_cat1.root", 
	    "datacards/PreApproval/controlPlots/MEM_New_rec_std_cat1.root", 
	    "MEM_best_1_eta_b2", "MEM_best_1_eta_b2", 
	    "ttbarPlusBBbar", "ttbar", 
	    "nominal", "nominal", 
	    "tt+bb","tt+jj", 
	    "SL Cat 1", "cat1_eta_b2_TTJetsHFbb_vs_TTJetsLF");

  return;

  for(int cat = 1; cat<7 ; cat++){
    
    cout << "Doing............." << cat << endl;

    if( cat>3 && cat<6 ) continue;
    
    string catNameL = "";
    string catNameH = "";
    if( cat<=3 ){
      catNameL = string(Form("SL-Cat. %d (L)", cat)) ;
      catNameH = string(Form("SL-Cat. %d (H)", cat)) ;
    }
    else{
      catNameL = "DL (L)" ;
      catNameH = "DL (H)" ;
    }

    plot_comp("datacards/"+dir+"/MEM_New_rec_std_sb.root" , "datacards/"+dir+"/MEM_New_rec_std_sb.root" , 
	      Form("MEM_cat%d_H",cat),    Form("MEM_cat%d_H",cat), 
	      "TTJetsHFbb", "TTH125", 
	      "nominal", "nominal", 
	      "t#bar{t}+bb", "t#bar{t}H",
	      catNameH.c_str(),
	      Form("cat%d_H_TTH125_vs_TTJetsHFbb", cat));

 
    plot_comp("datacards/"+dir+"/MEM_New_rec_std_sb.root" , "datacards/"+dir+"/MEM_New_rec_std_sb.root" , 
	      Form("MEM_cat%d_L",cat),    Form("MEM_cat%d_L",cat), 
	      "TTJetsHFbb", "TTH125", 
	      "nominal", "nominal", 
	      "t#bar{t}+bb", "t#bar{t}H",
	      catNameL.c_str(),
	      Form("cat%d_L_TTH125_vs_TTJetsHFbb", cat));


    plot_comp("datacards/"+dir+"/MEM_New_rec_std_sb.root" , "datacards/"+dir+"/MEM_New_rec_std_sb.root" , 
	      Form("MEM_cat%d_H",cat),    Form("MEM_cat%d_H",cat), 
	      "TTJetsLF", "TTH125", 
	      "nominal", "nominal", 
	      "t#bar{t}+jj", "t#bar{t}H",
	      catNameH.c_str(),
	      Form("cat%d_H_TTH125_vs_TTJetsLF", cat));

 
    plot_comp("datacards/"+dir+"/MEM_New_rec_std_sb.root" , "datacards/"+dir+"/MEM_New_rec_std_sb.root" , 
	      Form("MEM_cat%d_L",cat),    Form("MEM_cat%d_L",cat), 
	      "TTJetsLF", "TTH125", 
	      "nominal", "nominal", 
	      "t#bar{t}+jj", "t#bar{t}H",
	      catNameL.c_str(),
	      Form("cat%d_L_TTH125_vs_TTJetsLF", cat));




    plot_comp("datacards/"+dir+"/MEM_New_rec_std_sb.root" , "datacards/"+dir+"/MEM_New_rec_std_sb.root" , 
	      Form("MEM_cat%d_H",cat),    Form("MEM_cat%d_H",cat), 
	      "TTJetsHFb", "TTH125", 
	      "nominal", "nominal", 
	      "t#bar{t}+b", "t#bar{t}H",
	      catNameH.c_str(),
	      Form("cat%d_H_TTH125_vs_TTJetsHFb", cat));

 
    plot_comp("datacards/"+dir+"/MEM_New_rec_std_sb.root" , "datacards/"+dir+"/MEM_New_rec_std_sb.root" , 
	      Form("MEM_cat%d_L",cat),    Form("MEM_cat%d_L",cat), 
	      "TTJetsHFb", "TTH125", 
	      "nominal", "nominal", 
	      "t#bar{t}+b", "t#bar{t}H",
	      catNameL.c_str(),
	      Form("cat%d_L_TTH125_vs_TTJetsHFb", cat));






    plot_comp("datacards/"+dir+"/MEM_New_rec_std_sb.root" , "datacards/"+dir+"/MEM_New_rec_std_sb.root" , 
	      Form("MEM_cat%d_H",cat),    Form("MEM_cat%d_H",cat), 
	      "TTJetsHFb", "TTJetsHFbb", 
	      "nominal", "nominal", 
	      "t#bar{t}+b", "t#bar{t}+bb",
	      catNameH.c_str(),
	      Form("cat%d_H_TTJetsHFbb_vs_TTJetsHFb", cat));

 
    plot_comp("datacards/"+dir+"/MEM_New_rec_std_sb.root" , "datacards/"+dir+"/MEM_New_rec_std_sb.root" , 
	      Form("MEM_cat%d_L",cat),    Form("MEM_cat%d_L",cat), 
	      "TTJetsHFb", "TTJetsHFbb", 
	      "nominal", "nominal", 
	      "t#bar{t}+b", "t#bar{t}+bb",
	      catNameL.c_str(),
	      Form("cat%d_L_TTJetsHFbb_vs_TTJetsHFb", cat));





    plot_comp("datacards/"+dir+"/MEM_New_rec_std_sb.root" , "datacards/"+dir+"/MEM_New_rec_std_sb.root" , 
	      Form("MEM_cat%d_H",cat),    Form("MEM_cat%d_H",cat), 
	      "TTJetsLF", "TTJetsHFbb", 
	      "nominal", "nominal", 
	      "t#bar{t}+jj", "t#bar{t}+bb",
	      catNameH.c_str(),
	      Form("cat%d_H_TTJetsHFbb_vs_TTJetsLF", cat));

    plot_comp("datacards/"+dir+"/MEM_New_rec_std_sb.root" , "datacards/"+dir+"/MEM_New_rec_std_sb.root" , 
	      Form("MEM_cat%d_L",cat),    Form("MEM_cat%d_L",cat), 
	      "TTJetsLF", "TTJetsHFbb", 
	      "nominal", "nominal", 
	      "t#bar{t}+jj", "t#bar{t}+bb",
	      catNameL.c_str(),
	      Form("cat%d_L_TTJetsHFbb_vs_TTJetsLF", cat));


    plot_comp("datacards/"+dir+"/MEM_New_rec_std_sb.root" , "datacards/"+dir+"/MEM_New_rec_std_sb.root" , 
	      Form("MEM_cat%d_H",cat),    Form("MEM_cat%d_H",cat), 
	      "TTJetsLFcc", "TTJetsHFbb", 
	      "nominal", "nominal", 
	      "t#bar{t}+cc", "t#bar{t}+bb",
	      catNameH.c_str(),
	      Form("cat%d_H_TTJetsHFbb_vs_TTJetsLFcc", cat));

    plot_comp("datacards/"+dir+"/MEM_New_rec_std_sb.root" , "datacards/"+dir+"/MEM_New_rec_std_sb.root" , 
	      Form("MEM_cat%d_L",cat),    Form("MEM_cat%d_L",cat), 
	      "TTJetsLFcc", "TTJetsHFbb", 
	      "nominal", "nominal", 
	      "t#bar{t}+cc", "t#bar{t}+bb",
	      catNameL.c_str(),
	      Form("cat%d_L_TTJetsHFbb_vs_TTJetsLFcc", cat));

    plot_comp("datacards/"+dir+"/MEM_New_rec_std_sb.root" , "datacards/"+dir+"/MEM_New_rec_std_sb.root" , 
	      Form("MEM_cat%d_H",cat),    Form("MEM_cat%d_H",cat), 
	      "TTJetsLF", "TTJetsLFcc", 
	      "nominal", "nominal", 
	      "t#bar{t}+jj", "t#bar{t}+cc",
	      catNameH.c_str(),
	      Form("cat%d_H_TTJetsLFcc_vs_TTJetsLF", cat));

    plot_comp("datacards/"+dir+"/MEM_New_rec_std_sb.root" , "datacards/"+dir+"/MEM_New_rec_std_sb.root" , 
	      Form("MEM_cat%d_L",cat),    Form("MEM_cat%d_L",cat), 
	      "TTJetsLF", "TTJetsLFcc", 
	      "nominal", "nominal", 
	      "t#bar{t}+jj", "t#bar{t}+cc",
	      catNameL.c_str(),
	      Form("cat%d_L_TTJetsLFcc_vs_TTJetsLF", cat));


    plot_comp("datacards/"+dir+"/MEM_New_rec_std_sb_nb.root" , "datacards/"+dir+"/MEM_New_rec_std_sb_nb.root" , 
	      Form("MEM_cat%d",cat),    Form("MEM_cat%d",cat), 
	      "TTJetsHFbb", "TTH125", 
	      "nominal", "nominal", 
	      "t#bar{t}+bb", "t#bar{t}H",
	      catNameH.c_str(),
	      Form("cat%d_H_TTH125_vs_TTJetsHFbb_nb", cat),
	      "P_{s/b}^{1 perm}");
 
    plot_comp("datacards/"+dir+"/MEM_New_rec_std_sb.root" , "datacards/"+dir+"/MEM_New_rec_std_sb.root" , 
	      Form("MEM_cat%d_L",cat),    Form("MEM_cat%d_L",cat), 
	      "SingleT", "TTH125", 
	      "nominal", "nominal", 
	      "Single-t", "t#bar{t}H",
	      catNameL.c_str(),
	      Form("cat%d_L_TTH125_vs_SingleT", cat));

    plot_comp("datacards/"+dir+"/MEM_New_rec_std_sb.root" , "datacards/"+dir+"/MEM_New_rec_std_sb.root" , 
	      Form("MEM_cat%d_H",cat),    Form("MEM_cat%d_H",cat), 
	      "SingleT", "TTH125", 
	      "nominal", "nominal", 
	      "Single-t", "t#bar{t}H",
	      catNameH.c_str(),
	      Form("cat%d_H_TTH125_vs_SingleT", cat));

    plot_comp("datacards/"+dir+"/MEM_New_rec_std_sb.root" , "datacards/"+dir+"/MEM_New_rec_std_sb.root" , 
	      Form("MEM_cat%d_L",cat),    Form("MEM_cat%d_L",cat), 
	      "TTV", "TTH125", 
	      "nominal", "nominal", 
	      "t#bar{t}V", "t#bar{t}H",
	      catNameL.c_str(),
	      Form("cat%d_L_TTH125_vs_TTV", cat));

    plot_comp("datacards/"+dir+"/MEM_New_rec_std_sb.root" , "datacards/"+dir+"/MEM_New_rec_std_sb.root" , 
	      Form("MEM_cat%d_H",cat),    Form("MEM_cat%d_H",cat), 
	      "TTV", "TTH125", 
	      "nominal", "nominal", 
	      "t#bar{t}V", "t#bar{t}H",
	      catNameH.c_str(),
	      Form("cat%d_H_TTH125_vs_TTV", cat));

  

  }


}




void plot_comp_nice(TString file_1 = "",
		    TString file_2 = "",
		    TString cat_1  = "cat1",
		    TString cat_2  = "cat1",
		    TString proc_1 = "TTH125",
		    TString proc_2 = "TTH125",
		    TString syst_1 = "nominal",
		    TString syst_2 = "nominal",
		    TString header_1 = "",
		    TString header_2 = "",
		    TString header = "",
		    TString fname  = "",
		    TString titleX = "P_{s/b}",
		    int mergeFirstBins = 1
		    ){
  
  gStyle->SetOptStat(0);

  gStyle->SetOptStat(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleH(0.04);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleOffset(1.3,"y");

  TCanvas *c1 = new TCanvas("c1","",5,30,640,580);
  //c1->SetGrid(0,0);
  //c1->SetFillStyle(4000);
  //c1->SetFillColor(10);
  //c1->SetTicky();
  //c1->SetObjectStat(0);
  TPad *p1 = new TPad("p1","",0,0,1,1);
  p1->SetGrid(0,0);
  p1->SetFillStyle(4000);
  p1->SetFillColor(10);
  p1->SetTicky();
  p1->SetObjectStat(0);
  p1->SetLogy(0);
  p1->Draw();
  p1->cd();
  p1->SetTopMargin(0.05);

  TString cmsinfo = "";
  if(RUNONDATA==0) cmsinfo = "Simulation                                  #sqrt{s}=8 TeV, L=19.5 fb^{-1}";
  else cmsinfo = "CMS Preliminary                                                  #sqrt{s}=8 TeV";

  TPaveText *pt_title = new TPaveText(0.1, 0.952, 0.9, 1.0,"brNDC");
  pt_title->SetFillStyle(1001);
  pt_title->SetBorderSize(0);
  pt_title->SetFillColor(0);
  pt_title->SetTextFont(42); 
  pt_title->SetTextSize(0.04); 
  pt_title->AddText(cmsinfo);

  TLegend* leg = new TLegend(0.605,0.622,0.896,0.888,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.05); 


  TH1F* h_1 = 0;
  TFile* f_1 = TFile::Open(file_1);
  if(f_1==0 || f_1->IsZombie() ){
    cout << "Could not open file1=" << string(file_1.Data()) << endl;
    return;
  }
  TString syst_1_nickname = string(syst_1.Data()).find("nominal")!=string::npos ? "" : "_"+syst_1;
  cout << string((cat_1+"/"+proc_1+syst_1_nickname).Data()) << endl;
  h_1 = (TH1F*)f_1->Get(cat_1+"/"+proc_1+syst_1_nickname);

  if(h_1==0){
    cout << "Missing histo" << endl;
    return;
  }
  h_1->SetLineWidth(1);
  h_1->SetLineStyle(kSolid);
  h_1->SetLineColor(kBlack);


  TH1F* h_2 = 0;
  TFile* f_2 = TFile::Open(file_2);
  if(f_2==0 || f_2->IsZombie() ){
    cout << "Could not open file1=" << string(file_2.Data()) << endl;
    return;
  }
  TString syst_2_nickname = string(syst_2.Data()).find("nominal")!=string::npos ? "" : "_"+syst_2;
  cout << string((cat_2+"/"+proc_2+syst_2_nickname).Data()) << endl;
  h_2 = (TH1F*)f_2->Get(cat_2+"/"+proc_2+syst_2_nickname);
  if(h_2==0){
    cout << "Missing histo" << endl;
    return;
  }
  h_2->SetLineWidth(1);
  h_2->SetLineColor(kRed);

  //h_2->SetFillStyle(3004);
  //h_2->SetFillColor(kRed);

  if( mergeFirstBins ){

    int nBins = h_1->GetNbinsX();

    TH1F* h_1_new = new TH1F( "h_1_new", "", nBins-1, 0, 1);
    for(int b = 1; b <= nBins; b++){
      if(b<=2){
	float err1 = h_1->GetBinError(1);
	float err2 = h_1->GetBinError(2);
	float bin1 = h_1->GetBinContent(1);
	float bin2 = h_1->GetBinContent(2);
	h_1_new->SetBinContent(b, bin1+bin2 );
	h_1_new->SetBinError  (b, sqrt(err1*err1 + err2*err2) );
      }
      else{
	float err1 = h_1->GetBinError(b); 
	float bin1 = h_1->GetBinContent(b);
	h_1_new->SetBinContent(b-1, bin1 );
	h_1_new->SetBinError  (b-1, err1 );
      }
    }

    TH1F* h_2_new = new TH1F( "h_2_new", "", nBins-1, 0, 1);
    for(int b = 1; b <= nBins; b++){
      if(b<=2){
	float err1 = h_2->GetBinError(1);
	float err2 = h_2->GetBinError(2);
	float bin1 = h_2->GetBinContent(1);
	float bin2 = h_2->GetBinContent(2);
	h_2_new->SetBinContent(b, bin1+bin2 );
	h_2_new->SetBinError  (b, sqrt(err1*err1 + err2*err2) );
      }
      else{
	float err1 = h_2->GetBinError(b); 
	float bin1 = h_2->GetBinContent(b);
	h_2_new->SetBinContent(b-1, bin1 );
	h_2_new->SetBinError  (b-1, err1 );
      }
    }

    h_1 = h_1_new;
    h_2 = h_2_new;

    h_1->SetLineWidth(3);
    h_1->SetLineStyle(kDashed);
    h_1->SetLineColor(kBlue);    
    //h_1->SetFillStyle(3005);
    //h_1->SetFillColor(kBlue);
    h_2->SetLineWidth(3);
    h_2->SetLineColor(kRed);    
    //h_2->SetFillStyle(3004);
    //h_2->SetFillColor(kRed);
  }

  leg->AddEntry( h_1, header_1 , "L");
  leg->AddEntry( h_2, header_2 , "L");


  //h_1->SetTitle("CMS Simulation #sqrt{s}=8 TeV");
  h_1->SetTitleSize(0.05,"X");
  h_1->SetTitleSize(0.05,"Y");
  h_1->SetXTitle( titleX );
  h_1->GetXaxis()->SetTitleOffset(0.85);
  h_1->GetYaxis()->SetTitleOffset(0.88);
  h_1->SetYTitle("normalized to unity");

  h_1->Scale( 1./h_1->Integral() );
  h_2->Scale( 1./h_2->Integral() );

  //h_1->SetMaximum( TMath::Max( float(h_1->GetMaximum()), float(h_2->GetMaximum()) )*1.4 );h_1->SetMaximum( TMath::Max( float(h_1->GetMaximum()), float(h_2->GetMaximum()) )*1.4 );
  h_1->SetMaximum( 0.7 );
  h_1->SetMinimum(0);

  h_1->GetXaxis()->SetTitleFont(62);
  h_1->GetYaxis()->SetTitleFont(62);


  h_1->Draw("HIST");
  h_2->Draw("HISTSAME");
  
  TH1D* h_1_clone = (TH1D*)h_1->Clone("h_1_clone");
  h_1_clone  ->SetFillColor(kBlue);
  h_1_clone  ->SetFillStyle( 3001 );
  h_1_clone  ->SetLineWidth(1);
  h_1_clone  ->Draw("E2SAME");
  
  TH1D* h_2_clone = (TH1D*)h_2->Clone("h_2_clone");
  h_2_clone  ->SetFillColor(kRed);
  h_2_clone  ->SetFillStyle( 3001 );
  h_2_clone  ->SetLineWidth(1);
  h_2_clone  ->Draw("E2SAME");

  TH1D* h_3_clone = (TH1D*)h_1->Clone("h_3_clone");
  h_3_clone->Reset();
  h_3_clone  ->SetLineWidth(1);
  h_3_clone  ->SetLineStyle(kSolid);
  h_3_clone  ->SetLineColor(kBlack);
  h_3_clone  ->SetFillColor(kBlack);
  h_3_clone  ->SetFillStyle( 3001 );
  leg->AddEntry( h_3_clone, "MC stat." , "F");

  cout << "*********" << endl;
  cout << "h_1 := " << h_1->Integral() << ", mu = " << h_1->GetMean() << endl;
  cout << "h_2 := " << h_2->Integral() << ", mu = " << h_2->GetMean() << endl;
  cout << "*********" << endl;

  leg->SetHeader( header );
  leg->Draw();

  pt_title->Draw();


  if(1){
    c1->SaveAs("Plots/PAS/Plot_Comp_"+fname+"_nice.pdf");
    c1->SaveAs("Plots/PAS/Plot_Comp_"+fname+"_nice.png");
    delete c1;
  }

}

void plot_comp_niceAll(){

  plot_comp_nice("datacards/PreApproval/MEM_New_rec_std_sb.root", "datacards/PreApproval/MEM_New_rec_std_sb.root", 
		 "MEM_cat1_H", "MEM_cat1_H", 
		 "ttbarPlusBBbar", "ttH_hbb", 
		 "nominal", "nominal", 
		 "t#bar{t}+bb","t#bar{t}H", 
		 "SL Cat-1 (H)", "cat1_H_TTH125_vs_TTJetsHFbb", "P_{s/b}", 1);

  plot_comp_nice("datacards/PreApproval/MEM_New_rec_std_sb.root", "datacards/PreApproval/MEM_New_rec_std_sb.root", 
		 "MEM_cat2_H", "MEM_cat2_H", 
		 "ttbarPlusBBbar", "ttH_hbb", 
		 "nominal", "nominal", 
		 "t#bar{t}+bb","t#bar{t}H", 
		 "SL Cat-2 (H)", "cat2_H_TTH125_vs_TTJetsHFbb", "P_{s/b}", 1);

  plot_comp_nice("datacards/PreApproval/MEM_New_rec_std_sb.root", "datacards/PreApproval/MEM_New_rec_std_sb.root", 
		 "MEM_cat3_H", "MEM_cat3_H", 
		 "ttbarPlusBBbar", "ttH_hbb", 
		 "nominal", "nominal", 
		 "t#bar{t}+bb","t#bar{t}H", 
		 "SL Cat-3 (H)", "cat3_H_TTH125_vs_TTJetsHFbb", "P_{s/b}", 1);

  plot_comp_nice("datacards/PreApproval/MEM_New_rec_std_sb.root", "datacards/PreApproval/MEM_New_rec_std_sb.root", 
		 "MEM_cat6_H", "MEM_cat6_H", 
		 "ttbarPlusBBbar", "ttH_hbb", 
		 "nominal", "nominal", 
		 "t#bar{t}+bb","t#bar{t}H", 
		 "DL (H)", "cat6_H_TTH125_vs_TTJetsHFbb", "P_{s/b}", 1);
}




void plot_limit(string version = "_New", 
		int doExpOnly = 1, int doMLFit = 0,
		float minY = 0.8, float maxY=40, 
		int compare = 1,
		string label = "",
		int plotCombined = 1,
		int doSignalInjected = 0
		){

  gStyle->SetPaintTextFormat("g");
  gStyle->SetOptStat(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetTitleFillColor(0);
  gStyle->SetTitleBorderSize(0);
  gStyle->SetTitleH(0.04);
  gStyle->SetTitleStyle(0);
  gStyle->SetTitleOffset(1.3,"y");

  TCanvas *c1 = new TCanvas("c1","",5,30,640,580);
  /*
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  c1->SetLogy(!doMLFit);
  */

  TPad *p1 = new TPad("p1","",0,0,1,1);
  p1->SetGrid(0,0);
  p1->SetFillStyle(4000);
  p1->SetFillColor(10);
  p1->SetTicky();
  p1->SetObjectStat(0);
  p1->SetLogy(!doMLFit);
  p1->Draw();
  p1->cd();
  p1->SetTopMargin(0.05);

  TString cmsinfo = "";
  if(RUNONDATA==0) cmsinfo = "Simulation                                  #sqrt{s}=8 TeV, L=19.5 fb^{-1}";
  else cmsinfo = "CMS Preliminary                           #sqrt{s}=8 TeV, L=19.5 fb^{-1}";

  TPaveText *pt_title = new TPaveText(0.1, 0.952, 0.9, 1.0,"brNDC");
  pt_title->SetFillStyle(1001);
  pt_title->SetBorderSize(0);
  pt_title->SetFillColor(0);
  pt_title->SetTextFont(42); 
  pt_title->SetTextSize(0.04); 
  pt_title->AddText(cmsinfo);
  
  TLegend* leg = new TLegend(0.12,0.17,0.32, 0.36, NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->SetFillColor(0);


  float X[]        = {0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5, 13.5, 14.5 };
  float expY[]     = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float obsY[]     = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float expY1sL[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float expY1sH[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float expY2sL[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float expY2sH[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

  float expSIY[]     = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float obsSIY[]     = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float expSIY1sL[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float expSIY1sH[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float expSIY2sL[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float expSIY2sH[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

  float expXs[]  = {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5};
  float expYs[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};


  vector<string> categories;        vector<string> names;
  
  if( plotCombined==0 ){
    categories.push_back("MEM_cat1_H"+version);  names.push_back("SL Cat-1 (H)");
    categories.push_back("MEM_cat2_H"+version);  names.push_back("SL Cat-2 (H)");
    categories.push_back("MEM_cat3_H"+version);  names.push_back("SL Cat-3 (H)");
    categories.push_back("MEM_cat6_H"+version);  names.push_back("DL (H)"); 
    categories.push_back("MEM_cat1_L"+version);  names.push_back("SL Cat-1 (L)");
    categories.push_back("MEM_cat2_L"+version);  names.push_back("SL Cat-2 (L)");
    categories.push_back("MEM_cat3_L"+version);  names.push_back("SL Cat-3 (L)");
    categories.push_back("MEM_cat6_L"+version);  names.push_back("DL (L)");
    categories.push_back("MEM_cat1"+version);    names.push_back("SL Cat-1");
    categories.push_back("MEM_cat2"+version);    names.push_back("SL Cat-2");
    categories.push_back("MEM_cat3"+version);    names.push_back("SL Cat-3");
    categories.push_back("MEM_cat6"+version);    names.push_back("DL"); 
    categories.push_back("MEM_SL"+version);      names.push_back("SL comb.");
    categories.push_back("MEM_DL"+version);      names.push_back("DL comb."); 
    categories.push_back("MEM_COMB"+version);    names.push_back("All comb.");
  }
  
  else{
    categories.push_back("MEM_cat1"+version);  names.push_back("SL Cat-1");
    categories.push_back("MEM_cat2"+version);  names.push_back("SL Cat-2");
    categories.push_back("MEM_cat3"+version);  names.push_back("SL Cat-3");
    categories.push_back("MEM_cat6"+version);  names.push_back("DL"); 
    //categories.push_back("MEM_SL"+version);    names.push_back("SL comb.");
    //categories.push_back("MEM_DL"+version);    names.push_back("DL comb."); 
    categories.push_back("MEM_COMB"+version);  names.push_back("All comb.");
  }

  int nBins = categories.size();

  for( int b = 0; b < nBins; b++){

    TFile* f = doMLFit ? 
      TFile::Open(("datacards/PreApproval//higgsCombine"+categories[b]+".MaxLikelihoodFit.mH120.root").c_str()) :
      TFile::Open(("datacards/PreApproval//higgsCombine"+categories[b]+".Asymptotic.mH120.root").c_str());
    if( f==0 ) continue;
    
    TFile* f2 = 0;
    if( doSignalInjected ){
      f2 = TFile::Open(("datacards/PreApproval//higgsCombine"+categories[b]+"_MC-sb_signalInjected.Asymptotic.mH120.root").c_str());
      if( f2==0 ) continue;
    }

    Double_t r;
    TTree* limit = (TTree*)f->Get("limit");
    limit->SetBranchAddress("limit",&r);

    for(int k = 0 ; k< limit->GetEntries() ; k++){
      limit->GetEntry(k);      

      string val(Form("%.1f",r));
      r = atof(val.c_str());

      //cout << r << endl;

      if(!doMLFit){
	if(k==0) expY2sL[b] = r;
	if(k==1) expY1sL[b] = r;
	if(k==2) expY[b]    = r;
	if(k==3) expY1sH[b] = r;
	if(k==4) expY2sH[b] = r;
	if(k==5) obsY[b]    = r;
      }
      else{
	if(k==0) expY[b]    = r;       
	if(k==1) expY1sL[b] = r;
	if(k==2) expY1sH[b] = r;	
      }
    }



    ////////////////////////////////////////////////////////////////////////////////////////
    if( doSignalInjected ){

      limit = (TTree*)f2->Get("limit");
      limit->SetBranchAddress("limit",&r);
      
      for(int k = 0 ; k< limit->GetEntries() ; k++){
	limit->GetEntry(k);      
	
	string val(Form("%.1f",r));
	r = atof(val.c_str());
	
	if(k==0) expSIY2sL[b] = r;
	if(k==1) expSIY1sL[b] = r;
	if(k==2) expSIY[b]    = r;
	if(k==3) expSIY1sH[b] = r;
	if(k==4) expSIY2sH[b] = r;
	if(k==5) obsSIY[b]    = r;
     
      }
    }
    ////////////////////////////////////////////////////////////////////////////////////////

    std::cout.precision(2);
    if( !doSignalInjected )
      cout << names[b] << " & " << obsY[b] << " & " << expY2sL[b] << " & " <<  expY1sL[b] << " & "
	   <<  expY[b]  << " & " <<  expY1sH[b] << " & " <<  expY2sH[b]  << " \\\\" << endl;
    else
      cout << names[b]   << " & " 
	   << obsY[b]    << " & "
	   << obsSIY[b]  << " & "
	   << expY[b]    << " & "
	   << "[" << expY1sL[b] << ", " << expY1sH[b] << "]"   << " & "
	   << "[" << expY2sL[b] << ", " << expY2sH[b] << "]"
	   << " \\\\" << endl;
	
    expY1sH[b] = TMath::Abs(expY1sH[b]-expY[b]);
    expY1sL[b] = TMath::Abs(expY1sL[b]-expY[b]);
    expY2sH[b] = TMath::Abs(expY2sH[b]-expY[b]);
    expY2sL[b] = TMath::Abs(expY2sL[b]-expY[b]);


    if( doSignalInjected ){
      expSIY1sH[b] = TMath::Abs(expSIY1sH[b]-expSIY[b]);
      expSIY1sL[b] = TMath::Abs(expSIY1sL[b]-expSIY[b]);
      expSIY2sH[b] = TMath::Abs(expSIY2sH[b]-expSIY[b]);
      expSIY2sL[b] = TMath::Abs(expSIY2sL[b]-expSIY[b]);
    }

    f->Close();
    if( f2 ) f2->Close();

  }

  TMultiGraph *mg = new TMultiGraph();
  //mg->SetTitle("CMS Preliminary #sqrt{s}=8 TeV, L=19.5 fb^{-1}");


  TGraphAsymmErrors* observed  = new TGraphAsymmErrors(15, X, obsY, expXs ,expXs,  expYs,   expYs);
  TGraphAsymmErrors* expected  = new TGraphAsymmErrors(15, X, expY, expXs ,expXs,  expYs,   expYs);
  TGraphAsymmErrors* oneSigma  = new TGraphAsymmErrors(15, X, expY, expXs, expXs,  expY1sL, expY1sH);
  TGraphAsymmErrors* twoSigma  = new TGraphAsymmErrors(15, X, expY, expXs, expXs,  expY2sL, expY2sH);
  TGraphAsymmErrors* signalInj = new TGraphAsymmErrors(15, X, obsSIY, expXs ,expXs,  expYs,   expYs);

  oneSigma->SetMarkerColor(kBlack);
  oneSigma->SetMarkerStyle(kOpenCircle);
  oneSigma->SetFillColor(kGreen);
  oneSigma->SetFillStyle(1001);

  twoSigma->SetMarkerColor(kBlack);
  twoSigma->SetMarkerStyle(kOpenCircle);
  twoSigma->SetFillColor(kYellow);
  twoSigma->SetFillStyle(1001);

  expected->SetMarkerColor(kBlack);
  expected->SetMarkerStyle(kOpenCircle);
  expected->SetMarkerSize(1.0);
  expected->SetLineColor(kBlack);
  expected->SetLineWidth(2);

  observed->SetMarkerColor(kBlack);
  observed->SetMarkerStyle(kFullCircle);
  observed->SetMarkerSize(1.0);
  observed->SetLineColor(kBlack);
  observed->SetLineWidth(2);
 

  if(!doMLFit)
    mg->Add(twoSigma);
  mg->Add(oneSigma);
  mg->Add(expected);
  if( doExpOnly==0 && !doMLFit )
    mg->Add(observed);

  mg->Draw("a2");  
  expected->Draw("pSAME");

  if( doExpOnly==0  && !doMLFit)
    observed->Draw("pSAME");


  TH1F* hT = new TH1F("hT", "", nBins, 0, nBins);
  for(int k = 1; k <= hT->GetNbinsX(); k++){
    string val(Form("%.1f",doMLFit ? expY[k-1] : obsY[k-1] ));
    hT->SetBinContent(k,  atof(val.c_str())   );
  }
  hT->SetMarkerSize(2.);
  if( doExpOnly==1 || true)
    hT->Draw("TEXT0SAME");


  TH1F* hTSI = new TH1F("hTSI", "", nBins, 0, nBins);
  if( doSignalInjected ){
    /*
    for(int k = 1; k <= hTSI->GetNbinsX(); k++){
      string val(Form("%.1f", obsSIY[k-1] ));
      hTSI->SetBinContent(k,  atof(val.c_str())   );
      hTSI->SetBinError(k, 0.);
    }
    hTSI->SetMarkerSize(1.0);
    hTSI->SetMarkerStyle(kOpenCircle);
    hTSI->SetLineWidth(2);
    hTSI->SetLineStyle(kDashed);
    hTSI->SetLineColor(kRed);
    hTSI->SetMarkerColor(kRed);
    */

    signalInj->SetMarkerSize(1.0);
    signalInj->SetMarkerStyle(kOpenCircle);
    signalInj->SetLineWidth(2);
    signalInj->SetLineStyle(kDashed);
    signalInj->SetLineColor(kRed);
    signalInj->SetMarkerColor(kRed);
    
    signalInj->Draw("pSAME");
  }


  TF1 *line = new TF1("line","1",0,nBins);
  line->SetLineColor(kRed);
  line->SetLineWidth(2);

  line->Draw("SAME");


  TF1 *line0 = new TF1("line0","0",0,nBins);
  line0->SetLineColor(kBlack);
  line0->SetLineWidth(1);  line0->SetLineStyle(kDashed);
  if( doMLFit ) line0->Draw("SAME");


  TF1 *lineML = new TF1("lineML","2.4",0,nBins);
  lineML->SetLineColor(kBlue);
  lineML->SetLineStyle(kDashed);
  lineML->SetLineWidth(3);


  if(compare) lineML->Draw("SAME");

  TF1 *lineTTH = new TF1("lineTTH","4.1",0,nBins);
  lineTTH->SetLineColor(kMagenta);
  lineTTH->SetLineStyle(kDashed);
  lineTTH->SetLineWidth(3);


  if(compare) lineTTH->Draw("SAME");

  c1->cd();
  p1->cd();
  gPad->Modified();


  mg->GetXaxis()->Set(nBins,0,nBins);


  //mg->GetXaxis()->SetRange(-1,11);


  for( int b = 0; b < nBins; b++){
    mg->GetXaxis()->SetBinLabel(b+1, names[b].c_str() );
  }


  mg->GetXaxis()->SetTitleFont(62);
  mg->GetYaxis()->SetTitleFont(62);
  mg->GetXaxis()->SetTitleSize(0.05);
  mg->GetYaxis()->SetTitleSize(0.05);
  mg->GetYaxis()->SetTitleOffset(0.80);
  mg->GetXaxis()->SetLabelSize(0.06);
  mg->SetMinimum( minY );
  mg->SetMaximum( maxY );
  mg->GetXaxis()->SetTitle("");
  if(!doMLFit) 
    mg->GetYaxis()->SetTitle("95% CL upper limit on #mu = #sigma/#sigma_{SM}");
  else
    mg->GetYaxis()->SetTitle("Best fit of #mu = #sigma/#sigma_{SM}");



  if( doExpOnly==0 && !doMLFit ){
    //leg->AddEntry( expected, "Median exp.", "PL");
    leg->AddEntry( oneSigma, "Exp. 68%", "PFL");
    leg->AddEntry( twoSigma, "Exp. 95%", "PFL");
    if( doSignalInjected ){
      leg->AddEntry( signalInj /*hTSI*/ , "Median exp. (signal injected)", "PL");
    }
    leg->AddEntry( observed, "Observed", "PL");
  }
  else if(!doMLFit)
    leg->AddEntry( expected, "Expected (pre-fit)", "P");
  else{
    //leg->AddEntry( expected, "Observed", "P");
  }


  if(compare){
    leg->AddEntry( lineTTH, "HIG-13-019 (post-fit)", "L");
    leg->AddEntry( lineML,  "HIG-13-020 (post-fit)", "L");
  }

  if(!doMLFit)
    leg->Draw();

  pt_title->Draw();


  TPaveText *pt = new TPaveText(0.106811,0.155594,0.407121,0.286713,"brNDC");
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetFillColor(10);
  pt->SetTextSize(0.04);
  pt->SetTextAlign(11);
  if(!doMLFit){
    pt->AddText(Form("Comb: #mu < %.1f at 95%% CL", expY[categories.size()-1]))->SetTextColor(kRed);
    pt->AddText(Form("(SL: #mu < %.1f)", expY[categories.size()-3]))->SetTextColor(kBlack);
    pt->AddText(Form("(DL: #mu < %.1f)", expY[categories.size()-2]))->SetTextColor(kBlack);
  }

  TLine* del = new TLine( categories.size()-7, maxY, categories.size()-7, minY );
  del->SetLineWidth(3);
  del->SetLineStyle(kSolid);
  del->SetLineColor(kBlack);
  if(plotCombined==0) del->Draw("SAME");

  TLine* del_ = new TLine( categories.size()-3, maxY, categories.size()-3, minY );
  del_->SetLineWidth(3);
  del_->SetLineStyle(kSolid);
  del_->SetLineColor(kBlack);
  if(plotCombined==0) del_->Draw("SAME");

  TLine* del2 = new TLine(categories.size()-1, maxY, categories.size()-1, minY );
  del2->SetLineWidth(3);
  del2->SetLineStyle(kSolid);
  del2->SetLineColor(kBlack);
  del2->Draw("SAME");


  if( 1 ){
    c1->SaveAs("datacards/PreApproval/Limits"+TString(version.c_str())+"_"+TString(label.c_str())+".pdf");
  }

}


void plot_limitAll(){

  // exp limits

  //plot_limit("_New_rec_std_sbprefit",    1  , 0  , 0.8  , 25  , 0, "limit_comb_pre" , 1);
  //plot_limit("_New_rec_std_sb",          0  , 0  , 0.8  , 25  , 0, "limit_comb_post", 1);

  //return;

  plot_limit("_New_rec_std_sb",    0  , 0  , 0.8  , 25  , 0, "limit_comb", 1, 1);
  //plot_limit("_New_rec_std_sb",    0  , 0  , 0.8  , 80  , 0, "limit_all",  0, 1);

  plot_limit("_New_rec_std_sb",    0  , 1  , -10  , 10  , 0, "mlfit_comb", 1);
  //plot_limit("_New_rec_std_sb",    0  , 1  , -25  , 25  , 0, "mlfit_all",  0);

}



void plot_limit_comp(string version = "_New", 
		     int doExpOnly = 1, int doMLFit = 0,
		     float minY = 0.8, float maxY=40, 
		     int compare = 1){

  gStyle->SetPaintTextFormat("g");

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  

  float X[]        = {0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5 , 10.5, 11.5};
  float expY[]     = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float obsY[]     = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float expY1sL[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float expY1sH[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float expY2sL[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float expY2sH[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

  float expXs[]  = {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5};
  float expYs[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};


  vector<string> categories;        vector<string> names;
 
  categories.push_back("MEM_COMB_New_rec_std_sb-STATONLY");  names.push_back("Stat. only");
  categories.push_back("MEM_COMB_New_rec_std_sb-bin");       names.push_back("bin-by-bin");
  categories.push_back("MEM_COMB_New_rec_std_sb-CMS_res_j");  names.push_back("JER");
  categories.push_back("MEM_COMB_New_rec_std_sb-CMS_scale_j");  names.push_back("JES");
  categories.push_back("MEM_COMB_New_rec_std_sb-CSV");  names.push_back("CSV");
  categories.push_back("MEM_COMB_New_rec_std_sb-CMS_ttH_QCDscale_ttcc");  names.push_back("ttcc");
  categories.push_back("MEM_COMB_New_rec_std_sb-CMS_ttH_QCDscale_ttb");  names.push_back("ttb");
  categories.push_back("MEM_COMB_New_rec_std_sb-CMS_ttH_QCDscale_ttbb");  names.push_back("ttbb");
  categories.push_back("MEM_COMB_New_rec_std_sb-Q2");  names.push_back("Q^{2} scale");
  categories.push_back("MEM_COMB_New_rec_std_sb-CMS_ttH_topPtcorr");  names.push_back("top p_{T}");
  categories.push_back("MEM_COMB_New_rec_std_sb-QCDscale_ttH");  names.push_back("ttH scale");
  categories.push_back("MEM_COMB_New_rec_std_sb");  names.push_back("All systematics");

  int nBins = categories.size();

  for( int b = 0; b < nBins; b++){

    TFile* f = doMLFit ? 
      TFile::Open(("datacards/PreApproval//higgsCombine"+categories[b]+".MaxLikelihoodFit.mH120.root").c_str()) :
      TFile::Open(("datacards/PreApproval//higgsCombine"+categories[b]+".Asymptotic.mH120.root").c_str());
    if( f==0 ) continue;
    
    Double_t r;
    TTree* limit = (TTree*)f->Get("limit");
    limit->SetBranchAddress("limit",&r);

    for(int k = 0 ; k< limit->GetEntries() ; k++){
      limit->GetEntry(k);      

      string val(Form("%.1f",r));
      r = atof(val.c_str());

      cout << r << endl;

      if(!doMLFit){
	if(k==0) expY2sL[b] = r;
	if(k==1) expY1sL[b] = r;
	if(k==2) expY[b]    = r;
	if(k==3) expY1sH[b] = r;
	if(k==4) expY2sH[b] = r;
	if(k==5) obsY[b]    = r;
      }
      else{
	if(k==0) expY[b]    = r;
	if(k==1) expY1sL[b] = r;
	if(k==2) expY1sH[b] = r;	
      }
    }

    std::cout.precision(2);
    cout << names[b] << " & " << obsY[b] << " & " << expY2sL[b] << " & " <<  expY1sL[b] << " & "
	 <<  expY[b]  << " & " <<  expY1sH[b] << " & " <<  expY2sH[b]  << " \\\\" << endl;

    expY1sH[b] = TMath::Abs(expY1sH[b]-expY[b]);
    expY1sL[b] = TMath::Abs(expY1sL[b]-expY[b]);
    expY2sH[b] = TMath::Abs(expY2sH[b]-expY[b]);
    expY2sL[b] = TMath::Abs(expY2sL[b]-expY[b]);


  }


  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("CMS Preliminary #sqrt{s}=8 TeV, L=19.5 fb^{-1}");

  TGraphAsymmErrors* observed = new TGraphAsymmErrors(12, X, obsY, expXs ,expXs,  expYs,   expYs);
  TGraphAsymmErrors* expected = new TGraphAsymmErrors(12, X, expY, expXs ,expXs,  expYs,   expYs);
  TGraphAsymmErrors* oneSigma = new TGraphAsymmErrors(12, X, expY, expXs, expXs,  expY1sL, expY1sH);
  TGraphAsymmErrors* twoSigma = new TGraphAsymmErrors(12, X, expY, expXs, expXs,  expY2sL, expY2sH);


  oneSigma->SetMarkerColor(kBlack);
  oneSigma->SetMarkerStyle(kOpenCircle);
  oneSigma->SetFillColor(kGreen);
  oneSigma->SetFillStyle(1001);

  twoSigma->SetMarkerColor(kBlack);
  twoSigma->SetMarkerStyle(kOpenCircle);
  twoSigma->SetFillColor(kYellow);
  twoSigma->SetFillStyle(1001);

  expected->SetMarkerColor(kBlack);
  expected->SetMarkerStyle(kOpenCircle);
  expected->SetMarkerSize(1.0);
  expected->SetLineColor(kBlack);
  expected->SetLineWidth(2);

  observed->SetMarkerColor(kBlack);
  observed->SetMarkerStyle(kFullCircle);
  observed->SetMarkerSize(1.0);
  observed->SetLineColor(kBlack);
  observed->SetLineWidth(2);
 
  if(!doMLFit)
    mg->Add(twoSigma);
  mg->Add(oneSigma);
  mg->Add(expected);
  if( doExpOnly==0 && !doMLFit  )
    mg->Add(observed);

  mg->Draw("a2");  
  expected->Draw("pSAME");
  if( doExpOnly==0  && !doMLFit )
    observed->Draw("pSAME");
  

  TH1F* hT = new TH1F("hT", "", nBins, 0, nBins);
  for(int k = 1; k <= hT->GetNbinsX(); k++){
    int approx = int(expY[k-1]*10);
    hT->SetBinContent(k, approx/10.   );
  }
  hT->SetMarkerSize(2.);
  hT->Draw("TEXT0SAME");

  TF1 *line = new TF1("line","1",0,nBins);
  line->SetLineColor(kRed);
  line->SetLineWidth(2);

  c1->cd();
  gPad->Modified();
  mg->GetXaxis()->Set(nBins,0,nBins);

  //mg->GetXaxis()->SetRange(-1,11);

  for( int b = 0; b < nBins; b++){
    mg->GetXaxis()->SetBinLabel(b+1, names[b].c_str() );
  }

  TLine* del2 = new TLine(categories.size()-1, maxY, categories.size()-1, 0. );
  del2->SetLineWidth(3);
  del2->SetLineStyle(kSolid);
  del2->SetLineColor(kBlack);
  del2->Draw("SAME");

  mg->GetXaxis()->SetTitleSize(0.05);
  mg->GetYaxis()->SetTitleSize(0.05);
  mg->GetYaxis()->SetTitleOffset(0.80);
  mg->SetMinimum( minY);
  mg->SetMaximum( maxY );
  mg->GetXaxis()->SetTitle("");
  if(!doMLFit) 
    mg->GetYaxis()->SetTitle("95% CL upper limit on #mu = #sigma/#sigma_{SM}");
  else
    mg->GetYaxis()->SetTitle("ML fit of #mu = #sigma/#sigma_{SM}");

  if(1){
    c1->SaveAs("datacards/PreApproval/Limits_Sytstematics"+TString(version.c_str())+".pdf");
  }

}

void plot_limit_compAll(){

  // exp limits
  plot_limit_comp("_New_rec_std_sb",    1  , 0  , 0.  , 7 , 0);

}




void plot_syst(string input = "MEM_New_rec_std_sb",
	       string type = "MEM",
	       string cat  = "cat1_H",
	       string proc = "TTH125",
	       string syst     = "JES",
	       string systName = "JES",
	       string header = "SL-Cat. 1",
	       string fname  = "",
	       TString titleX = "P_{s/b}"
	       ){
  
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

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600/0.85);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TPad *p1 = new TPad("p1","",0,0.155,1,1);
  p1->SetGrid(0,0);
  p1->SetFillStyle(4000);
  p1->SetFillColor(10);
  p1->SetTicky();
  p1->SetObjectStat(0);
  p1->Draw();
  p1->cd();

  TLegend* leg = new TLegend(0.35913,0.695804,0.609907,0.884615,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04); 

  vector<string> samples;
  samples.push_back( proc );

  TH1F* hN = 0;
  TH1F* hU = 0;
  TH1F* hD = 0;

  TFile* f = TFile::Open(("datacards/PreApproval/"+input+".root").c_str());
  if(f==0 || f->IsZombie() ) return;

  for(unsigned int sample = 0; sample < samples.size(); sample++){

    cout << "Doing " << samples[sample] << endl;

  
    hN = ((TH1F*)f->Get((type+"_"+cat+"/"+samples[sample]).c_str()));
    if( hN==0 ){
      cout << (type+"_"+cat+"/"+samples[sample]) << " not found" << endl;
      continue;
    }
    else{
      cout << "hN= " << hN->Integral() << endl;
    }

    hU = ((TH1F*)f->Get((type+"_"+cat+"/"+samples[sample]+"_"+syst+"Up").c_str()));
    if( hU==0 ){
      cout << (type+"_"+cat+"/"+samples[sample]+"_"+syst+"Up") << " not found" << endl;
      continue;
    }
    else{
      cout << "hU= " << hU->Integral() << endl;
    }

    hD = ((TH1F*)f->Get((type+"_"+cat+"/"+samples[sample]+"_"+syst+"Down").c_str()));
    if( hD==0 ){
      cout << (type+"_"+cat+"/"+samples[sample]+"_"+syst+"Down") << " not found" << endl;
      continue;
    }
    else{
      cout << "hD= " << hD->Integral() << endl;
    }

  }

  if(hN==0 || hU==0 || hD==0){
    cout << "Return!" << endl;
    return;
  }

  hN->SetLineWidth(3);
  hN->SetLineColor(kBlack);
  hN->SetLineStyle(kSolid);

  hU->SetLineWidth(3);
  hU->SetLineColor(kRed);
  hU->SetLineStyle(kDashed);

  hD->SetLineWidth(3);
  hD->SetLineColor(kBlue);
  hD->SetLineStyle(kDashed);


  if(hU->GetMaximum()>hD->GetMaximum()){
    hU->SetMinimum(0.0);
    hU->SetMaximum( hU->GetMaximum()*1.4 );
    hU->GetYaxis()->SetTitle("Events");
    hU->GetXaxis()->SetTitle( titleX );
    hU->SetTitle("Simulation #sqrt{s}=8 TeV, L=19.5 fb^{-1}");
    hU->SetTitleSize  (0.04,"X");
    hU->SetTitleOffset(0.95,"X");
    hU->Draw("HIST");
    hN->Draw("HISTSAME");
    hD->Draw("HISTSAME");
  }
  else{
    hD->SetMinimum(0.0);
    hD->SetMaximum( hD->GetMaximum()*1.4 );
    hD->GetYaxis()->SetTitle("Events");
    hD->GetXaxis()->SetTitle( titleX );
    hD->SetTitle("Simulation #sqrt{s}=8 TeV, L=19.5 fb^{-1}");
    hD->SetTitleSize  (0.04,"X");
    hD->SetTitleOffset(0.95,"X");
    hD->Draw("HIST");
    hN->Draw("HISTSAME");
    hU->Draw("HISTSAME");
  }

  leg->SetHeader(header.c_str());
  leg->AddEntry( hU, (systName+" up").c_str(),      "L");
  leg->AddEntry( hN, "nominal", "L");
  leg->AddEntry( hD, (systName+" down").c_str(),    "L");

  leg->Draw();


  c1->cd();
  TPad *p2 = new TPad("p2","",0,0,1,0.175);
  p2->SetGrid(0,0);
  p2->SetFillStyle(4000);
  p2->SetFillColor(10);
  p2->SetTicky();
  p2->SetObjectStat(0);
  p2->Draw();
  p2->cd();
 
  TF1* one = new TF1("one","1", hN->GetXaxis()->GetXmin(),hN->GetXaxis()->GetXmax());
  one->SetLineColor(kBlack);

  TH1F* hU_copy = (TH1F*)hU->Clone("hU_copy");
  TH1F* hD_copy = (TH1F*)hD->Clone("hD_copy");
  
  hU_copy->Divide( hU_copy , hN, 1.0, 1.0);
  hD_copy->Divide( hD_copy , hN, 1.0, 1.0);

  hU_copy->GetYaxis()->SetRangeUser(0.6,1.4);
  hU_copy->SetTitle(0);
  hU_copy->GetXaxis()->SetTitle(0);
  hU_copy->GetXaxis()->SetLabelSize(0);
  hU_copy->GetXaxis()->SetNdivisions(0);
  hU_copy->GetXaxis()->SetTickLength(0.1);
  
  hU_copy->GetYaxis()->SetTitleSize(0.16);
  hU_copy->GetYaxis()->SetTitleOffset(0.28);
  hU_copy->GetYaxis()->SetTitle("Ratio     ");
  hU_copy->GetYaxis()->SetNdivisions(202);
  hU_copy->GetYaxis()->SetLabelSize(0.16);

  hU_copy->Draw("HIST");
  one->Draw("SAME");
  hD_copy->Draw("HISTSAME");

  if(1){
    c1->SaveAs(  ("Plots/PAS/"+fname+".pdf").c_str() );
    delete c1; delete leg;
    f->Close();
  }

}


void plot_systAll(){

  vector<string> sys;
  vector<string> sysName;

  sys.push_back("csv_sys_CSVHF");       sysName.push_back("CSV-HF");
  sys.push_back("csv_sys_CSVLF");       sysName.push_back("CSV-LF");
  sys.push_back("csv_sys_CSVHFStats1"); sysName.push_back("CSV-HF stat. 1");
  sys.push_back("csv_sys_CSVHFStats2"); sysName.push_back("CSV-HF stat. 2");
  sys.push_back("csv_sys_CSVLFStats1"); sysName.push_back("CSV-LF stat. 1");
  sys.push_back("csv_sys_CSVLFStats2"); sysName.push_back("CSV-LF stat. 2");
  sys.push_back("csv_sys_CSVCErr1");    sysName.push_back("CSV-C Err.1");
  sys.push_back("csv_sys_CSVCErr2");    sysName.push_back("CSV-C Err.2");
  sys.push_back("Q2Scale1p");           sysName.push_back("Q^{2}-scale (1p)");
  sys.push_back("Q2Scale2p");           sysName.push_back("Q^{2}-scale (2p)");
  sys.push_back("Q2Scale3p");           sysName.push_back("Q^{2}-scale (3p)");
  sys.push_back("Q2ScaleHFbb");         sysName.push_back("Q^{2}-scale (bb)");
  sys.push_back("Q2ScaleHFb");          sysName.push_back("Q^{2}-scale (b)");
  sys.push_back("Q2ScaleLFcc");         sysName.push_back("Q^{2}-scale (cc)");
  sys.push_back("Q2ScaleLF");           sysName.push_back("Q^{2}-scale (jj)");
  sys.push_back("TopPt");               sysName.push_back("top p_{T}");
  sys.push_back("JEC");                 sysName.push_back("JES");
  sys.push_back("JER");                 sysName.push_back("JER");
  

  vector<string> samples;
  vector<string> samplesName;
  samples.push_back("TTH125");        samplesName.push_back("t#bar{t}H");
  samples.push_back("TTJetsHFbb");    samplesName.push_back("t#bar{t}+bb");
  samples.push_back("TTJetsHFb");     samplesName.push_back("t#bar{t}+b");
  samples.push_back("TTJetsLFcc");    samplesName.push_back("t#bar{t}+cc");
  samples.push_back("TTJetsLF");      samplesName.push_back("t#bar{t}+jj");
  samples.push_back("TTV");           samplesName.push_back("t#bar{t}V");

  vector<string> cat;
  vector<string> catNames;
  cat.push_back("cat1_H");   catNames.push_back("SL-Cat.1 (H)");
  cat.push_back("cat1_L");   catNames.push_back("SL-Cat.1 (L)");
  cat.push_back("cat2_H");   catNames.push_back("SL-Cat.2 (H)");
  cat.push_back("cat2_L");   catNames.push_back("SL-Cat.2 (L)"); 
  cat.push_back("cat3_H");   catNames.push_back("SL-Cat.3 (H)");
  cat.push_back("cat3_L");   catNames.push_back("SL-Cat.3 (L)");
  cat.push_back("cat6_H");   catNames.push_back("DL (H)");
  cat.push_back("cat6_L");   catNames.push_back("DL (L)");

  for( int i = 0 ; i < sys.size() ; i++){
    for( int j = 0 ; j < samples.size() ; j++){
      for( int k = 0 ; k < cat.size() ; k++){
	plot_syst("MEM_New_rec_std_sb","MEM", cat[k], samples[j], sys[i], sysName[i], catNames[k]+", "+samplesName[j], "Plots_Syst_"+cat[k]+"_"+samples[j]+"_"+sys[i]);
      }
    }
  }

}



void plot_analysis_MW(TString cat   ="SL", 
		      string header="type==3",
		      TString title = "",
		      string fname= "type0",
		      float fact = 1, float fact2 = 1
		      ){

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


  TCanvas *c1 = new TCanvas("c1","",5,30,650,379);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TLegend* leg = new TLegend(0.25,0.65,0.65,0.88,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.06); 

  TString version   =  "rec_std" ;
  TString inputpath =  "files/byLLR/Apr23_2014/";

  //TTJetsSemiLept
  TFile* fS = 0;
  if(fname.find("SL")!=string::npos) 
    fS = TFile::Open(inputpath+"MEAnalysisNew_all_"+version+"_TTH125.root");
  else
    fS = TFile::Open(inputpath+"MEAnalysisNew_all_"+version+"_TTH125.root");
   
  TFile* fB = 0;
  if(fname.find("SL")!=string::npos) 
    fB = TFile::Open(inputpath+"MEAnalysisNew_all_"+version+"_TTJets.root");
  else
    fB = TFile::Open(inputpath+"MEAnalysisNew_all_"+version+"_TTJets.root");


  TTree* tS = (TTree*)fS->Get("tree");
  TTree* tB = (TTree*)fB->Get("tree");

  TH1F* hS = new TH1F("hS","Simulation #sqrt{s}=8 TeV "+title+"; P_{S/B}; units",20, 0,1);
  hS->SetLineColor(kRed);
  hS->SetLineWidth(3);
  hS->SetFillStyle(3005);
  hS->SetFillColor(kRed);
  hS->Sumw2();
  //if( fname=="DL" ) hS->Rebin(2);

  TH1F* hB = new TH1F("hB","Simulation #sqrt{s}=8 TeV "+title+"; P_{S/B}; units",20, 0,1);
  hB->SetLineColor(kBlue);
  hB->SetLineWidth(3);
  hB->SetFillStyle(3004);
  hB->SetFillColor(kBlue);
  hB->Sumw2();
  //if( fname=="DL" ) hB->Rebin(2);

  tS->Draw(Form("p_125_all_s/(p_125_all_s+p_125_all_b*%f)>>hS", fact), "weight"*TCut((header+" && nSimBs>=2").c_str()));
  tB->Draw(Form("p_125_all_s/(p_125_all_s+p_125_all_b*%f)>>hB", fact), "weight"*TCut((header+" && nSimBs>=2").c_str()));

  cout << "hS=" << hS->Integral() << endl;
  cout << "hB=" << hB->Integral() << endl;

  hS->SetTitleSize(0.05,"X");
  hS->SetTitleSize(0.05,"Y");
  hS->SetTitleOffset(0.88,"X");
  hS->SetTitleOffset(0.75,"Y");
  
  hB->SetTitleSize(0.05,"X");
  hB->SetTitleSize(0.05,"Y");
  hB->SetTitleOffset(0.88,"X");
  hB->SetTitleOffset(0.75,"Y");

  hS->Scale( 1./hS->Integral());
  hB->Scale( 1./hB->Integral());

  if(hS->GetMaximum()>=hB->GetMaximum()){  
    hS->SetMaximum(0.4);
    hS->SetMinimum(0.);
    hS->Draw("HISTE");
    hB->Draw("HISTESAME");
  }
  else{
    hB->SetMaximum(0.4);
    hB->SetMinimum(0.);
    hB->Draw("HISTE");
    hS->Draw("HISTESAME");
  }

  leg->AddEntry(hS,"t#bar{t}H");
  leg->AddEntry(hB,"t#bar{t}+jets");
  leg->Draw();

  TH1F* hRatio = (TH1F*) hS->Clone("hRatio");
  hRatio->Reset();
  hRatio->SetFillColor(0);
  hRatio->Sumw2();
  hRatio->Divide(hS,hB,1,1,"B");
  hRatio->SetMinimum(0.5);
  hRatio->SetMaximum(0.4);

  TPaveText *pt1 = new TPaveText(0.61, 0.43, 0.71, 0.49,"brNDC");
  pt1->SetFillStyle(1001);
  pt1->SetBorderSize(0);
  pt1->SetFillColor(kWhite);
  pt1->SetTextSize(0.09); 
  if( fname=="SL" )
    pt1->AddText("single-lepton channel");
  else
    pt1->AddText("di-lepton channel");
  pt1->Draw();


  if(1)
    c1->SaveAs(("Plots/PAS/Plot_"+fname+"_CompMadWeight.pdf").c_str());

}

void plot_analysis_MWAll(){

  plot_analysis_MW("SL", "type==0 && syst==0",                 "", "SL", 0.020);
  plot_analysis_MW("DL", "type==6 && syst==0 && btag_LR>0.99", "", "DL", 0.012);

}



void plot2DwoBBB( TString fname="mlfitMEM_COMB_New_rec_std_sb",
		  TString hname = "covariance_fit_s",
		  TString title = ""){

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

  TCanvas *c1 = new TCanvas("c1","",5,30,int(650*1.),int(650*1.));
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  c1->SetRightMargin(-0.08);
  c1->SetLeftMargin(0.25);
  c1->SetBottomMargin(0.25);

  /*
  TPad *p1 = new TPad("p1","",0.155,0.155,1,1);
  p1->SetGrid(0,0);
  p1->SetFillStyle(4000);
  p1->SetFillColor(10);
  p1->SetTicky();
  p1->SetObjectStat(0);
  p1->Draw();
  p1->cd();
  */

  TFile* f= TFile::Open("datacards/PreApproval/"+fname+".root");

  TH2D* h2 = (TH2D*)f->Get( hname );
  float binWidth = h2->GetXaxis()->GetBinWidth(1);
  int nBins      =  h2->GetXaxis()->GetNbins();

  int counter = 0;
  for( int bX = 1; bX <= nBins ; bX++ ){
    string labelX( h2->GetXaxis()->GetBinLabel(bX) );
    if( labelX.find("bin")==string::npos ){
      counter++;
    }        
  }
  cout << "Nusiances = " << counter << endl;

  TH2D* h2_copy = new TH2D("h2_copy", title, counter, 0, counter, counter, 0, counter  );
  int nBins2 = h2_copy->GetXaxis()->GetNbins();

  counter = 0;
  h2_copy->Reset();
  for( int bX = 1; bX <= nBins ; bX++ ){
    string labelX( h2->GetXaxis()->GetBinLabel(bX) );
    if( labelX.find("bin")==string::npos ){
      h2_copy->GetXaxis()->SetBinLabel( counter+1,      labelX.c_str() );
      h2_copy->GetYaxis()->SetBinLabel( nBins2-counter, labelX.c_str() );
      counter++;
    }        
  }
  h2_copy->GetXaxis()->SetLabelSize(0.03);
  h2_copy->GetYaxis()->SetLabelSize(0.03);
  h2_copy->GetXaxis()->LabelsOption("v");

  for( int bX = 1; bX <= nBins ; bX++ ){
    string labelX( h2->GetXaxis()->GetBinLabel(bX) );
    if( labelX.find("bin")!=string::npos ) continue;

    for( int bY = 1; bY <= nBins ; bY++ ){
      string labelY( h2->GetYaxis()->GetBinLabel(bY) );
      if( labelY.find("bin")!=string::npos ) continue;      
      
      float bin = h2->GetBinContent( bX, bY);
      cout << "Target = (" << labelX << "," << labelY << ") --> " << bin  << endl;

      for( int bX2 = 1; bX2 <= nBins2 ; bX2++ ){
	for( int bY2 = 1; bY2 <= nBins2 ; bY2++ ){
	  string labelX2( h2_copy->GetXaxis()->GetBinLabel(bX2) );
	  string labelY2( h2_copy->GetYaxis()->GetBinLabel(bY2) );
	  if( labelX2==labelX && labelY2==labelY ) 
	    h2_copy->SetBinContent( bX2, bY2, bin);
	}
      }
      
    }
  }
  
  h2_copy->SetTitle( title );
  h2_copy->Draw("COLZ");

  if(1){
    c1->SaveAs("datacards/PreApproval/Limits_Correlation_"+fname+"_"+hname+".pdf");
  }


}


void plot2DwoBBBAll(){
  plot2DwoBBB( "mlfitMEM_COMB_New_rec_std_sb", "covariance_fit_s", "Combined fit (s+b)"); 
  //plot2DwoBBB( "mlfitMEM_SL_New_rec_std_sb",   "covariance_fit_s", "SL fit (s+b)"); 
  //plot2DwoBBB( "mlfitMEM_DL_New_rec_std_sb",   "covariance_fit_s", "DL fit (s+b)"); 
  //plot2DwoBBB( "mlfitMEM_cat1_H_New_rec_std_sb", "covariance_fit_s", "SL Cat.1 (H) (s+b)"); 
  //plot2DwoBBB( "mlfitMEM_cat1_L_New_rec_std_sb", "covariance_fit_s", "SL Cat.1 (L) (s+b)"); 
  //plot2DwoBBB( "mlfitMEM_cat2_H_New_rec_std_sb", "covariance_fit_s", "SL Cat.2 (H) (s+b)"); 
  //plot2DwoBBB( "mlfitMEM_cat2_L_New_rec_std_sb", "covariance_fit_s", "SL Cat.2 (L) (s+b)"); 
  //plot2DwoBBB( "mlfitMEM_cat3_H_New_rec_std_sb", "covariance_fit_s", "SL Cat.3 (H) (s+b)"); 
  //plot2DwoBBB( "mlfitMEM_cat3_L_New_rec_std_sb", "covariance_fit_s", "SL Cat.3 (L) (s+b)"); 
  //plot2DwoBBB( "mlfitMEM_cat6_H_New_rec_std_sb", "covariance_fit_s", "DL (H) (s+b)"); 
  //plot2DwoBBB( "mlfitMEM_cat6_L_New_rec_std_sb", "covariance_fit_s", "DL (L) (s+b)"); 
  //plot2DwoBBB( "mlfitMEM_cat1_New_rec_std_sb", "covariance_fit_s", "SL Cat. 1 (s+b)"); 
  //plot2DwoBBB( "mlfitMEM_cat2_New_rec_std_sb", "covariance_fit_s", "SL Cat. 2 (s+b)"); 
  //plot2DwoBBB( "mlfitMEM_cat3_New_rec_std_sb", "covariance_fit_s", "SL Cat. 3 (s+b)"); 
  //plot2DwoBBB( "mlfitMEM_cat6_New_rec_std_sb", "covariance_fit_s", "DL (s+b)"); 
}
