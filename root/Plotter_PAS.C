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

#define RUNONDATA 0

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
  
  plot_param_3(1, "b-jets: quark #rightarrow jet",   "TEST");
  plot_param_3(1, "b-jets: gen-jet #rightarrow jet", "TEST_std_gen");

  plot_param_3(0, "udcs-jets: quark #rightarrow jet",   "TEST");
  plot_param_3(0, "udcs-jets: gen-jet #rightarrow jet", "TEST_std_gen");
    

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

  int nmax = 4000;
  if(flag>=8) nmax = 4000;

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
  if(flag<4) 
    leg->SetHeader("SL, N_{jet}#geq6");
  else if(flag>=4 && flag<8)   
    leg->SetHeader("SL, N_{jet}=5");
  else if(flag>=8 && flag<12)   
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

  if(oldFiles==0){
  
    if(flag<8){
      files.push_back(path+"MEAnalysisNew_CMVA_rec_std_TTH125.root");
      files.push_back(path+"MEAnalysisNew_CMVA_rec_std_TTJetsSemiLept.root");
    }
    else{
      files.push_back(path+"MEAnalysisNew_CMVA_rec_std_TTH125.root");
      files.push_back(path+"MEAnalysisNew_CMVA_rec_std_TTJetsFullLept.root");
    }
  
  }
  else{

    if(flag<8){
      files.push_back(path+"MEAnalysisNew_CSV_rec_std_TTH125.root");
      files.push_back(path+"MEAnalysisNew_CSV_rec_std_TTJetsSemiLept.root");
    }
    else{
      files.push_back(path+"MEAnalysisNew_CSV_rec_std_TTH125.root");
      files.push_back(path+"MEAnalysisNew_CSV_rec_std_TTJetsFullLept.root");
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
    int nMatchSimBs;
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

    int counter2M = 0;
    int counter3M = 0;
    int counter4M = 0;

    Long64_t nentries = t->GetEntries(); 
    int counter = 0;

    for (Long64_t i = 0; i < nentries ; i++){

      t->GetEntry(i);

      //if( oldFiles==0 ){
      if( flag<4             && type!=-1) continue;
      if( flag>=4 && flag<8  && type!=-2) continue;
      if( flag>=8 && flag<12 && type!=-3) continue;
      //}
      
      
      if(j==1 && flag%4==0     && nSimBs==0) continue; // tt+X
      if(j==1 && (flag-1)%4==0 && nSimBs >2) continue; // tt+LF
      if(j==1 && (flag-2)%4==0 && (nSimBs==2 || (nSimBs >2 && nMatchSimBs>=2) )) continue; // tt+b
      if(j==1 && (flag-3)%4==0 && (nSimBs==2 || (nSimBs >2 && nMatchSimBs<2) ))  continue; // tt+bb

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

    if( flag<=4 &&  hS->GetBinCenter(b)>0.959 &&   hS->GetBinCenter(b)<0.961 && xWP[0]<0 ){
      xWP[0] = intS;
      yWP[0] = intB;
      cout << "\e[1;32m" << intS << " --- " << intB << "\e[0m" << endl;
    }
    if( flag>4 && flag<=8 &&  hS->GetBinCenter(b)>0.969 &&   hS->GetBinCenter(b)<0.971 && xWP[0]<0 ){
      xWP[0] = intS;
      yWP[0] = intB;
      cout << "\e[1;32m" << intS << " --- " << intB << "\e[0m" << endl;
    }
    if(flag>8 && flag<=11  &&  hS->GetBinCenter(b)>0.859 &&   hS->GetBinCenter(b)<0.861 && xWP[0]<0 ){
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
  if( flag%4==0    ) title = "Simulation #sqrt{s}=8 TeV; selection efficiency in t#bar{t}H; selection efficiency in t#bar{t}+jets";
  if( (flag-1)%4==0) title = "Simulation #sqrt{s}=8 TeV; selection efficiency in t#bar{t}H; selection efficiency in t#bar{t}+jj";
  if( (flag-2)%4==0) title = "Simulation #sqrt{s}=8 TeV; selection efficiency in t#bar{t}H; selection efficiency in t#bar{t}+b ";
  if( (flag-3)%4==0) title = "Simulation #sqrt{s}=8 TeV; selection efficiency in t#bar{t}H; selection efficiency in t#bar{t}+b#bar{b}";

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
    h->SetMinimum(0.001);
    h->SetMaximum(0.20);
  }
  if(flag==5){
    h->SetMinimum(0.0002);
    h->SetMaximum(0.20);
  }
 if(flag==6){
    h->SetMinimum(0.001);
    h->SetMaximum(0.30);
  }
 if(flag==7){
    h->SetMinimum(0.01);
    h->SetMaximum(0.40);
  }

 if(flag==8){
    h->SetMinimum(0.001);
    h->SetMaximum(0.05);
  }
 if(flag==9){
    h->SetMinimum(0.0001);
    h->SetMaximum(0.05);
  }
 if(flag==10){
    h->SetMinimum(0.001);
    h->SetMaximum(0.20);
  }
 if(flag==11){
   h->SetMinimum(0.04);
   h->SetMaximum(0.40);
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
    if(flag==0 || flag==4 || flag==8)
      leg->AddEntry(hB, "t#bar{t}+jets","F");
    if(flag==1 || flag==5 || flag==9)
      leg->AddEntry(hB, "t#bar{t}+jj","F");
    if(flag==2 || flag==6 || flag==10)
      leg->AddEntry(hB, "t#bar{t}+b","F");
    if(flag==3 || flag==7 || flag==11)
      leg->AddEntry(hB, "t#bar{t}+b#bar{b}","F");

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
 
  float pos1[4] = {0.,0.,0.,0.};
  float pos2[4] = {0.,0.,0.,0.};
  float pos3[4] = {0.,0.,0.,0.};
  if( flag==0 ){
    pos1[0] = 0.61;  
    pos1[1] = 0.70;  
    pos1[2] = 0.83;  
    pos1[3] = 0.78;

    pos2[0] = 0.10;  
    pos2[1] = 0.31;  
    pos2[2] = 0.41;  
    pos2[3] = 0.43;

    pos3[0] = 0.31;  
    pos3[1] = 0.40;  
    pos3[2] = 0.48;  
    pos3[3] = 0.47;
  }
  if( flag==1 ){
    pos1[0] = 0.61;  
    pos1[1] = 0.70;  
    pos1[2] = 0.83;  
    pos1[3] = 0.78;

    pos2[0] = 0.10;  
    pos2[1] = 0.31;  
    pos2[2] = 0.41;  
    pos2[3] = 0.43;

    pos3[0] = 0.31;  
    pos3[1] = 0.40;  
    pos3[2] = 0.48;  
    pos3[3] = 0.47;
  }
  if( flag==2 ){
    pos1[0] = 0.61;  
    pos1[1] = 0.70;  
    pos1[2] = 0.83;  
    pos1[3] = 0.78;

    pos2[0] = 0.10;  
    pos2[1] = 0.31;  
    pos2[2] = 0.41;  
    pos2[3] = 0.43;

    pos3[0] = 0.31;  
    pos3[1] = 0.40;  
    pos3[2] = 0.48;  
    pos3[3] = 0.47;
  }
  if( flag==3 ){
    pos1[0] = 0.61;  
    pos1[1] = 0.70;  
    pos1[2] = 0.83;  
    pos1[3] = 0.78;

    pos2[0] = 0.10;  
    pos2[1] = 0.31;  
    pos2[2] = 0.41;  
    pos2[3] = 0.43;

    pos3[0] = 0.31;  
    pos3[1] = 0.40;  
    pos3[2] = 0.48;  
    pos3[3] = 0.47;
  }


  if( flag==4 ){
    pos1[0] = 0.49;  
    pos1[1] = 0.55;  
    pos1[2] = 0.71;  
    pos1[3] = 0.63;

    pos2[0] = 0.32;  
    pos2[1] = 0.09;  
    pos2[2] = 0.63;  
    pos2[3] = 0.20;

    pos3[0] = 0.19;  
    pos3[1] = 0.22;  
    pos3[2] = 0.36;  
    pos3[3] = 0.29;
  }
  if( flag==5 ){
    pos1[0] = 0.50;  
    pos1[1] = 0.64;  
    pos1[2] = 0.73;  
    pos1[3] = 0.72;

    pos2[0] = 0.33;  
    pos2[1] = 0.18;  
    pos2[2] = 0.64;  
    pos2[3] = 0.29;

    pos3[0] = 0.21;  
    pos3[1] = 0.34;  
    pos3[2] = 0.38;  
    pos3[3] = 0.41;
  }
  if( flag==6 ){
    pos1[0] = 0.50;  
    pos1[1] = 0.64;  
    pos1[2] = 0.73;  
    pos1[3] = 0.72;

    pos2[0] = 0.33;  
    pos2[1] = 0.18;  
    pos2[2] = 0.64;  
    pos2[3] = 0.29;

    pos3[0] = 0.21;  
    pos3[1] = 0.34;  
    pos3[2] = 0.38;  
    pos3[3] = 0.41;
  }
  if( flag==7 ){
    pos1[0] = 0.50;  
    pos1[1] = 0.64;  
    pos1[2] = 0.73;  
    pos1[3] = 0.72;

    pos2[0] = 0.33;  
    pos2[1] = 0.18;  
    pos2[2] = 0.64;  
    pos2[3] = 0.29;

    pos3[0] = 0.21;  
    pos3[1] = 0.34;  
    pos3[2] = 0.38;  
    pos3[3] = 0.41;
  }

  if( flag==8 ){
    pos1[0] = 0.50;  
    pos1[1] = 0.64;  
    pos1[2] = 0.73;  
    pos1[3] = 0.72;

    pos2[0] = 0.33;  
    pos2[1] = 0.18;  
    pos2[2] = 0.64;  
    pos2[3] = 0.29;

    pos3[0] = 0.21;  
    pos3[1] = 0.34;  
    pos3[2] = 0.38;  
    pos3[3] = 0.41;
  }
  if( flag==9 ){
    pos1[0] = 0.50;  
    pos1[1] = 0.64;  
    pos1[2] = 0.73;  
    pos1[3] = 0.72;

    pos2[0] = 0.33;  
    pos2[1] = 0.18;  
    pos2[2] = 0.64;  
    pos2[3] = 0.29;

    pos3[0] = 0.21;  
    pos3[1] = 0.34;  
    pos3[2] = 0.38;  
    pos3[3] = 0.41;
  }
 if( flag==10 ){
    pos1[0] = 0.37;  
    pos1[1] = 0.59;  
    pos1[2] = 0.60;  
    pos1[3] = 0.67;

    pos2[0] = 0.33;  
    pos2[1] = 0.18;  
    pos2[2] = 0.64;  
    pos2[3] = 0.29;

    pos3[0] = 0.21;  
    pos3[1] = 0.34;  
    pos3[2] = 0.38;  
    pos3[3] = 0.41;
  }
 if( flag==11 ){
    pos1[0] = 0.37;  
    pos1[1] = 0.59;  
    pos1[2] = 0.60;  
    pos1[3] = 0.67;

    pos2[0] = 0.33;  
    pos2[1] = 0.18;  
    pos2[2] = 0.64;  
    pos2[3] = 0.29;

    pos3[0] = 0.21;  
    pos3[1] = 0.34;  
    pos3[2] = 0.38;  
    pos3[3] = 0.41;
  }


  TPaveText *pt = new TPaveText(pos1[0], pos1[1], pos1[2], pos1[3],"brNDC");
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetFillColor(10);
  pt->SetTextSize(0.04);
  pt->SetTextAlign(11);
  pt->AddText("N_{CSVL}#geq4, N_{CSVM}#geq3")->SetTextColor(kBlue);
  //if(plotDistr==0) pt->Draw();

  TPaveText *pt2 = new TPaveText(pos2[0], pos2[1], pos2[2], pos2[3],"brNDC");
  pt2->SetFillStyle(0);
  pt2->SetBorderSize(0);
  pt2->SetFillColor(10);
  pt2->SetTextSize(0.04);
  pt2->SetTextAlign(11);
  pt2->AddText("N_{CSVM}#geq4, N_{CSVT}#geq2")->SetTextColor(kBlue);
  //if(plotDistr==0) pt2->Draw();

  TPaveText *pt3 = new TPaveText(pos3[0], pos3[1], pos3[2], pos3[3],"brNDC");
  pt3->SetFillStyle(0);
  pt3->SetBorderSize(0);
  pt3->SetFillColor(10);
  pt3->SetTextSize(0.04);
  pt3->SetTextAlign(11);
  pt3->AddText("N_{CSVM}#geq4")->SetTextColor(kBlue);
  //if(plotDistr==0) pt3->Draw();

  if(1){
    c1->SaveAs("Plots/PAS/Plot_BTag_ROC_"+extraname+".pdf");
  }

}

void plot_BTagAll(){
  
  for(int a = 0; a<12; a++){
    for(int b = 0; b<2; b++){
      plot_BTag(TString(Form("%d_%d",a,b)),a,b);
    }
  }
  return;
}



void plot_category(string dir  = "Apr07_2014",
		   string type = "MEM", 
		   string cat  = "cat1",
		   string tag = "New_rec_std_sb",
		   string header = "Cat 1",
		   string fname  = "",
		   int    drawseparation = 1,
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

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TLegend* leg = new TLegend(0.52,0.50,0.77,0.88,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04); 

  THStack* aStack = new THStack("aStack","CMS preliminary #sqrt{s}=8 TeV, L=19.04 fb^{-1}; P_{s/b} ; events ");

  vector<string> samples;
  samples.push_back("SingleT");
  samples.push_back("TTV");
  samples.push_back("EWK");
  samples.push_back("DiBoson");
  samples.push_back("TTJetsHFbb");
  samples.push_back("TTJetsHFb");
  samples.push_back("TTJetsLF");
  samples.push_back("TTH125");
  if(RUNONDATA)  samples.push_back("data_obs");


  TH1F* hS    = 0;
  TH1F* hData = 0;
  TH1F* hErr  = 0;

  TFile* f = TFile::Open(("datacards/"+dir+"/"+type+"_"+tag+".root").c_str());
  if(f==0 || f->IsZombie() ) return;


  ////////////////////////////////////////////////////////////////////////

  TH1F* h_csvUp                 = 0;
  TH1F* h_JECUp                 = 0;
  TH1F* h_JERUp                 = 0;
  TH1F* h_lumiUp                = 0;
  TH1F* h_pdf_ggUp              = 0;
  TH1F* h_Norm_TTVUp            = 0;
  TH1F* h_Norm_SingleTUp        = 0;
  TH1F* h_QCDscale_TTHUp        = 0;
  TH1F* h_QCDscale_TTbbUp       = 0;
  TH1F* h_QCDscale_TTbUp        = 0;
  TH1F* h_QCDscale_TTJetsHFUp   = 0;
  TH1F* h_QCDscale_TTJetsLFUp   = 0;

  map<string,TH1F*> h_Ups;


  f->cd( (type+"_"+cat).c_str());
 
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
    
    if( hname.find("TTH125")==string::npos && 
	hname.find("data")==string::npos && 
	hname.find("bin")==string::npos && 
	hname.find("Down")==string::npos &&
	hname.find("Up")!=string::npos ){

      string sysname = hname;
      sysname.erase(0, (sysname.find_first_of("_"))+1 );
      cout << sysname << endl;

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


    if( samples[sample].find("TTH125") != string::npos ){
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
    if( samples[sample].find("TTJetsHFbb") != string::npos ){
      leg->AddEntry(h, "t#bar{t} + bb", "F");
      //h->Scale( rescaleTTbb );
      h->SetLineColor( 16 );
      h->SetFillColor( 16 );
    }
    else if( samples[sample].find("TTJetsHFb") != string::npos ){
      leg->AddEntry(h, "t#bar{t} + b", "F");
      //h->Scale( rescaleTTb );
      h->SetLineColor( 17 );
      h->SetFillColor( 17 );
    }
    if( samples[sample].find("TTJetsLF") != string::npos ){
      leg->AddEntry(h, "t#bar{t} + jj", "F");
      //h->Scale( rescaleTTjj );
      h->SetLineColor( 18 );
      h->SetFillColor( 18 );
    }    
    if( samples[sample].find("SingleT") != string::npos ){
      h->SetLineColor( kMagenta );
      h->SetFillColor( kMagenta );
      leg->AddEntry(h, "Single top", "F");
    }
    if( samples[sample].find("EWK") != string::npos ){
      h->SetLineColor( kGreen );
      h->SetFillColor( kGreen );
      leg->AddEntry(h, "V+jets", "F");
    }
    if( samples[sample].find("DiBoson") != string::npos ){
      h->SetLineColor( kYellow );
      h->SetFillColor( kYellow );
      leg->AddEntry(h, "VV", "F");
    }
    if( samples[sample].find("TTV") != string::npos ){
      h->SetLineColor( 30 );
      h->SetFillColor( 30 );
      leg->AddEntry(h, "t#bar{t}V", "F");
    }
    if( samples[sample].find("data_obs") != string::npos && RUNONDATA){
      hData = (TH1F*)h->Clone("hData");
      //cout <<  hData->GetNbinsX() << endl;
      hData->Sumw2();
      hData->SetMarkerStyle(kFullCircle);
      hData->SetMarkerSize(1.5);
      hData->SetMarkerColor(kBlack);
      //hData->SetBinContent( hData->GetNbinsX()   , 0.);
      //hData->SetBinContent( hData->GetNbinsX()-1 , 0.);
      //hData->SetBinContent( hData->GetNbinsX()-2 , 0.);
      //hData->SetBinError( hData->GetNbinsX()   , 0.);
      //hData->SetBinError( hData->GetNbinsX()-1 , 0.);
      //hData->SetBinError( hData->GetNbinsX()-2 , 0.);
      leg->AddEntry(hData, "data", "P");
    }

    if( !(samples[sample].find("TTH125")!=string::npos || samples[sample].find("data_obs")!=string::npos) ){
      if( hErr==0 ){
	hErr = (TH1F*)h->Clone("hErr");
	hErr->Reset();
	hErr->SetFillStyle(0);
	hErr->SetLineColor(kBlack);
	//leg->AddEntry(hErr, "MC unc. (stat.)", "L");
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
    hErr->SetTitle("Simulation #sqrt{s}=8 TeV, L=19.6 fb^{-1}");
  else
    hErr->SetTitle("CMS Preliminary #sqrt{s}=8 TeV, L=19.6 fb^{-1}");
  hErr->SetTitleSize  (0.04,"X");
  hErr->SetTitleOffset(0.95,"X");
  float max =  TMath::Max(float(hErr->GetMaximum()*1.55),(hData!=0 ? hData->GetMaximum()*1.55 : -1.));
  hErr->GetYaxis()->SetRangeUser(0., max );
  hErr->SetLineColor(kBlack);

  hErr->Draw("HIST");
  TH1D* hErr_clone = (TH1D*)hErr->Clone("hErr_clone");
  hErr_clone  ->SetFillColor(kBlack);
  hErr_clone  ->SetFillStyle(3004);
  hErr_clone  ->SetLineWidth(1);
  hErr_clone  ->Draw("E2SAME");
  leg->AddEntry(hErr_clone, "Bkg unc.", "F");

  //hErr->Draw("HISTE1");
  aStack->Draw("HISTSAME");
  leg->SetHeader( header.c_str() );
  leg->AddEntry(hS, "signal x 10", "L");
  hS->Draw("HISTSAME");
  hErr_clone  ->Draw("E2SAME");


  if(hData!=0) hData->Draw("ESAME");
  leg->Draw();

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

  if(1){
    c1->SaveAs(  ("Plots/PAS/Plot_"+cat+"_"+fname+"_"+tag+".pdf").c_str() );
  }

  cout << "REMOVE: " << endl;
  for( map<string,TH1F*>::iterator it = h_Ups.begin() ; it != h_Ups.end() ; it++){
    if( (it->second) )  delete (it->second);
  }

}


void plotAll(){

  plot_category("Apr07_2014","MEM", "cat1_H", "New_rec_std_sb", "SL Cat-1 (H)", "",   1);
  plot_category("Apr07_2014","MEM", "cat2_H", "New_rec_std_sb", "SL Cat-2 (H)", "",   1);
  plot_category("Apr07_2014","MEM", "cat3_H", "New_rec_std_sb", "SL Cat-3 (H)", "",   1);
  plot_category("Apr07_2014","MEM", "cat6_H", "New_rec_std_sb", "DL (H)", "", 1);

  plot_category("Apr07_2014","MEM", "cat1_L", "New_rec_std_sb", "SL Cat-1 (L)", "",   0);
  plot_category("Apr07_2014","MEM", "cat2_L", "New_rec_std_sb", "SL Cat-2 (L)", "",   0);
  plot_category("Apr07_2014","MEM", "cat3_L", "New_rec_std_sb", "SL Cat-3 (L)", "",   0);
  plot_category("Apr07_2014","MEM", "cat6_L", "New_rec_std_sb", "DL (L)", "", 0);

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

  delete c1;
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
    cout << hnames[k] << endl;
    if( hnames[k].find("resp")!=string::npos )
      plot_param( TString(titles[k].c_str()) ,"quark energy (GeV)","GeV",TString(hnames[k]),  50 , 200, 50, 200, "parametrization", ""); 
    if( hnames[k].find("resol")!=string::npos )
      plot_param( TString(titles[k].c_str()) ,"quark energy (GeV)","GeV",TString(hnames[k]),  50 , 200, 0, 100, "parametrization", ""); 

    cout << endl;
  }


}


