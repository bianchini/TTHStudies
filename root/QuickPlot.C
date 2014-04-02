#include <cstdlib>
#include <iostream> 
#include <map>
#include <string>
#include <vector>

#include "TMath.h"
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
#include "TList.h"
#include "TCollection.h"
#include "TObject.h"
#include "TLegend.h"
#include "THStack.h"
#include "TCut.h"
#include "TLine.h"

#include "TVector3.h"
#include "TLorentzVector.h"






void plot_analysis( string header="syst==0 && type==0",
		    float fact = 1
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

  TLegend* leg = new TLegend(0.25,0.65,0.65,0.88,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.06); 

  TFile* fS = TFile::Open("files/byLLR/CSVcalibration_V3/MEAnalysisNew_all_CSVcalibration_rec_std_TTJets.root");
  TFile* fB = TFile::Open("MEAnalysisNew_test_analyticalME_v2_rec_std_TTJets.root");

  TTree* tS = (TTree*)fS->Get("tree");
  TTree* tB = (TTree*)fB->Get("tree");

  TH1F* hS = new TH1F("hS","Simulation #sqrt{s}=8 TeV; P_{S/B}(y); units",12, 0,1);
  hS->SetLineColor(kRed);
  hS->SetLineWidth(3);
  hS->SetFillStyle(3005);
  hS->SetFillColor(kRed);
  hS->Sumw2();
  TH1F* hB = new TH1F("hB","Simulation #sqrt{s}=8 TeV,; P_{S/B}(y); units",12, 0,1);
  hB->SetLineColor(kBlue);
  hB->SetLineWidth(3);
  hB->SetFillStyle(3004);
  hB->SetFillColor(kBlue);
  hB->Sumw2();

  tS->Draw(Form("p_125_all_s/(p_125_all_s+p_125_all_b*%f)>>hS", fact), "weight"*TCut((header+" && nSimBs>=2").c_str()));
  tB->Draw(Form("p_125_all_s/(p_125_all_s+p_125_all_b*%f)>>hB", fact), "weight"*TCut((header+" && nSimBs>=2").c_str()));

  cout << "hS=" << hS->GetEntries() << endl;
  cout << "hB=" << hB->GetEntries() << endl;

  hS->SetTitleSize(0.05,"X");
  hS->SetTitleSize(0.05,"Y");
  hS->SetTitleOffset(0.88,"X");
  hS->SetTitleOffset(0.75,"Y");
  
  hB->SetTitleSize(0.05,"X");
  hB->SetTitleSize(0.05,"Y");
  hB->SetTitleOffset(0.88,"X");
  hB->SetTitleOffset(0.75,"Y");

  if(hS->GetMaximum()>=hB->GetMaximum()){
    hS->Scale( 1./hS->Integral());
    hB->Scale( 1./hB->Integral());
    hS->SetMaximum(1);
    hS->SetMinimum(0.);
    hS->Draw("HISTE");
    hB->Draw("HISTESAME");
  }
  else{
    hS->Scale( 1./hS->Integral());
    hB->Scale( 1./hB->Integral());
    hB->SetMaximum(1);
    hB->SetMinimum(0.);
    hB->Draw("HISTE");
    hS->Draw("HISTESAME");
  }

  //leg->SetHeader( Form("#sigma(S)/#sigma(B)=%.2f 10^{-2}",1./fact ) );
  leg->AddEntry(hS,"OLD");
  leg->AddEntry(hB,"NEW");
  leg->Draw();

  TH1F* hRatio = (TH1F*) hS->Clone("hRatio");
  hRatio->Reset();
  hRatio->SetFillColor(0);
  hRatio->Sumw2();
  hRatio->Divide(hS,hB,1,1,"B");
  hRatio->SetMinimum(0.5);
  hRatio->SetMaximum(0.4);
  //hRatio->Draw("HISTE");
  //c1->SaveAs("SoB_"+cat+"_"+fname+".png");


  //c1->SaveAs("ttbb_vs_ttjj.png");


}

void plot_analysis_ttbb(TString cat   ="SL", 
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


  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
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

  //TFile* fB = TFile::Open("MEAnalysis_"+cat+"_VType2_nominal_v5_TTH125.root");

  TFile* fS = TFile::Open("MEAnalysis_"+cat+"_nominal_v7_TTJets.root");
  TFile* fB = TFile::Open("MEAnalysis_"+cat+"_nominal_v7_TTJets.root");

  //TFile* fS = TFile::Open("MEAnalysis_"+cat+"_VType2_nominal_nSvs_TTH125.root");
  //TFile* fS = TFile::Open("MEAnalysis_"+cat+"_VType2_nominal_nSvs_TTJetsSemiLept.root");
  //TFile* fS = TFile::Open("MEAnalysis_"+cat+"_nominal_v6_TTH125.root");

  //TFile* fS = TFile::Open("MEAnalysis_"+cat+"_VType2_nominal_nSvs_TTJetsSemiLept.root");
  //TFile* fB = TFile::Open("MEAnalysis_"+cat+"_VType2_nominal_nSvs_TTJetsSemiLept.root");
  //TFile* fS = TFile::Open("MEAnalysis_"+cat+"_VType2_nominal_v5_TTJetsSemiLept.root");
  
  //TFile* fS = TFile::Open("MEAnalysis_"+cat+"_nominal_v6_TTJetsSemiLept.root");
  //TFile* fB = TFile::Open("MEAnalysis_"+cat+"_nominal_v6_TTJetsFullLept.root");



  TTree* tS = (TTree*)fS->Get("tree");
  TTree* tB = (TTree*)fB->Get("tree");

  TH1F* hS = new TH1F("hS","Simulation #sqrt{s}=8 TeV, "+title+"; P_{S/B}(y); units",5, 0,1);
  hS->SetLineColor(kRed);
  hS->SetLineWidth(3);
  hS->SetFillStyle(3005);
  hS->SetFillColor(kRed);
  hS->Sumw2();
  TH1F* hB = new TH1F("hB","Simulation #sqrt{s}=8 TeV, "+title+"; P_{S/B}(y); units",5, 0,1);
  hB->SetLineColor(kBlue);
  hB->SetLineWidth(3);
  hB->SetFillStyle(3004);
  hB->SetFillColor(kBlue);
  hB->Sumw2();

  //tS->Draw(Form("p_125_all_rec/2.5e+16/(p_125_all_rec/2.5e+16+p_125_all_alt_rec/(2.7e+18)/%f)>>hS", fact));
  //tB->Draw(Form("p_125_all_rec/2.5e+16/(p_125_all_rec/2.5e+16+p_125_all_alt_rec/(2.7e+18)/%f)>>hB", fact));

  //tS->Draw(Form("p_125_all_s/(p_125_all_s+p_125_all_b*%f)>>hS", fact), TCut((header+" && nSimBs>2 && nMatchSimBs>=2").c_str()));
  //tB->Draw(Form("p_125_all_s/(p_125_all_s+p_125_all_b*%f)>>hB", fact), TCut((header+" && nSimBs>2 && nMatchSimBs<2").c_str()));

  //tS->Draw(Form("p_125_all_b_ttbb/((p_125_all_b_ttbb+%f*p_125_all_b_ttjj))>>hS", fact), TCut((header+" && nSimBs>2  && p_125_all_s_ttbb/(p_125_all_s_ttbb+0.02*((1-0.5)*p_125_all_b_ttbb+0.5*p_125_all_b_ttjj))<0.1").c_str()));
  //tB->Draw(Form("p_125_all_b_ttbb/((p_125_all_b_ttbb+%f*p_125_all_b_ttjj))>>hB", fact), TCut((header+" && nSimBs==2 && p_125_all_s_ttbb/(p_125_all_s_ttbb+0.02*((1-0.5)*p_125_all_b_ttbb+0.5*p_125_all_b_ttjj))<0.1").c_str()));

  if(fname.find("L")!=string::npos){
    tS->Draw(Form("p_125_all_s/(p_125_all_s+p_125_all_b*%f)>>hS", fact), TCut((header+" && nSimBs>2  && nMatchSimBs>=-999").c_str()));
    tB->Draw(Form("p_125_all_s/(p_125_all_s+p_125_all_b*%f)>>hB", fact), TCut((header+" && nSimBs==2 && nMatchSimBs>=-999").c_str()));
  }
  else{
    tS->Draw(Form("p_125_all_s_ttbb/(p_125_all_s_ttbb+%f*((1-%f)*p_125_all_b_ttbb+%f*p_125_all_b_ttjj))>>hS", fact, fact2, fact2), TCut((header+" && nSimBs>2 && nMatchSimBs>=-999").c_str()) );
    tB->Draw(Form("p_125_all_s_ttbb/(p_125_all_s_ttbb+%f*((1-%f)*p_125_all_b_ttbb+%f*p_125_all_b_ttjj))>>hB", fact, fact2, fact2), TCut((header+" && nSimBs==2 && nMatchSimBs>=-999").c_str()) );
  }

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

  if(hS->GetMaximum()>=hB->GetMaximum()){
    hS->Scale( 1./hS->Integral());
    hB->Scale( 1./hB->Integral());
    hS->SetMaximum(1);
    hS->SetMinimum(0.);
    hS->Draw("HISTE");
    hB->Draw("HISTESAME");
  }
  else{
    hS->Scale( 1./hS->Integral());
    hB->Scale( 1./hB->Integral());
    hB->SetMaximum(1);
    hB->SetMinimum(0.);
    hB->Draw("HISTE");
    hS->Draw("HISTESAME");
  }

  //leg->SetHeader( Form("#sigma(S)/#sigma(B)=%.2f 10^{-2}",1./fact ) );
  leg->AddEntry(hS,"t#bar{t}bb");
  leg->AddEntry(hB,"t#bar{t}jj");
  leg->Draw();

  TH1F* hRatio = (TH1F*) hS->Clone("hRatio");
  hRatio->Reset();
  hRatio->SetFillColor(0);
  hRatio->Sumw2();
  hRatio->Divide(hS,hB,1,1,"B");
  hRatio->SetMinimum(0.5);
  hRatio->SetMaximum(0.4);
  //hRatio->Draw("HISTE");
  //c1->SaveAs("SoB_"+cat+"_"+fname+".png");

  if(1)
    c1->SaveAs(("Plot_"+fname+"_AN_ttbb_vs_ttb.pdf").c_str());


}


void plot_analysis_bj(TString cat   ="SL", 
		      string header="type==3",
		      TString title = "",
		      string fname= "type0",
		      string selS = "nSimBs>2  && nMatchSimBs>=-999",
		      string selB = "nSimBs==2  && nMatchSimBs>=-999"
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

  TLegend* leg = new TLegend(0.25,0.65,0.65,0.88,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.06); 


  TFile* fS = TFile::Open("MEAnalysis_"+cat+"_nominal_v7_TTJets.root");
  TFile* fB = TFile::Open("MEAnalysis_"+cat+"_nominal_v7_TTJets.root");
  //TFile* fB = TFile::Open("MEAnalysis_"+cat+"_nominal_v6_TTH125.root");

  TTree* tS = (TTree*)fS->Get("tree");
  TTree* tB = (TTree*)fB->Get("tree");

  TH1F* hS = new TH1F("hS","Simulation #sqrt{s}=8 TeV, "+title+"; P_{b/j}(y); units",10, 0,1);
  hS->SetLineColor(kRed);
  hS->SetLineWidth(3);
  hS->SetFillStyle(3005);
  hS->SetFillColor(kRed);
  hS->Sumw2();
  TH1F* hB = new TH1F("hB","Simulation #sqrt{s}=8 TeV, "+title+"; P_{b/j}(y); units",10, 0,1);
  hB->SetLineColor(kBlue);
  hB->SetLineWidth(3);
  hB->SetFillStyle(3004);
  hB->SetFillColor(kBlue);
  hB->Sumw2();

  tS->Draw("p_125_all_b_ttbb/(p_125_all_b_ttbb+0.5*p_125_all_b_ttjj)>>hS", TCut((header+" && "+selS).c_str()));
  tB->Draw("p_125_all_b_ttbb/(p_125_all_b_ttbb+0.5*p_125_all_b_ttjj)>>hB", TCut((header+" && "+selB).c_str()));


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

  if(hS->GetMaximum()>=hB->GetMaximum()){
    hS->Scale( 1./hS->Integral());
    hB->Scale( 1./hB->Integral());
    hS->SetMaximum(1);
    hS->SetMinimum(0.);
    hS->Draw("HISTE");
    hB->Draw("HISTESAME");
  }
  else{
    hS->Scale( 1./hS->Integral());
    hB->Scale( 1./hB->Integral());
    hB->SetMaximum(1);
    hB->SetMinimum(0.);
    hB->Draw("HISTE");
    hS->Draw("HISTESAME");
  }

  //leg->SetHeader( Form("#sigma(S)/#sigma(B)=%.2f 10^{-2}",1./fact ) );
  leg->AddEntry(hS,"t#bar{t}bb");
  if(fname.find("ttjj")!=string::npos)  
    leg->AddEntry(hB,"t#bar{t}jj");
  else
    leg->AddEntry(hB,"t#bar{t}b");
  leg->Draw();


  if(1)
    c1->SaveAs(("Plot_"+fname+".pdf").c_str());

  return;

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


  //TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
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

  TString version   =  "v4_gen_std" ;
  TString inputpath = "gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/Trees/MEM/Dec06_2013/WS/";

  //TTJetsSemiLept
  TFile* fS = 0;
  if(fname.find("SL")!=string::npos) 
    fS = TFile::Open(inputpath+"MEAnalysisNew_SL_nominal_"+version+"_TTH125.root");
  else
    fS = TFile::Open(inputpath+"MEAnalysisNew_DL_nominal_"+version+"_TTH125.root");
   
  TFile* fB = 0;
  if(fname.find("SL")!=string::npos) 
    fB = TFile::Open(inputpath+"MEAnalysisNew_SL_nominal_"+version+"_TTJets.root");
  else
    fB = TFile::Open(inputpath+"MEAnalysisNew_DL_nominal_"+version+"_TTJets.root");


  TTree* tS = (TTree*)fS->Get("tree");
  TTree* tB = (TTree*)fB->Get("tree");

  TH1F* hS = new TH1F("hS","Simulation #sqrt{s}=8 TeV "+title+"; P_{S/B}(y); units",20, 0,1);
  hS->SetLineColor(kRed);
  hS->SetLineWidth(3);
  hS->SetFillStyle(3005);
  hS->SetFillColor(kRed);
  hS->Sumw2();
  //if( fname=="DL" ) hS->Rebin(2);

  TH1F* hB = new TH1F("hB","Simulation #sqrt{s}=8 TeV "+title+"; P_{S/B}(y); units",20, 0,1);
  hB->SetLineColor(kBlue);
  hB->SetLineWidth(3);
  hB->SetFillStyle(3004);
  hB->SetFillColor(kBlue);
  hB->Sumw2();
  //if( fname=="DL" ) hB->Rebin(2);

  tS->Draw(Form("p_125_all_s/(p_125_all_s+p_125_all_b*%f)>>hS", fact), TCut((header+" && nSimBs>=2").c_str()));
  tB->Draw(Form("p_125_all_s/(p_125_all_s+p_125_all_b*%f)>>hB", fact), TCut((header+" && nSimBs>=2").c_str()));

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

  //leg->SetHeader( Form("#sigma(S)/#sigma(B)=%.2f 10^{-2}",1./fact ) );
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
  //hRatio->Draw("HISTE");
  //c1->SaveAs("SoB_"+cat+"_"+fname+".png");


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
    c1->SaveAs(("Plots/WS/Plot_"+fname+"_CompMadWeight.png").c_str());

  //c1->SaveAs("ttbb_vs_ttjj.png");


}



void plot_syst(TString cat   = "SL",
	       TString syst  = "csv", 
	       string header = "type==0",
	       TString fname = "type0",
	       TString title = "",
	       int doRatio = 0,
	       TString sample = "v6_TTJetsSemiLept",
	       float fact = 0.02, float fact2 = 0.5
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

  TLegend* leg = new TLegend(0.25,0.65,0.65,0.88,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04); 


  TFile* fN  = TFile::Open("MEAnalysis_"+cat+"_nominal_"+sample+".root");
  TFile* fU  = TFile::Open("MEAnalysis_"+cat+"_"+syst+"Up_"+sample+".root");
  TFile* fD  = TFile::Open("MEAnalysis_"+cat+"_"+syst+"Down_"+sample+".root");


  TTree* tN = (TTree*)fN->Get("tree");
  TTree* tU = (TTree*)fU->Get("tree");
  TTree* tD = (TTree*)fD->Get("tree");

  TH1F* hN = new TH1F("hN","Simulation #sqrt{s}=8 TeV, "+title+"; S/(S+B); units",4, 0,1);
  hN->SetLineColor(kBlack);
  hN->SetLineWidth(3);
  hN->Sumw2();
  TH1F* hU = new TH1F("hU","Simulation #sqrt{s}=8 TeV, "+title+"; S/(S+B); units",4, 0,1);
  hU->SetLineColor(kBlue);
  hU->SetLineWidth(3);
  hU->Sumw2();
  TH1F* hD = new TH1F("hD","Simulation #sqrt{s}=8 TeV, "+title+"; S/(S+B); units",4, 0,1);
  hD->SetLineColor(kRed);
  hD->SetLineWidth(3);
  hD->Sumw2();

  tN->Draw(Form("p_125_all_s_ttbb/(p_125_all_s_ttbb+%f*((1-%f)*p_125_all_b_ttbb+%f*p_125_all_b_ttjj))>>hN", fact, fact2, fact2), TCut((header+" && nSimBs>=2").c_str()));
  tU->Draw(Form("p_125_all_s_ttbb/(p_125_all_s_ttbb+%f*((1-%f)*p_125_all_b_ttbb+%f*p_125_all_b_ttjj))>>hU", fact, fact2, fact2), TCut((header+" && nSimBs>=2").c_str()));
  tD->Draw(Form("p_125_all_s_ttbb/(p_125_all_s_ttbb+%f*((1-%f)*p_125_all_b_ttbb+%f*p_125_all_b_ttjj))>>hD", fact, fact2, fact2), TCut((header+" && nSimBs>=2").c_str()));

  hN->SetMaximum(  hN->GetMaximum()*1.8);
  hN->SetMinimum( 0 );

  hN->Draw("HISTE");
  hU->Draw("HISTESAME");
  hD->Draw("HISTESAME");
 

  //leg->SetHeader( Form("#sigma(S)/#sigma(B)=%.2f 10^{-2}",1./fact ) );
  leg->AddEntry(hN,"nominal");
  leg->AddEntry(hU,syst+" up");
  leg->AddEntry(hD,syst+" down");
  leg->Draw();

  TH1F* hRatioN = (TH1F*) hN->Clone("hRatioN");
  hRatioN->Reset();
  hRatioN->SetFillColor(0);
  hRatioN->Sumw2();
  hRatioN->Divide(hN,hN,1,1,"B");
  TH1F* hRatioU = (TH1F*) hN->Clone("hRatioU");
  hRatioU->Reset();
  hRatioU->SetFillColor(0);
  hRatioU->Sumw2();
  hRatioU->Divide(hU,hN,1,1,"B");
  TH1F* hRatioD = (TH1F*) hN->Clone("hRatioD");
  hRatioD->Reset();
  hRatioD->SetFillColor(0);
  hRatioD->Sumw2();
  hRatioD->Divide(hD,hN,1,1,"B");

  if(doRatio){
    hRatioN->SetMinimum(0.4);
    hRatioN->SetMaximum(1.6);
    hRatioN->Draw("HISTE");
    hRatioU->Draw("HISTESAME");
    hRatioD->Draw("HISTESAME");
  }

  //c1->SaveAs("SoB_"+cat+"_"+fname+".png");

}

void plot_time(
	       TString sample = "TTH125",
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

  TLegend* leg = new TLegend(0.55,0.65,0.80,0.88,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04); 


  TFile* fN_SL  = TFile::Open("MEAnalysis_SL_nominal_v6_"+sample+".root");
  TTree* tN_SL = (TTree*)fN_SL->Get("tree");

  TFile* fN_DL  = string(sample.Data()).find("TTJets")!=string::npos ?  TFile::Open("MEAnalysis_DL_nominal_v5_TTJetsFullLept.root") : TFile::Open("MEAnalysis_DL_nominal_v5_"+sample+".root");
  TTree* tN_DL = (TTree*)fN_DL->Get("tree");

  TH1F* hN_type0 = new TH1F("hN_type0","Simulation #sqrt{s}=8 TeV "+title+"; time/event [sec]; units",50, 0,200);
  hN_type0->SetLineColor(kRed);
  hN_type0->SetLineWidth(3);
  hN_type0->SetFillColor(kRed);
  hN_type0->SetFillStyle(3004);
  hN_type0->SetTitleSize(0.05,"X");
  hN_type0->SetTitleSize(0.05,"Y");
  hN_type0->SetTitleOffset(0.90,"X");
  hN_type0->SetTitleOffset(0.75,"Y");
  tN_SL->Draw("time>>hN_type0","type==0 || type==3");

  TH1F* hN_type1 = new TH1F("hN_type1","Simulation #sqrt{s}=8 TeV, "+title+"; time/event [sec]; units",50, 0,200);
  hN_type1->SetLineColor(kBlue);
  hN_type1->SetLineWidth(3);
  hN_type1->SetFillColor(kBlue);
  hN_type1->SetFillStyle(3004);
  tN_SL->Draw("time>>hN_type1","type==1");

  TH1F* hN_type2 = new TH1F("hN_type2","Simulation #sqrt{s}=8 TeV, "+title+"; time/event [sec]; units",50, 0,200);
  hN_type2->SetLineColor(kGreen);
  hN_type2->SetLineWidth(3);
  hN_type2->SetFillColor(kGreen);
  hN_type2->SetFillStyle(3004);
  tN_SL->Draw("time>>hN_type2","type==2");

  TH1F* hN_type6 = new TH1F("hN_type6","Simulation #sqrt{s}=8 TeV, "+title+"; time/event [sec]; units",50, 0,200);
  hN_type6->SetLineColor(kMagenta);
  hN_type6->SetLineWidth(3);
  hN_type6->SetFillColor(kMagenta);
  hN_type6->SetFillStyle(3005);
  tN_DL->Draw("time>>hN_type6","type==6 || type==7");

  //hN->SetMaximum(  hN->GetMaximum()*1.8);
  hN_type0->SetMinimum( 0 );
  hN_type0->DrawNormalized("HIST");
  hN_type1->DrawNormalized("HISTSAME");
  hN_type2->DrawNormalized("HISTSAME");
  hN_type6->DrawNormalized("HISTSAME");

  leg->AddEntry(hN_type0,"Cat. 1,5", "F");
  leg->AddEntry(hN_type1,"Cat. 2",   "F");
  leg->AddEntry(hN_type2,"Cat. 3,4", "F");
  leg->AddEntry(hN_type6,"Cat. 6,7", "F");
  leg->Draw();

  if(1){
    c1->SaveAs("Plots/Time_"+sample+".png");
    c1->SaveAs("Plots/Time_"+sample+".pdf");
  }

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

  TFile *f = TFile::Open("../bin/root/ControlPlotsTEST_std_gen.root","READ");
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
  hBestMass->SetTitle("");
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
  if(func) leg->AddEntry(func, param, "L");

  if(func){
    //hBestMass->Draw("PE");
    func->Draw();
  }
  else{
    hBestMass->SetFillColor(kRed);
    hBestMass->SetFillStyle(3004);
    hBestMass->Draw("HIST");
  }

  leg->Draw();

  if(func) cout << "f(50,100,150)=" << string(Form("(%.1f, %.1f, %.1f)",func->Eval(50), func->Eval(100), func->Eval(150))) << endl;

  if(0){
    c1->SaveAs("Plots/Plot_"+hname+"_tf_"+extraLabel+".png");
    c1->SaveAs("Plots/Plot_"+hname+"_tf_"+extraLabel+".pdf");
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
  hnames.push_back("resolHeavyBin0");
  hnames.push_back("respHeavyBin0");
  hnames.push_back("resolHeavyBin1");
  hnames.push_back("respHeavyBin1");
  hnames.push_back("respG1HeavyBin0");
  hnames.push_back("respG2HeavyBin0");
  hnames.push_back("resolG1HeavyBin0");
  hnames.push_back("resolG2HeavyBin0");
  hnames.push_back("respG1HeavyBin1");
  hnames.push_back("respG2HeavyBin1");
  hnames.push_back("resolG1HeavyBin1");
  hnames.push_back("resolG2HeavyBin1");


  for(int k = 0 ; k < hnames.size(); k++){
    cout << hnames[k] << endl;
    plot_param( "","","",TString(hnames[k])); 
    cout << endl;
  }


}


void plot_param_3(int isB = 1, TString header = "Light jet TF"){

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
 

  TH1F* h = new TH1F("h"," ; parton energy (GeV) ; transfer function (1/GeV)",1,5,225);
  h->GetYaxis()->SetRangeUser( 0, 0.07);
  h->GetYaxis()->SetTitleSize(0.05);
  h->GetXaxis()->SetTitleSize(0.05);

  h->GetYaxis()->SetTitleOffset(0.95);
  h->GetXaxis()->SetTitleOffset(0.90);

  TFile *f = TFile::Open("../bin/root/ControlPlotsTEST.root","READ");


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
  leg->AddEntry(l_Bin0_1 , "50 GeV" ,  "L");
  leg->AddEntry(l_Bin0_2 , "100 GeV" , "L");
  leg->AddEntry(l_Bin0_3 , "150 GeV" , "L");

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

  if(0){
  }

  return;
}





void plot_Slice(
		TString header = "Light jet TF",
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
  hBestMass->SetTitle("");
  hBestMass->SetXTitle( unitsX );
  //hBestMass->SetYTitle( "units" );
  hBestMass->SetTitleSize(0.04,"X");
  hBestMass->SetTitleSize(0.04,"Y");

  TLegend* leg = new TLegend(0.11,0.71,0.45,0.94,NULL,"brNDC");
  leg->SetHeader( header );
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.06); 
  //leg->AddEntry(hBestMass, unitsY);

  float normFunc = 0.;  
  if(func) {
    leg->AddEntry(func, param, "L");
    func->SetNpx(1000);
    normFunc = func->GetParameter(0);
    cout << "Bin: " << func->Integral(-800,800) << endl;
  }
  if(funcP){
    leg->AddEntry(funcP,  "2G parametrization", "L");
    funcP->SetNpx(1000);
    funcP->SetParameter(0, normFunc);
    cout << "Param 2G: " << funcP->Integral(-800,800) << endl;
    if( TMath::Abs(funcP->Integral(-800,800)/func->Integral(-800,800)-1)>1e-02  )
      funcP->SetParameter(0, funcP->GetParameter(0)*func->Integral(-800,800)/funcP->Integral(-800,800)  );
    cout << "RESCALE: Param 2G: " << funcP->Integral(-800,800) << endl;
  }
  if(funcP2){
    leg->AddEntry(funcP2, "1G parametrization", "L");
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
  if( funcP)  funcP->Draw("SAME");
  if( funcP2) funcP2->Draw("SAME");
  leg->Draw();

  if(1){
    string eventStr = event<10 ? string(Form("0%d",event)) : string(Form("%d",event));
    c1->SaveAs(Form("Plots/Plot_%s_%s_%d%s.png", flavor.c_str(), eventStr.c_str(), bin, extraLabel.Data()));
    //c1->SaveAs(Form("Plots/Plot_%s_%d_%d.pdf", flavor.c_str(), event, bin));
  }


  return;
}


void plot_Slice_Reg(
		    TString header = "Light jet TF",
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

  TFile *f     = TFile::Open("../bin/root/ControlPlotsTEST.root",    "READ");
  TFile *f_reg = TFile::Open("../bin/root/ControlPlotsTEST_reg.root","READ");

  TH1F* hBestMass = (TH1F*)f->Get(TString(Form("Bin%s%d_%d/hCorr%s_py",flavor.c_str(), event, bin,flavor.c_str())));

  TF1* func     = 0;
  TF1* func_reg = 0;
  if( flavor.find("Heavy")!=string::npos ){
    func     = (TF1*)f    ->Get(TString(Form("Bin%s%d_%d/resolDG_fit%d_%d",flavor.c_str(), event, bin, event, bin)));
    func_reg = (TF1*)f_reg->Get(TString(Form("Bin%s%d_%d/resolDG_fit%d_%d",flavor.c_str(), event, bin, event, bin)));
  }
  
  if(!func)     cout << "Could not find func"     << endl;
  if(!func_reg) cout << "Could not find func_reg" << endl;

  TLegend* leg = new TLegend(0.11,0.71,0.45,0.94,NULL,"brNDC");
  leg->SetHeader( header );
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.06); 

  float normFunc = 0.;  
  if(func) {
    leg->AddEntry(func, "Reco", "L");
    func->SetNpx(1000);
    normFunc = func->GetParameter(0);
    cout << "Bin: " << func->Integral(-800,800) << endl;
  }
  if(func_reg){
    leg->AddEntry(func_reg,  "Regression", "L");
    func_reg->SetNpx(1000);
    func_reg->SetParameter(0, normFunc);
    cout << "Param 2G: " << func_reg->Integral(-800,800) << endl;
    if( TMath::Abs(func_reg->Integral(-800,800)/func->Integral(-800,800)-1)>1e-02  )
      func_reg->SetParameter(0, func_reg->GetParameter(0)*func->Integral(-800,800)/func_reg->Integral(-800,800)  );
    cout << "RESCALE: Param 2G: " << func_reg->Integral(-800,800) << endl;
  }

  hBestMass->GetXaxis()->SetRangeUser(0,350);
  hBestMass->Scale(0.);
  hBestMass->Draw("P");
  if( func )     func    ->Draw("SAME");
  if( func_reg)  func_reg->Draw("SAME");
  leg->Draw();

  if(0){
    string eventStr = event<10 ? string(Form("0%d",event)) : string(Form("%d",event));
    c1->SaveAs(Form("Plots/Plot_%s_%s_%d%s.png", flavor.c_str(), eventStr.c_str(), bin, extraLabel.Data()));
    //c1->SaveAs(Form("Plots/Plot_%s_%d_%d.pdf", flavor.c_str(), event, bin));
  }


  return;
}




void draw(vector<float> param, TTree* t = 0, TString var = "", TH1F* h = 0, TCut cut = ""){

  h->Reset();
  if(param.size()==0){
    t->Draw(var+">>"+TString(h->GetName()),   "weight"*cut );
    return;
  }
  else{
    if( param.size()!=2 ){
      cout << "Error in draw()" << endl;
      return;
    }

    float cutval = h->GetBinLowEdge( 3 );
    TCut    cut1( Form("%s<%f", var.Data(), cutval) );    
    TString var2(Form("p_125_all_b_ttbb/((p_125_all_b_ttbb+%f*p_125_all_b_ttjj))",    param[0]));

    TCut    cut2(Form("p_125_all_b_ttbb/((p_125_all_b_ttbb+%f*p_125_all_b_ttjj))<%f", param[0],param[1]));    
    t->Draw(var2+">>"+TString(h->GetName()),   "weight"*(cut&&cut1&&cut2) );
    float binL    = h->Integral();
    float binLErr = h->GetEntries()>0 ? sqrt(h->GetEntries())*h->Integral()/h->GetEntries() : 0.;
    h->Reset();
    t->Draw(var2+">>"+TString(h->GetName()),   "weight"*(cut&&cut1&&(!cut2)) );
    float binH    = h->Integral();
    float binHErr = h->GetEntries()>0 ? sqrt(h->GetEntries())*h->Integral()/h->GetEntries() : 0.;
    h->Reset();

    t->Draw(var+">>"+TString(h->GetName()),   "weight"*cut );
    h->SetBinContent(1, binL);  h->SetBinError(1, binLErr);
    h->SetBinContent(2, binH);  h->SetBinError(2, binHErr);
  }

  return;
}


void plot_category(string cat = "SL_VType2_nominal", 
		   string  cut = "type==0",
		   float fact1 = 0.02,
		   float fact2 = 0.,
		   float lumiScale = 3.2,
		   string header = "",
		   string fname  = "",
		   int newvar = 1,
		   int nBins  = 6,
		   int splitFirstBin = 0
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

  THStack* aStack = new THStack("aStack","Simulation #sqrt{s}=8 TeV, L=19.5 fb^{-1};  L_{S}/(L_{S}+L_{B}) ; events ");
  //TH1F* hStack = (TH1F*)aStack->GetHistogram();
  //hStack->SetTitle(cutName.c_str());
  //hStack->SetXTitle(xTitle.c_str());
  //hStack->SetYTitle(yTitle.c_str());

 //  vector<string> samples;
//   samples.push_back("SingleT");
//   samples.push_back("TTV");
//   samples.push_back("EWK");
//   samples.push_back("DiBoson");
//   if( cat.find("VType2")!=string::npos || cat.find("VType3")!=string::npos ){
//     samples.push_back("TTJetsSemiLept");
//     samples.push_back("TTJetsSemiLept");
//   }
//   if( cat.find("VType0")!=string::npos || cat.find("VType1")!=string::npos ){
//     samples.push_back("TTJetsFullLept");
//     samples.push_back("TTJetsFullLept");
//   }
//   samples.push_back("TTH125");

  string version = cat.find("SL")!=string::npos ? "_v6" : "_v5";


  TString var("");
  if(newvar)
    var = TString(Form("p_125_all_s_ttbb/(p_125_all_s_ttbb+%f*((1-%f)*p_125_all_b_ttbb+%f*p_125_all_b_ttjj))", fact1, fact2, fact2));
  else
    var = TString(Form("p_125_all_s/(p_125_all_s+p_125_all_b*%f)", fact1));
  vector<float> param;
  param.clear();
  if(splitFirstBin!=0){
    param.push_back(1.0);
    param.push_back(0.5);
  }

  TArrayF bins (nBins+1);
  TArrayF bins2(nBins+2);
  cout << "Making histograms with " << nBins << " bins:" << endl;
  if( param.size()==0){
    for(int b = 0; b < nBins+1; b++)
      bins[b] = b*1.0/nBins;
  }
  else{
    for(int b = 0; b < nBins+2; b++){
      if(b<=2) bins2[b] = b*0.5/(nBins);
      else     bins2[b] = (b-1)*1.0/(nBins);
      cout <<  bins2[b] << ", ";
    }
    cout << endl;
  }



  vector<string> samples;
  if( cat.find("SL")!=string::npos ){
    samples.push_back("nominal"+version+"_SingleT");
    samples.push_back("nominal"+version+"_TTV");

    //samples.push_back("nominal"+version+"_EWK");
    //samples.push_back("nominal"+version+"_DiBoson");
    samples.push_back("nominal"+version+"_TTJets");
    samples.push_back("nominal"+version+"_TTJets");
    samples.push_back("nominal"+version+"_TTH125");

    //samples.push_back("VType2_nominal_wtag_TTJetsSemiLept");
    //samples.push_back("VType2_nominal_wtag_TTJetsSemiLept");
    //samples.push_back("VType2_nominal_wtag_TTH125");
  }
  if( cat.find("DL")!=string::npos ){
    samples.push_back("nominal"+version+"_SingleT");
    samples.push_back("nominal"+version+"_TTV");
    //samples.push_back("nominal"+version+"_EWK");
    //samples.push_back("nominal"+version+"_DiBoson");
    samples.push_back("nominal"+version+"_TTJetsFullLept");
    samples.push_back("nominal"+version+"_TTJetsFullLept");
    samples.push_back("nominal"+version+"_TTH125");
  }

  int countTTJets  = 0;
  int countSingleT = 0;
  int countTTV     = 0;
  int countEWK     = 0;
  int countDiBoson = 0;
  int countTTH     = 0;

  TH1F* hS = new TH1F("hS","Simulation 8 TeV",  param.size()==0 ? nBins : nBins+1 ,  param.size()==0 ? bins.GetArray() : bins2.GetArray() );
  hS->SetLineWidth(3);
  hS->SetLineStyle(kDashed);
  hS->SetLineColor(kRed);

  TH1F* hErr = new TH1F("hErr","Simulation 8 TeV L=19.5 fb^{-1}" , param.size()==0 ? nBins : nBins+1 ,  param.size()==0 ? bins.GetArray() : bins2.GetArray() );
  hErr->Sumw2();
  leg->AddEntry(hErr, "MC unc. (stat.)", "L");

  for(unsigned int sample = 0; sample < samples.size(); sample++){

    cout << "Doing " << samples[sample] << endl;
    TFile* f = TFile::Open(("MEAnalysis_"+cat+"_"+samples[sample]+".root").c_str());
    if(f==0 || f->IsZombie() ) continue;
    TTree* tree     = (TTree*)f->Get("tree");
    //TH1F*  hcounter = (TH1F*) f->Get("hcounter");

    TCut tcut(cut.c_str());

    //tcut = tcut && TCut("p_125_all_b_ttbb/((p_125_all_b_ttbb+2*p_125_all_b_ttjj))>0.18");

    if( samples[sample].find("TTJets") != string::npos && countTTJets==0){
      countTTJets++;
      tcut = tcut&&TCut("nSimBs>2");
    }
    else if( samples[sample].find("TTJets") != string::npos && countTTJets==1){
      countTTJets++;
      tcut = tcut&&TCut("nSimBs<=2");
    }

    TH1F* h = new TH1F(("h_"+samples[sample]).c_str(),"Simulation 8 TeV L=19.5 fb^{-1}; events; L_{S}/{L_{S}+L_{B}}",  param.size()==0 ? nBins : nBins+1 ,  param.size()==0 ? bins.GetArray() : bins2.GetArray() );
    h->Sumw2();
    draw( param, tree , var, h , tcut );

    //if(newvar)
    //tree->Draw(Form("p_125_all_s_ttbb/(p_125_all_s_ttbb+%f*((1-%f)*p_125_all_b_ttbb+%f*p_125_all_b_ttjj))>>h_%s", fact1, fact2, fact2, samples[sample].c_str()), "weight"*tcut);
    //else
    //tree->Draw(Form("p_125_all_s/(p_125_all_s+%f*p_125_all_b)>>h_%s", fact1, samples[sample].c_str()), "weight"*tcut);
    //h->Scale(1./hcounter->GetBinContent(1) * lumiScale);

    h->Scale(lumiScale);
    cout << "Events = " << h->Integral() << endl;

    //for(int b = 1; b<= h->GetNbinsX(); b++ ){
      //h->SetBinContent( b, h->GetBinContent(b)/h->GetBinWidth(b) );
      //h->SetBinError  ( b, h->GetBinError  (b)/h->GetBinWidth(b) );
      //if( param.size()>0 ){
      //if(b==1)
      //h->GetXaxis()->SetBinLabel(b, "L_{LH}<0.5" );
      //else if(b==2) 
      //  h->GetXaxis()->SetBinLabel(b, "L_{LH}>0.5" );
      //else
      //  h->GetXaxis()->SetBinLabel(b, Form("%.1f",  h->GetBinCenter(b)) );
      //h->GetXaxis()->SetLabelSize(0.10);
      //h->GetXaxis()->SetNdivisions(nBins*100 + 2*nBins,kTRUE);	  
      //}
    //}


    if(samples[sample].find("TTH125") == string::npos){
      hErr->Add( h, 1.0);
    }
   
    if( samples[sample].find("TTH125") != string::npos ){
      h->SetLineColor( kRed );
      h->SetFillColor( kRed );
      if(countTTH==0) leg->AddEntry(h, "t#bar{t}H", "F");
      countTTH++;
      hS->Add( h, 5.);
    }
    if( samples[sample].find("TTJets") != string::npos ){
      if(countTTJets==1){
	leg->AddEntry(h, "t#bar{t}+jets HF", "F");
	h->SetLineColor( 17 );
	h->SetFillColor( 17 );
      }
      else if(countTTJets==2){
	h->SetLineColor( 18 );
	//h->SetLineWidth( 0 );
	leg->AddEntry(h, "t#bar{t}+jets LF", "F");
	h->SetFillColor( 18 );
      }
    }
    if( samples[sample].find("SingleT") != string::npos ){
      h->SetLineColor( kMagenta );
      h->SetFillColor( kMagenta );
      if(countSingleT==0) leg->AddEntry(h, "Single top", "F");
      countSingleT++;
    }
    if( samples[sample].find("EWK") != string::npos ){
      h->SetLineColor( kGreen );
      h->SetFillColor( kGreen );
      if(countEWK==0) leg->AddEntry(h, "V+jets", "F");
      countEWK++;
    }
    if( samples[sample].find("DiBoson") != string::npos ){
      h->SetLineColor( kYellow );
      h->SetFillColor( kYellow );
      if(countDiBoson==0) leg->AddEntry(h, "VV", "F");
      countDiBoson++;
    }
    if( samples[sample].find("TTV") != string::npos ){
      h->SetLineColor( 30 );
      h->SetFillColor( 30 );
      if(countTTV==0) leg->AddEntry(h, "t#bar{t}V", "F");
      countTTV++;
    }

    aStack->Add( h );
  }

  hErr->GetYaxis()->SetTitle("Events");
  hErr->GetXaxis()->SetTitle("L_{S} = P_{tth} / ( P_{tth} + P_{ttbb} + P_{ttjj} )");
  hErr->SetTitleSize  (0.04,"X");
  hErr->SetTitleOffset(0.95,"X");
  float max =  hErr->GetMaximum()*1.35;
  hErr->GetYaxis()->SetRangeUser(0., max );
  hErr->SetLineColor(kBlack);
  hErr->Draw("HISTE1");
  aStack->Draw("HISTSAME");
  leg->SetHeader( header.c_str() );
  leg->AddEntry(hS, "signal x 5", "L");
  hS->Draw("HISTSAME");
  hErr->Draw("HISTE1SAME");
  leg->Draw();

  TLine* line = new TLine(hErr->GetBinLowEdge(3), max , hErr->GetBinLowEdge(3), 0.);
  line->SetLineWidth(4);
  line->SetLineStyle(kSolid);
  line->SetLineColor(kBlack);
  if( param.size()>0 ) line->Draw("SAME");

  TPaveText *pt1 = new TPaveText(0.101, 0.839161, 0.198142, 0.895105,"brNDC");
  pt1->SetFillStyle(1001);
  pt1->SetBorderSize(0);
  pt1->SetFillColor(kWhite);
  pt1->SetTextSize(0.03); 
  pt1->AddText("L_{bb}<0.5");

  TPaveText *pt2 = new TPaveText(0.191, 0.839161, 0.294118, 0.895105,"brNDC");
  pt2->SetFillStyle(1001);
  pt2->SetBorderSize(0);
  pt2->SetFillColor(kWhite);
  pt2->SetTextSize(0.03); 
  pt2->AddText("L_{bb}>0.5");
 
  if( param.size()>0 ){
    pt1->Draw();
    pt2->Draw();
  }

  cout << "Signal = " << hS->Integral()/5. << endl;

  c1->SaveAs(  (fname+version+".png").c_str() );
  c1->SaveAs(  (fname+version+".pdf").c_str() );
  

}



void plotAll(){


  
  plot_category( "SL", 
		 "type==0",
		 0.02, 0.50,
		 19.5/12.1,
		 "Cat1 (4b2j W-tag)",
		 "Plot_SL_Cat1",
		 1,
		 4,
		 1
		 );
  
  
  plot_category( "SL", 
		 "type==1",
		 0.02, 0.70, 
		 19.5/12.1,
		 "Cat2 (4b2j !W-tag)",
		 "Plot_SL_Cat2",
		 1,
		 5,
		 1
		 );


  plot_category( "SL", 
		 "type==2 && flag_type2>0",
		 0.02, 0.80,
		 19.5/12.1,
		 "Cat3 (4b1j W-tag)",
		 "Plot_SL_Cat3",
		 1,
		 5,
		 1
		 );
  

  plot_category( "SL", 
		 "type==2 && flag_type2<=0",
		 0.02, 0.70,
		 19.5/12.1,
		 "Cat4 (4b1j !W-tag)",
		 "Plot_SL_Cat4",
		 1,
		 5,
		 1
		 );


  plot_category( "SL", 
		 "type==3",
		 0.02, 0.50,
		 19.5/12.1,
		 "Cat5 (4b3j)",
		 "Plot_SL_Cat5",
		 1,
		 5,
		 1
		 );

  
  plot_category( "DL", 
		 "type==6",
		 0.02, 0,
		 19.5/12.1 * 2,
		 "Cat6 (4b), tight",
		 "Plot_DL_Cat6",
		 0,
		 5,
		 0
		 );

  plot_category( "DL", 
		 "type==7",
		 0.02, 0,
		 19.5/12.1 * 2,
		 "Cat7 (4b), loose",
		 "Plot_DL_Cat7",
		 0,
		 6,
		 0
		 );
  
  
  

}





void plotTF(){

  int bin;

  
  bin = 11;
  plot_Slice(Form("#splitline{udcs-quarks}{E_{q}#in[%.0f,%.0f], |#eta_{q}|#in[0.0,1.5]}", 30.+(bin-1.)*4. , 30.+bin*4.), "jet energy [GeV]", "", bin, "Light", 0, "" , "N(E,#mu,#sigma)");
  bin = 21;
  plot_Slice(Form("#splitline{udcs-quarks}{E_{q}#in[%.0f,%.0f], |#eta_{q}|#in[0.0,1.5]}", 30.+(bin-1.)*4. , 30.+bin*4.), "jet energy [GeV]", "", bin, "Light", 0, "" , "N(E,#mu,#sigma)");
  bin = 31;
  plot_Slice(Form("#splitline{udcs-quarks}{E_{q}#in[%.0f,%.0f], |#eta_{q}|#in[0.0,1.5]}", 30.+(bin-1.)*4. , 30.+bin*4.), "jet energy [GeV]", "", bin, "Light", 0, "" , "N(E,#mu,#sigma)");
  bin = 41;
  plot_Slice(Form("#splitline{udcs-quarks}{E_{q}#in[%.0f,%.0f], |#eta_{q}|#in[0.0,1.5]}", 30.+(bin-1.)*4. , 30.+bin*4.), "jet energy [GeV]", "", bin, "Light", 0, "" , "N(E,#mu,#sigma)");

  bin = 11;
  plot_Slice(Form("#splitline{udcs-quarks}{E_{q}#in[%.0f,%.0f], |#eta_{q}|#in[1.5,2.5]}", 30.+(bin-1.)*4. , 30.+bin*4.), "jet energy [GeV]", "", bin, "Light", 1, "" , "N(E,#mu,#sigma)");
  bin = 21;
  plot_Slice(Form("#splitline{udcs-quarks}{E_{q}#in[%.0f,%.0f], |#eta_{q}|#in[1.5,2.5]}", 30.+(bin-1.)*4. , 30.+bin*4.), "jet energy [GeV]", "", bin, "Light", 1, "" , "N(E,#mu,#sigma)");
  bin = 31;
  plot_Slice(Form("#splitline{udcs-quarks}{E_{q}#in[%.0f,%.0f], |#eta_{q}|#in[1.5,2.5]}", 30.+(bin-1.)*4. , 30.+bin*4.), "jet energy [GeV]", "", bin, "Light", 1, "" , "N(E,#mu,#sigma)");
  bin = 41;
  plot_Slice(Form("#splitline{udcs-quarks}{E_{q}#in[%.0f,%.0f], |#eta_{q}|#in[1.5,2.5]}", 30.+(bin-1.)*4. , 30.+bin*4.), "jet energy [GeV]", "", bin, "Light", 1, "" , "N(E,#mu,#sigma)");

  bin = 11;
  plot_Slice(Form("#splitline{b-quarks}{E_{b}#in[%.0f,%.0f], |#eta_{b}|#in[0.0,1.5]}", 30.+(bin-1.)*4. , 30.+bin*4.), "jet energy [GeV]", "", bin, "Heavy", 0, "" , "fN(E,#mu,#sigma)+(1-f)N(E,#mu',#sigma')");
  bin = 21;
  plot_Slice(Form("#splitline{b-quarks}{E_{b}#in[%.0f,%.0f], |#eta_{b}|#in[0.0,1.5]}", 30.+(bin-1.)*4. , 30.+bin*4.), "jet energy [GeV]", "", bin, "Heavy", 0, "" , "fN(E,#mu,#sigma)+(1-f)N(E,#mu',#sigma')");
  bin = 31;
  plot_Slice(Form("#splitline{b-quarks}{E_{b}#in[%.0f,%.0f], |#eta_{b}|#in[0.0,1.5]}", 30.+(bin-1.)*4. , 30.+bin*4.), "jet energy [GeV]", "", bin, "Heavy", 0, "" , "fN(E,#mu,#sigma)+(1-f)N(E,#mu',#sigma')");
  bin = 41;
  plot_Slice(Form("#splitline{b-quarks}{E_{b}#in[%.0f,%.0f], |#eta_{b}|#in[0.0,1.5]}", 30.+(bin-1.)*4. , 30.+bin*4.), "jet energy [GeV]", "", bin, "Heavy", 0, "" , "fN(E,#mu,#sigma)+(1-f)N(E,#mu',#sigma')");

  bin = 11;
  plot_Slice(Form("#splitline{b-quarks}{E_{b}#in[%.0f,%.0f], |#eta_{b}|#in[1.5,2.5]}", 30.+(bin-1.)*4. , 30.+bin*4.), "jet energy [GeV]", "", bin, "Heavy", 1, "" , "fN(E,#mu,#sigma)+(1-f)N(E,#mu',#sigma')");
  bin = 21;
  plot_Slice(Form("#splitline{b-quarks}{E_{b}#in[%.0f,%.0f], |#eta_{b}|#in[1.5,2.5]}", 30.+(bin-1.)*4. , 30.+bin*4.), "jet energy [GeV]", "", bin, "Heavy", 1, "" , "fN(E,#mu,#sigma)+(1-f)N(E,#mu',#sigma')");
  bin = 31;
  plot_Slice(Form("#splitline{b-quarks}{E_{b}#in[%.0f,%.0f], |#eta_{b}|#in[1.5,2.5]}", 30.+(bin-1.)*4. , 30.+bin*4.), "jet energy [GeV]", "", bin, "Heavy", 1, "" , "fN(E,#mu,#sigma)+(1-f)N(E,#mu',#sigma')");
  bin = 41;
  plot_Slice(Form("#splitline{b-quarks}{E_{b}#in[%.0f,%.0f], |#eta_{b}|#in[1.5,2.5]}", 30.+(bin-1.)*4. , 30.+bin*4.), "jet energy [GeV]", "", bin, "Heavy", 1, "" , "fN(E,#mu,#sigma)+(1-f)N(E,#mu',#sigma')");
 

  

  return;
  
  plot_param( "udcs-quarks,  |#eta_{q}|#in[0.0,1.5]", "parton energy [GeV]", "#sigma_{q} [GeV]", "resolLightBin0", 40, 250, 0, 50, "a #oplus b#sqrt{E} #oplus cE", "");
  plot_param( "udcs-quarks,  |#eta_{q}|#in[1.5,2.5]", "parton energy [GeV]", "#sigma_{q} [GeV]", "resolLightBin1", 40, 250, 0, 50, "a #oplus b#sqrt{E} #oplus cE", "");
  
  plot_param( "udcs-quarks,  |#eta_{q}|#in[0.0,1.5]", "parton energy [GeV]", "#mu_{q} [GeV]",    "respLightBin0",  40, 250, 40, 250, "mE", "");
  plot_param( "udcs-quarks,  |#eta_{q}|#in[1.5,2.5]", "parton energy [GeV]", "#mu_{q} [GeV]",    "respLightBin1",  40, 250, 40, 250, "mE", "");

  plot_param( "b-quarks,  |#eta_{b}|#in[0.0,1.5]",    "parton energy [GeV]", "#sigma_{b} [GeV]", "resolG1HeavyBin0", 40, 250, 0, 70, "a #oplus b#sqrt{E} #oplus cE", "");
  plot_param( "b-quarks,  |#eta_{b}|#in[1.5,2.5]",    "parton energy [GeV]", "#sigma_{b} [GeV]", "resolG1HeavyBin1", 40, 250, 0, 70, "a #oplus b#sqrt{E} #oplus cE", "");

  plot_param( "b-quarks,  |#eta_{b}|#in[0.0,1.5]",    "parton energy [GeV]", "#sigma'_{b} [GeV]", "resolG2HeavyBin0", 40, 250, 0, 70, "a #oplus b#sqrt{E} #oplus cE", "");
  plot_param( "b-quarks,  |#eta_{b}|#in[1.5,2.5]",    "parton energy [GeV]", "#sigma'_{b} [GeV]", "resolG2HeavyBin1", 40, 250, 0, 70, "a #oplus b#sqrt{E} #oplus cE", "");

  plot_param( "b-quarks,  |#eta_{b}|#in[0.0,1.5]", "parton energy [GeV]", "#mu_{b} [GeV]",    "respG1HeavyBin0",  40, 250, 40, 250, "n+mE", "");
  plot_param( "b-quarks,  |#eta_{b}|#in[1.5,2.5]", "parton energy [GeV]", "#mu_{b} [GeV]",    "respG1HeavyBin1",  40, 250, 40, 250, "n+mE", "");
 
  plot_param( "b-quarks,  |#eta_{b}|#in[0.0,1.5]", "parton energy [GeV]", "#mu'_{b} [GeV]",    "respG2HeavyBin0",  40, 250, 40, 250, "n+mE", "");
  plot_param( "b-quarks,  |#eta_{b}|#in[1.5,2.5]", "parton energy [GeV]", "#mu'_{b} [GeV]",    "respG2HeavyBin1",  40, 250, 40, 250, "n+mE", "");
 
  plot_param( "E_{x,y}^{miss}", "#sum E_{T} [GeV]", "#sigma_{E_{T}} [GeV]", "hWidthsPxResol", 500, 3000, 0, 50, "a #oplus b#sqrt{E}", "");

  plot_param( "b-quarks,  |#eta_{b}|#in[0.0,1.5]", "#xi", "df/d#xi", "csv_b_Bin0__csvReco", 0.679, 1.0, 0, 35, "", "");
  plot_param( "b-quarks,  |#eta_{b}|#in[1.5,2.5]", "#xi", "df/d#xi", "csv_b_Bin1__csvReco", 0.679, 1.0, 0, 35, "", "");

  plot_param( "uds-quarks,  |#eta_{q}|#in[0.0,1.5]", "#xi", "df/d#xi", "csv_l_Bin0__csvReco", 0.679, 1.0, 0, 35, "", "");
  plot_param( "uds-quarks,  |#eta_{q}|#in[1.5,2.5]", "#xi", "df/d#xi", "csv_l_Bin1__csvReco", 0.679, 1.0, 0, 35, "", "");

  return;

}


void plot_MadWeight_comp(){

  plot_analysis_MW("SL","type==1 || type==3","SL","SL",0.0145);
  plot_analysis_MW("DL","type>=6","DL","DL",0.012);
}


void plot_ttbb(){

  plot_analysis_ttbb("SL", "type==0", "SL Cat. 1",                  "SL_cat1_bbjj_SB", 0.02, 0.50);
  plot_analysis_ttbb("SL", "type==1", "SL Cat. 2",                  "SL_cat2_bbjj_SB", 0.02, 0.70);
  plot_analysis_ttbb("SL", "type==2 && flag_type2>0",  "SL Cat. 3", "SL_cat3_bbjj_SB", 0.02, 0.80);
  plot_analysis_ttbb("SL", "type==2 && flag_type2<=0", "SL Cat. 4", "SL_cat4_bbjj_SB", 0.02, 0.70);
  plot_analysis_ttbb("SL", "type==3", "SL Cat. 5",                  "SL_cat5_bbjj_SB", 0.02, 0.50);
  plot_analysis_ttbb("SL", "type<=3", "SL All Cat.",                "SL_Comb_bbjj_SB", 0.02, 0.50);

  plot_analysis_ttbb("DL", "type==6", "DL Cat. 6",                  "DL_cat6_bbjj", 0.02, 0.);
  plot_analysis_ttbb("DL", "type==7", "DL Cat. 7",                  "DL_cat7_bbjj", 0.02, 0.);
  plot_analysis_ttbb("DL", "type>=6", "DL All Cat.",                "DL_Comb_bbjj", 0.02, 0.);
  

}


void plot_bj(){

  //plot_analysis_bj("SL", "type<=3", "SL All Cat.", "SL_Comb_BJ_AN_ttbb_vs_ttH",
  //	   "nSimBs>2  && nMatchSimBs>=-999", "nSimBs>2" );


  //return;

  plot_analysis_bj("SL", "type<=3", "SL All Cat.", "SL_Comb_BJ_AN_ttbb_vs_ttjj",
		   "nSimBs>2  && nMatchSimBs>=-999", "nSimBs==2  && nMatchSimBs>=-999" );
  
  plot_analysis_bj("SL", "type<=3", "SL All Cat.", "SL_Comb_BJ_AN_ttbb_vs_ttb",
		   "nSimBs>2  && nMatchSimBs>=2", "nSimBs>2  && nMatchSimBs<2" );

}



void plotTFAll(float binW=(300.-30.)/58., TString extra = ""){

  for(int b = 0 ; b <2 ; b++){
    for(int s = 1; s <= 58; s++){
      plot_Slice(string(Form("b-jet %.1f<|#eta|<%.1f, %.0f<E<%.0f",     b==0 ? 0. : 1., b==0 ? 1.0 : 2.5, 30+(s-1)*binW,30+s*binW)),"jet energy (GeV)","f(#hat{E}|E)", s, "Heavy", b, extra, "2G fit");
      plot_Slice(string(Form("light-jet %.1f<|#eta|<%.1f, %.0f<E<%.0f", b==0 ? 0. : 1., b==0 ? 1.0 : 2.5, 30+(s-1)*binW,30+s*binW)),"jet energy (GeV)","f(#hat{E}|E)", s, "Light", b, extra, "1G fit");
    }
  }

}
