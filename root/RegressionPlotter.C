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

  int firstFilledBin=-99;
  for(int bb = 1; bb < hBestMass->GetNbinsX() && firstFilledBin<0; bb++){
    if(hBestMass->GetBinContent(bb)>0) firstFilledBin=bb;
  }
  float firstVal = firstFilledBin>0 ? hBestMass->GetBinCenter(firstFilledBin) : 0.;
  cout << "firstVal = " << firstVal << endl;

  float normFunc = 0.;  
  if(func) {
    func->SetLineColor(kBlack);
    func->SetLineStyle(kDashed);
    leg->AddEntry(func, param, "L");
    func->SetNpx(1000);
    normFunc = func->GetParameter(0);
    cout << "Bin: " << func->Integral(-800,800) << endl;
    func->SetParameter(0, func->GetParameter(0)*hBestMass->Integral("width")/func->Integral(-800,800) / (1- func->Integral(-800,firstVal)/func->Integral(-800, 800) )  );
    cout << "RESCALE: Bin 2G: " << func->Integral(-800,800) << endl;
  }
  if(funcP){
    funcP->SetLineWidth(2);
    leg->AddEntry(funcP,  "2G parametrization", "L");
    funcP->SetNpx(1000);
    funcP->SetParameter(0, normFunc);
    cout << "Param 2G: " << funcP->Integral(-800,800) << endl;
    funcP->SetParameter(0, funcP->GetParameter(0)*hBestMass->Integral("width")/funcP->Integral(-800,800)  / (1- funcP->Integral(-800,firstVal)/funcP->Integral(-800, 800) ) );
    cout << "RESCALE: Param 2G: " << funcP->Integral(-800,800) << endl;
  }
  if(funcP2){
    funcP2->SetLineWidth(2);
    leg->AddEntry(funcP2, "1G parametrization", "L");
    funcP2->SetNpx(1000);
    funcP2->SetParameter(0, normFunc);
    cout << "Param 1G: " << funcP2->Integral(-800,800) << endl;
    funcP2->SetParameter(0, funcP2->GetParameter(0)*hBestMass->Integral("width")/funcP2->Integral(-800,800) /  (1- funcP2->Integral(-800,firstVal)/funcP2->Integral(-800, 800) ) );
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
		    TString param = "DG"
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
    func     = (TF1*)f    ->Get(TString(Form("Bin%s%d_%d/resol%s_fit%d_%d",flavor.c_str(), event, bin, param.Data(), event, bin)));
    func_reg = (TF1*)f_reg->Get(TString(Form("Bin%s%d_%d/resol%s_fit%d_%d",flavor.c_str(), event, bin, param.Data(), event, bin)));
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
    func->SetLineWidth(2);
    func->SetLineColor(kBlue);
    leg->AddEntry(func, "Standard ("+param+")", "L");
    func->SetNpx(1000);
    normFunc = func->GetParameter(0);
    cout << "Bin: " << func->Integral(-800,800) << endl;
  }
  if(func_reg){
    func_reg->SetLineWidth(2);
    leg->AddEntry(func_reg,  "Regression ("+param+")", "L");
    func_reg->SetNpx(1000);
    func_reg->SetLineColor(kRed);
    func_reg->SetParameter(0, normFunc);
    cout << "Param 2G: " << func_reg->Integral(-800,800) << endl;
    if( TMath::Abs(func_reg->Integral(-800,800)/func->Integral(-800,800)-1)>1e-02  )
      func_reg->SetParameter(0, func_reg->GetParameter(0)*func->Integral(-800,800)/func_reg->Integral(-800,800)  );
    cout << "RESCALE: Param 2G: " << func_reg->Integral(-800,800) << endl;
  }

  hBestMass->GetXaxis()->SetRangeUser(0,350);
  hBestMass->Reset();
  if(func) hBestMass->GetYaxis()->SetRangeUser(0, func->GetMaximum()*1.3);
  hBestMass->Draw("");
  if( func )     func    ->Draw("SAME");
  if( func_reg)  func_reg->Draw("SAME");
  leg->Draw();

  if(1){
    string eventStr = event<10 ? string(Form("0%d",event)) : string(Form("%d",event));
    c1->SaveAs(Form("Plots/Plot_%s_%s_%d%s.png", flavor.c_str(), eventStr.c_str(), bin, (extraLabel+param).Data()));
    //c1->SaveAs(Form("Plots/Plot_%s_%d_%d.pdf", flavor.c_str(), event, bin));
  }


  return;
}




void plotTFAll(float binW=(300.-30.)/58.){

  TString extra = ""; // standard

  
  for(int b = 0 ; b <2 ; b++){
    for(int s = 1; s <= 58; s++){
      plot_Slice(string(Form("b-jet %.1f<|#eta|<%.1f, %.0f<E<%.0f",     b==0 ? 0. : 1., b==0 ? 1.0 : 2.5, 30+(s-1)*binW,30+s*binW)),"jet energy (GeV)","f(#hat{E}|E)", s, "Heavy", b, extra, "2G fit");
      plot_Slice(string(Form("light-jet %.1f<|#eta|<%.1f, %.0f<E<%.0f", b==0 ? 0. : 1., b==0 ? 1.0 : 2.5, 30+(s-1)*binW,30+s*binW)),"jet energy (GeV)","f(#hat{E}|E)", s, "Light", b, extra, "1G fit");
    }
  }
  

  
  extra = "_reg"; // regression

  for(int b = 0 ; b < 2 ; b++){
    for(int s = 1; s <= 58; s++){
      plot_Slice(string(Form("b-jet %.1f<|#eta|<%.1f, %.0f<E<%.0f",     b==0 ? 0. : 1., b==0 ? 1.0 : 2.5, 30+(s-1)*binW,30+s*binW)),"jet energy (GeV)","f(#hat{E}|E)", s, "Heavy", b, extra, "2G fit");
    }
  }
 

  
  extra = "_comp"; // comparsion

  for(int b = 0 ; b <2 ; b++){
    for(int s = 1; s <= 58; s++){
      plot_Slice_Reg(string(Form("b-jet %.1f<|#eta|<%.1f, %.0f<E<%.0f",     b==0 ? 0. : 1., b==0 ? 1.0 : 2.5, 30+(s-1)*binW,30+s*binW)),"jet energy (GeV)","f(#hat{E}|E)", s, "Heavy", b, extra , "SG");
      plot_Slice_Reg(string(Form("b-jet %.1f<|#eta|<%.1f, %.0f<E<%.0f",     b==0 ? 0. : 1., b==0 ? 1.0 : 2.5, 30+(s-1)*binW,30+s*binW)),"jet energy (GeV)","f(#hat{E}|E)", s, "Heavy", b, extra , "DG");
    }
  }
  
  

}


void plot_correlation( TString target = "pt_gen", TString varRec = "e_rec_reg", TString varGen = "e_gen", 
		       TString cut = "(TMath::Abs(eta_part)<1.0)",  
		       TString xTitle = "gen-jet energy (GeV)", TString yTitle = "#frac{reco-gen}{gen}",
		       TString header = "b-jet 0#leq#eta#leq1, (regression)",
		       TString outname = "pt_gen_Bin0_e_gen"){

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
  c1->SetGrid(0,1);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TFile* f = TFile::Open("../bin/root/treeProducerTEST"+target+".root");
  TTree* tree = (TTree*)f->Get("genJetHeavyTree");

  TH2F* h2 = new TH2F("h2", "Simulation #sqrt{s}=8 TeV; "+xTitle+";"+yTitle, 100, 30, 300, 120, -0.9, 0.9);

  tree->Draw("("+varRec+"-"+varGen+")/"+varGen+":"+varGen+">>h2", varGen+">20 && "+varRec+">20 &&"+cut, "");

  h2->GetYaxis()->CenterTitle();
  h2->Scale( 1./h2->GetMaximum());
  h2->Draw("COLZ");

  TLegend* leg = new TLegend(0.11,0.026,0.45,0.25,NULL,"brNDC");
  leg->SetHeader( header );
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.06); 
  leg->Draw();

  if(1){
    c1->SaveAs("Plots/Plot_Regression_correlation_"+outname+".png");
  }

}

void plot_correlationAll(){

  plot_correlation("pt_gen",  "e_rec_reg", "e_part", "abs(eta_part)<1.0", "gen-parton energy (GeV)", "#frac{reco-gen}{gen}", "b-jet 0#leq#eta#leq1, (regression)", "pt_gen_Bin0_e_rec_reg_vs_e_part");
  plot_correlation("pt_gen",  "e_rec",     "e_part", "abs(eta_part)<1.0", "gen-parton energy (GeV)", "#frac{reco-gen}{gen}", "b-jet 0#leq#eta#leq1, (stdandard)",  "pt_gen_Bin0_e_rec_vs_e_part");

  plot_correlation("pt_gen",  "e_rec_reg", "e_part", "abs(eta_part)>1.0", "gen-parton energy (GeV)", "#frac{reco-gen}{gen}", "b-jet 1#leq#eta#leq2.5, (regression)", "pt_gen_Bin1_e_rec_reg_vs_e_part");
  plot_correlation("pt_gen",  "e_rec",     "e_part", "abs(eta_part)>1.0", "gen-parton energy (GeV)", "#frac{reco-gen}{gen}", "b-jet 1#leq#eta#leq2.5, (stdandard)",  "pt_gen_Bin1_e_rec_vs_e_part");

  plot_correlation("pt_gen",  "e_rec_reg", "e_gen", "abs(eta_part)<1.0", "gen-jet energy (GeV)", "#frac{reco-gen}{gen}", "b-jet 0#leq#eta#leq1, (regression)", "pt_gen_Bin0_e_rec_reg_vs_e_gen");
  plot_correlation("pt_gen",  "e_rec",     "e_gen", "abs(eta_part)<1.0", "gen-jet energy (GeV)", "#frac{reco-gen}{gen}", "b-jet 0#leq#eta#leq1, (stdandard)",  "pt_gen_Bin0_e_rec_vs_e_gen");

  plot_correlation("pt_gen",  "e_rec_reg", "e_gen", "abs(eta_part)>1.0", "gen-jet energy (GeV)", "#frac{reco-gen}{gen}", "b-jet 1#leq#eta#leq2.5, (regression)", "pt_gen_Bin1_e_rec_reg_vs_e_gen");
  plot_correlation("pt_gen",  "e_rec",     "e_gen", "abs(eta_part)>1.0", "gen-jet energy (GeV)", "#frac{reco-gen}{gen}", "b-jet 1#leq#eta#leq2.5, (stdandard)",  "pt_gen_Bin1_e_rec_vs_e_gen");


  plot_correlation("pt_part",  "e_rec_reg", "e_part", "abs(eta_part)<1.0", "gen-parton energy (GeV)", "#frac{reco-gen}{gen}", "b-jet 0#leq#eta#leq1, (regression)", "pt_part_Bin0_e_rec_reg_vs_e_part");
  plot_correlation("pt_part",  "e_rec",     "e_part", "abs(eta_part)<1.0", "gen-parton energy (GeV)", "#frac{reco-gen}{gen}", "b-jet 0#leq#eta#leq1, (stdandard)",  "pt_part_Bin0_e_rec_vs_e_part");

  plot_correlation("pt_part",  "e_rec_reg", "e_part", "abs(eta_part)>1.0", "gen-parton energy (GeV)", "#frac{reco-gen}{gen}", "b-jet 1#leq#eta#leq2.5, (regression)", "pt_part_Bin1_e_rec_reg_vs_e_part");
  plot_correlation("pt_part",  "e_rec",     "e_part", "abs(eta_part)>1.0", "gen-parton energy (GeV)", "#frac{reco-gen}{gen}", "b-jet 1#leq#eta#leq2.5, (stdandard)",  "pt_part_Bin1_e_rec_vs_e_part");

  plot_correlation("pt_part",  "e_rec_reg", "e_gen", "abs(eta_part)<1.0", "gen-jet energy (GeV)", "#frac{reco-gen}{gen}", "b-jet 0#leq#eta#leq1, (regression)", "pt_part_Bin0_e_rec_reg_vs_e_gen");
  plot_correlation("pt_part",  "e_rec",     "e_gen", "abs(eta_part)<1.0", "gen-jet energy (GeV)", "#frac{reco-gen}{gen}", "b-jet 0#leq#eta#leq1, (stdandard)",  "pt_part_Bin0_e_rec_vs_e_gen");

  plot_correlation("pt_part",  "e_rec_reg", "e_gen", "abs(eta_part)>1.0", "gen-jet energy (GeV)", "#frac{reco-gen}{gen}", "b-jet 1#leq#eta#leq2.5, (regression)", "pt_part_Bin1_e_rec_reg_vs_e_gen");
  plot_correlation("pt_part",  "e_rec",     "e_gen", "abs(eta_part)>1.0", "gen-jet energy (GeV)", "#frac{reco-gen}{gen}", "b-jet 1#leq#eta#leq2.5, (stdandard)",  "pt_part_Bin1_e_rec_vs_e_gen");

}
