#include <cstdlib>
#include <iostream> 
#include <fstream>
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
#include "TGraphAsymmErrors.h"
#include "TGraphPainter.h"
#include "TMultiGraph.h"
#include "TArrayF.h"
#include "TLine.h"



void bbb( TH1F* hin, TH1F* hout_Down, TH1F* hout_Up, int bin){

  float bin_down = hin->GetBinContent( bin ) - hin->GetBinError( bin );
  float bin_up   = hin->GetBinContent( bin ) + hin->GetBinError( bin );

  hout_Down->Reset();
  hout_Down->Add(hin,1.0);
  hout_Down->SetBinContent(bin, TMath::Max(bin_down,0.01));

  hout_Up->Reset();
  hout_Up->Add(hin,1.0);
  hout_Up->SetBinContent  (bin, TMath::Max(bin_up,0.01));

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




void produce( TString fname = "SL_VType2", string cut = "type==0", TString category = "cat0", float fact1 = 0.02, float fact2 = 0., float lumiScale = 20./12.1,
	      int newvar = 1, int nBins=6, int splitFirstBin=0){


  string version = (string(fname.Data())).find("SL")!=string::npos ? "_v6" : "_v5";
  //string version = "";

  cout << "Doing version " << version << " and category " << category << endl;

  // input files path
  TString inputpath = "gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/Trees/MEM/Oct04_2013/";




  float scaleTTH      = 1.0;//4.7;
  float scaleTTJetsLF = 1.0;//3.9;
  float scaleTTJetsHF = 1.0;//3.9;

  cout << "*****************************************************" << endl;
  cout << "Attention!!! The following rescaling will be applied:" << endl;
  cout << "TTH:       " << scaleTTH << endl;
  cout << "TTJets LF: " << scaleTTJetsLF << endl;
  cout << "TTJets HF: " << scaleTTJetsHF << endl;
  cout << "*****************************************************" << endl;

  vector<float> param;
  param.clear();
  if(splitFirstBin){
    param.push_back(1.0);
    param.push_back(0.5);
  }

  TFile* fout = new TFile(fname+".root","UPDATE");
  TDirectory* dir =  fout->GetDirectory( fname+"_"+category); 
  if( !dir) dir = fout->mkdir( fname+"_"+category ) ;

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


  TH1F* h = new TH1F("h","Simulation #sqrt{s}=8 TeV, "+fname+"; S/(S+B); units",  param.size()==0 ? nBins : nBins+1 ,  param.size()==0 ? bins.GetArray() : bins2.GetArray());
 
  TString var("");
  if(newvar)
    var = TString(Form("p_125_all_s_ttbb/(p_125_all_s_ttbb+%f*((1-%f)*p_125_all_b_ttbb+%f*p_125_all_b_ttjj))", fact1, fact2, fact2));
  else
    var = TString(Form("p_125_all_s/(p_125_all_s+p_125_all_b*%f)", fact1));

  string basecut = cut;
  TString ttjets = (string(fname.Data())).find("SL")!=string::npos ? "TTJets" : "TTJetsFullLept";

  ////////////////////////////// TTH
  
  // NOMINAL
  TFile* f_TTH125 = TFile::Open(inputpath+"MEAnalysis_"+fname+"_nominal"+version+"_TTH125.root");
  TTree* t_TTH125 = (TTree*)f_TTH125->Get("tree");
  TH1F* h_TTH125 = (TH1F*)h->Clone("h_TTH125");
  h_TTH125->Reset();
  h_TTH125->Sumw2();
  draw( param, t_TTH125, var, h_TTH125, TCut(cut.c_str()));
  //t_TTH125->Draw(var+">>h_TTH125",   "weight"*TCut(cut.c_str()) );
  h_TTH125->Scale(lumiScale*scaleTTH);
  dir->cd();
  h_TTH125->Write("TTH125", TObject::kOverwrite);

  // BBB
  for(int bin = 1; bin<=h_TTH125->GetNbinsX(); bin++ ){
    TH1F* h_TTH125_b_up   = (TH1F*)h->Clone(Form("h_TTH125_%d_Up",bin));
    TH1F* h_TTH125_b_down = (TH1F*)h->Clone(Form("h_TTH125_%d_Down",bin));
    bbb( h_TTH125, h_TTH125_b_up, h_TTH125_b_down, bin);

    dir->cd();
    h_TTH125_b_up  ->Write(Form("TTH125_TTH125%sbin%dUp",  category.Data(), bin), TObject::kOverwrite);
    h_TTH125_b_down->Write(Form("TTH125_TTH125%sbin%dDown",category.Data(), bin), TObject::kOverwrite);
 
  }

  // JEC UP
  TFile* f_TTH125_JECUp = TFile::Open(inputpath+"MEAnalysis_"+fname+"_JECUp"+version+"_TTH125.root");
  TTree* t_TTH125_JECUp = (TTree*)f_TTH125_JECUp->Get("tree");
  TH1F* h_TTH125_JECUp = (TH1F*)h->Clone("h_TTH125_JECUp");
  h_TTH125_JECUp->Reset();
  h_TTH125_JECUp->Sumw2();
  draw( param, t_TTH125_JECUp, var, h_TTH125_JECUp, TCut(cut.c_str()));
  //t_TTH125_JECUp->Draw(var+">>h_TTH125_JECUp",   "weight"*TCut(cut.c_str()) );
  h_TTH125_JECUp->Scale(lumiScale*scaleTTH);
  dir->cd();
  h_TTH125_JECUp->Write("TTH125_JECUp", TObject::kOverwrite);

  // JEC DOWN
  TFile* f_TTH125_JECDown = TFile::Open(inputpath+"MEAnalysis_"+fname+"_JECDown"+version+"_TTH125.root");
  TTree* t_TTH125_JECDown = (TTree*)f_TTH125_JECDown->Get("tree");
  TH1F* h_TTH125_JECDown = (TH1F*)h->Clone("h_TTH125_JECDown");
  h_TTH125_JECDown->Reset();
  h_TTH125_JECDown->Sumw2();
  draw( param, t_TTH125_JECDown, var, h_TTH125_JECDown, TCut(cut.c_str()));
  //t_TTH125_JECDown->Draw(var+">>h_TTH125_JECDown",   "weight"*TCut(cut.c_str()) );
  h_TTH125_JECDown->Scale(lumiScale*scaleTTH);
  dir->cd();
  h_TTH125_JECDown->Write("TTH125_JECDown", TObject::kOverwrite);

  // csv UP
  TFile* f_TTH125_csvUp = TFile::Open(inputpath+"MEAnalysis_"+fname+"_csvUp"+version+"_TTH125.root");
  TTree* t_TTH125_csvUp = (TTree*)f_TTH125_csvUp->Get("tree");
  TH1F* h_TTH125_csvUp = (TH1F*)h->Clone("h_TTH125_csvUp");
  h_TTH125_csvUp->Reset();
  h_TTH125_csvUp->Sumw2();
  draw( param, t_TTH125_csvUp, var, h_TTH125_csvUp, TCut(cut.c_str()));
  //t_TTH125_csvUp->Draw(var+">>h_TTH125_csvUp",   "weight"*TCut(cut.c_str()) );
  h_TTH125_csvUp->Scale(lumiScale*scaleTTH);
  dir->cd();
  h_TTH125_csvUp->Write("TTH125_csvUp", TObject::kOverwrite);

  // csv DOWN
  TFile* f_TTH125_csvDown = TFile::Open(inputpath+"MEAnalysis_"+fname+"_csvDown"+version+"_TTH125.root");
  TTree* t_TTH125_csvDown = (TTree*)f_TTH125_csvDown->Get("tree");
  TH1F* h_TTH125_csvDown = (TH1F*)h->Clone("h_TTH125_csvDown");
  h_TTH125_csvDown->Reset();
  h_TTH125_csvDown->Sumw2();
  draw( param, t_TTH125_csvDown, var, h_TTH125_csvDown, TCut(cut.c_str()));
  //t_TTH125_csvDown->Draw(var+">>h_TTH125_csvDown",   "weight"*TCut(cut.c_str()) );
  h_TTH125_csvDown->Scale(lumiScale*scaleTTH);
  dir->cd();
  h_TTH125_csvDown->Write("TTH125_csvDown", TObject::kOverwrite);

  ////////////////////////////// TTJets: HF

  cut = basecut+"&&nSimBs>2";

  // NOMINAL"+VERSION+"
  TFile* f_TTJetsSemiLept_HF = TFile::Open(inputpath+"MEAnalysis_"+fname+"_nominal"+version+"_"+ttjets+".root");
  TTree* t_TTJetsSemiLept_HF = (TTree*)f_TTJetsSemiLept_HF->Get("tree");
  TH1F* h_TTJetsSemiLept_HF = (TH1F*)h->Clone("h_TTJetsSemiLept_HF");
  h_TTJetsSemiLept_HF->Reset();
  h_TTJetsSemiLept_HF->Sumw2();
  draw( param, t_TTJetsSemiLept_HF, var, h_TTJetsSemiLept_HF, TCut(cut.c_str()));
  //t_TTJetsSemiLept_HF->Draw(var+">>h_TTJetsSemiLept_HF",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_HF->Scale(lumiScale*scaleTTJetsHF);
  dir->cd();
  h_TTJetsSemiLept_HF->Write("TTJetsHF", TObject::kOverwrite);

  // BBB
  for(int bin = 1; bin<=h_TTJetsSemiLept_HF->GetNbinsX(); bin++ ){
    TH1F* h_TTJetsSemiLept_HF_b_up   = (TH1F*)h->Clone(Form("h_TTJetsSemiLept_HF_%d_Up",bin));
    TH1F* h_TTJetsSemiLept_HF_b_down = (TH1F*)h->Clone(Form("h_TTJetsSemiLept_HF_%d_Down",bin));
    bbb( h_TTJetsSemiLept_HF, h_TTJetsSemiLept_HF_b_up, h_TTJetsSemiLept_HF_b_down, bin);

    dir->cd();
    h_TTJetsSemiLept_HF_b_up  ->Write(Form("TTJetsHF_TTJetsHF%sbin%dUp",    category.Data(),  bin), TObject::kOverwrite);
    h_TTJetsSemiLept_HF_b_down->Write(Form("TTJetsHF_TTJetsHF%sbin%dDown",  category.Data(),bin), TObject::kOverwrite);
 
  }

  // JEC UP
  TFile* f_TTJetsSemiLept_HF_JECUp = TFile::Open(inputpath+"MEAnalysis_"+fname+"_JECUp"+version+"_"+ttjets+".root");
  TTree* t_TTJetsSemiLept_HF_JECUp = (TTree*)f_TTJetsSemiLept_HF_JECUp->Get("tree");
  TH1F* h_TTJetsSemiLept_HF_JECUp = (TH1F*)h->Clone("h_TTJetsSemiLept_HF_JECUp");
  h_TTJetsSemiLept_HF_JECUp->Reset();
  h_TTJetsSemiLept_HF_JECUp->Sumw2();
  draw( param, t_TTJetsSemiLept_HF_JECUp , var, h_TTJetsSemiLept_HF_JECUp, TCut(cut.c_str()));
  //t_TTJetsSemiLept_HF_JECUp->Draw(var+">>h_TTJetsSemiLept_HF_JECUp",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_HF_JECUp->Scale(lumiScale*scaleTTJetsHF);
  dir->cd();
  h_TTJetsSemiLept_HF_JECUp->Write("TTJetsHF_JECUp", TObject::kOverwrite);

  // JEC DOWN
  TFile* f_TTJetsSemiLept_HF_JECDown = TFile::Open(inputpath+"MEAnalysis_"+fname+"_JECDown"+version+"_"+ttjets+".root");
  TTree* t_TTJetsSemiLept_HF_JECDown = (TTree*)f_TTJetsSemiLept_HF_JECDown->Get("tree");
  TH1F* h_TTJetsSemiLept_HF_JECDown = (TH1F*)h->Clone("h_TTJetsSemiLept_HF_JECDown");
  h_TTJetsSemiLept_HF_JECDown->Reset();
  h_TTJetsSemiLept_HF_JECDown->Sumw2();
  draw( param, t_TTJetsSemiLept_HF_JECDown, var, h_TTJetsSemiLept_HF_JECDown, TCut(cut.c_str()));
  //t_TTJetsSemiLept_HF_JECDown->Draw(var+">>h_TTJetsSemiLept_HF_JECDown",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_HF_JECDown->Scale(lumiScale*scaleTTJetsHF);
  dir->cd();
  h_TTJetsSemiLept_HF_JECDown->Write("TTJetsHF_JECDown", TObject::kOverwrite);

  // csv UP
  TFile* f_TTJetsSemiLept_HF_csvUp = TFile::Open(inputpath+"MEAnalysis_"+fname+"_csvUp"+version+"_"+ttjets+".root");
  TTree* t_TTJetsSemiLept_HF_csvUp = (TTree*)f_TTJetsSemiLept_HF_csvUp->Get("tree");
  TH1F* h_TTJetsSemiLept_HF_csvUp = (TH1F*)h->Clone("h_TTJetsSemiLept_HF_csvUp");
  h_TTJetsSemiLept_HF_csvUp->Reset();
  h_TTJetsSemiLept_HF_csvUp->Sumw2();
  draw( param, t_TTJetsSemiLept_HF_csvUp, var, h_TTJetsSemiLept_HF_csvUp, TCut(cut.c_str()));
  //t_TTJetsSemiLept_HF_csvUp->Draw(var+">>h_TTJetsSemiLept_HF_csvUp",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_HF_csvUp->Scale(lumiScale*scaleTTJetsHF);
  dir->cd();
  h_TTJetsSemiLept_HF_csvUp->Write("TTJetsHF_csvUp", TObject::kOverwrite);

  // csv DOWN
  TFile* f_TTJetsSemiLept_HF_csvDown = TFile::Open(inputpath+"MEAnalysis_"+fname+"_csvDown"+version+"_"+ttjets+".root");
  TTree* t_TTJetsSemiLept_HF_csvDown = (TTree*)f_TTJetsSemiLept_HF_csvDown->Get("tree");
  TH1F* h_TTJetsSemiLept_HF_csvDown = (TH1F*)h->Clone("h_TTJetsSemiLept_HF_csvDown");
  h_TTJetsSemiLept_HF_csvDown->Reset();
  h_TTJetsSemiLept_HF_csvDown->Sumw2();
  draw( param, t_TTJetsSemiLept_HF_csvDown, var, h_TTJetsSemiLept_HF_csvDown, TCut(cut.c_str()));
  //t_TTJetsSemiLept_HF_csvDown->Draw(var+">>h_TTJetsSemiLept_HF_csvDown",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_HF_csvDown->Scale(lumiScale*scaleTTJetsHF);
  dir->cd();
  h_TTJetsSemiLept_HF_csvDown->Write("TTJetsHF_csvDown", TObject::kOverwrite);


  ////////////////////////////// TTJets: LF

  cut = basecut+"&&nSimBs==2";

  // NOMINAL"+VERSION+"
  TFile* f_TTJetsSemiLept_LF = TFile::Open(inputpath+"MEAnalysis_"+fname+"_nominal"+version+"_"+ttjets+".root");
  TTree* t_TTJetsSemiLept_LF = (TTree*)f_TTJetsSemiLept_LF->Get("tree");
  TH1F* h_TTJetsSemiLept_LF = (TH1F*)h->Clone("h_TTJetsSemiLept_LF");
  h_TTJetsSemiLept_LF->Reset();
  h_TTJetsSemiLept_LF->Sumw2();
  draw( param, t_TTJetsSemiLept_LF, var, h_TTJetsSemiLept_LF, TCut(cut.c_str()));
  //t_TTJetsSemiLept_LF->Draw(var+">>h_TTJetsSemiLept_LF",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_LF->Scale(lumiScale*scaleTTJetsLF);
  dir->cd();
  h_TTJetsSemiLept_LF->Write("TTJetsLF", TObject::kOverwrite);

  // BBB
  for(int bin = 1; bin<=h_TTJetsSemiLept_LF->GetNbinsX(); bin++ ){
    TH1F* h_TTJetsSemiLept_LF_b_up   = (TH1F*)h->Clone(Form("h_TTJetsSemiLept_LF_%d_Up",bin));
    TH1F* h_TTJetsSemiLept_LF_b_down = (TH1F*)h->Clone(Form("h_TTJetsSemiLept_LF_%d_Down",bin));
    bbb( h_TTJetsSemiLept_LF, h_TTJetsSemiLept_LF_b_up, h_TTJetsSemiLept_LF_b_down, bin);

    dir->cd();
    h_TTJetsSemiLept_LF_b_up  ->Write(Form("TTJetsLF_TTJetsLF%sbin%dUp",  category.Data(),  bin), TObject::kOverwrite);
    h_TTJetsSemiLept_LF_b_down->Write(Form("TTJetsLF_TTJetsLF%sbin%dDown",category.Data(),bin), TObject::kOverwrite);
 
  }

  // JEC UP
  TFile* f_TTJetsSemiLept_LF_JECUp = TFile::Open(inputpath+"MEAnalysis_"+fname+"_JECUp"+version+"_"+ttjets+".root");
  TTree* t_TTJetsSemiLept_LF_JECUp = (TTree*)f_TTJetsSemiLept_LF_JECUp->Get("tree");
  TH1F* h_TTJetsSemiLept_LF_JECUp = (TH1F*)h->Clone("h_TTJetsSemiLept_LF_JECUp");
  h_TTJetsSemiLept_LF_JECUp->Reset();
  h_TTJetsSemiLept_LF_JECUp->Sumw2();
  draw( param, t_TTJetsSemiLept_LF_JECUp, var, h_TTJetsSemiLept_LF_JECUp, TCut(cut.c_str()));
  //t_TTJetsSemiLept_LF_JECUp->Draw(var+">>h_TTJetsSemiLept_LF_JECUp",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_LF_JECUp->Scale(lumiScale*scaleTTJetsLF);
  dir->cd();
  h_TTJetsSemiLept_LF_JECUp->Write("TTJetsLF_JECUp", TObject::kOverwrite);

  // JEC DOWN
  TFile* f_TTJetsSemiLept_LF_JECDown = TFile::Open(inputpath+"MEAnalysis_"+fname+"_JECDown"+version+"_"+ttjets+".root");
  TTree* t_TTJetsSemiLept_LF_JECDown = (TTree*)f_TTJetsSemiLept_LF_JECDown->Get("tree");
  TH1F* h_TTJetsSemiLept_LF_JECDown = (TH1F*)h->Clone("h_TTJetsSemiLept_LF_JECDown");
  h_TTJetsSemiLept_LF_JECDown->Reset();
  h_TTJetsSemiLept_LF_JECDown->Sumw2();
  draw( param, t_TTJetsSemiLept_LF_JECDown , var, h_TTJetsSemiLept_LF_JECDown, TCut(cut.c_str()));
  //t_TTJetsSemiLept_LF_JECDown->Draw(var+">>h_TTJetsSemiLept_LF_JECDown",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_LF_JECDown->Scale(lumiScale*scaleTTJetsLF);
  dir->cd();
  h_TTJetsSemiLept_LF_JECDown->Write("TTJetsLF_JECDown", TObject::kOverwrite);

  // csv UP
  TFile* f_TTJetsSemiLept_LF_csvUp = TFile::Open(inputpath+"MEAnalysis_"+fname+"_csvUp"+version+"_"+ttjets+".root");
  TTree* t_TTJetsSemiLept_LF_csvUp = (TTree*)f_TTJetsSemiLept_LF_csvUp->Get("tree");
  TH1F* h_TTJetsSemiLept_LF_csvUp = (TH1F*)h->Clone("h_TTJetsSemiLept_LF_csvUp");
  h_TTJetsSemiLept_LF_csvUp->Reset();
  h_TTJetsSemiLept_LF_csvUp->Sumw2();
  draw( param, t_TTJetsSemiLept_LF_csvUp, var, h_TTJetsSemiLept_LF_csvUp, TCut(cut.c_str()));
  //t_TTJetsSemiLept_LF_csvUp->Draw(var+">>h_TTJetsSemiLept_LF_csvUp",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_LF_csvUp->Scale(lumiScale*scaleTTJetsLF);
  dir->cd();
  h_TTJetsSemiLept_LF_csvUp->Write("TTJetsLF_csvUp", TObject::kOverwrite);

  // csv DOWN
  TFile* f_TTJetsSemiLept_LF_csvDown = TFile::Open(inputpath+"MEAnalysis_"+fname+"_csvDown"+version+"_"+ttjets+".root");
  TTree* t_TTJetsSemiLept_LF_csvDown = (TTree*)f_TTJetsSemiLept_LF_csvDown->Get("tree");
  TH1F* h_TTJetsSemiLept_LF_csvDown = (TH1F*)h->Clone("h_TTJetsSemiLept_LF_csvDown");
  h_TTJetsSemiLept_LF_csvDown->Reset();
  h_TTJetsSemiLept_LF_csvDown->Sumw2();
  draw( param, t_TTJetsSemiLept_LF_csvDown, var, h_TTJetsSemiLept_LF_csvDown, TCut(cut.c_str()));
  //t_TTJetsSemiLept_LF_csvDown->Draw(var+">>h_TTJetsSemiLept_LF_csvDown",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_LF_csvDown->Scale(lumiScale*scaleTTJetsLF);
  dir->cd();
  h_TTJetsSemiLept_LF_csvDown->Write("TTJetsLF_csvDown", TObject::kOverwrite);



  ////////////////////////////// DIBOSON
  cut = basecut;
  
  // NOMINAL
  TFile* f_DiBoson = TFile::Open(inputpath+"MEAnalysis_"+fname+"_nominal"+version+"_DiBoson.root");
  TTree* t_DiBoson = (TTree*)f_DiBoson->Get("tree");
  TH1F* h_DiBoson = (TH1F*)h->Clone("h_DiBoson");
  h_DiBoson->Reset();
  h_DiBoson->Sumw2();
  draw( param, t_DiBoson , var, h_DiBoson, TCut(cut.c_str()));
  //t_DiBoson->Draw(var+">>h_DiBoson",   "weight"*TCut(cut.c_str()) );
  h_DiBoson->Scale(lumiScale);
  dir->cd();
  h_DiBoson->Write("DiBoson", TObject::kOverwrite);

  // BBB
  for(int bin = 1; bin<=h_DiBoson->GetNbinsX(); bin++ ){
    TH1F* h_DiBoson_b_up   = (TH1F*)h->Clone(Form("h_DiBoson_%d_Up",bin));
    TH1F* h_DiBoson_b_down = (TH1F*)h->Clone(Form("h_DiBoson_%d_Down",bin));
    bbb( h_DiBoson, h_DiBoson_b_up, h_DiBoson_b_down, bin);

    dir->cd();
    h_DiBoson_b_up  ->Write(Form("DiBoson_DiBoson%sbin%dUp",  category.Data(), bin), TObject::kOverwrite);
    h_DiBoson_b_down->Write(Form("DiBoson_DiBoson%sbin%dDown",category.Data(), bin), TObject::kOverwrite);
 
  }


  ////////////////////////////// SINGLE T
  
  // NOMINAL
  TFile* f_SingleT = TFile::Open(inputpath+"MEAnalysis_"+fname+"_nominal"+version+"_SingleT.root");
  TTree* t_SingleT = (TTree*)f_SingleT->Get("tree");
  TH1F* h_SingleT = (TH1F*)h->Clone("h_SingleT");
  h_SingleT->Reset();
  h_SingleT->Sumw2();
  draw( param, t_SingleT , var, h_SingleT, TCut(cut.c_str()));
  //t_SingleT->Draw(var+">>h_SingleT",   "weight"*TCut(cut.c_str()) );
  h_SingleT->Scale(lumiScale);
  dir->cd();
  h_SingleT->Write("SingleT", TObject::kOverwrite);

  // BBB
  for(int bin = 1; bin<=h_SingleT->GetNbinsX(); bin++ ){
    TH1F* h_SingleT_b_up   = (TH1F*)h->Clone(Form("h_SingleT_%d_Up",bin));
    TH1F* h_SingleT_b_down = (TH1F*)h->Clone(Form("h_SingleT_%d_Down",bin));
    bbb( h_SingleT, h_SingleT_b_up, h_SingleT_b_down, bin);

    dir->cd();
    h_SingleT_b_up  ->Write(Form("SingleT_SingleT%sbin%dUp",  category.Data(), bin), TObject::kOverwrite);
    h_SingleT_b_down->Write(Form("SingleT_SingleT%sbin%dDown",category.Data(), bin), TObject::kOverwrite);
 
  }

  // NOMINAL
  TFile* f_TTV = TFile::Open(inputpath+"MEAnalysis_"+fname+"_nominal"+version+"_TTV.root");
  TTree* t_TTV = (TTree*)f_TTV->Get("tree");
  TH1F* h_TTV = (TH1F*)h->Clone("h_TTV");
  h_TTV->Reset();
  h_TTV->Sumw2();
  draw( param, t_TTV , var, h_TTV, TCut(cut.c_str()));
  //t_TTV->Draw(var+">>h_TTV",   "weight"*TCut(cut.c_str()) );
  h_TTV->Scale(lumiScale);
  dir->cd();
  h_TTV->Write("TTV", TObject::kOverwrite);

  // BBB
  for(int bin = 1; bin<=h_TTV->GetNbinsX(); bin++ ){
    TH1F* h_TTV_b_up   = (TH1F*)h->Clone(Form("h_TTV_%d_Up",bin));
    TH1F* h_TTV_b_down = (TH1F*)h->Clone(Form("h_TTV_%d_Down",bin));
    bbb( h_TTV, h_TTV_b_up, h_TTV_b_down, bin);

    dir->cd();
    h_TTV_b_up  ->Write(Form("TTV_TTV%sbin%dUp",  category.Data(), bin), TObject::kOverwrite);
    h_TTV_b_down->Write(Form("TTV_TTV%sbin%dDown",category.Data(), bin), TObject::kOverwrite);
 
  }

  // NOMINAL
  TFile* f_EWK = TFile::Open(inputpath+"MEAnalysis_"+fname+"_nominal"+version+"_EWK.root");
  TTree* t_EWK = (TTree*)f_EWK->Get("tree");
  TH1F* h_EWK = (TH1F*)h->Clone("h_EWK");
  h_EWK->Reset();
  h_EWK->Sumw2();
  draw( param, t_EWK, var, h_EWK, TCut(cut.c_str()));
  //t_EWK->Draw(var+">>h_EWK",   "weight"*TCut(cut.c_str()) );
  h_EWK->Scale(lumiScale);
  dir->cd();
  h_EWK->Write("EWK", TObject::kOverwrite);

  // BBB
  for(int bin = 1; bin<=h_EWK->GetNbinsX(); bin++ ){
    TH1F* h_EWK_b_up   = (TH1F*)h->Clone(Form("h_EWK_%d_Up",bin));
    TH1F* h_EWK_b_down = (TH1F*)h->Clone(Form("h_EWK_%d_Down",bin));
    bbb( h_EWK, h_EWK_b_up, h_EWK_b_down, bin);

    dir->cd();
    h_EWK_b_up  ->Write(Form("EWK_EWK%sbin%dUp",  category.Data(), bin), TObject::kOverwrite);
    h_EWK_b_down->Write(Form("EWK_EWK%sbin%dDown",category.Data(), bin), TObject::kOverwrite);
 
  }




  cout << "TTH125:      "  << h_TTH125->Integral()         << endl;
  cout << "TTJets (HF): "  << h_TTJetsSemiLept_HF->Integral() << endl;
  cout << "TTJets (LF): "  << h_TTJetsSemiLept_LF->Integral() << endl;
  cout << "SingleT:     "  << h_SingleT->Integral() << endl;
  cout << "TTV:         "  << h_TTV->Integral() << endl;
  cout << "EWK:         "  << h_EWK->Integral() << endl;
  cout << "DiBoson:     "  << h_DiBoson->Integral() << endl;

  float rTTH125    = h_TTH125->Integral() ;
  float rTTJets_HF = h_TTJetsSemiLept_HF->Integral();
  float rTTJets_LF = h_TTJetsSemiLept_LF->Integral();
  float rSingleT   = h_SingleT->Integral() ;
  float rTTV       = h_TTV->Integral() ;
  float rEWK       = h_EWK->Integral() ;
  float rDiBoson   = h_DiBoson->Integral() ;

  float observation = int(rTTH125 + rTTJets_HF + rTTJets_LF  + rTTV + rSingleT /*+ rEWK + rDiBoson*/ );

  TH1F* h_data = (TH1F*)h->Clone("h_data");
  h_data->Reset();
  h_data->Add(h_TTH125,1.0);
  h_data->Add(h_TTJetsSemiLept_HF,1.0);
  h_data->Add(h_TTJetsSemiLept_LF,1.0);
  h_data->Add(h_TTV);
  h_data->Add(h_SingleT);
  //h_data->Add(h_EWK);
  //h_data->Add(h_DiBoson);

  h_data->Scale(int(observation)/h_data->Integral());
  dir->cd();
  h_data->Write("data_obs", TObject::kOverwrite);

  fout->Close();


  ofstream out(Form("%s.txt",(fname+"_"+category).Data()));
  out.precision(8);
  string longspace  = "              ";
  string shortspace = "          ";


  //string null(category.Data());
  string null("");


  string line1("imax 1"); out << line1 << endl;
  string line2("jmax *"); out << line2 << endl;
  string line3("kmax *"); out << line3 << endl;
  string line4(Form("shapes *  *    %s.root  $CHANNEL/$PROCESS $CHANNEL/$PROCESS_$SYSTEMATIC", fname.Data())); out << line4 << endl;
  out<<endl;
  string line5(Form("observation %.0f", observation )); out << line5 << endl;
  out<<endl;


  // PROC AND NORM
  string line("bin                         ");
  if( rTTH125>0 )    line += string(Form("%s    ", (fname+"_"+category).Data()));
  if( rTTJets_HF>0 ) line += string(Form("%s    ", (fname+"_"+category).Data()));
  if( rTTJets_LF>0 ) line += string(Form("%s    ", (fname+"_"+category).Data()));
  if( rTTV>0 )       line += string(Form("%s    ", (fname+"_"+category).Data()));
  if( rSingleT>0 )   line += string(Form("%s    ", (fname+"_"+category).Data()));
  out<<line;
  out<<endl; 
  line = "process                     ";
  if( rTTH125>0 )    line += "TTH125      ";
  if( rTTJets_HF>0 ) line += "TTJetsHF    ";
  if( rTTJets_LF>0 ) line += "TTJetsLF    ";
  if( rTTV>0 )       line += "TTV      ";
  if( rSingleT>0 )   line += "SingleT     ";
  out<<line;
  out<<endl;
  line = "process                       ";
  if( rTTH125>0 )    line += "0          ";
  if( rTTJets_HF>0 ) line += "1          ";
  if( rTTJets_LF>0 ) line += "2          ";
  if( rTTV>0 )       line += "3          ";
  if( rSingleT>0 )   line += "4          ";
  out<<line;
  out<<endl;
  line = "rate                        ";
  if( rTTH125>0 )    line += string(Form("%.4f     ", rTTH125 ));
  if( rTTJets_HF>0 ) line += string(Form("%.4f     ", rTTJets_HF ));
  if( rTTJets_LF>0 ) line += string(Form("%.4f     ", rTTJets_LF));
  if( rTTV>0 )       line += string(Form("%.4f     ", rTTV));
  if( rSingleT>0 )   line += string(Form("%.4f     ", rSingleT));
  out<<line;
  out<<endl;
  out<<endl;
  out<< "-----------------------------------------------------------------" << endl;


  // SYSTEMATICS
  line = "lumi                  lnN  ";
  if( rTTH125>0 )    line += "1.022      ";
  if( rTTJets_HF>0 ) line += "1.022      ";
  if( rTTJets_LF>0 ) line += "1.022      ";
  if( rTTV>0 )       line += "1.022      ";
  if( rSingleT>0 )   line += "1.022      ";
  out<<line;
  out<<endl;
  line = "csv                   shape  ";
  if( rTTH125>0 )    line += "1.0        ";
  if( rTTJets_HF>0 ) line += "1.0        ";
  if( rTTJets_LF>0 ) line += "1.0        ";
  if( rTTV>0 )       line += " -         ";
  if( rSingleT>0 )   line += " -         ";
  out<<line;
  out<<endl;
  line = "JEC                   shape  ";
  if( rTTH125>0 )    line += "1.0        ";
  if( rTTJets_HF>0 ) line += "1.0        ";
  if( rTTJets_LF>0 ) line += "1.0        ";
  if( rTTV>0 )       line += " -         ";
  if( rSingleT>0 )   line += " -         ";
  out<<line;
  out<<endl;

  for(int bin = 1; bin<=h_TTH125->GetNbinsX(); bin++ ){
    if(h_TTH125->GetBinContent(bin)>0){
      line = string(Form("TTH125%sbin%d        shape  ", category.Data(), bin));
      if( rTTH125>0 )    line += "1.0        ";
      if( rTTJets_HF>0 ) line += " -         ";
      if( rTTJets_LF>0 ) line += " -         ";
      if( rTTV>0 )       line += " -         ";
      if( rSingleT>0 )   line += " -         ";
      out<<line;
      out<<endl;
    }
  }  
  for(int bin = 1; bin<=h_TTJetsSemiLept_HF->GetNbinsX(); bin++ ){
    if(h_TTJetsSemiLept_HF->GetBinContent(bin)>0){
      line = string(Form("TTJetsHF%sbin%d      shape  ", category.Data(), bin));
      if( rTTH125>0 )    line += " -         ";
      if( rTTJets_HF>0 ) line += "1.0        ";
      if( rTTJets_LF>0 ) line += " -         ";
      if( rTTV>0 )       line += " -         ";
      if( rSingleT>0 )   line += " -         ";
      out<<line;
      out<<endl;
    }
  }  
  for(int bin = 1; bin<=h_TTJetsSemiLept_LF->GetNbinsX(); bin++ ){
    if(h_TTJetsSemiLept_LF->GetBinContent(bin)>0){
      line = string(Form("TTJetsLF%sbin%d      shape  ", category.Data(), bin));
      if( rTTH125>0 )    line += " -         ";
      if( rTTJets_HF>0 ) line += " -         ";
      if( rTTJets_LF>0 ) line += "1.0        ";
      if( rTTV>0 )       line += " -         ";
      if( rSingleT>0 )   line += " -         ";
      out<<line;
      out<<endl;
    }
  }  
  for(int bin = 1; bin<=h_TTV->GetNbinsX(); bin++ ){
    if(h_TTV->GetBinContent(bin)>0){
      line = string(Form("TTV%sbin%d           shape  ", category.Data(), bin));
      if( rTTH125>0 )    line += " -         ";
      if( rTTJets_HF>0 ) line += " -         ";
      if( rTTJets_LF>0 ) line += " -         ";
      if( rTTV>0 )       line += "1.0        ";
      if( rSingleT>0 )   line += " -         ";
      out<<line;
      out<<endl;
    }
  } 
  for(int bin = 1; bin<=h_SingleT->GetNbinsX(); bin++ ){
    if(h_SingleT->GetBinContent(bin)>0){
      line = string(Form("SingleT%sbin%d       shape  ", category.Data(), bin));
      if( rTTH125>0 )    line += " -         ";
      if( rTTJets_HF>0 ) line += " -         ";
      if( rTTJets_LF>0 ) line += " -         ";
      if( rTTV>0 )       line += " -         ";
      if( rSingleT>0 )   line += "1.0        ";
      out<<line;
      out<<endl;
    }
  } 

  line = "Norm_TTbb             lnN    ";
  if( rTTH125>0 )    line += " -         ";
  if( rTTJets_HF>0 ) line += "1.50       ";
  if( rTTJets_LF>0 ) line += " -         ";
  if( rTTV>0 )       line += " -         ";
  if( rSingleT>0 )   line += " -         ";
  out<<line;
  out<<endl;

  line = "Norm_TTV              lnN    ";
  if( rTTH125>0 )    line += " -         ";
  if( rTTJets_HF>0 ) line += " -         ";
  if( rTTJets_LF>0 ) line += " -         ";
  if( rTTV>0 )       line += "1.20       ";
  if( rSingleT>0 )   line += " -         ";
  out<<line;
  out<<endl;

  line = "Norm_SingleT          lnN    ";
  if( rTTH125>0 )    line += " -         ";
  if( rTTJets_HF>0 ) line += " -         ";
  if( rTTJets_LF>0 ) line += " -         ";
  if( rTTV>0 )       line += " -         ";
  if( rSingleT>0 )   line += "1.20       ";
  out<<line;
  out<<endl;

  line = "QCDscale_TTH          lnN    ";
  if( rTTH125>0 )    line += "1.12       ";
  if( rTTJets_HF>0 ) line += " -         ";
  if( rTTJets_LF>0 ) line += " -         ";
  if( rTTV>0 )       line += " -         ";
  if( rSingleT>0 )   line += " -         ";
  out<<line;
  out<<endl;

  line = string(Form("QCDscale%s_TTJetsHF     lnN    ", null.c_str()));
  if( rTTH125>0 )    line += " -         ";
  if( rTTJets_HF>0 ) line += "1.35       ";
  if( rTTJets_LF>0 ) line += " -         ";
  if( rTTV>0 )       line += " -         ";
  if( rSingleT>0 )   line += " -         ";
  out<<line;
  out<<endl;

  line = string(Form("QCDscale%s_TTJetsLF     lnN    ", null.c_str()));
  if( rTTH125>0 )    line += " -         ";
  if( rTTJets_HF>0 ) line += " -         ";
  if( rTTJets_LF>0 ) line += "1.35       ";
  if( rTTV>0 )       line += " -         ";
  if( rSingleT>0 )   line += " -         ";
  out<<line;
  out<<endl;

  line = "pdf_gg                lnN    ";
  if( rTTH125>0 )    line += "1.03       ";
  if( rTTJets_HF>0 ) line += "1.03       ";
  if( rTTJets_LF>0 ) line += "1.03       ";
  if( rTTV>0 )       line += "1.03       ";
  if( rSingleT>0 )   line += "1.03       ";
  out<<line;
  out<<endl;
  out<<endl;
  out<< "-----------------------------------------------------------------" << endl;

}



void produceAll( float LumiScale = 19.5/12.1){

  //produce("SL", "type==0",  "cat1",                    0.02, 0.50, LumiScale ,   1, 7);
  //produce("SL", "type==1",  "cat2",                    0.02, 0.7,  LumiScale ,   1, 7);
  //produce("SL", "type==2 && flag_type2>0",  "cat3",    0.02, 0.8,  LumiScale ,   1, 7);
  //produce("SL", "type==2 && flag_type2<0",  "cat4",    0.02, 0.7,  LumiScale ,   1, 7);
  //produce("SL", "type==3",  "cat5",                    0.02, 0.5,  LumiScale ,   1, 7);  
  //produce("DL", "type==6",  "cat6",                    0.02, 0.,   LumiScale*2 , 0, 6);
  //produce("DL", "type==7",  "cat7",                    0.02, 0.,   LumiScale*2 , 0, 6);

  produce("SL", "type==0",  "cat1",                       0.02, 0.50, LumiScale ,   1, 4, 1);
  //produce("SL", "type==0 && flag_type0<=1",  "cat1a",     0.02, 0.50, LumiScale ,   1, 4, 1);
  //produce("SL", "type==0 && flag_type0>1",   "cat1b",     0.02, 0.50, LumiScale ,   1, 4, 1);
  produce("SL", "type==1",  "cat2",                       0.02, 0.7,  LumiScale ,   1, 5, 1);
  produce("SL", "type==2 && flag_type2>0",   "cat3",      0.02, 0.8,  LumiScale ,   1, 5, 1);
  produce("SL", "type==2 && flag_type2<=0",  "cat4",      0.02, 0.7,  LumiScale ,   1, 5, 1);
  produce("SL", "type==3",  "cat5",                       0.02, 0.5,  LumiScale ,   1, 6, 1);  
  produce("DL", "type==6",  "cat6",                       0.02, 0.,   LumiScale*2 , 0, 5, 0);
  produce("DL", "type==7",  "cat7",                       0.02, 0.,   LumiScale*2 , 0, 6, 0);


}







void produce_2( TString fname = "SL_VType2", string cut = "type==0", TString category = "cat0", float fact1 = 0.02, float fact2 = 0., float lumiScale = 20./12.1,
		int newvar = 1, int nBins=6, int splitFirstBin=0){


  string version = (string(fname.Data())).find("SL")!=string::npos ? "_v6" : "_v5";
  //string version = "";

  string  versionTTJets  = (string(fname.Data())).find("SL")!=string::npos ? "_v7_" : "_v7_";

  cout << "Doing version " << version << " and category " << category << endl;


  // input files path
  TString inputpath = "gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/Trees/MEM/Oct04_2013/";


  float scaleTTH      = 1.0;//4.7;
  float scaleTTJetsLF = 1.0;//3.9;
  float scaleTTJetsHF = 1.0;//3.9;

  cout << "*****************************************************" << endl;
  cout << "Attention!!! The following rescaling will be applied:" << endl;
  cout << "TTH:       " << scaleTTH << endl;
  cout << "TTJets LF: " << scaleTTJetsLF << endl;
  cout << "TTJets HF: " << scaleTTJetsHF << endl;
  cout << "*****************************************************" << endl;

  vector<float> param;
  param.clear();
  if(splitFirstBin){
    param.push_back(1.0);
    param.push_back(0.5);
  }

  TFile* fout = new TFile(fname+"_2.root","UPDATE");
  TDirectory* dir =  fout->GetDirectory( fname+"_"+category); 
  if( !dir) dir = fout->mkdir( fname+"_"+category ) ;

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


  TH1F* h = new TH1F("h","Simulation #sqrt{s}=8 TeV, "+fname+"; S/(S+B); units",  param.size()==0 ? nBins : nBins+1 ,  param.size()==0 ? bins.GetArray() : bins2.GetArray());
 
  TString var("");
  if(newvar)
    var = TString(Form("p_125_all_s_ttbb/(p_125_all_s_ttbb+%f*((1-%f)*p_125_all_b_ttbb+%f*p_125_all_b_ttjj))", fact1, fact2, fact2));
  else
    var = TString(Form("p_125_all_s/(p_125_all_s+p_125_all_b*%f)", fact1));

  string basecut = cut;
  TString ttjets = (string(fname.Data())).find("SL")!=string::npos ? "TTJets" : "TTJets";

  ////////////////////////////// TTH
  
  // NOMINAL
  TFile* f_TTH125 = TFile::Open(inputpath+"MEAnalysis_"+fname+"_nominal"+version+"_TTH125.root");
  TTree* t_TTH125 = (TTree*)f_TTH125->Get("tree");
  TH1F* h_TTH125 = (TH1F*)h->Clone("h_TTH125");
  h_TTH125->Reset();
  h_TTH125->Sumw2();
  draw( param, t_TTH125, var, h_TTH125, TCut(cut.c_str()));
  //t_TTH125->Draw(var+">>h_TTH125",   "weight"*TCut(cut.c_str()) );
  h_TTH125->Scale(lumiScale*scaleTTH);
  dir->cd();
  h_TTH125->Write("TTH125", TObject::kOverwrite);

  // BBB
  for(int bin = 1; bin<=h_TTH125->GetNbinsX(); bin++ ){
    TH1F* h_TTH125_b_up   = (TH1F*)h->Clone(Form("h_TTH125_%d_Up",bin));
    TH1F* h_TTH125_b_down = (TH1F*)h->Clone(Form("h_TTH125_%d_Down",bin));
    bbb( h_TTH125, h_TTH125_b_up, h_TTH125_b_down, bin);

    dir->cd();
    h_TTH125_b_up  ->Write(Form("TTH125_TTH125%sbin%dUp",  category.Data(), bin), TObject::kOverwrite);
    h_TTH125_b_down->Write(Form("TTH125_TTH125%sbin%dDown",category.Data(), bin), TObject::kOverwrite);
 
  }

  // JEC UP
  TFile* f_TTH125_JECUp = TFile::Open(inputpath+"MEAnalysis_"+fname+"_JECUp"+version+"_TTH125.root");
  TTree* t_TTH125_JECUp = (TTree*)f_TTH125_JECUp->Get("tree");
  TH1F* h_TTH125_JECUp = (TH1F*)h->Clone("h_TTH125_JECUp");
  h_TTH125_JECUp->Reset();
  h_TTH125_JECUp->Sumw2();
  draw( param, t_TTH125_JECUp, var, h_TTH125_JECUp, TCut(cut.c_str()));
  //t_TTH125_JECUp->Draw(var+">>h_TTH125_JECUp",   "weight"*TCut(cut.c_str()) );
  h_TTH125_JECUp->Scale(lumiScale*scaleTTH);
  dir->cd();
  h_TTH125_JECUp->Write("TTH125_JECUp", TObject::kOverwrite);

  // JEC DOWN
  TFile* f_TTH125_JECDown = TFile::Open(inputpath+"MEAnalysis_"+fname+"_JECDown"+version+"_TTH125.root");
  TTree* t_TTH125_JECDown = (TTree*)f_TTH125_JECDown->Get("tree");
  TH1F* h_TTH125_JECDown = (TH1F*)h->Clone("h_TTH125_JECDown");
  h_TTH125_JECDown->Reset();
  h_TTH125_JECDown->Sumw2();
  draw( param, t_TTH125_JECDown, var, h_TTH125_JECDown, TCut(cut.c_str()));
  //t_TTH125_JECDown->Draw(var+">>h_TTH125_JECDown",   "weight"*TCut(cut.c_str()) );
  h_TTH125_JECDown->Scale(lumiScale*scaleTTH);
  dir->cd();
  h_TTH125_JECDown->Write("TTH125_JECDown", TObject::kOverwrite);

  // csv UP
  TFile* f_TTH125_csvUp = TFile::Open(inputpath+"MEAnalysis_"+fname+"_csvUp"+version+"_TTH125.root");
  TTree* t_TTH125_csvUp = (TTree*)f_TTH125_csvUp->Get("tree");
  TH1F* h_TTH125_csvUp = (TH1F*)h->Clone("h_TTH125_csvUp");
  h_TTH125_csvUp->Reset();
  h_TTH125_csvUp->Sumw2();
  draw( param, t_TTH125_csvUp, var, h_TTH125_csvUp, TCut(cut.c_str()));
  //t_TTH125_csvUp->Draw(var+">>h_TTH125_csvUp",   "weight"*TCut(cut.c_str()) );
  h_TTH125_csvUp->Scale(lumiScale*scaleTTH);
  dir->cd();
  h_TTH125_csvUp->Write("TTH125_csvUp", TObject::kOverwrite);

  // csv DOWN
  TFile* f_TTH125_csvDown = TFile::Open(inputpath+"MEAnalysis_"+fname+"_csvDown"+version+"_TTH125.root");
  TTree* t_TTH125_csvDown = (TTree*)f_TTH125_csvDown->Get("tree");
  TH1F* h_TTH125_csvDown = (TH1F*)h->Clone("h_TTH125_csvDown");
  h_TTH125_csvDown->Reset();
  h_TTH125_csvDown->Sumw2();
  draw( param, t_TTH125_csvDown, var, h_TTH125_csvDown, TCut(cut.c_str()));
  //t_TTH125_csvDown->Draw(var+">>h_TTH125_csvDown",   "weight"*TCut(cut.c_str()) );
  h_TTH125_csvDown->Scale(lumiScale*scaleTTH);
  dir->cd();
  h_TTH125_csvDown->Write("TTH125_csvDown", TObject::kOverwrite);




  ////////////////////////////// TTJets: HFbb


  cut = basecut+"&&nSimBs>2 && nMatchSimBs>=2";

  // NOMINAL"+VERSION+"
  TFile* f_TTJetsSemiLept_HFbb = TFile::Open(inputpath+"MEAnalysis_"+fname+"_nominal"+versionTTJets+ttjets+".root");
  TTree* t_TTJetsSemiLept_HFbb = (TTree*)f_TTJetsSemiLept_HFbb->Get("tree");
  TH1F* h_TTJetsSemiLept_HFbb = (TH1F*)h->Clone("h_TTJetsSemiLept_HFbb");
  h_TTJetsSemiLept_HFbb->Reset();
  h_TTJetsSemiLept_HFbb->Sumw2();
  draw( param, t_TTJetsSemiLept_HFbb, var, h_TTJetsSemiLept_HFbb, TCut(cut.c_str()));
  //t_TTJetsSemiLept_HFbb->Draw(var+">>h_TTJetsSemiLept_HFbb",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_HFbb->Scale(lumiScale*scaleTTJetsHF);
  dir->cd();
  h_TTJetsSemiLept_HFbb->Write("TTJetsHFbb", TObject::kOverwrite);

  // BBB
  for(int bin = 1; bin<=h_TTJetsSemiLept_HFbb->GetNbinsX(); bin++ ){
    TH1F* h_TTJetsSemiLept_HFbb_b_up   = (TH1F*)h->Clone(Form("h_TTJetsSemiLept_HFbb_%d_Up",bin));
    TH1F* h_TTJetsSemiLept_HFbb_b_down = (TH1F*)h->Clone(Form("h_TTJetsSemiLept_HFbb_%d_Down",bin));
    bbb( h_TTJetsSemiLept_HFbb, h_TTJetsSemiLept_HFbb_b_up, h_TTJetsSemiLept_HFbb_b_down, bin);

    dir->cd();
    h_TTJetsSemiLept_HFbb_b_up  ->Write(Form("TTJetsHFbb_TTJetsHFbb%sbin%dUp",    category.Data(),  bin), TObject::kOverwrite);
    h_TTJetsSemiLept_HFbb_b_down->Write(Form("TTJetsHFbb_TTJetsHFbb%sbin%dDown",  category.Data(),bin), TObject::kOverwrite);
 
  }

  // JEC UP
  TFile* f_TTJetsSemiLept_HFbb_JECUp = TFile::Open(inputpath+"MEAnalysis_"+fname+"_JECUp"+versionTTJets+ttjets+".root");
  TTree* t_TTJetsSemiLept_HFbb_JECUp = (TTree*)f_TTJetsSemiLept_HFbb_JECUp->Get("tree");
  TH1F* h_TTJetsSemiLept_HFbb_JECUp = (TH1F*)h->Clone("h_TTJetsSemiLept_HFbb_JECUp");
  h_TTJetsSemiLept_HFbb_JECUp->Reset();
  h_TTJetsSemiLept_HFbb_JECUp->Sumw2();
  draw( param, t_TTJetsSemiLept_HFbb_JECUp , var, h_TTJetsSemiLept_HFbb_JECUp, TCut(cut.c_str()));
  //t_TTJetsSemiLept_HFbb_JECUp->Draw(var+">>h_TTJetsSemiLept_HFbb_JECUp",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_HFbb_JECUp->Scale(lumiScale*scaleTTJetsHF);
  dir->cd();
  h_TTJetsSemiLept_HFbb_JECUp->Write("TTJetsHFbb_JECUp", TObject::kOverwrite);

  // JEC DOWN
  TFile* f_TTJetsSemiLept_HFbb_JECDown = TFile::Open(inputpath+"MEAnalysis_"+fname+"_JECDown"+versionTTJets+ttjets+".root");
  TTree* t_TTJetsSemiLept_HFbb_JECDown = (TTree*)f_TTJetsSemiLept_HFbb_JECDown->Get("tree");
  TH1F* h_TTJetsSemiLept_HFbb_JECDown = (TH1F*)h->Clone("h_TTJetsSemiLept_HFbb_JECDown");
  h_TTJetsSemiLept_HFbb_JECDown->Reset();
  h_TTJetsSemiLept_HFbb_JECDown->Sumw2();
  draw( param, t_TTJetsSemiLept_HFbb_JECDown, var, h_TTJetsSemiLept_HFbb_JECDown, TCut(cut.c_str()));
  //t_TTJetsSemiLept_HFbb_JECDown->Draw(var+">>h_TTJetsSemiLept_HFbb_JECDown",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_HFbb_JECDown->Scale(lumiScale*scaleTTJetsHF);
  dir->cd();
  h_TTJetsSemiLept_HFbb_JECDown->Write("TTJetsHFbb_JECDown", TObject::kOverwrite);

  // csv UP
  TFile* f_TTJetsSemiLept_HFbb_csvUp = TFile::Open(inputpath+"MEAnalysis_"+fname+"_csvUp"+versionTTJets+ttjets+".root");
  TTree* t_TTJetsSemiLept_HFbb_csvUp = (TTree*)f_TTJetsSemiLept_HFbb_csvUp->Get("tree");
  TH1F* h_TTJetsSemiLept_HFbb_csvUp = (TH1F*)h->Clone("h_TTJetsSemiLept_HFbb_csvUp");
  h_TTJetsSemiLept_HFbb_csvUp->Reset();
  h_TTJetsSemiLept_HFbb_csvUp->Sumw2();
  draw( param, t_TTJetsSemiLept_HFbb_csvUp, var, h_TTJetsSemiLept_HFbb_csvUp, TCut(cut.c_str()));
  //t_TTJetsSemiLept_HFbb_csvUp->Draw(var+">>h_TTJetsSemiLept_HFbb_csvUp",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_HFbb_csvUp->Scale(lumiScale*scaleTTJetsHF);
  dir->cd();
  h_TTJetsSemiLept_HFbb_csvUp->Write("TTJetsHFbb_csvUp", TObject::kOverwrite);

  // csv DOWN
  TFile* f_TTJetsSemiLept_HFbb_csvDown = TFile::Open(inputpath+"MEAnalysis_"+fname+"_csvDown"+versionTTJets+ttjets+".root");
  TTree* t_TTJetsSemiLept_HFbb_csvDown = (TTree*)f_TTJetsSemiLept_HFbb_csvDown->Get("tree");
  TH1F* h_TTJetsSemiLept_HFbb_csvDown = (TH1F*)h->Clone("h_TTJetsSemiLept_HFbb_csvDown");
  h_TTJetsSemiLept_HFbb_csvDown->Reset();
  h_TTJetsSemiLept_HFbb_csvDown->Sumw2();
  draw( param, t_TTJetsSemiLept_HFbb_csvDown, var, h_TTJetsSemiLept_HFbb_csvDown, TCut(cut.c_str()));
  //t_TTJetsSemiLept_HFbb_csvDown->Draw(var+">>h_TTJetsSemiLept_HFbb_csvDown",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_HFbb_csvDown->Scale(lumiScale*scaleTTJetsHF);
  dir->cd();
  h_TTJetsSemiLept_HFbb_csvDown->Write("TTJetsHFbb_csvDown", TObject::kOverwrite);




  ////////////////////////////// TTJets: HFb


  cut = basecut+"&&nSimBs>2 && nMatchSimBs<2";

  // NOMINAL"+VERSION+"
  TFile* f_TTJetsSemiLept_HFb = TFile::Open(inputpath+"MEAnalysis_"+fname+"_nominal"+versionTTJets+ttjets+".root");
  TTree* t_TTJetsSemiLept_HFb = (TTree*)f_TTJetsSemiLept_HFb->Get("tree");
  TH1F* h_TTJetsSemiLept_HFb = (TH1F*)h->Clone("h_TTJetsSemiLept_HFb");
  h_TTJetsSemiLept_HFb->Reset();
  h_TTJetsSemiLept_HFb->Sumw2();
  draw( param, t_TTJetsSemiLept_HFb, var, h_TTJetsSemiLept_HFb, TCut(cut.c_str()));
  //t_TTJetsSemiLept_HFb->Draw(var+">>h_TTJetsSemiLept_HFb",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_HFb->Scale(lumiScale*scaleTTJetsHF);
  dir->cd();
  h_TTJetsSemiLept_HFb->Write("TTJetsHFb", TObject::kOverwrite);

  // BBB
  for(int bin = 1; bin<=h_TTJetsSemiLept_HFb->GetNbinsX(); bin++ ){
    TH1F* h_TTJetsSemiLept_HFb_b_up   = (TH1F*)h->Clone(Form("h_TTJetsSemiLept_HFb_%d_Up",bin));
    TH1F* h_TTJetsSemiLept_HFb_b_down = (TH1F*)h->Clone(Form("h_TTJetsSemiLept_HFb_%d_Down",bin));
    bbb( h_TTJetsSemiLept_HFb, h_TTJetsSemiLept_HFb_b_up, h_TTJetsSemiLept_HFb_b_down, bin);

    dir->cd();
    h_TTJetsSemiLept_HFb_b_up  ->Write(Form("TTJetsHFb_TTJetsHFb%sbin%dUp",    category.Data(),  bin), TObject::kOverwrite);
    h_TTJetsSemiLept_HFb_b_down->Write(Form("TTJetsHFb_TTJetsHFb%sbin%dDown",  category.Data(),bin), TObject::kOverwrite);
 
  }

  // JEC UP
  TFile* f_TTJetsSemiLept_HFb_JECUp = TFile::Open(inputpath+"MEAnalysis_"+fname+"_JECUp"+versionTTJets+ttjets+".root");
  TTree* t_TTJetsSemiLept_HFb_JECUp = (TTree*)f_TTJetsSemiLept_HFb_JECUp->Get("tree");
  TH1F* h_TTJetsSemiLept_HFb_JECUp = (TH1F*)h->Clone("h_TTJetsSemiLept_HFb_JECUp");
  h_TTJetsSemiLept_HFb_JECUp->Reset();
  h_TTJetsSemiLept_HFb_JECUp->Sumw2();
  draw( param, t_TTJetsSemiLept_HFb_JECUp , var, h_TTJetsSemiLept_HFb_JECUp, TCut(cut.c_str()));
  //t_TTJetsSemiLept_HFb_JECUp->Draw(var+">>h_TTJetsSemiLept_HFb_JECUp",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_HFb_JECUp->Scale(lumiScale*scaleTTJetsHF);
  dir->cd();
  h_TTJetsSemiLept_HFb_JECUp->Write("TTJetsHFb_JECUp", TObject::kOverwrite);

  // JEC DOWN
  TFile* f_TTJetsSemiLept_HFb_JECDown = TFile::Open(inputpath+"MEAnalysis_"+fname+"_JECDown"+versionTTJets+ttjets+".root");
  TTree* t_TTJetsSemiLept_HFb_JECDown = (TTree*)f_TTJetsSemiLept_HFb_JECDown->Get("tree");
  TH1F* h_TTJetsSemiLept_HFb_JECDown = (TH1F*)h->Clone("h_TTJetsSemiLept_HFb_JECDown");
  h_TTJetsSemiLept_HFb_JECDown->Reset();
  h_TTJetsSemiLept_HFb_JECDown->Sumw2();
  draw( param, t_TTJetsSemiLept_HFb_JECDown, var, h_TTJetsSemiLept_HFb_JECDown, TCut(cut.c_str()));
  //t_TTJetsSemiLept_HFb_JECDown->Draw(var+">>h_TTJetsSemiLept_HFb_JECDown",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_HFb_JECDown->Scale(lumiScale*scaleTTJetsHF);
  dir->cd();
  h_TTJetsSemiLept_HFb_JECDown->Write("TTJetsHFb_JECDown", TObject::kOverwrite);

  // csv UP
  TFile* f_TTJetsSemiLept_HFb_csvUp = TFile::Open(inputpath+"MEAnalysis_"+fname+"_csvUp"+versionTTJets+ttjets+".root");
  TTree* t_TTJetsSemiLept_HFb_csvUp = (TTree*)f_TTJetsSemiLept_HFb_csvUp->Get("tree");
  TH1F* h_TTJetsSemiLept_HFb_csvUp = (TH1F*)h->Clone("h_TTJetsSemiLept_HFb_csvUp");
  h_TTJetsSemiLept_HFb_csvUp->Reset();
  h_TTJetsSemiLept_HFb_csvUp->Sumw2();
  draw( param, t_TTJetsSemiLept_HFb_csvUp, var, h_TTJetsSemiLept_HFb_csvUp, TCut(cut.c_str()));
  //t_TTJetsSemiLept_HFb_csvUp->Draw(var+">>h_TTJetsSemiLept_HFb_csvUp",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_HFb_csvUp->Scale(lumiScale*scaleTTJetsHF);
  dir->cd();
  h_TTJetsSemiLept_HFb_csvUp->Write("TTJetsHFb_csvUp", TObject::kOverwrite);

  // csv DOWN
  TFile* f_TTJetsSemiLept_HFb_csvDown = TFile::Open(inputpath+"MEAnalysis_"+fname+"_csvDown"+versionTTJets+ttjets+".root");
  TTree* t_TTJetsSemiLept_HFb_csvDown = (TTree*)f_TTJetsSemiLept_HFb_csvDown->Get("tree");
  TH1F* h_TTJetsSemiLept_HFb_csvDown = (TH1F*)h->Clone("h_TTJetsSemiLept_HFb_csvDown");
  h_TTJetsSemiLept_HFb_csvDown->Reset();
  h_TTJetsSemiLept_HFb_csvDown->Sumw2();
  draw( param, t_TTJetsSemiLept_HFb_csvDown, var, h_TTJetsSemiLept_HFb_csvDown, TCut(cut.c_str()));
  //t_TTJetsSemiLept_HFb_csvDown->Draw(var+">>h_TTJetsSemiLept_HFb_csvDown",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_HFb_csvDown->Scale(lumiScale*scaleTTJetsHF);
  dir->cd();
  h_TTJetsSemiLept_HFb_csvDown->Write("TTJetsHFb_csvDown", TObject::kOverwrite);



  ////////////////////////////// TTJets: LF

  cut = basecut+"&&nSimBs==2";

  // NOMINAL"+VERSION+"
  TFile* f_TTJetsSemiLept_LF = TFile::Open(inputpath+"MEAnalysis_"+fname+"_nominal"+versionTTJets+ttjets+".root");
  TTree* t_TTJetsSemiLept_LF = (TTree*)f_TTJetsSemiLept_LF->Get("tree");
  TH1F* h_TTJetsSemiLept_LF = (TH1F*)h->Clone("h_TTJetsSemiLept_LF");
  h_TTJetsSemiLept_LF->Reset();
  h_TTJetsSemiLept_LF->Sumw2();
  draw( param, t_TTJetsSemiLept_LF, var, h_TTJetsSemiLept_LF, TCut(cut.c_str()));
  //t_TTJetsSemiLept_LF->Draw(var+">>h_TTJetsSemiLept_LF",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_LF->Scale(lumiScale*scaleTTJetsLF);
  dir->cd();
  h_TTJetsSemiLept_LF->Write("TTJetsLF", TObject::kOverwrite);

  // BBB
  for(int bin = 1; bin<=h_TTJetsSemiLept_LF->GetNbinsX(); bin++ ){
    TH1F* h_TTJetsSemiLept_LF_b_up   = (TH1F*)h->Clone(Form("h_TTJetsSemiLept_LF_%d_Up",bin));
    TH1F* h_TTJetsSemiLept_LF_b_down = (TH1F*)h->Clone(Form("h_TTJetsSemiLept_LF_%d_Down",bin));
    bbb( h_TTJetsSemiLept_LF, h_TTJetsSemiLept_LF_b_up, h_TTJetsSemiLept_LF_b_down, bin);

    dir->cd();
    h_TTJetsSemiLept_LF_b_up  ->Write(Form("TTJetsLF_TTJetsLF%sbin%dUp",  category.Data(),  bin), TObject::kOverwrite);
    h_TTJetsSemiLept_LF_b_down->Write(Form("TTJetsLF_TTJetsLF%sbin%dDown",category.Data(),bin), TObject::kOverwrite);
 
  }

  // JEC UP
  TFile* f_TTJetsSemiLept_LF_JECUp = TFile::Open(inputpath+"MEAnalysis_"+fname+"_JECUp"+versionTTJets+ttjets+".root");
  TTree* t_TTJetsSemiLept_LF_JECUp = (TTree*)f_TTJetsSemiLept_LF_JECUp->Get("tree");
  TH1F* h_TTJetsSemiLept_LF_JECUp = (TH1F*)h->Clone("h_TTJetsSemiLept_LF_JECUp");
  h_TTJetsSemiLept_LF_JECUp->Reset();
  h_TTJetsSemiLept_LF_JECUp->Sumw2();
  draw( param, t_TTJetsSemiLept_LF_JECUp, var, h_TTJetsSemiLept_LF_JECUp, TCut(cut.c_str()));
  //t_TTJetsSemiLept_LF_JECUp->Draw(var+">>h_TTJetsSemiLept_LF_JECUp",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_LF_JECUp->Scale(lumiScale*scaleTTJetsLF);
  dir->cd();
  h_TTJetsSemiLept_LF_JECUp->Write("TTJetsLF_JECUp", TObject::kOverwrite);

  // JEC DOWN
  TFile* f_TTJetsSemiLept_LF_JECDown = TFile::Open(inputpath+"MEAnalysis_"+fname+"_JECDown"+versionTTJets+ttjets+".root");
  TTree* t_TTJetsSemiLept_LF_JECDown = (TTree*)f_TTJetsSemiLept_LF_JECDown->Get("tree");
  TH1F* h_TTJetsSemiLept_LF_JECDown = (TH1F*)h->Clone("h_TTJetsSemiLept_LF_JECDown");
  h_TTJetsSemiLept_LF_JECDown->Reset();
  h_TTJetsSemiLept_LF_JECDown->Sumw2();
  draw( param, t_TTJetsSemiLept_LF_JECDown , var, h_TTJetsSemiLept_LF_JECDown, TCut(cut.c_str()));
  //t_TTJetsSemiLept_LF_JECDown->Draw(var+">>h_TTJetsSemiLept_LF_JECDown",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_LF_JECDown->Scale(lumiScale*scaleTTJetsLF);
  dir->cd();
  h_TTJetsSemiLept_LF_JECDown->Write("TTJetsLF_JECDown", TObject::kOverwrite);

  // csv UP
  TFile* f_TTJetsSemiLept_LF_csvUp = TFile::Open(inputpath+"MEAnalysis_"+fname+"_csvUp"+versionTTJets+ttjets+".root");
  TTree* t_TTJetsSemiLept_LF_csvUp = (TTree*)f_TTJetsSemiLept_LF_csvUp->Get("tree");
  TH1F* h_TTJetsSemiLept_LF_csvUp = (TH1F*)h->Clone("h_TTJetsSemiLept_LF_csvUp");
  h_TTJetsSemiLept_LF_csvUp->Reset();
  h_TTJetsSemiLept_LF_csvUp->Sumw2();
  draw( param, t_TTJetsSemiLept_LF_csvUp, var, h_TTJetsSemiLept_LF_csvUp, TCut(cut.c_str()));
  //t_TTJetsSemiLept_LF_csvUp->Draw(var+">>h_TTJetsSemiLept_LF_csvUp",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_LF_csvUp->Scale(lumiScale*scaleTTJetsLF);
  dir->cd();
  h_TTJetsSemiLept_LF_csvUp->Write("TTJetsLF_csvUp", TObject::kOverwrite);

  // csv DOWN
  TFile* f_TTJetsSemiLept_LF_csvDown = TFile::Open(inputpath+"MEAnalysis_"+fname+"_csvDown"+versionTTJets+ttjets+".root");
  TTree* t_TTJetsSemiLept_LF_csvDown = (TTree*)f_TTJetsSemiLept_LF_csvDown->Get("tree");
  TH1F* h_TTJetsSemiLept_LF_csvDown = (TH1F*)h->Clone("h_TTJetsSemiLept_LF_csvDown");
  h_TTJetsSemiLept_LF_csvDown->Reset();
  h_TTJetsSemiLept_LF_csvDown->Sumw2();
  draw( param, t_TTJetsSemiLept_LF_csvDown, var, h_TTJetsSemiLept_LF_csvDown, TCut(cut.c_str()));
  //t_TTJetsSemiLept_LF_csvDown->Draw(var+">>h_TTJetsSemiLept_LF_csvDown",   "weight"*TCut(cut.c_str()) );
  h_TTJetsSemiLept_LF_csvDown->Scale(lumiScale*scaleTTJetsLF);
  dir->cd();
  h_TTJetsSemiLept_LF_csvDown->Write("TTJetsLF_csvDown", TObject::kOverwrite);



  ////////////////////////////// DIBOSON
  cut = basecut;
  
  // NOMINAL
  TFile* f_DiBoson = TFile::Open(inputpath+"MEAnalysis_"+fname+"_nominal"+version+"_DiBoson.root");
  TTree* t_DiBoson = (TTree*)f_DiBoson->Get("tree");
  TH1F* h_DiBoson = (TH1F*)h->Clone("h_DiBoson");
  h_DiBoson->Reset();
  h_DiBoson->Sumw2();
  draw( param, t_DiBoson , var, h_DiBoson, TCut(cut.c_str()));
  //t_DiBoson->Draw(var+">>h_DiBoson",   "weight"*TCut(cut.c_str()) );
  h_DiBoson->Scale(lumiScale);
  dir->cd();
  h_DiBoson->Write("DiBoson", TObject::kOverwrite);

  // BBB
  for(int bin = 1; bin<=h_DiBoson->GetNbinsX(); bin++ ){
    TH1F* h_DiBoson_b_up   = (TH1F*)h->Clone(Form("h_DiBoson_%d_Up",bin));
    TH1F* h_DiBoson_b_down = (TH1F*)h->Clone(Form("h_DiBoson_%d_Down",bin));
    bbb( h_DiBoson, h_DiBoson_b_up, h_DiBoson_b_down, bin);

    dir->cd();
    h_DiBoson_b_up  ->Write(Form("DiBoson_DiBoson%sbin%dUp",  category.Data(), bin), TObject::kOverwrite);
    h_DiBoson_b_down->Write(Form("DiBoson_DiBoson%sbin%dDown",category.Data(), bin), TObject::kOverwrite);
 
  }


  ////////////////////////////// SINGLE T
  
  // NOMINAL
  TFile* f_SingleT = TFile::Open(inputpath+"MEAnalysis_"+fname+"_nominal"+version+"_SingleT.root");
  TTree* t_SingleT = (TTree*)f_SingleT->Get("tree");
  TH1F* h_SingleT = (TH1F*)h->Clone("h_SingleT");
  h_SingleT->Reset();
  h_SingleT->Sumw2();
  draw( param, t_SingleT , var, h_SingleT, TCut(cut.c_str()));
  //t_SingleT->Draw(var+">>h_SingleT",   "weight"*TCut(cut.c_str()) );
  h_SingleT->Scale(lumiScale);
  dir->cd();
  h_SingleT->Write("SingleT", TObject::kOverwrite);

  // BBB
  for(int bin = 1; bin<=h_SingleT->GetNbinsX(); bin++ ){
    TH1F* h_SingleT_b_up   = (TH1F*)h->Clone(Form("h_SingleT_%d_Up",bin));
    TH1F* h_SingleT_b_down = (TH1F*)h->Clone(Form("h_SingleT_%d_Down",bin));
    bbb( h_SingleT, h_SingleT_b_up, h_SingleT_b_down, bin);

    dir->cd();
    h_SingleT_b_up  ->Write(Form("SingleT_SingleT%sbin%dUp",  category.Data(), bin), TObject::kOverwrite);
    h_SingleT_b_down->Write(Form("SingleT_SingleT%sbin%dDown",category.Data(), bin), TObject::kOverwrite);
 
  }

  // NOMINAL
  TFile* f_TTV = TFile::Open(inputpath+"MEAnalysis_"+fname+"_nominal"+version+"_TTV.root");
  TTree* t_TTV = (TTree*)f_TTV->Get("tree");
  TH1F* h_TTV = (TH1F*)h->Clone("h_TTV");
  h_TTV->Reset();
  h_TTV->Sumw2();
  draw( param, t_TTV , var, h_TTV, TCut(cut.c_str()));
  //t_TTV->Draw(var+">>h_TTV",   "weight"*TCut(cut.c_str()) );
  h_TTV->Scale(lumiScale);
  dir->cd();
  h_TTV->Write("TTV", TObject::kOverwrite);

  // BBB
  for(int bin = 1; bin<=h_TTV->GetNbinsX(); bin++ ){
    TH1F* h_TTV_b_up   = (TH1F*)h->Clone(Form("h_TTV_%d_Up",bin));
    TH1F* h_TTV_b_down = (TH1F*)h->Clone(Form("h_TTV_%d_Down",bin));
    bbb( h_TTV, h_TTV_b_up, h_TTV_b_down, bin);

    dir->cd();
    h_TTV_b_up  ->Write(Form("TTV_TTV%sbin%dUp",  category.Data(), bin), TObject::kOverwrite);
    h_TTV_b_down->Write(Form("TTV_TTV%sbin%dDown",category.Data(), bin), TObject::kOverwrite);
 
  }

  // NOMINAL
  TFile* f_EWK = TFile::Open(inputpath+"MEAnalysis_"+fname+"_nominal"+version+"_EWK.root");
  TTree* t_EWK = (TTree*)f_EWK->Get("tree");
  TH1F* h_EWK = (TH1F*)h->Clone("h_EWK");
  h_EWK->Reset();
  h_EWK->Sumw2();
  draw( param, t_EWK, var, h_EWK, TCut(cut.c_str()));
  //t_EWK->Draw(var+">>h_EWK",   "weight"*TCut(cut.c_str()) );
  h_EWK->Scale(lumiScale);
  dir->cd();
  h_EWK->Write("EWK", TObject::kOverwrite);

  // BBB
  for(int bin = 1; bin<=h_EWK->GetNbinsX(); bin++ ){
    TH1F* h_EWK_b_up   = (TH1F*)h->Clone(Form("h_EWK_%d_Up",bin));
    TH1F* h_EWK_b_down = (TH1F*)h->Clone(Form("h_EWK_%d_Down",bin));
    bbb( h_EWK, h_EWK_b_up, h_EWK_b_down, bin);

    dir->cd();
    h_EWK_b_up  ->Write(Form("EWK_EWK%sbin%dUp",  category.Data(), bin), TObject::kOverwrite);
    h_EWK_b_down->Write(Form("EWK_EWK%sbin%dDown",category.Data(), bin), TObject::kOverwrite);
 
  }




  cout << "TTH125:      "  << h_TTH125->Integral()         << endl;
  cout << "TTJets (HFbb): "<< h_TTJetsSemiLept_HFbb->Integral() << endl;
  cout << "TTJets (HFb): " << h_TTJetsSemiLept_HFb->Integral() << endl;
  cout << "TTJets (LF): "  << h_TTJetsSemiLept_LF->Integral() << endl;
  cout << "SingleT:     "  << h_SingleT->Integral() << endl;
  cout << "TTV:         "  << h_TTV->Integral() << endl;
  cout << "EWK:         "  << h_EWK->Integral() << endl;
  cout << "DiBoson:     "  << h_DiBoson->Integral() << endl;

  float rTTH125    = h_TTH125->Integral() ;
  float rTTJets_HFbb = h_TTJetsSemiLept_HFbb->Integral();
  float rTTJets_HFb  = h_TTJetsSemiLept_HFb->Integral();
  float rTTJets_LF = h_TTJetsSemiLept_LF->Integral();
  float rSingleT   = h_SingleT->Integral() ;
  float rTTV       = h_TTV->Integral() ;
  float rEWK       = h_EWK->Integral() ;
  float rDiBoson   = h_DiBoson->Integral() ;

  float observation = int(rTTH125 + rTTJets_HFbb + rTTJets_HFb + rTTJets_LF  + rTTV + rSingleT /*+ rEWK + rDiBoson*/ );

  TH1F* h_data = (TH1F*)h->Clone("h_data");
  h_data->Reset();
  h_data->Add(h_TTH125,1.0);
  h_data->Add(h_TTJetsSemiLept_HFbb,1.0);
  h_data->Add(h_TTJetsSemiLept_HFb,1.0);
  h_data->Add(h_TTJetsSemiLept_LF,1.0);
  h_data->Add(h_TTV);
  h_data->Add(h_SingleT);
  //h_data->Add(h_EWK);
  //h_data->Add(h_DiBoson);

  h_data->Scale(int(observation)/h_data->Integral());
  dir->cd();
  h_data->Write("data_obs", TObject::kOverwrite);

  fout->Close();


  ofstream out(Form("%s_2.txt",(fname+"_"+category).Data()));
  out.precision(8);
  string longspace  = "              ";
  string shortspace = "          ";


  //string null(category.Data());
  string null("");


  string line1("imax 1"); out << line1 << endl;
  string line2("jmax *"); out << line2 << endl;
  string line3("kmax *"); out << line3 << endl;
  string line4(Form("shapes *  *    %s_2.root  $CHANNEL/$PROCESS $CHANNEL/$PROCESS_$SYSTEMATIC", fname.Data())); out << line4 << endl;
  out<<endl;
  string line5(Form("observation %.0f", observation )); out << line5 << endl;
  out<<endl;


  // PROC AND NORM
  string line("bin                         ");
  if( rTTH125>0 )    line += string(Form("%s    ", (fname+"_"+category).Data()));
  if( rTTJets_HFbb>0 ) line += string(Form("%s    ", (fname+"_"+category).Data()));
  if( rTTJets_HFb>0 )  line += string(Form("%s    ", (fname+"_"+category).Data()));
  if( rTTJets_LF>0 ) line += string(Form("%s    ", (fname+"_"+category).Data()));
  if( rTTV>0 )       line += string(Form("%s    ", (fname+"_"+category).Data()));
  if( rSingleT>0 )   line += string(Form("%s    ", (fname+"_"+category).Data()));
  out<<line;
  out<<endl; 
  line = "process                     ";
  if( rTTH125>0 )    line += "TTH125      ";
  if( rTTJets_HFbb>0 ) line += "TTJetsHFbb    ";
  if( rTTJets_HFb>0 )  line += "TTJetsHFb    ";

  if( rTTJets_LF>0 ) line += "TTJetsLF    ";
  if( rTTV>0 )       line += "TTV      ";
  if( rSingleT>0 )   line += "SingleT     ";
  out<<line;
  out<<endl;
  line = "process                       ";
  if( rTTH125>0 )    line += "0          ";
  if( rTTJets_HFbb>0 ) line += "1          ";
  if( rTTJets_HFb>0 )  line += "2          ";
  if( rTTJets_LF>0 ) line += "3          ";
  if( rTTV>0 )       line += "4          ";
  if( rSingleT>0 )   line += "5          ";
  out<<line;
  out<<endl;
  line = "rate                        ";
  if( rTTH125>0 )    line += string(Form("%.4f     ", rTTH125 ));
  if( rTTJets_HFbb>0 ) line += string(Form("%.4f     ", rTTJets_HFbb ));
  if( rTTJets_HFb>0 )  line += string(Form("%.4f     ", rTTJets_HFb ));
  if( rTTJets_LF>0 ) line += string(Form("%.4f     ", rTTJets_LF));
  if( rTTV>0 )       line += string(Form("%.4f     ", rTTV));
  if( rSingleT>0 )   line += string(Form("%.4f     ", rSingleT));
  out<<line;
  out<<endl;
  out<<endl;
  out<< "-----------------------------------------------------------------" << endl;


  // SYSTEMATICS
  line = "lumi                  lnN  ";
  if( rTTH125>0 )    line += "1.022      ";
  if( rTTJets_HFbb>0 ) line += "1.022      ";
  if( rTTJets_HFb>0 )  line += "1.022      ";
  if( rTTJets_LF>0 ) line += "1.022      ";
  if( rTTV>0 )       line += "1.022      ";
  if( rSingleT>0 )   line += "1.022      ";
  out<<line;
  out<<endl;
  line = "csv                   shape  ";
  if( rTTH125>0 )    line += "1.0        ";
  if( rTTJets_HFbb>0 ) line += "1.0        ";
  if( rTTJets_HFb>0 )  line += "1.0        ";
  if( rTTJets_LF>0 ) line += "1.0        ";
  if( rTTV>0 )       line += " -         ";
  if( rSingleT>0 )   line += " -         ";
  out<<line;
  out<<endl;
  line = "JEC                   shape  ";
  if( rTTH125>0 )    line += "1.0        ";
  if( rTTJets_HFbb>0 ) line += "1.0        ";
  if( rTTJets_HFb>0 )  line += "1.0        ";
  if( rTTJets_LF>0 ) line += "1.0        ";
  if( rTTV>0 )       line += " -         ";
  if( rSingleT>0 )   line += " -         ";
  out<<line;
  out<<endl;

  for(int bin = 1; bin<=h_TTH125->GetNbinsX(); bin++ ){
    if(h_TTH125->GetBinContent(bin)>0){
      line = string(Form("TTH125%sbin%d        shape  ", category.Data(), bin));
      if( rTTH125>0 )    line += "1.0        ";
      if( rTTJets_HFbb>0 ) line += " -         ";
      if( rTTJets_HFb>0 )  line += " -         ";
      if( rTTJets_LF>0 ) line += " -         ";
      if( rTTV>0 )       line += " -         ";
      if( rSingleT>0 )   line += " -         ";
      out<<line;
      out<<endl;
    }
  }  
  for(int bin = 1; bin<=h_TTJetsSemiLept_HFbb->GetNbinsX(); bin++ ){
    if(h_TTJetsSemiLept_HFbb->GetBinContent(bin)>0){
      line = string(Form("TTJetsHFbb%sbin%d      shape  ", category.Data(), bin));
      if( rTTH125>0 )    line += " -         ";
      if( rTTJets_HFbb>0 ) line += "1.0        ";
      if( rTTJets_HFb>0 )  line += " -        ";
      if( rTTJets_LF>0 ) line += " -         ";
      if( rTTV>0 )       line += " -         ";
      if( rSingleT>0 )   line += " -         ";
      out<<line;
      out<<endl;
    }
  }  
  for(int bin = 1; bin<=h_TTJetsSemiLept_HFb->GetNbinsX(); bin++ ){
    if(h_TTJetsSemiLept_HFb->GetBinContent(bin)>0){
      line = string(Form("TTJetsHFb%sbin%d      shape  ", category.Data(), bin));
      if( rTTH125>0 )    line += " -         ";
      if( rTTJets_HFbb>0 ) line += " -        ";
      if( rTTJets_HFb>0 )  line += "1.0        ";
      if( rTTJets_LF>0 ) line += " -         ";
      if( rTTV>0 )       line += " -         ";
      if( rSingleT>0 )   line += " -         ";
      out<<line;
      out<<endl;
    }
  }  
  for(int bin = 1; bin<=h_TTJetsSemiLept_LF->GetNbinsX(); bin++ ){
    if(h_TTJetsSemiLept_LF->GetBinContent(bin)>0){
      line = string(Form("TTJetsLF%sbin%d      shape  ", category.Data(), bin));
      if( rTTH125>0 )    line += " -         ";
      if( rTTJets_HFbb>0 ) line += " -         ";
      if( rTTJets_HFb>0 )  line += " -         ";
      if( rTTJets_LF>0 ) line += "1.0        ";
      if( rTTV>0 )       line += " -         ";
      if( rSingleT>0 )   line += " -         ";
      out<<line;
      out<<endl;
    }
  }  
  for(int bin = 1; bin<=h_TTV->GetNbinsX(); bin++ ){
    if(h_TTV->GetBinContent(bin)>0){
      line = string(Form("TTV%sbin%d           shape  ", category.Data(), bin));
      if( rTTH125>0 )    line += " -         ";
      if( rTTJets_HFbb>0 ) line += " -         ";
      if( rTTJets_HFb>0 )  line += " -         ";
      if( rTTJets_LF>0 ) line += " -         ";
      if( rTTV>0 )       line += "1.0        ";
      if( rSingleT>0 )   line += " -         ";
      out<<line;
      out<<endl;
    }
  } 
  for(int bin = 1; bin<=h_SingleT->GetNbinsX(); bin++ ){
    if(h_SingleT->GetBinContent(bin)>0){
      line = string(Form("SingleT%sbin%d       shape  ", category.Data(), bin));
      if( rTTH125>0 )    line += " -         ";
      if( rTTJets_HFbb>0 ) line += " -         ";
      if( rTTJets_HFb>0 )  line += " -         ";
      if( rTTJets_LF>0 ) line += " -         ";
      if( rTTV>0 )       line += " -         ";
      if( rSingleT>0 )   line += "1.0        ";
      out<<line;
      out<<endl;
    }
  } 

  line = "Norm_TTbb             lnN    ";
  if( rTTH125>0 )    line += " -         ";
  if( rTTJets_HFbb>0 ) line += "1.50       ";
  if( rTTJets_HFb>0 )  line += " -       ";
  if( rTTJets_LF>0 ) line += " -         ";
  if( rTTV>0 )       line += " -         ";
  if( rSingleT>0 )   line += " -         ";
  out<<line;
  out<<endl;

  line = "Norm_TTb              lnN    ";
  if( rTTH125>0 )    line += " -         ";
  if( rTTJets_HFbb>0 ) line += " -       ";
  if( rTTJets_HFb>0 )  line += "1.50       ";
  if( rTTJets_LF>0 ) line += " -         ";
  if( rTTV>0 )       line += " -         ";
  if( rSingleT>0 )   line += " -         ";
  out<<line;
  out<<endl;

  line = "Norm_TTV              lnN    ";
  if( rTTH125>0 )    line += " -         ";
  if( rTTJets_HFbb>0 ) line += " -         ";
  if( rTTJets_HFb>0 )  line += " -         ";
  if( rTTJets_LF>0 ) line += " -         ";
  if( rTTV>0 )       line += "1.20       ";
  if( rSingleT>0 )   line += " -         ";
  out<<line;
  out<<endl;

  line = "Norm_SingleT          lnN    ";
  if( rTTH125>0 )    line += " -         ";
  if( rTTJets_HFbb>0 ) line += " -         ";
  if( rTTJets_HFb>0 )  line += " -         ";
  if( rTTJets_LF>0 ) line += " -         ";
  if( rTTV>0 )       line += " -         ";
  if( rSingleT>0 )   line += "1.20       ";
  out<<line;
  out<<endl;

  line = "QCDscale_TTH          lnN    ";
  if( rTTH125>0 )    line += "1.12       ";
  if( rTTJets_HFbb>0 ) line += " -         ";
  if( rTTJets_HFb>0 )  line += " -         ";
  if( rTTJets_LF>0 ) line += " -         ";
  if( rTTV>0 )       line += " -         ";
  if( rSingleT>0 )   line += " -         ";
  out<<line;
  out<<endl;

  line = string(Form("QCDscale%s_TTJetsHF     lnN    ", null.c_str()));
  if( rTTH125>0 )    line += " -         ";
  if( rTTJets_HFbb>0 ) line += "1.35       ";
  if( rTTJets_HFb>0 )  line += "1.35       ";
  if( rTTJets_LF>0 ) line += " -         ";
  if( rTTV>0 )       line += " -         ";
  if( rSingleT>0 )   line += " -         ";
  out<<line;
  out<<endl;

  line = string(Form("QCDscale%s_TTJetsLF     lnN    ", null.c_str()));
  if( rTTH125>0 )    line += " -         ";
  if( rTTJets_HFbb>0 ) line += " -         ";
  if( rTTJets_HFb>0 )  line += " -         ";
  if( rTTJets_LF>0 ) line += "1.35       ";
  if( rTTV>0 )       line += " -         ";
  if( rSingleT>0 )   line += " -         ";
  out<<line;
  out<<endl;

  line = "pdf_gg                lnN    ";
  if( rTTH125>0 )    line += "1.03       ";
  if( rTTJets_HFbb>0 ) line += "1.03       ";
  if( rTTJets_HFb>0 )  line += "1.03       ";
  if( rTTJets_LF>0 ) line += "1.03       ";
  if( rTTV>0 )       line += "1.03       ";
  if( rSingleT>0 )   line += "1.03       ";
  out<<line;
  out<<endl;
  out<<endl;
  out<< "-----------------------------------------------------------------" << endl;

}



void produceAll_2( float LumiScale = 19.5/12.1){

  produce_2("SL", "type==0",  "cat1",                       0.02, 0.50, LumiScale ,   1, 4, 1);
  produce_2("SL", "type==1",  "cat2",                       0.02, 0.7,  LumiScale ,   1, 5, 1);
  produce_2("SL", "type==2 && flag_type2>0",   "cat3",      0.02, 0.8,  LumiScale ,   1, 5, 1);
  produce_2("SL", "type==2 && flag_type2<=0",  "cat4",      0.02, 0.7,  LumiScale ,   1, 5, 1);
  produce_2("SL", "type==3",  "cat5",                       0.02, 0.5,  LumiScale ,   1, 6, 1);  
  produce_2("DL", "type==6",  "cat6",                       0.02, 0.,   LumiScale*2 , 0, 5, 0);
  produce_2("DL", "type==7",  "cat7",                       0.02, 0.,   LumiScale*2 , 0, 6, 0);


}

// cat1  
// 0.1  => 8.28
// 0.25 => 7.93
// 0.5  => 7.95
// 0.75 => 8.03
// 1.0  => 8.09

// cat2
// 0.10 => 12.2
// 0.25 => 11.3
// 0.50 => 11.3
// 1.50 => 12.0

// cat4
// 0.25 => 8.8
// 0.50 => 8.6
// 0.75 => 9.2

// cat5
// 0.10 => 10.8
// 0.25 => 11.3
// 0.5  => 10.8
// 0.75 => 11.3

// cat6
// 0.25 => 9.3
// 0.50 => 8.8
// 0.75 => 8.8
// 1.50 => 8.8


void makePlot(){

  gStyle->SetPaintTextFormat("g");

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
  c1->SetLogy(1);
  
  TLegend* leg = new TLegend(0.120743,0.153846,0.321981, 0.26049, NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->SetFillColor(0);


  float X[]        = {0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5 };
  float expY[]     = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float expY1sL[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float expY1sH[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float expY2sL[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float expY2sH[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

  float expXs[]  = {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5};
  float expYs[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};


  vector<string> categories;
  categories.push_back("SL_cat1");
  categories.push_back("SL_cat2");
  categories.push_back("SL_cat3");
  categories.push_back("SL_cat4");
  categories.push_back("SL_cat5");
  categories.push_back("SL");
  categories.push_back("DL_cat6");
  categories.push_back("DL_cat7");
  categories.push_back("DL");
  categories.push_back("COMB");

  int nBins = categories.size();

  for( int b = 0; b < nBins; b++){

    TFile* f = TFile::Open(("higgsCombine"+categories[b]+".Asymptotic.mH120.root").c_str());
    if( f==0 ) continue;
    
    Double_t r;
    TTree* limit = (TTree*)f->Get("limit");
    limit->SetBranchAddress("limit",&r);

    for(int k = 0 ; k< limit->GetEntries() ; k++){
      limit->GetEntry(k);
      if(k==0) expY2sL[b] = r;
      if(k==1) expY1sL[b] = r;
      if(k==2) expY[b]    = r;
      if(k==3) expY1sH[b] = r;
      if(k==4) expY2sH[b] = r;
    }


    cout << categories[b] << ": r<" <<  expY[b] << " @ 95% CL" << endl;

    expY1sH[b] = TMath::Abs(expY1sH[b]-expY[b]);
    expY1sL[b] = TMath::Abs(expY1sL[b]-expY[b]);
    expY2sH[b] = TMath::Abs(expY2sH[b]-expY[b]);
    expY2sL[b] = TMath::Abs(expY2sL[b]-expY[b]);


  }


  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("Simulation #sqrt{s}=8 TeV, L=19.5 fb^{-1}");

  TGraphAsymmErrors* expected = new TGraphAsymmErrors(10, X, expY, expXs ,expXs,  expYs,   expYs);
  TGraphAsymmErrors* oneSigma = new TGraphAsymmErrors(10, X, expY, expXs, expXs,  expY1sL, expY1sH);
  TGraphAsymmErrors* twoSigma = new TGraphAsymmErrors(10, X, expY, expXs, expXs,  expY2sL, expY2sH);


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
 
  mg->Add(twoSigma);
  mg->Add(oneSigma);
  mg->Add(expected);

  mg->Draw("a2");  
  expected->Draw("pSAME");

  TH1F* hT = new TH1F("hT", "", 10, 0, 10);
  for(int k = 1; k <= hT->GetNbinsX(); k++){
    int approx = int(expY[k-1]*10);
    hT->SetBinContent(k, approx/10.   );
  }
  hT->SetMarkerSize(2.);
  hT->Draw("TEXT0SAME");

  TF1 *line = new TF1("line","1",0,10);
  line->SetLineColor(kRed);
  line->SetLineWidth(2);

  line->Draw("SAME");

  TF1 *lineML = new TF1("lineML","2.4",0,10);
  lineML->SetLineColor(kBlue);
  lineML->SetLineStyle(kDashed);
  lineML->SetLineWidth(3);

  lineML->Draw("SAME");

  TF1 *lineTTH = new TF1("lineTTH","4.1",0,10);
  lineTTH->SetLineColor(kMagenta);
  lineTTH->SetLineStyle(kDashed);
  lineTTH->SetLineWidth(3);

  lineTTH->Draw("SAME");

  c1->cd();
  gPad->Modified();
  mg->GetXaxis()->Set(10,0,10);

  //mg->GetXaxis()->SetRange(-1,11);

  for( int b = 0; b < nBins; b++){
    mg->GetXaxis()->SetBinLabel(b+1, categories[b].c_str() );
  }

  mg->GetYaxis()->SetTitleOffset(0.97);
  mg->SetMinimum(0.8);
  mg->SetMaximum( 50);
  mg->GetXaxis()->SetTitle("");
  mg->GetYaxis()->SetTitle("95% CL upper limit on #mu = #sigma/#sigma_{SM}");

  leg->AddEntry( lineTTH, "HIG-13-019", "L");
  leg->AddEntry( lineML,  "HIG-13-020", "L");
  leg->Draw();

  TPaveText *pt = new TPaveText(0.106811,0.155594,0.407121,0.286713,"brNDC");
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetFillColor(10);
  pt->SetTextSize(0.04);
  pt->SetTextAlign(11);
  pt->AddText(Form("Comb: #mu < %.1f at 95%% CL", expY[categories.size()-1]))->SetTextColor(kRed);
  pt->AddText(Form("(SL: #mu < %.1f)", expY[categories.size()-3]))->SetTextColor(kBlack);
  pt->AddText(Form("(DL: #mu < %.1f)", expY[categories.size()-2]))->SetTextColor(kBlack);

  //pt->Draw();

  //c1->Modified();
  //c1->Draw();

  TLine* del = new TLine(6., 50., 6., 0. );
  del->SetLineWidth(3);
  del->SetLineStyle(kSolid);
  del->SetLineColor(kBlack);
  del->Draw("SAME");

  TLine* del_ = new TLine(5., 50., 5., 0. );
  del_->SetLineWidth(2);
  del_->SetLineStyle(kDotted);
  del_->SetLineColor(kBlack);
  del_->Draw("SAME");

  TLine* del2 = new TLine(9., 50., 9., 0. );
  del2->SetLineWidth(3);
  del2->SetLineStyle(kSolid);
  del2->SetLineColor(kBlack);
  del2->Draw("SAME");

  TLine* del2_ = new TLine(8., 50., 8., 0. );
  del2_->SetLineWidth(2);
  del2_->SetLineStyle(kDotted);
  del2_->SetLineColor(kBlack);
  del2_->Draw("SAME");

  if(1){
    c1->SaveAs("Limits.png");
    c1->SaveAs("Limits.pdf");
  }

}



void plot_category(string type = "SL", 
		   string cat  = "cat1",
		   string header = "Cat 1",
		   string fname  = ""
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

  vector<string> samples;
  if( type.find("SL")!=string::npos ){
    samples.push_back("SingleT");
    samples.push_back("TTV");
    samples.push_back("TTJetsHF");
    samples.push_back("TTJetsLF");
    samples.push_back("TTH125");
  }
  if( type.find("DL")!=string::npos ){
    samples.push_back("SingleT");
    samples.push_back("TTV");
    samples.push_back("TTJetsHF");
    samples.push_back("TTJetsLF");
    samples.push_back("TTH125");
  }

  TH1F* hS = 0;
  TH1F* hErr = 0;

  TFile* f = TFile::Open((type+".root").c_str());
  if(f==0 || f->IsZombie() ) return;

  for(unsigned int sample = 0; sample < samples.size(); sample++){

    cout << "Doing " << samples[sample] << endl;

  
    TH1F* h = (TH1F*)f->Get((type+"_"+cat+"/"+samples[sample]).c_str());
    if( h==0 ){
      cout << (type+"_"+cat+"/"+samples[sample]) << " not found" << endl;
      continue;
    }
    cout << "Events = $" << h->Integral() << " \\pm " << (h->GetEntries()>0 ? sqrt(h->GetEntries())*h->Integral()/h->GetEntries() : 0.0 ) << " $" << endl;

    if(samples[sample].find("TTH125") == string::npos){
      if( hErr==0 ){
	hErr = (TH1F*)h->Clone("hErr");
	hErr->Reset();
	leg->AddEntry(hErr, "MC unc. (stat.)", "L");
      }
      hErr->Add( h, 1.0);
    }
   
    if( samples[sample].find("TTH125") != string::npos ){
      hS = (TH1F*)h->Clone("hS");
      hS->Reset();
      hS->Add(h,5.0);
      hS->SetLineWidth(3);
      hS->SetLineStyle(kDashed);
      hS->SetLineColor(kRed);

      h->SetLineColor( kRed );
      h->SetFillColor( kRed );
      h->SetFillStyle( 3002 );
      leg->AddEntry(h, "t#bar{t}H", "F");
    }
    if( samples[sample].find("TTJetsHF") != string::npos ){
      leg->AddEntry(h, "t#bar{t}+jets HF", "F");
      h->SetLineColor( 17 );
      h->SetFillColor( 17 );
    }
    if( samples[sample].find("TTJetsLF") != string::npos ){
      leg->AddEntry(h, "t#bar{t}+jets LF", "F");
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


    aStack->Add( h );
  }

  if(hErr==0 || hS==0) return;

  hErr->GetYaxis()->SetTitle("Events");
  if(type.find("SL")!=string::npos) 
    //hErr->GetXaxis()->SetTitle("f_{ttH} / ( f_{ttH} + f_{ttbb} + f_{ttjj} )");
    hErr->GetXaxis()->SetTitle("P_{s/b}(y)");
  else
    //hErr->GetXaxis()->SetTitle("f_{ttH} / ( f_{ttH} + f_{ttbb} )");
    hErr->GetXaxis()->SetTitle("P_{s/b}(y)");
  hErr->SetTitle("Simulation #sqrt{s}=8 TeV, L=19.5 fb^{-1}");
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
  if(type.find("SL")!=string::npos) line->Draw("SAME");

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
 
  if(   type.find("SL")!=string::npos ){
    pt1->Draw();
    pt2->Draw();
  }

  cout << "Signal = " << hS->Integral()/5. << endl;

  if(1){
    c1->SaveAs(  (fname+"_AN"+".png").c_str() );
    c1->SaveAs(  (fname+"_AN"+".pdf").c_str() );
  }

}


void plotAll(){


  
  plot_category( "SL", 
		 "cat1",
		 "Cat1 (4b2j W-tag)",
		 "Plot_SL_Cat1"
		 );
  
  
  plot_category( "SL", 
		 "cat2",
		 "Cat2 (4b2j !W-tag)",
		 "Plot_SL_Cat2"
		 );


  plot_category( "SL", 
		 "cat3",
		 "Cat3 (4b1j cs-tag)",
		 "Plot_SL_Cat3"
		 );
  

  plot_category( "SL", 
		 "cat4",
		 "Cat4 (4b1j !cs-tag)",
		 "Plot_SL_Cat4"
		 );


  plot_category( "SL", 
		 "cat5",
		 "Cat5 (4b3j)",
		 "Plot_SL_Cat5"
		 );

  
  plot_category( "DL", 
		 "cat6",
		 "Cat6 (4b), tight",
		 "Plot_DL_Cat6"
		 );

  plot_category( "DL", 
		 "cat7",	
		 "Cat7 (4b), loose",
		 "Plot_DL_Cat7"
		 );
    
  

}




void plot_category_2(string type = "SL", 
		     string cat  = "cat1",
		     string header = "Cat 1",
		     string fname  = ""
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

  vector<string> samples;
  if( type.find("SL")!=string::npos ){
    samples.push_back("SingleT");
    samples.push_back("TTV");
    samples.push_back("TTJetsHFbb");
    samples.push_back("TTJetsHFb");
    samples.push_back("TTJetsLF");
    samples.push_back("TTH125");
  }
  if( type.find("DL")!=string::npos ){
    samples.push_back("SingleT");
    samples.push_back("TTV");
    samples.push_back("TTJetsHFbb");
    samples.push_back("TTJetsHFb");
    samples.push_back("TTJetsLF");
    samples.push_back("TTH125");
  }

  TH1F* hS = 0;
  TH1F* hErr = 0;

  TFile* f = TFile::Open((type+"_2.root").c_str());
  if(f==0 || f->IsZombie() ) return;

  for(unsigned int sample = 0; sample < samples.size(); sample++){

    cout << "Doing " << samples[sample] << endl;

  
    TH1F* h = (TH1F*)f->Get((type+"_"+cat+"/"+samples[sample]).c_str());
    if( h==0 ){
      cout << (type+"_"+cat+"/"+samples[sample]) << " not found" << endl;
      continue;
    }
    cout << "Events = $" << h->Integral() << " \\pm " << (h->GetEntries()>0 ? sqrt(h->GetEntries())*h->Integral()/h->GetEntries() : 0.0 ) << " $" << endl;

    if(samples[sample].find("TTH125") == string::npos){
      if( hErr==0 ){
	hErr = (TH1F*)h->Clone("hErr");
	hErr->Reset();
	leg->AddEntry(hErr, "MC unc. (stat.)", "L");
      }
      hErr->Add( h, 1.0);
    }
   
    if( samples[sample].find("TTH125") != string::npos ){
      hS = (TH1F*)h->Clone("hS");
      hS->Reset();
      hS->Add(h,5.0);
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
      h->SetLineColor( 16 );
      h->SetFillColor( 16 );
    }
    else if( samples[sample].find("TTJetsHFb") != string::npos ){
      leg->AddEntry(h, "t#bar{t} + b", "F");
      h->SetLineColor( 17 );
      h->SetFillColor( 17 );
    }
    if( samples[sample].find("TTJetsLF") != string::npos ){
      leg->AddEntry(h, "t#bar{t} + jj", "F");
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


    aStack->Add( h );
  }

  if(hErr==0 || hS==0) return;

  hErr->GetYaxis()->SetTitle("Events");
  if(type.find("SL")!=string::npos) 
    //hErr->GetXaxis()->SetTitle("f_{ttH} / ( f_{ttH} + f_{ttbb} + f_{ttjj} )");
    hErr->GetXaxis()->SetTitle("P_{s/b}(y)");
  else
    //hErr->GetXaxis()->SetTitle("f_{ttH} / ( f_{ttH} + f_{ttbb} )");
    hErr->GetXaxis()->SetTitle("P_{s/b}(y)");
  hErr->SetTitle("Simulation #sqrt{s}=8 TeV, L=19.5 fb^{-1}");
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
  if(type.find("SL")!=string::npos) line->Draw("SAME");

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
 
  if(   type.find("SL")!=string::npos ){
    pt1->Draw();
    pt2->Draw();
  }

  cout << "Signal = " << hS->Integral()/5. << endl;

  if(0){
    c1->SaveAs(  (fname+"_AN"+"_2.png").c_str() );
    c1->SaveAs(  (fname+"_AN"+"_2.pdf").c_str() );
  }

}


void plotAll_2(){


  
  plot_category_2( "SL", 
		   "cat1",
		   "Cat1 (4b2j W-tag)",
		   "Plot_SL_Cat1_2"
		   );
  
  
  plot_category_2( "SL", 
		   "cat2",
		   "Cat2 (4b2j !W-tag)",
		   "Plot_SL_Cat2_2"
		   );
  

  plot_category_2( "SL", 
		   "cat3",
		   "Cat3 (4b1j cs-tag)",
		   "Plot_SL_Cat3_2"
		   );
  

  plot_category_2( "SL", 
		   "cat4",
		   "Cat4 (4b1j !cs-tag)",
		   "Plot_SL_Cat4_2"
		   );
  

  plot_category_2( "SL", 
		   "cat5",
		   "Cat5 (4b3j)",
		   "Plot_SL_Cat5_2"
		   );
  
  
  plot_category_2( "DL", 
		   "cat6",
		   "Cat6 (4b), tight",
		   "Plot_DL_Cat6_2"
		   );
  
  plot_category_2( "DL", 
		   "cat7",	
		   "Cat7 (4b), loose",
		   "Plot_DL_Cat7_2"
		   );
    
  

}





void plot_syst(string type = "SL",
	       string cat  = "cat1",
	       string proc = "TTH125",
	       string syst = "csv",
	       string header = "Cat 1",
	       string fname  = ""
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

  vector<string> samples;
  samples.push_back( proc );

  TH1F* hN = 0;
  TH1F* hU = 0;
  TH1F* hD = 0;

  TFile* f = TFile::Open((type+".root").c_str());
  if(f==0 || f->IsZombie() ) return;

  for(unsigned int sample = 0; sample < samples.size(); sample++){

    cout << "Doing " << samples[sample] << endl;

  
    hN = (TH1F*)((TH1F*)f->Get((type+"_"+cat+"/"+samples[sample]).c_str()))->Clone("hN");
    if( hN==0 ){
      cout << (type+"_"+cat+"/"+samples[sample]) << " not found" << endl;
      continue;
    }
    else{
      cout << hN->Integral() << endl;
    }

    hU = (TH1F*)((TH1F*)f->Get((type+"_"+cat+"/"+samples[sample]+"_"+syst+"Up").c_str()))->Clone("hU");
    if( hU==0 ){
      cout << (type+"_"+cat+"/"+samples[sample]+"_"+syst+"Up") << " not found" << endl;
      continue;
    }
    else{
      cout << hU->Integral() << endl;
    }

    hD = (TH1F*)((TH1F*)f->Get((type+"_"+cat+"/"+samples[sample]+"_"+syst+"Down").c_str()))->Clone("hD");
    if( hD==0 ){
      cout << (type+"_"+cat+"/"+samples[sample]+"_"+syst+"Down") << " not found" << endl;
      continue;
    }
    else{
      cout << hD->Integral() << endl;
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
    hU->GetYaxis()->SetTitle("Events");
    if(type.find("SL")!=string::npos) 
      hU->GetXaxis()->SetTitle("P_{s/b}(y)");
    else
      hU->GetXaxis()->SetTitle("P_{s/b}(y)");
    hU->SetTitle("Simulation #sqrt{s}=8 TeV, L=19.5 fb^{-1}");
    hU->SetTitleSize  (0.04,"X");
    hU->SetTitleOffset(0.95,"X");
    hU->Draw("HISTE");
    hN->Draw("HISTSAME");
    hD->Draw("HISTESAME");
  }
  else{
    hD->SetMinimum(0.0);
    hD->GetYaxis()->SetTitle("Events");
    if(type.find("SL")!=string::npos) 
      hD->GetXaxis()->SetTitle("P_{s/b}(y)");
    else
      hD->GetXaxis()->SetTitle("P_{s/b}(y)");
    hD->SetTitle("Simulation #sqrt{s}=8 TeV, L=19.5 fb^{-1}");
    hD->SetTitleSize  (0.04,"X");
    hD->SetTitleOffset(0.95,"X");
    hD->Draw("HISTE");
    hN->Draw("HISTSAME");
    hU->Draw("HISTESAME");
  }

  leg->SetHeader(header.c_str());
  leg->AddEntry( hN, (syst+" nominal").c_str(), "L");
  leg->AddEntry( hU, (syst+" up").c_str(),      "L");
  leg->AddEntry( hD, (syst+" down").c_str(),    "L");

  leg->Draw();

  if(1){
    c1->SaveAs(  ("Plots/"+fname+"_AN"+".png").c_str() );
    c1->SaveAs(  ("Plots/"+fname+"_AN"+".pdf").c_str() );
  }

}



void plot_systAll(){

  plot_syst("SL", "cat1", "TTH125", "csv", "Cat 1, ttH", "Plots_SL_Cat1_TTH_csv");
  plot_syst("SL", "cat2", "TTH125", "csv", "Cat 2, ttH", "Plots_SL_Cat2_TTH_csv");
  plot_syst("SL", "cat3", "TTH125", "csv", "Cat 3, ttH", "Plots_SL_Cat3_TTH_csv");
  plot_syst("SL", "cat4", "TTH125", "csv", "Cat 4, ttH", "Plots_SL_Cat4_TTH_csv");
  plot_syst("SL", "cat5", "TTH125", "csv", "Cat 5, ttH", "Plots_SL_Cat5_TTH_csv");
  plot_syst("DL", "cat6", "TTH125", "csv", "Cat 6, ttH", "Plots_DL_Cat6_TTH_csv");
  plot_syst("DL", "cat7", "TTH125", "csv", "Cat 7, ttH", "Plots_DL_Cat7_TTH_csv");

  plot_syst("SL", "cat1", "TTJetsHF", "csv", "Cat 1, ttbb", "Plots_SL_Cat1_TTJetsHF_csv");
  plot_syst("SL", "cat2", "TTJetsHF", "csv", "Cat 2, ttbb", "Plots_SL_Cat2_TTJetsHF_csv");
  plot_syst("SL", "cat3", "TTJetsHF", "csv", "Cat 3, ttbb", "Plots_SL_Cat3_TTJetsHF_csv");
  plot_syst("SL", "cat4", "TTJetsHF", "csv", "Cat 4, ttbb", "Plots_SL_Cat4_TTJetsHF_csv");
  plot_syst("SL", "cat5", "TTJetsHF", "csv", "Cat 5, ttbb", "Plots_SL_Cat5_TTJetsHF_csv");
  plot_syst("DL", "cat6", "TTJetsHF", "csv", "Cat 6, ttbb", "Plots_DL_Cat6_TTJetsHF_csv");
  plot_syst("DL", "cat7", "TTJetsHF", "csv", "Cat 7, ttbb", "Plots_DL_Cat7_TTJetsHF_csv");

  plot_syst("SL", "cat1", "TTJetsLF", "csv", "Cat 1, ttjj", "Plots_SL_Cat1_TTJetsLF_csv");
  plot_syst("SL", "cat2", "TTJetsLF", "csv", "Cat 2, ttjj", "Plots_SL_Cat2_TTJetsLF_csv");
  plot_syst("SL", "cat3", "TTJetsLF", "csv", "Cat 3, ttjj", "Plots_SL_Cat3_TTJetsLF_csv");
  plot_syst("SL", "cat4", "TTJetsLF", "csv", "Cat 4, ttjj", "Plots_SL_Cat4_TTJetsLF_csv");
  plot_syst("SL", "cat5", "TTJetsLF", "csv", "Cat 5, ttjj", "Plots_SL_Cat5_TTJetsLF_csv");
  plot_syst("DL", "cat6", "TTJetsLF", "csv", "Cat 6, ttjj", "Plots_DL_Cat6_TTJetsLF_csv");
  plot_syst("DL", "cat7", "TTJetsLF", "csv", "Cat 7, ttjj", "Plots_DL_Cat7_TTJetsLF_csv");

  plot_syst("SL", "cat1", "TTH125", "JEC", "Cat 1, ttH", "Plots_SL_Cat1_TTH_JEC");
  plot_syst("SL", "cat2", "TTH125", "JEC", "Cat 2, ttH", "Plots_SL_Cat2_TTH_JEC");
  plot_syst("SL", "cat3", "TTH125", "JEC", "Cat 3, ttH", "Plots_SL_Cat3_TTH_JEC");
  plot_syst("SL", "cat4", "TTH125", "JEC", "Cat 4, ttH", "Plots_SL_Cat4_TTH_JEC");
  plot_syst("SL", "cat5", "TTH125", "JEC", "Cat 5, ttH", "Plots_SL_Cat5_TTH_JEC");
  plot_syst("DL", "cat6", "TTH125", "JEC", "Cat 6, ttH", "Plots_DL_Cat6_TTH_JEC");
  plot_syst("DL", "cat7", "TTH125", "JEC", "Cat 7, ttH", "Plots_DL_Cat7_TTH_JEC");

  plot_syst("SL", "cat1", "TTJetsHF", "JEC", "Cat 1, ttbb", "Plots_SL_Cat1_TTJetsHF_JEC");
  plot_syst("SL", "cat2", "TTJetsHF", "JEC", "Cat 2, ttbb", "Plots_SL_Cat2_TTJetsHF_JEC");
  plot_syst("SL", "cat3", "TTJetsHF", "JEC", "Cat 3, ttbb", "Plots_SL_Cat3_TTJetsHF_JEC");
  plot_syst("SL", "cat4", "TTJetsHF", "JEC", "Cat 4, ttbb", "Plots_SL_Cat4_TTJetsHF_JEC");
  plot_syst("SL", "cat5", "TTJetsHF", "JEC", "Cat 5, ttbb", "Plots_SL_Cat5_TTJetsHF_JEC");
  plot_syst("DL", "cat6", "TTJetsHF", "JEC", "Cat 6, ttbb", "Plots_DL_Cat6_TTJetsHF_JEC");
  plot_syst("DL", "cat7", "TTJetsHF", "JEC", "Cat 7, ttbb", "Plots_DL_Cat7_TTJetsHF_JEC");

  plot_syst("SL", "cat1", "TTJetsLF", "JEC", "Cat 1, ttjj", "Plots_SL_Cat1_TTJetsLF_JEC");
  plot_syst("SL", "cat2", "TTJetsLF", "JEC", "Cat 2, ttjj", "Plots_SL_Cat2_TTJetsLF_JEC");
  plot_syst("SL", "cat3", "TTJetsLF", "JEC", "Cat 3, ttjj", "Plots_SL_Cat3_TTJetsLF_JEC");
  plot_syst("SL", "cat4", "TTJetsLF", "JEC", "Cat 4, ttjj", "Plots_SL_Cat4_TTJetsLF_JEC");
  plot_syst("SL", "cat5", "TTJetsLF", "JEC", "Cat 5, ttjj", "Plots_SL_Cat5_TTJetsLF_JEC");
  plot_syst("DL", "cat6", "TTJetsLF", "JEC", "Cat 6, ttjj", "Plots_DL_Cat6_TTJetsLF_JEC");
  plot_syst("DL", "cat7", "TTJetsLF", "JEC", "Cat 7, ttjj", "Plots_DL_Cat7_TTJetsLF_JEC");

}
