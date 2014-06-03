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
#include "TString.h"

#define UPDATE 0

void smoothen(TString dir   = "datacards/May25_2014_195fb/",
	      TString input = "MEM_New_rec_std_sb",
	      TString cat  = "MEM_cat1_L/",
	      TString hist = "TTJetsLFcc",
	      TString syst = "JEC"
	      ){

  
  // input file
  TFile *f = 0;
  if(UPDATE) 
    f = TFile::Open(dir+"/"+input+".root","UPDATE");
  else 
    f = TFile::Open(dir+"/"+input+".root","READ");

  // nominal histogram
  TH1F* h   = (TH1F*)f->Get( cat + "/" + hist );

  // shifted Up
  TH1F* h_u = (TH1F*)f->Get( cat + "/" + hist + "_" + syst + "Up" );

  // shifted Down
  TH1F* h_d = (TH1F*)f->Get( cat + "/" + hist + "_" + syst + "Down");

  if( !h || !h_u || !h_d ){
    cout << "Miss histo" << endl;
    return;
  }

  // number of bins
  int   nBins = h->GetNbinsX();

  int mergeBins = 0;
  if( string(cat.Data()).find("_H")!=string::npos ) 
    mergeBins = 1;

  TH1F* h2   = 0;
  TH1F* h2_u = 0;
  TH1F* h2_d = 0;

  if( mergeBins ){

    // merge first two bins
    h2 = new TH1F("h2", "", nBins-1, 0, 1);
    for(int b = 1; b <= nBins; b++){
      
      if(b<=2){
	float err1 = h->GetBinError(1);
	float err2 = h->GetBinError(2);
	float bin1 = h->GetBinContent(1);
	float bin2 = h->GetBinContent(2);
	h2->SetBinContent(b, bin1+bin2 );
	h2->SetBinError  (b, sqrt(err1*err1 + err2*err2) );
      }
      else{
	float err1 = h->GetBinError(b);
	float bin1 = h->GetBinContent(b);
	h2->SetBinContent(b-1, bin1 );
	h2->SetBinError  (b-1, err1 );
      }
      
    }
    
    // merge first two bins
    h2_u = new TH1F("h2_u", "", nBins-1, 0, 1);
    for(int b = 1; b <= nBins; b++){
      
      if(b<=2){
	float err1 = h_u->GetBinError(1);
	float err2 = h_u->GetBinError(2);
	float bin1 = h_u->GetBinContent(1);
	float bin2 = h_u->GetBinContent(2);
	h2_u->SetBinContent(b, bin1+bin2 );
	h2_u->SetBinError  (b, sqrt(err1*err1 + err2*err2) );
      }
      else{
	float err1 = h_u->GetBinError(b);
	float bin1 = h_u->GetBinContent(b);
	h2_u->SetBinContent(b-1, bin1 );
      h2_u->SetBinError  (b-1, err1 );
      }
      
    }
    
    // merge first two bins
    h2_d = new TH1F("h2_d", "", nBins-1, 0, 1);
    for(int b = 1; b <= nBins; b++){
      
      if(b<=2){
	float err1 = h_d->GetBinError(1);
	float err2 = h_d->GetBinError(2);
	float bin1 = h_d->GetBinContent(1);
	float bin2 = h_d->GetBinContent(2);
	h2_d->SetBinContent(b, bin1+bin2 );
	h2_d->SetBinError  (b, sqrt(err1*err1 + err2*err2) );
      }
      else{
	float err1 = h_d->GetBinError(b);
	float bin1 = h_d->GetBinContent(b);
	h2_d->SetBinContent(b-1, bin1 );
	h2_d->SetBinError  (b-1, err1 );
      }
      
    }
    
  }

  if( h2==0 ){
    h2   = h;
    h2_u = h_u;
    h2_d = h_d;
  }

  // for nominal shape
  TF1* pol4   = new TF1("pol4",  "-([0]+[1]+[2]+[3]+[4])+[0]*x+[1]*x*x + [2]*x*x*x + [3]*x*x*x*x + [4]*x*x*x*x*x",0,1);
  pol4 ->SetLineColor(kBlack);
  h2->Fit  ( pol4, "", "", 0,1);

  // for shifted Up
  TF1* pol4_u   = new TF1("pol4_u",  "-([0]+[1]+[2]+[3]+[4])+[0]*x+[1]*x*x + [2]*x*x*x + [3]*x*x*x*x + [4]*x*x*x*x*x",0,1);
  pol4_u ->SetLineColor(kRed);
  h2_u->Fit( pol4_u, "", "", 0,1);

  // for shifted Down
  TF1* pol4_d   = new TF1("pol4_d",  "-([0]+[1]+[2]+[3]+[4])+[0]*x+[1]*x*x + [2]*x*x*x + [3]*x*x*x*x + [4]*x*x*x*x*x",0,1);
  pol4_d ->SetLineColor(kBlue);
  h2_d->Fit( pol4_d, "", "", 0,1);


  // the original one
  TH1F* h3_u = (TH1F*)h->Clone(  hist + "_" + syst + "Up"   ); 
  TH1F* h3_d = (TH1F*)h->Clone(  hist + "_" + syst + "Down" );
 

  for(int b = 1; b <= nBins ; b++){

    float nominal =  h->GetBinContent(b);
    float x       =  h->GetBinCenter(b);

    if( mergeBins && b<=2 ){
      x = h->GetBinLowEdge(2);      
    }

    h3_u->SetBinContent( b, nominal*(pol4_u->Eval(x)/pol4->Eval(x) ));
    h3_d->SetBinContent( b, nominal*(pol4_d->Eval(x)/pol4->Eval(x) ));

  }

  if(UPDATE==0){
    h3_u->SetLineColor(kRed);  h3_u ->SetLineWidth(2);
    h3_d->SetLineColor(kBlue); h3_d ->SetLineWidth(2);
    h_u->SetMarkerColor(kRed); h_u->SetMarkerStyle(kOpenCircle);
    h_d->SetMarkerColor(kBlue);h_d->SetMarkerStyle(kOpenCircle);

    delete pol4; delete pol4_u; delete pol4_d;  

    h->SetMaximum( h->GetMaximum()*1.5 ) ;
    h->Draw("HIST");
    h3_u->Draw("HISTSAME");
    h3_d->Draw("HISTSAME");
    h_u->Draw("PSAME");
    h_d->Draw("PSAME");

    gPad->SaveAs("Plots/Studies/Smoothen_"+cat+"_"+hist+"_"+syst+".png");
    f->Close();
  }
  else{
    f->cd(cat);
    h3_u->Write( hist + "_" + syst + "Up" ,   TObject::kOverwrite);
    h3_d->Write( hist + "_" + syst + "Down" , TObject::kOverwrite);
    delete pol4; delete pol4_u; delete pol4_d;  
    f->Close();
  }


}


void smoothenAll(){

  vector<TString> hists;
  hists.push_back("TTJetsLF");
  hists.push_back("TTJetsLFcc");
  hists.push_back("TTJetsHFbb");
  hists.push_back("TTJetsHFb");
  hists.push_back("TTH125");
  hists.push_back("TTV");

  vector<TString> cats;
  cats.push_back("MEM_cat1_H");
  //cats.push_back("MEM_cat1_L");
  cats.push_back("MEM_cat2_H");
  //cats.push_back("MEM_cat2_L");
  cats.push_back("MEM_cat3_H");
  //cats.push_back("MEM_cat3_L");
  cats.push_back("MEM_cat6_H");
  //cats.push_back("MEM_cat6_L");
  
  for(unsigned int i = 0 ; i < hists.size(); i++ ){
    for(unsigned int j = 0 ; j < cats.size(); j++ ){
      smoothen("datacards/May25_2014_195fb/", "MEM_New_rec_std_sb",  cats[j], hists[i] , "JER");
      smoothen("datacards/May25_2014_195fb/", "MEM_New_rec_std_sb",  cats[j], hists[i] , "JEC");
    }
  }

}
