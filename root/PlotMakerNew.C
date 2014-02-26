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

#include "TIterator.h"

#include "TRandom3.h"

typedef TMatrixT<double> TMatrixD;

#define RUNONDATA 0



pair<double,double> getMaxValue( TH1F* hMassProb){
 
  float est  = -99;
  float prob = -99;
  float a = -99;
  float b = -99;
  float c = -99;

  if( hMassProb->GetEntries()==0 ){
    est  =  0.;
    prob =  0.;
    return make_pair(est, prob);
  }

  
  int maxBin = hMassProb->GetMaximumBin();
  int bD = -999;
  int bC = -999;
  int bU = -999;

  int isRightMostBin = 1;
  for(int k=1; k<=hMassProb->GetNbinsX(); k++){
    if( hMassProb->GetBinContent(maxBin+k)>0 ) isRightMostBin=0;
  }
  int isLeftMostBin = 1;
  for(int k=1; k<=hMassProb->GetNbinsX(); k++){
    if( hMassProb->GetBinContent(maxBin-k)>0 ) isLeftMostBin=0;
  }

  if( (maxBin < hMassProb->GetNbinsX() && maxBin > 1) && !isRightMostBin && !isLeftMostBin ){
    //cout << "Center" << endl;
    bC = maxBin;
    for(int k=1; k<=hMassProb->GetNbinsX() && bD==-999; k++){
      if( hMassProb->GetBinContent(maxBin-k)>0 ){
	bD = maxBin-k;	
	////cout << "Found bin " << bD << " with " <<  hMassProb->GetBinContent(bD) << endl;
      }
    }
    for(int k=1; k<=hMassProb->GetNbinsX() && bU==-999; k++){
      if( hMassProb->GetBinContent(maxBin+k)>0 ){
	bU = maxBin+k;	
	////cout << "Found bin " << bU << " with " <<  hMassProb->GetBinContent(bU) << endl;	
      }
    }   
  }
  else if( (maxBin == hMassProb->GetNbinsX()) || isRightMostBin){
    //cout << "Right" << endl;
    bU = maxBin; 
    for(int k=1; k<=hMassProb->GetNbinsX() && bC==-999; k++){
      if( hMassProb->GetBinContent(maxBin-k)>0 ) bC = maxBin-k;	
    }
    for(int k=1; k<=hMassProb->GetNbinsX() && bD==-999; k++){
      if( hMassProb->GetBinContent(bC-k)>0 ) bD = bC-k;	
    }
  }
  else if( (maxBin == 1) || isLeftMostBin){
    //cout << "Left" << endl;
    bD = maxBin; 
    for(int k=1; k<=hMassProb->GetNbinsX() && bC==-999; k++){
      if( hMassProb->GetBinContent(maxBin+k)>0 ) bC = maxBin+k;	
    }
    for(int k=1; k<=hMassProb->GetNbinsX() && bU==-999; k++){
      if( hMassProb->GetBinContent(bC+k)>0 ) bU = bC+k;	
    }
  }
  else{
    //cout << "Unknown" << endl;
    est  =  hMassProb->GetBinCenter  (hMassProb->GetMaximumBin());
    prob =  hMassProb->GetBinContent (hMassProb->GetMaximumBin());
    //cout << "M=" << est << endl;
    return make_pair(est, prob);
  }

  if( bD==-999 || bU==-999 || bC==-999){
    //cout << "Problems" << endl;   
    //cout << bD << "," << bC << "," << bU << endl;
    est  =  hMassProb->GetBinCenter  (hMassProb->GetMaximumBin());
    prob =  hMassProb->GetBinContent (hMassProb->GetMaximumBin());
    //cout << "M=" << est << endl;
    return make_pair(est, prob);
  }

  double xD =  hMassProb->GetBinCenter (bD);
  double xC =  hMassProb->GetBinCenter (bC);
  double xU =  hMassProb->GetBinCenter (bU);
  double yD =  hMassProb->GetBinContent(bD);
  double yC =  hMassProb->GetBinContent(bC);
  double yU =  hMassProb->GetBinContent(bU);
  
  TMatrixD A(3,3);
  const double elements[9] = 
    { 1, xD,  xD*xD,
      1, xC,  xC*xC,
      1, xU,  xU*xU 
    };
  A.SetMatrixArray(elements, "C");
  TMatrixD AInv(3,3);
  double det;
  AInv = A.Invert(&det);
  
  TMatrixD Y(3,1);
  const double yPos[3] = 
    { yD, yC, yU
    };
  Y.SetMatrixArray(yPos, "C");
  
  TMatrixD C(3,1);
  const double dummy[3] = 
    { 1., 1., 1.
    };
  C.SetMatrixArray(dummy,"C");
  C.Mult(AInv,Y);
  
  a = C(2,0);
  b = C(1,0);
  c = C(0,0);
  
  est  = -b/2/a ;
  prob = a*est*est + b*est + c;

  if(est<0){
    est  =  hMassProb->GetBinCenter  (hMassProb->GetMaximumBin());
    prob =  hMassProb->GetBinContent (hMassProb->GetMaximumBin());
    return make_pair(est, prob);  
  }
  
  return make_pair(est, prob);
  
}






void plot_category(string type = "SL", 
		   string cat  = "cat1",
		   string tag = "New_rec_std_byLLR_bj",
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

  THStack* aStack = new THStack("aStack","CMS preliminary #sqrt{s}=8 TeV, L=19.04 fb^{-1};  L_{S}/(L_{S}+L_{B}) ; events ");

  vector<string> samples;
  samples.push_back("SingleT");
  samples.push_back("TTV");
  samples.push_back("TTJetsHFbb");
  samples.push_back("TTJetsHFb");
  samples.push_back("TTJetsLF");
  samples.push_back("TTH125");
  if(RUNONDATA)  samples.push_back("data_obs");


  TH1F* hS    = 0;
  TH1F* hData = 0;
  TH1F* hErr  = 0;

  TFile* f = TFile::Open(("datacards/Feb24_2014/"+type+"_"+tag+".root").c_str());
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
    
    if( hname.find("csvUp")!=string::npos && hname.find("TTH125")==string::npos && hname.find("data")==string::npos ){
      if( h_csvUp==0 ) h_csvUp = (TH1F*)hist->Clone("csvUp");
      else             h_csvUp->Add( hist );
      cout << "Merge csvUp" << endl; 
    }

    if( hname.find("JECUp")!=string::npos && hname.find("TTH125")==string::npos && hname.find("data")==string::npos ){
      if( h_JECUp==0 ) h_JECUp = (TH1F*)hist->Clone("JECUp");
      else             h_JECUp->Add( hist );	    	
      cout << "Merge JECUp" << endl; 
    }

    if( hname.find("JERUp")!=string::npos && hname.find("TTH125")==string::npos && hname.find("data")==string::npos ){
      if( h_JERUp==0 ) h_JERUp = (TH1F*)hist->Clone("JERUp");
      else             h_JERUp->Add( hist );
      cout << "Merge JERUp" << endl; 	    	
    }

    if( hname.find("lumiUp")!=string::npos && hname.find("TTH125")==string::npos && hname.find("data")==string::npos ){
      if( h_lumiUp==0 ) h_lumiUp = (TH1F*)hist->Clone("lumiUp");
      else              h_lumiUp->Add( hist );	    	
    }

    if( hname.find("pdf_ggUp")!=string::npos && hname.find("TTH125")==string::npos && hname.find("data")==string::npos ){
      if( h_pdf_ggUp==0 ) h_pdf_ggUp = (TH1F*)hist->Clone("pdf_ggUp");
      else                h_pdf_ggUp->Add( hist );	    	
      cout << "Merge lumiUp" << endl; 
    }

    if( hname.find("QCDscale_TTbbUp")!=string::npos && hname.find("TTH125")==string::npos && hname.find("data")==string::npos ){
      if( h_QCDscale_TTbbUp==0 ) h_QCDscale_TTbbUp = (TH1F*)hist->Clone("QCDscale_TTbbUp");
      else                h_QCDscale_TTbbUp->Add( hist );	    	
      cout << "Merge QCDscale_TTbbUp" << endl; 
    }

   if( hname.find("QCDscale_TTbUp")!=string::npos && hname.find("TTH125")==string::npos && hname.find("data")==string::npos ){
      if( h_QCDscale_TTbUp==0 ) h_QCDscale_TTbUp = (TH1F*)hist->Clone("QCDscale_TTbUp");
      else             h_QCDscale_TTbUp->Add( hist );	  
      cout << "Merge QCDscale_TTbUp" << endl;   	
    }

   if( hname.find("QCDscale_TTJetsHFUp")!=string::npos && hname.find("TTH125")==string::npos && hname.find("data")==string::npos ){
      if( h_QCDscale_TTJetsHFUp==0 ) h_QCDscale_TTJetsHFUp = (TH1F*)hist->Clone("QCDscale_TTJetsHFUp");
      else             h_QCDscale_TTJetsHFUp->Add( hist );
      cout << "Merge QCDscale_TTJetsHFUp" << endl; 	    	
    }
   
   if( hname.find("QCDscale_TTJetsLFUp")!=string::npos && hname.find("TTH125")==string::npos && hname.find("data")==string::npos ){
      if( h_QCDscale_TTJetsLFUp==0 ) h_QCDscale_TTJetsLFUp = (TH1F*)hist->Clone("QCDscale_TTJetsLFUp");
      else             h_QCDscale_TTJetsLFUp->Add( hist );
      cout << "Merge QCDscale_TTJetsLFUp" << endl; 	    		    	
    }

   if( hname.find("Norm_TTVUp")!=string::npos && hname.find("TTH125")==string::npos && hname.find("data")==string::npos ){
      if( h_Norm_TTVUp==0 ) h_Norm_TTVUp = (TH1F*)hist->Clone("Norm_TTVUp");
      else             h_Norm_TTVUp->Add( hist );
      cout << "Merge Norm_TTVUp" << endl; 	    		    	
    }

   if( hname.find("Norm_SingleTUp")!=string::npos && hname.find("TTH125")==string::npos && hname.find("data")==string::npos ){
      if( h_Norm_SingleTUp==0 ) h_Norm_SingleTUp = (TH1F*)hist->Clone("Norm_SingleTUp");
      else             h_Norm_SingleTUp->Add( hist );
      cout << "Merge Norm_SingleTUp" << endl; 	    		    		    	
    }

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
    if( h_csvUp ) syst_b2 += (hErr->GetBinContent(b)-h_csvUp->GetBinContent(b))*(hErr->GetBinContent(b)-h_csvUp->GetBinContent(b));
    if( h_JECUp ) syst_b2 += (hErr->GetBinContent(b)-h_JECUp->GetBinContent(b))*(hErr->GetBinContent(b)-h_JECUp->GetBinContent(b));
    if( h_JERUp ) syst_b2 += (hErr->GetBinContent(b)-h_JERUp->GetBinContent(b))*(hErr->GetBinContent(b)-h_JERUp->GetBinContent(b));
    if( h_lumiUp ) syst_b2 += (hErr->GetBinContent(b)-h_lumiUp->GetBinContent(b))*(hErr->GetBinContent(b)-h_lumiUp->GetBinContent(b));
    if( h_pdf_ggUp ) syst_b2 += (hErr->GetBinContent(b)-h_pdf_ggUp->GetBinContent(b))*(hErr->GetBinContent(b)-h_pdf_ggUp->GetBinContent(b));
    if( h_QCDscale_TTbbUp ) syst_b2 += (hErr->GetBinContent(b)-h_QCDscale_TTbbUp->GetBinContent(b))*(hErr->GetBinContent(b)-h_QCDscale_TTbbUp->GetBinContent(b));
    if( h_QCDscale_TTbUp ) syst_b2 += (hErr->GetBinContent(b)-h_QCDscale_TTbUp->GetBinContent(b))*(hErr->GetBinContent(b)-h_QCDscale_TTbUp->GetBinContent(b));
    if( h_QCDscale_TTJetsHFUp ) syst_b2 += (hErr->GetBinContent(b)-h_QCDscale_TTJetsHFUp->GetBinContent(b))*(hErr->GetBinContent(b)-h_QCDscale_TTJetsHFUp->GetBinContent(b));
    if( h_QCDscale_TTJetsLFUp ) syst_b2 += (hErr->GetBinContent(b)-h_QCDscale_TTJetsLFUp->GetBinContent(b))*(hErr->GetBinContent(b)-h_QCDscale_TTJetsLFUp->GetBinContent(b));
    if( h_Norm_TTVUp ) syst_b2 += (hErr->GetBinContent(b)-h_Norm_TTVUp->GetBinContent(b))*(hErr->GetBinContent(b)-h_Norm_TTVUp->GetBinContent(b));
    if( h_Norm_SingleTUp ) syst_b2 += (hErr->GetBinContent(b)-h_Norm_SingleTUp->GetBinContent(b))*(hErr->GetBinContent(b)-h_Norm_SingleTUp->GetBinContent(b));

    hErr->SetBinError(b, sqrt(hErr->GetBinError(b)*hErr->GetBinError(b) + syst_b2) );
  }


  hErr->GetYaxis()->SetTitle("Events");
  if( tag.find("_sb") == string::npos )
    hErr->GetXaxis()->SetTitle("P_{s/b}");
  if( tag.find("_bj") == string::npos )
    hErr->GetXaxis()->SetTitle("P_{b/j}");

  if(RUNONDATA==0)
    hErr->SetTitle("Simulation #sqrt{s}=8 TeV, L=19.04 fb^{-1}");
  else
    hErr->SetTitle("CMS Preliminary #sqrt{s}=8 TeV, L=19.04 fb^{-1}");
  hErr->SetTitleSize  (0.04,"X");
  hErr->SetTitleOffset(0.95,"X");
  float max =  TMath::Max(float(hErr->GetMaximum()*1.35),(hData!=0 ? hData->GetMaximum()*1.35 : -1.));
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
  leg->AddEntry(hS, "signal x 5", "L");
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
    c1->SaveAs(  ("Plots/Feb24_2014/Plot_"+cat+"_"+fname+"_"+tag+".png").c_str() );
  }

}


void plotAll(){


/////////////////////////////////////////////


  plot_category( "MEM", 
		 "cat1",
		 "New_rec_reg_byLLR_sb",
		 "Cat. 1 (4b2j W-tag)",
		 "",
		 0, 
		 1.0, 1.0, 1.0
		 );

  plot_category( "MEM", 
		 "cat2",
		 "New_rec_reg_byLLR_sb",
		 "Cat. 2 (4b2j !W-tag)",
		 "",
		 0, 
		 1.0, 1.0, 1.0
		 );


  plot_category( "MEM", 
		 "cat3",
		 "New_rec_reg_byLLR_sb",
		 "Cat. 3 (4b1j cs-tag)",
		 "",
		 0, 
		 1.0, 1.0, 1.0
		 );

  plot_category( "MEM", 
		 "cat4",
		 "New_rec_reg_byLLR_sb",
		 "Cat. 4 (4b1j !cs-tag)",
		 "",
		 0, 
		 1.0, 1.0, 1.0
		 );

  plot_category( "MEM", 
		 "cat5",
		 "New_rec_reg_byLLR_sb",
		 "Cat. 5 (4b3j)",
		 "",
		 0, 
		 1.0, 1.0, 1.0
		 );

  plot_category( "MEM", 
		 "cat6",
		 "New_rec_reg_byLLR_sb",
		 "Cat. 6 (4b)",
		 "",
		 0, 
		 1.0, 1.0, 1.0
		 );


  return;
  ////////////////////////////////////////
  
  plot_category( "MEM", 
		 "cat1",
		 "New_rec_std_byLLR_sb_nb",
		 "Cat. 1 (4b2j W-tag)",
		 "Scale",
		 0, 
		 1.6, 1.6, 1.0
		 );

  plot_category( "MEM", 
		 "cat2",
		 "New_rec_std_byLLR_sb_nb",
		 "Cat. 2 (4b2j !W-tag)",
		 "Scale",
		 0, 
		 1.6, 1.6, 1.0
		 );


  plot_category( "MEM", 
		 "cat3",
		 "New_rec_std_byLLR_sb_nb",
		 "Cat. 3 (4b1j cs-tag)",
		 "Scale",
		 0, 
		 1.6, 1.6, 1.0
		 );

  plot_category( "MEM", 
		 "cat4",
		 "New_rec_std_byLLR_sb_nb",
		 "Cat. 4 (4b1j !cs-tag)",
		 "Scale",
		 0, 
		 1.6, 1.6, 1.0
		 );

  plot_category( "MEM", 
		 "cat5",
		 "New_rec_std_byLLR_sb_nb",
		 "Cat. 5 (4b3j)",
		 "Scale",
		 0, 
		 1.6, 1.6, 1.0
		 );

  plot_category( "MEM", 
		 "cat6",
		 "New_rec_std_byLLR_sb_nb",
		 "Cat. 6 (4b)",
		 "Scale",
		 0, 
		 1.6, 1.6, 1.0
		 );


  /////////////////////////////////////////////

  plot_category( "MEM", 
		 "cat1",
		 "New_rec_std_byLLR_bj",
		 "Cat. 1 (4b2j W-tag)",
		 "Scale",
		 0, 
		 1.6, 1.6, 1.0
		 );

  plot_category( "MEM", 
		 "cat2",
		 "New_rec_std_byLLR_bj",
		 "Cat. 2 (4b2j !W-tag)",
		 "Scale",
		 0, 
		 1.6, 1.6, 1.0
		 );


  plot_category( "MEM", 
		 "cat3",
		 "New_rec_std_byLLR_bj",
		 "Cat. 3 (4b1j cs-tag)",
		 "Scale",
		 0, 
		 1.6, 1.6, 1.0
		 );

  plot_category( "MEM", 
		 "cat4",
		 "New_rec_std_byLLR_bj",
		 "Cat. 4 (4b1j !cs-tag)",
		 "Scale",
		 0, 
		 1.6, 1.6, 1.0
		 );

  plot_category( "MEM", 
		 "cat5",
		 "New_rec_std_byLLR_bj",
		 "Cat. 5 (4b3j)",
		 "Scale",
		 0, 
		 1.6, 1.6, 1.0
		 );

  plot_category( "MEM", 
		 "cat6",
		 "New_rec_std_byLLR_bj",
		 "Cat. 6 (4b)",
		 "Scale",
		 0, 
		 1.6, 1.6, 1.0
		 );


  /////////////////////////////////////////////


  plot_category( "MEM", 
		 "cat1",
		 "New_rec_std_byLLR_sb_nb",
		 "Cat. 1 (4b2j W-tag)",
		 "",
		 0, 
		 1.0, 1.0, 1.0
		 );

  plot_category( "MEM", 
		 "cat2",
		 "New_rec_std_byLLR_sb_nb",
		 "Cat. 2 (4b2j !W-tag)",
		 "",
		 0, 
		 1.0, 1.0, 1.0
		 );


  plot_category( "MEM", 
		 "cat3",
		 "New_rec_std_byLLR_sb_nb",
		 "Cat. 3 (4b1j cs-tag)",
		 "",
		 0, 
		 1.0, 1.0, 1.0
		 );

  plot_category( "MEM", 
		 "cat4",
		 "New_rec_std_byLLR_sb_nb",
		 "Cat. 4 (4b1j !cs-tag)",
		 "",
		 0, 
		 1.0, 1.0, 1.0
		 );

  plot_category( "MEM", 
		 "cat5",
		 "New_rec_std_byLLR_sb_nb",
		 "Cat. 5 (4b3j)",
		 "",
		 0, 
		 1.0, 1.0, 1.0
		 );

  plot_category( "MEM", 
		 "cat6",
		 "New_rec_std_byLLR_sb_nb",
		 "Cat. 6 (4b)",
		 "",
		 0, 
		 1.0, 1.0, 1.0
		 );


  /////////////////////////////////////////////


  plot_category( "MEM", 
		 "cat1",
		 "New_rec_std_byLLR_bj",
		 "Cat. 1 (4b2j W-tag)",
		 "",
		 0, 
		 1.0, 1.0, 1.0
		 );

  plot_category( "MEM", 
		 "cat2",
		 "New_rec_std_byLLR_bj",
		 "Cat. 2 (4b2j !W-tag)",
		 "",
		 0, 
		 1.0, 1.0, 1.0
		 );


  plot_category( "MEM", 
		 "cat3",
		 "New_rec_std_byLLR_bj",
		 "Cat. 3 (4b1j cs-tag)",
		 "",
		 0, 
		 1.0, 1.0, 1.0
		 );

  plot_category( "MEM", 
		 "cat4",
		 "New_rec_std_byLLR_bj",
		 "Cat. 4 (4b1j !cs-tag)",
		 "",
		 0, 
		 1.0, 1.0, 1.0
		 );

  plot_category( "MEM", 
		 "cat5",
		 "New_rec_std_byLLR_bj",
		 "Cat. 5 (4b3j)",
		 "",
		 0, 
		 1.0, 1.0, 1.0
		 );

  plot_category( "MEM", 
		 "cat6",
		 "New_rec_std_byLLR_bj",
		 "Cat. 6 (4b)",
		 "",
		 0, 
		 1.0, 1.0, 1.0
		 );



 
  
  //\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
 
 
  plot_category( "MEM", 
		 "cat1",
		 "New_rec_std_byCSV_sb_nb",
		 "Cat. 1 (4b2j W-tag)",
		 "Scale",
		 0, 
		 1.6, 1.6, 1.0
		 );

  plot_category( "MEM", 
		 "cat2",
		 "New_rec_std_byCSV_sb_nb",
		 "Cat. 2 (4b2j !W-tag)",
		 "Scale",
		 0, 
		 1.6, 1.6, 1.0
		 );


  plot_category( "MEM", 
		 "cat3",
		 "New_rec_std_byCSV_sb_nb",
		 "Cat. 3 (4b1j cs-tag)",
		 "Scale",
		 0, 
		 1.6, 1.6, 1.0
		 );

  plot_category( "MEM", 
		 "cat4",
		 "New_rec_std_byCSV_sb_nb",
		 "Cat. 4 (4b1j !cs-tag)",
		 "Scale",
		 0, 
		 1.6, 1.6, 1.0
		 );

  plot_category( "MEM", 
		 "cat5",
		 "New_rec_std_byCSV_sb_nb",
		 "Cat. 5 (4b3j)",
		 "Scale",
		 0, 
		 1.6, 1.6, 1.0
		 );

  plot_category( "MEM", 
		 "cat6",
		 "New_rec_std_byCSV_sb_nb",
		 "Cat. 6 (4b)",
		 "Scale",
		 0, 
		 1.6, 1.6, 1.0
		 );


  /////////////////////////////////////////////

  plot_category( "MEM", 
		 "cat1",
		 "New_rec_std_byCSV_bj",
		 "Cat. 1 (4b2j W-tag)",
		 "Scale",
		 0, 
		 1.6, 1.6, 1.0
		 );

  plot_category( "MEM", 
		 "cat2",
		 "New_rec_std_byCSV_bj",
		 "Cat. 2 (4b2j !W-tag)",
		 "Scale",
		 0, 
		 1.6, 1.6, 1.0
		 );


  plot_category( "MEM", 
		 "cat3",
		 "New_rec_std_byCSV_bj",
		 "Cat. 3 (4b1j cs-tag)",
		 "Scale",
		 0, 
		 1.6, 1.6, 1.0
		 );

  plot_category( "MEM", 
		 "cat4",
		 "New_rec_std_byCSV_bj",
		 "Cat. 4 (4b1j !cs-tag)",
		 "Scale",
		 0, 
		 1.6, 1.6, 1.0
		 );

  plot_category( "MEM", 
		 "cat5",
		 "New_rec_std_byCSV_bj",
		 "Cat. 5 (4b3j)",
		 "Scale",
		 0, 
		 1.6, 1.6, 1.0
		 );

  plot_category( "MEM", 
		 "cat6",
		 "New_rec_std_byCSV_bj",
		 "Cat. 6 (4b)",
		 "Scale",
		 0, 
		 1.6, 1.6, 1.0
		 );


  /////////////////////////////////////////////


  plot_category( "MEM", 
		 "cat1",
		 "New_rec_std_byCSV_sb_nb",
		 "Cat. 1 (4b2j W-tag)",
		 "",
		 0, 
		 1.0, 1.0, 1.0
		 );

  plot_category( "MEM", 
		 "cat2",
		 "New_rec_std_byCSV_sb_nb",
		 "Cat. 2 (4b2j !W-tag)",
		 "",
		 0, 
		 1.0, 1.0, 1.0
		 );


  plot_category( "MEM", 
		 "cat3",
		 "New_rec_std_byCSV_sb_nb",
		 "Cat. 3 (4b1j cs-tag)",
		 "",
		 0, 
		 1.0, 1.0, 1.0
		 );

  plot_category( "MEM", 
		 "cat4",
		 "New_rec_std_byCSV_sb_nb",
		 "Cat. 4 (4b1j !cs-tag)",
		 "",
		 0, 
		 1.0, 1.0, 1.0
		 );

  plot_category( "MEM", 
		 "cat5",
		 "New_rec_std_byCSV_sb_nb",
		 "Cat. 5 (4b3j)",
		 "",
		 0, 
		 1.0, 1.0, 1.0
		 );

  plot_category( "MEM", 
		 "cat6",
		 "New_rec_std_byCSV_sb_nb",
		 "Cat. 6 (4b)",
		 "",
		 0, 
		 1.0, 1.0, 1.0
		 );

  /////////////////////////////////////////////

  
  plot_category( "MEM", 
		 "cat1",
		 "New_rec_std_byCSV_bj",
		 "Cat. 1 (4b2j W-tag)",
		 "",
		 0, 
		 1.0, 1.0, 1.0
		 );

  plot_category( "MEM", 
		 "cat2",
		 "New_rec_std_byCSV_bj",
		 "Cat. 2 (4b2j !W-tag)",
		 "",
		 0, 
		 1.0, 1.0, 1.0
		 );


  plot_category( "MEM", 
		 "cat3",
		 "New_rec_std_byCSV_bj",
		 "Cat. 3 (4b1j cs-tag)",
		 "",
		 0, 
		 1.0, 1.0, 1.0
		 );

  plot_category( "MEM", 
		 "cat4",
		 "New_rec_std_byCSV_bj",
		 "Cat. 4 (4b1j !cs-tag)",
		 "",
		 0, 
		 1.0, 1.0, 1.0
		 );

  plot_category( "MEM", 
		 "cat5",
		 "New_rec_std_byCSV_bj",
		 "Cat. 5 (4b3j)",
		 "",
		 0, 
		 1.0, 1.0, 1.0
		 );

  plot_category( "MEM", 
		 "cat6",
		 "New_rec_std_byCSV_bj",
		 "Cat. 6 (4b)",
		 "",
		 0, 
		 1.0, 1.0, 1.0
		 );


}




void plot_limit(string version = "_New", 
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
  c1->SetLogy(!doMLFit);
  
  TLegend* leg = new TLegend(0.121,0.154,0.32, 0.29, NULL,"brNDC");
  leg->SetBorderSize(0);
  leg->SetTextSize(0.04);
  leg->SetFillColor(0);


  float X[]        = {0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5 };
  float expY[]     = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float obsY[]     = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float expY1sL[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float expY1sH[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float expY2sL[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  float expY2sH[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};

  float expXs[]  = {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5};
  float expYs[]  = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};


  vector<string> categories;        vector<string> names;
  categories.push_back("MEM_cat1"+version);  names.push_back("Cat 1");
  categories.push_back("MEM_cat2"+version);  names.push_back("Cat 2");
  categories.push_back("MEM_cat3"+version);  names.push_back("Cat 3");
  categories.push_back("MEM_cat4"+version);  names.push_back("Cat 4");
  categories.push_back("MEM_cat5"+version);  names.push_back("Cat 5");
  categories.push_back("MEM_cat6"+version);  names.push_back("Cat 6"); 
  categories.push_back("MEM_SL"+version);       names.push_back("SL comb.");
  //categories.push_back("DL_cat7");names.push_back();
  categories.push_back("MEM_DL"+version);       names.push_back("DL comb."); 
  categories.push_back("MEM_COMB"+version);     names.push_back("All comb.");

  int nBins = categories.size();

  for( int b = 0; b < nBins; b++){

    TFile* f = doMLFit ? 
      TFile::Open(("datacards/higgsCombine"+categories[b]+".MaxLikelihoodFit.mH120.root").c_str()) :
      TFile::Open(("datacards/higgsCombine"+categories[b]+".Asymptotic.mH120.root").c_str());
    if( f==0 ) continue;
    
    Double_t r;
    TTree* limit = (TTree*)f->Get("limit");
    limit->SetBranchAddress("limit",&r);

    for(int k = 0 ; k< limit->GetEntries() ; k++){
      limit->GetEntry(k);
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


    cout << categories[b] << ": r<" <<  expY[b] << " @ 95% CL" << endl;

    expY1sH[b] = TMath::Abs(expY1sH[b]-expY[b]);
    expY1sL[b] = TMath::Abs(expY1sL[b]-expY[b]);
    expY2sH[b] = TMath::Abs(expY2sH[b]-expY[b]);
    expY2sL[b] = TMath::Abs(expY2sL[b]-expY[b]);


  }


  TMultiGraph *mg = new TMultiGraph();
  mg->SetTitle("Simulation #sqrt{s}=8 TeV, L=19.04 fb^{-1}");

  TGraphAsymmErrors* observed = new TGraphAsymmErrors(10, X, obsY, expXs ,expXs,  expYs,   expYs);
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

  line->Draw("SAME");

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
  gPad->Modified();
  mg->GetXaxis()->Set(nBins,0,nBins);

  //mg->GetXaxis()->SetRange(-1,11);

  for( int b = 0; b < nBins; b++){
    mg->GetXaxis()->SetBinLabel(b+1, names[b].c_str() );
  }

  mg->GetYaxis()->SetTitleOffset(0.97);
  mg->SetMinimum( minY);
  mg->SetMaximum( maxY );
  mg->GetXaxis()->SetTitle("");
  if(!doMLFit) 
    mg->GetYaxis()->SetTitle("95% CL upper limit on #mu = #sigma/#sigma_{SM}");
  else
    mg->GetYaxis()->SetTitle("ML fit of #mu = #sigma/#sigma_{SM}");

  if( doExpOnly==0 && !doMLFit ){
    leg->AddEntry( observed, "MEM observed", "P");
    leg->AddEntry( expected, "MEM expected (post-fit)", "P");
  }
  else if(!doMLFit)
    leg->AddEntry( expected, "MEM expected (pre-fit)", "P");
  else 
    leg->AddEntry( expected, "observed", "P");

  if(compare){
    leg->AddEntry( lineTTH, "HIG-13-019 (post-fit)", "L");
    leg->AddEntry( lineML,  "HIG-13-020 (post-fit)", "L");
  }

  leg->Draw();

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

  //pt->Draw();

  //c1->Modified();
  //c1->Draw();

  TLine* del = new TLine(6., maxY, 6., 0. );
  del->SetLineWidth(3);
  del->SetLineStyle(kSolid);
  del->SetLineColor(kBlack);
  del->Draw("SAME");

  TLine* del_ = new TLine(5., maxY, 5., 0. );
  del_->SetLineWidth(2);
  del_->SetLineStyle(kDotted);
  del_->SetLineColor(kBlack);
  //del_->Draw("SAME");

  TLine* del2 = new TLine(8., maxY, 8., 0. );
  del2->SetLineWidth(3);
  del2->SetLineStyle(kSolid);
  del2->SetLineColor(kBlack);
  del2->Draw("SAME");

  TLine* del2_ = new TLine(7., maxY, 7., 0. );
  del2_->SetLineWidth(2);
  del2_->SetLineStyle(kDotted);
  del2_->SetLineColor(kBlack);
  //del2_->Draw("SAME");

  if(1){
    c1->SaveAs("Plots/Feb24_2014/Limits"+TString(version.c_str())+".png");
    //c1->SaveAs("datacards/Limits"+TString(version.c_str())+".pdf");
  }

}


void plot_limitAll(){

  // exp limits
  plot_limit("_New_rec_std_byCSV_sb",    1, 0, 0.8, 40, 1);
  plot_limit("_New_rec_std_byLLR_sb",    1, 0, 0.8, 40, 1);
  plot_limit("_New_rec_reg_byLLR_sb",    1, 0, 0.8, 40, 1);
  plot_limit("_New_rec_reg_byLLR_sb_scaleUp_TTbb",    1, 0, 0.8, 40, 1);
  plot_limit("_New_rec_std_byCSV_sb_nb", 1, 0, 0.8, 40, 1);
  plot_limit("_New_rec_std_byLLR_sb_nb", 1, 0, 0.8, 40, 1);


  // ML fit
  plot_limit("_New_rec_std_byLLR_bj",    0, 1, 0., 4, 0);
  plot_limit("_New_rec_std_byCSV_bj",    0, 1, 0., 4, 0);
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
    hU->SetTitle("Simulation #sqrt{s}=8 TeV, L=19.04 fb^{-1}");
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
    hD->SetTitle("Simulation #sqrt{s}=8 TeV, L=19.04 fb^{-1}");
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




void plot_BTag( TString extraname = "",	       
	       int flag = 0,
	       int plotDistr = 0){

  int nmax = 2000;
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
    leg->SetHeader("Single-lepton, N_{jet}=6");
  else if(flag>=4 && flag<8)   
    leg->SetHeader("Single-lepton, N_{jet}=5");
  else if(flag>=8 && flag<12)   
    leg->SetHeader("Double-lepton, N_{jet}=4");

  float kappa1 = 1;

  float min_bb = 0; 

  TString titleDistr("Simulation #sqrt{s}=8 TeV; likelihood ratio; units");

  TH1F* hS = new TH1F("hS",titleDistr, nmax, 0, 1); hS->SetLineColor(kRed) ; hS->SetLineWidth(2); hS->SetFillStyle(3004); hS->SetFillColor(kRed); 
  TH1F* hB = new TH1F("hB",titleDistr, nmax, 0, 1); hB->SetLineColor(kBlue); hB->SetLineWidth(2); hB->SetFillStyle(3005); hB->SetFillColor(kBlue); 
  
  float xM[3], yM[3];

  float xWP[1], yWP[1];
  xWP[0] = -99;
  yWP[0] = -99;

  int oldFiles = 0;
  vector<string> files;
  string path = oldFiles ? "~/public/BTagging/Nov04_2013/" : "~/public/BTagging/Nov21_2013/";
  string name =  "New" ;
  string jets =  "J" ;

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
    
    int counter2M = 0;
    int counter3M = 0;
    int counter4M = 0;

    Long64_t nentries = t->GetEntries(); 
    int counter = 0;

    for (Long64_t i = 0; i < nentries ; i++){

      t->GetEntry(i);
      
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

    if( flag<=4 &&  hS->GetBinCenter(b)>0.974 &&   hS->GetBinCenter(b)<0.976 && xWP[0]<0 ){
      xWP[0] = intS;
      yWP[0] = intB;
      cout << "\e[1;32m" << intS << " --- " << intB << "\e[0m" << endl;
    }
    if( flag>4 && flag<=8 &&  hS->GetBinCenter(b)>0.989 &&   hS->GetBinCenter(b)<0.991 && xWP[0]<0 ){
      xWP[0] = intS;
      yWP[0] = intB;
      cout << "\e[1;32m" << intS << " --- " << intB << "\e[0m" << endl;
    }
    if(flag>8 && flag<=11  &&  hS->GetBinCenter(b)>0.952 &&   hS->GetBinCenter(b)<0.954 && xWP[0]<0 ){
      xWP[0] = intS;
      yWP[0] = intB;
      cout << "\e[1;32m" << intS << " --- " << intB << "\e[0m" << endl;
    }

  }
  x[nmax]=0.00;
  y[nmax]=0.00;
  cout << "Eff.=" << x[nmax] << " -- vs -- FkR.=" << y[nmax] << endl;


  ///////////////////////////////////////

  TString title("Simulation #sqrt{s}=8 TeV; selection efficiency in t#bar{t}H ( #it{PITHYA} ); selection efficiency in t#bar{t}+jets ( #it{MadGraph} )");
  if( flag%4==0    ) title = "Simulation #sqrt{s}=8 TeV; selection efficiency in t#bar{t}H ( #it{PITHYA} ); selection efficiency in t#bar{t}+jets ( #it{MadGraph} )";
  if( (flag-1)%4==0) title = "Simulation #sqrt{s}=8 TeV; selection efficiency in t#bar{t}H ( #it{PITHYA} ); selection efficiency in t#bar{t}+jj ( #it{MadGraph} )";
  if( (flag-2)%4==0) title = "Simulation #sqrt{s}=8 TeV; selection efficiency in t#bar{t}H ( #it{PITHYA} ); selection efficiency in t#bar{t}+b ( #it{MadGraph} ) ";
  if( (flag-3)%4==0) title = "Simulation #sqrt{s}=8 TeV; selection efficiency in t#bar{t}H ( #it{PITHYA} ); selection efficiency in t#bar{t}+b#bar{b} ( #it{MadGraph} )";

  TH1F* h = new TH1F("h",title, 100, 0.0, 0.3);  
  if(flag==0){
    h->SetMinimum(0.001);
    h->SetMaximum(0.10);
  }
  if(flag==1){
    h->SetMinimum(0.0005);
    h->SetMaximum(0.10);
  }
  if(flag==2){
    h->SetMinimum(0.001);
    h->SetMaximum(0.20);
  }
  if(flag==3){
    h->SetMinimum(0.04);
    h->SetMaximum(0.30);
  }

  if(flag==4){
    h->SetMinimum(0.001);
    h->SetMaximum(0.10);
  }
  if(flag==5){
    h->SetMinimum(0.0002);
    h->SetMaximum(0.10);
  }
 if(flag==6){
    h->SetMinimum(0.001);
    h->SetMaximum(0.20);
  }
 if(flag==7){
    h->SetMinimum(0.005);
    h->SetMaximum(0.30);
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
   h->SetMaximum(0.30);
 }

  h->Draw();
  TGraph* hROC = new TGraph(nmax+1, x, y);
  hROC->SetLineWidth(2);
  hROC->SetLineColor(kRed);
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
   leg->AddEntry(hROCM, "cut-and-count: 4L3M, 4M, 4M2T","P");
   leg->AddEntry(hROC, "btag-LR","L");
   leg->AddEntry(hWP, "btag-LR WP","P");
   leg->Draw();
  }

  if(plotDistr==1){
    leg->AddEntry(hS, "t#bar{t}H","F");
    if(flag==0 || flag==3 || flag==6)
      leg->AddEntry(hB, "t#bar{t}+jets","F");
    if(flag==1 || flag==4 || flag==7)
      leg->AddEntry(hB, "t#bar{t}+jj","F");
    if(flag==2 || flag==5 || flag==8)
      leg->AddEntry(hB, "t#bar{t}+b#bar{b}","F");

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
    c1->SaveAs("Plots/Feb24_2014/Plot_BTag_ROC_"+extraname+".png");
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


void plot_match( TString fname = "all_rec_std", 
		string cut = "type==0", 
		TString category = "cat0",
		string header = "",
		int doSgn = 1){
  

  string basecut = cut;

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
  c1->SetGrid(1,1);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  TLegend* leg = new TLegend(0.15,0.75,0.40,0.89,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04); 
  leg->SetHeader(header.c_str());

  TFile* f_B = TFile::Open("files/byLLR/MEAnalysisNew_"+fname+"_TTJets.root");
  TTree* tB = (TTree*)f_B->Get("tree");
  TFile* f_S = TFile::Open("files/byLLR/MEAnalysisNew_"+fname+"_TTH125.root");
  TTree* tS = (TTree*)f_S->Get("tree");

  TTree* tree = doSgn ? tS : tB;
  int nmax    = doSgn ? 5 : 5;

  float all         = tree->GetEntries(TCut(cut.c_str()));
  float matchesT    = tree->GetEntries(TCut("matchesT==2")&&TCut(cut.c_str()));
  float matchesW    = tree->GetEntries(TCut("matchesW==2")&&TCut(cut.c_str()));
  float matchesTop  = tree->GetEntries(TCut("matchesW==2 && matchesT==2")&&TCut(cut.c_str()));
  float matchesW1   = tree->GetEntries(TCut("matchesW==1")&&TCut(cut.c_str()));
  float matchesTop1 = tree->GetEntries(TCut("matchesW==1 && matchesT==2")&&TCut(cut.c_str()));
  float matchesH    = tree->GetEntries(TCut("matchesH==2")&&TCut(cut.c_str()));
  float matchesAll  = tree->GetEntries(TCut("matchesH==2 && matchesW==2 && matchesT==2")&&TCut(cut.c_str()));
  float matchesAll1 = tree->GetEntries(TCut("matchesH==2 && matchesW==1 && matchesT==2")&&TCut(cut.c_str()));
  float matchesAll2 = tree->GetEntries(TCut("matchesH==2 && matchesT==2")&&TCut(cut.c_str()));

  TString titleDistr("Simulation #sqrt{s}=8 TeV; ; efficiency");
  TH1F* hS = new TH1F("hS",titleDistr, nmax, 0, 1);
  hS->Sumw2();

  hS->SetBinContent(1, matchesT/all);  
  hS->SetBinError(1, sqrt(matchesT/all*(1-matchesT/all)/all));    
  hS->GetXaxis()->SetBinLabel(1, "b's from top");

  if((string(category.Data())).find("cat1")!=string::npos || 
     (string(category.Data())).find("cat5")!=string::npos){

    hS->SetBinContent(2, matchesW/all);   
    hS->SetBinError  (2, sqrt(matchesW/all*(1-matchesW/all))/all); 
    hS->GetXaxis()->SetBinLabel(2, "q's from top");

    hS->SetBinContent(3, matchesTop/all); 
    hS->SetBinError  (3, sqrt(matchesTop/all*(1-matchesTop/all))/all); 
    hS->GetXaxis()->SetBinLabel(3, "top's");

    if(doSgn){
      hS->SetBinContent(5, matchesAll/all); 
      hS->SetBinError  (5, sqrt(matchesAll/all*(1-matchesAll/all))/all); 
      hS->GetXaxis()->SetBinLabel(5, "top+higgs");
    }
    else{
      hS->SetBinContent(5,0.); 
      hS->SetBinError  (5,0.); 
      hS->GetXaxis()->SetBinLabel(5, "top+higgs");
    }
  }
  else if((string(category.Data())).find("cat2")!=string::npos || 
	  (string(category.Data())).find("cat3")!=string::npos ||
	  (string(category.Data())).find("cat4")!=string::npos){

    hS->SetBinContent(2, matchesW1/all);   
    hS->SetBinError  (2, sqrt(matchesW1/all*(1-matchesW1/all))/all); 
    hS->GetXaxis()->SetBinLabel(2, "q from top");

    hS->SetBinContent(3, matchesTop1/all); 
    hS->SetBinError  (3, sqrt(matchesTop1/all*(1-matchesTop1/all))/all); 
    hS->GetXaxis()->SetBinLabel(3, "top's");

    if(doSgn){
      hS->SetBinContent(5, matchesAll1/all); 
      hS->SetBinError  (5, sqrt(matchesAll1/all*(1-matchesAll1/all))/all); 
      hS->GetXaxis()->SetBinLabel(5, "top+higgs");
    }
    else{
      hS->SetBinContent(5,0.); 
      hS->SetBinError  (5,0.); 
      hS->GetXaxis()->SetBinLabel(5, "top+higgs");
    }

  }
  else if((string(category.Data())).find("cat6")!=string::npos ){
    
    hS->SetBinContent(2, 0.);   
    hS->SetBinError  (2, 0.);
    hS->GetXaxis()->SetBinLabel(2, "q's from top");

    hS->SetBinContent(3, 0.); 
    hS->SetBinError  (3, 0.); 
    hS->GetXaxis()->SetBinLabel(3, "top's");

    if(doSgn){
      hS->SetBinContent(5, matchesAll2/all); 
      hS->SetBinError  (5, sqrt(matchesAll2/all*(1-matchesAll2/all))/all); 
      hS->GetXaxis()->SetBinLabel(5, "top+higgs");
    }
    else{
      hS->SetBinContent(5,0.); 
      hS->SetBinError  (5,0.); 
      hS->GetXaxis()->SetBinLabel(5, "top+higgs");
    }
  }
    hS->SetBinContent(4, matchesH/all); 
    hS->SetBinError  (4, sqrt(matchesH/all*(1-matchesH/all))/all); 
    hS->GetXaxis()->SetBinLabel(4, "b's from higgs");
  

  hS->SetMarkerStyle(kFullCircle);
  hS->SetMarkerColor(kRed);
  hS->SetMarkerSize(1.8);
  hS->SetMinimum(0.);
  hS->SetMaximum(1.3);
  hS->Draw("PE");

  TF1* line = new TF1("line","1",0,1);
  line->SetLineWidth(3);
  line->SetLineStyle(kDashed);
  line->SetLineColor(kBlue);
  line->Draw("SAME");

  TLine* line2 = new TLine(hS->GetBinLowEdge(nmax), 1.3  , hS->GetBinLowEdge(nmax), 0.);
  line2->SetLineWidth(4);
  line2->SetLineStyle(kSolid);
  line2->SetLineColor(kBlack);
  line2->Draw("SAME");

  //leg->AddEntry(hS, "MC matching efficiency","P");
  leg->Draw();

  if(0){
    c1->SaveAs("Plots/Feb24_2014/Plot_Match_"+category+".png");
    delete c1;
  }
}

void plot_matchAll(){

  plot_match("all_rec_std", "syst==0 && type==0", "cat1_s_byLLR", "Cat 1, t#bar{t}H",     1);
  plot_match("all_rec_std", "syst==0 && type==0", "cat1_b_byLLR", "Cat 1, t#bar{t}+jets", 0);

  plot_match("all_rec_std", "syst==0 && type==1", "cat2_s_byLLR", "Cat 2, t#bar{t}H",     1);
  plot_match("all_rec_std", "syst==0 && type==1", "cat2_b_byLLR", "Cat 2, t#bar{t}+jets", 0);

  plot_match("all_rec_std", "syst==0 && type==2 && flag_type2>0", "cat3_s_byLLR", "Cat 3, t#bar{t}H",     1);
  plot_match("all_rec_std", "syst==0 && type==2 && flag_type2>0", "cat3_b_byLLR", "Cat 3, t#bar{t}+jets", 0);

  plot_match("all_rec_std", "syst==0 && type==2 && flag_type2<=0", "cat4_s_byLLR", "Cat 4, t#bar{t}H",     1);
  plot_match("all_rec_std", "syst==0 && type==2 && flag_type2<=0", "cat4_b_byLLR", "Cat 4, t#bar{t}+jets", 0);

  plot_match("all_rec_std", "syst==0 && type==3", "cat5_s_byLLR", "Cat 5, t#bar{t}H",     1);
  plot_match("all_rec_std", "syst==0 && type==3", "cat5_b_byLLR", "Cat 5, t#bar{t}+jets", 0);

  plot_match("all_rec_std", "syst==0 && type==6", "cat6_s_byLLR", "Cat 6, t#bar{t}H",     1);
  plot_match("all_rec_std", "syst==0 && type==6", "cat6_b_byLLR", "Cat 6, t#bar{t}+jets", 0);


}



void plot_mass(// secondary name of the input trees 
	       TString fname = "SL", 

	       // selection cut
	       string cut = "type==0", 
	       
	       // a header for the output file name
	       TString category = "cat0",
	       
	       // a name tag to save the plot 
	       TString plotname = "SL6jets",
	       
	       // 0 = stack ; 1 = ttH/ttV/ttjets ; 2 = ttbb/ttb/ttjj 
	       int plotShapes = 0,

	       // luminosity normalization (if 1, samples normalized to 12.1 fb-1)
	       float lumiScale = 20./12.1,
	       
	       // histo binning
	       int nBins =20, float xLow=50, float xHigh=250,
	       
	       // normalize weights ? If 1, xec = M^{-fact1}
	       int normXSec=1, float fact1 = 5,

	       // which matching pattern is a full match ? Take the OR of two patterns
	       int matchPattern = 111111, int matchPattern2 = 999999,

	       // choose which perm to use (LO, NLO, NNLO, NNNLO)
	       int choosePerm = 999,

	       // 0 = mass scan ; 1 = MEM
	       int analysis = 0,

	       // parameters for the ME weight
	       float p0 = 1., float p1 = 1., float p2 = 1.
	       ){
  
 
  // setup the style
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

  // canvas
  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);

  // legend
  TLegend* leg = new TLegend(0.59,0.44,0.84,0.90,NULL,"brNDC");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetFillColor(10);
  leg->SetTextSize(0.04); 
  leg->SetHeader(category);


  TRandom3* ran = new TRandom3();


  // inclusive counter (for ttH)
  int counterH         = 0;

  // counter incr. if all quarks in acceptance (for ttH)
  int counterHQuarks   = 0;

  // counter incr. if all quarks in acceptanc & best comb. is matched (for ttH)
  int counterHMatch    = 0;

  // counter incr. if best comb. has the higgs's b's matched (for ttH)
  int counterHMatchAny = 0;
  
  // cross-section norm. (a simple power law)
  TF1* xsec = new TF1("xsec",Form("x^(-%f)", fact1 ), 20, 500);

  // histo labels
  TString histolabel = "Simulation #sqrt{s}=8 TeV, L=19.04 fb^{-1}; M_{H} estimator [GeV] ; events";

  // the stack and its MC stat. error plot
  THStack* aStack = new THStack("aStack",histolabel);
  TH1F* hErr      = 0;

  // the parent histogram
  TH1F* hMass     = new TH1F("hMass", histolabel,nBins, xLow, xHigh);

  // histograms for sepcific processes
  TH1F* hTTH      = new TH1F("hTTH",  histolabel,nBins, xLow, xHigh);
  TH1F* hTTV      = new TH1F("hTTV",  histolabel,nBins, xLow, xHigh);
  TH1F* hTT       = new TH1F("hTT",   histolabel,nBins, xLow, xHigh);
  TH1F* hTTbb     = new TH1F("hTTbb", histolabel,nBins, xLow, xHigh);
  TH1F* hTTbj     = new TH1F("hTTbj", histolabel,nBins, xLow, xHigh);
  TH1F* hTTjj     = new TH1F("hTTjj", histolabel,nBins, xLow, xHigh);

  // histogram filled ev-by-ev with mass-dependent probability
  TH1F* hTmp      = new TH1F("hTmp", "", 500,0,500);
  
  // create structure for errors
  hTTbb->Sumw2();
  hTTbj->Sumw2();
  hTTjj->Sumw2();
  hTTH->Sumw2();
  hTTV->Sumw2();

  // the input samples
  vector<string> samples;
  samples.push_back("TTV");
  samples.push_back("SingleT");
  samples.push_back("DiBoson");
  samples.push_back("TTJetsBB");
  samples.push_back("TTJetsBJ");
  samples.push_back("TTJetsJJ");
  samples.push_back("TTH125");
  //samples.push_back("EWK");

  // main name of the trees
  string name = "New_MHscan";

  // which scan ?
  int doMHscan = 0;

  if( !doMHscan ) name = "New_MTscan";
  
  // version
  string version =  "_v4_gen_std" ;
  cout << "Doing version " << version << " and category " << category << endl;
 
  // input files path
  string inputpath = "gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/Trees/MEM/Nov04_2013/mhscan/";
  inputpath = "./";


  // loop over the samples...
  for( unsigned int s = 0 ; s < samples.size(); s++){

    string sample = samples[s];

    // cut for the plot
    TCut sample_cut(cut.c_str());

    // get a temporary histogram for the process under study
    TH1F* hMass_s = (TH1F*)hMass->Clone(("hMass_"+sample).c_str());
    hMass_s->Reset();
    hMass_s->Sumw2();
    
    // split TTJets into final states
    int color;
    if( sample.find("TTJetsBB")!=string::npos ){
      sample = "TTJets";
      sample_cut = sample_cut && TCut("nSimBs>2 && nMatchSimBs>=2");
      color = 16;
    }
    if( sample.find("TTJetsBJ")!=string::npos ){
      sample = "TTJets";
      sample_cut = sample_cut && TCut("nSimBs>2 && nMatchSimBs<2");
      color = 17;
    }
    if( sample.find("TTJetsJJ")!=string::npos ){
      sample = "TTJets";
      sample_cut = sample_cut && TCut("nSimBs==2");
      color = 18;
    }

    // open the file
    TFile* f = TFile::Open((inputpath+"MEAnalysis"+name+"_"+string(fname.Data())+"_nominal"+version+"_"+sample+".root").c_str());

    if(f==0 || f->IsZombie()){
      cout << "Missing " << sample << " file" << endl;
      return;
    }
    else{
      cout << "Doing sample " << sample << endl;
    }   

    // copy the events form the tree passing the selection cut
    TTree* tFull = (TTree*)f->Get("tree");
    TFile* dummy = new TFile("dummy.root","RECREATE");
    TTree* t = (TTree*)tFull->CopyTree( sample_cut );

    // permutation ME probability vs Higgs mass
    float p_vsMH_s[999];
    float p_vsMT_b[999];

    // permutation btag probability
    float p_tt_bb [999];
    float p_tt_jj [999];

    // num. of permutations
    int nPermut_s;   
    int nPermut_b;

    // num. of mass points scanned
    int nMassPoints;

    // Higgs mass points scanned
    float mH_scan[999];
    float mT_scan[999];
  
    // total number of integrations ( nPermut_s*nMassPoints)
    int nTotInteg_s;
    int nTotInteg_b;
    

    // matching of each permutation
    int perm_to_gen_s[999];
    int perm_to_gen_b[999];

    // the normalization weight
    float weight;

    t->SetBranchAddress("p_vsMH_s",     p_vsMH_s);
    t->SetBranchAddress("p_vsMT_b",     p_vsMT_b);
    t->SetBranchAddress("p_tt_bb",      p_tt_bb);
    t->SetBranchAddress("p_tt_jj",      p_tt_jj);
    t->SetBranchAddress("nPermut_s",    &nPermut_s);
    t->SetBranchAddress("nPermut_b",    &nPermut_b);
    t->SetBranchAddress("nTotInteg_s",  &nTotInteg_s);
    t->SetBranchAddress("nTotInteg_b",  &nTotInteg_b);
    t->SetBranchAddress("nMassPoints",  &nMassPoints);
    t->SetBranchAddress("mH_scan",      mH_scan);
    t->SetBranchAddress("mT_scan",      mT_scan);
    t->SetBranchAddress("weight",       &weight);
    t->SetBranchAddress("perm_to_gen_s",perm_to_gen_s);
    t->SetBranchAddress("perm_to_gen_b",perm_to_gen_b);
    
    // total tree entries
    Long64_t nentries = t->GetEntries(); 
    cout << "Total entries: " << nentries << endl;

    // inclusive counter 
    int counter         = 0;
    
    // counter incr. if all quarks in acceptance
    int counterQuarks   = 0;
    
    // counter incr. if all quarks in acceptanc & best comb. is matched
    int counterMatch    = 0;
    
    // counter incr. if best comb. has the higgs's b's matched
    int counterMatchAny = 0;

    // loop over the events...
    for (Long64_t i = 0; i < nentries ; i++){

      // read the event...
      t->GetEntry(i);

      // increment the inclusive counter
      counter++;

      // choose what scan to perform
      int nPermut      = doMHscan ? nPermut_s     : nPermut_b;
      float* m_scan    = doMHscan ? mH_scan       : mT_scan;
      float* p_vsM     = doMHscan ? p_vsMH_s      : p_vsMT_b;
      int* perm_to_gen = doMHscan ? perm_to_gen_s : perm_to_gen_b;



      // which permutations to consider for the ME analysis
      vector<int> whichpermut;
      int stop    = 0;
      int numiter = 0;
      while(stop!=1 && choosePerm<0){
	
	// take a random permutation
	int newpermut =  int(ran->Uniform(-0.5, nPermut+0.5));

	// if empty, fill the vector with the first permutation and stop if only 1 is needed
	if( whichpermut.size()==0 ){
	  whichpermut.push_back( newpermut );
	  if( whichpermut.size() == abs(choosePerm) ) stop = 1;
	}

	// the vector is already filled, check for replicas
	else{
	  int isThere = 0;
	  for(int ppp = 0 ; ppp < whichpermut.size() ; ppp++)
	    if( newpermut == whichpermut[ppp] ) isThere=1;

	  // is the permutation is new, push it back
	  if(isThere==0){
	    whichpermut.push_back( newpermut );
	    if( whichpermut.size() == abs(choosePerm) ) stop = 1;
	  }
	}

	// count how many iterations have been done (protetion agains infinite loop)
	numiter++;
	if(numiter>2000){
	  cout << "Exceeded max number of tries... exit" << endl;
	  stop = 1;
	}
      }

      //cout << whichpermut.size() << endl;
      //if( whichpermut.size()>0 )  cout << " p = " << whichpermut[0] << endl;

      // mass scan analysis
      if(analysis==0){

	// reset the temporary histogram
	hTmp->Reset();
	
	// reset the per-permutation total probability (summed over MH)
	float perm_prob [nPermut];
	float perm_match[nPermut];      
	for( int perm_it = 0 ; perm_it < nPermut ; perm_it++){
	  perm_prob [perm_it] = 0.;
	  perm_match[perm_it] = 0;
	}

	// a flag for events with all quarks in the acceptance
	int quarks = 0;
	
	// loop over the Higgs mass grid
	for(int mH_it = 0; mH_it<nMassPoints; mH_it++){
	  
	  float mH =  m_scan[mH_it];
	  
	  // change values here to retsrict the search range
	  if( mH<0 || mH>300) continue;
	  
	  // for each mass point, loop over pewrmutations
	  for( int perm_it = 0 ; perm_it<nPermut ; perm_it++){
	    
	    // the ME probability 
	    float ME_prob = p_vsM[mH_it*nPermut + perm_it];
	    
	    // the btag probability
	    float bb_prob =  p_tt_bb[perm_it] ;
	    
	    // normalize the weights
	    float norm    = normXSec ? 1./xsec->Eval( mH ) : 1.0;
	    
	    // matching status of this particular permutation
	    int match     = perm_to_gen[perm_it];

	    // the total probability
	    double p =  ME_prob /** bb_prob*/ * norm;
	    
	    // add this into the per-permutation probability
	    perm_prob [perm_it] += p;
	    perm_match[perm_it] = match;
	    
	    // flag events where at least one permutation matches the pattern
	    if( match == matchPattern ) quarks++;
	  
	  //hTmp->Fill( mH, p );
	  }
	}
    
	// MAX, NMAX, NNMAX, NNNMAX permutations ranked by prob. 
	double maxP       = 0.;
	double nmaxP      = 0.;
	double nnmaxP     = 0.;
	double nnnmaxP    = 0.;
	
	int   maxM       = 0;
	int   nmaxM      = 0;
	int   nnmaxM     = 0;
	int   nnnmaxM    = 0;
	
	int   maxPerm    = -99;
	int   nmaxPerm   = -99;
	int   nnmaxPerm  = -99;
	int   nnnmaxPerm = -99;

      
	// loop over the mass integrated probability and rank the permutations
	for( int perm_it = 0 ; perm_it < nPermut ; perm_it++){
	  
	  // LO
	  if(perm_prob [perm_it] > maxP ){
	    maxP    = perm_prob [perm_it];
	    maxM    = perm_match[perm_it];
	    maxPerm = perm_it;
	  }
	  // NLO
	  if(perm_prob [perm_it] > nmaxP && perm_prob [perm_it]< maxP){
	    nmaxP    = perm_prob [perm_it];
	    nmaxM    = perm_match[perm_it];
	    nmaxPerm = perm_it;
	  }
	  // NNLO
	  if(perm_prob [perm_it] > nnmaxP && perm_prob [perm_it]< nmaxP){
	    nnmaxP    = perm_prob [perm_it];
	    nnmaxM    = perm_match[perm_it];
	    nnmaxPerm = perm_it;
	  }
	  // NNNLO
	  if(perm_prob [perm_it] > nnnmaxP && perm_prob [perm_it]< nnmaxP){
	    nnnmaxP    = perm_prob [perm_it];
	    nnnmaxM    = perm_match[perm_it];
	    nnnmaxPerm = perm_it;
	  }
	}
	
	// choose which permutation to plot (default is maxPerm)
	if( (choosePerm==0 || choosePerm==999) && maxPerm != -99 ){
	  /*  do not do anything ... */
	}

	// if no permutation with prob>0 exist, take the first
	// (it will give 0 prob, and getMaxValue() will return M=0)
	else if( (choosePerm==0 || choosePerm==999) && maxPerm == -99 ){
	  maxPerm = 0;
	}

	
	else if(choosePerm==1){
	  if( nmaxPerm != -99 ){
	    maxPerm = nmaxPerm;
	    maxM    = nmaxM;
	  }
	  else{
	    /*  do not do anything ... */
	  }
	}
	if(choosePerm==2){
	  if( nnmaxPerm != -99 ){
	    maxPerm = nnmaxPerm;
	    maxM    = nnmaxM;
	  }
	  else if( nmaxPerm != -99 ){
	    maxPerm = nmaxPerm;
	    maxM    = nmaxM;
	  }
	  else{
	    /*  do not do anything ... */
	  }
	}
	if(choosePerm==3){
	  if( nnnmaxPerm != -99 ){
	    maxPerm = nnnmaxPerm;
	    maxM    = nnnmaxM;
	  }
	  else if( nnmaxPerm != -99 ){
	    maxPerm = nnmaxPerm;
	    maxM    = nnmaxM;
	  }
	  else if( nmaxPerm != -99 ){
	    maxPerm = nmaxPerm;
	    maxM    = nmaxM;
	  }
	  else{
	    /*  do not do anything ... */
	  }
	}

	// loop over the mass points and fill hTmp with the mass dependent probability for maxPerm ...
	for(int mH_it = 0; mH_it<nMassPoints; mH_it++){
	  
	  // the mass point value
	  float mH =  m_scan[mH_it];
	  
	  if( mH<0 || mH>300) continue;
	  
	  // again, get the probability, this time for maxPerm
	  float ME_prob =  p_vsM[mH_it*nPermut + maxPerm];
	  float bb_prob =  p_tt_bb[maxPerm] ;
	  float norm    = normXSec ? 1./xsec->Eval( mH ) : 1.0;	
	  
	  // probability for maxPerm under mH hypothesis
	  double p =  ME_prob /** bb_prob */ * norm;       

	  // fill hTmp at mH
	  hTmp->Fill( mH, p );
	}
	
	// get the maximum of hTmp by quadratic interpolation
	pair<double,double> bestMass = getMaxValue(hTmp);

	// the observables
	double mass = bestMass.first;
	double prob = bestMass.second;
	
	// if there is at leats one permutation matched, increment the counters
	if(quarks>0){
	  
	  // at least one permutation matched
	  counterQuarks++;
	  
	  // the best permutation matched the pattern ( OR between two )
	  if( maxM == matchPattern || maxM == matchPattern2 || (maxM%100 == 11 ))  counterMatch++;	
	}
	
	// if the matched permutation assigns correctly at least the Higgs b-quarks...
	if( maxM%100 == 11 )  counterMatchAny++;	
	
	// fill the histogram with the mass     
	hMass_s->Fill(mass, weight*lumiScale);
      }
	
      // ME analysis
      else{
	
	// this number are the total probability under the ttH, ttbb, and ttjj hypotheses
	double p_125_all_s_ttbb = 0.;
	double p_125_all_b_ttbb = 0.;
	double p_125_all_b_ttjj = 0.;
	
	// consider only the mass scan at the nominal Higgs mass
	for(int mH_it = 0; mH_it<nMassPoints; mH_it++){	
	  
	  // skip other mass values...
	  float mH =  m_scan[mH_it];
	  if(mH>126 || mH<124) continue;	  
	  
	  // once we got the good Higgs mass, loop over permutations...
	  for( int perm_it = 0 ; perm_it<nPermut ; perm_it++){	
	    float ME_prob     = p_vsM[mH_it*nPermut + perm_it];
	    float bb_prob     = p_tt_bb[perm_it] ;

	    // if choosePerm<0, choose randomly abs(choosePerm) permutations
	    if( choosePerm<0 ){

	      // check if perm_it is among those selected
	      int isThere = 0;
	      for(int pp = 0 ; pp < whichpermut.size(); pp++){
		if(perm_it == whichpermut[pp]) isThere++;
	      }

	      // if not, skip it
	      if(isThere==0) continue;
	    }

	    // if choosePerm is positive but not 999, choose exacrtly the choosePerm'th permutation
	    else if( choosePerm!=999 ){
	      if(perm_it != choosePerm ) continue;
	    }

	    // else choosePerm==999, so consider perm_it anyway
	    else{}

	    // add to the total probability...
	    p_125_all_s_ttbb += ME_prob/**bb_prob*/;

	  }
	}
	
	// consider only the mass scan at the nominal top mass
	for(int mT_it = 0; mT_it<nMassPoints; mT_it++){
	  
	  // skip other mass values...
	  float mT =  mT_scan[mT_it];
	  if(mT>175 || mT<173) continue;
	  
	  // once we got the good Top mass, loop over permutations...
	  for( int perm_it = 0 ; perm_it<nPermut_b ; perm_it++){	
	    float ME_prob     = p_vsMT_b[mT_it*nPermut_b + perm_it];
	    float bb_prob     = p_tt_bb[perm_it] ;
	    float jj_prob     = p_tt_jj[perm_it] ;

	    // if choosePerm<0, choose randomly abs(choosePerm) permutations
	    if( choosePerm<0 ){

	      // check if perm_it is among those selected
	      int isThere = 0;
	      for(int pp = 0 ; pp < whichpermut.size(); pp++){
		if(perm_it == whichpermut[pp]) isThere++;
	      }

	      // if not, skip it
	      if(isThere==0) continue;
	    } 

	    // if choosePerm is positive but not 999, choose exacrtly the choosePerm'th permutation	    
	    else if( choosePerm!=999 ){
	      if(perm_it != choosePerm ) continue;
	    }

	    // else choosePerm==999, so consider perm_it anyway	    
	    else{}

	    // add to the total probability...	    
	    p_125_all_b_ttbb += ME_prob/**bb_prob*/;
	    p_125_all_b_ttjj += ME_prob/**jj_prob*/;

	  }
	}

	// prepare the sb likelihood ratio ingredients
	double wDen = p_125_all_s_ttbb + p0*( p1*p_125_all_b_ttbb + p2*p_125_all_b_ttjj );
	double wNum = p_125_all_s_ttbb;
	
	// the sb likelihood ratio
	double w = wDen>0 ? wNum/wDen : 0.;
	
	// fill the histogram     
	hMass_s->Fill( w, weight);	

      }
    }

    // underflow bin added to the first bin
    int firstBin              =  1;
    float firstBinContent     =  hMass_s->GetBinContent( firstBin ); 
    float underflowBinContent =  hMass_s->GetBinContent( firstBin-1 ); 
    hMass_s->SetBinContent( firstBin, firstBinContent+underflowBinContent); 

    // overflow bin added to the last bin
    int lastBin              =  hMass_s->GetNbinsX();
    float lastBinContent     =  hMass_s->GetBinContent( lastBin ); 
    float overflowBinContent =  hMass_s->GetBinContent( lastBin+1 ); 
    hMass_s->SetBinContent( lastBin, lastBinContent+overflowBinContent );    

    // a recapitulation...
    cout << "Total: " << counter << endl;
    cout << "All quarks in acceptance: " << counterQuarks << endl;
    cout << "Max permutation is correct: " << counterMatch << endl;
    cout << " ==> eff.      = "      << (counter>0 ? float(counterMatchAny)/float(counter) : 0. ) << endl;
    cout << " ==> acc.      = " << (counter>0 ? float(counterQuarks)/float(counter) : 0.) << endl;
    cout << " ==> eff.| acc = " << (counterQuarks>0 ? float(counterMatch)/float(counterQuarks) : 0. ) << endl;
    cout << "Overflow="  << overflowBinContent << endl;
    cout << "Underflow=" << underflowBinContent << endl;
    
    // save a snapshot of the matching status for ttH
    if( sample.find("TTH")!=string::npos ){
       counterH         = counter;
       counterHQuarks   = counterQuarks;
       counterHMatch    = counterMatch;
       counterHMatchAny = counterMatchAny;
    }
    
    // fill hErr with all the MC but signal
    if( hErr==0 ){
      hErr = (TH1F*)hMass_s->Clone("hErr");
      hErr->Reset();
      leg->AddEntry(hErr, "MC unc. (stat.)", "L");
    }
    if( sample.find("TTH")==string::npos ) hErr->Add( hMass_s, 1.0);

    // deal with tt+jets
    if( sample.find("TTJets")!=string::npos ){
      hMass_s->SetFillColor(color);
      if(color==16){
	leg->AddEntry(hMass_s, "t#bar{t} + bb", "F");
	hTTbb->Add(hMass_s, 1.0);
      }
      if(color==17){
	leg->AddEntry(hMass_s, "t#bar{t} + b", "F");
	hTTbj->Add(hMass_s, 1.0);
      }
      if(color==18){
	leg->AddEntry(hMass_s, "t#bar{t} + jj", "F");
	hTTjj->Add(hMass_s, 1.0);
      }
      hTT->Add(hMass_s, 1.0);
    }

    // deal with ttH
    if( sample.find("TTH")!=string::npos ){
      hMass_s->SetFillColor(kRed);
      hTTH->Add(hMass_s, 1.0);
      leg->AddEntry(hMass_s, "t#bar{t}H", "F");
    }

    // deal with ttV
    if( sample.find("TTV")!=string::npos ){
      hMass_s->SetFillColor(kBlue);
      hTTV->Add(hMass_s, 1.0);
      leg->AddEntry(hMass_s, "t#bar{t}V", "F");
    }

    // deal with single top 
    if( sample.find("SingleT")!=string::npos ){
      hMass_s->SetFillColor(kMagenta);
      leg->AddEntry(hMass_s, "Single top", "F");
    }

    // deal with VV
    if( sample.find("DiBoson")!=string::npos ){
      hMass_s->SetFillColor(kYellow);
      leg->AddEntry(hMass_s, "VV", "F");
    }

    // deal with EWK
    if( sample.find("EWK")!=string::npos ){
      hMass_s->SetFillColor(kGreen);
      leg->AddEntry(hMass_s, "V+jets", "F");
    }

    // add to the stack
    cout << "Adding " << hMass_s->Integral() << " weighted events to the stack" << endl;
    aStack->Add( hMass_s );

    // close the file
    f->Close();

  }


  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  // drawing options 

  hErr->GetYaxis()->SetTitle("Events");
  if(doMHscan) 
    hErr->GetXaxis()->SetTitle("M_{H} estimator [GeV]");
  else
    hErr->GetXaxis()->SetTitle("M_{T} estimator [GeV]");

  hErr->SetTitle("Simulation #sqrt{s}=8 TeV, L=19.04 fb^{-1}");
  hErr->SetTitleSize  (0.04,"X");
  hErr->SetTitleOffset(0.95,"X");
  float max =  hErr->GetMaximum()*1.45;
  hErr->GetYaxis()->SetRangeUser(0., max );
  hErr->SetLineColor(kBlack);
  hErr->Draw("HISTE1");
  aStack->Draw("HISTSAME");
  hTTH->SetLineWidth(3);
  hTTH->SetLineColor(kRed);
  hTTH->SetLineStyle(kDashed);
  hTTH->SetFillColor(kRed);
  hTTH->SetFillStyle(3004);
  hTTH->Scale(5.0);
  hTTH->Draw("HISTSAME");
  hTTV->SetLineWidth(3);
  hTTV->SetLineColor(kBlue);
  hTTV->SetLineStyle(kDashed);
  hTTV->SetFillColor(kBlue);
  hTTV->SetFillStyle(3005);
  hTTV->Scale(5.0);
  hTTV->Draw("HISTSAME");
  hErr->Draw("HISTE1SAME");
  leg->AddEntry(hTTH, "t#bar{t}H x 5", "L");
  leg->AddEntry(hTTV, "t#bar{t}V x 5", "L");

  TPaveText *pt = new TPaveText(0.11,0.76,0.41,0.89,"brNDC");
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetFillColor(10);
  pt->SetTextSize(0.03);
  pt->SetTextAlign(11);
  pt->AddText(Form("All quarks reconstructed.: %.0f%%",    (counterH>0 ? float(counterHQuarks)/float(counterH)*100 : 0.) ))->SetTextColor(kBlack);
  pt->AddText(Form("(%.0f%% w/ Higgs-b's matched)",       (counterHQuarks>0 ? float(counterHMatch)/float(counterHQuarks)*100 : 0. )))->SetTextColor(kRed);
  pt->AddText(Form("Higgs-b's matched : %.0f%%",           (counterH>0 ? float(counterHMatchAny)/float(counterH)*100 : 0. )))->SetTextColor(kBlue);

  if(analysis==0)  pt->Draw();

  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////

  // if 1, then plot the normalized shapes of tt+jest, tt, and ttH
  if(plotShapes==1){

    hTT->SetLineWidth(3);
    hTT->SetLineColor(kBlack);
    hTT->SetLineStyle(kSolid);

    hTTH->Scale(1./hTTH->Integral());
    hTTV->Scale(1./hTTV->Integral());
    hTT ->Scale(1./hTT ->Integral());

    hTTH->GetYaxis()->SetRangeUser(0., TMath::Max(hTTV->GetMaximum(), hTTH->GetMaximum())*1.45);
    //hTTH->GetYaxis()->SetRangeUser(0., doMHscan ? 0.35 : 0.50);
    hTTH->Draw("HISTE");
    hTTV->Draw("HISTESAME");
    //hTT->Draw("HISTESAME");
    leg->Clear();
    leg->SetHeader(category);
    leg->AddEntry(hTTH, "t#bar{t}H",     "F");
    leg->AddEntry(hTTV, "t#bar{t}V",     "F");
    //leg->AddEntry(hTT,  "t#bar{t}+jets", "L");
    leg->Draw();
  }
  // if 2, then plot the normalized shapes of tt+bb, tt+b, and tt+jj
  else if(plotShapes==2){

    hTTbb->SetLineWidth(3);
    hTTbb->SetLineColor(kRed);
    hTTbb->SetLineStyle(kSolid);

    hTTbj->SetLineWidth(3);
    hTTbj->SetLineColor(kBlue);
    hTTbj->SetLineStyle(kSolid);
    
    hTTjj->SetLineWidth(3);
    hTTjj->SetLineColor(kMagenta);
    hTTjj->SetLineStyle(kSolid);

    hTTbb->Scale(1./hTTbb->Integral());
    hTTbj->Scale(1./hTTbj->Integral());
    hTTjj->Scale(1./hTTjj->Integral());

    hTTbb->GetYaxis()->SetRangeUser(0., TMath::Max( hTTbb->GetMaximum(), hTTjj->GetMaximum())*1.45);
    hTTbb->Draw("HISTE");
    hTTbj->Draw("HISTESAME");
    hTTjj->Draw("HISTESAME");
    leg->Clear();
    leg->SetHeader(category);
    leg->AddEntry(hTTbb, "t#bar{t}+bb", "F");
    leg->AddEntry(hTTbj, "t#bar{t}+b",  "F");
    leg->AddEntry(hTTjj, "t#bar{t}+jj", "F");
    leg->Draw();
  }
  else{
    leg->Draw();
  }

  // save the plots
  if(0){

    TString what = doMHscan ? "MHscan" : "MTscan";

    TString whichVar;
    if     (analysis==0) whichVar = what+"_Mass_"; 
    else if(analysis==1) whichVar = what+"_MEM_"; 
    else                 whichVar = what+"_UNKNOW_"; 
    
    TString whichPlot;
    if     (plotShapes==0) whichPlot = "Spectrum_"; 
    else if(plotShapes==1) whichPlot = "Shapes_ttH_"; 
    else if(plotShapes==2) whichPlot = "Shapes_ttjets_"; 
    else                   whichPlot = "UNKNOW_";

    TString whichPerm;
    if     ( choosePerm == 999 ) whichPerm = "AllPerm_";
    else if( choosePerm < 0   )  whichPerm = Form("%dRandom_", abs(choosePerm));
    else                         whichPerm = Form("%dthPerm_",  choosePerm );
    
    c1->SaveAs("Plots/Plot_"+whichVar+whichPlot+whichPerm+fname+"_"+plotname+".png");
  }

  // delete stuff
  delete ran;

  // remove the dummy file and return
  cout << "Remove dummy file" << endl;
  gSystem->Exec("rm dummy.root");
  return;
}



void plot_massAll( int analysis = 0 ){

  float xMin = 100;   //25
  float xMax = 250;  //300
  int nBins  = 9;    //9
  
  if( analysis == 0 ){

    // the mass analysis => shapes and specta for the 1st rank
    for(int i = 0; i < 1; i++){
      //plot_mass("SL", "type==0", "SL + 6 jets", "SL6jets", i, 19.5/12.1 ,   nBins, xMin, xMax, 0, 1, 111111 , 111111,  999, 0);
      //plot_mass("SL", "type==2", "SL + 5 jets", "SL5jets", i, 19.5/12.1 ,   nBins, xMin, xMax, 0, 1, 111111 , 111111,  999, 0);
      //plot_mass("SL", "type==3", "SL + >6 jets","SL7jets", i, 19.5/12.1 ,   nBins, xMin, xMax, 0, 1, 111111 , 111111,  999, 0);
      plot_mass("DL", "type==6", "DL + 4 jets", "DL4jets", i, 2*19.5/12.1 ,   nBins, xMin, xMax, 0, 1,   100111 , 11,      999, 0);
    }

    return;
    
    // the mass analysis => shapes and specta for the 2nd rank
    for(int i = 0; i < 3; i++){
      plot_mass("SL", "type==0", "SL + 6 jets", "SL6jets", i, 19.5/12.1 ,   nBins, xMin, xMax, 0, 1, 111111 , 111111,  1, 0);
      plot_mass("SL", "type==2", "SL + 5 jets", "SL5jets", i, 19.5/12.1 ,   nBins, xMin, xMax, 0, 1, 111111 , 111111,  1, 0);
      plot_mass("SL", "type==3", "SL + >6 jets","SL7jets", i, 19.5/12.1 ,   nBins, xMin, xMax, 0, 1, 111111 , 111111,  1, 0);
      plot_mass("DL", "type==6", "DL + 4 jets", "DL4jets", i, 2*19.5/12.1 , nBins, xMin, xMax, 0, 1, 100111 , 11,      1, 0);
    }
    
    // the mass analysis => shapes and specta for the 3rd rank
    for(int i = 0; i < 3; i++){
      plot_mass("SL", "type==0", "SL + 6 jets", "SL6jets", i, 19.5/12.1 ,   nBins, xMin, xMax, 0, 1, 111111 , 111111,  2, 0);
      plot_mass("SL", "type==2", "SL + 5 jets", "SL5jets", i, 19.5/12.1 ,   nBins, xMin, xMax, 0, 1, 111111 , 111111,  2, 0);
      plot_mass("SL", "type==3", "SL + >6 jets","SL7jets", i, 19.5/12.1 ,   nBins, xMin, xMax, 0, 1, 111111 , 111111,  2, 0);
      plot_mass("DL", "type==6", "DL + 4 jets", "DL4jets", i, 2*19.5/12.1 , nBins, xMin, xMax, 0, 1, 100111 , 11,      2, 0);
    }
    
    // the mass analysis => shapes and specta for the 4th rank
    for(int i = 0; i < 3; i++){
      plot_mass("SL", "type==0", "SL + 6 jets", "SL6jets", i, 19.5/12.1 ,   nBins, xMin, xMax, 0, 1, 111111 , 111111,  3, 0);
      plot_mass("SL", "type==2", "SL + 5 jets", "SL5jets", i, 19.5/12.1 ,   nBins, xMin, xMax, 0, 1, 111111 , 111111,  3, 0);
      plot_mass("SL", "type==3", "SL + >6 jets","SL7jets", i, 19.5/12.1 ,   nBins, xMin, xMax, 0, 1, 111111 , 111111,  3, 0);
      plot_mass("DL", "type==6", "DL + 4 jets", "DL4jets", i, 2*19.5/12.1 , nBins, xMin, xMax, 0, 1, 100111 , 11,      3, 0);
    }

  }

 if( analysis == 1 ){

   /*
   // the mass analysis => shapes and specta for sum over permutations
   for(int i = 0; i < 3; i++){
     plot_mass("SL", "type==0", "Cat. 1", "cat1",                                   i,  19.5/12.1 ,   8 , 0,  1,  0, 1, 111111 , 111111,  999, 1, 1.0, 0.00611416, 45.1513);
     //plot_mass("SL", "type==1", "Cat. 2", "cat2",                                 i,  19.5/12.1 ,   8 , 0,  1,  0, 1, 111111 , 111111,  999, 1, 1.7, 0.00633982, 7.12076);
     plot_mass("SL", "type==2 && flag_type2>0",  "Cat. 3", "cat3",                  i,  19.5/12.1 ,   8 , 0,  1,  0, 1, 111111 , 111111,  999, 1, 2.2, 0.00203578, 6.99398);
     plot_mass("SL", "type==2 && flag_type2<=0", "Cat. 4", "cat4",                  i,  19.5/12.1 ,   8 , 0,  1,  0, 1, 111111 , 111111,  999, 1, 2.0, 0.00346887, 8.96759);
     plot_mass("SL", "type==3 && flag_type3>0 && p_125_all_s>0", "Cat. 5", "cat5",  i,  19.5/12.1 ,   8 , 0,  1,  0, 1, 111111 , 111111,  999, 1, 5.5, 0.00163639, 6.92078);
     plot_mass("DL", "type==6", "Cat. 6", "cat6",                                   i,2*19.5/12.1 ,   8 , 0,  1,  0, 1, 111111 , 111111,  999, 1, 1.5, 0.00398151, 5.95277);
   }
   */

   // the mass analysis => shapes and specta for sum over permutations
   for(int d = -4; d<-3 ; d++){
     for(int i = 0; i < 3; i++){
       plot_mass("SL", "type==0", "Cat. 1", "cat1",                                   i,  19.5/12.1 ,   8 , 0,  1,  0, 1, 111111 , 111111,  d, 1, 1.0, 0.00611416, 45.1513);
       //plot_mass("SL", "type==1", "Cat. 2", "cat2",                                 i,  19.5/12.1 ,   8 , 0,  1,  0, 1, 111111 , 111111,  d, 1, 1.7, 0.00633982, 7.12076);
       plot_mass("SL", "type==2 && flag_type2>0",  "Cat. 3", "cat3",                  i,  19.5/12.1 ,   8 , 0,  1,  0, 1, 111111 , 111111,  d, 1, 2.2, 0.00203578, 6.99398);
       plot_mass("SL", "type==2 && flag_type2<=0", "Cat. 4", "cat4",                  i,  19.5/12.1 ,   8 , 0,  1,  0, 1, 111111 , 111111,  d, 1, 2.0, 0.00346887, 8.96759);
       plot_mass("SL", "type==3 && flag_type3>0 && p_125_all_s>0", "Cat. 5", "cat5",  i,  19.5/12.1 ,   8 , 0,  1,  0, 1, 111111 , 111111,  d, 1, 5.5, 0.00163639, 6.92078);
       plot_mass("DL", "type==6", "Cat. 6", "cat6",                                   i,2*19.5/12.1 ,   8 , 0,  1,  0, 1, 111111 , 111111,  d, 1, 1.5, 0.00398151, 5.95277);
     }
   }


 }



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
	       TString fname  = ""
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
  h_1->SetLineStyle(kDashed);
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
  leg->AddEntry( h_2, header_2 , "L");

  h_1->SetXTitle("discriminant");
  h_1->SetYTitle("units");

  h_1->SetMaximum( TMath::Max( float(h_1->GetMaximum()), float(h_2->GetMaximum()) )*1.3 );
  h_1->SetMinimum(0);
  h_1->DrawNormalized("HISTE");
  h_2->DrawNormalized("HISTESAME");

  cout << "*********" << endl;
  cout << "h_1 := " << h_1->Integral() << ", mu = " << h_1->GetMean() << endl;
  cout << "h_2 := " << h_2->Integral() << ", mu = " << h_2->GetMean() << endl;
  cout << "*********" << endl;

  leg->SetHeader( header );
  leg->Draw();

  if(1){
    c1->SaveAs("Plots/Feb24_2014/Plot_Comp_"+fname+".png");
    delete c1;
  }

}


void plot_compAll(){


  for(int cat = 1; cat<7 ; cat++){

    cout << "Doing............." << cat << endl;

    plot_comp("datacards/Feb24_2014/MEM_New_rec_std_byLLR_sb.root" , "datacards/Feb24_2014/MEM_New_rec_reg_byLLR_sb.root" , 
	      Form("MEM_cat%d",cat),    Form("MEM_cat%d",cat), 
	      "TTH125", "TTH125", 
	      "nominal", "nominal", 
	      "standard", "regression",
	      Form("ttH, Cat. %d", cat),
	      Form("TTH125_cat%d_std_vs_reg", cat));

    plot_comp("datacards/Feb24_2014/MEM_New_rec_std_byLLR_sb.root" , "datacards/Feb24_2014/MEM_New_rec_reg_byLLR_sb.root" , 
	      Form("MEM_cat%d",cat),    Form("MEM_cat%d",cat), 
	      "TTJetsHFbb", "TTJetsHFbb", 
	      "nominal", "nominal", 
	      "standard", "regression",
	      Form("tt+bb, Cat. %d", cat),
	      Form("TTJetsHFbb_cat%d_std_vs_reg", cat));

    plot_comp("datacards/Feb24_2014/MEM_New_rec_std_byLLR_sb.root" , "datacards/Feb24_2014/MEM_New_rec_reg_byLLR_sb.root" , 
	      Form("MEM_cat%d",cat),    Form("MEM_cat%d",cat), 
	      "TTJetsLF", "TTJetsLF", 
	      "nominal", "nominal", 
	      "standard", "regression",
	      Form("tt+jj, Cat. %d", cat),
	      Form("TTJetsLF_cat%d_std_vs_reg", cat));

    plot_comp("datacards/Feb24_2014/MEM_New_rec_reg_byLLR_sb.root" , "datacards/Feb24_2014/MEM_New_rec_reg_byLLR_sb.root" , 
	      Form("MEM_cat%d",cat),    Form("MEM_cat%d",cat), 
	      "TTJetsHFbb", "TTH125", 
	      "nominal", "nominal", 
	      "t#bar{t}+b#bar{b}", "ttH",
	      Form("Cat. %d", cat),
	      Form("cat%d_TTH125_vs_TTJetsHFbb", cat));

    plot_comp("datacards/Feb24_2014/MEM_New_rec_reg_byLLR_sb_nb.root" , "datacards/Feb24_2014/MEM_New_rec_reg_byLLR_sb_nb.root" , 
	      Form("MEM_cat%d",cat),    Form("MEM_cat%d",cat), 
	      "TTJetsHFbb", "TTH125", 
	      "nominal", "nominal", 
	      "t#bar{t}+b#bar{b}", "ttH",
	      Form("Cat. %d", cat),
	      Form("cat%d_TTH125_vs_TTJetsHFbb_nb", cat));

    plot_comp("datacards/Feb24_2014/MEM_New_rec_reg_byLLR_bj.root" , "datacards/Feb24_2014/MEM_New_rec_reg_byLLR_bj.root" , 
	      Form("MEM_cat%d",cat),    Form("MEM_cat%d",cat), 
	      "TTJetsHFbb", "TTH125", 
	      "nominal", "nominal", 
	      "t#bar{t}+b#bar{b}", "ttH",
	      Form("Cat. %d", cat),
	      Form("cat%d_TTH125_vs_TTJetsHFbb_bj", cat));

    plot_comp("datacards/Feb24_2014/MEM_New_rec_reg_byLLR_sb.root" , "datacards/Feb24_2014/MEM_New_rec_reg_byLLR_sb.root" , 
	      Form("MEM_cat%d",cat),    Form("MEM_cat%d",cat), 
	      "TTJetsHFbb", "TTJetsLF", 
	      "nominal", "nominal", 
	      "t#bar{t}+b#bar{b}", "t#bar{t}+jj",
	      Form("Cat. %d", cat),
	      Form("cat%d_TTJetsHFbb_vs_TTJetsLF", cat));

    plot_comp("datacards/Feb24_2014/MEM_New_rec_reg_byLLR_bj.root" , "datacards/Feb24_2014/MEM_New_rec_reg_byLLR_bj.root" , 
	      Form("MEM_cat%d",cat),    Form("MEM_cat%d",cat), 
	      "TTJetsHFbb", "TTJetsLF", 
	      "nominal", "nominal", 
	      "t#bar{t}+b#bar{b}", "t#bar{t}+jj",
	      Form("Cat. %d", cat),
	      Form("cat%d_TTJetsHFbb_vs_TTJetsLF_bj", cat));
    
    plot_comp("datacards/Feb24_2014/MEM_New_rec_reg_byLLR_sb.root" , "datacards/Feb24_2014/MEM_New_rec_reg_byLLR_sb.root" , 
	      Form("MEM_cat%d",cat),    Form("MEM_cat%d",cat), 
	      "TTV", "TTH125", 
	      "nominal", "nominal", 
	      "t#bar{t}V", "ttH",
	      Form("Cat. %d", cat),
	      Form("cat%d_TTH125_vs_TTV", cat));

    plot_comp("datacards/Feb24_2014/MEM_New_rec_reg_byLLR_sb.root" , "datacards/Feb24_2014/MEM_New_rec_reg_byLLR_sb.root" , 
	      Form("MEM_cat%d",cat),    Form("MEM_cat%d",cat), 
	      "SingleT", "TTH125", 
	      "nominal", "nominal", 
	      "Single-top", "ttH",
	      Form("Cat. %d", cat),
	      Form("cat%d_TTH125_vs_SingleT", cat));
    
    /*
    plot_comp(Form("SL_VType%d_New_std.root", cat<6 ? 2 : 0), Form("SL_VType%d_New_reg.root",  cat<6 ? 2 : 0), 
	      Form("SL_VType%d_cat%d",  cat<6 ? 2 : 0, cat),  Form("SL_VType%d_cat%d",  cat<6 ? 2 : 0, cat), 
	      "TTH125", "TTH125", 
	      "nominal", "nominal", 
	      "Nov21 standard", "Nov21 regression",
	      Form("ttH, Cat. %d", cat),
	      Form("TTH_cat%d_std_vs_reg", cat));

    plot_comp(Form("~/public/Nov04_2013/mem/%sL_New.root", cat<6 ? "S" : "D" ) , Form("SL_VType%d_New_std.root",  cat<6 ? 2 : 0), 
	      Form("%sL_cat%d", cat<6 ? "S" : "D", cat),  Form("SL_VType%d_cat%d",  cat<6 ? 2 : 0, cat), 
	      "TTH125", "TTH125", 
	      "nominal", "nominal", 
	      "Nov04", "Nov21",
	      Form("ttH, Cat. %d", cat),
	      Form("TTH_cat%d_old_vs_new", cat));


    plot_comp(Form("SL_VType%d_New_std.root", cat<6 ? 2 : 0), Form("SL_VType%d_New_reg.root",  cat<6 ? 2 : 0), 
	      Form("SL_VType%d_cat%d",  cat<6 ? 2 : 0, cat),  Form("SL_VType%d_cat%d",  cat<6 ? 2 : 0, cat), 
	      "TTJetsHFbb", "TTJetsHFbb", 
	      "nominal", "nominal", 
	      "Nov21 standard", "Nov21 regression",
	      Form("ttbb, Cat. %d", cat),
	      Form("TTbb_cat%d_std_vs_reg", cat));

    plot_comp(Form("~/public/Nov04_2013/mem/%sL_New.root", cat<6 ? "S" : "D" ) , Form("SL_VType%d_New_std.root",  cat<6 ? 2 : 0), 
	      Form("%sL_cat%d", cat<6 ? "S" : "D", cat),  Form("SL_VType%d_cat%d",  cat<6 ? 2 : 0, cat), 
	      "TTJetsHFbb", "TTJetsHFbb", 
	      "nominal", "nominal", 
	      "Nov04", "Nov21",
	      Form("ttbb, Cat. %d", cat),
	      Form("TTbb_cat%d_old_vs_new", cat)); 

    plot_comp(Form("SL_VType%d_New_std.root", cat<6 ? 2 : 0), Form("SL_VType%d_New_reg.root",  cat<6 ? 2 : 0), 
	      Form("SL_VType%d_cat%d",  cat<6 ? 2 : 0, cat),  Form("SL_VType%d_cat%d",  cat<6 ? 2 : 0, cat), 
	      "TTJetsLF", "TTJetsLF", 
	      "nominal", "nominal", 
	      "Nov21 standard", "Nov21 regression",
	      Form("ttjj, Cat. %d", cat),
	      Form("TTjj_cat%d_std_vs_reg", cat));

    plot_comp(Form("~/public/Nov04_2013/mem/%sL_New.root", cat<6 ? "S" : "D" ) , Form("SL_VType%d_New_std.root",  cat<6 ? 2 : 0), 
	      Form("%sL_cat%d", cat<6 ? "S" : "D", cat),  Form("SL_VType%d_cat%d",  cat<6 ? 2 : 0, cat), 
	      "TTJetsLF", "TTJetsLF", 
	      "nominal", "nominal", 
	      "Nov04", "Nov21",
	      Form("ttjj, Cat. %d", cat),
	      Form("TTjj_cat%d_old_vs_new", cat));    
    */

  }


}
