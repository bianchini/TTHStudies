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

#include "TRandom3.h"

typedef TMatrixT<double> TMatrixD;

#define RUNONDATA 0


string version   =  "_v4_gen_std" ;
string inputpath = "gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/Trees/MEM/Dec06_2013/WS/";
string tag = "New_v4_gen_std";

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
		   string header = "Cat 1",
		   string fname  = "",
		   string extraname = "_MEM",
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

  THStack* aStack = new THStack("aStack","Simulation #sqrt{s}=8 TeV, L=19.5 fb^{-1};  L_{S}/(L_{S}+L_{B}) ; events ");

  vector<string> samples;
  if( type.find("SL")!=string::npos ){
    samples.push_back("SingleT");
    samples.push_back("TTV");
    samples.push_back("TTJetsHFbb");
    samples.push_back("TTJetsHFb");
    samples.push_back("TTJetsLF");
    samples.push_back("TTH125");
    if(RUNONDATA)  samples.push_back("data_obs");
  }
  if( type.find("DL")!=string::npos ){
    samples.push_back("SingleT");
    samples.push_back("TTV");
    samples.push_back("TTJetsHFbb");
    samples.push_back("TTJetsHFb");
    samples.push_back("TTJetsLF");
    samples.push_back("TTH125");
    if(RUNONDATA)  samples.push_back("data_obs");
  }

  TH1F* hS    = 0;
  TH1F* hData = 0;
  TH1F* hErr  = 0;

  TFile* f = TFile::Open(("datacards/WS/"+type+"_"+tag+extraname+".root").c_str());
  if(f==0 || f->IsZombie() ) return;

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
      h->Scale( rescaleTTbb );
      h->SetLineColor( 16 );
      h->SetFillColor( 16 );
    }
    else if( samples[sample].find("TTJetsHFb") != string::npos ){
      leg->AddEntry(h, "t#bar{t} + b", "F");
      h->Scale( rescaleTTb );
      h->SetLineColor( 17 );
      h->SetFillColor( 17 );
    }
    if( samples[sample].find("TTJetsLF") != string::npos ){
      leg->AddEntry(h, "t#bar{t} + jj", "F");
      h->Scale( rescaleTTjj );
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
      hData->SetMarkerSize(1.8);
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
	leg->AddEntry(hErr, "MC unc. (stat.)", "L");
      }
      hErr->Add( h, 1.0);
    }

    if( samples[sample].find("data_obs") == string::npos )
      aStack->Add( h );
  }

  if(hErr==0 || hS==0) return;

  hErr->GetYaxis()->SetTitle("Events");
  if(fname.find("MEM")!=string::npos) 
    hErr->GetXaxis()->SetTitle("P_{s/b}(y)");
  if(fname.find("bbjj")!=string::npos) 
    hErr->GetXaxis()->SetTitle("P_{b/j}(y)");
  hErr->SetTitle("Simulation #sqrt{s}=8 TeV, L=19.5 fb^{-1}");
  hErr->SetTitleSize  (0.04,"X");
  hErr->SetTitleOffset(0.95,"X");
  float max =  TMath::Max(float(hErr->GetMaximum()*1.35),(hData!=0 ? hData->GetMaximum()*1.35 : -1.));
  hErr->GetYaxis()->SetRangeUser(0., max );
  hErr->SetLineColor(kBlack);
  hErr->Draw("HISTE1");
  aStack->Draw("HISTSAME");
  leg->SetHeader( header.c_str() );
  leg->AddEntry(hS, "signal x 5", "L");
  hS->Draw("HISTSAME");
  hErr->Draw("HISTE1SAME");
  if(hData!=0) hData->Draw("ESAME");
  leg->Draw();

  TLine* line = new TLine(hErr->GetBinLowEdge(3), max , hErr->GetBinLowEdge(3), 0.);
  line->SetLineWidth(4);
  line->SetLineStyle(kSolid);
  line->SetLineColor(kBlack);
  //if(type.find("SL")!=string::npos) 
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
 
  TPaveText *pt3 = new TPaveText(0.52,0.42,0.81,0.48,"brNDC");
  pt3->SetFillStyle(1001);
  pt3->SetBorderSize(0);
  pt3->SetFillColor(kWhite);
  pt3->SetTextSize(0.042); 
  pt3->SetTextColor(kRed); 
  if(fname.find("MEM_nb")!=string::npos) 
    pt3->AddText("ttH vs. ttbb, one permut.");
  else if(fname.find("MEM")!=string::npos) 
    pt3->AddText("ttH vs. ttbb, all permut.");
  else if(fname.find("bbjj_nb")!=string::npos) 
    pt3->AddText("ttbb vs. ttjj, one permut.");
  else if(fname.find("bbjj")!=string::npos) 
    pt3->AddText("ttbb vs. ttjj, all permut.");


  if(drawseparation) pt1->Draw();
  if(drawseparation) pt2->Draw();
  pt3->Draw();

  cout << "Signal = " << hS->Integral()/5. << endl;

  if(1){
    c1->SaveAs(  ("Plots/WS/Plot_"+cat+"_"+fname+"_"+tag+".png").c_str() );
  }

}


void plotAll(string extraname = "_MEM"){


  
  plot_category( "SL", 
		 "cat1",
		 "SL, 6j (W-tag)",
		 extraname,
		 extraname,
		 0
		 );
  
  
  plot_category( "SL", 
		 "cat2",
		 "SL, 6j (!W-tag)",
		 extraname,
		 extraname,
		 0
		 ); 
  
  plot_category( "SL", 
		 "cat4",
		 "SL, 5j",
		 extraname,
		 extraname,
		 0
		 );
  
  plot_category( "DL", 
		 "cat6",
		 "DL, #geq4j",
		 extraname,
		 extraname,
		 0
		 );  

}




void plot_limit(string version = "_New_v4_gen_std_MEM", int scan = 0, string cat = "cat1", string title = "SL, 6j (W-tag)"){

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


  vector<string> categories;        vector<string> names;
  categories.push_back("SL_cat1"+version);  names.push_back("SL, 6j (W-tag)");
  categories.push_back("SL_cat2"+version);  names.push_back("SL, 6j (!W-tag)");
  categories.push_back("SL_cat4"+version);  names.push_back("SL, 5j");
  categories.push_back("SL"+version);       names.push_back("SL comb.");
  categories.push_back("DL_cat6"+version);  names.push_back("DL #geq4j"); 
  categories.push_back("COMB"+version);     names.push_back("All comb.");

  if(scan){
    string analysis = cat.find("cat6")!=string::npos ? "DL" : "SL";
    categories.clear(); names.clear();
    categories.push_back(analysis+"_"+cat+""+version+"_0");  names.push_back("0.8%");
    categories.push_back(analysis+"_"+cat+""+version+"_1");  names.push_back("1.2%");
    categories.push_back(analysis+"_"+cat+""+version+"_2");  names.push_back("1.6%");
    categories.push_back(analysis+"_"+cat+""+version+"_3");  names.push_back("2.0%");
    categories.push_back(analysis+"_"+cat+""+version+"_4");  names.push_back("2.4%");
    categories.push_back(analysis+"_"+cat+""+version+"_5");  names.push_back("2.8%");
    categories.push_back(analysis+"_"+cat+""+version+"_6");  names.push_back("3.2%");
  }



  int nBins = categories.size();

  for( int b = 0; b < nBins; b++){

    TFile* f = TFile::Open(("datacards/WS/higgsCombine"+categories[b]+".Asymptotic.mH120.root").c_str());
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

  TH1F* hT = new TH1F("hT", "", nBins, 0, nBins);
  for(int k = 1; k <= hT->GetNbinsX(); k++){
    int approx = int(expY[k-1]*10);
    hT->SetBinContent(k, approx/10.   );
  }
  hT->SetMarkerSize(2.);
  hT->Draw("TEXT0SAME");

  TF1 *line = new TF1("line","1",0,nBins);
  line->SetLineColor(kRed);
  line->SetLineWidth(3);

  line->Draw("SAME");

  TF1 *lineML = new TF1("lineML","2.4",0,nBins);
  lineML->SetLineColor(kBlue);
  lineML->SetLineStyle(kDashed);
  lineML->SetLineWidth(3);

  //lineML->Draw("SAME");

  TF1 *lineTTH = new TF1("lineTTH","4.1",0,nBins);
  lineTTH->SetLineColor(kMagenta);
  lineTTH->SetLineStyle(kDashed);
  lineTTH->SetLineWidth(3);

  //lineTTH->Draw("SAME");

  c1->cd();
  gPad->Modified();
  mg->GetXaxis()->Set(nBins,0,nBins);

  //mg->GetXaxis()->SetRange(-1,11);

  for( int b = 0; b < nBins; b++){
    mg->GetXaxis()->SetBinLabel(b+1, names[b].c_str() );
  }
  mg->GetXaxis()->SetLabelSize(0.05);
  mg->GetYaxis()->SetLabelSize(0.04);
  mg->GetYaxis()->SetTitleSize(0.05);


  float maxY = 50;

  mg->GetYaxis()->SetTitleOffset(0.82);
  mg->SetMinimum(0.8);
  mg->SetMaximum( maxY );
  mg->GetXaxis()->SetTitle("");
  mg->GetYaxis()->SetTitle("95% CL upper limit on #mu = #sigma/#sigma_{SM}");

  //leg->AddEntry( lineTTH, "HIG-13-019 (post-fit)", "L");
  //leg->AddEntry( lineML,  "HIG-13-020 (post-fit)", "L");
  leg->Draw();

  TPaveText *pt = new TPaveText(0.106811,0.155594,0.407121,0.286713,"brNDC");
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetFillColor(10);
  pt->SetTextSize(0.04);
  pt->SetTextAlign(11);
  if(!scan){
    pt->AddText(Form("Comb: #mu < %.1f at 95%% CL", expY[categories.size()-1]))->SetTextColor(kRed);
    pt->AddText(Form("(SL: #mu < %.1f)", expY[categories.size()-3]))->SetTextColor(kBlack);
    pt->AddText(Form("(DL: #mu < %.1f)", expY[categories.size()-2]))->SetTextColor(kBlack);
  }
  else{
    pt->AddText( title.c_str() )->SetTextColor(kRed);
    pt->AddText( "Upper limit as a function of  #sigma_{ttH}/#sigma_{ttbb}" );
  }

  pt->Draw();

  //c1->Modified();
  //c1->Draw();

  TLine* del = new TLine(3., maxY, 3., 0. );
  del->SetLineWidth(3);
  del->SetLineStyle(kSolid);
  del->SetLineColor(kBlack);
  if(!scan) del->Draw("SAME");

  TLine* del_ = new TLine(4., maxY, 4., 0. );
  del_->SetLineWidth(2);
  del_->SetLineStyle(kDotted);
  del_->SetLineColor(kBlack);
  //del_->Draw("SAME");

  TLine* del2 = new TLine(5., maxY, 5., 0. );
  del2->SetLineWidth(3);
  del2->SetLineStyle(kSolid);
  del2->SetLineColor(kBlack);
  if(!scan) del2->Draw("SAME");

  TLine* del2_ = new TLine(7., maxY, 7., 0. );
  del2_->SetLineWidth(2);
  del2_->SetLineStyle(kDotted);
  del2_->SetLineColor(kBlack);
  //del2_->Draw("SAME");

  if(1){
    c1->SaveAs("Plots/WS/Limits"+TString((version+cat).c_str())+".png");
    //c1->SaveAs("datacards/Limits"+TString(version.c_str())+".pdf");
  }

}



void plot_limitAll(){

  plot_limit("_New_v4_gen_std_MEM", 0, "");
  plot_limit("_New_v4_gen_std_MEM", 1, "cat1", "SL, 6j (W-tag)");
  plot_limit("_New_v4_gen_std_MEM", 1, "cat2", "SL, 6j (!W-tag)");
  plot_limit("_New_v4_gen_std_MEM", 1, "cat4", "SL, 5j");
  plot_limit("_New_v4_gen_std_MEM", 1, "cat6", "DL, #geq4j");

}


void plot_match( TString fname = "SL", 
		string cut = "type==0", 
		TString category = "cat0",
		string header = "",
		int doSgn = 1){
  
  cout << "Doing version " << version << " and category " << category << endl;

  string basecut = cut;
  TString ttjets =  "_TTJets" ;


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

  TFile* f_B = TFile::Open(inputpath+"MEAnalysisNew_"+fname+"_nominal"+version+ttjets+".root");
  TTree* tB = (TTree*)f_B->Get("tree");
  TFile* f_S = TFile::Open(inputpath+"MEAnalysisNew_"+fname+"_nominal"+version+"_TTH125.root");
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
    c1->SaveAs("Plots/WS/Plot_Match_"+category+".png");
  }

}

void plot_matchAll(){

  plot_match("SL", "type==0", "cat1_s", "Cat 1, t#bar{t}H",     1);
  plot_match("SL", "type==0", "cat1_b", "Cat 1, t#bar{t}+jets", 0);

  plot_match("SL", "type==1", "cat2_s", "Cat 2, t#bar{t}H",     1);
  plot_match("SL", "type==1", "cat2_b", "Cat 2, t#bar{t}+jets", 0);

  plot_match("SL", "type==2", "cat3_s", "Cat 3, t#bar{t}H",     1);
  plot_match("SL", "type==2", "cat3_b", "Cat 3, t#bar{t}+jets", 0);

  plot_match("DL", "type==6", "cat6_s", "Cat 6, t#bar{t}H",     1);
  plot_match("DL", "type==6", "cat6_b", "Cat 6, t#bar{t}+jets", 0);

}



void plot_mass(// secondary name of the input trees 
	       TString fname = "SL", 

	       // MT or MH ?
	       int doMHscan = 1,

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

  // inclusive counter (for ttjj)
  int counterT         = 0;

  // counter incr. if all quarks in acceptance (for ttjj)
  int counterTQuarks   = 0;

  // counter incr. if all quarks in acceptanc & best comb. is matched (for ttjj)
  int counterTMatch    = 0;

  // counter incr. if best comb. has the higgs's b's matched (for ttjj)
  int counterTMatchAny = 0;

  
  // cross-section norm. (a simple power law)
  TF1* xsec = new TF1("xsec",Form("x^(-%f)", fact1 ), 20, 500);

  // histo labels
  TString histolabel = "Simulation #sqrt{s}=8 TeV, L=19.5 fb^{-1}; M_{H} estimator [GeV] ; events";

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
  if( !doMHscan ) name = "New_MTscan";
  
  // version
  cout << "Doing version " << version << " and category " << category << endl;
 
 
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
	  if( (doMHscan && maxM%100 == 11 ) || (!doMHscan && maxM%100==11) )  counterMatch++;	
	}
	
	// if the matched permutation assigns correctly at least the Higgs b-quarks...
	if( (doMHscan && maxM%100 == 11) || (!doMHscan && maxM%100==11) )  counterMatchAny++;	
	
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

    // save a snapshot of the matching status for ttjj
    if( sample.find("TTJets")!=string::npos && color==16){
       counterT         = counter;
       counterTQuarks   = counterQuarks;
       counterTMatch    = counterMatch;
       counterTMatchAny = counterMatchAny;
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

  hErr->SetTitle("Simulation #sqrt{s}=8 TeV, L=19.5 fb^{-1}");
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

  float recoEff    = doMHscan ? (counterH>0 ? float(counterHQuarks)/float(counterH)*100 : 0.) :             (counterT>0 ? float(counterTQuarks)/float(counterT)*100 : 0.);
  float recoEffAcc = doMHscan ? (counterHQuarks>0 ? float(counterHMatch)/float(counterHQuarks)*100 : 0. ) : (counterTQuarks>0 ? float(counterTMatch)/float(counterTQuarks)*100 : 0. );
  float matchEff   = doMHscan ? (counterH>0 ? float(counterHMatchAny)/float(counterH)*100 : 0. ) :          (counterT>0 ? float(counterTMatchAny)/float(counterT)*100 : 0. );

  // always refer to ttH:
  recoEff    = (counterH>0 ? float(counterHQuarks)/float(counterH)*100 : 0.);
  recoEffAcc = (counterHQuarks>0 ? float(counterHMatch)/float(counterHQuarks)*100 : 0. );
  matchEff   = (counterH>0 ? float(counterHMatchAny)/float(counterH)*100 : 0. );

  pt->AddText(Form("All quarks reconstructed: %.0f%%",     recoEff ))->SetTextColor(kBlack);
  //if(doMHscan) 
  pt->AddText(Form("(%.0f%% w/ Higgs-b's matched)",       recoEffAcc))->SetTextColor(kRed);
  //else
  //pt->AddText(Form("(%.0f%% w/ top quarks matched)",      recoEffAcc))->SetTextColor(kRed);
  pt->AddText(Form("Higgs-b's matched : %.0f%%",           matchEff))->SetTextColor(kBlue);

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
  if(1){

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
    
    c1->SaveAs("Plots/WS/Plot_"+whichVar+whichPlot+whichPerm+fname+"_"+plotname+".png");
  }

  // delete stuff
  delete ran;

  // remove the dummy file and return
  cout << "Remove dummy file" << endl;
  gSystem->Exec("rm dummy.root");
  return;
}



void plot_massAll( int doMHscan = 1 ){

  float xMin = doMHscan ?  25 :  55; 
  float xMax = doMHscan ? 300 : 295;
  int nBins  = doMHscan ?  10 :  17;
    
  for(int i = 0; i < 2; i++){
    plot_mass("SL", doMHscan, "type==0", "SL, 6j",     "SL6jets", i,   19.5/12.1 , nBins,                             xMin, xMax, 0, 1, doMHscan ? 111111 : 111111,  doMHscan ? 111111 : 11111,   999, 0);
    plot_mass("SL", doMHscan, "type==2", "SL, 5j",     "SL5jets", i,   19.5/12.1 , int(nBins*1.4),                    xMin, xMax, 0, 1, doMHscan ? 111111 : 111111 , doMHscan ? 111111 : 111111,  999, 0);
    plot_mass("DL", doMHscan, "type==6", "DL, #geq4j", "DL4jets", i, 2*19.5/12.1 , int(nBins/(doMHscan ? 1.2 : 1.3)), xMin, xMax, 0, 1, doMHscan ? 100111 : 100111,  doMHscan ?     11 : 11,      999, 0);
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
  leg->SetTextSize(0.045); 


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

  h_1->SetTitle("Simulation #sqrt{s}=8 TeV");
  h_1->GetYaxis()->SetTitle("normalized to 1");
  h_1->GetYaxis()->SetTitleSize(0.05);
  h_1->GetYaxis()->SetTitleOffset(1.1);
  h_1->GetXaxis()->SetTitle("discriminant");
  h_1->GetXaxis()->SetTitleSize(0.05);
  h_1->GetXaxis()->SetTitleOffset(0.9);


  h_1->Scale(1./h_1->Integral());
  h_2->Scale(1./h_2->Integral());

  h_1->SetMaximum( TMath::Max( float(h_1->GetMaximum()/h_1->Integral()), float(h_2->GetMaximum()/h_2->Integral()) )*1.3 );
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
    c1->SaveAs("Plots/WS/Plot_Comp_"+fname+".png");
  }

}


void plot_compAll(int i = 0){

  vector<int> cats;
  cats.push_back(1);
  cats.push_back(2);
  cats.push_back(4);
  cats.push_back(6);

  for(int j = 0; j < cats.size(); j++){

    string analysis = "SL";
    string cat;

    if( cats[j]==1 )
      cat = "SL, 6j (W-tag)";

    if( cats[j]==2 )
      cat = "SL, 6j (!W-tag)";
      
    if( cats[j]==4 )
      cat = "SL, 5j";

    if( cats[j]==6 ) {
      cat = "DL, #geq4j";
      analysis = "DL";
    }

    plot_comp( "datacards/WS/"+analysis+"_New_v4_gen_std_MEM.root", "datacards/WS/"+analysis+"_New_v4_gen_std_MEM.root", 
	       string(Form("%s_cat%d",analysis.c_str(), cats[j])), string(Form("%s_cat%d",analysis.c_str(), cats[j])), 
	       "TTH125",    "TTJetsHFbb", 
	       "nominal", "nominal", 
	       "t#bar{t}H", "t#bar{t}b#bar{b}", 
	       cat+", S/B disc.", 
	       string(Form("cat%d_ttH_vs_ttbb",cats[j])) );
    
    plot_comp( "datacards/WS/"+analysis+"_New_v4_gen_std_MEM.root", "datacards/WS/"+analysis+"_New_v4_gen_std_MEM.root", 
	       string(Form("%s_cat%d",analysis.c_str(), cats[j])), string(Form("%s_cat%d",analysis.c_str(), cats[j])), 
	       "TTH125",    "TTV",        
	       "nominal", "nominal", 
	       "t#bar{t}H", "t#bar{t}V", 
	       cat+", S/B disc.", 	       
	       string(Form("cat%d_ttH_vs_ttV",cats[j]))  );
    
    plot_comp( "datacards/WS/"+analysis+"_New_v4_gen_std_MEM.root", "datacards/WS/"+analysis+"_New_v4_gen_std_MEM.root", 
	       string(Form("%s_cat%d",analysis.c_str(), cats[j])), string(Form("%s_cat%d",analysis.c_str(), cats[j])), 
	       "TTJetsHFb", "TTJetsHFbb", 
	       "nominal", "nominal", 
	       "t#bar{t}b", "t#bar{t}b#bar{b}", 
	       cat+", S/B disc.", 
	       string(Form("cat%d_ttb_vs_ttbb",cats[j]))  );
    
    plot_comp( "datacards/WS/"+analysis+"_New_v4_gen_std_MEM.root", "datacards/WS/"+analysis+"_New_v4_gen_std_MEM.root", 
	       string(Form("%s_cat%d",analysis.c_str(), cats[j])), string(Form("%s_cat%d",analysis.c_str(), cats[j])), 
	       "TTJetsLF",  "TTJetsHFbb", 
	       "nominal", "nominal", 
	       "t#bar{t}jj", "t#bar{t}b#bar{b}", 
	       cat+",S/B disc.", 
	       string(Form("cat%d_ttjj_vs_ttbb",cats[j]))  );

    plot_comp( "datacards/WS/"+analysis+"_New_v4_gen_std_MEM_nb.root", "datacards/WS/"+analysis+"_New_v4_gen_std_MEM_nb.root", 
	       string(Form("%s_cat%d",analysis.c_str(), cats[j])), string(Form("%s_cat%d",analysis.c_str(), cats[j])), 
	       "TTH125",    "TTJetsHFbb", 
	       "nominal", "nominal", 
	       "t#bar{t}H", "t#bar{t}b#bar{b}", 
	       cat+", S/B disc. (ran)", 
	       string(Form("cat%d_nb_ttH_vs_ttbb",cats[j])) );

    plot_comp( "datacards/WS/"+analysis+"_New_v4_gen_std_bbjj.root", "datacards/WS/"+analysis+"_New_v4_gen_std_bbjj.root", 
	       string(Form("%s_cat%d",analysis.c_str(), cats[j])), string(Form("%s_cat%d",analysis.c_str(), cats[j])), 
	       "TTH125",    "TTJetsHFbb", 
	       "nominal", "nominal", 
	       "t#bar{t}H", "t#bar{t}b#bar{b}", 
	       cat+", bb/jj disc.", 
	       string(Form("cat%d_bbjj_ttH_vs_ttbb",cats[j])) );

    plot_comp( "datacards/WS/"+analysis+"_New_v4_gen_std_bbjj.root", "datacards/WS/"+analysis+"_New_v4_gen_std_bbjj.root", 
	       string(Form("%s_cat%d",analysis.c_str(), cats[j])), string(Form("%s_cat%d",analysis.c_str(), cats[j])), 
	       "TTJetsHFbb",    "TTJetsLF", 
	       "nominal", "nominal", 
	       "t#bar{t}b#bar{b}", "t#bar{t}jj", 
	       cat+", bb/jj disc.", 
	       string(Form("cat%d_bbjj_ttbb_vs_ttjj",cats[j])) );

    plot_comp( "datacards/WS/"+analysis+"_New_v4_gen_std_bbjj.root", "datacards/WS/"+analysis+"_New_v4_gen_std_bbjj.root", 
	       string(Form("%s_cat%d",analysis.c_str(), cats[j])), string(Form("%s_cat%d",analysis.c_str(), cats[j])), 
	       "TTJetsHFbb",    "TTJetsHFb", 
	       "nominal", "nominal", 
	       "t#bar{t}b#bar{b}", "t#bar{t}b", 
	       cat+", bb/jj disc.", 
	       string(Form("cat%d_bbjj_ttbb_vs_ttb",cats[j])) );

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
    c1->SaveAs(("Plots/WS/Plot_"+fname+"_CompMadWeight.png").c_str());

}

void plot_analysis_MWAll(){

  plot_analysis_MW("SL", "(type<=3)", "", "SL", 0.020);
  plot_analysis_MW("DL", "(type<=6)", "", "DL", 0.020);

}



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
 

  TH1F* h = new TH1F("h"," ; parton energy (GeV) ; transfer function (1/GeV)",1,5,225);
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
    c1->SaveAs("Plots/WS/Plot_TF_"+CP+"_"+TString(Form("%d",isB))+".png");
  }

  return;
}

void plot_param_3All(){
  
  plot_param_3(1, "Heavy-jets: quark #rightarrow reco-jet",   "TEST");
  plot_param_3(1, "Heavy-jets: gen-jet #rightarrow reco-jet", "TEST_std_gen");

  plot_param_3(0, "Light-jets: quark #rightarrow reco-jet",   "TEST");
  plot_param_3(0, "Light-jets: gen-jet #rightarrow reco-jet", "TEST_std_gen");
    

}
