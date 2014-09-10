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
#include "TFormula.h"
#include "TTree.h"
#include "TTreeFormula.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TROOT.h"
#include "TRandom3.h"

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
#include "TChain.h"

#include "TIterator.h"

#define SAVEPLOTS  1
#define SAVEHISTOS 0
#define VERBOSE    1
#define ASSUMEFLAT 0

typedef struct 
{
  int run;
  int lumi;
  int event;
  int json;
} EventInfo;


void get_migration_matrix(TString name = "fullSysts",
			  TString systname    = "CMS_scale_j", 
			  TString cutname = "cat1", 
			  TString cut = "( type==0 || (type==3 && flag_type3>0))   && btag_LR>=0.", 
			  double fact = 0.0189667*1.2 ){

  cout << "Now doing " << string(cutname.Data()) << endl ;
  
  vector<TString> files;
  //files.push_back( "./files/byLLR/Apr23_2014/MEAnalysisNew_all_fullSysts_rec_std_TTH125.root" );
  //files.push_back( "./files/byLLR/Apr23_2014/MEAnalysisNew_all_rec_std_TTV.root" );
  files.push_back( "./files/byLLR/Apr23_2014/MEAnalysisNew_all_fullSysts_rec_std_TTJets.root");

  TFile* fout = TFile::Open("2Dtemplates_matrix_"+name+".root","UPDATE");

  const int bins = 5;
  float limits[ bins+1 ] = {0., 0.166667, 0.333333, 0.5, 0.666667, 1. };
  TH2F* histosUp   = new TH2F(systname+"_"+cutname+"_Up",   "; obs ; obs_Up",   bins, limits, bins, limits ); 
  TH2F* histosDown = new TH2F(systname+"_"+cutname+"_Down", "; obs ; obs_Down", bins, limits, bins, limits ); 

  histosUp  ->Sumw2();
  histosDown->Sumw2();

  Long64_t pass = 0;
  Long64_t filled = 0;

  TFile* f = 0;
  
  for(int file = 0 ; file<files.size() ; file++){

    f = TFile::Open( files[file], "READ");
    TTree* t = (TTree*) f->Get("tree");

    TTreeFormula* treeformula = new TTreeFormula("cat_selection", cut , t );
    
    float p_125_all_s;
    float p_125_all_b;
    int type;
    int syst;
    float weight;
    float PUweight;
    float trigger;
    float weightTopPt;
    float weightEle;
    float weightCSV[19];
    
    t->SetBranchAddress("p_125_all_s",        &p_125_all_s );
    t->SetBranchAddress("p_125_all_b",        &p_125_all_b );
    t->SetBranchAddress("type",               &type); 
    t->SetBranchAddress("syst",               &syst); 
    
    // event information
    EventInfo EVENT;
    t->SetBranchAddress("EVENT",        &EVENT);   
    t->SetBranchAddress("weight",       &weight);
    t->SetBranchAddress("PUweight",     &PUweight);
    t->SetBranchAddress("trigger",      &trigger);
    t->SetBranchAddress("weightTopPt",  &weightTopPt);
    t->SetBranchAddress("weightCSV",     weightCSV  );   
    t->SetBranchAddress("weightEle",    &weightEle);

    for (Long64_t i = 0; i < t->GetEntries() ; i++){
      
      t->GetEntry(i);
      
      t->LoadTree(i);
      if( treeformula-> EvalInstance() == 0){
	continue; 
      }
      
      float fillN = -99;
      float fillU = -99;
      float fillD = -99;
      float total_weight = 1.0;
      
      if(syst==0){

	pass++;
	
	int eventi = EVENT.event;
	int typei  = type;
	
	for (Long64_t j = i; j < i+5 && j <t->GetEntries() ; j++){
	  t->GetEntry(j);
	  
	  t->LoadTree(i);
	  if( treeformula-> EvalInstance() == 0){
	    continue; 
	  }

	  int eventj = EVENT.event;
	  int typej  = type;
	  
	  float obs =   1./(1+fact*p_125_all_b/p_125_all_s) ;
	  
	  if( j==i   && syst==0 && eventj==eventi /*&& typej==typei*/){	  	    
	    fillN   = obs;
	    total_weight =  1.0; /*weight*weightEle*PUweight*trigger*weightCSV[0]*/;
	  }

	  if( systname=="CMS_scale_j" ){
	    if( j==i+1 && syst==3 && eventj==eventi /*&& typej!=typei*/){
	      fillU   = obs;
	    }
	    if( j==i+2 && syst==4 && eventj==eventi /*&& typej==typei*/){
	      fillD   = obs;
	    }
	  }
	  if( systname=="CMS_res_j"){
	    if( j==i+3 && syst==5 && eventj==eventi /*&& typej!=typei*/){
	      fillU   = obs;
	    }
	    if( j==i+4 && syst==6 && eventj==eventi /*&& typej==typei*/){
	      fillD   = obs;
	    }
	  }

	}
      }
	
      if( fillN>=0. && 
	  fillU>=0. &&
	  fillD>=0. ) {
	filled++;
	histosUp   ->Fill( fillN, fillU, total_weight);
	histosDown ->Fill( fillN, fillD, total_weight);
      }
      
      
    }
   
  }
  
  cout << "  >>> Passing:  " << pass   << " events" << endl;
  cout << "  >>> Filled:  " << filled << " events" << endl;
  
  fout->cd();
  for(int x = 1 ; x <= histosUp->GetXaxis()->GetNbins() ; x++){
    float normUp_x   = (histosUp  ->ProjectionY("_py",x,x))->Integral();
    float normDown_x = (histosDown->ProjectionY("_py",x,x))->Integral();
    for(int y = 1 ; y <= histosUp->GetYaxis()->GetNbins() ; y++){
      histosUp  ->SetBinContent( x,y, histosUp  ->GetBinContent( x,y )/normUp_x );
      histosDown->SetBinContent( x,y, histosDown->GetBinContent( x,y )/normDown_x );
   
      histosUp  ->SetBinError( x,y, histosUp  ->GetBinError( x,y )/normUp_x );
      histosDown->SetBinError( x,y, histosDown->GetBinError( x,y )/normDown_x );

      histosUp  ->Write(histosUp->GetName(),  TObject::kOverwrite);
      histosDown->Write(histosDown->GetName(),TObject::kOverwrite);
    }
  }
  fout->Close();

}

void get_migration_matrix_all(TString name = "fullSysts", TString syst = "CMS_scale_j", int average = 0){

  get_migration_matrix(name,  syst,  "cat1_H",  "( type==0 || (type==3 && flag_type3>0))   && btag_LR>=0.995",                0.0189667 *1.2 );
  get_migration_matrix(name,  syst,  "cat1_L",  "( type==0 || (type==3 && flag_type3>0)) && btag_LR<0.995 && btag_LR>=0.960", 0.017771  *1.2 );
  get_migration_matrix(name,  syst,  "cat2_H",  "(( type==1 || (type==3 && flag_type3<=0) ) && btag_LR>=0.9925) || (type==2 && flag_type2<=999 && btag_LR>=0.995)", 0.0198814 *1.2 );
  get_migration_matrix(name,  syst,  "cat2_L",  "(( type==1 || (type==3 && flag_type3<=0) ) && btag_LR<0.9925 && btag_LR>=0.960) || (type==2 && flag_type2<=999 && btag_LR<0.995 && btag_LR>=0.970)", 0.0185351 *1.2 );
  get_migration_matrix(name,  syst,  "cat6_H",  "type==6 && btag_LR>=0.925",                  0.0119848 *1.2 );
  get_migration_matrix(name,  syst,  "cat6_L",  "type==6 && btag_LR<0.925 && btag_LR>=0.850", 0.00746372*1.2 );

  if(average){
    TFile* f2D = TFile::Open( "2Dtemplates_matrix_"+name+".root", "UPDATE");
    vector<TString> cats;
    cats.push_back( "cat1");
    cats.push_back( "cat2");
    cats.push_back( "cat6");

    for( unsigned int c = 0 ; c < cats.size() ; c++){
      TH2F* pdf_Up_0 = (TH2F*)f2D->Get(syst+"_"+cats[c]+"_H_Up");
      TH2F* pdf_Up_1 = (TH2F*)f2D->Get(syst+"_"+cats[c]+"_L_Up");
      TH2F* pdf_Up   = (TH2F*) pdf_Up_0->Clone(syst+"_"+cats[c]+"_Up");
      pdf_Up->Add(pdf_Up_1);
      pdf_Up->Scale(0.5);      
      pdf_Up->Write(pdf_Up->GetName(), TObject::kOverwrite);

      TH2F* pdf_Down_0 = (TH2F*)f2D->Get(syst+"_"+cats[c]+"_H_Down");
      TH2F* pdf_Down_1 = (TH2F*)f2D->Get(syst+"_"+cats[c]+"_L_Down");
      TH2F* pdf_Down   = (TH2F*) pdf_Down_0->Clone(syst+"_"+cats[c]+"_Down");
      pdf_Down->Add(pdf_Down_1);
      pdf_Down->Scale(0.5);      
      pdf_Down->Write(pdf_Down->GetName(), TObject::kOverwrite);
    }
    f2D->Close();
  }
}


void apply_matrix(TString name = "fullSysts",
		  int rescale = 1, int ratio = 1, 
		  TString cat_    = "MEM_cat2_H",
		  TString dirpdf_ = "cat2",
		  TString syst_   = "CMS_res_j",
		  string bkg_     = "ttbarPlusBBbar",
		  float n_max_sigma = 1.0, 
		  float max_rel_err = 0.10){

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);


  TFile* f2D = TFile::Open( "2Dtemplates_matrix_"+name+".root", "READ"); 	  
 
  //TString  path = "./datacards/PostApproval/MyTest_new_merged_smooth/";
  TString  path = "./datacards/PostApproval/2D_OptFactors/";

  vector<TString> files;
  files.push_back("MEM_New_rec_std_sb.root");

  vector<TString> dirs;
  dirs.push_back(cat_);
  
  vector<TString> dirspdf;
  dirspdf.push_back(dirpdf_);

  vector<TString> systs;
  systs.push_back(syst_);

  vector<string> proc;
  //proc.push_back("ttbarV");
  //proc.push_back("ttH_hbb");
  //proc.push_back("singlet");
  //proc.push_back("ttbar");
  //proc.push_back("ttbarPlusB");
  proc.push_back(bkg_);
  //proc.push_back("ttbarPlusCCbar");


  TH1F* h0    = 0;
  TH1F* hUp   = 0;
  TH1F* hDown = 0;

  TH1F* hUp0   = 0;
  TH1F* hDown0 = 0;

  TFile* file = 0;

 for(unsigned int f = 0; f < files.size(); f++){

   file = TFile::Open( path+files[f], "UPDATE" );
    
    for( unsigned int d = 0; d < dirs.size(); d++ ){
      file->cd( dirs[d] ); 
      
      TIter iter( gDirectory->GetListOfKeys()  );
      TKey *key;
  
      // merge all bins
      while ( (key = (TKey*)iter.Next()) ){   
	
	TH1F* hist   =  (TH1F*)key->ReadObj(); 
	if(hist==0){
	  cout << "Null" << endl;
	  continue;
	}	

	string hname = string(hist->GetName());
	int nBins    = hist->GetNbinsX();
	float integ  = hist->Integral();

	if( hname.find("Up")!=string::npos || hname.find("Down")!=string::npos || hname.find("data_obs")!=string::npos)
	  continue;

	int isAmong = 0;
	for( unsigned int p = 0; p < proc.size(); p++){
	  if( hname == proc[p] ) isAmong++;
	}
	if(isAmong<1) continue;

	if(VERBOSE) cout << "Histo " << hname << " with " << nBins << endl;		
	//h0 = hist;

	for( unsigned int s = 0; s < systs.size(); s++  ){


	  TH2F* pdf_Up = (TH2F*)f2D->Get(systs[s]+"_"+dirspdf[d]+"_Up");
	  if( pdf_Up ==0 ){
	    cout << "Cannot find histo up " << string( (systs[s]+"_"+dirspdf[d]+"_Up").Data()) << endl;
	    continue;
	  }
	  TH2F* pdf_Down = (TH2F*)f2D->Get(systs[s]+"_"+dirspdf[d]+"_Down");
	  if( pdf_Down ==0 ){
	    cout << "Cannot find histo up " << string( (systs[s]+"_"+dirspdf[d]+"_Down").Data()) << endl;
	    continue;
	  }
	  
	  hUp0   = (TH1F*)gDirectory->Get( (hname+"_"+systs[s]+"Up"));
	  hDown0 = (TH1F*)gDirectory->Get( (hname+"_"+systs[s]+"Down"));
	  if( hUp0==0 || hDown0==0 ){
	    cout << "Cannot find" << endl;
	    continue;
	  }
	  double ratioUp   = hist->Integral()>0 ? hUp0->Integral()/hist->Integral() : 1.;
	  double ratioDown = hist->Integral()>0 ? hDown0->Integral()/hist->Integral() : 1.;

	  if(VERBOSE) cout << "Ratio Up: " << ratioUp << " --- Ratio Down: " << ratioDown << endl;

	  h0    = new TH1F((hname+"_"+systs[s]+"_copy"),     "", nBins, (*(hist->GetXaxis()->GetXbins())).GetArray() );
	  hUp   = new TH1F((hname+"_"+systs[s]+"Up_copy"),   "", nBins, (*(hist->GetXaxis()->GetXbins())).GetArray() );
	  hDown = new TH1F((hname+"_"+systs[s]+"Down_copy"), "", nBins, (*(hist->GetXaxis()->GetXbins())).GetArray() );

	  hUp->SetLineColor(kRed);    hUp->SetLineWidth(2);
	  hDown->SetLineColor(kBlue); hDown->SetLineWidth(2);

	  double normUp   = 0.;
	  double normDown = 0.;	

	  for(int b = 1 ; b <= nBins ; b++){

	    double int_b = hist->GetBinContent(b);
	    double p_b   = int_b/integ;

	    double c_Up   = 0.;
	    double c_Down = 0.;

	    double cErr2_Up   = 0.;
	    double cErr2_Down = 0.;

	    float xLow  = hist->GetXaxis()->GetBinLowEdge(b);
	    float xHigh = hist->GetXaxis()->GetBinUpEdge(b);
	    if(b>(nBins/2)){
	      xLow -=1.;
	      xHigh-=1.;
	    }
	    if(VERBOSE) cout << b << "th bin: [" << xLow << "," << xHigh << "], p=" << p_b ;

	    double norm0Up   = 0.;
	    double norm0Down = 0.;

	    for(int b0 = 1 ; b0 <=nBins ; b0++){

	      double p_b0       = hist->GetBinContent(b0)/integ;
	      double pErr2_b0   = hist->GetBinError(b0)*hist->GetBinError(b0)/integ/integ;

	      if(b<=(nBins/2) && b0>(nBins/2))  continue;
	      if(b>(nBins/2)  && b0<=(nBins/2)) continue;
	      if( b<=(nBins/2) ){
		c_Up   += p_b0*pdf_Up  ->GetBinContent(b0,b);
		c_Down += p_b0*pdf_Down->GetBinContent(b0,b);
		norm0Up   += pdf_Up  ->GetBinContent(b0,b);
		norm0Down += pdf_Down->GetBinContent(b0,b);
		cErr2_Up   +=  pErr2_b0*pdf_Up  ->GetBinContent(b0,b)*pdf_Up  ->GetBinContent(b0,b);
		cErr2_Down +=  pErr2_b0*pdf_Down->GetBinContent(b0,b)*pdf_Down->GetBinContent(b0,b);
	      }
	      else{
		c_Up   += p_b0*pdf_Up  ->GetBinContent(b0-(nBins/2),b-(nBins/2));
		c_Down += p_b0*pdf_Down->GetBinContent(b0-(nBins/2),b-(nBins/2));
		norm0Up   += pdf_Up  ->GetBinContent(b0-(nBins/2),b-(nBins/2));
		norm0Down += pdf_Down->GetBinContent(b0-(nBins/2),b-(nBins/2));
		cErr2_Up   +=  pErr2_b0*pdf_Up  ->GetBinContent(b0-(nBins/2),b-(nBins/2))*pdf_Up  ->GetBinContent(b0-(nBins/2),b-(nBins/2));
		cErr2_Down +=  pErr2_b0*pdf_Down->GetBinContent(b0-(nBins/2),b-(nBins/2))*pdf_Down->GetBinContent(b0-(nBins/2),b-(nBins/2));
	      }
		
	    }

	    if(VERBOSE) cout << ": (" << c_Up << ", " << c_Down << ")" << ". From matrix" << ": (" << norm0Up << ", " << norm0Down << ")" << endl;
	    normUp   += c_Up   ;
	    normDown += c_Down ;

	    h0     ->SetBinContent(b, int_b);
	    hUp    ->SetBinContent(b, rescale ? (ratioUp   -1 + (p_b>0 ? c_Up/p_b  : 1.)  )  *int_b : (p_b>0 ? c_Up/p_b : 1.0)  *int_b );
	    hDown  ->SetBinContent(b, rescale ? (ratioDown -1 + (p_b>0 ? c_Down/p_b: 1.)  )  *int_b : (p_b>0 ? c_Down/p_b: 1.0) *int_b );

	    hUp    ->SetBinError(b, sqrt(cErr2_Up)  *int_b );
	    hDown  ->SetBinError(b, sqrt(cErr2_Down)*int_b );
	   
	  }// bins	

	  if(VERBOSE){
	    cout << "Normalization: Up:"         << normUp << ", Down:" << normDown << endl;
	    cout << "Normalization (Histo): Up:" << hUp->Integral()/hist->Integral() << ", Down:" << hDown->Integral()/hist->Integral() << endl;
	  }


	  ////////////////
	  float maxDelta = 0.;
	  int maxBin     = 0;
	  float biasUp    = hUp->Integral()/h0->Integral(); 
	  float biasDown  = hDown->Integral()/h0->Integral(); 
	  for(int b = 1 ; b <=hUp->GetNbinsX() ; b++){
	    float bin       = h0->GetBinContent(b);
	    float binUp     = hUp->GetBinContent(b);
	    float binDown   = hDown->GetBinContent(b);
	    float binErr    = h0->GetBinError(b);  
	    float deltaUp   = binErr>0 ? TMath::Abs(binUp-bin*biasUp)/sqrt(binErr*binErr) : 0.;
	    float deltaDown = binErr>0 ? TMath::Abs(binDown-bin*biasDown)/sqrt(binErr*binErr) : 0.;
	    if(  deltaUp>maxDelta || deltaDown>maxDelta){
	      maxDelta = TMath::Max( deltaUp, deltaDown);
	      maxBin = b;
	    }
	  }
	  //if( ASSUMEFLAT && !(maxDelta > n_max_sigma && h0->GetBinError(maxBin)/h0->GetBinContent(maxBin) < max_rel_err) ){
	  if( ASSUMEFLAT && !((cat_=="MEM_cat2_H" || cat_=="MEM_cat2_L") && bkg_!="singlet" ) ){
	    cout << "Assume flat systematics, because max deviation from flat is " << maxDelta << " s.d." << endl;
	    for(int b = 1 ; b <=hUp->GetNbinsX() ; b++){
	      hUp  ->SetBinContent(b, h0->GetBinContent(b)*biasUp);
	      hUp  ->SetBinError  (b, h0->GetBinError  (b)*biasUp);
	      hDown  ->SetBinContent(b, h0->GetBinContent(b)*biasDown);
	      hDown  ->SetBinError  (b, h0->GetBinError  (b)*biasDown);
	    }
	  }
	  /////////


	  if(SAVEHISTOS){
	    file->cd( dirs[d] ); 
	    hUp   ->Write(hname+"_"+systs[s]+"Up",    TObject::kOverwrite);
	    hDown ->Write(hname+"_"+systs[s]+"Down",  TObject::kOverwrite);
	  }

	  h0 = hist;

	}	

      }

    } // directories


 }// files




 hUp0   ->SetLineColor(kRed);
 hDown0->SetLineColor(kBlue);
 hUp0   ->SetMarkerStyle(kFullCircle);
 hDown0 ->SetMarkerStyle(kFullCircle);

 hUp->SetStats(0);
 hUp->SetTitle(cat_+", "+syst_+" --> "+bkg_);
 h0->SetStats(0);
 h0->SetTitle(cat_+", "+syst_+" --> "+bkg_);

 h0->SetMinimum(0.);
 h0->SetMaximum(h0->GetMaximum()*1.35);
 h0   ->Draw();
 hUp  ->Draw("SAME");
 hDown->Draw("SAME");
 hUp0  ->Draw("PSAME");
 hDown0->Draw("PSAME");

 if(ratio){

   hUp->Divide(h0);
   hDown->Divide(h0);

   //hUp0->Sumw2();
   //hDown0->Sumw2();
   //h0->Sumw2();

   hUp0->Divide(hUp0, h0, 1,1);
   hDown0->Divide(hDown0, h0, 1,1);

   //h0->Divide(h0, h0, 1,1,"B");
   for(int b0 = 1 ; b0 <=h0->GetNbinsX() ; b0++){
     float binc = h0->GetBinContent(b0);
     float bine = h0->GetBinError(b0);
     h0->SetBinContent(b0, 1.);
     h0->SetBinError(b0, bine/binc);
   }

   hUp->SetMaximum(2.0);
   hUp->SetMinimum(0.5);

   //hUp  ->SetFillColor(kRed);
   //hUp  ->SetFillStyle(3004);
   //hDown->SetFillColor(kBlue);
   //hDown->SetFillStyle(3005);

   //hUp  ->Draw("E3");
   //hDown->Draw("E3SAME");
   hUp  ->Draw("HIST");
   hDown->Draw("HISTSAME");


   hUp0  ->Draw("ESAME");
   hDown0->Draw("ESAME");

   h0->SetFillColor(kGreen);
   h0->SetFillStyle(3002);
   h0->Draw("E3SAME");
 }

 if(SAVEPLOTS){
   c1->SaveAs(path+"/Comp_2D_"+cat_+"_"+syst_+"_"+bkg_+".pdf");
   file->Close();
 }
 
 f2D->Close();

 if(SAVEHISTOS)  file->Close();

 delete c1;

 return;
}

void saveAll(TString name = "fullSysts"){
  
 vector<TString> dirs;
 dirs.push_back("MEM_cat1_H");
 dirs.push_back("MEM_cat1_L");
 dirs.push_back("MEM_cat2_H");
 dirs.push_back("MEM_cat2_L");
 dirs.push_back("MEM_cat3_H");
 dirs.push_back("MEM_cat3_L");
 dirs.push_back("MEM_cat6_H");
 dirs.push_back("MEM_cat6_L");
  
 vector<TString> dirspdf;
 dirspdf.push_back("cat1_H");
 dirspdf.push_back("cat1_L");
 dirspdf.push_back("cat2_H");
 dirspdf.push_back("cat2_L");
 dirspdf.push_back("cat2_H");
 dirspdf.push_back("cat2_L");
 dirspdf.push_back("cat6");
 dirspdf.push_back("cat6");
 
 
 vector<TString> systs;
 systs.push_back( "CMS_scale_j" );
 systs.push_back( "CMS_res_j" );
 
 vector<string> proc;
 proc.push_back("ttbarV");
 proc.push_back("ttH_hbb");
 proc.push_back("singlet");
 proc.push_back("ttbar");
 proc.push_back("ttbarPlusB");
 proc.push_back("ttbarPlusBBbar");
 proc.push_back("ttbarPlusCCbar");
 
  for( unsigned int d = 0; d < dirs.size(); d++ ){
    for( unsigned int s = 0; s < systs.size(); s++  ){
      for( unsigned int p = 0; p < proc.size(); p++){
	apply_matrix(name, 1,1,dirs[d], dirspdf[d], systs[s], proc[p]);
      }
    }
  }

}


void get_shape_cumulative( TString syst_   = "CMS_res_j", 
			   TString cumName = "bkg", 
			   TString cat = "cat1_H",  
			   int ratio = 1){

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);
 
  string syststr(syst_.Data());
  string cumNamestr(cumName.Data());
  string catstr(cat.Data());

  TString  path = "./datacards/PostApproval/MyTest_new_fullSysts/";

  vector<TString> files;
  files.push_back("MEM_New_rec_std_sb.root");

  vector<TString> dirs;
  if(catstr=="cat1_H"){
    dirs.push_back("MEM_cat1_H");
    dirs.push_back("MEM_cat1_L");  
    dirs.push_back("MEM_cat2_H");
    dirs.push_back("MEM_cat2_L");
  }
  if(catstr=="cat1_L"){
    dirs.push_back("MEM_cat1_H");
    dirs.push_back("MEM_cat1_L");  
    dirs.push_back("MEM_cat2_H");
    dirs.push_back("MEM_cat2_L");
  }
  if(catstr=="cat2_H"){
    dirs.push_back("MEM_cat1_H");
    dirs.push_back("MEM_cat1_L");  
    dirs.push_back("MEM_cat2_H");
    dirs.push_back("MEM_cat2_L");
  }
  if(catstr=="cat2_L"){
    dirs.push_back("MEM_cat1_H");
    dirs.push_back("MEM_cat1_L");  
    dirs.push_back("MEM_cat2_H");
    dirs.push_back("MEM_cat2_L");
  }
  if(catstr=="cat6_H"){
    dirs.push_back("MEM_cat6_H");
    dirs.push_back("MEM_cat6_L");
  }
  if(catstr=="cat6_L"){
    dirs.push_back("MEM_cat6_H");
    dirs.push_back("MEM_cat6_L");
  }
  
  TH1F* h0    = 0;
  TH1F* hUp   = 0;
  TH1F* hDown = 0;

  float Nominal   = -1.;
  float Up        = -1.;
  float Down      = -1.;

  float NominalCat   = 0.;
  float UpCat        = 0.;
  float DownCat      = 0.;

  vector<string> proc;
  if(cumNamestr=="ttH_hbb" && catstr.find("cat6")==string::npos ){
    proc.push_back("ttH_hbb");
    //proc.push_back("ttbarV");
    Nominal = 0;
    Up      = 0;
    Down    = 0;
  }
  if(cumNamestr=="ttH_hbb" && catstr.find("cat6")!=string::npos ){
    proc.push_back("ttbar");
    proc.push_back("ttbarPlusB");
    proc.push_back("ttbarPlusBBbar");
    proc.push_back("ttbarPlusCCbar");
    Nominal = 0;
    Up      = 0;
    Down    = 0;
  }
  if(cumNamestr=="ttbarV"){
    proc.push_back("ttbar");
    proc.push_back("ttbarPlusB");
    proc.push_back("ttbarPlusBBbar");
    proc.push_back("ttbarPlusCCbar");
    Nominal = 0;
    Up      = 0;
    Down    = 0;
  }
  if(cumNamestr=="singlet"){
    proc.push_back("ttbar");
    proc.push_back("ttbarPlusB");
    proc.push_back("ttbarPlusBBbar");
    proc.push_back("ttbarPlusCCbar");
    Nominal = 0;
    Up      = 0;
    Down    = 0;
  }
  if(cumNamestr=="ttbar"){
    proc.push_back("ttbar");
    proc.push_back("ttbarPlusB");
    proc.push_back("ttbarPlusBBbar");
    proc.push_back("ttbarPlusCCbar");
  }

  TFile* file = 0;

  for(unsigned int f = 0; f < files.size(); f++){

    file = TFile::Open( path+files[f], "READ" );

    for( unsigned int d = 0; d < dirs.size(); d++ ){
      file->cd( dirs[d] ); 
      
      TIter iter( gDirectory->GetListOfKeys()  );
      TKey *key;

      // merge all bins
      while ( (key = (TKey*)iter.Next()) ){   
	
	TH1F* hist   =  (TH1F*)key->ReadObj(); 
	if(hist==0){
	  cout << "Null" << endl;
	  continue;
	}	
	
	string hname = string(hist->GetName());
	int nBins    = hist->GetNbinsX();
	float integ  = hist->Integral();

	int isAmong = 0;
	string process = "";
	for( unsigned int p = 0; p < proc.size(); p++){
	  if( hname.find(proc[p])!=string::npos ){
	    isAmong++;
	    process = proc[p];
	  }
	}

	if( dirs[d]==("MEM_"+cat) ){
	  if(  hname ==  cumNamestr                     && Nominal> 0 ) Nominal += hist->Integral();
	  if(  hname == (cumNamestr+"_"+syststr+"Up")   &&    Up  > 0 ) Up      += hist->Integral();
	  if(  hname == (cumNamestr+"_"+syststr+"Down") &&  Down  > 0 ) Down    += hist->Integral();
	}


	if(isAmong<1) continue;

	//if(VERBOSE) cout << "Histo " << hname << " with " << nBins << ".... " ;

	if( hname == process  ){

	  if(h0==0){
	    h0  = new TH1F("cumulative_"+cumName, "", nBins, (*(hist->GetXaxis()->GetXbins())).GetArray() );
	    h0->Add( hist, 1.0);
	  }
	  else
	    h0->Add( hist, 1.0);

	  if( dirs[d]==("MEM_"+cat) ) NominalCat +=  hist->Integral();
	  //cout << hname << ": filled h0 with " << hist->Integral() << endl;
	  continue;
	}

	if( hname == (process+"_"+syststr+"Up")  ){

	  if(hUp==0){
	    hUp  = new TH1F(("cumulative_"+cumName+"_"+syststr+"Up"), "", nBins, (*(hist->GetXaxis()->GetXbins())).GetArray() );
	    hUp->Add( hist, 1.0);
	  }
	  else
	    hUp->Add( hist, 1.0);
	  if( dirs[d]==("MEM_"+cat) ) UpCat +=  hist->Integral();
	  //cout << hname << ": filled hUp with " << hist->Integral() << endl;
	  continue;
	}

	if( hname == (process+"_"+syststr+"Down")  ){

	  if( process == cumNamestr && Down > 0 ) Down = hist->Integral();

	  if(hDown==0){
	    hDown  = new TH1F(("cumulative_"+cumName+"_"+syststr+"Down"), "", nBins, (*(hist->GetXaxis()->GetXbins())).GetArray() );
	    hDown->Add( hist, 1.0);
	  }
	  else
	    hDown->Add( hist, 1.0);
	  if( dirs[d]==("MEM_"+cat) ) DownCat +=  hist->Integral();
	  //cout << hname << ": filled hDown with " << hist->Integral() << endl;
	  continue;
	}

	//cout << endl;

      }
      
    } // dirs

    
  } // files

  for(int b = 1 ; b <=h0->GetNbinsX() ; b++){
    float bin0      = h0->GetBinContent(b);
    float binUp     = hUp->GetBinContent(b);
    float binDown   = hDown->GetBinContent(b);
    float ratUp     = bin0>0 ? binUp/bin0 : 1.;
    if(ratUp<1)   ratUp   = 1./ratUp;
    float ratDown   = bin0>0 ? binDown/bin0 : 1.;
    if(ratDown<1) ratDown = 1./ratDown;
    double rat = TMath::Max( ratUp, ratDown );

    if( bin0>0. ){
      hUp  ->SetBinContent(b, bin0* (binUp>binDown ? rat    : 1./rat)     );
      hDown->SetBinContent(b, bin0* (binUp>binDown ? 1./rat : rat)        );
    }
    else{
      hUp    ->SetBinContent(b, bin0);
      hDown  ->SetBinContent(b, bin0);
    }
  }

  if(Nominal>0){
    cout << "Ratio Up: " << Up/Nominal << ", ratio Down: " << Down/Nominal << endl;
    hUp   ->Scale( Up/Nominal  *h0->Integral()/hUp->Integral() );
    hDown ->Scale( Down/Nominal*h0->Integral()/hDown->Integral() );
  }
  else{
    cout << "Ratio Up: " << UpCat/NominalCat << ", ratio Down: " << DownCat/NominalCat << endl;
    hUp   ->Scale( UpCat/NominalCat  *h0->Integral()/hUp->Integral() );
    hDown ->Scale( DownCat/NominalCat*h0->Integral()/hDown->Integral() );
  }

  
  TFile* out = TFile::Open(path+"Map_systs.root", "UPDATE");
  out->cd();
  TH1F* hUp_clone   = (TH1F*)h0->Clone( hUp->GetName() );
  TH1F* hDown_clone = (TH1F*)h0->Clone( hDown->GetName() );
  hUp_clone  ->Divide( hUp,   h0, 1, 1);
  hDown_clone->Divide( hDown, h0, 1, 1);  
  //h0    ->Write(cat+"_"+cumName,                   TObject::kOverwrite);
  hUp_clone   ->Write(cat+"_"+cumName+"_"+syst_+"Up",    TObject::kOverwrite);
  hDown_clone ->Write(cat+"_"+cumName+"_"+syst_+"Down",  TObject::kOverwrite);
  out->Close();
  

  cout << "Normalization (Histo): " << h0->Integral() << ", Up:" << hUp->Integral() << ", Down:" << hDown->Integral() << endl;
  hUp->SetLineColor(kRed);    hUp->SetLineWidth(2);
  hDown->SetLineColor(kBlue); hDown->SetLineWidth(2);

  h0->SetMinimum(0.);
  h0->SetMaximum(h0->GetMaximum()*1.35);
  h0   ->Draw();
  hUp  ->Draw("SAME");
  hDown->Draw("SAME");

  
  


  //hUp  ->Scale(h0->Integral()/hUp->Integral());
  //hDown->Scale(h0->Integral()/hDown->Integral());

  if(ratio){

    hUp->Divide(h0);
    hDown->Divide(h0);

    for(int b0 = 1 ; b0 <=h0->GetNbinsX() ; b0++){
      float binc = h0->GetBinContent(b0);
      float bine = h0->GetBinError(b0);
      h0->SetBinContent(b0, 1.);
      h0->SetBinError(b0, bine/binc);
    }

    hUp->SetMaximum(1.5);
    hUp->SetMinimum(0.5);

    //hUp  ->SetFillColor(kRed);
    //hUp  ->SetFillStyle(3004);
    //hDown->SetFillColor(kBlue);
    //hDown->SetFillStyle(3005);    
    //hUp  ->Draw("E3");
    //hDown->Draw("E3SAME");

    hUp  ->Draw("HIST");
    hDown->Draw("HISTSAME");

    h0->SetFillColor(kGreen);
    h0->SetFillStyle(3002);
    h0->Draw("E3SAME");
  }
  
  if(SAVEPLOTS){
    c1->SaveAs(path+"/Cum_"+syst_+"_"+cumName+"_"+cat+".pdf");
    file->Close();
  }
  
  if(SAVEHISTOS)  file->Close();
  
  delete c1;
  
  return;

}


void tes_all(){

  // ttbar
  get_shape_cumulative("CMS_scale_j", "ttbar", "cat1_H", 1);
  get_shape_cumulative("CMS_scale_j", "ttbar", "cat1_L", 1);

  get_shape_cumulative("CMS_scale_j", "ttbar", "cat2_H", 1);
  get_shape_cumulative("CMS_scale_j", "ttbar", "cat2_L", 1);
  
  get_shape_cumulative("CMS_scale_j", "ttbar", "cat6_H", 1);
  get_shape_cumulative("CMS_scale_j", "ttbar", "cat6_L", 1);

  get_shape_cumulative("CMS_res_j", "ttbar", "cat1_H", 1);
  get_shape_cumulative("CMS_res_j", "ttbar", "cat1_L", 1);

  get_shape_cumulative("CMS_res_j", "ttbar", "cat2_H", 1);
  get_shape_cumulative("CMS_res_j", "ttbar", "cat2_L", 1);
  
  get_shape_cumulative("CMS_res_j", "ttbar", "cat6_H", 1);
  get_shape_cumulative("CMS_res_j", "ttbar", "cat6_L", 1);

  // ttH
  get_shape_cumulative("CMS_scale_j", "ttH_hbb", "cat1_H", 1);
  get_shape_cumulative("CMS_scale_j", "ttH_hbb", "cat1_L", 1);

  get_shape_cumulative("CMS_scale_j", "ttH_hbb", "cat2_H", 1);
  get_shape_cumulative("CMS_scale_j", "ttH_hbb", "cat2_L", 1);
  
  get_shape_cumulative("CMS_scale_j", "ttH_hbb", "cat6_H", 1);
  get_shape_cumulative("CMS_scale_j", "ttH_hbb", "cat6_L", 1);

  get_shape_cumulative("CMS_res_j", "ttH_hbb", "cat1_H", 1);
  get_shape_cumulative("CMS_res_j", "ttH_hbb", "cat1_L", 1);

  get_shape_cumulative("CMS_res_j", "ttH_hbb", "cat2_H", 1);
  get_shape_cumulative("CMS_res_j", "ttH_hbb", "cat2_L", 1);
  
  get_shape_cumulative("CMS_res_j", "ttH_hbb", "cat6_H", 1);
  get_shape_cumulative("CMS_res_j", "ttH_hbb", "cat6_L", 1);

  // ttbarV
  get_shape_cumulative("CMS_scale_j", "ttbarV", "cat1_H", 1);
  get_shape_cumulative("CMS_scale_j", "ttbarV", "cat1_L", 1);

  get_shape_cumulative("CMS_scale_j", "ttbarV", "cat2_H", 1);
  get_shape_cumulative("CMS_scale_j", "ttbarV", "cat2_L", 1);
  
  get_shape_cumulative("CMS_scale_j", "ttbarV", "cat6_H", 1);
  get_shape_cumulative("CMS_scale_j", "ttbarV", "cat6_L", 1);

  get_shape_cumulative("CMS_res_j", "ttbarV", "cat1_H", 1);
  get_shape_cumulative("CMS_res_j", "ttbarV", "cat1_L", 1);

  get_shape_cumulative("CMS_res_j", "ttbarV", "cat2_H", 1);
  get_shape_cumulative("CMS_res_j", "ttbarV", "cat2_L", 1);
  
  get_shape_cumulative("CMS_res_j", "ttbarV", "cat6_H", 1);
  get_shape_cumulative("CMS_res_j", "ttbarV", "cat6_L", 1);

  // singlet
  get_shape_cumulative("CMS_scale_j", "singlet", "cat1_H", 1);
  get_shape_cumulative("CMS_scale_j", "singlet", "cat1_L", 1);

  get_shape_cumulative("CMS_scale_j", "singlet", "cat2_H", 1);
  get_shape_cumulative("CMS_scale_j", "singlet", "cat2_L", 1);
  
  get_shape_cumulative("CMS_scale_j", "singlet", "cat6_H", 1);
  get_shape_cumulative("CMS_scale_j", "singlet", "cat6_L", 1);

  get_shape_cumulative("CMS_res_j", "singlet", "cat1_H", 1);
  get_shape_cumulative("CMS_res_j", "singlet", "cat1_L", 1);

  get_shape_cumulative("CMS_res_j", "singlet", "cat2_H", 1);
  get_shape_cumulative("CMS_res_j", "singlet", "cat2_L", 1);
  
  get_shape_cumulative("CMS_res_j", "singlet", "cat6_H", 1);
  get_shape_cumulative("CMS_res_j", "singlet", "cat6_L", 1);
  

}
