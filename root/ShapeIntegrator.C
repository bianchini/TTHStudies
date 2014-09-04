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

#define SAVE 1

typedef struct 
{
  int run;
  int lumi;
  int event;
  int json;
} EventInfo;


void get_migration_matrix(TString systname    = "CMS_scale_j", 
			  TString cutname = "cat1", 
			  TString cut = "( type==0 || (type==3 && flag_type3>0))   && btag_LR>=0.", 
			  double fact = 0.0189667*1.2 ){

  cout << "Now doing " << string(cutname.Data()) << endl ;
  
  vector<TString> files;
  files.push_back( "../../root/files/byLLR/Apr23_2014/MEAnalysisNew_all_rec_std_TTH125.root" );
  files.push_back( "../../root/files/byLLR/Apr23_2014/MEAnalysisNew_all_rec_std_TTV.root" );
  files.push_back( "../../root/files/byLLR/Apr23_2014/MEAnalysisNew_all_rec_std_TTJets_nC.root");

  TFile* fout = TFile::Open("2Dtemplates_matrix.root","UPDATE");

  const int bins = 5;
  float limits[ bins+1 ] = {0., 0.166667, 0.333333, 0.5, 0.666667, 1. };
  TH2F* histosUp   = new TH2F(systname+"_"+cutname+"_Up",   "; obs ; obs_Up",   bins, limits, bins, limits ); 
  TH2F* histosDown = new TH2F(systname+"_"+cutname+"_Down", "; obs ; obs_Down", bins, limits, bins, limits ); 

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
    
    
    t->SetBranchAddress("p_125_all_s",        &p_125_all_s );
    t->SetBranchAddress("p_125_all_b",        &p_125_all_b );
    t->SetBranchAddress("type",               &type); 
    t->SetBranchAddress("syst",               &syst); 
    
    // event information
    EventInfo EVENT;
    t->SetBranchAddress("EVENT",        &EVENT);   
    
    for (Long64_t i = 0; i < t->GetEntries() ; i++){
      
      t->GetEntry(i);
      
      t->LoadTree(i);
      if( treeformula-> EvalInstance() == 0){
	continue; 
      }
      
      float fillN = -99;
      float fillU = -99;
      float fillD = -99;
      
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
	histosUp   ->Fill( fillN, fillU);
	histosDown ->Fill( fillN, fillD);
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
      histosUp  ->Write(histosUp->GetName(),  TObject::kOverwrite);
      histosDown->Write(histosDown->GetName(),TObject::kOverwrite);
    }
  }
  fout->Close();

}

void get_migration_matrix_all(TString syst = "CMS_scale_j", int average = 0){

  get_migration_matrix( syst,  "cat1_H",  "( type==0 || (type==3 && flag_type3>0))   && btag_LR>=0.995",                0.0189667 *1.2 );
  get_migration_matrix( syst,  "cat1_L",  "( type==0 || (type==3 && flag_type3>0)) && btag_LR<0.995 && btag_LR>=0.960", 0.017771  *1.2 );
  get_migration_matrix( syst,  "cat2_H",  "(( type==1 || (type==3 && flag_type3<=0) ) && btag_LR>=0.9925) || (type==2 && flag_type2<=999 && btag_LR>=0.995)", 0.0198814 *1.2 );
  get_migration_matrix( syst,  "cat2_L",  "(( type==1 || (type==3 && flag_type3<=0) ) && btag_LR<0.9925 && btag_LR>=0.960) || (type==2 && flag_type2<=999 && btag_LR<0.995 && btag_LR>=0.970)", 0.0185351 *1.2 );
  get_migration_matrix( syst,  "cat6_H",  "type==6 && btag_LR>=0.925",                  0.0119848 *1.2 );
  get_migration_matrix( syst,  "cat6_L",  "type==6 && btag_LR<0.925 && btag_LR>=0.850", 0.00746372*1.2 );

  if(average){
    TFile* f2D = TFile::Open( "2Dtemplates_matrix.root", "UPDATE");
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


void apply_matrix(int rescale = 1, int ratio = 1, 
		  TString cat_    = "MEM_cat2_H",
		  TString dirpdf_ = "cat2",
		  TString syst_   = "CMS_res_j",
		  string bkg_    = "ttbarPlusBBbar"){

  TCanvas *c1 = new TCanvas("c1","",5,30,650,600);
  c1->SetGrid(0,0);
  c1->SetFillStyle(4000);
  c1->SetFillColor(10);
  c1->SetTicky();
  c1->SetObjectStat(0);


  TFile* f2D = TFile::Open( "2Dtemplates_matrix.root", "READ"); 	  
 
  TString  path = "../../root/datacards/PostApproval/MyTest_new_merged/";

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

	cout << "Histo " << hname << " with " << nBins << endl;		
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

	  cout << "Ratio Up: " << ratioUp << " --- Ratio Down: " << ratioDown << endl;

	  h0    = new TH1F((hname+"_"+systs[s]+"_copy"),     "", nBins, (*(hist->GetXaxis()->GetXbins())).GetArray() );
	  hUp   = new TH1F((hname+"_"+systs[s]+"Up_copy"),   "", nBins, (*(hist->GetXaxis()->GetXbins())).GetArray() );
	  hDown = new TH1F((hname+"_"+systs[s]+"Down_copy"), "", nBins, (*(hist->GetXaxis()->GetXbins())).GetArray() );

	  hUp->SetLineColor(kRed);    hUp->SetLineWidth(2);
	  hDown->SetLineColor(kBlue); hDown->SetLineWidth(2);

	  for(int b = 1 ; b <= nBins ; b++){

	    float p_b   = hist->GetBinContent(b)/integ;
	    float int_b = hist->GetBinContent(b);

	    float c_Up   = 0.;
	    float c_Down = 0.;

	    float xLow  = hist->GetXaxis()->GetBinLowEdge(b);
	    float xHigh = hist->GetXaxis()->GetBinUpEdge(b);
	    if(b>(nBins/2)){
	      xLow -=1.;
	      xHigh-=1.;
	    }
	    cout << b << "th bin: [" << xLow << "," << xHigh << "]" << endl;

	    for(int b0 = 1 ; b0 <=nBins ; b0++){
	      if(b<=(nBins/2) && b0>(nBins/2))  continue;
	      if(b>(nBins/2)  && b0<=(nBins/2)) continue;
	      if( b<=(nBins/2) ){
		c_Up   += p_b*pdf_Up  ->GetBinContent(b0,b);
		c_Down += p_b*pdf_Down->GetBinContent(b0,b);
	      }
	      else{
		c_Up   += p_b*pdf_Up  ->GetBinContent(b0-(nBins/2),b-(nBins/2));
		c_Down += p_b*pdf_Down->GetBinContent(b0-(nBins/2),b-(nBins/2));
	      }
		
	    }

	    h0   ->SetBinContent(b, int_b);
	    hUp    ->SetBinContent(b, rescale ? (ratioUp  +(c_Up/p_b-1))  *int_b : (c_Up/p_b)  *int_b );
	    hDown  ->SetBinContent(b, rescale ? (ratioDown+(c_Down/p_b-1))*int_b : (c_Down/p_b)*int_b );
	   
	  }// bins	 

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
   hUp  ->Draw("HIST");
   hDown->Draw("HISTSAME");
   hUp0  ->Draw("ESAME");
   hDown0->Draw("ESAME");

   h0->SetFillColor(kGreen);
   h0->SetFillStyle(3002);
   h0->Draw("E3SAME");
 }

 if(SAVE){
   c1->SaveAs("Comp_2D_"+cat_+"_"+syst_+"_"+bkg_+".pdf");
   file->Close();
 }
 
 f2D->Close();

 return;
}

void saveAll(){
  
 vector<TString> dirs;
 dirs.push_back("MEM_cat1_H");
 dirs.push_back("MEM_cat1_L");
 dirs.push_back("MEM_cat2_H");
 dirs.push_back("MEM_cat2_L");
 dirs.push_back("MEM_cat6_H");
 dirs.push_back("MEM_cat6_L");
  
 vector<TString> dirspdf;
 dirspdf.push_back("cat1_H");
 dirspdf.push_back("cat1_L");
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
 proc.push_back("ttbarPlusCCbar");
 
  for( unsigned int d = 0; d < dirs.size(); d++ ){
    for( unsigned int s = 0; s < systs.size(); s++  ){
      for( unsigned int p = 0; p < proc.size(); p++){
	apply_matrix(1,1,dirs[d], dirspdf[d], systs[s], proc[p]);
      }
    }
  }

}


/*



void get_templates( TString cutname = "cat6", TString cut = "(type==6)"){

  cout << "Now doing " << string(cutname.Data()) << endl ;

  TFile* f  = TFile::Open("../../root/files/byLLR/Apr23_2014/MEAnalysisNew_all_rec_std_TTH125.root");
  TTree* t  = (TTree*)f->Get("tree");

  TFile* fout = TFile::Open("2Dtemplates.root","UPDATE");

  const int bins = 5;
  float limits[ bins+1 ] = {0., 0.166667, 0.333333, 0.5, 0.666667, 1. };
  vector<TH2F*> histos2D;
  vector<TH1F*> histos1D_Up;
  vector<TH1F*> histos1D_Down;
  for(int i = 0 ; i < bins ; i++){
    TH2F* h2 = new TH2F(cutname+TString(Form("_2D_%d",i)), "; Up ; Down", 40, 0, 4 , 40, 0, 4 ); 
    histos2D.push_back( h2 );

    TH1F* hUp = new TH1F(cutname+TString(Form("_Up_%d",i)), "; Up", 40, 0, 4); 
    histos1D_Up.push_back( hUp );
    TH1F* hDown = new TH1F(cutname+TString(Form("_Down_%d",i)), "; Down", 40, 0, 4); 
    histos1D_Down.push_back( hDown );
  }
  
  TTreeFormula* treeformula = new TTreeFormula("cat_selection", cut , t );

  float p_125_all_s;
  float p_125_all_b;
  int type;
  int syst;

  t->SetBranchAddress("p_125_all_s",        &p_125_all_s );
  t->SetBranchAddress("p_125_all_b",        &p_125_all_b );
  t->SetBranchAddress("type",               &type); 
  t->SetBranchAddress("syst",               &syst); 

  // event information
  EventInfo EVENT;
  t->SetBranchAddress("EVENT",        &EVENT);

  Long64_t pass = 0;

  for (Long64_t i = 0; i < t->GetEntries() ; i++){

    t->GetEntry(i);

    t->LoadTree(i);
    if( treeformula-> EvalInstance() == 0){
      continue; 
    }
    pass++;

    float fillN = -99;
    float fillU = -99;
    float fillD = -99;

    if(syst==0){

      int eventi = EVENT.event;
      int typei  = type;

      int whichbin = -99;

      for (Long64_t j = i; j < i+5 && j <t->GetEntries() ; j++){
	t->GetEntry(j);

	int eventj = EVENT.event;
	int typej  = type;

	float o =   1./(1+0.0189667*1.2*p_125_all_b/p_125_all_s) ;

	if( syst==0 && eventj==eventi && typej==typei){

	  double logo = o;//TMath::Finite(o) ?  TMath::Log(o) : +998. ;
	  //if( logo<-999.) logo = -998.;

	  for(int b = 0 ; b < bins ; b++){
	    if( logo>= limits[b] && logo<limits[b+1]) whichbin = b;
	  }
	  if(whichbin<0){
	    cout << logo << ", " << o << ", " << p_125_all_b << ", " << p_125_all_s << endl;;
	    whichbin = bins-1;
	  }
	  
	  fillN   = o;
	}
	if( syst==3 && eventj==eventi && typej==typei){
	  fillU   = o;
	}
	if( syst==4 && eventj==eventi && typej==typei){
	  fillD   = o;
	}
      }
      
      if( fillN>=0. && 
	  fillU>=0. &&
	  fillD>=0. ) {
	if(whichbin>=0){
	  //if( TMath::Log(fillN)>10 ){
	  //histos[whichbin]->Fill( 1.0, 1.0 );
	  //cout << "o = " << fillN << ", o+ = " << fillU << ", o- = " << fillD << endl;
	  //}
	  //else
	  histos2D[whichbin]     ->Fill( (fillU)/fillN, (fillD)/fillN );
	  histos1D_Up[whichbin]  ->Fill( (fillU)/fillN );
	  histos1D_Down[whichbin]->Fill( (fillD)/fillN );
	}
	else{
	  cout << "... cannot ind this bin..." << endl;
	}
      }

    }

  }
  f->Close();

  cout << "  >>> Passing:  " << pass << " events" << endl;

  fout->cd();
  for(int i = 0 ; i < bins ; i++){
    TH2F* h2 = histos2D[i];
    cout << string(h2->GetName()) << " has " << h2->GetEntries() << " entries" << endl;
    h2->Write(h2->GetName(),TObject::kOverwrite);

    TH1F* hUp = histos1D_Up[i];
    hUp->Write(hUp->GetName(),TObject::kOverwrite);

    TH1F* hDown = histos1D_Down[i];
    hDown->Write(hDown->GetName(),TObject::kOverwrite);
  }
  fout->Close();

}


void get_templates_all(){

  get_templates("cat1", "( type==0 || (type==3 && flag_type3>0))   && btag_LR>=0.");
  //get_templates("cat2", "(( type==1 || (type==3 && flag_type3<=0) ) && btag_LR>=0.) || (type==2 && flag_type2<=999 && btag_LR>=0.)");
  //get_templates("cat6", "type==6 && btag_LR>=0."); 
}


void apply(int nmax = 10000){

  const int bins = 5;
  float limits[ bins+1 ] = {0., 0.166667, 0.333333, 0.5, 0.666667, 1. };


  TFile* f2D = TFile::Open( "2Dtemplates.root", "READ");
  
  TRandom3* ran = new TRandom3();

  TString  path = "../../root/datacards/PostApproval/2D_OptFactors/";

  vector<TString> files;
  files.push_back("MEM_New_rec_std_sb.root");

  vector<TString> dirs;
  dirs.push_back("MEM_cat1_H");
  
  vector<TString> dirspdf;
  dirspdf.push_back("cat1");

  vector<double> forms;
  forms.push_back(0.0189667*1.2);

  vector<TString> systs;
  systs.push_back("CMS_scale_j");

  vector<string> proc;
  //proc.push_back("ttbarV");
  //proc.push_back("ttH_hbb");
  //proc.push_back("singlet");
  //proc.push_back("ttbar");
  //proc.push_back("ttbarPlusB");
  proc.push_back("ttbarPlusBBbar");
  //proc.push_back("ttbarPlusCCbar");


  TH1F* h0    = 0;
  TH1F* hUp   = 0;
  TH1F* hDown = 0;

 for(unsigned int f = 0; f < files.size(); f++){

    TFile* file = TFile::Open( path+files[f], "UPDATE" );
    
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

	cout << "Histo " << hname << endl;		
	//h0 = hist;

	for( unsigned int s = 0; s < systs.size(); s++  ){
	  TH1F* hUp0   = (TH1F*)gDirectory->Get( (hname+"_"+systs[s]+"Up"));
	  TH1F* hDown0 = (TH1F*)gDirectory->Get( (hname+"_"+systs[s]+"Down"));
	  if( hUp0==0 || hDown0==0 ){
	    cout << "Cannot find" << endl;
	    continue;
	  }
	  double ratioUp   = hist->Integral()>0 ? hUp0->Integral()/hist->Integral() : 1.;
	  double ratioDown = hist->Integral()>0 ? hDown0->Integral()/hist->Integral() : 1.;

	  h0    = new TH1F((hname+"_"+systs[s]+"_copy"),     "", nBins, (*(hist->GetXaxis()->GetXbins())).GetArray() );
	  hUp   = new TH1F((hname+"_"+systs[s]+"Up_copy"),   "", nBins, (*(hist->GetXaxis()->GetXbins())).GetArray() );
	  hDown = new TH1F((hname+"_"+systs[s]+"Down_copy"), "", nBins, (*(hist->GetXaxis()->GetXbins())).GetArray() );

	  hUp->SetLineColor(kRed); hUp->SetLineWidth(2);
	  hDown->SetLineColor(kBlue); hDown->SetLineWidth(2);

	  for(int b = 1 ; b <=nBins ; b++){

	    float p_b = hist->GetBinContent(b)/integ;
	    int n_b = int(p_b*nmax);

	    n_b = int(nmax/nBins);

	    cout << "Number of iterations for bin b=" << n_b << endl;
	    float xLow  = hist->GetXaxis()->GetBinLowEdge(b);
	    float xHigh = hist->GetXaxis()->GetBinUpEdge(b);
	    if(b>(nBins/2)){
	      xLow -=1.;
	      xHigh-=1.;
	    }
	    cout << "[" << xLow << "," << xHigh << "]" << endl;

	    int it=0;
	    while(it<n_b){	      
	      double ran_w = ran->Uniform(xLow, xHigh);
	      double ran_o = ran_w;//1./forms[d]*(1./ran_w-1);

	      if(b==1) cout << hist->FindBin(ran_w) << ": " ;

	      if(b<=(nBins/2))
		h0->Fill( ran_w );
	      else
		h0->Fill( ran_w+1. );

	      double logo = ran_o;//TMath::Finite(ran_o) ?  TMath::Log(ran_o) : +998. ;
	      //if( logo<-999.) logo = -998.;
	      
	      int whichbin = -99;
	      for(int l = 0 ; l < bins ; l++){
		if( logo>= limits[l] && logo<limits[l+1]) whichbin = l;
	      }

	      //whichbin = 0;
	      
	      // doing for up
	      TH2F* pdf = (TH2F*)f2D->Get(dirspdf[d]+"_2D_"+TString(Form("%d",whichbin)));
	      if( pdf ==0 ){
		cout << "Cannot find histo up" << string((dirspdf[d]+"_2D_"+TString(Form("%d",whichbin))).Data()) << endl;
		continue;
	      }	      
	      
	      pdf->Reset();
	      for(int k=0;k<1000;k++){
		pdf->Fill( ran->Gaus(1,0.2), ran->Gaus(1,0.2));
	      }

	      double smear_o_up ;  
	      double smear_o_down ; 
	      pdf->GetRandom2( smear_o_up, smear_o_down);
	      smear_o_up   *= ran_o;
	      smear_o_down *= ran_o;

	      smear_o_up   = TMath::Min( smear_o_up, 1.0);
	      smear_o_down = TMath::Min( smear_o_down, 1.0);

	      smear_o_up   = TMath::Max( smear_o_up, 0.);
	      smear_o_down = TMath::Max( smear_o_down, 0.);

	      double smear_w_up = smear_o_up;//1./(1+forms[d]*smear_o_up);

	      if(b==1) cout << hist->FindBin(smear_w_up) << ", ";
	      //smear_w_up = ran_w;

	      if(b<=(nBins/2))
		hUp->Fill( smear_w_up );
	      else
		hUp->Fill( smear_w_up+1. );

	      // doing for down
	      //TH1F* pdf_down = (TH1F*)f2D->Get(dirspdf[d]+"_Down_"+TString(Form("%d",whichbin)));
	      //if( pdf_down ==0 ){
	      //cout << "Cannot find histo down" << string((dirspdf[d]+"_Down_"+TString(Form("%d",whichbin))).Data()) << endl;
	      //continue;
	      //}	      
	      //double smear_o_down = pdf_down->GetRandom()*ran_o;
	      double smear_w_down = smear_o_down;//1./(1+forms[d]*smear_o_down);

	      //smear_w_down = ran_w;
	      if(b==1) cout << hist->FindBin(smear_w_down) << endl;

	      if(b<=(nBins/2))
		hDown->Fill( smear_w_down );
	      else
		hDown->Fill( smear_w_down+1. );


	      it++;
	    }//itaeartions

	  }// bins
	  

	  //for(int bb=1; bb<=hist->GetNbinsX();bb++){
	  //hUp  ->SetBinContent( bb, hist->GetBinContent(bb)*ratioUp);
	  //hDown->SetBinContent( bb, hist->GetBinContent(bb)*ratioDown);
	  //}
	}
      } 

    } // directories

 }// files


 h0   ->DrawNormalized();
 hUp  ->DrawNormalized("SAME");
 hDown->DrawNormalized("SAME");

 cout << h0->GetEntries() << endl;
 cout << hUp->GetEntries() << endl;
 cout << hDown->GetEntries() << endl;


}

*/
