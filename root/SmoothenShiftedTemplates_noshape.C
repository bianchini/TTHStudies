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


void merge(){

  //TString  path = "datacards/PostApproval/2D_noshape/";
  //TString  path = "datacards/PostApproval/MyTest_new_merged/";
  TString  path = "datacards/PostApproval/2D_OptFactors/";

  vector<TString> files;
  //files.push_back("MEM_New_rec_std_2D.root");
  files.push_back("MEM_New_rec_std_sb.root");

  vector<TString> dirs;
  dirs.push_back("MEM_cat1_H");
  dirs.push_back("MEM_cat2_H");
  dirs.push_back("MEM_cat3_H");
  dirs.push_back("MEM_cat6_H");
  dirs.push_back("MEM_cat1_L");
  dirs.push_back("MEM_cat2_L");
  dirs.push_back("MEM_cat3_L");
  dirs.push_back("MEM_cat6_L");

  vector<TString> systs;
  systs.push_back("CMS_scale_j");
  systs.push_back("CMS_res_j");

  vector<string> proc;
  proc.push_back("ttbarV");
  proc.push_back("ttH_hbb");
  //proc.push_back("singlet");
  proc.push_back("ttbar");
  proc.push_back("ttbarPlusB");
  proc.push_back("ttbarPlusBBbar");
  proc.push_back("ttbarPlusCCbar");

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

	if( hname.find("Up")!=string::npos || hname.find("Down")!=string::npos || hname.find("data_obs")!=string::npos)
	  continue;

	int isAmong = 0;
	for( unsigned int p = 0; p < proc.size(); p++){
	  if( hname == proc[p] ) isAmong++;
	}
	if(isAmong<1) continue;

	cout << "Histo " << hname ;	

	for( unsigned int s = 0; s < systs.size(); s++  ){
	  TH1F* hUp   = (TH1F*)gDirectory->Get( (hname+"_"+systs[s]+"Up"));
	  TH1F* hDown = (TH1F*)gDirectory->Get( (hname+"_"+systs[s]+"Down"));
	  if( hUp==0 || hDown==0 ){
	    cout << "Cannot find" << endl;
	    continue;
	  }
	  double ratioUp   = hist->Integral()>0 ? hUp->Integral()/hist->Integral() : 1.;
	  double ratioDown = hist->Integral()>0 ? hDown->Integral()/hist->Integral() : 1.;

	  for(int bb=1; bb<=hist->GetNbinsX();bb++){
	    hUp  ->SetBinContent( bb, hist->GetBinContent(bb)*ratioUp);
	    hDown->SetBinContent( bb, hist->GetBinContent(bb)*ratioDown);
	  }

	  hUp->Write(   hUp->GetName() ,   TObject::kOverwrite);
	  hDown->Write( hDown->GetName() , TObject::kOverwrite);
	  cout << "...saved"  << endl;
	}
      }           

    }// dirs

    file->Close();
  }
  

}
