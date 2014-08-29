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

  TString  path = "datacards/Approval/zero_bins/fromNominal/";

  vector<TString> files;
  files.push_back("MEM_New_rec_std_sb.root");

  vector<TString> dirs;
  dirs.push_back("MEM_cat1_H");
  dirs.push_back("MEM_cat2_H");
  dirs.push_back("MEM_cat3_H");
  dirs.push_back("MEM_cat6_H");


  for(unsigned int f = 0; f < files.size(); f++){

    TFile* file = TFile::Open( path+files[f], "UPDATE" );
    
    for( unsigned int d = 0; d < dirs.size(); d++ ){
      file->cd( dirs[d] ); 

      if(gDirectory->FindObject("h2")!=0){
	gDirectory->Remove(gDirectory->FindObject("h2"));
      }     
      
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
	if(hname == "h2") continue;

	// number of bins
	int   nBins = hist->GetNbinsX();

	// merge first two bins
	TH1F* h2 = new TH1F("h2", "", nBins-1, 0, 1);
	
	cout << "Doing..." << hname ;

	bool isBBB = hname.find("Hbin")!=string::npos || hname.find("Lbin")!=string::npos ;

	for(int b = 1; b <= nBins; b++){
      
	  if(b<=2){
	    float err1 = hist->GetBinError(1);
	    float err2 = hist->GetBinError(2);
	    float bin1 = hist->GetBinContent(1);
	    float bin2 = hist->GetBinContent(2);
	    if( !isBBB ){	      
	      h2->SetBinContent(b, bin1+bin2 );
	      h2->SetBinError  (b, sqrt(err1*err1 + err2*err2) );
	    }
	    if( isBBB && (hname.find("bin1")!=string::npos || hname.find("bin2")!=string::npos) ){	      	      
	      if(  hname.find("Up")  !=string::npos ){
		if( hname.find("bin1")!=string::npos ) bin1 -= err1;
		if( hname.find("bin2")!=string::npos ) bin2 -= err2;		
		h2->SetBinContent(b, bin1+bin2 + sqrt(err1*err1 + err2*err2) );
	      }
	      if(  hname.find("Down")!=string::npos ){
		if( hname.find("bin1")!=string::npos ) bin1 += err1;
		if( hname.find("bin2")!=string::npos ) bin2 += err2;
		h2->SetBinContent(b, bin1+bin2 - sqrt(err1*err1 + err2*err2) );
	      }
	    }
	  }
	  else{
	    float err1 = hist->GetBinError(b);
	    float bin1 = hist->GetBinContent(b);
	    h2->SetBinContent(b-1, bin1 );
	    h2->SetBinError  (b-1, err1 );
	  }
	  
	}
    
	h2->Write( hname.c_str() , TObject::kOverwrite); 
	cout << "...saved"  << endl;
	delete h2;
      }           

    }// dirs

    file->Close();
  }
  

}
