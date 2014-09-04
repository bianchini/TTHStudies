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


void merge( int binL = 1, int binH = 2){


  string binNameL( Form("bin%d", binL) );
  string binNameH( Form("bin%d", binH) );

  TString  path = "datacards/PostApproval/MyTest_new_merged//";

  vector<TString> files;
  files.push_back("MEM_New_rec_std_sb.root");

  vector<TString> dirs;
  dirs.push_back("MEM_cat1_H");
  dirs.push_back("MEM_cat2_H");
  //dirs.push_back("MEM_cat3_H");
  dirs.push_back("MEM_cat6_H");
  dirs.push_back("MEM_cat1_L");
  dirs.push_back("MEM_cat2_L");
  //dirs.push_back("MEM_cat3_L");
  dirs.push_back("MEM_cat6_L");


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

	cout << "List: " << (gDirectory->GetListOfKeys())->GetSize() << endl;

	TH1F* hist   =  (TH1F*)key->ReadObj(); 
	if(hist==0){
	  cout << "Null" << endl;
	  continue;
	}	
	string hname = string(hist->GetName());
	if(hname == "h2"){
	  cout << "Found h2.. continue" << endl;
	  continue;
	}

	cout << "Doing..." << hname ;

	// number of bins
	int   nBins = hist->GetNbinsX();

	// merge first two bins
	TArrayD newBins(nBins);
	const TArrayD* oldBins = hist->GetXaxis()->GetXbins();
	for( int k = 0; k <= nBins ; k++ ){
	  if(k<=(binL-1)) newBins[k]   = (*oldBins)[k];
	  if(k==binL)     continue;
	  if(k>=binH)     newBins[k-1] = (*oldBins)[k];
	}


	TH1F* h2 = new TH1F("h2", "", nBins-1, newBins.GetArray());
	//cout << "Bins: " << endl;
	for(int w = 1 ; w <= nBins-1; w++){
	  //cout << " [" << w << "] --> [" << h2->GetXaxis()->GetBinLowEdge(w) << "," <<  h2->GetXaxis()->GetBinUpEdge(w) << "]" << endl;
	}
	//continue;


	bool isBBB = 
	  (hname.find("Hbin")!=string::npos     || hname.find("Lbin")!=string::npos) &&
	  (hname.find( binNameL )!=string::npos || hname.find( binNameH )!=string::npos) ;

	if(!isBBB){

	  for(int b = 1; b <= nBins; b++){      
	    
	    if( b == binL || b == binH){
	      float err1 = hist->GetBinError(binL);
	      float err2 = hist->GetBinError(binH);
	      float bin1 = hist->GetBinContent(binL);
	      float bin2 = hist->GetBinContent(binH);
	      h2->SetBinContent(b, bin1+bin2 );
	      h2->SetBinError  (b, sqrt(err1*err1 + err2*err2) );	    
	    }
	    else{
	      float err1 = hist->GetBinError(b);
	      float bin1 = hist->GetBinContent(b);
	      
	      int binToFill = b;
	      if( b<binL ) binToFill = b;
	      if( b>binH ) binToFill = b-1;
	      h2->SetBinContent(binToFill, bin1 );
	      h2->SetBinError  (binToFill, err1 );
	    }	    
	  }

	}
	else{

	  for(int b = 1; b <= nBins; b++){      
	    
	    if( b == binL || b == binH){
	      float err1 = hist->GetBinError(binL);
	      float err2 = hist->GetBinError(binH);
	      float bin1 = hist->GetBinContent(binL);
	      float bin2 = hist->GetBinContent(binH);

	      if(  hname.find("Up")  !=string::npos ){
		if( hname.find( binNameL )!=string::npos ) bin1 -= err1;
		if( hname.find( binNameH )!=string::npos ) bin2 -= err2;		
		h2->SetBinContent(b, bin1+bin2 + sqrt(err1*err1 + err2*err2) );
	      }
	      if(  hname.find("Down")!=string::npos ){
		if( hname.find( binNameL )!=string::npos ) bin1 += err1;
		if( hname.find( binNameH )!=string::npos ) bin2 += err2;
		h2->SetBinContent(b, bin1+bin2 - sqrt(err1*err1 + err2*err2) );
	      }

	    }
	    else{
	      float err1 = hist->GetBinError(b);
	      float bin1 = hist->GetBinContent(b);
	      
	      int binToFill = b;
	      if( b<binL ) binToFill = b;
	      if( b>binH ) binToFill = b-1;
	      h2->SetBinContent(binToFill, bin1 );
	      h2->SetBinError  (binToFill, err1 );
	    }	    
	  }

	}
    
	h2->SetName( hname.c_str() );
	h2->Write( hname.c_str() , TObject::kWriteDelete ); 
	cout << "...saved"  << " --> list: " << (gDirectory->GetListOfKeys())->GetSize() << endl;
	delete h2;
      }           


    }// dirs

    file->Close();
  }
  

}


void mergeAll(){

  merge(11,12);
  merge( 5, 6);

}
