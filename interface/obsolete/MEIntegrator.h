#ifndef MEINTEGRATOR_H
#define MEINTEGRATOR_H

#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include "TH1F.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TF1.h"
#include "TF2.h"
#include "TLegend.h"
#include "TList.h"
#include "THStack.h"
#include "TCut.h"
#include "TArrayF.h"
#include "TObjArray.h"
#include "TVector3.h"
#include "TStyle.h"
#include "TGraph.h"
#include "TKey.h"
#include "TMultiGraph.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "Math/GenVector/LorentzVector.h"
#include "TLorentzVector.h"
#include "RooWorkspace.h"
#include "TH1F.h"

#include <string>
#include <map>
#include <limits>  

#define DEBUG 0

using namespace RooFit;
using namespace std;

#define MTOP 175.
#define MW    81.
#define Mb    5.


class MEIntegrator {

 public:

  MEIntegrator( string , int );
  ~MEIntegrator();

  double Eval(const double* ) const;  
  void   setInputLV( TLorentzVector , TLorentzVector );
  void   SetPar(int);
  void   setJets( vector<TLorentzVector>* );
  void   setBtag( std::vector<float>* );
  void   createMash();
  double        probability(const double*) const;
  unsigned int  findMatch(double, double) const;
  void   saveJetParam( string );
  void   cachePdf( string , string , int );
  TH1F*  getCachedPdf( string ) const;
  TH2F*  getMash( );
  TH1F*  getDebugHisto( );
  void   debug();


 private:
  
  RooWorkspace *w_;
  std::map<string, double> jetParam_; 
  std::map<string, TH1F*> variables_;
  vector<TLorentzVector> jets_;
  vector<float> bTagging_;
  TLorentzVector higgsLV_;
  TLorentzVector topLepLV_;
  int par_;
  float pStar_;
  float EbStar_;
  float EWStar_;
  float EuStar_;
  TH2F* mash_;

  TH1F* debugHisto1_;

};


MEIntegrator::MEIntegrator( string fileName , int param  ) {

  cout << "Begin constructor" << endl;

  par_      = param;
  higgsLV_. SetPxPyPzE(0.,0.,0.,0.);
  topLepLV_.SetPxPyPzE(0.,0.,0.,0.);
  jets_.    clear();
  bTagging_.clear();

  pStar_  = TMath::Sqrt( (MTOP*MTOP-(MW+Mb)*(MW+Mb) )*( MTOP*MTOP-(MW-Mb)*(MW-Mb) ) )/2./MTOP;
  EbStar_ =  (MTOP*MTOP - MW*MW + Mb*Mb)/2./MTOP;
  EWStar_ =  (MTOP*MTOP + MW*MW - Mb*Mb)/2./MTOP;
  EuStar_ =  MW/2;

  mash_        = new TH2F("mash","",500,-2.5,2.5, 628, -TMath::Pi(), TMath::Pi());
  debugHisto1_ = new TH1F("debugHisto1","w1 pt", 100,0,400);

  TFile* file = TFile::Open(fileName.c_str(),"READ");
  w_ = (RooWorkspace*)file->Get("transferFuntions");

  // jets
  RooArgSet allVars = w_->allVars();
  TIterator* iter = allVars.createIterator();
  RooRealVar* var = 0;
  while( (var = (RooRealVar*)(*iter)() ) ){
    jetParam_[ string(var->GetName()) ] = var->getVal();
  }

  cachePdf( "pdfCsvLight",  "csvReco", 100);
  cachePdf( "pdfCsvHeavy",  "csvReco", 100);
  cachePdf( "pdfV1",             "X1", 100);
  cachePdf( "pdfV2",             "X2", 100);
  cachePdf( "pdfPtTTH",       "PtTTH", 100);
  cachePdf( "pdfPtHStar",   "PtHStar", 100);
  cachePdf( "pdfBetaW",       "BetaW", 100);
  cachePdf( "pdfGammaW",     "GammaW", 100);


  cout << "End constructor" << endl;
  //file->Close();

}


MEIntegrator::~MEIntegrator(){

  for(std::map<string, TH1F*>::iterator it = variables_.begin(); it!=variables_.end(); it++){
    if(it->second) 
      delete (it->second);
  }

  delete mash_; delete debugHisto1_;

}


double MEIntegrator::Eval(const double* x) const {
  double prob = probability(x);      
  if ( TMath::IsNaN(prob) ) prob = 0.;
  return prob;
}

void MEIntegrator::setInputLV( TLorentzVector lv1, TLorentzVector lv2){ 
  higgsLV_  = lv1;
  topLepLV_ = lv2;
}

void MEIntegrator::SetPar(int p){ 
  par_ = p; 
}
 

void MEIntegrator::setJets( std::vector<TLorentzVector>* jets){
  jets_.clear();
  for(unsigned int k = 0 ; k<jets->size() ; k++)
    jets_.push_back( (*jets)[k] );
}

void MEIntegrator::setBtag( std::vector<float>* bTagging ){
  bTagging_.clear();
  for(unsigned int k = 0 ; k<bTagging->size() ; k++)
    bTagging_.push_back( (*bTagging)[k] );
}


void   MEIntegrator::createMash(){

  mash_->Reset();

  for(unsigned int j = 0; j < jets_.size() ; j++){

    float eta = jets_[j].Eta();
    float phi = jets_[j].Phi();

    for(int binX = 1; binX<mash_->GetNbinsX(); binX++ ){
      float binCenterX = mash_->GetXaxis()->GetBinCenter(binX);
      for(int binY = 1; binY<mash_->GetNbinsY(); binY++ ){
	float binCenterY = mash_->GetYaxis()->GetBinCenter(binY);

	float dist = TMath::Sqrt( (eta-binCenterX)*(eta-binCenterX) + (phi-binCenterY)*(phi-binCenterY) );
	if( j<4 ){ // jets to veto
	  if( dist<0.50 || mash_->GetBinContent(binX,binY)==-1 ) 
	    mash_->SetBinContent(binX,binY,  -1);
	  else
	    mash_->SetBinContent(binX,binY,   0);
	}
	else{ // jets candidate W had
	  if( mash_->GetBinContent(binX,binY)==0 && dist<0.30 ) 
	    mash_->SetBinContent(binX,binY,   j);
	}

      }
    }

  }


}


TH2F*  MEIntegrator::getMash( ){
  return mash_;
}

TH1F* MEIntegrator::getDebugHisto(){
  return debugHisto1_;
}

void MEIntegrator::saveJetParam( string paramName){

  RooRealVar* var = w_->var(paramName.c_str());
  jetParam_[ paramName ] = var->getVal();

}

void MEIntegrator::cachePdf( string pdfName, string varName, int nBins){
  
  RooAbsPdf* pdf = w_->pdf( pdfName.c_str() );
  if(pdf){
    RooRealVar* var = w_->var( varName.c_str() );

    TH1F* hCache = new TH1F("hCache","", nBins, var->getMin(), var->getMax() );

    for(int j = 1; j < hCache->GetNbinsX(); j++){
      float lvalue = ( TMath::Abs(var->getMax() - var->getMin() )/nBins)*(j-0.5) +  var->getMin();
      var->setVal( lvalue );
      hCache->SetBinContent(j, pdf->getVal( RooArgSet(*var) ) );
    }

    variables_[ pdfName ] = hCache;

  }

}


TH1F* MEIntegrator::getCachedPdf( string pdfName ) const{

  return (variables_.find(pdfName))->second;

}

 



double MEIntegrator::probability(const double* x) const{

  double prob = 1.0;
  
  double v1       = x[0];
  double v2       = x[1];
  double PtTTH    = x[2];
  double PhiTTH   = x[3];
  double BetaW    = x[4];
  double DeltaW   = x[5];
  double GammaW   = x[6];
  double EpsilonW = x[7];

  if(DEBUG){
    cout << "v1    = " << v1 << endl; 
    cout << "v2    = " << v2 << endl; 
    cout << "PtTTH = " << PtTTH << endl; 
    cout << "PhiTTH = " <<PhiTTH  << endl; 
    cout << "BetaW = " << BetaW  << endl; 
    cout << "DeltaW = " << DeltaW  << endl; 
    cout << "GammaW = " <<  GammaW << endl; 
    cout << "EpsilonW = " <<  EpsilonW << endl; 
  }

  TLorentzVector higgsLV  = higgsLV_;
  TLorentzVector topLepLV = topLepLV_;
  TLorentzVector PTOT;
  PTOT.SetPxPyPzE( PtTTH*TMath::Cos( PhiTTH ), PtTTH*TMath::Sin( PhiTTH ), 8000*v2, 
		   TMath::Sqrt( PtTTH*PtTTH + 8000*v2*8000*v2  + 8000*8000*v1*v1 )  );

  TVector3 boostToCMS(PTOT.Px()/PTOT.E(), PTOT.Py()/PTOT.E(), PTOT.Pz()/PTOT.E() );

  TLorentzVector topHadLV = PTOT - higgsLV - topLepLV ;
  if( topHadLV.E()<MTOP ){
    prob =   std::numeric_limits<double>::min();
    return prob; 
  }

  topHadLV.SetE( TMath::Sqrt( (PTOT - higgsLV - topLepLV).P()*(PTOT - higgsLV - topLepLV).P() + MTOP*MTOP ) );
  if(DEBUG){
    cout << "Top had in lab " << topHadLV.Pt() << "," << topHadLV.Eta() << ", " << topHadLV.Phi() << endl;
  }

  higgsLV. Boost( -boostToCMS );
  topLepLV.Boost( -boostToCMS );

  double PtHStar   = higgsLV.P()/PTOT.M();
  double DphiHStar = TMath::Cos((higgsLV.Vect()).Angle(topLepLV.Vect()));


  //////////////////////// 
  TVector3 boostToTopHadCMS(topHadLV.Px()/topHadLV.E(), topHadLV.Py()/topHadLV.E(), topHadLV.Pz()/topHadLV.E() );
  TLorentzVector topHadLVInCMS = topHadLV;
  topHadLVInCMS.Boost( -boostToTopHadCMS );

  TVector3 versor1 = boostToTopHadCMS.Unit();
  TVector3 versor2 = (versor1.Orthogonal()).Unit();
  TVector3 versor3 = versor1.Cross(versor2);

  TVector3 dirB = BetaW*versor1 + TMath::Sqrt(1-BetaW*BetaW)*TMath::Cos( DeltaW )*versor2 + TMath::Sqrt(1-BetaW*BetaW)*TMath::Sin( DeltaW )*versor3;
  //TLorentzVector bLV, wLV;
  //bLV.SetPxPyPzE(  dirB.X()*pStar_,  dirB.Y()*pStar_,  dirB.Z()*pStar_ , EbStar_ );
  //wLV.SetPxPyPzE( -dirB.X()*pStar_, -dirB.Y()*pStar_, -dirB.Z()*pStar_ , EWStar_ );
  TLorentzVector bLV( (dirB*pStar_),    EbStar_ );
  TLorentzVector wLV( (dirB*(-pStar_)), EWStar_ );

  TLorentzVector bInLab = bLV;
  bInLab.Boost( boostToTopHadCMS );

  if(DEBUG){
    cout << "B had in CMS " << bLV.Pt() << "," << bLV.Eta() << ", " << bLV.Phi() << endl;
    cout << "B had in Lab " << bInLab.Pt() << "," << bInLab.Eta() << ", " << bInLab.Phi() << endl;
  }

  
  TVector3 boostToWHadCMS(wLV.Px()/wLV.E(), wLV.Py()/wLV.E(), wLV.Pz()/wLV.E() );
  TLorentzVector WHadLVInCMS = wLV;
  WHadLVInCMS.Boost( -boostToWHadCMS );

  versor1 = boostToWHadCMS.Unit();
  versor2 = (versor1.Orthogonal()).Unit();
  versor3 = versor1.Cross(versor2);

  TVector3 dirW = GammaW*versor1 + TMath::Sqrt(1-GammaW*GammaW)*TMath::Cos( EpsilonW )*versor2 + TMath::Sqrt(1-GammaW*GammaW)*TMath::Sin( EpsilonW )*versor3;
  //TLorentzVector w1LV, w2LV;
  //w1LV.SetPxPyPzE(  dirW.X()*EuStar_,  dirW.Y()*EuStar_,  dirW.Z()*EuStar_ , EuStar_ );
  //w2LV.SetPxPyPzE( -dirW.X()*EuStar_, -dirW.Y()*EuStar_, -dirW.Z()*EuStar_ , EuStar_ );
  TLorentzVector w1LV( (dirW*EuStar_),    EuStar_);
  TLorentzVector w2LV( (dirW*(-EuStar_)), EuStar_);

  TLorentzVector w1InLab = w1LV;
  TLorentzVector w2InLab = w2LV;
  w1InLab.Boost( boostToWHadCMS );
  w2InLab.Boost( boostToWHadCMS );
  w1InLab.Boost( boostToTopHadCMS );
  w2InLab.Boost( boostToTopHadCMS );
  //////////////////////// 

  if(DEBUG){
    //debugHisto1_->Fill( w1InLab.Pt() );
    cout << "W1 in Lab " << w1InLab.Pt() << ", " << w1InLab.Eta() << ", " << w1InLab.Phi() << endl;
    cout << "W2 in Lab " << w2InLab.Pt() << ", " << w2InLab.Eta() << ", " << w2InLab.Phi() << endl;

    cout << "Check Pt: " <<   "Top had in lab " << topHadLV.Pt() << " from components: " << (w1InLab+w2InLab+bInLab).Pt()   << endl;
    cout << "Check Eta: " <<  "Top had in lab " << topHadLV.Eta() << " from components: " << (w1InLab+w2InLab+bInLab).Eta() << endl;
    cout << "Check BetaW..." << endl;
    TLorentzVector bHelp = bInLab;
    bHelp.Boost( -boostToTopHadCMS );
    cout << BetaW << " <=> " << TMath::Cos( (bHelp.Vect()).Angle(boostToTopHadCMS) ) << endl;
    cout << "Check GammaW..." << endl;
    TLorentzVector w1Help = w1InLab;
    TLorentzVector wHelp  = w1InLab+w2InLab;
    wHelp.Boost( -boostToTopHadCMS );
    w1Help.Boost( -boostToTopHadCMS );
    w1Help.Boost( -wHelp.BoostVector() );
    cout << GammaW << " <=> " << TMath::Cos( (w1Help.Vect()).Angle(wHelp.BoostVector()) ) << endl;
  }


  unsigned int j1 = findMatch(w1InLab.Eta(), w1InLab.Phi());
  unsigned int j2 = findMatch(w2InLab.Eta(), w2InLab.Phi());
  unsigned int bb = findMatch(bInLab.Eta(),  bInLab.Phi());


  double AccJ1 =  TMath::Abs(w1InLab.Eta())<1.0 ?  
    TMath::Erf( (jetParam_.find("param0AccLightBin0"))->second * w1InLab.Pt()+(jetParam_.find("param1AccLightBin0"))->second )+(jetParam_.find("param2AccLightBin0"))->second :
    TMath::Erf( (jetParam_.find("param0AccLightBin1"))->second * w1InLab.Pt()+(jetParam_.find("param1AccLightBin1"))->second )+(jetParam_.find("param2AccLightBin1"))->second ;
  double AccJ2 =  TMath::Abs(w2InLab.Eta())<1.0 ?  
    TMath::Erf( (jetParam_.find("param0AccLightBin0"))->second * w2InLab.Pt()+(jetParam_.find("param1AccLightBin0"))->second )+(jetParam_.find("param2AccLightBin0"))->second :
    TMath::Erf( (jetParam_.find("param0AccLightBin1"))->second * w2InLab.Pt()+(jetParam_.find("param1AccLightBin1"))->second )+(jetParam_.find("param2AccLightBin1"))->second ;
  double AccB =  TMath::Abs(bInLab.Eta())<1.0 ?  
    TMath::Erf( (jetParam_.find("param0AccHeavyBin0"))->second * bInLab.Pt()+(jetParam_.find("param1AccHeavyBin0"))->second )+(jetParam_.find("param2AccHeavyBin0"))->second :
    TMath::Erf( (jetParam_.find("param0AccHeavyBin1"))->second * bInLab.Pt()+(jetParam_.find("param1AccHeavyBin1"))->second )+(jetParam_.find("param2AccHeavyBin1"))->second ;
   
  double ResolJ1 =  999.;
  if( w1InLab.Pt()>0){
    ResolJ1 = TMath::Abs(w1InLab.Eta())<1.0 ?  
      w1InLab.Pt()*TMath::Sqrt( (jetParam_.find("param0resolLightBin0"))->second*(jetParam_.find("param0resolLightBin0"))->second/w1InLab.Pt() +  
				(jetParam_.find("param1resolLightBin0"))->second*(jetParam_.find("param1resolLightBin0"))->second/w1InLab.Pt()/w1InLab.Pt()) :
      w1InLab.Pt()*TMath::Sqrt( (jetParam_.find("param0resolLightBin1"))->second*(jetParam_.find("param0resolLightBin1"))->second/w1InLab.Pt() +  
				(jetParam_.find("param1resolLightBin1"))->second*(jetParam_.find("param1resolLightBin1"))->second/w1InLab.Pt()/w1InLab.Pt()) ;
  }
  double ResolJ2 =  999.;
  if( w2InLab.Pt()>0){
    ResolJ2 = TMath::Abs(w2InLab.Eta())<1.0 ?  
      w2InLab.Pt()*TMath::Sqrt( (jetParam_.find("param0resolLightBin0"))->second*(jetParam_.find("param0resolLightBin0"))->second/w2InLab.Pt() +  
				(jetParam_.find("param1resolLightBin0"))->second*(jetParam_.find("param1resolLightBin0"))->second/w2InLab.Pt()/w2InLab.Pt()) :
      w2InLab.Pt()*TMath::Sqrt( (jetParam_.find("param0resolLightBin1"))->second*(jetParam_.find("param0resolLightBin1"))->second/w2InLab.Pt() +  
				(jetParam_.find("param1resolLightBin1"))->second*(jetParam_.find("param1resolLightBin1"))->second/w2InLab.Pt()/w2InLab.Pt()) ;
  }
  double ResolB =  999.;
  if( bInLab.Pt()>0){
    ResolB = TMath::Abs(bInLab.Eta())<1.0 ?  
      bInLab.Pt()*TMath::Sqrt( (jetParam_.find("param0resolHeavyBin0"))->second*(jetParam_.find("param0resolHeavyBin0"))->second/bInLab.Pt() +  
			       (jetParam_.find("param1resolHeavyBin0"))->second*(jetParam_.find("param1resolHeavyBin0"))->second/bInLab.Pt()/bInLab.Pt()) :
      bInLab.Pt()*TMath::Sqrt( (jetParam_.find("param0resolHeavyBin1"))->second*(jetParam_.find("param0resolHeavyBin1"))->second/bInLab.Pt() +  
			       (jetParam_.find("param1resolHeavyBin1"))->second*(jetParam_.find("param1resolHeavyBin1"))->second/bInLab.Pt()/bInLab.Pt()) ;
  }

  if(DEBUG){
    cout << "AccJ1 " << AccJ1 << endl;
    cout << "AccJ2 " << AccJ2 << endl;
    cout << "AccB "  << AccB << endl;

    cout << "ResolJ1 " << ResolJ1 << endl;
    cout << "ResolJ2 " << ResolJ2 << endl;
    cout << "ResolB "  << ResolB << endl;

  }

  TH1F* bTagLight = this->getCachedPdf(string("pdfCsvLight"));
  TH1F* bTagHeavy = this->getCachedPdf(string("pdfCsvHeavy"));

  double TFJ1 = 1.0;
  if( j1==98) 
    TFJ1 = 0.0;
  else if(j1==99)
    TFJ1 = TMath::Max(1-AccJ1, 0.0);
  else
    TFJ1 = TMath::Min( AccJ1,  1.0)*TMath::Gaus( w1InLab.Pt(), ResolJ1 );
  TFJ1 *= bTagLight->Interpolate(  bTagging_[j1] );

 double TFJ2 = 1.0;
  if( j2==98) 
    TFJ2 = 0.0;
  else if(j2==99)
    TFJ2 = TMath::Max(1-AccJ2, 0.0);
  else
    TFJ2 = TMath::Min( AccJ2,  1.0)*TMath::Gaus( w2InLab.Pt(), ResolJ2 );
  TFJ2 *= bTagLight->Interpolate(  bTagging_[j2] );

  double TFB = 1.0;
  if( bb==98) 
    TFB = 0.0;
  else if(bb==99)
    TFB = TMath::Max(1-AccB, 0.0);
  else
    TFB = TMath::Min( AccB,  1.0)*TMath::Gaus( bInLab.Pt(), ResolB );
  TFB *= bTagHeavy->Interpolate(  bTagging_[bb] );

  prob *= TFJ1;
  prob *= TFJ2;
  prob *= TFB;

  TH1F* pdfV1      = this->getCachedPdf(string("pdfV1"));
  TH1F* pdfV2      = this->getCachedPdf(string("pdfV2"));
  TH1F* pdfPtHStar = this->getCachedPdf(string("pdfPtHStar"));
  TH1F* pdfPtTTH   = this->getCachedPdf(string("pdfPtTTH"));
  TH1F* pdfBetaW   = this->getCachedPdf(string("pdfBetaW"));
  TH1F* pdfGammaW  = this->getCachedPdf(string("pdfGammaW"));


  double pdf = 
    pdfV1->Interpolate( v1 ) *
    pdfV2->Interpolate( v2 ) *
    pdfPtTTH->Interpolate( PtTTH ) *
    (1./2./TMath::Pi())    ;

  double cms = 
    pdfPtHStar->Interpolate( PtHStar ) ;

  double decay = 
    pdfBetaW->Interpolate(  BetaW  ) * 
    (1./2./TMath::Pi())  * 
    pdfGammaW->Interpolate( GammaW ) * 
    (1./2./TMath::Pi())  ;

  prob *= pdf;
  prob *= cms;
  prob *= decay;

 if(DEBUG){
   cout << "TFJ1 " << TFJ1 << endl;
   cout << "TFJ2 " << TFJ2 << endl;
   cout << "TFB  " << TFB << endl;
   debugHisto1_->Fill( w1InLab.Pt() , decay*pdf*cms);
 }

 prob = pdf*cms*decay;

  ////////////////////////
  return prob;
}


unsigned int MEIntegrator::findMatch(double eta, double phi) const{

  unsigned int match = 99;
  int bin   = mash_->FindBin(eta, phi);
  float res = mash_->GetBinContent(bin);
  if( res>-0.5 && res<0.5 )
    match = 99; // no match
  else if(  res<-0.5 )
    match = 98; // veto
  else
    match = res;
 

  return match;
}




void MEIntegrator::debug(){

  cout << "*** debug start *** " << endl;

  for(std::map<string, double>::iterator it = jetParam_.begin(); it!=jetParam_.end(); it++){
    cout << it->first << " => " << it->second << endl;
  }
  for(std::map<string, TH1F*>::iterator it = variables_.begin(); it!=variables_.end(); it++){
    if(it->second) 
      cout << it->first << " => " << (it->second)->Integral() << endl;
    else
      cout << it->first << " => null pointer" << endl;
  }

  cout << "*** debug end *** " << endl;

}


#endif
