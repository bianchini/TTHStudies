#include "Bianchi/TTHStudies/interface/DiJetCandidateFinder.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <algorithm>

using namespace std;


void DiJetCandidateFinder::run(int VType , std::vector<JetByPt>& jetsCount,  std::vector<JetByPt>& jetsRecover, std::vector<LV>& leptons){

  
  switch (VType){
  case 0:
    this->runLL(jetsCount, jetsRecover, leptons);
    break;
  case 1:
    this->runLL(jetsCount, jetsRecover, leptons);
    break;
  case 2:
    this->runLJ(jetsCount, jetsRecover, leptons);
    break;
  case 3:
    this->runLJ(jetsCount, jetsRecover, leptons);
    break;
  default:
    cout << "VType " << VType << " unsupported for DiJetCandidateFinder::run" << endl;
    break;
   }

  return;

}

int DiJetCandidateFinder::jetCountMult(){
  return numJetCount_;
}

int DiJetCandidateFinder::jetRecoverMult(){
  return numJetRecov_;
}

int DiJetCandidateFinder::eventFlag(){
  return eventFlag_;
}

int DiJetCandidateFinder::errorFlag(){
  return error_;
}


void DiJetCandidateFinder::runLL(std::vector<JetByPt>& jetsCount, std::vector<JetByPt>& jetsRecover, std::vector<LV>& leptons){


}


void DiJetCandidateFinder::runLJ(std::vector<JetByPt>& jetsCount, std::vector<JetByPt>& jetsRecover, std::vector<LV>& leptons){

  if(leptons.size()!=2){
    error_ = 1;
    cout << "Two leptons are needed to run DiJetCandidateFinder::runLJ: the charged letpon and the MET... return" << endl;
    return;
  }

  std::vector<Permutation> permutations;
  int myMatch = DiJetCandidateFinder::kNoMatch;
  Permutation output;

  std::pair<int, unsigned int> scanResults;
  std::vector<unsigned int> visited;

  numJetCount_ =   jetsCount.size();
  numJetRecov_ =   jetsRecover.size();
  

  if( numJetCount_ >= 4 ){
    eventFlag_ = DiJetCandidateFinder::kFirstRank;
    if(verbose_)  cout << "The event contains " << numJetCount_ << " jets, proceed with findMatch..." << endl;
    permutations = findMatchLJ(jetsCount, leptons, visited);
    scanResults  = scanLJ( permutations );
    myMatch      = scanResults.first;
  }

  if(myMatch == kTTbarFound){
    if(verbose_)  cout << "A valid match has been found, proceed with sortUnMasked..." << endl;
    Permutation permutation = permutations[scanResults.second];
    permutations.clear();  
    permutations.push_back(permutation);
    sortUnMasked( output, permutations  , jetsCount);
  }
  else{

    if( numJetCount_ < 4 ){
      if(verbose_) cout << "The event contains " << numJetCount_ << " jets, try to increase acceptance with additional " << numJetRecov_ << " jets" << endl;
      eventFlag_ = DiJetCandidateFinder::kThirdRank;
    }
    else{
      if(verbose_) cout << "Try to increase acceptance with additional " << numJetRecov_ << " jets..." << endl;
      eventFlag_ = DiJetCandidateFinder::kSecondRank;
    }

    jetsCount.insert( jetsCount.end(), jetsRecover.begin(), jetsRecover.end() );
    permutations = findMatchLJ(jetsCount, leptons, visited);
    scanResults  = scanLJ( permutations );
    myMatch      = scanResults.first;

    if(myMatch == kTTbarFound){
      if(verbose_)  cout << "A valid match has been found after increasing acceptance, proceed with sortUnMasked..." << endl;
      Permutation permutation = permutations[scanResults.second];
      permutations.clear();  
      permutations.push_back(permutation);
      sortUnMasked( output, permutations  , jetsCount);  
    }
    else{ 
      if(verbose_)  cout << "No valid match has been found after increasing acceptance, proceed with split..." << endl;

      eventFlag_ = DiJetCandidateFinder::kFourthRank;
      splitJets( jetsCount );
      permutations = findMatchLJ(jetsCount, leptons, visited);
      scanResults  = scanLJ( permutations );
      myMatch      = scanResults.first;
    
      if(myMatch == kTTbarFound){
	if(verbose_)  cout << "A valid match has been found after splitting,  proceed with sortUnMasked..." << endl;
	Permutation permutation = permutations[scanResults.second];
	permutations.clear();  
	permutations.push_back(permutation);
	sortUnMasked( output, permutations  , jetsCount);
      }
      else{
	if(verbose_)  cout << "No valid match has been found after splitting, mask matched jets ..." << endl;
	eventFlag_ = DiJetCandidateFinder::kFifthRank;
	sortUnMasked(output, permutations , jetsCount );   // do something here
      }
    }
  }

  return;

}


std::vector<Permutation>  DiJetCandidateFinder::findMatchLJ(std::vector<JetByPt>& jets, 
										  std::vector<LV>& leptons,
										  std::vector<unsigned int>& visited){

  std::vector<Permutation > output;

  unsigned int indices[jets.size()]; // indices sorted as [0,1,2,...]
  for(unsigned int k = 0; k<jets.size() ; k++)
    indices[k]=k; 


  int counter = 0;
  if(verbose_) std::cout << "DiJetCandidateFinder::findMatchLJ... start permutations" << endl;
  do {

    // indices are permuted now...
    counter++;
    if(verbose_) cout << "#" << counter << ":"<< endl;

    Permutation permutation;
    std::map<unsigned int, int> flagMask;
    for(unsigned int k = 0; k<jets.size() ; k++) // set all bits to zero
      flagMask[k] = 0;

    unsigned int indexWdau1Cand = indices[0];
    unsigned int indexWdau2Cand = indices[1];
    unsigned int indexbCand     = indices[2];
    unsigned int indexbbarCand  = indices[3];

      
    if( !bookKeeper( indexWdau1Cand, indexWdau2Cand, indexbCand, indexbbarCand, visited) ) continue;
    
    if(verbose_){
      cout << "W(1): jets[" << indexWdau1Cand << "]" << endl;
      cout << "W(2): jets[" << indexWdau2Cand << "]" << endl;
      cout << "b:    jets[" << indexbCand << "]" << endl;
      cout << "bbar: jets[" << indexbbarCand << "]" << endl;
    }

    LV Wdau1Cand( jets[indexWdau1Cand].pt, jets[indexWdau1Cand].eta, jets[indexWdau1Cand].phi, jets[indexWdau1Cand].mass);
    LV Wdau2Cand( jets[indexWdau2Cand].pt, jets[indexWdau2Cand].eta, jets[indexWdau2Cand].phi, jets[indexWdau2Cand].mass);
    LV bCand(     jets[indexbCand].pt,     jets[indexbCand].eta,     jets[indexbCand].phi,     jets[indexbCand].mass);
    LV bbarCand(  jets[indexbbarCand].pt,  jets[indexbbarCand].eta,  jets[indexbbarCand].phi,  jets[indexbbarCand].mass);
    LV chargedLep = leptons[0];
    LV met        = leptons[1];
    std::vector<LV> neutrinos  = buildNeutrinos(chargedLep, met);

    double likelihood     = 1.0;
    double WHadLikelihood = 1.0;
    double bLikelihood    = 1.0;
    double bbarLikelihood = 1.0;

    if(func1D_.find("csvLF")!=func1D_.end()){
      WHadLikelihood *= func1D_["csvLF"]->Eval( jets[indexWdau1Cand].csv );
      WHadLikelihood *= func1D_["csvLF"]->Eval( jets[indexWdau2Cand].csv );
    }
    if(func1D_.find("csvHF")!=func1D_.end()){
      bLikelihood    *= func1D_["csvHF"]->Eval( jets[indexbCand].csv ); 
      bbarLikelihood *= func1D_["csvHF"]->Eval( jets[indexbbarCand].csv ); 
    }
    if(func1D_.find("wHadMass")!=func1D_.end()){
      WHadLikelihood *= func1D_["wHadMass"]->Eval( (Wdau1Cand+Wdau2Cand).M() ); // W had mass
    }
    if(func1D_.find("topLepMass")!=func1D_.end()){
      double topLepLike = TMath::Max(func1D_["topLepMass"]->Eval( (bbarCand+neutrinos[0]).M() ), // W lep mass
				     func1D_["topLepMass"]->Eval( (bbarCand+neutrinos[1]).M() ));
      bbarLikelihood *= topLepLike;
    }
    if(func2D_.find("topHadMass")!=func2D_.end()){
      bLikelihood  *= func2D_["topWHadMass"]->Eval( (Wdau1Cand+Wdau2Cand).M(), (Wdau1Cand+Wdau2Cand+bCand).M() ); // top vs W mass
    }

    likelihood *= (WHadLikelihood*bLikelihood*bbarLikelihood);

    switch(eventFlag_){ // control to which depth the algo has run

    default: // For the moment, use the same cut values

      bool WHadMatch   = false;
      bool TopHadMatch = false;
      bool TopLepMatch = false;
      bool FullMatch   = false;

      if( WHadLikelihood> 99999 ){  // mask these jets as W candidates
	flagMask[ indexWdau1Cand ] = 1;
	flagMask[ indexWdau2Cand ] = 1;
	WHadMatch = true;
      }
      if( bLikelihood > 99999 ){    // mask these jets as b candidates
	flagMask[ indexbCand ]     = 2;
	TopHadMatch = true;
      }
      if( bbarLikelihood > 99999 ){ // mask these jets as bbar candidates
	flagMask[ indexbbarCand ]  = 3;
	TopLepMatch = true;
      }
      if( likelihood > 9999 ){
	flagMask[ indexWdau1Cand ] = 1;
	flagMask[ indexWdau2Cand ] = 1;
	flagMask[ indexbCand ]     = 2;
	flagMask[ indexbbarCand ]  = 3;
	FullMatch = true;
      }

      int matchType;
      if(      WHadMatch && !TopHadMatch && !TopLepMatch && !FullMatch) matchType = DiJetCandidateFinder::kOnlyWHadFound;
      else if( WHadMatch &&  TopHadMatch && !TopLepMatch && !FullMatch) matchType = DiJetCandidateFinder::kOnlyTFound;
      else if(!WHadMatch && !TopHadMatch &&  TopLepMatch && !FullMatch) matchType = DiJetCandidateFinder::kOnlyTbarFound;
      else if( WHadMatch && !TopHadMatch &&  TopLepMatch && !FullMatch) matchType = DiJetCandidateFinder::kOnlyWHadTbarFound;
      else if( FullMatch )                                              matchType = DiJetCandidateFinder::kTTbarFound;
      else                                                              matchType = DiJetCandidateFinder::kNoMatch;

      permutation.flagMask       = flagMask;
      permutation.likelihood     = likelihood;
      permutation.WHadLikelihood = WHadLikelihood;
      permutation.bLikelihood    = bLikelihood;
      permutation.bbarLikelihood = bbarLikelihood;
      permutation.tagging        = eventFlag_;
      permutation.matchType      = matchType;

      if( WHadMatch || TopHadMatch ||  TopLepMatch || FullMatch){ // prune combinations...
	output.push_back(permutation);
	if(verbose_) cout << " ==> " << string(Form("L=%f : WHad=%f : b=%f : bbar=%f", likelihood, WHadLikelihood, bLikelihood, bbarLikelihood)) 
			  << " ... Accetp!" << endl;
      }

      break;
    }

  } while ( std::next_permutation( indices, indices+jets.size()  ) );



  return output;
}


bool DiJetCandidateFinder::bookKeeper(unsigned int index1, unsigned int index2, 
				      unsigned int index3, unsigned int index4, 
				      vector<unsigned int>& visited){

  unsigned int combination1 = 1000*index1 + 100*index2 + 10*index3 + 1*index3;
  unsigned int combination2 = 1000*index2 + 100*index1 + 10*index3 + 1*index3;
  bool isThere = false;

  for(unsigned int k = 0; k < visited.size() ; k++ )
    if( combination1 == visited[k] || combination2 == visited[k] ) isThere = true;

  if(!isThere){
    visited.push_back( combination1 );
    visited.push_back( combination2 );
    return true;
  }
  else{
    return false;
  }

}

std::pair<int, unsigned int> DiJetCandidateFinder::scanLJ(std::vector<Permutation> permutations){

  std::pair<int, unsigned int> out;
  double bestLikelihood = -999; unsigned int pos = 0;
  for(unsigned int j = 0 ; j < permutations.size() ; j++){
    if( permutations[j].likelihood>bestLikelihood){
      bestLikelihood = permutations[j].likelihood;
      pos = j;
    }
  }

  if( bestLikelihood> 9999){
    out = make_pair( DiJetCandidateFinder::kTTbarFound, pos);
    return out;
  }
  else{
    out = make_pair( DiJetCandidateFinder::kNoMatch, 0);
    return out;
  }

}

void DiJetCandidateFinder::splitJets(std::vector<JetByPt>& jetsCount){
  return;
}


void DiJetCandidateFinder::sortUnMasked(Permutation& out, vector<Permutation> permutations, std::vector<JetByPt>& jets){  

  //Permutation out;
  
  std::pair<unsigned int, unsigned int> diJetIndices;
  
  std::map<unsigned int, int> flagMask;
  for(unsigned int k = 0; k<jets.size() ; k++) // set all bits to zero
    flagMask[k] = 0;

  out.reset();  
  out.flagMask = flagMask;

  if( permutations.size() == 0) return;


  // many combinations... give preference:
  
  bool WHadMatch   = false;     unsigned int posWHadMatch=0;
  bool TopHadMatch = false;     unsigned int posTopHadMatch=0;
  bool TopLepMatch = false;     unsigned int posTopLepMatch=0;
  bool WHadTopLepMatch = false; unsigned int posWHadTopLepMatch=0;
  bool FullMatch   = false;     unsigned int posFullMatch=0;
  
  for(unsigned int k = 0; k<permutations.size() ; k++){
    if( permutations[k].matchType == DiJetCandidateFinder::kTTbarFound){
      posFullMatch = k;
      FullMatch    = true;
    }
    if( permutations[k].matchType == DiJetCandidateFinder::kOnlyTFound){
      posTopHadMatch = k;
      TopHadMatch    = true;
    }
    if( permutations[k].matchType == DiJetCandidateFinder::kOnlyTbarFound){
      posTopLepMatch = k;
      TopLepMatch    = true;
    }
    if( permutations[k].matchType == DiJetCandidateFinder::kOnlyWHadFound){
      posWHadMatch   = k;
      WHadMatch      = true;
    }
    if( permutations[k].matchType == DiJetCandidateFinder::kOnlyWHadTbarFound){
      posWHadTopLepMatch = k;
      WHadTopLepMatch = true;
    }
  }
  
  if( FullMatch ){          // if one combination has full match, use it first
    flagMask = (permutations[posFullMatch]).flagMask;
    out      = permutations[posFullMatch]; 
    }
  else if(WHadTopLepMatch){ // else, use the one with WHad and TopLep only match
    flagMask = (permutations[posWHadTopLepMatch]).flagMask;
    out      = permutations[posWHadTopLepMatch]; 
  }
  else if(TopHadMatch){    // else, use the one with TopHad only match
    flagMask = (permutations[posTopHadMatch]).flagMask;
    out      = permutations[posTopHadMatch]; 
  }
  else if(TopLepMatch){    // else, use the one with TopLep only match
    flagMask = (permutations[posTopLepMatch]).flagMask;
      out      = permutations[posTopLepMatch]; 
  }
  else if(WHadMatch){      // else, use the one with WHad only match
    flagMask = (permutations[posWHadMatch]).flagMask;
    out      = permutations[posWHadMatch]; 
  }else{}
  
  
  // choose pair with highest csv score and Mjj>100; if no pair with Mjj>100, choose highest csv one

  float bestCsvSum    = -999;
  bool atLeastOnePair = false;

  for(unsigned int i = 0; i < jets.size()-1; i++){
    if( (flagMask)[i] != 0 ) continue;
    LV j1( jets[i].pt, jets[i].eta, jets[i].phi, jets[i].mass);
    for(unsigned int j = i+1; j < jets.size(); j++){
      if( (flagMask)[j] != 0 ) continue; 
      LV j2( jets[j].pt, jets[j].eta, jets[j].phi, jets[j].mass);
      float mass = (j1+j2).M();
      if(mass>100. && (jets[i].csv + jets[j].csv)>bestCsvSum ){
	bestCsvSum = (jets[i].csv + jets[j].csv);
	diJetIndices = std::make_pair( i, j);
	atLeastOnePair = true;
      }
    }
  }

  if(!atLeastOnePair){
    for(unsigned int i = 0; i < jets.size()-1; i++){
      if( (flagMask)[i] != 0 ) continue;
      LV j1( jets[i].pt, jets[i].eta, jets[i].phi, jets[i].mass);
      for(unsigned int j = i+1; j < jets.size(); j++){
	if( (flagMask)[j] != 0 ) continue; 
	LV j2( jets[j].pt, jets[j].eta, jets[j].phi, jets[j].mass);
	if( (jets[i].csv + jets[j].csv)>bestCsvSum ){
	  bestCsvSum = (jets[i].csv + jets[j].csv);
	  diJetIndices =  std::make_pair( i, j);
	}
      }
    }
  }

  out.flagMask[ diJetIndices.first ]  = 4; // flag as "Higgs"
  out.flagMask[ diJetIndices.second ] = 4; // flag as "Higgs"

  //out = make_pair( 0, 0);

  return;
}


std::vector<LV> DiJetCandidateFinder::buildNeutrinos(LV lep, LV met){

  vector<LV> output;
  output.push_back( met);

  return output;

}
