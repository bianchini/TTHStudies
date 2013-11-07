#ifndef DIJETCANDIDATEFINDER_H
#define DIJETCANDIDATEFINDER_H

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

#include <string>
#include <utility>
#include <map>
#include <vector>

using namespace std;

typedef  ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LV;
typedef struct
  {
    int index;
    float pt;
    float eta;
    float phi;
    float mass;
    float csv;
    float topB;
    float topW;
    float atopB;
    float atopW;
    float higgsB;
    float flavor;
    float unc;
    void reset(){
      index = -99; pt = -99; eta = -99; phi = -99; mass = -99; csv = -99; topB = -99; topW = -99;  atopW = -99;  atopB = -99;
      higgsB = -99; flavor = -99; unc = -99;
    }
    void set(int index_,  float pt_, float eta_, float phi_, float mass_,float csv_, float topB_, float topW_, float atopB_,
	     float atopW_, float higgsB_, float flavor_, float unc_){
      index = index_; 
      pt = pt_; 
      eta = eta_; 
      phi = phi_; 
      mass = mass_; 
      csv = csv_; 
      topB = topB_; 
      topW = topW_;  
      atopW = atopW_;  
      atopB = atopB_;
      higgsB = higgsB_; 
      flavor = flavor_; 
      unc = unc_;
    }
  } JetByPt;


typedef struct {
  int tagging;
  int matchType;
  std::map<unsigned int, int> flagMask;
  double likelihood;
  double WHadLikelihood;
  double bLikelihood;
  double bbarLikelihood;
  void reset(){
    tagging   = -99;
    matchType  = 0;
    likelihood     = -99;
    WHadLikelihood = -99;
    bLikelihood    = -99;
    bbarLikelihood = -99;
  }
} Permutation;



class DiJetCandidateFinder{

 public:
  
  enum MatchTypeLJ {
    kNoMatch = 0,
    kTTbarFound,
    kOnlyTFound,
    kOnlyTbarFound,
    kOnlyWHadFound,
    kOnlyWHadTbarFound,
    kFound
  };
  
  enum LevelOfTagging {
    kFirstRank = 0, // Nj>=4, TTbar candidates found within "jetsCount"
    kSecondRank,    // Nj>=4, found a TTbar match after increasing acceptance
    kThirdRank,     // Nj<4,  found a TTbar match after increasing acceptance
    kFourthRank,    //        found a TTbar match after aplitting jets
    kFifthRank      //        no TTbar match found, mask matched jets
  };
  
  
  DiJetCandidateFinder(bool verbose, std::map<string, TF1*> func1D, std::map<string, TF2*> func2D):  
    verbose_(verbose), eventFlag_(-99), numJetCount_(-99), numJetRecov_(-99), error_(0) {
    
    for(std::map<string, TF1*>::iterator it = func1D.begin(); it!=func1D.end(); it++){
      func1D_[it->first] = it->second;
    }
    for(std::map<string, TF2*>::iterator it = func2D.begin(); it!=func2D.end(); it++){
      func2D_[it->first] = it->second;
    }
    
  }

    void run(int, std::vector<JetByPt>&, std::vector<JetByPt>&, std::vector<LV>& );
    int jetCountMult();
    int jetRecoverMult();
    int eventFlag();
    int errorFlag();
 
    struct sorterByPt {
      bool operator() (float i,float j) const { return (i>j);}
    };
    
    struct sorterByJetPt {
      bool operator() (JetByPt i, JetByPt j) const { return (i.pt>j.pt);}
    };
    
    void runLJ(std::vector<JetByPt>&, std::vector<JetByPt>&, std::vector<LV>&);
    void runLL(std::vector<JetByPt>&, std::vector<JetByPt>&, std::vector<LV>&);
    
    vector<Permutation> findMatchLJ(std::vector<JetByPt>&, std::vector<LV>&, std::vector<unsigned int>&);
    std::pair<int, unsigned int> scanLJ(std::vector<Permutation>); 
    void splitJets(std::vector<JetByPt>&);
    
    std::vector<LV> buildNeutrinos(LV, LV);
    void sortUnMasked(Permutation&, std::vector<Permutation>, std::vector<JetByPt>&); 
    bool bookKeeper(unsigned int, unsigned int, unsigned int, unsigned int, vector<unsigned int>&);
   
 private:
    
    bool verbose_;
    int eventFlag_;
    int numJetCount_;
    int numJetRecov_;
    int error_;
    
    std::map<string, TF1*> func1D_;
    std::map<string, TF2*> func2D_;
    
};

#endif
