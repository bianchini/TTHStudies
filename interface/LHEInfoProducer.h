#ifndef Bianchi_TTHStudies_LHEInfoProducer_h
#define Bianchi_TTHStudies_LHEInfoProducer_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"


class LHEInfoProducer : public edm::EDAnalyzer{

 public:

  explicit LHEInfoProducer(const edm::ParameterSet &);
  ~LHEInfoProducer();

  void beginJob() ;
  void analyze(const edm::Event&  iEvent, const edm::EventSetup& iSetup);
  void endJob() ;

 private:

  TTree* tree_;
  float lheHT_, lheVpt_, lheVM_;
  int lheNj_, lheNjLF_, lheNjC_, lheNjB_, nSimBs_;
  unsigned long run_,event_,lumi_;
  float SimBs_pt[100];
  float SimBs_eta[100];
  float SimBs_phi[100];
  float lheBs_pt[100];
  float lheBs_eta[100];
  float lheBs_phi[100];

};


#endif
