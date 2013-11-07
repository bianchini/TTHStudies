#include "Bianchi/TTHStudies/interface/LHEInfoProducer.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/Common/interface/Handle.h"
#include "ZSV/BAnalysis/interface/SimBHadron.h"
#include "TLorentzVector.h"
#include <vector>
#include <utility>

using namespace std;

LHEInfoProducer::LHEInfoProducer(const edm::ParameterSet & iConfig){
}


void LHEInfoProducer::beginJob(){

  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("tree","tree");

  tree_->Branch("lheHT",   &lheHT_,  "lheHT/F");
  tree_->Branch("lheVpt",  &lheVpt_, "lheVpt/F");
  tree_->Branch("lheVM",   &lheVM_,  "lheVM/F");
  tree_->Branch("lheNj",   &lheNj_,  "lheNj/I");
  tree_->Branch("lheNjLF", &lheNjLF_,"lheNjLF/I");
  tree_->Branch("lheNjC",  &lheNjC_, "lheNjC/I");
  tree_->Branch("lheNjB",  &lheNjB_, "lheNjB/I");
  tree_->Branch("nSimBs",  &nSimBs_, "nSimBs/I");

  tree_->Branch("SimBs_pt",   SimBs_pt,  "SimBs_pt[nSimBs]/F");
  tree_->Branch("SimBs_eta",  SimBs_eta, "SimBs_eta[nSimBs]/F");
  tree_->Branch("SimBs_phi",  SimBs_phi, "SimBs_phi[nSimBs]/F");

  tree_->Branch("lheBs_pt",   lheBs_pt,  "lheBs_pt[lheNjB]/F");
  tree_->Branch("lheBs_eta",  lheBs_eta, "lheBs_eta[lheNjB]/F");
  tree_->Branch("lheBs_phi",  lheBs_phi, "lheBs_phi[lheNjB]/F");


  tree_->Branch("run",    &run_,    "run/l");
  tree_->Branch("event",  &event_,  "event/l");
  tree_->Branch("lumi",   &lumi_,   "lumi/l");


}

LHEInfoProducer::~LHEInfoProducer(){}

void LHEInfoProducer::analyze(const edm::Event & iEvent, const edm::EventSetup & iSetup){

  lheHT_  = -99.;
  lheVM_  = -99.;
  lheVpt_ = -99.;
  lheNj_  = -99;
  lheNjLF_= -99;
  lheNjC_ = -99;
  lheNjB_ = -99;
  nSimBs_ = -99;
  for(int k = 0; k < 100; k++){
    SimBs_pt[k]  = -99.;
    SimBs_eta[k] = -99.;
    SimBs_phi[k] = -99.;
    lheBs_pt[k]  = -99.;
    lheBs_eta[k] = -99.;
    lheBs_phi[k] = -99.;
  }


  run_   = iEvent.run();
  event_ = (iEvent.eventAuxiliary()).event();
  lumi_  = iEvent.luminosityBlock();


  edm::Handle<SimBHadronCollection> SBHC;
  iEvent.getByLabel("bhadrons", SBHC);
  if(SBHC.isValid() && (SBHC.product()!=0) ){
    const SimBHadronCollection *sbhc = SBHC.product();
    nSimBs_ = sbhc->size();
    for(int k = 0; k < nSimBs_; k++){
      if(k>=100) continue;
      SimBs_pt[k]  = (*sbhc)[k].pt();
      SimBs_eta[k] = (*sbhc)[k].eta();
      SimBs_phi[k] = (*sbhc)[k].phi();
    }
  }


  edm::Handle<LHEEventProduct> evt;
  iEvent.getByLabel("source",evt);


  int it = 0;
  if( evt.isValid() ){


    bool lCheck     = false;
    bool lbarCheck  = false;
    bool vlCheck    = false;
    bool vlbarCheck =false;

    lheHT_  =0.;
    lheVM_  =0.;
    lheNj_  =0 ;
    lheNjLF_=0 ;
    lheNjC_ =0 ;
    lheNjB_ =0 ;

    TLorentzVector l,lbar,vl,vlbar, tmpB;
    TLorentzVector V_tlv(0.,0.,0.,0.);

    const lhef::HEPEUP hepeup_ = evt->hepeup();
    const std::vector<lhef::HEPEUP::FiveVector> pup_ = hepeup_.PUP; // px, py, pz, E, M

    for(unsigned int i=0; i<pup_.size(); ++i){

      int id     = hepeup_.IDUP[i];
      int status = hepeup_.ISTUP[i];
      int idabs  = TMath::Abs(id); 

      if( status == 1 && ( ( idabs == 21 ) || (idabs > 0 && idabs < 7) ) ){ // gluons and quarks
	lheHT_ += TMath::Sqrt( TMath::Power(hepeup_.PUP[i][0],2) + TMath::Power(hepeup_.PUP[i][1],2) );
	lheNj_++;

	if( ( idabs == 21 ) || (idabs > 0 && idabs < 7) ) lheNjLF_++;
	if( idabs == 4 ) lheNjC_++;
	if( idabs == 5 ){
	  lheNjB_++;
	  if(it<100){
	    tmpB.SetPxPyPzE( hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]);
	    lheBs_pt[it]  = tmpB.Pt();
	    lheBs_eta[it] = tmpB.Eta();
	    lheBs_phi[it] = tmpB.Phi();
	  }
	  it++;
	}
      }	    

      if(id==11){  l.SetPxPyPzE(    hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); lCheck=true;}
      if(id==-11){ lbar.SetPxPyPzE( hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); lbarCheck=true;}
      if(id==12){  vl.SetPxPyPzE(   hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); vlCheck=true;}
      if(id==-12){ vlbar.SetPxPyPzE(hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); vlbarCheck=true;}
      
      if(id==13){  l.SetPxPyPzE(    hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); lCheck=true;}
      if(id==-13){ lbar.SetPxPyPzE( hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); lbarCheck=true;}
      if(id==14){  vl.SetPxPyPzE(   hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); vlCheck=true;}
      if(id==-14){ vlbar.SetPxPyPzE(hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); vlbarCheck=true;}
	    
      if(id==15){  l.SetPxPyPzE(    hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); lCheck=true;}
      if(id==-15){ lbar.SetPxPyPzE( hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); lbarCheck=true;}
      if(id==16){  vl.SetPxPyPzE(   hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); vlCheck=true;}
      if(id==-16){ vlbar.SetPxPyPzE(hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); vlbarCheck=true;}
	    
    }

    if( lCheck && lbarCheck )   V_tlv = l + lbar;   // ZtoLL
    if( vlCheck && vlbarCheck ) V_tlv = vl + vlbar; // ZtoNuNu
    if( lCheck && vlbarCheck )  V_tlv = l + vlbar;  // WToLNu
    if( lbarCheck && vlCheck )  V_tlv = lbar + vl;  // WToLNu       

    lheVpt_ = V_tlv.Pt();
    lheVM_  = V_tlv.M();
    
  }


  tree_->Fill();


}


void LHEInfoProducer::endJob(){}


#include "FWCore/Framework/interface/MakerMacros.h"
 
DEFINE_FWK_MODULE(LHEInfoProducer);
