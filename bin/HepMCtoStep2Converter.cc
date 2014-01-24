#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include <cstdlib>
#include <iostream> 
#include <fstream>
#include <map>
#include <string>
#include <algorithm>

#include "TSystem.h"
#include "TROOT.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVectorD.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include "DataFormats/FWLite/interface/LuminosityBlock.h"
#include "DataFormats/FWLite/interface/Run.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"

#include "HepMC/IO_GenEvent.h"
#include "HepMC/GenEvent.h"
#include "fastjet/ClusterSequence.hh"

#include <math.h>
#include <algorithm>
#include <list>

// dR distance for reco-gen matching
#define GENJETDR  0.3

using namespace std;

typedef struct 
{
  float et; 
  float sumet;
  float sig;
  float phi;
} metInfo;


typedef struct 
{
  int run;
  int lumi;
  int event;
  int json;
} EventInfo;



typedef struct
{
  float bmass;
  float bpt;
  float beta;
  float bphi;
  float bstatus;
  float wdau1mass;
  float wdau1pt;
  float wdau1eta;
  float wdau1phi;
  float wdau1id;
  float wdau2mass;
  float wdau2pt;
  float wdau2eta;
  float wdau2phi;
  float wdau2id;
} genTopInfo;

typedef struct 
{
  float mass; 
  float pt;
  float eta;
  float phi;
  float status;
  float charge;
  float momid;
} genParticleInfo;




float dR( TLorentzVector reco, TLorentzVector gen){
  return TMath::Sqrt( (reco.Eta()-gen.Eta())*(reco.Eta()-gen.Eta()) +  TMath::ACos( TMath::Cos(reco.Phi()-gen.Phi()) )*TMath::ACos( TMath::Cos(reco.Phi()-gen.Phi()) ) )  ;
}

class IsStateFinal {
public:
  bool operator()( const HepMC::GenParticle* p ) {
    if ( !p->end_vertex() && p->status()==1 ) return 1;
    return 0;
  }
};


class IsLquark {
public:
  bool operator()( const HepMC::GenParticle* p ) {
    if ( ( abs( p->pdg_id() )<5 || p->pdg_id()==21 ) && p->status()==1 && !p->end_vertex() ) return 1; // check here !!!
    return 0;
  }
};

class IsCquark {
public:
  bool operator()( const HepMC::GenParticle* p ) {
    if ( abs(p->pdg_id())==4 && p->status()==1  && !p->end_vertex() ) return 1; // check here !!!
    return 0;
  }
};

class IsBquark {
public:
  bool operator()( const HepMC::GenParticle* p ) {
    if ( abs(p->pdg_id())==5 && p->status()==1  && !p->end_vertex() ) return 1; // check here !!!
    return 0;
  }
};

class IsHiggs {
public:
  bool operator()( const HepMC::GenParticle* p ) {
    if ( abs(p->pdg_id())==25 && p->status()==1  && !p->end_vertex() ) return 1; // check here !!!
    return 0;
  }
};

class IsTquark {
public:
  bool operator()( const HepMC::GenParticle* p ) {
    if ( abs(p->pdg_id())==6 && p->status()==1  && !p->end_vertex() ) return 1; // check here !!!
    return 0;
  }
};



class IsSLEvent {
public:
  bool operator()( const HepMC::GenEvent* evt ) {
    int counterlep = 0;
    for ( HepMC::GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); ++p ){
      int pdgid = (*p)->pdg_id();
      float pt  = (*p)->momentum().perp();
      float eta = TMath::Abs( (*p)->momentum().eta() );
      if(  (pdgid == 11 || pdgid == 13 ) && pt > 25.&& eta < 2.5 ) counterlep++;	  
    }
    if( counterlep==1 ) return 1;
    return 0;
  }
};

class IsDLEvent {
public:
  bool operator()( const HepMC::GenEvent* evt ) {
    int counterlep  = 0;
    int countcharge = 0;
    for ( HepMC::GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); ++p ){
      int pdgid = (*p)->pdg_id();
      float pt  = (*p)->momentum().perp();
      float eta = TMath::Abs( (*p)->momentum().eta() );
      if(  (pdgid == 11 || pdgid == 13 ) && pt > 25.&& eta < 2.5 ){
	counterlep++;	  
	countcharge += abs(pdgid)/pdgid;
      }
    }
    if( counterlep==2 && countcharge==0 ) return 1;
    return 0;
  }
};
 


int main(int argc, const char* argv[])
{

  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
  /* @@@@@@@@@@@@@@@@@@@@@@@@ FWLITE @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */

 
  std::cout << "HEPMC --> Step2 converter" << std::endl;
  gROOT->SetBatch(true);
 
  gSystem->Load("libFWCoreFWLite");
  gSystem->Load("libDataFormatsFWLite");

  AutoLibraryLoader::enable();


  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
  /* @@@@@@@@@@@@@@@@@@@@@@@ CONFIGURATION @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */


  PythonProcessDesc builder(argv[1]);
  const edm::ParameterSet& in = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteInput");

  // SAMPLES
  std::string         outFileName    ( in.getParameter<std::string>                 ("outFileName" ) );
  vector<std::string> pathToFile     ( in.getParameter< std::vector<std::string> >  ("pathToFile"  ) );
  
  // FLAGS
  bool   verbose             ( in.getParameter<bool>         ("verbose" ) );
  bool   filter              ( in.getParameter<bool>         ("filter" ) );


  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
  /* @@@@@@@@@@@@@@@@@@@@@@@@@ INITIALIZE @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */

  // a random engine (needed for GEN+SMEAR analysis)
  TRandom3* ran = 0;

  // a clock to monitor the integration CPU time
  TStopwatch* clock = new TStopwatch();

  // a clock to monitor the total job
  TStopwatch* clock2 = new TStopwatch();
  clock2->Start();

  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
  /* @@@@@@@@@@@@@@@@@@@@@@@ OUTPUT @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */


  // clean output file (if any)
  gSystem->Exec(("rm "+outFileName).c_str());

  // output file
  TFile* fout_tmp = TFile::Open(outFileName.c_str(),"UPDATE");

  // total event counter for normalization
  TH1F*  hcounter = new TH1F("CountWithPU","",1,0,1);
  
  // output tree
  TTree* tree  = new TTree("tree","");


  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
  /* @@@@@@@@@@@@@@@@@@@@@@@@@ OPEN FILES @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
 


    
  // variables to be used from the input files
  genParticleInfo genB, genBbar;
  genTopInfo      genTop, genTbar;
  metInfo         METtype1p2corr;
  EventInfo       EVENT;
  int             nvlep, nSimBs, nhJets, naJets, nPVs, Vtype;
  float           PUweight, PUweightP, PUweightM;
  float           lheNj;
  float           weightTrig2012;
  UChar_t         triggerFlags[70];
  Int_t   vLepton_type      [2];
  Float_t vLepton_mass      [2];
  Float_t vLepton_pt        [2];
  Float_t vLepton_eta       [2];
  Float_t vLepton_phi       [2];
  Float_t vLepton_genPt     [2];
  Float_t vLepton_genEta    [2];
  Float_t vLepton_genPhi    [2];
  Float_t vLepton_charge    [2];
  Float_t vLepton_pfCorrIso [2];
  Float_t hJet_pt           [999];
  Float_t hJet_eta          [999];
  Float_t hJet_phi          [999];
  Float_t hJet_e            [999];
  Float_t hJet_flavour      [999];
  Float_t hJet_puJetIdL     [999];
  Float_t hJet_csv_nominal  [999];
  Float_t hJet_csv_upBC     [999];
  Float_t hJet_csv_downBC   [999];
  Float_t hJet_csv_upL      [999];
  Float_t hJet_csv_downL    [999];
  Float_t hJet_JECUnc       [999];
  Float_t hJet_genPt        [999];
  Float_t hJet_genEta       [999];
  Float_t hJet_genPhi       [999];
  Float_t aJet_pt           [999];
  Float_t aJet_eta          [999];
  Float_t aJet_phi          [999];
  Float_t aJet_e            [999];
  Float_t aJet_flavour      [999];
  Float_t aJet_puJetIdL     [999];
  Float_t aJet_csv_nominal  [999];
  Float_t aJet_csv_upBC     [999];
  Float_t aJet_csv_downBC   [999];
  Float_t aJet_csv_upL      [999];
  Float_t aJet_csv_downL    [999];
  Float_t aJet_JECUnc       [999];
  Float_t aJet_genPt        [999];
  Float_t aJet_genEta       [999];
  Float_t aJet_genPhi       [999];
  float SimBsmass           [999];
  float SimBspt             [999];
  float SimBseta            [999];
  float SimBsphi            [999];
  
  tree->Branch("EVENT",            &EVENT,            "run/I:lumi/I:event/I:json/I");
  tree->Branch("PUweight",         &PUweight,         "PUweight/F");
  tree->Branch("PUweightP",        &PUweightP,        "PUweightP/F");
  tree->Branch("PUweightM",        &PUweightM,        "PUweightM/F");
  tree->Branch("lheNj",            &lheNj,            "lheNj/F"); 
  tree->Branch("weightTrig2012",   &weightTrig2012,   "weightTrig2012/F"); 
  tree->Branch("triggerFlags",     triggerFlags,      "triggerFlags[70]/b"); 
  tree->Branch("Vtype",            &Vtype,            "Vtype/I");        
  tree->Branch("nhJets",           &nhJets,           "nhJets/I");
  tree->Branch("naJets",           &naJets,           "naJets/I");
  tree->Branch("nSimBs",           &nSimBs,           "nSimBs/I");
  tree->Branch("nvlep",            &nvlep,            "nvlep/I");
  tree->Branch("nPVs",             &nPVs,             "nPVs/I");
  tree->Branch("genB",             &genB,             "mass/F:pt/F:eta:phi/F:status/F:charge:momid/F");
  tree->Branch("genBbar",          &genBbar,          "mass/F:pt/F:eta:phi/F:status/F:charge:momid/F");
  tree->Branch("genTop",           &genTop,           "bmass/F:bpt/F:beta:bphi/F:bstatus/F:wdau1mass/F:wdau1pt/F:wdau1eta:wdau1phi/F:wdau1id/F:wdau2mass/F:wdau2pt/F:wdau2eta:wdau2phi/F:wdau2id/F");
  tree->Branch("genTbar",          &genTbar,          "bmass/F:bpt/F:beta:bphi/F:bstatus/F:wdau1mass/F:wdau1pt/F:wdau1eta:wdau1phi/F:wdau1id/F:wdau2mass/F:wdau2pt/F:wdau2eta:wdau2phi/F:wdau2id/F");
  tree->Branch("METtype1p2corr",   &METtype1p2corr,   "et/F:sumet:sig/F:phi/F");
  tree->Branch("vLepton_charge",   vLepton_charge,    "vLepton_charge[nvlep]/F");
  tree->Branch("vLepton_mass"  ,   vLepton_mass,      "vLepton_mass[nvlep]/F");
  tree->Branch("vLepton_pt"    ,   vLepton_pt,        "vLepton_pt[nvlep]/F");
  tree->Branch("vLepton_eta"   ,   vLepton_eta,       "vLepton_eta[nvlep]/F");
  tree->Branch("vLepton_phi"   ,   vLepton_phi,       "vLepton_phi[nvlep]/F");
  tree->Branch("vLepton_genPt" ,   vLepton_genPt,     "vLepton_genPt[nvlep]/F");
  tree->Branch("vLepton_genEta",   vLepton_genEta,    "vLepton_genEta[nvlep]/F");
  tree->Branch("vLepton_genPhi",   vLepton_genPhi,    "vLepton_genPhi[nvlep]/F");
  tree->Branch("vLepton_pfCorrIso",vLepton_pfCorrIso, "vLepton_pfCorrIso[nvlep]/F");
  tree->Branch("vLepton_type",     vLepton_type,      "vLepton_type[nvlep]/F");

  tree->Branch("hJet_pt",          hJet_pt,         "hJet_pt[nhJets]/F");    
  tree->Branch("hJet_eta",         hJet_eta,        "hJet_eta[nhJets]/F");    
  tree->Branch("hJet_phi",         hJet_phi,        "hJet_phi[nhJets]/F");    
  tree->Branch("hJet_e",           hJet_e,          "hJet_e[nhJets]/F"); 
  tree->Branch("hJet_flavour",     hJet_flavour,    "hJet_flavour[nhJets]/F");    
  tree->Branch("hJet_puJetIdL",    hJet_puJetIdL,   "hJet_puJetIdL[nhJets]/F");
  tree->Branch("hJet_csv_nominal", hJet_csv_nominal,"hJet_csv_nominal[nhJets]/F");
  tree->Branch("hJet_csv_upBC",    hJet_csv_upBC,   "hJet_csv_upBC[nhJets]/F");
  tree->Branch("hJet_csv_downBC",  hJet_csv_downBC, "hJet_csv_downBC[nhJets]/F");
  tree->Branch("hJet_csv_upL",     hJet_csv_upL,    "hJet_csv_upL[nhJets]/F");
  tree->Branch("hJet_csv_downL",   hJet_csv_downL,  "hJet_csv_downL[nhJets]/F");
  tree->Branch("hJet_JECUnc",      hJet_JECUnc,     "hJet_JECUnc[nhJets]/F");
  tree->Branch("hJet_genPt",       hJet_genPt,      "hJet_genPt[nhJets]/F");
  tree->Branch("hJet_genEta",      hJet_genEta,     "hJet_genEta[nhJets]/F");
  tree->Branch("hJet_genPhi",      hJet_genPhi,     "hJet_genPhi[nhJets]/F");

  tree->Branch("aJet_pt",          aJet_pt,         "aJet_pt[naJets]/F");    
  tree->Branch("aJet_eta",         aJet_eta,        "aJet_eta[naJets]/F");    
  tree->Branch("aJet_phi",         aJet_phi,        "aJet_phi[naJets]/F");    
  tree->Branch("aJet_e",           aJet_e,          "aJet_e[naJets]/F"); 
  tree->Branch("aJet_flavour",     aJet_flavour,    "aJet_flavour[naJets]/F");    
  tree->Branch("aJet_puJetIdL",    aJet_puJetIdL,   "aJet_puJetIdL[naJets]/F");
  tree->Branch("aJet_csv_nominal", aJet_csv_nominal,"aJet_csv_nominal[naJets]/F");
  tree->Branch("aJet_csv_upBC",    aJet_csv_upBC,   "aJet_csv_upBC[naJets]/F");
  tree->Branch("aJet_csv_downBC",  aJet_csv_downBC, "aJet_csv_downBC[naJets]/F");
  tree->Branch("aJet_csv_upL",     aJet_csv_upL,    "aJet_csv_upL[naJets]/F");
  tree->Branch("aJet_csv_downL",   aJet_csv_downL,  "aJet_csv_downL[naJets]/F");
  tree->Branch("aJet_JECUnc",      aJet_JECUnc,     "aJet_JECUnc[naJets]/F");
  tree->Branch("aJet_genPt",       aJet_genPt,      "aJet_genPt[naJets]/F");
  tree->Branch("aJet_genEta",      aJet_genEta,     "aJet_genEta[naJets]/F");
  tree->Branch("aJet_genPhi",      aJet_genPhi,     "aJet_genPhi[naJets]/F");


  tree->Branch("SimBs_mass",       SimBsmass, "SimBsmass[nSimBs]/F");
  tree->Branch("SimBs_pt",         SimBspt,   "SimBspt[nSimBs]/F");
  tree->Branch("SimBs_eta",        SimBseta,  "SimBseta[nSimBs]/F");
  tree->Branch("SimBs_phi",        SimBsphi,  "SimBsphi[nSimBs]/F");
 
  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
  /* @@@@@@@@@@@@@@@@@@@@@@@@@ EVENT LOOP @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */      

  // declare an instance of the event selection predicate
  IsSLEvent    is_SL_event;
  IsDLEvent    is_DL_event;
  IsStateFinal is_fs_particle;
  
  IsLquark is_l_quark;
  IsBquark is_b_quark;
  IsCquark is_c_quark;
  IsTquark is_t_quark;
  IsHiggs  is_higgs;
  

  int icount = 0;
  for(unsigned int file = 0; file < pathToFile.size() ; file++){

    HepMC::IO_GenEvent ascii_in( pathToFile[file] ,std::ios::in);

    if ( ascii_in.rdstate() == std::ios::failbit ) {
      cout << "File " <<   pathToFile[file] << " not found" << endl;
      continue;
    }
      
    
    HepMC::GenEvent* evt = ascii_in.read_next_event();

    while ( evt ) {

      icount++;
      
      if ( icount%200==1 ) 
	std::cout << "Processing event " << icount
		  << "  [event number " << evt->event_number() << "]"
		  << std::endl;
      
      
      if ( ( filter && ( is_SL_event(evt) || is_DL_event(evt) ) ) || !filter ) {
	
	
	// get all particles
	std::list<HepMC::GenParticle*> tquarks;
	std::list<HepMC::GenParticle*> bquarks;
	std::list<HepMC::GenParticle*> cquarks;
	std::list<HepMC::GenParticle*> lquarks;
	std::list<HepMC::GenParticle*> higgs;
	std::list<HepMC::GenParticle*> finalstateparticles;
	
	for ( HepMC::GenEvent::particle_iterator p = evt->particles_begin();
	      p != evt->particles_end(); ++p ) {
	  
	  if ( is_fs_particle(*p) && 
	       !is_t_quark(*p) && 
	       !is_higgs(*p) ) 
	    finalstateparticles.push_back(*p);

	  if ( is_t_quark(*p) )     tquarks.push_back(*p);
	  if ( is_b_quark(*p) )     bquarks.push_back(*p);
	  if ( is_c_quark(*p) )     cquarks.push_back(*p);
	  if ( is_l_quark(*p) )     lquarks.push_back(*p);
	  if ( is_higgs(*p)   )     higgs.push_back(*p);
	}
	
	if( verbose ){
	  cout << "Event " << icount << ", found:" << endl;
	  cout << " -- " << finalstateparticles.size() << " FS particles" << endl; 
	  cout << " -- " << tquarks.size() << " t quarks" << endl;
	  cout << " -- " << bquarks.size() << " b quarks" << endl;
	  cout << " -- " << cquarks.size() << " c quarks" << endl;
	  cout << " -- " << lquarks.size() << " l quarks" << endl;
	  cout << " -- " << higgs.size() << " higgs" << endl;
	}
	
      vector<TLorentzVector> leptons ;
      vector<int> leptonscharge;

      vector<TLorentzVector> jets ;
      vector<int> jetscharge;

      TLorentzVector invisible(0.,0.,0.,0) ;
      
      // get the leptons and MET
      vector<fastjet::PseudoJet> input_particles;
      for ( std::list<HepMC::GenParticle*>::iterator p_fs = finalstateparticles.begin() ; 
	    p_fs!=finalstateparticles.end() ; p_fs++  ){
	
	int pdgid  = (*p_fs)->pdg_id();
	float pt   = (*p_fs)->momentum().perp();
	float eta  = (*p_fs)->momentum().eta();
	float phi  = (*p_fs)->momentum().phi();
	float m    = (*p_fs)->momentum().m();
	float px   = (*p_fs)->momentum().px();
	float py   = (*p_fs)->momentum().py();
	float pz   = (*p_fs)->momentum().pz();
	float e    = (*p_fs)->momentum().e();

	input_particles.push_back(fastjet::PseudoJet(px,py,pz,e)); 

	if( (abs(pdgid)==11 || abs(pdgid)==13) && pt>30 && TMath::Abs(eta)<2.5 ){
	  TLorentzVector lep;
	  lep.SetPtEtaPhiM( pt, eta, phi, m);
	  leptons.push_back( lep );
	  leptonscharge.push_back( pdgid );
	}
	if( (abs(pdgid)==12 || abs(pdgid)==14 || abs(pdgid)==16) ){
	  TLorentzVector nu;
	  nu.SetPtEtaPhiM( pt, eta, phi, m);
	  invisible += nu;
	}
      }

      // jet clustering algo
      float sumEt = 0.;
      double R = 0.5;
      fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
      
      // run fastjet
      fastjet::ClusterSequence clust_seq(input_particles, jet_def);
 
      // get the resulting jets ordered in pt     
      double ptmin = 5.0;
      vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets(ptmin));

      if( verbose ) cout << "Found " << inclusive_jets.size() << " jets" << endl;
      for(unsigned int j = 0 ; j < inclusive_jets.size() ; j++){

	// count if the jet overlaps with a lepton
	int overlap = 0;

	// the pseudo jet from fastjet
	fastjet::PseudoJet pjet_j = inclusive_jets[j] ;	
	TLorentzVector jet_j( pjet_j.px(),  pjet_j.py(),  pjet_j.pz(),  pjet_j.e() );

	if( verbose ) cout << "Jet #" << j << ": (" << jet_j.Pt() << "," << jet_j.Eta() << "," << jet_j.Phi() << "," << jet_j.M() << ")" << endl; 

	// loop over selected leptons
	for(unsigned int l = 0 ; l < leptons.size() ; l++){
	  if( dR( leptons[l],jet_j)<0.50 ) overlap++;
	}

	// if the jet does NOT overlap, then...
	if(overlap==0){
	  
	  // save it
	  jets.push_back( jet_j );

	  // check for matching
	  int closestPdg  = -99;
	  float closestdR = 999.;

	  for ( std::list<HepMC::GenParticle*>::iterator p = bquarks.begin() ; p != bquarks.end() ; ++p ){
	    TLorentzVector part( ((*p)->momentum()).px(),  ((*p)->momentum()).py(),  ((*p)->momentum()).pz(),  ((*p)->momentum()).e() );
	    float dist = dR(part,jet_j) ;
	    if( dist < GENJETDR && dist < closestdR ){
	      closestdR  = dist;
	      closestPdg = (*p)->pdg_id();
	    }
	  }
	  for ( std::list<HepMC::GenParticle*>::iterator p = cquarks.begin() ; p!=cquarks.end() ; p++ ){
	    TLorentzVector part( ((*p)->momentum()).px(),  ((*p)->momentum()).py(),  ((*p)->momentum()).pz(),  ((*p)->momentum()).e() );
	    float dist = dR(part,jet_j) ;
	    if( dist < GENJETDR && dist < closestdR ){
	      closestdR  = dist;
	      closestPdg = (*p)->pdg_id();
	    }
	  }
	  for ( std::list<HepMC::GenParticle*>::iterator p = lquarks.begin() ; p!=lquarks.end() ; p++ ){
	    TLorentzVector part( ((*p)->momentum()).px(),  ((*p)->momentum()).py(),  ((*p)->momentum()).pz(),  ((*p)->momentum()).e() );
	    float dist = dR(part,jet_j) ;
	    if( dist < GENJETDR && dist < closestdR ){
	      closestdR  = dist;
	      closestPdg = (*p)->pdg_id();
	    }
	  }

	  if( verbose ) cout << closestPdg << endl;
	  jetscharge.push_back( closestPdg );


	}
      } // incluisve_jets


      /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
      /* @@@@@@@@@@@@@@@@@@@@@@@@@ FILL !!!! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */      

      nvlep = leptons.size();

      if( !(nvlep==1 || nvlep==2) && filter){
	ascii_in >> evt;
	continue;
      }

      for(unsigned int l = 0 ; l < (unsigned int)nvlep; l++ ){
	sumEt += leptons[l].Pt();
	vLepton_type[l]      = abs(leptonscharge[l]) ;
	vLepton_mass[l]      = leptons[l].M();
	vLepton_pt[l]        = leptons[l].Pt();
	vLepton_eta[l]       = leptons[l].Eta();
	vLepton_phi[l]       = leptons[l].Phi();
	vLepton_genPt[l]     = leptons[l].Pt();
	vLepton_genEta[l]    = leptons[l].Eta();
	vLepton_genPhi[l]    = leptons[l].Phi();
	vLepton_charge[l]    = abs( leptonscharge[l] )/leptonscharge[l];
	vLepton_pfCorrIso[l] = 0.;
      }


      if(nvlep==2 && abs(leptonscharge[0])==13 && abs(leptonscharge[1])==13 )
	Vtype = 0;
      else if(nvlep==2 && abs(leptonscharge[0])==11 && abs(leptonscharge[1])==11 )
	Vtype = 1;
      else if(nvlep==2 && ( (abs(leptonscharge[0])==13 && abs(leptonscharge[1])==11) ||
			    (abs(leptonscharge[1])==13 && abs(leptonscharge[0])==11)
			    ))	
	Vtype = 4;
      else if(nvlep==1 && abs(leptonscharge[0])==13 )
	Vtype = 2;
      else if(nvlep==1 && abs(leptonscharge[0])==11 )
	Vtype = 3;
      else
	Vtype = -99;


      if(jets.size()<2 /*&& filter*/){
	ascii_in >> evt;
	continue;
      }
      
      nhJets = 2;
      naJets = jets.size()-2;

      for(unsigned int j = 0 ; j < jets.size() ; j++){

	sumEt += jets[j].Pt();

	// fill hJets
	if( j<2 ){
	  hJet_pt [j]     = jets[j].Pt();
	  hJet_eta[j]     = jets[j].Eta();
	  hJet_phi[j]     = jets[j].Phi();
	  hJet_e  [j]     = jets[j].E();
	  hJet_flavour[j] = jetscharge[j];

	  hJet_puJetIdL[j]   = 1.0;
	  hJet_csv_nominal[j]= 1.0;
	  hJet_csv_upBC[j]   = 1.0;
	  hJet_csv_downBC[j] = 1.0;
	  hJet_csv_upL[j]    = 1.0;
	  hJet_csv_downL[j]  = 1.0;
	  hJet_JECUnc[j]     = 1.0;
	  hJet_genPt[j]      = jets[j].Pt();
	  hJet_genEta [j]    = jets[j].Eta();
	  hJet_genPhi [j]    = jets[j].Phi();
	}
	// fill aJets
	else{
	  aJet_pt [j-2]     = jets[j].Pt();
	  aJet_eta[j-2]     = jets[j].Eta();
	  aJet_phi[j-2]     = jets[j].Phi();
	  aJet_e  [j-2]     = jets[j].E();
	  aJet_flavour[j-2] = jetscharge[j];

	  aJet_puJetIdL[j-2]   = 1.0;
	  aJet_csv_nominal[j-2]= 1.0;
	  aJet_csv_upBC[j-2]   = 1.0;
	  aJet_csv_downBC[j-2] = 1.0;
	  aJet_csv_upL[j-2]    = 1.0;
	  aJet_csv_downL[j-2]  = 1.0;
	  aJet_JECUnc[j-2]     = 1.0;
	  aJet_genPt[j-2]      = jets[j].Pt();
	  aJet_genEta [j-2]    = jets[j].Eta();
	  aJet_genPhi [j-2]    = jets[j].Phi();
	}
      }


      nSimBs = bquarks.size();
      int b_iter = 0;
      for ( std::list<HepMC::GenParticle*>::iterator p = bquarks.begin() ; p!=bquarks.end() ; p++ ){
	SimBsmass[b_iter] = ((*p)->momentum()).m();
	SimBspt  [b_iter] = ((*p)->momentum()).perp();
	SimBseta [b_iter] = ((*p)->momentum()).eta();
	SimBsphi [b_iter] = ((*p)->momentum()).phi();
	b_iter++;
      }

      EVENT.event = evt->event_number();
      EVENT.lumi  = 1;
      EVENT.run   = 1;
      EVENT.json  = 1;

      PUweight  = 1.0;
      PUweightP = 1.0;
      PUweightM = 1.0;
      lheNj     = 0.;
      nPVs      = 1;
      weightTrig2012 = 1.0;

      for(int k = 0; k < 70 ; k++){
	triggerFlags[k] = 1;
      }

      METtype1p2corr.et     = invisible.Pt();
      METtype1p2corr.phi    = invisible.Phi();
      METtype1p2corr.sig    = 1.0;
      METtype1p2corr.sumet  = sumEt;


      for ( std::list<HepMC::GenParticle*>::iterator h = higgs.begin() ; h!=higgs.end() ; h++ ){
	if ( (*h)->end_vertex() ) {
	  for ( HepMC::GenVertex::particle_iterator des = (*h)->end_vertex()->particles_begin(HepMC::descendants);
		des != (*h)->end_vertex()->particles_end(HepMC::descendants); ++des ) {
	    if( abs((*des)->pdg_id())==+5 && (*des)->status()==2 ){
	      genB.mass   = (*des)->momentum().m();
	      genB.pt     = (*des)->momentum().perp();
	      genB.eta    = (*des)->momentum().eta();
	      genB.phi    = (*des)->momentum().phi();
	      genB.status = (*des)->pdg_id();
	      genB.charge = +2/3;
	      genB.momid  = (*h)->pdg_id();
	    }
	    if( abs((*des)->pdg_id())==-5 && (*des)->status()==2 ){
	      genBbar.mass   = (*des)->momentum().m();
	      genBbar.pt     = (*des)->momentum().perp();
	      genBbar.eta    = (*des)->momentum().eta();
	      genBbar.phi    = (*des)->momentum().phi();
	      genBbar.status = (*des)->pdg_id();
	      genBbar.charge = -1/3;
	      genBbar.momid  = (*h)->pdg_id();
	    }
	  }

	}
      }

      for ( std::list<HepMC::GenParticle*>::iterator t = tquarks.begin() ; t!=tquarks.end() ; t++ ){
	
	int pdg = (*t)->pdg_id();	

	genTopInfo* thetop;
	if     ( pdg==+6 ) thetop = &genTop;
	else if( pdg==-6)  thetop = &genTbar;
	else{
	  cout << "Error: top" << endl;
	  ascii_in >> evt;
	  continue;
	}

	if ( (*t)->end_vertex() ) {
	  for ( HepMC::GenVertex::particle_iterator des = (*t)->end_vertex()->particles_begin(HepMC::descendants);
		des != (*t)->end_vertex()->particles_end(HepMC::descendants); ++des ) {

	    // first, deal w/ the b-quark
	    if( abs((*des)->pdg_id())==5 && (*des)->status()==1 ){
	      thetop->bmass   = (*des)->momentum().m();
	      thetop->bpt     = (*des)->momentum().perp();
	      thetop->beta    = (*des)->momentum().eta();
	      thetop->bphi    = (*des)->momentum().phi();
	      thetop->bstatus = (*des)->pdg_id();
	    }

	    // then, deal w/ the W boson
	    if( abs((*des)->pdg_id())==24 && (*des)->status()==1 && (*des)->end_vertex()!=0 ){
	      int countdau = 0;
	      for ( HepMC::GenVertex::particle_iterator des2 = (*des)->end_vertex()->particles_begin(HepMC::descendants);
		    des2 != (*des)->end_vertex()->particles_end(HepMC::descendants); ++des2 ) {
		int pdg2 = (*des2)->pdg_id();
		int isWdaughter = abs(pdg2)<=5 || ( abs(pdg2)>=11 && abs(pdg2)<=16 ) ;
		if(isWdaughter && countdau==0){
		  thetop->wdau1mass   = (*des2)->momentum().m();
		  thetop->wdau1pt     = (*des2)->momentum().perp();
		  thetop->wdau1eta    = (*des2)->momentum().eta();
		  thetop->wdau1phi    = (*des2)->momentum().phi();
		  thetop->wdau1id     = (*des2)->pdg_id();
		}
		else if(isWdaughter && countdau==1){
		  thetop->wdau2mass   = (*des2)->momentum().m();
		  thetop->wdau2pt     = (*des2)->momentum().perp();
		  thetop->wdau2eta    = (*des2)->momentum().eta();
		  thetop->wdau2phi    = (*des2)->momentum().phi();
		  thetop->wdau2id     = (*des2)->pdg_id();
		}
		else{}
		countdau++;		
	      }
	    }


	  }
	}

      }
      
      tree->Fill();
      
      } // is good event

      delete evt;
      ascii_in >> evt;
      
    } // event loop
    
  } // input file loop


  // save the tree and the counting histo in the ROOT file
  fout_tmp->cd();
  hcounter->SetBinContent(1, icount);
  hcounter->Write("",TObject::kOverwrite );
  tree->Write("",TObject::kOverwrite );
  fout_tmp->Close();

  // find out the total time needed for the job
  clock2->Stop();
  float totalTime = clock2->RealTime();
  cout << "*******************" << endl;
  cout << "Job done in " << totalTime/60. << " min." << endl;

  // delete MEIntegrator and other allocated objects
  cout << "Delete clocks and random engine..." << endl;
  delete clock; delete clock2;
  if(ran!=0)    delete ran;
  cout << "Finished!!!" << endl;
  cout << "*******************" << endl;
  
  // Done!!!
  return 0;
}
