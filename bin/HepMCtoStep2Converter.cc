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
#define GENJETDR  0.5
#define NJETSMAX   99

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
  int signalID;
  float eventscale;
  float alphaQCD;
  float alphaQED;
  float xSec;
  float xSecErr;
  int id1;
  int id2;
  float x1;
  float x2;
  float scalePDF;
  float pdf1;
  float pdf2;
  float weights[8];
} GenEventInfo;

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

class IsInvisible {
public:
  bool operator()( const HepMC::GenParticle* p ) {
    if ( ( abs( p->pdg_id() )==12 || abs(p->pdg_id())==14 || abs(p->pdg_id())==16 ) && p->status()==1 && !p->end_vertex() ) return 1; // check here !!!
    return 0;
  }
};


class IsLquark {
public:
  bool operator()( const HepMC::GenParticle* p ) {
    if ( ( abs( p->pdg_id() )<4 || p->pdg_id()==21 ) && p->status()==1 && !p->end_vertex() ) return 1; // check here !!!
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

class IsLhadron {
public:
  bool operator()( const HepMC::GenParticle* p ) {
    int pdg = abs(p->pdg_id());
    if ( pdg!=22 && ( (pdg%1000/100>0 && pdg%1000/100<4) || (pdg%10000/1000>0 && pdg%10000/1000<4) ) && p->status()==1 ) return 1; // check here !!!
    return 0;
  }
};


class IsChadron {
public:
  bool operator()( const HepMC::GenParticle* p ) {

    int pdg = abs(p->pdg_id());

    int hasCharm   = (pdg%1000/100==4 || pdg%10000/1000==4);
    int isTransient = p->status()==2;
    int decaysToNonCharm = 1;
    if(  p->end_vertex() ){
      for ( HepMC::GenVertex::particle_iterator des = p->end_vertex()->particles_begin(HepMC::children);
	    des != p->end_vertex()->particles_end(HepMC::children); ++des ){
	int pdg_des = abs((*des)->pdg_id());
	if(pdg_des%1000/100==4 || pdg_des%10000/1000==4) decaysToNonCharm = 0; 
      }
    }
    if (  hasCharm && isTransient && decaysToNonCharm ) return 1; // check here !!!
    return 0;
  }
};


class IsBhadron {
public:
  bool operator()( const HepMC::GenParticle* p ) {

    int pdg = abs(p->pdg_id());

    int hasBeauty   = (pdg%1000/100==5 || pdg%10000/1000==5);
    int isTransient = p->status()==2;
    int decaysToNonBeauty = 1;
    if(  p->end_vertex() ){
      for ( HepMC::GenVertex::particle_iterator des = p->end_vertex()->particles_begin(HepMC::children);
	    des != p->end_vertex()->particles_end(HepMC::children); ++des ){
	int pdg_des = abs((*des)->pdg_id());
	if(pdg_des%1000/100==5 || pdg_des%10000/1000==5) decaysToNonBeauty = 0; 
      }
    }
    if (  hasBeauty && isTransient && decaysToNonBeauty ) return 1; // check here !!!
    return 0;
  }
};


class IsHiggs_decayed {
public:
  bool operator()( const HepMC::GenParticle* p ) {
    if ( abs(p->pdg_id())==25 && p->status()==3  && p->end_vertex() ) return 1; // check here !!!
    return 0;
  }
};

class IsTquark_decayed {
public:
  bool operator()( const HepMC::GenParticle* p ) {
    if ( abs(p->pdg_id())==6 && p->status()==3  && p->end_vertex() ) return 1; // check here !!!
    return 0;
  }
};

class IsWboson_decayed {
public:
  bool operator()( const HepMC::GenParticle* p ) {
    if ( abs(p->pdg_id())==24 && p->status()==11  && p->end_vertex() ) return 1; // check here !!!
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

class IsWboson {
public:
  bool operator()( const HepMC::GenParticle* p ) {
    if ( abs(p->pdg_id())==24 && p->status()==1  && !p->end_vertex() ) return 1; // check here !!!
    return 0;
  }
};

class IsZboson {
public:
  bool operator()( const HepMC::GenParticle* p ) {
    if ( (abs(p->pdg_id())==23 || p->pdg_id()==22)  && p->status()==1  && !p->end_vertex() ) return 1; // check here !!!
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
      if( (*p)->status()==1  && !(*p)->end_vertex() && (pdgid == 11 || pdgid == 13 ) && pt > 25.&& eta < 2.5 ) counterlep++;	  
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
      if( (*p)->status()==1  && !(*p)->end_vertex() &&  (pdgid == 11 || pdgid == 13 ) && pt > 25.&& eta < 2.5 ){
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

 
  std::cout << "HepMC --> ROOT Step2 converter" << std::endl;
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
  bool   verbose             ( in.getParameter<bool>         ("verbose")   );
  bool   topDecay            ( in.getParameter<bool>         ("topDecay")  );
  bool   higgsDecay          ( in.getParameter<bool>         ("higgsDecay"));
  bool   fragmentation       ( in.getParameter<bool>         ("fragmentation"));
  bool   shower              ( in.getParameter<bool>         ("shower"));

  bool   filter              ( in.getParameter<bool>         ("filter")    );

  double jetRadius           ( in.getParameter<double>       ("jetRadius") );
  double ptMin               ( in.getParameter<double>       ("ptMin")     );

  double ptCut               ( in.getParameter<double>       ("ptCut")     );
  double etaCut              ( in.getParameter<double>       ("etaCut")    );
  int    nJetsMin            ( in.getParameter<int>          ("nJetsMin")  );

  double rIsoLep             ( in.getParameter<double>       ("rIsoLep"));
  double dRIsoLep            ( in.getParameter<double>       ("dRIsoLep"));
  double overlapLep          ( in.getParameter<double>       ("overlapLep"));
  double ptCutLep            ( in.getParameter<double>       ("ptCutLep")     );
  double etaCutLep           ( in.getParameter<double>       ("etaCutLep")    );
 
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
  GenEventInfo    GENEVENT;

  int             nvlep, nSimBs /*, nhJets*/, naJets /*, nPVs*/, Vtype;
  //float           PUweight, PUweightP, PUweightM;
  //float           lheNj;
  //float           weightTrig2012;
  //UChar_t         triggerFlags[70];
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
  Float_t aJet_pt           [NJETSMAX];
  Float_t aJet_eta          [NJETSMAX];
  Float_t aJet_phi          [NJETSMAX];
  Float_t aJet_e            [NJETSMAX];
  Float_t aJet_flavour      [NJETSMAX];

  Float_t aJet_nBs          [NJETSMAX];
  Float_t aJet_nCs          [NJETSMAX];
  Float_t aJet_nLs          [NJETSMAX];
   
  Float_t aJet_puJetIdL     [NJETSMAX];
  Float_t aJet_csv_nominal  [NJETSMAX];
  Float_t aJet_csv_upBC     [NJETSMAX];
  Float_t aJet_csv_downBC   [NJETSMAX];
  Float_t aJet_csv_upL      [NJETSMAX];
  Float_t aJet_csv_downL    [NJETSMAX];
  Float_t aJet_JECUnc       [NJETSMAX];
  Float_t aJet_genPt        [NJETSMAX];
  Float_t aJet_genEta       [NJETSMAX];
  Float_t aJet_genPhi       [NJETSMAX];
  float SimBsmass           [NJETSMAX];
  float SimBspt             [NJETSMAX];
  float SimBseta            [NJETSMAX];
  float SimBsphi            [NJETSMAX];
  
  tree->Branch("EVENT",            &EVENT,            "run/I:lumi/I:event/I:json/I");
  tree->Branch("GENEVENT",         &GENEVENT,         "signalID/I:eventscale/F:alphaQCD/F:alphaQED/F:xSec/F:xSecErr/F:id1/I:id2/I:x1/F:x2/F:scalePDF/F:pdf1/F:pdf2/F:weights[8]/F");
 
  //tree->Branch("PUweight",         &PUweight,         "PUweight/F");
  //tree->Branch("PUweightP",        &PUweightP,        "PUweightP/F");
  //tree->Branch("PUweightM",        &PUweightM,        "PUweightM/F");
  //tree->Branch("lheNj",            &lheNj,            "lheNj/F"); 
  //tree->Branch("weightTrig2012",   &weightTrig2012,   "weightTrig2012/F"); 
  //tree->Branch("triggerFlags",     triggerFlags,      "triggerFlags[70]/b"); 
  tree->Branch("Vtype",            &Vtype,            "Vtype/I");        
  //tree->Branch("nhJets",           &nhJets,           "nhJets/I");
  tree->Branch("naJets",           &naJets,           "naJets/I");
  tree->Branch("nSimBs",           &nSimBs,           "nSimBs/I");
  tree->Branch("nvlep",            &nvlep,            "nvlep/I");
  //tree->Branch("nPVs",             &nPVs,             "nPVs/I");
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
  tree->Branch("vLepton_type",     vLepton_type,      "vLepton_type[nvlep]/I");

  tree->Branch("aJet_pt",          aJet_pt,         "aJet_pt[naJets]/F");    
  tree->Branch("aJet_eta",         aJet_eta,        "aJet_eta[naJets]/F");    
  tree->Branch("aJet_phi",         aJet_phi,        "aJet_phi[naJets]/F");    
  tree->Branch("aJet_e",           aJet_e,          "aJet_e[naJets]/F"); 
  tree->Branch("aJet_flavour",     aJet_flavour,    "aJet_flavour[naJets]/F");    

  tree->Branch("aJet_nBs",         aJet_nBs,        "aJet_nBs[naJets]/F");    
  tree->Branch("aJet_nCs",         aJet_nCs,        "aJet_nCs[naJets]/F");    
  tree->Branch("aJet_nLs",         aJet_nLs,        "aJet_nLs[naJets]/F");    

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
  
  IsInvisible is_invisible;

  IsLquark is_l_quark;
  IsBquark is_b_quark;
  IsCquark is_c_quark;

  IsLhadron is_l_hadron;
  IsBhadron is_b_hadron;
  IsChadron is_c_hadron;

  IsTquark is_t_quark;
  IsWboson is_w_boson;
  IsHiggs  is_h_boson;

  IsTquark_decayed is_t_quark_decayed;
  IsWboson_decayed is_w_boson_decayed;
  IsHiggs_decayed  is_h_boson_decayed;
  
  // total event counter
  int icount = 0;

  // total accepted event counter  
  int pcount = 0;

  // loop over input files
  for(unsigned int file = 0; file < pathToFile.size() ; file++){

    // read the event colection from the input file
    HepMC::IO_GenEvent ascii_in( pathToFile[file] ,std::ios::in);

    // if no valid input hepmc file found, go to the next file
    if ( ascii_in.rdstate() == std::ios::failbit ) {
      cout << "File " <<   pathToFile[file] << " not found" << endl;
      continue;
    }
      
    // read an event...
    HepMC::GenEvent* evt = ascii_in.read_next_event();    

    // start the event loop
    while ( evt ) {

      // increment event counter
      icount++;
      
      // to monitor the job
      if ( icount%200==1 ) 
	std::cout << "Processing event " << icount
		  << "  [event number " << evt->event_number() << "]"
		  << std::endl;
      
      // if the event passes the selection cut, or no cut at all is applied...
      if ( ( filter && ( is_SL_event(evt) || is_DL_event(evt) ) ) || !filter ) {			

	if( verbose ){
	  cout << endl;
	  cout << "*********************" << endl;
	  cout << "Event " << icount << ", found:" << endl;
	}

	// gen-level infos
	HepMC::PdfInfo* pdfInfo              = evt->pdf_info();
	HepMC::GenCrossSection* crosssection = evt->cross_section();

	GENEVENT.signalID   = evt->signal_process_id();
	GENEVENT.eventscale = evt->event_scale();
	GENEVENT.alphaQCD   = evt->alphaQCD();
	GENEVENT.alphaQED   = evt->alphaQED();

	if( pdfInfo ){
	  GENEVENT.id1      = pdfInfo->id1();
	  GENEVENT.id2      = pdfInfo->id2();
	  GENEVENT.x1       = pdfInfo->x1();
	  GENEVENT.x2       = pdfInfo->x2();
	  GENEVENT.pdf1     = pdfInfo->pdf1();
	  GENEVENT.pdf2     = pdfInfo->pdf2();
	  GENEVENT.scalePDF = pdfInfo->scalePDF();
	}
	else{
	  GENEVENT.id1      = -99;
	  GENEVENT.id2      = -99;
	  GENEVENT.x1       = -99;
	  GENEVENT.x2       = -99;
	  GENEVENT.pdf1     = -99;
	  GENEVENT.pdf2     = -99;
	  GENEVENT.scalePDF = -99;
	}

	if( crosssection ){
	  GENEVENT.xSec    = crosssection->cross_section();
	  GENEVENT.xSecErr = crosssection->cross_section_error();
	}
	else{
	  GENEVENT.xSec    = -99;
	  GENEVENT.xSecErr = -99;
	}
	
	for(int w = 0; w < 8 ; w++)
	  GENEVENT.weights[w] = -99;
	unsigned int count_wit = 0;
	for( HepMC::WeightContainer::iterator wit = (evt->weights()).begin() ; wit != (evt->weights()).end() ; wit++){
	  if(  count_wit>7 ) continue;
	  GENEVENT.weights[ count_wit ] = (*wit);
	  count_wit++;
	}


	// get all particles
	std::list<HepMC::GenParticle*> tquarks;
	std::list<HepMC::GenParticle*> bquarks;
	std::list<HepMC::GenParticle*> cquarks;
	std::list<HepMC::GenParticle*> lquarks;
	std::list<HepMC::GenParticle*> higgs;
	std::list<HepMC::GenParticle*> finalstateparticles;
	
	// loop over all particles in the event record
	for ( HepMC::GenEvent::particle_iterator p = evt->particles_begin();
	      p != evt->particles_end(); ++p ) {
	  
	  // push back final-state particles that are neither top nor higgs
	  // [ is_fs_particle() selects only status=1 particles, but if the top and higgs are undecayed
	  // they will ahve such status ]
	  if ( is_fs_particle(*p) && 
	       !is_t_quark(*p) && 
	       !is_h_boson(*p) && 
	       !is_w_boson(*p) ) 
	    finalstateparticles.push_back(*p);

	  if(!fragmentation){
	    if ( is_b_quark(*p) )     bquarks.push_back(*p);
	    if ( is_c_quark(*p) )     cquarks.push_back(*p);
	    if ( is_l_quark(*p) )     lquarks.push_back(*p);
	  }
	  else{
	    if      ( is_b_hadron(*p) ){
	      bquarks.push_back(*p);
	      if( verbose ) cout << "B-hadron " << (*p)->pdg_id() << ", " << (*p)->status() << endl;
	    }
	    else if ( is_c_hadron(*p) ){
	      cquarks.push_back(*p);
	      if( verbose ) cout << "C-hadron " << (*p)->pdg_id() << ", " << (*p)->status() << endl;
	    }
	    else if ( is_l_hadron(*p) ){
	      lquarks.push_back(*p);
	      if( verbose ) cout << "L-hadron " << (*p)->pdg_id() << ", " << (*p)->status() << endl;
	    }
	  }

	  if ( (!topDecay   && is_t_quark(*p)) || (topDecay   && is_t_quark_decayed(*p)) )  tquarks.push_back(*p);
	  if ( (!higgsDecay && is_h_boson(*p)) || (higgsDecay && is_h_boson_decayed(*p)) )  higgs.push_back(*p);

	}
	
	if( verbose ){
	  cout << " -- " << finalstateparticles.size() << " FS particles" << endl; 
	  cout << " -- " << tquarks.size() << " t quarks" << endl;
	  cout << " -- " << bquarks.size() << " b quarks" << endl;
	  cout << " -- " << cquarks.size() << " c quarks" << endl;
	  cout << " -- " << lquarks.size() << " l quarks" << endl;
	  cout << " -- " << higgs.size() << " higgs" << endl;
	}
	
	// container of leptons 4-vectors
	vector<TLorentzVector> leptons ;
	
	// container of leptons pdgid
	vector<int> leptonscharge;

	// container of leptons isolation
	vector<float> leptonsiso;
	
	// container of jets 4-vectors
	vector<TLorentzVector> jets ;

	// container of jets closest particle pdgid
	vector<int> jetscharge;
	
	// container of jets multiplicity of b-quarks
	vector<int> jetsnBs;

	// container of jets multiplicity of c-quarks
	vector<int> jetsnCs;

	// container of jets multiplicity of udcsg-quarks
	vector<int> jetsnLs;

	// total neutrinos momentum
	TLorentzVector invisible(0.,0.,0.,0) ;
      
	// get the particles for jets
	vector<fastjet::PseudoJet> input_particles;

	vector<TLorentzVector> neutrinos_from_hadron_decay;

	for ( std::list<HepMC::GenParticle*>::iterator p_fs = finalstateparticles.begin() ; 
	      p_fs!=finalstateparticles.end() ; p_fs++  ){
	  
	  float px   = (*p_fs)->momentum().px();
	  float py   = (*p_fs)->momentum().py();
	  float pz   = (*p_fs)->momentum().pz();
	  float e    = (*p_fs)->momentum().e();
	  
	  // don't throw in neutrinos !
	  if( ! is_invisible(*p_fs) )
	    input_particles.push_back(fastjet::PseudoJet(px,py,pz,e)); 	  	
	  else{
	    if( verbose ) cout << "jet inputs: this is invisible (pdgid = " << (*p_fs)->pdg_id() << ") : don't consider it!" << endl;

	    int fromHadronDecay = 1;
	    if( (*p_fs)->production_vertex() ){
	      for ( HepMC::GenVertex::particle_iterator par = (*p_fs)->production_vertex()->particles_begin(HepMC::parents);
		    par != (*p_fs)->production_vertex()->particles_end(HepMC::parents); ++par ){
		if( verbose ) cout << "parents of neutrinos: " << (*par)->pdg_id() << ", " << (*par)->status() << endl;		 
		if(  abs((*par)->pdg_id())>=11 && abs((*par)->pdg_id())<=16 ) fromHadronDecay=0;
	      }
	    }
	    if( fromHadronDecay ) neutrinos_from_hadron_decay.push_back( TLorentzVector(px,py,pz,e) );
	    	    
	  }
	}



	// get the leptons & MET
	for ( std::list<HepMC::GenParticle*>::iterator p_fs = finalstateparticles.begin() ; 
	      p_fs!=finalstateparticles.end() ; p_fs++  ){
	  
	  int pdgid  = (*p_fs)->pdg_id();
	  float pt   = (*p_fs)->momentum().perp();
	  float eta  = (*p_fs)->momentum().eta();
	  float phi  = (*p_fs)->momentum().phi();
	  float m    = (*p_fs)->momentum().m();
	  
	  // if it is a muon or electron passing the acceptance cuts
	  if( (abs(pdgid)==11 || abs(pdgid)==13) && pt>ptCutLep && TMath::Abs(eta)<etaCutLep ){

	    TLorentzVector lep;
	    lep.SetPtEtaPhiM( pt, eta, phi, m);
	    
	    // calculate the relative isolation
	    float rIso = 0.;
	    for ( std::list<HepMC::GenParticle*>::iterator p_fs2 = finalstateparticles.begin() ; 
		  p_fs2!=finalstateparticles.end() ; p_fs2++  ){
	      
	      float pt2   = (*p_fs2)->momentum().perp();
	      float eta2  = (*p_fs2)->momentum().eta();
	      float phi2  = (*p_fs2)->momentum().phi();
	      float m2    = (*p_fs2)->momentum().m();
	      
	      TLorentzVector part;
	      part.SetPtEtaPhiM( pt2, eta2, phi2, m2);
	      
	      // if the particle is close (but not overlapping)
	      if( dR( lep, part)<1e-04 && TMath::Abs(lep.E()-part.E())/part.E()<1e-04 ) continue;

	      // if the particole is in the cone, add its pt to the sum
	      if( dR( lep , part)<dRIsoLep )  rIso += part.Pt();
	    }

	    // make relative isolation
	    rIso /= (lep.Pt()>0 ? lep.Pt() : 1.0);

	    // save only isolated leptons
	    if( rIso < rIsoLep ){
	      leptons.push_back( lep );
	      leptonscharge.push_back( pdgid );
	      leptonsiso.push_back( rIso );
	    }
	  }

	  else if( (abs(pdgid)==12 || abs(pdgid)==14 || abs(pdgid)==16) ){
	    TLorentzVector nu;
	    nu.SetPtEtaPhiM( pt, eta, phi, m);
	    invisible += nu;
	  }

	}


	// the sum Et 
	float sumEt    = 0.;
	int jetCounter = 0;

	// jet clustering algo
	double R    = jetRadius;
	fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, R);
	
	// run fastjet
	fastjet::ClusterSequence clust_seq(input_particles, jet_def);
 
	// get the resulting jets ordered in pt     
	vector<fastjet::PseudoJet> inclusive_jets = sorted_by_pt(clust_seq.inclusive_jets( ptMin ));
	
	if( verbose ) cout << "Found " << inclusive_jets.size() << " jets" << endl;

	// loop over jet collection
	for(unsigned int j = 0 ; j < inclusive_jets.size() ; j++){
	  
	  // count if the jet overlaps with a lepton
	  int overlap = 0;
	  
	  // the pseudo jet from fastjet
	  fastjet::PseudoJet pjet_j = inclusive_jets[j] ;	
	  TLorentzVector jet_j( pjet_j.px(),  pjet_j.py(),  pjet_j.pz(),  pjet_j.e() );

	  // count jets above threshold
	  if(  jet_j.Pt()> ptCut && TMath::Abs( jet_j.Eta() )< etaCut ) jetCounter++;

	  if( verbose ) cout << "Jet #" << j << ": (" << jet_j.Pt() << "," << jet_j.Eta() << "," << jet_j.Phi() << "," << jet_j.M() << ")" << endl; 

	  // loop over selected leptons
	  float deltaR = -99;
	  for(unsigned int l = 0 ; l < leptons.size() ; l++){
	    deltaR = dR( leptons[l],jet_j);
	    if( deltaR < overlapLep ) overlap++;
	  }

	  // if the jet does NOT overlap, then...
	  if(overlap==0){
	  
	    // save it
	    jets.push_back( jet_j );

	    // count multiplicity of b,c, or light quarks
	    int numBs = 0;
	    int numCs = 0;
	    int numLs = 0;

	    // the fast-jet constituents
	    vector<fastjet::PseudoJet> constituents_j = clust_seq.constituents (pjet_j);

	    if(verbose) cout << "Pseudo-jet " << j << " has " << constituents_j.size() << " constituents" << endl;	  

	    // if doing a parton-level analysis, check for constituents
	    if(!fragmentation){

	      // loop over the constituents
	      for(unsigned int k = 0 ; k < constituents_j.size() ; k++){
		
		// take the kth constituent...
		fastjet::PseudoJet const_k = constituents_j[k] ;	
		TLorentzVector p_k( const_k.px(),  const_k.py(),  const_k.pz(),  const_k.e() );
		
		if( verbose ) cout << " > const #" << k << ": (" << p_k.Pt() << "," << p_k.Eta() << "," << p_k.Phi() << "," << p_k.M() << ")" << endl; 
		
		// make sure each constituent is unambiguously identified
		int matches = 0;
				
		// check of which type it is (use angular matching)
		for ( std::list<HepMC::GenParticle*>::iterator p = lquarks.begin() ; p!=lquarks.end() ; p++ ){
		  TLorentzVector part( ((*p)->momentum()).px(),  
				       ((*p)->momentum()).py(),  
				       ((*p)->momentum()).pz(),  
				       ((*p)->momentum()).e() );
		  float dist = dR(part, p_k) ;
		  if( dist<1e-04 && TMath::Abs(part.E()-p_k.E())/p_k.E()<1e-04 ){
		    numLs++;
		    matches++;
		    if(verbose) cout << "  ... matches with ligh-quark " << ": (" << part.Pt() << "," << part.Eta() << "," << part.Phi() << "," << part.M() << ")" << endl; 
		  }
		  
		}
		
		for ( std::list<HepMC::GenParticle*>::iterator p = cquarks.begin() ; p!=cquarks.end() ; p++ ){
		  TLorentzVector part( ((*p)->momentum()).px(),  
				       ((*p)->momentum()).py(),  
				       ((*p)->momentum()).pz(),  
				       ((*p)->momentum()).e() );
		  float dist = dR(part, p_k) ;
		  if( dist<1e-04 && TMath::Abs(part.E()-p_k.E())/p_k.E()<1e-04  ){
		    numCs++;
		    matches++;
		    if(verbose) cout << "  ... matches with c-quark " << ": (" << part.Pt() << "," << part.Eta() << "," << part.Phi() << "," << part.M() << ")" << endl; 
		  }
		}
		
		for ( std::list<HepMC::GenParticle*>::iterator p = bquarks.begin() ; p!=bquarks.end() ; p++ ){
		  TLorentzVector part( ((*p)->momentum()).px(),  
				       ((*p)->momentum()).py(),  
				       ((*p)->momentum()).pz(),  
				       ((*p)->momentum()).e() );
		  float dist = dR(part, p_k) ;
		  if( dist<1e-04 && TMath::Abs(part.E()-p_k.E())/p_k.E()<1e-04  ){
		    numBs++;
		    matches++;
		    if(verbose) cout << "  ... matches with b-quark " << ": (" << part.Pt() << "," << part.Eta() << "," << part.Phi() << "," << part.M() << ")" << endl; 
		  }
		}
		
		if(matches>1) cout << "*** Ambiguous match of jet constituent ***" << endl;
		
	      }
 
	    }
	    
	  
	    // check for closest particle matching
	    int closestPdg  = -99;
	    float closestdR = 999.;

	    int matchesBhad = 0;
	    int matchesChad = 0;
	    int matchesLhad = 0;
	    
	    for ( std::list<HepMC::GenParticle*>::iterator p = lquarks.begin() ; p!=lquarks.end() ; p++ ){
	      TLorentzVector part( ((*p)->momentum()).px(),  
				   ((*p)->momentum()).py(),  
				   ((*p)->momentum()).pz(),  
				   ((*p)->momentum()).e() );
	      float dist = dR(part,jet_j) ;
	      if( dist < GENJETDR ){
		matchesLhad++;
		if(dist < closestdR){
		  closestdR  = dist;
		  closestPdg = (*p)->pdg_id();
		}
	      }
	    }
	    
	    for ( std::list<HepMC::GenParticle*>::iterator p = cquarks.begin() ; p!=cquarks.end() ; p++ ){
	      TLorentzVector part( ((*p)->momentum()).px(),  
				   ((*p)->momentum()).py(),  
				   ((*p)->momentum()).pz(),  
				   ((*p)->momentum()).e() );
	      float dist = dR(part,jet_j) ;
	      if( dist < GENJETDR ){
		matchesChad++;
		if(dist < closestdR ){
		  closestdR  = dist;
		  closestPdg = (*p)->pdg_id();
		}
	      }
	    }
	    
	    for ( std::list<HepMC::GenParticle*>::iterator p = bquarks.begin() ; p != bquarks.end() ; ++p ){
	      TLorentzVector part( ((*p)->momentum()).px(),  
				   ((*p)->momentum()).py(),  
				   ((*p)->momentum()).pz(),  
				   ((*p)->momentum()).e() );
	      float dist = dR(part,jet_j) ;
	      if( dist < GENJETDR ){
		matchesBhad++;
		if(dist < closestdR){
		  closestdR  = dist;
		  closestPdg = (*p)->pdg_id();
		}
	      }
	    }	         	 
	    
	    // save the results
	    if( verbose ) cout << closestPdg << endl;
	    
	    if(!fragmentation){
	      if     ( numBs>0 ) closestPdg = 5;
	      else if( numCs>0 ) closestPdg = 4;
	      else if( numLs>0 ) closestPdg = 21;
	      else closestPdg = -99;
	    }
	    else{
	      if     ( matchesBhad>0 ) closestPdg = 5;
	      else if( matchesChad>0 ) closestPdg = 4;
	      else if( matchesLhad>0 ) closestPdg = 21;
	      else closestPdg = -99;
	    }
	    
	    jetscharge.push_back( closestPdg );
	    jetsnBs.push_back( numBs );
	    jetsnCs.push_back( numCs );
	    jetsnLs.push_back( numLs );
	    
	  }
	  else{
	    if( verbose ) cout << "Pseudo-jet matches within " << deltaR << " with a lepton: skip it!" << endl;
	  }

	} // incluisve_jets
	
	if( verbose ) cout << "Total number of jets after cleaning: " << jets.size() << endl;
	
	/* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */
	/* @@@@@@@@@@@@@@@@@@@@@@@@@ FILL !!!! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@  */      
	
	// number of leptons
	nvlep = leptons.size();

	// if not enough leptons, continue (unless filter=False)
	if( !(nvlep==1 || nvlep==2) && filter){
	  delete evt;
	  ascii_in >> evt;
	  continue;
	}
	
	for(unsigned int l = 0 ; l < (unsigned int)nvlep; l++ ){
	  sumEt                += leptons[l].Pt();
	  vLepton_type[l]      = abs( leptonscharge[l] ) ;
	  vLepton_mass[l]      = leptons[l].M();
	  vLepton_pt[l]        = leptons[l].Pt();
	  vLepton_eta[l]       = leptons[l].Eta();
	  vLepton_phi[l]       = leptons[l].Phi();
	  vLepton_genPt[l]     = leptons[l].Pt();
	  vLepton_genEta[l]    = leptons[l].Eta();
	  vLepton_genPhi[l]    = leptons[l].Phi();
	  vLepton_charge[l]    = abs( leptonscharge[l] )/leptonscharge[l];	 
	  vLepton_pfCorrIso[l] = leptonsiso[l];
	}
	

	if(nvlep==2 && abs(leptonscharge[0])==13 && abs(leptonscharge[1])==13 && leptonscharge[0]*leptonscharge[1]<0)
	  Vtype = 0;
	else if(nvlep==2 && abs(leptonscharge[0])==11 && abs(leptonscharge[1])==11 && leptonscharge[0]*leptonscharge[1]<0 )
	  Vtype = 1;
	else if(nvlep==2 && ( (abs(leptonscharge[0])==13 && abs(leptonscharge[1])==11) ||
			      (abs(leptonscharge[1])==13 && abs(leptonscharge[0])==11)
			      ) && leptonscharge[0]*leptonscharge[1]<0 )	
	  Vtype = 4;
	else if(nvlep==1 && abs(leptonscharge[0])==13 )
	  Vtype = 2;
	else if(nvlep==1 && abs(leptonscharge[0])==11 )
	  Vtype = 3;
	else
	  Vtype = -99;
	
	// if not enough jets, continue (unless filter=False)
	if( jetCounter < nJetsMin && filter ){
	  delete evt;
	  ascii_in >> evt;
	  continue;
	}
      
	//nhJets = 2;
	naJets = jets.size();
	
	// loop over jets
	for(unsigned int j = 0 ; j < (unsigned int)TMath::Min( naJets , NJETSMAX ) /*jets.size()*/ ; j++){
	  
	  sumEt += jets[j].Pt();
	  		  
	  // fill aJets
	  aJet_pt [j]     = jets[j].Pt();
	  aJet_eta[j]     = jets[j].Eta();
	  aJet_phi[j]     = jets[j].Phi();
	  aJet_e  [j]     = jets[j].E();
	  aJet_flavour[j] = jetscharge[j];

	  aJet_nBs[j] = jetsnBs[j];
	  aJet_nCs[j] = jetsnCs[j];
	  aJet_nLs[j] = jetsnLs[j];
	
	  aJet_puJetIdL[j]   = 1.0;
	  aJet_csv_nominal[j]= 1.0;
	  aJet_csv_upBC[j]   = 1.0;
	  aJet_csv_downBC[j] = 1.0;
	  aJet_csv_upL[j]    = 1.0;
	  aJet_csv_downL[j]  = 1.0;
	  aJet_JECUnc[j]     = 1.0;

	  TLorentzVector extrap4(0.,0.,0.,0.);
	  for(unsigned nus = 0 ; nus < neutrinos_from_hadron_decay.size() ; nus++){
	    if( dR(jets[j],neutrinos_from_hadron_decay[nus])<jetRadius ){
	      extrap4 += neutrinos_from_hadron_decay[nus];
	      if( verbose ) cout << "Neutrino (" 
				 << neutrinos_from_hadron_decay[nus].E() << "," 
				 << neutrinos_from_hadron_decay[nus].Eta() << "," 
				 << neutrinos_from_hadron_decay[nus].Phi() << ") is matched to jet " 
				 << j << "(" << jets[j].E() << "," << jets[j].Eta() << "," << jets[j].Phi() 
				 << "): add its 4-vector..." << endl;
	    }
	  }
	  if( verbose ) cout << "Total extrap4 energy = " << extrap4.E() << endl;

	  aJet_genPt[j]      = (jets[j]+extrap4).Pt();
	  aJet_genEta [j]    = (jets[j]+extrap4).Eta();
	  aJet_genPhi [j]    = (jets[j]+extrap4).Phi();
	  
	}
	
	
	// number of b-quarks
	nSimBs = bquarks.size();
	int b_iter = 0;
	for ( std::list<HepMC::GenParticle*>::iterator p = bquarks.begin() ; p!=bquarks.end() ; p++ ){
	  SimBsmass[b_iter] = ((*p)->momentum()).m();
	  SimBspt  [b_iter] = ((*p)->momentum()).perp();
	  SimBseta [b_iter] = ((*p)->momentum()).eta();
	  SimBsphi [b_iter] = ((*p)->momentum()).phi();
	  b_iter++;
	}

	// event information
	EVENT.event = evt->event_number();
	EVENT.lumi  = 1;
	EVENT.run   = 1;
	EVENT.json  = 1;

	//PUweight  = 1.0;
	//PUweightP = 1.0;
	//PUweightM = 1.0;
	//lheNj     = 0.;
	//nPVs      = 1;
	//weightTrig2012 = 1.0;
	//for(int k = 0; k < 70 ; k++){
	//triggerFlags[k] = 1;
	//}
	
	// the neutrinos 4-vector
	METtype1p2corr.et     = invisible.Pt();
	METtype1p2corr.phi    = invisible.Phi();
	METtype1p2corr.sig    = 1.0;
	METtype1p2corr.sumet  = sumEt;
	


	// fill Higgs decay kinematics (if any Higgs decay found)
	if( higgs.size()<1 ){
	  genB.mass      = -99;
	  genB.pt        = -99;
	  genB.eta       = -99;
	  genB.phi       = -99;
	  genB.status    = -99;
	  genB.charge    = -99;
	  genB.momid     = -99;
	  genBbar.mass   = -99;
	  genBbar.pt     = -99;
	  genBbar.eta    = -99;
	  genBbar.phi    = -99;
	  genBbar.status = -99;
	  genBbar.charge = -99;
	  genBbar.momid  = -99;
	}	
	for ( std::list<HepMC::GenParticle*>::iterator h = higgs.begin() ; h!=higgs.end() ; h++ ){

	  if ( (*h)->end_vertex() ) {

	    if(verbose ) cout << "Children of the higgs " << (*h)->pdg_id()  << ":" << endl;
	    
	    for ( HepMC::GenVertex::particle_iterator des = (*h)->end_vertex()->particles_begin(HepMC::children);
		  des != (*h)->end_vertex()->particles_end(HepMC::children); ++des ) {

	      if(verbose ) cout << " > ID " << (*des)->pdg_id() << ", status " << (*des)->status()<< ", mass " << (*des)->momentum().m() << endl;

	      int goodstatus = shower ? (*des)->status()==11 : (*des)->status()==1;

	      if( (*des)->pdg_id()==+5 && goodstatus ){
		genB.mass   = (*des)->momentum().m();
		genB.pt     = (*des)->momentum().perp();
		genB.eta    = (*des)->momentum().eta();
		genB.phi    = (*des)->momentum().phi();
		genB.status = (*des)->pdg_id();
		genB.charge = +2/3;
		genB.momid  = (*h)->pdg_id();
	      }
	      if( (*des)->pdg_id()==-5 && goodstatus ){
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
	  else{
	    genB.mass      = ((*h)->momentum()).m();
	    genB.pt        = ((*h)->momentum()).perp();
	    genB.eta       = ((*h)->momentum()).eta();
	    genB.phi       = ((*h)->momentum()).phi();
	    genB.status    = -99;
	    genB.charge    = -99;
	    genB.momid     = -99;
	    genBbar.mass   = 0;
	    genBbar.pt     = 0;
	    genBbar.eta    = 0;
	    genBbar.phi    = 0;
	    genBbar.status = 0;
	    genBbar.charge = 0;
	    genBbar.momid  = 0;
	  }
	}
	
	// fill top decay kinematics (if any too decay found)
	if( tquarks.size()<2 ){
	  genTop.bmass     = -99;
	  genTop.bpt       = -99;
	  genTop.beta      = -99;
	  genTop.bphi      = -99;
	  genTop.bstatus   = -99;
	  genTop.wdau1mass = -99;
	  genTop.wdau1pt   = -99;
	  genTop.wdau1eta  = -99;
	  genTop.wdau1phi  = -99;
	  genTop.wdau1id   = -99;
	  genTop.wdau2mass = -99;
	  genTop.wdau2pt   = -99;
	  genTop.wdau2eta  = -99;
	  genTop.wdau2phi  = -99;
	  genTop.wdau2id   = -99;
	  genTbar.bmass    = -99;
	  genTbar.bpt      = -99;
	  genTbar.beta     = -99;
	  genTbar.bphi     = -99;
	  genTbar.bstatus  = -99;
	  genTbar.wdau1mass= -99;
	  genTbar.wdau1pt  = -99;
	  genTbar.wdau1eta = -99;
	  genTbar.wdau1phi = -99;
	  genTbar.wdau1id  = -99;
	  genTbar.wdau2mass= -99;
	  genTbar.wdau2pt  = -99;
	  genTbar.wdau2eta = -99;
	  genTbar.wdau2phi = -99;
	  genTbar.wdau2id  = -99;	  
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

	    if(verbose ) cout << "Children of the top " <<  pdg << ":" << endl;

	    for ( HepMC::GenVertex::particle_iterator des = (*t)->end_vertex()->particles_begin(HepMC::children);
		  des != (*t)->end_vertex()->particles_end(HepMC::children); ++des ) {
	      
	      if(verbose ) cout << " > ID " << (*des)->pdg_id() << ", status " << (*des)->status()<< ", mass " << (*des)->momentum().m() << endl;

	      int goodstatus = shower ? (*des)->status()==11 : (*des)->status()==1;

	      // first, deal w/ the b-quark
	      if( abs((*des)->pdg_id())==5 && goodstatus ){
		thetop->bmass   = (*des)->momentum().m();
		thetop->bpt     = (*des)->momentum().perp();
		thetop->beta    = (*des)->momentum().eta();
		thetop->bphi    = (*des)->momentum().phi();
		thetop->bstatus = (*des)->pdg_id();
	      }

	      // then, deal w/ the W boson
	      if( is_w_boson_decayed( (*des) ) ){

		if(verbose ) cout << "    Children of the W :" << endl;

		int countdau = 0;
		for ( HepMC::GenVertex::particle_iterator des2 = (*des)->end_vertex()->particles_begin(HepMC::children);
		      des2 != (*des)->end_vertex()->particles_end(HepMC::children); ++des2 ) {

		  int goodstatus2 = shower ? (*des2)->status()==11 : (*des2)->status()==1;
		  int pdg2    = (*des2)->pdg_id();
		  int isWdaughter = goodstatus2 && (abs(pdg2)<5 || ( abs(pdg2)>=11 && abs(pdg2)<=16 )) ;

		  if(verbose ) cout << "     > ID " << (*des2)->pdg_id() << ", status " << (*des2)->status() << ", mass " << (*des2)->momentum().m() << endl;

		  if(isWdaughter && countdau==0){
		    thetop->wdau1mass   = (*des2)->momentum().m()>0 ? (*des2)->momentum().m() : 0.;
		    thetop->wdau1pt     = (*des2)->momentum().perp();
		    thetop->wdau1eta    = (*des2)->momentum().eta();
		    thetop->wdau1phi    = (*des2)->momentum().phi();
		    thetop->wdau1id     = (*des2)->pdg_id();
		  }
		  else if(isWdaughter && countdau==1){
		    thetop->wdau2mass   = (*des2)->momentum().m()>0 ? (*des2)->momentum().m() : 0.;
		    thetop->wdau2pt     = (*des2)->momentum().perp();
		    thetop->wdau2eta    = (*des2)->momentum().eta();
		    thetop->wdau2phi    = (*des2)->momentum().phi();
		    thetop->wdau2id     = (*des2)->pdg_id();
		  }
		  else{}
		  countdau++;		
		}
	      }
	    } // des
	  }
	  else{
	    thetop->bmass     = 0;
	    thetop->bpt       = 0;
	    thetop->beta      = 0;
	    thetop->bphi      = 0;
	    thetop->bstatus   = -99;
	    thetop->wdau1mass = ((*t)->momentum()).m();
	    thetop->wdau1pt   = ((*t)->momentum()).perp();
	    thetop->wdau1eta  = ((*t)->momentum()).eta();
	    thetop->wdau1phi  = ((*t)->momentum()).phi();
	    thetop->wdau1id   = -99;
	    thetop->wdau2mass = 0;
	    thetop->wdau2pt   = 0;
	    thetop->wdau2eta  = 0;
	    thetop->wdau2phi  = 0;
	    thetop->wdau2id   = -99;  
	  }

	}
      

	// fill the tree
	pcount++;
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
  cout << string(Form("Filter efficiency: %.2f %%", float(pcount)/float(icount)*100.)) << endl;

  // delete MEIntegrator and other allocated objects
  cout << "Delete clocks and random engine..." << endl;
  delete clock; delete clock2;
  if(ran!=0)    delete ran;
  cout << "Finished!!!" << endl;
  cout << "*******************" << endl;
  
  // Done!!!
  return 0;
}
