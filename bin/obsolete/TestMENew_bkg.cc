#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include <cstdlib>
#include <iostream> 
#include <fstream>
#include <map>
#include <string>

#include "TSystem.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraph.h"
#include "Math/GenVector/LorentzVector.h"
#include "TLorentzVector.h"
#include "TVectorD.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TMatrixT.h"

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GSLMCIntegrator.h"
#include "Math/AllIntegrationTypes.h"

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


#include "RooWorkspace.h"
#include "RooAbsData.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooCBShape.h"
#include "RooExponential.h"
#include "RooLandau.h"
#include "RooUniform.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooFitResult.h"
#include "RooKeysPdf.h"
#include "RooBifurGauss.h"
#include "RooVoigtian.h"
#include "RooSimultaneous.h"
#include "RooCategory.h"
#include "RooGaussian.h"
#include "RooAbsCategory.h"
#include "RooNumConvPdf.h"
#include "RooFFTConvPdf.h"
#include "RooTruthModel.h"
#include "RooDataHist.h"
#include "RooGenericPdf.h"
#include "RooPolynomial.h"
#include "RooBernstein.h"
#include "RooNDKeysPdf.h"
#include "RooChiSquarePdf.h"

#include "Bianchi/TTHStudies/interface/MEIntegratorNew.h"
#include "Bianchi/TTHStudies/interface/Samples.h"
#include "Bianchi/TTHStudies/interface/CalcME.h"

#define GENJETDR 0.3
#define VERBOSE  false

using namespace std;
using namespace RooFit;

typedef struct 
{
  float et; 
  float sumet;
  float sig;
  float phi;
} metInfo;

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




bool appendEventLHCO(int counter, 
		     ofstream& out, 
		     TLorentzVector leptonLV,
		     int charge,
		     TLorentzVector neutrinoLV,
		     vector<TLorentzVector> myJetsFilt, 
		     std::map<unsigned int, JetByPt> mapFilt,
		     float& closestMbb, float& closestMjj, float& closestMbjj
		     ){

  int algo1 = 0;
  int algo2 = 1;

  string line_inc;
  line_inc = "#    typ     eta       phi       pt     jmass     ntrk     btag    had/em   dummy  dummy";
  out << line_inc << endl;

  string line_num;
  line_num = "0    " +  string(Form("%d",(counter-1))) + "          1";
  out << line_num << endl;

  double phi = leptonLV.Phi() + TMath::Pi();
  out << "1    2    " << string(Form("%.3f", leptonLV.Eta())) 
      << "    " 
      << string(Form("%.3f", phi )) 
      << "    " 
      << string(Form("%.3f", leptonLV.Pt())) 
      << "    "
      << string(Form("%.3f", 0. ))
      << "    "
      << charge
      << "    "
      << 0.0 << "    " << 0.0 << "    " << 0.0 << "    " << 0.0
      << endl;


  float mindMass  = 9999.;
  unsigned int w1 = 999;
  unsigned int w2 = 999;
  for( unsigned int j1 = 0; j1 < myJetsFilt.size()-1 ; j1++ ){
    if( !(mapFilt[j1].csv < 0.679 &&  mapFilt[j1].pt > 30  && TMath::Abs(mapFilt[j1].eta)<2.5) ) continue;
    for( unsigned int j2 = j1+1; j2 < myJetsFilt.size() ; j2++ ){
      if( !(mapFilt[j2].csv < 0.679 &&  mapFilt[j2].pt > 30  && TMath::Abs(mapFilt[j2].eta)<2.5) ) continue;
      float mass  = (myJetsFilt[j1] + myJetsFilt[j2]).M();
      float dMass = TMath::Abs( mass - 81.);
      if( dMass < mindMass ){
	mindMass = dMass;
	w1 = j1; w2 = j2;
	closestMjj = mass;
      }
    }
  }

  if(w1==999 || w2==999) return false;

  vector<unsigned int> bjets;

  int counterBtag   = 0;
  int counterBUntag = 0;
  int counterJet    = 0;
  for( unsigned int j = 0; j < myJetsFilt.size() ; j++ ){

    if(counterBtag<4 && mapFilt[j].csv > 0.679 &&  mapFilt[j].pt > 30 && TMath::Abs(mapFilt[j].eta)<2.5 ){

      bjets.push_back( j );

      double phi = myJetsFilt[j].Phi() + TMath::Pi();
      out << string(Form("%d",counterJet+2)) << "    4    " << string(Form("%.3f", myJetsFilt[j].Eta())) 
	  << "    " 
	  << string(Form("%.3f", phi )) 
	  << "    " 
	  << string(Form("%.3f", myJetsFilt[j].Pt())) 
	  << "    "
	  << string(Form("%.3f", myJetsFilt[j].M()))  
	  << "    "
	  << 0.0
	  << "    "
	  << 2.0 << "    " << 0.0 << "    " << 0.0 << "    " << 0.0
	  << endl;

      counterJet++;
      counterBtag++;
    }

    if(algo1){
      if(counterBUntag<2 && mapFilt[j].csv < 0.679 &&  mapFilt[j].pt > 30  && TMath::Abs(mapFilt[j].eta)<2.5){

	double phi = myJetsFilt[j].Phi() + TMath::Pi();
	out << string(Form("%d",counterJet+2)) << "    4    " << string(Form("%.3f", myJetsFilt[j].Eta())) 
	    << "    " 
	    << string(Form("%.3f", phi )) 
	    << "    " 
	    << string(Form("%.3f", myJetsFilt[j].Pt())) 
	    << "    "
	    << string(Form("%.3f", myJetsFilt[j].M()))  
	    << "    "
	    << 0.0
	    << "    "
	    << 0.0 << "    " << 0.0 << "    " << 0.0 << "    " << 0.0
	    << endl;
	
	counterJet++;
	counterBUntag++;
      }
    }
  }

  if(algo2){
    double phi = myJetsFilt[w1].Phi() + TMath::Pi();
    out << string(Form("%d",counterJet+2)) << "    4    " << string(Form("%.3f", myJetsFilt[w1].Eta())) 
	<< "    " 
	<< string(Form("%.3f", phi )) 
	<< "    " 
	<< string(Form("%.3f", myJetsFilt[w1].Pt())) 
	<< "    "
	<< string(Form("%.3f", myJetsFilt[w1].M()))  
	<< "    "
	<< 0.0
	<< "    "
	<< 0.0 << "    " << 0.0 << "    " << 0.0 << "    " << 0.0
	<< endl;
    
    counterJet++;
    counterBUntag++;

    phi = myJetsFilt[w2].Phi() + TMath::Pi();
    out << string(Form("%d",counterJet+2)) << "    4    " << string(Form("%.3f", myJetsFilt[w2].Eta())) 
	<< "    " 
	<< string(Form("%.3f", phi )) 
	<< "    " 
	<< string(Form("%.3f", myJetsFilt[w2].Pt())) 
	<< "    "
	<< string(Form("%.3f", myJetsFilt[w2].M()))  
	<< "    "
	<< 0.0
	<< "    "
	<< 0.0 << "    " << 0.0 << "    " << 0.0 << "    " << 0.0
	<< endl;
    
    counterJet++;
    counterBUntag++;
  }

  phi = neutrinoLV.Phi() + TMath::Pi();
  out <<  "8    6    " 
      << string(Form("%.3f", 0. )) 
      << "    " 
      << string(Form("%.3f", phi )) 
      << "    " 
      << string(Form("%.3f", neutrinoLV.Pt())) 
      << "    "
      << 0.0 
      << "    "
      << 0.0
      << "    "
      << 0.0 << "    " << 0.0 << "    " << 0.0 << "    " << 0.0
      << endl;

  float minMbb = 999.;
  for(unsigned int bb1 = 0; bb1 != bjets.size()-1 ; bb1++){
    for(unsigned int bb2 = bb1+1; bb2 != bjets.size(); bb2++){
      float dMass = TMath::Abs( 125. - (myJetsFilt[bb1] + myJetsFilt[bb2]).M() );
      if( dMass < minMbb ){
	minMbb = dMass;
	closestMbb = (myJetsFilt[bb1] + myJetsFilt[bb2]).M() ;
      }
    }
  }

  float minMbjj = 999.;
  for(unsigned int bb1 = 0; bb1 != bjets.size() ; bb1++){
    float dMass = TMath::Abs( (myJetsFilt[bb1]+myJetsFilt[w1]+myJetsFilt[w2]).M() -175. );
    if ( dMass < minMbjj ){
      minMbjj = dMass;
      closestMbjj = (myJetsFilt[bb1]+myJetsFilt[w1]+myJetsFilt[w2]).M();
    }
  }


  if( !(counterJet==6 && counterBtag==4 && counterBUntag==2) )
    return false;
  
  return true;
}


int main(int argc, const char* argv[])
{

  std::cout << "TestMENew_bkg" << std::endl;
  gROOT->SetBatch(true);
 
  gSystem->Load("libFWCoreFWLite");
  gSystem->Load("libDataFormatsFWLite");

  AutoLibraryLoader::enable();

  PythonProcessDesc builder(argv[1]);
  const edm::ParameterSet& in = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteInput");

  const edm::VParameterSet& samples = in.getParameter<edm::VParameterSet>("samples") ;
  std::string outFileName  ( in.getParameter<std::string>  ("outFileName" ) );
  std::string pathToFile   ( in.getParameter<std::string>  ("pathToFile" ) );
  std::string ordering     ( in.getParameter<std::string>  ("ordering" ) );
  std::string pathToTF     ( in.getParameter<std::string>  ("pathToTF"   ) );
  std::string madweight    ( in.getParameter<std::string>  ("madweight"   ) );
  int vegasPoints          ( in.getParameter<int>   ("vegasPoints"   ) );
  bool verbose             ( in.getParameter<bool>  ("verbose") );
  double lumi              ( in.getParameter<double>("lumi") );

  float met       ( in.getParameter<double> ("met")        );
  float pertW1    ( in.getParameter<double> ("pertW1")     );
  float pertW2    ( in.getParameter<double> ("pertW2")     );
  float pertBHad  ( in.getParameter<double> ("pertBHad")   );
  float pertBLep  ( in.getParameter<double> ("pertBLep")   );
  float enlargeE1 ( in.getParameter<double> ("enlargeE1")  );
  float enlargeEh1( in.getParameter<double> ("enlargeEh1") );
  float enlargePt ( in.getParameter<double> ("enlargePt")  );

  float scaleH  ( in.getParameter<double> ("scaleH")   );
  float scaleL  ( in.getParameter<double> ("scaleL")   );
  float scaleMET( in.getParameter<double> ("scaleMET") );

  int   useME     ( in.getParameter<int>    ("useME")      );
  int   useJac    ( in.getParameter<int>    ("useJac")     );
  int   useMET    ( in.getParameter<int>    ("useMET")     );
  int   useTF     ( in.getParameter<int>    ("useTF")      );
  int   usePDF    ( in.getParameter<int>    ("usePDF")     );
  int   shiftMomenta    ( in.getParameter<int>    ("shiftMomenta")     );
  int   testMassScan    ( in.getParameter<int>    ("testMassScan")     );
  int   testPermutations( in.getParameter<int>    ("testPermutations") );
  int   compAcceptance  ( in.getUntrackedParameter<int>    ("compAcceptance", 0) );
  int   printP4         ( in.getParameter<int>    ("printP4")          );

  vector<int> evLimits( in.getParameter<vector<int> >  ("evLimits") );
  int evLow = evLimits[0];
  int evHigh = evLimits[1];

  int   mode         ( in.getUntrackedParameter<int>    ("mode", 0)          );

  ofstream out("./root/input.lhco");
  out.precision(8);

  gSystem->Exec(("rm "+outFileName).c_str());

  TStopwatch* clock = new TStopwatch();
  TRandom3*   ran   = new TRandom3();

  int par = mode==0 ? 4 : 6;
  MEIntegratorNew* meIntegrator = new MEIntegratorNew( pathToTF , par, int(verbose));
  if(mode == 0)
    meIntegrator->setIntType( MEIntegratorNew::SL2wj );
  else if(mode == 1)
    meIntegrator->setIntType( MEIntegratorNew::SL1wj );
  else{
    cout << "Unsupported mode... exit" << endl;
    return 0;
  }
  meIntegrator->setWeightNorm( 1 ); // divide weights by xsec

  //const int nMassPoints = 30;
  //double mH[nMassPoints] = {95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240};
  //const int nMassPoints  = 18;                                                                                                                        
  //double mH[nMassPoints] = {90, 95, 100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175}; 
 

  // 21 bins 80 - 180
  const int nMassPoints  = 39;//31;                                                                                                                        
  double mH[nMassPoints] = 
    { /*50,  55,  */60,  65,  70,  75,  
      80 , 85,  90, 95, 100, 105, 110, 
      115, 120, 125, 130, 135, 140, 145, 150, 155, 
      160, 165, 170, 175, 180 ,185, 190, 195, 200, 205, 210, 215, 220, 225, 230, 235, 240, 245, 250 
			     // 235, 240, 245, 250 
    }; 


  TH1F*  hMass         = new TH1F("hMass","",nMassPoints, mH[0]-2.5, mH[nMassPoints-1]+2.5);
  TH1F*  hBestMass     = new TH1F("hBestMass","",40, 50, 250);
  TH1F*  hMassProb     = new TH1F("hMassProb","",nMassPoints, mH[0]-2.5, mH[nMassPoints-1]+2.5);
  TH1F*  hBestMassProb = new TH1F("hBestMassProb","",41, 47.5, 252.5);

  TTree* tProb         = new TTree("tree","");
  float prob_;
  float mass_;
  float initCpu_,initRea_,instCpu_, instRea_,evalCpu_,evalRea_;
  tProb->Branch("prob",&prob_,"prob/F");
  tProb->Branch("mass",&mass_,"mass/F");
  tProb->Branch("initCpu",&initCpu_,"initCpu/F");
  tProb->Branch("initRea",&initRea_,"initRea/F");
  tProb->Branch("instCpu",&instCpu_,"instCpu/F");
  tProb->Branch("instRea",&instRea_,"instRea/F");
  tProb->Branch("evalCpu",&evalCpu_,"evalCpu/F");
  tProb->Branch("evalRea",&evalRea_,"evalRea/F");

  TTree* tMW           = new TTree("tree_madweight","");
  int evt_;
  float weight_sgn_, weightErr_sgn_;
  float weight_bkg_, weightErr_bkg_;
  float closestMbb_, closestMjj_, closestMbjj_;
  float hPt_;

  tMW->Branch("evt",          &evt_,            "evt/I");
  tMW->Branch("weight_sgn",   &weight_sgn_,   "weight_sgn/F");
  tMW->Branch("weightErr_sgn",&weightErr_sgn_,"weightErr_sgn/F");
  tMW->Branch("weight_bkg",   &weight_bkg_,   "weight_bkg/F");
  tMW->Branch("weightErr_bkg",&weightErr_bkg_,"weightErr_bkg/F");
  tMW->Branch("hPt",          &hPt_,          "hPt/F");

  tMW->Branch("closestMbb", &closestMbb_,"closestMbb/F");
  tMW->Branch("closestMjj", &closestMjj_,"closestMjj/F");
  tMW->Branch("closestMbjj",&closestMbjj_,"closestMbjj/F");

  evt_           = -99;
  weight_sgn_    = -99.;
  weightErr_sgn_ = -99.;
  weight_bkg_    = -99.;
  weightErr_bkg_ = -99.;
  closestMbb_    = -99.;
  closestMjj_    = -99.; 
  closestMbjj_   = -99.;
  hPt_           = -99.;
  ////////////////////////////////////////////////////////


  bool openAllFiles  = false;
  Samples* mySamples = new Samples(openAllFiles, pathToFile, ordering, samples, lumi, verbose);
  vector<string> mySampleFiles;

  if(mySamples->IsOk()){

    cout << "Ok!" << endl;
    mySampleFiles = mySamples->Files();

    for( unsigned int i = 0 ; i < mySampleFiles.size(); i++){
      string sampleName       = mySampleFiles[i];

      if(verbose){
	cout << mySampleFiles[i] << " ==> " << mySamples->GetXSec(sampleName) 
	     << " pb,"
	     << " ==> weight = "            << mySamples->GetWeight(sampleName) << endl;
      }
    }
  }
  else{
    cout << "Problems... leaving" << endl;
    return 0;
  }


  ////////////////////////////////////////////////////////


  for(unsigned int sample = 0 ; sample < mySampleFiles.size(); sample++){
    
    string currentName       = mySampleFiles[sample];

    mySamples->OpenFile( currentName );
    cout << "Opening file " << currentName << endl;
    TTree* currentTree       = mySamples->GetTree( currentName, "tree");
    cout << "Done!!" << endl;
    
    //genParticleInfo genB, genBbar;
    //genTopInfo genTop, genTbar;
    JetByPt jet1, jet2, jet3, jet4, jet5, jet6, jet7, jet8, jet9, jet10;
    int nvlep, nSimBs;
    Float_t vLepton_mass[999];
    Float_t vLepton_pt  [999];
    Float_t vLepton_eta [999];
    Float_t vLepton_phi [999];
    Float_t vLepton_charge [999];
    metInfo METtype1p2corr;
    genParticleInfo genB, genBbar;

    currentTree->SetBranchAddress("jet1",   &jet1);
    currentTree->SetBranchAddress("jet2",   &jet2);
    currentTree->SetBranchAddress("jet3",   &jet3);
    currentTree->SetBranchAddress("jet4",   &jet4);
    currentTree->SetBranchAddress("jet5",   &jet5);
    currentTree->SetBranchAddress("jet6",   &jet6);
    currentTree->SetBranchAddress("jet7",   &jet7);
    currentTree->SetBranchAddress("jet8",   &jet8);
    currentTree->SetBranchAddress("jet9",   &jet9);
    currentTree->SetBranchAddress("jet10",  &jet10);
    currentTree->SetBranchAddress("vLepton_mass"  ,vLepton_mass);
    currentTree->SetBranchAddress("vLepton_pt"    ,vLepton_pt);
    currentTree->SetBranchAddress("vLepton_eta"   ,vLepton_eta);
    currentTree->SetBranchAddress("vLepton_phi"   ,vLepton_phi);
    currentTree->SetBranchAddress("vLepton_charge",vLepton_charge);
    currentTree->SetBranchAddress("METtype1p2corr",&METtype1p2corr);
    currentTree->SetBranchAddress("nvlep",  &nvlep);
    currentTree->SetBranchAddress("genB",   &genB);
    currentTree->SetBranchAddress("genBbar",&genBbar);
    currentTree->SetBranchAddress("nSimBs", &nSimBs);


    int counter    = 0;
    int counterTot = 0;
    Long64_t nentries = currentTree->GetEntries();
    for (Long64_t i = 0; i < nentries ; i++){

      counterTot++;

      if( (counter+1)>evHigh ) continue; // skip remaining events
      
      if(i%5000==0) cout << i << endl;
      currentTree->GetEntry(i);
      
      if(compAcceptance) counterTot -= ((nSimBs>2) ? 0 : 1);
      //continue;
      
      TLorentzVector genBLV   ;
      TLorentzVector genBbarLV;
      
      if(genB.mass>0 && genB.momid==25){
	genBLV.SetPtEtaPhiM(genB.pt,genB.eta ,genB.phi, genB.mass );
      }
      if(genBbar.mass>0 && genBbar.momid==25){
	genBbarLV.SetPtEtaPhiM(genBbar.pt,genBbar.eta ,genBbar.phi, genBbar.mass );
      }

      vector<TLorentzVector> myJets;
      std::map<unsigned int, JetByPt> map;
      vector<TLorentzVector> myJetsFilt;
      std::map<unsigned int, JetByPt> mapFilt;
      
      map[0] = jet1;  map[1] = jet2; map[2] = jet3; map[3] = jet4; 
      map[4] = jet5;  map[5] = jet6; map[6] = jet7; map[7] = jet8; 
      map[8] = jet9;  map[9] = jet10; 
    

      for(unsigned int m = 0; m <10; m++){
	TLorentzVector jetLV; 
	if(map[m].pt>20){
	  jetLV.SetPtEtaPhiM( map[m].pt, map[m].eta, map[m].phi, map[m].mass );
	  myJets.push_back(jetLV);
	}
      }

      unsigned int iter = 0;    
      for(unsigned int k = 0; k<myJets.size() ; k++){
	
	if(myJets[k].Pt()>30 && TMath::Abs(myJets[k].Eta()) < 2.5 ){
	  
	  myJetsFilt.push_back( myJets[k] );
	  
	  switch(k){
	  case 0:
	    mapFilt[iter] = jet1;
	    break;
	  case 1:
	    mapFilt[iter] = jet2;
	    break;
	  case 2:
	    mapFilt[iter] = jet3;
	    break;
	  case 3:
	    mapFilt[iter] = jet4;
	    break;
	  case 4:
	    mapFilt[iter] = jet5;
	    break;
	  case 5:
	    mapFilt[iter] = jet6;
	    break;
	  case 6:
	    mapFilt[iter] = jet7;
	    break;
	  case 7:
	    mapFilt[iter] = jet8;
	    break;
	  case 8:
	    mapFilt[iter] = jet9;
	    break;
	  case 9:
	    mapFilt[iter] = jet10;
	    break;
	  default:
	    break;
	  }
	  
	  iter++;
	}
      }
      


      int numJets20 =0;     
      int numJets30 =0; 
      int numJets60 =0; 
      int numJets20BtagL=0; 
      int numJets20BtagM=0; 
      int numJets20BtagT=0;
      int numJets30BtagL=0; 
      int numJets30BtagM=0; 
      int numJets30BtagT=0;
      int numJets60BtagL=0; 
      int numJets60BtagM=0; 
      int numJets60BtagT=0;

      for(unsigned int k = 0; k < myJetsFilt.size(); k++){  
	float csv_k = mapFilt[k].csv>0 ? mapFilt[k].csv : 0.0 ;
	float pt_k  = mapFilt[k].pt>0  ? mapFilt[k].pt  : 0.0 ;
	
	if( pt_k>20 && csv_k>0.244 ) numJets20BtagL++;
	if( pt_k>20 && csv_k>0.679 ) numJets20BtagM++;
	if( pt_k>20 && csv_k>0.898 ) numJets20BtagT++;
	if( pt_k>30 && csv_k>0.244 ) numJets30BtagL++;
	if( pt_k>30 && csv_k>0.679 ) numJets30BtagM++;
	if( pt_k>30 && csv_k>0.898 ) numJets30BtagT++;
	if( pt_k>60 && csv_k>0.244 ) numJets60BtagL++;
	if( pt_k>60 && csv_k>0.679 ) numJets60BtagM++;
	if( pt_k>60 && csv_k>0.898 ) numJets60BtagT++;
	if( pt_k>20     ) numJets20++;
	if( pt_k>30     ) numJets30++;
	if( pt_k>60     ) numJets60++;
       
      }
      
      if( myJetsFilt.size()<6 || numJets30<6 || numJets30BtagM<4) continue;   // RELAXED
      bool properEvent = true;

      properEvent = nvlep==1;
      if(!properEvent) continue;

      TLorentzVector leptonLV;
      leptonLV.SetPtEtaPhiM(vLepton_pt[0],vLepton_eta[0],vLepton_phi[0],vLepton_mass[0]);
      properEvent = leptonLV.Pt()>20;
     
      if(!properEvent) continue;
      counter++;

      if(compAcceptance){
	counter    -= ((nSimBs>2) ? 0 : 1);
	continue;
      }

      float nuPx = METtype1p2corr.et*TMath::Cos(METtype1p2corr.phi);
      float nuPy = METtype1p2corr.et*TMath::Sin(METtype1p2corr.phi);
      float nuE  = TMath::Sqrt(nuPx*nuPx+nuPy*nuPy);
      TLorentzVector neutrinoLV;
      neutrinoLV.SetPxPyPzE(nuPx,nuPy,0. ,nuE);


      if( !(counter>=evLow && counter<=evHigh) ) continue;
      cout << "Processing event # " << counter << endl;

      prob_    = 0.;
      mass_    = -99;
      initCpu_ = -99;
      initRea_ = -99;
      instCpu_ = -99;
      evalCpu_ = -99;
      evalCpu_ = -99;
      evalRea_ = -99;


      /*
      appendEventLHCO( counter, out,  leptonLV, vLepton_charge[0], neutrinoLV,  myJetsFilt, mapFilt, closestMbb_, closestMjj_, closestMbjj_);

      bool match_sgn = false;
      bool match_bkg = false;

      ifstream in_sgn;
      ifstream in_bkg;
      char* c = new char[1000];
      in_sgn.open("./root/"+madweight+"_sgn.out");
      in_bkg.open("./root/"+madweight+"_bkg.out");

      while (in_sgn.good()){
	in_sgn.getline(c,1000,'\n');
	if (in_sgn.good()){	  
	  string line(c);
	  if(line.find("#")!=string::npos) continue;
	  size_t first   = line.find_first_of(" ");
	  size_t second  = line.find_first_of(" ", first+1);
	  size_t third   = line.find_first_of(" ", second+1);

	  string number_str = line.substr(0, first);
	  int number_int    = atoi(number_str.c_str());

	  if( (counter-1) != number_int ) continue; 

	  match_sgn = true;

	  string trig_str   = line.substr(first+1, second);
	  int trig_int      = atoi(trig_str.c_str());

	  string w_str    = line.substr(second+1, third);
	  float  w_float  = atof(w_str.c_str());

	  string we_str    = line.substr(third+1, string::npos);
	  float  we_float  = atof(we_str.c_str());

	  //cout << number_int << " " << trig_int << " " << w_float << " " << we_float << endl;
	  evt_           = number_int;
	  weight_sgn_    = w_float;
	  weightErr_sgn_ = we_float;

	}
      }
      while (in_bkg.good()){
	in_bkg.getline(c,1000,'\n');
	if (in_bkg.good()){	  
	  string line(c);
	  if(line.find("#")!=string::npos) continue;
	  size_t first   = line.find_first_of(" ");
	  size_t second  = line.find_first_of(" ", first+1);
	  size_t third   = line.find_first_of(" ", second+1);

	  string number_str = line.substr(0, first);
	  int number_int    = atoi(number_str.c_str());

	  if( (counter-1) != number_int ) continue; 

	  match_bkg = true;
 
	  string trig_str   = line.substr(first+1, second);
	  int trig_int      = atoi(trig_str.c_str());

	  string w_str    = line.substr(second+1, third);
	  float  w_float  = atof(w_str.c_str());

	  string we_str    = line.substr(third+1, string::npos);
	  float  we_float  = atof(we_str.c_str());

	  //cout << number_int << " " << trig_int << " " << w_float << " " << we_float << endl;
	  evt_           = number_int;
	  weight_bkg_    = w_float;
	  weightErr_bkg_ = we_float;
	}
      }
      hPt_ =  (genBLV+genBbarLV).Pt();
      if(match_sgn && match_bkg)  tMW->Fill();
      delete c;
      */





      //////////////////////////////////////////////////////////
      // test permutations
      if(testPermutations) {


	hMassProb->Reset();

	for(int m = 0; m < nMassPoints ; m++){
	  cout << "Processing " << mH[m] << endl;
	  meIntegrator->setMass( mH[m] );

	  for(unsigned int pp1=0; pp1!=myJetsFilt.size(); pp1++){

	    for(unsigned int pp2=0; pp2!=myJetsFilt.size(); pp2++){
	      if(pp2==pp1) continue;
	    
	      for(unsigned int pp3=0; pp3!=myJetsFilt.size()-1; pp3++){
		if(pp3==pp1 || pp3==pp2) continue;
		
		for(unsigned int pp4=pp3+1; pp4!=myJetsFilt.size(); pp4++){
		  if(pp4==pp1 || pp4==pp2 || pp4==pp3) continue;

		  for(unsigned int pp5=0; pp5!=myJetsFilt.size()-1; pp5++){
		    if(pp5==pp1 || pp5==pp2 || pp5==pp3 || pp5==pp4) continue;

		    for(unsigned int pp6=pp5+1; pp6!=myJetsFilt.size(); pp6++){
		      if(pp6==pp1 || pp6==pp2 || pp6==pp3 || pp6==pp4 || pp6==pp5) continue;


		      if( !(  mapFilt[pp1].csv > /*0.244*/0.679 &&
			      mapFilt[pp2].csv > /*0.244*/0.679 &&
			      mapFilt[pp5].csv > /*0.244*/0.679 &&
			      mapFilt[pp6].csv > /*0.244*/0.679 &&
			      mapFilt[pp3].csv < 999   &&
			      mapFilt[pp4].csv < 999 
			      ) ) continue;

		      float massPair = (myJetsFilt[pp5]+myJetsFilt[pp6]).M(); 
		      if( TMath::Abs( massPair - mH[m]) > 2.0*scaleH*massPair ){
			continue; 
		      }

		      
		      
		      vector<TLorentzVector> jets;
		      jets.push_back(  leptonLV   );
		      jets.push_back(  neutrinoLV );
		      jets.push_back( myJetsFilt[pp1] ); // top lep b
		      jets.push_back( myJetsFilt[pp3] ); // W q
		      jets.push_back( myJetsFilt[pp4] ); // W q'
		      jets.push_back( myJetsFilt[pp2] ); // top had b
		      jets.push_back( myJetsFilt[pp5] ); // h1
		      jets.push_back( myJetsFilt[pp6] ); // h2
		    
		    
		      meIntegrator->setJets(&jets);
		      meIntegrator->initVersors(1);
		      meIntegrator->initTF();
		      meIntegrator->setTopFlags( vLepton_charge[0]==1 ? +1 : -1 , vLepton_charge[0]==1 ? -1 : +1 );
		      meIntegrator->setUseME (useME);
		      meIntegrator->setUseJac(useJac);
		      meIntegrator->setUseMET(useMET);
		      meIntegrator->setUseTF (useTF);
		      meIntegrator->setUsePDF(usePDF);
		      
		      double E1low   = (meIntegrator->getCachedTF("tfWjet1"))->GetXaxis()->GetXmin() * (1-enlargeE1);
		      double E1high  = (meIntegrator->getCachedTF("tfWjet1"))->GetXaxis()->GetXmax() * (1+enlargeE1);
		      double Ptlow   = (meIntegrator->getCachedTF("tfMetPt"))->GetXaxis()->GetXmin() * (1-enlargePt);
		      double Pthigh  = (meIntegrator->getCachedTF("tfMetPt"))->GetXaxis()->GetXmax() * (1+enlargePt);
		      Ptlow  = -1; 
		      Pthigh = +1; 
		      double Philow   = -(meIntegrator->getCachedTF("tfMetPhi"))->GetYaxis()->GetXmax();
		      double Phihigh  =  (meIntegrator->getCachedTF("tfMetPhi"))->GetYaxis()->GetXmax();
		      double Eh1low   = (meIntegrator->getCachedTF("tfHiggs1"))->GetXaxis()->GetXmin() * (1-enlargeEh1);
		      double Eh1high  = (meIntegrator->getCachedTF("tfHiggs1"))->GetXaxis()->GetXmax() * (1+enlargeEh1);	  
		      
		      double xLmode0[4] = {  E1low,  Ptlow,   Philow,   Eh1low};
		      double xUmode0[4] = {  E1high, Pthigh , Phihigh,  Eh1high};

		      double xLmode1[6] = {  E1low,  -1, -TMath::Pi(), Ptlow,   Philow,   Eh1low};
		      double xUmode1[6] = {  E1high, +1, +TMath::Pi(), Pthigh , Phihigh,  Eh1high};

		      
		      ROOT::Math::Functor toIntegrate(meIntegrator, &MEIntegratorNew::Eval, par);
		      ROOT::Math::GSLMCIntegrator ig2( ROOT::Math::IntegrationMultiDim::kVEGAS , 1.e-12, 1.e-5, vegasPoints);
		      ig2.SetFunction(toIntegrate);
		      meIntegrator->SetPar(par);
		      
		      double p = 0.;
		      if(mode==0)
			p = ig2.Integral(xLmode0, xUmode0);
		      else if(mode==1)
			p = ig2.Integral(xLmode1, xUmode1);
		      else{ }		   
		      
		      float old = hMassProb->GetBinContent( hMassProb->FindBin(mH[m]) );
		      hMassProb->SetBinContent( hMassProb->FindBin(mH[m]), old+p );
		      
		      meIntegrator->deleteTF();	
		    }
		    
		  }
		}
	      }
	    }
	  }
	}// mass points

	float est = -99;
	float a = -99;
	float b = -99;
	float c = -99;

	int maxBin = hMassProb->GetMaximumBin();

	if( maxBin<  hMassProb->GetNbinsX() && maxBin>1 ){

	  double xD =  hMassProb->GetBinCenter (hMassProb->GetMaximumBin()-1);
	  double xC =  hMassProb->GetBinCenter (hMassProb->GetMaximumBin());
	  double xU =  hMassProb->GetBinCenter (hMassProb->GetMaximumBin()+1);
	  double yD =  hMassProb->GetBinContent(hMassProb->GetMaximumBin()-1);
	  double yC =  hMassProb->GetBinContent(hMassProb->GetMaximumBin());
	  double yU =  hMassProb->GetBinContent(hMassProb->GetMaximumBin()+1);
	  
	  TMatrixD A(3,3);
	  const double elements[9] = 
	    { 1, xD,  xD*xD,
	      1, xC,  xC*xC,
	      1, xU,  xU*xU 
	    };
	  A.SetMatrixArray(elements, "C");
	  //A.Print("");
	  TMatrixD AInv(3,3);
	  double det;
	  AInv = A.Invert(&det);
	  //AInv.Print("");
	  
	  TMatrixD Y(3,1);
	  const double yPos[3] = 
	    { yD, yC, yU
	    };
	  Y.SetMatrixArray(yPos, "C");
	  //Y.Print("");
	  
	  TMatrixD C(3,1);
	  const double dummy[3] = 
	    { 1., 1., 1.
	    };
	  C.SetMatrixArray(dummy,"C");
	  C.Mult(AInv,Y);
	  //C.Print("");
	  
	  a = C(2,0);
	  b = C(1,0);
	  c = C(0,0);
	  //cout << "a = " << a << endl;
	  //cout << "b = " << b << endl;
	  //cout << "c = " << c << endl;
	  //cout << "Estim. = " << -b/2/a << endl; 
	  est = -b/2/a ;
	}
	else if(maxBin== hMassProb->GetNbinsX()){

	  double xD =  hMassProb->GetBinCenter (hMassProb->GetMaximumBin()-2);
	  double xC =  hMassProb->GetBinCenter (hMassProb->GetMaximumBin()-1);
	  double xU =  hMassProb->GetBinCenter (hMassProb->GetMaximumBin());
	  double yD =  hMassProb->GetBinContent(hMassProb->GetMaximumBin()-2);
	  double yC =  hMassProb->GetBinContent(hMassProb->GetMaximumBin()-1);
	  double yU =  hMassProb->GetBinContent(hMassProb->GetMaximumBin());
	  
	  TMatrixD A(3,3);
	  const double elements[9] = 
	    { 1, xD,  xD*xD,
	      1, xC,  xC*xC,
	      1, xU,  xU*xU 
	    };
	  A.SetMatrixArray(elements, "C");
	  //A.Print("");
	  TMatrixD AInv(3,3);
	  double det;
	  AInv = A.Invert(&det);
	  //AInv.Print("");
	  
	  TMatrixD Y(3,1);
	  const double yPos[3] = 
	    { yD, yC, yU
	    };
	  Y.SetMatrixArray(yPos, "C");
	  //Y.Print("");
	  
	  TMatrixD C(3,1);
	  const double dummy[3] = 
	    { 1., 1., 1.
	    };
	  C.SetMatrixArray(dummy,"C");
	  C.Mult(AInv,Y);
	  //C.Print("");
	  
	  float a = C(2,0);
	  float b = C(1,0);
	  float c = C(0,0);
	  //cout << "a = " << a << endl;
	  //cout << "b = " << b << endl;
	  //cout << "c = " << c << endl;
	  //cout << "Estim. = " << -b/2/a << endl; 
	  est = -b/2/a ;	  
	}
	else if(maxBin==1){

	  double xD =  hMassProb->GetBinCenter (hMassProb->GetMaximumBin());
	  double xC =  hMassProb->GetBinCenter (hMassProb->GetMaximumBin()+1);
	  double xU =  hMassProb->GetBinCenter (hMassProb->GetMaximumBin()+2);
	  double yD =  hMassProb->GetBinContent(hMassProb->GetMaximumBin());
	  double yC =  hMassProb->GetBinContent(hMassProb->GetMaximumBin()+1);
	  double yU =  hMassProb->GetBinContent(hMassProb->GetMaximumBin()+2);
	  
	  TMatrixD A(3,3);
	  const double elements[9] = 
	    { 1, xD,  xD*xD,
	      1, xC,  xC*xC,
	      1, xU,  xU*xU 
	    };
	  A.SetMatrixArray(elements, "C");
	  //A.Print("");
	  TMatrixD AInv(3,3);
	  double det;
	  AInv = A.Invert(&det);
	  //AInv.Print("");
	  
	  TMatrixD Y(3,1);
	  const double yPos[3] = 
	    { yD, yC, yU
	    };
	  Y.SetMatrixArray(yPos, "C");
	  //Y.Print("");
	  
	  TMatrixD C(3,1);
	  const double dummy[3] = 
	    { 1., 1., 1.
	    };
	  C.SetMatrixArray(dummy,"C");
	  C.Mult(AInv,Y);
	  //C.Print("");
	  
	  a = C(2,0);
	  b = C(1,0);
	  c = C(0,0);
	  //cout << "a = " << a << endl;
	  //cout << "b = " << b << endl;
	  //cout << "c = " << c << endl;
	  //cout << "Estim. = " << -b/2/a << endl; 
	  est = -b/2/a ;
	}
	else{
	  est =  hMassProb->GetBinCenter (hMassProb->GetMaximumBin());
	}



	hBestMassProb->Fill( est );
	mass_   = est;
	prob_   = a*est*est + b*est + c;
	tProb->Fill();

	TFile* fout_tmp = TFile::Open(outFileName.c_str(),"UPDATE");
	TDirectory *dir = fout_tmp->mkdir(Form("Event_%d", counter));
	if(dir)
	  dir->cd();
	else
	  fout_tmp->cd(Form("Event_%d", counter));
	hMassProb->Write("hMassProb",TObject::kOverwrite);
	fout_tmp->Close();
	
      }
      //////////////////////////////////////////////////////////
    
    }
   
    cout << "Sample " << currentName << ": acceptance " << counter << "/" << counterTot << " = " << float(counter)/counterTot << endl;
 
  }
  

  TFile* fout_tmp = TFile::Open(outFileName.c_str(),"UPDATE");
  hBestMass->Write("hBestMass",TObject::kOverwrite);
  hBestMassProb->Write("hBestMassProb",TObject::kOverwrite);
  tProb->Write("",TObject::kOverwrite );
  tMW->Write("",TObject::kOverwrite );
  fout_tmp->Close();

  
  cout << "Delete meIntegrator..." << endl;
  delete meIntegrator;
  delete clock; delete ran;
  cout << "Finished!!!" << endl;
  return 0;
  

}
