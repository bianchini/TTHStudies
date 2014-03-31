#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include <iostream> 

#include "TSystem.h"
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TGraph.h"
#include "Math/GenVector/LorentzVector.h"
#include "TLorentzVector.h"

#include "TopQuarkAnalysis/TopKinFitter/interface/TtSemiLepKinFitter.h"
#include "TopQuarkAnalysis/TopKinFitter/plugins/TtSemiLepKinFitProducer.h"

using namespace std;

#define GENJETDR 0.5

typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > LV;
typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LV2;
   
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
  float et; 
  float sumet;
  float sig;
  float phi;
} metInfo;

Float_t deltaR(LV j1, LV j2){
  return TMath::Sqrt(TMath::Power(j1.eta()-j2.eta(),2) + TMath::Power(j1.phi()-j2.phi(),2));
}

Int_t match(LV j1, LV j2){
  return deltaR(j1,j2) < GENJETDR;
}

// -------------------------------------------------------------------------
// -------------------------------------------------------------------------

void fit(Int_t nbEventsMax = 10000000, 
	 Int_t nbEventsPrint = 0, 
	 Int_t doOnlyRightCombination = 0,
	 Int_t printOnlyRightCombination = 0,
	 Int_t maxNbIter = 200,
	 Int_t nbBins = 100)
{

  TH1D* hPRightCombinationAll = new TH1D("hPRightCombinationAll",";;",nbBins,0,1);
  TH1D* hPWrongCombinationAll = new TH1D("hPWrongCombinationAll",";;",nbBins,0,1);
  TH1D* hLogPRightCombinationAll = new TH1D("hLogPRightCombinationAll",";;",nbBins,-300,10);
  TH1D* hLogPWrongCombinationAll = new TH1D("hLogPWrongCombinationAll",";;",nbBins,-300,10);
  Float_t totalRightAll = 0;
  Float_t totalWrongAll = 0;

  TH1D* hPRightCombinationConverged = new TH1D("hPRightCombinationConverged",";;",nbBins,0,1);
  TH1D* hPWrongCombinationConverged = new TH1D("hPWrongCombinationConverged",";;",nbBins,0,1);
  TH1D* hLogPRightCombinationConverged = new TH1D("hLogPRightCombinationConverged",";;",nbBins,-300,10);
  TH1D* hLogPWrongCombinationConverged = new TH1D("hLogPWrongCombinationConverged",";;",nbBins,-300,10);
  Float_t totalRightConverged = 0;
  Float_t totalWrongConverged = 0;

  TH1D* hFracCombConverged = new TH1D("hFracCombConverged",";;",100,0,1);

  TH1D* hFloatingTopMass = new TH1D("hFloatingTopMass",";top mass;",100,0,400);

  // -----------------------------------------------------------------------

  TFile* fBkgd = TFile::Open("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store//user/bianchi/HBB_EDMNtuple/AllHDiJetPt_VType2/DiJetPt_TTJets_SemiLeptMGDecays_8TeV-madgraph-part_VType2.root"); 
  TTree* tBkgd = (TTree*)fBkgd->Get("tree");

  Float_t lumTarget = 10; // fb-1
  Float_t sigmaBkgd = 103.0; // pb 
  Float_t countBkgd = min(((TH1D*)fBkgd->Get("Count"))->Integral(),(Double_t)nbEventsMax);
  Float_t lumBkgd = 0.001*countBkgd/sigmaBkgd; // fb-1
  Float_t weightLumBkgd = lumTarget/lumBkgd;

  Float_t hJet_e[2];
  Float_t aJet_e[99];
  genTopInfo genTop;
  genTopInfo genTbar;
  Int_t index1;
  Int_t index2;
  Int_t index3;
  Int_t index4;
  Int_t index5;
  Int_t index6;
  Float_t pt1;
  Float_t pt2;
  Float_t pt3;
  Float_t pt4;
  Float_t pt5;
  Float_t pt6;
  Float_t eta1;
  Float_t eta2;
  Float_t eta3;
  Float_t eta4;
  Float_t eta5;
  Float_t eta6;
  Float_t phi1;
  Float_t phi2;
  Float_t phi3;
  Float_t phi4;
  Float_t phi5;
  Float_t phi6;
  Int_t numJets30;
  Int_t numJets30bTag;
  Float_t weightTrig2012;
  Float_t PUweight;
  Float_t vLepton_mass[99];
  Float_t vLepton_pt[99];
  Float_t vLepton_eta[99];
  Float_t vLepton_phi[99];
  metInfo METtype1p2corr;

  tBkgd->SetBranchAddress("hJet_e",hJet_e);
  tBkgd->SetBranchAddress("aJet_e",aJet_e);
  tBkgd->SetBranchAddress("genTop",&genTop);
  tBkgd->SetBranchAddress("genTbar",&genTbar);
  tBkgd->SetBranchAddress("index1",&index1);
  tBkgd->SetBranchAddress("index2",&index2);
  tBkgd->SetBranchAddress("index3",&index3);
  tBkgd->SetBranchAddress("index4",&index4);
  tBkgd->SetBranchAddress("index5",&index5);
  tBkgd->SetBranchAddress("index6",&index6);
  tBkgd->SetBranchAddress("pt1",&pt1);
  tBkgd->SetBranchAddress("pt2",&pt2);
  tBkgd->SetBranchAddress("pt3",&pt3);
  tBkgd->SetBranchAddress("pt4",&pt4);
  tBkgd->SetBranchAddress("pt5",&pt5);
  tBkgd->SetBranchAddress("pt6",&pt6);
  tBkgd->SetBranchAddress("eta1",&eta1);
  tBkgd->SetBranchAddress("eta2",&eta2);
  tBkgd->SetBranchAddress("eta3",&eta3);
  tBkgd->SetBranchAddress("eta4",&eta4);
  tBkgd->SetBranchAddress("eta5",&eta5);
  tBkgd->SetBranchAddress("eta6",&eta6);
  tBkgd->SetBranchAddress("phi1",&phi1);
  tBkgd->SetBranchAddress("phi2",&phi2);
  tBkgd->SetBranchAddress("phi3",&phi3);
  tBkgd->SetBranchAddress("phi4",&phi4);
  tBkgd->SetBranchAddress("phi5",&phi5);
  tBkgd->SetBranchAddress("phi6",&phi6);
  tBkgd->SetBranchAddress("numJets30",&numJets30);
  tBkgd->SetBranchAddress("numJets30bTag",&numJets30bTag);
  tBkgd->SetBranchAddress("weightTrig2012",&weightTrig2012);
  tBkgd->SetBranchAddress("PUweight",&PUweight);
  tBkgd->SetBranchAddress("vLepton_mass",vLepton_mass);
  tBkgd->SetBranchAddress("vLepton_pt",vLepton_pt);
  tBkgd->SetBranchAddress("vLepton_eta",vLepton_eta);
  tBkgd->SetBranchAddress("vLepton_phi",vLepton_phi);
  tBkgd->SetBranchAddress("METtype1p2corr",&METtype1p2corr);

  Double_t weight;
  Int_t index[6];
  Double_t pt[6], eta[6], phi[6], energy[6], mass[6];
  Int_t neutrinoSide = 0;
  LV* jetLV[6];
  LV* chargedLeptonLV;
  Double_t nuPx, nuPy, nuE;
  Int_t nbReqJets = 0;
  Int_t nbGen = 0;
  Int_t nbAcc1 = 0;
  Int_t nbAcc2 = 0;
  Int_t nbMatchedJets = 0;
  Int_t nbMatchedLept = 0;
  Int_t countCombinations;
  Int_t isCorrectCombination;
  Int_t status;
  Int_t nbConverged = 0, nbFail = 0;
  Int_t nbConvergedWeighted = 0, nbFailWeighted = 0;
  Int_t countConverged, countFail;
  Double_t p, logp;

  // -----------------------------------------------------------------------

  Int_t nentriesBkgd = (Int_t)tBkgd->GetEntries();

  if (nbEventsMax<nentriesBkgd){
    cout<<endl<<" Warning : only "<<nbEventsMax<<" events out of "<<nentriesBkgd<<" will be considered."<<endl<<endl;
    nentriesBkgd = nbEventsMax;
  }

  for (Int_t i=0; i<nentriesBkgd; i++) {

    tBkgd->GetEntry(i);
    if (i%100==0) cout<<" "<<i<<" / "<<nentriesBkgd<<endl;

    weight = weightLumBkgd*weightTrig2012*PUweight;

    if (!(numJets30>=6 && numJets30bTag>=1)) continue; // (only 3 & 1 by construction)
    nbReqJets++;

    LV topHadrBLV(0.,0.,0.,0.);
    LV topHadrW1LV(0.,0.,0.,0.); 
    LV topHadrW2LV(0.,0.,0.,0.);
    LV topLeptBLV(0.,0.,0.,0.);
    LV topLeptW1LV(0.,0.,0.,0.);
    LV topLeptW2LV(0.,0.,0.,0.);

    if(genTop.bmass>0 && genTbar.bmass>0){

      if(    -5<=genTop.wdau1id  && genTop.wdau1id <=5 
	 &&  -5<=genTop.wdau2id  && genTop.wdau2id <=5
	 && ( (-16<=genTbar.wdau1id && genTbar.wdau1id<=-11)
	   || ( 11<=genTbar.wdau1id && genTbar.wdau1id<=16 ))
	 && ( (-16<=genTbar.wdau2id && genTbar.wdau2id<=-11)
	   || ( 11<=genTbar.wdau2id && genTbar.wdau2id<=16 ))
	){

	topHadrBLV.SetPt(   genTop.bpt );
	topHadrBLV.SetEta(  genTop.beta );
	topHadrBLV.SetPhi(  genTop.bphi );
	topHadrBLV.SetM(    genTop.bmass );
	topHadrW1LV.SetPt(  genTop.wdau1pt );
	topHadrW1LV.SetEta( genTop.wdau1eta );
	topHadrW1LV.SetPhi( genTop.wdau1phi );
	topHadrW1LV.SetM(   genTop.wdau1mass );
	topHadrW2LV.SetPt(  genTop.wdau2pt );
	topHadrW2LV.SetEta( genTop.wdau2eta );
	topHadrW2LV.SetPhi( genTop.wdau2phi );
	topHadrW2LV.SetM(   genTop.wdau2mass );
	topLeptBLV.SetPt(   genTbar.bpt );
	topLeptBLV.SetEta(  genTbar.beta );
	topLeptBLV.SetPhi(  genTbar.bphi );
	topLeptBLV.SetM(    genTbar.bmass );
	topLeptW1LV.SetPt(  genTbar.wdau1pt );
	topLeptW1LV.SetEta( genTbar.wdau1eta );
	topLeptW1LV.SetPhi( genTbar.wdau1phi );
	topLeptW1LV.SetM(   genTbar.wdau1mass );
	topLeptW2LV.SetPt(  genTbar.wdau2pt );
	topLeptW2LV.SetEta( genTbar.wdau2eta );
	topLeptW2LV.SetPhi( genTbar.wdau2phi );
	topLeptW2LV.SetM(   genTbar.wdau2mass );

	if (genTbar.wdau1id==12
	    || genTbar.wdau1id==14
	    || genTbar.wdau1id==16
	    || genTbar.wdau1id==-12
	    || genTbar.wdau1id==-14
	    || genTbar.wdau1id==-16){
	  neutrinoSide = 1;
	} else if (genTbar.wdau2id==12
	    || genTbar.wdau2id==14
	    || genTbar.wdau2id==16
	    || genTbar.wdau2id==-12
	    || genTbar.wdau2id==-14
	    || genTbar.wdau2id==-16){
	  neutrinoSide = 2;
	} else continue;

      } else if(    -5<=genTbar.wdau1id && genTbar.wdau1id<=5 
		&&  -5<=genTbar.wdau2id && genTbar.wdau2id<=5
		&& ( (-16<=genTop.wdau1id && genTop.wdau1id<=-11)
	          || ( 11<=genTop.wdau1id && genTop.wdau1id<=16 ))
	        && ( (-16<=genTop.wdau2id && genTop.wdau2id<=-11)
	          || ( 11<=genTop.wdau2id && genTop.wdau2id<=16 ))
	){
	
	topHadrBLV.SetPt(   genTbar.bpt );
	topHadrBLV.SetEta(  genTbar.beta );
	topHadrBLV.SetPhi(  genTbar.bphi );
	topHadrBLV.SetM(    genTbar.bmass );
	topHadrW1LV.SetPt(  genTbar.wdau1pt );
	topHadrW1LV.SetEta( genTbar.wdau1eta );
	topHadrW1LV.SetPhi( genTbar.wdau1phi );
	topHadrW1LV.SetM(   genTbar.wdau1mass );
	topHadrW2LV.SetPt(  genTbar.wdau2pt );
	topHadrW2LV.SetEta( genTbar.wdau2eta );
	topHadrW2LV.SetPhi( genTbar.wdau2phi );
	topHadrW2LV.SetM(   genTbar.wdau2mass );
	topLeptBLV.SetPt(   genTop.bpt );
	topLeptBLV.SetEta(  genTop.beta );
	topLeptBLV.SetPhi(  genTop.bphi );
	topLeptBLV.SetM(    genTop.bmass );
	topLeptW1LV.SetPt(  genTop.wdau1pt );
	topLeptW1LV.SetEta( genTop.wdau1eta );
	topLeptW1LV.SetPhi( genTop.wdau1phi );
	topLeptW1LV.SetM(   genTop.wdau1mass );
	topLeptW2LV.SetPt(  genTop.wdau2pt );
	topLeptW2LV.SetEta( genTop.wdau2eta );
	topLeptW2LV.SetPhi( genTop.wdau2phi );
	topLeptW2LV.SetM(   genTop.wdau2mass );

	if (genTop.wdau1id==12
	    || genTop.wdau1id==14
	    || genTop.wdau1id==16
	    || genTop.wdau1id==-12
	    || genTop.wdau1id==-14
	    || genTop.wdau1id==-16){
	  neutrinoSide = 1;
	} else if (genTop.wdau2id==12
	    || genTop.wdau2id==14
	    || genTop.wdau2id==16
	    || genTop.wdau2id==-12
	    || genTop.wdau2id==-14
	    || genTop.wdau2id==-16){
	  neutrinoSide = 2;
	} else continue;

      } else continue;
    } else continue;
    nbGen++;

    if (!(abs(topHadrBLV.eta())<5  && topHadrBLV.pt()>30
	  && abs(topHadrW1LV.eta())<5 && topHadrW1LV.pt()>30
	  && abs(topHadrW2LV.eta())<5 && topHadrW2LV.pt()>30
	  && abs(topLeptBLV.eta())<5  && topLeptBLV.pt()>30  
	  )) continue;
    nbAcc1++;

    if (!(abs(topLeptW1LV.eta())<5 && topLeptW1LV.pt()>20
	  && abs(topLeptW2LV.eta())<5 && topLeptW2LV.pt()>20   
	  )) continue;
    nbAcc2++;

    index[0] = index1;
    index[1] = index2;
    index[2] = index3;
    index[3] = index4;
    index[4] = index5;
    index[5] = index6;
    pt[0] = pt1;
    pt[1] = pt2;
    pt[2] = pt3;
    pt[3] = pt4;
    pt[4] = pt5;
    pt[5] = pt6;
    eta[0] = eta1;
    eta[1] = eta2;
    eta[2] = eta3;
    eta[3] = eta4;
    eta[4] = eta5;
    eta[5] = eta6;
    phi[0] = phi1;
    phi[1] = phi2;
    phi[2] = phi3;
    phi[3] = phi4;
    phi[4] = phi5;
    phi[5] = phi6;
    for (Int_t j=0; j<6; j++) {
      energy[j] = (index[j]>=0) ? hJet_e[index[j]] : aJet_e[-index[j]-1];
      mass[j] = TMath::Sqrt(energy[j]*energy[j] - TMath::Power(pt[j]*TMath::CosH(eta[j]),2));
      jetLV[j] = new LV(pt[j],eta[j],phi[j],mass[j]); 
    }
    chargedLeptonLV = new LV(vLepton_pt[0],vLepton_eta[0],vLepton_phi[0],vLepton_mass[0]);
    
    Int_t indexOfJetMatchedTo[4]= {-1,-1,-1,-1};
    Int_t nbJetsMatchedTo[4] = {0,0,0,0};
    // 0 : b/bbar from top -> bjj
    // 1 : 1st quark from top -> bjj
    // 2 : 2nd quark from top -> bjj
    // 3 : b/bbar from top -> blnu
    for (Int_t j=0; j<6; j++) {
      if (match(*jetLV[j],topHadrBLV)
	  && !match(*jetLV[j],topHadrW1LV)
	  && !match(*jetLV[j],topHadrW2LV)
	  && !match(*jetLV[j],topLeptBLV)) {
	nbJetsMatchedTo[0]++;
	indexOfJetMatchedTo[0]=j;	    
      } else if (match(*jetLV[j],topHadrW1LV)
		 && !match(*jetLV[j],topHadrBLV)
		 && !match(*jetLV[j],topHadrW2LV)
		 && !match(*jetLV[j],topLeptBLV)) { 
	nbJetsMatchedTo[1]++;
	indexOfJetMatchedTo[1]=j;	 
      } else if (match(*jetLV[j],topHadrW2LV)
		 && !match(*jetLV[j],topHadrBLV)
		 && !match(*jetLV[j],topHadrW1LV)
		 && !match(*jetLV[j],topLeptBLV)) {
	nbJetsMatchedTo[2]++;
	indexOfJetMatchedTo[2]=j;	
      } else if (match(*jetLV[j],topLeptBLV)
		 && !match(*jetLV[j],topHadrBLV)
		 && !match(*jetLV[j],topHadrW1LV)
		 && !match(*jetLV[j],topHadrW2LV)) { 
	nbJetsMatchedTo[3]++;
	indexOfJetMatchedTo[3]=j;	 
      } 
    }

    if(!(nbJetsMatchedTo[0]==1 
	 && nbJetsMatchedTo[1]==1 
	 && nbJetsMatchedTo[2]==1 
	 && nbJetsMatchedTo[3]==1 
	 )) continue;
    nbMatchedJets++;

    if (neutrinoSide==2){ 
      if(!match(topLeptW1LV,*chargedLeptonLV)) continue;
    } else if (neutrinoSide==1){
      if(!match(topLeptW2LV,*chargedLeptonLV)) continue;
    } else continue;
    nbMatchedLept++;  

    // -----------------------------------------------------------------------
    
    std::vector< TtSemiLepKinFitter::Constraint > constraints; 
    std::vector< edm::ParameterSet > udscResolutions;
    std::vector< edm::ParameterSet > bResolutions;
    std::vector< edm::ParameterSet > lepResolutions;
    std::vector< edm::ParameterSet > metResolutions;
    std::vector< double > jetEnergyResolutionScaleFactors;
    std::vector< double > jetEnergyResolutionEtaBinning;
    
    edm::ParameterSet paramSetB_bin1;
    edm::ParameterSet paramSetB_bin2;
    edm::ParameterSet paramSetB_bin3;
    edm::ParameterSet paramSetUdsc_bin1;
    edm::ParameterSet paramSetUdsc_bin2;
    edm::ParameterSet paramSetUdsc_bin3;
    edm::ParameterSet paramSetLep_bin1;
    edm::ParameterSet paramSetLep_bin2;
    edm::ParameterSet paramSetLep_bin3;
    edm::ParameterSet paramSetMet;  
    
    paramSetB_bin1   .addParameter("bin",string("abs(eta)>=0.  && abs(eta)<1.0"));
    paramSetB_bin2   .addParameter("bin",string("abs(eta)>=1.0 && abs(eta)<2.0"));
    paramSetB_bin3   .addParameter("bin",string("abs(eta)>=2.0 && abs(eta)<5.0")); // /!\ 2.5
    paramSetUdsc_bin1.addParameter("bin",string("abs(eta)>=0.  && abs(eta)<1.0"));
    paramSetUdsc_bin2.addParameter("bin",string("abs(eta)>=1.0 && abs(eta)<2.0"));
    paramSetUdsc_bin3.addParameter("bin",string("abs(eta)>=2.0 && abs(eta)<5.0")); // /!\ 2.5
    paramSetLep_bin1 .addParameter("bin",string("abs(eta)>=0.  && abs(eta)<1.0"));
    paramSetLep_bin2 .addParameter("bin",string("abs(eta)>=1.0 && abs(eta)<2.0"));
    paramSetLep_bin3 .addParameter("bin",string("abs(eta)>=2.0 && abs(eta)<5.0")); // /!\ 2.5
    
    // et : fits from fitResolutionPt_ttjets.C, or getResolutionPt_ttjets.C (for MET)
    paramSetB_bin1   .addParameter("et",string("sqrt(4.11576e+01 + 1.38054e+00*1.38054e+00*et + 0.00000e+00*0.00000e+00*et*et)"));
    paramSetB_bin2   .addParameter("et",string("sqrt(3.75956e+01 + 1.46691e+00*1.46691e+00*et + 0.00000e+00*0.00000e+00*et*et)"));
    paramSetB_bin3   .addParameter("et",string("sqrt(5.40734e+01 + 1.19942e+00*1.19942e+00*et + 1.36122e-06*1.36122e-06*et*et)"));
    paramSetUdsc_bin1.addParameter("et",string("sqrt(8.51479e+00 + 9.85110e-01*9.85110e-01*et + 5.50086e-02*5.50086e-02*et*et)"));
    paramSetUdsc_bin2.addParameter("et",string("sqrt(9.93662e+00 + 1.13533e+00*1.13533e+00*et + 4.88416e-02*4.88416e-02*et*et)"));
    paramSetUdsc_bin3.addParameter("et",string("sqrt(1.50064e+01 + 9.67855e-01*9.67855e-01*et + 1.15671e-02*1.15671e-02*et*et)"));
    paramSetLep_bin1 .addParameter("et",string("sqrt(1.44169e-04*et*et + 4.15168e-08*et*et*et*et)"));
    paramSetLep_bin2 .addParameter("et",string("sqrt(3.77878e-04*et*et + 5.65095e-08*et*et*et*et)"));
    paramSetLep_bin3 .addParameter("et",string("sqrt(9.02986e-04*et*et + 3.24309e-07*et*et*et*et)"));
    paramSetMet      .addParameter("et",string(Form("sqrt(1.12544e+02 + 3.00905e-01*%.5f)",METtype1p2corr.sumet)));
    //paramSetMet      .addParameter("et",string("20"));
    
    // eta : fits from fitResolutionEta_ttjets.C ; csts from getResolutionEta_ttjets.C
    paramSetB_bin1   .addParameter("eta",string("2.38983e-02 + 2.23730e-03*abs(eta) + (-1.00830e-03)*eta*eta + 1.04913e-03*abs(eta*eta*eta)"));
    paramSetB_bin2   .addParameter("eta",string("2.38983e-02 + 2.23730e-03*abs(eta) + (-1.00830e-03)*eta*eta + 1.04913e-03*abs(eta*eta*eta)"));
    paramSetB_bin3   .addParameter("eta",string("2.38983e-02 + 2.23730e-03*abs(eta) + (-1.00830e-03)*eta*eta + 1.04913e-03*abs(eta*eta*eta)"));
    paramSetUdsc_bin1.addParameter("eta",string("2.63342e-02 + 3.48643e-03*abs(eta) + (-2.62017e-03)*eta*eta + 1.45496e-03*abs(eta*eta*eta)"));
    paramSetUdsc_bin2.addParameter("eta",string("2.63342e-02 + 3.48643e-03*abs(eta) + (-2.62017e-03)*eta*eta + 1.45496e-03*abs(eta*eta*eta)"));
    paramSetUdsc_bin3.addParameter("eta",string("2.63342e-02 + 3.48643e-03*abs(eta) + (-2.62017e-03)*eta*eta + 1.45496e-03*abs(eta*eta*eta)"));
    paramSetLep_bin1 .addParameter("eta",string("0.0003"));
    paramSetLep_bin2 .addParameter("eta",string("0.0003"));
    paramSetLep_bin3 .addParameter("eta",string("0.0003"));
    paramSetMet      .addParameter("eta",string("999")); // no prior knowledge 
    
    // phi : csts from getResolutionPhi_ttjets.C ; fit from fitResolutionPhi_ttjets.C
    paramSetB_bin1   .addParameter("phi",string("0.0290914"));
    paramSetB_bin2   .addParameter("phi",string("0.0290914"));
    paramSetB_bin3   .addParameter("phi",string("0.0290914"));
    paramSetUdsc_bin1.addParameter("phi",string("0.0323073"));
    paramSetUdsc_bin2.addParameter("phi",string("0.0323073"));
    paramSetUdsc_bin3.addParameter("phi",string("0.0323073"));
    paramSetLep_bin1 .addParameter("phi",string("0.000148728")); 
    paramSetLep_bin2 .addParameter("phi",string("0.000148728")); 
    paramSetLep_bin3 .addParameter("phi",string("0.000148728"));  
    paramSetMet      .addParameter("phi",string(Form("sqrt(6.54808e+01 + 7.35296e-01*%.5f)/%.5f",METtype1p2corr.sumet,METtype1p2corr.et)));
    //paramSetMet      .addParameter("phi",string("0.746123"));
    
    bResolutions   .push_back(paramSetB_bin1);
    bResolutions   .push_back(paramSetB_bin2);
    bResolutions   .push_back(paramSetB_bin3);
    udscResolutions.push_back(paramSetUdsc_bin1);
    udscResolutions.push_back(paramSetUdsc_bin2);
    udscResolutions.push_back(paramSetUdsc_bin3);
    lepResolutions .push_back(paramSetLep_bin1);
    lepResolutions .push_back(paramSetLep_bin2);
    lepResolutions .push_back(paramSetLep_bin3);
    metResolutions .push_back(paramSetMet);
    
    constraints.push_back(TtSemiLepKinFitter::kWHadMass);
    constraints.push_back(TtSemiLepKinFitter::kWLepMass);
    constraints.push_back(TtSemiLepKinFitter::kTopHadMass);
    constraints.push_back(TtSemiLepKinFitter::kTopLepMass);
    //constraints.push_back(TtSemiLepKinFitter::kEqualTopMasses);
    //constraints.push_back(TtSemiLepKinFitter::kNeutrinoMass); // causes crash
    jetEnergyResolutionScaleFactors.push_back(1.);
    jetEnergyResolutionEtaBinning.push_back(0.);
    jetEnergyResolutionEtaBinning.push_back(-1.);
    
    TtSemiLepKinFitter *fitter = new TtSemiLepKinFitter(TtSemiLepKinFitter::kEtEtaPhi,TtSemiLepKinFitter::kEtEtaPhi,TtSemiLepKinFitter::kEtEtaPhi,maxNbIter,5e-5,1e-4,constraints,80.4,173.,&udscResolutions,&bResolutions,&lepResolutions,&metResolutions,&jetEnergyResolutionScaleFactors,&jetEnergyResolutionEtaBinning);

    // -----------------------------------------------------------------------

    if (i<nbEventsPrint && !doOnlyRightCombination && !printOnlyRightCombination){
      cout<<" event no. "<<i<<"  ----------------------------------------------"<<endl<<endl;
    }

    countCombinations = 0;
    countConverged = 0;
    countFail = 0;

    for (Int_t j0=0; j0<6; j0++){
    for (Int_t j1=0; j1<6; j1++){
    for (Int_t j2=0; j2<6; j2++){
    for (Int_t j3=0; j3<6; j3++){

      if (j0!=j1 && j0!=j2 && j0!=j3 && j1!=j2 && j1!=j3 && j2!=j3 
	  && j1<j2
	  ){  // 180 combinations per event
    
	countCombinations++;

	isCorrectCombination = (
				j0==indexOfJetMatchedTo[0]
				&& j1==indexOfJetMatchedTo[1]
				&& j2==indexOfJetMatchedTo[2]
				&& j3==indexOfJetMatchedTo[3]
				)||(
				    j0==indexOfJetMatchedTo[0]
				    && j1==indexOfJetMatchedTo[2]
				    && j2==indexOfJetMatchedTo[1]
				    && j3==indexOfJetMatchedTo[3]
				    );

	if (doOnlyRightCombination && !isCorrectCombination) continue;

	TLorentzVector p4HadB(jetLV[j0]->Px(),jetLV[j0]->Py(),jetLV[j0]->Pz(),jetLV[j0]->E());
	TLorentzVector p4HadP(jetLV[j1]->Px(),jetLV[j1]->Py(),jetLV[j1]->Pz(),jetLV[j1]->E());
	TLorentzVector p4HadQ(jetLV[j2]->Px(),jetLV[j2]->Py(),jetLV[j2]->Pz(),jetLV[j2]->E());
	TLorentzVector p4LepB(jetLV[j3]->Px(),jetLV[j3]->Py(),jetLV[j3]->Pz(),jetLV[j3]->E());
	TLorentzVector p4Lepton(chargedLeptonLV->Px(),chargedLeptonLV->Py(),chargedLeptonLV->Pz(),chargedLeptonLV->E());
	nuPx = METtype1p2corr.et*TMath::Cos(METtype1p2corr.phi);
	nuPy = METtype1p2corr.et*TMath::Sin(METtype1p2corr.phi);
	nuE = TMath::Sqrt(nuPx*nuPx+nuPy*nuPy);
	TLorentzVector p4Neutrino(nuPx,nuPy,0,nuE);

	status = fitter->fit(p4HadP, p4HadQ, p4HadB, p4LepB, p4Lepton, p4Neutrino, -1, CovarianceMatrix::kMuon);

	pat::Particle resHadP = fitter->fittedHadP();
	pat::Particle resHadQ = fitter->fittedHadQ();
	pat::Particle resHadB = fitter->fittedHadB();
	pat::Particle resLepB = fitter->fittedLepB();
	pat::Particle resLepton = fitter->fittedLepton();
	pat::Particle resNeutrino = fitter->fittedNeutrino();
	TLorentzVector lvHadP(resHadP.px(),resHadP.py(),resHadP.pz(),resHadP.energy());
	TLorentzVector lvHadQ(resHadQ.px(),resHadQ.py(),resHadQ.pz(),resHadQ.energy());
	TLorentzVector lvHadB(resHadB.px(),resHadB.py(),resHadB.pz(),resHadB.energy());
	TLorentzVector lvLepB(resLepB.px(),resLepB.py(),resLepB.pz(),resLepB.energy());
	TLorentzVector lvLepton(resLepton.px(),resLepton.py(),resLepton.pz(),resLepton.energy());
	TLorentzVector lvNeutrino(resNeutrino.px(),resNeutrino.py(),resNeutrino.pz(),resNeutrino.energy());
	TLorentzVector lvHadW   = lvHadP + lvHadQ ;
	TLorentzVector lvHadTop = lvHadW + lvHadB ;
	TLorentzVector lvLepW   = lvLepton + lvNeutrino ;
	TLorentzVector lvLepTop = lvLepW + lvLepB ;

	if (status==0) {
	  countConverged++;
	} else {
	  countFail++;
	}

	p = fitter->fitProb();
    
	if (i<nbEventsPrint && (!printOnlyRightCombination || isCorrectCombination) ){
	  cout<<"  event "<<i<<", ";
	  cout<<"combination "<<countCombinations<<"/180 ("<<j0<<j1<<j2<<j3;
	  cout<<" - "<<(isCorrectCombination?"right":"wrong")<<")"<<endl;
	  cout<<"   |  fit status :    "<<status<<endl;
	  cout<<"   |  (test) :        "<<resHadP.px()<<" "<<resHadQ.px()<<" "<<resHadB.px()<<" "<<resLepB.px()<<" "<<resLepton.px()<<" "<<resNeutrino.px()<<endl;
	  cout<<"   |  hadW mass :     "<<lvHadW.M()<<endl;
	  cout<<"   |  hadTop mass :   "<<lvHadTop.M()<<endl;
	  cout<<"   |  lepW mass :     "<<lvLepW.M()<<endl;
	  cout<<"   |  lepTop mass :   "<<lvLepTop.M()<<endl;
	  cout<<"   |  probability :   "<<p<<endl;
	  cout<<"   |  final nu mass : "<<lvNeutrino.M()<<endl;
	  cout<<"   |  final nu pz :   "<<lvNeutrino.Pz()<<endl;
	  cout<<endl;
	}

	logp = log(p);
	//logp = (p<=0)?-150:log(p);
	
	if(isCorrectCombination){

	  hPRightCombinationAll->Fill(p,weight);
	  hLogPRightCombinationAll->Fill(logp,weight);
	  totalRightAll += weight; 

	  if (status==0) {
	    hPRightCombinationConverged->Fill(p,weight);
	    hLogPRightCombinationConverged->Fill(logp,weight);
	    totalRightConverged += weight; 
	  }

	  hFloatingTopMass->Fill(lvHadTop.M(),weight); // identical with lvLepTop, if constraint kEqualTopMasses is set

	} else {

	  hPWrongCombinationAll->Fill(p,weight);
	  hLogPWrongCombinationAll->Fill(logp,weight);
	  totalWrongAll += weight;

	  if (status==0) {
	    hPWrongCombinationConverged->Fill(p,weight);
	    hLogPWrongCombinationConverged->Fill(logp,weight);
	    totalWrongConverged += weight;
	  }

	}

      }

    }}}}

    hFracCombConverged->Fill((Float_t)countConverged/(countConverged+countFail),weight);
    nbConverged += countConverged;
    nbFail += countFail;
    nbConvergedWeighted += countConverged*weight;
    nbFailWeighted += countFail*weight;

  } 

  // -----------------------------------------------------------------------

  cout<<" impact of the cuts : "<<endl;
  cout<<"  || nentries :      "<<nentriesBkgd<<endl;
  cout<<"  || nbReqJets :     "<<nbReqJets<<endl;
  cout<<"  || nbGen :         "<<nbGen<<endl;
  cout<<"  || nbAcc1 :        "<<nbAcc1<<endl;
  cout<<"  || nbAcc2 :        "<<nbAcc2<<endl;
  cout<<"  || nbMatchedJets : "<<nbMatchedJets<<endl;
  cout<<"  || nbMatchedLept : "<<nbMatchedLept<<endl;
  cout<<endl;

  if (doOnlyRightCombination) {
      cout<<" Warning : only the right combinations were considered."<<endl;
      cout<<endl;
  }

  cout<<" fraction of combinations where fitter converged successfully : "<<endl;
  cout<<"   "<<(Float_t)nbConverged/(nbConverged+nbFail)<<endl;
  cout<<"   or "<<(Float_t)nbConvergedWeighted/(nbConvergedWeighted+nbFailWeighted)<<" (weighted)"<<endl;
  cout<<endl;

  cout<<" hPRightCombinationAll : "<<endl;
  cout<<"   expected :  "<<totalRightAll<<endl;
  cout<<"   integral :  "<<hPRightCombinationAll->Integral()<<endl;
  cout<<"   underflow/overflow : "<<hPRightCombinationAll->GetBinContent(0);
  cout<<" / "<<hPRightCombinationAll->GetBinContent(nbBins+1)<<endl;
  cout<<"   mean : "<<hPRightCombinationAll->GetMean()<<endl;
  cout<<"   RMS :  "<<hPRightCombinationAll->GetRMS()<<endl;
  cout<<" hPWrongCombinationAll : "<<endl;
  cout<<"   expected :  "<<totalWrongAll<<endl;
  cout<<"   integral :  "<<hPWrongCombinationAll->Integral()<<endl;
  cout<<"   underflow/overflow : "<<hPWrongCombinationAll->GetBinContent(0);
  cout<<" / "<<hPWrongCombinationAll->GetBinContent(nbBins+1)<<endl;
  cout<<"   mean : "<<hPWrongCombinationAll->GetMean()<<endl;
  cout<<"   RMS :  "<<hPWrongCombinationAll->GetRMS()<<endl;
  if (!doOnlyRightCombination) {
    cout<<" <P_right> / <P_wrong> = "<<hPRightCombinationAll->GetMean()/hPWrongCombinationAll->GetMean()<<endl;
  }
  cout<<endl;

  cout<<" hPRightCombinationConverged : "<<endl;
  cout<<"   expected :  "<<totalRightConverged<<endl;
  cout<<"   integral :  "<<hPRightCombinationConverged->Integral()<<endl;
  cout<<"   underflow/overflow : "<<hPRightCombinationConverged->GetBinContent(0);
  cout<<" / "<<hPRightCombinationConverged->GetBinContent(nbBins+1)<<endl;
  cout<<"   mean : "<<hPRightCombinationConverged->GetMean()<<endl;
  cout<<"   RMS :  "<<hPRightCombinationConverged->GetRMS()<<endl;
  cout<<" hPWrongCombinationConverged : "<<endl;
  cout<<"   expected :  "<<totalWrongConverged<<endl;
  cout<<"   integral :  "<<hPWrongCombinationConverged->Integral()<<endl;
  cout<<"   underflow/overflow : "<<hPWrongCombinationConverged->GetBinContent(0);
  cout<<" / "<<hPWrongCombinationConverged->GetBinContent(nbBins+1)<<endl;
  cout<<"   mean : "<<hPWrongCombinationConverged->GetMean()<<endl;
  cout<<"   RMS :  "<<hPWrongCombinationConverged->GetRMS()<<endl;
  if (!doOnlyRightCombination) {
    cout<<" <P_right> / <P_wrong> = "<<hPRightCombinationConverged->GetMean()/hPWrongCombinationConverged->GetMean()<<endl;
  }
  cout<<endl;

  Float_t vRightAll[nbBins];
  Float_t vWrongAll[nbBins];
  Float_t intRightAll = hPRightCombinationAll->Integral(1,nbBins);
  Float_t intWrongAll = hPWrongCombinationAll->Integral(1,nbBins);
  for (Int_t i=1; i<=nbBins; i++){
    vRightAll[i-1] = hPRightCombinationAll->Integral(i,nbBins)/intRightAll;
    vWrongAll[i-1] = hPWrongCombinationAll->Integral(i,nbBins)/intWrongAll;
  }
  TGraph* gRocAll = new TGraph(nbBins,vRightAll,vWrongAll); 
  gRocAll->GetXaxis()->SetTitle("right combination efficiency");
  gRocAll->GetYaxis()->SetTitle("wrong combination efficiency");
  gRocAll->GetXaxis()->SetLimits(0.,1.);
  gRocAll->SetMinimum(0);
  gRocAll->SetMaximum(1);

  Float_t vRightConverged[nbBins];
  Float_t vWrongConverged[nbBins];
  Float_t intRightConverged = hPRightCombinationConverged->Integral(1,nbBins);
  Float_t intWrongConverged = hPWrongCombinationConverged->Integral(1,nbBins);
  for (Int_t i=1; i<=nbBins; i++){
    vRightConverged[i-1] = hPRightCombinationConverged->Integral(i,nbBins)/intRightConverged;
    vWrongConverged[i-1] = hPWrongCombinationConverged->Integral(i,nbBins)/intWrongConverged;
  }
  TGraph* gRocConverged = new TGraph(nbBins,vRightConverged,vWrongConverged); 
  gRocConverged->GetXaxis()->SetTitle("right combination efficiency");
  gRocConverged->GetYaxis()->SetTitle("wrong combination efficiency");
  gRocConverged->GetXaxis()->SetLimits(0.,1.);
  gRocConverged->SetMinimum(0);
  gRocConverged->SetMaximum(1);

  TFile fOut("compareProb_ttjets.root","recreate");
  hPRightCombinationAll->Write();
  hPWrongCombinationAll->Write();
  hLogPRightCombinationAll->Write();
  hLogPWrongCombinationAll->Write();
  gRocAll->Write("gRocAll");
  hPRightCombinationConverged->Write();
  hPWrongCombinationConverged->Write();
  hLogPRightCombinationConverged->Write();
  hLogPWrongCombinationConverged->Write();
  gRocConverged->Write("gRocConverged");
  if (!doOnlyRightCombination) hFracCombConverged->Write();
  hFloatingTopMass->Write();
  fOut.Close();

}



// -------------------------------------------------------------------------
// -------------------------------------------------------------------------

int main(int argc, const char* argv[])
{

  gROOT->SetBatch(true);
 
  gSystem->Load("libFWCoreFWLite");
  AutoLibraryLoader::enable();

  fit(
      argc>=2 ? atoi(argv[1]) : 10000000 ,
      argc>=3 ? atoi(argv[2]) : 0 ,
      argc>=4 ? atoi(argv[3]) : 0 ,
      argc>=5 ? atoi(argv[4]) : 0 ,
      argc>=6 ? atoi(argv[5]) : 200
      );

  return 0;

}
