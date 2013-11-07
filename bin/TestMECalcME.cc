#include "FWCore/FWLite/interface/AutoLibraryLoader.h"

#include <cstdlib>
#include <iostream> 
#include <fstream>
#include <map>
#include <string>

#include "TInterpreter.h"
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

#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/GSLMCIntegrator.h"

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

#include "Bianchi/TTHStudies/interface/Samples.h"
#include "Bianchi/TTHStudies/interface/CalcME.h"


#define GENJETDR 0.3
#define VERBOSE  false

using namespace std;
using namespace RooFit;


//extern float bar(float x);
extern "C" {
void bar_(float const *,float*);
}
float bar(float x) { float y; bar_(&x,&y); return y;}


extern "C" {
  void parameters_init_(double *Mass_E, double *Mass_M, double *Mass_L, double *Mass_U, double *Mass_D, double *Mass_S, double *Mass_C, double *Width_C, double *Mass_B, double *Width_B, double *Mass_T, double *Width_T, double *Mass_W, double *Width_W, double *Mass_Z, double *Width_Z, double *Mass_H, double *Width_H, double *Coupl_Alpha_QED, double *Coupl_Alpha_QCD, int *last_switch, int *amp_switch, int *amp_switch_rescue, int *use_coli_cache, int *check_Ward_tree, int *check_Ward_loop, int *out_symmetry, int *leading_colour);
  //void ppeexa03callme2born_(double *ME2, double cP[], int *no_prc, int *no_mod);
}

extern "C" {
  void set_permutation_pphtt_ttxhgg_1_( int perm[] );
  void amp2_pphtt_ttxhgg_1_( double[][4], double* );
}

extern "C" {
  void pphttxcallme2born_( double*, double[20], double*, double* );
}

double EH  = 120.;
double PH  = sqrt(EH*EH - 120.*120);
double ET  = (1000.- EH)/2;
double PT  = sqrt(ET*ET - 172*172);

/*
double cPP_[20] = {
		500.,   0.,   0.,    500., // g1
		500.,   0.,   0.,   -500., // g2
	EH,     0.,   0.,      0.,  // h
	ET,     0.,   PT,      0., // t
	ET,     0.,  -PT,      0.  // t~
};
*/

double cPP_[20] = {
  127.609, 0., 0., 127.609,
  1144.76, 0., 0., -1144.76,
  698.44, -127.365, 176.845, -652.628,
  371.474, 76.137, -170.77, -271.016,
  127.609+1144.76-698.44-371.474 , +127.365-76.137, -176.845+170.77 ,+652.628+271.016+127.609-1144.76
};

double Mt = 172; double Mh = 80;

double M_e   = 0.0005;
double M_mu  = 0.105;
double M_tau = 1.7;
double M_u   = 0.002;
double M_d   = 0.005;
double M_s   = 0.095;
double M_c   = 1.3;
double M_b   = 4.8;
double M_t   = 175;
double M_W   = 80.4;
double M_Z   = 91;
double M_H   = 125;
double Gamma_c = 0.01; // ???
double Gamma_b = 0.01; // ???
double Gamma_t = 2; 
double Gamma_W = 2; 
double Gamma_Z = 2.5; 
double Gamma_H = 0.01; // ??? 
double alpha_e = 1./137;
double alpha_S = 0.11;


double ccM_e     = M_e;
double ccM_mu    = M_mu;
double ccM_tau   = M_tau;
double ccM_u     = M_u;
double ccM_d     = M_d;
double ccM_s     = M_s;
double ccM_c     = M_c;
double ccGamma_c = Gamma_c;
double ccM_b     = M_b;
double ccGamma_b = Gamma_b;
double ccM_t     = M_t;
double ccGamma_t = Gamma_t;
double ccM_W     = M_W;
double ccGamma_W = Gamma_W;
double ccM_Z     = M_Z;
double ccGamma_Z = Gamma_Z;
double ccM_H     = M_H;
double ccGamma_H = Gamma_H;
double ccalpha_e = alpha_e;
double ccalpha_S = alpha_S;

int last_switch       = 1;
int amp_switch        = 1;
int amp_switch_rescue = 7; // 1 or 7 !!!
int use_coli_cache    = 1;
int check_Ward_tree   = 0;
int check_Ward_loop   = 0;
int out_symmetry      = 1;
int leading_colour    = 0;


typedef struct 
{
  float et; 
  float sumet;
  float sig;
  float phi;
} metInfo;


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


//namespace LHAPDF {
//    void   initPDFSet(int nset, const std::string& filename, int member=0);
//    int    numberPDF  (int nset);
//    void   usePDFMember(int nset, int member);
//    double xfx(int nset, double x, double Q, int fl);
//    double getXmin(int nset, int member);
//    double getXmax(int nset, int member);
//    double getQ2min(int nset, int member);
//    double getQ2max(int nset, int member);
//    void   extrapolate(bool extrapolate=true);
//}



int main(int argc, const char* argv[])
{

  std::cout << "TestMENew" << std::endl;
  gROOT->SetBatch(true);
 
  gSystem->Load("libFWCoreFWLite");
  gSystem->Load("libDataFormatsFWLite");

  //gSystem->Load("libbar.so");
  //cout << bar(10.) << endl;

  cout << "Start Init openloops " << endl;
  //parameters_init_(&ccM_e, &ccM_mu, &ccM_tau, &ccM_u, &ccM_d, &ccM_s, &ccM_c, &ccGamma_c, &ccM_b, &ccGamma_b, &ccM_t, &ccGamma_t, &ccM_W, &ccGamma_W, &ccM_Z, &ccGamma_Z, &ccM_H, &ccGamma_H, &ccalpha_e, &ccalpha_S, &last_switch, &amp_switch, &amp_switch_rescue, &use_coli_cache, &check_Ward_tree, &check_Ward_loop, &out_symmetry, &leading_colour);
  //cout << "Param init done " << endl;
  //int vec[] = {4, 5, 3, 2, 1};
  //set_permutation_pphtt_ttxhgg_1_( vec );
  //cout << "set permut done " << endl;
  //double P_scatt[6][4]; double M2;
  //P_scatt[0][0] = 100.; P_scatt[0][1] = 0.;P_scatt[0][2] = 0.;P_scatt[0][3] =  100.;
  //P_scatt[1][0] = 100.; P_scatt[1][1] = 0.;P_scatt[1][2] = 0.;P_scatt[1][3] = -100.;
  //P_scatt[2][0] = 100.; P_scatt[2][1] = 0.;P_scatt[2][2] = 0.;P_scatt[2][3] = -100.;
  //P_scatt[3][0] = 100.; P_scatt[3][1] = 0.;P_scatt[3][2] = 0.;P_scatt[3][3] = -100.;
  //P_scatt[4][0] = 100.; P_scatt[4][1] = 0.;P_scatt[4][2] = 0.;P_scatt[4][3] = -100.;
  //P_scatt[5][0] = 100.; P_scatt[5][1] = 0.;P_scatt[5][2] = 0.;P_scatt[5][3] = -100.;
  //amp2_pphtt_ttxhgg_1_( P_scatt, &M2 );
  //cout << "M2 = " << M2 << endl;

  TFile* f = new TFile("root/MEopenloop.root","RECREATE");
  TTree* treeME = new TTree("treeME", "");
  float me2;
  treeME->Branch("me2",&me2,"me2/F");

  //double M2_;
  //pphttxcallme2born_(&M2_, cPP_, &Mt, &Mh);
  //cout << "M2 = " << M2_ << endl;
  cout << "End Init openloops " << endl;

  //return 0;

  AutoLibraryLoader::enable();

  PythonProcessDesc builder(argv[1]);
  const edm::ParameterSet& in = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteInput");

  const edm::VParameterSet& samples = in.getParameter<edm::VParameterSet>("samples") ;
  std::string outFileName  ( in.getParameter<std::string>  ("outFileName" ) );
  std::string pathToFile   ( in.getParameter<std::string>  ("pathToFile" ) );
  std::string ordering     ( in.getParameter<std::string>  ("ordering" ) );
  bool verbose             ( in.getParameter<bool>  ("verbose") );
  double lumi              ( in.getParameter<double>("lumi") );

  float met       ( in.getParameter<double> ("met") );

  vector<int> evLimits( in.getParameter<vector<int> >  ("evLimits") );
  int evLow = evLimits[0];
  int evHigh = evLimits[1];



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


  for(unsigned int i = 0 ; i < mySampleFiles.size(); i++){
    
    string currentName       = mySampleFiles[i];

    mySamples->OpenFile( currentName );
    cout << "Opening file " << currentName << endl;
    TTree* currentTree       = mySamples->GetTree( currentName, "tree");
    cout << "Done!!" << endl;
    
    genParticleInfo genB, genBbar;
    genTopInfo genTop, genTbar;
    metInfo METtype1p2corr;
 
    currentTree->SetBranchAddress("genB",   &genB);
    currentTree->SetBranchAddress("genBbar",&genBbar);
    currentTree->SetBranchAddress("genTop", &genTop);
    currentTree->SetBranchAddress("genTbar",&genTbar);
    currentTree->SetBranchAddress("METtype1p2corr",&METtype1p2corr);
 
    int counter = 0;
    Long64_t nentries = currentTree->GetEntries();
    for (Long64_t i = 0; i < nentries ; i++){
      
      if(i%5000==0) cout << i << endl;
      currentTree->GetEntry(i);
      
      //continue;
      
      TLorentzVector topBLV   ;
      TLorentzVector topW1LV  ; 
      TLorentzVector topW2LV  ;
      TLorentzVector atopBLV  ; 
      TLorentzVector atopW1LV ; 
      TLorentzVector atopW2LV ;
      TLorentzVector genBLV   ;
      TLorentzVector genBbarLV;
      if(genTop.bmass>0){
	topBLV.SetPtEtaPhiM(  genTop.bpt,genTop.beta,genTop.bphi,genTop.bmass );
	topW1LV.SetPtEtaPhiM( genTop.wdau1pt,genTop.wdau1eta, genTop.wdau1phi,genTop.wdau1mass);
	topW2LV.SetPtEtaPhiM( genTop.wdau2pt,genTop.wdau2eta, genTop.wdau2phi,genTop.wdau2mass);
      }
      if(genTbar.bmass>0){
	atopBLV.SetPtEtaPhiM(  genTbar.bpt,genTbar.beta,genTbar.bphi,genTbar.bmass );
	atopW1LV.SetPtEtaPhiM( genTbar.wdau1pt,genTbar.wdau1eta, genTbar.wdau1phi,genTbar.wdau1mass);
	atopW2LV.SetPtEtaPhiM( genTbar.wdau2pt,genTbar.wdau2eta,genTbar.wdau2phi,genTbar.wdau2mass);
      }
      if(genB.mass>0 && genB.momid==25){
	genBLV.SetPtEtaPhiM(genB.pt,genB.eta ,genB.phi, genB.mass );
      }
      if(genBbar.mass>0 && genBbar.momid==25){
	genBbarLV.SetPtEtaPhiM(genBbar.pt,genBbar.eta ,genBbar.phi, genBbar.mass );
      }
  
      TLorentzVector TOPLEP;
      TLorentzVector TOPHAD;
      TLorentzVector TOPHADW1;
      TLorentzVector TOPHADW2;
      TLorentzVector TOPHADB;
      TLorentzVector TOPLEPW1;
      TLorentzVector TOPLEPW2;
      TLorentzVector TOPLEPB;
      TLorentzVector HIGGS;

      float neutCosTheta;

      bool properEvent = (genBLV.Pt()>0 && genBbarLV.Pt()>0 && topBLV.Pt()>0 && topW1LV.Pt()>0 && topW2LV.Pt()>0 && atopBLV.Pt()>0 && atopW1LV.Pt()>0 && atopW2LV.Pt()>0);

      HIGGS.SetPxPyPzE( (genBLV+genBbarLV).Px(), (genBLV+genBbarLV).Py(),(genBLV+genBbarLV).Pz(),(genBLV+genBbarLV).E());
      if( abs(genTop.wdau1id)>6 && abs(genTbar.wdau1id)<6){
	TOPLEP.SetPxPyPzE( (topBLV+topW1LV+topW2LV).Px(), (topBLV+topW1LV+topW2LV).Py(), (topBLV+topW1LV+topW2LV).Pz(), (topBLV+topW1LV+topW2LV).E() );
	TOPHAD.SetPxPyPzE( (atopBLV+atopW1LV+atopW2LV).Px(), (atopBLV+atopW1LV+atopW2LV).Py(), (atopBLV+atopW1LV+atopW2LV).Pz(), (atopBLV+atopW1LV+atopW2LV).E() );
	TOPHADW1.SetPxPyPzE( atopW1LV.Px(), atopW1LV.Py(), atopW1LV.Pz(), atopW1LV.E());
	TOPHADW2.SetPxPyPzE( atopW2LV.Px(), atopW2LV.Py(), atopW2LV.Pz(), atopW2LV.E());
	TOPHADB.SetPxPyPzE( atopBLV.Px(),  atopBLV.Py(),   atopBLV.Pz(),  atopBLV.E());
	if( abs(genTop.wdau1id)==11 || abs(genTop.wdau1id)==13 || abs(genTop.wdau1id)==15 ){
	  TOPLEPW1.SetPxPyPzE  ( topW1LV.Px(), topW1LV.Py(), topW1LV.Pz(), topW1LV.E());
	  TOPLEPW2.SetPtEtaPhiM( topW2LV.Pt(), 0.0, topW2LV.Phi(), 0.0);
	  neutCosTheta = TMath::Cos( topW2LV.Vect().Theta() );
	}
	else{
	  TOPLEPW2.SetPtEtaPhiM( topW1LV.Pt(), 0.0, topW1LV.Phi(), 0.0);
	  TOPLEPW1.SetPxPyPzE  ( topW2LV.Px(), topW2LV.Py(), topW2LV.Pz(), topW2LV.E());
	  neutCosTheta = TMath::Cos( topW1LV.Vect().Theta() );
	}
	TOPLEPB.SetPxPyPzE(  topBLV.Px(),  topBLV.Py(),   topBLV.Pz(), topBLV.E());
      }
      else if(abs(genTop.wdau1id)<6 && abs(genTbar.wdau1id)>6){
	TOPHAD.SetPxPyPzE( (topBLV+topW1LV+topW2LV).Px(), (topBLV+topW1LV+topW2LV).Py(), (topBLV+topW1LV+topW2LV).Pz(), (topBLV+topW1LV+topW2LV).E() );
	TOPLEP.SetPxPyPzE( (atopBLV+atopW1LV+atopW2LV).Px(), (atopBLV+atopW1LV+atopW2LV).Py(), (atopBLV+atopW1LV+atopW2LV).Pz(), (atopBLV+atopW1LV+atopW2LV).E() );
	TOPHADW1.SetPxPyPzE( topW1LV.Px(), topW1LV.Py(), topW1LV.Pz(), topW1LV.E());
	TOPHADW2.SetPxPyPzE( topW2LV.Px(), topW2LV.Py(), topW2LV.Pz(), topW2LV.E());
	TOPHADB.SetPxPyPzE( topBLV.Px(),  topBLV.Py(),   topBLV.Pz(),  topBLV.E());
	if( abs(genTbar.wdau1id)==11 || abs(genTbar.wdau1id)==13 || abs(genTbar.wdau1id)==15 ){
	  TOPLEPW1.SetPxPyPzE  ( atopW1LV.Px(), atopW1LV.Py(), atopW1LV.Pz(), atopW1LV.E());
	  TOPLEPW2.SetPtEtaPhiM( atopW2LV.Pt(), 0.0, atopW2LV.Phi(), 0.0);
	  neutCosTheta = TMath::Cos( atopW2LV.Vect().Theta() );
	}
	else{
	  TOPLEPW2.SetPtEtaPhiM( atopW1LV.Pt(), 0.0, atopW1LV.Phi(), 0.0);
	  TOPLEPW1.SetPxPyPzE  ( atopW2LV.Px(), atopW2LV.Py(), atopW2LV.Pz(), atopW2LV.E());
	  neutCosTheta = TMath::Cos( atopW1LV.Vect().Theta() );
	}
	TOPLEPB.SetPxPyPzE( atopBLV.Px(),  atopBLV.Py(),   atopBLV.Pz(),  atopBLV.E());
      }      
      else{
	properEvent=false;
      }

      properEvent = ( TOPLEPW1.Pt() >20 && TMath::Abs(TOPLEPW1.Eta()) <2.1 &&
		      TOPLEPB.Pt()  >30 && TMath::Abs(TOPLEPB.Eta())  <2.5 &&
		      TOPHADW1.Pt() >30 && TMath::Abs(TOPHADW1.Eta()) <2.5 &&
		      TOPHADW1.Pt() >30 && TMath::Abs(TOPHADW2.Eta()) <2.5 &&
		      TOPHADB.Pt()  >30 && TMath::Abs(TOPHADB.Eta())  <2.5 &&
		      genBLV.Pt()   >30 && TMath::Abs(genBLV.Eta())   <2.5 &&
		      genBbarLV.Pt()>30 && TMath::Abs(genBbarLV.Eta())<2.5
		      );

      if(!properEvent) continue;
      counter++;

      if( !(counter>=evLow && counter<=evHigh) ) continue;
      cout << "Processing event # " << counter << endl;

      TLorentzVector v3 = (TOPLEPW1+TOPLEPW2+TOPLEPB);           
      TLorentzVector v4 = (TOPHADW1+TOPHADW2+TOPHADB);           
      TLorentzVector v5 = (genBLV+genBbarLV);                   
      TLorentzVector vSum = v3+v4+v5;

      v3.SetPtEtaPhiM( (TOPLEPW1+TOPLEPW2+TOPLEPB).Pt(), (TOPLEPW1+TOPLEPW2+TOPLEPB).Eta(), (TOPLEPW1+TOPLEPW2+TOPLEPB).Phi(), Mt );
      v4.SetPtEtaPhiM( (TOPHADW1+TOPHADW2+TOPHADB).Pt(), (TOPHADW1+TOPHADW2+TOPHADB).Eta(), (TOPHADW1+TOPHADW2+TOPHADB).Phi(), Mt );
      v5.SetPtEtaPhiM( (genBLV+genBbarLV).Pt(), (genBLV+genBbarLV).Eta(), (genBLV+genBbarLV).Phi(), Mh );

      //cout << "Before : " << endl;
      //cout << "Pt = " << vSum.Pt() << ", Pz = " << vSum.Pz() << endl;
      TVector3 boostPt( vSum.Px()/vSum.E(), vSum.Py()/vSum.E(), 0.0 );
      vSum.Boost( -boostPt );
      //cout << "After : " << endl;
      //cout << "Pt = " << vSum.Pt() << ", Pz = " << vSum.Pz() << endl;

      v3.Boost( -boostPt );
      v4.Boost( -boostPt );
      v5.Boost( -boostPt );

      double v5Px = -(v3.Px() + v4.Px());
      double v5Py = -(v3.Py() + v4.Py());
      double v5Pz = v5.Pz();
      v5.SetPxPyPzE( v5Px, v5Py, v5Pz, sqrt(v5Px*v5Px + v5Py*v5Py + v5Pz*v5Pz + Mh*Mh) );

      double E  = v3.E() +v4.E() + v5.E();
      double Pz = v3.Pz()+v4.Pz()+ v5.Pz();
      TLorentzVector v1 = TLorentzVector(0.,0.,  (E+Pz)/2., (E+Pz)/2.); 
      TLorentzVector v2 = TLorentzVector(0.,0., -(E-Pz)/2., (E-Pz)/2.); 
      //v1.SetPxPyPzE(0.,0.,  (E+Pz)/2., (E+Pz)/2.);
      //v2.SetPxPyPzE(0.,0., -(E-Pz)/2., (E-Pz)/2.);
    
      //cout <<
      //v1.E() << ", " <<  v1.Px() << ", " <<  v1.Py() << ", " <<  v1.Pz() << " => " << v1.M() << endl << 
      //v2.E() << ", " <<  v2.Px() << ", " <<  v2.Py() << ", " <<  v2.Pz() << " => " << v2.M() << endl <<
      //v5.E() << ", " <<  v5.Px() << ", " <<  v5.Py() << ", " <<  v5.Pz() << " => " << v5.M()/120. - 1 << endl <<
      //v3.E() << ", " <<  v3.Px() << ", " <<  v3.Py() << ", " <<  v3.Pz() << " => " << v3.M()/172. -1 << endl <<
      //v4.E() << ", " <<  v4.Px() << ", " <<  v4.Py() << ", " <<  v4.Pz() << " => " << v4.M()/172. -1 << endl;

      //TLorentzVector tot = v1 + v2 - v3 - v4 - v5;
      //cout << "balance = " << tot.E() << ", " << tot.Px() << ", " << tot.Py() << ", " << tot.Pz() << endl;

      double MhDefault = 120.;
      double M2;
      double ccP[20] = {
	v1.E(), v1.Px(), v1.Py(), v1.Pz(),
	v2.E(), v2.Px(), v2.Py(), v2.Pz(),
	v5.E(), v5.Px(), v5.Py(), v5.Pz(),
	v3.E(), v3.Px(), v3.Py(), v3.Pz(),
	v4.E(), v4.Px(), v4.Py(), v4.Pz()
      };
      pphttxcallme2born_(&M2, ccP, &Mt, &Mh);
      //cout << "M2 = " << M2 << endl;
      me2 = M2;
      treeME->Fill();
      continue;


      CalcME* calc = new CalcME(
				v1,v2,v3,v4,v5,
				175., met);
      cout << calc->meSquared() << endl;
      delete calc;

      TVector3 boost = genBLV.BoostVector();
      v1.Boost(-boost);
      v2.Boost(-boost);
      v3.Boost(-boost); 
      v4.Boost(-boost);
      v5.Boost(-boost); 

      CalcME* calc2 = new CalcME(
				 v1,v2,v3,v4,v5,
				 175., met);
      cout << calc2->meSquared() << endl;
      delete calc2;
      cout << "<=> " << 4*v3.Dot(v4) + 4*175*175 << endl;
      cout << "t1 : Px,Py,Pz,E = " << v3.Px() << ", " << v3.Py() << ", " << v3.Pz() << ", " << v3.E() << endl;
      cout << "t2 : Px,Py,Pz,E = " << v4.Px() << ", " << v4.Py() << ", " << v4.Pz() << ", " << v4.E() << endl;
 
    }
    
  }
  
  f->cd();
  treeME->Write("",TObject::kOverwrite);
  f->Close();

  cout << "Finished!!!" << endl;
  return 0;
  

}
