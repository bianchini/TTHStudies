#ifndef HELPERFUNCTIONS_H // header guards
#define HELPERFUNCTIONS_H


#include <cassert>
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "TMath.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TString.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#endif

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
  TLorentzVector p4;
  float csv;
} JetObservable;


struct JetObservableListerByPt
{
  bool operator()( const JetObservable &lx, const JetObservable& rx ) const {
    return (lx.p4).Pt() > (rx.p4).Pt();
  }
};

struct JetObservableListerByCSV
{
  bool operator()( const JetObservable &lx, const JetObservable& rx ) const {
    return (lx.csv) > (rx.csv);
  }
};


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


//https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution
float resolutionBias(float eta, int shift, int alreadyCorrect = 0)
{
  if(eta < 0.5){
    if(shift ==  0) 
      return alreadyCorrect ?  0. : 0.052 ;
    if(shift == +1) 
      return 0.115 - resolutionBias(eta, 0 , !alreadyCorrect);  
    if(shift == -1) 
      return 0.010 - resolutionBias(eta, 0 , !alreadyCorrect);  
  }

  if(eta < 1.1){
    if(shift ==  0) 
      return alreadyCorrect ?  0. : 0.057 ;
    if(shift == +1) 
      return 0.114 - resolutionBias(eta, 0 , !alreadyCorrect);  
    if(shift == -1) 
      return 0.001 - resolutionBias(eta, 0 , !alreadyCorrect);  
  }

  if(eta < 1.7){
    if(shift ==  0) 
      return alreadyCorrect ?  0. : 0.096 ;
    if(shift == +1) 
      return 0.161 - resolutionBias(eta, 0 , !alreadyCorrect);  
    if(shift == -1) 
      return 0.032 - resolutionBias(eta, 0 , !alreadyCorrect);  
  }

  if(eta < 2.3){
    if(shift ==  0) 
      return alreadyCorrect ?  0. : 0.134 ;
    if(shift == +1) 
      return 0.228 - resolutionBias(eta, 0 , !alreadyCorrect);  
    if(shift == -1) 
      return 0.042 - resolutionBias(eta, 0 , !alreadyCorrect);  
  }

  if(eta < 5.0){
    if(shift ==  0) 
      return alreadyCorrect ?  0. : 0.288 ;
    if(shift == +1) 
      return 0.488 - resolutionBias(eta, 0 , !alreadyCorrect);  
    if(shift == -1) 
      return 0.089 - resolutionBias(eta, 0 , !alreadyCorrect);  
  }

  return 0.;
}

std::vector<std::string> GetInputExpressionsReg() {
    std::vector<std::string> values;
    values.resize(15, "");
    values[0] = "breg_rawptJER := smear_pt_res(hJet_ptRaw, hJet_genPt, hJet_eta)";
    values[1] = "breg_pt := hJet_pt";
    values[2] = "breg_et := evalEt(hJet_pt, hJet_eta, hJet_phi, hJet_e)";
    values[3] = "breg_mt := evalMt(hJet_pt, hJet_eta, hJet_phi, hJet_e)";
    values[4] = "breg_leadtrackpt := max(0,hJet_ptLeadTrack)";
    values[5] = "breg_vtx3dL := max(0,hJet_vtx3dL)";
    values[6] = "breg_vtx3deL := max(0,hJet_vtx3deL)";
    values[7] = "breg_vtxMass := max(0,hJet_vtxMass)";
    values[8] = "breg_vtxPt := max(0,hJet_vtxPt)";
    values[9] = "breg_cef := hJet_cef";
    values[10] = "breg_ntot := hJet_nconstituents";
    values[11] = "breg_eJEC := hJet_JECUnc";
    values[12] = "breg_softlepptrel := max(0,hJet_SoftLeptptRel*(hJet_SoftLeptIdlooseMu==1 || hJet_SoftLeptId95==1))";
    values[13] = "breg_softleppt := max(0,hJet_SoftLeptPt*(hJet_SoftLeptIdlooseMu==1 || hJet_SoftLeptId95==1))";
    values[14] = "breg_softlepdR := max(0,hJet_SoftLeptdR*(hJet_SoftLeptIdlooseMu==1 || hJet_SoftLeptId95==1))";
    return values;
}

double quadsum(double a, double b)
{
    return TMath::Sqrt(a*a + b*b);
}

float smear_pt_res(float pt, float genpt, float eta)
{
    eta = fabs(eta);

    if (genpt>15. && (fabs(pt - genpt) / pt)<0.5) {  // limit the effect to the core

        double res    = 1.0;

        if (eta <= 0.5) {
            res    = 1.052;
        } else if (0.5 < eta && eta <= 1.1) {
            res    = 1.057;
        } else if (1.1 < eta && eta <= 1.7) {
            res    = 1.096;
        } else if (1.7 < eta && eta <= 2.3) {
            res    = 1.134;
        } else {
            res    = 1.288;
        }

        float deltapt = (pt - genpt) * res;
        return TMath::Max(float(0.), genpt + deltapt);
    }
    return pt;
}

double evalEt(double pt, double eta, double phi, double e)
{
    TLorentzVector j;
    j.SetPtEtaPhiE(pt, eta, phi, e);
    return j.Et();
}

double evalMt(double pt, double eta, double phi, double e)
{
    TLorentzVector j;
    j.SetPtEtaPhiE(pt, eta, phi, e);
    return j.Mt();
}


TMVA::Reader* getTMVAReader( std::string regMethod , Float_t readerVars[15]){

  TMVA::Tools::Instance();  //< This loads the library
  TMVA::Reader *reader = new TMVA::Reader("!Color:!Silent");

  /// Get the variables
  const std::vector<std::string>& inputExpressionsReg = GetInputExpressionsReg();
  const UInt_t nvars = inputExpressionsReg.size();

  for (UInt_t iexpr = 0; iexpr < nvars; iexpr++) {
    const TString& expr = inputExpressionsReg.at(iexpr);
    reader->AddVariable(expr, &readerVars[iexpr]);
  }
  TString weightdir  = "../../../VHbbUF/weights/";
  TString weightfile = weightdir + "TMVARegression_" + regMethod + ".testweights.xml";
  reader->BookMVA(regMethod + " method", weightfile);

  return reader;
}

float getRegressionEnergy(std::string regMethod, TMVA::Reader *reader, Float_t readerVars[15], TTree* tree, Long64_t entry , int jetpos, int verbose){

  float Jet_pt[99],         Jet_genPt[99],         Jet_ptRaw[99],        Jet_eta[99],  Jet_phi[99],    Jet_e[99], Jet_JECUnc[99];
  float Jet_vtx3dL[99],     Jet_vtx3deL[99],       Jet_vtxMass[99],      Jet_vtxPt[99];
  float Jet_cef[99],        Jet_ptLeadTrack[99],   Jet_nconstituents[99];
  float Jet_SoftLeptPt[99], Jet_SoftLeptptRel[99], Jet_SoftLeptdR[99];
  int Jet_SoftLeptIdlooseMu[99], Jet_SoftLeptId95[99];

  TString jetColl = jetpos>=0 ? "h" : "a";
  int j      = jetpos>=0 ? jetpos : -jetpos-1;

  tree->SetBranchAddress(jetColl+"Jet_pt",      Jet_pt);
  tree->SetBranchAddress(jetColl+"Jet_genPt",   Jet_genPt);
  tree->SetBranchAddress(jetColl+"Jet_ptRaw",   Jet_ptRaw);
  tree->SetBranchAddress(jetColl+"Jet_eta",     Jet_eta);
  tree->SetBranchAddress(jetColl+"Jet_phi",     Jet_phi);
  tree->SetBranchAddress(jetColl+"Jet_e",       Jet_e);
  tree->SetBranchAddress(jetColl+"Jet_JECUnc",  Jet_JECUnc);
  tree->SetBranchAddress(jetColl+"Jet_vtx3dL",  Jet_vtx3dL );
  tree->SetBranchAddress(jetColl+"Jet_vtx3deL", Jet_vtx3deL);
  tree->SetBranchAddress(jetColl+"Jet_vtxMass", Jet_vtxMass);
  tree->SetBranchAddress(jetColl+"Jet_vtxPt",   Jet_vtxPt);
  tree->SetBranchAddress(jetColl+"Jet_cef",           Jet_cef);
  tree->SetBranchAddress(jetColl+"Jet_ptLeadTrack",   Jet_ptLeadTrack  );
  tree->SetBranchAddress(jetColl+"Jet_nconstituents", Jet_nconstituents);
  tree->SetBranchAddress(jetColl+"Jet_SoftLeptPt",        Jet_SoftLeptPt);
  tree->SetBranchAddress(jetColl+"Jet_SoftLeptptRel",     Jet_SoftLeptptRel);
  tree->SetBranchAddress(jetColl+"Jet_SoftLeptdR",        Jet_SoftLeptdR);
  tree->SetBranchAddress(jetColl+"Jet_SoftLeptIdlooseMu", Jet_SoftLeptIdlooseMu );
  tree->SetBranchAddress(jetColl+"Jet_SoftLeptId95",      Jet_SoftLeptId95);
 

  tree->GetEntry(entry);

  // fill the input variables
  readerVars[0]  = smear_pt_res( Jet_ptRaw[j], Jet_genPt[j], Jet_eta[j]);
  readerVars[1]  = Jet_pt[j];
  readerVars[2]  = evalEt( Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_e[j]);
  readerVars[3]  = evalMt( Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_e[j]);
  readerVars[4]  = TMath::Max(float(0.), Jet_ptLeadTrack[j] );
  readerVars[5]  = TMath::Max(float(0.), Jet_vtx3dL[j] );
  readerVars[6]  = TMath::Max(float(0.), Jet_vtx3deL[j] );
  readerVars[7]  = TMath::Max(float(0.), Jet_vtxMass[j] );
  readerVars[8]  = TMath::Max(float(0.), Jet_vtxPt[j] );;
  readerVars[9]  = Jet_cef[j];
  readerVars[10] = Jet_nconstituents[j];
  readerVars[11] = Jet_JECUnc[j];
  readerVars[12] = TMath::Max(float(0.), Jet_SoftLeptptRel[j]*(Jet_SoftLeptIdlooseMu[j]==1 || Jet_SoftLeptId95[j]==1));
  readerVars[13] = TMath::Max(float(0.), Jet_SoftLeptPt[j]*   (Jet_SoftLeptIdlooseMu[j]==1 || Jet_SoftLeptId95[j]==1));
  readerVars[14] = TMath::Max(float(0.), Jet_SoftLeptdR[j]*   (Jet_SoftLeptIdlooseMu[j]==1 || Jet_SoftLeptId95[j]==1));

  float output = (reader->EvaluateRegression(regMethod+" method"))[0];

  if(verbose){
    cout << ">>> Regression" << endl;
    cout << "   --> input:" << endl;
    for(int k = 0 ; k < 15 ; k++){
      cout << "        -> " << "readerVars[" << k << "] = " <<  readerVars[k] << endl;
    }
    cout << "    --> output = " << output << endl;
  }

  return output;

}



#endif
