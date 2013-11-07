import FWCore.ParameterSet.Types as CfgTypes
import FWCore.ParameterSet.Config as cms

VType = "_VType2"

xsecTT_FH = 106.9
xsecTT_SL = 103.0
xsecTT_FL = 24.8

process = cms.Process("TopHinFitterSL")

process.fwliteInput = cms.PSet(

    pathToFile    = cms.string("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store//user/bianchi/HBB_EDMNtuple/AllHDiJetPt"+VType+"/v2/"),
    ordering      = cms.string("DiJetPt_"),
    lumi          = cms.double(12.1),
    verbose             = cms.bool(False),
    computeCMSVariables = cms.bool(True),
    
    samples       = cms.VPSet(
    
    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('TTH_HToBB_M-125_8TeV-pythia6'+VType),
    nickName = cms.string('TTH125'),
    color    = cms.int32(2),
    xSec     = cms.double(0.1302*0.569)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTJets_SemiLeptMGDecays_8TeV-madgraph-part'+VType),
    nickName = cms.string('TTJetsSemiLept'),
    color    = cms.int32(41),
    xSec     = cms.double(xsecTT_SL),
    ),
    
    ),

    likelihoodFile = cms.string("../test/Macro/likelihoods_bkgd.root"),
    #likelihoods    = cms.vstring("KIN","BTAG","HEL"),
    likelihoods    = cms.vstring("KIN","HMASS", "HEL"),

    udscPtCut = cms.double(30),
    bPtCut    = cms.double(30),

    metResolutionsCoeff = cms.vdouble(1.12544e+02, 3.00905e-01, 6.54808e+01, 7.35296e-01),
    helicityCoeff       = cms.vdouble(0.0315287, -0.00783723, -0.0290744, 0.00656129,      # leptonic pol3
                                      0.0296338, 0.00, -0.0121257, 0.00, -0.0182667),      # hadronic pol4
    chi2Coeff           = cms.vdouble(3.09651e+00),
    sgnMassCoeff        = cms.vdouble(1.15716e+02,1.98236e+01),
    
    udscResolutions = cms.VPSet(

    cms.PSet(
    bin = cms.string("abs(eta)>=0.  && abs(eta)<1.0"),
    et  = cms.string("sqrt(8.51479e+00 + 9.85110e-01*9.85110e-01*et + 5.50086e-02*5.50086e-02*et*et)"),
    eta = cms.string("2.63342e-02 + 3.48643e-03*abs(eta) + (-2.62017e-03)*eta*eta + 1.45496e-03*abs(eta*eta*eta)"),
    phi = cms.string("0.0323073")
    ),

    cms.PSet(
    bin = cms.string("abs(eta)>=1.0 && abs(eta)<2.0"),
    et  = cms.string("sqrt(9.93662e+00 + 1.13533e+00*1.13533e+00*et + 4.88416e-02*4.88416e-02*et*et)"),
    eta = cms.string("2.63342e-02 + 3.48643e-03*abs(eta) + (-2.62017e-03)*eta*eta + 1.45496e-03*abs(eta*eta*eta)"),
    phi = cms.string("0.0323073")
    ),

    cms.PSet(
    bin = cms.string("abs(eta)>=2.0 && abs(eta)<5.0"),
    et  = cms.string("sqrt(1.50064e+01 + 9.67855e-01*9.67855e-01*et + 1.15671e-02*1.15671e-02*et*et)"),
    eta = cms.string("2.63342e-02 + 3.48643e-03*abs(eta) + (-2.62017e-03)*eta*eta + 1.45496e-03*abs(eta*eta*eta)"),
    phi = cms.string("0.0323073")
    ),
  
   
    ),

    bResolutions = cms.VPSet(

    cms.PSet(
    bin = cms.string("abs(eta)>=0.  && abs(eta)<1.0"),
    et  = cms.string("sqrt(4.11576e+01 + 1.38054e+00*1.38054e+00*et + 0.00000e+00*0.00000e+00*et*et)"),
    eta = cms.string("2.38983e-02 + 2.23730e-03*abs(eta) + (-1.00830e-03)*eta*eta + 1.04913e-03*abs(eta*eta*eta)"),
    phi = cms.string("0.0290914")
    ),

    cms.PSet(
    bin = cms.string("abs(eta)>=1.0 && abs(eta)<2.0"),
    et  = cms.string("sqrt(3.75956e+01 + 1.46691e+00*1.46691e+00*et + 0.00000e+00*0.00000e+00*et*et)"),
    eta = cms.string("2.38983e-02 + 2.23730e-03*abs(eta) + (-1.00830e-03)*eta*eta + 1.04913e-03*abs(eta*eta*eta)"),
    phi = cms.string("0.0290914")
    ),

    cms.PSet(
    bin = cms.string("abs(eta)>=2.0 && abs(eta)<5.0"),
    et  = cms.string("sqrt(5.40734e+01 + 1.19942e+00*1.19942e+00*et + 1.36122e-06*1.36122e-06*et*et)"),
    eta = cms.string("2.38983e-02 + 2.23730e-03*abs(eta) + (-1.00830e-03)*eta*eta + 1.04913e-03*abs(eta*eta*eta)"),
    phi = cms.string("0.0290914")
    ),
    
    ),

    
    lepResolutions = cms.VPSet(

    cms.PSet(
    bin = cms.string("abs(eta)>=0.  && abs(eta)<1.0"),
    et  = cms.string("sqrt(1.44169e-04*et*et + 4.15168e-08*et*et*et*et)"),
    eta = cms.string("0.0003"),
    phi = cms.string("0.000148728")
    ),

    cms.PSet(
    bin = cms.string("abs(eta)>=1.0 && abs(eta)<2.0"),
    et  = cms.string("sqrt(3.77878e-04*et*et + 5.65095e-08*et*et*et*et)"),
    eta = cms.string("0.0003"),
    phi = cms.string("0.000148728")
    ),

    cms.PSet(
    bin = cms.string("abs(eta)>=2.0 && abs(eta)<2.5"),
    et  = cms.string("sqrt(9.02986e-04*et*et + 3.24309e-07*et*et*et*et)"),
    eta = cms.string("0.0003"),
    phi = cms.string("0.000148728")
    ),
 
    
    ),


    )
