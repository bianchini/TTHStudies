import FWCore.ParameterSet.Types  as CfgTypes
import FWCore.ParameterSet.Config as cms

VType     = "_VType2"
xsecTT_SL = 103.0

process = cms.Process("TestMENew")

process.fwliteInput = cms.PSet(

    outFileName   = cms.string("./root/TestMENew.root"),
    madweight     = cms.string("weights_ttjets"),
    pathToTF      = cms.string("./root/transferFunctions_partonE.root"),
    pathToFile    = cms.string("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store//user/bianchi/HBB_EDMNtuple/AllHDiJetPt"+VType+"/v2/"),
    ordering      = cms.string("DiJetPt_"),
    lumi          = cms.double(12.1),

    samples       = cms.VPSet(


    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTH_HToBB_M-110_8TeV-pythia6'+VType),
    nickName = cms.string('TTH110'),
    color    = cms.int32(2),
    xSec     = cms.double(0.1887*0.744)
    ),

    
    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTH_HToBB_M-115_8TeV-pythia6'+VType),
    nickName = cms.string('TTH115'),
    color    = cms.int32(2),
    xSec     = cms.double(0.1663*0.703)
    ),


    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTH_HToBB_M-120_8TeV-pythia6'+VType),
    nickName = cms.string('TTH120'),
    color    = cms.int32(2),
    xSec     = cms.double(0.1470*0.648)
    ),


    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('TTH_HToBB_M-125_8TeV-pythia6'+VType),
    nickName = cms.string('TTH125'),
    color    = cms.int32(2),
    xSec     = cms.double(0.1302*0.569)
    ),


    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTH_HToBB_M-130_8TeV-pythia6'+VType),
    nickName = cms.string('TTH130'),
    color    = cms.int32(2),
    xSec     = cms.double(0.1157*0.494)
    ),

    
    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTH_HToBB_M-135_8TeV-pythia6'+VType),
    nickName = cms.string('TTH135'),
    color    = cms.int32(2),
    xSec     = cms.double(0.1031*0.404)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTWJets_8TeV-madgraph'+VType),
    nickName = cms.string('TTW'),
    color    = cms.int32(18),
    xSec     = cms.double(0.232),
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTZJets_8TeV-madgraph'+VType),
    nickName = cms.string('TTZ'),
    color    = cms.int32(18),
    xSec     = cms.double(0.2057),
    ),


    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTJets_SemiLeptMGDecays_8TeV-madgraph-part'+VType),
    nickName = cms.string('TTJetsSemiLept'),
    color    = cms.int32(41),
    xSec     = cms.double(xsecTT_SL),
    ),
    
    
    ),

    vegasPoints   = cms.int32(15000000), #6000000
    mode          = cms.untracked.int32(2),

    functions     = cms.vstring('1.39e+23*x^(-2.81e+00)',
                                '7.84e+17*TMath::Landau(x,6.17e+01,1.61e+01)',
                                'x>=12 ? x^(-2.010e-01)*exp((-1.5785e-02)*x) : 4.184e-02*x'
                                ),
    
    useME         = cms.int32(1),
    useJac        = cms.int32(1),
    useMET        = cms.int32(0),
    useTF         = cms.int32(0),
    usePDF        = cms.int32(1),

    shiftMomenta      = cms.int32(1),
    testMassScan      = cms.int32(1),
    testPermutations  = cms.int32(1),
    printP4           = cms.int32(1),
    
    compAcceptance    = cms.untracked.int32(0),
    
    verbose       = cms.bool(False),
    met           = cms.double(120.),
    masses        = cms.vdouble(120,125),
    evLimits      = cms.vint32(1,1),

    pertBLep      = cms.double(1.0),
    pertW1        = cms.double(1.0),
    pertW2        = cms.double(1.0),
    pertBHad      = cms.double(1.0),

    enlargeE1     = cms.double(0.),
    enlargeEh1    = cms.double(0.),
    enlargePt     = cms.double(0.),

    scaleH        = cms.double(0.18),
    scaleL        = cms.double(0.15),
    scaleMET      = cms.double(20),

    )
