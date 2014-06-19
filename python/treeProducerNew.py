import FWCore.ParameterSet.Types as CfgTypes
import FWCore.ParameterSet.Config as cms


xsecTT_FH = 106.9
xsecTT_SL = 103.0
xsecTT_FL = 24.8

process = cms.Process("TreeProducer")

process.fwliteInput = cms.PSet(

    pathToFile    = cms.string("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store//user/bianchi/HBB_EDMNtuple/AllHDiJetPt_V4/"),
    ordering      = cms.string("DiJetPt_"),
    lumi          = cms.double(12.1),
    verbose       = cms.bool(False),
    evalReg       = cms.bool(False),
    maxnum        = cms.int32(100000),

    samples       = cms.VPSet(
    
    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('TTH_HToBB_M-120_8TeV-pythia6'),
    nickName = cms.string('TTH120'),
    color    = cms.int32(2),
    xSec     = cms.double(0.1470*0.648)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTJets_SemiLeptMGDecays_8TeV-madgraph'),
    nickName = cms.string('TTJetsSemiLept'),
    color    = cms.int32(41),
    xSec     = cms.double(xsecTT_SL),
    ),
    
    ),



)
