import FWCore.ParameterSet.Types as CfgTypes
import FWCore.ParameterSet.Config as cms

VType = "_VType2"

xsecTT_FH = 106.9
xsecTT_SL = 103.0
xsecTT_FL = 24.8

process = cms.Process("TreeProducer")

process.fwliteInput = cms.PSet(

    pathToFile    = cms.string("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store//user/bianchi/HBB_EDMNtuple/AllHDiJetPt"+VType+"/v2/"),
    ordering      = cms.string("DiJetPt_"),
    lumi          = cms.double(12.1),
    verbose       = cms.bool(False),
    computeCMSVariables = cms.bool(True),
    printP4       = cms.bool(True),
    toPrint       = cms.int32(10),

    cut = cms.string(""),
    
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



)
