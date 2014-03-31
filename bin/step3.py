import FWCore.ParameterSet.Types as CfgTypes
import FWCore.PythonUtilities.LumiList as LumiList

import FWCore.ParameterSet.Config as cms



process = cms.Process("Step3")


process.fwliteInput = cms.PSet(

    # LOCAL 
    #pathToFile    = cms.string("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store//user/bianchi/HBB_EDMNtuple/ZllHDiJetPt/")
    #pathToFile    = cms.string("/scratch/bianchi/HBB_EDMNtuple/All.H.DiJetPt/"),
    #ordering      = cms.string("ZllH.DiJetPt.Oct22."),

    # BATCH
    pathToFile    = cms.string("gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt/"),
    outPath       = cms.string("gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_Step2/"),
    ordering      = cms.string("DiJetPt_"),

    lumi          = cms.double(12.1),
    verbose       = cms.bool(False),
    eventsToDebug = cms.int32(-1),

    lumisToProcess      = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange()),
    

    samples       = cms.VPSet(

    cms.PSet(
    skip     = cms.bool(True),
    name     = cms.string('DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball'),
    nickName = cms.string('DYJets'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71),
    update   = cms.bool(True)
    ),

    cms.PSet(
    skip     = cms.bool(False),
    name     = cms.string('SingleMuRun2012'),
    nickName = cms.string('SingleMu'),
    color    = cms.int32(18),
    xSec     = cms.double(-1),
    update   = cms.bool(False)
    ),

    
    ),


    )


#JSONfile = 'Cert_190456-203002_8TeV_PromptReco_Collisions12_JSON.txt'
#JSONfile = 'Cert_190456-196531_8TeV_13Jul2012ReReco_ert_190782-190949_8TeV_06Aug2012ReReco_Cert_198022-198523_8TeV_24Aug2012ReReco_Cert_198941-203002_8TeV_PromptReco_Collisions12_JSON.txt'
JSONfile = 'Cert_190456-203002_8TeV_PromptReco_Collisions12_JSON.txt'
lumiList = LumiList.LumiList (filename = JSONfile).getCMSSWString().split(',')
process.fwliteInput.lumisToProcess.extend(lumiList)
