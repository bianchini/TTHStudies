import FWCore.ParameterSet.Types as CfgTypes
import FWCore.ParameterSet.Config as cms

process = cms.Process("TFuncSL")

process.fwliteInput = cms.PSet(

    pathToFile    = cms.string("./root/treeProducer_new.root"),
    outFileName   = cms.string("transferFunctions_partonE_recoil.root"),
    verbose       = cms.bool(False),
    
    
    )
