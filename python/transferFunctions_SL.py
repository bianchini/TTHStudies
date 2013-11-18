import FWCore.ParameterSet.Types as CfgTypes
import FWCore.ParameterSet.Config as cms

process = cms.Process("TFuncSL")

process.fwliteInput = cms.PSet(

    pathToFile    = cms.string("./root/treeProducerTESTpt_gen.root"),
    outFileName   = cms.string("transferFunctionsTESTpt_gen.root"),
    verbose       = cms.bool(False),
    doOnlyJetTF   = cms.untracked.int32(1),

    relWeightBin0 = cms.untracked.double(0.65),
    relWeightBin1 = cms.untracked.double(0.65),
    relShiftBin0  = cms.untracked.double(0.92),
    
    )
