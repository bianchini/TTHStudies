import FWCore.ParameterSet.Types as CfgTypes
import FWCore.ParameterSet.Config as cms

process = cms.Process("HEPMC")

process.fwliteInput = cms.PSet(

    pathToFile    = cms.string("tryfull.hepmc2g"),
    outFileName   = cms.string("TEST.root"),
    verbose       = cms.bool(True),
    filter        = cms.bool(False),
    
    
    )
