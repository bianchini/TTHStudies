import FWCore.ParameterSet.Types as CfgTypes
import FWCore.ParameterSet.Config as cms

process = cms.Process("HEPMC")

process.fwliteInput = cms.PSet(

    pathToFile    = cms.vstring(
    './root/hep_1BJpu_999468319.hepmc2g',
    './root/hep_1BJpu_993036056.hepmc2g'),

    outFileName   = cms.string("TEST.root"),
    verbose       = cms.bool(False),
    filter        = cms.bool(False),        

    )
