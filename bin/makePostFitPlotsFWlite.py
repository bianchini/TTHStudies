import FWCore.ParameterSet.Config as cms

process = cms.Process("makePostFitPlotsFWlite")

process.fwliteInput = cms.PSet(

    path2Workspace  = cms.string("../../../HiggsAnalysis/CombinedLimit/test/mytest.root"),
    path2FitResults = cms.string("../../../HiggsAnalysis/CombinedLimit/test/mlfitMEM_COMB_New_rec_std_sb_wplots.root "),
    outputName      = cms.string("Test.root"),
    dirName         = cms.string("test"),

)
