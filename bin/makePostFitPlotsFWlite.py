import FWCore.ParameterSet.Config as cms

process = cms.Process("makePostFitPlotsFWlite")

process.fwliteInput = cms.PSet(

    path2Workspace  = cms.string("NONE"),
    path2Datacard   = cms.string("MEM_New_rec_std_cat6_H_gt05.txt"),
    path2FitResults = cms.string("../../../HiggsAnalysis/CombinedLimit/test/mlfitMEM_COMB_New_rec_std_sb_wplots.root "),
    outputName      = cms.string("PostFit_bj.root"),
    dirName         = cms.string("cat6_H_gt05"),

)
