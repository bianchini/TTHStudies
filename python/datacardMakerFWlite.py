import FWCore.ParameterSet.Config as cms

from Bianchi.TTHStudies.mem_categories_cff import *

process = cms.Process("datacardMaker")

process.fwliteInput = cms.PSet(

    name = cms.string("New"),
    version= cms.string("_CSVcalibration_rec_std"),
    extraname= cms.string(""),
    fname= cms.string("MEM"),
    inputpath= cms.string("../root/"),
    directory= cms.string("Mar25_2014"),
    cut= cms.string("(numJets>=6 && numBTagM==2 && Vtype==2)"),
    category= cms.string("lepton_pt"),
    doMEM= cms.int32(4),
    fact1= cms.double(-99),
    fact2= cms.double(-99),
    factbb= cms.double(-99),
    lumiScale= cms.double(19.6/12.1),
    nBins= cms.int32(15),
    splitFirstBin= cms.int32(0),
    binvec= cms.vdouble(30., 40, 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180.),   
    #samples= cms.vstring( "TTV", "SingleT", "DiBoson", "TTJetsBB",
    #                      "TTJetsBJ", "TTJetsJJ", "TTH125", "EWK",
    #                      "Run2012_SingleMu", "Run2012_SingleElectron"),
    samples= cms.vstring("SingleT"),
    nparts= cms.int32(2),
    part  = cms.int32(1),
    )



#process.fwliteInput = cat1_sb
