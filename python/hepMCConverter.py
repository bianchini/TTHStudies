import FWCore.ParameterSet.Types as CfgTypes
import FWCore.ParameterSet.Config as cms

process = cms.Process("HEPMC")

process.fwliteInput = cms.PSet(

    # input file names
    pathToFile    = cms.vstring(
    #'./root/S_1.hepmc2g',
    '/shome/bianchi/Generators/Sherpa_run/ttH_LO/sample_LOPS_ttH_uweighted.hepmc2g'
    ),

    # output file name
    outFileName   = cms.string("./root/DiJetPt_TTH_HToBB_M-125_sherpa_uweighted.root"),

    # print out intermediate steps
    verbose       = cms.bool(False),

    # select only events passing lepton cut
    filter        = cms.bool(False),

    # are the tops and Higgs decayed ?
    higgsDecay    = cms.bool(True),
    topDecay      = cms.bool(True),

    # filter parameters on jet multiplicity
    ptCut         = cms.double(5.0),
    etaCut        = cms.double(5.0),
    nJetsMin      = cms.int32(0),
    
    # clean jets that overlap with leptons
    overlapLep    = cms.double(0.5),

    # lepton definition
    rIsoLep       = cms.double(0.10),
    dRIsoLep      = cms.double(0.5),
    ptCutLep      = cms.double(20),
    etaCutLep     = cms.double(2.5),

    # FastJet parameters
    jetRadius     = cms.double(0.4),

    # save onlt jets with pt in excess of...
    ptMin         = cms.double(20.0),

    )
