import FWCore.ParameterSet.Types as CfgTypes
import FWCore.ParameterSet.Config as cms

process = cms.Process("TF")

process.fwliteInput = cms.PSet(

    # the tree with jet kinematics and gen event information (pt_gen -> target)
    pathToFile    = cms.string("./root/treeProducerTESTpt_gen.root"),

    #### the name of the out file with the RooWorkspace ###
    #
    # --> the actual name will be outFileName+'TEST'+argv[2]+'.root'
    #     if argv[2] = _reg, the regressed energy will be used
    #
    ####    
    outFileName   = cms.string("transferFunctions"),

    # verbosity flag
    verbose       = cms.bool(False),

    # if 1, compute only jet TF
    doOnlyJetTF   = cms.untracked.int32(0),

    # parameters of the 2G parametrization for b-jets
    relWeightBin0 = cms.untracked.double(0.65),
    relWeightBin1 = cms.untracked.double(0.65),
    relShiftBin0  = cms.untracked.double(0.92),
    
    )
