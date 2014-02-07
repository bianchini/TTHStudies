import FWCore.ParameterSet.Types as CfgTypes
import FWCore.ParameterSet.Config as cms


inputDir = '/scratch/bianchi/HepMC/ttH_LO_UE_HAD/'
gen      = 'LOPSUEHAD'
proc     = 'ttH'


process = cms.Process("HEPMC")

process.fwliteInput = cms.PSet(

    # input file names
    pathToFile    = cms.vstring(
    #'./root/S_1.hepmc2g',
    inputDir+'/sample_'+gen+'_'+proc+'_unweighted.hepmc2g',
    #inputDir+'/sample_'+gen+'_'+proc+'_unweighted.1.hepmc2g',
    #inputDir+'/sample_'+gen+'_'+proc+'_unweighted.2.hepmc2g',
    #inputDir+'/sample_'+gen+'_'+proc+'_unweighted.3.hepmc2g',
    #inputDir+'/sample_'+gen+'_'+proc+'_unweighted.4.hepmc2g',
    #inputDir+'/sample_'+gen+'_'+proc+'_unweighted.5.hepmc2g',
    #inputDir+'/sample_'+gen+'_'+proc+'_unweighted.6.hepmc2g',
    #inputDir+'/sample_'+gen+'_'+proc+'_unweighted.7.hepmc2g',
    #inputDir+'/sample_'+gen+'_'+proc+'_unweighted.8.hepmc2g',
    #inputDir+'/sample_'+gen+'_'+proc+'_unweighted.9.hepmc2g',
    ),

    # output file name
    outFileName   = cms.string("/scratch/bianchi/HBB_EDMNtuple/Sherpa_run/DiJetPt_TTH125_sherpa_"+gen+"_unweighted.root"),
    #outFileName   = cms.string("./root/TEST.root"),

    # print out intermediate steps
    verbose       = cms.bool(True),

    # select only events passing lepton cut
    filter        = cms.bool(False),

    # are the tops and Higgs decayed ?
    higgsDecay    = cms.bool(True),
    topDecay      = cms.bool(True),

    # was the sample generated with parton shower ?
    shower        = cms.bool(True),
    
    # was the sample generated with fragmentation ?
    fragmentation = cms.bool(True),

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
    jetRadius     = cms.double(0.5),

    # save only jets with pt in excess of...
    ptMin         = cms.double(10.0),

    )
