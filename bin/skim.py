import FWCore.ParameterSet.Types as CfgTypes
import FWCore.ParameterSet.Config as cms


# script that skims the Step2 trees and save the skimmed trees (and histos) into new files

##########################################
#
#  INPUT:     pathToFile/ordering+(samples.name).root
#
#    |
#    |
#    |
#    V
#
#  OUTPUT:    outPath_(skims.name)/newDir/ordering+(samples.name)_(skims.name).root
#
##########################################


# process name
process = cms.Process("Step3")

# main part
process.fwliteInput = cms.PSet(

    # PFN of the input directory (files in there need to be opneded via TFile::Open()
    pathToFile    = cms.string("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store//user/bianchi/HBB_EDMNtuple/AllHDiJetPt_V2/"),

    # PFN of the output directory (no '/' at the end, an extra name tag is appended depending on 'name')
    outPath       = cms.string("gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_V2"),

    # the name tag at the beginning of the files
    ordering      = cms.string("DiJetPt_"),

    # an extra sub-directory inside 'outPath' (optional)
    newDir        = cms.string(""),

    # the target luminosity (needed to normalize MCs)
    lumi          = cms.double(18.9),

    # verbosity (to debug the opening of files)
    verbose       = cms.bool(False),

    # set of skims to be run
    skims         = cms.VPSet(

    cms.PSet(
    name = cms.string("VType0"),
    cut  = cms.string("Vtype==0 ")
    ),

    cms.PSet(
    name = cms.string("VType1"),
    cut  = cms.string("Vtype==1")
    ),

    cms.PSet(
    name = cms.string("VType2"),
    cut  = cms.string("Vtype==2")
    ),

    cms.PSet(
    name = cms.string("VType3"),
    cut  = cms.string("Vtype==3")
    ),


    ),
    
    # the samples to process
    samples       = cms.VPSet(

    cms.PSet(
    skip     = cms.bool(True),
    name     = cms.string('DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball'),
    nickName = cms.string('DYJets'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71),
    update   = cms.bool(False)
    ),
    
    
    ),
    

    )
