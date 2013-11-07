import FWCore.ParameterSet.Types as CfgTypes
import FWCore.ParameterSet.Config as cms



process = cms.Process("Step3")


process.fwliteInput = cms.PSet(

    pathToFile    = cms.string("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store//user/bianchi/HBB_EDMNtuple/AllHDiJetPt_Step2/"),
    #pathToFile    = cms.string("/scratch/bianchi/HBB_EDMNtuple/All.H.DiJetPt/"),
    outPath       = cms.string("gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/"),
    #outPath       = cms.string("/scratch/bianchi/HBB_EDMNtuple/All.H.DiJetPt/"),
    ordering      = cms.string("DiJetPt_"),
    newDir        = cms.string("v2"),
    #ordering      = cms.string("ZllH.DiJetPt.Oct22."),
    lumi          = cms.double(12.1),
    verbose       = cms.bool(False),

    skims         = cms.VPSet(

    cms.PSet(
    name = cms.string("VType0"),
    cut  = cms.string("numJets30bTag>=1 && Vtype==0 && numJets30>=2 && jet1.pt>30 && jet2.pt>30")
    ),

    cms.PSet(
    name = cms.string("VType1"),
    cut  = cms.string("numJets30bTag>=1 && Vtype==1 && numJets30>=2 && jet1.pt>30 && jet2.pt>30")
    ),

    cms.PSet(
    name = cms.string("VType2"),
    cut  = cms.string("numJets30bTag>=1 && Vtype==2 && numJets30>=3 && jet1.pt>30 && jet2.pt>30 && jet3.pt>30")
    ),

    cms.PSet(
    name = cms.string("VType3"),
    cut  = cms.string("numJets30bTag>=1 && Vtype==3 && numJets30>=3 && jet1.pt>30 && jet2.pt>30 && jet3.pt>30")
    ),


    ),
    

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
