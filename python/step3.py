import FWCore.ParameterSet.Types as CfgTypes
import FWCore.ParameterSet.Config as cms



process = cms.Process("Step3")


process.fwliteInput = cms.PSet(

    #pathToFile    = cms.string("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/HBB_EDMNtuple/V42/Oct22/env/sys/MVAout/"),
    #pathToFile    = cms.string("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store//user/bianchi/HBB_EDMNtuple/ZllHDiJetPt/")
    pathToFile    = cms.string("/scratch/bianchi/HBB_EDMNtuple/All.H.DiJetPt/"),
    #ordering      = cms.string("ZllH.DiJetPt.Oct22."),
    ordering      = cms.string("DiJetPt_"),
    lumi          = cms.double(12.1),
    verbose       = cms.bool(False),
    

    samples       = cms.VPSet(


    cms.PSet(
    skip     = cms.bool(False),
    name     = cms.string('DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball'),
    nickName = cms.string('DYJets'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71),
    update   = cms.bool(True)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('DY1JetsToLL_M-50_TuneZ2Star_8TeV-madgraph'),
    nickName = cms.string('DYJets1'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71),
    ),

    cms.PSet(
    skip     = cms.bool(True),
    name     = cms.string('DY2JetsToLL_M-50_TuneZ2Star_8TeV-madgraph'),
    nickName = cms.string('DYJets2'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('DY3JetsToLL_M-50_TuneZ2Star_8TeV-madgraph'),
    nickName = cms.string('DYJets3'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('DY4JetsToLL_M-50_TuneZ2Star_8TeV-madgraph'),
    nickName = cms.string('DYJets4'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph'),
    nickName = cms.string('DYJetsHT1'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph'),
    nickName = cms.string('DYJetsHT2'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71)
    ),

    ##cms.PSet(
    ##name     = cms.string('DYJetsToLL_PtZ-180_TuneZ2star_8TeV-madgraph'),
    ##nickName = cms.string('DYJetsPt4'),
    ##color    = cms.int32(18),
    ##xSec     = cms.double(3503.71)
    ##),
    
    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('DYJetsToLL_PtZ-100_TuneZ2star_8TeV-madgraph'),
    nickName = cms.string('DYJetsPt3'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71)
    ),
    
    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('DYJetsToLL_PtZ-70To100_TuneZ2star_8TeV-madgraph-tarball'),
    nickName = cms.string('DYJetsPt2'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('DYJetsToLL_PtZ-50To70_TuneZ2star_8TeV-madgraph-tarball'),
    nickName = cms.string('DYJetsPt1'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTJets_Merged'),
    nickName = cms.string('TTJets'),
    color    = cms.int32(5),
    xSec     = cms.double(234),
    update   = cms.bool(True),
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('TTJets_HadronicMGDecays_8TeV-madgraph-part'),
    nickName = cms.string('TTJetsFullHad'),
    color    = cms.int32(5),
    xSec     = cms.double(133.62),
    update   = cms.bool(True),
    ),
    
    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('TTJets_FullLeptMGDecays_8TeV-madgraph-part'),
    nickName = cms.string('TTJetsFullLept'),
    color    = cms.int32(5),
    xSec     = cms.double(24.56),
    update   = cms.bool(True),
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTJets_SemiLeptMGDecays_8TeV-madgraph-part'),
    nickName = cms.string('TTJetsSemiLept'),
    color    = cms.int32(5),
    xSec     = cms.double(75.82),
    update   = cms.bool(True),
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola'),
    nickName = cms.string('TtW'),
    color    = cms.int32(6),
    xSec     = cms.double(11.1),
    update   = cms.bool(True)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('T_t-channel_TuneZ2star_8TeV-powheg-tauola'),
    nickName = cms.string('Tt'),
    color    = cms.int32(6),
    xSec     = cms.double(56.4),
    update   = cms.bool(True)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('T_s-channel_TuneZ2star_8TeV-powheg-tauola'),
    nickName = cms.string('Ts'),
    color    = cms.int32(6),
    xSec     = cms.double(3.79),
    update   = cms.bool(True)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola'),
    nickName = cms.string('TbartW'),
    color    = cms.int32(6),
    xSec     = cms.double(11.1),
    update   = cms.bool(True)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola'),
    nickName = cms.string('Tbart'),
    color    = cms.int32(6),
    xSec     = cms.double(30.7),
    update   = cms.bool(True)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola'),
    nickName = cms.string('Tbars'),
    color    = cms.int32(6),
    xSec     = cms.double(1.76),
    update   = cms.bool(True),
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('WW_TuneZ2star_8TeV_pythia6_tauola'),
    nickName = cms.string('WW'),
    color    = cms.int32(4),
    xSec     = cms.double(56.75),
    update   = cms.bool(True)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('WZ_TuneZ2star_8TeV_pythia6_tauola'),
    nickName = cms.string('WZ'),
    color    = cms.int32(3),
    xSec     = cms.double(33.85),
    update   = cms.bool(True)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('ZZ_TuneZ2star_8TeV_pythia6_tauola'),
    nickName = cms.string('ZZ'),
    color    = cms.int32(8),
    xSec     = cms.double(8.297),
    update   = cms.bool(True)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('ZH_ZToLL_HToBB_M-110_8TeV-powheg-herwigpp'),
    nickName = cms.string('ZH110'),
    color    = cms.int32(2),
    xSec     = cms.double(0.04414),
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('ZH_ZToLL_HToBB_M-115_8TeV-powheg-herwigpp'),
    nickName = cms.string('ZH115'),
    color    = cms.int32(2),
    xSec     = cms.double(0.036375),
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('ZH_ZToLL_HToBB_M-120_8TeV-powheg-herwigpp'),
    nickName = cms.string('ZH120'),
    color    = cms.int32(2),
    xSec     = cms.double(0.0293327854)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('ZH_ZToLL_HToBB_M-125_8TeV-powheg-herwigpp'),
    nickName = cms.string('ZH125'),
    color    = cms.int32(2),
    xSec     = cms.double(0.0229727058),
    update   = cms.bool(True)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('ZH_ZToLL_HToBB_M-130_8TeV-powheg-herwigpp'),
    nickName = cms.string('ZH130'),
    color    = cms.int32(2),
    xSec     = cms.double(0.017288657)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('ZH_ZToLL_HToBB_M-135_8TeV-powheg-herwigpp'),
    nickName = cms.string('ZH135'),
    color    = cms.int32(2),
    xSec     = cms.double(0.01250888)
    ),


    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('WH_WToLNu_HToBB_M-110_8TeV-powheg-herwigpp'),
    nickName = cms.string('WH110'),
    color    = cms.int32(2),
    xSec     = cms.double(0.78864),
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('WH_WToLNu_HToBB_M-115_8TeV-powheg-herwigpp'),
    nickName = cms.string('WH115'),
    color    = cms.int32(2),
    xSec     = cms.double(0.644299499),
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('WH_WToLNu_HToBB_M-120_8TeV-powheg-herwigpp'),
    nickName = cms.string('WH120'),
    color    = cms.int32(2),
    xSec     = cms.double(0.5161968)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('WH_WToLNu_HToBB_M-125_8TeV-powheg-herwigpp'),
    nickName = cms.string('WH125'),
    color    = cms.int32(2),
    xSec     = cms.double(0.4019381),
    update   = cms.bool(True)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('WH_WToLNu_HToBB_M-130_8TeV-powheg-herwigpp'),
    nickName = cms.string('WH130'),
    color    = cms.int32(2),
    xSec     = cms.double(0.3010930)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('WH_WToLNu_HToBB_M-135_8TeV-powheg-herwigpp'),
    nickName = cms.string('WH135'),
    color    = cms.int32(2),
    xSec     = cms.double(0.230628)
    ),








    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('DataZmm'),
    nickName = cms.string('DataZmm'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    update   = cms.bool(True)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('DataZee'),
    nickName = cms.string('DataZee'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    update   = cms.bool(True)
    ),
    
    
    ),


    )
