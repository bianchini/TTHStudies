import FWCore.ParameterSet.Types as CfgTypes
import FWCore.ParameterSet.Config as cms

stitchDY = False

# trigger matching? ID? ISO?

mmBasic = "Vtype==0 && H.HiggsFlag==1 && vLepton_charge[0]*vLepton_charge[1]<0"
mmKin   = "vLepton_pt[0]>30 && vLepton_pt[1]>20 && abs(vLepton_eta[0])<2.1 && abs(vLepton_eta[1])<2.4"
mmId    = "vLepton_pfCorrIso[0]<0.10 && vLepton_pfCorrIso[1]<0.20"
#mmId   =
#"(vLepton_chargedHadIso[0]+TMath::Max((vLepton_neutralHadIso[0])+(vLepton_photonIso[0])-0.5*(vLepton_chargedPUIso[0]) ,0.0))/vLepton_pt[0]<0.1 &&" +
#"(vLepton_chargedHadIso[1]+TMath::Max((vLepton_neutralHadIso[1])+(vLepton_photonIso[1])-0.5*(vLepton_chargedPUIso[1]) ,0.0))/vLepton_pt[1]<0.2" 



eeBasic = "Vtype==1 && H.HiggsFlag==1 && vLepton_charge[0]*vLepton_charge[1]<0"
eeKin   = "vLepton_pt[0]>30 && vLepton_pt[1]>20 && abs(vLepton_eta[0])<2.5 && abs(vLepton_eta[1])<2.5"
eeId    = "vLepton_pfCorrIso[0]<0.10 && vLepton_pfCorrIso[1]<0.20"
#eeId   =
#"(vLepton_chargedHadIso[0]+TMath::Max((vLepton_neutralHadIso[0])+(vLepton_photonIso[0])-0.5*(vLepton_chargedPUIso[0]) ,0.0))/vLepton_pt[0]<0.1 &&" +
#"(vLepton_chargedHadIso[1]+TMath::Max((vLepton_neutralHadIso[1])+(vLepton_photonIso[1])-0.5*(vLepton_chargedPUIso[1]) ,0.0))/vLepton_pt[1]<0.2" 



process = cms.Process("FWLitePlots")

process.fwliteInput = cms.PSet(

    #pathToFile    = cms.string("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store/HBB_EDMNtuple/V42/Oct22/env/sys/MVAout/"),
    pathToFile    = cms.string("."),
    ordering      = cms.string("ZllH.DiJetPt.Oct22."),
    lumi          = cms.double(12.1),


    samples       = cms.VPSet(

    cms.PSet(
    skip     = cms.bool(True),
    name     = cms.string('DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball'),
    nickName = cms.string('DYJets'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71)
    ),

    cms.PSet(
    skip     = cms.bool(not stitchDY),
    name     = cms.string('DY1JetsToLL_M-50_TuneZ2Star_8TeV-madgraph'),
    nickName = cms.string('DYJets1'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71)
    ),

    cms.PSet(
    skip     = cms.bool(not stitchDY),
    name     = cms.string('DY2JetsToLL_M-50_TuneZ2Star_8TeV-madgraph'),
    nickName = cms.string('DYJets2'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71)
    ),

    cms.PSet(
    skip     = cms.bool(not stitchDY),  
    name     = cms.string('DY3JetsToLL_M-50_TuneZ2Star_8TeV-madgraph'),
    nickName = cms.string('DYJets3'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71)
    ),

    cms.PSet(
    skip     = cms.bool(not stitchDY),  
    name     = cms.string('DY4JetsToLL_M-50_TuneZ2Star_8TeV-madgraph'),
    nickName = cms.string('DYJets4'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71)
    ),

    cms.PSet(
    skip     = cms.bool(not stitchDY),  
    name     = cms.string('DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph'),
    nickName = cms.string('DYJetsHT1'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71)
    ),

    cms.PSet(
    skip     = cms.bool(not stitchDY),  
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
    skip     = cms.bool(not stitchDY),  
    name     = cms.string('DYJetsToLL_PtZ-100_TuneZ2star_8TeV-madgraph'),
    nickName = cms.string('DYJetsPt3'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71)
    ),
    
    cms.PSet(
    skip     = cms.bool(not stitchDY),  
    name     = cms.string('DYJetsToLL_PtZ-70To100_TuneZ2star_8TeV-madgraph-tarball'),
    nickName = cms.string('DYJetsPt2'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71)
    ),

    cms.PSet(
    skip     = cms.bool(not stitchDY),  
    name     = cms.string('DYJetsToLL_PtZ-50To70_TuneZ2star_8TeV-madgraph-tarball'),
    nickName = cms.string('DYJetsPt1'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('TTJets_Merged'),
    nickName = cms.string('TTJets'),
    color    = cms.int32(5),
    xSec     = cms.double(234)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola'),
    nickName = cms.string('TtW'),
    color    = cms.int32(6),
    xSec     = cms.double(11.1)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('T_t-channel_TuneZ2star_8TeV-powheg-tauola'),
    nickName = cms.string('Tt'),
    color    = cms.int32(6),
    xSec     = cms.double(56.4)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('T_s-channel_TuneZ2star_8TeV-powheg-tauola'),
    nickName = cms.string('Ts'),
    color    = cms.int32(6),
    xSec     = cms.double(3.79)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola'),
    nickName = cms.string('TbartW'),
    color    = cms.int32(6),
    xSec     = cms.double(11.1)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola'),
    nickName = cms.string('Tbart'),
    color    = cms.int32(6),
    xSec     = cms.double(30.7)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola'),
    nickName = cms.string('Tbars'),
    color    = cms.int32(6),
    xSec     = cms.double(1.76)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('WW_TuneZ2star_8TeV_pythia6_tauola'),
    nickName = cms.string('WW'),
    color    = cms.int32(4),
    xSec     = cms.double(56.75)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('WZ_TuneZ2star_8TeV_pythia6_tauola'),
    nickName = cms.string('WZ'),
    color    = cms.int32(3),
    xSec     = cms.double(33.85)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('ZZ_TuneZ2star_8TeV_pythia6_tauola'),
    nickName = cms.string('ZZ'),
    color    = cms.int32(8),
    xSec     = cms.double(8.297)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('ZH_ZToLL_HToBB_M-110_8TeV-powheg-herwigpp'),
    nickName = cms.string('ZH110'),
    color    = cms.int32(2),
    xSec     = cms.double(0.04414)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('ZH_ZToLL_HToBB_M-115_8TeV-powheg-herwigpp'),
    nickName = cms.string('ZH115'),
    color    = cms.int32(2),
    xSec     = cms.double(0.036375)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('ZH_ZToLL_HToBB_M-120_8TeV-powheg-herwigpp'),
    nickName = cms.string('ZH120'),
    color    = cms.int32(2),
    xSec     = cms.double(0.0293327854)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('ZH_ZToLL_HToBB_M-125_8TeV-powheg-herwigpp'),
    nickName = cms.string('ZH125'),
    color    = cms.int32(2),
    xSec     = cms.double(0.0229727058)
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
    name     = cms.string('DataZmm'),
    nickName = cms.string('DataZmm'),
    color    = cms.int32(1),
    xSec     = cms.double(-1)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('DataZee'),
    nickName = cms.string('DataZee'),
    color    = cms.int32(1),
    xSec     = cms.double(-1)
    ),
    
    
    ),



    plots      = cms.VPSet(

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(10),
    xHigh     = cms.double(250),
    nBins     = cms.int32(120),
    variable  = cms.string("vLepton_pt[0]"),
    xTitle    = cms.string("p_{T} lead #mu"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("Zmm_leadMuPt"),
    cut       = cms.string(mmBasic+" && "+mmKin+" && "+mmId+" && "+"pt1>30 && pt2>30"),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(10),
    xHigh     = cms.double(250),
    nBins     = cms.int32(120),
    variable  = cms.string("vLepton_pt[1]"),
    xTitle    = cms.string("p_{T} trail #mu"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("Zmm_trailMuPt"),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(0),
    xHigh     = cms.double(30),
    nBins     = cms.int32(30),
    variable  = cms.string("nPVs"),
    xTitle    = cms.string("vertices"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("Zmm_vertex"),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(0),
    xHigh     = cms.double(200),
    nBins     = cms.int32(50),
    variable  = cms.string("METtype1p2corr"),
    xTitle    = cms.string("MET"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("Zmm_met"),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(60),
    xHigh     = cms.double(200),
    nBins     = cms.int32(70),
    variable  = cms.string("V.mass"),
    xTitle    = cms.string("di#mu mass"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("Zmm_diMuMass"),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(-3),
    xHigh     = cms.double(3),
    nBins     = cms.int32(60),
    variable  = cms.string("vLepton_eta[0]"),
    xTitle    = cms.string("#eta lead #mu"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("Zmm_leadMuEta"),
    logy      = cms.int32(0),
    ),

    
    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(-3),
    xHigh     = cms.double(3),
    nBins     = cms.int32(60),
    variable  = cms.string("vLepton_eta[1]"),
    xTitle    = cms.string("#eta trail #mu"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("Zmm_trailMuEta"),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(20),
    xHigh     = cms.double(420),
    nBins     = cms.int32(100),
    variable  = cms.string("hJet_pt[0]"),
    xTitle    = cms.string("p_{T} lead jet"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("Zmm_leadJetPt"),
    logy      = cms.int32(1),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(20),
    xHigh     = cms.double(420),
    nBins     = cms.int32(100),
    variable  = cms.string("hJet_pt[1]"),
    xTitle    = cms.string("p_{T} trail jet"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("Zmm_trailJetPt"),
    logy      = cms.int32(1),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(0),
    xHigh     = cms.double(10),
    nBins     = cms.int32(10),
    variable  = cms.string("numJets"),
    xTitle    = cms.string("jet mult"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("Zmm_numJet"),
    logy      = cms.int32(1),
    ),

    cms.PSet(
    skip      = cms.bool(True),
    xLow      = cms.double(0),
    xHigh     = cms.double(10),
    nBins     = cms.int32(10),
    variable  = cms.string("numBJets"),
    xTitle    = cms.string("b-jet mult"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("Zmm_numBJet"),
    logy      = cms.int32(1),
    ),

    # cms.PSet(
    #xLow      = cms.double(),
    #xHigh     = cms.double(),
    #nBins     = cms.int32(),
    #xTitle    = cms.string(),
    #yTitle    = cms.string(),
    #histoName = cms.string(),
    #logy      = cms.int32(),
    #),

    ),


    fileList   = cms.vstring(
    'DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball',
    'DYJetsToLL_PtZ-50To70_TuneZ2star_8TeV-madgraph-tarball',
    'DYJetsToLL_PtZ-70To100_TuneZ2star_8TeV-madgraph-tarball',
    'DYJetsToLL_PtZ-100_TuneZ2star_8TeV-madgraph',
    'DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph',
    'DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph',
    'DY1JetsToLL_M-50_TuneZ2Star_8TeV-madgraph',
    'DY2JetsToLL_M-50_TuneZ2Star_8TeV-madgraph',
    'DY3JetsToLL_M-50_TuneZ2Star_8TeV-madgraph',
    'DY4JetsToLL_M-50_TuneZ2Star_8TeV-madgraph',
    'ZJetsToLL_Pt-100_8TeV-herwigpp',
    'TTJets_Merged',    
    'TToTlepWhad_tW-channel-DR_8TeV-powheg-tauola',
    'TToThadWlep_tW-channel-DR_8TeV-powheg-tauola',
    'TToDilepton_tW-channel-DR_8TeV-powheg-tauola',
    'TToLeptons_t-channel_8TeV-powheg-tauola',
    'TToLeptons_s-channel_8TeV-powheg-tauola',
    'T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola',
    'T_t-channel_TuneZ2star_8TeV-powheg-tauola',
    'T_s-channel_TuneZ2star_8TeV-powheg-tauola',
    'TBarToDilepton_tW-channel-DR_8TeV-powheg-tauola',
    'TBarToLeptons_t-channel_8TeV-powheg-tauola',
    'TBarToLeptons_s-channel_8TeV-powheg-tauola',   
    'Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola',
    'Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola',
    'Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola',
    'WW_TuneZ2star_8TeV_pythia6_tauola',
    'WZ_TuneZ2star_8TeV_pythia6_tauola',
    'WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola',
    'ZZ_TuneZ2star_8TeV_pythia6_tauola',
    'ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola',
    'ZH_ZToLL_HToBB_M-110_8TeV-powheg-herwigpp',
    'ZH_ZToLL_HToBB_M-115_8TeV-powheg-herwigpp',
    'ZH_ZToLL_HToBB_M-120_8TeV-powheg-herwigpp',
    'ZH_ZToLL_HToBB_M-125_8TeV-powheg-herwigpp',
    'ZH_ZToLL_HToBB_M-130_8TeV-powheg-herwigpp',
    'ZH_ZToLL_HToBB_M-135_8TeV-powheg-herwigpp',
    'DataZmm',
    'DataZee',
    ),
    
    xsection      = cms.vdouble()

    )
