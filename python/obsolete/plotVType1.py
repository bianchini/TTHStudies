import FWCore.ParameterSet.Types as CfgTypes
import FWCore.ParameterSet.Config as cms

from Bianchi.TTHStudies.jetsBasic_cff import *
from Bianchi.TTHStudies.aNNInputs_cff import *

VType = "_VType1"

mmBasic      = "(Vtype==1 && H.HiggsFlag==1 && vLepton_charge[0]*vLepton_charge[1]<0 && V.mass>20)"
mmKin        = "(vLepton_pt[0]>30 && vLepton_pt[1]>20 && abs(vLepton_eta[0])<2.1 && abs(vLepton_eta[1])<2.4)"
mmId         = "(vLepton_pfCorrIso[0]<0.10 && vLepton_pfCorrIso[1]<0.20)"
mmJetBasic   = "(numJets30>=2 && pt1>30 && pt2>30 && numJets30bTag>=1)"

common    = mmBasic+" && "+mmKin+" && "+mmId+" && "+mmJetBasic
cat2j2b   = "numJets30==2 && numJets30bTag==2"
cat3j3b   = "numJets30>=3 && numJets30bTag>=3"

dataCut      = "((EVENT.json == 1 || EVENT.run < 196532) && (triggerFlags[5]>0 || triggerFlags[6]>0))"


ttjetsLF = \
         "nSimBs==2 && nC-nCTop==0"

ttjetsC  = \
        "nSimBs==2 && nC-nCTop>0"

ttjetsB  = \
        "nSimBs>2"

zjetsLF = \
        "nSimBs==0 && nC==0"

zjetsC  = \
       "nSimBs==0 && nC>0"

zjetsB  = \
       "nSimBs>0"

xsecTT_FH = 106.9
xsecTT_SL = 103.0
xsecTT_FL = 24.8

process = cms.Process("FWLitePlots")


process.fwliteInput = cms.PSet(

    pathToFile    = cms.string("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store//user/bianchi/HBB_EDMNtuple/AllHDiJetPt"+VType),
    #pathToFile    = cms.string("/scratch/bianchi/HBB_EDMNtuple/Zll.H.DiJetPt/"),
    ordering      = cms.string("DiJetPt_"),
    #ordering      = cms.string("ZllH.DiJetPt.Oct22."),
    lumi          = cms.double(12.1),
    debug         = cms.bool(False),
    


    samples       = cms.VPSet(

    cms.PSet(
    skip     = cms.bool(False),
    name     = cms.string('DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph'+VType),
    nickName = cms.string('DYJets10to50'),
    color    = cms.int32(18),
    xSec     = cms.double(12765.)
    ),
 
    cms.PSet(
    skip     = cms.bool(False),
    name     = cms.string('DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball'+VType),
    nickName = cms.string('DYJets50'),
    color    = cms.int32(18),
    xSec     = cms.double(3503.71)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball'+VType),
    nickName = cms.string('WJets'),
    color    = cms.int32(29),
    xSec     = cms.double(37509.0),
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTJets_HadronicMGDecays_8TeV-madgraph-part'+VType),
    nickName = cms.string('TTJetsFullHad'),
    color    = cms.int32(41),
    xSec     = cms.double(xsecTT_FH),
    ),
    
    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('TTJets_HadronicMGDecays_8TeV-madgraph-part'+VType),
    nickName = cms.string('TTJetsFullHad_LF'),
    color    = cms.int32(41),
    xSec     = cms.double(xsecTT_FH),
    cut      = cms.string(ttjetsLF),
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('TTJets_HadronicMGDecays_8TeV-madgraph-part'+VType),
    nickName = cms.string('TTJetsFullHad_C'),
    color    = cms.int32(44),
    xSec     = cms.double(xsecTT_FH),
    cut      = cms.string(ttjetsC),
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('TTJets_HadronicMGDecays_8TeV-madgraph-part'+VType),
    nickName = cms.string('TTJetsFullHad_B'),
    color    = cms.int32(46),
    xSec     = cms.double(xsecTT_FH),
    cut      = cms.string(ttjetsB),
    ),
    

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTJets_FullLeptMGDecays_8TeV-madgraph-part'+VType),
    nickName = cms.string('TTJetsFullLept'),
    color    = cms.int32(41),
    xSec     = cms.double(xsecTT_FL),
    ),
    
    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('TTJets_FullLeptMGDecays_8TeV-madgraph-part'+VType),
    nickName = cms.string('TTJetsFullLept_LF'),
    color    = cms.int32(41),
    xSec     = cms.double(xsecTT_FL),
    cut      = cms.string(ttjetsLF),
    ),
    
    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('TTJets_FullLeptMGDecays_8TeV-madgraph-part'+VType),
    nickName = cms.string('TTJetsFullLept_C'),
    color    = cms.int32(44),
    xSec     = cms.double(xsecTT_FL),
    cut      = cms.string(ttjetsC)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('TTJets_FullLeptMGDecays_8TeV-madgraph-part'+VType),
    nickName = cms.string('TTJetsFullLept_B'),
    color    = cms.int32(46),
    xSec     = cms.double(xsecTT_FL),
    cut      = cms.string(ttjetsB)
    ),
    
    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTJets_SemiLeptMGDecays_8TeV-madgraph-part'+VType),
    nickName = cms.string('TTJetsSemiLept'),
    color    = cms.int32(41),
    xSec     = cms.double(xsecTT_SL),
    ),
    
    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('TTJets_SemiLeptMGDecays_8TeV-madgraph-part'+VType),
    nickName = cms.string('TTJetsSemiLept_LF'),
    color    = cms.int32(41),
    xSec     = cms.double(xsecTT_SL),
    cut      = cms.string(ttjetsLF),
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('TTJets_SemiLeptMGDecays_8TeV-madgraph-part'+VType),
    nickName = cms.string('TTJetsSemiLept_C'),
    color    = cms.int32(44),
    xSec     = cms.double(xsecTT_SL),
    cut      = cms.string(ttjetsC),
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('TTJets_SemiLeptMGDecays_8TeV-madgraph-part'+VType),
    nickName = cms.string('TTJetsSemiLept_B'),
    color    = cms.int32(46),
    xSec     = cms.double(xsecTT_SL),
    cut      = cms.string(ttjetsB),
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('TTWJets_8TeV-madgraph'+VType),
    nickName = cms.string('TTW'),
    color    = cms.int32(18),
    xSec     = cms.double(0.232),
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('TTZJets_8TeV-madgraph'+VType),
    nickName = cms.string('TTZ'),
    color    = cms.int32(18),
    xSec     = cms.double(0.2057),
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola'+VType),
    nickName = cms.string('TtW'),
    color    = cms.int32(6),
    xSec     = cms.double(11.1)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('T_t-channel_TuneZ2star_8TeV-powheg-tauola'+VType),
    nickName = cms.string('Tt'),
    color    = cms.int32(6),
    xSec     = cms.double(56.4)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('T_s-channel_TuneZ2star_8TeV-powheg-tauola'+VType),
    nickName = cms.string('Ts'),
    color    = cms.int32(6),
    xSec     = cms.double(3.79)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola'+VType),
    nickName = cms.string('TbartW'),
    color    = cms.int32(6),
    xSec     = cms.double(11.1)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola'+VType),
    nickName = cms.string('Tbart'),
    color    = cms.int32(6),
    xSec     = cms.double(30.7)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola'+VType),
    nickName = cms.string('Tbars'),
    color    = cms.int32(6),
    xSec     = cms.double(1.76)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('WW_TuneZ2star_8TeV_pythia6_tauola'+VType),
    nickName = cms.string('WW'),
    color    = cms.int32(4),
    xSec     = cms.double(56.75)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('WZ_TuneZ2star_8TeV_pythia6_tauola'+VType),
    nickName = cms.string('WZ'),
    color    = cms.int32(4),
    xSec     = cms.double(33.85)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('ZZ_TuneZ2star_8TeV_pythia6_tauola'+VType),
    nickName = cms.string('ZZ'),
    color    = cms.int32(4),
    xSec     = cms.double(8.297)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('ZH_ZToLL_HToBB_M-110_8TeV-powheg-herwigpp'+VType),
    nickName = cms.string('ZH110'),
    color    = cms.int32(2),
    xSec     = cms.double(0.04414)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('ZH_ZToLL_HToBB_M-115_8TeV-powheg-herwigpp'+VType),
    nickName = cms.string('ZH115'),
    color    = cms.int32(2),
    xSec     = cms.double(0.036375)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('ZH_ZToLL_HToBB_M-120_8TeV-powheg-herwigpp'+VType),
    nickName = cms.string('ZH120'),
    color    = cms.int32(2),
    xSec     = cms.double(0.0293327854)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('ZH_ZToLL_HToBB_M-125_8TeV-powheg-herwigpp'+VType),
    nickName = cms.string('ZH125'),
    color    = cms.int32(2),
    xSec     = cms.double(0.0229727058)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('ZH_ZToLL_HToBB_M-130_8TeV-powheg-herwigpp'+VType),
    nickName = cms.string('ZH130'),
    color    = cms.int32(2),
    xSec     = cms.double(0.017288657)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('ZH_ZToLL_HToBB_M-135_8TeV-powheg-herwigpp'+VType),
    nickName = cms.string('ZH135'),
    color    = cms.int32(2),
    xSec     = cms.double(0.01250888)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('WH_WToLNu_HToBB_M-115_8TeV-powheg-herwigpp'+VType),
    nickName = cms.string('WH115'),
    color    = cms.int32(2),
    xSec     = cms.double(0.644299499),
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('WH_WToLNu_HToBB_M-120_8TeV-powheg-herwigpp'+VType),
    nickName = cms.string('WH120'),
    color    = cms.int32(2),
    xSec     = cms.double(0.5161968)
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('WH_WToLNu_HToBB_M-125_8TeV-powheg-herwigpp'+VType),
    nickName = cms.string('WH125'),
    color    = cms.int32(2),
    xSec     = cms.double(0.4019381),
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('WH_WToLNu_HToBB_M-130_8TeV-powheg-herwigpp'+VType),
    nickName = cms.string('WH130'),
    color    = cms.int32(2),
    xSec     = cms.double(0.3010930)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('WH_WToLNu_HToBB_M-135_8TeV-powheg-herwigpp'+VType),
    nickName = cms.string('WH135'),
    color    = cms.int32(2),
    xSec     = cms.double(0.230628)
    ),


    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTH_HToBB_M-110_8TeV-pythia6'+VType),
    nickName = cms.string('TTH110'),
    color    = cms.int32(2),
    xSec     = cms.double(0.1887*0.744)
    ),

    
    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTH_HToBB_M-115_8TeV-pythia6'+VType),
    nickName = cms.string('TTH115'),
    color    = cms.int32(2),
    xSec     = cms.double(0.1663*0.703)
    ),


    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTH_HToBB_M-120_8TeV-pythia6'+VType),
    nickName = cms.string('TTH120'),
    color    = cms.int32(2),
    xSec     = cms.double(0.1470*0.648)
    ),

    
    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('TTH_HToBB_M-125_8TeV-pythia6'+VType),
    nickName = cms.string('TTH125'),
    color    = cms.int32(2),
    xSec     = cms.double(0.1302*0.569)
    ),

    
    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTH_HToBB_M-130_8TeV-pythia6'+VType),
    nickName = cms.string('TTH130'),
    color    = cms.int32(2),
    xSec     = cms.double(0.1157*0.494)
    ),

    
    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTH_HToBB_M-135_8TeV-pythia6'+VType),
    nickName = cms.string('TTH135'),
    color    = cms.int32(2),
    xSec     = cms.double(0.1031*0.404)
    ),


    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('SingleElectronRun2012AAug06EdmV42'+VType),
    nickName = cms.string('DataSingleElectron_1'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    cut      = cms.string(dataCut),
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('SingleElectronRun2012AJul13EdmV42b'+VType),
    nickName = cms.string('DataSingleElectron_2'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    cut      = cms.string(dataCut),
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('SingleElectronRun2012BJul13EdmV42'+VType),
    nickName = cms.string('DataSingleElectron_3'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    cut      = cms.string(dataCut),
    ),

     cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('SingleElectronRun2012CAug24RerecoEdmV42'+VType),
    nickName = cms.string('DataSingleElectron_4'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    cut      = cms.string(dataCut),
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('SingleElectronRun2012CPromptv2EdmV42'+VType),
    nickName = cms.string('DataSingleElectron_5'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    cut      = cms.string(dataCut),
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('SingleElectronRun2012CPromptV2TopUpEdmV42'+VType),
    nickName = cms.string('DataSingleElectron_6'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    cut      = cms.string(dataCut),
    ),
    

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('DoubleElectron_Run2012A-recover-06Aug2012-v1_ProcV2'+VType),
    nickName = cms.string('DataDoubleElectron_1'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    cut      = cms.string(dataCut),
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('DoubleElectron_Run2012A-13Jul2012-v1_ProcFIXED'+VType),
    nickName = cms.string('DataDoubleElectron_2'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    cut      = cms.string(dataCut),
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('DoubleElectron_Run2012B-13Jul2012-v1_ProcFIXED'+VType),
    nickName = cms.string('DataDoubleElectron_3'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    cut      = cms.string(dataCut),
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('DoubleElectronRun2012CAug24RerecoEdmV42'+VType),
    nickName = cms.string('DataDoubleElectron_4'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    cut      = cms.string(dataCut),
    ),
    
    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('DoubleElectron_Run2012C-PromptReco-v2_HBB_EDMNtupleV42_ProcV1'+VType),
    nickName = cms.string('DataDoubleElectron_5'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    cut      = cms.string(dataCut),
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('DoubleElectron_Run2012C-PromptReco-v2_HBB_EDMNtupleV42_ProcV2'+VType),
    nickName = cms.string('DataDoubleElectron_6'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    cut      = cms.string(dataCut),
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
    xSec     = cms.double(-1),
    cut      = cms.string(dataCut)
    ),
    
    
    ),



    plots      = cms.VPSet(

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(10),
    xHigh     = cms.double(250),
    nBins     = cms.int32(120),
    variable  = cms.string("vLepton_pt[0]"),
    xTitle    = cms.string("p_{T} lead e"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_leadElePt"),
    cut       = cms.string(common),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(10),
    xHigh     = cms.double(250),
    nBins     = cms.int32(120),
    variable  = cms.string("vLepton_pt[1]"),
    xTitle    = cms.string("p_{T} trail e"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_trailElePt"),
    cut       = cms.string(common),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(30),
    nBins     = cms.int32(30),
    variable  = cms.string("nPVs"),
    xTitle    = cms.string("vertices"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_vertex"),
    cut       = cms.string(common),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(200),
    nBins     = cms.int32(50),
    variable  = cms.string("METtype1p2corr.et"),
    xTitle    = cms.string("MET"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_met"),
    cut       = cms.string(common),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(400),
    nBins     = cms.int32(100),
    variable  = cms.string("V.mass"),
    xTitle    = cms.string("di-e mass"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_diEleMass"),
    cut       = cms.string(common),
    logy      = cms.int32(1),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(-3),
    xHigh     = cms.double(3),
    nBins     = cms.int32(60),
    variable  = cms.string("vLepton_eta[0]"),
    xTitle    = cms.string("#eta lead e"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_leadEleEta"),
    cut       = cms.string(common),
    logy      = cms.int32(0),
    ),

    
    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(-3),
    xHigh     = cms.double(3),
    nBins     = cms.int32(60),
    variable  = cms.string("vLepton_eta[1]"),
    xTitle    = cms.string("#eta trail e"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType1_trailEleEta"),
    cut       = cms.string(common),
    logy      = cms.int32(0),
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

    )


for plot in plots_jetsBasic:
    newplot = plot.clone()
    newplot.histoName = cms.string("VType1_"+(plot.histoName).value()+"_preselection")
    newplot.cut       = cms.string(common)
    newplot.cutName   = cms.string("ee")
    oldVPset = process.fwliteInput.plots
    oldVPset.append(newplot)

for plot in plots_aNNInputs:
    oldVPset = process.fwliteInput.plots

    newplot = plot.clone()
    newplot.histoName = cms.string("VType1_"+(plot.histoName).value()+"_2j2b")
    newplot.cut       = cms.string(common+" && "+cat2j2b)
    newplot.cutName   = cms.string("ee, 2j2b")
    oldVPset.append(newplot)

    newplot2 = plot.clone()
    newplot2.histoName = cms.string("VType1_"+(plot.histoName).value()+"_3j3b")
    newplot2.cut       = cms.string(common+" && "+cat3j3b)
    newplot2.nBins     = cms.int32(10)
    newplot2.cutName   = cms.string("ee, #geq3j#geq3b")
    oldVPset.append(newplot2)


