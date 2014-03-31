import FWCore.ParameterSet.Types as CfgTypes
import FWCore.ParameterSet.Config as cms

from Bianchi.TTHStudies.jetsBasic_cff import *
from Bianchi.TTHStudies.aNNInputs_cff import *

VType = "_VType2"

mmBasic      = "(V.pt>0 && V.Mt>40 && H.pt>0 && Vtype==2 && H.HiggsFlag==1)"
mmKin        = "(vLepton_pt[0]>30 && abs(vLepton_eta[0])<2.1)"
mmId         = "(vLepton_pfCorrIso[0]<0.10)"
mmJetBasic   = "(numJets30>=3 && numJets30bTag>=2)"

common       =  mmJetBasic + " && " + mmBasic+" && "+mmKin+" && "+mmId
#cat4j2b      = "numJets40==4 && numJets30bTag==2"
#cat4j3b      = "numJets40==4 && numJets30bTag==3"
#cat4j4b      = "numJets40==4 && numJets30bTag==4"
#cat5j2b      = "numJets40==5 && numJets30bTag==2"
#cat5j3b      = "numJets40==5 && numJets30bTag==3"
#cat5j4b      = "numJets40==5 && numJets30bTag>=4"
#cat6j2b      = "numJets40>=6 && numJets30bTag==2"
#cat6j3b      = "numJets40>=6 && numJets30bTag==3"
#cat6j4b      = "numJets40>=6 && numJets30bTag>=4"

cat4j2b      = "pt1>40 && pt2>40 && pt3>40 && pt4>30 && pt5<30 && numJets30bTag==2"
cat4j3b      = "pt1>40 && pt2>40 && pt3>40 && pt4>30 && pt5<30 && numJets30bTag==3"
cat4j4b      = "pt1>40 && pt2>40 && pt3>40 && pt4>30 && pt5<30 && numJets30bTag==4"
cat5j2b      = "pt1>40 && pt2>40 && pt3>40 && pt4>30 && pt5>30 && pt6<30 && numJets30bTag==2"
cat5j3b      = "pt1>40 && pt2>40 && pt3>40 && pt4>30 && pt5>30 && pt6<30 && numJets30bTag==3"
cat5j4b      = "pt1>40 && pt2>40 && pt3>40 && pt4>30 && pt5>30 && pt6<30 && numJets30bTag>=4"
cat6j2b      = "pt1>40 && pt2>40 && pt3>40 && pt4>30 && pt5>30 && pt6>30 && numJets30bTag==2"
cat6j3b      = "pt1>40 && pt2>40 && pt3>40 && pt4>30 && pt5>30 && pt6>30 && numJets30bTag==3"
cat6j4b      = "pt1>40 && pt2>40 && pt3>40 && pt4>30 && pt5>30 && pt6>30 && numJets30bTag>=4"



dataCut      = "((EVENT.json == 1 || EVENT.run < 196532) && (triggerFlags[14] || triggerFlags[15]  || triggerFlags[21]  || triggerFlags[22] || triggerFlags[23]))"

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
    #pathToFile    = cms.string("/scratch/bianchi/HBB_EDMNtuple/All.H.DiJetPt/"),
    ordering      = cms.string("DiJetPt_"),
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
    color    = cms.int32(19),
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
    name     = cms.string('WJetsToLNu_PtW-70To100_TuneZ2star_8TeV-madgraph'+VType),
    nickName = cms.string('WJetsPt70100'),
    color    = cms.int32(29),
    xSec     = cms.double(37509.0*0.0152),
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('WJetsToLNu_PtW-100_TuneZ2star_8TeV-madgraph'+VType),
    nickName = cms.string('WJetsPt100'),
    color    = cms.int32(29),
    xSec     = cms.double(37509.0*0.00883),
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
    name     = cms.string('WH_WToLNu_HToBB_M-110_8TeV-powheg-herwigpp'+VType),
    nickName = cms.string('WH110'),
    color    = cms.int32(2),
    xSec     = cms.double(0.78864),
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
    skip     = cms.bool(False),  
    name     = cms.string('SingleMuRun2012AAug06EdmV42'+VType),
    nickName = cms.string('DataSingleMu_1'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    cut      = cms.string(dataCut),
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('SingleMuRun2012AJul13EdmV42'+VType),
    nickName = cms.string('DataSingleMu_2'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    cut      = cms.string(dataCut),
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('SingleMuRun2012BJul13EdmV42'+VType),
    nickName = cms.string('DataSingleMu_3'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    cut      = cms.string(dataCut),
    ),

     cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('SingleMuRun2012CAug24RerecoEdmV42'+VType),
    nickName = cms.string('DataSingleMu_4'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    cut      = cms.string(dataCut),
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('SingleMuRun2012CPromptv2EdmV42'+VType),
    nickName = cms.string('DataSingleMu_5'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    cut      = cms.string(dataCut),
    ),

    cms.PSet(
    skip     = cms.bool(False),  
    name     = cms.string('SingleMuRun2012CPromptV2TopUpEdmV42'+VType),
    nickName = cms.string('DataSingleMu_6'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    cut      = cms.string(dataCut),
    ),

    
    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('SingleMuRun2012'),
    nickName = cms.string('DataSingleMu'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    cut      = cms.string(dataCut),
    ),

    
    ),



    plots      = cms.VPSet(

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(10),
    xHigh     = cms.double(250),
    nBins     = cms.int32(80),
    variable  = cms.string("vLepton_pt[0]"),
    xTitle    = cms.string("p_{T} lead #mu"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType2_leadMuPt"),
    cut       = cms.string(common),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(0.1),
    nBins     = cms.int32(40),
    variable  = cms.string("vLepton_pfCorrIso[0]"),
    xTitle    = cms.string("iso"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType2_iso"),
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
    histoName = cms.string("VType2_vertex"),
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
    histoName = cms.string("VType2_met"),
    cut       = cms.string(common),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(-3),
    xHigh     = cms.double(3),
    nBins     = cms.int32(60),
    variable  = cms.string("vLepton_eta[0]"),
    xTitle    = cms.string("#eta lead #mu"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("VType2_leadMuEta"),
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
    newplot.skip      = cms.bool(True) 
    newplot.histoName = cms.string("VType2_"+(plot.histoName).value()+"_preselection")
    newplot.cut       = cms.string(common)
    newplot.cutName   = cms.string("Single #mu")
    oldVPset = process.fwliteInput.plots
    oldVPset.append(newplot)

for plot in plots_aNNInputs:

    oldVPset = process.fwliteInput.plots

    newplot = plot.clone()
    newplot.skip      = cms.bool(False) 
    newplot.histoName  = cms.string("VType2_"+(plot.histoName).value()+"_4j2b")
    newplot.cut        = cms.string(common+" && "+cat4j2b)
    newplot.cutName    = cms.string("Single #mu, 4j2b")
    oldVPset.append(newplot)

    newplot2 = plot.clone()
    newplot2.skip      = cms.bool(False) 
    newplot2.histoName = cms.string("VType2_"+(plot.histoName).value()+"_4j3b")
    newplot2.cut       = cms.string(common+" && "+cat4j3b)
    newplot2.cutName   = cms.string("Single #mu, 4j3b")
    oldVPset.append(newplot2)

    newplot3 = plot.clone()
    newplot3.skip      = cms.bool(False) 
    newplot3.histoName = cms.string("VType2_"+(plot.histoName).value()+"_4j4b")
    newplot3.cut       = cms.string(common+" && "+cat4j4b)
    newplot3.cutName   = cms.string("Single #mu, 4j4b")
    oldVPset = process.fwliteInput.plots
    oldVPset.append(newplot3)

    newplot4 = plot.clone()
    newplot4.skip      = cms.bool(False) 
    newplot4.histoName = cms.string("VType2_"+(plot.histoName).value()+"_5j2b")
    newplot4.cut       = cms.string(common+" && "+cat5j2b)
    newplot4.cutName   = cms.string("Single #mu, 5j2b")
    oldVPset = process.fwliteInput.plots
    oldVPset.append(newplot4)

    newplot5 = plot.clone()
    newplot5.skip      = cms.bool(False) 
    newplot5.histoName = cms.string("VType2_"+(plot.histoName).value()+"_5j3b")
    newplot5.cut       = cms.string(common+" && "+cat5j3b)
    newplot5.cutName   = cms.string("Single #mu, 5j3b")
    oldVPset.append(newplot5)

    newplot6 = plot.clone()
    newplot6.skip      = cms.bool(False) 
    newplot6.histoName = cms.string("VType2_"+(plot.histoName).value()+"_5j4b")
    newplot6.cut       = cms.string(common+" && "+cat5j4b)
    newplot6.cutName   = cms.string("Single #mu, 5j4b")
    oldVPset = process.fwliteInput.plots
    oldVPset.append(newplot6)

    newplot7 = plot.clone()
    newplot7.skip      = cms.bool(False) 
    newplot7.histoName = cms.string("VType2_"+(plot.histoName).value()+"_6j2b")
    newplot7.cut       = cms.string(common+" && "+cat6j2b)
    newplot7.cutName   = cms.string("Single #mu, 6j2b")
    oldVPset = process.fwliteInput.plots
    oldVPset.append(newplot7)

    newplot8 = plot.clone()
    newplot8.skip      = cms.bool(False) 
    newplot8.histoName = cms.string("VType2_"+(plot.histoName).value()+"_6j3b")
    newplot8.cut       = cms.string(common+" && "+cat6j3b)
    newplot8.cutName   = cms.string("Single #mu, 6j3b")
    oldVPset.append(newplot8)

    newplot9 = plot.clone()
    newplot9.skip      = cms.bool(False) 
    newplot9.histoName = cms.string("VType2_"+(plot.histoName).value()+"_6j4b")
    newplot9.cut       = cms.string(common+" && "+cat6j4b)
    newplot9.cutName   = cms.string("Single #mu, 6j4b")
    oldVPset = process.fwliteInput.plots
    oldVPset.append(newplot9)
