import FWCore.ParameterSet.Types  as CfgTypes
import FWCore.ParameterSet.Config as cms

VType     = "_VType1"

xsecTT_SL = 103.0
xsecTT_FL = 24.8

process = cms.Process("MEAnalysisNew")

process.fwliteInput = cms.PSet(

    outFileName   = cms.string("./root/MEAnalysisNew.root"),
    pathToTF      = cms.string("./root/transferFunctionsNew_partonE_new.root"),
    pathToCP      = cms.string("./root/ControlPlotsNew_new.root"),
    pathToFile    = cms.string("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store//user/bianchi/HBB_EDMNtuple/AllHDiJetPt"+VType+"/v2/"),
    ordering      = cms.string("DiJetPt_"),
    lumi          = cms.double(12.1),

    samples       = cms.VPSet(

    cms.PSet(
    skip     = cms.bool(True),
    name     = cms.string('DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph'+VType),
    nickName = cms.string('DYJets10to50'),
    color    = cms.int32(18),
    xSec     = cms.double(12765.)
    ),
 
    cms.PSet(
    skip     = cms.bool(True),
    name     = cms.string('DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball'+VType),
    nickName = cms.string('DYJets50'),
    color    = cms.int32(19),
    xSec     = cms.double(3503.71)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball'+VType),
    nickName = cms.string('WJets'),
    color    = cms.int32(29),
    xSec     = cms.double(37509.0),
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola'+VType),
    nickName = cms.string('TtW'),
    color    = cms.int32(6),
    xSec     = cms.double(11.1)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('T_t-channel_TuneZ2star_8TeV-powheg-tauola'+VType),
    nickName = cms.string('Tt'),
    color    = cms.int32(6),
    xSec     = cms.double(56.4)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('T_s-channel_TuneZ2star_8TeV-powheg-tauola'+VType),
    nickName = cms.string('Ts'),
    color    = cms.int32(6),
    xSec     = cms.double(3.79)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola'+VType),
    nickName = cms.string('TbartW'),
    color    = cms.int32(6),
    xSec     = cms.double(11.1)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola'+VType),
    nickName = cms.string('Tbart'),
    color    = cms.int32(6),
    xSec     = cms.double(30.7)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola'+VType),
    nickName = cms.string('Tbars'),
    color    = cms.int32(6),
    xSec     = cms.double(1.76)
    ),


    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('WW_TuneZ2star_8TeV_pythia6_tauola'+VType),
    nickName = cms.string('WW'),
    color    = cms.int32(4),
    xSec     = cms.double(56.75)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('WZ_TuneZ2star_8TeV_pythia6_tauola'+VType),
    nickName = cms.string('WZ'),
    color    = cms.int32(4),
    xSec     = cms.double(33.85)
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('ZZ_TuneZ2star_8TeV_pythia6_tauola'+VType),
    nickName = cms.string('ZZ'),
    color    = cms.int32(4),
    xSec     = cms.double(8.297)
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
    name     = cms.string('TTWJets_8TeV-madgraph'+VType),
    nickName = cms.string('TTW'),
    color    = cms.int32(18),
    xSec     = cms.double(0.232),
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTZJets_8TeV-madgraph'+VType),
    nickName = cms.string('TTZ'),
    color    = cms.int32(18),
    xSec     = cms.double(0.2057),
    ),


    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTJets_SemiLeptMGDecays_8TeV-madgraph-part'+VType),
    nickName = cms.string('TTJetsSemiLept'),
    color    = cms.int32(41),
    xSec     = cms.double(xsecTT_SL),
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTJets_FullLeptMGDecays_8TeV-madgraph-part'+VType),
    nickName = cms.string('TTJetsFullLept'),
    color    = cms.int32(41),
    xSec     = cms.double(xsecTT_FL),
    ),
    
    
    ),


    #SLNoBLep 10000
    #SLNoBHad 10000
    #SLNoHiggs 8000
    #SL2wj     2000
    #SL1wj     4000
    #DL       10000
    
    hypo          = cms.untracked.int32(0),
    SoB           = cms.untracked.int32(1),
    maxChi2       = cms.untracked.double(2.5), # 2.5
    norm          = cms.untracked.int32(0),

     functions     = cms.vstring(
    '8.95351e+18*TMath::Landau(x, 5.67600e+01,1.01258e+01)',                # incl
    '2.95547e+17*TMath::Landau(x,7.61581e+01 ,1.89245e+01)',                # SL2wj
    '2.98474e+17*TMath::Landau(x,7.40196e+01 ,1.80142e+01)',                # SL1wj
    '6.28300e+16*TMath::Landau(x,8.03060e+01 ,1.81679e+01)',                # SLNoBHad
    'x>150?2.44515e+27*x^(-5.35628e+00):1.24208e+18*exp((-3.63162e-02)*x)', # SLNoHiggs
    'x>=12 ? x^(-2.010e-01)*exp((-1.5785e-02)*x) : 4.184e-02*x'),           # tth Pt
 
    switchoffOL   = cms.untracked.int32(0), ###### CHECK HERE
    speedup       = cms.untracked.int32(0), ###### CHECK HERE

    doTypeBTag6   = cms.untracked.int32(0),  #SL 6 jets
    doTypeBTag5   = cms.untracked.int32(0),  #SL 5 jets
    doTypeBTag4   = cms.untracked.int32(0),  #DL 4 jets
    
    
    doType0       = cms.untracked.int32(0),  #SL(4,2)  w/  W-tag
    doType1       = cms.untracked.int32(0),  #SL(4,2)  w/o W-tag
    doType2       = cms.untracked.int32(0),  #SL(4,1)
    doType3       = cms.untracked.int32(0),  #SL(4,3) 
    doType4       = cms.untracked.int32(0),  #SL(3,2)
    doType6       = cms.untracked.int32(0),  #DL(4,X)
    doType7       = cms.untracked.int32(0),  #DL(3M+1L,X)

    doType0ByBTagShape = cms.untracked.int32(0),
    doType1ByBTagShape = cms.untracked.int32(0),
    doType2ByBTagShape = cms.untracked.int32(0),
    doType3ByBTagShape = cms.untracked.int32(0),
    doType6ByBTagShape = cms.untracked.int32(1),

    useME         = cms.int32(1),
    useJac        = cms.int32(1),
    useMET        = cms.int32(1),
    useTF         = cms.int32(1),
    usePDF        = cms.int32(1),


    doubleGaussianB  = cms.untracked.int32(1),
    useBtag          = cms.untracked.int32(1),
    selectByBTagShape= cms.untracked.int32(1),
    
    printout     = cms.int32(0),
    debug        = cms.int32(0),   
    verbose      = cms.bool(False),

    MH           = cms.untracked.double(125.00),
    MT           = cms.untracked.double(174.30),
    MW           = cms.untracked.double( 80.19),

    MwL          = cms.untracked.double(60),
    MwH          = cms.untracked.double(100),
    MhL          = cms.untracked.double(110),
    MhH          = cms.untracked.double(140),

    btag_prob_cut_6jets = cms.untracked.double(0.988),
    btag_prob_cut_5jets = cms.untracked.double(0.992),
    btag_prob_cut_4jets = cms.untracked.double(0.992),
    
    massesH      = cms.vdouble(125.),
    #massesH      = cms.vdouble(105.,110.,115.,120.,125.,130.,135.,140.),
    massesT      = cms.vdouble(174.3),
    #massesT      = cms.vdouble(145, 155, 165, 174, 185, 195, 205),

    fixNumEvJob    = cms.untracked.int32(0),
    evLimits       = cms.vint32(0,-1),

    doJERbias  = cms.untracked.int32(0),   
    doCSVup    = cms.untracked.int32(0),
    doCSVdown  = cms.untracked.int32(0),
    doJECup    = cms.untracked.int32(0),
    doJECdown  = cms.untracked.int32(0),
    doJERup    = cms.untracked.int32(0),
    doJERdown  = cms.untracked.int32(0),

    )
