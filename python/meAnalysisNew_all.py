import FWCore.ParameterSet.Types  as CfgTypes
import FWCore.ParameterSet.Config as cms

VType     = ""

xsecTT_SL = 103.0
xsecTT_FL = 24.8

process = cms.Process("MEAnalysisNew")

process.fwliteInput = cms.PSet(

    outFileName   = cms.string("./root/MEAnalysisNewTEST.root"),
    pathToTF      = cms.string("./root/transferFunctionsTEST.root"),
    pathToCP      = cms.string("./root/ControlPlotsTEST.root"),
    pathToCP_smear= cms.string("./root/ControlPlotsTEST_std_gen.root"),

    pathToFile    = cms.string("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store//user/bianchi/HBB_EDMNtuple/AllHDiJetPt_V2"+VType+"/"),
    #pathToFile    = cms.string("./"),
    ordering      = cms.string("DiJetPt_"),
    lumi          = cms.untracked.double(12.1),

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
    name     = cms.string('TTH_HToBB_M-125_8TeV-pythia6_v2'+VType),
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
    name     = cms.string('TTZJets_8TeV-madgraph_v2'+VType),
    nickName = cms.string('TTZ'),
    color    = cms.int32(18),
    xSec     = cms.double(0.2057),
    ),


    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTJets_SemiLeptMGDecays_8TeV-madgraph'+VType),
    nickName = cms.string('TTJetsSemiLept'),
    color    = cms.int32(41),
    xSec     = cms.double(xsecTT_SL),
    ),

    cms.PSet(
    skip     = cms.bool(True),  
    name     = cms.string('TTJets_FullLeptMGDecays_8TeV-madgraph'+VType),
    nickName = cms.string('TTJetsFullLept'),
    color    = cms.int32(41),
    xSec     = cms.double(xsecTT_FL),
    ),

    cms.PSet(
    skip     = cms.bool(True), 
    name     = cms.string('DoubleElectronRun2012C-EcalRecover_11Dec2012-v1_v2'+VType),
    nickName = cms.string('Run2012_DoubleElectronRun2012C-EcalRecover_11Dec2012-v1_v2'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    ),
    cms.PSet(
    skip     = cms.bool(True), 
    name     = cms.string('DoubleElectronRun2012CAug24RerecoEdmV42'+VType),
    nickName = cms.string('Run2012_DoubleElectronRun2012CAug24RerecoEdmV42'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    ),
    cms.PSet(
    skip     = cms.bool(True), 
    name     = cms.string('DoubleElectronRun2012D'+VType),
    nickName = cms.string('Run2012_DoubleElectronRun2012D'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    ),
    cms.PSet(
    skip     = cms.bool(True), 
    name     = cms.string('DoubleElectron_Run2012A-13Jul2012-v1_ProcFIXED'+VType),
    nickName = cms.string('Run2012_DoubleElectron_Run2012A-13Jul2012-v1_ProcFIXED'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    ),
    cms.PSet(
    skip     = cms.bool(True), 
    name     = cms.string('DoubleElectron_Run2012A-recover-06Aug2012-v1_ProcV2'+VType),
    nickName = cms.string('Run2012_DoubleElectron_Run2012A-recover-06Aug2012-v1_ProcV2'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    ),
    cms.PSet(
    skip     = cms.bool(True), 
    name     = cms.string('DoubleElectron_Run2012B-13Jul2012-v1_ProcFIXED'+VType),
    nickName = cms.string('Run2012_DoubleElectron_Run2012B-13Jul2012-v1_ProcFIXED'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    ),
    cms.PSet(
    skip     = cms.bool(True), 
    name     = cms.string('DoubleElectron_Run2012C-PromptReco-v2_HBB_EDMNtupleV42_ProcV1'+VType),
    nickName = cms.string('Run2012_DoubleElectron_Run2012C-PromptReco-v2_HBB_EDMNtupleV42_ProcV1'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    ),
    cms.PSet(
    skip     = cms.bool(True), 
    name     = cms.string('DoubleElectron_Run2012C-PromptReco-v2_HBB_EDMNtupleV42_ProcV2'+VType),
    nickName = cms.string('Run2012_DoubleElectron_Run2012C-PromptReco-v2_HBB_EDMNtupleV42_ProcV2'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    ),



    cms.PSet(
    skip     = cms.bool(True),  # 148139
    name     = cms.string('SingleElectronRun2012AAug06EdmV42'+VType),
    nickName = cms.string('Run2012_SingleElectronRun2012AAug06EdmV42'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    ),
    cms.PSet(
    skip     = cms.bool(True), # 1551019
    name     = cms.string('SingleElectronRun2012AJul13EdmV42b'+VType),
    nickName = cms.string('Run2012_SingleElectronRun2012AJul13EdmV42b'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    ),
    cms.PSet(
    skip     = cms.bool(True), # 9351330
    name     = cms.string('SingleElectronRun2012BJul13EdmV42'+VType),
    nickName = cms.string('Run2012_SingleElectronRun2012BJul13EdmV42'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    ),
    cms.PSet(
    skip     = cms.bool(True), # 263593
    name     = cms.string('SingleElectronRun2012C-EcalRecover_11Dec2012-v1_v2'+VType),
    nickName = cms.string('Run2012_SingleElectronRun2012C-EcalRecover_11Dec2012-v1_v2'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    ),
    cms.PSet(
    skip     = cms.bool(True), # 1064158
    name     = cms.string('SingleElectronRun2012CAug24RerecoEdmV42'+VType),
    nickName = cms.string('Run2012_SingleElectronRun2012CAug24RerecoEdmV42'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    ),
    cms.PSet(
    skip     = cms.bool(True), # 9768094
    name     = cms.string('SingleElectronRun2012CPromptv2EdmV42'+VType),
    nickName = cms.string('Run2012_SingleElectronRun2012CPromptv2EdmV42'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    ),
    cms.PSet(
    skip     = cms.bool(True), # 3491407
    name     = cms.string('SingleElectronRun2012CPromptV2TopUpEdmV42'+VType),
    nickName = cms.string('Run2012_SingleElectronRun2012CPromptV2TopUpEdmV42'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    ),
    cms.PSet(
    skip     = cms.bool(True), # 16178887
    name     = cms.string('SingleElectronRun2012D-PromptReco-v1_v3'+VType),
    nickName = cms.string('Run2012_SingleElectronRun2012D-PromptReco-v1_v3'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    ),



    cms.PSet(
    skip     = cms.bool(True), # 90889
    name     = cms.string('SingleMuRun2012AAug06EdmV42'+VType),
    nickName = cms.string('Run2012_SingleMuRun2012AAug06EdmV42'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    ),
    cms.PSet(
    skip     = cms.bool(True), # 916855
    name     = cms.string('SingleMuRun2012AJul13EdmV42'+VType),
    nickName = cms.string('Run2012_SingleMuRun2012AJul13EdmV42'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    ),
    cms.PSet(
    skip     = cms.bool(True), # 6121904
    name     = cms.string('SingleMuRun2012BJul13EdmV42'+VType),
    nickName = cms.string('Run2012_SingleMuRun2012BJul13EdmV42'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    ),
    cms.PSet(
    skip     = cms.bool(True), 
    name     = cms.string('SingleMuRun2012C-EcalRecover_11Dec2012-v1_v2'+VType),
    nickName = cms.string('Run2012_SingleMuRun2012C-EcalRecover_11Dec2012-v1_v2'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    ),
    cms.PSet(
    skip     = cms.bool(True), 
    name     = cms.string('SingleMuRun2012CAug24RerecoEdmV42'+VType),
    nickName = cms.string('Run2012_SingleMuRun2012CAug24RerecoEdmV42'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    ),
    cms.PSet(
    skip     = cms.bool(True), 
    name     = cms.string('SingleMuRun2012CPromptV2TopUpEdmV42'+VType),
    nickName = cms.string('Run2012_SingleMuRun2012CPromptV2TopUpEdmV42'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    ),
    cms.PSet(
    skip     = cms.bool(True), 
    name     = cms.string('SingleMuRun2012CPromptv2EdmV42'+VType),
    nickName = cms.string('Run2012_SingleMuRun2012CPromptv2EdmV42'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    ),
    cms.PSet(
    skip     = cms.bool(True), # 11860310
    name     = cms.string('SingleMuRun2012D-PromptReco-v1'+VType),
    nickName = cms.string('Run2012_SingleMuRun2012D-PromptReco-v1'),
    color    = cms.int32(1),
    xSec     = cms.double(-1),
    ),
   
    
    ),


    #SLNoBLep 10000
    #SLNoBHad 10000
    #SLNoHiggs 8000
    #SL2wj     2000
    #SL1wj     4000
    #DL       10000

    # run both S and B hypotheses
    SoB             = cms.untracked.int32(1),

    # in case SoB=0, run only this hypothesis
    hypo            = cms.untracked.int32(0),

    # optimization0: re-run integral if bad chi2
    integralOption0 = cms.untracked.int32(1),
    
    # max chi2 for optimization0
    maxChi2         = cms.untracked.double(2.5),

    # optimization1: skip combinations loosely compatible with H/t/W mass
    integralOption1 = cms.untracked.int32(0),

    # optimization2: re-run the integral using last VEGAS grid only
    integralOption2 = cms.untracked.int32(0),

    # number of iterations for option2
    integralOption2_stage     = cms.untracked.int32(1),
    
    # number of iterations for option2
    integralOption2_niter     = cms.untracked.int32(1),

    # increase in evaluations for option2
    integralOption2_nevalfact = cms.untracked.double(1.),

    norm            = cms.untracked.int32(0),

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
    
    
    doType0       = cms.untracked.int32(1),  #SL(4,2)  w/  W-tag
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
    doType6ByBTagShape = cms.untracked.int32(0),

    useME         = cms.untracked.int32(1),
    useJac        = cms.untracked.int32(1),
    useMET        = cms.untracked.int32(1),
    useTF         = cms.untracked.int32(1),
    usePDF        = cms.untracked.int32(1),


    doubleGaussianB  = cms.untracked.int32(1),
    useBtag          = cms.untracked.int32(1),
    selectByBTagShape= cms.untracked.int32(0),
    useRegression    = cms.untracked.int32(0),
    
    printout     = cms.untracked.int32(1),
    debug        = cms.untracked.int32(0),   
    verbose      = cms.bool(False),

    MH           = cms.untracked.double(125.00),
    MT           = cms.untracked.double(174.30),
    MW           = cms.untracked.double( 80.19),

    MwL          = cms.untracked.double(60),
    MwH          = cms.untracked.double(100),
    MhL          = cms.untracked.double(110),
    MhH          = cms.untracked.double(140),

    btag_prob_cut_6jets = cms.untracked.double( 0.96675 ), #0.96675
    btag_prob_cut_5jets = cms.untracked.double( 0.98225 ),
    btag_prob_cut_4jets = cms.untracked.double( 0.95295 ),
    
    massesH      = cms.vdouble(125.),
    #massesH      = cms.vdouble(105.,110.,115.,120.,125.,130.,135.,140.),
    massesT      = cms.vdouble(174.3),
    #massesT      = cms.vdouble(145, 155, 165, 174, 185, 195, 205),

    fixNumEvJob    = cms.untracked.int32(1),
    evLimits       = cms.vint32(0,1),

    doJERbias  = cms.untracked.int32(0),   
    doCSVup    = cms.untracked.int32(0),
    doCSVdown  = cms.untracked.int32(0),
    doJECup    = cms.untracked.int32(0),
    doJECdown  = cms.untracked.int32(0),
    doJERup    = cms.untracked.int32(0),
    doJERdown  = cms.untracked.int32(0),

    systematics= cms.vint32(0,1,2,3,4),


    doGenLevelAnalysis  = cms.untracked.int32(0),
    ntuplizeAll         = cms.untracked.int32(0),
    
    )

#print (process.fwliteInput.doType0)
