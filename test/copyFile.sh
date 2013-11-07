#! /bin/sh

VType="_VType2"

#pathI=" srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store//user/bianchi/HBB_EDMNtuple/ZllHDiJetPt/"
pathI=" srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store//user/bianchi/HBB_EDMNtuple/AllHDiJetPt"$VType"/"
#fileN="ZllH.DiJetPt.Oct22."
fileN="DiJetPt_"
#pathO=" file:////scratch/bianchi/HBB_EDMNtuple/Zll.H.DiJetPt/"
pathO=" file:////scratch/bianchi/HBB_EDMNtuple/All.H.DiJetPt/"

command1="srmcp -2 "
command2=" srmrm "

declare -a arr=(  
    
    
    #WW_TuneZ2star_8TeV_pythia6_tauola
    #WZ_TuneZ2star_8TeV_pythia6_tauola
    #ZZ_TuneZ2star_8TeV_pythia6_tauola
    #DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball
    #T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola
    #T_s-channel_TuneZ2star_8TeV-powheg-tauola
    #T_t-channel_TuneZ2star_8TeV-powheg-tauola
    #Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola
    #Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola
    #Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola
    #TTJets_HadronicMGDecays_8TeV-madgraph-part
    #TTJets_FullLeptMGDecays_8TeV-madgraph-part
    #TTJets_SemiLeptMGDecays_8TeV-madgraph-part
    #WJetsToLNu_PtW-180_TuneZ2star_8TeV-madgraph-tarball-part
    #WJetsToLNu_PtW-100_TuneZ2star_8TeV-madgraph
    #WJetsToLNu_PtW-70To100_TuneZ2star_8TeV-madgraph
    #ZH_ZToLL_HToBB_M-125_8TeV-powheg-herwigpp
    #WH_WToLNu_HToBB_M-125_8TeV-powheg-herwigpp
    #DataZee
    #DataZmm
    #DoubleElectronRun2012CAug24RerecoEdmV42
    #DoubleElectron_Run2012A-13Jul2012-v1_ProcFIXED
    #DoubleElectron_Run2012A-13Jul2012-v1_ProcV2
    #DoubleElectron_Run2012A-recover-06Aug2012-v1_ProcV2
    #DoubleElectron_Run2012B-13Jul2012-v1_ProcV2
    #DoubleElectron_Run2012C-PromptReco-v2_HBB_EDMNtupleV42_ProcV1
    #DoubleElectron_Run2012C-PromptReco-v2_HBB_EDMNtupleV42_ProcV2
    #SingleElectronRun2012AAug06EdmV42
    #SingleElectronRun2012AJul13EdmV42b
    #SingleElectronRun2012BJul13EdmV42
    #SingleElectronRun2012CAug24RerecoEdmV42
    #SingleElectronRun2012CPromptV2TopUpEdmV42
    #SingleElectronRun2012CPromptv2EdmV42
    SingleMuRun2012AAug06EdmV42
    SingleMuRun2012AJul13EdmV42
    SingleMuRun2012CAug24RerecoEdmV42
    SingleMuRun2012CPromptv2EdmV42
    SingleMuRun2012BJul13EdmV42
    SingleMuRun2012CPromptV2TopUpEdmV42
    #DoubleElectron_Run2012B-13Jul2012-v1_ProcFIXED
    #TTJets_Merged
    #ZH_ZToLL_HToBB_M-110_8TeV-powheg-herwigpp
    #DYJetsToLL_PtZ-100_TuneZ2star_8TeV-madgraph
    #TToLeptons_t-channel_8TeV-powheg-tauola 
    #DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph
    #DY3JetsToLL_M-50_TuneZ2Star_8TeV-madgraph
    #ZH_ZToLL_HToBB_M-120_8TeV-powheg-herwigpp
    #TBarToLeptons_s-channel_8TeV-powheg-tauola
    #TBarToLeptons_t-channel_8TeV-powheg-tauola
    #ZH_ZToLL_HToBB_M-130_8TeV-powheg-herwigpp
    #DY4JetsToLL_M-50_TuneZ2Star_8TeV-madgraph
    #DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph
    #ZJetsToLL_Pt-100_8TeV-herwigpp
    #ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola
    #DY1JetsToLL_M-50_TuneZ2Star_8TeV-madgraph
    #ZH_ZToLL_HToBB_M-115_8TeV-powheg-herwigpp
    #TToThadWlep_tW-channel-DR_8TeV-powheg-tauola
    #TBarToDilepton_tW-channel-DR_8TeV-powheg-tauola
    #ZH_ZToLL_HToBB_M-135_8TeV-powheg-herwigpp
    #DY2JetsToLL_M-50_TuneZ2Star_8TeV-madgraph
    #TToTlepWhad_tW-channel-DR_8TeV-powheg-tauola
    #DYJetsToLL_PtZ-50To70_TuneZ2star_8TeV-madgraph-tarball
    #TToLeptons_s-channel_8TeV-powheg-tauola
    #TToDilepton_tW-channel-DR_8TeV-powheg-tauola
    #DYJetsToLL_PtZ-70To100_TuneZ2star_8TeV-madgraph-tarball
    #WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola
)
 

for i in ${arr[@]}
do
   #echo $command2 $pathO$i
   #echo `$command2 $pathO$i`

   echo $command1 $pathI$fileN$i$VType".root" $pathO$fileN$i$VType".root"
   echo `$command1 $pathI$fileN$i$VType".root " $pathO$fileN$i`
   #mv $fileN$i /scratch/bianchi/HBB_EDMNtuple/All.H.DiJetPt
done

