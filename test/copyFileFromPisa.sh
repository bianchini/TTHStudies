#! /bin/sh


pathI="  srm://stormfe1.pi.infn.it:8444/srm/managerv2?SFN=/cms/store/user/arizzi/Ntuple_Step1V42_Step2Tag_EDMV42_Step2_V4a_MC/"
pathO="  srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt/"
fileN="DiJetPt_"

command=" lcg-cp -v -b -D srmv2 "
#command=" srmcp -2 -globus_tcp_port_range 20000,25000  --streams_num=1  "

declare -a arr=(  

##Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola.root
##Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola.root
##Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola.root
##T_t-channel_TuneZ2star_8TeV-powheg-tauola.root
##T_s-channel_TuneZ2star_8TeV-powheg-tauola.root
##T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola.root
##TTJets_FullLeptMGDecays_8TeV-madgraph-part.root
#TTJets_SemiLeptMGDecays_8TeV-madgraph-part.root
##TTJets_HadronicMGDecays_8TeV-madgraph-part.root
#WH_WToLNu_HToBB_M-110_8TeV-powheg-herwigpp.root
#WH_WToLNu_HToBB_M-115_8TeV-powheg-herwigpp.root
#WH_WToLNu_HToBB_M-120_8TeV-powheg-herwigpp.root
#WH_WToLNu_HToBB_M-125_8TeV-powheg-herwigpp.root
#WH_WToLNu_HToBB_M-130_8TeV-powheg-herwigpp.root
#WH_WToLNu_HToBB_M-135_8TeV-powheg-herwigpp.root
DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball.root
WJetsToLNu_PtW-180_TuneZ2star_8TeV-madgraph-tarball-part.root
WJetsToLNu_PtW-100_TuneZ2star_8TeV-madgraph.root
WJetsToLNu_PtW-70To100_TuneZ2star_8TeV-madgraph.root
#WW_TuneZ2star_8TeV_pythia6_tauola.root
#WZ_TuneZ2star_8TeV_pythia6_tauola.root
ZZ_TuneZ2star_8TeV_pythia6_tauola.root
#WZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola.root
#ZZJetsTo2L2Q_TuneZ2star_8TeV-madgraph-tauola.root
#ZH_ZToLL_HToBB_M-110_8TeV-powheg-herwigpp3.root
#ZH_ZToLL_HToBB_M-115_8TeV-powheg-herwigpp.root
#ZH_ZToLL_HToBB_M-125_8TeV-powheg-herwigpp.root
#ZH_ZToLL_HToBB_M-130_8TeV-powheg-herwigpp3.root
#ZH_ZToLL_HToBB_M-135_8TeV-powheg-herwigpp.root
#ZH_ZToLL_HToBB_M-120_8TeV-powheg-herwigpp.root

#DoubleElectronRun2012CAug24RerecoEdmV42.root
#DoubleElectron_Run2012A-13Jul2012-v1_ProcFIXED.root
#DoubleElectron_Run2012A-13Jul2012-v1_ProcV2.root
#DoubleElectron_Run2012A-recover-06Aug2012-v1_ProcV2.root
#DoubleElectron_Run2012B-13Jul2012-v1_ProcV2.root
#DoubleElectron_Run2012C-PromptReco-v2_HBB_EDMNtupleV42_ProcV1.root
#DoubleElectron_Run2012C-PromptReco-v2_HBB_EDMNtupleV42_ProcV2.root
#SingleElectronRun2012AAug06EdmV42.root
#SingleElectronRun2012AJul13EdmV42b.root
#SingleElectronRun2012BJul13EdmV42.root
#SingleElectronRun2012CAug24RerecoEdmV42.root
#SingleElectronRun2012CPromptV2TopUpEdmV42.root
#SingleMuRun2012AAug06EdmV42.root
#SingleMuRun2012AJul13EdmV42.root
#SingleMuRun2012CAug24RerecoEdmV42.root
#SingleElectronRun2012CPromptv2EdmV42.root
#SingleMuRun2012CPromptv2EdmV42.root
#DoubleElectron_Run2012B-13Jul2012-v1_ProcFIXED.root
#SingleMuRun2012BJul13EdmV42.root
#SingleMuRun2012CPromptV2TopUpEdmV42.root

#DYJetsToLL_PtZ-70To100_TuneZ2star_8TeV-madgraph-tarball.root
#DY2JetsToLL_M-50_TuneZ2Star_8TeV-madgraph.root
#DYJetsToLL_HT-400ToInf_TuneZ2Star_8TeV-madgraph.root
#TBarToDilepton_tW-channel-DR_8TeV-powheg-tauola.root
#TBarToLeptons_s-channel_8TeV-powheg-tauola.root
#TToLeptons_s-channel_8TeV-powheg-tauola.root
#TToLeptons_t-channel_8TeV-powheg-tauola.root
#TToThadWlep_tW-channel-DR_8TeV-powheg-tauola.root
#TToTlepWhad_tW-channel-DR_8TeV-powheg-tauola.root
#ZH_ZToNuNu_HToBB_M-110_8TeV-powheg-herwigpp.root
#ZH_ZToNuNu_HToBB_M-115_8TeV-powheg-herwigpp.root
#ZH_ZToNuNu_HToBB_M-120_8TeV-powheg-herwigpp.root
#ZH_ZToNuNu_HToBB_M-125_8TeV-powheg-herwigpp.root
#ZH_ZToNuNu_HToBB_M-130_8TeV-powheg-herwigpp.root
#ZH_ZToNuNu_HToBB_M-135_8TeV-powheg-herwigpp.root
#ZJetsToNuNu_100_HT_200_TuneZ2Star_8TeV_madgraph.root
#ZJetsToNuNu_200_HT_400_TuneZ2Star_8TeV_madgraph.root
#ZJetsToNuNu_400_HT_inf_TuneZ2Star_8TeV_madgraph.root
#ZJetsToNuNu_50_HT_100_TuneZ2Star_8TeV_madgraph.root
#DY1JetsToLL_M-50_TuneZ2Star_8TeV-madgraph.root
#DYJetsToLL_PtZ-100_TuneZ2star_8TeV-madgraph.root
#TBarToLeptons_t-channel_8TeV-powheg-tauola.root
#TT_CT10_TuneZ2star_8TeV-powheg-tauola-7M.root
#TToDilepton_tW-channel-DR_8TeV-powheg-tauola.root
#DY4JetsToLL_M-50_TuneZ2Star_8TeV-madgraph.root
#TT_CT10_TuneZ2star_8TeV-powheg-tauola.root
#DYJetsToLL_HT-200To400_TuneZ2Star_8TeV-madgraph.root
#DYJetsToLL_PtZ-50To70_TuneZ2star_8TeV-madgraph-tarball.root
#ZJetsToNuNu_PtZ-100_8TeV-madgraph.root
#ZJetsToNuNu_Pt-100_8TeV-herwigpp.root
#ZJetsToLL_Pt-100_8TeV-herwigpp.root
#DY1JetsToLL_M-50_TuneZ2Star_8TeV-madgraph_procV2_mergeV1V2.root
#DY4JetsToLL_M-50_TuneZ2Star_8TeV-madgraph_mergeV1V2.root
#WJetsToLNu_PtW-100_8TeV-herwigpp.root
#QCD_HT-100To250_TuneZ2star_8TeV-madgraph-pythia.root
#QCD_HT-250To500_TuneZ2star_8TeV-madgraph-pythia6.root
#QCD_Pt-1000to1400_TuneZ2star_8TeV_pythia6.root
#QCD_Pt-120to170_TuneZ2star_8TeV_pythia6.root
#QCD_Pt-1400to1800_TuneZ2star_8TeV_pythia6.root
#QCD_Pt-170to300_TuneZ2star_8TeV_pythia6.root
#QCD_Pt-1800_TuneZ2star_8TeV_pythia6.root
#QCD_Pt-300to470_TuneZ2star_8TeV_pythia6.root
#QCD_Pt-470to600_TuneZ2star_8TeV_pythia6.root
#QCD_Pt-50to80_TuneZ2star_8TeV_pythia6.root
#QCD_Pt-600to800_TuneZ2star_8TeV_pythia6.root
#QCD_Pt-800to1000_TuneZ2star_8TeV_pythia6.root
#QCD_Pt-80to120_TuneZ2star_8TeV_pythia6.root
#DYJetsToLL_PtZ-180_TuneZ2star_8TeV-madgraph-tarball-part.root
#DYJetsToLL_PtZ-100_TuneZ2star_8TeV_ext-madgraph-tarball.root
#DYJetsToLL_PtZ-70To100_TuneZ2star_8TeV_ext-madgraph-tarball.root
#DYJetsToLL_PtZ-50To70_TuneZ2star_8TeV_ext-madgraph-tarball.root
#QCD_Pt-120to170_MuEnrichedPt5_TuneZ2star_8TeV_pythia6.root
#QCD_Pt-50to80_MuEnrichedPt5_TuneZ2star_8TeV_pythia6.root
#QCD_Pt-80to120_MuEnrichedPt5_TuneZ2star_8TeV_pythia6.root
#QCD_Pt-150_bEnriched_TuneZ2star_8TeV-pythia6-evtgen.root
#QCD_Pt-30To50_bEnriched_TuneZ2star_8TeV-pythia6-evtgen.root
#QCD_Pt-50To150_bEnriched_TuneZ2star_8TeV-pythia6-evtgen.root
#TT_8TeV-mcatnlo-part.root
#QCD_Pt-1000_MuEnrichedPt5_TuneZ2star_8TeV_pythia6.root
#QCD_Pt-170to300_MuEnrichedPt5_TuneZ2star_8TeV_pythia6.root
#QCD_Pt-300to470_MuEnrichedPt5_TuneZ2star_8TeV_pythia6.root
#QCD_Pt-470to600_MuEnrichedPt5_TuneZ2star_8TeV_pythia6.root
#QCD_Pt-600to800_MuEnrichedPt5_TuneZ2star_8TeV_pythia6.root
#QCD_Pt-800to1000_MuEnrichedPt5_TuneZ2star_8TeV_pythia6.root
#QCD_Pt_170_250_BCtoE_TuneZ2star_8TeV_pythia6.root
#QCD_Pt_250_350_BCtoE_TuneZ2star_8TeV_pythia6.root
#QCD_Pt_350_BCtoE_TuneZ2star_8TeV_pythia6.root
#QCD_Pt_80_170_BCtoE_TuneZ2star_8TeV_pythia6.root
#TT_8TeV-mcatnlo-withGENWEIGHT.root
#ZJetsToNuNu_Pt-100_8TeV-herwigppSummer12_DR53Xv2.root

)


for i in ${arr[@]}
do
  echo $command $pathI$fileN$i $pathO$fileN$i
  echo `$command $pathI$fileN$i $pathO$fileN$i`
done