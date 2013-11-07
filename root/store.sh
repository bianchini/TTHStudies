#! /bin/sh

#srmcp -2  file:///MEAnalysis_SL_nominal_TTJetsSemiLept.root  srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/bianchi/Trees/MEM/Sept09_2013/MEAnalysis_SL_nominal_TTJetsSemiLept.root

NAME=""
VERSION=v7

declare -a arr=(  

#MEAnalysis${NAME}_SL_VType0_nominal_${VERSION}_DiBoson.root
#MEAnalysis${NAME}_SL_VType0_nominal_${VERSION}_EWK.root
#MEAnalysis${NAME}_SL_VType0_nominal_${VERSION}_SingleT.root
#MEAnalysis${NAME}_SL_VType0_nominal_${VERSION}_TTH125.root
#MEAnalysis${NAME}_SL_VType0_nominal_${VERSION}_TTJets.root
#MEAnalysis${NAME}_SL_VType0_nominal_${VERSION}_TTV.root
#MEAnalysis${NAME}_SL_VType0_csvDown_${VERSION}_TTH125.root
#MEAnalysis${NAME}_SL_VType0_csvDown_${VERSION}_TTJets.root
#MEAnalysis${NAME}_SL_VType0_csvUp_${VERSION}_TTH125.root
#MEAnalysis${NAME}_SL_VType0_csvUp_${VERSION}_TTJets.root
#MEAnalysis${NAME}_SL_VType0_JECDown_${VERSION}_TTH125.root
#MEAnalysis${NAME}_SL_VType0_JECDown_${VERSION}_TTJets.root
#MEAnalysis${NAME}_SL_VType0_JECUp_${VERSION}_TTH125.root
#MEAnalysis${NAME}_SL_VType0_JECUp_${VERSION}_TTJets.root

#MEAnalysis${NAME}_SL_VType1_nominal_${VERSION}_DiBoson.root
#MEAnalysis${NAME}_SL_VType1_nominal_${VERSION}_EWK.root
#MEAnalysis${NAME}_SL_VType1_nominal_${VERSION}_SingleT.root
#MEAnalysis${NAME}_SL_VType1_nominal_${VERSION}_TTH125.root
#MEAnalysis${NAME}_SL_VType1_nominal_${VERSION}_TTJets.root
#MEAnalysis${NAME}_SL_VType1_nominal_${VERSION}_TTV.root
#MEAnalysis${NAME}_SL_VType1_csvDown_${VERSION}_TTH125.root
#MEAnalysis${NAME}_SL_VType1_csvDown_${VERSION}_TTJets.root
#MEAnalysis${NAME}_SL_VType1_csvUp_${VERSION}_TTH125.root
#MEAnalysis${NAME}_SL_VType1_csvUp_${VERSION}_TTJets.root
#MEAnalysis${NAME}_SL_VType1_JECDown_${VERSION}_TTH125.root
#MEAnalysis${NAME}_SL_VType1_JECDown_${VERSION}_TTJets.root
#MEAnalysis${NAME}_SL_VType1_JECUp_${VERSION}_TTH125.root
#MEAnalysis${NAME}_SL_VType1_JECUp_${VERSION}_TTJets.root

#MEAnalysis${NAME}_DL_nominal_${VERSION}_DiBoson.root
#MEAnalysis${NAME}_DL_nominal_${VERSION}_EWK.root
#MEAnalysis${NAME}_DL_nominal_${VERSION}_SingleT.root
#MEAnalysis${NAME}_DL_nominal_${VERSION}_TTH125.root
MEAnalysis${NAME}_DL_nominal_${VERSION}_TTJets.root
#MEAnalysis${NAME}_DL_nominal_${VERSION}_TTV.root
#MEAnalysis${NAME}_DL_csvDown_${VERSION}_TTH125.root
MEAnalysis${NAME}_DL_csvDown_${VERSION}_TTJets.root
#MEAnalysis${NAME}_DL_csvUp_${VERSION}_TTH125.root
MEAnalysis${NAME}_DL_csvUp_${VERSION}_TTJets.root
#MEAnalysis${NAME}_DL_JECDown_${VERSION}_TTH125.root
MEAnalysis${NAME}_DL_JECDown_${VERSION}_TTJets.root
#MEAnalysis${NAME}_DL_JECUp_${VERSION}_TTH125.root
MEAnalysis${NAME}_DL_JECUp_${VERSION}_TTJets.root

#MEAnalysis${NAME}_SL_VType2_nominal_${VERSION}_DiBoson.root
#MEAnalysis${NAME}_SL_VType2_nominal_${VERSION}_EWK.root
#MEAnalysis${NAME}_SL_VType2_nominal_${VERSION}_SingleT.root
#MEAnalysis${NAME}_SL_VType2_nominal_${VERSION}_TTH125.root
#MEAnalysis${NAME}_SL_VType2_nominal_${VERSION}_TTJets.root
#MEAnalysis${NAME}_SL_VType2_nominal_${VERSION}_TTV.root
#MEAnalysis${NAME}_SL_VType2_csvDown_${VERSION}_TTH125.root
#MEAnalysis${NAME}_SL_VType2_csvDown_${VERSION}_TTJets.root
#MEAnalysis${NAME}_SL_VType2_csvUp_${VERSION}_TTH125.root
#MEAnalysis${NAME}_SL_VType2_csvUp_${VERSION}_TTJets.root
#MEAnalysis${NAME}_SL_VType2_JECDown_${VERSION}_TTH125.root
#MEAnalysis${NAME}_SL_VType2_JECDown_${VERSION}_TTJets.root
#MEAnalysis${NAME}_SL_VType2_JECUp_${VERSION}_TTH125.root
#MEAnalysis${NAME}_SL_VType2_JECUp_${VERSION}_TTJets.root

#MEAnalysis${NAME}_SL_VType3_nominal_${VERSION}_DiBoson.root
#MEAnalysis${NAME}_SL_VType3_nominal_${VERSION}_EWK.root
#MEAnalysis${NAME}_SL_VType3_nominal_${VERSION}_SingleT.root
#MEAnalysis${NAME}_SL_VType3_nominal_${VERSION}_TTH125.root
#MEAnalysis${NAME}_SL_VType3_nominal_${VERSION}_TTJets.root
#MEAnalysis${NAME}_SL_VType3_nominal_${VERSION}_TTV.root
#MEAnalysis${NAME}_SL_VType3_csvDown_${VERSION}_TTH125.root
#MEAnalysis${NAME}_SL_VType3_csvDown_${VERSION}_TTJets.root
#MEAnalysis${NAME}_SL_VType3_csvUp_${VERSION}_TTH125.root
#MEAnalysis${NAME}_SL_VType3_csvUp_${VERSION}_TTJets.root
#MEAnalysis${NAME}_SL_VType3_JECDown_${VERSION}_TTH125.root
#MEAnalysis${NAME}_SL_VType3_JECDown_${VERSION}_TTJets.root
#MEAnalysis${NAME}_SL_VType3_JECUp_${VERSION}_TTH125.root
#MEAnalysis${NAME}_SL_VType3_JECUp_${VERSION}_TTJets.root

#MEAnalysis${NAME}_SL_nominal_${VERSION}_DiBoson.root
#MEAnalysis${NAME}_SL_nominal_${VERSION}_EWK.root
#MEAnalysis${NAME}_SL_nominal_${VERSION}_SingleT.root
#MEAnalysis${NAME}_SL_nominal_${VERSION}_TTH125.root
#MEAnalysis${NAME}_SL_nominal_${VERSION}_TTJets.root
#MEAnalysis${NAME}_SL_nominal_${VERSION}_TTV.root
#MEAnalysis${NAME}_SL_csvDown_${VERSION}_TTH125.root
#MEAnalysis${NAME}_SL_csvDown_${VERSION}_TTJets.root
#MEAnalysis${NAME}_SL_csvUp_${VERSION}_TTH125.root
#MEAnalysis${NAME}_SL_csvUp_${VERSION}_TTJets.root
#MEAnalysis${NAME}_SL_JECDown_${VERSION}_TTH125.root
#MEAnalysis${NAME}_SL_JECDown_${VERSION}_TTJets.root
#MEAnalysis${NAME}_SL_JECUp_${VERSION}_TTH125.root
#MEAnalysis${NAME}_SL_JECUp_${VERSION}_TTJets.root

#transferFunctionsNew_partonE_new.root
#ControlPlotsNew_new.root
)

for i in ${arr[@]}
do
  #echo "Doing file " $i
  ls $i
  srmcp -2  file:///$i srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/bianchi/Trees/MEM/$1/$i
done