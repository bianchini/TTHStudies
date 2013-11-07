#! /bin/sh

# ANALYSIS:
#  0 if only nominal 
#  1 if only systematics
#  2 if both 
DOSYS=1

# secondary name of the trees
NAME=New_MHscan

ls -ltr

declare -a arr=(  
    TTH125
    TTJets
    EWK
    SingleT
    DiBoson
    TTV
)

for i in ${arr[@]}
do
  if [ $DOSYS == 0 ] || [ $DOSYS == 2 ] ; then
      hadd -f MEAnalysis${NAME}_DL_nominal_$1_$i.root MEAnalysis${NAME}_SL_VType0_nominal_$1_$i.root MEAnalysis${NAME}_SL_VType1_nominal_$1_$i.root
      hadd -f MEAnalysis${NAME}_SL_nominal_$1_$i.root MEAnalysis${NAME}_SL_VType2_nominal_$1_$i.root MEAnalysis${NAME}_SL_VType3_nominal_$1_$i.root
  fi
  if ([ "$i" == "TTH125" ] || [ "$i" == "TTJets" ]) && ( [ $DOSYS == 1 ] || [ $DOSYS == 2 ] ); then
      hadd -f MEAnalysis${NAME}_SL_csvUp_$1_$i.root   MEAnalysis${NAME}_SL_VType2_csvUp_$1_$i.root   MEAnalysis${NAME}_SL_VType3_csvUp_$1_$i.root
      hadd -f MEAnalysis${NAME}_SL_csvDown_$1_$i.root MEAnalysis${NAME}_SL_VType2_csvDown_$1_$i.root MEAnalysis${NAME}_SL_VType3_csvDown_$1_$i.root
      hadd -f MEAnalysis${NAME}_SL_JECUp_$1_$i.root   MEAnalysis${NAME}_SL_VType2_JECUp_$1_$i.root   MEAnalysis${NAME}_SL_VType3_JECUp_$1_$i.root
      hadd -f MEAnalysis${NAME}_SL_JECDown_$1_$i.root MEAnalysis${NAME}_SL_VType2_JECDown_$1_$i.root MEAnalysis${NAME}_SL_VType3_JECDown_$1_$i.root
      hadd -f MEAnalysis${NAME}_DL_csvUp_$1_$i.root   MEAnalysis${NAME}_SL_VType0_csvUp_$1_$i.root   MEAnalysis${NAME}_SL_VType1_csvUp_$1_$i.root
      hadd -f MEAnalysis${NAME}_DL_csvDown_$1_$i.root MEAnalysis${NAME}_SL_VType0_csvDown_$1_$i.root MEAnalysis${NAME}_SL_VType1_csvDown_$1_$i.root
      hadd -f MEAnalysis${NAME}_DL_JECUp_$1_$i.root   MEAnalysis${NAME}_SL_VType0_JECUp_$1_$i.root   MEAnalysis${NAME}_SL_VType1_JECUp_$1_$i.root
      hadd -f MEAnalysis${NAME}_DL_JECDown_$1_$i.root MEAnalysis${NAME}_SL_VType0_JECDown_$1_$i.root MEAnalysis${NAME}_SL_VType1_JECDown_$1_$i.root
  fi
done