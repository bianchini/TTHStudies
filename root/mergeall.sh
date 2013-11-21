#! /bin/sh

# ANALYSIS:
#  0 if only nominal 
#  1 if only systematics
#  2 if both 
DOSYS=0

# secondary name of the trees
NAME=New

# just for debug
ls -ltr

if [ $DOSYS == 0 ] || [ $DOSYS == 2 ] ; then
    ./hadd.sh  SL_$1_nominal_$2_TTJetsSemiLept  $NAME
    ./hadd.sh  SL_$1_nominal_$2_TTJetsFullLept  $NAME
    ./hadd.sh  SL_$1_nominal_$2_TTH125          $NAME
    ./hadd.sh  SL_$1_nominal_$2_TTZ             $NAME
    ./hadd.sh  SL_$1_nominal_$2_DYJets10to50    $NAME
    ./hadd.sh  SL_$1_nominal_$2_DYJets50        $NAME
    ./hadd.sh  SL_$1_nominal_$2_TTW             $NAME
    ./hadd.sh  SL_$1_nominal_$2_Tt              $NAME
    ./hadd.sh  SL_$1_nominal_$2_WJets           $NAME
    ./hadd.sh  SL_$1_nominal_$2_Ts              $NAME
    ./hadd.sh  SL_$1_nominal_$2_ZZ              $NAME
    ./hadd.sh  SL_$1_nominal_$2_Tbars           $NAME
    ./hadd.sh  SL_$1_nominal_$2_WW              $NAME
    ./hadd.sh  SL_$1_nominal_$2_WZ              $NAME
    ./hadd.sh  SL_$1_nominal_$2_Tbart           $NAME
    ./hadd.sh  SL_$1_nominal_$2_TbartW          $NAME
    ./hadd.sh  SL_$1_nominal_$2_TtW             $NAME
fi
if [ $DOSYS == 1 ] || [ $DOSYS == 2 ] ; then
    ./hadd.sh  SL_$1_csvDown_$2_TTH125          $NAME
    ./hadd.sh  SL_$1_csvUp_$2_TTH125            $NAME
    ./hadd.sh  SL_$1_JECDown_$2_TTH125          $NAME
    ./hadd.sh  SL_$1_JECUp_$2_TTH125            $NAME
    ./hadd.sh  SL_$1_csvDown_$2_TTJetsSemiLept  $NAME
    ./hadd.sh  SL_$1_csvUp_$2_TTJetsSemiLept    $NAME
    ./hadd.sh  SL_$1_JECDown_$2_TTJetsSemiLept  $NAME
    ./hadd.sh  SL_$1_JECUp_$2_TTJetsSemiLept    $NAME
    ./hadd.sh  SL_$1_csvDown_$2_TTJetsFullLept  $NAME
    ./hadd.sh  SL_$1_csvUp_$2_TTJetsFullLept    $NAME
    ./hadd.sh  SL_$1_JECDown_$2_TTJetsFullLept  $NAME
    ./hadd.sh  SL_$1_JECUp_$2_TTJetsFullLept    $NAME
fi

if [ $DOSYS == 0 ] || [ $DOSYS == 2 ] ; then

    if ls MEAnalysis${NAME}_SL_$1_nominal_$2_DYJets10to50.root MEAnalysis${NAME}_SL_$1_nominal_$2_DYJets50.root MEAnalysis${NAME}_SL_$1_nominal_$2_WJets.root  &> /dev/null ; then
	hadd -f MEAnalysis${NAME}_SL_$1_nominal_$2_EWK.root  MEAnalysis${NAME}_SL_$1_nominal_$2_DYJets10to50.root MEAnalysis${NAME}_SL_$1_nominal_$2_DYJets50.root MEAnalysis${NAME}_SL_$1_nominal_$2_WJets.root
	rm   MEAnalysis${NAME}_SL_$1_nominal_$2_DYJets10to50.root MEAnalysis${NAME}_SL_$1_nominal_$2_DYJets50.root MEAnalysis${NAME}_SL_$1_nominal_$2_WJets.root
    fi

    if ls MEAnalysis${NAME}_SL_$1_nominal_$2_Tt.root MEAnalysis${NAME}_SL_$1_nominal_$2_Ts.root MEAnalysis${NAME}_SL_$1_nominal_$2_Tbars.root MEAnalysis${NAME}_SL_$1_nominal_$2_Tbart.root MEAnalysis${NAME}_SL_$1_nominal_$2_TbartW.root MEAnalysis${NAME}_SL_$1_nominal_$2_TtW.root   &> /dev/null ; then
	hadd -f  MEAnalysis${NAME}_SL_$1_nominal_$2_SingleT.root MEAnalysis${NAME}_SL_$1_nominal_$2_Tt.root MEAnalysis${NAME}_SL_$1_nominal_$2_Ts.root MEAnalysis${NAME}_SL_$1_nominal_$2_Tbars.root MEAnalysis${NAME}_SL_$1_nominal_$2_Tbart.root MEAnalysis${NAME}_SL_$1_nominal_$2_TbartW.root MEAnalysis${NAME}_SL_$1_nominal_$2_TtW.root
	rm MEAnalysis${NAME}_SL_$1_nominal_$2_Tt.root MEAnalysis${NAME}_SL_$1_nominal_$2_Ts.root MEAnalysis${NAME}_SL_$1_nominal_$2_Tbars.root MEAnalysis${NAME}_SL_$1_nominal_$2_Tbart.root MEAnalysis${NAME}_SL_$1_nominal_$2_TbartW.root MEAnalysis${NAME}_SL_$1_nominal_$2_TtW.root
    fi

    if ls MEAnalysis${NAME}_SL_$1_nominal_$2_ZZ.root MEAnalysis${NAME}_SL_$1_nominal_$2_WW.root MEAnalysis${NAME}_SL_$1_nominal_$2_WZ.root &> /dev/null ; then
	hadd -f  MEAnalysis${NAME}_SL_$1_nominal_$2_DiBoson.root MEAnalysis${NAME}_SL_$1_nominal_$2_ZZ.root MEAnalysis${NAME}_SL_$1_nominal_$2_WW.root MEAnalysis${NAME}_SL_$1_nominal_$2_WZ.root
	rm  MEAnalysis${NAME}_SL_$1_nominal_$2_ZZ.root MEAnalysis${NAME}_SL_$1_nominal_$2_WW.root MEAnalysis${NAME}_SL_$1_nominal_$2_WZ.root
    fi

    if ls MEAnalysis${NAME}_SL_$1_nominal_$2_TTZ.root MEAnalysis${NAME}_SL_$1_nominal_$2_TTW.root &> /dev/null ; then
	hadd -f  MEAnalysis${NAME}_SL_$1_nominal_$2_TTV.root MEAnalysis${NAME}_SL_$1_nominal_$2_TTZ.root MEAnalysis${NAME}_SL_$1_nominal_$2_TTW.root
	rm MEAnalysis${NAME}_SL_$1_nominal_$2_TTZ.root MEAnalysis${NAME}_SL_$1_nominal_$2_TTW.root
    fi

    if ls MEAnalysis${NAME}_SL_$1_nominal_$2_TTJetsFullLept.root MEAnalysis${NAME}_SL_$1_nominal_$2_TTJetsSemiLept.root  &> /dev/null ; then
	hadd -f MEAnalysis${NAME}_SL_$1_nominal_$2_TTJets.root    MEAnalysis${NAME}_SL_$1_nominal_$2_TTJetsFullLept.root MEAnalysis${NAME}_SL_$1_nominal_$2_TTJetsSemiLept.root 
	rm MEAnalysis${NAME}_SL_$1_nominal_$2_TTJetsFullLept.root MEAnalysis${NAME}_SL_$1_nominal_$2_TTJetsSemiLept.root
    fi
fi

if [ $DOSYS == 1 ] || [ $DOSYS == 2 ] ; then
    if ls  MEAnalysis${NAME}_SL_$1_csvDown_$2_TTJetsFullLept.root MEAnalysis${NAME}_SL_$1_csvDown_$2_TTJetsSemiLept.root  &> /dev/null ; then 
	hadd -f MEAnalysis${NAME}_SL_$1_csvDown_$2_TTJets.root    MEAnalysis${NAME}_SL_$1_csvDown_$2_TTJetsFullLept.root MEAnalysis${NAME}_SL_$1_csvDown_$2_TTJetsSemiLept.root 
	rm MEAnalysis${NAME}_SL_$1_csvDown_$2_TTJetsFullLept.root MEAnalysis${NAME}_SL_$1_csvDown_$2_TTJetsSemiLept.root
    fi
    if ls MEAnalysis${NAME}_SL_$1_csvUp_$2_TTJetsFullLept.root   MEAnalysis${NAME}_SL_$1_csvUp_$2_TTJetsSemiLept.root  &> /dev/null ; then 
	hadd -f MEAnalysis${NAME}_SL_$1_csvUp_$2_TTJets.root      MEAnalysis${NAME}_SL_$1_csvUp_$2_TTJetsFullLept.root   MEAnalysis${NAME}_SL_$1_csvUp_$2_TTJetsSemiLept.root 
	rm MEAnalysis${NAME}_SL_$1_csvUp_$2_TTJetsFullLept.root   MEAnalysis${NAME}_SL_$1_csvUp_$2_TTJetsSemiLept.root 
    fi
    if ls MEAnalysis${NAME}_SL_$1_JECDown_$2_TTJetsFullLept.root MEAnalysis${NAME}_SL_$1_JECDown_$2_TTJetsSemiLept.root   &> /dev/null ; then 
	hadd -f MEAnalysis${NAME}_SL_$1_JECDown_$2_TTJets.root    MEAnalysis${NAME}_SL_$1_JECDown_$2_TTJetsFullLept.root MEAnalysis${NAME}_SL_$1_JECDown_$2_TTJetsSemiLept.root 
	rm MEAnalysis${NAME}_SL_$1_JECDown_$2_TTJetsFullLept.root MEAnalysis${NAME}_SL_$1_JECDown_$2_TTJetsSemiLept.root    
    fi
    if ls MEAnalysis${NAME}_SL_$1_JECUp_$2_TTJetsFullLept.root   MEAnalysis${NAME}_SL_$1_JECUp_$2_TTJetsSemiLept.root  &> /dev/null ; then 
	hadd -f MEAnalysis${NAME}_SL_$1_JECUp_$2_TTJets.root      MEAnalysis${NAME}_SL_$1_JECUp_$2_TTJetsFullLept.root   MEAnalysis${NAME}_SL_$1_JECUp_$2_TTJetsSemiLept.root 
	rm MEAnalysis${NAME}_SL_$1_JECUp_$2_TTJetsFullLept.root   MEAnalysis${NAME}_SL_$1_JECUp_$2_TTJetsSemiLept.root 
    fi
fi