#! /bin/sh

# ANALYSIS:
#  0 if only nominal 
#  1 if only systematics
#  2 if both 
DOSYS=1

# secondary name of the trees
NAME=New_MHscan

# just for debug
ls -ltr

if [ $DOSYS == 0 ] || [ $DOSYS == 2 ] ; then
    if [ "$1" == "VType2" ] ||  [ "$1" == "VType3" ] ; then
	./hadd.sh  SL_$1_nominal_$2_TTJetsSemiLept
    fi
    ./hadd.sh  SL_$1_nominal_$2_TTJetsFullLept
    ./hadd.sh  SL_$1_nominal_$2_TTH125
    ./hadd.sh  SL_$1_nominal_$2_TTZ
    ./hadd.sh  SL_$1_nominal_$2_DYJets10to50
    ./hadd.sh  SL_$1_nominal_$2_DYJets50
    ./hadd.sh  SL_$1_nominal_$2_TTW
    ./hadd.sh  SL_$1_nominal_$2_Tt
    ./hadd.sh  SL_$1_nominal_$2_WJets
    ./hadd.sh  SL_$1_nominal_$2_Ts
    ./hadd.sh  SL_$1_nominal_$2_ZZ
    ./hadd.sh  SL_$1_nominal_$2_Tbars
    ./hadd.sh  SL_$1_nominal_$2_WW
    ./hadd.sh  SL_$1_nominal_$2_WZ
    ./hadd.sh  SL_$1_nominal_$2_Tbart
    ./hadd.sh  SL_$1_nominal_$2_TbartW
    ./hadd.sh  SL_$1_nominal_$2_TtW
fi
if [ $DOSYS == 1 ] || [ $DOSYS == 2 ] ; then
    ./hadd.sh  SL_$1_csvDown_$2_TTH125
    ./hadd.sh  SL_$1_csvUp_$2_TTH125
    ./hadd.sh  SL_$1_JECDown_$2_TTH125
    ./hadd.sh  SL_$1_JECUp_$2_TTH125
    if [ "$1" == "VType2" ] ||  [ "$1" == "VType3" ] ; then
	./hadd.sh  SL_$1_csvDown_$2_TTJetsSemiLept
	./hadd.sh  SL_$1_csvUp_$2_TTJetsSemiLept
	./hadd.sh  SL_$1_JECDown_$2_TTJetsSemiLept
	./hadd.sh  SL_$1_JECUp_$2_TTJetsSemiLept
    fi
    ./hadd.sh  SL_$1_csvDown_$2_TTJetsFullLept
    ./hadd.sh  SL_$1_csvUp_$2_TTJetsFullLept
    ./hadd.sh  SL_$1_JECDown_$2_TTJetsFullLept
    ./hadd.sh  SL_$1_JECUp_$2_TTJetsFullLept
fi

if [ $DOSYS == 0 ] || [ $DOSYS == 2 ] ; then

    hadd -f MEAnalysis${NAME}_SL_$1_nominal_$2_EWK.root  MEAnalysis${NAME}_SL_$1_nominal_$2_DYJets10to50.root MEAnalysis${NAME}_SL_$1_nominal_$2_DYJets50.root MEAnalysis${NAME}_SL_$1_nominal_$2_WJets.root
    rm   MEAnalysis${NAME}_SL_$1_nominal_$2_DYJets10to50.root MEAnalysis${NAME}_SL_$1_nominal_$2_DYJets50.root MEAnalysis${NAME}_SL_$1_nominal_$2_WJets.root

    hadd -f  MEAnalysis${NAME}_SL_$1_nominal_$2_SingleT.root MEAnalysis${NAME}_SL_$1_nominal_$2_Tt.root MEAnalysis${NAME}_SL_$1_nominal_$2_Ts.root MEAnalysis${NAME}_SL_$1_nominal_$2_Tbars.root MEAnalysis${NAME}_SL_$1_nominal_$2_Tbart.root MEAnalysis${NAME}_SL_$1_nominal_$2_TbartW.root MEAnalysis${NAME}_SL_$1_nominal_$2_TtW.root
    rm MEAnalysis${NAME}_SL_$1_nominal_$2_Tt.root MEAnalysis${NAME}_SL_$1_nominal_$2_Ts.root MEAnalysis${NAME}_SL_$1_nominal_$2_Tbars.root MEAnalysis${NAME}_SL_$1_nominal_$2_Tbart.root MEAnalysis${NAME}_SL_$1_nominal_$2_TbartW.root MEAnalysis${NAME}_SL_$1_nominal_$2_TtW.root
    
    hadd -f  MEAnalysis${NAME}_SL_$1_nominal_$2_DiBoson.root MEAnalysis${NAME}_SL_$1_nominal_$2_ZZ.root MEAnalysis${NAME}_SL_$1_nominal_$2_WW.root MEAnalysis${NAME}_SL_$1_nominal_$2_WZ.root
    rm  MEAnalysis${NAME}_SL_$1_nominal_$2_ZZ.root MEAnalysis${NAME}_SL_$1_nominal_$2_WW.root MEAnalysis${NAME}_SL_$1_nominal_$2_WZ.root
    
    hadd -f  MEAnalysis${NAME}_SL_$1_nominal_$2_TTV.root MEAnalysis${NAME}_SL_$1_nominal_$2_TTZ.root MEAnalysis${NAME}_SL_$1_nominal_$2_TTW.root
    rm MEAnalysis${NAME}_SL_$1_nominal_$2_TTZ.root MEAnalysis${NAME}_SL_$1_nominal_$2_TTW.root
    
    if [ "$1" == "VType2" ] ||  [ "$1" == "VType3" ] ; then
	hadd -f MEAnalysis${NAME}_SL_$1_nominal_$2_TTJets.root    MEAnalysis${NAME}_SL_$1_nominal_$2_TTJetsFullLept.root MEAnalysis${NAME}_SL_$1_nominal_$2_TTJetsSemiLept.root 
	rm MEAnalysis${NAME}_SL_$1_nominal_$2_TTJetsFullLept.root MEAnalysis${NAME}_SL_$1_nominal_$2_TTJetsSemiLept.root
    else
	hadd -f MEAnalysis${NAME}_SL_$1_nominal_$2_TTJets.root MEAnalysis${NAME}_SL_$1_nominal_$2_TTJetsFullLept.root
	rm MEAnalysis${NAME}_SL_$1_nominal_$2_TTJetsFullLept.root
    fi
fi

if [ $DOSYS == 1 ] || [ $DOSYS == 2 ] ; then
    if [ "$1" == "VType2" ] ||  [ "$1" == "VType3" ] ; then
	hadd -f MEAnalysis${NAME}_SL_$1_csvDown_$2_TTJets.root    MEAnalysis${NAME}_SL_$1_csvDown_$2_TTJetsFullLept.root MEAnalysis${NAME}_SL_$1_csvDown_$2_TTJetsSemiLept.root 
	rm MEAnalysis${NAME}_SL_$1_csvDown_$2_TTJetsFullLept.root MEAnalysis${NAME}_SL_$1_csvDown_$2_TTJetsSemiLept.root
	hadd -f MEAnalysis${NAME}_SL_$1_csvUp_$2_TTJets.root      MEAnalysis${NAME}_SL_$1_csvUp_$2_TTJetsFullLept.root   MEAnalysis${NAME}_SL_$1_csvUp_$2_TTJetsSemiLept.root 
	rm MEAnalysis${NAME}_SL_$1_csvUp_$2_TTJetsFullLept.root   MEAnalysis${NAME}_SL_$1_csvUp_$2_TTJetsSemiLept.root 
	hadd -f MEAnalysis${NAME}_SL_$1_JECDown_$2_TTJets.root    MEAnalysis${NAME}_SL_$1_JECDown_$2_TTJetsFullLept.root MEAnalysis${NAME}_SL_$1_JECDown_$2_TTJetsSemiLept.root 
	rm MEAnalysis${NAME}_SL_$1_JECDown_$2_TTJetsFullLept.root MEAnalysis${NAME}_SL_$1_JECDown_$2_TTJetsSemiLept.root    
	hadd -f MEAnalysis${NAME}_SL_$1_JECUp_$2_TTJets.root      MEAnalysis${NAME}_SL_$1_JECUp_$2_TTJetsFullLept.root   MEAnalysis${NAME}_SL_$1_JECUp_$2_TTJetsSemiLept.root 
	rm MEAnalysis${NAME}_SL_$1_JECUp_$2_TTJetsFullLept.root   MEAnalysis${NAME}_SL_$1_JECUp_$2_TTJetsSemiLept.root 
    else
	hadd -f MEAnalysis${NAME}_SL_$1_csvDown_$2_TTJets.root MEAnalysis${NAME}_SL_$1_csvDown_$2_TTJetsFullLept.root 
	rm MEAnalysis${NAME}_SL_$1_csvDown_$2_TTJetsFullLept.root 
	hadd -f MEAnalysis${NAME}_SL_$1_csvUp_$2_TTJets.root   MEAnalysis${NAME}_SL_$1_csvUp_$2_TTJetsFullLept.root   
	rm MEAnalysis${NAME}_SL_$1_csvUp_$2_TTJetsFullLept.root   
	hadd -f MEAnalysis${NAME}_SL_$1_JECDown_$2_TTJets.root MEAnalysis${NAME}_SL_$1_JECDown_$2_TTJetsFullLept.root 
	rm MEAnalysis${NAME}_SL_$1_JECDown_$2_TTJetsFullLept.root 
	hadd -f MEAnalysis${NAME}_SL_$1_JECUp_$2_TTJets.root   MEAnalysis${NAME}_SL_$1_JECUp_$2_TTJetsFullLept.root   
	rm MEAnalysis${NAME}_SL_$1_JECUp_$2_TTJetsFullLept.root   
    fi
fi