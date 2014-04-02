#! /bin/sh


if [ -z $ROOTSYS ]; then
 echo "ROOTSYS is not defined: source ROOT, or hadd won't work!"
 exit
fi


# ANALYSIS:
#  0 if only nominal 
#  1 if only systematics
#  2 if both 
DOSYS=0

# name tag for nominal analysis
NOMINAL=all

# secondary name of the trees
NAME=New

# just for debug
ls -ltr

if [ $DOSYS == 0 ] || [ $DOSYS == 2 ] ; then
    ./hadd.sh  $1${NOMINAL}_$2_TTJetsSemiLept  $NAME
    ./hadd.sh  $1${NOMINAL}_$2_TTJetsFullHad   $NAME
    ./hadd.sh  $1${NOMINAL}_$2_TTJetsFullLept  $NAME
    ./hadd.sh  $1${NOMINAL}_$2_TTH125          $NAME
    ./hadd.sh  $1${NOMINAL}_$2_TTZ             $NAME
    ./hadd.sh  $1${NOMINAL}_$2_DYJets10to50    $NAME
    ./hadd.sh  $1${NOMINAL}_$2_DYJets50        $NAME
    ./hadd.sh  $1${NOMINAL}_$2_TTW             $NAME
    ./hadd.sh  $1${NOMINAL}_$2_Tt              $NAME
    ./hadd.sh  $1${NOMINAL}_$2_WJets           $NAME
    ./hadd.sh  $1${NOMINAL}_$2_Ts              $NAME
    ./hadd.sh  $1${NOMINAL}_$2_ZZ              $NAME
    ./hadd.sh  $1${NOMINAL}_$2_Tbars           $NAME
    ./hadd.sh  $1${NOMINAL}_$2_WW              $NAME
    ./hadd.sh  $1${NOMINAL}_$2_WZ              $NAME
    ./hadd.sh  $1${NOMINAL}_$2_Tbart           $NAME
    ./hadd.sh  $1${NOMINAL}_$2_TbartW          $NAME
    ./hadd.sh  $1${NOMINAL}_$2_TtW             $NAME
fi
if [ $DOSYS == 1 ] || [ $DOSYS == 2 ] ; then
    ./hadd.sh  $1_csvDown_$2_TTH125          $NAME
    ./hadd.sh  $1_csvUp_$2_TTH125            $NAME
    ./hadd.sh  $1_JECDown_$2_TTH125          $NAME
    ./hadd.sh  $1_JECUp_$2_TTH125            $NAME
    ./hadd.sh  $1_csvDown_$2_TTJetsSemiLept  $NAME
    ./hadd.sh  $1_csvUp_$2_TTJetsSemiLept    $NAME
    ./hadd.sh  $1_JECDown_$2_TTJetsSemiLept  $NAME
    ./hadd.sh  $1_JECUp_$2_TTJetsSemiLept    $NAME
    ./hadd.sh  $1_csvDown_$2_TTJetsFullLept  $NAME
    ./hadd.sh  $1_csvUp_$2_TTJetsFullLept    $NAME
    ./hadd.sh  $1_JECDown_$2_TTJetsFullLept  $NAME
    ./hadd.sh  $1_JECUp_$2_TTJetsFullLept    $NAME
fi

if [ $DOSYS == 0 ] || [ $DOSYS == 2 ] ; then

    if ls MEAnalysis${NAME}_$1${NOMINAL}_$2_DYJets10to50.root MEAnalysis${NAME}_$1${NOMINAL}_$2_DYJets50.root MEAnalysis${NAME}_$1${NOMINAL}_$2_WJets.root  &> /dev/null ; then
	hadd -f MEAnalysis${NAME}_$1${NOMINAL}_$2_EWK.root  MEAnalysis${NAME}_$1${NOMINAL}_$2_DYJets10to50.root MEAnalysis${NAME}_$1${NOMINAL}_$2_DYJets50.root MEAnalysis${NAME}_$1${NOMINAL}_$2_WJets.root
	if ls  MEAnalysis${NAME}_$1${NOMINAL}_$2_EWK.root  &> /dev/null ; then
	    rm   MEAnalysis${NAME}_$1${NOMINAL}_$2_DYJets10to50.root MEAnalysis${NAME}_$1${NOMINAL}_$2_DYJets50.root MEAnalysis${NAME}_$1${NOMINAL}_$2_WJets.root
	fi
    fi

    if ls MEAnalysis${NAME}_$1${NOMINAL}_$2_Tt.root MEAnalysis${NAME}_$1${NOMINAL}_$2_Ts.root MEAnalysis${NAME}_$1${NOMINAL}_$2_Tbars.root MEAnalysis${NAME}_$1${NOMINAL}_$2_Tbart.root MEAnalysis${NAME}_$1${NOMINAL}_$2_TbartW.root MEAnalysis${NAME}_$1${NOMINAL}_$2_TtW.root   &> /dev/null ; then
	hadd -f  MEAnalysis${NAME}_$1${NOMINAL}_$2_SingleT.root MEAnalysis${NAME}_$1${NOMINAL}_$2_Tt.root MEAnalysis${NAME}_$1${NOMINAL}_$2_Ts.root MEAnalysis${NAME}_$1${NOMINAL}_$2_Tbars.root MEAnalysis${NAME}_$1${NOMINAL}_$2_Tbart.root MEAnalysis${NAME}_$1${NOMINAL}_$2_TbartW.root MEAnalysis${NAME}_$1${NOMINAL}_$2_TtW.root
	if  ls MEAnalysis${NAME}_$1${NOMINAL}_$2_SingleT.root  &> /dev/null ; then
	    rm MEAnalysis${NAME}_$1${NOMINAL}_$2_Tt.root MEAnalysis${NAME}_$1${NOMINAL}_$2_Ts.root MEAnalysis${NAME}_$1${NOMINAL}_$2_Tbars.root MEAnalysis${NAME}_$1${NOMINAL}_$2_Tbart.root MEAnalysis${NAME}_$1${NOMINAL}_$2_TbartW.root MEAnalysis${NAME}_$1${NOMINAL}_$2_TtW.root
	fi
    fi

    if ls MEAnalysis${NAME}_$1${NOMINAL}_$2_ZZ.root MEAnalysis${NAME}_$1${NOMINAL}_$2_WW.root MEAnalysis${NAME}_$1${NOMINAL}_$2_WZ.root &> /dev/null ; then
	hadd -f  MEAnalysis${NAME}_$1${NOMINAL}_$2_DiBoson.root MEAnalysis${NAME}_$1${NOMINAL}_$2_ZZ.root MEAnalysis${NAME}_$1${NOMINAL}_$2_WW.root MEAnalysis${NAME}_$1${NOMINAL}_$2_WZ.root
	if ls  MEAnalysis${NAME}_$1${NOMINAL}_$2_DiBoson.root  &> /dev/null ; then
	    rm  MEAnalysis${NAME}_$1${NOMINAL}_$2_ZZ.root MEAnalysis${NAME}_$1${NOMINAL}_$2_WW.root MEAnalysis${NAME}_$1${NOMINAL}_$2_WZ.root
	fi
    fi

    if ls MEAnalysis${NAME}_$1${NOMINAL}_$2_TTZ.root MEAnalysis${NAME}_$1${NOMINAL}_$2_TTW.root &> /dev/null ; then
	hadd -f  MEAnalysis${NAME}_$1${NOMINAL}_$2_TTV.root MEAnalysis${NAME}_$1${NOMINAL}_$2_TTZ.root MEAnalysis${NAME}_$1${NOMINAL}_$2_TTW.root
	if ls  MEAnalysis${NAME}_$1${NOMINAL}_$2_TTV.root &> /dev/null ; then
	    rm MEAnalysis${NAME}_$1${NOMINAL}_$2_TTZ.root MEAnalysis${NAME}_$1${NOMINAL}_$2_TTW.root
	fi
    fi

    if ls MEAnalysis${NAME}_$1${NOMINAL}_$2_TTJetsFullLept.root MEAnalysis${NAME}_$1${NOMINAL}_$2_TTJetsFullHad.root MEAnalysis${NAME}_$1${NOMINAL}_$2_TTJetsSemiLept.root  &> /dev/null ; then
	hadd -f MEAnalysis${NAME}_$1${NOMINAL}_$2_TTJets.root  MEAnalysis${NAME}_$1${NOMINAL}_$2_TTJetsFullLept.root  MEAnalysis${NAME}_$1${NOMINAL}_$2_TTJetsFullHad.root MEAnalysis${NAME}_$1${NOMINAL}_$2_TTJetsSemiLept.root 
	if ls MEAnalysis${NAME}_$1${NOMINAL}_$2_TTJets.root  &> /dev/null ; then
	    rm MEAnalysis${NAME}_$1${NOMINAL}_$2_TTJetsFullLept.root MEAnalysis${NAME}_$1${NOMINAL}_$2_TTJetsFullHad.root MEAnalysis${NAME}_$1${NOMINAL}_$2_TTJetsSemiLept.root
	fi
    fi
fi

if [ $DOSYS == 1 ] || [ $DOSYS == 2 ] ; then
    if ls  MEAnalysis${NAME}_$1_csvDown_$2_TTJetsFullLept.root MEAnalysis${NAME}_$1_csvDown_$2_TTJetsSemiLept.root  &> /dev/null ; then 
	hadd -f MEAnalysis${NAME}_$1_csvDown_$2_TTJets.root    MEAnalysis${NAME}_$1_csvDown_$2_TTJetsFullLept.root MEAnalysis${NAME}_$1_csvDown_$2_TTJetsSemiLept.root 
	if ls MEAnalysis${NAME}_$1_csvDown_$2_TTJets.root  &> /dev/null ; then
	    rm MEAnalysis${NAME}_$1_csvDown_$2_TTJetsFullLept.root MEAnalysis${NAME}_$1_csvDown_$2_TTJetsSemiLept.root
	fi
    fi
    if ls MEAnalysis${NAME}_$1_csvUp_$2_TTJetsFullLept.root   MEAnalysis${NAME}_$1_csvUp_$2_TTJetsSemiLept.root  &> /dev/null ; then 
	hadd -f MEAnalysis${NAME}_$1_csvUp_$2_TTJets.root      MEAnalysis${NAME}_$1_csvUp_$2_TTJetsFullLept.root   MEAnalysis${NAME}_$1_csvUp_$2_TTJetsSemiLept.root 
	if ls MEAnalysis${NAME}_$1_csvUp_$2_TTJets.root  &> /dev/null ; then
	    rm MEAnalysis${NAME}_$1_csvUp_$2_TTJetsFullLept.root   MEAnalysis${NAME}_$1_csvUp_$2_TTJetsSemiLept.root 
	fi
    fi
    if ls MEAnalysis${NAME}_$1_JECDown_$2_TTJetsFullLept.root MEAnalysis${NAME}_$1_JECDown_$2_TTJetsSemiLept.root   &> /dev/null ; then 
	hadd -f MEAnalysis${NAME}_$1_JECDown_$2_TTJets.root    MEAnalysis${NAME}_$1_JECDown_$2_TTJetsFullLept.root MEAnalysis${NAME}_$1_JECDown_$2_TTJetsSemiLept.root 
	if ls MEAnalysis${NAME}_$1_JECDown_$2_TTJets.root   &> /dev/null ; then
	    rm MEAnalysis${NAME}_$1_JECDown_$2_TTJetsFullLept.root MEAnalysis${NAME}_$1_JECDown_$2_TTJetsSemiLept.root    
	fi
    fi
    if ls MEAnalysis${NAME}_$1_JECUp_$2_TTJetsFullLept.root   MEAnalysis${NAME}_$1_JECUp_$2_TTJetsSemiLept.root  &> /dev/null ; then 
	hadd -f MEAnalysis${NAME}_$1_JECUp_$2_TTJets.root      MEAnalysis${NAME}_$1_JECUp_$2_TTJetsFullLept.root   MEAnalysis${NAME}_$1_JECUp_$2_TTJetsSemiLept.root 
	if ls MEAnalysis${NAME}_$1_JECUp_$2_TTJets.root    &> /dev/null ; then
	    rm MEAnalysis${NAME}_$1_JECUp_$2_TTJetsFullLept.root   MEAnalysis${NAME}_$1_JECUp_$2_TTJetsSemiLept.root 
	fi
    fi
fi



if ls ./MEAnalysis${NAME}_$1${NOMINAL}_$2_*SingleMu*_p* &> /dev/null ; then
    ls MEAnalysis${NAME}_$1${NOMINAL}_$2_*SingleMu*_p*.root
    hadd -f MEAnalysis${NAME}_$1${NOMINAL}_$2_Run2012_SingleMu.root MEAnalysis${NAME}_$1${NOMINAL}_$2_*SingleMu*_p*.root 
    #echo "Hadd"
    if ls  MEAnalysis${NAME}_$1${NOMINAL}_$2_Run2012_SingleMu.root  &> /dev/null ; then
	#echo "Remove"
	rm MEAnalysis${NAME}_$1${NOMINAL}_$2_*SingleMu*_p*.root
    fi
else
    echo "MEAnalysis${NAME}_*SingleMu* not found"
fi

if ls ./MEAnalysis${NAME}_$1${NOMINAL}_$2_*SingleElectron*_p* &> /dev/null ; then
    ls MEAnalysis${NAME}_$1${NOMINAL}_$2_*SingleElectron*_p*.root
    hadd -f MEAnalysis${NAME}_$1${NOMINAL}_$2_Run2012_SingleElectron.root MEAnalysis${NAME}_$1${NOMINAL}_$2_*SingleElectron*_p*.root 
    #echo "Hadd"
    if ls  MEAnalysis${NAME}_$1${NOMINAL}_$2_Run2012_SingleElectron.root  &> /dev/null ; then
	#echo "Remove"
	rm MEAnalysis${NAME}_$1${NOMINAL}_$2_*SingleElectron*_p*.root
    fi
else
    echo "MEAnalysis${NAME}_*SingleElectron* not found"
fi

if ls ./MEAnalysis${NAME}_$1${NOMINAL}_$2_*DoubleElectron*_p* &> /dev/null ; then
    ls MEAnalysis${NAME}_$1${NOMINAL}_$2_*DoubleElectron*_p*.root
    hadd -f MEAnalysis${NAME}_$1${NOMINAL}_$2_Run2012_DoubleElectron.root MEAnalysis${NAME}_$1${NOMINAL}_$2_*DoubleElectron*_p*.root 
    #echo "Hadd"
    if ls  MEAnalysis${NAME}_$1${NOMINAL}_$2_Run2012_DoubleElectron.root  &> /dev/null ; then
	#echo "Remove"
	rm MEAnalysis${NAME}_$1${NOMINAL}_$2_*DoubleElectron*_p*.root
    fi
else
    echo "MEAnalysis${NAME}_*DoubleElectron* not found"
fi