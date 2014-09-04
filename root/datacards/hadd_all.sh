
INDIR=$1 #Specify directory with input root trees

# files in INDIR should be named as basename_PROCESS_jobNr.root

if [ -z $ROOTSYS ]; then
 echo "ROOTSYS is not defined: source ROOT, or hadd won't work!"
 exit
fi

OUTDIR=$INDIR"_merged"
echo Saving merged root files to $OUTDIR

if [ ! -d "$OUTDIR" ]; then
    mkdir $OUTDIR
fi

PROCESSES=(
"TTH125"
"TTJetsBB"
"TTJetsJJ"
"TTJetsBJ"
"TTJetsCC"
"TTV"
"SingleT"
"EWK"
"DiBoson"
"Run2012_SingleElectron"
"Run2012_DoubleElectron"
"Run2012_SingleMu"
#'bEnriched'
#'BCtoE'

)

CUTS=(
#"SL_g5jg3t" #tight control region
#"SL_g5jg2t" #loose control retion

#"SL_g5jg2t_eta15"
#"SL_g5jg3t"

#----- btag LR--------
#"SL_5j"
#"SL_6j"
#"SL_5jg2t" 
#"SL_g6jg2t"
#----------------------
#"SL_g4jg2t"

#"SL_g4j"
#"SL_4j"

#"SL_g4jg2t"
#"SL_5jg2t"

#----- BDT SL table ----
#"SL_g6jg2t"
#"SL_4j3t"
#"SL_5j3t"
#"SL_g6j3t"
#"SL_4j4t"
#"SL_5jg4t"
#"SL_g6jg4t"
#---------------------
#"SL_g5jg3t"
#"SL_g6jg3t"
#---------- DL --------------
#"DL_g2jg2t"

#"DL_g4j_z"
#"DL_g4j"
#"DL_g2jg2t_z"
#"DL_g2jg2t"

#----- BDT DL table ----
#"DL_g4j2t"
#"DL_3j2t"
#"DL_g3jg3t"

#"DL_3j3t"
#"DL_g43t"
#"DL_g4g4t"
#---------------------
#"DL_g3jg3t_mm"
#"DL_g3jg3t_ee"
#"DL_g3jg3t_em"
#"DL_g4j"

#"DL_cat6_HP_mm"
#"DL_cat6_HP_ee"
#"DL_cat6_HP_em"
#"DL_cat6_LP_mm"
#"DL_cat6_LP_ee"
#"DL_cat6_LP_em"

#"SL_cat1_HP"
#"SL_cat2_HP"
#"SL_cat3_HP"
#"DL_cat4_HP"

#"SL_cat1_LP"
#"SL_cat2_LP"
#"SL_cat3_LP"
#"DL_cat4_LP"
#----- cat 2 tests------
"SL_cat2_HP"
#"SL_cat2_HP_muon"
#"SL_cat2_HP_electron"
#"SL_cat2_loose"

)

VARS=(
#"MET_pt"
#"MTln"
#"Mll"
#"MET_sumEt"

#"nPVs"

#"btag_LR"
#"jetsAboveCut"
#"numJets"
#"numBTagM"

#"bjet_pt"
#"bjet_eta"
#"leadjet_pt"
#"leadjet_eta"

#"muon_pt"
#"muon_eta"
#"muon_rIso"

#"Vtype"

#"electron_pt"
#"electron_eta"
#"electron_rIso"
#"electron_dxy"

#---- cat 2 tests ---
#"btag_LR"
"p_sb"
"p_bj"

)


for CUT in ${CUTS[@]}
do
  for VAR in ${VARS[@]}
    do
    for PROC in ${PROCESSES[@]} 
      do
      echo "Start processing: proc = "$PROC", var = "$VAR
      if ls ${INDIR}/*"${VAR}"_${CUT}*_${PROC}_*.root &> /dev/null ; then # check whether input files exist
	  INFILENAMES=(`ls ${INDIR}/*"${VAR}"_${CUT}*_${PROC}_*.root`)
      else
	  echo for proc = $PROC and var = $VAR -- input files not found
	  continue
      fi
      
      INFILENAME=$(basename `echo ${INFILENAMES[0]} | cut -d'.' --complement -f2-`) #drop .root extension and path
      INFILENAME=${INFILENAME%${PROC}*}"$PROC" # drop _xx after process name

      if [ ${PROC:0:3} != "QCD" ]; then
	  echo Merging to output file: $INFILENAME".root"
	  hadd -f ${OUTDIR}/${INFILENAME}".root" "${INDIR}/${INFILENAME}"*.root
      fi

    done
#    if ls ${INDIR}/*${VAR}*${CUT}*QCD*BCtoE*.root &> /dev/null ; then
#	OUTFILE=${OUTDIR}/${VAR}_${CUT}_"QCD_BCtoE.root"
#	echo "Merging QCD_BCtoE to "$OUTFILE
#	hadd -f $OUTFILE ${INDIR}/*${VAR}*${CUT}*QCD*BCtoE*.root
#    fi

#    if ls ${INDIR}/*${VAR}*${CUT}*QCD*bEnriched*.root &> /dev/null ; then
#	OUTFILE=${OUTDIR}/${VAR}_${CUT}_"QCD_bEnriched.root"
#	echo "Merging QCD_bEnriched to "$OUTFILE
#        hadd -f $OUTFILE ${INDIR}/*${VAR}*${CUT}*QCD*bEnriched*.root
#    fi

  done
done


