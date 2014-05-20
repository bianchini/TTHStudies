
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
"TTV"
"SingleT"
"EWK"
"DiBoson"
"Run2012_SingleElectron"
"Run2012_DoubleElectron"
"Run2012_SingleMu"
)

CUTS=(
#"SL_g6jg3t" #tight control region
"SL_g5jg2t" #loose control retion
#"SL_g5jg2t_eta15"
#"SL_5j"
#"SL_g6j"
#"SL_g4jg2t"

#"SL_g6j2t"
#"SL_4j3t"
#"SL_5j3t"
#"SL_g6j3t"
#"SL_4j4t"
#"SL_5jg4t"
#"SL_g6jg4t"

#"SL_g5jg3t"
#"SL_g6jg3t"
#"DL_g2jg2t"
#"DL_g4j"

"DL_g2jg2t"
#"DL_3j2t"
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

)

VARS=(
"MET_pt"
#"btag_LR"
#"numJets"
#"numBTagM"
#"Mll"
#"MTln"
#"muon_pt"
#"muon_eta"
#"muon_rIso"

#"electron_pt"
#"electron_eta"
#"electron_rIso"
)

for PROC in ${PROCESSES[@]} 
do
  for CUT in ${CUTS[@]}
  do
    for VAR in ${VARS[@]}
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
      echo Merging to output file: $INFILENAME".root"

      hadd -f ${OUTDIR}/${INFILENAME}".root" "${INDIR}/${INFILENAME}"*.root
    done
  done
done


