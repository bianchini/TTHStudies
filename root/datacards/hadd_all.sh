
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
"Run2012_SingleMu"
)

CUTS=(
"SL_g6j2t"
"SL_4j3t"
"SL_5j3t"
"SL_g6j3t"
"SL_4j4t"
"SL_5jg4t"
"SL_g6jg4t"
#"DL_4j2t"
#"DL_4j4t"
)

for PROC in ${PROCESSES[@]} 
do
  for CUT in ${CUTS[@]}
  do
    echo "Start processing: "$PROC
    if ls ${INDIR}/*${CUT}*_${PROC}_*.root &> /dev/null ; then # check whether input files exist
	INFILENAMES=(`ls ${INDIR}/*${CUT}*_${PROC}_*.root`)
    else
	echo $PROC input files not found
	continue
    fi

    INFILENAME=$(basename `echo ${INFILENAMES[0]} | cut -d'.' --complement -f2-`) #drop .root extension and path
    INFILENAME=${INFILENAME%${PROC}*}"$PROC" # drop _xx after process name
    echo Merging to output file: $INFILENAME".root"

    hadd -f ${OUTDIR}/${INFILENAME}".root" ${INDIR}/${INFILENAME}*.root
  done
done

