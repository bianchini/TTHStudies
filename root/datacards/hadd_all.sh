
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

for PROC in ${PROCESSES[@]} 
do
  echo "Start processing: "$PROC
  if ls ${INDIR}/*_${PROC}_*.root &> /dev/null ; then # check whether input files exist
      INFILENAMES=(`ls ${INDIR}/*_${PROC}_*.root`)
  else
      echo $PROC input files not found
      continue
  fi

  INFILENAME=$(basename `echo ${INFILENAMES[0]} | cut -d'.' --complement -f2-`) #drop .root extension and path
  INFILENAME=${INFILENAME%${PROC}*}"$PROC" # drop _xx after process name
  echo Merging to output file: $INFILENAME".root"

  hadd -f ${OUTDIR}/${INFILENAME}".root" ${INDIR}/${INFILENAME}*.root

done

