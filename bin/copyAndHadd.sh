#! /bin/sh

if [ -z $ROOTSYS ]; then
 echo "ROOTSYS is not defined: source ROOT, or hadd won't work!"
 exit
fi

INDIR='/pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_V3_tmp/'
OUTDIR='/pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_V3/'
ACCESS='srm://t3se01.psi.ch:8443/srm/managerv2?SFN='


data=$(ls ${INDIR}*${1}*)

if !(ls /scratch/bianchi/ &> /dev/null) ; then
    mkdir /scratch/bianchi/
fi

for name in ${data}
do
  newname=`echo $name | sed 's/\/pnfs\/psi.ch\/cms\/trivcat\/store\/user\/bianchi\/HBB_EDMNtuple\/AllHDiJetPt_V3_tmp\///'`
  echo "Copy: " ${ACCESS}"/"${INDIR}"/"${newname}
  echo " -->  " /scratch/bianchi/${newname}
  echo ""
  COPY="lcg-cp -b -D srmv2 ${ACCESS}"/"${INDIR}"/"${newname} /scratch/bianchi/${newname}"
  #echo ${COPY}
  eval "$COPY"
  #srmcp -2 ${ACCESS}"/"${INDIR}"/"${newname} file:////scratch/bianchi/${newname}
done

echo "hadd -f /scratch/bianchi/DiJetPt_$1.root /scratch/bianchi/DiJetPt_$1_*.root"
hadd -f /scratch/bianchi/DiJetPt_$1.root /scratch/bianchi/DiJetPt_$1_*.root

echo "rm /scratch/bianchi/DiJetPt_$1_*.root"
rm /scratch/bianchi/DiJetPt_$1_*.root

echo "srmrm "${ACCESS}"/"${OUTDIR}"/"DiJetPt_$1.root
srmrm ${ACCESS}"/"${OUTDIR}"/"DiJetPt_$1.root

echo "srmcp -2 file:////scratch/bianchi/DiJetPt_$1.root "${ACCESS}"/"${OUTDIR}"/"DiJetPt_$1.root 
srmcp -2 file:////scratch/bianchi/DiJetPt_$1.root ${ACCESS}"/"${OUTDIR}"/"DiJetPt_$1.root

echo "rm /scratch/bianchi/DiJetPt_$1.root"
rm /scratch/bianchi/DiJetPt_$1.root
