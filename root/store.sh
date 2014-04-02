#! /bin/sh


VERSION="all_CSVcalibration_rec_std"

data=$(ls *${VERSION}*)
outdir="srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/bianchi/Trees/MEM/"

for name in ${data}
do
  echo "Copy: " $name 
  echo " -->  " ${outdir}"/"$1"/"${name}
  srmcp -2  file:///$name  ${outdir}"/"$1"/"${name}
done