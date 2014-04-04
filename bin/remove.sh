#! /bin/sh

INDIR="obsolete/"

data=$(ls ${INDIR}*${1}*)

for name in ${data}
do
  newname=`echo $name | sed 's/obsolete\///'`
  echo "cp ${INDIR}/${newname} ."
  cp ${INDIR}/${newname} .
  echo "git rm ${newname}"
  git rm ${newname}
done
