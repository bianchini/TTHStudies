#! /bin/sh

#tail -n 30 slurm-* | grep "Job done in" | wc -l

files=$(ls slurm-*)

for file in ${files}
do
  is=`tail -n 50 $file | grep "Job done in"`
  if [ -z "$is" ]; then
      echo $file
      echo `cat $file | grep "cannot remove"` >> log.txt
  fi
done