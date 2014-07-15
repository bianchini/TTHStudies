
#! /bin/sh

#tail -n 30 slurm-* | grep "Job done in" | wc -l

files=$(ls job*.o*)

for file in ${files}
do
<<<<<<< HEAD
  is=`tail -n 50 $file | grep "Job done in"`
  #is=`tail -n 50 $file | grep "Missing histogram for a process. Return."`
  #is=`tail -n 50 $file | grep "Recap"`
=======
  #is=`tail -n 50 $file | grep "Job done in"`
  is=`tail -n 50 $file | grep "Recap: "`
>>>>>>> b326668b738ca563b1461f3928d2c0a86ab1b39a
  if [ -z "$is" ]; then
      echo $file
      echo `cat $file | grep "cannot remove"` >> log.txt
  fi
done