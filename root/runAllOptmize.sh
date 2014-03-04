#! /bin/sh


max=6
for i in `seq 0 $max`
  do
  echo "Doing $1_$2_$3_$i"
  combine -M Asymptotic -t -1 -n $1_$2_$3_$i  -d $1_$2_$3_$i.txt 
done
