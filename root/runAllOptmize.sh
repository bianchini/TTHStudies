#! /bin/sh


OPTIMIZE1=1
OPTIMIZE2=0

if [ $OPTIMIZE1==1 ] ; then
    max=8
    for i in `seq 0 $max`
      do
      echo "Doing $1_$2_$3_$i"
      combineCards.py Name1=$1_$2_$3-${i}_0.txt  Name2=$1_$2_$3-${i}_1.txt > $1_$2_$3-${i}.txt
      combine -M Asymptotic -t -1 -n $1_$2_$3_$i  -d $1_$2_$3-$i.txt 
    done
fi

if [ $OPTIMIZE2==1 ] ; then
    max=8
    for i in `seq 0 $max`
      do
      echo "Doing $1_$2_$3_$i"
      combine -M Asymptotic -t -1 -n $1_$2_$3_$i  -d $1_$2_$3_$i.txt 
    done
fi
