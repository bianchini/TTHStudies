#! /bin/sh

if [ -z $ROOTSYS ]; then
 echo "ROOTSYS is not defined: source ROOT, or hadd won't work!"
 exit
fi

NAME=$2

if ls ./MEAnalysis${NAME}_$1_p* &> /dev/null ; then
    ls MEAnalysis${NAME}_$1_p*.root
    hadd -f MEAnalysis${NAME}_$1.root MEAnalysis${NAME}_$1_p*.root 
    if ls MEAnalysis${NAME}_$1.root &> /dev/null ; then
	rm MEAnalysis${NAME}_$1_p*.root
    fi
else
    echo "MEAnalysis${NAME}_$1 not found"
fi