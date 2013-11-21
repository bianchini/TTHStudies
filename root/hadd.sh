#! /bin/sh

NAME=$2

if ls ./MEAnalysis${NAME}_$1_p* &> /dev/null ; then
    ls MEAnalysis${NAME}_$1_p*.root
    hadd -f MEAnalysis${NAME}_$1.root MEAnalysis${NAME}_$1_p*.root 
    rm  MEAnalysis${NAME}_$1_p*.root
else
    echo "MEAnalysis${NAME}_$1 not found"
fi