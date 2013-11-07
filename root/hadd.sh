#! /bin/sh

NAME=New_MHscan

ls MEAnalysis${NAME}_$1_p*.root
hadd -f MEAnalysis${NAME}_$1.root MEAnalysis${NAME}_$1_p*.root 
rm  MEAnalysis${NAME}_$1_p*.root