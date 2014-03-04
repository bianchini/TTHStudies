#! /bin/sh

OPTION='Asymptotic -t -1'
#OPTION='MaxLikelihoodFit '

combineCards.py Name1=MEM_cat1_$1.txt Name2=MEM_cat2_$1.txt Name3=MEM_cat3_$1.txt Name4=MEM_cat4_$1.txt Name5=MEM_cat5_$1.txt> MEM_SL_$1.txt

combineCards.py Name1=MEM_cat6_$1.txt > MEM_DL_$1.txt 

combineCards.py Name1=MEM_SL_$1.txt Name2=MEM_DL_$1.txt > MEM_COMB_$1.txt


combine -M ${OPTION} -n MEM_COMB_$1  -d MEM_COMB_$1.txt 
combine -M ${OPTION} -n MEM_SL_$1    -d MEM_SL_$1.txt 
combine -M ${OPTION} -n MEM_DL_$1    -d MEM_DL_$1.txt 

combine -M ${OPTION} -n MEM_cat1_$1  -d MEM_cat1_$1.txt 
combine -M ${OPTION} -n MEM_cat2_$1  -d MEM_cat2_$1.txt 
combine -M ${OPTION} -n MEM_cat3_$1  -d MEM_cat3_$1.txt 
combine -M ${OPTION} -n MEM_cat4_$1  -d MEM_cat4_$1.txt 
combine -M ${OPTION} -n MEM_cat5_$1  -d MEM_cat5_$1.txt 
combine -M ${OPTION} -n MEM_cat6_$1  -d MEM_cat6_$1.txt 
