#! /bin/sh

# first combine datacards
combineCards.py Name1=SL_cat1_New.txt Name2=SL_cat2_New.txt Name3=SL_cat3_New.txt Name4=SL_cat4_New.txt Name5=SL_cat5_New.txt > SL_New.txt
combineCards.py Name1=DL_cat6_New.txt > DL_New.txt 
combineCards.py Name1=SL_New.txt Name2=DL_New.txt > COMB_New.txt

# run the combined limit 
combine -M Asymptotic -t -1 -n COMB     -d COMB_New.txt 
combine -M Asymptotic -t -1 -n SL       -d SL_New.txt 
combine -M Asymptotic -t -1 -n DL       -d DL_New.txt 

# run the individual limits
combine -M Asymptotic -t -1 -n SL_cat1  -d SL_cat1_New.txt 
combine -M Asymptotic -t -1 -n SL_cat2  -d SL_cat2_New.txt 
combine -M Asymptotic -t -1 -n SL_cat3  -d SL_cat3_New.txt 
combine -M Asymptotic -t -1 -n SL_cat4  -d SL_cat4_New.txt 
combine -M Asymptotic -t -1 -n SL_cat5  -d SL_cat5_New.txt 
combine -M Asymptotic -t -1 -n DL_cat6  -d DL_cat6_New.txt 
