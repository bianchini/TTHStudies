#! /bin/sh


combineCards.py Name1=SL_cat1_$1.txt Name2=SL_cat2_$1.txt Name3=SL_cat3_$1.txt Name4=SL_cat4_$1.txt Name5=SL_cat5_$1.txt> SL_$1.txt
combineCards.py Name1=DL_cat6_$1.txt > DL_$1.txt 
combineCards.py Name1=SL_$1.txt Name2=DL_$1.txt > COMB_$1.txt

#combineCards.py Name1=SL_cat1_$1.txt Name2=SL_cat2_$1.txt Name3=DL_cat4_$1.txt > SL_$1.txt    
#combine -M Asymptotic -t -1 -n SL_cat1_$1 -d SL_cat1_$1.txt
#combine -M Asymptotic -t -1 -n SL_cat2_$1 -d SL_cat2_$1.txt
#combine -M Asymptotic -t -1 -n DL_cat4_$1 -d DL_cat4_$1.txt
#combine -M Asymptotic -t -1 -n SL_$1      -d SL_$1.txt     

combine -M Asymptotic -t -1 -n COMB_$1     -d COMB_$1.txt 
combine -M Asymptotic -t -1 -n SL_$1       -d SL_$1.txt 
combine -M Asymptotic -t -1 -n DL_$1       -d DL_$1.txt 

combine -M Asymptotic -t -1 -n SL_cat1_$1  -d SL_cat1_$1.txt 
combine -M Asymptotic -t -1 -n SL_cat2_$1  -d SL_cat2_$1.txt 
combine -M Asymptotic -t -1 -n SL_cat3_$1  -d SL_cat3_$1.txt 
combine -M Asymptotic -t -1 -n SL_cat4_$1  -d SL_cat4_$1.txt 
combine -M Asymptotic -t -1 -n SL_cat5_$1  -d SL_cat5_$1.txt 
combine -M Asymptotic -t -1 -n DL_cat6_$1  -d DL_cat6_$1.txt 
