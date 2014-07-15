#! /bin/sh

data=$(ls MEAnalysisNew_*csvVH*)

for name in ${data}
do
    new="$(echo $name | sed 's/csvVH_rec_std/ntuplizeAll_csvVH_rec_std/')"
    echo $name "-->" $new
    mv $name $new
done