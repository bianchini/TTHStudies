#! /bin/sh

data=$(ls MEAnalysisNew_SL_VType0*)

for name in ${data}
do
    new="$(echo $name | sed 's/SL_VType0/DL/')"
    echo $name "-->" $new
    mv $name $new
done