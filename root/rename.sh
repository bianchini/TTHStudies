#! /bin/sh

data=$(ls MEAnalysisNew_MHscan_SL_VType2*v4*)

for name in ${data}
do
    new="$(echo $name | sed 's/SL_VType2/SL/')"
    echo $name "-->" $new
    mv $name $new
done