#! /bin/sh

data=$(ls /pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_V2/ | grep Run)
count=1

for name in ${data}
do
  count=$((count+1))
  new="$(echo $name | sed 's/DiJetPt_//')"
  new2="$(echo $new | sed 's/.root//')"
  echo "'"$new2"',"
  #echo "    cms.PSet("
  #echo "    skip     = cms.bool(True), " 
  #echo "    name     = cms.string('${new2}'+VType),"
  #echo "    nickName = cms.string('Run2012_${new2}'),"
  #echo "    color    = cms.int32(41),"
  #echo "    xSec     = cms.double(-1),"
  #echo "    ),"
done