#!/usr/bin/env python


import commands
import re
import os
import FWCore.ParameterSet.Config as cms

import sys
sys.path.append('./')

from ntuple  import process 

def processAllBatch(jobName, outName, split):

    process.fwliteInput.fileNames = ()
    
    input = open('fileList_'+jobName+'.txt','r')
    counter = 1
    for line in input:
        foo = line.strip('\n')
        if counter>=split[0] and counter<=split[1]:
            process.fwliteInput.fileNames.append('dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/'+foo)
        counter = counter+1

    input.close()
    process.fwliteOutput.fileName = cms.string(outName+'_'+str(split[0])+'-'+str(split[1])+'.root')

    out = open('ntuple_'+jobName+'_'+str(split[0])+'-'+str(split[1])+'.py','w')
    out.write(process.dumpPython())
    out.close()

    scriptName = 'myStep2_'+jobName+'_'+str(split[0])+'-'+str(split[1])+'.sh'
    
    f = open(scriptName,'w')
    f.write('#!/bin/bash\n\n')
    f.write('cd $HOME/CMSSW_5_3_3_patch2/src/VHbbAnalysis/VHbbDataFormats/bin\n')
    f.write('source /swshare/psit3/etc/profile.d/cms_ui_env.sh\n')
    f.write('export SCRAM_ARCH="slc5_amd64_gcc462"\n')
    f.write('source $VO_CMS_SW_DIR/cmsset_default.sh\n')
    f.write('eval `scramv1 runtime -sh`\n')
    f.write('export LD_PRELOAD="libglobus_gssapi_gsi_gcc64pthr.so.0":${LD_PRELOAD}\n')
    f.write('export LD_LIBRARY_PATH=/swshare/glite/globus/lib/:/swshare/glite/d-cache/dcap/lib64/:$LD_LIBRARY_PATH\n')
    f.write('\n\n')
    f.write('\n\n')
    f.write('Ntupler ntuple_'+jobName+'_'+str(split[0])+'-'+str(split[1])+'.py')
    f.close()

    os.system('chmod +x '+scriptName)
    submitToQueue = 'qsub -V -cwd -l h_vmem=6G -q all.q -N '+jobName+'_'+str(split[0])+'-'+str(split[1])+' '+scriptName 
    print submitToQueue
    os.system(submitToQueue)


###########################################################################
###########################################################################


#processAllBatch("WJets", "DiJetPt_WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball", [1, 500])
for k in range(22):
    #processAllBatch("WJets", "DiJetPt_WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball", [k*200+1,(k+1)*200])
    processAllBatch("DYJets-10to50", "DiJetPt_DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph", [k*200+1,(k+1)*200])


#processAllBatch("TTH110", "DiJetPt_TTH_HToBB_M-110_8TeV-pythia6")
#processAllBatch("TTH115", "DiJetPt_TTH_HToBB_M-115_8TeV-pythia6")
#processAllBatch("TTH120", "DiJetPt_TTH_HToBB_M-120_8TeV-pythia6")
#processAllBatch("TTH125", "DiJetPt_TTH_HToBB_M-125_8TeV-pythia6")
#processAllBatch("TTH130", "DiJetPt_TTH_HToBB_M-130_8TeV-pythia6")
#processAllBatch("TTH135", "DiJetPt_TTH_HToBB_M-135_8TeV-pythia6")
#processAllBatch("TTW",    "DiJetPt_TTWJets_8TeV-madgraph")
#processAllBatch("TTZ",    "DiJetPt_TTZJets_8TeV-madgraph")

