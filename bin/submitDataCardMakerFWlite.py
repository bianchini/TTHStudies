#!/usr/bin/env python


import commands
import re
import os
import ROOT

import sys
sys.path.append('./')

import FWCore.ParameterSet.Config as cms

from Bianchi.TTHStudies.mem_categories_cff import *

###########################################
###########################################
 

def submitDataCardMakerFWlite_Limits(category):

    print "Overload datacardMakerFWlite.py..."
    os.system('cp ../python/datacardMakerFWlite.py ./')

    from datacardMakerFWlite import process

    if   category == "cat1_sb":
        process.fwliteInput = cat1_sb
    elif category == "cat2_sb":
        process.fwliteInput = cat2_sb
    elif category == "cat3_sb":
        process.fwliteInput = cat3_sb
    elif category == "cat6_sb":
        process.fwliteInput = cat6_sb
    elif category == "cat1_sb_nb":
        process.fwliteInput = cat1_sb_nb
    elif category == "cat2_sb_nb":
        process.fwliteInput = cat2_sb_nb
    elif category == "cat3_sb_nb":
        process.fwliteInput = cat3_sb_nb
    elif category == "cat6_sb_nb":
        process.fwliteInput = cat6_sb_nb
    elif category == "cat1_bj":
        process.fwliteInput = cat1_bj
    elif category == "cat2_bj":
        process.fwliteInput = cat2_bj
    elif category == "cat3_bj":
        process.fwliteInput = cat3_bj
    elif category == "cat6_bj":
        process.fwliteInput = cat6_bj
    else:
        print "Cannot find this category... exit"
        return

    process.fwliteInput.extraname = cms.string( process.fwliteInput.extraname.value()+"_"+category)
        
    print "Creating the shell file for the batch..."
    scriptName = 'job_'+category+'.sh'
    jobName    = 'job_'+category

    out = open(jobName+'.py','w')
    out.write(process.dumpPython())
   
    f = open(scriptName,'w')
    f.write('#!/bin/bash\n\n')
    f.write('cd ${CMSSW_BASE}/src/Bianchi/TTHStudies/bin/\n')
    f.write('export SCRAM_ARCH="slc5_amd64_gcc462"\n')
    f.write('source $VO_CMS_SW_DIR/cmsset_default.sh\n')
    f.write('eval `scramv1 runtime -sh`\n')
    f.write('\n\n')
    f.write('\n\n')
    f.write('DataCardMakerNewFWlite ./'+jobName+'.py\n')
    f.write('rm /scratch/bianchi/dummy*.root\n')
    f.close()
    os.system('chmod +x '+scriptName)

    submitToQueue = 'qsub -V -cwd -l h_vmem=6G -q all.q -N job'+category+' '+scriptName 
    print submitToQueue
    #os.system(submitToQueue)
    
    print "\n@@@@@ END JOB @@@@@@@@@@@@@@@"

###########################################
###########################################
 

def submitDataCardMakerFWlite(category, cut, script, samples, extraname, nparts, part):

    print "Overload datacardMakerFWlite.py..."
    os.system('cp ../python/datacardMakerFWlite.py ./')

    from datacardMakerFWlite import process

    process.fwliteInput.samples   = samples    
    process.fwliteInput.category  = cms.string( category )
    process.fwliteInput.cut       = cms.string( cut )
    process.fwliteInput.extraname = cms.string( extraname )
    process.fwliteInput.nparts    = cms.int32(nparts)
    process.fwliteInput.part      = cms.int32(part)
    
    print "Creating the shell file for the batch..."
    scriptName = 'job_'+script+'.sh'
    jobName    = 'job_'+script

    out = open(jobName+'.py','w')
    out.write(process.dumpPython())
   
    f = open(scriptName,'w')
    f.write('#!/bin/bash\n\n')
    f.write('cd ${CMSSW_BASE}/src/Bianchi/TTHStudies/bin/\n')
    f.write('export SCRAM_ARCH="slc5_amd64_gcc462"\n')
    f.write('source $VO_CMS_SW_DIR/cmsset_default.sh\n')
    f.write('eval `scramv1 runtime -sh`\n')
    f.write('\n\n')
    f.write('\n\n')
    f.write('DataCardMakerNewFWlite ./'+jobName+'.py\n')
    f.write('rm /scratch/bianchi/dummy*.root\n')
    f.close()
    os.system('chmod +x '+scriptName)

    submitToQueue = 'qsub -V -cwd -l h_vmem=6G -q all.q -N job'+category+' '+scriptName 
    print submitToQueue
    os.system(submitToQueue)
    
    print "\n@@@@@ END JOB @@@@@@@@@@@@@@@"
        
###########################################
###########################################



def submitDataCardMakerFWlite_all( category , cut, selection):



    sampless = [ [["TTV"],     1],
                 [["SingleT"], 1],
                 [["DiBoson"], 1],
                 [["TTJetsBB"],22],
                 [["TTJetsBJ"],22],
                 [["TTJetsJJ"],50],
                 [["TTH125"],  1],
                 [["EWK"],     1],
                 #[["Run2012_SingleMu", "Run2012_SingleElectron"],1 ]
                 [["Run2012_SingleElectron"],1 ]
                 ]

    counter = 0
    for samples in sampless:
        toBeRun = cms.vstring( samples[0] )
        outfile = ""
        for sample in samples[0]:
            if sample != (samples[0])[len(samples[0])-1]:
                outfile += (sample+"-")
            else:
                outfile += (sample)
        for split in range( samples[1] ):
            counter += 1
            script  = category+"_p"+str(counter)
            submitDataCardMakerFWlite( category, cut, script, toBeRun, "_"+selection+"_"+outfile+"_"+str(split), samples[1], split)
      
      
###########################################
###########################################

#submitDataCardMakerFWlite_all( "lepton_pt", "(numJets>=6 && numBTagM==2 && Vtype==3)", "6j2t_ele" )

submitDataCardMakerFWlite_Limits("cat1_sb")
#submitDataCardMakerFWlite_Limits("cat2_sb")
#submitDataCardMakerFWlite_Limits("cat3_sb")
#submitDataCardMakerFWlite_Limits("cat6_sb")


