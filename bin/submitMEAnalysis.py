#!/usr/bin/env python


import commands
import re
import os

import sys
sys.path.append('./')

import FWCore.ParameterSet.Config as cms


doTest     = False
doAnalysis = True

doSL = True
doDL = False

type = "DUMMY"
if doDL:
    type = "_VType1"
if doSL:
    type = "_VType2"

masses = cms.vdouble(125)

analysis = "nominal_v5"
#analysis = "csvUp_v5"
#analysis = "csvDown_v5"
#analysis = "JECUp_v5"
#analysis = "JECDown_v5"

doCSVup    =  cms.untracked.int32(0)
doCSVdown  =  cms.untracked.int32(0)
doJECup    =  cms.untracked.int32(0)
doJECdown  =  cms.untracked.int32(0)

def submitMEAnalysis(script,
                     sample,
                     norm,
                     useME, useJac, useMET, useTF, usePDF,
                     met,
                     masses,
                     evLow,evHigh
                     ):

    print "Overload meAnalysis.py..."
    os.system('cp ../python/meAnalysis.py ./')

    from meAnalysis import process, VType

    if type!=VType:
        print "Attention! Mismatch between VType... Please, check and run again."
        return

    print "Creating the shell file for the batch..."
    scriptName = 'job_'+script+'.sh'
    jobName    = 'job_'+script

    for sam in process.fwliteInput.samples:
        if sam.nickName != sample:
            sam.skip = cms.bool(True)
        else:
            sam.skip = cms.bool(False)
            
    process.fwliteInput.outFileName      = cms.string('../root/MEAnalysis_'+script+'.root')
    process.fwliteInput.norm             = cms.untracked.int32(norm)
    process.fwliteInput.useME            = cms.int32(useME)
    process.fwliteInput.useJac           = cms.int32(useJac)
    process.fwliteInput.useMET           = cms.int32(useMET)
    process.fwliteInput.useTF            = cms.int32(useTF)
    process.fwliteInput.usePDF           = cms.int32(usePDF)

    process.fwliteInput.masses           = masses
    process.fwliteInput.evLimits         = cms.vint32(evLow,evHigh)


    process.fwliteInput.doCSVup          = doCSVup
    process.fwliteInput.doCSVdown        = doCSVdown
    process.fwliteInput.doJECup          = doJECup
    process.fwliteInput.doJECdown        = doJECdown


    out = open(jobName+'.py','w')
    out.write(process.dumpPython())
   
    f = open(scriptName,'w')
    f.write('#!/bin/bash\n\n')
    f.write('cd /shome/bianchi/CMSSW_5_3_3_patch2/src/Bianchi/TTHStudies/bin/\n')
    f.write('source /swshare/psit3/etc/profile.d/cms_ui_env.sh\n')
    f.write('export SCRAM_ARCH="slc5_amd64_gcc462"\n')
    f.write('source $VO_CMS_SW_DIR/cmsset_default.sh\n')
    f.write('eval `scramv1 runtime -sh`\n')
    f.write('export LD_PRELOAD="libglobus_gssapi_gsi_gcc64pthr.so.0":${LD_PRELOAD}\n')
    f.write('\n\n')
    f.write('\n\n')
    f.write('MEAnalysis ./'+jobName+'.py\n')
    f.close()
    os.system('chmod +x '+scriptName)

    submitToQueue = 'qsub -V -cwd -l h_vmem=6G -q all.q -N '+script+' '+scriptName 
    print submitToQueue
    os.system(submitToQueue)
    
    print "\n@@@@@ END JOB @@@@@@@@@@@@@@@"
        
    #os.system('rm testME_tmp.py')
  
###########################################
###########################################



if doTest:
    counter = 0
    for i in range(1):
        counter = counter + 1
        submitMEAnalysis('SL4X_type2new_rec_SoB_allpermut_bkg_p'+str(counter), "TTH", 0, 1,1,1,1,1,  125,  masses,  i*20+1, (i+1)*20 )



if doAnalysis:

    # DYJets10to50
    sample  = 'DYJets10to50'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysis('SL'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1,  125,  masses,  0, -1 )

    # DYJets50
    sample  = 'DYJets50'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysis('SL'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1,  125,  masses,  0, -1 )

    # WJets 
    sample  = 'WJets'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysis('SL'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1,  125,  masses,  0, -1 )
    
    # TtW
    sample  = 'TtW'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysis('SL'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1,  125,  masses,  0, -1 )
        
   # Tt
    sample  = 'Tt'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysis('SL'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1,  125,  masses,  0, -1 )

   # Ts
    sample  = 'Ts'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysis('SL'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1,  125,  masses,  0, -1 )

   # TbartW
    sample  = 'TbartW'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysis('SL'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1,  125,  masses,  0, -1 )

   # Tbart
    sample  = 'Tbart'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysis('SL'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1,  125,  masses,  0, -1 )

   # Tbars
    sample  = 'Tbars'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysis('SL'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1,  125,  masses,  0, -1 )

   # WW
    sample  = 'WW'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysis('SL'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1,  125,  masses,  0, -1 )

   # WZ
    sample  = 'WZ'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysis('SL'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1,  125,  masses,  0, -1 )

   # ZZ
    sample  = 'ZZ'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysis('SL'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1,  125,  masses,  0, -1 )
    
    # TTH125
    sample  = 'TTH125'
    counter = 0
    num_of_jobs = 1
    evs_per_job = 100
    if doSL:
         num_of_jobs = 110
         evs_per_job = 500
    else:
         num_of_jobs = 8
         evs_per_job = 300
    for i in range(num_of_jobs):
        counter = counter + 1
        #submitMEAnalysis('SL'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1,  125,  masses,   i*evs_per_job+1, (i+1)*evs_per_job )

    # TTZ
    sample  = 'TTZ'
    counter = 0
    num_of_jobs = 1
    evs_per_job = 1000
    if doSL:
         num_of_jobs = 7
         evs_per_job = 3000
    else:
         num_of_jobs = 3
         evs_per_job = 500
    for i in range(num_of_jobs):
        counter = counter + 1
        if counter!=7:
            continue
        submitMEAnalysis('SL'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1,  125,  masses,  i*evs_per_job+1, (i+1)*evs_per_job )

   # TTW
    sample  = 'TTW'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysis('SL'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1,  125,  masses,  0, -1 )

   # TTJetsSemiLept
    sample  = 'TTJetsSemiLept'
    counter = 0
    for i in range(150):
        counter = counter + 1
        #if doSL:
        #    submitMEAnalysis('SL'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1,  125,  masses,  i*25000+1, (i+1)*25000 )
            

    # TTJetsFullLept
    sample  = 'TTJetsFullLept'
    counter = 0
    for i in range(12):
        counter = counter + 1
        if doDL:
            submitMEAnalysis('SL'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1,  125,  masses,  i*20000+1, (i+1)*20000 )
