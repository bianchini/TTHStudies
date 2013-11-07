#!/usr/bin/env python


import commands
import re
import os

import sys
sys.path.append('./')

import FWCore.ParameterSet.Config as cms


masses = cms.vdouble(125)




def submitMEAnalysis(type,
                     script,
                     sample,
                     norm,
                     useME, useJac, useMET, useTF, usePDF,
                     doCSVup, doCSVdown, doJECup, doJECdown,
                     met,
                     masses,
                     evLow,evHigh
                     ):

    print "Overload meAnalysis_",type,".py..."
    os.system('cp ../python/meAnalysis_'+type+'.py ./')

    if   re.search("VType0",   type )!=None:
        from meAnalysis_VType0 import process, VType
    elif re.search("VType1",   type )!=None:
        from meAnalysis_VType1 import process, VType
    elif re.search("VType2",   type )!=None:
        from meAnalysis_VType2 import process, VType
    elif re.search("VType3",   type )!=None:
        from meAnalysis_VType3 import process, VType
    else:
        print "Unsupported analysis..."
        return

    if VType!=("_"+type):
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


    process.fwliteInput.doCSVup          = cms.untracked.int32(doCSVup)
    process.fwliteInput.doCSVdown        = cms.untracked.int32(doCSVdown)
    process.fwliteInput.doJECup          = cms.untracked.int32(doJECup)
    process.fwliteInput.doJECdown        = cms.untracked.int32(doJECdown)


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



def submitFullMEAnalysis( type, analysis ):

    print "Running full anamlysis for ", type, " (", analysis, ")"


    doSL = False
    doDL = False

    if   (re.search("VType0",   type )!=None or re.search("VType1",   type )!=None):
        doDL = True
    elif (re.search("VType2",   type )!=None or re.search("VType3",   type )!=None):
        doSL = True
    else:
        print "Exit"
        return

    doCSVup  = 0
    doCSVdown= 0
    doJECup  = 0
    doJECdown= 0
    
    if re.search("csvUp",   analysis )!=None:
        doCSVup   = 1
    elif re.search("csvDown", analysis )!=None:
        doCSVdown = 1
    elif re.search("JECUp",   analysis )!=None:
        doJECup   = 1
    elif re.search("JECDown", analysis )!=None:
        doJECdown = 1
    else:
        print "Doing nominal analysis"

   # TTJetsSemiLept
    sample  = 'TTJetsSemiLept'
    counter = 0
    for i in range(150):
        counter = counter + 1
        #if doSL:
        #    submitMEAnalysis(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1, doCSVup, doCSVdown, doJECup, doJECdown,  125,  masses,  i*25000+1, (i+1)*25000 )
            

    # TTJetsFullLept
    sample  = 'TTJetsFullLept'
    counter = 0
    for i in range(12):
        counter = counter + 1
        #if doDL:
        #    submitMEAnalysis(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1, doCSVup, doCSVdown, doJECup, doJECdown,  125,  masses,  i*20000+1, (i+1)*20000 )

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
        #submitMEAnalysis(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1, doCSVup, doCSVdown, doJECup, doJECdown,  125,  masses,   i*evs_per_job+1, (i+1)*evs_per_job )


    if re.search("nominal",   analysis )==None:
        return


    # DYJets10to50
    sample  = 'DYJets10to50'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysis(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1, doCSVup, doCSVdown, doJECup, doJECdown,  125,  masses,  0, -1 )

    # DYJets50
    sample  = 'DYJets50'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysis(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1, doCSVup, doCSVdown, doJECup, doJECdown,  125,  masses,  0, -1 )

    # WJets 
    sample  = 'WJets'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysis(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1, doCSVup, doCSVdown, doJECup, doJECdown,  125,  masses,  0, -1 )
    
    # TtW
    sample  = 'TtW'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysis(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1, doCSVup, doCSVdown, doJECup, doJECdown,  125,  masses,  0, -1 )
        
   # Tt
    sample  = 'Tt'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysis(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1, doCSVup, doCSVdown, doJECup, doJECdown,  125,  masses,  0, -1 )

   # Ts
    sample  = 'Ts'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysis(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1, doCSVup, doCSVdown, doJECup, doJECdown,  125,  masses,  0, -1 )

   # TbartW
    sample  = 'TbartW'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysis(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1, doCSVup, doCSVdown, doJECup, doJECdown,  125,  masses,  0, -1 )

   # Tbart
    sample  = 'Tbart'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysis(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1, doCSVup, doCSVdown, doJECup, doJECdown,  125,  masses,  0, -1 )

   # Tbars
    sample  = 'Tbars'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysis(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1, doCSVup, doCSVdown, doJECup, doJECdown,  125,  masses,  0, -1 )

   # WW
    sample  = 'WW'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysis(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1, doCSVup, doCSVdown, doJECup, doJECdown,  125,  masses,  0, -1 )

   # WZ
    sample  = 'WZ'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysis(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1, doCSVup, doCSVdown, doJECup, doJECdown,  125,  masses,  0, -1 )

   # ZZ
    sample  = 'ZZ'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysis(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1, doCSVup, doCSVdown, doJECup, doJECdown,  125,  masses,  0, -1 )
    
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
        submitMEAnalysis(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1, doCSVup, doCSVdown, doJECup, doJECdown,  125,  masses,  i*evs_per_job+1, (i+1)*evs_per_job )

   # TTW
    sample  = 'TTW'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysis(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample , 0, 1,1,1,1,1, doCSVup, doCSVdown, doJECup, doJECdown,  125,  masses,  0, -1 )

###########################################
###########################################


analyses = ['nominal_v5',
            #'csvUp_v5',
            #'csvDown_v5',
            #'JECUp_v5',
            #'JECDown_v5',
            ]

types = ['VType0',
         'VType1',
         #'VType2',
         'VType3',
         ]

for type in types:
    for analysis in analyses:
        submitFullMEAnalysis(type, analysis)
