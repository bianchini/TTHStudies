#!/usr/bin/env python


import commands
import re
import os

import sys
sys.path.append('./')

import FWCore.ParameterSet.Config as cms

#### mass scan
massesH     = cms.vdouble(125)
#massesH     = cms.vdouble(55., 65., 75., 85., 95., 105., 115., 125., 135., 145., 155., 165., 185., 205., 225., 250., 275., 300.)
massesT     = cms.vdouble(174.3)
MH          = 125.00
MT          = 174.30

# 0 = all events between evLow and evHigh will be considered
# 1 = process (evHigh-evLow) events of the desired type starting from the beginning
fixNumEvJob = 0

#### flags
norm        = 0
useME       = 1
useJac      = 1
useMET      = 1
useTF       = 1
usePDF      = 1

# to print intermediate steps
printout    = 1

# cut values to select events
btag_prob_cut_6jets = 0.988
btag_prob_cut_5jets = 0.992
btag_prob_cut_4jets = 0.992

# regression
useRegression = 0


def submitMEAnalysisNew(type,
                        script,
                        sample,
                        doCSVup, doCSVdown, doJECup, doJECdown,
                        evLow,evHigh
                        ):

    print "Overload meAnalysisNew_",type,".py..."
    os.system('cp ../python/meAnalysisNew_'+type+'.py ./')

    if   re.search("VType0",   type )!=None:
        from meAnalysisNew_VType0 import process, VType
    elif re.search("VType1",   type )!=None:
        from meAnalysisNew_VType1 import process, VType
    elif re.search("VType2",   type )!=None:
        from meAnalysisNew_VType2 import process, VType
    elif re.search("VType3",   type )!=None:
        from meAnalysisNew_VType3 import process, VType
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
            
    process.fwliteInput.outFileName      = cms.string('../root/MEAnalysisNew_'+script+'.root')
    process.fwliteInput.pathToFile       = cms.string("dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store//user/bianchi/HBB_EDMNtuple/AllHDiJetPt_V2/")
    if useRegression:
        process.fwliteInput.pathToTF         = cms.string("./root/transferFunctionsTEST_reg.root"),
        process.fwliteInput.pathToCP         = cms.string("./root/ControlPlotsTEST_reg.root"),
    else:
        process.fwliteInput.pathToTF         = cms.string("./root/transferFunctionsTEST.root"),
        process.fwliteInput.pathToCP         = cms.string("./root/ControlPlotsTEST_reg.root"),
    process.fwliteInput.norm             = cms.untracked.int32(norm)
    process.fwliteInput.useME            = cms.int32(useME)
    process.fwliteInput.useJac           = cms.int32(useJac)
    process.fwliteInput.useMET           = cms.int32(useMET)
    process.fwliteInput.useTF            = cms.int32(useTF)
    process.fwliteInput.usePDF           = cms.int32(usePDF)

    process.fwliteInput.btag_prob_cut_6jets = cms.untracked.double(btag_prob_cut_6jets)
    process.fwliteInput.btag_prob_cut_5jets = cms.untracked.double(btag_prob_cut_5jets)
    process.fwliteInput.btag_prob_cut_4jets = cms.untracked.double(btag_prob_cut_4jets)

    process.fwliteInput.useRegression    = cms.untracked.int32(useRegression)

    process.fwliteInput.massesH          = massesH
    process.fwliteInput.massesT          = massesT
    process.fwliteInput.MH               = cms.untracked.double(MH)
    process.fwliteInput.MT               = cms.untracked.double(MT)
    process.fwliteInput.fixNumEvJob      = cms.untracked.int32(fixNumEvJob)
    process.fwliteInput.evLimits         = cms.vint32(evLow,evHigh)

    process.fwliteInput.printout         = cms.int32(printout)

    process.fwliteInput.doCSVup          = cms.untracked.int32(doCSVup)
    process.fwliteInput.doCSVdown        = cms.untracked.int32(doCSVdown)
    process.fwliteInput.doJECup          = cms.untracked.int32(doJECup)
    process.fwliteInput.doJECdown        = cms.untracked.int32(doJECdown)


    out = open(jobName+'.py','w')
    out.write(process.dumpPython())
   
    f = open(scriptName,'w')
    f.write('#!/bin/bash\n\n')
    f.write('cd ${CMSSW_BASE}/src/Bianchi/TTHStudies/bin/\n')
    f.write('source /swshare/psit3/etc/profile.d/cms_ui_env.sh\n')
    f.write('export SCRAM_ARCH="slc5_amd64_gcc462"\n')
    f.write('source $VO_CMS_SW_DIR/cmsset_default.sh\n')
    f.write('eval `scramv1 runtime -sh`\n')
    f.write('export LD_PRELOAD="libglobus_gssapi_gsi_gcc64pthr.so.0":${LD_PRELOAD}\n')
    f.write('\n\n')
    f.write('\n\n')
    f.write('MEAnalysisNew ./'+jobName+'.py\n')
    f.close()
    os.system('chmod +x '+scriptName)

    submitToQueue = 'qsub -V -cwd -l h_vmem=2G -q all.q -N job'+sample+' '+scriptName 
    print submitToQueue
    os.system(submitToQueue)
    
    print "\n@@@@@ END JOB @@@@@@@@@@@@@@@"
        
    #os.system('rm testME_tmp.py')
  
###########################################
###########################################



def submitFullMEAnalysisNew( type, analysis ):

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
    
    if   re.search("csvUp",   analysis )!=None:
        doCSVup   = 1
    elif re.search("csvDown", analysis )!=None:
        doCSVdown = 1
    elif re.search("JECUp",   analysis )!=None:
        doJECup   = 1
    elif re.search("JECDown", analysis )!=None:
        doJECdown = 1
    else:
        print "Doing nominal analysis"

   # TTJetsSemiLept --> 16749255
    sample  = 'TTJetsSemiLept'
    counter = 0
    for i in range(300):   # ---> ~ 40/job
        counter = counter + 1
        if doSL:
            submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown,i*56000+1, (i+1)*56000 )
            

    # TTJetsFullLept
    sample  = 'TTJetsFullLept'
    counter = 0
    for i in range(12):
        counter = counter + 1
        if doDL:
            submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown,i*20000+1, (i+1)*20000 )

    # TTJetsFullLept
    sample  = 'TTJetsFullLept'
    counter = 0
    for i in range(12):
        counter = counter + 1
        if doSL:
            submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown,i*30000+1, (i+1)*30000 )

    # TTH125 --> 194808
    sample  = 'TTH125'
    counter = 0
    num_of_jobs = 1
    evs_per_job = 100
    if doSL:
         num_of_jobs =  195
         evs_per_job = 1000  # ---> ~ 40/job
    else:
         num_of_jobs =    10
         evs_per_job = 21000 # ---> ~ 40/job
    for i in range(num_of_jobs):
        counter = counter + 1
        submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown, i*evs_per_job+1, (i+1)*evs_per_job )


    if re.search("nominal",   analysis )==None:
        return


    # DYJets10to50
    sample  = 'DYJets10to50'
    counter = 0
    for i in range(1):
        counter = counter + 1
        submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown,0, -1 )

    # DYJets50
    sample  = 'DYJets50'
    counter = 0
    for i in range(1):
        counter = counter + 1
        submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown,0, -1 )

    # WJets 
    sample  = 'WJets'
    counter = 0
    for i in range(1):
        counter = counter + 1
        submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown,0, -1 )
    
    # TtW
    sample  = 'TtW'
    counter = 0
    for i in range(1):
        counter = counter + 1
        submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown,0, -1 )
        
   # Tt
    sample  = 'Tt'
    counter = 0
    for i in range(1):
        counter = counter + 1
        submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown,0, -1 )

   # Ts
    sample  = 'Ts'
    counter = 0
    for i in range(1):
        counter = counter + 1
        submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown,0, -1 )

   # TbartW
    sample  = 'TbartW'
    counter = 0
    for i in range(1):
        counter = counter + 1
        submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown,0, -1 )

   # Tbart
    sample  = 'Tbart'
    counter = 0
    for i in range(1):
        counter = counter + 1
        submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown,0, -1 )

   # Tbars
    sample  = 'Tbars'
    counter = 0
    for i in range(1):
        counter = counter + 1
        submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown,0, -1 )

   # WW
    sample  = 'WW'
    counter = 0
    for i in range(1):
        counter = counter + 1
        submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown,0, -1 )

   # WZ
    sample  = 'WZ'
    counter = 0
    for i in range(1):
        counter = counter + 1
        submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown,0, -1 )

   # ZZ
    sample  = 'ZZ'
    counter = 0
    for i in range(1):
        counter = counter + 1
        submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown,0, -1 )
    
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
        submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown,i*evs_per_job+1, (i+1)*evs_per_job )

   # TTW
    sample  = 'TTW'
    counter = 0
    for i in range(1):
        counter = counter + 1
        submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown,0, -1 )

###########################################
###########################################


analyses = ['nominal_v2',
            #'csvUp_v1',
            #'csvDown_v1',
            #'JECUp_v1',
            #'JECDown_v1',
            ]

types = ['VType0',
         'VType1',
         'VType2',
         'VType3',
         ]

for type in types:
    for analysis in analyses:
        submitFullMEAnalysisNew(type, analysis)
