#!/usr/bin/env python


import commands
import re
import os

import sys
sys.path.append('./')

import FWCore.ParameterSet.Config as cms


def submitTest(script,
               mode,
               masses,
               useME, useJac, useMET, useTF, usePDF,
               shiftMomenta,testMassScan,testPermutations,
               scaleH,scaleL,scaleMET,
               evLow,evHigh):

    print "Overload testME.py..."
    os.system('cp ../python/testME.py ./')

    from testME import process
    
    print "Creating the shell file for the batch..."
    scriptName = 'job_'+script+'.sh'
    jobName    = 'job_'+script

    process.fwliteInput.outFileName      = cms.string('./root/TestMENew_'+script+'.root')
    process.fwliteInput.mode             = cms.untracked.int32(mode)
    process.fwliteInput.masses           = masses
    process.fwliteInput.useME            = cms.int32(useME)
    process.fwliteInput.useJac           = cms.int32(useJac)
    process.fwliteInput.useMET           = cms.int32(useMET)
    process.fwliteInput.useTF            = cms.int32(useTF)
    process.fwliteInput.usePDF           = cms.int32(usePDF)
    process.fwliteInput.shiftMomenta     = cms.int32(shiftMomenta)
    process.fwliteInput.testMassScan     = cms.int32(testMassScan)
    process.fwliteInput.testPermutations = cms.int32(testPermutations)
    process.fwliteInput.scaleH           = cms.double(scaleH)
    process.fwliteInput.scaleL           = cms.double(scaleL)
    process.fwliteInput.scaleMET         = cms.double(scaleMET)          
    process.fwliteInput.evLimits         = cms.vint32(evLow,evHigh)

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
    f.write('TestMENew ./'+jobName+'.py\n')
    f.close()
    os.system('chmod +x '+scriptName)

    submitToQueue = 'qsub -V -cwd -l h_vmem=6G -q all.q -N '+jobName+' '+scriptName 
    print submitToQueue
    os.system(submitToQueue)
    
    print "\n@@@@@ END JOB @@@@@@@@@@@@@@@"
        
    #os.system('rm testME_tmp.py')
  
###########################################
###########################################

masses = cms.vdouble(120)

#for i in range(40):
#    #print "(%d,%d)" % (i*50+1, (i+1)*50)
#    submitTest('rec_openloops_ttjets_p'+str(i+1), 0, masses,   1,1,1,1,1,    1,1,1,   0.15,0.18, 20,    i*100+1, (i+1)*100)

#submitTest('rec_default_largerLR_range', 0, masses,   1,1,1,1,1,    1,1,1,   0.15,0.18, 20,    0,1000)
#submitTest('rec_METup', 0, masses,     1,1,1,1,1,    1,1,1,   0.15,0.18, 40,    0,1000)
#submitTest('rec_ScaleBup', 0, masses,  1,1,1,1,1,    1,1,1,   0.15,0.27, 40,    0,1000)
#submitTest('rec_ScaleLup', 0, masses,  1,1,1,1,1,    1,1,1,   0.22,0.18, 40,    0,1000)
#submitTest('rec_noJac'   ,   1,0,1,1,1,    1,1,1,   0.15,0.18, 20,    0,1000)

#submitTest('rec_default_noH', 0, masses,   1,1,1,1,1,    1,1,1,   0.15,0.18, 20,    0,1000)

#submitTest('gen_default_largerLR_range', 0, masses,   1,1,1,1,1,    0,1,1,  0.15,0.18, 20,      0,1000)
#submitTest('gen_noME'   ,   0,1,1,1,1,    0,1,1,  0.15,0.18, 20,      0,1000)
#submitTest('gen_noJac'  ,   1,0,1,1,1,    0,1,1,  0.15,0.18, 20,      0,1000)
#submitTest('gen_noMET', 0, masses,     1,1,0,1,1,    0,1,1,  0.15,0.18, 20,      0,1000)
#submitTest('gen_noTF', 0, masses,      1,1,1,0,1,    0,1,1,  0.15,0.18, 20,      0,1000)
#submitTest('gen_noPDF', 0, masses,     1,1,1,1,0,    0,1,1,  0.15,0.18, 20,      0,1000)

#submitTest('gen_openloops_SL1wj_p1', 0, masses,     1,1,1,1,1,   0,1,1,  0.15,0.18, 20,  1,20)
#submitTest('gen_openloops_SL1wj_p2', 0, masses,     1,1,1,1,1,   0,1,1,  0.15,0.18, 20, 21,40)
#submitTest('gen_openloops_SL1wj_p3', 0, masses,     1,1,1,1,1,   0,1,1,  0.15,0.18, 20, 41,60)
#submitTest('gen_openloops_SL1wj_p4', 0, masses,     1,1,1,1,1,   0,1,1,  0.15,0.18, 20, 61,80)
#submitTest('gen_openloops_SL1wj_p5', 0, masses,     1,1,1,1,1,   0,1,1,  0.15,0.18, 20, 81,100)

#submitTest('rec_openloops_SL1wj_p1', 0, masses,     1,1,1,1,1,   1,1,1,  0.15,0.18, 20,  0,20)
#submitTest('rec_openloops_SL1wj_p2', 0, masses,     1,1,1,1,1,   1,1,1,  0.15,0.18, 20, 21,40)
#submitTest('rec_openloops_SL1wj_p3', 0, masses,     1,1,1,1,1,   1,1,1,  0.15,0.18, 20, 41,60)
#submitTest('rec_openloops_SL1wj_p4', 0, masses,     1,1,1,1,1,   1,1,1,  0.15,0.18, 20, 61,80)
#submitTest('rec_openloops_SL1wj_p5', 0, masses,     1,1,1,1,1,   1,1,1,  0.15,0.18, 20, 81,100)


#submitTest('rec_noME_highstat'   ,0, 0,1,1,1,1,1,1,1,201,1000)
#submitTest('rec_noPDF_highstat'  ,0, 1,1,1,1,0,1,1,1,201,1000)


counter = 0;
#for i in [60,70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210]:
for i in [60,65,70,75,80,85,90,95,100,105,110,115,120, 125,130, 135,140, 145,150,155, 160, 165,170,175, 180,185, 190, 195,200, 205,210]:
    counter = counter + 1
    #print "(%d,%d)" % (i*50+1, (i+1)*50)
    mymasses = cms.vdouble(i)
    #print "(%d,%d)" % (i, i+5)
    #print 'rec_XSec_p'+str(counter)
    #submitTest('rec_XSec_p'+str(counter), 2, mymasses,   1,1,1,1,1,    1,1,1,   0.15,0.18, 20,  1, 1)
    submitTest('rec_Acc_p'+str(counter),  3, mymasses,   1,1,1,1,1,    1,1,1,   0.15,0.18, 20,  1, 1)
