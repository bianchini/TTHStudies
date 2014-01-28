#!/usr/bin/env python


import commands
import re
import os

import FWCore.ParameterSet.Config as cms

import sys
sys.path.append('./')


indir  = '/pnfs/psi.ch/cms/trivcat/store/user/bianchi/HepMC/SherpaOpenLoops/Jan24_2014/default_1/'
outdir = '/scratch/bianchi/HBB_EDMNtuple/SherpaOpenLoops/'

output   = commands.getoutput("ls -1 "+indir)
outFiles = re.split(r'\n',output)

files_per_job = 10

l         = []
counter   = 0
totcounter= 0
job       = 0

print "\nProcessing files in %s" % indir
print "Total number of files %s" % len(outFiles)

for file in outFiles:

    #if job>0:
    #    continue
    
    l.append( file )

    counter    += 1
    totcounter += 1
    
    if counter==files_per_job or counter==len(outFiles):

        job += 1

        print "\n@@@@@ START JOB @@@@@@@@@@@@@@@"
        
        print "Overload hepMCConverter.py ..."
        os.system('cp ../python/hepMCConverter.py ./')

        print "Creating the shell file for the batch..."
        scriptName = 'job_hepMC_default_1_p'+str(job)+'.sh'
        jobName    = 'job_hepMC_default_1_p'+str(job)

        from hepMCConverter import process

        l_in = []
        for file_i in l:
            l_in.append( '/scratch/bianchi/'+file_i)

        process.fwliteInput.pathToFile  = cms.vstring( l_in )
        process.fwliteInput.outFileName = cms.string(outdir+"hepMC_default_1_p"+str(job)+".root") 

        out = open(jobName+'.py','w')
        out.write(process.dumpPython())     

        f = open(scriptName,'w')
        f.write('#!/bin/bash\n\n')
        f.write('cd ${CMSSW_BASE}/src/Bianchi/TTHStudies/bin/\n')
        f.write('source /swshare/psit3/etc/profile.d/cms_ui_env.sh\n')
        f.write('export SCRAM_ARCH="slc5_amd64_gcc462"\n')
        f.write('source $VO_CMS_SW_DIR/cmsset_default.sh\n')
        f.write('eval `scramv1 runtime -sh`\n')
        f.write('\n\n')
        f.write('if !(ls /scratch/bianchi/ &> /dev/null) ; then\n')
        f.write('    mkdir /scratch/bianchi/\n')
        f.write('fi\n')
        f.write('ls -ltr /scratch/bianchi/\n')
        for file_i in l:
            f.write('srmcp -2 srm://t3se01.psi.ch:8443/srm/managerv2?SFN='+indir+file_i+' file:////scratch/bianchi/'+file_i+'\n')
        f.write('ls -ltr /scratch/bianchi/\n')
        f.write('HepMCtoStep2Converter ./'+jobName+'.py\n')
        f.write('echo "Remove input files..."\n')   
        for file_i in l:
            f.write('rm /scratch/bianchi/'+file_i+'\n')
        f.write('echo "Done:"\n')
        f.write('ls -ltr /scratch/bianchi/\n')
        f.close()
        os.system('chmod +x '+scriptName)

        submitToQueue = 'qsub -V -cwd -l h_vmem=1G -q all.q -N job_hepMC_p'+str(job)+' '+scriptName 
        print submitToQueue
        os.system(submitToQueue)
    
        print "@@@@@ END JOB @@@@@@@@@@@@@@@"

        counter = 0
        l[:] = []

        
print "\n==> Submitted %s jobs over a total of %s files" % (job,totcounter)
