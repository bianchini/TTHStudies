#!/usr/bin/env python

import commands
import re
import os

import FWCore.ParameterSet.Config as cms

import sys
sys.path.append('./')


def submit( sample = 'default' ):

    indir  = '/hdfs/local/bianchi/Sherpa/'+sample

    outdir = '../root/files/Sherpa/'

    output   = commands.getoutput("ls -1 "+indir)
    outFiles = re.split(r'\n',output)

    files_per_job = 5

    l         = []
    counter   = 0
    totcounter= 0
    job       = 0
    
    print "\nProcessing files in %s" % indir
    print "Total number of files %s" % len(outFiles)
    
    for file in outFiles:
        
        if re.search(".hepmc2g", file)==None:
            print file, " is not a valid HepMC file"
            continue

        l.append( file )

        counter    += 1
        totcounter += 1
        
        if counter==files_per_job or totcounter==len(outFiles):

            job += 1

            print "\n@@@@@ START JOB @@@@@@@@@@@@@@@"
        
            print "Overload hepMCConverter.py ..."
            os.system('cp ../python/hepMCConverter.py ./')

            print "Creating the shell file for the batch..."
            scriptName = 'job_hepMC_'+sample+'_p'+str(job)+'.sh'
            jobName    = 'job_hepMC_'+sample+'_p'+str(job)

            from hepMCConverter import process

            l_in = []
            for file_i in l:
                l_in.append( indir+'/'+file_i)

            process.fwliteInput.pathToFile  = cms.vstring( l_in )
            process.fwliteInput.outFileName = cms.string(outdir+"hepMC_"+sample+"_p"+str(job)+".root") 

            out = open(jobName+'.py','w')
            out.write(process.dumpPython())     

            f = open(scriptName,'w')
            f.write('#!/bin/bash\n\n')
            f.write('cd ${CMSSW_BASE}/src/Bianchi/TTHStudies/bin/\n')
            f.write('export SCRAM_ARCH="slc5_amd64_gcc462"\n')
            f.write('source $VO_CMS_SW_DIR/cmsset_default.sh\n')
            f.write('eval `scramv1 runtime -sh`\n')
            f.write('\n\n')
            f.write('HepMCtoStep2Converter ./'+jobName+'.py\n')

            f.close()
            os.system('chmod +x '+scriptName)
            
            submitToQueue = 'qsub -q main -N job_hepMC_p'+str(job)+' '+scriptName 
            print submitToQueue
            os.system(submitToQueue)
    
            print "@@@@@ END JOB @@@@@@@@@@@@@@@"

            counter = 0
            l[:] = []

        
            print "\n==> Submitted %s jobs over a total of %s files" % (job,totcounter)

#######################################################
#######################################################


samples = [ #'default',
            #'defaultX05',
            #'defaultX2',
            #'glo_soft',
            #'R_Mbb',
            #'Q_CMMPS',
            #'NNPDF',
            #'MSTW',
            #'MPI_up',
            #'MPI_down',
            #'CSS_KIN',
            'LO_full'
            ]

for sample in samples:
    submit( sample )
