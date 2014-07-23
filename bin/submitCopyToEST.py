#!/usr/bin/env python

import commands
import re
import os

import sys
sys.path.append('./')

samples = [
    ['full', ['1', '3', '11', '51', '101', '201', '301', '401', '501', '601', '701', '801', '901', '1001', '1101', '1201', '1301', '1401', '1501', '1601', '1701', '1801', '2001'] ],
    #['default', ['1', '101', '201'] ],
    #['defaultX05', ['1', '101', '201'] ],
    #['defaultX2',  [ '101', '201'] ],
    #['glo_soft',   ['1', '101', '201'] ],
    #['R_Mbb',      ['1', '101', '201'] ],
    #['Q_CMMPS',    ['1', '101', '201'] ],
    #['NNPDF',      ['1', '101', '201'] ],
    #['MSTW',       ['1', '101', '201'] ],
    #['MPI_up',     ['1', '101', '201'] ],
    #['MPI_down',   ['1', '101', '201'] ],
    #['CSS_KIN',    ['1', '101', '201'] ],
    #['LO_full',  [''] ],
    ]

for sample in samples:

    print "Copy sample ", sample[0], " with " , len( sample[1] ), " elements" 

    for job in range(len( sample[1] )):    

        separator = '_'
        if (sample[1])[job]=='':
            separator = ''
        
        scriptName = 'job_copyHepMC_'+sample[0]+'_p'+str(job)+'.sh'
        jobName    = 'job_copyHepMC_'+sample[0]+'_p'+str(job)
        
        f = open(scriptName,'w')
        f.write('#!/bin/bash\n')
        f.write('cd /hdfs/local/bianchi/Sherpa\n')
        f.write('if !(ls ./'+sample[0]+' &> /dev/null) ; then\n')
        f.write('    mkdir ./'+sample[0]+'\n')
        f.write('fi\n')
        f.write('cd ./'+sample[0]+'\n')
        f.write('wget http://www.physik.uzh.ch/data/ttHsim/ttbb_CTEQ10_inclusive/Samples_'+sample[0]+'/'+sample[0]+separator+(sample[1])[job]+'.tar.gz\n')
        f.write('tar -xvf '+sample[0]+separator+(sample[1])[job]+'.tar.gz\n')
        f.write('rm '+sample[0]+separator+(sample[1])[job]+'.tar.gz\n')

        f.close()
        os.system('chmod +x '+scriptName)

        submitToQueue = 'qsub -q main -N job_copyHepMC_p'+str(job)+' '+scriptName 
        print submitToQueue
        os.system(submitToQueue)
    
        print "@@@@@ END JOB @@@@@@@@@@@@@@@"
