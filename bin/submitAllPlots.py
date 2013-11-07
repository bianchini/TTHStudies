#!/usr/bin/env python


import commands
import re
import os

import sys
sys.path.append('./')



def submitPlots(VType):

    scriptName = 'batchPlotZmm_'+VType+'.sh'
    jobName    = 'batchPlotZmm_'+VType

    f = open(scriptName,'w')
    f.write('#!/bin/bash\n\n')
    f.write('cd /shome/bianchi/CMSSW_5_3_3_patch2/src/Bianchi/TTHStudies/bin\n')
    f.write('source /swshare/psit3/etc/profile.d/cms_ui_env.sh\n')
    f.write('export SCRAM_ARCH="slc5_amd64_gcc462"\n')
    f.write('source $VO_CMS_SW_DIR/cmsset_default.sh\n')
    f.write('eval `scramv1 runtime -sh`\n')
    f.write('export LD_PRELOAD="libglobus_gssapi_gsi_gcc64pthr.so.0":${LD_PRELOAD}\n')
    f.write('export LD_LIBRARY_PATH=/swshare/glite/globus/lib/:/swshare/glite/d-cache/dcap/lib64/:$LD_LIBRARY_PATH\n')
    f.write('\n\n')
    f.write('\n\n')
    f.write("PlotZmm ../python/plot"+VType+'.py')
    f.close()

    os.system('chmod +x '+scriptName)
    submitToQueue = 'qsub -V -cwd -l h_vmem=6G -q all.q -N '+jobName+' '+scriptName 
    print submitToQueue
    os.system(submitToQueue)


#################################################
#################################################

submitPlots("VType0")
submitPlots("VType1")
submitPlots("VType2")
submitPlots("VType3")
