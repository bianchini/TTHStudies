#!/usr/bin/env python


import commands
import re
import os
import string
import FWCore.ParameterSet.Config as cms

import sys
sys.path.append('./')

from ntuple  import process 

pushToT3 = True
debug    = False

outdir = '/scratch/bianchi/'
indir  = '/pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_V3_tmp'

#outdir = './'
#indir  = '/pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_V3'

def processAllBatch(jobName, isPisa, outName, split):

    process.fwliteInput.fileNames = ()

    if string.find( jobName, 'Run2012' )>0:
        process.Analyzer.isMC  = cms.bool(False)
    
    input = open('fileList_'+jobName+'.txt','r')
    counter     = 1
    counterJobs = 0
    for line in input:
        foo = line.strip('\n')
        if counter>=split[0] and counter<=split[1]:
            if isPisa:
                process.fwliteInput.fileNames.append('root://xrootd.ba.infn.it//'+foo)
            else:
                process.fwliteInput.fileNames.append('dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/'+foo)
            counterJobs += 1
        counter += 1

    input.close()

    outFileName = outName+'_'+str(split[0])+'-'+str(split[1])+'.root'
    
    process.fwliteOutput.fileName = cms.string(outdir+'/'+outFileName)

    out = open('myStep2_'+jobName+'_'+str(split[0])+'-'+str(split[1])+'.py','w')
    out.write(process.dumpPython())
    out.close()

    scriptName = 'myStep2_'+jobName+'_'+str(split[0])+'-'+str(split[1])+'.sh'
    
    f = open(scriptName,'w')
    f.write('#!/bin/bash\n\n')
    f.write('cd /shome/bianchi/CMSSW_5_3_3_patch2_New/src/VHbbAnalysis/VHbbDataFormats/bin\n')
    f.write('export SCRAM_ARCH="slc5_amd64_gcc462"\n')
    f.write('source $VO_CMS_SW_DIR/cmsset_default.sh\n')
    f.write('eval `scramv1 runtime -sh`\n')
    f.write('\n')
    f.write('\n')
    if pushToT3:
        f.write('if !(ls /scratch/bianchi/ &> /dev/null) ; then\n')
        f.write('    mkdir /scratch/bianchi/\n')
        f.write('fi\n')
        f.write('ls -ltr /scratch/bianchi/\n\n')

    mainexec = 'Ntupler myStep2_'+jobName+'_'+str(split[0])+'-'+str(split[1])+'.py'
    print mainexec
    f.write(mainexec+'\n\n')

    if pushToT3:
        push = 'lcg-cp -b -D srmv2 '+outdir+'/'+outFileName+' srm://t3se01.psi.ch:8443/srm/managerv2?SFN='+indir+'/'+outFileName
        print push
        f.write(push+'\n\n')
        remove = 'rm '+outdir+'/'+outFileName
        f.write(remove+'\n')

    f.close()

    os.system('chmod +x '+scriptName)
    submitToQueue = 'qsub -V -cwd -l h_vmem=2G -q all.q -N my'+jobName+'_'+str(split[0])+'-'+str(split[1])+' '+scriptName 

    if not debug:
        os.system(submitToQueue)
        print submitToQueue


    return float(counterJobs)/float(counter-1)

###########################################################################
###########################################################################


#processAllBatch("WJets", "DiJetPt_WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball", [1, 500])
#processAllBatch("TTH110", "DiJetPt_TTH_HToBB_M-110_8TeV-pythia6", [1, 500])
#processAllBatch("TTH115", "DiJetPt_TTH_HToBB_M-115_8TeV-pythia6", [1, 500])
#processAllBatch("TTH120", "DiJetPt_TTH_HToBB_M-120_8TeV-pythia6", [1, 500])
#processAllBatch("TTH130", "DiJetPt_TTH_HToBB_M-130_8TeV-pythia6", [1, 500])
#processAllBatch("TTH135", "DiJetPt_TTH_HToBB_M-135_8TeV-pythia6", [1, 500])

total = 0.
for k in range(8):
    #print "\nProcessing job....", k
    #total += processAllBatch("T_t-channel", 1, "DiJetPt_T_t-channel_TuneZ2star_8TeV-powheg-tauola", [k*5  +1,(k+1)*5 ])
print '\n**************************************\nFraction of processed sample: %s\n**************************************' % total

total = 0.
for k in range(70):
    #print "Processing....", k
    #total += processAllBatch("DYJets10to50", 0, "DiJetPt_DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph", [k*40+1,(k+1)*40])
print '\n**************************************\nFraction of processed sample: %s\n**************************************' % total

total = 0.
for k in range(30):
    #print "Processing....", k
    #total += processAllBatch("DYJets50", 1, "DiJetPt_DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph", [k*10+1,(k+1)*10])
print '\n**************************************\nFraction of processed sample: %s\n**************************************' % total


total = 0.
for k in range(15):
    #print "Processing....", k
    #total += processAllBatch("WJetsToLNu", 1, "DiJetPt_WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball", [k*10+1,(k+1)*10])
print '\n**************************************\nFraction of processed sample: %s\n**************************************' % total


total = 0.
for k in range(2):
    #print "Processing....", k
    #total += processAllBatch("T_s-channel", 1, "DiJetPt_T_s-channel_TuneZ2star_8TeV-powheg-tauola", [k*10  +1,(k+1)*10 ])
print '\n**************************************\nFraction of processed sample: %s\n**************************************' % total


total = 0.
for k in range(4):
    #print "Processing....", k
    #total += processAllBatch("T_t-channel", 1, "DiJetPt_T_t-channel_TuneZ2star_8TeV-powheg-tauola", [k*10  +1,(k+1)*10 ])
print '\n**************************************\nFraction of processed sample: %s\n**************************************' % total


total = 0.
for k in range(3):
    #print "Processing....", k
    #total += processAllBatch("T_tW-channel-DR", 1, "DiJetPt_T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola", [k*10  +1,(k+1)*10 ])
print '\n**************************************\nFraction of processed sample: %s\n**************************************' % total


total = 0.
for k in range(2):
    #print "Processing....", k
    #total += processAllBatch("Tbar_s-channel", 1, "DiJetPt_Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola", [k*10  +1,(k+1)*10 ])
print '\n**************************************\nFraction of processed sample: %s\n**************************************' % total


total = 0.
for k in range(4):
    #print "Processing....", k
    #total += processAllBatch("Tbar_t-channel", 1, "DiJetPt_Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola", [k*10  +1,(k+1)*10 ])
print '\n**************************************\nFraction of processed sample: %s\n**************************************' % total


total = 0.
for k in range(3):
    #print "Processing....", k
    #total += processAllBatch("Tbar_tW-channel-DR", 1, "DiJetPt_Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola", [k*10  +1,(k+1)*10 ])
print '\n**************************************\nFraction of processed sample: %s\n**************************************' % total



total = 0.
for k in range(14*5):
    print "Processing....", k
    total += processAllBatch("WW", 1, "DiJetPt_WW_TuneZ2star_8TeV_pythia6_tauola", [k*10  +1,(k+1)*10 ])
print '\n**************************************\nFraction of processed sample: %s\n**************************************' % total


total = 0.
for k in range(14*5):
    print "Processing....", k
    total += processAllBatch("WZ", 1, "DiJetPt_WZ_TuneZ2star_8TeV_pythia6_tauola", [k*10  +1,(k+1)*10 ])
print '\n**************************************\nFraction of processed sample: %s\n**************************************' % total


total = 0.
for k in range(14*5):
    print "Processing....", k
    total += processAllBatch("ZZ", 1, "DiJetPt_ZZ_TuneZ2star_8TeV_pythia6_tauola", [k*10  +1,(k+1)*10 ])
print '\n**************************************\nFraction of processed sample: %s\n**************************************' % total


total = 0.
for k in range(34*5):
    print "Processing....", k
    total += processAllBatch("TTJets_SemiLeptMGDecays", 1, "DiJetPt_TTJets_SemiLeptMGDecays_8TeV-madgraph", [k*10  +1,(k+1)*10 ])
print '\n**************************************\nFraction of processed sample: %s\n**************************************' % total


total = 0.
for k in range(16*5):
    print "Processing....", k
    total += processAllBatch("TTJets_FullLeptMGDecays", 1, "DiJetPt_TTJets_FullLeptMGDecays_8TeV-madgraph", [k*10 +1,(k+1)*10 ])
print '\n**************************************\nFraction of processed sample: %s\n**************************************' % total


total = 0.
for k in range(20*5):
    print "Processing....", k
    total += processAllBatch("TTJets_HadronicMGDecays",  1, "DiJetPt_TTJets_HadronicMGDecays_8TeV-madgraph", [k*10  +1,(k+1)*10 ])
print '\n**************************************\nFraction of processed sample: %s\n**************************************' % total


total = 0.
for k in range(8):
    print "Processing....", k
    total += processAllBatch("TTJets_MassiveBinDECAY", 1, "DiJetPt_TTJets_MassiveBinDECAY_8TeV-madgraph", [k*5  +1,(k+1)*5 ])
print '\n**************************************\nFraction of processed sample: %s\n**************************************' % total


#for k in range(4):
#    print "Processing....", k
#    processAllBatch("SingleMuRun2012AJul13", 1, "DiJetPt_SingleMuRun2012AJul13", [k*50  +1,(k+1)*50 ])

#for k in range(1):
#    print "Processing....", k
#    processAllBatch("SingleMuRun2012AAug06", 1, "DiJetPt_SingleMuRun2012AAug06", [k*50  +1,(k+1)*50 ])

#for k in range(32):
#    print "Processing....", k
#    processAllBatch("SingleMuRun2012BJul13", 1, "DiJetPt_SingleMuRun2012BJul13", [k*50  +1,(k+1)*50 ])

#for k in range(4):
#    print "Processing....", k
#    processAllBatch("SingleMuRun2012CAug24Rereco", 1, "DiJetPt_SingleMuRun2012CAug24Rereco", [k*50  +1,(k+1)*50 ])

#for k in range(34):
#    print "Processing....", k
#    processAllBatch("SingleMuRun2012CPromptv2", 1, "DiJetPt_SingleMuRun2012CPromptv2", [k*50  +1,(k+1)*50 ])

#for k in range(12):
#    print "Processing....", k
#    processAllBatch("SingleMuRun2012CPromptV2TopUp", 1, "DiJetPt_SingleMuRun2012CPromptV2TopUp", [k*50  +1,(k+1)*50 ])

#for k in range(30):
#    print "Processing....", k
#    processAllBatch("SingleMuRun2012D-PromptReco-v1", 1, "DiJetPt_SingleMuRun2012D-PromptReco-v1", [k*50  +1,(k+1)*50 ])

#for k in range(1):
#    print "Processing....", k
#    processAllBatch("SingleMuRun2012C-EcalRecover_11Dec2012-v1_v2", 1, "DiJetPt_SingleMuRun2012C-EcalRecover_11Dec2012-v1_v2", [k*51  +1,(k+1)*51 ])




#for k in range(14*4):
#    print "Processing....", k
#    processAllBatch("DYJets10to50", 0, "DiJetPt_DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph", [k*50+1,(k+1)*50])

#processAllBatch("TTH125", 0,"DiJetPt_TTH_HToBB_M-125_8TeV-pythia6", [1, 500])
#processAllBatch("TTW",    0,"DiJetPt_TTWJets_8TeV-madgraph",        [1, 500])
#processAllBatch("TTZ",    0,"DiJetPt_TTZJets_8TeV-madgraph",        [1, 500])

