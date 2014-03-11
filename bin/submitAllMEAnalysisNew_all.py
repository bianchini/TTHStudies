#!/usr/bin/env python


import commands
import re
import os

import sys
sys.path.append('./')

import FWCore.ParameterSet.Config as cms

########### mass scan ###########
massesH     = cms.vdouble(125)
#massesH     = cms.vdouble(45., 55., 65., 75., 85., 95., 105., 115., 125., 135., 145., 155., 165., 185., 205., 225., 250., 275., 300.)
massesT     = cms.vdouble(174.3)
#massesT     = cms.vdouble(115., 125., 135., 145., 155., 165., 175., 185., 195., 205.,215., 225.)
#################################

# cetral mass values
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

# speed up the job not doing VEGAS integration
speedup     = 0

# cut values to select events
btag_prob_cut_6jets = 0.96675 # <--- 0.988
btag_prob_cut_5jets = 0.98225 # <--- 0.992
btag_prob_cut_4jets = 0.95295 # <--- 0.85 <--- 0.95295 <--- 0.992

# regression
useRegression = 0

# use gen-jets or reco-jets ???
doGenLevelAnalysis = 0

# smear jets by TF_smear
smearJets = 0

# btag-thresholds
csv_WP_L = 0.244
csv_WP_M = 0.679 
csv_WP_T = 0.898

# select by btag_LR
selectByBTagShape  = 1

# recover the <4 btag bin
recoverTopBTagBin  = 1

# test SLw1j hypothesis on type3 events
testSLw1jType3     = 1

# test hypo SLw1j on up to nMaxJetsSLw1jType3 untagged jets
nMaxJetsSLw1jType3 = 4

# an extra name for the output files
extraoutname = ""
if len(massesH)>1:
    extraoutname = "MHscan_"
if len(massesT)>1:
    extraoutname = "MTscan_"

# set integral optimizations
integralOption0 = 0 # chi2-optimization
integralOption1 = 0 # permutation pruning
integralOption2 = 1 # integration speed-up

# ntuplize all events ?
ntuplizeAll = 0

# systematics
systematics = cms.vint32(0,1,2,3,4,5,6)



def submitMEAnalysisNew_all(script,
                            sample,
                            evLow,evHigh):

    print "Overload meAnalysisNew_all.py..."
    os.system('cp ../python/meAnalysisNew_all.py ./')

    from meAnalysisNew_all import process
    
    print "Creating the shell file for the batch..."
    scriptName = 'job_'+script+'.sh'
    jobName    = 'job_'+script

    for sam in process.fwliteInput.samples:
        if sam.nickName != sample:
            sam.skip = cms.bool(True)
        else:
            sam.skip = cms.bool(False)
            
    process.fwliteInput.outFileName      = cms.string('../root/MEAnalysisNew_'+extraoutname+script+'.root')
    process.fwliteInput.pathToFile       = cms.string('dcap://t3se01.psi.ch:22125//pnfs/psi.ch/cms/trivcat/store//user/bianchi/HBB_EDMNtuple/AllHDiJetPt_V2/')

    if useRegression:
        process.fwliteInput.pathToTF         = cms.string('./root/transferFunctionsTEST_reg.root')
        process.fwliteInput.pathToCP         = cms.string('./root/ControlPlotsTEST_reg.root')
        process.fwliteInput.pathToCP_smear   = cms.string("./root/ControlPlotsTEST_reg_gen.root")
    else:
        process.fwliteInput.pathToTF         = cms.string('./root/transferFunctionsTEST.root')
        process.fwliteInput.pathToCP         = cms.string('./root/ControlPlotsTEST.root')
        process.fwliteInput.pathToCP_smear   = cms.string("./root/ControlPlotsTEST_std_gen.root")

    process.fwliteInput.norm             = cms.untracked.int32(norm)
    process.fwliteInput.useME            = cms.untracked.int32(useME)
    process.fwliteInput.useJac           = cms.untracked.int32(useJac)
    process.fwliteInput.useMET           = cms.untracked.int32(useMET)
    process.fwliteInput.useTF            = cms.untracked.int32(useTF)
    process.fwliteInput.usePDF           = cms.untracked.int32(usePDF)

    process.fwliteInput.integralOption0  = cms.untracked.int32( integralOption0 )
    process.fwliteInput.integralOption1  = cms.untracked.int32( integralOption1 )
    process.fwliteInput.integralOption2  = cms.untracked.int32( integralOption2 )

    process.fwliteInput.selectByBTagShape  = cms.untracked.int32(selectByBTagShape)
    process.fwliteInput.recoverTopBTagBin  = cms.untracked.int32(recoverTopBTagBin)

    process.fwliteInput.testSLw1jType3     = cms.untracked.int32(testSLw1jType3)
    process.fwliteInput.nMaxJetsSLw1jType3 = cms.untracked.int32(nMaxJetsSLw1jType3)

    process.fwliteInput.doType0            = cms.untracked.int32(not selectByBTagShape)
    process.fwliteInput.doType1            = cms.untracked.int32(not selectByBTagShape)
    process.fwliteInput.doType2            = cms.untracked.int32(not selectByBTagShape)
    process.fwliteInput.doType3            = cms.untracked.int32(not selectByBTagShape)
    process.fwliteInput.doType6            = cms.untracked.int32(not selectByBTagShape)
    process.fwliteInput.doType7            = cms.untracked.int32(0)

    #process.fwliteInput.doType0ByBTagShape = cms.untracked.int32(    selectByBTagShape)
    #process.fwliteInput.doType1ByBTagShape = cms.untracked.int32(    selectByBTagShape)
    #process.fwliteInput.doType2ByBTagShape = cms.untracked.int32(    selectByBTagShape)        
    #process.fwliteInput.doType3ByBTagShape = cms.untracked.int32(    selectByBTagShape)
    #process.fwliteInput.doType6ByBTagShape = cms.untracked.int32(    selectByBTagShape)
    process.fwliteInput.doType0ByBTagShape = cms.untracked.int32(    selectByBTagShape)
    process.fwliteInput.doType1ByBTagShape = cms.untracked.int32(    selectByBTagShape)
    process.fwliteInput.doType2ByBTagShape = cms.untracked.int32(    0)        
    process.fwliteInput.doType3ByBTagShape = cms.untracked.int32(    selectByBTagShape)
    process.fwliteInput.doType6ByBTagShape = cms.untracked.int32(    0)
  
    process.fwliteInput.btag_prob_cut_6jets = cms.untracked.double(btag_prob_cut_6jets)
    process.fwliteInput.btag_prob_cut_5jets = cms.untracked.double(btag_prob_cut_5jets)
    process.fwliteInput.btag_prob_cut_4jets = cms.untracked.double(btag_prob_cut_4jets)

    process.fwliteInput.csv_WP_L            =  cms.untracked.double(csv_WP_L)
    process.fwliteInput.csv_WP_M            =  cms.untracked.double(csv_WP_M)
    process.fwliteInput.csv_WP_T            =  cms.untracked.double(csv_WP_T)

    process.fwliteInput.useRegression       = cms.untracked.int32(useRegression)

    process.fwliteInput.massesH             = massesH
    process.fwliteInput.massesT             = massesT
    process.fwliteInput.MH                  = cms.untracked.double(MH)
    process.fwliteInput.MT                  = cms.untracked.double(MT)
    process.fwliteInput.fixNumEvJob         = cms.untracked.int32(fixNumEvJob)
    process.fwliteInput.evLimits            = cms.vint32(evLow,evHigh)

    process.fwliteInput.printout            = cms.untracked.int32(printout)

    process.fwliteInput.doGenLevelAnalysis  = cms.untracked.int32(doGenLevelAnalysis)
    process.fwliteInput.smearJets           = cms.untracked.int32(smearJets)

    process.fwliteInput.speedup             = cms.untracked.int32(speedup)
    process.fwliteInput.ntuplizeAll         = cms.untracked.int32(ntuplizeAll)

    process.fwliteInput.systematics         = systematics
 
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
    f.write('MEAnalysisNew_all ./'+jobName+'.py\n')
    f.close()
    os.system('chmod +x '+scriptName)

    submitToQueue = 'qsub -V -cwd -l h_vmem=2G -q all.q -N job'+sample+' '+scriptName 
    print submitToQueue
    os.system(submitToQueue)
    
    print "\n@@@@@ END JOB @@@@@@@@@@@@@@@"
        
    #os.system('rm testME_tmp.py')
  
###########################################
###########################################



def submitFullMEAnalysisNew_all( analysis ):

    print "Running full analysis for ",  analysis

    ###################################################### TTH
    ######################################################

    # TTH125 --> 194808
    sample  = 'TTH125'
    counter = 0
    num_of_jobs =   1
    evs_per_job =   1
    
    num_of_jobs = 270#180
    evs_per_job = 725#1085 

    for i in range(num_of_jobs):
        counter = counter + 1
        submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample,  i*evs_per_job+1, (i+1)*evs_per_job )

    #return

    
    ###################################################### TTJets SL
    ######################################################
        
    # TTJetsSemiLept --> 16749255
    sample  = 'TTJetsSemiLept'
    counter = 0
    num_of_jobs =   1
    evs_per_job =   1 

    num_of_jobs =       150
    evs_per_job =    112000  

    for i in range(num_of_jobs):
        counter = counter + 1
        submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample, i*evs_per_job+1, (i+1)*evs_per_job )

    #return

    ###################################################### TTJets FL
    ######################################################

    # TTJetsFullLept --> 8932897
    sample  = 'TTJetsFullLept'
    counter = 0
    num_of_jobs =   1
    evs_per_job =   1

    num_of_jobs =     46
    evs_per_job = 202000  

    for i in range(num_of_jobs):
        counter = counter + 1
        submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample,  i*evs_per_job+1, (i+1)*evs_per_job )

 
    if re.search("all",   analysis )==None and re.search("nominal",   analysis )==None:
        return


    ###################################################### EWK
    ######################################################

    # DYJets10to50
    sample  = 'DYJets10to50'
    counter = 0
    for i in range(1):
        counter = counter + 1
        submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample, 0, -1 )

    # DYJets50
    sample  = 'DYJets50'
    counter = 0
    for i in range(1):
        counter = counter + 1
        submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample, 0, -1 )

    # WJets 
    sample  = 'WJets'
    counter = 0
    for i in range(1):
        counter = counter + 1
        submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample, 0, -1 )

    ###################################################### Single-top
    ######################################################

    # TtW
    sample  = 'TtW'
    counter = 0
    for i in range(1):
        counter = counter + 1
        submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample, 0, -1 )
        
   # Tt
    sample  = 'Tt'
    counter = 0
    for i in range(1):
        counter = counter + 1
        submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample, 0, -1 )

   # Ts
    sample  = 'Ts'
    counter = 0
    for i in range(1):
        counter = counter + 1
        submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample, 0, -1 )

    # TbartW
    sample  = 'TbartW'
    counter = 0
    for i in range(1):
        counter = counter + 1
        submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample, 0, -1 )

    # Tbart
    sample  = 'Tbart'
    counter = 0
    for i in range(1):
        counter = counter + 1
        submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample, 0, -1 )

    # Tbars
    sample  = 'Tbars'
    counter = 0
    for i in range(1):
        counter = counter + 1
        submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample, 0, -1 )


    ###################################################### Di-boson
    ######################################################


   # WW
    sample  = 'WW'
    counter = 0
    for i in range(1):
        counter = counter + 1
        submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample, 0, -1 )

   # WZ
    sample  = 'WZ'
    counter = 0
    for i in range(1):
        counter = counter + 1
        submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample, 0, -1 )

   # ZZ
    sample  = 'ZZ'
    counter = 0
    for i in range(1):
        counter = counter + 1
        submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample, 0, -1 )


    ###################################################### TTZ
    ######################################################

    # TTZ   --> 112517
    sample  = 'TTZ'
    counter = 0
    num_of_jobs =    1
    evs_per_job =    1

    num_of_jobs =   15
    evs_per_job = 8000     # ---> ~ 40/job

    for i in range(num_of_jobs):
        counter = counter + 1
        submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample, i*evs_per_job+1, (i+1)*evs_per_job )

    # TTW
    sample  = 'TTW'
    counter = 0
    for i in range(1):
        counter = counter + 1
        submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample, 0, -1 )


    ###################################################### Data
    ######################################################

    datasamples = ['DoubleElectron_Run2012A-13Jul2012-v1_ProcFIXED',
                   'DoubleElectron_Run2012A-recover-06Aug2012-v1_ProcV2',
                   'DoubleElectron_Run2012B-13Jul2012-v1_ProcFIXED',
                   'DoubleElectron_Run2012C-PromptReco-v2_HBB_EDMNtupleV42_ProcV1',
                   'DoubleElectron_Run2012C-PromptReco-v2_HBB_EDMNtupleV42_ProcV2',
                   'DoubleElectronRun2012C-EcalRecover_11Dec2012-v1_v2',
                   'DoubleElectronRun2012CAug24RerecoEdmV42',
                   'DoubleElectronRun2012D',

                   'SingleElectronRun2012AAug06EdmV42',
                   'SingleElectronRun2012AJul13EdmV42b',
                   'SingleElectronRun2012BJul13EdmV42',
                   'SingleElectronRun2012C-EcalRecover_11Dec2012-v1_v2',
                   'SingleElectronRun2012CAug24RerecoEdmV42',
                   'SingleElectronRun2012CPromptv2EdmV42',
                   'SingleElectronRun2012CPromptV2TopUpEdmV42',
                   'SingleElectronRun2012D-PromptReco-v1_v3',

                   'SingleMuRun2012AAug06EdmV42',
                   'SingleMuRun2012AJul13EdmV42',
                   'SingleMuRun2012BJul13EdmV42',
                   'SingleMuRun2012C-EcalRecover_11Dec2012-v1_v2',
                   'SingleMuRun2012CAug24RerecoEdmV42',
                   'SingleMuRun2012CPromptv2EdmV42',
                   'SingleMuRun2012CPromptV2TopUpEdmV42',
                   'SingleMuRun2012D-PromptReco-v1'
                   ]

    for datasample in datasamples:

        continue

        sample  = 'Run2012_'+datasample

        counter     =    0
        num_of_jobs =    1
        evs_per_job =    1
                      
        if (re.search("Run2012AAug06EdmV42",datasample)!=None):  # 148139
            for i in range(1):
                counter = counter + 1
                submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample, 0, -1 )

        elif (re.search("Run2012AJul13EdmV42b", datasample)!=None):  # 1551019
            for i in range(1):
                counter = counter + 1
                submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample, 0, -1 )

        elif (re.search("Run2012BJul13EdmV42", datasample)!=None):  # 9351330
            num_of_jobs =        4
            evs_per_job =  2340000
            for i in range(num_of_jobs):
                counter = counter + 1
                submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample,  i*evs_per_job+1, (i+1)*evs_per_job  )

        elif (re.search("Run2012C-EcalRecover_11Dec2012-v1_v2", datasample)!=None):  # 263593
            for i in range(1):
                counter = counter + 1
                submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample, 0, -1 )

        elif (re.search("Run2012CAug24RerecoEdmV42", datasample)!=None):  # 1064158
            for i in range(1):
                counter = counter + 1
                submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample, 0, -1 )

        elif (re.search("Run2012CPromptv2EdmV42", datasample)!=None):  # 9768094
            num_of_jobs =        4
            evs_per_job =  2443000
            for i in range(num_of_jobs):
                counter = counter + 1
                submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample,  i*evs_per_job+1, (i+1)*evs_per_job  )

        elif (re.search("Run2012CPromptV2TopUpEdmV42", datasample)!=None):  # 3491407
            for i in range(1):
                counter = counter + 1
                submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample, 0, -1 )

        elif (re.search("Run2012D-PromptReco-v1", datasample)!=None):  # 16178887
            num_of_jobs =        6*6
            evs_per_job =  2800000/6
            for i in range(num_of_jobs):
                counter = counter + 1
                submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample,  i*evs_per_job+1, (i+1)*evs_per_job  )
        else:
            for i in range(1):
                counter = counter + 1
                submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample, 0, -1 )

###########################################
###########################################


analyses = ['all_testNewType3']

for analysis in analyses:
    if doGenLevelAnalysis:
        analysis = analysis+'_gen'
    else:
        analysis = analysis+'_rec'   
    if useRegression:
        analysis = analysis+'_reg'
    else:
        analysis = analysis+'_std'
    submitFullMEAnalysisNew_all( analysis)
