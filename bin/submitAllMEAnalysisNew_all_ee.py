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
printout    = 0

# speed up the job not doing VEGAS integration
speedup     = 1

# cut values to select events
btag_prob_cut_6jets = 0.96675 # <--- 0.988
btag_prob_cut_5jets = 0.98225 # <--- 0.992
btag_prob_cut_4jets = 0.95295 # <--- 0.85 <--- 0.95295 <--- 0.992

# regression
useRegression = 0

#allows additional btags
recoverTopBTagBin = cms.untracked.int32(1)

# use gen-jets or reco-jets ???
doGenLevelAnalysis = 0

# btag-thresholds
csv_WP_L = 0.244
csv_WP_M = 0.679 
csv_WP_T = 0.898

# select bt btag_LR
selectByBTagShape = 0

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
ntuplizeAll = 1

# systematics
#systematics = cms.vint32([0,1,2,3,4,5,6])
systematics = cms.vint32(0)

def submitMEAnalysisNew_all(script,
                            sample,
                            evLow,evHigh, isZmm = False):

    print "Overload meAnalysisNew_all_ee.py..."
    os.system('cp ../python/meAnalysisNew_all_ee.py ./')

    from meAnalysisNew_all_ee import process
    
    print "Creating the shell file for the batch..."
    scriptName = 'job_'+script+'.sh'
    jobName    = 'job_'+script

    for sam in process.fwliteInput.samples:
        if sam.nickName != sample:
            sam.skip = cms.bool(True)
        else:
            sam.skip = cms.bool(False)
            
    process.fwliteInput.outFileName      = cms.string('../root/MEAnalysisNew_'+extraoutname+script+'.root')
    process.fwliteInput.pathToFile       = cms.string('/hdfs/cms/store/user/liis/TTH_Ntuples_v3/')
#    process.fwliteInput.pathToFile       = cms.string('/home/liis/TTH_Ntuples/CMSSW_5_3_3/src/VHbbAnalysis/VHbbDataFormats/bin/Ntuples_withPDF/Ntuples_merged')

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

    process.fwliteInput.selectByBTagShape = cms.untracked.int32(selectByBTagShape)
    
    process.fwliteInput.doType6ByBTagShape = cms.untracked.int32(    selectByBTagShape)
    process.fwliteInput.doType6            = cms.untracked.int32(not selectByBTagShape)
    process.fwliteInput.doType7            = cms.untracked.int32(not selectByBTagShape)

    process.fwliteInput.doType0ByBTagShape = cms.untracked.int32(    selectByBTagShape)

    if len(massesH) == 1 and len(massesT) == 1:
        process.fwliteInput.doType1ByBTagShape = cms.untracked.int32(    selectByBTagShape)
    else:
        process.fwliteInput.doType1ByBTagShape = cms.untracked.int32(0)        

    process.fwliteInput.doType2ByBTagShape = cms.untracked.int32(    selectByBTagShape)
        
    if len(massesH) == 1 and len(massesT) == 1:
        process.fwliteInput.doType3ByBTagShape = cms.untracked.int32(    selectByBTagShape)
    else:
        process.fwliteInput.doType3ByBTagShape = cms.untracked.int32(0)

    process.fwliteInput.doType0            = cms.untracked.int32(not selectByBTagShape)

    if len(massesH) == 1 and len(massesT) == 1:
        process.fwliteInput.doType1            = cms.untracked.int32(not selectByBTagShape)

    else:
        process.fwliteInput.doType1         = cms.untracked.int32(not selectByBTagShape)

    process.fwliteInput.doType2             = cms.untracked.int32(not selectByBTagShape)

    if len(massesH) == 1 and len(massesT) == 1:
        process.fwliteInput.doType3         = cms.untracked.int32(not selectByBTagShape)
    else:
        process.fwliteInput.doType3         = cms.untracked.int32(not selectByBTagShape)
    
    process.fwliteInput.btag_prob_cut_6jets = cms.untracked.double(btag_prob_cut_6jets)
    process.fwliteInput.btag_prob_cut_5jets = cms.untracked.double(btag_prob_cut_5jets)
    process.fwliteInput.btag_prob_cut_4jets = cms.untracked.double(btag_prob_cut_4jets)

    process.fwliteInput.csv_WP_L            =  cms.untracked.double(csv_WP_L)
    process.fwliteInput.csv_WP_M            =  cms.untracked.double(csv_WP_M)
    process.fwliteInput.csv_WP_T            =  cms.untracked.double(csv_WP_T)

    process.fwliteInput.useCSVcalibration   = cms.untracked.int32(1)
    process.fwliteInput.useRegression       = cms.untracked.int32(useRegression)
    process.fwliteInput.recoverTopBTagBin   = cms.untracked.int32(1)

    process.fwliteInput.massesH             = massesH
    process.fwliteInput.massesT             = massesT
    process.fwliteInput.MH                  = cms.untracked.double(MH)
    process.fwliteInput.MT                  = cms.untracked.double(MT)
    process.fwliteInput.fixNumEvJob         = cms.untracked.int32(fixNumEvJob)
    process.fwliteInput.evLimits            = cms.vint32(evLow,evHigh)

    process.fwliteInput.printout            = cms.untracked.int32(printout)

    process.fwliteInput.doGenLevelAnalysis  = cms.untracked.int32(doGenLevelAnalysis)
    process.fwliteInput.speedup             = cms.untracked.int32(speedup)
    process.fwliteInput.ntuplizeAll         = cms.untracked.int32(ntuplizeAll)

    process.fwliteInput.systematics         = systematics
 
    out = open(jobName+'.py','w')
    out.write(process.dumpPython())
   
    f = open(scriptName,'w')
    f.write('#!/bin/bash\n\n')
    f.write('cd ${CMSSW_BASE}/src/Bianchi/TTHStudies/bin/\n')
    f.write('source /swshare/psit3/etc/profile.d/cms_ui_env.sh\n')
    f.write('export SCRAM_ARCH="slc5_amd64_gcc462"\n')
    f.write('source $VO_CMS_SW_DIR/cmsset_default.sh\n')
    f.write('eval `scramv1 runtime -sh`\n')
    #f.write('export LD_PRELOAD="libglobus_gssapi_gsi_gcc64pthr.so.0":${LD_PRELOAD}\n')
    f.write('\n\n')
    f.write('\n\n')
    f.write('MEAnalysisNew_all ./'+jobName+'.py\n')
    f.close()
    os.system('chmod +x '+scriptName)

#    submitToQueue = 'qsub -N job'+sample+' '+scriptName #EE
    submitToQueue = submitToQueue = 'qsub -N job'+sample+' '+scriptName #T2_EE
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
    if speedup:
        num_of_jobs =   1
        evs_per_job =   1
    else:
        num_of_jobs =  180
        evs_per_job = 1085 

    for i in range(num_of_jobs):
        counter = counter + 1
        if speedup:
            submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample, 0, -1 )
        else:
            submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample,  i*evs_per_job+1, (i+1)*evs_per_job )

    #return

    
    ###################################################### TTJets SL
    ######################################################
        
    # TTJetsSemiLept --> 16749255
    sample  = 'TTJetsSemiLept'
    counter = 0
    if speedup:
        num_of_jobs =   4
        evs_per_job =   4200000 
    else:
        num_of_jobs =       150
        evs_per_job =    112000  

    for i in range(num_of_jobs):
        counter = counter + 1
        if speedup:
            submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample, i*evs_per_job+1, (i+1)*evs_per_job )
        else:
            submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample, i*evs_per_job+1, (i+1)*evs_per_job )

#    return
    ###################################################### TTJets FL
    ######################################################

    # TTJetsFullHad --> 8932897
    sample  = 'TTJetsFullLept'
    counter = 0
    if speedup:
        num_of_jobs =   1
        evs_per_job =   1
    else:
        num_of_jobs =     46
        evs_per_job = 202000  

    for i in range(num_of_jobs):
        counter = counter + 1
        if speedup:
            submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample, 0, -1 )
        else:
            submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample,  i*evs_per_job+1, (i+1)*evs_per_job )


    ###################################################### TTJets FH ######################################################
            
    # TTJetsFullLept --> 8932897
    sample  = 'TTJetsFullHad'
    submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(1), sample, 0, -1 )
    
    ###################################################### EWK
    ######################################################

    # DYJets10to50
    sample  = 'DYJets10to50'
    submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(1), sample, 0, -1 )

    # DYJets50
    sample  = 'DYJets50'
    counter = 0
    for i in range(1):
        counter = counter + 1
        submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample, 0, -1 )

 #   return

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
    if speedup:
        num_of_jobs =    1
        evs_per_job =    -1
    else:
        num_of_jobs =   15
        evs_per_job = 8000     # ---> ~ 40/job

    for i in range(num_of_jobs):
        counter = counter + 1
        if speedup:
            submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample, 0, -1 )
        else:
            submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample, i*evs_per_job+1, (i+1)*evs_per_job )

    # TTW
    sample  = 'TTW'
    counter = 0
    for i in range(1):
        counter = counter + 1
        submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample, 0, -1 )


    ###################################################### Data
    ######################################################

    datasamples = [
#        'DoubleElectron_Run2012A-13Jul2012-v1_ProcFIXED',
#        'DoubleElectron_Run2012A-recover-06Aug2012-v1_ProcV2',
#        'DoubleElectron_Run2012B-13Jul2012-v1_ProcFIXED',
#        'DoubleElectron_Run2012C-PromptReco-v2_HBB_EDMNtupleV42_ProcV1',
#        'DoubleElectron_Run2012C-PromptReco-v2_HBB_EDMNtupleV42_ProcV2',
#        'DoubleElectronRun2012C-EcalRecover_11Dec2012-v1_v2',
#        'DoubleElectronRun2012CAug24RerecoEdmV42',
#        'DoubleElectronRun2012D',

#        'SingleElectronRun2012AAug06EdmV42',
#        'SingleElectronRun2012AJul13EdmV42b',
#        'SingleElectronRun2012BJul13EdmV42',
#        'SingleElectronRun2012C-EcalRecover_11Dec2012-v1_v2',
#        'SingleElectronRun2012CAug24RerecoEdmV42',
#        'SingleElectronRun2012CPromptv2EdmV42',
#        'SingleElectronRun2012CPromptV2TopUpEdmV42',
#        'SingleElectronRun2012D-PromptReco-v1_v3',
        
        'SingleMuRun2012AAug06',
        'SingleMuRun2012AJul13',
        'SingleMuRun2012BJul13',
        'SingleMuRun2012C-EcalRecover_11Dec2012-v1_v2',
        'SingleMuRun2012CAug24Rereco',
        'SingleMuRun2012CPromptv2',
        'SingleMuRun2012CPromptV2TopUp',
        'SingleMuRun2012D-PromptReco-v1'
        ]

    for datasample in datasamples:

        sample  = 'Run2012_'+datasample

        counter     =    0
        num_of_jobs =    1
        evs_per_job =    1

        for i in range(1):
            counter = counter + 1
            submitMEAnalysisNew_all(analysis+'_'+sample+'_p'+str(counter), sample, 0, -1 )

#-----------------------------------------------------------------------------------------------------



###########################################
###########################################


analyses = ['nominal_0-0-1']

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
