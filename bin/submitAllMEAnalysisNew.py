#!/usr/bin/env python


import commands
import re
import os

import sys
sys.path.append('./')

import FWCore.ParameterSet.Config as cms

# mass scan

massesH     = cms.vdouble(125)
#massesH     = cms.vdouble(45., 55., 65., 75., 85., 95., 105., 115., 125., 135., 145., 155., 165., 185., 205., 225., 250., 275., 300.)
massesT     = cms.vdouble(174.3)
#massesT     = cms.vdouble(115., 125., 135., 145., 155., 165., 175., 185., 195., 205.,215., 225.)

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

# cut values to select events
btag_prob_cut_6jets = 0.96675 # <--- 0.988
btag_prob_cut_5jets = 0.98225 # <--- 0.992
btag_prob_cut_4jets = 0.95295 # <--- 0.85 <--- 0.95295 <--- 0.992

# regression
useRegression = 0

# use gen-jets or reco-jets ???
doGenLevelAnalysis = 1

# btag-thresholds
csv_WP_L = 0.244
csv_WP_M = 0.600   #0.679
csv_WP_T = 0.898

# select bt btag_LR
selectByBTagShape = 0

# an extra name for the output files
extraoutname = ""
if len(massesH)>1:
    extraoutname = "MHscan_"
if len(massesT)>1:
    extraoutname = "MTscan_"

# run integral optimization
integralOption2 = 1


def submitMEAnalysisNew(type,
                        script,
                        sample,
                        doCSVup, doCSVdown, doJECup, doJECdown, doJERup, doJERdown,
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

    #if VType!=("_"+type):
    #    print "Attention! Mismatch between VType... Please, check and run again."
    #    return

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
    process.fwliteInput.useME            = cms.int32(useME)
    process.fwliteInput.useJac           = cms.int32(useJac)
    process.fwliteInput.useMET           = cms.int32(useMET)
    process.fwliteInput.useTF            = cms.int32(useTF)
    process.fwliteInput.usePDF           = cms.int32(usePDF)

    process.fwliteInput.integralOption2  = cms.untracked.int32( integralOption2 )

    process.fwliteInput.selectByBTagShape = cms.untracked.int32(selectByBTagShape)
    
    if re.search("VType0",   type )!=None or re.search("VType1",   type )!=None:
        process.fwliteInput.doType6ByBTagShape = cms.untracked.int32(    selectByBTagShape)
        process.fwliteInput.doType6            = cms.untracked.int32(not selectByBTagShape)
        process.fwliteInput.doType7            = cms.untracked.int32(0)
    else:
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
            #process.fwliteInput.doType1            = cms.untracked.int32(not selectByBTagShape)
            process.fwliteInput.doType1            = cms.untracked.int32(0)
        else:
            process.fwliteInput.doType1            = cms.untracked.int32(0)

        process.fwliteInput.doType2            = cms.untracked.int32(not selectByBTagShape)

        if len(massesH) == 1 and len(massesT) == 1:
            process.fwliteInput.doType3            = cms.untracked.int32(not selectByBTagShape)
        else:
            process.fwliteInput.doType3            = cms.untracked.int32(0)

    
    process.fwliteInput.btag_prob_cut_6jets = cms.untracked.double(btag_prob_cut_6jets)
    process.fwliteInput.btag_prob_cut_5jets = cms.untracked.double(btag_prob_cut_5jets)
    process.fwliteInput.btag_prob_cut_4jets = cms.untracked.double(btag_prob_cut_4jets)

    process.fwliteInput.csv_WP_L            =  cms.untracked.double(csv_WP_L)
    process.fwliteInput.csv_WP_M            =  cms.untracked.double(csv_WP_M)
    process.fwliteInput.csv_WP_T            =  cms.untracked.double(csv_WP_T)

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
    process.fwliteInput.doJERup          = cms.untracked.int32(doJERup)
    process.fwliteInput.doJERdown        = cms.untracked.int32(doJERdown)

    process.fwliteInput.doGenLevelAnalysis  = cms.untracked.int32(doGenLevelAnalysis)
    process.fwliteInput.speedup             = cms.untracked.int32(0)
    process.fwliteInput.ntuplizeAll         = cms.untracked.int32(0)

 
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
    doJERup  = 0
    doJERdown= 0
    
    if   re.search("csvUp",   analysis )!=None:
        doCSVup   = 1
    elif re.search("csvDown", analysis )!=None:
        doCSVdown = 1
    elif re.search("JECUp",   analysis )!=None:
        doJECup   = 1
    elif re.search("JECDown", analysis )!=None:
        doJECdown = 1
    elif re.search("JERUp",   analysis )!=None:
        doJERup   = 1
    elif re.search("JERDown", analysis )!=None:
        doJERdown = 1
    else:
        print "Doing nominal analysis"


    ###################################################### TTJets SL
    ######################################################
        
    # TTJetsSemiLept --> 16749255
    sample  = 'TTJetsSemiLept'
    counter = 0
    num_of_jobs =   1
    evs_per_job =   1 
    if doSL:
        num_of_jobs =       150
        evs_per_job =    112000  # ---> ~ 40*2/job
    else:
        num_of_jobs =         4*2
        evs_per_job =   4250000/2  # ---> <~ 40/job
    for i in range(num_of_jobs):
        counter = counter + 1
        #submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown, doJERup, doJERdown,i*evs_per_job+1, (i+1)*evs_per_job )

    ###################################################### TTJets FL
    ######################################################

    # TTJetsFullLept --> 8932897
    sample  = 'TTJetsFullLept'
    counter = 0
    num_of_jobs =   1
    evs_per_job =   1
    if doSL:
         num_of_jobs =     22
         evs_per_job = 419000  # ---> ~ 40*2/job
    else:
         num_of_jobs =     23*2
         evs_per_job = 404000/2  # ---> ~ 40/job
    for i in range(num_of_jobs):
        counter = counter + 1
        #submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown, doJERup, doJERdown, i*evs_per_job+1, (i+1)*evs_per_job )

    ###################################################### TTH
    ######################################################

    # TTH125 --> 194808
    sample  = 'TTH125'
    counter = 0
    num_of_jobs =   1
    evs_per_job =   1
    if doSL:
         num_of_jobs =  180
         evs_per_job = 1085  # ---> ~ 40/job
    else:
         num_of_jobs =    10*2
         evs_per_job = 21000/2 # ---> ~ 40/job
    for i in range(num_of_jobs):
        counter = counter + 1
        submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown, doJERup, doJERdown, i*evs_per_job+1, (i+1)*evs_per_job )


    if re.search("nominal",   analysis )==None:
        return


    ###################################################### EWK
    ######################################################

    # DYJets10to50
    sample  = 'DYJets10to50'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown, doJERup, doJERdown,0, -1 )

    # DYJets50
    sample  = 'DYJets50'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown, doJERup, doJERdown,0, -1 )

    # WJets 
    sample  = 'WJets'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown, doJERup, doJERdown,0, -1 )

    ###################################################### Single-top
    ######################################################

    # TtW
    sample  = 'TtW'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown, doJERup, doJERdown,0, -1 )
        
   # Tt
    sample  = 'Tt'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown, doJERup, doJERdown,0, -1 )

   # Ts
    sample  = 'Ts'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown, doJERup, doJERdown,0, -1 )

   # TbartW
    sample  = 'TbartW'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown, doJERup, doJERdown,0, -1 )

   # Tbart
    sample  = 'Tbart'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown, doJERup, doJERdown,0, -1 )

   # Tbars
    sample  = 'Tbars'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown, doJERup, doJERdown,0, -1 )


    ###################################################### Di-boson
    ######################################################


   # WW
    sample  = 'WW'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown, doJERup, doJERdown,0, -1 )

   # WZ
    sample  = 'WZ'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown, doJERup, doJERdown,0, -1 )

   # ZZ
    sample  = 'ZZ'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown, doJERup, doJERdown,0, -1 )


    ###################################################### TTZ
    ######################################################

    # TTZ   --> 112517
    sample  = 'TTZ'
    counter = 0
    num_of_jobs =    1
    evs_per_job =    1
    if doSL:
         num_of_jobs =   15
         evs_per_job = 8000     # ---> ~ 40/job
    else:
         num_of_jobs =     2*2
         evs_per_job = 60000/2    # ---> ~ 40/job
    for i in range(num_of_jobs):
        counter = counter + 1
        #submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown, doJERup, doJERdown,i*evs_per_job+1, (i+1)*evs_per_job )

    # TTW
    sample  = 'TTW'
    counter = 0
    for i in range(1):
        counter = counter + 1
        #submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown, doJERup, doJERdown,0, -1 )


    ###################################################### Data
    ######################################################

    datasamples = [#'DoubleElectron_Run2012A-13Jul2012-v1_ProcFIXED',
                   #'DoubleElectron_Run2012A-recover-06Aug2012-v1_ProcV2',
                   #'DoubleElectron_Run2012B-13Jul2012-v1_ProcFIXED',
                   #'DoubleElectron_Run2012C-PromptReco-v2_HBB_EDMNtupleV42_ProcV1',
                   #'DoubleElectron_Run2012C-PromptReco-v2_HBB_EDMNtupleV42_ProcV2',
                   #'DoubleElectronRun2012C-EcalRecover_11Dec2012-v1_v2',
                   #'DoubleElectronRun2012CAug24RerecoEdmV42',
                   #'DoubleElectronRun2012D',

                   #'SingleElectronRun2012AAug06EdmV42',
                   #'SingleElectronRun2012AJul13EdmV42b',
                   #'SingleElectronRun2012BJul13EdmV42',
                   #'SingleElectronRun2012C-EcalRecover_11Dec2012-v1_v2',
                   #'SingleElectronRun2012CAug24RerecoEdmV42',
                   #'SingleElectronRun2012CPromptv2EdmV42',
                   #'SingleElectronRun2012CPromptV2TopUpEdmV42',
                   #'SingleElectronRun2012D-PromptReco-v1_v3',

                   #'SingleMuRun2012AAug06EdmV42',
                   #'SingleMuRun2012AJul13EdmV42',
                   #'SingleMuRun2012BJul13EdmV42',
                   #'SingleMuRun2012C-EcalRecover_11Dec2012-v1_v2',
                   #'SingleMuRun2012CAug24RerecoEdmV42',
                   #'SingleMuRun2012CPromptv2EdmV42',
                   #'SingleMuRun2012CPromptV2TopUpEdmV42',
                   'SingleMuRun2012D-PromptReco-v1']

    for datasample in datasamples:

        continue
        
        sample  = 'Run2012_'+datasample

        counter     =    0
        num_of_jobs =    1
        evs_per_job =    1

        #print sample
        
        if (doSL and re.search('Double',datasample)!=None):
            continue


        if (doDL):
            for i in range(1):
                counter = counter + 1
                submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown, doJERup, doJERdown,0, -1 )

                      
        if (doSL and re.search("Run2012AAug06EdmV42",datasample)!=None):  # 148139
            for i in range(1):
                counter = counter + 1
                submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown, doJERup, doJERdown,0, -1 )

        if (doSL and re.search("Run2012AJul13EdmV42b", datasample)!=None):  # 1551019
            for i in range(1):
                counter = counter + 1
                submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown, doJERup, doJERdown,0, -1 )

        if (doSL and re.search("Run2012BJul13EdmV42", datasample)!=None):  # 9351330
            num_of_jobs =        4
            evs_per_job =  2340000
            for i in range(num_of_jobs):
                counter = counter + 1
                submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown, doJERup, doJERdown, i*evs_per_job+1, (i+1)*evs_per_job  )

        if (doSL and re.search("Run2012C-EcalRecover_11Dec2012-v1_v2", datasample)!=None):  # 263593
            for i in range(1):
                counter = counter + 1
                submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown, doJERup, doJERdown,0, -1 )

        if (doSL and re.search("Run2012CAug24RerecoEdmV42", datasample)!=None):  # 1064158
            for i in range(1):
                counter = counter + 1
                submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown, doJERup, doJERdown,0, -1 )

        if (doSL and re.search("Run2012CPromptv2EdmV42", datasample)!=None):  # 9768094
            num_of_jobs =        4
            evs_per_job =  2443000
            for i in range(num_of_jobs):
                counter = counter + 1
                submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown, doJERup, doJERdown, i*evs_per_job+1, (i+1)*evs_per_job  )

        if (doSL and re.search("Run2012CPromptV2TopUpEdmV42", datasample)!=None):  # 3491407
            for i in range(1):
                counter = counter + 1
                submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown, doJERup, doJERdown,0, -1 )

        if (doSL and re.search("Run2012D-PromptReco-v1", datasample)!=None):  # 16178887
            num_of_jobs =        6*6
            evs_per_job =  2800000/6
            for i in range(num_of_jobs):
                counter = counter + 1
                submitMEAnalysisNew(type,'SL_'+type+'_'+analysis+'_'+sample+'_p'+str(counter), sample, doCSVup, doCSVdown, doJECup, doJECdown, doJERup, doJERdown, i*evs_per_job+1, (i+1)*evs_per_job  )


###########################################
###########################################


analyses = ['nominal_test',
            #'csvUp_v2',
            #'csvDown_v2',
            #'JECUp_v2',
            #'JECDown_v2',
            #'JERUp',
            #'JERDown',
            ]

types = ['VType0',
         ####'VType1',
         'VType2',
         ####'VType3',
         ]

for type in types:
    for analysis in analyses:
        if doGenLevelAnalysis:
            analysis = analysis+'_gen'
        else:
            analysis = analysis+'_rec'   
        if useRegression:
            analysis = analysis+'_reg'
        else:
            analysis = analysis+'_std'
        submitFullMEAnalysisNew(type, analysis)
