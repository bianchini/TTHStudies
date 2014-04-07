#!/usr/bin/env python


import commands
import re
import os
import ROOT

import sys
sys.path.append('./')

import FWCore.ParameterSet.Config as cms

from Bianchi.TTHStudies.mem_categories_cff import *


###########################################
###########################################

PSI  = 0
qsub = ''
if PSI==1:
    qsub            = 'qsub -V -cwd -l h_vmem=2G -q all.q'
else:
    qsub            = 'qsub -q main'


###########################################
###########################################

    
def addZllVeto   ( process ):
    old_cut     = process.cut.value()
    process.cut =  cms.string("("+old_cut+") && (Vtype==4 || (Vtype!=4 && TMath::Abs(Mll-91.2)>8.))")
        
def addDiJetPtCut( process ):
    old_cut     = process.cut.value()
    process.cut =  cms.string("("+old_cut+") && hJetAmong>=2")
        

###########################################
###########################################
 

def submitDataCardMakerFWlite_Limits(category):

    print "Overload datacardMakerFWlite.py..."
    os.system('cp ../python/datacardMakerFWlite.py ./')

    from datacardMakerFWlite import process

    if   category == "cat1_sb":
        process.fwliteInput = cat1_sb.clone()
    elif category == "cat2_sb":
        process.fwliteInput = cat2_sb.clone()
    elif category == "cat3_sb":
        process.fwliteInput = cat3_sb.clone()
    elif category == "cat6_sb":
        process.fwliteInput = cat6_sb.clone()
    elif category == "cat1_sb_nb":
        process.fwliteInput = cat1_sb_nb.clone()
    elif category == "cat2_sb_nb":
        process.fwliteInput = cat2_sb_nb.clone()
    elif category == "cat3_sb_nb":
        process.fwliteInput = cat3_sb_nb.clone()
    elif category == "cat6_sb_nb":
        process.fwliteInput = cat6_sb_nb.clone()
    elif category == "cat1_bj":
        process.fwliteInput = cat1_bj.clone()
    elif category == "cat2_bj":
        process.fwliteInput = cat2_bj.clone()
    elif category == "cat3_bj":
        process.fwliteInput = cat3_bj.clone()
    elif category == "cat6_bj":
        process.fwliteInput = cat6_bj.clone()

    elif category == "cat1_sb_L":
        process.fwliteInput = cat1_sb_L.clone()
    elif category == "cat2_sb_L":
        process.fwliteInput = cat2_sb_L.clone()
    elif category == "cat3_sb_L":
        process.fwliteInput = cat3_sb_L.clone()
    elif category == "cat6_sb_L":
        process.fwliteInput = cat6_sb_L.clone()

    elif category == "cat1_sb_H":
        process.fwliteInput = cat1_sb_H.clone()
    elif category == "cat2_sb_H":
        process.fwliteInput = cat2_sb_H.clone()
    elif category == "cat3_sb_H":
        process.fwliteInput = cat3_sb_H.clone()
    elif category == "cat6_sb_H":
        process.fwliteInput = cat6_sb_H.clone()
        
    else:
        print "Cannot find this category... exit"
        return

    if ADDZLLVETO and re.search("cat6", category)!=None:
        addZllVeto( process.fwliteInput )

    print "Creating the shell file for the batch..."
    scriptName = 'job_'+category+'.sh'
    jobName    = 'job_'+category

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
    f.write('DataCardMakerNewFWlite ./'+jobName+'.py '+category+'\n')
    f.write('rm /scratch/bianchi/dummy*.root\n')
    f.close()
    os.system('chmod +x '+scriptName)

    submitToQueue = qsub+' -N job'+category+' '+scriptName 
    print submitToQueue
    #os.system(submitToQueue)
    
    print "\n@@@@@ END JOB @@@@@@@@@@@@@@@"

###########################################
###########################################
 

def submitDataCardMakerFWlite_Limits_Optimization(category, extracut, trial, fact1=-99):

    print "Overload datacardMakerFWlite.py..."
    os.system('cp ../python/datacardMakerFWlite.py ./')

    from datacardMakerFWlite import process


    if   category == "cat1_sb":
        process.fwliteInput = cat1_sb.clone()
    elif category == "cat2_sb":
        process.fwliteInput = cat2_sb.clone()
    elif category == "cat3_sb":
        process.fwliteInput = cat3_sb.clone()
    elif category == "cat6_sb":
        process.fwliteInput = cat6_sb.clone()
    elif category == "cat1_sb_nb":
        process.fwliteInput = cat1_sb_nb.clone()
    elif category == "cat2_sb_nb":
        process.fwliteInput = cat2_sb_nb.clone()
    elif category == "cat3_sb_nb":
        process.fwliteInput = cat3_sb_nb.clone()
    elif category == "cat6_sb_nb":
        process.fwliteInput = cat6_sb_nb.clone()
    elif category == "cat1_bj":
        process.fwliteInput = cat1_bj.clone()
    elif category == "cat2_bj":
        process.fwliteInput = cat2_bj.clone()
    elif category == "cat3_bj":
        process.fwliteInput = cat3_bj.clone()
    elif category == "cat6_bj":
        process.fwliteInput = cat6_bj.clone()

    elif category == "cat1_sb_L":
        process.fwliteInput = cat1_sb_L.clone()
    elif category == "cat2_sb_L":
        process.fwliteInput = cat2_sb_L.clone()
    elif category == "cat3_sb_L":
        process.fwliteInput = cat3_sb_L.clone()
    elif category == "cat6_sb_L":
        process.fwliteInput = cat6_sb_L.clone()

    elif category == "cat1_sb_H":
        process.fwliteInput = cat1_sb_H.clone()
    elif category == "cat2_sb_H":
        process.fwliteInput = cat2_sb_H.clone()
    elif category == "cat3_sb_H":
        process.fwliteInput = cat3_sb_H.clone()
    elif category == "cat6_sb_H":
        process.fwliteInput = cat6_sb_H.clone()
        
    else:
        print "Cannot find this category... exit"
        return     

    oldcut = process.fwliteInput.cut.value()
    newcut = oldcut+" && "+extracut
    process.fwliteInput.cut       = cms.string(newcut)
    process.fwliteInput.extraname = cms.string("_sb_"+trial)

    if fact1>=0:
         process.fwliteInput.fact1 = cms.double(fact1)

    if ADDZLLVETO and re.search("cat6", category)!=None:
        addZllVeto( process.fwliteInput )

    print "Creating the shell file for the batch..."
    scriptName = 'job_'+category+'_'+trial+'.sh'
    jobName    = 'job_'+category+'_'+trial

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
    f.write('DataCardMakerNewFWlite ./'+jobName+'.py\n')
    f.write('rm /scratch/bianchi/dummy*.root\n')
    f.close()
    os.system('chmod +x '+scriptName)

    submitToQueue = qsub+' -N job'+category+' '+scriptName 
    print submitToQueue
    os.system(submitToQueue)
    
    print "\n@@@@@ END JOB @@@@@@@@@@@@@@@"


###########################################
###########################################
 

def submitDataCardMakerFWlite(category, cut, script, samples, extraname, nparts, part, binvec):

    print "Overload datacardMakerFWlite.py..."
    os.system('cp ../python/datacardMakerFWlite.py ./')

    from datacardMakerFWlite import process

    process.fwliteInput.samples   = samples    
    process.fwliteInput.category  = cms.string( category )
    process.fwliteInput.cut       = cms.string( cut )
    process.fwliteInput.extraname = cms.string( extraname )
    process.fwliteInput.nparts    = cms.int32(nparts)
    process.fwliteInput.part      = cms.int32(part)
    process.fwliteInput.binvec    = cms.vdouble(binvec)
    process.fwliteInput.nBins     = cms.int32(len(binvec)-1)

    if ADDZLLVETO and re.search("cat6", category)!=None:
        addZllVeto( process.fwliteInput )
        
    print "Creating the shell file for the batch..."
    scriptName = 'job_'+script+'.sh'
    jobName    = 'job_'+script

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
    f.write('DataCardMakerNewFWlite ./'+jobName+'.py\n')
    f.write('rm /scratch/bianchi/dummy*.root\n')
    f.close()
    os.system('chmod +x '+scriptName)

    submitToQueue = qsub+' -N job'+category+' '+scriptName 
    print submitToQueue
    os.system(submitToQueue)
    
    print "\n@@@@@ END JOB @@@@@@@@@@@@@@@"
        
###########################################
###########################################



def submitDataCardMakerFWlite_all( category , cut, selection, binvec):



    sampless = [ #[["TTV"],     1],
                 #[["SingleT"], 1],
                 #[["DiBoson"], 1],
                 [["TTJetsBB"], 20],
                 #[["TTJetsBJ"],22],
                 #[["TTJetsJJ"],50],
                 #[["TTH125"],  1],
                 #[["EWK"],     1],
                 #[["Run2012_SingleMu", "Run2012_SingleElectron"],1 ]
                 #[["Run2012_SingleElectron"],1 ]
                 ]

    counter = 0
    for samples in sampless:
        toBeRun = cms.vstring( samples[0] )
        outfile = ""
        for sample in samples[0]:
            if sample != (samples[0])[len(samples[0])-1]:
                outfile += (sample+"-")
            else:
                outfile += (sample)
        for split in range( samples[1] ):
            counter += 1
            script  = category+"_p"+str(counter)
            submitDataCardMakerFWlite( category, cut, script, toBeRun, "_"+selection+"_"+outfile+"_"+str(split), samples[1], split, binvec)
      
      
###########################################
###########################################
binvec = cms.vdouble(30., 40, 50., 60., 70., 80., 90., 100., 110., 120., 130., 140., 150., 160., 170., 180.)
#submitDataCardMakerFWlite_all( "lepton_pt", "(numJets>=6 && numBTagM==3)", "test" , binvec)


#submitDataCardMakerFWlite_Limits("cat1_sb_L")
#submitDataCardMakerFWlite_Limits("cat2_sb_L")
#submitDataCardMakerFWlite_Limits("cat3_sb_L")
#submitDataCardMakerFWlite_Limits("cat1_sb_H")
#submitDataCardMakerFWlite_Limits("cat2_sb_H")
#submitDataCardMakerFWlite_Limits("cat3_sb_H")
#submitDataCardMakerFWlite_Limits("cat6_sb_H")

#submitDataCardMakerFWlite_Limits("cat1_sb")
#submitDataCardMakerFWlite_Limits("cat2_sb")
#submitDataCardMakerFWlite_Limits("cat3_sb")
#submitDataCardMakerFWlite_Limits("cat6_sb")

#submitDataCardMakerFWlite_Limits("cat1_sb_nb")
#submitDataCardMakerFWlite_Limits("cat2_sb_nb")
#submitDataCardMakerFWlite_Limits("cat3_sb_nb")
#submitDataCardMakerFWlite_Limits("cat6_sb_nb")

#submitDataCardMakerFWlite_Limits("cat1_bj")
#submitDataCardMakerFWlite_Limits("cat2_bj")
#submitDataCardMakerFWlite_Limits("cat3_bj")
#submitDataCardMakerFWlite_Limits("cat6_bj")

########################################### optimize btagLR lower cut

cuts = [0.85, 0.875, 0.90, 0.925 , 0.950, 0.975, 0.980, 0.990 ]

trial = 0
for cut in cuts:
    #submitDataCardMakerFWlite_Limits_Optimization("cat1_sb",  ("btag_LR>=%f" % cut), "cat1_"+str(trial) )
    #submitDataCardMakerFWlite_Limits_Optimization("cat2_sb",  ("btag_LR>=%f" % cut), "cat2_"+str(trial) )
    #submitDataCardMakerFWlite_Limits_Optimization("cat3_sb",  ("btag_LR>=%f" % cut), "cat3_"+str(trial) )
    #submitDataCardMakerFWlite_Limits_Optimization("cat6_sb",  ("btag_LR>=%f" % cut), "cat6_"+str(trial) )
    trial += 1


########################################### optimize btagLR splitting


cuts =  [0.80, 0.825, 0.850, 0.875, 0.900]

trial = 0
for cut in range(len(cuts)):
    for trial in range(2):
      if trial==0:
          #submitDataCardMakerFWlite_Limits_Optimization("cat1_sb",  ("btag_LR>=%f" % cuts[cut] ), "cat1-"+str(cut)+"_"+str(trial) )
          #submitDataCardMakerFWlite_Limits_Optimization("cat2_sb",  ("btag_LR>=%f" % cuts[cut] ), "cat2-"+str(cut)+"_"+str(trial) )
          #submitDataCardMakerFWlite_Limits_Optimization("cat3_sb",  ("btag_LR>=%f" % cuts[cut] ), "cat3-"+str(cut)+"_"+str(trial) )
          #submitDataCardMakerFWlite_Limits_Optimization("cat6_sb",  ("btag_LR>=%f" % cuts[cut] ), "cat6-"+str(cut)+"_"+str(trial) )

          #submitDataCardMakerFWlite_Limits_Optimization("cat1_sb",  ("btag_LR>=%f" % 0.995 ), "cat1-"+str(cut)+"_"+str(trial) )
          #submitDataCardMakerFWlite_Limits_Optimization("cat2_sb",  ("btag_LR>=%f" % 0.9925), "cat2-"+str(cut)+"_"+str(trial) )
          #submitDataCardMakerFWlite_Limits_Optimization("cat3_sb",  ("btag_LR>=%f" % 0.995),  "cat3-"+str(cut)+"_"+str(trial) )
          #submitDataCardMakerFWlite_Limits_Optimization("cat6_sb",  ("btag_LR>=%f" % 0.925),  "cat6-"+str(cut)+"_"+str(trial) )
          trial += 0
      else:
          #submitDataCardMakerFWlite_Limits_Optimization("cat1_sb",  ("btag_LR<%f"  % cuts[cut] ), "cat1-"+str(cut)+"_"+str(trial) )
          #submitDataCardMakerFWlite_Limits_Optimization("cat2_sb",  ("btag_LR<%f"  % cuts[cut] ), "cat2-"+str(cut)+"_"+str(trial) )
          #submitDataCardMakerFWlite_Limits_Optimization("cat3_sb",  ("btag_LR<%f"  % cuts[cut] ), "cat3-"+str(cut)+"_"+str(trial) )
          #submitDataCardMakerFWlite_Limits_Optimization("cat6_sb",  ("btag_LR<%f"  % cuts[cut] ), "cat6-"+str(cut)+"_"+str(trial) )

          #submitDataCardMakerFWlite_Limits_Optimization("cat1_sb",  ("btag_LR<%f && btag_LR>=%f" % (0.995, cuts[cut])  ), "cat1-"+str(cut)+"_"+str(trial) )
          #submitDataCardMakerFWlite_Limits_Optimization("cat2_sb",  ("btag_LR<%f && btag_LR>=%f" % (0.9925,cuts[cut])  ), "cat2-"+str(cut)+"_"+str(trial) )
          #submitDataCardMakerFWlite_Limits_Optimization("cat3_sb",  ("btag_LR<%f && btag_LR>=%f" % (0.995, cuts[cut])  ), "cat3-"+str(cut)+"_"+str(trial) )
          #submitDataCardMakerFWlite_Limits_Optimization("cat6_sb",  ("btag_LR<%f && btag_LR>=%f" % (0.925, cuts[cut])  ), "cat6-"+str(cut)+"_"+str(trial) )
          trial += 0


########################################### optimize S/B normalization


cuts =  [1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8]

trial = 0
for cut in cuts:
    #submitDataCardMakerFWlite_Limits_Optimization("cat1_sb_L",  ("btag_LR>=%f" % 0.), "cat1_"+str(trial), cut )
    #submitDataCardMakerFWlite_Limits_Optimization("cat2_sb_L",  ("btag_LR>=%f" % 0.), "cat2_"+str(trial), cut )
    #submitDataCardMakerFWlite_Limits_Optimization("cat3_sb_L",  ("btag_LR>=%f" % 0.), "cat3_"+str(trial), cut )
    submitDataCardMakerFWlite_Limits_Optimization("cat6_sb_L",  ("btag_LR<%f" % 0.925), "cat6_"+str(trial), cut )
    trial += 1



