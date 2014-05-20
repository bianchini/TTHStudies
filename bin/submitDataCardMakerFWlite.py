#!/usr/bin/env python

import commands
import re
import os
import ROOT
import subprocess
from time import sleep

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
    process.cut =  cms.string("("+old_cut+") && (Vtype==2 || Vtype==3 || Vtype==4 || (Vtype<=1 && TMath::Abs(Mll-91.2)>8. && Mll>15.))")

def addPixelVeto ( process ):
    old_cut     = process.cut.value()
    process.cut =  cms.string("("+old_cut+") && (EVENT.run < 207883 || EVENT.run > 208307)")
    
def addDiJetPtCut( process ):
    old_cut     = process.cut.value()
    process.cut =  cms.string("("+old_cut+") && hJetAmong>=2")
        
def addJetPt40Cut ( process, numJets ) :
    old_cut = process.cut.value()
    process.cut = cms.string( "(" + old_cut + ") && jetsAboveCut >= " + str(numJets) ) 

###########################################
###########################################
 

def submitDataCardMakerFWlite_Limits(category):

    if len(sys.argv)<2:
        print "Specify how to run the job: [batch,local]"
        return
    
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

    if ADDZLLVETO:
        addZllVeto( process.fwliteInput )

    if ADDPIXELVETO:
        addPixelVeto ( process.fwliteInput )

    if ADDDIJETPTCUT:
        addDiJetPtCut( process.fwliteInput )

    scriptName = 'job_'+category+'.sh'
    jobName    = 'job_'+category


    out = open(jobName+'.py','w')
    out.write(process.dumpPython())
    out.close()
    
    print "Creating the shell file for the batch..."
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

    if sys.argv[1]== "batch":
        print submitToQueue
        os.system(submitToQueue)
    elif sys.argv[1]== "local":
        subprocess.call( ['DataCardMakerNewFWlite', jobName+'.py', category] )
    else:
        print "Unsupported job type"
        return
    
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

    if ADDZLLVETO:
        addZllVeto( process.fwliteInput )

    if ADDPIXELVETO:
        addPixelVeto( process.fwliteInput )

    if ADDDIJETPTCUT:
        addDiJetPtCut( process.fwliteInput )

    
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
 

def submitDataCardMakerFWlite(varname, category, cut, script, samples, extraname, nparts, part, binvec, analysis, inputpath="", version="", outdir=""):

    if len(sys.argv)<2:
        print "Specify how to run the job: [batch,local]"
        return
    
    print "Overload datacardMakerFWlite.py..."
    os.system('cp ../python/datacardMakerFWlite.py ./')

    from datacardMakerFWlite import process

    process.fwliteInput = cat.clone()

    process.fwliteInput.directory = cms.string( process.fwliteInput.directory.value()+"/controlPlots")
    process.fwliteInput.samples   = samples    
    process.fwliteInput.category  = cms.string( category )
    process.fwliteInput.varname   = cms.string( varname )
    process.fwliteInput.cut       = cms.string( cut )
    process.fwliteInput.extraname = cms.string( extraname )
    process.fwliteInput.nparts    = cms.int32(nparts)
    process.fwliteInput.part      = cms.int32(part)
    process.fwliteInput.binvec    = cms.vdouble(binvec)
    process.fwliteInput.nBins     = cms.int32(len(binvec)-1)
    if inputpath:
        process.fwliteInput.inputpath = cms.string(inputpath) #otherwise default
    if outdir:
        process.fwliteInput.directory = cms.string(outdir)
    if version:
        process.fwliteInput.version = cms.string(version)
    
    process.fwliteInput.analysis  = cms.untracked.int32(analysis)

    if ADDZLLVETO:
        addZllVeto( process.fwliteInput )

    if ADDPIXELVETO:
        addPixelVeto( process.fwliteInput )

    if ADDDIJETPTCUT:
        addDiJetPtCut( process.fwliteInput )

    if ADDJETPT40CUT:
        addJetPt40Cut( process.fwliteInput, NR_PT40_JETS )
        
    print "Creating the shell file for the batch..."
    scriptName = 'job_'+script+'.sh'
    jobName    = 'job_'+script

    out = open(jobName+'.py','w')
    out.write(process.dumpPython())
    out.close()
   
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

    submitToQueue = qsub+' -N job_'+script+' '+scriptName

    if sys.argv[1]== "batch":
        print submitToQueue
        os.system(submitToQueue)
    elif sys.argv[1]== "local":
        subprocess.call( ['DataCardMakerNewFWlite', jobName+'.py', category] )
    else:
        print "Unsupported job type"
        return
    
    print "\n@@@@@ END JOB @@@@@@@@@@@@@@@"
        
###########################################
###########################################

sampless_ALL = [
    [["TTV"],     1],
    [["SingleT"], 1],
    [["DiBoson"], 1],
    [["TTJetsBB"],1],
    [["TTJetsBJ"],1],
    [["TTJetsJJ"],1],
    [["TTH125"],  1],
    [["EWK"],     1],
    [["Run2012_SingleMu", "Run2012_SingleElectron"],1 ],
    ]


sampless_SL = [
    [["TTV", "SingleT", "DiBoson","EWK" ],          1],
    [["TTJetsBB", "TTJetsBJ", "TTJetsJJ"],          1],
    [["TTH125"],                                    1],
    [["Run2012_SingleMu", "Run2012_SingleElectron"],1],
    ]

sampless_DL = [
    [["TTV", "SingleT", "DiBoson","EWK" ],          1],
    [["TTJetsBB", "TTJetsBJ", "TTJetsJJ"],          1],
    [["TTH125"],                                    1],
    [["Run2012_SingleMu", "Run2012_DoubleElectron"],1],
    ]

def submitDataCardMakerFWlite_all( varname, category, cut, selection, binvec, analysis, sampless=[], inputpath="", version="", outdir=""):
    """
    Optional arguments (otherwise use default values from ../python/mem_categories.py):
    inputpath -- path for inputfiles
    outdir -- directory, where the output root files are saved in ../root/datacards/
    version -- production version of MEM ntuples (e.g. '_ntuplizeAll_v3_rec_std')
    """
    counter = 0

    if sampless == []:
        if analysis==0:
            sampless = sampless_SL
        else:
            sampless = sampless_DL

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
            script  = selection + "_" + category + "_p"+str(counter)

            submitDataCardMakerFWlite( varname, category, cut, script, toBeRun, "_"+selection+"_"+outfile+"_"+str(split), samples[1], split, binvec, analysis, inputpath, version, outdir)
            sleep(1)

      
###########################################
###########################################



