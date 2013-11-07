#!/usr/bin/env python


import commands
import re
import os

import sys
sys.path.append('./')

from skim  import process as skimP
from step3 import process as step3P
import FWCore.ParameterSet.Config as cms

###########################################
###########################################

pathI   = "srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store//user/bianchi/HBB_EDMNtuple/AllHDiJetPt/"
pathI2  = "srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store//user/bianchi/HBB_EDMNtuple/AllHDiJetPt_Step2/"
pathI3  = "gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt/"
pathI4  = "gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/AllHDiJetPt_Step2/"
pathI5  = "gsidcap://t3se01.psi.ch:22128//pnfs/psi.ch/cms/trivcat/store/user/bianchi/HBB_EDMNtuple/" 

pathO   = "/scratch/bianchi/HBB_EDMNtuple/All.H.DiJetPt/"
fileN   = "DiJetPt_"
command ="srmcp -2 "


def processAll(samples):

    for sample in samples:

        output = commands.getoutput("ls -ltr "+pathO+fileN+sample[0]+".root")
        outFiles = re.split(r'\n',output)
        for name in outFiles:
            if re.search("file",name)!=None:
                print command+pathI+fileN+sample[0]+".root" +" file:///"+pathO+fileN+sample[0]+".root"
                os.system(command+pathI+fileN+sample[0]+".root" +" file:///"+pathO+fileN+sample[0]+".root")
                print "Step3 step3.py "+sample[0]+" "+sample[1]
                os.system("Step3 step3.py "+sample[0]+" "+sample[1])
                print "Skim  skim.py "+sample[0]+" "+sample[1]
                os.system("Skim  skim.py "+sample[0]+" "+sample[1])
                print "rm "+pathO+"*"+sample[0]+"_VType*.root"
                os.system("rm "+pathO+"*"+sample[0]+"_VType*.root")
            else:
                print "File is there... skimming"
                print "Skim  skim.py "+sample[0]+" "+sample[1]
                os.system("Skim  skim.py "+sample[0]+" "+sample[1])
                print "rm "+pathO+"*"+sample[0]+"_VType*.root"
                os.system("rm "+pathO+"*"+sample[0]+"_VType*.root")

###########################################
###########################################

def processAllDebugBatch(samples):

    step3P.fwliteInput.pathToFile = cms.string(pathI)
    step3P.fwliteInput.outPath    = cms.string(pathI2)

    skimP.fwliteInput.pathToFile  = cms.string(pathI)
    skimP.fwliteInput.outPath     = cms.string(pathI2)

    out = open('step3_tmp.py','w')
    out.write(step3P.dumpPython())
    out = open('skim_tmp.py','w')
    out.write(skimP.dumpPython())
    
    for sample in samples:

        print "Step3Clone step3_tmp.py "+sample[0]+" "+sample[1]
        #os.system("Step3Clone step3.py "+sample[0]+" "+sample[1])
        print "Skim  skim_tmp.py "+sample[0]+" "+sample[1]
        #os.system("Skim  skim.py "+sample[0]+" "+sample[1])

    os.system('rm step3_tmp.py')
    os.system('rm skim_tmp.py')

###########################################
###########################################

def processAllBatch(samples, skims, redoStep3):

    print "Processing: "
    i = 0
    while i < len(samples):
        print samples[i][1]+", "
        i = i + 1 
    print "Skims: "+str(len(skims))

    if redoStep3:
        print "Step3 will be redone. Files will be stored in: \n"
        print pathI2
        
    step3P.fwliteInput.pathToFile = cms.string(pathI3)
    step3P.fwliteInput.outPath    = cms.string(pathI4)

    skimP.fwliteInput.pathToFile  = cms.string(pathI4)
    skimP.fwliteInput.outPath     = cms.string(pathI5)

    out = open('step3_tmp.py','w')
    out.write(step3P.dumpPython())
    out = open('skim_tmp.py','w')
    out.write(skimP.dumpPython())

    readSkimsFromFile = False

    if len(skims)==0:
        print "Skims read from file:"
        for k in skimP.fwliteInput.skims:
            print k.name
        skims = [['','']]
        readSkimsFromFile = True

    print "\n@@@@@ BEGIN JOBS @@@@@@@@@@@@@@@\n"
        
    for skim in skims:
        for sample in samples:

            print "Creating the shell file for the batch..."
            scriptName = 'job_'+sample[1]
            if not readSkimsFromFile:
                scriptName = scriptName+'_'+skim[0]
            scriptName = scriptName+'.sh'
                
            f = open(scriptName,'w')
            f.write('#!/bin/bash\n\n')
            f.write('cd /shome/bianchi/CMSSW_5_3_3_patch2/src/Bianchi/TTHStudies/bin/\n')
            f.write('source /swshare/psit3/etc/profile.d/cms_ui_env.sh\n')
            f.write('export SCRAM_ARCH="slc5_amd64_gcc462"\n')
            f.write('source $VO_CMS_SW_DIR/cmsset_default.sh\n')
            f.write('eval `scramv1 runtime -sh`\n')
            f.write('export LD_PRELOAD="libglobus_gssapi_gsi_gcc64pthr.so.0":${LD_PRELOAD}\n')
            f.write('\n\n')
            f.write('\n\n')
            
            #print "srmls "+pathI2+fileN+sample[0]+".root"
            output = commands.getoutput("srmls "+pathI2+fileN+sample[0]+".root")
            outFiles = re.split(r'\n',output)
            skipStep3 = True
            for name in outFiles:
                if re.search("FAILURE",name)!=None:
                    skipStep3 = False
            skipStep3 = (skipStep3 and not redoStep3)
            if skipStep3:
                print "File "+fileN+sample[0]+".root is already at destination, proceed with Skim..."
            else:
                print "File "+fileN+sample[0]+".root is NOT at destination, proceed with Step3Clone + Skim..."
                f.write('Step3Clone step3_tmp.py '+sample[0]+' '+sample[1]+'\n')
            if readSkimsFromFile:
                f.write('Skim  skim_tmp.py '+sample[0]+' '+sample[1]+'\n')
            else:
                f.write('Skim  skim_tmp.py '+sample[0]+' '+sample[1]+' '+skim[0]+' "'+skim[1]+'"\n')
            f.close()
            jobName = 'job_'+sample[1]
            if not readSkimsFromFile:
                jobName = jobName+'_'+skim[0]

            os.system('chmod +x '+scriptName)
            submitToQueue = 'qsub -V -cwd -l h_vmem=6G -q all.q -N '+jobName+' '+scriptName 
            print submitToQueue
            os.system(submitToQueue)

        print "\n@@@@@ END JOBS @@@@@@@@@@@@@@@"

###########################################
###########################################


samples = [ 
    #['WW_TuneZ2star_8TeV_pythia6_tauola', 'WW'],
    #['WZ_TuneZ2star_8TeV_pythia6_tauola', 'WZ'],
    #['ZZ_TuneZ2star_8TeV_pythia6_tauola', 'ZZ'],
    #['DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball', 'DYJets'],
    #['T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola', 'TtW'],
    #['T_s-channel_TuneZ2star_8TeV-powheg-tauola', 'Ts'],
    #['T_t-channel_TuneZ2star_8TeV-powheg-tauola', 'Tt'],
    #['Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola', 'Tbars'],
    #['Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola', 'Tbart'],
    #['Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola', 'TbartW'],
    #['TTJets_HadronicMGDecays_8TeV-madgraph-part', 'TTJetsFullHad'],
    ['TTJets_FullLeptMGDecays_8TeV-madgraph-part', 'TTJetsFullLept'],
    #['TTJets_SemiLeptMGDecays_8TeV-madgraph-part', 'TTJetsSemiLept'],
    #['WJetsToLNu_PtW-100_TuneZ2star_8TeV-madgraph',     'WJets100'],
    #['WJetsToLNu_PtW-70To100_TuneZ2star_8TeV-madgraph', 'WJets70100'],
    #['ZH_ZToLL_HToBB_M-125_8TeV-powheg-herwigpp',  'ZH25'],
    #['WH_WToLNu_HToBB_M-125_8TeV-powheg-herwigpp', 'WH25'],
    #['SingleElectronRun2012AJul13EdmV42b',  'SingleElectron_1'],
    #['SingleElectronRun2012AAug06EdmV42',   'SingleElectron_2'],
    #['SingleElectronRun2012BJul13EdmV42',   'SingleElectron_3'],
    #['SingleElectronRun2012CAug24RerecoEdmV42', 'SingleElectron_4'],
    #['SingleElectronRun2012CPromptv2EdmV42', 'SingleElectron_5'],
    #['SingleElectronRun2012CPromptV2TopUpEdmV42', 'SingleElectron_6'],
    #['SingleMuRun2012AJul13EdmV42', 'SingleMu_1'],
    #['SingleMuRun2012AAug06EdmV42', 'SingleMu_2'],
    #['SingleMuRun2012BJul13EdmV42', 'SingleMu_3'],
    #['SingleMuRun2012CAug24RerecoEdmV42',  'SingleMu_4'],
    #['SingleMuRun2012CPromptv2EdmV42',     'SingleMu_5'],
    #['SingleMuRun2012CPromptV2TopUpEdmV42','SingleMu_6']
    #['DoubleElectron_Run2012A-13Jul2012-v1_ProcFIXED',   'DoubleEle_1'],
    #['DoubleElectron_Run2012A-recover-06Aug2012-v1_ProcV2','DoubleEle_2'],
    #['DoubleElectron_Run2012B-13Jul2012-v1_ProcFIXED','DoubleEle_3'],
    #['DoubleElectronRun2012CAug24RerecoEdmV42','DoubleEle_4'],
    #['DoubleElectron_Run2012C-PromptReco-v2_HBB_EDMNtupleV42_ProcV1','DoubleEle_5'],
    #['DoubleElectron_Run2012C-PromptReco-v2_HBB_EDMNtupleV42_ProcV2','DoubleEle_6'],
    #['TTH_HToBB_M-110_8TeV-pythia6','TTH110'],
    #['TTH_HToBB_M-115_8TeV-pythia6','TTH115'],
    #['TTH_HToBB_M-120_8TeV-pythia6','TTH120'],
    #['TTH_HToBB_M-125_8TeV-pythia6','TTH125'],
    #['TTH_HToBB_M-130_8TeV-pythia6','TTH130'],
    #['TTH_HToBB_M-135_8TeV-pythia6','TTH135'],
    #['TTWJets_8TeV-madgraph','TTW'],
    #['TTZJets_8TeV-madgraph','TTZ'],
    #['WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball','WJets'],
    #['DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph','DYJets10to50']
    ]

skims = [
    #['VType0','numJets30bTag>=2&&Vtype==0&&numJets30>=2&&jet1.pt>30&&jet2.pt>30'],
    ['VType1','numJets30bTag>=2&&Vtype==1&&numJets30>=2&&jet1.pt>30&&jet2.pt>30'],
    #['VType2','numJets30bTag>=2&&Vtype==2&&numJets30>=3&&jet1.pt>30&&jet2.pt>30&&jet3.pt>30'],
    #['VType3','numJets30bTag>=2&&Vtype==3&&numJets30>=3&&jet1.pt>30&&jet2.pt>30&&jet3.pt>30']
    ]
#skims = []

#processAll(samples)
processAllBatch(samples, skims, False)
#processAllDebugBatch(samples)
