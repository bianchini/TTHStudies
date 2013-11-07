#!/usr/bin/env python


import commands
import re
import os

import sys
sys.path.append('./')


def createList(sampleName, path, step, split):

     out = open('fileListToCopy_'+sampleName+'_'+str(split[0])+'-'+str(split[1])+'.txt','w')

     counter = 1
     for k in range(10):
          output = commands.getoutput("lcg-ls -o "+str(k*step)+" -c "+str(step)+" -b -D srmv2 "+path)
          outFiles = re.split(r'\n',output)
          for name in outFiles:
               if re.search(".root",name)!=None:
                    #foo = name.replace("/pnfs/lcg.cscs.ch/cms/trivcat","")
                    foo = name.replace("/pnfs/psi.ch/cms/trivcat", "")
                    if counter>=split[0] and counter<=split[1]:
                         out.write(foo+'\n')
                    counter = counter + 1
     out.close()


#################################################
#################################################


def useDataReplica(sampleName, dest, split):
    print "data_replica.py  --from-site  T2_CH_CSCS --to-site T3_CH_PSI text/fileListToCopy_"+sampleName+'_'+str(split[0])+'-'+str(split[1])+".txt "+dest
    os.system("data_replica.py  --from-site  T2_CH_CSCS --to-site T3_CH_PSI text/fileListToCopy_"+sampleName+'_'+str(split[0])+'-'+str(split[1])+".txt "+dest)


#################################################
#################################################

def createFileListAndCopy(sampleName, path, step, split, dest):

     out = open('fileListToCopy_'+sampleName+'_'+str(split[0])+'-'+str(split[1])+'.txt','w')

     counter = 0
     counted = 0;
     for k in range(10):
          output = commands.getoutput("lcg-ls -o "+str(k*step)+" -c "+str(step)+" -b -D srmv2 "+path)
          outFiles = re.split(r'\n',output)
          for name in outFiles:
               if re.search(".root",name)!=None:
                    counter = counter + 1
                    if not(counter>=split[0] and counter<=split[1]):
                         continue

                    toBeMasked = ""
                    skip = False
                    words = re.split(r'/',name)
                    for word in words:
                         if re.search(".root",word)!=None:
                              toBeMasked = word
                              #print "Trying to match "+toBeMasked
                              maskFile = open('fileListToMask_'+sampleName+'.txt','r')
                              for line in maskFile:
                                   #print line
                                   words2 = re.split(r'/',line)
                                   #print words2
                                   for word2 in words2:
                                        if re.search(".root",word2)!=None:
                                             #print "Comparing with "+word2
                                             if re.search(toBeMasked,word2)!=None:
                                                  print word2+" will be masked"
                                                  skip = True
                         
                    foo = name.replace("/pnfs/lcg.cscs.ch/cms/trivcat","")
                    #foo = name.replace("/pnfs/psi.ch/cms/trivcat", "")

                    if (not skip) and counter>=split[0] and counter<=split[1]:
                         out.write(foo+'\n')
                         counted = counted + 1

     out.close()

     if counted<1:
          print "Not submitting this job..."
          return 

     jobName    = 'dr_'+sampleName+'_'+str(split[0])+'-'+str(split[1])
     scriptName = 'batchDatReplica_'+sampleName+'_'+str(split[0])+'-'+str(split[1])+'.sh'
    
     f = open(scriptName,'w')
     f.write('#!/bin/bash\n\n')
     f.write('cd /shome/bianchi/CMSSW_5_3_3_patch2/src/Bianchi/TTHStudies/test\n')
     f.write('source /swshare/psit3/etc/profile.d/cms_ui_env.sh\n')
     f.write('export SCRAM_ARCH="slc5_amd64_gcc462"\n')
     f.write('source $VO_CMS_SW_DIR/cmsset_default.sh\n')
     f.write('eval `scramv1 runtime -sh`\n')
     f.write('export LD_PRELOAD="libglobus_gssapi_gsi_gcc64pthr.so.0":${LD_PRELOAD}\n')
     f.write('export LD_LIBRARY_PATH=/swshare/glite/globus/lib/:/swshare/glite/d-cache/dcap/lib64/:$LD_LIBRARY_PATH\n')
     f.write('\n\n')
     f.write('\n\n')
     f.write("data_replica.py  --from-site  T2_CH_CSCS --to-site T3_CH_PSI fileListToCopy_"+sampleName+'_'+str(split[0])+'-'+str(split[1])+".txt "+dest)
     f.close()

     os.system('chmod +x '+scriptName)
     submitToQueue = 'qsub -V -cwd -l h_vmem=2G -q all.q -N '+jobName+' '+scriptName 
     print submitToQueue
     os.system(submitToQueue)



#################################################
#################################################




# T2

#createList("DYJets10-50","srm://storage01.lcg.cscs.ch:8443/srm/managerv2?SFN=/pnfs/lcg.cscs.ch/cms/trivcat/store/user/bianchi/VHbbNTuples/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph", 500, [   1, 500])
#createList("DYJets10-50","srm://storage01.lcg.cscs.ch:8443/srm/managerv2?SFN=/pnfs/lcg.cscs.ch/cms/trivcat/store/user/bianchi/VHbbNTuples/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph", 500, [ 501,1000])
#createList("DYJets10-50","srm://storage01.lcg.cscs.ch:8443/srm/managerv2?SFN=/pnfs/lcg.cscs.ch/cms/trivcat/store/user/bianchi/VHbbNTuples/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph", 500, [1001,2000])


#createFileListAndCopy("DYJets10to50","srm://storage01.lcg.cscs.ch:8443/srm/managerv2?SFN=/pnfs/lcg.cscs.ch/cms/trivcat/store/user/bianchi/VHbbNTuples/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph", 500, [   1, 200], "/store/user/bianchi//HBB_EDMNtuple/AllHDiJetPt_Step1/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph")
#createFileListAndCopy("DYJets10to50","srm://storage01.lcg.cscs.ch:8443/srm/managerv2?SFN=/pnfs/lcg.cscs.ch/cms/trivcat/store/user/bianchi/VHbbNTuples/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph", 500, [ 201, 400], "/store/user/bianchi//HBB_EDMNtuple/AllHDiJetPt_Step1/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph")
#createFileListAndCopy("DYJets10to50","srm://storage01.lcg.cscs.ch:8443/srm/managerv2?SFN=/pnfs/lcg.cscs.ch/cms/trivcat/store/user/bianchi/VHbbNTuples/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph", 500, [ 401, 600], "/store/user/bianchi//HBB_EDMNtuple/AllHDiJetPt_Step1/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph")
#createFileListAndCopy("DYJets10to50","srm://storage01.lcg.cscs.ch:8443/srm/managerv2?SFN=/pnfs/lcg.cscs.ch/cms/trivcat/store/user/bianchi/VHbbNTuples/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph", 500, [ 601, 800], "/store/user/bianchi//HBB_EDMNtuple/AllHDiJetPt_Step1/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph")
#createFileListAndCopy("DYJets10to50","srm://storage01.lcg.cscs.ch:8443/srm/managerv2?SFN=/pnfs/lcg.cscs.ch/cms/trivcat/store/user/bianchi/VHbbNTuples/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph", 500, [ 801,1000], "/store/user/bianchi//HBB_EDMNtuple/AllHDiJetPt_Step1/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph")
#createFileListAndCopy("DYJets10to50","srm://storage01.lcg.cscs.ch:8443/srm/managerv2?SFN=/pnfs/lcg.cscs.ch/cms/trivcat/store/user/bianchi/VHbbNTuples/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph", 500, [1001,1200], "/store/user/bianchi//HBB_EDMNtuple/AllHDiJetPt_Step1/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph")
#createFileListAndCopy("DYJets10to50","srm://storage01.lcg.cscs.ch:8443/srm/managerv2?SFN=/pnfs/lcg.cscs.ch/cms/trivcat/store/user/bianchi/VHbbNTuples/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph", 500, [1201,1400], "/store/user/bianchi//HBB_EDMNtuple/AllHDiJetPt_Step1/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph")
#createFileListAndCopy("DYJets10to50","srm://storage01.lcg.cscs.ch:8443/srm/managerv2?SFN=/pnfs/lcg.cscs.ch/cms/trivcat/store/user/bianchi/VHbbNTuples/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph", 500, [1401,1600], "/store/user/bianchi//HBB_EDMNtuple/AllHDiJetPt_Step1/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph")
#createFileListAndCopy("DYJets10to50","srm://storage01.lcg.cscs.ch:8443/srm/managerv2?SFN=/pnfs/lcg.cscs.ch/cms/trivcat/store/user/bianchi/VHbbNTuples/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph", 500, [1601,1800], "/store/user/bianchi//HBB_EDMNtuple/AllHDiJetPt_Step1/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph")
#createFileListAndCopy("DYJets10to50","srm://storage01.lcg.cscs.ch:8443/srm/managerv2?SFN=/pnfs/lcg.cscs.ch/cms/trivcat/store/user/bianchi/VHbbNTuples/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph", 500, [1801,2000], "/store/user/bianchi//HBB_EDMNtuple/AllHDiJetPt_Step1/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph")
#createFileListAndCopy("DYJets10to50","srm://storage01.lcg.cscs.ch:8443/srm/managerv2?SFN=/pnfs/lcg.cscs.ch/cms/trivcat/store/user/bianchi/VHbbNTuples/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph", 500, [2001,2200], "/store/user/bianchi//HBB_EDMNtuple/AllHDiJetPt_Step1/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph")
#createFileListAndCopy("DYJets10to50","srm://storage01.lcg.cscs.ch:8443/srm/managerv2?SFN=/pnfs/lcg.cscs.ch/cms/trivcat/store/user/bianchi/VHbbNTuples/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph", 500, [2201,2400], "/store/user/bianchi//HBB_EDMNtuple/AllHDiJetPt_Step1/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph")
#createFileListAndCopy("DYJets10to50","srm://storage01.lcg.cscs.ch:8443/srm/managerv2?SFN=/pnfs/lcg.cscs.ch/cms/trivcat/store/user/bianchi/VHbbNTuples/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph", 500, [2401,2600], "/store/user/bianchi//HBB_EDMNtuple/AllHDiJetPt_Step1/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph")
#createFileListAndCopy("DYJets10to50","srm://storage01.lcg.cscs.ch:8443/srm/managerv2?SFN=/pnfs/lcg.cscs.ch/cms/trivcat/store/user/bianchi/VHbbNTuples/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph", 500, [2601,2800], "/store/user/bianchi//HBB_EDMNtuple/AllHDiJetPt_Step1/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph")


# T3
#createList("WJets", "srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/bianchi//HBB_EDMNtuple/AllHDiJetPt_Step1/WJetsToLNu_TuneZ2Star_8TeV-madgraph-tarball",500)
createList("DYJets10to50", "srm://t3se01.psi.ch:8443/srm/managerv2?SFN=/pnfs/psi.ch/cms/trivcat/store/user/bianchi//HBB_EDMNtuple/AllHDiJetPt_Step1/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph",500, [1, 5000])



#useDataReplica("DYJets10-50", "/store/user/bianchi//HBB_EDMNtuple/AllHDiJetPt_Step1/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph", [ 1, 500])
#useDataReplica("DYJets10-50", "/store/user/bianchi//HBB_EDMNtuple/AllHDiJetPt_Step1/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph", [ 501, 1000])
#useDataReplica("DYJets10-50", "/store/user/bianchi//HBB_EDMNtuple/AllHDiJetPt_Step1/DYJetsToLL_M-10To50_TuneZ2Star_8TeV-madgraph", [ 1001, 2000])
