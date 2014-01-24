#!/usr/bin/python
import time
import commands
import re
import os
import sys

outdir = '/pnfs/psi.ch/cms/trivcat/store/user/bianchi/HepMC/SherpaOpenLoops/Jan24_2014/default_1/'

print "Start : %s" % time.ctime()

# how many files have been copied
counter        = 0

# number of hepmc files in the PWD
numoflastfiles = 0

# wait time before new cicle (subject to feedback)
penalty = 2 #sec

# number of times no hepmc files are found (the source hangs, or it is over)
numberoffailures = 0

for i in range(2000):

    print 'Itearation number %s' % i
    print '--------------------'

    numoffiles     = 0
    
    onlyfiles = [ f for f in os.listdir('./') if os.path.isfile(os.path.join('./',f)) ]

    for f in onlyfiles:
        if re.search("hepmc2g", f)==None:
            continue
        try:
            with open( f ):
                print ''
                print 'File: ', f, ' opened'            
                numoffiles += 1
                if 'END_EVENT_LISTING' in open(f).read():
                    print '..... file properly closed: PROCEED'
                    if os.path.isfile( os.path.join( outdir,f) ):
                        print "..... is already at destination: SKIP"
                    else:
                        print "..... is not at destination: COPY"
                        start = time.time()
                        #time.sleep(2.3)
                        os.system('lcg-cp -b -D srmv2 '+f+' srm://t3se01.psi.ch:8443/srm/managerv2?SFN='+outdir+f)
                        counter += 1
                        os.system('rm '+f)
                        elapsed = (time.time() - start)
                        print '..... Done in %s sec' % elapsed
                else:
                    print '..... file incomplete: WAIT FOR NEXT ITERATION'
                    
        except IOError:
            print 'Not found'

    print ''
    print '--------------------'

    if numoffiles==0:
        numberoffailures += 1
    if numberoffailures>10:
        break

    # feedback
    if  numoffiles > numoflastfiles:
        penalty  -= 0.2
    elif  numoffiles < numoflastfiles:
        penalty  += 0.2
        
    if penalty <= 0.:        
        penalty = 0.01

    numoflastfiles = numoffiles

    print 'Iteration number %s done: wait for %s sec...' % (i,penalty)
    time.sleep( penalty )

print "Copied %s files" % counter
print "End : %s" % time.ctime()
