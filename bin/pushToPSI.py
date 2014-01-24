#!/usr/bin/python
import time
import commands
import re
import os
import sys

outdir = 'bianchi@t3ui07.psi.ch:/scratch/bianchi/'

print "Start : %s" % time.ctime()

for i in range(1):

    print 'Itearation number %s' % i
    print '--------------------'

    onlyfiles = [ f for f in os.listdir('/Volumes/Macintosh HD/Data') if os.path.isfile(os.path.join('/Volumes/Macintosh HD/Data',f)) ]

    for f in onlyfiles:
        if re.search("hepmc2g", f)==None:
            continue
        try:
            with open( os.path.join('/Volumes/Macintosh HD/Data',f) ):
                print ''
                print 'File: ', f, ' opened'

                if 'END_EVENT_LISTING' in open(os.path.join('/Volumes/Macintosh HD/Data',f)).read():
                    print '..... file properly closed: PROCEED'                               
                    start = time.time()
                    copycommand = 'scp /Volumes/Macintosh\ HD/Data/'+f+' '+outdir+f
                    print copycommand
                    os.system('scp /Volumes/Macintosh\ HD/Data/'+f+' '+outdir+f)
                    #os.system('rm '+f)
                    elapsed = (time.time() - start)
                    print '..... Done in %s sec: now wait 0 more seconds' % elapsed
                    time.sleep( 0.01 )
                else:
                    print '..... file incomplete: WAIT FOR NEXT ITERATION'
                    
        except IOError:
            print 'Not found'

    print ''
    print '--------------------'
    print 'Iteration number %s done'


print "End : %s" % time.ctime()
