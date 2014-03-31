#!/usr/bin/env python

import os

for queue in range(60):
    print "Queue %s" % queue
    os.system('qrsh -q debug.q -l hostname=t3wn'+str(queue)+' ls /scratch/$USER')
    os.system('qrsh -q debug.q -l hostname=t3wn'+str(queue)+' find /scratch/$USER/ -user $USER -exec rm -rf {} \;')

