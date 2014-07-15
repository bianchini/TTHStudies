#!/usr/bin/env python


import commands
import re
import os

SHIFT=2
extra='MH90_nC_'

file = open('log.txt', 'r')

for line in file.readlines():
    if re.search("MEAnalysis", line)!=None:
        split = line.split("_")
        #print split
        sample = split[4+SHIFT]
        job    = (split[5+SHIFT].split("."))[0]
        print sample,job
        command = 'qsub -q main job_all_'+extra+'rec_std_'+sample+'_'+job+'.sh'
        print command
        os.system(command)
