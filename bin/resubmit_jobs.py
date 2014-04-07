#!/usr/bin/env python


import commands
import re
import os


file = open('log.txt')

for line in file.readlines():
    if re.search("MEAnalysis", line)!=None:
        split = line.split("_")
        sample = split[4]
        job    = (split[5].split("."))[0]
        print sample,job
        command = 'qsub -q main job_all_rec_std_'+sample+'_'+job+'.sh'
        print command
        os.system(command)
