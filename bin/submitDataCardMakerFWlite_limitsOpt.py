#!/usr/bin/env python

import FWCore.ParameterSet.Config as cms

import sys
sys.path.append('./')

from submitDataCardMakerFWlite import *


    
########################################### optimize btagLR lower cut

cuts = [0.85, 0.875, 0.90, 0.925 , 0.950, 0.975, 0.980, 0.990 ]

trial = 0
for cut in cuts:
     #   submitDataCardMakerFWlite_Limits_Optimization("cat1_sb",  ("btag_LR>=%f" % cut), "cat1_"+str(trial) )
     #   submitDataCardMakerFWlite_Limits_Optimization("cat2_sb",  ("btag_LR>=%f" % cut), "cat2_"+str(trial) )
     #   submitDataCardMakerFWlite_Limits_Optimization("cat3_sb",  ("btag_LR>=%f" % cut), "cat3_"+str(trial) )
     #   submitDataCardMakerFWlite_Limits_Optimization("cat6_sb",  ("btag_LR>=%f" % cut), "cat6_"+str(trial) )
     trial += 1


########################################### optimize btagLR splitting


cuts =  [0.80, 0.825, 0.850, 0.875, 0.900]

trial = 0
for cut in range(len(cuts)):
    for trial in range(2):
        if trial==0:
            #  submitDataCardMakerFWlite_Limits_Optimization("cat1_sb",  ("btag_LR>=%f" % cuts[cut] ), "cat1-"+str(cut)+"_"+str(trial) )
            #  submitDataCardMakerFWlite_Limits_Optimization("cat2_sb",  ("btag_LR>=%f" % cuts[cut] ), "cat2-"+str(cut)+"_"+str(trial) )
            #  submitDataCardMakerFWlite_Limits_Optimization("cat3_sb",  ("btag_LR>=%f" % cuts[cut] ), "cat3-"+str(cut)+"_"+str(trial) )
            #  submitDataCardMakerFWlite_Limits_Optimization("cat6_sb",  ("btag_LR>=%f" % cuts[cut] ), "cat6-"+str(cut)+"_"+str(trial) )
            
            #  submitDataCardMakerFWlite_Limits_Optimization("cat1_sb",  ("btag_LR>=%f" % 0.995 ), "cat1-"+str(cut)+"_"+str(trial) )
            #  submitDataCardMakerFWlite_Limits_Optimization("cat2_sb",  ("btag_LR>=%f" % 0.9925), "cat2-"+str(cut)+"_"+str(trial) )
            #  submitDataCardMakerFWlite_Limits_Optimization("cat3_sb",  ("btag_LR>=%f" % 0.995),  "cat3-"+str(cut)+"_"+str(trial) )
            #  submitDataCardMakerFWlite_Limits_Optimization("cat6_sb",  ("btag_LR>=%f" % 0.925),  "cat6-"+str(cut)+"_"+str(trial) )
            trial += 0
        else:
            #  submitDataCardMakerFWlite_Limits_Optimization("cat1_sb",  ("btag_LR<%f"  % cuts[cut] ), "cat1-"+str(cut)+"_"+str(trial) )
            #  submitDataCardMakerFWlite_Limits_Optimization("cat2_sb",  ("btag_LR<%f"  % cuts[cut] ), "cat2-"+str(cut)+"_"+str(trial) )
            #  submitDataCardMakerFWlite_Limits_Optimization("cat3_sb",  ("btag_LR<%f"  % cuts[cut] ), "cat3-"+str(cut)+"_"+str(trial) )
            #  submitDataCardMakerFWlite_Limits_Optimization("cat6_sb",  ("btag_LR<%f"  % cuts[cut] ), "cat6-"+str(cut)+"_"+str(trial) )
            
            #  submitDataCardMakerFWlite_Limits_Optimization("cat1_sb",  ("btag_LR<%f && btag_LR>=%f" % (0.995, cuts[cut])  ), "cat1-"+str(cut)+"_"+str(trial) )
            #  submitDataCardMakerFWlite_Limits_Optimization("cat2_sb",  ("btag_LR<%f && btag_LR>=%f" % (0.9925,cuts[cut])  ), "cat2-"+str(cut)+"_"+str(trial) )
            #  submitDataCardMakerFWlite_Limits_Optimization("cat3_sb",  ("btag_LR<%f && btag_LR>=%f" % (0.995, cuts[cut])  ), "cat3-"+str(cut)+"_"+str(trial) )
            #  submitDataCardMakerFWlite_Limits_Optimization("cat6_sb",  ("btag_LR<%f && btag_LR>=%f" % (0.925, cuts[cut])  ), "cat6-"+str(cut)+"_"+str(trial) )
            trial += 0
            
            
########################################### optimize S/B normalization

cuts =  [1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.8]
    
trial = 0
for cut in cuts:
    #submitDataCardMakerFWlite_Limits_Optimization("cat1_sb_L",  ("btag_LR>=%f" % 0.), "cat1_"+str(trial), cut )
    #submitDataCardMakerFWlite_Limits_Optimization("cat2_sb_L",  ("btag_LR>=%f" % 0.), "cat2_"+str(trial), cut )
    #submitDataCardMakerFWlite_Limits_Optimization("cat3_sb_L",  ("btag_LR>=%f" % 0.), "cat3_"+str(trial), cut )
    #submitDataCardMakerFWlite_Limits_Optimization("cat6_sb_L",  ("btag_LR<%f" % 0.925), "cat6_"+str(trial), cut )
    trial += 1
    


