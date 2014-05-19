#!/usr/bin/env python

import FWCore.ParameterSet.Config as cms

import sys
sys.path.append('./')

from submitDataCardMakerFWlite import *


submitDataCardMakerFWlite_Limits("cat1_sb_L")
submitDataCardMakerFWlite_Limits("cat2_sb_L")
submitDataCardMakerFWlite_Limits("cat3_sb_L")
submitDataCardMakerFWlite_Limits("cat6_sb_L")
submitDataCardMakerFWlite_Limits("cat1_sb_H")
submitDataCardMakerFWlite_Limits("cat2_sb_H")
submitDataCardMakerFWlite_Limits("cat3_sb_H")
submitDataCardMakerFWlite_Limits("cat6_sb_H")

submitDataCardMakerFWlite_Limits("cat1_sb")
submitDataCardMakerFWlite_Limits("cat2_sb")
submitDataCardMakerFWlite_Limits("cat3_sb")
submitDataCardMakerFWlite_Limits("cat6_sb")
    
submitDataCardMakerFWlite_Limits("cat1_sb_nb")
submitDataCardMakerFWlite_Limits("cat2_sb_nb")
submitDataCardMakerFWlite_Limits("cat3_sb_nb")
submitDataCardMakerFWlite_Limits("cat6_sb_nb")
submitDataCardMakerFWlite_Limits("cat1_bj")
submitDataCardMakerFWlite_Limits("cat2_bj")
submitDataCardMakerFWlite_Limits("cat3_bj")
submitDataCardMakerFWlite_Limits("cat6_bj")
