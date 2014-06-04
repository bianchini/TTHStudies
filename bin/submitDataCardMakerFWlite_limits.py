#!/usr/bin/env python

import FWCore.ParameterSet.Config as cms

import sys
sys.path.append('./')

from submitDataCardMakerFWlite import *


submitDataCardMakerFWlite_Limits("cat1_sb_L", "_rec_std",125)
submitDataCardMakerFWlite_Limits("cat2_sb_L", "_rec_std",125)
submitDataCardMakerFWlite_Limits("cat3_sb_L", "_rec_std",125)
submitDataCardMakerFWlite_Limits("cat6_sb_L", "_rec_std",125)
submitDataCardMakerFWlite_Limits("cat1_sb_H", "_rec_std",125)
submitDataCardMakerFWlite_Limits("cat2_sb_H", "_rec_std",125)
submitDataCardMakerFWlite_Limits("cat3_sb_H", "_rec_std",125)
submitDataCardMakerFWlite_Limits("cat6_sb_H", "_rec_std",125)

#submitDataCardMakerFWlite_Limits("cat1_sb")
#submitDataCardMakerFWlite_Limits("cat2_sb")
#submitDataCardMakerFWlite_Limits("cat3_sb")
#submitDataCardMakerFWlite_Limits("cat6_sb")

    
submitDataCardMakerFWlite_Limits("cat1_sb_nb", "_rec_std",125)
submitDataCardMakerFWlite_Limits("cat2_sb_nb", "_rec_std",125)
submitDataCardMakerFWlite_Limits("cat3_sb_nb", "_rec_std",125)
submitDataCardMakerFWlite_Limits("cat6_sb_nb", "_rec_std",125)
submitDataCardMakerFWlite_Limits("cat1_bj",    "_rec_std",125)
submitDataCardMakerFWlite_Limits("cat2_bj",    "_rec_std",125)
submitDataCardMakerFWlite_Limits("cat3_bj",    "_rec_std",125)
submitDataCardMakerFWlite_Limits("cat6_bj",    "_rec_std",125)

submitDataCardMakerFWlite_Limits("cat1_sb_L", "_MH90_rec_std",90)
submitDataCardMakerFWlite_Limits("cat2_sb_L", "_MH90_rec_std",90)
submitDataCardMakerFWlite_Limits("cat3_sb_L", "_MH90_rec_std",90)
submitDataCardMakerFWlite_Limits("cat6_sb_L", "_MH90_rec_std",90)
submitDataCardMakerFWlite_Limits("cat1_sb_H", "_MH90_rec_std",90)
submitDataCardMakerFWlite_Limits("cat2_sb_H", "_MH90_rec_std",90)
submitDataCardMakerFWlite_Limits("cat3_sb_H", "_MH90_rec_std",90)
submitDataCardMakerFWlite_Limits("cat6_sb_H", "_MH90_rec_std",90)

submitDataCardMakerFWlite_Limits("cat1_sb_L", "_MH160_rec_std",160)
submitDataCardMakerFWlite_Limits("cat2_sb_L", "_MH160_rec_std",160)
submitDataCardMakerFWlite_Limits("cat3_sb_L", "_MH160_rec_std",160)
submitDataCardMakerFWlite_Limits("cat6_sb_L", "_MH160_rec_std",160)
submitDataCardMakerFWlite_Limits("cat1_sb_H", "_MH160_rec_std",160)
submitDataCardMakerFWlite_Limits("cat2_sb_H", "_MH160_rec_std",160)
submitDataCardMakerFWlite_Limits("cat3_sb_H", "_MH160_rec_std",160)
submitDataCardMakerFWlite_Limits("cat6_sb_H", "_MH160_rec_std",160)
