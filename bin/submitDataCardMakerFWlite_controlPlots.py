#!/usr/bin/env python

import FWCore.ParameterSet.Config as cms

import sys
sys.path.append('./')

from submitDataCardMakerFWlite import *



###### w0 vs w1

binvec  = cms.vdouble()

binvec.append(0. )
for b in range(6):
    binvec.append(20. + 4*b)
binvecY = cms.vdouble()
for b in range(44):
    binvecY.append(0. + 1*b)    
submitDataCardMakerFWlite_all( "log(p_125_all_b):log(p_125_all_s)", "logPbvslogPs", cut_cat1_H , "cat1_H" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_b):log(p_125_all_s)", "logPbvslogPs", cut_cat1_L , "cat1_L" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_b):log(p_125_all_s)", "logPbvslogPs", cut_cat2_H , "cat2_H" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_b):log(p_125_all_s)", "logPbvslogPs", cut_cat2_L , "cat2_L" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_b):log(p_125_all_s)", "logPbvslogPs", cut_cat3_H , "cat3_H" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_b):log(p_125_all_s)", "logPbvslogPs", cut_cat3_L , "cat3_L" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_b):log(p_125_all_s)", "logPbvslogPs", cut_cat6_H , "cat6_H" , binvec, binvecY, 1)
submitDataCardMakerFWlite_all( "log(p_125_all_b):log(p_125_all_s)", "logPbvslogPs", cut_cat6_L , "cat6_L" , binvec, binvecY, 1)



###### Pbb vs Pjj

binvec = cms.vdouble()
binvec.append(-5.)
for b in range(3):
    binvec.append(2.5 + 2.5*b)
binvec.append(15.)
binvecY = cms.vdouble()
for b in range(56):
    binvecY.append(-10. + 0.5*b)
submitDataCardMakerFWlite_all( "log(p_125_all_b_ttbb/p_125_all_b):log(p_125_all_b_ttjj/p_125_all_b)", "logPbbvslogPjj", cut_cat1_H , "cat1_H" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_b_ttbb/p_125_all_b):log(p_125_all_b_ttjj/p_125_all_b)", "logPbbvslogPjj", cut_cat2_H , "cat2_H" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_b_ttbb/p_125_all_b):log(p_125_all_b_ttjj/p_125_all_b)", "logPbbvslogPjj", cut_cat3_H , "cat3_H" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_b_ttbb/p_125_all_b):log(p_125_all_b_ttjj/p_125_all_b)", "logPbbvslogPjj", cut_cat6_H , "cat6_H" , binvec, binvecY, 1)
binvec = cms.vdouble()
for b in range(14):
    binvec.append(-10. + 2*b)
submitDataCardMakerFWlite_all( "log(p_125_all_b_ttbb/p_125_all_b):log(p_125_all_b_ttjj/p_125_all_b)", "logPbbvslogPjj", cut_cat1_L , "cat1_L" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_b_ttbb/p_125_all_b):log(p_125_all_b_ttjj/p_125_all_b)", "logPbbvslogPjj", cut_cat2_L , "cat2_L" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_b_ttbb/p_125_all_b):log(p_125_all_b_ttjj/p_125_all_b)", "logPbbvslogPjj", cut_cat3_L , "cat3_L" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_b_ttbb/p_125_all_b):log(p_125_all_b_ttjj/p_125_all_b)", "logPbbvslogPjj", cut_cat6_L , "cat6_L" , binvec, binvecY, 1)



###### w0 vs Pbb

binvec  = cms.vdouble()
binvec.append(0.)
for b in range(6):
    binvec.append(15 + 4*b)
binvecY = cms.vdouble()
for b in range(56):
    binvecY.append(-10. + 0.5*b)
submitDataCardMakerFWlite_all( "log(p_125_all_s):log(p_125_all_s_ttbb/p_125_all_s)", "logPsvslogPbb", cut_cat1_H , "cat1_H" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_s):log(p_125_all_s_ttbb/p_125_all_s)", "logPsvslogPbb", cut_cat1_L , "cat1_L" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_s):log(p_125_all_s_ttbb/p_125_all_s)", "logPsvslogPbb", cut_cat2_H , "cat2_H" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_s):log(p_125_all_s_ttbb/p_125_all_s)", "logPsvslogPbb", cut_cat2_L , "cat2_L" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_s):log(p_125_all_s_ttbb/p_125_all_s)", "logPsvslogPbb", cut_cat3_H , "cat3_H" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_s):log(p_125_all_s_ttbb/p_125_all_s)", "logPsvslogPbb", cut_cat3_L , "cat3_L" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_s):log(p_125_all_s_ttbb/p_125_all_s)", "logPsvslogPbb", cut_cat6_H , "cat6_H" , binvec, binvecY, 1)
submitDataCardMakerFWlite_all( "log(p_125_all_s):log(p_125_all_s_ttbb/p_125_all_s)", "logPsvslogPbb", cut_cat6_L , "cat6_L" , binvec, binvecY, 1)


###### w0 vs Pjj

binvec  = cms.vdouble()
binvec.append(0.)
for b in range(6):
    binvec.append(15 + 4*b)
binvecY = cms.vdouble()
for b in range(56):
    binvecY.append(-10. + 0.5*b)
submitDataCardMakerFWlite_all( "log(p_125_all_s):log(p_125_all_b_ttjj/p_125_all_b)", "logPsvslogPjj", cut_cat1_H , "cat1_H" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_s):log(p_125_all_b_ttjj/p_125_all_b)", "logPsvslogPjj", cut_cat1_L , "cat1_L" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_s):log(p_125_all_b_ttjj/p_125_all_b)", "logPsvslogPjj", cut_cat2_H , "cat2_H" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_s):log(p_125_all_b_ttjj/p_125_all_b)", "logPsvslogPjj", cut_cat2_L , "cat2_L" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_s):log(p_125_all_b_ttjj/p_125_all_b)", "logPsvslogPjj", cut_cat3_H , "cat3_H" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_s):log(p_125_all_b_ttjj/p_125_all_b)", "logPsvslogPjj", cut_cat3_L , "cat3_L" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_s):log(p_125_all_b_ttjj/p_125_all_b)", "logPsvslogPjj", cut_cat6_H , "cat6_H" , binvec, binvecY, 1)
submitDataCardMakerFWlite_all( "log(p_125_all_s):log(p_125_all_b_ttjj/p_125_all_b)", "logPsvslogPjj", cut_cat6_L , "cat6_L" , binvec, binvecY, 1)



###### w1 vs Pbb


binvec  = cms.vdouble()
binvec.append(0.)
for b in range(6):
    binvec.append(15 + 4*b)
binvecY = cms.vdouble()
for b in range(56):
    binvecY.append(-10. + 0.5*b)
submitDataCardMakerFWlite_all( "log(p_125_all_b):log(p_125_all_b_ttbb/p_125_all_b)", "logPbvslogPbb", cut_cat1_H , "cat1_H" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_b):log(p_125_all_b_ttbb/p_125_all_b)", "logPbvslogPbb", cut_cat1_L , "cat1_L" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_b):log(p_125_all_b_ttbb/p_125_all_b)", "logPbvslogPbb", cut_cat2_H , "cat2_H" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_b):log(p_125_all_b_ttbb/p_125_all_b)", "logPbvslogPbb", cut_cat2_L , "cat2_L" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_b):log(p_125_all_b_ttbb/p_125_all_b)", "logPbvslogPbb", cut_cat3_H , "cat3_H" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_b):log(p_125_all_b_ttbb/p_125_all_b)", "logPbvslogPbb", cut_cat3_L , "cat3_L" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_b):log(p_125_all_b_ttbb/p_125_all_b)", "logPbvslogPbb", cut_cat6_H , "cat6_H" , binvec, binvecY, 1)
submitDataCardMakerFWlite_all( "log(p_125_all_b):log(p_125_all_b_ttbb/p_125_all_b)", "logPbvslogPbb", cut_cat6_L , "cat6_L" , binvec, binvecY, 1)



###### w1 vs Pjj

binvec  = cms.vdouble()
binvec.append(0.)
for b in range(6):
    binvec.append(15 + 4*b)
binvecY = cms.vdouble()
for b in range(56):
    binvecY.append(-10. + 0.5*b)
submitDataCardMakerFWlite_all( "log(p_125_all_b):log(p_125_all_b_ttjj/p_125_all_b)", "logPbvslogPjj", cut_cat1_H , "cat1_H" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_b):log(p_125_all_b_ttjj/p_125_all_b)", "logPbvslogPjj", cut_cat1_L , "cat1_L" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_b):log(p_125_all_b_ttjj/p_125_all_b)", "logPbvslogPjj", cut_cat2_H , "cat2_H" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_b):log(p_125_all_b_ttjj/p_125_all_b)", "logPbvslogPjj", cut_cat2_L , "cat2_L" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_b):log(p_125_all_b_ttjj/p_125_all_b)", "logPbvslogPjj", cut_cat3_H , "cat3_H" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_b):log(p_125_all_b_ttjj/p_125_all_b)", "logPbvslogPjj", cut_cat3_L , "cat3_L" , binvec, binvecY, 0)
submitDataCardMakerFWlite_all( "log(p_125_all_b):log(p_125_all_b_ttjj/p_125_all_b)", "logPbvslogPjj", cut_cat6_H , "cat6_H" , binvec, binvecY, 1)
submitDataCardMakerFWlite_all( "log(p_125_all_b):log(p_125_all_b_ttjj/p_125_all_b)", "logPbvslogPjj", cut_cat6_L , "cat6_L" , binvec, binvecY, 1)


    
binvec = cms.vdouble()
for b in range(17):
    binvec.append(50. + 8*b)
submitDataCardMakerFWlite_all( "best_0_mass_WHad", "best_0_mass_WHad",     cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_mass_WHad", "best_1_mass_WHad",     cut_cat1 , "cat1" , binvec, binvec, 0)

binvec = cms.vdouble()
for b in range(12):
    binvec.append(50. + 20*b)
submitDataCardMakerFWlite_all( "best_0_mass_TopHad", "best_0_mass_TopHad", cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_mass_TopHad", "best_1_mass_TopHad", cut_cat1 , "cat1" , binvec, binvec, 0)

binvec = cms.vdouble()
for b in range(12):
    binvec.append(50. + 20*b)
submitDataCardMakerFWlite_all( "best_0_mass_H", "best_0_mass_H", cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_mass_H", "best_0_mass_H", cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_mass_H", "best_0_mass_H", cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_mass_H", "best_0_mass_H", cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_1_mass_H", "best_1_mass_H", cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_mass_H", "best_1_mass_H", cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_mass_H", "best_1_mass_H", cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_mass_H", "best_1_mass_H", cut_cat6 , "cat6" , binvec, binvec, +1)

binvec = cms.vdouble()
for b in range(13):
    binvec.append(50. + 30*b)
submitDataCardMakerFWlite_all( "best_0_mass_TopLep1",   "best_0_mass_TopLep1",   cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_mass_TopLep1",   "best_0_mass_TopLep1",   cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_mass_TopLep1",   "best_0_mass_TopLep1",   cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_mass_TopLep1",   "best_0_mass_TopLep1",   cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_0_mass_TopLep2",   "best_0_mass_TopLep2",   cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_1_mass_TopLep1",   "best_1_mass_TopLep1",   cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_mass_TopLep1",   "best_1_mass_TopLep1",   cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_mass_TopLep1",   "best_1_mass_TopLep1",   cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_mass_TopLep1",   "best_1_mass_TopLep1",   cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_1_mass_TopLep2",   "best_1_mass_TopLep2",   cut_cat6 , "cat6" , binvec, binvec, +1)

binvec = cms.vdouble()
for b in range(12):
    binvec.append(0. + 20*b)
submitDataCardMakerFWlite_all( "best_0_e_MET",            "best_0_e_MET",            cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_e_MET",            "best_0_e_MET",            cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_e_MET",            "best_0_e_MET",            cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_e_MET",            "best_0_e_MET",            cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_1_e_MET",            "best_1_e_MET",            cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_e_MET",            "best_1_e_MET",            cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_e_MET",            "best_1_e_MET",            cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_e_MET",            "best_1_e_MET",            cut_cat6 , "cat6" , binvec, binvec, +1)

binvec = cms.vdouble()
for b in range(12):
    binvec.append(30. + 20*b)
submitDataCardMakerFWlite_all( "best_0_pt_WLep1",            "best_0_pt_WLep1",            cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_WLep1",            "best_0_pt_WLep1",            cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_WLep1",            "best_0_pt_WLep1",            cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_WLep1",            "best_0_pt_WLep1",            cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_0_pt_WLep2",            "best_0_pt_WLep2",            cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_0_pt_bLep",             "best_0_pt_bLep",             cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_bLep",             "best_0_pt_bLep",             cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_bLep",             "best_0_pt_bLep",             cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_bLep",             "best_0_pt_bLep",             cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_0_pt_WHad1",            "best_0_pt_WHad1",            cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_WHad1",            "best_0_pt_WHad1",            cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_WHad1",            "best_0_pt_WHad1",            cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_WHad1",            "best_0_pt_WHad1",            cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_0_pt_WHad2",            "best_0_pt_WHad2",            cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_WHad2",            "best_0_pt_WHad2",            cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_WHad2",            "best_0_pt_WHad2",            cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_WHad2",            "best_0_pt_WHad2",            cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_0_pt_bHad",             "best_0_pt_bHad",             cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_bHad",             "best_0_pt_bHad",             cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_bHad",             "best_0_pt_bHad",             cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_bHad",             "best_0_pt_bHad",             cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_0_pt_b1",               "best_0_pt_b1",             cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_b1",               "best_0_pt_b1",             cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_b1",               "best_0_pt_b1",             cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_b1",               "best_0_pt_b1",             cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_0_pt_b2",               "best_0_pt_b2",             cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_b2",               "best_0_pt_b2",             cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_b2",               "best_0_pt_b2",             cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_b2",               "best_0_pt_b2",             cut_cat6 , "cat6" , binvec, binvec, +1)

submitDataCardMakerFWlite_all( "best_1_pt_WLep1",            "best_1_pt_WLep1",            cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_WLep1",            "best_1_pt_WLep1",            cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_WLep1",            "best_1_pt_WLep1",            cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_WLep1",            "best_1_pt_WLep1",            cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_1_pt_WLep2",            "best_1_pt_WLep2",            cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_1_pt_bLep",             "best_1_pt_bLep",             cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_bLep",             "best_1_pt_bLep",             cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_bLep",             "best_1_pt_bLep",             cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_bLep",             "best_1_pt_bLep",             cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_1_pt_WHad1",            "best_1_pt_WHad1",            cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_WHad1",            "best_1_pt_WHad1",            cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_WHad1",            "best_1_pt_WHad1",            cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_WHad1",            "best_1_pt_WHad1",            cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_1_pt_WHad2",            "best_1_pt_WHad2",            cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_WHad2",            "best_1_pt_WHad2",            cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_WHad2",            "best_1_pt_WHad2",            cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_WHad2",            "best_1_pt_WHad2",            cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_1_pt_bHad",             "best_1_pt_bHad",             cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_bHad",             "best_1_pt_bHad",             cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_bHad",             "best_1_pt_bHad",             cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_bHad",             "best_1_pt_bHad",             cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_1_pt_b1",               "best_1_pt_b1",             cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_b1",               "best_1_pt_b1",             cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_b1",               "best_1_pt_b1",             cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_b1",               "best_1_pt_b1",             cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_1_pt_b2",               "best_1_pt_b2",             cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_b2",               "best_1_pt_b2",             cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_b2",               "best_1_pt_b2",             cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_b2",               "best_1_pt_b2",             cut_cat6 , "cat6" , binvec, binvec, +1)
 
binvec = cms.vdouble()
for b in range(11):
    binvec.append(-2.5 + 0.5*b)
submitDataCardMakerFWlite_all( "best_0_eta_WLep1",            "best_0_eta_WLep1",            cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_WLep1",            "best_0_eta_WLep1",            cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_WLep1",            "best_0_eta_WLep1",            cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_WLep1",            "best_0_eta_WLep1",            cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_0_eta_WLep2",            "best_0_eta_WLep2",            cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_0_eta_bLep",             "best_0_eta_bLep",             cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_bLep",             "best_0_eta_bLep",             cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_bLep",             "best_0_eta_bLep",             cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_bLep",             "best_0_eta_bLep",             cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_0_eta_WHad1",            "best_0_eta_WHad1",            cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_WHad1",            "best_0_eta_WHad1",            cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_WHad1",            "best_0_eta_WHad1",            cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_WHad1",            "best_0_eta_WHad1",            cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_0_eta_WHad2",            "best_0_eta_WHad2",            cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_WHad2",            "best_0_eta_WHad2",            cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_WHad2",            "best_0_eta_WHad2",            cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_WHad2",            "best_0_eta_WHad2",            cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_0_eta_bHad",             "best_0_eta_bHad",             cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_bHad",             "best_0_eta_bHad",             cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_bHad",             "best_0_eta_bHad",             cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_bHad",             "best_0_eta_bHad",             cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_0_eta_b1",               "best_0_eta_b1",               cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_b1",               "best_0_eta_b1",               cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_b1",               "best_0_eta_b1",               cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_b1",               "best_0_eta_b1",               cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_0_eta_b2",               "best_0_eta_b2",               cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_b2",               "best_0_eta_b2",               cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_b2",               "best_0_eta_b2",               cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_b2",               "best_0_eta_b2",               cut_cat6 , "cat6" , binvec, binvec, +1)

submitDataCardMakerFWlite_all( "best_1_eta_WLep1",            "best_1_eta_WLep1",            cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_WLep1",            "best_1_eta_WLep1",            cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_WLep1",            "best_1_eta_WLep1",            cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_WLep1",            "best_1_eta_WLep1",            cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_1_eta_WLep2",            "best_1_eta_WLep2",            cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_1_eta_bLep",             "best_1_eta_bLep",             cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_bLep",             "best_1_eta_bLep",             cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_bLep",             "best_1_eta_bLep",             cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_bLep",             "best_1_eta_bLep",             cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_1_eta_WHad1",            "best_1_eta_WHad1",            cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_WHad1",            "best_1_eta_WHad1",            cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_WHad1",            "best_1_eta_WHad1",            cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_WHad1",            "best_1_eta_WHad1",            cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_1_eta_WHad2",            "best_1_eta_WHad2",            cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_WHad2",            "best_1_eta_WHad2",            cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_WHad2",            "best_1_eta_WHad2",            cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_WHad2",            "best_1_eta_WHad2",            cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_1_eta_bHad",             "best_1_eta_bHad",             cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_bHad",             "best_1_eta_bHad",             cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_bHad",             "best_1_eta_bHad",             cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_bHad",             "best_1_eta_bHad",             cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_1_eta_b1",               "best_1_eta_b1",               cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_b1",               "best_1_eta_b1",               cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_b1",               "best_1_eta_b1",               cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_b1",               "best_1_eta_b1",               cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_1_eta_b2",               "best_1_eta_b2",               cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_b2",               "best_1_eta_b2",               cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_b2",               "best_1_eta_b2",               cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_b2",               "best_1_eta_b2",               cut_cat6 , "cat6" , binvec, binvec, +1)
 


binvec = cms.vdouble()
for b in range(11):
    binvec.append(0.3141592*b)
submitDataCardMakerFWlite_all( "best_0_dphi_MET_WLep1", "best_0_dphi_MET_WLep1", cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_dphi_MET_WLep1", "best_0_dphi_MET_WLep1", cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_dphi_MET_WLep1", "best_0_dphi_MET_WLep1", cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_dphi_MET_WLep1", "best_0_dphi_MET_WLep1", cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_0_dphi_MET_WLep2", "best_0_dphi_MET_WLep2", cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_0_dphi_b1_b2", "best_0_dphi_b1_b2", cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_dphi_b1_b2", "best_0_dphi_b1_b2", cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_dphi_b1_b2", "best_0_dphi_b1_b2", cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_0_dphi_b1_b2", "best_0_dphi_b1_b2", cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_1_dphi_MET_WLep1", "best_1_dphi_MET_WLep1", cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_dphi_MET_WLep1", "best_1_dphi_MET_WLep1", cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_dphi_MET_WLep1", "best_1_dphi_MET_WLep1", cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_dphi_MET_WLep1", "best_1_dphi_MET_WLep1", cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_1_dphi_MET_WLep2", "best_1_dphi_MET_WLep2", cut_cat6 , "cat6" , binvec, binvec, +1)
submitDataCardMakerFWlite_all( "best_1_dphi_b1_b2", "best_1_dphi_b1_b2", cut_cat1 , "cat1" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_dphi_b1_b2", "best_1_dphi_b1_b2", cut_cat2 , "cat2" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_dphi_b1_b2", "best_1_dphi_b1_b2", cut_cat3 , "cat3" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "best_1_dphi_b1_b2", "best_1_dphi_b1_b2", cut_cat6 , "cat6" , binvec, binvec, +1)

binvec = cms.vdouble()
for b in range(11):
    binvec.append(10. + 3*b)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b)", "logPb", cut_cat1_H , "cat1_H" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b)", "logPb", cut_cat1_L , "cat1_L" , binvec, binvec, 0)


binvec = cms.vdouble()
for b in range(11):
    binvec.append(20. + 3*b)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b)", "logPb", cut_cat2_H , "cat2_H" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b)", "logPb", cut_cat3_H , "cat3_H" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b)", "logPb", cut_cat6_H , "cat6_H" , binvec, binvec, 1)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b)", "logPb", cut_cat2_L , "cat2_L" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b)", "logPb", cut_cat3_L , "cat3_L" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b)", "logPb", cut_cat6_L , "cat6_L" , binvec, binvec, 1)

binvec = cms.vdouble()
for b in range(12):
    binvec.append(0. + 4*b)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_s)", "logPs", cut_cat1_H , "cat1_H" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_s)", "logPs", cut_cat2_H , "cat2_H" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_s)", "logPs", cut_cat3_H , "cat3_H" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_s)", "logPs", cut_cat6_H , "cat6_H" , binvec, binvec, 1)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_s)", "logPs", cut_cat1_L , "cat1_L" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_s)", "logPs", cut_cat2_L , "cat2_L" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_s)", "logPs", cut_cat3_L , "cat3_L" , binvec, binvec, 0)    
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_s)", "logPs", cut_cat6_L , "cat6_L" , binvec, binvec, 1)
    
    

binvec = cms.vdouble()
for b in range(13):
    binvec.append(-4. + 2*b)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b_ttbb/p_125_all_b)", "logPbb", cut_cat1_H , "cat1_H" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b_ttbb/p_125_all_b)", "logPbb", cut_cat2_H , "cat2_H" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b_ttbb/p_125_all_b)", "logPbb", cut_cat3_H , "cat3_H" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b_ttbb/p_125_all_b)", "logPbb", cut_cat6_H , "cat6_H" , binvec, binvec, 1)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b_ttbb/p_125_all_b)", "logPbb", cut_cat1_L , "cat1_L" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b_ttbb/p_125_all_b)", "logPbb", cut_cat2_L , "cat2_L" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b_ttbb/p_125_all_b)", "logPbb", cut_cat3_L , "cat3_L" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b_ttbb/p_125_all_b)", "logPbb", cut_cat6_L , "cat6_L" , binvec, binvec, 1)

binvec = cms.vdouble()
for b in range(11):
    binvec.append(-10. + 2*b)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b_ttjj/p_125_all_b)", "logPjj", cut_cat1_H , "cat1_H" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b_ttjj/p_125_all_b)", "logPjj", cut_cat2_H , "cat2_H" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b_ttjj/p_125_all_b)", "logPjj", cut_cat3_H , "cat3_H" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b_ttjj/p_125_all_b)", "logPjj", cut_cat6_H , "cat6_H" , binvec, binvec, 1)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b_ttjj/p_125_all_b)", "logPjj", cut_cat1_L , "cat1_L" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b_ttjj/p_125_all_b)", "logPjj", cut_cat2_L , "cat2_L" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b_ttjj/p_125_all_b)", "logPjj", cut_cat3_L , "cat3_L" , binvec, binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b_ttjj/p_125_all_b)", "logPjj", cut_cat6_L , "cat6_L" , binvec, binvec, 1)

    


