#!/usr/bin/env python

import FWCore.ParameterSet.Config as cms

import sys
sys.path.append('./')

from submitDataCardMakerFWlite import *


binvec = cms.vdouble()
for b in range(17):
    binvec.append(50. + 8*b)
submitDataCardMakerFWlite_all( "best_0_mass_WHad", "best_0_mass_WHad",     cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_mass_WHad", "best_1_mass_WHad",     cut_cat1 , "cat1" , binvec, 0)

binvec = cms.vdouble()
for b in range(12):
    binvec.append(50. + 20*b)
submitDataCardMakerFWlite_all( "best_0_mass_TopHad", "best_0_mass_TopHad", cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_mass_TopHad", "best_1_mass_TopHad", cut_cat1 , "cat1" , binvec, 0)

binvec = cms.vdouble()
for b in range(12):
    binvec.append(50. + 20*b)
submitDataCardMakerFWlite_all( "best_0_mass_H", "best_0_mass_H", cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_mass_H", "best_0_mass_H", cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_mass_H", "best_0_mass_H", cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_mass_H", "best_0_mass_H", cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_1_mass_H", "best_1_mass_H", cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_mass_H", "best_1_mass_H", cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_mass_H", "best_1_mass_H", cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_mass_H", "best_1_mass_H", cut_cat6 , "cat6" , binvec, +1)

binvec = cms.vdouble()
for b in range(13):
    binvec.append(50. + 30*b)
submitDataCardMakerFWlite_all( "best_0_mass_TopLep1",   "best_0_mass_TopLep1",   cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_mass_TopLep1",   "best_0_mass_TopLep1",   cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_mass_TopLep1",   "best_0_mass_TopLep1",   cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_mass_TopLep1",   "best_0_mass_TopLep1",   cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_0_mass_TopLep2",   "best_0_mass_TopLep2",   cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_1_mass_TopLep1",   "best_1_mass_TopLep1",   cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_mass_TopLep1",   "best_1_mass_TopLep1",   cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_mass_TopLep1",   "best_1_mass_TopLep1",   cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_mass_TopLep1",   "best_1_mass_TopLep1",   cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_1_mass_TopLep2",   "best_1_mass_TopLep2",   cut_cat6 , "cat6" , binvec, +1)

binvec = cms.vdouble()
for b in range(12):
    binvec.append(0. + 20*b)
submitDataCardMakerFWlite_all( "best_0_e_MET",            "best_0_e_MET",            cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_e_MET",            "best_0_e_MET",            cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_e_MET",            "best_0_e_MET",            cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_e_MET",            "best_0_e_MET",            cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_1_e_MET",            "best_1_e_MET",            cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_e_MET",            "best_1_e_MET",            cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_e_MET",            "best_1_e_MET",            cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_e_MET",            "best_1_e_MET",            cut_cat6 , "cat6" , binvec, +1)

binvec = cms.vdouble()
for b in range(12):
    binvec.append(30. + 20*b)
submitDataCardMakerFWlite_all( "best_0_pt_WLep1",            "best_0_pt_WLep1",            cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_WLep1",            "best_0_pt_WLep1",            cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_WLep1",            "best_0_pt_WLep1",            cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_WLep1",            "best_0_pt_WLep1",            cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_0_pt_WLep2",            "best_0_pt_WLep2",            cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_0_pt_bLep",             "best_0_pt_bLep",             cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_bLep",             "best_0_pt_bLep",             cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_bLep",             "best_0_pt_bLep",             cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_bLep",             "best_0_pt_bLep",             cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_0_pt_WHad1",            "best_0_pt_WHad1",            cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_WHad1",            "best_0_pt_WHad1",            cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_WHad1",            "best_0_pt_WHad1",            cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_WHad1",            "best_0_pt_WHad1",            cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_0_pt_WHad2",            "best_0_pt_WHad2",            cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_WHad2",            "best_0_pt_WHad2",            cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_WHad2",            "best_0_pt_WHad2",            cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_WHad2",            "best_0_pt_WHad2",            cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_0_pt_bHad",             "best_0_pt_bHad",             cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_bHad",             "best_0_pt_bHad",             cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_bHad",             "best_0_pt_bHad",             cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_bHad",             "best_0_pt_bHad",             cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_0_pt_b1",               "best_0_pt_b1",             cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_b1",               "best_0_pt_b1",             cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_b1",               "best_0_pt_b1",             cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_b1",               "best_0_pt_b1",             cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_0_pt_b2",               "best_0_pt_b2",             cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_b2",               "best_0_pt_b2",             cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_b2",               "best_0_pt_b2",             cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_pt_b2",               "best_0_pt_b2",             cut_cat6 , "cat6" , binvec, +1)

submitDataCardMakerFWlite_all( "best_1_pt_WLep1",            "best_1_pt_WLep1",            cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_WLep1",            "best_1_pt_WLep1",            cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_WLep1",            "best_1_pt_WLep1",            cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_WLep1",            "best_1_pt_WLep1",            cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_1_pt_WLep2",            "best_1_pt_WLep2",            cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_1_pt_bLep",             "best_1_pt_bLep",             cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_bLep",             "best_1_pt_bLep",             cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_bLep",             "best_1_pt_bLep",             cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_bLep",             "best_1_pt_bLep",             cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_1_pt_WHad1",            "best_1_pt_WHad1",            cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_WHad1",            "best_1_pt_WHad1",            cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_WHad1",            "best_1_pt_WHad1",            cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_WHad1",            "best_1_pt_WHad1",            cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_1_pt_WHad2",            "best_1_pt_WHad2",            cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_WHad2",            "best_1_pt_WHad2",            cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_WHad2",            "best_1_pt_WHad2",            cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_WHad2",            "best_1_pt_WHad2",            cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_1_pt_bHad",             "best_1_pt_bHad",             cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_bHad",             "best_1_pt_bHad",             cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_bHad",             "best_1_pt_bHad",             cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_bHad",             "best_1_pt_bHad",             cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_1_pt_b1",               "best_1_pt_b1",             cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_b1",               "best_1_pt_b1",             cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_b1",               "best_1_pt_b1",             cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_b1",               "best_1_pt_b1",             cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_1_pt_b2",               "best_1_pt_b2",             cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_b2",               "best_1_pt_b2",             cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_b2",               "best_1_pt_b2",             cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_pt_b2",               "best_1_pt_b2",             cut_cat6 , "cat6" , binvec, +1)
 
binvec = cms.vdouble()
for b in range(11):
    binvec.append(-2.5 + 0.5*b)
submitDataCardMakerFWlite_all( "best_0_eta_WLep1",            "best_0_eta_WLep1",            cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_WLep1",            "best_0_eta_WLep1",            cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_WLep1",            "best_0_eta_WLep1",            cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_WLep1",            "best_0_eta_WLep1",            cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_0_eta_WLep2",            "best_0_eta_WLep2",            cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_0_eta_bLep",             "best_0_eta_bLep",             cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_bLep",             "best_0_eta_bLep",             cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_bLep",             "best_0_eta_bLep",             cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_bLep",             "best_0_eta_bLep",             cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_0_eta_WHad1",            "best_0_eta_WHad1",            cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_WHad1",            "best_0_eta_WHad1",            cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_WHad1",            "best_0_eta_WHad1",            cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_WHad1",            "best_0_eta_WHad1",            cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_0_eta_WHad2",            "best_0_eta_WHad2",            cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_WHad2",            "best_0_eta_WHad2",            cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_WHad2",            "best_0_eta_WHad2",            cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_WHad2",            "best_0_eta_WHad2",            cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_0_eta_bHad",             "best_0_eta_bHad",             cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_bHad",             "best_0_eta_bHad",             cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_bHad",             "best_0_eta_bHad",             cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_bHad",             "best_0_eta_bHad",             cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_0_eta_b1",               "best_0_eta_b1",               cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_b1",               "best_0_eta_b1",               cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_b1",               "best_0_eta_b1",               cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_b1",               "best_0_eta_b1",               cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_0_eta_b2",               "best_0_eta_b2",               cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_b2",               "best_0_eta_b2",               cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_b2",               "best_0_eta_b2",               cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_eta_b2",               "best_0_eta_b2",               cut_cat6 , "cat6" , binvec, +1)

submitDataCardMakerFWlite_all( "best_1_eta_WLep1",            "best_1_eta_WLep1",            cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_WLep1",            "best_1_eta_WLep1",            cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_WLep1",            "best_1_eta_WLep1",            cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_WLep1",            "best_1_eta_WLep1",            cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_1_eta_WLep2",            "best_1_eta_WLep2",            cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_1_eta_bLep",             "best_1_eta_bLep",             cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_bLep",             "best_1_eta_bLep",             cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_bLep",             "best_1_eta_bLep",             cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_bLep",             "best_1_eta_bLep",             cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_1_eta_WHad1",            "best_1_eta_WHad1",            cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_WHad1",            "best_1_eta_WHad1",            cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_WHad1",            "best_1_eta_WHad1",            cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_WHad1",            "best_1_eta_WHad1",            cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_1_eta_WHad2",            "best_1_eta_WHad2",            cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_WHad2",            "best_1_eta_WHad2",            cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_WHad2",            "best_1_eta_WHad2",            cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_WHad2",            "best_1_eta_WHad2",            cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_1_eta_bHad",             "best_1_eta_bHad",             cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_bHad",             "best_1_eta_bHad",             cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_bHad",             "best_1_eta_bHad",             cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_bHad",             "best_1_eta_bHad",             cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_1_eta_b1",               "best_1_eta_b1",               cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_b1",               "best_1_eta_b1",               cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_b1",               "best_1_eta_b1",               cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_b1",               "best_1_eta_b1",               cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_1_eta_b2",               "best_1_eta_b2",               cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_b2",               "best_1_eta_b2",               cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_b2",               "best_1_eta_b2",               cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_eta_b2",               "best_1_eta_b2",               cut_cat6 , "cat6" , binvec, +1)
 


binvec = cms.vdouble()
for b in range(11):
    binvec.append(0.3141592*b)
submitDataCardMakerFWlite_all( "best_0_dphi_MET_WLep1", "best_0_dphi_MET_WLep1", cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_dphi_MET_WLep1", "best_0_dphi_MET_WLep1", cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_dphi_MET_WLep1", "best_0_dphi_MET_WLep1", cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_dphi_MET_WLep1", "best_0_dphi_MET_WLep1", cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_0_dphi_MET_WLep2", "best_0_dphi_MET_WLep2", cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_0_dphi_b1_b2", "best_0_dphi_b1_b2", cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_dphi_b1_b2", "best_0_dphi_b1_b2", cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_dphi_b1_b2", "best_0_dphi_b1_b2", cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_0_dphi_b1_b2", "best_0_dphi_b1_b2", cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_1_dphi_MET_WLep1", "best_1_dphi_MET_WLep1", cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_dphi_MET_WLep1", "best_1_dphi_MET_WLep1", cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_dphi_MET_WLep1", "best_1_dphi_MET_WLep1", cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_dphi_MET_WLep1", "best_1_dphi_MET_WLep1", cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_1_dphi_MET_WLep2", "best_1_dphi_MET_WLep2", cut_cat6 , "cat6" , binvec, +1)
submitDataCardMakerFWlite_all( "best_1_dphi_b1_b2", "best_1_dphi_b1_b2", cut_cat1 , "cat1" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_dphi_b1_b2", "best_1_dphi_b1_b2", cut_cat2 , "cat2" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_dphi_b1_b2", "best_1_dphi_b1_b2", cut_cat3 , "cat3" , binvec, 0)
submitDataCardMakerFWlite_all( "best_1_dphi_b1_b2", "best_1_dphi_b1_b2", cut_cat6 , "cat6" , binvec, +1)

binvec = cms.vdouble()
for b in range(11):
    binvec.append(10. + 3*b)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b)", "logPb", cut_cat1_H , "cat1_H" , binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b)", "logPb", cut_cat1_L , "cat1_L" , binvec, 0)


binvec = cms.vdouble()
for b in range(11):
    binvec.append(20. + 3*b)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b)", "logPb", cut_cat2_H , "cat2_H" , binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b)", "logPb", cut_cat3_H , "cat3_H" , binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b)", "logPb", cut_cat6_H , "cat6_H" , binvec, 1)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b)", "logPb", cut_cat2_L , "cat2_L" , binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b)", "logPb", cut_cat3_L , "cat3_L" , binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b)", "logPb", cut_cat6_L , "cat6_L" , binvec, 1)

binvec = cms.vdouble()
for b in range(12):
    binvec.append(0. + 4*b)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_s)", "logPs", cut_cat1_H , "cat1_H" , binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_s)", "logPs", cut_cat2_H , "cat2_H" , binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_s)", "logPs", cut_cat3_H , "cat3_H" , binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_s)", "logPs", cut_cat6_H , "cat6_H" , binvec, 1)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_s)", "logPs", cut_cat1_L , "cat1_L" , binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_s)", "logPs", cut_cat2_L , "cat2_L" , binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_s)", "logPs", cut_cat3_L , "cat3_L" , binvec, 0)    
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_s)", "logPs", cut_cat6_L , "cat6_L" , binvec, 1)
    
    

binvec = cms.vdouble()
for b in range(13):
    binvec.append(-4. + 2*b)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b_ttbb/p_125_all_b)", "logPbb", cut_cat1_H , "cat1_H" , binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b_ttbb/p_125_all_b)", "logPbb", cut_cat2_H , "cat2_H" , binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b_ttbb/p_125_all_b)", "logPbb", cut_cat3_H , "cat3_H" , binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b_ttbb/p_125_all_b)", "logPbb", cut_cat6_H , "cat6_H" , binvec, 1)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b_ttbb/p_125_all_b)", "logPbb", cut_cat1_L , "cat1_L" , binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b_ttbb/p_125_all_b)", "logPbb", cut_cat2_L , "cat2_L" , binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b_ttbb/p_125_all_b)", "logPbb", cut_cat3_L , "cat3_L" , binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b_ttbb/p_125_all_b)", "logPbb", cut_cat6_L , "cat6_L" , binvec, 1)

binvec = cms.vdouble()
for b in range(11):
    binvec.append(-10. + 2*b)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b_ttjj/p_125_all_b)", "logPjj", cut_cat1_H , "cat1_H" , binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b_ttjj/p_125_all_b)", "logPjj", cut_cat2_H , "cat2_H" , binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b_ttjj/p_125_all_b)", "logPjj", cut_cat3_H , "cat3_H" , binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b_ttjj/p_125_all_b)", "logPjj", cut_cat6_H , "cat6_H" , binvec, 1)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b_ttjj/p_125_all_b)", "logPjj", cut_cat1_L , "cat1_L" , binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b_ttjj/p_125_all_b)", "logPjj", cut_cat2_L , "cat2_L" , binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b_ttjj/p_125_all_b)", "logPjj", cut_cat3_L , "cat3_L" , binvec, 0)
submitDataCardMakerFWlite_all( "TMath::Log(p_125_all_b_ttjj/p_125_all_b)", "logPjj", cut_cat6_L , "cat6_L" , binvec, 1)

    

