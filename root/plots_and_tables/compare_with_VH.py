import ROOT
import tdrstyle
tdrstyle.tdrstyle()
from collections import OrderedDict as dict
import argparse
import math

cut_cat1   = "( type==0 || (type==3 && flag_type3>0))   && btag_LR>=0."
cut_cat2   = "( type==1 || (type==3 && flag_type3<=0) ) && btag_LR>=0."
cut_cat3   = "type==2 && flag_type2<=999 && btag_LR>=0."
cut_cat6   = "type==6 && btag_LR>=0."

cut_cat1_H = "( type==0 || (type==3 && flag_type3>0)) && btag_LR>=0.995"
cut_cat1_L = "( type==0 || (type==3 && flag_type3>0)) && btag_LR<0.995 && btag_LR>=0.960"

cut_cat2_H = "( type==1 || (type==3 && flag_type3<=0) ) && btag_LR>=0.9925"
cut_cat2_L = "( type==1 || (type==3 && flag_type3<=0) ) && btag_LR<0.9925 && btag_LR>=0.960"

cut_cat3_H = "type==2 && flag_type2<=999 && btag_LR>=0.995"
cut_cat3_L = "type==2 && flag_type2<=999 && btag_LR<0.995 && btag_LR>=0.970"

cut_cat6_H = "type==6 && btag_LR>=0.925"
cut_cat6_L = "type==6 && btag_LR<0.925 && btag_LR>=0.850"


indir = "../datacards/June03_PAS/data_trees/"
root_files = dict()
root_files["TTH"] = ROOT.TFile(indir + "MEAnalysisNew_all_rec_std_TTH125.root")


trees = dict()
for file in root_files:
    trees[file] = root_files[file].Get("tree")

#trees["TTH"].Draw("MET_pt")

final_selection_SL = "(" + cut_cat1_H + ") || (" + cut_cat2_H + ") || (" + cut_cat3_H + ") || (" + cut_cat1_L + ") || ( " + cut_cat2_L + ") || (" + cut_cat3_L + ")"

final_selection_DL = "( (" + cut_cat6_H + ") || (" + cut_cat6_L + ") ) && ( (Vtype<=1 && abs(Mll-91.2)>8. && Mll>15.) )"

ptV = ROOT.TH1F("ptV", "ptV", 100,0,400)
met = ROOT.TH1F("met", "met", 100,0,250)
mll = ROOT.TH1F("mll", "mll", 50,0,200)

#trees["TTH"].Draw("(sqrt((sin(MET_phi)*MET_pt + sin(lepton_phi)*lepton_pt)**2 + (cos(MET_phi)*MET_pt + cos(lepton_phi)*lepton_pt)**2))>>ptV", final_selection_SL)

trees["TTH"].Draw("(sqrt( (sin(lepton_phi[0])*lepton_pt[0] + sin(lepton_phi[1])*lepton_pt[1])**2 + (cos(lepton_phi[0])*lepton_pt[0] + cos(lepton_phi[1])*lepton_pt[1])**2 ) )>>ptV", final_selection_DL)


#trees["TTH"].Draw("MET_pt>>met", final_selection_SL)
#trees["TTH"].Draw("Mll>>mll", final_selection_DL)
