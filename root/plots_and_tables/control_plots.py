import ROOT
import tdrstyle
tdrstyle.tdrstyle()
from collections import OrderedDict as dict
import argparse
from histlib import colors, stackplot, get_jet_count_hist
from systematics import systematics_list, find_sum_sys, get_tot_sys

ROOT.gROOT.SetBatch(ROOT.kTRUE) #dont show graphics (messes things up)

parser = argparse.ArgumentParser()
parser.add_argument('--mode', dest='mode',  choices=["DL", "SL"], required=True, help="specify DL or SL analysis")
args = parser.parse_args()

#mode = "SL"

#inpath = "../datacards/Apr15_2014_3jets40_merged/"
#inpath = "../datacards/Apr24_2014_control_plots_merged/"
inpath = "../datacards/Apr28_checks/controlPlots_merged/"

version = "MEM_New_ntuplizeAll_v3_rec_std_"
#version = "MEM_New_rec_std_"

vars = {
    "btag_LR": "b_{LR}",
#    "MET_pt": "MET",
#    "electron_eta": "electron #eta",
#    "electron_pt": "electron #pt",
#    "electron_rIso": "electron relIso",

#    "muon_eta": "muon #eta",
#    "muon_pt": "muon #pt",
#    "muon_rIso": "muon relIso",

#    "Mll": "m(l^{+}l^{-})",
#    "MTln": "m_{T}(l #nu)",

#    "numJets": "nr jets",
#    "numBTagM": "nr b-tag (CSV medium)",
    }


if args.mode == "SL":
    regs = [
        "SL_5j_ewk",
        "SL_6j_ewk",
        "SL_g6jg2t_ewk", 
        "SL_5jg2t_ewk",
#        "SL_g5jg2tmt60",
#        "SL_g6j",
#        "SL_5j",
#        "SL_g4jg2t",
#        "SL_g5jg2t",
#        "SL_g5jg3t",
#        "SL_tot_3j40",
#        "SL_tot_mt60_3j40"
#        "SL_tot_3t_3tag"
        ]
if args.mode == "DL":
    regs = [
        "DL_4j_z",
#        "DL_4j",
#        "DL_g4j",
#        "DL_g2jg2t",
#        "DL_g2jg2t",
        #"DL_g3jg2t",
        ]

proc_mc = dict() #filename: [histname, pretty name]
proc_mc["TTH125"] = ["TTH125", "$t\\bar{t}H(\\rightarrow b\\bar{b})$ (125)" ]
proc_mc["SingleT"] = ["SingleT", "single-$t$"]
proc_mc["TTV"] = ["TTV", " $t\\bar{t} + V$ "]
proc_mc["DiBoson"] = ["DiBoson", "diboson"]
proc_mc["EWK"] = ["EWK", "$V$ + jets"]
proc_mc["TTJetsJJ"] = ["TTJetsLF", " $t\\bar{t} + jj$"]
proc_mc["TTJetsBJ"] = ["TTJetsHFb", " $t\\bar{t} + bj$"]
proc_mc["TTJetsBB"] = ["TTJetsHFbb", " $t\\bar{t} + b\\bar{b}$ "]

proc_data = dict()

if args.mode == "SL":
    proc_data["Run2012_SingleElectron"] = ["data_obs", ""]
    proc_data["Run2012_SingleMu"] = ["data_obs", ""]
if args.mode == "DL":
    proc_data["Run2012_SingleMu"] = ["data_obs", ""]
    proc_data["Run2012_DoubleElectron"] = ["data_obs", ""]

for reg in regs:
    for var in vars:
        hist = "MEM_" + var

        mc=dict()
        mc_up=dict()
        mc_down=dict()
        jet_count = dict()
        jet_count_up = dict()
        jet_count_down = dict()

        inputfiles={} # needed ?

        for proc in proc_mc:
            infile = inpath + version +var + "_" + reg+ "_" +proc + ".root"
    
            try:
                f=open(infile)
                f.close()
            except IOError:
                print "File " + infile + " doesn't exist -- continue"
                break

            inputfiles[proc] =ROOT.TFile(infile)

            
#            if var[:8] == "electron" or var[:4] == "muon":
#                hist= "MEM_lepton_" + var.split("_")[1]
#                print "opening histogram: " + hist
#            else:
#                histvar=var

            mc[proc] = inputfiles[proc].Get(hist + "/" + proc_mc[proc][0])
            print mc[proc].Integral()
            mc_up[proc] = find_sum_sys(proc, proc_mc[proc][0], systematics_list, inputfiles[proc], hist, "Up") # sys_err_up per process
            mc_down[proc] = find_sum_sys(proc, proc_mc[proc][0], systematics_list, inputfiles[proc], hist, "Down")

            jet_labels = {} # xaxis-range, label name (for jet count hists)
            if args.mode == "SL":
                jet_labels["numJets"] = ( [4,10], "jets")
                jet_labels["numBTagM"] =( [2,6], "b-tags")
            if args.mode == "DL":
                jet_labels["numJets"] = ( [2,6], "jets")
                jet_labels["numBTagM"] =( [0,6], "b-tags")

            if var == "numJets" or var == "numBTagM":
                jet_count[proc] = get_jet_count_hist(mc[proc], jet_labels[var][0], jet_labels[var][1])
                jet_count_up[proc] = get_jet_count_hist(mc_up[proc], jet_labels[var][0], jet_labels[var][1])
                jet_count_down[proc] = get_jet_count_hist(mc_down[proc], jet_labels[var][0], jet_labels[var][1])

                print "jet_count = " + str(jet_count[proc].GetBinContent(1))
                mc[proc].Delete()
                mc_up[proc].Delete()
                mc_down[proc].Delete()
                mc[proc] = jet_count[proc]
                mc_up[proc] = jet_count_up[proc]
                mc_down[proc] = jet_count_down[proc]

        data=dict()
        for proc in proc_data:
            infile = inpath + version +var + "_" + reg + "_" +proc + ".root"
            print "opening file: " + infile

            try:
                f=open(infile)
                f.close()
            except IOError:
                print "File " + infile + "doesn't exist -- continue"
                break

            inputfiles[proc] =ROOT.TFile(infile)

 #           if var[:8] == "electron" or var[:4] == "muon":
 #               hist= "MEM_lepton_" + var.split("_")[1]

            print "data hist = " + hist + "/" + proc_data[proc][0]
            data[proc] = inputfiles[proc].Get(hist + "/" + proc_data[proc][0])
    
        if len(data) == 0 or len(mc)==0:
            print "No input files found, skip"
            continue

        dataSum = data["Run2012_SingleMu"].Clone("dataSum")
        if args.mode == "SL":
            dataSum.Add(data["Run2012_SingleElectron"])
        elif args.mode == "DL":
            dataSum.Add(data["Run2012_DoubleElectron"])

        if var == "numJets" or var == "numBTagM":
            jet_count["dataSum"] = get_jet_count_hist(dataSum, jet_labels[var][0], jet_labels[var][1])
            dataSum.Delete()
            dataSum = jet_count["dataSum"]
             

#        print mc

        signal = mc["TTH125"].Clone("signal")
        signal.SetLineColor(ROOT.kRed-3)
        signal.SetLineStyle(ROOT.kDashed)
        signal.SetLineWidth(3)
        signal.SetFillStyle(0)

        stackplot(dataSum, mc, mc_up, mc_down, signal, var, vars[var], reg)

