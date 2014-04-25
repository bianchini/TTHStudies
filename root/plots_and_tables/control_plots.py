import ROOT
import tdrstyle
tdrstyle.tdrstyle()
from collections import OrderedDict as dict
import argparse
from histlib import colors, stackplot
from systematics import systematics_list, find_sum_sys, get_tot_sys

ROOT.gROOT.SetBatch(ROOT.kTRUE) #dont show graphics (messes things up)

parser = argparse.ArgumentParser()
parser.add_argument('--mode', dest='mode',  choices=["DL", "SL"], required=True, help="specify DL or SL analysis")
args = parser.parse_args()

#mode = "SL"

#inpath = "../datacards/Apr15_2014_3jets40_merged/"
inpath = "../datacards/Apr24_2014_control_plots_merged/"

version = "MEM_New_ntuplizeAll_v3_rec_std_"
#var = "btag_LR"
vars = {
#    "btag_LR": "btag CSV LR",
#    "electron_eta": "electron #eta",
#    "electron_pt": "electron #pt",
#    "electron_rIso": "electron relIso",
#    "lepton_eta": "muon #eta",
#    "lepton_pt": "muon p_{T}",
#    "lepton_rIso": "muon relIso",

#    "lepton_eta": "electron #eta",
#    "lepton_pt": "electron p_{T}",
#    "lepton_rIso": "electron relIso",

#    "lepton_eta": "muon #eta",
#    "lepton_pt": "muon p_{T}",
#    "lepton_rIso": "muon relIso",  

    "numJets": "nr jets",
#    "numBTagM": "nr b-tag (CSV medium)",
    }


if args.mode == "SL":
    regs = [
#        "SL_g6j", 
#        "SL_5j",
        "SL_g4jg2t",
#        "SL_g5jg2t",
 #       "SL_g6jg3t",
        ]
if args.mode == "DL":
    regs = [
        "DL_g4j",
        "DL_g4jg2t",
        "DL_g2jg2t",
        #"DL_g3jg2t",
        ]

proc_mc = dict() #filename: [histname, pretty name]
proc_mc["TTH125"] = ["TTH125", "$t\\bar{t}H(\\rightarrow b\\bar{b})$ (125)" ]
proc_mc["EWK"] = ["EWK", "$V$ + jets"]
proc_mc["SingleT"] = ["SingleT", "single-$t$"]
proc_mc["DiBoson"] = ["DiBoson", "diboson"]
proc_mc["TTV"] = ["TTV", " $t\\bar{t} + V$ "]
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
#                histvar= "lepton_" + var.split("_")[1]
#            else:
#                histvar=var

            print "mc_hist = " + var + "/" + proc_mc[proc][0]
            mc[proc] = inputfiles[proc].Get(hist + "/" + proc_mc[proc][0])
            mc_up[proc] = find_sum_sys(proc, proc_mc[proc][0], systematics_list, inputfiles[proc], var, "Up") # sys_err_up per process
            mc_down[proc] = find_sum_sys(proc, proc_mc[proc][0], systematics_list, inputfiles[proc], var, "Down")


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

#            if var[:8] == "electron" or var[:4] == "muon":
#                histvar= "lepton_" + var.split("_")[1]
#            else:
#            histvar=var
            print "data hist = " + hist + "/" + proc_data[proc][0]

            data[proc] = inputfiles[proc].Get(hist + "/" + proc_data[proc][0])
    
        print data

        if len(data) == 0 or len(mc)==0:
            print "No input files found, skip"
            continue

        dataSum = data["Run2012_SingleMu"].Clone("dataSum")
        if args.mode == "SL":
            dataSum.Add(data["Run2012_SingleElectron"])
        elif args.mode == "DL":
            dataSum.Add(data["Run2012_DoubleElectron"])

        signal = mc["TTH125"].Clone("signal")
        signal.SetLineColor(ROOT.kRed-3)
        signal.SetLineWidth(2)
        signal.SetFillStyle(0)
        signal.Scale(50)

        stackplot(dataSum, mc, mc_up, mc_down, signal, var, vars[var], reg)

