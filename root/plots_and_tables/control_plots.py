import ROOT
import tdrstyle
tdrstyle.tdrstyle()
from collections import OrderedDict as dict
import argparse
from histlib import colors, stackplot, get_jet_count_hist
from systematics import systematics_list, find_sum_sys, get_tot_sys

ROOT.gROOT.SetBatch(ROOT.kTRUE) #dont show graphics (messes things up)

parser = argparse.ArgumentParser()
parser.add_argument('--mode', dest='mode',  choices=["DL", "SL"], required=True, help="specify *DL* or *SL* analysis")
parser.add_argument('--skipSys', dest='skipSys', action="store_true")#, default=True, required=False)
args = parser.parse_args()
plot_style="paper"

rebin = 1 # FIXME, doesnt work atm

if args.skipSys:
    print "Omitting systematic uncertainties"
else:
    print "Running with MC systematics"

#inpath = "../datacards/May26_PAS/control_plots_merged/"
#inpath = "../datacards/June03_PAS/control_plots_merged/"
#inpath = "../datacards/Sep14_cat2_studies/"
inpath = "../datacards/Sep_PAPER/btag_LR_5j1t_10bin/"

version = "MEM_New_ntuplizeAll_v3_rec_std_"
#version = "MEM_New_rec_std_"

vars = { # x-axis title, x-axis range
    "btag_LR": [" \\mathscr{F}", [0,1] ],

    "electron_eta": ["electron #eta", [-2.5, 2.5] ],
    "electron_pt": ["electron p_{T}", [30, 250] ],
    "electron_rIso": ["electron r_{iso}", [0, 0.115] ], #needed for showing 0.12 upper bound
    "electron_dxy": ["electrond dxy", [0, 0.025] ],

    "muon_eta": ["muon #eta", [-2.5, 2.5] ],
    "muon_pt": ["muon p_{T}", [30, 250] ],
    "muon_rIso": ["muon r_{iso}", [0, 0.115] ],
    "muon_dxy": ["muon dxy", [0, 0.025] ],

    "Mll": ["m(l^{+}l^{-})", [30,350] ],
    "MTln": ["m_{T}(l #nu)", [30,450] ],
    "MET_pt": ["MET", [0,250] ],
    "MET_sumEt": ["MET_sumEt", [350,3500]],

    "numJets": ["jet multiplicity", [4, 10] ],
    "numBTagM": ["multiplicity of b-tagged jets (CSVM)", [0, 5] ],
    "jetsAboveCut": ["nr jets with p_{T} > 40", [0, 10] ],

    "nPVs": ["# primary vertices", [0,40]],

    "bjet_pt": ["b-jet p_{T}", [30, 250]],
    "bjet_eta": ["b-jet #eta", [-2.5, 2.5]],

    "leadjet_pt": ["leading jet p_{T}", [30, 340]],
    "leadjet_eta": ["leading jet #eta", [-2.5, 2.5]],

    "p_bj": ["P_{b/j}", [0,1,10] ],
    "p_sb": ["P_{s/b}", [0,1,10] ],
    }

if args.mode == "DL":
    vars["electron_pt"][1] = [20, 250]
    vars["muon_pt"][1] = [20, 250]
    vars["numJets"][1] = [0, 10]
    vars["numBTagM"][1] = [1, 4]


if args.mode == "SL":
    regs = {
#        "SL_5j": ["btag_LR"],
#        "SL_6j": ["btag_LR"],
#        "SL_g6jg2t": ["btag_LR"],
#        "SL_5jg2t": ["btag_LR"],
        "SL_5jg1t": ["btag_LR"],
        "SL_g6jg1t": ["btag_LR"],

#        "SL_g4jg2t": [ "numBTagM", "numJets"],

#        "SL_g5jg3t": [
#            "electron_pt",
#            "electron_eta",
#            "electron_rIso",

#            "muon_eta",
#            "muon_pt",
#            "muon_rIso",

#            "MET_pt",
            #"MET_sumEt",
 #           "MTln",
  #          "nPVs",

#            "leadjet_pt",
#            "leadjet_eta",
#            "bjet_pt",
#            "bjet_eta",
#]
        #-------- test cat2 -------
 #       "SL_cat1_HP": [
 #           "btag_LR",
 #           "p_sb",
 #           "p_bj"
 #           ],
 #       "SL_cat1_HP_muon": [
 #           "btag_LR"
 #           ],
 #       "SL_cat1_HP_electron": [
 #           "btag_LR"
 #           ],

#        "SL_cat2_HP": [
#            "numBTagM",
#            "btag_LR", 
#            "p_sb", 
#            "p_bj"
#            ],
        }

"""
        "SL_cat2_HP_firstBin": [
            "numBTagM",
            "numJets",
            "btag_LR",
            "p_sb",
            "p_bj",
            "electron_pt",                                                                                                                                   
            "electron_eta",                                                                                                                                  
            "electron_rIso",
            "electron_dxy",
            "muon_eta", 
            "muon_pt",                                                                                                                                       
            "muon_rIso",
            "muon_dxy",
            "MET_pt",
            "MTln",
            "leadjet_pt",          
            "leadjet_eta",                                                                                                                                   
            "bjet_pt",                                                                                                                                       
            "bjet_eta"
            ],
"""

#        "SL_cat2_HP_muon": [
#            "btag_LR",
#            "p_sb",
#            "p_bj"
#            ],
#        "SL_cat2_HP_electron": [
#            "btag_LR",
#            "p_sb",
#            "p_bj" 
#            ]

#        "SL_g5jg2t_eta15": ["electron_dxy", "electron_eta", "electron_pt", "electron_rIso", "MET_pt", "MTln", "btag_LR", "jetsAboveCut", "numBTagM"], #investigating QCD in high eta electron events
#        "SL_g5jg2t": ["MET_pt"],
#        }


if args.mode == "DL":
    regs = {
        "DL_g2jg2t": [
#            "electron_pt",
#            "electron_eta",
#            "electron_rIso",

#            "muon_eta",
#            "muon_pt",
#            "muon_rIso",

#            "bjet_pt",
#            "bjet_eta",
#            "leadjet_pt",
#            "leadjet_eta",

#            "MET_pt",
#            "MET_sumEt",
#            "Mll",

            "numJets",
            ],
#        "DL_g2jg2t": ["Mll_z"],

#        "DL_g4j_z": ["btag_LR"],
#        "DL_g4j": ["btag_LR"],
        "DL_g4jg1t": ["btag_LR"],
#        "DL_g4j": ["btag_LR", "numBTagM"],
#        "DL_g4j": ["numBTagM"],
        }

do_QCD=False
proc_mc = dict() #filename: [histname, pretty name]
proc_mc["TTH125"] = ["ttH_hbb", "$t\\bar{t}H(\\rightarrow b\\bar{b})$ (125)" ]
if do_QCD:
    proc_mc["QCD_BCtoE"] = ["_TopPtUp", "QCD"]

proc_mc["DiBoson"] = ["diboson", "diboson"]
proc_mc["TTV"] = ["ttbarV", " $t\\bar{t} + V$ "]
proc_mc["SingleT"] = ["singlet", "single-$t$"]
proc_mc["EWK"] = ["ewk", "$V$ + jets"]

proc_mc["TTJetsBB"] = ["ttbarPlusBBbar", " $t\\bar{t} + b\\bar{b}$ "]
proc_mc["TTJetsBJ"] = ["ttbarPlusB", " $t\\bar{t} + bj$"]
proc_mc["TTJetsCC"] = ["ttbarPlusCCbar", " $t\\bar{t} + cc$"]
proc_mc["TTJetsJJ"] = ["ttbar", " $t\\bar{t} + jj$"]




proc_data = dict()

if args.mode == "SL":
    proc_data["Run2012_SingleElectron"] = ["data_obs", ""]
    proc_data["Run2012_SingleMu"] = ["data_obs", ""]
if args.mode == "DL":
    proc_data["Run2012_SingleMu"] = ["data_obs", ""]
    proc_data["Run2012_DoubleElectron"] = ["data_obs", ""]

for reg in regs:
    for var in regs[reg]:
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
            mc[proc].Rebin(rebin)

            try:
                test=mc[proc].Integral()
            except AttributeError:
                print "Failed to open histogram: " + hist + "/" + proc_mc[proc][0] + " from file " + str(inputfiles[proc])
                continue
            if proc != "QCD_BCtoE" and not args.skipSys:
                mc_up[proc] = find_sum_sys(proc, proc_mc[proc][0], systematics_list, inputfiles[proc], hist, "Up") # sys_err_up per process
                mc_up[proc].Rebin(rebin)
                mc_down[proc] = find_sum_sys(proc, proc_mc[proc][0], systematics_list, inputfiles[proc], hist, "Down")
                mc_down[proc].Rebin(rebin)
            else:
                mc_up[proc] = mc[proc]
                mc_up[proc].Rebin(rebin)
                mc_down[proc] = mc[proc]
                mc_down[proc].Rebin(rebin)

            jet_labels = {} # xaxis-range, label name (for jet count hists)
            if args.mode == "SL":
                jet_labels["numJets"] = ( [4,10], "")# "jets")
                jet_labels["numBTagM"] =( [2,5], "")#"b-tags")
            if args.mode == "DL":
                jet_labels["numJets"] = ( [2,8], "")#"jets")
                jet_labels["numBTagM"] =( [1,4], "")#"b-tags")

            if var == "numJets" or var == "numBTagM":
                jet_count[proc] = get_jet_count_hist(mc[proc], jet_labels[var][0], jet_labels[var][1])
                
                if not args.skipSys:
                    jet_count_up[proc] = get_jet_count_hist(mc_up[proc], jet_labels[var][0], jet_labels[var][1])
                    jet_count_down[proc] = get_jet_count_hist(mc_down[proc], jet_labels[var][0], jet_labels[var][1])

                print "jet_count = " + str(jet_count[proc].GetBinContent(1))
                mc[proc].Delete()                
                mc[proc] = jet_count[proc]
                if not args.skipSys:
                    mc_up[proc].Delete()
                    mc_down[proc].Delete()
                    
                    mc_up[proc] = jet_count_up[proc]
                    mc_down[proc] = jet_count_down[proc]
                else:
                    mc_up[proc] = mc[proc]
                    mc_down[proc] = mc[proc]
                    

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
            data[proc].Rebin(rebin)

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


        signal = mc["TTH125"].Clone("signal")

        stackplot(dataSum, mc, mc_up, mc_down, signal, var, vars[var][0], vars[var][1], reg, outdir="plots_paper", plot_style=plot_style)
