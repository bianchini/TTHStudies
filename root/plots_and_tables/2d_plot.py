import ROOT
import tdrstyle
tdrstyle.tdrstyle()
from collections import OrderedDict as dict
import argparse
from histlib import colors, stackplot, get_jet_count_hist
from systematics import systematics_list, find_sum_sys, get_tot_sys

#ROOT.gROOT.SetBatch(ROOT.kTRUE) #dont show graphics (messes things up)

parser = argparse.ArgumentParser()
parser.add_argument('--mode', dest='mode',  choices=["DL", "SL"], required=True, help="specify *DL* or *SL* analysis")
parser.add_argument('--skipSys', dest='skipSys', action="store_true")#, default=True, required=False)
args = parser.parse_args()

if args.skipSys:
    print "Omitting systematic uncertainties"
else:
    print "Running with MC systematics"

#inpath = "../datacards/May26_PAS/control_plots_merged/"
#inpath = "../datacards/June03_PAS/control_plots_merged/"
inpath = "../datacards/Sep14_cat2_studies/"

#version = "MEM_New_ntuplizeAll_v3_rec_std_"
version = "MEM_New_rec_std_"

vars = { # x-axis title, x-axis range
    "p_sb_vs_p_bj": ["P_{s/b}", [0,1,10], "P_{b/j}",[0,1,10] ],
    }

if args.mode == "SL":
    regs = {
        #-------- test cat2 -------
        "SL_cat2_HP": [
            "p_sb_vs_p_bj", 
            ],
        }
else:
    regs ={}


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

        inputfiles={} 

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

            try:
                test=mc[proc].Integral()
            except AttributeError:
                print "Failed to open histogram: " + hist + "/" + proc_mc[proc][0] + " from file " + str(inputfiles[proc])
                continue
            if proc != "QCD_BCtoE" and not args.skipSys:
                mc_up[proc] = find_sum_sys(proc, proc_mc[proc][0], systematics_list, inputfiles[proc], hist, "Up") # sys_err_up per process
                mc_down[proc] = find_sum_sys(proc, proc_mc[proc][0], systematics_list, inputfiles[proc], hist, "Down")
            else:
                mc_up[proc] = mc[proc]
                mc_down[proc] = mc[proc]

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

        signal.Draw("MEM_p_sb_vs_p_bj/ttH_hbb")

        h_sumMC = mc["TTV"].Clone("h_sumMC") #FIXME                                                                                                                         
        
        for proc in mc:
            if not proc=="TTH125" and not proc=="TTV":
                h_sumMC.Add(mc[proc])

        

        

#        stackplot(dataSum, mc, mc_up, mc_down, signal, var, vars[var][0], vars[var][1], reg, outdir="plots_cat2_tests/firstBin")
