import ROOT
import tdrstyle
tdrstyle.tdrstyle()
from collections import OrderedDict as dict
import argparse
import math
from systematics import systematics_list_profile as systematics_list
from systematics import get_profile_sys

infilepath = "../datacards/June03_PAS/2D_corr/"
infilename = "MEM_New_rec_std_cat1_H.root"

infilename_full = infilepath + infilename

proc_hist = [
    "ttH_hbb",
    "singlet",
    "ttbarV",
    "diboson",
    "ewk",
    "ttbar",
    "ttbarPlusCCbar",
    "ttbarPlusB",
    "ttbarPlusBBbar"
    ]


try:
    f=open(infilename_full)
    f.close()
except IOError:
    print "File " + infilename_full + " doesn't exist -- continue"
    # break

infile = ROOT.TFile(infilename_full)

histname = "MEM_logPbvslogPs"

data = infile.Get(histname + "/data_obs")


mc=dict()
for proc in proc_hist:
    mc[proc] = infile.Get(histname + "/" +proc)
    print proc + ": " +  str(mc[proc].Integral())

#    mc_up[proc] = find_sum_sys(proc, proc, systematics_list, infile, histname, "Up")

mc_sys = {}
for proc in proc_hist:
    mc_sys[proc] = {}
    for sys in systematics_list:
        for var in ["Up", "Down"]:
            print "getting sys hist: " + histname + "/" + proc +"_" + sys + var
            mc_sys[proc][sys + var] = infile.Get(histname + "/" + proc +"_" + sys + var)
            print proc + "_" + sys + ": " + str(mc_sys[proc][sys+var].Integral())

mc_sys_up = {}
for sys in systematics_list:
    mc_sys_up[sys] = mc_sys["ttH_hbb"][sys+"Up"].Clone("mc_sys_up" + sys)
    for proc in proc_hist:
        if proc != "ttH_hbb":
            mc_sys_up[sys].Add(mc_sys[proc][sys+"Up"])
    print "sysUp: " + sys + " -- " + str(mc_sys_up[sys].Integral())


mc_sys_down = {}
for sys in systematics_list:
    mc_sys_down[sys] = mc_sys["ttH_hbb"][sys+"Down"].Clone("mc_sys_down" + sys)
    for proc in proc_hist:
        if proc != "ttH_hbb":
            mc_sys_down[sys].Add(mc_sys[proc][sys+"Down"])
    print "sysDown: " + sys + " -- " + str(mc_sys_down[sys].Integral())


mc_nominal = mc["ttH_hbb"].Clone("mc_nominal")
for proc in proc_hist:
    if proc != "ttH_hbb":
        mc_nominal.Add(mc[proc])


#mc_nominal.Draw()
mc_up = get_profile_sys(mc_nominal, mc_sys_up)
mc_down = get_profile_sys(mc_nominal, mc_sys_down)

for ibinx in range(mc_nominal.GetNbinsX() + 1):
    for ibiny in range(mc_nominal.GetNbinsY() + 1):
        err_sys = max(mc_up.GetBinContent(ibinx+1, ibiny+1), mc_down.GetBinContent(ibinx+1, ibiny+1))
        err_stat = mc_nominal.GetBinError(ibinx+1, ibiny+1)

        mc_nominal.SetBinError(ibinx+1, ibiny+1, math.sqrt( err_stat**2 + err_sys**2) )

prof_data = data.ProfileY()
prof_mc = mc_nominal.ProfileY()

prof_data.GetXaxis().SetRangeUser(0,40)
prof_data.GetYaxis().SetRangeUser(0,40)
prof_data.GetXaxis().SetTitle("log(Ps)")
prof_data.GetYaxis().SetTitle("log(Pb)")

prof_data.SetMarkerColor(1)
prof_data.SetMarkerStyle(1)
prof_data.SetMarkerSize(1)
prof_data.Draw()

prof_mc.SetMarkerColor(ROOT.kRed)
prof_mc.SetLineColor(ROOT.kRed)
prof_mc.SetMarkerStyle(1)
prof_mc.Draw("epsame")

mc_nominal.SetLineColor(ROOT.kRed)

#data.Draw("lego")
#mc_nominal.Draw("legosame")
