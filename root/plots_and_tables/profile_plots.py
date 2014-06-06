import ROOT
import tdrstyle
tdrstyle.tdrstyle()
from collections import OrderedDict as dict
import argparse
import math
from systematics import systematics_list_profile as systematics_list
from systematics import get_profile_sys, get_profile_sys2
ROOT.gROOT.SetBatch(ROOT.kTRUE) #dont show graphics (messes things up)                                                       


infilepath = "../datacards/June03_PAS/2D_corr/controlPlots/"

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


regs = { 
    "cat1_H": "SL cat. 1 (H)", 
    "cat1_L": "SL cat. 1 (L)", 
    "cat2_H": "SL cat. 2 (H)",
    "cat2_L": "SL cat. 2 (L)",
    "cat3_H": "SL cat. 3 (H)",
    "cat3_L": "SL cat. 3 (L)",
    "cat6_H": "DL (H)",
    "cat6_L": "DL (L)", 
         }

do_sys_v2 = False # sum systematics of profiles, false for adding systematics to 2D histogram bin-by-bin

for reg in regs:
    infilename = "MEM_New_rec_std_" + reg + ".root"
    infilename_full = infilepath + infilename

    try:
        f=open(infilename_full)
        f.close()
    except IOError:
        print "File " + infilename_full + " doesn't exist -- continue"
        break;

    infile = ROOT.TFile(infilename_full)

    histnames = { #histname: [xaxis, yaxis, [xmin, xmax], [ymin,ymax], doProfileY]
        "MEM_logPbvslogPs": ["log(#omega_{1})", "log(#omega_{0})", [0, 40], [0,45], False],
        "MEM_logPbbvslogPjj": ["log(L_{bb}^{b-tag})", "log(L_{jj}^{b-tag})", [-5, 20], [-10, 10], False],

        "MEM_logPbvslogPbb": ["log(#omega_{1})", "log(L_{bb}^{b-tag})", [0,40], [-5,15], False],
        "MEM_logPbvslogPjj": ["log(#omega_{1})","log(L_{jj}^{b-tag})", [0,40], [-10, 5], False],

        "MEM_logPsvslogPbb": ["log(#omega_{0})","log(L_{bb}^{b-tag})", [0,40], [-5, 15], False],
        "MEM_logPsvslogPjj": ["log(#omega_{0})","log(L_{jj}^{b-tag})", [0,40], [-10, 5], False],
        }

    
    for histname in histnames:
        doProfileY = histnames[histname][4]

        data = infile.Get(histname + "/data_obs")


        mc=dict()
        print histname 
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


        mc_up = get_profile_sys(mc_nominal, mc_sys_up)
        mc_down = get_profile_sys(mc_nominal, mc_sys_down)

        for ibinx in range(mc_nominal.GetNbinsX() + 1):
            for ibiny in range(mc_nominal.GetNbinsY() + 1):
                err_sys = max(mc_up.GetBinContent(ibinx+1, ibiny+1), mc_down.GetBinContent(ibinx+1, ibiny+1))
                err_stat = mc_nominal.GetBinError(ibinx+1, ibiny+1)

                mc_nominal.SetBinError(ibinx+1, ibiny+1, math.sqrt( err_stat**2 + err_sys**2) )

        if doProfileY:
            prof_data = data.ProfileY()
            prof_mc = mc_nominal.ProfileY()
            prof_data.GetXaxis().SetTitle(histnames[histname][1])
            prof_data.GetYaxis().SetTitle(histnames[histname][0])

        else:
            prof_data = data.ProfileX()
            prof_mc = mc_nominal.ProfileX()
            prof_data.GetXaxis().SetTitle(histnames[histname][0])
            prof_data.GetYaxis().SetTitle(histnames[histname][1])

        prof_data.GetXaxis().SetRangeUser(histnames[histname][2][0], histnames[histname][2][1])
        prof_data.GetYaxis().SetRangeUser(histnames[histname][3][0], histnames[histname][3][1])
        prof_data.GetYaxis().SetTitleSize(0.0375)
        prof_data.GetXaxis().SetTitleOffset(1.2)
        prof_data.GetYaxis().SetLabelSize(0.0375)
        
        prof_data.SetMarkerColor(1)
        prof_data.SetLineColor(1)
        prof_data.SetMarkerStyle(20)
        prof_data.SetMarkerSize(1)

        prof_mc.SetMarkerColor(ROOT.kRed)
        prof_mc.SetLineColor(ROOT.kRed)
        prof_mc.SetFillColor(14)
        prof_mc.SetFillStyle(3004)
        prof_mc.SetLineColor(ROOT.kWhite)
        prof_mc.SetMarkerStyle(20)



        if do_sys_v2:
            xrange =  [mc_nominal.GetXaxis().GetXmin(), mc_nominal.GetXaxis().GetXmax()]
            error_band_mc = get_profile_sys2(mc_nominal, mc_sys_up, xrange, doProfileY)
        else:
            error_band_mc = prof_mc.Clone()
    
        error_band_mc.SetFillColor(14)
        error_band_mc.SetFillStyle(3004)
        error_band_mc.SetLineColor(ROOT.kWhite)
        error_band_mc.SetMarkerColor(ROOT.kRed)
        error_band_mc.SetMarkerStyle(20)

        for ibinx in range(prof_data.GetNbinsX()+1):
            if prof_data.GetBinError(ibinx) == 0:
                prof_data.SetBinContent(ibinx, -999)
            print prof_mc.GetBinError(ibinx)
            if prof_mc.GetBinError(ibinx) == 0:
                prof_mc.SetBinContent(ibinx, -999)
                error_band_mc.SetBinContent(ibinx, -999)

        c = ROOT.TCanvas("c" + reg ,"c" + reg, 800, 800)


        prof_data.Draw()
    #    prof_mc.Draw("e2same")
        error_band_mc.Draw("e2same")

        legend1 = ROOT.TLegend(0.7, 0.75, 0.95, 0.88, "", "brNDC")
        legend1.SetBorderSize(0)
        legend1.SetFillColor(0)
        legend1.SetTextSize(0.03)
        legend1.AddEntry(prof_data, "Data", "p")
        legend1.AddEntry(error_band_mc, "Expectation", "p")
        legend1.AddEntry(error_band_mc, "MC stat+sys", "f")
        legend1.Draw()
        
        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.03)
        latex.SetTextAlign(31)
        latex.SetTextAlign(11)
        
        cut = "CMS Preliminary"
        
        std_txt = cut + " #sqrt{s}=8 TeV, L=19.5 fb^{-1}"                                                       
        cat_txt = regs[reg]
        
        latex.DrawLatex(0.15, 0.97, std_txt)
        latex.DrawLatex(0.71, 0.89, cat_txt)
        

        outfilename = "profile_plots_v2/profile_" + histname + "_" + reg 
        if do_sys_v2:
            outfilename = outfilename + "_sysV2"

        c.SaveAs( outfilename + ".png")
        c.Close()

#    c[var].SaveAs("plots/" + var + "_" + reg + ".pdf")                                                                      
#c[var].Close()

#mc_nominal.SetLineColor(ROOT.kRed)
#data.Draw("lego")
#mc_nominal.Draw("legosame")
