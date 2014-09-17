import ROOT
import CMS_lumi
import tdrstyle
#tdrstyle.setTDRStyle()
tdrstyle.tdrstyle()
from numpy import arange

from collections import OrderedDict as dict
import argparse
import math

from array import array
from histlib import get_ratio, colors, style_hist, style_axes, get_poisson_err, get_poisson_ratio, add_cms_info
ROOT.gROOT.SetBatch(ROOT.kTRUE) #dont show graphics (messes things up)  

plot_style="new" # "pas" for old styling
CMS_lumi.writeExtraText = 0
CMS_lumi.lumi_8TeV = "19.5 fb^{-1}"
CMS_lumi.lumiTextSize = 0.75
CMS_lumi.cmsTextSize = 0.9

#infilepath = "../datacards/June03_PAS/control_final_fit/" # for pas
infilepath = "../datacards/Sep_PAPER/"

#mc_file = "mlfitMEM_COMB_New_rec_std_sb_wplots.root"
mc_file = "mlfitMEM_COMB_New_rec_std_sb.root"
data_file = "MEM_New_rec_std_sb.root"


infilename_mc = infilepath + mc_file
infilename_data = infilepath + data_file

print "Opening input file for MC: " + infilename_mc
try:
    f=open(infilename_mc)
    f.close()
except IOError:
    print "File " + infilename_mc + " doesn't exist -- continue"

print "Opening input file for data: " + infilename_data
try:
    f=open(infilename_data)
    f.close()
except IOError:
    print "File " + infilename_data + " doesn't exist -- continue"

infile_mc = ROOT.TFile(infilename_mc)
infile_data =ROOT.TFile(infilename_data)

prh = {} # dictionary of 8 purity histograms for each category 
step=0.6
#binning = [-6.3, -6.3+step, -6.3+2*step, -6.3+3*step, -6.3+4*step, -6.3+5*step,  -6.3+6*step, -1.3] #variable bin size

binvec_1 = arange(-6.3, -2.7, 0.45)
binning = binvec_1.tolist()
print binvec_1
binning.append(-2.7)
binning.append(-1.3)
#binning = [-6.3, -6.3+0.5*step, -6.3+step, -6.3+1.5*step, -6.3+2*step, -6.3+2.5*step, -6.3+3*step, -6.3+3.5*step, -6.3+4*step, -6.3+4.5*step, -6.3+5*step, -6.3+6*step, -1.3] #variable bin size

#purity_hist = ROOT.TH1F("purity","purity", 8, -7., -1.)
purity_hist = ROOT.TH1F("purity","purity", len(binning)-1,  array('d', binning) )
purity_hist_lin = purity_hist.Clone("purity_hist_lin")
purity_hist_bkg = purity_hist.Clone("purity_hist_bkg")
purity_hist_signal = purity_hist.Clone("purity_hist_signal")

data_purity_hist = purity_hist.Clone("data_purity_hist")

histname_map = { # different naming convention in data file
    "11": "MEM_cat1_H",
    "12": "MEM_cat1_L",
    "13": "MEM_cat2_H",
    "14": "MEM_cat2_L",
    "15": "MEM_cat3_H",
    "16": "MEM_cat3_L",
    "21": "MEM_cat6_H",
    "22": "MEM_cat6_L",
    }

pur_list = []
for i in range(1,3): # SL/DL
    for j in range(1,7): # category

        if i == 1 or (i == 2 and j <=2): # if SL or DL and cat 1,2 
            bkg_hist_name = "shapes_fit_s/Name" + str(i) + "_Name" + str(j) + "/total_background"
            signal_hist_name = "shapes_prefit/Name" + str(i) + "_Name" + str(j) + "/total_signal"
            signal_hist_name_pf = "shapes_fit_s/Name" + str(i) + "_Name" + str(j) + "/total_signal" #post fit signal

            data_hist_name = histname_map[str(i)+str(j)] + "/data_obs" 

            bkg = infile_mc.Get(bkg_hist_name)
            signal = infile_mc.Get(signal_hist_name)
            signal_pf = infile_mc.Get(signal_hist_name_pf)
            data = infile_data.Get(data_hist_name)

            for ibin in range(1,bkg.GetNbinsX()+1):
                if bkg.GetBinContent(ibin) != 0:
                    pur_bin = math.log(signal.GetBinContent(ibin)/bkg.GetBinContent(ibin)) #purity at bin in log scale
                    if pur_bin > -2.7:
                        print "SL/DL = " + str(i) + ", cat = " + str(j) + ", in bin " + str(ibin) + ", purity = " + str(pur_bin) + ", nr signal = " + str(signal_pf.GetBinContent(ibin)) + ", nr bkg = " + str(bkg.GetBinContent(ibin)) + ", nr data = " + str(data.GetBinContent(ibin)) 
                    pur_list.append(pur_bin)
                else:
                    print "Warning no bkg present in the bin"
                    pur_bin = 0 #log1

                purity_hist_bkg.Fill(pur_bin, bkg.GetBinContent(ibin) )
                purity_hist_signal.Fill(pur_bin, signal.GetBinContent(ibin)) # SM xs signal
#                purity_hist_signal.Fill(pur_bin, signal_pf.GetBinContent(ibin)) # fitted normalization signal
                data_purity_hist.Fill(pur_bin, data.GetBinContent(ibin))

data_purity_hist_poisson = get_poisson_err(data_purity_hist)
data_purity_hist_poisson = style_hist(data_purity_hist_poisson, is_data=True)

purity_hist = purity_hist_bkg.Clone("purity_hist")
purity_hist.Add(purity_hist_signal)

print "Purity values: "
print sorted(pur_list)

#----------- create signal + bkg stack --------------------
st = ROOT.THStack("st","")

purity_hist_bkg.SetFillColor(14)
purity_hist_bkg.SetLineColor(14)
purity_hist_signal.SetFillColor(colors["TTH125"])
purity_hist_signal.SetLineColor(colors["TTH125"])
st.Add(purity_hist_bkg)
st.Add(purity_hist_signal)

purity_lin = purity_hist.Clone("purity_lin")
#------------------ purity by bin -----------------------
#for ibin in range(1, purity_hist.GetNbinsX()+1):
#    purity_lin.SetBinContent(ibin,purity_hist_signal.GetBinContent(ibin)/purity_hist_bkg.GetBinContent(ibin) )


#-------------- get error band and ratio -------------------
error_band_main = purity_hist.Clone("error_band_main") # separate styling for error band

data_purity_hist_for_ratio = purity_hist.Clone("data_purity_hist_for_ratio")
for ibin in range(1, purity_hist.GetNbinsX()+1): # For data-only error-bars on ratio plot
    purity_hist_bkg.SetBinError(ibin, 0)
    data_purity_hist_for_ratio.SetBinError(ibin,0)

data_mc_ratio = get_ratio(data_purity_hist, purity_hist_bkg, ymin=0.5, ymax = 1.6, ratio_ytitle="Data/MC")
data_mc_ratio_sb = get_ratio(purity_hist, purity_hist_bkg, ymin=0.5, ymax = 1.6, ratio_ytitle="")
data_mc_ratio_poisson = get_poisson_ratio(data_purity_hist_poisson, data_purity_hist, purity_hist_bkg)

error_band = get_ratio(data_purity_hist_for_ratio, purity_hist, ymin=0., ymax = 2., ratio_ytitle= "")

#--------------------- style --------------------
data_mc_ratio_poisson = style_hist(data_mc_ratio_poisson, is_data=True)

purity_hist = style_axes(purity_hist, xTitle="log(S/B)", yTitle="Events") 
purity_hist.GetXaxis().SetNdivisions(8)

purity_hist.SetStats(False)
purity_hist.SetLineWidth(1)
purity_hist.SetMaximum(2*max(purity_hist.GetMaximum(), data_purity_hist.GetMaximum()) )
purity_hist.SetMinimum(20)

purity_lin.SetLineWidth(2)
purity_lin.SetLineColor(colors["TTH125"])

data_mc_ratio = style_hist(data_mc_ratio, is_data=True)
data_mc_ratio.SetTitle("log(S/B)")

error_band_main = style_hist(error_band_main, is_error_band=True)

#error_band_main.SetLineColor(ROOT.kBlack)
#error_band_main.SetFillColor(ROOT.kBlack)
#error_band_main.SetMarkerStyle(1)
#error_band_main.SetFillStyle(3004)

#------------------- plot -----------------------
c = ROOT.TCanvas("pur" , "pur", 800, 1000)

p1 = ROOT.TPad("p1", "p1", 0, 0.25, 1, 1)
p1.SetBottomMargin(0)

p1.Draw()
p1.SetTicks(0, 0);
p1.SetFillStyle(0);
p1.SetLogy()
p1.cd()

purity_hist.Draw()
st.Draw("histsame")
purity_hist.Draw("histsame")
error_band_main.Draw("e2same")
#data_purity_hist_poisson.SetMarkerColor(ROOT.kRed)
data_purity_hist_poisson.SetMarkerSize(0.8)
data_purity_hist_poisson.Draw("epsame")


legend1 = ROOT.TLegend(0.65, 0.75, 0.95, 0.9, "", "brNDC")
legend1.SetBorderSize(0)
legend1.SetFillColor(0)
legend1.SetTextSize(0.0375)
legend1.AddEntry(data_purity_hist, "Data", "lpe")
#legend1.AddEntry(purity_hist, "Expectation", "l")
legend1.AddEntry(purity_hist_signal, "Signal (#mu = 1)", "f")
legend1.AddEntry(purity_hist_bkg, "Background", "f")
legend1.AddEntry(error_band_main, "Bkg. Unc.", "f")
legend1.Draw()

#latex = ROOT.TLatex()
#latex.SetNDC()
#latex.SetTextSize(0.0375)
#latex.SetTextAlign(31)
#latex.SetTextAlign(11)
        
#cut = "CMS Preliminary"

#std_txt = cut + " #sqrt{s}=8 TeV, L=19.5 fb^{-1}"                                                       
#cat_txt = regs[reg]

#latex.DrawLatex(0.15, 0.96, std_txt)
#latex.DrawLatex(0.71, 0.89, cat_txt)

if plot_style == "pas":
    add_cms_info(lumi=19.5, com=8)
else: # according to new standard
    CMS_lumi.CMS_lumi(p1, 2, 11) # 2 -- 8TeV only, 3 -- 8 + 7

c.cd()

p2 = ROOT.TPad("p2","p2", 0, 0.02, 1, 0.18)
p2.SetTopMargin(0.0)
#p2.SetGrid();
p2.SetFillStyle(0);
p2.Draw()
p2.cd()



#data_mc_ratio_sb.SetLineWidth(2)
#data_mc_ratio_sb.SetLineColor(colors["TTH125"])
data_mc_ratio_sb = style_axes(data_mc_ratio_sb, yTitle = "Data/MC", is_ratio=True)
data_mc_ratio_sb = style_hist(data_mc_ratio_sb, is_signal=True)
error_band = style_hist(error_band, is_error_band=True)
one = error_band.Clone("one")
one = style_hist(one, line=True)


data_mc_ratio_sb.Draw("hist")
one.Draw("histsame")
error_band.Draw("e2same")
data_mc_ratio_poisson.SetMarkerSize(0.8)
data_mc_ratio_poisson.Draw("epsame")

#purity_lin.Draw("histsame") # to add purity on the ratio band

outfile = "plots/purity_2D.png"
print "saving output to: " + outfile
c.SaveAs(outfile)
c.Close()
