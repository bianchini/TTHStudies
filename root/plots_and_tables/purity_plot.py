import ROOT
import tdrstyle
tdrstyle.tdrstyle()
from collections import OrderedDict as dict
import argparse
import math

from array import array
from histlib import get_ratio
ROOT.gROOT.SetBatch(ROOT.kTRUE) #dont show graphics (messes things up)  

infilepath = "../datacards/June03_PAS/control_final_fit/"

mc_file = "mlfitMEM_COMB_New_rec_std_sb_wplots.root"
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
binning = [-6.3, -6.3+step, -6.3+2*step, -6.3+3*step, -6.3+4*step, -6.3+5*step,  -6.3+6*step, -1.3] #variable bin size

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

                purity_hist.Fill(pur_bin, bkg.GetBinContent(ibin) + signal_pf.GetBinContent(ibin))
  
                purity_hist_bkg.Fill(pur_bin, bkg.GetBinContent(ibin) )
                purity_hist_signal.Fill(pur_bin, signal_pf.GetBinContent(ibin))
                data_purity_hist.Fill(pur_bin, data.GetBinContent(ibin))


print "Purity values: "
print sorted(pur_list)


#----------- create signal + bkg stack --------------------
st = ROOT.THStack("st","")

purity_hist_bkg.SetFillColor(14)
purity_hist_bkg.SetLineColor(14)
purity_hist_signal.SetFillColor(ROOT.kRed)
purity_hist_signal.SetLineColor(ROOT.kRed)
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

data_mc_ratio = get_ratio(data_purity_hist, purity_hist_bkg, ymin=0.5, ymax = 1.5, ratio_ytitle="Data/MC")
data_mc_ratio_sb = get_ratio(purity_hist, purity_hist_bkg, ymin=0.5, ymax = 1.5, ratio_ytitle="")

error_band = get_ratio(data_purity_hist_for_ratio, purity_hist, ymin=0., ymax = 2., ratio_ytitle= "")

#--------------------- style --------------------
purity_hist.GetXaxis().SetTitle("log(S/B)")
purity_hist.GetXaxis().SetTitleSize(0.0375)
purity_hist.GetXaxis().SetNdivisions(8)
purity_hist.GetXaxis().SetLabelSize(0.0375)
purity_hist.GetYaxis().SetLabelSize(0.0375)
purity_hist.SetStats(False)
purity_hist.SetLineWidth(2)
purity_hist.SetMaximum(2*max(purity_hist.GetMaximum(), data_purity_hist.GetMaximum()) )
purity_hist.SetMinimum(20)

purity_lin.SetLineWidth(2)
purity_lin.SetLineColor(ROOT.kRed)

data_mc_ratio.SetTitle("log(S/B)")
data_mc_ratio.SetMarkerColor(1)
data_mc_ratio.SetMarkerStyle(20)
data_mc_ratio.SetMarkerSize(1)

error_band_main.SetLineColor(ROOT.kBlack)
error_band_main.SetFillColor(ROOT.kBlack)
error_band_main.SetMarkerStyle(1)
error_band_main.SetFillStyle(3004)

#------------------- plot -----------------------
c = ROOT.TCanvas("pur" , "pur", 800, 1000)

p1 = ROOT.TPad("p1", "p1", 0, 0.25, 1, 1)
p1.SetBottomMargin(0)

p1.Draw()
p1.SetTicks(1, 1);
p1.SetFillStyle(0);
p1.SetLogy()
p1.cd()

purity_hist.Draw()
st.Draw("histsame")
purity_hist.Draw("histsame")
error_band_main.Draw("e2same")
data_purity_hist.Draw("epsame")


legend1 = ROOT.TLegend(0.6, 0.75, 0.95, 0.88, "", "brNDC")
legend1.SetBorderSize(0)
legend1.SetFillColor(0)
legend1.SetTextSize(0.0375)
legend1.AddEntry(data_purity_hist, "Data", "p")
legend1.AddEntry(purity_hist, "Expectation", "l")
legend1.AddEntry(purity_hist_signal, "Signal (#mu = 0.67)", "f")
legend1.AddEntry(purity_hist_bkg, "Background", "f")
legend1.Draw()

latex = ROOT.TLatex()
latex.SetNDC()
latex.SetTextSize(0.0375)
latex.SetTextAlign(31)
latex.SetTextAlign(11)
        
cut = "CMS Preliminary"

std_txt = cut + " #sqrt{s}=8 TeV, L=19.5 fb^{-1}"                                                       
#cat_txt = regs[reg]

latex.DrawLatex(0.15, 0.96, std_txt)
#latex.DrawLatex(0.71, 0.89, cat_txt)

c.cd()

p2 = ROOT.TPad("p2","p2", 0, 0.02, 1, 0.18)
p2.SetTopMargin(0.0)
p2.SetGrid();
p2.SetFillStyle(0);
p2.Draw()
p2.cd()

data_mc_ratio.Draw("ep")
data_mc_ratio_sb.SetLineWidth(2)
data_mc_ratio_sb.SetLineColor(ROOT.kRed)
data_mc_ratio_sb.Draw("histsame")

error_band.SetLineColor(ROOT.kBlack)
error_band.SetLineWidth(2)
error_band.DrawCopy("histsame")
error_band.SetFillColor(ROOT.kBlack)
error_band.SetMarkerStyle(1)
error_band.SetFillStyle(3004)
error_band.Draw("e2same")

#purity_lin.Draw("histsame") # to add purity on the ratio band

c.SaveAs("plots/purity.png")
c.Close()
