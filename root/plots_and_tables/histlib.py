
import ROOT, sys, re
from systematics import get_tot_sys
from collections import OrderedDict as dict

proc_names = {
    "TTJetsJJ": "t#bar{t} + jj",
    "TTJetsCC": "t#bar{t} + cc",
    "TTJetsBJ": "t#bar{t} + b",
    "TTJetsBB": "t#bar{t} + bb",
    "TTV": "t#bar{t} + W,Z",
    "DiBoson": "VV",
    "TTH125": "t#bar{t}H (125)",
    "SingleT": "Single t",
    "EWK": "EWK",
    "QCD_BCtoE": "QCD",
}

colors = {
    "TTJetsJJ": ROOT.kRed-7,
    "TTJetsBJ": ROOT.kRed-2,
    "TTJetsBB": ROOT.kRed+3,
    "TTJetsCC": ROOT.kRed+1,
    "EWK": ROOT.kAzure+2,
    "DiBoson":ROOT.kBlue-9,
    "TTH125": ROOT.kBlue+2,
    "TTV": ROOT.kBlue-10,
    "SingleT": ROOT.kMagenta-2,
    "QCD_BCtoE": 15
    }

def style_axes(hist, xTitle="", yTitle="", is_ratio=False, is_jet_count=False):
  xAxis = hist.GetXaxis()
  yAxis = hist.GetYaxis()

  labelFont = 62

  if is_ratio:
    yLabelSize = 0.17

    xTitleOffset = 0.2
    yTitleOffset = 0.27

    xTitleSize = 0.2
    yTitleSize = 0.17

    xLabelSize = 0.0 #bdt 0.12

  else:
    yLabelSize = 0.038

    xTitleOffset = 0.97
    yTitleOffset = 1.4
    xTitleSize = 0.04
    yTitleSize = 0.04

    if is_jet_count:
      xLabelSize = 0.055
    else:
      xLabelSize = 0.038

  xAxis.SetTitle(xTitle)
  yAxis.SetTitle(yTitle)

  xAxis.SetTitleFont(labelFont)
  yAxis.SetTitleFont(labelFont)

  xAxis.SetLabelFont(labelFont)
  yAxis.SetLabelFont(labelFont)

  xAxis.SetLabelSize(xLabelSize)
  xAxis.SetTitleOffset(xTitleOffset)
  xAxis.SetTitleSize(xTitleSize)

  yAxis.SetLabelSize(yLabelSize)
  yAxis.SetTitleOffset(yTitleOffset)
  yAxis.SetTitleSize(yTitleSize)

  if is_ratio:
    xAxis.SetNdivisions(0)
    yAxis.SetNdivisions(3)

    yAxis.CenterTitle()

  return hist

def style_hist(hist, color=0, is_data = False, is_signal=False, is_error_band=False, line=False):
  if is_data:
    hist.SetMarkerStyle(20);
    hist.SetMarkerSize(1.5);
    hist.SetMarkerColor(ROOT.kBlack);
    hist.SetLineColor(ROOT.kBlack);
    hist.SetLineWidth(2);
#    hist.SetBinErrorOption(ROOT.TH1F.kPoisson)
    print hist.GetBinErrorOption()
    for ibin in range(hist.GetNbinsX()):
      print "err up =" + str(hist.GetBinErrorLow(ibin))
      print "err down = " + str(hist.GetBinErrorUp(ibin))

  elif is_signal:
    hist.SetLineColor(ROOT.kBlue+2)
#    signal.SetLineStyle(ROOT.kDashed)
    hist.SetLineWidth(4)
    hist.SetFillStyle(0)

  elif is_error_band:
    hist.SetMarkerSize(0)
    if color != 0:
      hist.SetFillColor(color)
      hist.SetFillStyle(1001)
    else:
      hist.SetFillColor(ROOT.kBlack)
      hist.SetFillStyle(3654) #3004

  elif line:
    hist.SetLineColor(ROOT.kBlack)
    hist.SetLineStyle(ROOT.kDashed)
    hist.SetFillStyle(0)
    hist.SetLineWidth(1)

  else: #regular MC stack
    hist.SetLineColor(color)
    hist.SetFillColor(color)
    hist.SetFillStyle(1001)


  return hist

def style_legend(legend):
    legend.SetBorderSize(0)
    legend.SetFillColor(0)
    legend.SetTextSize(0.04)

    return legend

def get_ratio(hist1, hist2, ymin=0., ymax=2, is_band = False, ratio_ytitle = ""):
    """
    hist1 -- numerator
    hist2 -- denominator
    """
    hist_ratio = hist1.Clone()
    hist_ratio.Divide(hist2)

    hist_ratio.SetStats(False)
    hist_ratio.SetMarkerStyle(20)
    hist_ratio.SetMarkerSize(0.35)

    if not is_band:
        hist_ratio.SetMarkerColor(ROOT.kBlack)
        hist_ratio.SetLineColor(ROOT.kBlack)
    hist_ratio.SetMaximum(ymax)
    hist_ratio.SetMinimum(ymin)

    return hist_ratio

def get_error_band(err_up, err_down, nominal_in, band_only=True):
    """
    compose histogram with errors assigned to each bin. Use max(err_up, err_down)
    band_only -- draw at errorband, nominal is a line at 1
    """
    if band_only:
        nominal = err_up.Clone("nominal")

        for ibin in range(nominal.GetNbinsX() + 1):
            nominal.SetBinContent(ibin+1, 1)
            nominal.SetBinError(ibin+1, max(err_up.GetBinContent(ibin+1)-1, err_down.GetBinContent(ibin+1)-1, 0) )

    else:
        nominal = nominal_in.Clone("nominal")

        for ibin in range(nominal.GetNbinsX() + 1):
            nominal.SetBinContent(ibin+1, nominal_in.GetBinContent(ibin+1))
            nominal.SetBinError(ibin+1, max(err_up.GetBinContent(ibin+1), err_down.GetBinContent(ibin+1), 0) )
            print "err = " + str(nominal.GetBinError(ibin+1))

    return nominal

def stackplot(dataSum, mc, mc_up, mc_down, signal, var, varname="", var_range=[-500,500], reg=""):
    """
    data -- sum of data histograms
    mc -- dictionary of mc histograms
    signal -- signal histogram
    var -- variable name
    varname -- for axis label
    """

    signal_scale = 50

    signal = style_hist(signal,is_signal=True)
    signal.Scale(signal_scale)

    for proc in mc:
        mc[proc] = style_hist(mc[proc], color = colors[proc])
    dataSum = style_hist(dataSum, is_data=True)

    sum = ROOT.THStack("sum","")

    for proc, mc_hist in mc.iteritems():
        if not proc == "TTH125": #FIXME
            sum.Add(mc_hist)

    h_sumMC = mc["TTV"].Clone("h_sumMC") #FIXME

    for proc in mc:
        if not proc=="TTH125" and not proc=="TTV":
            h_sumMC.Add(mc[proc])

    c={}
    c[var] = ROOT.TCanvas("c" + var ,"c" + var, 800, 1000)

    p1 = {}
    p1[var] = ROOT.TPad("p1", "p1", 0, 0.25, 1, 1)
    p1[var].SetBottomMargin(0)

    p1[var].Draw()
    p1[var].SetTicks(0, 0);
    p1[var].SetFillStyle(0);


    h_sumMC.SetTitle("")
    h_sumMC.SetStats(False)
    h_sumMC.SetLineWidth(1)
    h_sumMC.SetMaximum(1.3*max(h_sumMC.GetMaximum(), dataSum.GetMaximum()) )
    h_sumMC.SetMinimum(0.)
    h_sumMC.SetLineColor(ROOT.kBlack)
    h_sumMC.SetFillStyle(0)



    if var != "numJets" and var != "numBTagM":
        h_sumMC.GetXaxis().SetRangeUser( var_range[0], var_range[1])
        h_sumMC = style_axes(h_sumMC, xTitle=varname, yTitle="Events")
    else:
        h_sumMC = style_axes(h_sumMC, xTitle=varname, yTitle="Events", is_jet_count=True)



    if var == "numJets" or var == "numBTagM" or var == "cat_count" or var == "btag_LR" or var == "MTln" or var == "Mll" or var == "MET_pt" or var == "jetsAboveCut": #for logscale
        if var == "numJets" or var == "numBTagM":
            ymin_log = 1
        else:
            ymin_log = 1

        h_sumMC.SetMinimum(ymin_log)
        sum.SetMinimum(ymin_log)
        mc["TTH125"].SetMinimum(ymin_log)
        signal.SetMinimum(ymin_log)
        dataSum.SetMinimum(ymin_log)

        p1[var].SetLogy()

        if var == "numBTagM":
            h_sumMC.SetMaximum(10*ROOT.TMath.Max(h_sumMC.GetMaximum(), dataSum.GetMaximum()) )
        else:
            h_sumMC.SetMaximum(10*ROOT.TMath.Max(h_sumMC.GetMaximum(), dataSum.GetMaximum()) )


    h_sumMCup = get_tot_sys(mc_up)
    h_sumMCdown = get_tot_sys(mc_down)

    error_band_mc = get_error_band(h_sumMCup, h_sumMCdown, h_sumMC, 0)
    error_band_mc = style_hist(error_band_mc, is_error_band=True)


    p1[var].cd()

    h_sumMC.Draw("hist")
    sum.Draw("histsame")
    h_sumMC.Draw("histsame")

    error_band_mc.Draw("e2same")
    signal.Draw("histsame")
    dataSum.Draw("epsame")

    #-------------------- legend ----------------------------

    legend1 = ROOT.TLegend(0.45, 0.8, 0.66, 0.92, "", "brNDC")
    legend1 = style_legend(legend1)

    legend1.AddEntry(dataSum, "Data", "lpe")
    legend1.AddEntry(error_band_mc, "Bkg. Unc.", "f")
#    legend1.AddEntry(h_sumMC, "Expectation", "l")
    legend1.AddEntry(signal, "t#bar{t}H (125) x " + str(signal_scale) , "l")
    legend1.Draw()

    legend2 = ROOT.TLegend(0.73, 0.6, 0.92, 0.92, "", "brNDC")
#        legend2 = ROOT.TLegend(0.72, 0.57, 0.95, 0.79, "", "brNDC")
    legend2 = style_legend(legend2)

    mc_rev = mc.items()
    mc_rev.reverse()
    lmc = dict(mc_rev)

    for lname, lh in lmc.iteritems():
        if not (lname == "TTH125"):
            legend2.AddEntry(lh, proc_names[lname], "f")

    legend2.Draw()

    c[var].cd()

    p2 = {}
    p2[var] = ROOT.TPad("p2","p2", 0, 0.02, 1, 0.18)
    p2[var].SetTopMargin(0.0)
    p2[var].SetGrid();
    p2[var].SetFillStyle(0);
    p2[var].Draw()
    p2[var].cd()


    #-----Draw and style Data/MC points----
    hist_ratio = get_ratio(dataSum, h_sumMC, ymin=0.5, ymax=1.5, ratio_ytitle="Data/MC")
    hist_ratio = style_hist(hist_ratio, is_data=True)
    if var != "numJets" and var != "numBTagM":
        hist_ratio.GetXaxis().SetRangeUser(var_range[0], var_range[1])
        hist_ratio = style_axes(hist_ratio, yTitle="Data/MC", is_ratio=True)
    else:
        hist_ratio = style_axes(hist_ratio, yTitle="Data/MC", is_ratio=True, is_jet_count=True)

    #------Draw and style error band----
    ratio_up = get_ratio(h_sumMC+h_sumMCup, h_sumMC, is_band = True)
    ratio_down = get_ratio(h_sumMC-h_sumMCdown, h_sumMC, is_band = True)
    error_band = get_error_band(ratio_up, ratio_down, 1)
    error_band = style_hist(error_band, color=ROOT.kGreen, is_error_band=True)

    one = error_band.Clone("one")
    one = style_hist(one, line=True)

    hist_ratio.Draw("pe1")
    error_band.DrawCopy("e2same")
    one.Draw("histsame")
    hist_ratio.Draw("pe1same")

    #-----------------------------------

    c[var].cd()

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.038)
    latex.SetTextAlign(31)
    latex.SetTextAlign(11)

    cut = "CMS Preliminary"

    std_txt = cut + "                  #sqrt{s}=8 TeV, L=19.5 fb^{-1}" # (" + reg + ")"

    textlabel = std_txt
    latex.DrawLatex(0.16, 0.972, textlabel)

    c[var].SaveAs("plots/control_" + var + "_" + reg + ".png")
#    c[var].SaveAs("plots/" + var + "_" + reg + ".pdf")
    c[var].Close()

    print "dataSum = " + str(dataSum.Integral())
    print "sum MC = " + str(h_sumMC.Integral())


def get_jet_count_hist(jet_count_init, jet_range, label_name):
    print "jet range = " + str(jet_range)
    jet_count=jet_count_init.Clone("jet_count")
    jet_count.GetXaxis().SetRange( jet_range[0]+1, jet_range[1]+1)

    for i in range(jet_range[0], jet_range[1]+1):
        jet_count.GetXaxis().SetBinLabel(i+1, str(i) + " " + label_name)

    return jet_count
