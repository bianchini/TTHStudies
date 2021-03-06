import ROOT, sys, re
import CMS_lumi
from systematics import get_tot_sys
from collections import OrderedDict as dict

proc_names = {
    "TTJetsJJ": "t#bar{t} + lf",
    "TTJetsCC": "t#bar{t} + c#bar{c}",
    "TTJetsBJ": "t#bar{t} + b",
    "TTJetsBB": "t#bar{t} + b#bar{b}",
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
    "SingleT": ROOT.kMagenta,
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

def style_hist(hist, color=0, is_data = False, is_signal=False, is_error_band=False, line=False, yRange=[0.5,1.5]):
  if is_data:
    hist.SetMarkerStyle(20)
    hist.SetMarkerColor(ROOT.kBlack)
    hist.SetLineColor(ROOT.kBlack)

    hist.SetMarkerSize(1.5)
    hist.SetLineWidth(2)
#    hist.SetBinErrorOption(ROOT.TH1F.kPoisson)

  elif is_signal:
    hist.SetLineColor(colors["TTH125"])
    hist.SetFillStyle(0)
#    signal.SetLineStyle(ROOT.kDashed)
    
    hist.SetLineWidth(4)
#    hist.SetLineWidth(10)

  elif is_error_band:
    hist.SetMarkerSize(0)
    hist.SetMinimum(yRange[0])
    hist.SetMaximum(yRange[1])
    if color != 0:
      hist.SetFillColor(color)
      hist.SetFillStyle(1001)
    else:
      hist.SetFillColor(ROOT.kBlack)
      hist.SetFillStyle(3654) #3004

  elif line:
    hist.SetLineColor(ROOT.kBlack)
#    hist.SetFillStyle(0)

    hist.SetLineStyle(ROOT.kDashed)
#    hist.SetLineStyle(9)
    hist.SetLineWidth(2) 

#    hist.SetLineWidth(1)

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

def add_cms_info(lumi=19.5, com=8):
  latex = ROOT.TLatex()
  latex.SetNDC()
  latex.SetTextSize(0.038)
  latex.SetTextAlign(31)
  latex.SetTextAlign(11)
  
  cut = "CMS Preliminary"

  std_txt = cut + "                  #sqrt{s}=" + str(com) + " TeV, L=" + str(lumi) + " fb^{-1}" # (" + reg + ")"                      

  textlabel = std_txt
  latex.DrawLatex(0.16, 0.96, textlabel)

def get_poisson_err(hist):
    poissonErr = ROOT.TGraphAsymmErrors(hist)

    alpha = 1-0.6827
    for i in range(0,poissonErr.GetN()):
      N = poissonErr.GetY()[i]

      if N != 0:
        L = ROOT.Math.gamma_quantile(alpha/2,N,1.)
      else:
        L = 0

      U = ROOT.Math.gamma_quantile_c(alpha/2,N+1,1.)

      poissonErr.SetPointEYlow(i, N-L)
      poissonErr.SetPointEYhigh(i, U-N)
      
#      xPoint = hist.GetBinCenter(i+1)
#      if hist.GetBinContent(i+1) < 2:
#        poissonErr.SetPoint(i, xPoint, 0)
#      print "low = " + str(N-L)
#      print "high = " + str(U-N)

    return poissonErr

def get_poisson_ratio(hist1_poisson,hist1,hist2):
    gr = ROOT.TGraphAsymmErrors( hist1_poisson.GetN() )

    for iBin in range(0,hist1.GetNbinsX()):
      xPoint = hist1.GetBinCenter(iBin+1)
      xWidth = 0.5*hist1.GetBinWidth(iBin+1)

      yG = hist1_poisson.GetY()[iBin]
      yG_low = hist1_poisson.GetEYlow()[iBin]
      yG_high = hist1_poisson.GetEYhigh()[iBin]
      y1 = hist1.GetBinContent(iBin+1)
      y2 = hist2.GetBinContent(iBin+1)

      if y2 > 0:
        yG_ratio = yG/y2
        yG_ratio_low = yG_low/y2
        yG_ratio_high = yG_high/y2

      if y1 > 0:
        gr.SetPoint(iBin, xPoint, yG_ratio)
        gr.SetPointEYlow(iBin, yG_ratio_low)
        gr.SetPointEYhigh(iBin, yG_ratio_high)

        gr.SetPointEXlow(iBin, xWidth)
        gr.SetPointEXhigh(iBin, xWidth)

    return gr

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

def stackplot(dataSum, mc, mc_up, mc_down, signal, var, varname="", var_range=[-500,500], reg="", outdir="plots", plot_style = "pas", legend_header=""):
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
    dataSumPoisson = get_poisson_err(dataSum)

    sum = ROOT.THStack("sum","")

    for proc, mc_hist in mc.iteritems():
        if not proc == "TTH125": #FIXME
            sum.Add(mc_hist)

    h_sumMC = mc["TTV"].Clone("h_sumMC") #FIXME

    for proc in mc:
        if not proc=="TTH125" and not proc=="TTV":
            h_sumMC.Add(mc[proc])

    c={}
#    c[var] = ROOT.TCanvas("c" + var ,"c" + var, 3*800, 3*1000)
    c[var] = ROOT.TCanvas("c" + var, "c" + var,  5, 30, 640, int(580/0.85) ) # Daniel sizing
    
    p1 = {}
#    p1[var] = ROOT.TPad("p1", "p1", 0, 0.25, 1, 1)

    p1[var] = ROOT.TPad("p1", "p1", 0, 0.155, 1, 1) #Daniel
#    p1[var].SetBottomMargin(0)


    p1[var].SetGrid(0,0)
    p1[var].SetFillStyle(4000)
    p1[var].SetFillColor(10)
    p1[var].SetTicky()
    p1[var].SetTicks(0,0)
    p1[var].SetObjectStat(0)
    p1[var].Draw()
    p1[var].cd()

#    p1[var].SetTicks(0, 0)
#    p1[var].SetFillStyle(0)
#    p1.SetObjectStat(0)
#    p1[var].Draw()

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



    if var == "Mll" or var == "numJets" or var == "numBTagM" or var == "cat_count" or var == "btag_LR" or var == "MTln" or var == "Mll" or var == "MET_pt" or var == "jetsAboveCut": #for logscale
        if var == "numJets" or var == "numBTagM":
          ymin_log = 1
        elif var == "btag_LR":
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
        elif var == "btag_LR":
          h_sumMC.SetMaximum(70*ROOT.TMath.Max(h_sumMC.GetMaximum(), dataSum.GetMaximum()) )
        else:
          h_sumMC.SetMaximum(20*ROOT.TMath.Max(h_sumMC.GetMaximum(), dataSum.GetMaximum()) )


    h_sumMCup = get_tot_sys(mc_up)
    h_sumMCdown = get_tot_sys(mc_down)

    error_band_mc = get_error_band(h_sumMCup, h_sumMCdown, h_sumMC, 0)
    error_band_mc = style_hist(error_band_mc, is_error_band=True)


#    p1[var].cd()

    h_sumMC.Draw("hist")
    sum.Draw("histsame")
    h_sumMC.Draw("histsame")

    error_band_mc.Draw("e2same")
    signal.Draw("histsame")
#    dataSum.Draw("epsame")
    dataSum.Draw("pe1same")
#    dataSumPoisson.Draw("epsame")

    #-------------------- legend ----------------------------

    legend1 = ROOT.TLegend(0.45, 0.76, 0.66, 0.92, "", "brNDC")
    legend1 = style_legend(legend1)
    legend1.SetHeader(legend_header)

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

    if plot_style == "pas":
      print "Labels old/costum style"
      add_cms_info(lumi=19.5, com=8)
    else: # according to new standard
      print "Labels following latest CMS recommendations"
      plot_style="new" # "pas" for old styling
      CMS_lumi.writeExtraText = 0
      CMS_lumi.lumi_8TeV = "19.5 fb^{-1}"
      CMS_lumi.lumiTextSize = 0.75
      CMS_lumi.cmsTextSize = 0.9
      CMS_lumi.CMS_lumi(p1[var], 2, 11) # 2 -- 8TeV only, 

    c[var].cd()

    p2 = {}
#    p2[var] = ROOT.TPad("p2","p2", 0, 0.02, 1, 0.18)
    p2[var] = ROOT.TPad("p2","p2", 0, 0, 1, 0.175)
#    p2[var].SetTopMargin(0.0)
#    p2[var].SetGrid();
    p2[var].SetFillStyle(0);
    p2[var].Draw()
    p2[var].cd()


    #-----Draw and style Data/MC points----
    hist_ratio = get_ratio(dataSum, h_sumMC, ymin=0.5, ymax=1.5, ratio_ytitle="Data/Bkg")
    hist_ratio = style_hist(hist_ratio, is_data=True)
    if var != "numJets" and var != "numBTagM":
        hist_ratio.GetXaxis().SetRangeUser(var_range[0], var_range[1])
        hist_ratio = style_axes(hist_ratio, yTitle="Data/Bkg", is_ratio=True)
    else:
        hist_ratio = style_axes(hist_ratio, yTitle="Data/Bkg", is_ratio=True, is_jet_count=True)

    hist_ratio_poisson = get_poisson_ratio(dataSumPoisson, dataSum, h_sumMC)
    hist_ratio_poisson = style_hist(hist_ratio_poisson, is_data=True, is_error_band=True)

    #------Draw and style error band----
    ratio_up = get_ratio(h_sumMC+h_sumMCup, h_sumMC, is_band = True)
    ratio_down = get_ratio(h_sumMC-h_sumMCdown, h_sumMC, is_band = True)
    error_band = get_error_band(ratio_up, ratio_down, 1)
    error_band = style_hist(error_band, color=ROOT.kGreen, is_error_band=True)
    error_band = style_axes(error_band, is_ratio=True)

    hist_ratio.Draw("p")
    error_band.Draw("e2same")
    hist_ratio_poisson.Draw("pe1same")

    c[var].Update()
#    one = ROOT.TLine(c[var].GetUxmin(),1,c[var].GetUxmax(),1)
    print "xmin = " + str(c[var].GetUxmin()) + ", xmax = " + str(c[var].GetUxmax())
    print var_range[0], var_range[1] #hist_ratio.GetXaxis().GetXmin()

    one = ROOT.TLine(var_range[0],1,var_range[1],1)
    one = style_hist(one, line=True)

    one.Draw("same")
#    hist_ratio.Draw("pe1same")

    #-----------------------------------

    c[var].cd()


#    latex = ROOT.TLatex()
#    latex.SetNDC()
#    latex.SetTextSize(0.038)
#    latex.SetTextAlign(31)
#    latex.SetTextAlign(11)

#    cut = "CMS Preliminary"

#    std_txt = cut + "                  #sqrt{s}=8 TeV, L=19.5 fb^{-1}" # (" + reg + ")"

#    textlabel = std_txt
#    latex.DrawLatex(0.16, 0.972, textlabel)

    outfile = outdir + "/control_" + var + "_" + reg + ".png"
    print "saving outpit file to: " + outfile
    c[var].SaveAs( outfile ) 
    c[var].SaveAs( outdir + "/control_" + var + "_" + reg + ".root" )
    c[var].SaveAs( outdir + "/control_" + var + "_" + reg + ".pdf")
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
