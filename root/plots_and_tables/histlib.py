
import ROOT, sys, re
from systematics import get_tot_sys
from collections import OrderedDict as dict

proc_names = {
    "TTJetsJJ": "t#bar{t} + jj",
    "TTJetsBJ": "t#bar{t} + b",
    "TTJetsBB": "t#bar{t} + bb",
    "TTV": "t#bar{t}V",
    "DiBoson": "VV",
    "TTH125": "t#bar{t}H (125)",
    "SingleT": "Single top",
    "EWK": "EWK",
}

colors = {
    "TTJetsJJ": 18,
    "TTJetsBJ": 17,
    "TTJetsBB": 16,
    "EWK": ROOT.kGreen+3,
    "DiBoson": ROOT.kYellow+1,
    "TTH125": ROOT.kRed,
    "TTV": 30,
    "SingleT": ROOT.kMagenta
    }

def get_ratio(hist1, hist2, ymin=0., ymax=2, is_band = False, ratio_ytitle = ""):
    """
    hist1 -- numerator
    hist2 -- denominator
    """
    hist_ratio = hist1.Clone()
    hist_ratio.Divide(hist2)

    if is_band:
        print "considering band histo"
        for ibin in range(hist_ratio.GetNbinsX()+1):
            if hist_ratio.GetBinContent(ibin+1) == 0:
                hist_ratio.SetBinContent(ibin+1, 1)

    hist_ratio.SetStats(False)
    hist_ratio.SetMarkerStyle(20)
    hist_ratio.SetMarkerSize(0.35)
    if is_band:
        hist_ratio.SetLineColor(ROOT.kBlack) 
    else:
        hist_ratio.SetMarkerColor(ROOT.kBlack)
        hist_ratio.SetLineColor(ROOT.kBlack)
    hist_ratio.SetMaximum(ymax)
    hist_ratio.SetMinimum(ymin)
    
    xAxis = hist_ratio.GetXaxis()
    yAxis = hist_ratio.GetYaxis()

    yAxis.CenterTitle()
    yAxis.SetTitle(ratio_ytitle)
    yAxis.SetTitleOffset(0.2)
    yAxis.SetTitleSize(0.18)
    yAxis.SetLabelSize(0.15)
    yAxis.SetNdivisions(3)
    
    xAxis.SetLabelSize(0.01)
    xAxis.SetTitleSize(0.15)
    xAxis.SetTitleOffset(0.5)
    xAxis.SetTitle("")
                                                                             
    return hist_ratio

def get_error_band(err_up, err_down, nominal, band_only=True):
    """
    compose histogram with errors assigned to each bin. Use max(err_up, err_down)
    band_only -- draw at errorband, nominal is a line at 1
    """
    if band_only:
        nominal = err_up.Clone("nominal")

    for ibin in range(nominal.GetNbinsX() + 1):
        nominal.SetBinContent(ibin+1, 1)
        nominal.SetBinError(ibin+1, max(err_up.GetBinContent(ibin+1)-1, err_down.GetBinContent(ibin+1)-1 ) )

    return nominal

def stackplot(dataSum, mc, mc_up, mc_down, signal, var, varname="", reg=""):
    """
    data -- sum of data histograms
    mc -- dictionary of mc histograms
    signal -- signal histogram
    var -- variable name
    varname -- for axis label
    """

    for proc in mc:
        mc[proc].SetLineColor(colors[proc])
        mc[proc].SetFillColor(colors[proc])
        mc[proc].SetFillStyle(1001)

    dataSum.SetMarkerColor(1)
    dataSum.SetMarkerStyle(20)
    dataSum.SetMarkerSize(1)

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
    p1[var].SetTicks(1, 1);

    p1[var].SetFillStyle(0);
 
    
#    if var == "numJets" or var == "numBTagM": # Or Hist- == "btag_LR_5j" or hist == "btag_LR_6j" or hist == "btag_LR_4j":

    h_sumMC.SetTitle("")
    h_sumMC.SetStats(False)
    h_sumMC.SetLineWidth(2)
    h_sumMC.SetMaximum(1.3*max(h_sumMC.GetMaximum(), dataSum.GetMaximum()) )
    h_sumMC.SetMinimum(0.)
    h_sumMC.SetLineColor(ROOT.kBlack)
    h_sumMC.SetFillStyle(0)
    h_sumMC.GetXaxis().SetTitle(varname)
    h_sumMC.GetYaxis().SetTitle("")

    if var == "numJets" or var == "numBTagM" or var == "cat_count": # or var == "btag_LR" or var == "MTln" or var == "Mll": #for logscale
        h_sumMC.SetMinimum(1)
        sum.SetMinimum(1)
        mc["TTH125"].SetMinimum(1)
        signal.SetMinimum(1)
        dataSum.SetMinimum(1)
        p1[var].SetLogy()
#        if var == "btag_LR":
#            h_sumMC.GetXaxis().SetRange( 1, 20)
        if var == "numBTagM":
            h_sumMC.SetMaximum(50*ROOT.TMath.Max(h_sumMC.GetMaximum(), dataSum.GetMaximum()) )
        else:
            h_sumMC.SetMaximum(30*ROOT.TMath.Max(h_sumMC.GetMaximum(), dataSum.GetMaximum()) )

#            h_sumMC.GetXaxis().SetRange(1,5)
    
    p1[var].cd()
       
    h_sumMC.Draw("hist")
    sum.Draw("histsame")
    h_sumMC.Draw("histsame")
    signal.Draw("histsame")
    dataSum.Draw("epsame")

    #-------------------- legend ----------------------------

    legend1 = ROOT.TLegend(0.51, 0.8, 0.72, 0.92, "", "brNDC")
    legend1.SetBorderSize(0)
    legend1.SetFillColor(0)
    legend1.AddEntry(dataSum, "Data", "p")
    legend1.AddEntry(h_sumMC, "Expectation", "l")
    legend1.AddEntry(signal, "TTH125 x " + str(50) , "l")
    legend1.Draw()
    
#    if var == "btag_LR":
    legend2 = ROOT.TLegend(0.72, 0.7, 0.95, 0.92, "", "brNDC")


#    else:
#        legend1 = 
#        legend2 = ROOT.TLegend(0.72, 0.57, 0.95, 0.79, "", "brNDC")
        
    legend2.SetBorderSize(0)
    legend2.SetFillColor(0)
    
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
    hist_ratio.SetMarkerColor(1)
    hist_ratio.SetMarkerStyle(20)
    hist_ratio.SetMarkerSize(1)
    hist_ratio.Draw("ep")
    #------Draw and style error band----
    h_sumMCup = get_tot_sys(mc_up)
    h_sumMCdown = get_tot_sys(mc_down)
    ratio_up = get_ratio(h_sumMC+h_sumMCup, h_sumMC, is_band = True)
    ratio_down = get_ratio(h_sumMC-h_sumMCdown, h_sumMC, is_band = True)

    error_band = get_error_band(ratio_up, ratio_down, 1)
    error_band.DrawCopy("histsame")
    error_band.SetFillColor(ROOT.kBlack)
    error_band.SetFillStyle(3004)
    error_band.Draw("e2same")
    #-----------------------------------

    c[var].cd()

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.03)
    latex.SetTextAlign(31)
    latex.SetTextAlign(11)

    cut = "CMS Preliminary"
        
    std_txt = cut + " #sqrt{s}=8 TeV, L=19.04 fb^{-1}" # (" + reg + ")"
    
    textlabel = std_txt
    latex.DrawLatex(0.15, 0.975, textlabel)

    c[var].SaveAs("plots/control_" + var + "_" + reg + ".png")
#    c[var].SaveAs("plots/" + var + "_" + reg + ".pdf")
    c[var].Close()

    print "dataSum = " + str(dataSum.Integral())
    print "sum MC = " + str(h_sumMC.Integral())

    
def get_jet_count_hist(jet_count_init, jet_range, label_name):
    jet_count=jet_count_init.Clone("jet_count")
    jet_count.GetXaxis().SetRange( jet_range[0]+1, jet_range[1])

    for i in range(jet_range[0], jet_range[1]):
        jet_count.GetXaxis().SetBinLabel(i+1, str(i) + " " + label_name)

    return jet_count

