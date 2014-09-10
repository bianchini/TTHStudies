import ROOT

indir = "./" #location of input histograms
infile = ROOT.TFile(indir + "MEAnalysisNewTEST.root")
tree=infile.Get("tree")


h_gen=ROOT.TH1F("h_gen","h_gen",87,0,700)
h_rec=ROOT.TH1F("h_rec","h_rec",87,0,700)
h_recrec=ROOT.TH1F("h_recrec","h_recrec",87,0,700)

tree.Draw("p4H[0]>>h_gen", "p4H[0]>0")
tree.Draw("ptH_gen_matched>>h_rec", "ptH_gen_matched>0")
tree.Draw("ptH_matched>>h_recrec", "ptH_matched>0")


h_num = h_rec.Clone("h_num")
h_den = h_gen.Clone("h_den")


#
#h_den.Scale(1./h_den.Integral())
#h_num.Scale(1./h_num.Integral())

#h_num.GetNbinsX()

gr = ROOT.TGraphAsymmErrors( h_num.GetNbinsX() )
#gr.SetName( str(i_draw ) )
gr.BayesDivide( h_num, h_den)
gr.GetXaxis().SetTitle("Higgs p_{T}")
gr.GetYaxis().SetTitle("Reconstruction efficiency")
gr.GetYaxis().SetTitleOffset(1.5)
gr.SetTitle("")
gr.SetMarkerColor(ROOT.kBlue)
gr.SetLineColor(ROOT.kBlue)
gr.SetMarkerStyle(20)
gr.SetMarkerSize(1)
#gr.SetMaximum(0.00001)


c = ROOT.TCanvas("c1", "c1", 600, 600)
#c.SetGrid()
c.SetLogy()

h_den.SetMarkerColor(ROOT.kBlue)
h_den.SetLineColor(ROOT.kBlue)
h_den.SetMarkerStyle(20)

h_den.SetTitle("")
h_den.GetXaxis().SetTitle("sim. Higgs p_{T}")
h_den.GetXaxis().SetLabelSize(0.032)
h_den.SetStats(False)

h_num.SetMarkerColor(ROOT.kRed)
h_num.SetLineColor(ROOT.kRed)
h_num.SetMarkerStyle(20)

legend1 = ROOT.TLegend(0.45, 0.8, 0.66, 0.92, "", "brNDC")
legend1.AddEntry(h_den, "Loose preselection", "lpe")
legend1.AddEntry(h_num, "Reconstructed by ME", "lpe")

h_den.Draw("ep")
h_num.Draw("epsame")
h_recrec.Draw("epsame")
#legend1.Draw()

#gr.Draw("AP")



