import ROOT
from collections import OrderedDict as dict

indir = "../datacards/June03_PAS/data_trees/"

root_files = dict()
root_files["SingleMuon"] = ROOT.TFile(indir + "MEAnalysisNew_all_rec_std_Run2012_SingleMu.root")
root_files["SingleElectron"] = ROOT.TFile(indir + "MEAnalysisNew_all_rec_std_Run2012_SingleElectron.root")
root_files["DiElectron"] = ROOT.TFile(indir + "MEAnalysisNew_all_rec_std_Run2012_DoubleElectron.root")



trees = {}
for file in root_files:
    trees[file] = root_files[file].Get("tree")


cut_bin_67 = dict() # value at lower edge of bin
cut_bin_67["bin7"] = 0.835
cut_bin_67["bin6"] = 0.667

cut_bin_7 = dict()
cut_bin_7["bin7"] = 0.835

cuts = dict()
cuts["SL_cat1"] = "( type==0 || (type==3 && flag_type3>0)) && btag_LR>=0.995"
cuts["SL_cat2"] =  "( type==1 || (type==3 && flag_type3<=0) ) && btag_LR>=0.9925"
cuts["SL_cat3"] = "type==2 && flag_type2<=999 && btag_LR>=0.995"
cuts["DL"] = "type==6 && btag_LR>=0.925"

sEl_trigger = "(Vtype==3 && ( triggerFlags[44]>0 ) )"
sMu_trigger = "(Vtype==2 && ( triggerFlags[22]>0 || triggerFlags[23]>0 || triggerFlags[47]>0 ))"
dMu_trigger = "( (Vtype==0 || Vtype==4) && ( triggerFlags[22]>0 || triggerFlags[23]>0  ))"
dEl_trigger = "(Vtype==1 && ( triggerFlags[6]>0 ) )"
Zll_veto = "( (Vtype<=1 && abs(Mll-91.2)>8. && Mll>15.) || Vtype > 1)"

scan = "Vtype:btag_LR: numJets: numBTagM: jet_csv[2]: jet_csv[5]: jet_csv[6]: jet_csv[7]: jet_csv[3]: jet_csv[4]: jet_eta[2]: jet_eta[5]: jet_eta[6]: jet_eta[7]: jet_eta[3]: jet_eta[4]: lepton_pt[0]: lepton_eta[0]: MET_pt: Mll: MTln"

precision=3

for cut in cuts:
    print "---------------------------------------------"
    print "Starting: " + cut
    if cut == "DL" or cut == "SL_cat1":
        ps_cuts = cut_bin_67
    else:
        ps_cuts = cut_bin_7
        
    if cut == "SL_cat1":
        discriminant = "p_125_all_s_ttbb/(p_125_all_s_ttbb+1.200000*(0.005916*p_125_all_b_ttbb+85.645287*p_125_all_b_ttjj))"
    if cut == "SL_cat2":
        discriminant = "p_125_all_s_ttbb/(p_125_all_s_ttbb+0.600000*(0.013812*p_125_all_b_ttbb+68.346115*p_125_all_b_ttjj))"
    if cut == "SL_cat3":
        discriminant = "p_125_all_s_ttbb/(p_125_all_s_ttbb+0.600000*(0.012640*p_125_all_b_ttbb+64.279221*p_125_all_b_ttjj))"
    if cut == "DL":
        discriminant = "p_125_all_s_ttbb/(p_125_all_s_ttbb+0.600000*(0.014056*p_125_all_b_ttbb+6.107426*p_125_all_b_ttjj))"


    for ps_cut in ps_cuts:
        if ps_cut == "bin7": # last bin
            cutstring = cuts[cut] + "&& (" + discriminant + " > " + str(ps_cuts[ps_cut]) + ")"
        else: # second last bin, add higher threshold
            cutstring = cuts[cut] + "&& (" + discriminant + " > " + str(ps_cuts[ps_cut]) + "&&" + discriminant + " < " + str(ps_cuts["bin7"]) + ")"

        for tree in trees:
            if tree == "SingleElectron" and cut != "DL":
                print "Starting: " + tree + " in " + ps_cut
                trees[tree].Scan(scan, cutstring + " &&" + sEl_trigger, "precision=" + str(precision) )
            elif tree == "SingleMuon" and cut != "DL":
                print "Starting: " + tree + " in " + ps_cut
                trees[tree].Scan(scan, cutstring+ " &&" + sMu_trigger, "precision=" + str(precision) )
            elif tree == "SingleMuon" and cut == "DL":
                print "Starting: " + tree + " in " + ps_cut
                trees[tree].Scan(scan, cutstring + " &&" + dMu_trigger + "&&" +Zll_veto, "precision=" + str(precision) )
            elif tree =="DiElectron" and cut =="DL":
                print "Starting: " + tree + " in " + ps_cut
                trees[tree].Scan(scan, cutstring + " &&" + dEl_trigger + "&&" + Zll_veto, "precision=" + str(precision) )
