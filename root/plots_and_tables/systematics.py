import math
import ROOT
csv_sys_list = ["CMS_ttH_CSVLF", "CMS_ttH_CSVHFStats1", "CMS_ttH_CSVHFStats2", "CMS_ttH_CSVCErr1", "CMS_ttH_CSVCErr2",                                                                                                                         
                "CMS_ttH_CSVHF",  "CMS_ttH_CSVLFStats1", "CMS_ttH_CSVLFStats2"]

q2_sys_list = ["Q2scale_ttH_ttbar1p", "Q2scale_ttH_ttbar2p", "Q2scale_ttH_ttbar3p", "Q2scale_ttH_ttbar_bb", "Q2scale_ttH_ttbar_b", "Q2scale_ttH_ttbar_cc"]

pdf_sys_list = ["pdf_gg", "pdf_qq", "pdf_qg"]

extras = ["CMS_ttH_topPtcorr"] 

normalization = ["CMS_ttH_eff_lep", "lumi_8TeV", "Norm_TTbb", "Norm_TTb"]

systematics_list_profile = ["CMS_scale_j", "CMS_res_j"] + csv_sys_list + q2_sys_list + extras 
systematics_list = ["CMS_scale_j", "CMS_res_j"] + csv_sys_list + q2_sys_list + extras #+ pdf_sys_list + normalization 

#systematics_list = ["JEC", "JER", "lepton_id", "lumi"] + csv_sys_list + q2_sys_list + pdf_sys_list + extras

def find_sum_sys(proc, procname, systematics_list, infile, hist, shift):
#    print "Getting systematics from" + str(infile)
    nominal = infile.Get(hist + "/" + procname)

    nominal_m = nominal.Clone("nominal_m")
    nominal_m.Scale(-1)
#    print nominal.Integral()

    if procname == "TTH125":
        qcdScale = "QCDscale_TTH"
    else:
        qcdScale = "QCDscale_" + procname
    systematics_list.append(qcdScale)
    
    for idx, sys in enumerate(systematics_list):
#        print "Getting systematic variation " + sys
        sys_var_up = infile.Get(hist + "/" + procname + "_" + sys + shift)

        try:
            sys_var_up.Add(nominal_m)
        except AttributeError:
            print "WARNING: systematic variation '" + sys + shift + "' can not be read for " + procname
            continue

        sys_var_up2 = sys_var_up*sys_var_up
        
        if idx == 0:
            sum_sys_up2 = sys_var_up2.Clone("sum_sys_up2")
        else:
            sum_sys_up2.Add(sys_var_up2)


                           
    sys_tot_up = sum_sys_up2.Clone("sys_tot_up")
    for ibin in range(sum_sys_up2.GetNbinsX()+1): #loop over all bins and set each bin value to sqrt(err2) (cant sqrt(TH1F))
        isys_sum = sum_sys_up2.GetBinContent(ibin+1)
        sys_tot_up.SetBinContent(ibin+1, math.sqrt(isys_sum) )

    if qcdScale in systematics_list: # ask only for the sample being processed
        systematics_list.remove(qcdScale)

 #   print "Sumsys : " + proc + ": " + str(sys_tot_up.Integral())
    return sys_tot_up

def get_tot_sys( sys_proc):
    """
    sys_proc -- dictionary of histograms of systematic variations per mc sample
    """

    idx = 0
    for proc in sys_proc:
        sys_proc2 = sys_proc[proc]*sys_proc[proc]
        if idx == 0:
            tot_sys2 = sys_proc2.Clone("tot_sys2")
        else:
            tot_sys2.Add(sys_proc2)
        idx+=1

    tot_sys = tot_sys2.Clone("tot_sys")
    for ibin in range(tot_sys2.GetNbinsX()+1): #loop over all bins and set each bin value to sqrt(err2) 
        tot_sys_bin = tot_sys2.GetBinContent(ibin+1)
        tot_sys.SetBinContent(ibin+1, math.sqrt(tot_sys_bin) )

#    print "Sumsys (tot processes): " + str(tot_sys.Integral())

    return tot_sys

def get_profile_sys( mc_nominal, mc_sys):
    """
    mc_nominal -- nominal histogram
    mc_syst -- dictionary of systematic variations of mc_nominal
    """
#    print mc_nominal.Integral()
    mc_nominal_m = mc_nominal.Clone("mc_nominal_m")
    mc_nominal_m.Scale(-1)
#    print mc_nominal_m.Integral()

    for idx, sys in enumerate(mc_sys):
        mc_diff = mc_sys[sys].Clone("mc_diff")
        mc_diff.Add(mc_nominal_m)
        mc_diff2= mc_diff*mc_diff

        if idx == 0:
            sum_sys2 = mc_diff2.Clone("sum_sys2")
        else:
            sum_sys2.Add(mc_diff2)

        mc_diff.Clear()
        mc_diff2.Clear()
        

#    print "sys2 = " + str(sum_sys2.Integral())
    sum_sys = sum_sys2.Clone("sum_sys")
    for ibinx in range(sum_sys2.GetNbinsX()+1):
        for ibiny in range(sum_sys2.GetNbinsY()+1):
            isys_sum = sum_sys2.GetBinContent(ibinx+1, ibiny+1)
            sum_sys.SetBinContent(ibinx+1, ibiny+1, math.sqrt(isys_sum))

 #   print "final sys = " + str(sum_sys.Integral())
 
    return sum_sys

def get_profile_sys2( mc_nominal, mc_sys, xrange, doProfileY):
    print "---------get_profile2 ----------"
    
    mc_nominal_2 = mc_nominal.Clone("mc_nominal_2")

    if doProfileY:
        prof_nominal_2 = mc_nominal_2.ProfileY()
    else:
        prof_nominal_2 = mc_nominal_2.ProfileX()


    prof_sys = {}
    for sys in mc_sys:
        if doProfileY:
            prof_sys[sys] = mc_sys[sys].ProfileY()
        else:
            prof_sys[sys] = mc_sys[sys].ProfileX()

#    print range
    sys_tot = 0

    sys_band = ROOT.TH1F("sys_band","sys_band", prof_nominal_2.GetNbinsX(), xrange[0], xrange[1])
    for ibinx in range(prof_nominal_2.GetNbinsX()+1):
        print ibinx
        mean_nominal = prof_nominal_2.GetBinContent(ibinx)
        print "nominal = " + str( mean_nominal )

        for idx, sys in enumerate(mc_sys):

            mean_sys = prof_sys[sys].GetBinContent(ibinx) 
            var_sys2 = ( mean_nominal - mean_sys)**2
            
            if idx == 0:
                sum_sys2 = var_sys2
            else:
                sum_sys2 += var_sys2

            print sys + ": " + str( mean_sys )
        
        sys_tot = math.sqrt(sum_sys2)
        print "sys_tot = " + str( sys_tot )
        
#        prof_nominal_2.SetBinError(ibinx, math.sqrt( prof_nominal_2.GetBinError(ibinx)**2 + sys_tot**2) )
#        prof_nominal_2.SetBinError(ibinx, math.sqrt( sys_tot**2) )

#        print "stat error in bin= " + str(prof_nominal_2.GetBinError(ibinx))

 #       if ( prof_nominal_2.GetBinError(ibinx) != 0): # happens when only 1 evt in bin
        sys_band.SetBinContent(ibinx, prof_nominal_2.GetBinContent(ibinx))
        sys_band.SetBinError(ibinx, math.sqrt( prof_nominal_2.GetBinError(ibinx)**2 + sys_tot**2) )
#        else:
#            sys_band.SetBinContent(ibinx, 0)

    print "------------end getprofile2 ------------"
#    return prof_nominal_2
    return sys_band

        
