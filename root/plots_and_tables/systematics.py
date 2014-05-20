import math
csv_sys_list = ["csv_sys_CSVLF", "csv_sys_CSVHFStats1", "csv_sys_CSVHFStats2", "csv_sys_CSVCErr1", "csv_sys_CSVCErr2",
                "csv_sys_CSVHF",  "csv_sys_CSVLFStats1", "csv_sys_CSVLFStats2"]

q2_sys_list = ["Q2Scale1p", "Q2Scale2p", "Q2Scale3p", "Q2ScaleHFb", "Q2ScaleLF", "Q2ScaleHFbb"]#, "pdf_gg"]#, "pdf_qq", "pdf_qg"]

pdf_sys_list = ["pdf_gg", "pdf_qq", "pdf_qg"]

extras = ["TopPt", "Norm_TTbb", "Norm_TTb"]

systematics_list = ["JEC", "JER", "lepton_id", "lumi"] + csv_sys_list + q2_sys_list + pdf_sys_list + extras

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

    print "Sumsys (tot processes): " + str(tot_sys.Integral())

    return tot_sys
