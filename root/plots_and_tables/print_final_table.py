import ROOT
from collections import OrderedDict as dict
import argparse
import math

from systematics import systematics_list, find_sum_sys
parser = argparse.ArgumentParser()

args = parser.parse_args()

#inpath = "../datacards/Apr28_checks/controlPlots_merged/"
inpath = "../datacards/May26_PAS/tables/"
standalone = True

Lumi = 19.5

version = "MEM_New_ntuplizeAll_v3_rec_std_"
#version = "MEM_New_rec_std_"

#var = "jetsAboveCut"
#var = "numBTagM"
var = "Vtype"
hist = "MEM_" + var

cuts = dict()

cuts["SL_cat1_HP"] = "Cat 1 HP"
cuts["SL_cat2_HP"] = "Cat 2 HP"
cuts["SL_cat3_HP"] = "Cat 3 HP"
#cuts["DL_cat4_HP"] = "Cat 4 HP"

cuts["SL_cat1_LP"] = "Cat 1 LP"
cuts["SL_cat2_LP"] = "Cat 2 LP"
cuts["SL_cat3_LP"] = "Cat 3 LP"
#cuts["DL_cat4_LP"] = "Cat 4 LP"

#if args.mode == "DL":
#    cuts["DL_3j2t"] = "3j 2t"
#    cuts["DL_g4j2t"] = "$\ge$4j 2t"
#    cuts["DL_g3jg3t"] = "$\ge$ 3t"


processes = dict() #filename: [histname, pretty name]
processes["TTH125"] = ["TTH125", "$t\\bar{t}H(\\rightarrow b\\bar{b})$" ]
processes["TTJetsJJ"] = ["TTJetsLF", " $t\\bar{t} + jj$"]
processes["TTJetsCC"] = ["TTJetsLFcc", "$t\\bar{t} + cc$"]
processes["TTJetsBJ"] = ["TTJetsHFb", " $t\\bar{t} + bj$"]
processes["TTJetsBB"] = ["TTJetsHFbb", " $t\\bar{t} + b\\bar{b}$ "]
processes["TTV"] = ["TTV", " $t\\bar{t} + V$ "]
processes["SingleT"] = ["SingleT", "single-$t$"]
processes["EWK"] = ["EWK", "$V$ + jets"]
processes["DiBoson"] = ["DiBoson", "diboson"]

processes["Run2012_SingleElectron"] = ["data_obs", ""]
processes["Run2012_SingleMu"] = ["data_obs", ""]

#if args.mode == "DL":
#processes["Run2012_DoubleElectron"] = ["data_obs", ""]


# ------------------ Get Sum Bkg and Sum Data-------------
sumBkg=dict()
sumBkgErr2=dict()
sumBkgSys2=dict()

signal = dict()
sumData=dict()
    
for reg in cuts:
    sumBkg[reg] = 0
    sumBkgErr2[reg] = 0
    sumBkgSys2[reg] = 0
    sumData[reg] = 0
    signal[reg] = 0
    for proc in processes:

        infile = inpath + version + var + "_" + reg + "_" + proc + ".root"
        inputfile = ROOT.TFile(infile)
        h = inputfile.Get(hist + "/" + processes[proc][0])
        errVal = ROOT.Double(0)
        nr_evts = h.IntegralAndError(0, h.GetNbinsX(), errVal) # buggy with underflow bins
        nr_evts = h.Integral()

        if proc == "TTH125":
            signal[reg] = nr_evts
            continue

        if (reg == "DL_cat4_HP" or reg == "DL_cat4_LP") and (proc == "Run2012_SingleMu" or proc == "Run2012_DoubleElectron"):
#            print "add DL: " + reg + ", proc = " + proc
            sumData[reg] += nr_evts
        elif not (reg == "DL_cat4_HP" or reg == "DL_cat4_LP") and (proc == "Run2012_SingleMu" or proc == "Run2012_SingleElectron"):
 #           print "add SL: " + reg + ", proc = " + proc
            sumData[reg] += nr_evts

        elif proc == "Run2012_SingleMu" or proc == "Run2012_SingleElectron" or proc == "Run2012_DoubleElectron":
            continue
            
        else:
            sumBkg[reg] += nr_evts
            sumBkgErr2[reg] += errVal**2
            sumBkgSys2[reg] += ((find_sum_sys(proc, processes[proc][0], systematics_list, inputfile, hist, "Up")).Integral())**2
            
#-----------------------initialize table-----------------------------------
if standalone:
    print "\documentclass{article}"
    print "\usepackage[landscape]{geometry}"
    print "\\begin{document}"


print """
\\begin{table*}[htbp]
\\begin{center}"""
print "\\footnotesize{"


print "\\begin{tabular}{|",
for it in range(len(cuts) + 1):
    print "c|",
print "}"
print "\\hline"

for cutlabel in cuts.values():
    print " & " + cutlabel,
print '\\\ \\hline'

#-------------------Fill yields in table -------------------------------------
for proc in processes:
    if proc == "Run2012_SingleMu" or proc == "Run2012_SingleElectron" or proc == "Run2012_DoubleElectron":
        continue
 
    print processes[proc][1],

    for reg in cuts:
        infile = inpath + version + var + "_" + reg + "_" + proc + ".root"
        inputfile = ROOT.TFile(infile)
        h = inputfile.Get(hist + "/" + processes[proc][0])
        errVal = ROOT.Double(0)
        nr_evts = h.IntegralAndError(0, h.GetNbinsX(), errVal)
        nr_evts = h.Integral()

        mc_up = find_sum_sys(proc, processes[proc][0], systematics_list, inputfile, hist, "Up")
        sys = mc_up.Integral()

        tbl_str = str( round(nr_evts, 1) ) + " $\pm$ " + str( round( math.sqrt(errVal**2 + sys**2), 1) )
        
        print " & ",
        print tbl_str,
        
    print "\\\\"
    if proc == "TTH125":
        print "\\hline"

print "\\hline"
    
print "\\textbf {Total bkg}",
for reg in sumBkg:
    print " & \\textbf{" + str( round(sumBkg[reg], 1) ) + " $\pm$ " + str(round( math.sqrt(sumBkgSys2[reg] + sumBkgErr2[reg]), 1) ) + "}",

print "\\\\"
print "\\hline"
print "\\textbf{Data}",
for reg in sumData:
    print " & \\textbf{" + str( round(sumData[reg], 1) ) + " $\pm$ " + str(round( sumData[reg]**0.5) ) + "}",
print "\\\\"
print "\\hline"

print "Data/MC",
for reg in sumData:
    print " & " + str( round(sumData[reg]/(sumBkg[reg] + signal[reg]), 2 ) ),
print "\\\\"
print "\\hline"

print "$S/B$",
for reg in sumBkg:
    print " & " + str( round( signal[reg]/sumBkg[reg]*100 , 2) ) + "\\%",
print "\\\\"
print "\\hline"

print " $S/\sqrt{B}$ ",
for reg in sumBkg:
    print " & "+ str( round( signal[reg]/sumBkg[reg]**0.5 , 2) ),
print "\\\\"
print "\\hline"

#--------------------------- end table ---------------------------
print """
        \end{tabular}
"""        
print  "\caption{Final event yields for SL categories, L = " + str(Lumi) + " fb$^{-1}$. }"
print "\label{tab:cutflow_SL}"
print "}"

print """   
     \end{center}
        \end{table*}

        """
if standalone:
    print "\end{document}"
