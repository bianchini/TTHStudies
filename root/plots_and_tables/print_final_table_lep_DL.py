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


cuts["DL_cat6_HP_mm"] = "$\mu\mu$"
cuts["DL_cat6_HP_em"] = "$e\mu$"
cuts["DL_cat6_HP_ee"] = "ee"

cuts["DL_cat6_LP_mm"] = "$\mu\mu$"
cuts["DL_cat6_LP_em"] = "$e\mu$"
cuts["DL_cat6_LP_ee"] = "ee"

cutlabels = ["di-lepton HP", "di-lepton LP"]

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

#processes["Run2012_SingleElectron"] = ["data_obs", ""]
processes["Run2012_SingleMu"] = ["data_obs", ""]
processes["Run2012_DoubleElectron"] = ["data_obs", ""]


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

        if (reg[:7] == "DL_cat6") and (proc == "Run2012_SingleMu" or proc == "Run2012_DoubleElectron"):
            print "add DL: " + reg + ", proc = " + proc
            sumData[reg] += nr_evts
        elif not (reg[:7] == "DL_cat6") and (proc == "Run2012_SingleMu" or proc == "Run2012_SingleElectron"):
            print "add SL: " + reg + ", proc = " + proc
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

for cutlabel in cutlabels:
    print " & \multicolumn{3}{|c|}{" + cutlabel + "}",
print '\\\\'
print "\\cline{2-" + str(len(cuts) + 1) + "}"

for cutlabel in cutlabels:
    print "& $\mu\mu$ &  $e\mu$ & $ee$", 
print '\\\ \\hline'  

#-------------------Fill yields in table -------------------------------------
for proc in processes:
    if proc == "Run2012_SingleMu" or proc == "Run2012_SingleElectron" or proc == "Run2012_DoubleElectron":
        continue
 
    print processes[proc][1],
    
    if proc == "TTH125":
        sumSignal_HP = 0
        sumSignal_LP = 0
        sumSignalErr2_HP = 0
        sumSignalErr2_LP = 0
    for reg in cuts:
        infile = inpath + version + var + "_" + reg + "_" + proc + ".root"
        inputfile = ROOT.TFile(infile)
        h = inputfile.Get(hist + "/" + processes[proc][0])
        errVal = ROOT.Double(0)
        nr_evts = h.IntegralAndError(0, h.GetNbinsX(), errVal)
        nr_evts = h.Integral()

        mc_up = find_sum_sys(proc, processes[proc][0], systematics_list, inputfile, hist, "Up")
        sys = mc_up.Integral()
        toterr = math.sqrt(errVal**2 + sys**2)

        if proc == "TTH125":
            round_prec = 2
            if reg[8:][:2] == "HP":
                sumSignal_HP+=nr_evts
                sumSignalErr2_HP+=toterr**2
            if reg[8:][:2] == "LP":
                sumSignal_LP+=nr_evts
                sumSignalErr2_LP+=toterr**2
        else:
            round_prec = 1
        tbl_str = str( round(nr_evts, round_prec) ) + " $\pm$ " + str( round( toterr, 2) )
        
        print " & ",
        print tbl_str,

    print "\\\\"
    if proc == "TTH125":
        print "\\cline{2-" + str(len(cuts) + 1) + "}"
        for idx in range( len(cutlabels)):
            if idx == 0:
                sumSignal = sumSignal_HP
                sumSignalErr = math.sqrt(sumSignalErr2_HP)
            elif idx ==1:
                sumSignal = sumSignal_LP
                sumSignalErr = math.sqrt(sumSignalErr2_LP)
            else:
                sumSignal = 0
                sumSignalErr = 0

            print " & \multicolumn{3}{|c|}{" + str(round(sumSignal, round_prec)) + " $\pm$ " + str( round( sumSignalErr, 2) ) + "}",
        print "\\\\",
        print "\\hline" 

print "\\hline"
    
print "\\textbf {Total bkg}",

sumBkgTot_HP=0
sumBkgTot_LP=0
sumBkgErrTot2_HP=0
sumBkgErrTot2_LP=0
for reg in sumBkg:
    sumBkgErr = math.sqrt(sumBkgSys2[reg] + sumBkgErr2[reg])
    print " & \\textbf{" + str( round(sumBkg[reg], 1) ) + " $\pm$ " + str(round( sumBkgErr, 1) ) + "}",

    if reg[8:][:2] == "HP":
        sumBkgTot_HP += sumBkg[reg]
        sumBkgErrTot2_HP += sumBkgErr**2
    elif reg[8:][:2] == "LP":
        sumBkgTot_LP += sumBkg[reg]
        sumBkgErrTot2_LP += sumBkgErr**2

print "\\\\"
print "\\cline{2-" + str(len(cuts) + 1) + "}"

for idx in range( len(cutlabels)):
    if idx == 0:
        sumBkgTot = sumBkgTot_HP
        sumBkgErrTot = math.sqrt(sumBkgErrTot2_HP)
    elif idx ==1:
        sumBkgTot = sumBkgTot_LP
        sumBkgErrTot = math.sqrt(sumBkgErrTot2_LP)
    else:
        sumBkgTot = 0
        sumBkgErrTot = 0
    print " & \multicolumn{3}{|c|}{\\textbf{" + str(round(sumBkgTot, 1)) + " $\pm$ " + str( round( sumBkgErrTot, 1) ) + "} }",
print"\\\\",
print "\\hline"


sumDataTot_HP=0
sumDataTot_LP=0
print "\\textbf{Data}",
for reg in sumData:
    print " & \\textbf{" + str( round(sumData[reg], 1) ) + " $\pm$ " + str(round( sumData[reg]**0.5) ) + "}",
    if reg[8:][:2] == "HP":
        sumDataTot_HP += sumData[reg]
    elif reg[8:][:2] == "LP":
        sumDataTot_LP += sumData[reg]
print "\\\\"
print "\\cline{2-" + str(len(cuts) + 1) + "}"

for idx in range( len(cutlabels)):
    if idx == 0:
        sumDataTot = sumDataTot_HP
    elif idx ==1:
        sumDataTot = sumDataTot_LP
    else:
        sumDataTot = 0
    print " & \multicolumn{3}{|c|}{\\textbf{" + str(round(sumDataTot, 1)) + " $\pm$ " + str( round( sumDataTot**0.5, 1) ) + "} }",


print "\\\\"
print "\\hline"
print "\\hline"
print "Data/MC",
for reg in sumData:
    print " & " + str( round(sumData[reg]/(sumBkg[reg] + signal[reg]), 2 ) ),

print "\\\\"
print "\\cline{2-" + str(len(cuts) + 1) + "}"
for idx in range( len(cutlabels)):
    if idx == 0:
        print " & \multicolumn{3}{|c|}{ " + str( round(sumDataTot_HP/(sumBkgTot_HP + sumSignal_HP ), 2 ) ) + "}",
    elif idx ==1:
        print " & \multicolumn{3}{|c|}{ " + str( round(sumDataTot_LP/(sumBkgTot_LP + sumSignal_LP ), 2 ) ) + "}",

print "\\\\"
print "\\hline"
print "\\hline"

print "$S/B$",
for reg in sumBkg:
    print " & " + str( round( signal[reg]/sumBkg[reg]*100 , 2) ) + "\\%",

print "\\\\"
print "\\cline{2-" + str(len(cuts) + 1) + "}"

for idx in range( len(cutlabels)):
    if idx == 0:
        print " & \multicolumn{3}{|c|}{ " + str( round(sumSignal_HP/(sumBkgTot_HP)*100, 2 ) ) + "\\%}",
    elif idx ==1:
        print " & \multicolumn{3}{|c|}{ " + str( round(sumSignal_LP/(sumBkgTot_LP)*100, 2 ) ) + "\\%}",

print "\\\\"
print "\\hline"
print "\\hline"

print " $S/\sqrt{B}$ ",
for reg in sumBkg:
    print " & "+ str( round( signal[reg]/sumBkg[reg]**0.5 , 2) ),

print "\\\\"
print "\\cline{2-" + str(len(cuts) + 1) + "}"

for idx in range( len(cutlabels)):
    if idx == 0:
        print " & \multicolumn{3}{|c|}{ " + str( round(sumSignal_HP/(sumBkgTot_HP**0.5), 2 ) ) + "}",
    elif idx ==1:
        print " & \multicolumn{3}{|c|}{ " + str( round(sumSignal_LP/(sumBkgTot_LP**0.5), 2 ) ) + "}",

print "\\\\"
print "\\hline"

#--------------------------- end table ---------------------------
print """
        \end{tabular}
"""        
print "\caption{Final event yields for DL categories, L = " + str(Lumi) + " fb$^{-1}$. }"
print "\label{tab:cutflow_DL}"
print "}"

print """   
     \end{center}
        \end{table*}

        """
if standalone:
    print "\end{document}"
