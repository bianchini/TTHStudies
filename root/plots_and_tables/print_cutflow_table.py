import ROOT
from collections import OrderedDict as dict

inpath = "../datacards/Apr07_2014_merged/"
do_bdt_comparison = True
standalone = True

jetCut = "jet $p_{T} > 30$ GeV"
version = "MEM_New_ntuplizeAll_v2_rec_std_"

var = "btag_LR"
hist = "MEM_" + var

cuts = dict()
cuts["SL_g6j2t"] =  "$\ge$6j 2t"
cuts["SL_4j3t"] =  "4j 3t"
cuts["SL_5j3t"] = "5j 3t"
cuts["SL_g6j3t"] = "$\ge$6j 3t"
cuts["SL_4j4t"] = "4j 4t"
cuts["SL_5jg4t"] = "5j $\ge$ 4t"
cuts["SL_g6jg4t"] = "$\ge$6j $\ge$ 4t"

processes = dict() #filename: [histname, pretty name]
processes["TTH125"] = ["TTH125", "" ]
processes["TTJetsJJ"] = ["TTJetsLF", " $t\\bar{t} + jj$"]
processes["TTJetsBJ"] = ["TTJetsHFb", " $t\\bar{t} + bj$"]
processes["TTJetsBB"] = ["TTJetsHFbb", " $t\\bar{t} + b\\bar{b}$ "]
processes["TTV"] = ["TTV", " $t\\bar{t} + V$ "]
processes["SingleT"] = ["SingleT", "single-$t$"]
processes["EWK"] = ["EWK", "$V$ + jets"]
processes["DiBoson"] = ["DiBoson", "diboson"]

processes["Run2012_SingleElectron"] = ["data_obs", ""]
processes["Run2012_SingleMu"] = ["data_obs", ""]

# ------------------ Get Sum Bkg and Sum Data-------------
sumBkg=dict()
sumBkgErr2=dict()

signal = dict()
sumData=dict()
    
for reg in cuts:
    sumBkg[reg] = 0
    sumBkgErr2[reg] = 0
    sumData[reg] = 0
    signal[reg] = 0
    for proc in processes:

        infile = inpath + version + var + "_" + reg + "_" + proc + ".root"
        inputfile = ROOT.TFile(infile)
        h = inputfile.Get(hist + "/" + processes[proc][0])
        errVal = ROOT.Double(0)
        nr_evts = h.IntegralAndError(0, h.GetNbinsX(), errVal)

        if proc == "TTH125":
            signal[reg] = nr_evts
            continue

        if proc == "Run2012_SingleMu" or proc == "Run2012_SingleElectron":
            sumData[reg] +=nr_evts
        else:
            sumBkg[reg] += nr_evts
            sumBkgErr2[reg] += errVal**2


bdtData=[10724, 5667, 3983, 2426, 122, 219, 260]
bdtSignal = [33.4, 14., 21.1, 23.1, 1.8, 5.2, 8.3]

#-----------------------initialize table-----------------------------------
if standalone:
    print "\documentclass{article}"
    print "\usepackage[landscape]{geometry}"
    print "\\begin{document}"

print """
\\begin{table*}[htbp]
\\begin{center}"""

print "\label{tab:cutflow}"
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
    if proc == "Run2012_SingleMu" or proc == "Run2012_SingleElectron":
        continue
 
    print processes[proc][1],

    for reg in cuts:
        infile = inpath + version + var + "_" + reg + "_" + proc + ".root"
        inputfile = ROOT.TFile(infile)
        h = inputfile.Get(hist + "/" + processes[proc][0])
        errVal = ROOT.Double(0)
        nr_evts = h.IntegralAndError(0, h.GetNbinsX(), errVal)

        tbl_str = str( round(nr_evts, 1) ) + " $\pm$ " + str( round( errVal, 1) )
        
        print " & ",
        print tbl_str,
        
    print "\\\\"
    if proc == "TTH125":
        print "\\hline"

print "\\hline"
    
print "Total bkg",
for reg in sumBkg:
    print " & " + str( round(sumBkg[reg], 1) ) + " $\pm$ " + str(round( sumBkgErr2[reg]**0.5, 1) ),

print "\\\\"
print "\\hline"
print "Data",
for reg in sumData:
    print " & " + str( round(sumData[reg], 1) ) + " $\pm$ " + str(round( sumData[reg]**0.5) ),
print "\\\\"
print "\\hline"

#---------------------------------------------------------------------
if do_bdt_comparison:
    print "\\hline"
    print "BDT Data",

    for it_bdt in range(len(sumData)):
        print " & " + str( bdtData[it_bdt] ),

    print "\\\\"
    print "MEM/BDT",

    it_bdt=0
    for reg in sumData:
        print " & " + str( round(sumData[reg]/bdtData[it_bdt],1) ),
        it_bdt+=1;
    print "\\\\"
    print "\\hline"
    
    print "\\hline"
    
    print "BDT signal",

    for it_bdt in range(len(sumData)):
        print " & " + str( bdtSignal[it_bdt] ),

    print "\\\\"
    print "MEM/BDT",
    it_bdt=0
    for reg in sumData:
        print " & " + str( round( signal[reg]/bdtSignal[it_bdt],1) ),
        it_bdt+=1;

    print "\\\\"
    print "\\hline"
    

#--------------------------- end table ---------------------------
print """
        \end{tabular}
"""        
print "\caption{Cut flow (" + jetCut + ") }"

print """   
     \end{center}
        \end{table*}

        """
if standalone:
    print "\end{document}"
