import ROOT
from collections import OrderedDict as dict
import argparse
import math

from systematics import systematics_list, find_sum_sys
parser = argparse.ArgumentParser()

args = parser.parse_args()

do_hp=True
do_lp=False

# ------------ Open input files ----------------------
inpath = "../datacards/Sep_PAPER/"
infilename_mc = "mlfitMEM_COMB_New_rec_std_sbMu1.root"
infilename_data = "MEM_New_rec_std_sb.root"

print "Opening input file for MC: " + infilename_mc
try:
    f=open(inpath + infilename_mc)
    f.close()
except IOError:
    print "File " + inpath+infilename_mc + " doesn't exist!"
inputfile_mc = ROOT.TFile(inpath+infilename_mc)

print "Opening input file for data: " + infilename_data
try:
    f=open(inpath + infilename_data)
    f.close()
except IOError:
    print "File " + inpath+infilename_data + " doesn't exist!"
inputfile_data = ROOT.TFile(inpath+infilename_data)

standalone = True
Lumi = 19.5

# --------------------------------------------------
cuts_HP = dict()
cuts_HP["SL_cat1_HP"] = "Cat 1 (H)"
cuts_HP["SL_cat2_HP"] = "Cat 2 (H)"
cuts_HP["SL_cat3_HP"] = "Cat 3 (H)"
cuts_HP["DL_cat4_HP"] = "DL (H)"

cuts_LP = dict()
cuts_LP["SL_cat1_LP"] = "Cat 1 (L)"
cuts_LP["SL_cat2_LP"] = "Cat 1 (L)"
cuts_LP["SL_cat3_LP"] = "Cat 3 (L)"
cuts_LP["DL_cat4_LP"] = "DL (L)"

cuts = dict() 
cuts["SL_cat1_HP"] = "Cat 1 (H)"
cuts["SL_cat1_LP"] = "Cat 1 (L)"
cuts["SL_cat2_HP"] = "Cat 2 (H)"
cuts["SL_cat2_LP"] = "Cat 2 (L)"
cuts["SL_cat3_HP"] = "Cat 3 (H)"
cuts["SL_cat3_LP"] = "Cat 3 (L)"
cuts["DL_cat4_HP"] = "DL (H)"
cuts["DL_cat4_LP"] = "DL (L)"

if do_hp:
    cuts = cuts_HP
if do_lp:
    cuts = cuts_LP

histname_map = { # different naming convention in data file
    "11": ["MEM_cat1_H","Cat 1 HP"],
    "12": ["MEM_cat1_L","Cat 1 LP"],
    "13": ["MEM_cat2_H","Cat 2 HP"],
    "14": ["MEM_cat2_L","Cat 2 LP"],
    "15": ["MEM_cat3_H","Cat 3 HP"],
    "16": ["MEM_cat3_L","Cat 3 LP"],
    "21": ["MEM_cat6_H","DL HP"],
    "22": ["MEM_cat6_L","DL LP"],
    }

processes = dict() #filename: [histname, pretty name]
processes["TTH125"] = ["total_signal", "$t\\bar{t}H(\\rightarrow b\\bar{b})$" ]
processes["TTJetsJJ"] = ["ttbar", " $t\\bar{t} + jj$"]
processes["TTJetsCC"] = ["ttbarPlusCCbar", "$t\\bar{t} + cc$"]
processes["TTJetsBJ"] = ["ttbarPlusB", " $t\\bar{t} + bj$"]
processes["TTJetsBB"] = ["ttbarPlusBBbar", " $t\\bar{t} + b\\bar{b}$ "]
processes["TTV"] = ["ttbarV", " $t\\bar{t} + V$ "]
processes["SingleT"] = ["singlet", "single-$t$"]
processes["EWK"] = ["ewk", "$V$ + jets"]
#processes["DiBoson"] = ["diboson", "diboson"]

# ------------------ Get Sum Bkg and Sum Data-------------
sumBkg=dict()
sumBkgErr2=dict()
sumBkgSys2=dict()

signal = dict()
sumData=dict()
    

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


for cut in cuts:
    cutlabel = cuts[cut]
    print " & " + cutlabel,
print '\\\ \\hline'

#-------------------Fill yields in table -------------------------------------

nr_signal = dict()
for proc in processes:
    print processes[proc][1],

    for i in range(1,3): # SL/DL
        for j in range(1,7): # category
            if (i == 1 or (i == 2 and j <=2) ) and ( (not do_hp and not do_lp ) or (do_hp and j%2 != 0) or (do_lp and j%2 == 0) ): # if SL or DL and cat 1,2 + filter for SL/DL/All
                
                histname = "shapes_fit_s/Name" + str(i) + "_Name" + str(j) + "/" + processes[proc][0]
                h = inputfile_mc.Get(histname)

                errVal = ROOT.Double(0)
                try:
                    nr_evts = h.IntegralAndError(0, h.GetNbinsX(), errVal)
                except AttributeError:
                    if proc=="EWK":
                        nr_evts = 0
                    else:
                        print "Histogram: " + histname + " cant be opened" 
                if proc != "EWK": # if nr_evts=0 the histogram is missing
                    nr_evts = h.Integral()

                tbl_str = str( round(nr_evts, 1) ) + " $\pm$ " + str( round( math.sqrt(errVal**2), 1) )
        
                print " & ",
                print tbl_str,

                if proc == "TTH125":
                    nr_signal[ str(i)+str(j) ] = nr_evts
                    
    print "\\\\"
    if proc == "TTH125":
        print "\\hline"

print "\\hline"
    
nr_bkg=dict()
print "\\textbf {Total bkg}",
for i in range(1,3): # SL/DL                                                                                                              
    for j in range(1,7): # category
        if (i == 1 or (i == 2 and j <=2) ) and ( (not do_hp and not do_lp ) or (do_hp and j%2 != 0) or (do_lp and j%2 == 0) ):
            h = inputfile_mc.Get("shapes_fit_s/Name" + str(i) + "_Name" + str(j) + "/total_background" )
            errVal = ROOT.Double(0)
            nr_bkg[str(i)+str(j)] = h.IntegralAndError(0, h.GetNbinsX(), errVal)
            nr_bkg[str(i)+str(j)] = h.Integral()
            tbl_str = "\\textbf{" + str( round(nr_bkg[str(i)+str(j)], 1) ) + " $\pm$ " + str( round( math.sqrt(errVal**2), 1) ) + "} "

            print " & ",
            print tbl_str,

print "\\\\"
print "\\hline"

print "\\textbf{Data}",

nr_data=dict()
for i in range(1,3): # SL/DL
    for j in range(1,7): # category
        if (i == 1 or (i == 2 and j <=2) ) and ( (not do_hp and not do_lp ) or (do_hp and j%2 != 0) or (do_lp and j%2 == 0) ):
            h = inputfile_data.Get(histname_map[str(i)+str(j)][0]+ "/data_obs" )
            errVal = ROOT.Double(0)
            nr_data[str(i)+str(j)]= h.IntegralAndError(0, h.GetNbinsX(), errVal)
            nr_data[str(i)+str(j)] = h.Integral()
            tbl_str = "\\textbf{ " + str( round( nr_data[str(i)+str(j)], 1) ) + " $\pm$ " + str( round( math.sqrt(errVal**2), 1) ) + "} "

            print " & ",
            print tbl_str,

print "\\\\"
print "\\hline"

#print "Data/MC",
#for i in range(1,3): # SL/DL
#    for j in range(1,7): # category  
#        if (i == 1 or (i == 2 and j <=2) )
#            print " & " + str( round( nr_data[str(i)+str(j)]/( nr_bkg[str(i)+str(j)] + nr_signal[str(i)+str(j)] ), 2 ) ),
#print "\\\\"
#print "\\hline"

print "$S/B$",
for reg in nr_signal:
    print " & " + str( round( nr_signal[reg]/nr_bkg[reg]*100 , 2) ),
print "\\\\"
print "\\hline"

print " $S/\sqrt{B}$ ",
for reg in nr_signal:
    print " & "+ str( round( nr_signal[reg]/nr_bkg[reg]**0.5 , 2) ),
print "\\\\"
print "\\hline"

#--------------------------- end table ---------------------------
print """
        \end{tabular}
"""        
if do_hp:
    print  "\caption{Post-fit event yields for high-purity categories}" #, L = " + str(Lumi) + " fb$^{-1}$. }"
    print "\label{tab:cutflow_HP}"
if do_lp:
    print  "\caption{Post-fit event yields for low-purity categories}" #, L = " + str(Lumi) + " fb$^{-1}$. }"
    print "\label{tab:cutflow_LP}"

print "}"

print """   
     \end{center}
        \end{table*}

        """
if standalone:
    print "\end{document}"
