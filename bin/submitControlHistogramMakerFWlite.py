import FWCore.ParameterSet.Config as cms
from submitDataCardMakerFWlite import submitDataCardMakerFWlite_all
from numpy import arange

inpath="/home/bianchi/CMSSW_5_3_3_patch2/src/Bianchi/TTHStudies/root/"


outdir = "May23_PAS/control_plots"
prod_ver = "_ntuplizeAll_v3_rec_std"

do_qcd = False

samples_SL = [
    [["TTV"],      1],
    [["SingleT"],  1],
    [["DiBoson"],  1],
    [["TTJetsBB"], 8],
    [["TTJetsBJ"], 8],
    [["TTJetsJJ"], 8],
    [["TTJetsCC"], 8],
    [["TTH125"],   1],
    [["EWK"],      2],
    [["Run2012_SingleMu"], 1 ],
    [["Run2012_SingleElectron"], 1 ],
    ]

samples_DL = [
    [["TTV"],      1],
    [["SingleT"],  1],
    [["DiBoson"],  1],
    [["TTJetsBB"], 8],
    [["TTJetsBJ"], 8],
    [["TTJetsJJ"], 8],
    [["TTJetsCC"], 8],
    [["TTH125"],   1],
    [["EWK"],      5],
    [["Run2012_SingleMu"], 1 ],
    [["Run2012_DoubleElectron"], 1 ]
    ]


samples_QCD = [
    [["bEnriched"],            1],
    [["BCtoE"],                1],

    #[["QCD_Pt-50To150_bEnriched"],           1],
    #[["QCD_Pt-150_bEnriched"],               1],

    #[['QCD_Pt_80_170_BCtoE'],                1],
    #[['QCD_Pt_170_250_BCtoE'],               1],
    #[['QCD_Pt_250_350_BCtoE'],               1],
    #[['QCD_Pt_350_BCtoE'],                   1],

    ]

if do_qcd:
    inpath="/home/liis/TTH/CMSSW_5_3_3_patch2/src/Bianchi/TTHStudies/root/" #QCD path
    samples_SL = samples_QCD
    samples_DL = samples_QCD

cuts_SL = {
#    "SL_g6jg3t": "(Vtype==2 || Vtype==3) && numJets>=6 && numBTagM>=3", #tight preselection
#    "SL_g5jg3t": "(Vtype==2 || Vtype==3) && numJets>=5 && numBTagM>=3", #tight preselection   
#    "SL_g5jg2t_eta15": "(Vtype==2 || Vtype==3) && numJets>=5 && numBTagM>=2 && (lepton_eta[0] >= 1.5 || lepton_eta[0] <= -1.5)", #loose preselection
#------------------standard preselection------------------
#    "SL_g5jg2t": "(Vtype==2 || Vtype==3) && numJets>=5 && numBTagM>=2", #loose preselectiona
#-------------------btag LR ------------------------------
#    "SL_6j": "(Vtype==2 || Vtype==3) && numJets>=6",
#    "SL_5j": "(Vtype==2 || Vtype==3) && numJets==5",
#    "SL_5jg2t": "(Vtype==2 || Vtype==3) && numJets==5 && numBTagM>=2", #loose preselectiona
#    "SL_g6jg2t": "(Vtype==2 || Vtype==3) && numJets>=6 && numBTagM>=2", #loose preselection
#-----------------jet multiplicities------------------------
#    "SL_g4jg2t": "(Vtype==2 || Vtype==3) && numJets>=4 && numBTagM>=2",
#    "SL_g4j": "(Vtype==2 || Vtype==3) && numJets>=4", 
#-----------------------------------------------------------
   
#    "SL_g5jg3t": "(Vtype==2 || Vtype==3) && numJets>=5 && numBTagM>=3", #loose preselection


    
#    "SL_g6j2t":  "(Vtype==2 || Vtype==3) && numJets>=6 && numBTagM==2",
#    "SL_5j3t":   "(Vtype==2 || Vtype==3) && numJets==5 && numBTagM==3",
#    "SL_g6jg4t": "(Vtype==2 || Vtype==3) && numJets>=6 && numBTagM>=4",
#    "SL_5jg4t":  "(Vtype==2 || Vtype==3) && numJets==5 && numBTagM>=4",
#    "SL_4j4t":   "(Vtype==2 || Vtype==3) && numJets==4 && numBTagM==4",
#    "SL_4j3t":   "(Vtype==2 || Vtype==3) && numJets==4 && numBTagM==3",
#    "SL_6j": "(Vtype==2 || Vtype==3) && numJets>=6",
#    "SL_5j": "(Vtype==2 || Vtype==3) && numJets==5",

    "SL_cat1_HP": "( type==0 || (type==3 && flag_type3>0)) && btag_LR>=0.995",
    "SL_cat2_HP": "( type==1 || (type==3 && flag_type3<=0) ) && btag_LR>=0.9925",
    "SL_cat3_HP": "type==2 && flag_type2<=999 && btag_LR>=0.995",

    "SL_cat1_LP": "( type==0 || (type==3 && flag_type3>0)) && btag_LR<0.995 && btag_LR>=0.960",
    "SL_cat2_LP": "( type==1 || (type==3 && flag_type3<=0) ) && btag_LR<0.9925 && btag_LR>=0.960",
    "SL_cat3_LP": "type==2 && flag_type2<=999 && btag_LR<0.995 && btag_LR>=0.970",

}

#analysis = 0 SL
#analysis = 1 DL 

cuts_DL = {
#    "DL_g4j": "(Vtype==0 || Vtype==1 || Vtype==4) && numJets>=4", #btagLR, numBtag
#    "DL_g2jg2t": "(Vtype==0 || Vtype==1 || Vtype==4) && numJets>=2 && numBTagM >= 2",

#    "DL_g3jg2t": "(Vtype==0 || Vtype==1 || Vtype==4) && numJets>=3 && numBTagM >= 2",
#    "DL_g4j2t": "(Vtype==0 || Vtype==1 || Vtype==4) && numJets>=4 && numBTagM == 2",
#    "DL_3j2t": "(Vtype==0 || Vtype==1 || Vtype==4) && numJets==3 && numBTagM ==2",
#    "DL_g3jg3t": "(Vtype==0 || Vtype==1 || Vtype==4) && numJets>=3 && numBTagM >=3",
#    "DL_g4j": "(Vtype==0 || Vtype==1 || Vtype==4) && numJets>=4",

#    "DL_g3jg3t_ee": "(Vtype==1) && numJets>=3 && numBTagM >=3",
#    "DL_g3jg3t_mm": "(Vtype==0) && numJets>=3 && numBTagM >=3",
#    "DL_g3jg3t_em": "(Vtype==4) && numJets>=3 && numBTagM >=3",

    "DL_cat6_HP_mm": "Vtype==0 && type==6 && btag_LR>=0.925",
    "DL_cat6_HP_ee": "Vtype==1 && type==6 && btag_LR>=0.925",
    "DL_cat6_HP_em": "Vtype==4 && type==6 && btag_LR>=0.925",

    "DL_cat6_LP_mm": "Vtype==0 && type==6 && btag_LR<0.925 && btag_LR>=0.850",
    "DL_cat6_LP_ee": "Vtype==1 && type==6 && btag_LR<0.925 && btag_LR>=0.850",
    "DL_cat6_LP_em": "Vtype==4 && type==6 && btag_LR<0.925 && btag_LR>=0.850",

#    "DL_cat4_HP": "type==6 && btag_LR>=0.925",
#    "DL_cat4_LP": "type==6 && btag_LR<0.925 && btag_LR>=0.850",

    }


variables = {
#    "MTln": [30, 500, 30],
#    "Mll": [30, 500, 30],
#    "MET_pt": [0, 500, 30],
#    "MET_sumEt": [300, 3500, 30],

#    "nPVs": [0, 40, 40],

#    "lepton_pt": [20, 250, 30] ,
#    "lepton_eta": [-2.5, 2.5, 15], #nbins/2
#    "lepton_rIso": [0, 0.12, 30],
#    "lepton_dxy": [0, 0.025, 30],

    "Vtype":[0,6,6],
#    "numBTagM": [0, 10, 10],
#    "numJets": [0,12,12],

#    "btag_LR": [0, 1, 30], # (minbin, maxbin, nbins)
#    "jetsAboveCut": [0, 10, 10],
#    "numBTagM": [0, 7, 7], #Fixme -- need to change to allow int values

#    "jet_pt": [30, 350, 30],
#    "jet_eta": [-2.5, 2.5, 15],

    }

do_muon = False
do_electron = False

print "Read input files from: " + inpath
print "Version: " + prod_ver
print "Write output to: " + outdir

for var in variables:
    minbin = variables[var][0]
    maxbin = variables[var][1]
    nbins = variables[var][2]
    step = maxbin*1.0/nbins

    varname = var
    
    if var[:6] == "lepton" and do_muon:
        varname = "muon" + var[6:]
        var = var + "[0]"
    elif var[:6] == "lepton" and do_electron:
        varname = "electron" + var[6:]
        var = var + "[0]"
    elif var[:3] == "jet":
        var = var + "[0]"
        
    if var[:6] == "lepton":
        print "Plotting: " + varname + " with " + var

    binvec = cms.vdouble( arange(minbin, maxbin + step, step) )
    print "binvec = " + str(binvec)

    print "Submitting SL jobs... "
    for cut in cuts_SL:
        if do_muon: # and var[:6] == "lepton":
            cuts_SL[cut] = "(" + cuts_SL[cut] + ") && lepton_type[0]==13"
        elif do_electron: # and var[:6] == "lepton":
            cuts_SL[cut] = "(" + cuts_SL[cut] + ") && lepton_type[0]==11"

        print "Submit jobs for var: " + varname + ", cut = " + cut
        print "-----------------------------------------------"
        submitDataCardMakerFWlite_all( var, varname, cuts_SL[cut], varname + "_" + cut , binvec, 0, sampless=samples_SL, inputpath=inpath, version=prod_ver, outdir=outdir)

    print "Submitting DL jobs... "
    for cut in cuts_DL:
        
        if do_muon: # and var[:6] == "lepton":
            cuts_DL[cut] = "(" + cuts_DL[cut] + ") && lepton_type[0]==13"
        elif do_electron: # and var[:6] == "lepton":
            cuts_DL[cut] = "(" + cuts_DL[cut] + ") && lepton_type[0]==11"
                                        
        print "Submit jobs for var: " + varname + ", cut = " + cut
        print "-----------------------------------------------"
        submitDataCardMakerFWlite_all( var, varname, cuts_DL[cut], varname + "_" + cut , binvec, 1, sampless=samples_DL, inputpath=inpath, version=prod_ver, outdir=outdir)

    print " ...done"


