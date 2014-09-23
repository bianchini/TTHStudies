import FWCore.ParameterSet.Config as cms
from submitDataCardMakerFWlite import submitDataCardMakerFWlite_all
from numpy import arange

inpath="/home/bianchi/CMSSW_5_3_3_patch2/src/Bianchi/TTHStudies/root/"
#inpath="/home/bianchi/CMSSW_5_3_3_patch2/src/Bianchi/TTHStudies/root/files/byLLR/Apr23_2014/"

#inpath="../root/ME_trees/"
#inpath="../root/trees/"
prod_ver = "_ntuplizeAll_v3_rec_std"
#prod_ver = "_rec_std"

outdir = "Sep2014_cat2_test" # in ../root/datacards

do_qcd = False

samples_SL = [
    [["TTV"],      1],
    [["SingleT"],  1],
    [["DiBoson"],  1],
    [["TTJetsBB"], 15], #15
    [["TTJetsBJ"], 15], #15
    [["TTJetsJJ"], 15], #15
    [["TTJetsCC"], 15], #15
    [["TTH125"],   1],
    [["EWK"],      2], #2
    [["Run2012_SingleMu"], 1 ],
    [["Run2012_SingleElectron"], 1 ],
    ]

samples_DL = [
    [["TTV"],      1],
    [["SingleT"],  1],
    [["DiBoson"],  1],
    [["TTJetsBB"], 15],
    [["TTJetsBJ"], 15],
    [["TTJetsJJ"], 15],
    [["TTJetsCC"], 15],
    [["TTH125"],   1],
    [["EWK"],      5],
    [["Run2012_SingleMu"], 1 ],
    [["Run2012_DoubleElectron"], 1 ]
    ]


samples_QCD = [
    [["bEnriched"],            1],
    [["BCtoE"],                1],
    ]

if do_qcd:
    inpath="/home/liis/TTH/CMSSW_5_3_3_patch2/src/Bianchi/TTHStudies/root/" #QCD path
    samples_SL = samples_QCD
    samples_DL = samples_QCD

K = 0.0192521
fk = 1.8

Z = 3237.19
fz = 5.5

p_sb_cut = "p_125_all_s/(p_125_all_s + " + str(K*fk) + "*p_125_all_b)"
p_bj_cut = "p_125_all_b_ttbb/(p_125_all_b_ttbb + " + str(Z*fz) + "*p_125_all_b_ttjj)"

cuts_SL = {
#------------------standard preselection------------------
#    "SL_g5jg3t": "(Vtype==2 || Vtype==3) && numJets>=5 && numBTagM>=3", 

#-------------------btag LR ------------------------------
#    "SL_6j": "(Vtype==2 || Vtype==3) && numJets>=6",
#    "SL_5j": "(Vtype==2 || Vtype==3) && numJets==5",
#    "SL_5jg2t": "(Vtype==2 || Vtype==3) && numJets==5 && numBTagM>=2", #loose preselectiona
#    "SL_g6jg2t": "(Vtype==2 || Vtype==3) && numJets>=6 && numBTagM>=2", #loose preselection
    "SL_5jg1t": "(Vtype==2 || Vtype==3) && numJets==5 && numBTagM>=1", #loose preselectiona
    "SL_g6jg1t": "(Vtype==2 || Vtype==3) && numJets>=6 && numBTagM>=1", #loose preselection    
#----------------- jet multiplicities ------------------------
#    "SL_g4jg2t": "(Vtype==2 || Vtype==3) && numJets>=4 && numBTagM>=2",
#    "SL_g4j": "(Vtype==2 || Vtype==3) && numJets>=4", 
#    "SL_4j": "(Vtype==2 || Vtype==3) && numJets>=4",
#------------------ BDT table ---------------------
#    "SL_g6j2t":  "(Vtype==2 || Vtype==3) && numJets>=6 && numBTagM==2",
#    "SL_4j3t":   "(Vtype==2 || Vtype==3) && numJets==4 && numBTagM==3",
#    "SL_5j3t":   "(Vtype==2 || Vtype==3) && numJets==5 && numBTagM==3",
#    "SL_g6j3t": "(Vtype==2 || Vtype==3) && numJets>=6 && numBTagM==3",
#    "SL_4j4t":   "(Vtype==2 || Vtype==3) && numJets==4 && numBTagM==4",
#    "SL_5jg4t":  "(Vtype==2 || Vtype==3) && numJets==5 && numBTagM>=4",
#    "SL_g6jg4t": "(Vtype==2 || Vtype==3) && numJets>=6 && numBTagM>=4",
#------------------ final categories --------------------------
#    "SL_cat1_HP": "( type==0 || (type==3 && flag_type3>0)) && btag_LR>=0.9925", #0.995
#    "SL_cat2_HP": "( type==1 || (type==3 && flag_type3<=0) ) && btag_LR>=0.9925",
#    "SL_cat3_HP": "type==2 && flag_type2<=999 && btag_LR>=0.995",

#    "SL_cat1_LP": "( type==0 || (type==3 && flag_type3>0)) && btag_LR<0.995 && btag_LR>=0.960",
#    "SL_cat2_LP": "( type==1 || (type==3 && flag_type3<=0) ) && btag_LR<0.9925 && btag_LR>=0.960",
#    "SL_cat3_LP": "type==2 && flag_type2<=999 && btag_LR<0.995 && btag_LR>=0.970",
#--------------------------- All ----------------------------
#    "final_SL": "(type==0 || type==1 || type==2 || type == 3) && btag_LR>0.5",

#------------------- Test SL_cat2 excess ----------------------
    #    "SL_cat2_loose": "( type==1 || (type==3 && flag_type3<=0) ) && btag_LR>=0.960" #this is limited by computing ME   
#    "SL_cat2_HP": "( type==1 || (type==3 && flag_type3<=0) ) && btag_LR>=0.9925",
#    "SL_cat2_HP_firstBin": "( type==1 || (type==3 && flag_type3<=0) ) && btag_LR>=0.9925 && " + p_sb_cut + "< 0.2 && " + p_bj_cut + " < 0.2 "
}

cuts_DL = {
    #--------------------------- standard preselection ------------------------
#    "DL_g2jg2t": "(Vtype==0 || Vtype==1 || Vtype==4) && numJets>=2 && numBTagM >= 2",
    
#    "DL_g4j": "(Vtype==0 || Vtype==1 || Vtype==4) && numJets>=4", #btagLR, numBtag
#    "DL_g4j_z": "(Vtype==0 || Vtype==1 || Vtype==4) && numJets>=4", #btagLR, numBtag
#    "DL_g2jg2t_z": "(Vtype==0 || Vtype==1 || Vtype==4) && numJets>=2 && numBTagM >= 2",

    "DL_g4jg1t": "(Vtype==0 || Vtype==1 || Vtype==4) && numJets>=4 && numBTagM >= 1",

    #------------------ BDT table ---------------------------
#    "DL_g4j2t": "(Vtype==0 || Vtype==1 || Vtype==4) && numJets>=4 && numBTagM == 2",
#    "DL_3j2t": "(Vtype==0 || Vtype==1 || Vtype==4) && numJets==3 && numBTagM ==2",
#    "DL_g3jg3t": "(Vtype==0 || Vtype==1 || Vtype==4) && numJets>=3 && numBTagM >=3",

    # ----------------------- ATLAS table ---------------------
#    "DL_3j3t": "(Vtype==0 || Vtype==1 || Vtype==4) && numJets==3 && numBTagM ==3",
#    "DL_g43t": "(Vtype==0 || Vtype==1 || Vtype==4) && numJets>=4 && numBTagM ==3",
#    "DL_g4g4t": "(Vtype==0 || Vtype==1 || Vtype==4) && numJets>=4 && numBTagM >=4",

    #------------------ final table --------------------------
#    "DL_cat6_HP_mm": "Vtype==0 && type==6 && btag_LR>=0.925",
#    "DL_cat6_HP_ee": "Vtype==1 && type==6 && btag_LR>=0.925",
#    "DL_cat6_HP_em": "Vtype==4 && type==6 && btag_LR>=0.925",

#    "DL_cat6_LP_mm": "Vtype==0 && type==6 && btag_LR<0.925 && btag_LR>=0.850",
#    "DL_cat6_LP_ee": "Vtype==1 && type==6 && btag_LR<0.925 && btag_LR>=0.850",
#    "DL_cat6_LP_em": "Vtype==4 && type==6 && btag_LR<0.925 && btag_LR>=0.850",

#    "DL_cat4_HP": "type==6 && btag_LR>=0.925",
#    "DL_cat4_LP": "type==6 && btag_LR<0.925 && btag_LR>=0.850",
#    "final_DL": "type==6 && btag_LR>0.5",

    }


variables = {
#    "MTln": [30, 500, 30],
#    "Mll": [30, 500, 30],
#    "MET_pt": [0, 500, 30],
#    "MET_sumEt": [300, 3500, 30],

#    "nPVs": [0, 40, 40],

#    "lepton_pt": [20, 250, 20], #30
#    "lepton_eta": [-2.5, 2.5, 10], #nbins/2 15
#    "lepton_rIso": [0, 0.12, 20],
#    "lepton_dxy": [0, 0.025, 30],

#    "Vtype":[0,6,6],
#    "numBTagM": [0, 10, 10],
#    "numJets": [0,12,12],

#    "btag_LR": [0.9925, 1., 10], # (minbin, maxbin, nbins)
#    "jetsAboveCut": [0, 10, 10],

#    "bjet_pt": [30, 350, 30],
#    "bjet_eta": [-2.5, 2.5, 15],
#    "leadjet_pt": [30, 350, 30],
#    "leadjet_eta": [-2.5, 2.5, 15],

#    "best_0_pt_H": [0, 700, 30],
#    "best_0_pt_Top": [0, 700, 30],

    #---------- for cat2 checks ----------
#    "numBTagM": [0, 10, 10],
#    "numJets": [0,12,12],

    "btag_LR": [0., 1., 30], # (minbin, maxbin, nbins)
#    "p_sb": [0, 1, 10],
#    "p_bj": [0, 1, 10], 

#    "MET_pt": [0, 300, 10],
#    "MTln": [30, 500, 10],
#    "bjet_pt": [30, 350, 10],
#    "bjet_eta": [-2.5, 2.5, 10],
#    "leadjet_pt": [30, 350, 10],
#    "leadjet_eta": [-2.5, 2.5, 10],
    
#    "lepton_pt": [20, 250, 10], #30
#    "lepton_eta": [-2.5, 2.5, 10], #nbins/2 15
#    "lepton_rIso": [0, 0.12, 10],
#    "lepton_dxy": [0, 0.025, 10],               
    }

discriminant_mapping = {
    "p_sb": p_sb_cut,
    "p_bj": p_bj_cut,
   } 


do_muon = False
do_electron = False

print "Read input files from: " + inpath
print "Version: " + prod_ver
print "Write output to: " + outdir

# 2D histogram for p_sb vs p_bj
#binvec_2dDiscr = arange(0., 1.0001, 0.1)
#submitDataCardMakerFWlite_all(discriminant_mapping["p_sb"] + ":" + discriminant_mapping["p_bj"], "p_sb_vs_p_bj", cuts_SL["SL_cat2_HP"], "p_sb_vs_p_bj_SL_cat2_HP", binvec_2dDiscr, binvec_2dDiscr, 0, sampless=samples_SL, inputpath=inpath, version=prod_ver, outdir=outdir)

for var in variables:
    minbin = variables[var][0]
    maxbin = variables[var][1]
    nbins = variables[var][2]
    step = (maxbin-minbin)*1.0/nbins

    if var[:6] == "lepton" and do_muon:
        varname = "muon" + var[6:]
        var = var + "[0]"
    elif var[:6] == "lepton" and do_electron:
        varname = "electron" + var[6:]
        var = var + "[0]"
    else:
        varname = var

    if var == "p_sb" or var == "p_bj": # replace the discriminant by the corresponding expression
        var = discriminant_mapping[varname]
        print "Replacing discriminant " + varname + " with expression: " + var

    print "Plotting: " + varname + " with " + var

#    print "arange(" + str(minbin) + "," + str(maxbin+step) + "," + str(step) + ")"
    binvec = cms.vdouble( arange(minbin, maxbin + 0.00001, step) )
    print "binvec = " + str(binvec)

    print "Submitting SL jobs... "

    for cut in cuts_SL:
        cutname = cut
        if do_muon:
            cuts_SL[cut] = "(" + cuts_SL[cut] + ") && lepton_type[0]==13"
            if not var[:6] == "lepton":
                cutname = cut + "_muon"
                        
        elif do_electron:
            cuts_SL[cut] = "(" + cuts_SL[cut] + ") && lepton_type[0]==11"
            if not var[:6] == "lepton":
                cutname = cut + "_electron"
                    
        print "Submit jobs for var: " + varname + ", cut = " + cut
        print "-----------------------------------------------"
        submitDataCardMakerFWlite_all( var, varname, cuts_SL[cut], varname + "_" + cutname , binvec, binvec, 0, sampless=samples_SL, inputpath=inpath, version=prod_ver, outdir=outdir)

    print "Submitting DL jobs... "
    for cut in cuts_DL:
        
        if do_muon: # and var[:6] == "lepton":
            cuts_DL[cut] = "(" + cuts_DL[cut] + ") && lepton_type[0]==13"
        elif do_electron: # and var[:6] == "lepton":
            cuts_DL[cut] = "(" + cuts_DL[cut] + ") && lepton_type[0]==11"
                                        
        print "Submit jobs for var: " + varname + ", cut = " + cut
        print "-----------------------------------------------"
        submitDataCardMakerFWlite_all( var, varname, cuts_DL[cut], varname + "_" + cut , binvec, binvec, 1, sampless=samples_DL, inputpath=inpath, version=prod_ver, outdir=outdir)

    print " ...done"


