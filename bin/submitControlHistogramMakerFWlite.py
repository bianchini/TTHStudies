import FWCore.ParameterSet.Config as cms
from submitDataCardMakerFWlite import submitDataCardMakerFWlite_all
from numpy import arange

inpath="/home/bianchi/CMSSW_5_3_3_patch2/src/Bianchi/TTHStudies/root/"

samples_mc = [
    [["TTV"],      5],
    [["SingleT"],  1],
    [["DiBoson"],  5],
    [["TTJetsBB"], 20],
    [["TTJetsBJ"], 20],
    [["TTJetsJJ"], 20],
    [["TTH125"],   5],
    [["EWK"],      10],
    ]

samples_data_SL = [
    [["Run2012_SingleMu"],10 ],
    [["Run2012_SingleElectron"],10 ]
    ]

samples_data_DL = [
    [["Run2012_SingleMu"],10 ],
    [["Run2012_DoubleElectron"],10 ]
    ]

cuts_SL = {
    "SL_5j3t":   "(Vtype==2 || Vtype==3) && numJets==5 && numBTagM==3",
    "SL_g6j2t":  "(Vtype==2 || Vtype==3) && numJets>=6 && numBTagM==2",
    "SL_g6jg4t": "(Vtype==2 || Vtype==3) && numJets>=6 && numBTagM>=4",
    "SL_g6j3t":  "(Vtype==2 || Vtype==3) && numJets>=6 && numBTagM==3",
    "SL_5jg4t":  "(Vtype==2 || Vtype==3) && numJets==5 && numBTagM>=4",
    "SL_4j4t":   "(Vtype==2 || Vtype==3) && numJets==4 && numBTagM==4",
    "SL_4j3t":   "(Vtype==2 || Vtype==3) && numJets==4 && numBTagM==3",
#    "SL_6j": "(Vtype==2 || Vtype==3) && numJets>=6",
#    "SL_5j": "(Vtype==2 || Vtype==3) && numJets==5",
}

#analysis = 0 SL
#analysis = 1 DL 

cuts_DL = {
#    "DL_4j": "(Vtype==0 || Vtype==1 || Vtype==4) && numJets>=4",
#    "DL_4j2t": "(Vtype==0 || Vtype==1 || Vtype==4) && numJets>=4 && numBTagM >=2",
#    "DL_4j4t": "(Vtype==0 || Vtype==1 || Vtype==4) && numJets>=4 && numBTagM >=4",
    }


variables = {
    "btag_LR": [0, 1, 10], # (minbin, maxbin, nbins)
#    "jetsAboveCut": [0, 10, 10],
#    "numBTagM": [0, 7, 7], #Fixme -- need to change to allow int values
    }


for var in variables:
    minbin = variables[var][0]
    maxbin = variables[var][1]
    nbins = variables[var][2]
    step = maxbin*1.0/nbins

    binvec = cms.vdouble( arange(minbin, maxbin + step, step) )
    print "binvec = " + str(binvec)

    print "Submitting SL jobs... "
    for cut in cuts_SL:
        print "Submit jobs for var: " + var + ", cut = " + cut
        print "-----------------------------------------------"
        submitDataCardMakerFWlite_all( "mc_" + var + "_" + cut, cuts_SL[cut], var + "_" + cut , binvec, 0, samples_mc, inputpath=inpath)
        submitDataCardMakerFWlite_all( "data_" + var + "_" + cut, cuts_SL[cut], var + "_" + cut , binvec, 0, samples_data_SL, inputpath=inpath)

    print "Submitting DL jobs... "
    for cut in cuts_DL:
        print "Submit jobs for var: " + var + ", cut = " + cut
        print "-----------------------------------------------"
        submitDataCardMakerFWlite_all( "mc_" + var + "_" + cut , cuts_DL[cut], var + "_" + cut , binvec, 1, samples_mc, inputpath=inpath)
        submitDataCardMakerFWlite_all( "data_" + var + "_" + cut, cuts_DL[cut], var + "_" + cut , binvec, 1, samples_data_DL, inputpath=inpath)
    print " ...done"
