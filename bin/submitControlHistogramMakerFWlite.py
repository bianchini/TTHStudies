import FWCore.ParameterSet.Config as cms
from submitDataCardMakerFWlite import submitDataCardMakerFWlite_all

inpath="/home/bianchi/CMSSW_5_3_3_patch2/src/Bianchi/TTHStudies/root/"
samples = [
    #[["TTV"],     5],
    [["SingleT"], 1],
    #[["DiBoson"], 5],
    #[["TTJetsBB"],20],
    #[["TTJetsBJ"],20],
    #[["TTJetsJJ"],20],
    #[["TTH125"],   5],
    #[["EWK"],     10],
    #[["Run2012_SingleMu", "Run2012_SingleElectron"],10 ],
    #[["Run2012_SingleMu", "Run2012_DoubleElectron"],10 ]
    ]

binvec = cms.vdouble(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
submitDataCardMakerFWlite_all( "btag_LR", "(Vtype<=1 || Vtype==4) && numJets>=4", "test" , binvec, +1, samples, inputpath=inpath)

