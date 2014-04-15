import FWCore.ParameterSet.Config as cms


ADDZLLVETO    = 1

ADDPIXELVETO  = 1

ADDDIJETPTCUT = 1

ADDJETPT40CUT = 0

# baseline
cat = cms.PSet(
    
    name      = cms.string("New"),
    version   = cms.string("_ntuplizeAll_v2_rec_std"),
    extraname = cms.string(""),
    fname     = cms.string("MEM"),
    inputpath = cms.string("../root/"),
    directory = cms.string("Apr07_2014"),
    cut       = cms.string(""),
    category  = cms.string(""),
    doMEM     = cms.int32(4),
    fact1     = cms.double(0),
    fact2     = cms.double(0),
    factbb    = cms.double(0),
    lumiScale = cms.double(19.04/12.1),
    nBins     = cms.int32(6),
    splitFirstBin = cms.int32(0),
    binvec        = cms.vdouble(),
    samples       = cms.vstring( "TTV", "SingleT", "DiBoson", "TTJetsBB",
                                 "TTJetsBJ", "TTJetsJJ", "TTH125", "EWK",
                                 "Run2012_SingleMu", "Run2012_SingleElectron"),
    nparts= cms.int32(1),
    part  = cms.int32(0),

    # 0 = SL, 1 = DL
    analysis = cms.untracked.int32(-1)

    )


#################### ttbb vs ttjj discrimination

cat1_bj = cat.clone(
    extraname = cms.string("_bj"),
    cut       = cms.string("( type==0 || (type==3 && flag_type3>0)) && btag_LR>=0."), #0.975
    category  = cms.string("cat1"),
    doMEM     = cms.int32(3),
    factbb    = cms.double(0.15)
    )

cat2_bj = cat.clone(
    extraname = cms.string("_bj"),
    cut       = cms.string("( type==1 || (type==3 && flag_type3<=0) ) && btag_LR>=0."),#0.975
    category  = cms.string("cat2"),
    doMEM     = cms.int32(3),
    factbb    = cms.double(0.20)
    )

cat3_bj = cat.clone(
    extraname = cms.string("_bj"),
    cut       = cms.string("type==2 && flag_type2<=999 && btag_LR>=0."), #0.990
    category  = cms.string("cat3"),
    doMEM     = cms.int32(3),
    factbb    = cms.double(0.50)
    )

cat6_bj = cat.clone(
    extraname = cms.string("_bj"),
    cut       = cms.string("type==6 && btag_LR>=0."),#0.980
    category  = cms.string("cat6"),
    factbb    = cms.double(0.15),
    doMEM     = cms.int32(3),
    samples   = cms.vstring( "TTV", "SingleT", "DiBoson", "TTJetsBB",
                             "TTJetsBJ", "TTJetsJJ", "TTH125", "EWK",
                             "Run2012_SingleMu", "Run2012_DoubleElectron")
    ) 

#################### ttbb vs ttH separation (no bias)

cat1_sb_nb =  cat1_bj.clone(
    cut       = cms.string("( type==0 || (type==3 && flag_type3>0)) && btag_LR>=0.975"),
    extraname     = cms.string("_sb_nb"),
    doMEM         = cms.int32(-2),
    fact1         = cms.double(1.2),
    splitFirstBin = cms.int32(1),
    )

cat2_sb_nb =  cat2_bj.clone(
    cut       = cms.string("( type==1 || (type==3 && flag_type3<=0) ) && btag_LR>=0.975"),
    extraname = cms.string("_sb_nb"),
    doMEM     = cms.int32(-2),
    fact1     = cms.double(0.6),
    splitFirstBin = cms.int32(1),
    )

cat3_sb_nb =  cat3_bj.clone(
    cut       = cms.string("type==2 && flag_type2<=999 && btag_LR>=0.990"),
    extraname = cms.string("_sb_nb"),   
    doMEM     = cms.int32(-2),
    fact1     = cms.double(0.6),
    splitFirstBin = cms.int32(1),   
    )

cat6_sb_nb =  cat6_bj.clone(
    cut       = cms.string("type==6 && btag_LR>=0.925"),
    extraname = cms.string("_sb_nb"),
    doMEM     = cms.int32(-2),
    fact1     = cms.double(0.6),
    splitFirstBin = cms.int32(1),    
    )


#################### ttbb vs ttH separation

cat1_sb =  cat1_sb_nb.clone(
    cut       = cms.string("( type==0 || (type==3 && flag_type3>0)) && btag_LR>=0.0"),
    extraname = cms.string("_sb"),
    doMEM     = cms.int32(2),
    )

cat2_sb =  cat2_sb_nb.clone(
    cut       = cms.string("( type==1 || (type==3 && flag_type3<=0) ) && btag_LR>=0.0"),
    extraname = cms.string("_sb"),
    doMEM     = cms.int32(2),
    )

cat3_sb =  cat3_sb_nb.clone(
    cut       = cms.string("type==2 && flag_type2<=999 && btag_LR>=0.0"),
    extraname = cms.string("_sb"),
    doMEM     = cms.int32(2),
    )

cat6_sb =  cat6_sb_nb.clone(
    cut       = cms.string("type==6 && btag_LR>=0.0"),
    extraname = cms.string("_sb"),
    doMEM     = cms.int32(2),
    )

#################### ttbb vs ttH separation (no bias) !!! new !!!

cat1_sb_H =  cat1_bj.clone(
    cut           = cms.string("( type==0 || (type==3 && flag_type3>0)) && btag_LR>=0.995"),
    category      = cms.string("cat1_H"),
    extraname     = cms.string("_sb"),
    doMEM         = cms.int32(2),
    fact1         = cms.double(1.2),
    splitFirstBin = cms.int32(1),
    )


cat1_sb_L =  cat1_bj.clone(
    cut           = cms.string("( type==0 || (type==3 && flag_type3>0)) && btag_LR<0.995 && btag_LR>=0.960"),
    category      = cms.string("cat1_L"),
    extraname     = cms.string("_sb"),
    doMEM         = cms.int32(2),
    fact1         = cms.double(1.2),
    splitFirstBin = cms.int32(0),
    )


cat2_sb_H =  cat2_bj.clone(
    cut       = cms.string("( type==1 || (type==3 && flag_type3<=0) ) && btag_LR>=0.9925"),
    category  = cms.string("cat2_H"), 
    extraname = cms.string("_sb"),
    doMEM     = cms.int32(2),
    fact1     = cms.double(0.6),
    splitFirstBin = cms.int32(1),
    )

cat2_sb_L =  cat2_bj.clone(
    cut       = cms.string("( type==1 || (type==3 && flag_type3<=0) ) && btag_LR<0.9925 && btag_LR>=0.960"),
    category  = cms.string("cat2_L"), 
    extraname = cms.string("_sb"),
    doMEM     = cms.int32(2),
    fact1     = cms.double(1.8),
    splitFirstBin = cms.int32(0),
    )


cat3_sb_H =  cat3_bj.clone(
    cut       = cms.string("type==2 && flag_type2<=999 && btag_LR>=0.995"), 
    category  = cms.string("cat3_H"), 
    extraname = cms.string("_sb"),
    doMEM     = cms.int32(2),
    fact1     = cms.double(0.6),
    splitFirstBin = cms.int32(1),   
    )

cat3_sb_L =  cat3_bj.clone(
    cut       = cms.string("type==2 && flag_type2<=999 && btag_LR<0.995 && btag_LR>=0.970"), 
    category  = cms.string("cat3_L"), 
    extraname = cms.string("_sb"),
    doMEM     = cms.int32(2),
    fact1     = cms.double(2.0),
    splitFirstBin = cms.int32(0),   
    )

cat6_sb_H =  cat6_bj.clone(
    cut       = cms.string("type==6 && btag_LR>=0.925"),
    category  = cms.string("cat6_H"), 
    extraname = cms.string("_sb"),
    doMEM     = cms.int32(2),
    fact1     = cms.double(0.6),
    splitFirstBin = cms.int32(1),
    nBins     = cms.int32(6),
    )

cat6_sb_L =  cat6_bj.clone(
    cut       = cms.string("type==6 && btag_LR<0.925 && btag_LR>=0.850"),
    category  = cms.string("cat6_L"), 
    extraname = cms.string("_sb"),
    doMEM     = cms.int32(2),
    fact1     = cms.double(2.0),
    splitFirstBin = cms.int32(0),
    nBins     = cms.int32(6),
    )
