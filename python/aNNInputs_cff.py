import FWCore.ParameterSet.Types as CfgTypes
import FWCore.ParameterSet.Config as cms



plots_aNNInputs = cms.VPSet(


    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0.),
    xHigh     = cms.double(0.40),
    nBins     = cms.int32(20),
    variable  = cms.string("aplanarity"),
    xTitle    = cms.string("aplanarity"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("aplanarity"),
    cut       = cms.string(""),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(500),
    nBins     = cms.int32(20),
    variable  = cms.string("aveMbTag"),
    xTitle    = cms.string("avg mass(tag,tag)"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("aveMbTag"),
    cut       = cms.string(""),
    logy      = cms.int32(0),
    ),

     cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(500),
    nBins     = cms.int32(20),
    variable  = cms.string("aveMunTag"),
    xTitle    = cms.string("avg mass(untag,untag)"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("aveMunTag"),
    cut       = cms.string(""),
    logy      = cms.int32(0),
    ),

     cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(200),
    xHigh     = cms.double(2000),
    nBins     = cms.int32(20),
    variable  = cms.string("massAll"),
    xTitle    = cms.string("M_{INV}(jets,lepton,MET)"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("massAll"),
    cut       = cms.string(""),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(200),
    nBins     = cms.int32(20),
    variable  = cms.string("pt1"),
    xTitle    = cms.string("jet 1 p_{T}"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("pt1"),
    cut       = cms.string(""),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(200),
    nBins     = cms.int32(20),
    variable  = cms.string("pt2"),
    xTitle    = cms.string("jet 2 p_{T}"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("pt2"),
    cut       = cms.string(""),
    logy      = cms.int32(0),
    ),


    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(200),
    nBins     = cms.int32(20),
    variable  = cms.string("pt3"),
    xTitle    = cms.string("jet 3 p_{T}"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("pt3"),
    cut       = cms.string(""),
    logy      = cms.int32(0),
    ),

     cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(200),
    nBins     = cms.int32(20),
    variable  = cms.string("pt4"),
    xTitle    = cms.string("jet 4 p_{T}"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("pt4"),
    cut       = cms.string(""),
    logy      = cms.int32(0),
    ),

     cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(350),
    nBins     = cms.int32(20),
    variable  = cms.string("vLepton_pt[0]"),
    xTitle    = cms.string("Tight Lepton p_{T}"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("leptPt"),
    cut       = cms.string(""),
    logy      = cms.int32(0),
    ),

     cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(300),
    nBins     = cms.int32(20),
    variable  = cms.string("METtype1p2corr.et"),
    xTitle    = cms.string("MET"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("met"),
    cut       = cms.string(""),
    logy      = cms.int32(0),
    ),

     cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0.50),
    xHigh     = cms.double(4.0),
    nBins     = cms.int32(20),
    variable  = cms.string("minDeltaRBtag"),
    xTitle    = cms.string("min #Delta R(tag,tag)"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("minDeltaRBtag"),
    cut       = cms.string(""),
    logy      = cms.int32(0),
    ),

     cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(1000),
    nBins     = cms.int32(20),
    variable  = cms.string("M3"),
    xTitle    = cms.string("M3(1 tag)"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("M3"),
    cut       = cms.string(""),
    logy      = cms.int32(0),
    ),

     cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(300),
    nBins     = cms.int32(20),
    variable  = cms.string("massLJ"),
    xTitle    = cms.string("Mass(lep,tag)"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("massLJ"),
    cut       = cms.string(""),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(2000),
    nBins     = cms.int32(40),
    variable  = cms.string("sumAllPt"),
    xTitle    = cms.string("Sum of p_{T}(Lepton,Jets,MET)"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("sumAllPt"),
    cut       = cms.string(""),
    logy      = cms.int32(0),
    ),

     cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(-0.20),
    xHigh     = cms.double(1.0),
    nBins     = cms.int32(20),
    variable  = cms.string("h10"),
    xTitle    = cms.string("H1"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("H1"),
    cut       = cms.string(""),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(-0.20),
    xHigh     = cms.double(1.0),
    nBins     = cms.int32(20),
    variable  = cms.string("h20"),
    xTitle    = cms.string("H2"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("H2"),
    cut       = cms.string(""),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(-0.20),
    xHigh     = cms.double(1.0),
    nBins     = cms.int32(20),
    variable  = cms.string("h30"),
    xTitle    = cms.string("H3"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("H3"),
    cut       = cms.string(""),
    logy      = cms.int32(0),
    ),

      cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(-0.20),
    xHigh     = cms.double(1.0),
    nBins     = cms.int32(20),
    variable  = cms.string("h40"),
    xTitle    = cms.string("H4"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("H4"),
    cut       = cms.string(""),
    logy      = cms.int32(0),
    ),


    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(1),
    nBins     = cms.int32(20),
    variable  = cms.string("sphericity"),
    xTitle    = cms.string("sphericity"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("sphericity"),
    cut       = cms.string(""),
    logy      = cms.int32(0),
    ),


    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0.70),
    xHigh     = cms.double(1.05),
    nBins     = cms.int32(20),
    variable  = cms.string("aveCsv"),
    xTitle    = cms.string("avg disc. b Tag"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("aveCsv"),
    cut       = cms.string(""),
    logy      = cms.int32(0),
    ),


    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(0.025),
    nBins     = cms.int32(20),
    variable  = cms.string("varCsv"),
    xTitle    = cms.string("dev from avg disc. b Tag"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("varCsv"),
    cut       = cms.string(""),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0.70),
    xHigh     = cms.double(1.05),
    nBins     = cms.int32(20),
    variable  = cms.string("firstBtag"),
    xTitle    = cms.string("highest CSV"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("firstBtag"),
    cut       = cms.string(""),
    logy      = cms.int32(0),
    ),


    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0.70),
    xHigh     = cms.double(1.05),
    nBins     = cms.int32(20),
    variable  = cms.string("secondBtag"),
    xTitle    = cms.string("second-highest CSV"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("secondBtag"),
    cut       = cms.string(""),
    logy      = cms.int32(0),
    ),


    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0.70),
    xHigh     = cms.double(1.05),
    nBins     = cms.int32(20),
    variable  = cms.string("minCsv"),
    xTitle    = cms.string("lowest CSV"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("lowestCSV"),
    cut       = cms.string(""),
    logy      = cms.int32(0),
    ),


    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0.5),
    xHigh     = cms.double(4.0),
    nBins     = cms.int32(20),
    variable  = cms.string("aveDeltaRbTag"),
    xTitle    = cms.string("avg #Delta R(tag,tag)"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("aveDeltaRbTag"),
    cut       = cms.string(""),
    logy      = cms.int32(0),
    ),


    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(500),
    nBins     = cms.int32(10),
    variable  = cms.string("closestJJbTagMass"),
    xTitle    = cms.string("mass(closest b Tags)"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("closestJJbTagMass"),
    cut       = cms.string(""),
    logy      = cms.int32(0),
    ),

    cms.PSet(
    skip      = cms.bool(False),
    xLow      = cms.double(0),
    xHigh     = cms.double(1000),
    nBins     = cms.int32(10),
    variable  = cms.string("bestHiggsMass"),
    xTitle    = cms.string("Best Higgs Mass, Chi Square Fit"),
    yTitle    = cms.string("Events"),
    histoName = cms.string("bestHiggsMass"),
    cut       = cms.string(""),
    logy      = cms.int32(0),
    ),

  


    )
