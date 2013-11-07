import FWCore.ParameterSet.Config as cms

process = cms.Process("ANALYZE")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source(
    "PoolSource",
    fileNames = cms.untracked.vstring(
    'file:../../../VHbbAnalysis/HbbAnalyzer/test/test53.root'
    )
    )


process.bhadrons = cms.EDProducer(
    "MCBHadronProducer",
    quarkId = cms.uint32(5)
    )

process.LHEInfoProducer = cms.EDAnalyzer("LHEInfoProducer")

process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
    )

process.p = cms.Path( process.bhadrons * process.LHEInfoProducer )

process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("LHEInfos.root")
    )

#process.out = cms.OutputModule(
#    "PoolOutputModule",
#    outputCommands = cms.untracked.vstring( 'drop *'),
#    fileName = cms.untracked.string('LHEEvent.root'),
#)

#process.out.SelectEvents = cms.untracked.PSet(
#    SelectEvents = cms.vstring('p')
#    )
#process.outpath = cms.EndPath(process.out)
