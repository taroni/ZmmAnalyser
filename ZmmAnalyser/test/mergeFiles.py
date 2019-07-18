import FWCore.ParameterSet.Config as cms
process = cms.Process("mergingFiles")

process.load("FWCore.MessageService.MessageLogger_cfi")

lines=[
'file:/tmp/taroni/ZmmReceraCNoRec0.root',
'file:/tmp/taroni/ZmmReceraCNoRec1.root',
'file:/tmp/taroni/ZmmReceraCNoRec2.root',
'file:/tmp/taroni/ZmmReceraCNoRec3.root'
]
readFiles=cms.untracked.vstring()

readFiles.extend(lines)


process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
inputFiles = readFiles
outputFile = "/tmp/taroni/Zmm_noChange.root"
process.source = cms.Source ("PoolSource", fileNames = inputFiles )                             

process.Out = cms.OutputModule("PoolOutputModule",
        fileName = cms.untracked.string (outputFile)
)
process.end = cms.EndPath(process.Out)
