
import FWCore.ParameterSet.Config as cms
process = cms.Process("TestElectrons")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

readFiles=cms.untracked.vstring()
count=0
start=50*0
end=50+50*0

readFiles.extend([
'file:/eos/cms/store/user/taroni/DoubleMuAOD/PrasannaEvtList/pickevents_merged.root',
'file:/eos/cms/store/user/taroni/DoubleMuAOD/PrasannaEvtList/pickevents_merged001.root'
])

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10) )
inputFiles = readFiles
outputFile = "zMuMu_prasanna_noTagNoRecov.root"
process.source = cms.Source ("PoolSource", fileNames = inputFiles )                             
process.ntupler = cms.EDAnalyzer(
    'ZmmAnalyserXtalList',
    muons    = cms.untracked.InputTag("muons"),
    #inputRecHitsEB = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    inputRecHitsEB = cms.InputTag("ecalRecHit","reducedEcalRecHitsEB"),
    METInput = cms.InputTag("pfMet")

    )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string( outputFile )
                                   )

process.p = cms.Path(process.ntupler)


