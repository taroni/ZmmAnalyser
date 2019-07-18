
import FWCore.ParameterSet.Config as cms
process = cms.Process("TestElectrons")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")

from Configuration.AlCa.GlobalTag import GlobalTag
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
process.GlobalTag.toGet = cms.VPSet(
  cms.PSet(record = cms.string('EcalChannelStatusRcd'),
           tag = cms.string('EcalChannelStatus_isolated_deadchannels_data'),
           connect = cms.string("frontier://FrontierPrep/CMS_CONDITIONS")
           ),
)
readFiles=cms.untracked.vstring()
count=0
start=50*0
end=50+50*0

#lines=[]
#for line in open('pierreOutputFiles.txt'):
#    if (count >= start and count < end):
#        lines.append('file:%s' %(line.rstrip('\n')))
#    count+=1

readFiles.extend([
#'file:/eos/cms/store/user/taroni/DoubleMu_sum8gt0_pierre.root'
'file:/eos/cms/store/user/taroni/DoubleMu_sum8gt20_pierre.root'
])

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
inputFiles = readFiles
outputFile = "zMuMu_PierreTag_sum8gt20.root"
process.source = cms.Source ("PoolSource", fileNames = inputFiles )                             
process.ntupler = cms.EDAnalyzer(
    'ZmmAnalyser',
    muons    = cms.untracked.InputTag("muons"),
    inputRecHitsEB = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    METInput = cms.InputTag("pfMet")

    )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string( outputFile )
                                   )

process.p = cms.Path(process.ntupler)


