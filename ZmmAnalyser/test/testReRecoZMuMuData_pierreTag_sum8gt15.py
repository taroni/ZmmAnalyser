# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step3 --conditions 101X_dataRun2_Prompt_v9 -n -1 --era Run2_2018 --eventcontent RAWRECO --data -s RAW2DIGI,RECO --datatier RAW-RECO --python testReRecoZSkimData_fromRawReco.py --filein /store/data/Run2018B/EGamma/RAW-RECO/ZElectron-PromptReco-v1/000/317/864/00000/8A838407-EC71-E811-8525-FA163E54B47A.root --fileout file:step3_2018.root --scenario pp --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('RERECO',eras.Run2_2018)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)
process.MessageLogger.cerr.FwkReport.reportEvery = 10
# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(#'/store/data/Run2018B/EGamma/RAW-RECO/ZElectron-PromptReco-v1/000/317/182/00000/08C8D1B6-7966-E811-B307-FA163EAF8008.root',
'/store/data/Run2018C/DoubleMuon/RAW/v1/000/320/191/00000/F86BB976-838F-E811-A261-FA163EC28637.root',
'/store/data/Run2018C/DoubleMuon/RAW/v1/000/320/191/00000/F2983C9D-858F-E811-B444-FA163E50AC90.root',
'/store/data/Run2018C/DoubleMuon/RAW/v1/000/320/191/00000/E25C0DEC-848F-E811-8C79-FA163EFD5EF6.root',
'/store/data/Run2018C/DoubleMuon/RAW/v1/000/320/191/00000/DAA328DB-818F-E811-8E8A-FA163E7B5F86.root',
'/store/data/Run2018C/DoubleMuon/RAW/v1/000/320/191/00000/D861536B-838F-E811-BA02-FA163EF52E44.root',
'/store/data/Run2018C/DoubleMuon/RAW/v1/000/320/191/00000/D60E2668-AD90-E811-9B1E-FA163EE14499.root',
'/store/data/Run2018C/DoubleMuon/RAW/v1/000/320/191/00000/CEC320D3-828F-E811-AB62-FA163EA6A331.root',
'/store/data/Run2018C/DoubleMuon/RAW/v1/000/320/191/00000/C8AABDA0-828F-E811-8F12-FA163E95B0AD.root',
'/store/data/Run2018C/DoubleMuon/RAW/v1/000/320/191/00000/BC7AB5E2-848F-E811-AB05-FA163E9D506A.root',
'/store/data/Run2018C/DoubleMuon/RAW/v1/000/320/191/00000/B4AE9FA1-828F-E811-A8A7-FA163EAE298E.root',
'/store/data/Run2018C/DoubleMuon/RAW/v1/000/320/191/00000/B499E828-818F-E811-90A8-FA163E10C293.root',
'/store/data/Run2018C/DoubleMuon/RAW/v1/000/320/191/00000/AE499799-858F-E811-8B0F-FA163EE1C1A2.root'
#        '/store/data/Run2018B/DoubleMuon/RAW/v1/000/319/310/00000/8E350449-8381-E811-BF47-FA163EE23257.root'
#        '/store/data/Run2018C/EGamma/RAW-RECO/ZElectron-PromptReco-v3/000/320/369/00000/E2C382D4-7192-E811-945E-FA163E6A4783.root'
),
    secondaryFileNames = cms.untracked.vstring()
)
#events= cms.untracked.VEventRange()
#for line in open('zMuMu2018B.txt'):
#    events.append(line.rstrip("\n"))
#process.source.eventsToProcess = cms.untracked.VEventRange(events)

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '101X_dataRun2_Prompt_v9', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data_promptlike', '')

process.GlobalTag.toGet = cms.VPSet(
  cms.PSet(record = cms.string('EcalChannelStatusRcd'),
           tag = cms.string('EcalChannelStatus_isolated_deadchannels_data'),
           connect = cms.string("frontier://FrontierPrep/CMS_CONDITIONS")
           ),
)

process.options = cms.untracked.PSet(

)
process.options.numberOfThreads=cms.untracked.uint32(4)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step3 nevts:-1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

from Configuration.EventContent.EventContent_cff import RAWRECOEventContent
process.skimContent = process.RAWRECOEventContent.clone()
process.load("DPGAnalysis.Skims.filterRecHitsRecovery_cfi")
process.recoveryfilter = cms.Path(process.recHitRecoveryFilter)
from RecoLocalCalo.EcalRecProducers.ecalRecHit_cfi import ecalRecHit
process.ecalRecHit.singleChannelRecoveryThreshold=0.7
process.ecalRecHit.sum8ChannelRecoveryThreshold=15.

# Output definition

process.RAWRECOoutput = cms.OutputModule("PoolOutputModule",
                             
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('RAW-RECO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('file:step3_2018_sum8gt15.root'),
    outputCommands = process.RAWRECOEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring('recoveryfilter'))
)
process.RAWRECOoutput.outputCommands.append('drop *_*_*_RECO')


process.ntupler = cms.EDAnalyzer(
    'ZmmAnalyser',
    muons    = cms.untracked.InputTag("muons"),
    inputRecHitsEB = cms.InputTag("ecalRecHit","EcalRecHitsEB"),
    METInput = cms.InputTag("pfMet")
    )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("ZmmTree_2018Cv2_sum8gt15.root")
                                   )

process.load("DPGAnalysis/Skims/ZMuSkim_cff") 
process.p = cms.Path(process.recoverySequence*process.diMuonSelSeq*process.ntupler)

# Path and EndPath definitions

process.raw2digi_step = cms.Path(process.RawToDigi)
process.reconstruction_step = cms.Path(process.reconstruction)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWRECOoutput_step = cms.EndPath(process.RAWRECOoutput)

##process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.reconstruction_step,  process.recoveryfilter, process.p, process.RAWRECOoutput_step, process.endjob_step)
#from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
#associatePatAlgosToolsTask(process)


### Customisation from command line
##
###Have logErrorHarvester wait for the same EDProducers to finish as those providing data for the OutputModule
##from FWCore.Modules.logErrorHarvester_cff import customiseLogErrorHarvesterUsingOutputCommands
##process = customiseLogErrorHarvesterUsingOutputCommands(process)
##
### Add early deletion of temporary data products to reduce peak memory need
##from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
##process = customiseEarlyDelete(process)
### End adding early deletion
