#from CRABClient.UserUtilities import config
#config = config()
#the two following lines are due to https://hypernews.cern.ch/HyperNews/CMS/get/computing-tools/3731/1/1/1/1/1/2/1/4.html
from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'filteredZMuMu_pierreTag_Cv3_sum8gt15'
config.General.transferLogs = False

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
# Name of the CMSSW configuration file
config.JobType.psetName = 'testReRecoZMuMuData_pierreTag_sum8gt15.py'
#config.JobType.maxMemoryMB=2500
config.JobType.maxMemoryMB=16000
#config.JobType.maxJobRuntimeMin=10080
config.JobType.numCores=4

##config.JobType.inputFiles=['ZeeIntEvents.txt']

config.section_("Data")
config.Data.inputDataset = '/DoubleMuon/Run2018C-v1/RAW'
config.Data.splitting = 'Automatic'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob=40
config.Data.totalUnits=-1
config.Data.publication = True
#config.Data.outLFNDirBase='/store/group/dpg_ecal/comm_ecal/taroni/'

# This string is used to construct the output dataset name
config.Data.outputDatasetTag = 'ZMuMu_pierreTag_Cv3_sum8gt15'

# These values only make sense for processing data
#    Select input data based on a lumi mask

##ADD THE DCS ONLY LUMIMASK
#    Select input data based on run-ranges
#config.Data.runRange = '300742-302029'

# Where the output files will be transmitted to
config.section_("Site")
#config.Site.storageSite = 'T2_CH_CERN'
config.Site.storageSite = 'T3_US_NotreDame'
#config.Site.whitelist=['T2_US_MIT']
config.section_("User")
config.section_("Debug")
