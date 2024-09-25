# CRAB3 config template for flashgg
# More options available on the twiki :
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCRAB3Tutorial

from WMCore.Configuration import Configuration

config = Configuration()

config.section_('General')
config.General.requestName       = 'DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8'
config.General.transferLogs      = True
config.General.transferOutputs   = True

config.section_('JobType')
config.JobType.pluginName        = 'Analysis'

# Name of the CMSSW configuration file
config.JobType.psetName          = 'EleDumper_cfg.py'
config.JobType.inputFiles        = ["EleDumper_cfg.py","EleDumper_cfi.py"]
config.JobType.outputFiles       = ["output.root"]
#config.JobType.disableAutomaticOutputCollection = True
## Memory, cores, cmssw
config.JobType.allowUndistributedCMSSW = True
config.JobType.priority          = 30
config.JobType.maxMemoryMB       = 2500

config.section_('Data')
# This string determines the primary dataset of the newly-produced outputs.
config.Data.inputDataset         = '/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1/MINIAODSIM'
#config.Data.inputDataset         = '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1/MINIAODSIM'
#config.Data.inputDataset         = '/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer19UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v2/MINIAODSIM'
config.Data.inputDBS             = 'global'     
config.Data.splitting            = 'FileBased'
config.Data.unitsPerJob          = 1
config.Data.publication          = False
config.Data.ignoreLocality       = True
config.Data.outputDatasetTag     = 'RunIISummer20UL18MiniAOD-106X_upgrade2018_realistic_v11_L1v1-v1'

# This string is used to construct the output dataset name
#config.Data.publishDataName = 'CRAB3-tutorial'
config.Data.outLFNDirBase        =  '/store/group/dpg_ecal/alca_ecalcalib/bmarzocc/Regression_Studies/'

config.section_('Site')
# Where the output files will be transmitted to
config.Site.storageSite         = 'T2_CH_CERN'
config.Site.whitelist           = ['T2_CH_CERN']


## config.Data.allowNonValidInputDataset=True
