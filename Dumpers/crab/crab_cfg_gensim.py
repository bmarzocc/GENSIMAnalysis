

# CRAB3 config template for flashgg
# More options available on the twiki :
# https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCRAB3Tutorial

from WMCore.Configuration import Configuration

config = Configuration()

config.section_('General')
config.General.requestName       = 'ElectronGun_E_0-1000_eta1_phi1_106X_upgrade2018_realistic_v4_noBeamspot'
config.General.transferLogs      = True
config.General.transferOutputs   = True

config.section_('JobType')
config.JobType.pluginName        = 'PrivateMC'

# Name of the CMSSW configuration file
#config.JobType.psetName          = 'Photons_closeEcal_E1-1000_simAnalysis_cfi.py'
config.JobType.psetName          = 'Electrons_closeEcal_E1-1000_simAnalysis_cfi.py'
config.JobType.disableAutomaticOutputCollection = True
config.JobType.outputFiles       = ["output.root"]
## Memory, cores, cmssw
config.JobType.allowUndistributedCMSSW = True
config.JobType.maxMemoryMB = 5500
config.JobType.numCores    = 4

config.section_('Data')
# This string determines the primary dataset of the newly-produced outputs.

config.Data.splitting            = 'EventBased'
config.Data.unitsPerJob          = 100
config.Data.totalUnits           = 1000000
config.Data.publishDBS           = 'phys03'
config.Data.publication          = False
#config.Data.ignoreLocality       = True
config.Data.outputDatasetTag     = 'GEN-SIM'
config.Data.outputPrimaryDataset = 'ElectronGun_E_0-1000_eta1_phi1_106X_upgrade2018_realistic_v4_noBeamspot'

# This string is used to construct the output dataset name
#config.Data.publishDataName = 'CRAB3-tutorial'
config.Data.outLFNDirBase        =  '/store/group/phys_higgs/cmshgg/bmarzocc/FNUF/GEN-SIM-GUNS'

config.section_('Site')
# Where the output files will be transmitted to
config.Site.storageSite         = 'T2_CH_CERN'
#config.Site.whitelist           = ['T2_CH_CERN']


## config.Data.allowNonValidInputDataset=True
