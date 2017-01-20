from CRABClient.UserUtilities import config
config = config()

config.General.workArea = 'crabskim255031/ZB'
#config.General.workArea = '2016/z1'
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'dataZB031_config.py'
#config.JobType.psetName = 'python/data_cpnfig.py'

config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 4
# config.Data.outLFNDirBase = '/store/user/sakumar/2016/05/17/minBias'
config.Data.outLFNDirBase = '/store/user/wei/multiplicity/data/255031/ZB'
config.Data.publication = False
config.Site.storageSite = 'T2_CH_CERN'
config.Data.lumiMask = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Reprocessing/Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_LOWPU_JSON_MuonPhys.txt'
#config.Data.runRange = '255031-255031'
#config.Site.whitelist = ['T3_UK_London_UCL']
#config.Site.whitelist = ['T2_KR_KNU']
#config.Site.whitelist = ['T2_DE_DESY']
#config.Site.whitelist = ['T1_IT_CNAF', 'T1_US_FNAL']
#config.section_("General")
#config.section_("JobType")
#config.section_("MC")
#config.section_("Site")

#config.Data.inputDBS = 'global'
#config.Data.splitting = 'FileBased'
# config.Data.inputDataset = '/MinBias_TuneZ2star_13TeV-pythia6/RunIIFall15DR76-25nsNoPUFSQ_castor_76X_mcRun2_asymptotic_v12-v1/GEN-SIM-RECODEBUG'
#config.Data.inputDataset = '/ZeroBias2/Run2015D-PromptReco-v4/RECO'
#config.Data.inputDataset = '/MinBias_TuneEE5C_13TeV-herwigpp/RunIIWinter15GS-castor_MCRUN2_71_V0-v1/GEN-SIM'
config.Data.inputDataset = '/ZeroBias3/Run2015C_25ns-16Dec2015-v1/AOD'
