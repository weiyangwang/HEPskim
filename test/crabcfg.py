from CRABClient.UserUtilities import config
config = config()

config.General.workArea = '2016/minBias'
#config.General.workArea = '2016/z1'
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../python/mc_config.py'
#config.JobType.psetName = 'python/data_cpnfig.py'

config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/sakumar/2016/05/17/minBias'
#config.Data.outLFNDirBase = '/store/user/sakumar/2016/04/15/z1'
config.Data.publication = False
config.Site.storageSite = 'T2_IN_TIFR'
#config.Data.lumiMask = '/afs/cern.ch/work/s/sakumar/private/2016/01/1/CMSSW_7_6_3/src/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON_v2.txt'
#config.Site.whitelist = ['T2_KR_KNU']
#config.Site.whitelist = ['T2_DE_DESY']
#config.Site.whitelist = ['T1_IT_CNAF', 'T1_US_FNAL']
#config.section_("General")
#config.section_("JobType")
#config.section_("MC")
#config.section_("Site")

#config.Data.inputDBS = 'global'
#config.Data.splitting = 'FileBased'
config.Data.inputDataset = '/MinBias_TuneZ2star_13TeV-pythia6/RunIIFall15DR76-25nsNoPUFSQ_castor_76X_mcRun2_asymptotic_v12-v1/GEN-SIM-RECODEBUG'
#config.Data.inputDataset = '/ZeroBias2/Run2015D-PromptReco-v4/RECO'
#config.Data.inputDataset = '/MinBias_TuneEE5C_13TeV-herwigpp/RunIIWinter15GS-castor_MCRUN2_71_V0-v1/GEN-SIM'
#config.Data.inputDataset = '/ZeroBias1/Run2015B-PromptReco-v1/RECO'
