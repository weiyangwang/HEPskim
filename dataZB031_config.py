import FWCore.ParameterSet.Config as cms
process = cms.Process("anaTree")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000000) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#data set with HM triggers '/HighMultiplicity/Run2015B-PromptReco-v1/RECO'
# process.source = cms.Source("PoolSource",
#                            fileNames = cms.untracked.vstring(
#	'/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/02CB06F2-3A2C-E511-BF63-02163E01399F.root',
#	'/store/data/Run2015B/HighMultiplicity/RECO/PromptReco-v1/000/251/096/00000/1C89895E-5626-E511-A518-02163E0119E7.root'
#	'/store/data/Run2015B/HighMultiplicity/RECO/PromptReco-v1/000/251/721/00000/00BA646E-D82B-E511-A218-02163E01439E.root'
#        )
#                            )
readFiles = []
readFiles.extend( [
 '/store/data/Run2015C_25ns/ZeroBias1/AOD/16Dec2015-v1/20000/00230D52-A4B5-E511-BC87-0025904C6508.root'
] ) 
process.source = cms.Source ("PoolSource",fileNames = cms.untracked.vstring(readFiles))
process.load("Configuration.Geometry.GeometryRecoDB_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '74X_dataRun2_Prompt_v4')

process.tree = cms.EDAnalyzer('HEPskim')
process.TFileService = cms.Service("TFileService", fileName = cms.string('data_outputTree.root'))
process.p = cms.Path(process.tree)
