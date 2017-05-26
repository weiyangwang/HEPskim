import FWCore.ParameterSet.Config as cms
process = cms.Process("anaTree")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#data set with HM triggers '/HighMultiplicity/Run2015B-PromptReco-v1/RECO'
#process.source = cms.Source("PoolSource",
#                            fileNames = cms.untracked.vstring(
#	'/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/02CB06F2-3A2C-E511-BF63-02163E01399F.root',
#	'/store/data/Run2015B/HighMultiplicity/RECO/PromptReco-v1/000/251/096/00000/1C89895E-5626-E511-A518-02163E0119E7.root'
#	'/store/data/Run2015B/HighMultiplicity/RECO/PromptReco-v1/000/251/721/00000/00BA646E-D82B-E511-A218-02163E01439E.root'
#        )
#                            )
readFiles = []
readFiles.extend( [
#'/store/data/Run2015C_25ns/L1MinimumBiasHF2/AOD/16Dec2015-v1/10000/005CA025-82B0-E511-B1EA-0CC47A4D76AA.root',
#'/store/data/Run2015C_25ns/L1MinimumBiasHF3/AOD/16Dec2015-v1/10000/F6E4AA0D-5BB6-E511-9019-00266CFF0A88.root',
#'/store/data/Run2015C_25ns/L1MinimumBiasHF4/AOD/16Dec2015-v1/20000/F47DD467-DCB5-E511-8013-00266CF9AE10.root',
#'/store/data/Run2015C_25ns/L1MinimumBiasHF5/AOD/16Dec2015-v1/10000/FC116F30-B5B0-E511-98F8-0025905A60D6.root',
'/store/data/Run2015C_25ns/L1MinimumBiasHF6/AOD/16Dec2015-v1/20000/F4570B14-31B3-E511-84E6-C4346BC76CD8.root',
'/store/data/Run2015C_25ns/L1MinimumBiasHF7/AOD/16Dec2015-v1/10000/FA8FE0B8-48B6-E511-84AF-001E67397B11.root',
'/store/data/Run2015C_25ns/L1MinimumBiasHF8/AOD/16Dec2015-v1/10000/FC7D3EB8-8FB3-E511-94D1-002590593872.root'
] ) 
process.source = cms.Source ("PoolSource",fileNames = cms.untracked.vstring(readFiles))
process.load("Configuration.Geometry.GeometryRecoDB_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '74X_dataRun2_Prompt_v4')

process.tree = cms.EDAnalyzer('HEPskim')
process.TFileService = cms.Service("TFileService", fileName = cms.string('/afs/cern.ch/work/q/qileong/public/MB031test.root'))
process.p = cms.Path(process.tree)