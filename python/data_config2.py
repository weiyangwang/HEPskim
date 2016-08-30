import FWCore.ParameterSet.Config as cms
process = cms.Process("anaTree")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/02CB06F2-3A2C-E511-BF63-02163E01399F.root'
        )
                            )
process.load("Configuration.Geometry.GeometryRecoDB_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '74X_dataRun2_Prompt_v4')

process.tree = cms.EDAnalyzer('HEPskim')
process.TFileService = cms.Service("TFileService", fileName = cms.string('data_outputTree.root'))
process.p = cms.Path(process.tree)
