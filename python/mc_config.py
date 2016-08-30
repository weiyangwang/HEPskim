import FWCore.ParameterSet.Config as cms
process = cms.Process("anaTree")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('file:/afs/cern.ch/work/s/sakumar/public/2016/03/0E4C6224-773E-E511-A268-002590A370DC.root')
                            )

process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '76X_mcRun2_asymptotic_v12', '') 
#process.GlobalTag = GlobalTag(process.GlobalTag, '74X_mcRun2_asymptotic_v1')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')

process.mytree = cms.EDAnalyzer('CTSkim')
process.TFileService = cms.Service("TFileService", fileName = cms.string('mc_skimTree.root'))
process.p = cms.Path(process.mytree)
