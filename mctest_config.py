import FWCore.ParameterSet.Config as cms
process = cms.Process("anaTree")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(10000) )
#process.SkipEvent = cms.untracked.vstring('ProductNotFound')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
readFiles = cms.untracked.vstring([])
process.source = cms.Source ("PoolSource",fileNames = readFiles)
readFiles.extend( [
'/store/mc/RunIISpring15DR74/MinBias_TuneCUETP8M1_13TeV-pythia8/AODSIM/NoPU_MCRUN2_74_V8-v3/50000/00852C0E-5608-E511-8AB8-BC305B390A32.root',
'/store/mc/RunIISpring15DR74/MinBias_TuneCUETP8M1_13TeV-pythia8/AODSIM/NoPU_MCRUN2_74_V8-v3/50000/02B0E069-5508-E511-A455-0025904C66E4.root',
'/store/mc/RunIISpring15DR74/MinBias_TuneCUETP8M1_13TeV-pythia8/AODSIM/NoPU_MCRUN2_74_V8-v3/50000/049FC3E4-5608-E511-A1DE-20CF3019DF19.root'
] );

process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '76X_mcRun2_asymptotic_v12', '') 
#process.GlobalTag = GlobalTag(process.GlobalTag, '74X_mcRun2_asymptotic_v1')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc')

process.mytree = cms.EDAnalyzer('HEPmcskim')
process.TFileService = cms.Service("TFileService", fileName = cms.string('mc_skimTree.root'))
process.p = cms.Path(process.mytree)
