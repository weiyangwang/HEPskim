import FWCore.ParameterSet.Config as cms
process = cms.Process("anaTree")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
        #'/store/data/Run2015B/ZeroBias1/RECO/PromptReco-v1/000/251/025/00000/FC945C7F-3F26-E511-BC95-02163E011E0A.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/02CB06F2-3A2C-E511-BF63-02163E01399F.root'#,
"""       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/044A9089-392C-E511-9632-02163E0126D2.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/064C105E-3B2C-E511-8084-02163E013770.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/08A2BAD3-3B2C-E511-A477-02163E01241A.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/0E3B9426-392C-E511-B747-02163E01445F.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/1225F1F3-3E2C-E511-A49D-02163E011BB9.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/1849628E-3B2C-E511-8A5A-02163E014761.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/18951EC1-3C2C-E511-9710-02163E011DA4.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/1C94AD29-382C-E511-A0E7-02163E011A5A.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/260746F0-382C-E511-BD9C-02163E011C0F.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/26CEF7DF-392C-E511-A825-02163E012BD2.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/2A554F68-3D2C-E511-A502-02163E0139D5.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/2A7E382A-372C-E511-B03D-02163E013896.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/343A5406-382C-E511-ADBC-02163E012601.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/34ABEFD9-392C-E511-8FC8-02163E011CDB.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/3659AC4B-382C-E511-B364-02163E01208E.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/3893B68A-3F2C-E511-95E9-02163E0133DE.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/3A2E1480-3B2C-E511-B93E-02163E0144F6.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/3AE8EA9E-3B2C-E511-84A4-02163E012B30.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/4239BA97-382C-E511-9F2B-02163E0141E9.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/42F50C8E-382C-E511-9A59-02163E01366D.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/4A60C1BF-392C-E511-A26E-02163E01366D.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/4C8F154E-3A2C-E511-9C78-02163E014437.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/5A2EB254-442C-E511-8990-02163E011DFF.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/5A3C290A-372C-E511-814A-02163E0129A3.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/5AEAF61C-3A2C-E511-8CC3-02163E014543.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/602DAF41-3C2C-E511-80BD-02163E0133DE.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/607EE5F2-3A2C-E511-9360-02163E01364E.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/62CCDC31-412C-E511-98BC-02163E011824.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/68AA6620-3A2C-E511-B966-02163E014186.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/720DA60C-3D2C-E511-AA4A-02163E01258B.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/74E6AE8C-382C-E511-8F40-02163E01439E.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/7CE8FFD0-392C-E511-ACD4-02163E0133CA.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/80753004-362C-E511-96B6-02163E01439E.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/8CA01F9D-3A2C-E511-AC1F-02163E0139A2.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/8CBC80F7-3A2C-E511-9D16-02163E0144DA.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/8EB8B2F6-422C-E511-BBA4-02163E011C0F.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/9285A7B9-462C-E511-8419-02163E011A5A.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/928BCAC5-392C-E511-B943-02163E011824.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/9A667922-3C2C-E511-BA04-02163E013481.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/9E4CF297-362C-E511-8D12-02163E014437.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/9E8C0A21-3C2C-E511-8EBE-02163E0140DC.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/9EA838D9-392C-E511-A826-02163E0136E6.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/A8BBA289-392C-E511-B468-02163E0126D2.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/AC27542A-402C-E511-841C-02163E0133D1.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/B25CCE04-382C-E511-8597-02163E011A29.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/B6FC956C-372C-E511-AD67-02163E0133F9.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/C07DBCEB-382C-E511-B551-02163E013432.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/C0A101B6-3C2C-E511-8BA9-02163E013896.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/CC609034-372C-E511-B2B0-02163E012283.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/CCDC7842-382C-E511-A110-02163E014761.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/D075DD33-3B2C-E511-808A-02163E014729.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/D25F8C69-3D2C-E511-9B98-02163E012377.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/DE928C82-352C-E511-B7BF-02163E014272.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/E6591081-392C-E511-A25A-02163E012BD2.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/EC3FB2DF-382C-E511-BA51-02163E01359E.root',
       '/store/data/Run2015B/ZeroBias8/RECO/PromptReco-v1/000/251/721/00000/FAC95B15-3A2C-E511-9525-02163E01477B.root',"""
        )
                            )
process.load("Configuration.Geometry.GeometryRecoDB_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '74X_dataRun2_Prompt_v4')

process.tree = cms.EDAnalyzer('CTSkim')
process.TFileService = cms.Service("TFileService", fileName = cms.string('data_outputTree.root'))
process.p = cms.Path(process.tree)
