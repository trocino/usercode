import FWCore.ParameterSet.Config as cms

process = cms.Process("RATE")

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
#    'file:/tmp/trocino/loader/CMSSW_5_2_6_hltpatch1/src/outputA.root',
     'file:/tmp/trocino/noloader/CMSSW_5_2_6_hltpatch1/src/outputA.root',
#     'file:/tmp/trocino/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/195/937/0000/94B218E4-C4B4-E111-B12B-90E6BA19A22E.root',
#     'file:/tmp/trocino/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/195/937/0000/CC324419-C5B4-E111-A61B-90E6BA0D09AD.root',
#     'file:/tmp/trocino/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/195/937/0000/782633E5-C4B4-E111-A6F6-20CF3056171C.root',
#     'file:/tmp/trocino/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/195/950/0000/048A0693-6CB5-E111-9F9F-90E6BAE8CC0C.root'
    ),
                            secondaryFileNames = cms.untracked.vstring()
                            )

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = "GR_R_52_V7::All"

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('efficiency_Mu40_eta2p1_noloader.root'),
                                   closeFileFast = cms.untracked.bool(False)
                                   )

process.MuonTriggerEfficiencyAnalyzer =cms.EDAnalyzer("MuonTriggerEfficiencyAnalyzer",
                                                      vertexes = cms.InputTag("offlinePrimaryVertices"),
                                                      muons = cms.InputTag("muons"),
                                                      #triggerProcess = cms.string("HLT"),
                                                      triggerProcess = cms.string("TEST"),
                                                      tagTriggerName = cms.string("HLT_IsoMu24_eta2p1_v"),
                                                      #triggerName = cms.string("HLT_IsoMu24_eta2p1_v"),
                                                      triggerName = cms.string("HLT_Mu40_eta2p1_v"),
                                                      maxNumberMuons = cms.untracked.uint32(999)
                                                      )
process.HLTMuEffPath = cms.Path(process.MuonTriggerEfficiencyAnalyzer)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
