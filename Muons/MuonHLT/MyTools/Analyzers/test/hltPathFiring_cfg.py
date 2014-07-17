import FWCore.ParameterSet.Config as cms

process = cms.Process("RATE")

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring("file:/afs/cern.ch/user/r/rewang/public/toDaniele/duplicatedEvents/Run_DoubleMu_199336_826_714405303_1_1_5NU.root"), 
                            #fileNames = cms.untracked.vstring("file:/afs/cern.ch/user/r/rewang/public/toDaniele/duplicatedEvents/Run_SingleMu_199336_826_714405303_1_1_40y.root"), 
                            secondaryFileNames = cms.untracked.vstring()
                            )

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = "GR_R_52_V7::All"
#process.GlobalTag.globaltag = "GR_P_V32::All"
process.GlobalTag.globaltag = "GR_P_V41_AN2::All"

process.HLTPathFiring =cms.EDAnalyzer("HLTPathFiring",
                                      hltProcessName = cms.untracked.string("HLT"),         # from collection: edmTriggerResults_TriggerResults__HLT
                                      lumiProcessName = cms.untracked.string("RECO"),       # from collection: LumiScalerss_scalersRawToDigi__RECO
                                      lumiLabel = cms.untracked.string("scalersRawToDigi"), # from collection: LumiScalerss_scalersRawToDigi__RECO
                                      verbose = cms.untracked.bool(False), 
                                      printOutput = cms.untracked.bool(True), 
                                      printAtEndRun = cms.untracked.bool(True), 
                                      ListOfPaths = cms.untracked.vstring("HLT_IsoMu24_eta2p1_v", "HLT_Mu17_Mu8_v"), # List of paths for which you want to compute rate
                                   )
process.HLTPathFiringPath = cms.Path(process.HLTPathFiring)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
