import FWCore.ParameterSet.Config as cms

process = cms.Process("HLTMENU")

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    'file:/tmp/trocino/store/mc/Summer12_DR53X/ZZJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola/AODSIM/PU_S10_START53_V7A-v3/0000/1CD7D49D-58E0-E111-A6D3-E41F1318099C.root'
    ),
                            secondaryFileNames = cms.untracked.vstring()
                            )

#process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.FrameworkJobReport.default.limit = cms.untracked.int32(0)
#process.MessageLogger.FrameworkJobReport.FwkJob.limit = cms.untracked.int32(0)

process.MessageLogger = cms.Service("MessageLogger",
                                    debugModules = cms.untracked.vstring('*'),
                                    suppressDebug = cms.untracked.vstring('FwkJob', 'FwkReport'),
                                    suppressInfo = cms.untracked.vstring('FwkJob', 'FwkReport'),
                                    suppressWarning = cms.untracked.vstring('FwkJob', 'FwkReport'),
                                    cout = cms.untracked.PSet(threshold = cms.untracked.string('ERROR')),
                                    destinations = cms.untracked.vstring('cout')
                                    )

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = "GR_R_52_V7::All"
process.GlobalTag.globaltag = "GR_P_V32::All"

process.hltMenuVersionProvider =cms.EDAnalyzer("HLTMenuVersionProvider",
                                               hltProcessName = cms.untracked.string("HLT") 
                                               )
process.HLTMenuVersionProvider = cms.Path(process.hltMenuVersionProvider)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
