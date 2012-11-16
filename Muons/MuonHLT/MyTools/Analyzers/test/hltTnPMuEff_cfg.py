import FWCore.ParameterSet.Config as cms

process = cms.Process("RATE")

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/0C1A907B-ACBB-E111-8A7F-E0CB4E553642.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/0C7FC3E2-ACBB-E111-8EE4-20CF3056171D.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/140C0D69-ACBB-E111-8B9A-485B39897214.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/149DC466-ACBB-E111-B275-90E6BA442F31.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/18CF47E4-ACBB-E111-BB54-E0CB4EA0A906.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/1C151988-ACBB-E111-8347-20CF3027A596.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/1ED20D6A-ADBB-E111-B4A0-20CF3027A57B.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/22FB0ADC-ACBB-E111-A7C5-90E6BA0D09D4.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/2A38ADB3-ACBB-E111-B46F-BCAEC50971E3.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/2E5C5E71-ACBB-E111-8DAB-20CF3019DF0F.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/329AC3BA-ACBB-E111-9B8E-485B39800C14.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/38EAD6EB-ADBB-E111-A515-485B39897264.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/4C0D0F73-ADBB-E111-85A6-90E6BA442F1F.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/4EE8D967-ACBB-E111-9775-BCAEC54B304A.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/507F1B19-ADBB-E111-AFA0-20CF3027A639.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/509EB7AE-ACBB-E111-8A84-E0CB4E1A1192.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/52B67E3B-ADBB-E111-8BAC-20CF305B050D.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/52E64DF1-ACBB-E111-B8C0-E0CB4E29C4E8.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/52EADB22-ADBB-E111-B7EC-20CF305B04DA.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/56446158-ADBB-E111-966E-90E6BA19A231.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/5A30B98E-ACBB-E111-8EA3-20CF30724B0A.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/5A604536-ADBB-E111-8FA5-E0CB4E1A1177.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/5C2820BA-ACBB-E111-A5F2-485B39800BB3.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/5CACABAA-ACBB-E111-905F-E0CB4E29C51D.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/5E02324E-ADBB-E111-B0F1-20CF30561712.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/5E2A75FC-ACBB-E111-AA48-20CF3027A5E5.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/6476BF6A-ACBB-E111-B7CB-20CF305B052B.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/6A296FA4-AFBB-E111-83DE-0030487CDAC2.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/72EF9AA3-ACBB-E111-9EF0-20CF305B0508.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/78C5C1B9-ACBB-E111-91BA-20CF3027A635.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/7C420A91-ACBB-E111-89F5-20CF3027A5E5.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/7CC1C757-AEBB-E111-9056-0030487CDB2C.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/7CDC9957-ADBB-E111-9AF3-90E6BA442EEB.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/8486E2C1-ACBB-E111-85A8-20CF3027A578.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/8C4000C9-ACBB-E111-8D57-20CF3056171D.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/9075C8F4-ACBB-E111-B95B-20CF300E9ED3.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/90ECD1B0-ACBB-E111-830D-00261834B520.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/96ED355E-ADBB-E111-B2A5-90E6BAE8CC18.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/9A6CD8A7-ACBB-E111-8BAE-90E6BA19A212.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/A6B184BD-ADBB-E111-9726-485B39800C23.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/A829BC02-ADBB-E111-A8D9-E0CB4E19F972.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/A8A661EA-ACBB-E111-AC66-20CF3027A5C0.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/AC0DD790-ADBB-E111-B1EB-20CF3058709A.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/B468ED03-ADBB-E111-903D-BCAEC50971E3.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/B61DE437-ADBB-E111-AAC0-20CF305B050D.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/BC50DC1F-ADBB-E111-86B3-001EC9D8A8D0.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/C0561D9D-ADBB-E111-958D-00261834B5B1.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/C4E215D9-ACBB-E111-A133-00261834B5A0.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/C81A1E16-AEBB-E111-A6F3-E0CB4E29C517.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/CEDC8C30-ADBB-E111-9EDB-20CF3027A639.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/D0CF0D6A-ADBB-E111-BE67-20CF3027A57B.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/D2832B83-ACBB-E111-90FF-20CF305B050D.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/D2AF5B4E-ADBB-E111-B258-90E6BAE8CC08.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/D6BE12A2-ADBB-E111-A010-E0CB4E19F9BB.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/D884AB5A-ADBB-E111-83A0-E0CB4E1A1144.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/DE2FC9BE-ACBB-E111-B14A-90E6BA442F16.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/DE46BEE9-ACBB-E111-9DDA-20CF305B0585.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/DE51BBF5-ACBB-E111-967B-E0CB4E1A1172.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/E2199A5E-ADBB-E111-9B15-E0CB4E1A119E.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/E2EE4278-ACBB-E111-9713-E0CB4E553637.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/EAC6598F-ADBB-E111-B762-00261834B548.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/EE7D38BF-ACBB-E111-9891-485B39897259.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/F86B8999-ACBB-E111-BEC4-20CF300E9ECF.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/FA1A1C3C-ADBB-E111-8437-20CF3027A5AF.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/FCD3BF91-ADBB-E111-8A6C-00261834B5B0.root',
    '/store/data/Run2012B/SingleMu/RAW-RECO/ZMu-PromptSkim-v2/000/196/531/0000/FEE4268C-ACBB-E111-A998-485B39800BBE.root'
    ),
                            secondaryFileNames = cms.untracked.vstring()
                            )

process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = "GR_R_52_V7::All"

PATHNAME = "Mu40"
PROBEDEN = "hltL2fL1sMu16L1f0L2Filtered16Q"
PROBENUM = "hltL3fL1sMu16L1f0L2f16QL3Filtered40Q"

process.MuonTriggerEfficiencyAnalyzer =cms.EDAnalyzer("MuonTriggerEfficiencyAnalyzer",
                                                      vertexes = cms.InputTag("offlinePrimaryVertices"),
                                                      muons = cms.InputTag("muons"),
                                                      triggerProcess = cms.string("HLT"), #"TEST"
                                                      tagTriggerName = cms.string("HLT_IsoMu24_v"),
                                                      triggerName = cms.string("HLT_"+PATHNAME+"_v"), # hltL1sMu16, hltL2fL1sMu16L1f0L2Filtered16Q, hltL3fL1sMu16L1f0L2f16QL3Filtered40Q 
                                                      probeFilterDen = cms.string(PROBEDEN),
                                                      probeFilterNum = cms.string(PROBENUM),
                                                      maxNumberMuons = cms.untracked.uint32(999999)
                                                      )

OUTFILE = 'efficiency_'+PATHNAME

if PROBENUM != '': OUTFILE = OUTFILE+'_'+PROBENUM
else:              OUTFILE = OUTFILE+'_fullpath'

if PROBEDEN != '': OUTFILE = OUTFILE+'_over_'+PROBEDEN

OUTFILE = OUTFILE+'.root'

print "Output file: "+OUTFILE

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string(OUTFILE),
                                   closeFileFast = cms.untracked.bool(False)
                                   )

process.HLTMuEffPath = cms.Path(process.MuonTriggerEfficiencyAnalyzer)

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(10000))
