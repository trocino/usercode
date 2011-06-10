import FWCore.ParameterSet.Config as cms

process = cms.Process("HtoZZto2l2nu")

from CMGTools.HtoZZ2l2nu.localPatTuples_cff import *
process.source = cms.Source("PoolSource",
#                            fileNames = getLocalSourceFor('GluGluToHToZZTo2L2NuM400'),
                            inputCommands = cms.untracked.vstring('keep *',
                                                                  'drop *_MEtoEDMConverter_*_*'),
                            fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_10_1_ryb.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_11_1_XfK.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_12_1_8G1.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_13_1_j7M.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_14_1_YAn.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_15_1_mlw.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_16_1_sz2.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_17_1_Y05.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_18_1_r5m.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_19_1_QQP.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_1_1_Thy.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_20_1_IAN.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_21_1_6lb.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_22_1_sql.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_23_1_PVJ.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_24_1_4JJ.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_25_1_DtW.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_26_1_2wb.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_27_1_gV8.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_28_1_7hZ.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_29_1_rEL.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_2_1_r4O.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_30_1_Ygy.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_31_1_JDf.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_32_1_tP3.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_33_1_JMb.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_34_1_ZWe.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_35_1_K3X.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_36_1_5Ly.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_37_1_sas.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_38_1_tmk.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_39_1_bOn.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_3_1_Mr1.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_40_1_b1a.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_41_1_w0b.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_42_1_txg.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_43_1_Uod.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_44_1_8P0.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_4_1_kSW.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_5_1_LYI.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_6_1_yVo.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_7_1_udi.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_8_1_xlC.root', 
                                                              'rfio:/castor/cern.ch/cms/store/cmst3/user/psilva/Dileptons/ZZtoAnything/ZZtoAnything_Spring11/patTuple_9_1_LA6.root')
                            )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.load('CMGTools.HtoZZ2l2nu.NormalizationCounter_cfi')
process.load('CMGTools.HtoZZ2l2nu.CleanEventProducer_cfi')
process.load('CMGTools.HtoZZ2l2nu.CleanEventFilter_cfi')
process.load('CMGTools.HtoZZ2l2nu.CleanEventAnalyzer_cfi')


from CMGTools.HtoZZ2l2nu.StandardSelections_cfi import *
from SimGeneral.MixingModule.mix_E7TeV_FlatDist10_2011EarlyData_inTimeOnly_cfi import mix
process.zzllvvAnalyzer = cms.EDAnalyzer("ZZllvvAnalyzer",
                                        fileName = cms.untracked.string('/tmp/trocino/ZZllvvAnalyzer_ZZllnn_TEMPLOUTNAME.root'),
                                        GeneratedPU = mix.input.nbPileupEvents.probValue,
                                        PuDistrFile = cms.untracked.string('root://castorcms//castor/cern.ch/user/t/trocino/ZZllnn/PU-file/pudist_160404-163869.root'),
                                        Vertices = BaseVertexSelection.clone(),
                                        source = cms.untracked.InputTag("cleanEvent"),
                                        zmmInput = cms.untracked.InputTag("zMMCand"),
                                        luminosity = cms.untracked.double(191.091),
                                        xSection = cms.untracked.double(7.67),
                                        branchingRatio = cms.untracked.double(1.0),
                                        isMC = cms.untracked.bool(True),
                                        debug = cms.untracked.bool(False),
                                        RecoilLongWeight = cms.untracked.double(TEMPLRECOILLONG),
                                        RecoilPerpWeight = cms.untracked.double(TEMPLRECOILPERP),
                                        SigmaPtLongWeight = cms.untracked.double(TEMPLSIGMAPTLONG),
                                        SigmaPtPerpWeight = cms.untracked.double(TEMPLSIGMAPTPERP),
                                        PerpComponentWeight = cms.untracked.double(TEMPLPERPCOMP),
                                        RedMETMinCut = cms.untracked.double(TEMPLREDMETCUT),
                                        )


process.llnnFilter = cms.EDFilter("ZZllnnFilter")


process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('/tmp/trocino/evHyp.root'),
                               outputCommands = cms.untracked.vstring('drop *',
                                                                      'keep *_MEtoEDMConverter_*_*',
                                                                      'keep edmMergeableCounter_*_*_*',#FIXME: useful?
                                                                      'keep *_prunedGen_*_*',
                                                                      'keep *_genEventScale_*_*',
                                                                      'keep GenRunInfoProduct_*_*_*',
                                                                      'keep *_genMetTrue_*_*',
                                                                      'keep *_selectedPat*_*_*',
                                                                      'keep patMETs_*_*_*',
                                                                      'keep double*_*_rho_*',
                                                                      #'keep *_tcMet_*_*',
                                                                      'keep *_pfMet_*_*',
                                                                      'keep *_cleanEvent_*_*'),
                               SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') )
                               )


process.baseSeq = cms.Sequence(process.loadNormalizationCounters*process.cleanEvent*process.cleanEventFilter*process.zzllvvAnalyzer)
process.llnnFilt = cms.Sequence(process.llnnFilter)
process.llnnAntifilt = cms.Sequence(~process.llnnFilter)
process.p = cms.Path(process.llnnFilt+process.baseSeq)

#process.e = cms.EndPath(process.saveNormalizationCounters*process.out)

# message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 500
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
    )

process.TFileService = cms.Service("TFileService", fileName = cms.string("/tmp/trocino/histo_iso.root"))
