import FWCore.ParameterSet.Config as cms

process = cms.Process("HtoZZto2l2nu")

from CMGTools.HtoZZ2l2nu.localPatTuples_cff import *
process.source = cms.Source("PoolSource",
#                            fileNames = getLocalSourceFor('GluGluToHToZZTo2L2NuM400'),
                            inputCommands = cms.untracked.vstring('keep *',
                                                                  'drop *_MEtoEDMConverter_*_*'),
                            fileNames = cms.untracked.vstring('rfio:/castor/cern.ch/cms/store/cmst3/user/cerminar/ZZllvv_sel3/DYJetsToLL_PU2010/patTuple_99_1_mox.root')
                            )

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )

process.load('CMGTools.HtoZZ2l2nu.NormalizationCounter_cfi')
process.load('CMGTools.HtoZZ2l2nu.CleanEventProducer_cfi')
process.load('CMGTools.HtoZZ2l2nu.CleanEventFilter_cfi')
process.load('CMGTools.HtoZZ2l2nu.CleanEventAnalyzer_cfi')


from CMGTools.HtoZZ2l2nu.StandardSelections_cfi import *
from SimGeneral.MixingModule.mix_E7TeV_FlatDist10_2011EarlyData_inTimeOnly_cfi import mix
process.zzllvvAnalyzer = cms.EDAnalyzer("ZZllvvAnalyzer",
                                        GeneratedPU = mix.input.nbPileupEvents.probValue,
                                        PuDistrFile = cms.untracked.string('root://castorcms//castor/cern.ch/user/t/trocino/ZZllnn/PU-file/pudist_160404-163869.root'),
                                        fileName = cms.untracked.string('ZZllvvAnalyzer.root'),
                                        source = cms.untracked.InputTag("cleanEvent"),
                                        isMC = cms.untracked.bool(False),
                                        xSection = cms.untracked.double(1.0),
                                        branchingRatio = cms.untracked.double(1.0),
                                        luminosity = cms.untracked.double(1.0),
                                        debug = cms.untracked.bool(False),
                                        #Generator = BaseGeneratorSelection.clone(),
                                        Vertices = BaseVertexSelection.clone(),
                                        #Muons = BaseMuonsSelection.clone(),
                                        #Electrons = BaseElectronsSelection.clone(),
                                        #Dileptons = BaseDileptonSelection.clone(),
                                        #Jets = BaseJetSelection.clone(),
                                        #MET = BaseMetSelection.clone(),
                                        FlavorCombination = cms.untracked.int32(0),    ## 0=ee/mm (default), 1=mm, 2=ee, 3=em; 4=any combination
                                        ChargeCombination = cms.untracked.int32(-1),   ## 0=any, -1=opposite (default), 1=same
                                        RecoilLongWeight = cms.untracked.double(2.0),
                                        RecoilPerpWeight = cms.untracked.double(2.0),
                                        SigmaPtLongWeight = cms.untracked.double(2.8),
                                        SigmaPtPerpWeight = cms.untracked.double(2.8),
                                        PerpComponentWeight = cms.untracked.double(1.0),
                                        RedMETMinCut = cms.untracked.double(35.0),
                                        )

# flavorLabel=""
# if process.zzllvvAnalyzer.FlavorCombination == 0:
#     flavorLabel="eemm"
# elif process.zzllvvAnalyzer.FlavorCombination == 1:
#     flavorLabel="mm"
# elif process.zzllvvAnalyzer.FlavorCombination == 2:
#     flavorLabel="ee"
# elif process.zzllvvAnalyzer.FlavorCombination == 3:
#     flavorLabel="em"
# else:
#     flavorLabel="anyFlav"

# chargeLabel=""
# if process.zzllvvAnalyzer.ChargeCombination == -1:
#     chargeLabel="opp"
# elif process.zzllvvAnalyzer.ChargeCombination == 0:
#     chargeLabel="any"
# elif process.zzllvvAnalyzer.ChargeCombination == 1:
#     chargeLabel="same"
# else:
#     chargeLabel="wrong"

# process.zzllvvAnalyzer.fileName = cms.untracked.string('ZZllvvAnalyzer_'+flavorLabel+'_'+chargeLabel+'Charge.root')

process.llnnFilter = cms.EDFilter("ZZllnnFilter")


process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('evHyp.root'),
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
process.p = cms.Path(process.baseSeq)

#process.e = cms.EndPath(process.saveNormalizationCounters*process.out)

# message logger
process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 500
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
    )

process.TFileService = cms.Service("TFileService", fileName = cms.string("histo_iso.root"))
