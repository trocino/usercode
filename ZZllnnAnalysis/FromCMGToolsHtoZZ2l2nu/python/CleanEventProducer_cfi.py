import FWCore.ParameterSet.Config as cms

from CMGTools.HtoZZ2l2nu.StandardSelections_cfi import *
cleanEvent = cms.EDProducer("DileptonPlusMETEventProducer",
                            Generator = BaseGeneratorSelection.clone(),
                            Vertices = BaseVertexSelection.clone(),
                            Muons = BaseMuonsSelection.clone(),
                            Electrons = BaseElectronsSelection.clone(),
                            Dileptons = BaseDileptonSelection.clone(),
                            Tracks = BaseTrackSelection.clone(),
                            Jets = BaseJetSelection.clone(),
                            MET = BaseMetSelection.clone()
                            )
