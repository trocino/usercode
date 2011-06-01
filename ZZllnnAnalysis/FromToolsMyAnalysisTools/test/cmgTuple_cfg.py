from PhysicsTools.PatAlgos.patTemplate_cfg import *
import FWCore.ParameterSet.Config as cms

#---------------------------------------------

#---------------------------------------------


process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
)
process.source.fileNames = cms.untracked.vstring(
        '/store/cmst3/user/wreece/VJets/3_8_7/020211/Zmumu/PF2PAT_MuonFilter_99_1_8aQ.root'
)
process.setName_('CMG')
process.out.fileName = cms.untracked.string('file:/tmp/cerminar/cmgTuple.root')

isData = False

### GENERATOR LEVEL
if( not isData ):
    process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
    process.genHiggsFinalState = cms.EDProducer(
        "GenParticlePruner",
        src = cms.InputTag("genParticles"),
        select = cms.vstring(
        "drop  *",
        "keep++ pdgId={h0}",
        "drop pdgId={h0} & status=2",
        "drop pdgId={Z0} & status=2",
        "drop pdgId={W+} & status=2",
        "drop pdgId={W-} & status=2"
        )
        )
    process.load('CMGTools.Common.factories.cmgBaseGenCandidates_cfi')
    process.cmgHiggsFinalState = process.cmgGenCandidates.clone()
    process.cmgHiggsFinalState.cfg.inputCollection = 'genHiggsFinalState'
    process.genSequence = cms.Sequence( process.genHiggsFinalState+process.cmgHiggsFinalState)



# process.load('CMGTools.Common.cutsummary_cff')
# process.load('CMGTools.Common.muon_cff')
# process.load('CMGTools.Common.diMuon_cff')
# process.load('CMGTools.Common.jet_cff')
# process.load('CMGTools.Common.runInfoAccounting_cfi')

#config the pat trigger matching
process.load('CMGTools.ZmumuJetsTutorial.patTriggerMatching_cff')
#process.patMuonTriggerResults.processName=cms.string(trigMenu)
process.patMuonTriggerMatch.pathNames = cms.vstring('HLT_Mu9','HLT_Mu11','HLT_Mu15_v1')


### DI-MUON selection
process.load('CMGTools.Common.muon_cff')
process.load('CMGTools.Common.diMuon_cff')
process.cmgDiMuon.cuts.zmumu.pt = cms.string('pt() >= 0')
process.cmgDiMuon.cuts.zmumu.mass = cms.string('mass() >= 50')
process.cmgDiMuon.cuts.zmumu.leg1_quality = cms.string('leg1().getSelection(\"cuts_vbtfmuon_isGlobal\")')
process.cmgDiMuon.cuts.zmumu.leg1_kinematics.eta = cms.string('abs(leg1().eta()) < 2.4')
process.cmgMuon.cfg.inputCollection = cms.InputTag("triggeredPatMuons")

### JET selection
process.load('CMGTools.Common.jet_cff')
myjets=cms.InputTag('patJetsCachedReClusterZPUL1Corr')
process.cmgPFJet.cfg.inputCollection = myjets
process.cmgPFJet.cfg.btagType = cms.string('trackCountingHighEffBJetTags')
process.cmgPFBaseJet.cfg.inputCollection = myjets
process.cmgPFBaseJet.cfg.btagType = cms.string('trackCountingHighEffBJetTags')


#### MET
process.load('CMGTools.Common.met_cff')
#mymet=cms.InputTag('patMETsReClusterZPU')
mymet=cms.InputTag('patMETsPFlow')
process.cmgMHTPFJets.cfg.inputCollection = myjets
process.cmgMHTPFJets30.cfg.inputCollection = myjets
process.cmgMETPFCandidates.cfg.inputCollection = mymet
process.metSequence = cms.Sequence( process.cmgMHTPFJets+
                                    process.cmgMHTPFJets30+
                                    process.cmgMETPFCandidates)

### produce a di-muon+MET object
process.load('CMGTools.Common.factories.cmgDiMuonPlusMET_cfi')
process.selectedDiMuonMETSequence = cms.Sequence( process.cmgDiMuonPlusMET*process.selectedDiMuonPlusMETFilter )



from CMGTools.Common.EventContent.particleFlow_cff import particleFlow as particleFlowEventContent  
from CMGTools.Common.EventContent.particleFlow_cff import particleFlowBase as particleFlowEventContentBase
from CMGTools.Common.EventContent.particleFlow_cff import particleFlowMHT as particleFlowMHTContentBase
process.out.outputCommands = cms.untracked.vstring( 'drop *')
process.out.outputCommands.extend( particleFlowEventContent ) 
process.out.outputCommands.extend( particleFlowEventContentBase )
process.out.outputCommands.extend( particleFlowMHTContentBase )
process.out.outputCommands.append( 'keep *_triggeredPatMuons_*_*' ) # we keep the pat::Muon
process.out.outputCommands.append( 'keep *_patMETsPFlow_*_*' ) # MET
#process.out.outputCommands.append( 'keep *_genParticles_*_*' ) # GEN Particles
process.out.outputCommands.append( 'keep *_genHiggsFinalState_*_*' )
process.out.outputCommands.append( 'keep *_cmgHiggsFinalState_*_*' )
process.out.outputCommands.append( 'keep GenRunInfoProduct_generator_*_*' )
process.out.outputCommands.append( 'keep double*_*_MEtoEDMConverterRun_*' )
process.out.outputCommands.append( 'drop patMETs_*_*_*' )


#output file for histograms etc
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("histograms.root"))

### Produce a summary of cuts on the DiMuon
##process.zmumusummary = process.cutsummary.clone(inputCollection = cms.InputTag("cmgDiMuon"))

### Run the trigger matching selections
##from CMGTools.ZmumuJetsTutorial.selections.muontrigger_cfi import muontrigger
##process.cmgMuon.cuts.muontrigger = muontrigger.clone()
###

###The jet counting histograms
##process.load('CMGTools.ZmumuJetsTutorial.histograms.jetcounting_cfi')
##process.load('CMGTools.ZmumuJetsTutorial.skims.selEventsPFJet_cff')
##process.jetCountingClean = process.jetCounting.clone(inputCollection = cms.InputTag("jetIDJet"))
##process.jetCountingSequence = cms.Sequence(process.jetCounting+process.jetIDJet+process.jetCountingClean)
### The Z skimming
#process.load('CMGTools.ZmumuJetsTutorial.skims.selEventsZ_cff')
###

### the sequence
process.analysisSequence = cms.Sequence(
    process.patMuonTrigger +
    process.muonSequence +
    process.diMuonSequence +
    process.pfJetSequence +
    process.metSequence +
    process.selectedDiMuonMETSequence )

process.p = cms.Path(process.analysisSequence)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.suppressWarning=cms.untracked.vstring('patMuonTriggerResults')
