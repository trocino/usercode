import FWCore.ParameterSet.Config as cms

process = cms.Process("ME0MUON")

process.source = cms.Source("PoolSource", 
                            fileNames = cms.untracked.vstring( TEMPLFILELIST ), 
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck') 
                            ) 
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1)) 


process.load('Configuration.Geometry.GeometryExtended2023D4Reco_cff') 
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff') 
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAlong_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorOpposite_cfi")
process.load("TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi")
process.load("DQMServices/Core/DQMStore_cfg")
process.load("Validation.RecoMuon.associators_cff")
process.load('Validation.RecoMuon.selectors_cff')
process.load('Validation.RecoMuon.MuonTrackValidator_cfi')
process.load('Validation.RecoMuon.RecoMuonValidator_cfi')

# Some of these are useless, probably 
from Validation.RecoMuon.selectors_cff import *
from Validation.RecoMuon.associators_cff import *

# Configurations for MuonTrackValidators
import Validation.RecoMuon.MuonTrackValidator_cfi

# Configurations for RecoMuonValidators
from RecoMuon.TrackingTools.MuonServiceProxy_cff import *
from Validation.RecoMuon.RecoMuonValidator_cfi import *

import SimMuon.MCTruth.MuonAssociatorByHits_cfi
me0MuonAsso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
    UseTracker = True,
    UseMuon = False,
    useGEMs = cms.bool(False),
    EfficiencyCut_track = cms.double(0.0),
    PurityCut_track = cms.double(0.0),
    pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
    stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
    tracksTag = cms.InputTag(""),
    tpTag = cms.InputTag("")
    )

from Configuration.AlCa.GlobalTag import GlobalTag 
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase2_realistic', '') 

## TrackingParticle selector and producer 
aTkPprod = cms.EDProducer("TPMuonTrackProducer", 
                          tpCollectionTag  = cms.InputTag("mix:MergedTrackTruth")) 

############### 
##           ## 
## TEMPLATES ## 
##           ## 
############### 

segmentcoll = [ 
    # ## Strip: 768; Eta partitions: 4, 6, 8, 16 
    # "s768p004", 
    # "s768p006", 
    # "s768p008", 
    # "s768p012", 
    # "s768p016", 
    # #"s768p032", 
    # #"s768p064", 
    # ## Strip: 512; Eta partitions: 4, 6, 8, 16 
    # "s512p004", 
    "s512p006" #, 
    # "s512p008", 
    # "s512p012", 
    # "s512p016", 
    # #"s512p032", 
    # #"s512p064", 
    # ## Strip: 384; Eta partitions: 4, 6, 8, 16 
    # "s384p004", 
    # "s384p006", 
    # "s384p008", 
    # "s384p012", 
    # "s384p016", 
    # #"s384p032", 
    # #"s384p064", 
    # ## Strip: 256; Eta partitions: 4, 6, 8, 16 
    # "s256p004", 
    # "s256p006", 
    # "s256p008", 
    # "s256p012", 
    # "s256p016" #, 
    # #"s256p032", 
    # #"s256p064" 
    ] 

paramcoll = [
    ## 
    ## Default:  [3.0, 4.0, -1.0, -1.0, -1.0]
    ## 
    [-1.0,  4.0, -1.0, -1.0, -1.0], 
    [ 3.0,  4.0, -1.0, -1.0, -1.0], 
    [ 2.0,  4.0, -1.0, -1.0, -1.0], 
    [ 1.0,  4.0, -1.0, -1.0, -1.0], 
    [ 0.5,  4.0, -1.0, -1.0, -1.0], 
    [ 0.1,  4.0, -1.0, -1.0, -1.0], 

    [ 3.0, -1.0, -1.0, -1.0, -1.0], 
    [ 3.0,  3.0, -1.0, -1.0, -1.0], 
    [ 3.0,  2.0, -1.0, -1.0, -1.0], 
    [ 3.0,  1.0, -1.0, -1.0, -1.0], 
    [ 3.0,  0.5, -1.0, -1.0, -1.0], 

    [ 3.0,  4.0, 10.0, -1.0, -1.0], 
    [ 3.0,  4.0,  5.0, -1.0, -1.0], 
    [ 3.0,  4.0,  4.0, -1.0, -1.0], 
    [ 3.0,  4.0,  3.0, -1.0, -1.0], 
    [ 3.0,  4.0,  2.0, -1.0, -1.0], 
    [ 3.0,  4.0,  1.0, -1.0, -1.0], 
    [ 3.0,  4.0,  0.5, -1.0, -1.0], 

    [ 3.0,  4.0, -1.0, 20.0, -1.0], 
    [ 3.0,  4.0, -1.0, 10.0, -1.0], 
    [ 3.0,  4.0, -1.0,  5.0, -1.0], 
    [ 3.0,  4.0, -1.0,  4.0, -1.0], 
    [ 3.0,  4.0, -1.0,  3.0, -1.0], 
    [ 3.0,  4.0, -1.0,  2.0, -1.0], 
    [ 3.0,  4.0, -1.0,  1.0, -1.0], 

    [ 3.0,  4.0, -1.0, -1.0, 1.60], 
    [ 3.0,  4.0, -1.0, -1.0, 0.80], 
    [ 3.0,  4.0, -1.0, -1.0, 0.40], 
    [ 3.0,  4.0, -1.0, -1.0, 0.20], 
    [ 3.0,  4.0, -1.0, -1.0, 0.10] 
    ] 

## Muon selector and producer 
aMe0prod = cms.EDProducer("ME0MuonTrackProducer", 
                          muonsTag         = cms.InputTag("me0SegmentMatching"), 
                          me0segmentsTag   = cms.InputTag("me0Segments"), 
                          trackType        = cms.string("innerTrack"), 
                          maxPullX         = cms.untracked.double(      3.0 ), 
                          maxDiffX         = cms.untracked.double(      4.0 ), 
                          maxPullY         = cms.untracked.double(     -1.0 ), 
                          maxDiffY         = cms.untracked.double(     -1.0 ), 
                          maxDiffPhiDir    = cms.untracked.double(     -1.0 ), 
                          minP             = cms.untracked.double( 999999.0 ), 
                          minPt            = cms.untracked.double(      0.0 ) 
                          ) 

aMe0asso = SimMuon.MCTruth.MuonAssociatorByHits_cfi.muonAssociatorByHits.clone(
    UseTracker = True,
    UseMuon = False,
    useGEMs = cms.bool(False),
    EfficiencyCut_track = cms.double(0.0),
    PurityCut_track = cms.double(0.0),
    pixelSimLinkSrc = cms.InputTag("simSiPixelDigis", "Pixel"),
    stripSimLinkSrc = cms.InputTag("simSiStripDigis", "Tracker"),
    tracksTag = cms.InputTag("aMe0prod"),
    tpTag = cms.InputTag("tpcPro")
    )

aMe0anal = cms.EDAnalyzer("ME0MuonGeneralAnalyzer", 
                          verbose          = cms.untracked.bool(False), 
                          printOutput      = cms.untracked.bool(False), 
                          associatorNames  = cms.vstring("aMe0asso"), 
                          minP             = cms.untracked.double(999999.0), 
                          minPt            = cms.untracked.double(3.0), 
                          puInfoTag        = cms.InputTag("addPileupInfo"), 
                          me0segmentsTag   = cms.InputTag("me0Segments"), 
                          tpCollectionTag  = cms.InputTag("tpcPro"), 
                          muonTag          = cms.InputTag("me0SegmentMatching"), 
                          tracksTag        = cms.InputTag("aMe0prod")
                          ) 


for i in range(len(segmentcoll)): 
    for j in range(len(paramcoll)): 
        paramlab = "PullX" + ((str(paramcoll[j][0]) + "DiffX" + str(paramcoll[j][1]) + "PullY" + str(paramcoll[j][2]) + "DiffY" + str(paramcoll[j][3]) + "DiffPhi" + str(paramcoll[j][4])).replace(".", "p")).replace("-", "m") 

        thismuonprod = "muoPro" + segmentcoll[i] + paramlab 
        thisassocier = "assPro" + segmentcoll[i] + paramlab 
        thisanalyzer = "me0Ana" + segmentcoll[i] + paramlab 

        thissegmlabel = segmentcoll[i] + "Segments"        
        thismuonlabel = segmentcoll[i] + "SegmentMatching" 


        thisMe0prod = aMe0prod.clone(muonsTag         = cms.InputTag(thismuonlabel), 
                                     me0segmentsTag   = cms.InputTag(thissegmlabel), 
                                     maxPullX         = cms.untracked.double( paramcoll[j][0] ), 
                                     maxDiffX         = cms.untracked.double( paramcoll[j][1] ), 
                                     maxPullY         = cms.untracked.double( paramcoll[j][2] ), 
                                     maxDiffY         = cms.untracked.double( paramcoll[j][3] ), 
                                     maxDiffPhiDir    = cms.untracked.double( paramcoll[j][4] ) 
                                     ) 

        thisMe0asso = aMe0asso.clone(tracksTag        = cms.InputTag(thismuonprod)) 

        thisMe0anal = aMe0anal.clone(associatorNames  = cms.vstring(thisassocier), 
                                     me0segmentsTag   = cms.InputTag(thissegmlabel), 
                                     muonTag          = cms.InputTag(thismuonlabel), 
                                     tracksTag        = cms.InputTag(thismuonprod)) 

        setattr(process,       "tpcPro"    , aTkPprod    ) 
        setattr(process,       thismuonprod, thisMe0prod ) 
        setattr(process,       thisassocier, thisMe0asso ) 
        setattr(process,       thisanalyzer, thisMe0anal ) 
        setattr(process, "RUN"+thisanalyzer, cms.Path(aTkPprod*thisMe0prod*thisMe0asso*thisMe0anal)) 

process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string("TEMPLOUTFILE"), 
                                   closeFileFast = cms.untracked.bool(False)) 

