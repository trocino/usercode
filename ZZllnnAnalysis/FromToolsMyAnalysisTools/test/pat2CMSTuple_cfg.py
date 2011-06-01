from PhysicsTools.PatAlgos.patTemplate_cfg import *
import FWCore.ParameterSet.Config as cms

#---------------------------------------------

#---------------------------------------------
### rename this process
process.setName_('CMG')


process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
)
process.source.fileNames = cms.untracked.vstring(
#    '/store/cmst3/user/wreece/VJets/3_8_7/020211/ZJetsPU/PF2PAT_MuonFilter_99_1_G8o.root'
        '/store/cmst3/user/psilva/DYJetsToLL/patTuple_97_3_6Px.root'
)

process.out.fileName = cms.untracked.string('file:/tmp/cerminar/cmgTuple.root')


process.load('CMGTools.HtoZZ2l2nu.StandardSelection_cff')

isData = False
doElectron = False
doMuon = True
doJetMet = True

# put together the sequence
if isData:
    process.e = cms.EndPath(process.out)
    if doElectron:
        process.schedule = cms.Schedule(process.leptonPath,process.eePath,process.mumuPath,process.emuPath,process.jetmetPath,process.e)
    else:
        process.schedule = cms.Schedule(process.leptonPath,process.mumuPath,process.jetmetPath,process.e)
else:
    process.e = cms.EndPath(process.saveHistosInRunInfo*process.out)
    if doElectron:
        process.schedule = cms.Schedule(process.genPath,process.leptonPath,process.eePath,process.mumuPath,process.emuPath,process.jetmetPath,process.e)
    else:
        process.schedule = cms.Schedule(process.genPath,process.leptonPath,process.mumuPath,process.jetmetPath,process.e)



print '*** Scheduling the following sequence: '
print process.schedule

### output file for histograms etc
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('histograms.root'))


process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.out.outputCommands = process.outputConfig.outputCommandsHtoZZ

if doElectron:
    process.out.SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('eePath','mumuPath','emuPath') )
else:
    process.out.SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('mumuPath') )

process.options.SkipEvent = cms.untracked.vstring('ProductNotFound')

