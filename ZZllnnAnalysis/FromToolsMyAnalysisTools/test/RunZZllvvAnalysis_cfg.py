import FWCore.ParameterSet.Config as cms
import commands
import os
import FWCore.ParameterSet.VarParsing as VarParsing
options = VarParsing.VarParsing()
options.register('selection',
                "", #default value
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.string,
                "Selection pass")
options.register('sample',
                "", #default value
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.string,
                "Sample")

options.parseArguments()


#---------------------------------------------
baseInputDir = '/castor/cern.ch/cms/store/cmst3/user/cerminar/Analysis/ZZllvv_v01/'
baseOutputDir = '/data/Analysis/ZZllvv_v01/'
#---------------------------------------------

process = cms.Process("ANALYSIS")

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring()#firstRun = cms.untracked.uint32(10)
                            )


#process.source.fileNames.append('/store/cmst3/user/cerminar/Analysis/ZZllvv_v01/sel0/ZJetsPU/cmgTuple_0.root')


process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(-1)
)



process.zzllvvAnalyzer = cms.EDAnalyzer("ZZllvvAnalyzer",
                                        fileName = cms.untracked.string('ZZllvvAnalyzer.root'))





process.analysisSequence = cms.Sequence(process.zzllvvAnalyzer)

process.p = cms.Path(process.analysisSequence)
#process.MessageLogger.cerr.FwkReport.reportEvery = 100



if options.selection and options.sample:

    baseInputDir = '/castor/cern.ch/cms/store/cmst3/user/cerminar/Analysis/ZZllvv_v01/'
    inputDir = baseInputDir + options.selection + "/" + options.sample + "/"



    castorDir_cmd = "rfdir " + inputDir
    castorDir_out = commands.getstatusoutput(castorDir_cmd)
    if castorDir_out[0] != 0:
        print castorDir_out[1]



    storeDir = inputDir.split("cern.ch/cms")[1]
    #storeDir = "rfio://" + inputDir
    for castorFileLine in castorDir_out[1].split("\n"):
        castorFile = castorFileLine.split()[8]
        if "root" in castorFile and not "histo" in castorFile:
            process.source.fileNames.append(storeDir + castorFile)

    outputFile = baseOutputDir + options.selection + "/histos_" + options.sample + ".root"
    process.zzllvvAnalyzer.fileName = outputFile

else :
    process.source.fileNames.append('/store/cmst3/user/cerminar/Analysis/ZZllvv_v01/sel0/ZJetsPU/cmgTuple_0.root')
