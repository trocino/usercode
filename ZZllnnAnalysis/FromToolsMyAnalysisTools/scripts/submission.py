#!/usr/bin/env python

import os
import sys
import commands
import shutil
import time
from optparse import OptionParser, Option, OptionValueError
from ConfigParser import ConfigParser
from Tools.MyAnalysisTools.color_tools import *
import FWCore.ParameterSet.Config as cms


def createJobSetups(inputCfg, inputDir, outputDir, outputBaseName, isMC, xSect, bRatio, luminosity, JobName, nFilesPerJob, queue):

    JobName += "/"

    pwd = os.environ["PWD"]

    submissionArea = pwd + "/" + JobName


    outputBaseName = outputBaseName

    if not os.path.exists(JobName):
        print " directory " + JobName + " doesn't exist: creating it"
        os.mkdir(JobName)


    #shutil.copy(inputCfg, "input_cfg.py")
    shutil.copy(inputCfg, JobName + "input_cfg.py")

    os.chdir(submissionArea)
    sys.path.append(submissionArea)
    # list the files in the dir
    castorDir_cmd = "rfdir " + inputDir
    if castorDir_cmd.find(",") != -1:
        castorDir_cmd = castorDir_cmd.replace(",",";rfdir ")
    castorDir_out = commands.getstatusoutput(castorDir_cmd)
    if castorDir_out[0] != 0:
        print castorDir_out[1]
        sys.exit(1)

    # check the output dir
    outCastorDir_cmd = "rfdir " + outputDir
    outCastorDir_out = commands.getstatusoutput(outCastorDir_cmd)
    if outCastorDir_out[0] != 0:
        print outCastorDir_out[1]
        sys.exit(1)



    castorFileList = []
    #storeDir = inputDir.split("cern.ch/cms")[1]
    storeDir = "rfio:" + inputDir
    storeDir = storeDir.rstrip('/')+'/'
    for castorFileLine in castorDir_out[1].split("\n"):
        castorFile = castorFileLine.split()[8]
        if "root" in castorFile and not "histo" in castorFile:

            #print castorFile
            castorFileList.append(storeDir + castorFile)

    print "Input dir: " + inputDir
    print "# fo files: " + str(len(castorFileList))

    toNFiles = len(castorFileList)
    if len(castorFileList) < nFilesPerJob:
        nFilesPerJob = len(castorFileList)




    from input_cfg import process
    process.source.fileNames = cms.untracked.vstring()
    process.maxEvents.input = -1
    #if process.zzllvvAnalyzer != None and process.zzllvvAnalyzer.weight != None:
    #    process.zzllvvAnalyzer.weight = weight
    if "zzllvvAnalyzer" in process.analyzerNames().split():
        process.zzllvvAnalyzer.isMC = isMC
        process.zzllvvAnalyzer.xSection = xSect
        process.zzllvvAnalyzer.branchingRatio = bRatio
        process.zzllvvAnalyzer.luminosity = luminosity
    if JobName.find('only2l2n') != -1:
        process.p = cms.Path(process.llnnFilt * process.baseSeq)
    elif JobName.find('allBut2l2n') != -1:
        process.p = cms.Path(process.llnnAntifilt * process.baseSeq)
    else:
        process.p = cms.Path(process.baseSeq)
    
    # do the manipulatio on output and input files
    indexPart = 0
    indexTot  = 0
    indexJob = 0
    for inputFile in castorFileList:

        process.source.fileNames.append(inputFile)
        indexPart+=1
        indexTot+=1


        #print inputFile
        if indexPart == nFilesPerJob or indexTot == toNFiles:
            print "Writing cfg file for job # " + str(indexJob) + "...."
            outputFileName = outputBaseName + "_" + str(indexJob) + ".root"
            outputFileNameTmp =  outputBaseName + ".root"
            try:
                process.out.fileName = outputFileNameTmp
            except AttributeError:
                print "no output module \"out\" was found..."
            # write the previous cfg
            cfgfilename = "expanded_" + str(indexJob) + "_cfg.py"
            logfilename = "log_" + str(indexJob) + ".log"
            # dump it
            expanded = process.dumpPython()
            expandedFile = file(cfgfilename,"w")
            expandedFile.write(expanded)
            expandedFile.close()
            print "Writing submission script for job # " + str(indexJob) + "...."
            scriptname = "SubmissionJob_" +  str(indexJob) + ".csh"
            scriptfile = open(scriptname, 'w')
            scriptfile.write("#!/bin/tcsh\n")
            scriptfile.write("#BSUB -j " + JobName + "\n")
            scriptfile.write("#BSUB -q " + queue + "\n")
            scriptfile.write("setenv runningDir $PWD\n")
            scriptfile.write("cd " +  submissionArea + "\n")
            scriptfile.write("setenv SCRAM_ARCH slc5_amd64_gcc434 \n")
            scriptfile.write("eval `scram runtime -csh`\n")
            scriptfile.write("cp " + cfgfilename + " $runningDir\n")
            scriptfile.write("cd $runningDir\n")
            scriptfile.write("cmsRun " + cfgfilename + " >& " + logfilename + "\n")
            scriptfile.write("rfcp " + outputFileNameTmp + " " + outputDir + outputFileName + "\n")
            scriptfile.write("rfcp " + logfilename + " " + outputDir + "\n")
            scriptfile.write("#rfcp histograms.root " + outputDir + "/histograms_" +  str(indexJob) + ".root\n")
            scriptfile.write("\n")
            scriptfile.close()
            os.chmod(scriptname, 0777)

            indexJob += 1
            # reset the linst of input files    
            process.source.fileNames = cms.untracked.vstring()
        
            indexPart = 0

    os.chdir(pwd)

if __name__     ==  "__main__":
    # --- set the command line options
    parser = OptionParser()

    parser.add_option("-q", "--queue", dest="queue",
                      help="queue", type="str", metavar="<queue>")
    parser.add_option("-s", "--split", dest="split",
                      help="# files per job", type="int", metavar="<split>", default=100)
    parser.add_option("-c", "--cfg", dest="config",
                      help="configuration file", type="str", metavar="<config>")
    parser.add_option("-i", "--input-dir", dest="inputdir",
                      help="input directory", type="str", metavar="<input dir>")
    parser.add_option("-o", "--output-dir", dest="outputdir",
                      help="output directory", type="str", metavar="<output dir>")
    parser.add_option("-j", "--job-name", dest="jobname",
                      help="job name", type="str", metavar="<job name>")
    parser.add_option("-f", "--file-basename", dest="basename",
                      help="file basename", type="str", metavar="<file basename>")

    parser.add_option("--submit", action="store_true",dest="submit",default=False, help="submit the jobs")
    parser.add_option("--copy", action="store_true",dest="copy",default=False, help="copy the output of the jobs from castor to a local dir (argument)")
    parser.add_option("--status", action="store_true",dest="status",default=False, help="check the status of the jobs")

    parser.add_option("--file", dest="file",
                      help="submission config file", type="str", metavar="<file>")


    (options, args) = parser.parse_args()

    if options.status:

        if options.file != None:
            print 'Reading configuration file from ',options.file
            # read a global configuration file
            cfgfile = ConfigParser()
            cfgfile.optionxform = str

            cfgfile.read([options.file ])
            
            # get the releases currently managed
            listOfJobs = cfgfile.get('General','jobsToSubmit').split(',')
            flag = cfgfile.get('General','selFlag')
            configFile = cfgfile.get('General','configFile')
            filesPerJob = int(cfgfile.get('General','filesPerJob'))
            outputDirBase = cfgfile.get('General','outputDirBase')
            fileBaseName = cfgfile.get('General','fileBaseName')
            queue = cfgfile.get('General','queue')

            for job in listOfJobs:
                print "--- Task name: " + blue(job) + "  (" + flag + "_" + job + ")"
                status = ok("OK")
                ls_cmd = "ls " + flag + "_" + job
                ls_out = commands.getstatusoutput(ls_cmd)
                nJobs = 0
                if ls_out[0] == 0:
                    lslines = ls_out[1].split("\n")
                    for lsline in lslines:
                        if "SubmissionJob" in lsline:
                            nJobs += 1
                else:
                    print ls_out[1]
                print "    # of jobs: " + str(nJobs)
                outdir = outputDirBase + "/" + flag + "/" + job + "/"

                bjobs_cmd = "bjobs -J " + flag + "_" + job
                bjobs_out = commands.getstatusoutput(bjobs_cmd)
                if bjobs_out[0] == 0:
                    if not "not found" in bjobs_out[1]:
                        print bjobs_out[1]
                else:
                    print bjobs_out[1]
                pending = False
                if len(bjobs_out[1].split("\n")) > 1:
                    status = warning("PENDING")
                    pending = True

                print "    Output dir: " + outdir
                rfdir_cmd = "rfdir " + outdir
                nOutFile = 0
                outCastorDir_out = commands.getstatusoutput(rfdir_cmd)
                if outCastorDir_out[0] == 0:
                    castorLines = outCastorDir_out[1].split("\n")
                    if len(castorLines) != 0:
                        for castorFileLine in castorLines:
                            if 'cerminar' in castorFileLine and "root" in castorFileLine:
                                print "        - " + castorFileLine.split()[8]
                                nOutFile += 1
                else:
                    print outCastorDir_out[1]
                print "    # of output files: " + str(nOutFile)
                if nOutFile != nJobs and not pending:
                    status = error("ERROR")
                print "    Status: " + status

                
            sys.exit(0)
                
    if options.copy:
        if options.file == None:
            print "no cfg file specified!"
            sys.exit(1)

        
        # read a global configuration file
        cfgfile = ConfigParser()
        cfgfile.optionxform = str

        print 'Reading configuration file from ',options.file
        cfgfile.read([options.file ])

        # get the releases currently managed
        listOfJobs = cfgfile.get('General','jobsToSubmit').split(',')
        flag = cfgfile.get('General','selFlag')
        configFile = cfgfile.get('General','configFile')
        filesPerJob = int(cfgfile.get('General','filesPerJob'))
        outputDirBase = cfgfile.get('General','outputDirBase')
        fileBaseName = cfgfile.get('General','fileBaseName')
        queue = cfgfile.get('General','queue')

        for job in listOfJobs:
            inputDir = cfgfile.get(job,'inputDir')
            rfdir_cmd = "rfdir " + outputDirBase + "/" + flag + "/" + job + "/"
            outCastorDir_out = commands.getstatusoutput(rfdir_cmd)
            if outCastorDir_out[0] != 0:
                print outCastorDir_out[1]
                sys.exit(1)
            
            filesToCopy = []
            castorFileList = []
            castorLines = outCastorDir_out[1].split("\n")
            if len(castorLines) != 0:
                for castorFileLine in castorLines:
                    if 'root' in castorFileLine:
                        castorFile = castorFileLine.split()[8]
                        if "root" in castorFile and fileBaseName in castorFile:
                            filesToCopy.append(outputDirBase + "/" + flag + "/" + job + "/" + castorFile)
                            
            else:
                print "dir empty..."

            # do the actual copy
            unique = False
            if len(filesToCopy) == 1:
                unique = True

            for filename in filesToCopy:
                cp_cmd = "xrdcp root://castorcms/" + filename + " " + args[0]
                outputFile = ""
                if unique:
                    outputFile = fileBaseName + "_" + job + ".root"
                else:
                    outputFile = fileBaseName + "_" + job + filename.split(fileBaseName)[1]
                cp_cmd += "/" + outputFile
                
                if os.path.exists(args[0] + "/" + outputFile):
                      print "*** Warning, the file: " + args[0] + "/" + outputFile + " already exists!"
                      confirm = raw_input('Overwrite? (y/N)')
                      confirm = confirm.lower() #convert to lowercase
                      if confirm != 'y':
                          continue
                #print cp_cmd
                cp_out = commands.getstatusoutput(cp_cmd)
                print cp_out[1]

        sys.exit(0)
            
    if options.submit:
        #
        pwd = os.environ["PWD"]
        
        if len(args) == 0:
            if options.file == None:
                print "no workflow specified!"
                sys.exit(1)
            else: 
                # read a global configuration file
                cfgfile = ConfigParser()
                cfgfile.optionxform = str

                print 'Reading configuration file from ',options.file
                cfgfile.read([options.file ])

                # get the releases currently managed
                listOfJobs = cfgfile.get('General','jobsToSubmit').split(',')
                flag = cfgfile.get('General','selFlag')
                aqueue = cfgfile.get('General','queue')
                if options.queue != None:
                    aqueue = options.queue
                for job in listOfJobs:
                    args.append(flag + "_" + job)

        if aqueue == None:
            print "no queue specified!"
            sys.exit(1)


        counterSub = 0
        for job in args:
            if not os.path.exists(job):
                print "Dir: " + job + " doesn't exist!"
            else:
                os.chdir(job)
                fileList = os.listdir(".")
                for filename  in fileList:
                    if "SubmissionJob" in filename:
                        print "Submitting: " + filename + "..."
                        submit_cmd = "bsub -q " + aqueue + " -J " + job + " " + filename
                        submit_out = commands.getstatusoutput(submit_cmd)
                        print submit_out[1]
                        time.sleep(30)
                        counterSub += 1
                        if counterSub == 10:
                            print "will wait 500 sec....sorry!"
                            time.sleep(1200)
                            counterSub = 0
                os.chdir(pwd)
                
    elif options.file != None:
        # read a global configuration file
        cfgfile = ConfigParser()
        cfgfile.optionxform = str

        print 'Reading configuration file from ',options.file
        cfgfile.read([options.file ])

        # get the releases currently managed
        listOfJobs = cfgfile.get('General','jobsToSubmit').split(',')
        flag = cfgfile.get('General','selFlag')
        configFile = cfgfile.get('General','configFile')
        filesPerJob = int(cfgfile.get('General','filesPerJob'))
        outputDirBase = cfgfile.get('General','outputDirBase')
        fileBaseName = cfgfile.get('General','fileBaseName')
        queue = cfgfile.get('General','queue')

        for job in listOfJobs:
            inputDir   = cfgfile.get(job,'inputDir')
            isMC       = bool(int(cfgfile.get(job,'isMCsample')))
            xSect      = cfgfile.get(job,'procXsect')
            bRatio     = cfgfile.get(job,'procBranchRatio')
            luminosity = cfgfile.get(job,'luminosity')
            outputDir  = outputDirBase + "/" + flag + "/" + job + "/"
            mkdir_cmd  = "rfmkdir -p " + outputDir
            mkdir_out  =  commands.getstatusoutput(mkdir_cmd)
            print mkdir_out[1]
            createJobSetups(configFile, inputDir, outputDir, fileBaseName, isMC, xSect, bRatio, luminosity, flag + "_" + job, filesPerJob , queue)
    else:

        if options.queue == None:
            print "no queue specified!"
            sys.exit(1)
        if options.config == None:
            print "no cfg specified!"
            sys.exit(1)
        if options.inputdir == None:
            print "no input dir. specified!"
            sys.exit(1)
        if options.outputdir == None:
            print "no output dir. specified!"
            sys.exit(1)
        if options.jobname == None:
            print "no job name specified!"
            sys.exit(1)
        if options.basename == None:
            print "no fine basename specified!"
            sys.exit(1)

        # -----------------------------------------------------------------------
        inputCfg = options.config
        inputDir = options.inputdir
        outputDir = options.outputdir
        outputBaseName = options.basename
        JobName  = options.jobname
        nFilesPerJob = options.split
        queue = options.queue
        # -----------------------------------------------------------------------

        createJobSetups(inputCfg, inputDir, outputDir, outputBaseName, 0, 1.0, 1.0, 1.0, JobName, nFilesPerJob, queue)

sys.exit(0)
    
