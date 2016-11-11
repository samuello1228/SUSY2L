#!/usr/bin/env python
import ROOT
import logging
import shutil
import os

#gErrorIgnoreLevel = kPrint, kInfo, kWarning, kError, kBreak, kSysError, kFatal;
ROOT.gErrorIgnoreLevel = ROOT.kInfo
logging.basicConfig(level=logging.INFO)
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-f', "--inputFiles", help="all input file", default=None)
parser.add_option('-l', "--inputList", help="list of input file", default=None)
parser.add_option('-t', "--inputTag", help="tag of input file", default="test")
parser.add_option('-d', "--inputDir", help="dir of input file", default=None)
parser.add_option('-o', "--outputDir", help="direcotry for output files", default="test")
parser.add_option('-n', "--nevents", type=int, help="number of events to process for all the datasets", default=-1)
parser.add_option("-w", "--overwrite", action='store_true', default=False, help="overwrite previous submitDir")
parser.add_option('-a', "--isMC", type=int, help="type of input sample", default=0)
parser.add_option("--driver", help="select where to run", choices=("direct", "prooflite", "grid", "condor", "LSF"), default="direct")
parser.add_option("--inputDS", help="tag of input data set for grid driver", default=None)
parser.add_option("--outputTag", help="tag of output for grid driver", default='test')
parser.add_option("--shortName", help="shortname", default=None)
parser.add_option("--grl", help="good run list", default=None)
parser.add_option("--nFilesPerNode", type=int, help="number of files per node in LSF batch", default=15)
parser.add_option("-m", "--manualMatch", action='store_true', help="do the manual match or not", default=False)
parser.add_option("--MCTC", action='store_true', help="Whether to use MCTruthClassifier", default=True)
parser.add_option("--test", help="test run", action='store_true', default=False)

(options, args) = parser.parse_args()

import atexit
@atexit.register
def quite_exit():
    ROOT.gSystem.Exit(0)


logging.info("loading packages")
ROOT.gROOT.Macro("$ROOTCOREDIR/scripts/load_packages.C")

if options.overwrite:
    shutil.rmtree(options.outputDir, True)

#Set up the job for xAOD access:
ROOT.xAOD.Init().ignore();

# create a new sample handler to describe the data files we use
logging.info("creating new sample handler")
sh_all = ROOT.SH.SampleHandler()

if options.inputFiles:
    sample = ROOT.SH.SampleLocal(options.inputTag)
    for file in options.inputFiles.split(','): sample.add(file)
    sh_all.add (sample);
elif options.inputList:
    ROOT.SH.readFileList(sh_all, options.inputTag, options.inputList);
elif options.inputDir:
    dir = options.inputDir
    if dir[-1] != '/': dir+='/'
    sample = ROOT.SH.SampleLocal(options.inputTag)
    files = [dir+f for f in os.listdir(options.inputDir) if f.find('.root')!=-1]
    print files
    for file in files: sample.add(file)
    sh_all.add (sample)
elif options.inputDS:
    ROOT.SH.scanDQ2(sh_all, options.inputDS)

# print out the samples we found
logging.info("%d different datasets found scanning all directories", len(sh_all))

# set the name of the tree in our files
sh_all.setMetaString("nc_tree", "CollectionTree")

# this is the basic description of our job
logging.info("creating new job")
job = ROOT.EL.Job()
job.sampleHandler(sh_all)
 
# add our algorithm to the job
logging.info("creating algorithms")
alg = ROOT.TwoElectrons()
alg.outputName = "myOutput"
if options.MCTC:
    alg.matchMode = 2
elif options.manualMatch:
    alg.matchMode = 0
else:
    alg.matchMode = 1

# alg.m_isMC = options.isMC
if options.grl:
    alg.m_grlFile = options.grl
#output
output = ROOT.EL.OutputStream("myOutput")
job.outputAdd(output)
ntuple = ROOT.EL.NTupleSvc("myOutput")
job.algsAdd(ntuple);

logging.info("adding algorithms")
job.algsAdd(alg)

# job options
job.options().setDouble(ROOT.EL.Job.optMaxEvents, options.nevents)

# make the driver we want to use:
# this one works by running the algorithm directly
logging.info("creating driver")
driver = None
if (options.driver == "direct"):
        logging.info("running on direct")
        driver = ROOT.EL.DirectDriver()
        logging.info("submit job")
        driver.submit(job, options.outputDir)
elif (options.driver == "prooflite"):
        logging.info("running on prooflite")
        driver = ROOT.EL.ProofDriver()
        logging.info("submit job")
        driver.submit(job, options.outputDir)
elif (options.driver == "condor"):
        logging.info("running on condor")
        driver = ROOT.EL.CondorDriver()
        job.options().setString(ROOT.EL.Job.optCondorConf, "GetEnv = True");
        job.options().setDouble(ROOT.EL.Job.optFilesPerWorker, 10);
        logging.info("submit job")
#         driver.submit(job, options.outputDir)
        driver.submitOnly(job, options.outputDir)
elif (options.driver == "LSF"):
        logging.info("running on LSF batch")
        driver = ROOT.EL.LSFDriver()
        driver.options().setString (ROOT.EL.Job.optSubmitFlags, "-L /bin/bash");
        driver.shellInit = "export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase && source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh"
        job.options().setDouble(ROOT.EL.Job.optFilesPerWorker, options.nFilesPerNode);
        logging.info("submit job to LSF batch")
        driver.submitOnly(job, options.outputDir)
elif (options.driver == "grid"):
        logging.info("running on Grid")
        driver = ROOT.EL.PrunDriver()
        outname= "user."+os.environ["RUCIO_ACCOUNT"]+"."+ (options.shortName or options.outputTag+".%in:name[2]%.%in:name[3]%")
#         outname= "user."+os.environ["RUCIO_ACCOUNT"]+"."+ (options.shortName or options.outputTag+".%in:name[4]%")
        driver.options().setString("nc_outputSampleName", outname)
        if options.test:
            driver.options().setDouble(ROOT.EL.Job.optGridNFiles, 4)
            driver.options().setDouble(ROOT.EL.Job.optGridNFilesPerJob, 2)
        logging.info("submit job")
        driver.submitOnly(job, options.outputDir)
