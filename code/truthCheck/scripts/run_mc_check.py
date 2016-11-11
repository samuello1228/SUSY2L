#!/usr/bin/env python
import ROOT
import logging
import shutil
import os, re

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
parser.add_option("--nFilesPerNode", type=int, help="number of files per node in LSF batch", default=15)
parser.add_option("--grl", help="good run list", default=None)
parser.add_option("--test", help="test run", action='store_true', default=False)
parser.add_option("--samplesDir", help="samples dir", default=None)
parser.add_option("--samplePattern", help="sample pattern", default='(.*)')

(options, args) = parser.parse_args()

import atexit
@atexit.register
def quite_exit():
    ROOT.gSystem.Exit(0)

#gErrorIgnoreLevel = kPrint, kInfo, kWarning, kError, kBreak, kSysError, kFatal;
ROOT.gErrorIgnoreLevel = ROOT.kInfo

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
    dir1 = options.inputDir
    if dir1[-1] != '/': dir1+='/'
    sample = ROOT.SH.SampleLocal(options.inputTag)
    files = [dir1+f for f in os.listdir(options.inputDir) if f.find('.root')!=-1]
    print files
    for file in files: sample.add(file)
    sh_all.add(sample)
elif options.inputDS:
    ROOT.SH.scanDQ2(sh_all, options.inputDS)
elif options.samplesDir:
    dir0 = options.samplesDir
    if dir0[-1] != '/': dir0+='/'
    dirs = [dir0+d for d in os.listdir(dir0) if os.path.isdir(dir0+d)]
    for d in dirs:
        m1 = re.match(options.samplePattern, d)
        if not m1:
            print d, 'excluded'
            continue
        tag = m1.group(1)
        print tag, d
        sample = ROOT.SH.SampleLocal(tag)
        files = [d+'/'+f for f in os.listdir(d) if f.find('.root')!=-1]
        print files
        for f in files: sample.add(f)
        sh_all.add(sample)
    if options.test:
        exit(0)
    

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

alg = ROOT.mcChecker()
# alg.outputTreeName = "evt2l"
# alg.nLepCutExactly = 2
# alg.nLepCutMin = 2
# alg.nJetCutExactly = -1
# alg.nJetCutMin = -1
# 
# alg.ElPtCut = 20000
# alg.Eld0SigCut = 5
# alg.Elz0Cut = 0.5
# 
# alg.MuPtCut = 10000
# alg.Mud0SigCut = 3
# alg.Muz0Cut = 0.5
# 
# alg.JetPtCut = 20000
# alg.JetEtaCut = 2.8
# alg.JetJvtCut = 0.64
# 
# #trigger
# muonTrig = ["HLT_mu26_imedium", "HLT_mu24_imedium", "HLT_mu24_iloose_L1MU15", "HLT_mu20_iloose_L1MU15", "HLT_mu50", "HLT_mu60_0eta105_msonly"]
# dimuonTrig = ["HLT_2mu14","HLT_2mu10","HLT_mu24_mu8noL1","HLT_mu22_mu8noL1","HLT_mu20_mu8noL1","HLT_mu18_mu8noL1"]
# trimuonTrig = ["HLT_mu24_2mu4noL1","HLT_mu22_2mu4noL1","HLT_mu20_2mu4noL1","HLT_mu18_2mu4noL1","HLT_3mu6","HLT_3mu6_msonly"]
# 
# electronTrigCut = ["HLT_e26_tight_iloose","HLT_e60_medium","HLT_e24_tight_iloose","HLT_e24_medium_iloose_L1EM20VH","HLT_e24_medium_iloose_L1EM18VH"]
# electronTrigLH = ["HLT_e26_lhtight_iloose","HLT_e60_lhmedium","HLT_e24_lhtight_iloose","HLT_e24_lhmedium_iloose_L1EM20VH","HLT_e24_lhmedium_iloose_L1EM18VH"]
# 
# dielectronTrigCut = ["HLT_2e17_loose","HLT_2e15_loose_L12EM13VH","HLT_2e12_loose_L12EM10VH"]
# dielectronTrigLH = ["HLT_2e17_lhloose","HLT_2e15_lhloose_L12EM13VH","HLT_2e12_lhloose_L12EM10VH"]
# 
# trielectronTrigCut = ["HLT_e17_medium_2e9_medium"]
# trielectronTrigLH = ["HLT_e17_lhmedium_2e9_lhmedium"]
# 
# elemuonTrigCut = ["HLT_e17_loose_mu14","HLT_e7_medium_mu24","HLT_e26_medium_L1EM22VHI_mu8noL1","HLT_e24_medium_L1EM20VHI_mu8noL1"]
# elemuonTrigLH = ["HLT_e17_lhloose_mu14","HLT_e7_lhmedium_mu24","HLT_e26_lhmedium_L1EM22VHI_mu8noL1","HLT_e24_lhmedium_L1EM20VHI_mu8noL1"]
# 
# trielemuonTrigCut = ["HLT_e12_loose_2mu10","HLT_2e12_loose_mu10"]
# trielemuonTrigLH = ["HLT_e12_lhloose_2mu10","HLT_2e12_lhloose_mu10"]
# 
# for i in muonTrig: alg.trigNames.push_back(i)
# for i in dimuonTrig: alg.trigNames.push_back(i)
# for i in trimuonTrig: alg.trigNames.push_back(i)
# 
# for i in electronTrigCut: alg.trigNames.push_back(i)
# for i in electronTrigLH: alg.trigNames.push_back(i)
# 
# for i in dielectronTrigCut: alg.trigNames.push_back(i)
# for i in dielectronTrigLH: alg.trigNames.push_back(i)
# 
# for i in trielectronTrigCut: alg.trigNames.push_back(i)
# for i in trielectronTrigLH: alg.trigNames.push_back(i)
# 
# for i in elemuonTrigCut: alg.trigNames.push_back(i)
# #for i in elemuonTrigLH: alg.trigNames.push_back(i)
# 
# for i in trielemuonTrigCut: alg.trigNames.push_back(i)
# for i in trielemuonTrigLH: alg.trigNames.push_back(i)
# 
# alg.m_isMC = options.isMC
# if options.grl:
#     alg.m_grlFile = options.grl
#output
alg.outputName = "myOutput"
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
        driver.shellInit = "export HOME=$PWD && export PATH=$PATH:$PWD &&  export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase";
        #&& source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh && cd .. && rcSetup && cd $HOME";
        #job.options().setString(ROOT.EL.Job.optCondorConf, "GetEnv = True");
        if options.nFilesPerNode > 0: job.options().setDouble(ROOT.EL.Job.optFilesPerWorker, options.nFilesPerNode);
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
#         driver.options().setDouble("nc_disableAutoRetry", 1)
        logging.info("submit job")
        driver.submitOnly(job, options.outputDir)
