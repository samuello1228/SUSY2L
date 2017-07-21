#!/usr/bin/env python
import ROOT
import logging
import shutil
import os, re

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
parser.add_option('-s', "--doSys", type=int, help="do systematics", default=0)
parser.add_option("--driver", help="select where to run", choices=("direct", "prooflite", "grid", "condor", "LSF"), default="direct")
parser.add_option("--inputDS", help="tag of input data set for grid driver", default=None)
parser.add_option("--outputTag", help="tag of output for grid driver", default='test')
parser.add_option("--shortName", help="shortname", default=None)
parser.add_option("--nFilesPerNode", type=int, help="number of files per node in LSF batch", default=15)
parser.add_option("--grl", help="good run list", default='multiLepSearch/GRL/physics_25ns_20.7.xml')
parser.add_option("--conf", help="selection configuration file", default='multiLepSearch/sel_conf/SUSYTools_multilepAna.txt')
parser.add_option("--dataPRW", help='data pileup reweighting file list', default= 'multiLepSearch/prw_Data/ilumicalc_histograms_None_297730-311481_OflLumi-13TeV-005.root')
parser.add_option("--mcPRW", help='mc pileup reweighting file list', default='dev/SUSYTools/merged_prw_mc15c_latest.root')
parser.add_option("--test", help="test run", action='store_true', default=False)
parser.add_option("--fast", help="Fast submit for grid jobs.", action='store_true', default=False)
parser.add_option("--samplesDir", help="samples dir", default=None)
parser.add_option("--samplePattern", help="sample pattern", default='(.*)')
parser.add_option("--sampleList", help="sample list", default=None)
parser.add_option("--study", help="name of study",choices=("ss", "ssSlim", "3l", "fakes"),default ="3l")
parser.add_option("--mcMatch", help="MC truth match algorithm", choices=("MCTC", "dR", "TruthLink", "reverseTruthLink", "tryAll"), default="dR")
# parser.add_option("--isShortJob", action='store_true', default=False, help="use condor_submit_short")
parser.add_option("--ChargeID", type="int", help="Use ChargeIDSelector", default=1)
parser.add_option("--extraOptions", help="extra options", default=None)


(options, args) = parser.parse_args()

#gErrorIgnoreLevel = kPrint, kInfo, kWarning, kError, kBreak, kSysError, kFatal;
ROOT.gErrorIgnoreLevel = ROOT.kInfo
logging.basicConfig(level=logging.INFO)

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
    sh_all.add(sample);
elif options.inputList:
    if options.driver == 'grid':
        mergeSamples = []
        for ds in options.inputList.split(','):
            with open(ds) as fin1:
                for line in fin1.readlines():
                    if line[0] == '#':
                        print line, 'skipped'
                        continue
                    line = line.split()[0]
                    if options.fast: mergeSamples.append(line.rstrip())
                    else:
                        if line.find('*') == -1: ROOT.SH.addGrid(sh_all, line.rstrip())
                        else: ROOT.SH.scanRucio(sh_all, line.rstrip())
        if mergeSamples:
            ## create a single sample
            sample = ROOT.SH.SampleGrid("gridMergedSample")
            sample.meta().setString(ROOT.SH.MetaFields.gridName,','.join(mergeSamples))
            sample.meta().setString(ROOT.SH.MetaFields.gridFilter, ROOT.SH.MetaFields.gridFilter_default)
            sh_all.add(sample)
    else: 
        ROOT.SH.readFileList(sh_all, options.inputTag, options.inputList);
elif options.inputDir:
    dir0 = options.inputDir
    if dir0[-1] != '/': dir0+='/'
    sample = ROOT.SH.SampleLocal(options.inputTag)
    files = [dir0+f for f in os.listdir(options.inputDir) if f.find('.root')!=-1]
    print dir0, len(files), 'files'
    for file in files: sample.add(file)
    sh_all.add (sample)
elif options.inputDS:
    if options.fast:
        sample = ROOT.SH.SampleGrid("gridSamples")
        sample.meta().setString(ROOT.SH.MetaFields.gridName, options.inputDS)
        sample.meta().setString(ROOT.SH.MetaFields.gridFilter, ROOT.SH.MetaFields.gridFilter_default)
        sh_all.add(sample)
    else:
        if options.inputDS.find('*') == -1:
            ROOT.SH.addGrid(sh_all, options.inputDS)
        else: ROOT.SH.scanRucio(sh_all, options.inputDS)
elif options.samplesDir:
    dir0 = options.samplesDir
    if dir0[-1] != '/': dir0+='/'
    if options.sampleList:
        for samp1 in options.sampleList.split(','):
            with open(samp1) as fin1:
                for line in fin1.readlines():
                    if line[0] == '#':
                        print line, 'skipped'
                        continue
                    s = line.rstrip().split(':')
                    sample = ROOT.SH.SampleLocal(s[0])
                    d = dir0+s[-1]
                    print s[0],d
                    for f in filter(lambda x: x.find('.root')!=-1, os.listdir(d)): sample.add(d+'/'+f)
                    sh_all.add(sample)
    else:
        dirs = [dir0+d for d in os.listdir(dir0) if os.path.isdir(dir0+d)]
        for d in dirs:
            m1 = re.match(options.samplePattern, d)
            if not m1:
                print d, 'excluded'
                continue
#             tag = m1.group(1)
            tag = '_'.join(m1.groups())
            print tag, d,
            sample = ROOT.SH.SampleLocal(tag)
            files = [d+'/'+f for f in os.listdir(d) if f.find('.root')!=-1]
            print len(files), "files"
            for f in files: sample.add(f)
            sh_all.add(sample)
    if options.test:
        exit(0)


# print out the samples we found
logging.info("%d different datasets to be processed", len(sh_all))

# set the name of the tree in our files
sh_all.setMetaString("nc_tree", "CollectionTree")

# this is the basic description of our job
logging.info("creating new job")
job = ROOT.EL.Job()
job.sampleHandler(sh_all)

# add our algorithm to the job
logging.info("creating algorithms")

alg = ROOT.ssEvtSelection()
alg.CF_outputName = "myOutput"

if (options.study == "3l"):
    alg.CF_outputTreeName = "evt3l"
    muonTrig = ["HLT_mu26_imedium", "HLT_mu24_imedium", "HLT_mu24_iloose_L1MU15", "HLT_mu20_iloose_L1MU15", "HLT_mu50", "HLT_mu60_0eta105_msonly"]
    dimuonTrig = ["HLT_2mu14","HLT_2mu10","HLT_mu24_mu8noL1","HLT_mu22_mu8noL1","HLT_mu20_mu8noL1","HLT_mu18_mu8noL1"]
    for i in muonTrig: alg.CF_trigNames.push_back(i)
    for i in dimuonTrig: alg.CF_trigNames.push_back(i)
    alg.study = "3l"

elif(options.study == "ss" or options.study == "ssSlim" or options.study == "fakes" ):
    alg.CF_outputTreeName = "evt2l"

    #trigger
    #https://svnweb.cern.ch/trac/atlasoff/browser/PhysicsAnalysis/DerivationFramework/DerivationFrameworkSUSY/trunk/share/SUSY2.py
    #https://twiki.cern.ch/twiki/bin/viewauth/Atlas/LowestUnprescaled
    
    electronTrig = ["HLT_e24_lhmedium_L1EM20VH","HLT_e24_lhtight_nod0_ivarloose","HLT_e26_lhtight_nod0_ivarloose","HLT_e60_lhmedium_nod0","HLT_e60_medium","HLT_e140_lhloose_nod0","HLT_e300_etcut"]
    for i in electronTrig: alg.CF_trigNames.push_back(i)

    dielectronTrig = ["HLT_2e12_lhloose_L12EM10VH","HLT_2e15_lhvloose_nod0_L12EM13VH","HLT_2e17_lhvloose_nod0"]
    for i in dielectronTrig: alg.CF_trigNames.push_back(i)

    muonTrig = ["HLT_mu20_iloose_L1MU15","HLT_mu24_iloose","HLT_mu24_ivarloose","HLT_mu24_imedium","HLT_mu24_ivarmedium","HLT_mu26_imedium","HLT_mu26_ivarmedium","HLT_mu40","HLT_mu50"]
    for i in muonTrig: alg.CF_trigNames.push_back(i)

    dimuonTrig = ["HLT_2mu10","HLT_mu18_mu8noL1","HLT_mu18_2mu4noL1","HLT_2mu10_nomucomb","HLT_2mu14","HLT_2mu14_nomucomb","HLT_mu20_mu8noL1","HLT_mu22_mu8noL1"]
    for i in dimuonTrig: alg.CF_trigNames.push_back(i)

    elemuonTrig = ["HLT_e17_lhloose_mu14","HLT_e24_lhmedium_L1EM20VHI_mu8noL1","HLT_e7_lhmedium_mu24","HLT_e17_lhloose_nod0_mu14","HLT_e24_lhmedium_nod0_L1EM20VHI_mu8noL1","HLT_e26_lhmedium_nod0_L1EM22VHI_mu8noL1","HLT_e7_lhmedium_nod0_mu24"]
    for i in elemuonTrig: alg.CF_trigNames.push_back(i)
 
    alg.study = options.study

# elif(options.study == "fakes"):
#     alg.CF_outputTreeName = "evt2l"

#     # 2015 triggers
#     electronTrig = ["HLT_e24_lhmedium_L1EM20VH", "HLT_e60_lhmedium"]
#     for i in electronTrig: alg.CF_trigNames.push_back(i)

#     muonTrig = ["HLT_mu20_iloose_L1MU15", "HLT_mu50"]
#     for i in muonTrig: alg.CF_trigNames.push_back(i)

#     # 2016 triggers
#     electronTrig = ["HLT_e24_lhmedium_nod0_L1EM20VH", "HLT_e24_lhtight_nod0_ivarloose", "HLT_e26_lhtight_nod0_ivarloose", "HLT_e60_medium" ]
#     for i in electronTrig: alg.CF_trigNames.push_back(i)

#     muonTrig = ["HLT_mu24_ivarmedium", "HLT_mu50"]
#     for i in muonTrig: alg.CF_trigNames.push_back(i)

#     alg.study = options.study

alg.mcTruthMatch = options.mcMatch if (options.study!="fakes") else "MCTC"
alg.useChargeIDSelector = options.ChargeID

# alg.CF_nLepCutExactly = 2
# alg.CF_nLepCutMin = 2
# alg.CF_nJetCutExactly = -1
# alg.CF_nJetCutMin = -1

#alg.CF_ElPtCut = 20000
alg.CF_ElPtCut = -1
# alg.CF_Eld0SigCut = 5
# alg.CF_Elz0Cut = 0.5

#alg.CF_MuPtCut = 20000
alg.CF_MuPtCut = -1
# alg.CF_Mud0SigCut = 3
# alg.CF_Muz0Cut = 0.5

# alg.CF_JetPtCut = 20000
# alg.CF_JetEtaCut = 2.8
# alg.CF_JetJvtCut = 0.64
alg.CF_derivationName = 'SUSY2'
alg.CF_mT2_m0 = 0

alg.CF_isMC = options.isMC
alg.CF_ConfigFile = options.conf
alg.CF_vxTrkPtMin = -1
if options.grl:
    for x in options.grl.split(','): alg.CF_grlFiles.push_back(x)

alg.doSys = options.doSys

if options.dataPRW:
    for x in options.dataPRW.split(','): alg.CF_PRW_lcalcFiles.push_back(x)
if options.mcPRW:
    for x in options.mcPRW.split(','): alg.CF_PRW_confFiles.push_back(x)

#output -- done inside the alg
# output = ROOT.EL.OutputStream("myOutput")
# job.outputAdd(output)
## not needed
# ntuple = ROOT.EL.NTupleSvc("myOutput")
# job.algsAdd(ntuple);

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
#         if options.isShortJob: driver.shellInit += ' && export PATH=~/.tmp_bin/condor_submit:$PATH'
        #&& source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh && cd .. && rcSetup && cd $HOME";
        #job.options().setString(ROOT.EL.Job.optCondorConf, "GetEnv = True");
        job.options().setDouble(ROOT.EL.Job.optFilesPerWorker, options.nFilesPerNode);
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
        outname= "user."+os.environ["RUCIO_ACCOUNT"]+"."+ (options.shortName or options.outputTag+".%in:name[1]%.%in:name[2]%.%in:name[3]%")
#         outname= "user."+os.environ["RUCIO_ACCOUNT"]+"."+ (options.shortName or options.outputTag+".%in:name[4]%")
        if options.extraOptions:
            # "--allowTaskDuplication"
            driver.options().setString("nc_EventLoop_SubmitFlags", options.extraOptions);
        if options.test:
            driver.options().setDouble(ROOT.EL.Job.optGridNFiles, 4)
            driver.options().setDouble(ROOT.EL.Job.optGridNFilesPerJob, 2)
        if options.fast:
            job.options().setString(ROOT.EL.Job.optSubmitFlags, "--addNthFieldOfInDSToLFN=2,3 --useContElementBoundary");
            outname= "user."+os.environ["RUCIO_ACCOUNT"]+"."+ (options.shortName or options.outputTag)
        driver.options().setString("nc_outputSampleName", outname)
#         driver.options().setDouble("nc_disableAutoRetry", 1)
        logging.info("submit job")
        driver.submitOnly(job, options.outputDir)
