"""
 **********************************************************************************
 * Project: HistFitter - A ROOT-based package for statistical data analysis       *
 * Package: HistFitter                                                            *
 *                                                                                *
 * Description:                                                                   *
 *      Simple example configuration using a single bin and raw numbers           * 
 *                                                                                *
 * Authors:                                                                       *
 *      HistFitter group, CERN, Geneva                                            *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in the file          *
 * LICENSE.                                                                       *
 **********************************************************************************
"""

################################################################
## In principle all you have to setup is defined in this file ##
################################################################

## This configuration performs a trivial one-bin fit
## Only two systematics are considered:
##   -JES (Tree-based)
##   -Alpgen Kt scale (weight-based)
##


from configManager import configMgr
from ROOT import kBlack,kWhite,kGray,kRed,kPink,kMagenta,kViolet,kBlue,kAzure,kCyan,kTeal,kGreen,kSpring,kYellow,kOrange
from configWriter import fitConfig,Measurement,Channel,Sample
from systematic import Systematic
from math import sqrt

from ROOT import gROOT
#gROOT.LoadMacro("./macros/AtlasStyle.C")
import ROOT
#ROOT.SetAtlasStyle()

import re
import sys

from multiLepSearch import ssUtil

from optparse import OptionParser
myUserArgs= configMgr.userArg.split(' ')
myInputParser = OptionParser()
myInputParser.add_option('-b', "--bkgList" , help="list of BKG input ROOT file", default=None)
myInputParser.add_option('-s', "--sigList" , help="list of SIG input ROOT file", default=None)
myInputParser.add_option('-d', "--dataList", help="list of DATA input ROOT file", default=None)
myInputParser.add_option('-v', "--bdtCutVal", type="float", help="BDT cut value to use", default=0.1)
myInputParser.add_option('-c', "--channel"  , help="Channel (ISR or nonISR)", default="nonISR")
myInputParser.add_option('-m', "--dm"       , help="use BDT trained with for dm", default="50")
myInputParser.add_option('-l', "--bkgLumi"     , type="float", help="lumi of data driven background", default=33.2572)
(options, args) = myInputParser.parse_args(myUserArgs)

#-------------------------------
# Parameters for hypothesis test
#-------------------------------
#configMgr.doHypoTest=False
configMgr.nTOYs=5000
configMgr.calculatorType=2 ## 2=asymptotic calculator, 0=frequentist calculator
configMgr.testStatType=3 # one-sided test profile statistic (ATLAS standard)
#configMgr.testStatType=2
configMgr.nPoints=20


#-------------------------------------
# Now we start to build the data model
#-------------------------------------

# First define HistFactory attributes
configMgr.analysisName = "My2LSSAnalysis_%s_dm%s_Cut%s" % (options.channel, options.dm, str(options.bdtCutVal).replace(".","d"))
configMgr.histCacheFile = "data/"+configMgr.analysisName+".root"
configMgr.outputFileName = "results/"+configMgr.analysisName+"_Output.root"

# Scaling calculated by outputLumi / inputLumi
configMgr.inputLumi  = 0.001 # Luminosity of input TTree after weighting
#configMgr.inputLumi  = 3.248 # Luminosity of input TTree after weighting
#configMgr.outputLumi = 3.248 # Luminosity required for output histograms
configMgr.outputLumi = 33.2572 # Luminosity required for output histograms
#configMgr.outputLumi = 0.1 # Luminosity required for output histograms
configMgr.setLumiUnits("fb-1")
configMgr.blindSR=True       #    !!!!!
configMgr.blindCR=False       #    !!!!!
configMgr.useSignalInBlindedData=False  #if true the fake data that replace real data in blinded mode will include both sig and bkg

# Set the files to read from
bgdFiles = []
if configMgr.readFromTree:
    pass
else:
    bgdFiles = [configMgr.histCacheFile]
    pass
configMgr.setFileList(bgdFiles)

tarCutVals = [options.bdtCutVal]
tarDms     = [options.dm]
tarChs     = [options.channel]

for cutVal in tarCutVals:
  cutValStr = str(cutVal).replace('.','d')
  for dm in tarDms:
    for ch in tarChs:
      #SR
      cutName =  "SR_%s_dm%s_Cut%s" % (ch, dm, cutValStr)
      bdtCut  = "BDT_%s_dm%s > %e"  % (ch, dm, cutVal)
      basicCut =  ssUtil.getCut(useISR= (ch=="ISR")  )
      configMgr.cutsDict[cutName] = "(%s) && (%s)" % (basicCut , bdtCut)
      #CR, reverse BDT region
      cutName =  "CR_%s_dm%s" % (ch, dm )
      bdtCut  = "BDT_%s_dm%s < %e"  % (ch, dm, 0.0)
      basicCut =  ssUtil.getCut(useISR= (ch=="ISR")  )
      configMgr.cutsDict[cutName] = "(%s) && (%s)" % (basicCut , bdtCut)

configMgr.cutsDict["CR_Zpeak"] = "(!(%s)) && (%s) && (%s)" % (ssUtil.zMassCut, ssUtil.trigCut, ssUtil.sigLepSSWithDataBkgCut)

zeeMassCut   = "((int(abs(leps.ID[0])/1000)==11) && (int(abs(leps.ID[1])/1000)==11) && fabs(l12.m - 91.1876)<=10)"
zuuMassCut   = "((int(abs(leps.ID[0])/1000)==13) && (int(abs(leps.ID[1])/1000)==13) && fabs(l12.m - 91.1876)<=10)"
#configMgr.cutsDict["CR_Zeepeak"] = "(%s) && (%s) && (%s)" % (zeeMassCut, ssUtil.trigCut, ssUtil.sigLepSSWithDataBkgCut)
#configMgr.cutsDict["CR_Zeepeak"] = "&&".join(["(%s)"%cut for cut in [zeeMassCut, ssUtil.trigCut]])
#configMgr.cutsDict["CR_Zuupeak"] = "&&".join(["(%s)"%cut for cut in [zuuMassCut, ssUtil.trigCut]])
configMgr.cutsDict["CR_Zeepeak"] = "&&".join(["(%s)"%cut for cut in [zeeMassCut, ssUtil.trigCut, ssUtil.exact2LepCut, ssUtil.lepptCut, ssUtil.l12mCut, ssUtil.sigLepSSWithDataBkgCut]])
configMgr.cutsDict["CR_Zuupeak"] = "&&".join(["(%s)"%cut for cut in [zuuMassCut, ssUtil.trigCut, ssUtil.exact2LepCut, ssUtil.lepptCut, ssUtil.l12mCut, ssUtil.sigLepSSWithDataBkgCut]])




# Tuples of nominal weights without and with b-jet selection
#configMgr.weights = ("genWeight","weight","leptonWeight","triggerWeight","truthWptWeight","bTagWeight2Jet")
configMgr.weights = ("ElSF","MuSF", "BtagSF","weight","pwt", "(isMC? mcEvtW : 1.0)", "((%s)&&(%s))" % (ssUtil.sigLepCut, ssUtil.ssCut))

#--------------------
# List of systematics
#--------------------

reconSys = []
# weight-based
# None, as we still have a tree for weight-based sys

#EL_EFF_ID_UWeights      = ("genWeight","weight","EL_EFF_ID_TotalCorrUncertainty__1up"  ,"MuSF")
#EL_EFF_ID_DWeights      = ("genWeight","weight","EL_EFF_ID_TotalCorrUncertainty__1down","MuSF")
#EL_EFF_ID_Sys           = Systematic("EL_EFF_ID",configMgr.weights,EL_EFF_ID_UWeights,EL_EFF_ID_DWeights,"weight","overallSys")

# tree-based
def formTreeSys(sysName, upSuffix, dnSuffix, sysType):
  #syntax: Systematics( sysName, nominalTreeSuffix, 1sigmaPpSuffix, 1sigmaDownSuffix, "tree"/"weight", sysType)
  return Systematic( sysName, "", "_"+sysName+upSuffix, "_"+sysName+dnSuffix, "tree", sysType)

reconSys = [
       #     formTreeSys("EG_RESOLUTION_ALL"                              , "__1up", "__1down", "overallSys"), 
       #     formTreeSys("EG_SCALE_ALL"                                   , "__1up", "__1down", "overallSys"), 
       #     formTreeSys("EL_EFF_ID_TOTAL_1NPCOR_PLUS_UNCOR"              , "__1up", "__1down", "overallSys"), 
       #     formTreeSys("EL_EFF_Iso_TOTAL_1NPCOR_PLUS_UNCOR"             , "__1up", "__1down", "overallSys"), 
       #     formTreeSys("EL_EFF_Reco_TOTAL_1NPCOR_PLUS_UNCOR"            , "__1up", "__1down", "overallSys"), 
       #     formTreeSys("EL_EFF_TriggerEff_TOTAL_1NPCOR_PLUS_UNCOR"      , "__1up", "__1down", "overallSys"), 
       #     formTreeSys("EL_EFF_Trigger_TOTAL_1NPCOR_PLUS_UNCOR"         , "__1up", "__1down", "overallSys"), 
       #     formTreeSys("FT_EFF_B_systematics"                           , "__1up", "__1down", "overallSys"), 
       #     formTreeSys("FT_EFF_C_systematics"                           , "__1up", "__1down", "overallSys"), 
       #     formTreeSys("FT_EFF_Light_systematics"                       , "__1up", "__1down", "overallSys"), 
       #     formTreeSys("FT_EFF_extrapolation"                           , "__1up", "__1down", "overallSys"), 
       #     formTreeSys("FT_EFF_extrapolation_from_charm"                , "__1up", "__1down", "overallSys"), 
       #    #formTreeSys("JET_EtaIntercalibration_NonClosure"             , "__1up", "__1down", "overallSys"), 
       #    #formTreeSys("JET_GroupedNP_1"                                , "__1up", "__1down", "overallSys"), 
       #    #formTreeSys("JET_GroupedNP_2"                                , "__1up", "__1down", "overallSys"), 
       #    #formTreeSys("JET_GroupedNP_3"                                , "__1up", "__1down", "overallSys"), 
       #    #formTreeSys("JET_JER_SINGLE_NP"                              , "__1up", "__1up"  , "histoSysOneSide"), #!!
       #    #formTreeSys("JET_Rtrk_Baseline_Kin"                          , "__1up", "__1down", "overallSys"), 
       #    #formTreeSys("JET_Rtrk_Baseline_Sub"                          , "__1up", "__1down", "overallSys"), 
       #    #formTreeSys("JET_Rtrk_Modelling_Kin"                         , "__1up", "__1down", "overallSys"), 
       #    #formTreeSys("JET_Rtrk_Modelling_Sub"                         , "__1up", "__1down", "overallSys"), 
       #    #formTreeSys("JET_Rtrk_TotalStat_Kin"                         , "__1up", "__1down", "overallSys"), 
       #    #formTreeSys("JET_Rtrk_TotalStat_Sub"                         , "__1up", "__1down", "overallSys"), 
       #    #formTreeSys("JET_Rtrk_Tracking_Kin"                          , "__1up", "__1down", "overallSys"), 
       #    #formTreeSys("JET_Rtrk_Tracking_Sub"                          , "__1up", "__1down", "overallSys"), 
       #     formTreeSys("JvtEfficiency"                                  ,    "Up",    "Down", "overallSys"),      #!!
       #     formTreeSys("MET_SoftTrk_ResoPara"                           ,      "",        "", "histoSysOneSide"), #!!
       #     formTreeSys("MET_SoftTrk_ResoPerp"                           ,      "",        "", "histoSysOneSide"), #!!
       #     formTreeSys("MET_SoftTrk_Scale"                              ,    "Up",    "Down", "overallSys"),      #!!
       #     formTreeSys("MUONS_ID"                                       , "__1up", "__1down", "overallSys"), 
       #     formTreeSys("MUONS_MS"                                       , "__1up", "__1down", "overallSys"), 
       #     formTreeSys("MUONS_SCALE"                                    , "__1up", "__1down", "overallSys"), 
       #     formTreeSys("MUON_EFF_STAT"                                  , "__1up", "__1down", "overallSys"), 
       #     formTreeSys("MUON_EFF_STAT_LOWPT"                            , "__1up", "__1down", "overallSys"), 
       #     formTreeSys("MUON_EFF_SYS"                                   , "__1up", "__1down", "overallSys"), 
       #     formTreeSys("MUON_EFF_SYS_LOWPT"                             , "__1up", "__1down", "overallSys"), 
       #     formTreeSys("MUON_EFF_TrigStatUncertainty"                   , "__1up", "__1down", "overallSys"), 
       #     formTreeSys("MUON_ISO_STAT"                                  , "__1up", "__1down", "overallSys"), 
       #     formTreeSys("MUON_ISO_SYS"                                   , "__1up", "__1down", "overallSys"), 
       #     formTreeSys("PH_EFF_ID_Uncertainty"                          , "__1up", "__1down", "overallSys"), 
       #     formTreeSys("PH_EFF_TRKISO_Uncertainty"                      , "__1up", "__1down", "overallSys"), 
       #     formTreeSys("PH_Iso_DDonoff"                                 ,      "",        "", "histoSysOneSide"), #!!
       #     formTreeSys("PRW_DATASF"                                     , "__1up", "__1down", "overallSys"), 
       #    #formTreeSys("TAUS_TRUEELECTRON_EFF_ELEOLR_TOTAL"             , "__1up", "__1down", "overallSys"), 
       #    #formTreeSys("TAUS_TRUEHADTAU_EFF_ELEOLR_TOTAL"               , "__1up", "__1down", "overallSys"), 
       #    #formTreeSys("TAUS_TRUEHADTAU_EFF_JETID_HIGHPT"               , "__1up", "__1down", "overallSys"), 
       #    #formTreeSys("TAUS_TRUEHADTAU_EFF_JETID_TOTAL"                , "__1up", "__1down", "overallSys"), 
       #    #formTreeSys("TAUS_TRUEHADTAU_EFF_RECO_HIGHPT"                , "__1up", "__1down", "overallSys"), 
       #    #formTreeSys("TAUS_TRUEHADTAU_EFF_RECO_TOTAL"                 , "__1up", "__1down", "overallSys"), 
       #    #formTreeSys("TAUS_TRUEHADTAU_EFF_TRIGGER_STATDATA2015"       , "__1up", "__1down", "overallSys"), 
       #    #formTreeSys("TAUS_TRUEHADTAU_EFF_TRIGGER_STATMC2015"         , "__1up", "__1down", "overallSys"), 
       #    #formTreeSys("TAUS_TRUEHADTAU_EFF_TRIGGER_SYST2015"           , "__1up", "__1down", "overallSys"), 
       #    #formTreeSys("TAUS_TRUEHADTAU_EFF_TRIGGER_TOTAL2016"          , "__1up", "__1down", "overallSys"), 
       #    #formTreeSys("TAUS_TRUEHADTAU_SME_TES_DETECTOR"               , "__1up", "__1down", "overallSys"), 
       #    #formTreeSys("TAUS_TRUEHADTAU_SME_TES_INSITU"                 , "__1up", "__1down", "overallSys"), 
       #    #formTreeSys("TAUS_TRUEHADTAU_SME_TES_MODEL"                  , "__1up", "__1down", "overallSys"), 
           ]

qFlipSys = formTreeSys("CFLIP_SYS"                  , "__1up", "__1down", "overallSys")
fLepSys  = formTreeSys("FAKE_LEP_E_SYS"             , "__1up", "__1down", "overallSys")
fLepSys  = formTreeSys("FAKE_LEP_U_SYS"             , "__1up", "__1down", "overallSys")

#the input Lumi for mc is 0.001 but data driven bkg is obtained from 3.248fb-1 data
#so need a scale factor to convert data driven bkg back to 0.001 inputLumi...
dataDrivenBkgScale = "%e" % (configMgr.inputLumi/options.bkgLumi)

configMgr.nomName = ""

#-------------------------------------------
# List of samples and their plotting colours
#-------------------------------------------
wgammaSample   = Sample("wgamma" ,kGreen-9)
dibosonSample  = Sample("diboson",kAzure+1)
tribosonSample = Sample("triboson",kAzure+2)
topXSample     = Sample("topX" ,kRed+1)
higgsSample    = Sample("higgs" ,kRed+2)

qFlipSample   = Sample("qFlip" ,kBlue )
fLepSample    = Sample("fLep" ,kYellow )

#cflipSample.setNormFactor("mu_CFlip",1.,0.7,1.3,True) #no sys for CFlip yet, use 30% uncertain non-constant normFac
#fakelepSample.setNormFactor("mu_FakeLep",1.,0.7,1.3,True) #no sys for CFlip yet, use 30% uncertain non-constant normFac

dataSample = Sample("data",kBlack)
dataSample.setData()
#dataSample.buildHisto([1.0],"SR","cuts",0.5)
#dataSample.setStatConfig(False)

sampleDict = ssUtil.loadSampleDict("%s/../multiLepSearch/script/MCBgIDs_2LSS.txt" % os.environ["ROOTCOREBIN"])

wgammaFiles   = []
dibosonFiles  = []
tribosonFiles = []
topXFiles     = []
higgsFiles    = []
dataFiles     = []

f = open(options.bkgList).readlines()
f+= open(options.dataList).readlines()
for aLine in f:
  #remove and ignore comments in line
  aLine = aLine.split("#")[0] # content after # are considered commnets
  if len(aLine)==0: continue

  #read the line
  fName = aLine.split()[0]

  sampleType = None

  #infer type from datasetID
  match = re.search(".[0-9]{6}.", fName)
  if match:
    sampleID = int(match.group()[1:-1])
    sampleInfo = sampleDict.get( sampleID, None)
    if sampleInfo: sampleType = sampleInfo["phyType"]

  #infer type from dataset name
  if not sampleType:
    sampleType = ssUtil.guessSampleType(fName)

  if   sampleType==  "wgamma":   wgammaFiles.append(fName)
  elif sampleType== "diboson":  dibosonFiles.append(fName)
  elif sampleType=="triboson": tribosonFiles.append(fName)
  elif sampleType==    "topX":     topXFiles.append(fName)
  elif sampleType==   "higgs":    higgsFiles.append(fName)
  elif sampleType==    "data":     dataFiles.append(fName)
  else: print "Skipped unrecognized sample : %s" % fName

wgammaSample.setFileList(wgammaFiles)
dibosonSample.setFileList(dibosonFiles)
tribosonSample.setFileList(tribosonFiles)
topXSample.setFileList(topXFiles)
higgsSample.setFileList(higgsFiles)

qFlipSample.setFileList(dataFiles)
fLepSample.setFileList(dataFiles)
dataSample.setFileList(dataFiles)

for aSample in [wgammaSample, dibosonSample, tribosonSample, topXSample, higgsSample]:
  aSample.setNormByTheory(True)
  aSample.setStatConfig(True)
  for aSys in reconSys: aSample.addSystematic(aSys)

qFlipSample.addSystematic(qFlipSys)
qFlipWeights = ("ElSF","MuSF", "BtagSF","weight","pwt", "qFwt", dataDrivenBkgScale)
qFlipSample.setWeights(qFlipWeights) #overrides configMgr.weights
qFlipSample.setNormByTheory(False)
qFlipSample.setStatConfig(True)

fLepSample.addSystematic(fLepSys)
fLepWeights = ("ElSF","MuSF", "BtagSF","weight","pwt", "fLwt", dataDrivenBkgScale)
fLepSample.setWeights(fLepWeights) #overrides configMgr.weights
fLepSample.setNormByTheory(False)
fLepSample.setStatConfig(True)

dataWeights = ("((%s)&&(%s))" % (ssUtil.sigLepCut, ssUtil.ssCut), dataDrivenBkgScale)
dataSample.setWeights(dataWeights) #need to edit python/configManager.setWeightsCutsVariable() in HistFitter for this to take effect

#**************
# Bkg only fit
#**************

if myFitType==FitType.Background:
   print "background only fit not implemented"


#**************
# Discovery fit
#**************

if myFitType==FitType.Discovery:
   print "discovery fit not implemented"

#**************
# Exclusion fit
#**************
if myFitType==FitType.Exclusion:

    from fitConfig import fitConfig
    exclusionFitConfig = configMgr.addFitConfig("ExclusionTemplate")
    meas=exclusionFitConfig.addMeasurement(name="NormalMeasurement",lumi=1.0,lumiErr=0.029)
    meas.addPOI("mu_SIG")
    
    # Set stuff common to all fitConfigs
    exclusionFitConfig.addSamples([wgammaSample,dibosonSample,tribosonSample,topXSample,higgsSample,qFlipSample,fLepSample,dataSample])
    #exclusionFitConfig.addSamples([wgammaSample,dibosonSample,tribosonSample,topXSample,higgsSample,dataSample])
    #exclusionFitConfig.addSamples([dibosonSample,ttbarSample,cflipSample,fakelepSample,dataSample])
    exclusionFitConfig.setTreeName("BDT_PP1_evt2l")
    
    cutName = "CR_%s_dm%s"   % (options.channel, options.dm )
    varName = "BDT_%s_dm%s"  % (options.channel, options.dm )
    crBDTBin  = exclusionFitConfig.addChannel(   varName ,[cutName], 20, -0.5,   0.5)
    exclusionFitConfig.addBkgConstrainChannels([crBDTBin])

    #crZpeakBin = exclusionFitConfig.addChannel("l12.m" ,["CR_Zpeak"], 20, 80,   100)
    #exclusionFitConfig.addBkgConstrainChannels([crZpeakBin])
    crZeepeakBin = exclusionFitConfig.addChannel("l12.m" ,["CR_Zeepeak"], 20, 80,   100)
    crZuupeakBin = exclusionFitConfig.addChannel("l12.m" ,["CR_Zuupeak"], 20, 80,   100)
    exclusionFitConfig.addBkgConstrainChannels([crZeepeakBin, crZuupeakBin])

    #vrl12mBin = exclusionFitConfig.addChannel("l12.m"          ,["CR"],100, 50.0, 150.0)
    #exclusionFitConfig.setValidationChannels([vrl12mBin])


    f = open(options.sigList)
    for aLine in f:
      #remove and ignore comments in line
      aLine = aLine.split("#")[0] # content after # are considered commnets
      if len(aLine)==0: continue
    
      #read the line
      fName = aLine.split()[0]

      sigName = ssUtil.guessSampleType(fName)
      if (not sigName):
        print "Cannot infer signal type from file name, skipping " + sigfName
	continue

      print "Add fitConfig for " + sigName
      for cutVal in tarCutVals:
        cutValStr = str(cutVal).replace('.','d')
        for dm in tarDms:
          for ch in tarChs:
            cutName  = "SR_%s_dm%s_Cut%s" % (ch, dm, cutValStr)
            aSigExclusionFitConfig = configMgr.addFitConfigClone( exclusionFitConfig, "Exclusion_%s_%s_dm%s_%s"% (sigName, ch, dm, cutValStr))
    
            srBin    = aSigExclusionFitConfig.addChannel("cuts",[cutName],1,0.5,1.5)
            aSigExclusionFitConfig.addSignalChannels([srBin])
            #aSigExclusionFitConfig.setSignalChannels([srnonISRBin, srISRBin])

            sigSample = Sample( sigName, kCyan)
            sigSample.setFileList([fName])

            sigSample.setNormByTheory()
            sigSample.setNormFactor("mu_SIG",1.,0.,5.)                    

            #set PRW off and overrides configMgr.weights
            #sigWeights = ("ElSF","MuSF", "BtagSF","weight","1.0", "mcEvtW")
            #sigSample.setWeights(sigWeights) 

            aSigExclusionFitConfig.addSamples(sigSample)
            aSigExclusionFitConfig.setSignalSample(sigSample)

            break

    configMgr.removeFitConfig("ExclusionTemplate")
