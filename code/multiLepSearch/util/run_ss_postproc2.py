#!/usr/bin/env python
import ROOT
import logging
import shutil
import os, re
from multiLepSearch import ssUtil
from optparse import OptionParser

parser = OptionParser()
parser.add_option('-f', "--inputFile", help="input file", default=None)
parser.add_option('-t', "--outTreePrefix", help="prefix for output tree", default="auto")
parser.add_option('-m', "--isMC" , action="store_true", dest="isMC" , help="is MC or not", default=False)
parser.add_option('-s', "--isSig", action="store_true", dest="isSig",  help="is Signal sample or not", default=False)

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

from ROOT.SUSY import CrossSectionDB
xsecDB = CrossSectionDB("%s/data/SUSYTools/mc15_13TeV/"%os.environ["ROOTCOREBIN"] )

sampleDict = ssUtil.loadSampleDict("%s/../multiLepSearch/script/MCBgIDs_2LSS.txt" % os.environ["ROOTCOREBIN"])

def regBDT(alg, BDTName, varSet, weightFile):
  if varSet not in ["ISR", ""]: 
    print "unknown variable set ", varSet
    return 
  alg.BDTNameList.push_back( BDTName )
  alg.varSetList.push_back( varSet )
  alg.weightFileList.push_back( weightFile )

# "Main"
alg = ROOT.ssEvtPostProc2()
  
for (aName, aFormula) in ssUtil.basicBDTVars:
  alg.basicBDTVarNameList.push_back(aName)
  alg.basicBDTVarFormulaList.push_back(aFormula)

for (aName, aFormula) in ssUtil.isrBDTVars:
  alg.isrBDTVarNameList.push_back(aName)
  alg.isrBDTVarFormulaList.push_back(aFormula)

regBDT(alg, "BDT_nonISR_dm20"   , ""   , "/afs/cern.ch/work/y/ychan/SUSY_BDT/v7d11Ana/weights/TMVAClassification_BDT_nonISR_dm20_BDTD.weights.xml")
regBDT(alg, "BDT_nonISR_dm50"   , ""   , "/afs/cern.ch/work/y/ychan/SUSY_BDT/v7d11Ana/weights/TMVAClassification_BDT_nonISR_dm50_BDTD.weights.xml")
regBDT(alg, "BDT_nonISR_dm100"  , ""   , "/afs/cern.ch/work/y/ychan/SUSY_BDT/v7d11Ana/weights/TMVAClassification_BDT_nonISR_dm100_BDTD.weights.xml")
regBDT(alg, "BDT_nonISR_dm100up", ""   , "/afs/cern.ch/work/y/ychan/SUSY_BDT/v7d11Ana/weights/TMVAClassification_BDT_nonISR_dm100up_BDTD.weights.xml")
regBDT(alg,    "BDT_ISR_dm20"   , "ISR", "/afs/cern.ch/work/y/ychan/SUSY_BDT/v7d11Ana/weights/TMVAClassification_BDT_ISR_dm20_BDTD.weights.xml")
regBDT(alg,    "BDT_ISR_dm50"   , "ISR", "/afs/cern.ch/work/y/ychan/SUSY_BDT/v7d11Ana/weights/TMVAClassification_BDT_ISR_dm50_BDTD.weights.xml")
regBDT(alg,    "BDT_ISR_dm100"  , "ISR", "/afs/cern.ch/work/y/ychan/SUSY_BDT/v7d11Ana/weights/TMVAClassification_BDT_ISR_dm100_BDTD.weights.xml")
regBDT(alg,    "BDT_ISR_dm100up", "ISR", "/afs/cern.ch/work/y/ychan/SUSY_BDT/v7d11Ana/weights/TMVAClassification_BDT_ISR_dm100up_BDTD.weights.xml")

inputFiles = []
if options.inputFile.endswith(".root"):
  inputFiles.append( options.inputFile )
else:
  f = open(options.inputFile)
  for aLine in f:
    aLine = aLine.split("#")[0]
    if aLine=="": continue
    aFile = aLine.split()[0]
    inputFiles.append( aFile )

sumWdict = {}
if options.isMC:
  sumWdict.update( ssUtil.getSumWfromNTUPList( options.inputFile ) )

for aFile in inputFiles:

  alg.inputFile = aFile

  sampleID = -1
  if options.isMC:
    match = re.search(".[0-9]{6}.", aFile)
    if not match:
      print "Cannot infer datasetID from inputName %s , skipped" %aFile
      continue
    sampleID = int(match.group()[1:-1])

  #if options.outTreePrefix == "auto":
  #  sType = None
  #  sampleInfo = sampleDict.get( sampleID, None)
  #  if sampleInfo: sType = sampleInfo["phyType"]
  #  if not sType: sType =  ssUtil.guessSampleType(aFile)

  #  if not sType: 
  #     print "Cannot guess sample type from inputName %s skipped" % aFile
  #     continue
  #  alg.outTreePrefix = sType
  #else:
  #  alg.outTreePrefix = options.outTreePrefix

  alg.mcEvtW = 1.0
  if options.isMC:
    mcSumW = sumWdict.get(sampleID, -1)
    xSECxEff = xsecDB.xsectTimesEff(sampleID)
    #if options.isSig: xSECxEff = xsecDB.xsectTimesEff(sampleID, 125) + xsecDB.xsectTimesEff(sampleID, 127)
    if options.isSig: xSECxEff = xsecDB.xsectTimesEff(sampleID, 125) #125 is 125+127 at the DB (for Slep0d95 only)
    if xSECxEff<0:
      print "No xSEC available for inputName %s skipped" % aFile
      continue
    outputLumi = 1.0 #pb-1
    alg.mcEvtW = xSECxEff * outputLumi / mcSumW
    print mcSumW, xSECxEff, sampleID
    
  print "postproc2-ing %s %s" % (alg.outTreePrefix, aFile)
  alg.execute()
