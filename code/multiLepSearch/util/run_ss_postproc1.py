#!/usr/bin/env python
import ROOT
import logging
import shutil
import os, re

from optparse import OptionParser

parser = OptionParser()
parser.add_option('-f', "--inputFile", help="a .root file or a .txt containing a input file list", default=None)
parser.add_option('-m', "--isMC" , action="store_true", dest="isMC" , help="is MC or not", default=False)

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


alg = ROOT.ssEvtPostProc1()
alg.isData = 0 if options.isMC else 1

qFlipTool = ROOT.ChargeFlipBkgTool("myChargeFlipTool")
qFlipTool.setProperty("InputRatesFileName" , "$ROOTCOREBIN/data/multiLepSearch/root_files/chargeMisID_Zee_data_signal_wSys.root")
#qFlipTool.setProperty("InputRatesHistoName", "hFlipProb")
alg.mChargeFlipBkgTool = qFlipTool

fakeLepTool = ROOT.FakeLepBkgTool("myFakeLepTool")
#fakeLepTool.setProperty("Method", "Matrix")
#fakeLepTool.setProperty("InputFileName"    , "$ROOTCOREBIN/data/multiLepSearch/root_files/RealFakeLepEff_dummy.root")
#fakeLepTool.setProperty("RealeEffHistoName", "RealeEff")
#fakeLepTool.setProperty("RealuEffHistoName", "RealuEff")
#fakeLepTool.setProperty("FakeeEffHistoName", "FakeeEff")
#fakeLepTool.setProperty("FakeuEffHistoName", "FakeuEff")

fakeLepTool.setProperty("Method", "FakeFactor")
fakeLepTool.setProperty("InputFileName"    , "$ROOTCOREBIN/data/multiLepSearch/root_files/fakefactor_2D_Data16.root")
fakeLepTool.setProperty("eFakeFactorHistoName", "h_ff_ele")
#fakeLepTool.setProperty("eFakeFactorHistoName", "h_ff_ele_v2") #this histo has problem of bin error being just the sqrt of bin content
fakeLepTool.setProperty("uFakeFactorHistoName", "h_ff_mu")

alg.mFakeLepBkgTool = fakeLepTool

#alg.inputFile = "/media/DiracHDD/ytchan/SUSYData/NTUP/v7.0.00297730.physics_Main_PP.root"
#alg.execute()

#alg.inputFile = "/media/DiracHDD/ytchan/SUSYData/NTUP/v7.0.00298595.physics_Main_PP.root"
#alg.execute()
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

for aFile in inputFiles:
  alg.inputFile = aFile
  alg.execute()
