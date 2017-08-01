#!/usr/bin/env python
#this file is based heavily on $ROOTSYS/tutorials/TMVAClassification.C

import ROOT
ROOT.gROOT.Macro("$ROOTCOREDIR/scripts/load_packages.C")
from ROOT import TMVA

import os, re

from optparse import OptionParser
from multiLepSearch import ssUtil

parser = OptionParser()
parser.add_option('-b', "--bkgList" , help="list of BKG input ROOT file", default=None)
parser.add_option('-s', "--sigList" , help="list of SIG input ROOT file", default=None)
parser.add_option('-d', "--dataList", help="list of DATA input ROOT file", default=None)
parser.add_option('-t', "--nominalTreeName", help="name of the nominal tree in input file", default="evt2l")
parser.add_option('-l', "--inputLumi", help="total Luminosity(fb-1) in dataList samples. Affects weight ratio between MC & data driven bkg", 
      type="float", default=33257.2)
parser.add_option('-o', "--outFile", help="name of output file", default="testoutput.root")
parser.add_option('-c', "--channel", help="set channel, options are: 0=ee, 1=emu, 2=emu; +0=noISR, +10=ISR",type="int", default=None)
parser.add_option("--NodeSize", help="minimum node size (in %) of the trees", type="int", default=5)
parser.add_option("--NTrees", help="number of trees in forest", type="int", default=400)
parser.add_option("--Depth", help="tree depth", type="int", default=2)

(options, args) = parser.parse_args()


from ROOT.SUSY import CrossSectionDB
xsecDB = CrossSectionDB("%s/data/SUSYTools/mc15_13TeV/"%os.environ["ROOTCOREBIN"] )

#This loads the TMVA library
TMVA.Tools.Instance()

outputFile = ROOT.TFile( options.outFile, "RECREATE" );
outputTag = options.outFile.split("/")[-1].split(".")[0]
# factory = TMVA.Factory( "TMVAClassification_" + outputTag , outputFile,
#                        "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

factory = TMVA.Factory( "TMVAClassification_" + outputTag , outputFile,
                        "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G:AnalysisType=Classification" );

for (aVar, aFormula) in ssUtil.basicBDTVars:
  factory.AddVariable( "%s := %s" % (aVar, aFormula) ,  'F' )

if (int((options.channel%100)/10) != 0):
  for (aVar, aFormula) in ssUtil.isrBDTVars:
    factory.AddVariable( "%s := %s" % (aVar, aFormula) ,  'F' )

sumWdict = {}
sumWdict.update( ssUtil.getSumWfromNTUPList( options.sigList ) )
sumWdict.update( ssUtil.getSumWfromNTUPList( options.bkgList ) )

openedInFileList = []

def addTreesToTMVA(fileList, isMC, isSig):
  fileList = open(fileList,"r")
  for aLine in fileList:
    elements = aLine.split("#")[0]
    if elements=="": continue
    elements = elements.split()
    infname = elements[0]
    inFile = ROOT.TFile.Open( infname )
    openedInFileList.append(inFile)
  
    treeWeight = 1.0
    print infname

    if isMC:

      match = re.search(".[0-9]{6}.", infname)
      if not match:
        print "Cannot infer datasetID from filename %s , skipped" %ntupName
        continue

      sampleID = int(match.group()[1:-1])
      mcSumW = sumWdict.get(sampleID, -1)
      xSECxEff = -1.
      if isSig:
        #xSECxEff = xsecDB.xsectTimesEff(sampleID, 125) + xsecDB.xsectTimesEff(sampleID, 127)
        xSECxEff = xsecDB.xsectTimesEff(sampleID, 125) #for Slep0d95 I set 125 to be actually 125+127 in SUSYTools/data
      else:
        xSECxEff = xsecDB.xsectTimesEff(sampleID)
      print "mc :", sampleID, mcSumW, xSECxEff

      if isSig:
        treeWeight = 1.0 * options.inputLumi / mcSumW #treat diff SUSY scenario with equal weight
      else:
        treeWeight = xSECxEff * options.inputLumi / mcSumW
  
    if treeWeight<=0 : 
      print "Encounter <=0 weight sample %s , skipped" % infname
      continue
  
    inTree  = inFile.Get( options.nominalTreeName )
    if inTree.GetEntries()==0: continue
    
    if isSig:
      factory.AddSignalTree( inTree, treeWeight )
    else:
      factory.AddBackgroundTree( inTree, treeWeight )

  fileList.close()

addTreesToTMVA( options. sigList, isMC=True , isSig=True)
addTreesToTMVA( options. bkgList, isMC=True , isSig=False)
addTreesToTMVA( options.dataList, isMC=False, isSig=False)

# event-wise weights
factory.SetSignalWeightExpression    ( "ElSF*MuSF*BtagSF*weight*pwt" )
factory.SetBackgroundWeightExpression( "ElSF*MuSF*BtagSF*weight*pwt*(isMC? 1.0 : (qFwt+fLwt))" )

bkgCut = ROOT.TCut( ssUtil.getCut(int((options.channel)%100) ) ) # SS events or fake events from data only 
sigCut = ROOT.TCut( ssUtil.getCut( int(options.channel + (200 if int(options.channel/100)==0 else 0) ) ) ) # Allow OS signal events for training

factory.PrepareTrainingAndTestTree( sigCut, bkgCut,
    "nTrain_Signal=0:nTrain_Background=0:nTest_Background=0:SplitMode=Random:NormMode=EqualNumEvents:!V" )
    #"nTrain_Signal=0:nTrain_Background=2000:SplitMode=Random:NormMode=EqualNumEvents:!V" )

#Here we just use the BDT, see $ROOTSYS/tutorials/tmva/TMVAClassification.C for other available machine learning methods in TMVA
methodOpt = "!H:V:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=N"
methodOpt = "%s:MaxDepth=%d:NTrees=%d:MinNodeSize=%d%%" % (methodOpt, options.Depth, options.NTrees, options.NodeSize)
factory.BookMethod( TMVA.Types.kBDT, "BDTD", methodOpt )

factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()    

outputFile.Close()

print "=== wrote root file %s\n" % options.outFile
print "=== TMVAClassification is done!\n"
