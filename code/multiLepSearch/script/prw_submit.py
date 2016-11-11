#cmd to run in shell before runing this script:
#
#setupATLAS
#voms-proxy-init -voms atlas
#asetup AthAnalysisBase,2.4.16,here  (change 2.4.16 with latest athena ver in /afs/cern.ch/atlas/software/builds/AthAnalysisBase)
#lsetup panda
#lsetup pyami

import os
import sys
from subprocess import call

print sys.argv
if len(sys.argv)<2:
  print "usage: python prw_submit.py -g datasetNameList1.txt dataSetNameList2.txt ..."
  print "options"
  print "  -g : submit grid job"
  print "       if not present then just print the commands without actually submiting jobs"
  exit(0)

verTag = "v2"

#look for -g flag (FIXME: should use optParse?)
isDrill = True
for arg in sys.argv:
  if arg=="-g": isDrill=False

#init pyAMI
import pyAMI.client
import pyAMI.atlas.api as AtlasAPI
client = pyAMI.client.Client('atlas')
AtlasAPI.init()


for arg in sys.argv[1:]:
  if arg=="-g": continue
  f = open(arg)
  for aLine in f:
    aLine = aLine.split("#")[0]
    if aLine=="": continue

    dataSetName = aLine.split()[0]
    unskimmedDataSet = None

    ancestors = AtlasAPI.get_dataset_prov(client, dataSetName)
    ancestors = ancestors["node"]
    for aAncestor in ancestors:
      if aAncestor["dataType"] == "AOD" and aAncestor["events"]!="0":
        unskimmedDataSetName = aAncestor["logicalDatasetName"]
        print aAncestor
        break

    if unskimmedDataSetName:
      datasetID = unskimmedDataSetName.split(".")[1]
      outDataSetName = "user.%s.prw.%s.%s" % ( os.environ["USER"], verTag, datasetID)
      cmd = "pathena PileupReweighting/generatePRW_jobOptions.py --inDS %s --outDS %s"
      cmd = cmd % ( unskimmedDataSetName, outDataSetName)

      print cmd
      if not isDrill: call([cmd], shell=True)

if isDrill: 
  print "Drill mode printing commands only. Add -g flag to actually submit grid jobs."
