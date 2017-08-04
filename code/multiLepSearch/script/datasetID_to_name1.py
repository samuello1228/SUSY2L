#!/usr/bin/env python

#Usage:
#before runing this script, you have to setup pyAMi env by:
#source /afs/cern.ch/atlas/software/tools/pyAMI/setup.sh
#
#then:
#./datasetID_to_name.py inputIDs.csv
#See MCFakesID_WH.csv for format

import sys, csv

#use pyAMI to query datasets, see: https://ami.in2p3.fr/pyAMI/pyAMI5_atlas_api.html
import pyAMI.client
import pyAMI.atlas.api as AtlasAPI

def Data15_criteria(tags):
  """
  return true if the tag satisfy Data15 athena 20.7 reprocessing tags as in
  https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/SUSYxAODDerivationsr20#Physics_containers_for_data_peri
  """

  #This DAOD tag is recommended for athena20.7
  tarTags = ["p2667"]
  if not any( [ t in tags for t in tarTags] ): return False

  return True

def Data16_criteria(tags):
  """
  return true if the tag satisfy Data16 athena 20.7 reprocessing tags as in
  https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/SUSYxAODDerivationsr20#Physics_containers_for_data_peri
  """

  #This DAOD tag is recommended for athena20.7
  #tarTags = ["p2667", "p2689"]
  tarTags = ["p2880"]
  #tarTags = ["p2950"]
  if not any( [ t in tags for t in tarTags] ): return False

  return True


def MC15C_criteria(tags):
  """
  return true if the tag satisfy MC15C recommended tags as in
  https://twiki.cern.ch/twiki/bin/view/AtlasProtected/AtlasProductionGroupMC15c

  and DAOD related pTags as in
  https://twiki.cern.ch/twiki/bin/view/AtlasProtected/SUSYxAODDerivationsr20#Recommended_p_tags
  https://twiki.cern.ch/twiki/bin/view/AtlasProtected/DerivationProductionTeam
  https://twiki.cern.ch/twiki/bin/view/AtlasProtected/CentralMC15ProductionList
  """

  #One of the following must be present
  tarTags = ["r7725","r7772","a818","a821","a766"]
  #tarTags = ["a818","a821"] #tags for AFII fast sim?
  if not any( [ t in tags for t in tarTags] ): return False

  #This merge tag is needed
  if "r7676" not in tags: return False

  #This DAOD tag is recommended for athena20.7
  #tarTags = ["p2666","p2688"]
  #tarTags = ["p2879"]
  tarTags = ["p2949"]
  if not any( [ t in tags for t in tarTags] ): return False

  return True

#==============================================================
#user specified options
#==============================================================
tagCriteria = MC15C_criteria
#tagCriteria = Data15_criteria
#tagCriteria = Data16_criteria
inFile = open(sys.argv[1], 'rU')

#==============================================================
#"Main" program
#==============================================================
client = pyAMI.client.Client('atlas')
AtlasAPI.init()

#outputPath = "MCBG_sample_list.txt"
#outputPath = "Data_sample_list.txt"
# outputPath = "MCSig_sample_list.txt"
outputPath="MCFakes_sample_list.txt"
outputFile = open(outputPath,"w")

reader = csv.DictReader(inFile)
for aLine in reader:
  if aLine["#DSID"][0]=="#": continue

  #read the line
  runID   = int(aLine["#DSID"])
  phyProc   = aLine["Name"]
  generator = aLine["Generator"]
  customPTag = aLine["p2949 available?"]

  #find dataset with tags matching selected criteria
  candidates = AtlasAPI.list_datasets(client, dataset_number=str(runID) , type="DAOD_SUSY2", fields="events")

  matched = False
  for aCandidate in candidates:
    name = aCandidate['ldn']
    tags = name.split('.')[-1]
    if tagCriteria(tags):
      if matched: print "#MultiMatch ", 
      print name, aCandidate['events']
      outputFile.write(name)
      outputFile.write("\n")
      matched = True
  if not matched:
      print "#NoMatch ", runID, phyProc
outputFile.close()
