#Usage:
#before runing this script, you have to setup pyAMi env by:
#source /afs/cern.ch/atlas/software/tools/pyAMI/setup.sh
#
#then:
#python dataseID_to_name.py inputIDs.txt

import sys

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
  tarTags = ["p2667", "p2689"]
  if not any( [ t in tags for t in tarTags] ): return False

  return True


def MC15C_criteria(tags):
  """
  return true if the tag satisfy MC15C recommended tags as in
  https://twiki.cern.ch/twiki/bin/view/AtlasProtected/AtlasProductionGroupMC15c

  and DAOD related pTags as in
  https://twiki.cern.ch/twiki/bin/view/AtlasProtected/SUSYxAODDerivationsr20#Recommended_p_tags
  """

  #One of the following must be present
  tarTags = ["r7725","r7772","a818","a821"]
  #tarTags = ["a818","a821"] #tags for AFII fast sim?
  if not any( [ t in tags for t in tarTags] ): return False

  #This merge tag is needed
  if "r7676" not in tags: return False

  #This DAOD tag is recommended for athena20.7
  tarTags = ["p2666","p2688"]
  if not any( [ t in tags for t in tarTags] ): return False

  return True

#==============================================================
#user specified options
#==============================================================
tagCriteria = MC15C_criteria
#tagCriteria = Data15_criteria
#tagCriteria = Data16_criteria
inFile = open(sys.argv[1])

#==============================================================
#"Main" program
#==============================================================
client = pyAMI.client.Client('atlas')
AtlasAPI.init()

outputPath = "MCBG_sample_list.txt"
outputFile = open(outputPath,"w")

for aLine in inFile:
  #remove and ignore comments in line
  aLine = aLine.split("#")[0] # content after # are considered commnets
  if len(aLine)==0: continue

  #read the line
  elements = aLine.split()
  if len(elements)!=4: 
    print "#cannont understand line:"
    print "#" + aLine
    continue

  phyType   = elements[0]
  phyProc   = elements[1]
  generator = elements[2]
  runIDs    = elements[3].split(",")

  #interpret and expand runIDs range
  tmp = []
  for aStr in runIDs:
    if '-' in aStr:
      (startID, endID) = map(int, aStr.split('-') )
      tmp += [i for i in range(startID, endID+1)]
    else:
      tmp.append( int(aStr) )
  runIDs = tmp

  #find dataset with tags matching selected criteria
  for aRunID in runIDs:
    if phyType.lower() == "data": 
      candidates = AtlasAPI.list_datasets(client, run_number=str(aRunID) , stream="physics_Main", type="DAOD_SUSY2", fields="events")
    else:
      candidates = AtlasAPI.list_datasets(client, dataset_number=str(aRunID) , type="DAOD_SUSY2", fields="events")

    matched = False
    for aCandidate in candidates:
      name = aCandidate['ldn']
      tags = name.split('.')[-1]
      if tagCriteria(tags):
        if matched: print "#MultiMatch ", 
        print name, phyType, phyProc, aCandidate['events']
        outputFile.write(name)
        outputFile.write("\n")
        matched = True
    if not matched:
        print "#NoMatch ", aRunID, phyType, phyProc

outputFile.close()
