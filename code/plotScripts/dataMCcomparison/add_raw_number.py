#source /afs/cern.ch/atlas/software/tools/pyAMI/setup.sh
import pyAMI.atlas.api as AtlasAPI
import pyAMI.client

client = pyAMI.client.Client('atlas')
AtlasAPI.init()

MCSample=[]

#Signal
#MCSample.append("../../multiLepSearch/script/MCSig_sample_list.txt")

#MCBG
MCSample.append("../../multiLepSearch/script/MCBGZjetsSherpa_sample_list.txt")
MCSample.append("../../multiLepSearch/script/MCBGWjetsSherpa_sample_list.txt")
MCSample.append("../../multiLepSearch/script/MCBGDYSherpa_sample_list.txt")
MCSample.append("../../multiLepSearch/script/MCBGmultitop_sample_list.txt")
MCSample.append("../../multiLepSearch/script/MCBGmultitop_fast_sample_list.txt")
MCSample.append("../../multiLepSearch/script/MCBGVVVSherpa_sample_list.txt")
MCSample.append("../../multiLepSearch/script/MCBGsingletop_sample_list.txt")
MCSample.append("../../multiLepSearch/script/MCBGttbar_sample_list.txt")
MCSample.append("../../multiLepSearch/script/MCBGttV_sample_list.txt")
MCSample.append("../../multiLepSearch/script/MCBGVVSherpa_sample_list.txt")
MCSample.append("../../multiLepSearch/script/MCBGVgammaSherpa_sample_list.txt")
MCSample.append("../../multiLepSearch/script/MCBGhiggs_sample_list.txt")

#fill dictionary: raw number for each sample id
dict_AOD={}
dict_SUSY={}

for sample in MCSample:
  with open(sample) as sampleFiles:
    for aLine in sampleFiles:
      #remove and ignore comments in line
      aLine = aLine.split("#")[0]
      if len(aLine)==0: continue

      #read the line
      name = aLine.rstrip()
      elements = name.split(".")

      #find number of event of SUSY2
      info = AtlasAPI.get_dataset_info(client, name)[0]
      SUSYnumber = info['totalEvents']

      #find AOD sample
      prov = AtlasAPI.get_dataset_prov(client, name)
      node = prov['node']

      AODname = node[1]['logicalDatasetName']
      print AODname

      #find number of event of AOD
      info = AtlasAPI.get_dataset_info(client, AODname)[0]
      AODnumber = info['totalEvents']

      #find cross section and filter efficiency
      #xSec = info['crossSection']
      #unit = info['crossSection_unit']
      #eff = info['GenFiltEff_mean']
      #print elements[1], xSec, unit, eff

      #output file
      print elements[1], AODnumber, SUSYnumber
      dict_AOD[elements[1]] = AODnumber
      dict_SUSY[elements[1]] = SUSYnumber

#For Data driven BG
dict_AOD['999998'] = 0
dict_AOD['999999'] = 0
dict_SUSY['999998'] = 0
dict_SUSY['999999'] = 0


#append raw number in output file
#inputFile="SigSample.txt"
#outputFile="SigSample_add.txt"

inputFile="BGSample.txt"
outputFile="BGSample_add.txt"
print "reading", inputFile

with open(inputFile) as File:
  with open(outputFile, 'w') as outFile:
    for aLine in File:
      #remove and ignore comments in line
      aLine = aLine.split("#")[0]
      if len(aLine)==0: continue

      #find sample id
      elements = aLine.rstrip().split()
      sampleID = elements[0]

      #output file
      outStr = aLine.rstrip() + " " + str(dict_AOD[sampleID]) + " " + str(dict_SUSY[sampleID])
      print outStr
      outFile.write(outStr+"\n")
