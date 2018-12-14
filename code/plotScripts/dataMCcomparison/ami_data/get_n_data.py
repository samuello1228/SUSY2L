#source /afs/cern.ch/atlas/software/tools/pyAMI/setup.sh
import pyAMI.atlas.api as AtlasAPI
import pyAMI.client

client = pyAMI.client.Client('atlas')
AtlasAPI.init()

outputFile="DataSample_ami.txt"

with open(outputFile, 'w') as outFile:
  with open("../../../multiLepSearch/script/Data_sample_list.txt") as sampleFiles:
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

      dataset = info['contained_dataset']
      dataset = dataset.split(';')
      del dataset[-1]

      AODnumber = 0
      for sample in dataset:

        #find AOD sample
        prov = AtlasAPI.get_dataset_prov(client, sample)
        node = prov['node']
  
        AODname = node[1]['logicalDatasetName']
        print AODname
  
        #find number of event of AOD
        info = AtlasAPI.get_dataset_info(client, AODname)[0]
        AODnumber += int(info['totalEvents'])
  
        #runNumber = info['runNumber']

      #output file
      outStr = elements[0] + "." + elements[1] + " " + str(AODnumber) + " " + str(SUSYnumber)
      print outStr
      outFile.write(outStr+"\n")
