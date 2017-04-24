import pyAMI.atlas.api as AtlasAPI
import pyAMI.client

client = pyAMI.client.Client('atlas')
AtlasAPI.init()

sigFilesTxt="MCSig_sample_list.txt"
print "reading", sigFilesTxt

with open(sigFilesTxt) as sigFiles:
  with open("SigSample.txt", 'w') as outFile:
    for aLine in sigFiles:
      #remove and ignore comments in line
      aLine = aLine.split("#")[0]
      if len(aLine)==0: continue

      #read the line
      name = aLine.rstrip()
      elements = aLine.split(".")
      mass = elements[2].split("_")

      #find evgen.EVNT sample
      prov = AtlasAPI.get_dataset_prov(client, name)
      node = prov['node']
      evgen = node[10]
      name = evgen['logicalDatasetName']
      print name

      #find cross section and filter efficiency
      info = AtlasAPI.get_dataset_info(client, name)[0]
      xSec = info['crossSection']
      unit = info['crossSection_unit']
      eff = info['GenFiltEff_mean']
      print elements[1], xSec, unit, eff

      #output file
      outStr = elements[1] + " " + elements[2] + " " + mass[4] + " " + mass[5] + " " + str(float(xSec)*1000) + " " + eff
      print outStr, "\n"
      outFile.write(outStr+"\n")
