#!/usr/bin/env python
import pyAMI.atlas.api as AtlasAPI
import pyAMI.client

client = pyAMI.client.Client('atlas')
AtlasAPI.init()

def save1():
  ds = AtlasAPI.list_datasets(client, patterns = ['mc15_13TeV.%0p95_2L5%.e5668%'], type = 'EVNT')
  print ds
  with open("SigXsec.txt", 'w') as outFile:
    for d in ds:
      aLine = d['ldn']
      print aLine
      #remove and ignore comments in line
      aLine = aLine.split("#")[0]
      if len(aLine)==0: continue

      #read the line
      name = aLine.rstrip()
      elements = aLine.split(".")
      mass = elements[2].split("_")

#       #find evgen.EVNT sample
#       prov = AtlasAPI.get_dataset_prov(client, name)
#       print prov
#       node = prov['node']
#       evgen = node[10]
#       name = evgen['logicalDatasetName']
#       print name
      name = aLine

      #find cross section and filter efficiency
      info = AtlasAPI.get_dataset_info(client, name)[0]
      print info
      xSec = float(info['crossSection'])
      xSecLow = float(info['crossSection_min'])
      xSecHigh = float(info['crossSection_max'])
      relU = max(abs(xSecHigh-xSec),abs(xSec-xSecLow))/xSec
      unit = info['crossSection_unit']
      eff = info['GenFiltEff_mean']
      print elements[1], xSec, unit, eff

      #output file
#       outStr = elements[1] + " " + elements[2] + " " + mass[4] + " " + mass[5] + " " + str(float(xSec)*1000) + " " + eff
      outStr = elements[1] + " " + elements[2] + " " + str(xSec*1000) + " " + eff + " " + str(relU)
      print outStr, "\n"
      outFile.write(outStr+"\n")
save1()
