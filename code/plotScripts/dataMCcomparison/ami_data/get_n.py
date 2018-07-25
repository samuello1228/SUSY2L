#source /afs/cern.ch/atlas/software/tools/pyAMI/setup.sh
import pyAMI.atlas.api as AtlasAPI
import pyAMI.client

client = pyAMI.client.Client('atlas')
AtlasAPI.init()

MCSample=[]

#Signal
MCSample.append("../../../multiLepSearch/script/MCSig_sample_list.txt")
outputFile="SigSample_ami.txt"

#MCBG
#MCSample.append("../../../multiLepSearch/script/MCBGZjetsSherpa_sample_list.txt")
#MCSample.append("../../../multiLepSearch/script/MCBGWjetsSherpa_sample_list.txt")
#MCSample.append("../../../multiLepSearch/script/MCBGDYSherpa_sample_list.txt")
#MCSample.append("../../../multiLepSearch/script/MCBGmultitop_sample_list.txt")
#MCSample.append("../../../multiLepSearch/script/MCBGmultitop_fast_sample_list.txt")
#MCSample.append("../../../multiLepSearch/script/MCBGVVVSherpa_sample_list.txt")
#MCSample.append("../../../multiLepSearch/script/MCBGsingletop_sample_list.txt")
#MCSample.append("../../../multiLepSearch/script/MCBGttbar_sample_list.txt")
#MCSample.append("../../../multiLepSearch/script/MCBGttV_sample_list.txt")
#MCSample.append("../../../multiLepSearch/script/MCBGVVSherpa_sample_list.txt")
#MCSample.append("../../../multiLepSearch/script/MCBGVgammaSherpa_sample_list.txt")
#MCSample.append("../../../multiLepSearch/script/MCBGhiggs_sample_list.txt")
#outputFile="BGSample_ami.txt"


with open(outputFile, 'w') as outFile:
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
        outStr = elements[1] + " " + str(AODnumber) + " " + str(SUSYnumber)
        print outStr
        outFile.write(outStr+"\n")

