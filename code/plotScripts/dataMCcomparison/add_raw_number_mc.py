MCSample=[]

#Signal
#inputFile="SigSample_XS.txt"
#outputFile="SigSample_add.txt"
#input_ami="ami_data/SigSample_ami.txt"
#MCSample.append("../../multiLepSearch/script/MCSig_sample_list.txt")

#MCBG
inputFile="BGSample_XS.txt"
outputFile="BGSample_add.txt"
input_ami="ami_data/BGSample_ami.txt"
MCSample.append("../../multiLepSearch/script/MCBGZjetsSherpa_sample_list.txt")
MCSample.append("../../multiLepSearch/script/MCBGWjetsSherpa_sample_list.txt")
MCSample.append("../../multiLepSearch/script/MCBGDYSherpa_sample_list.txt")
MCSample.append("../../multiLepSearch/script/MCBGmultitop_sample_list.txt")
MCSample.append("../../multiLepSearch/script/MCBGmultitop_fast_sample_list.txt")
MCSample.append("../../multiLepSearch/script/MCBGVVVSherpa_sample_list.txt")
MCSample.append("../../multiLepSearch/script/MCBGsingletop_sample_list.txt")
MCSample.append("../../multiLepSearch/script/MCBGttbar_sample_list.txt")
MCSample.append("../../multiLepSearch/script/MCBGttV_sample_list.txt")
MCSample.append("../../multiLepSearch/script/MCBGVVSherpa_CT10_sample_list.txt")
MCSample.append("../../multiLepSearch/script/MCBGVVSherpa_221_sample_list.txt")
MCSample.append("../../multiLepSearch/script/MCBGVgammaSherpa_sample_list.txt")
MCSample.append("../../multiLepSearch/script/MCBGhiggs_Pythia8EvtGen_sample_list.txt")
MCSample.append("../../multiLepSearch/script/MCBGhiggs_HerwigppEvtGen_sample_list.txt")
MCSample.append("../../multiLepSearch/script/MCBGRare_sample_list.txt")

#fill dictionary: sample name for each sample id
dict_sampleName={}
#fill dictionary: rucio raw number for each sample id
dict_SUSY_rucio={}

for sample in MCSample:
  with open(sample) as sampleFiles:
    for aLine in sampleFiles:
      #remove and ignore comments in line
      aLine = aLine.split("#")[0]
      if len(aLine)==0: continue

      #read the line
      name = aLine.rstrip()
      elements = name.split(".")

      #read rucio data
      rucioFilePath = "rucio_data/" + name + ".txt"
      with open(rucioFilePath) as rucioFile:
        for aLine2 in rucioFile:
            name2 = aLine2.rstrip()

            if "Total events" in name2:
              number_str=""
              for digit in name2:
                if digit.isdigit():
                  number_str += digit
              print name2, number_str
              dict_SUSY_rucio[elements[1]] = number_str

            #if "Total size" in name2:
              #print name2

      #fill dictionary: sample name
      print elements[1], name
      dict_sampleName[elements[1]] = name

#For Data driven BG
dict_SUSY_rucio['999998'] = 0
dict_SUSY_rucio['999999'] = 0

#fill dictionary: ami raw number for each sample id
dict_AOD_ami={}
dict_SUSY_ami={}

with open(input_ami) as sampleFiles:
  for aLine in sampleFiles:
    #remove and ignore comments in line
    aLine = aLine.split("#")[0]
    if len(aLine)==0: continue

    #read the line
    name = aLine.rstrip()
    elements = name.split(" ")

    #fill dictionary
    print elements[0], elements[1], elements[2]
    dict_AOD_ami[elements[0]] = elements[1]
    dict_SUSY_ami[elements[0]] = elements[2]

#For Data driven BG
dict_AOD_ami['999998'] = 0
dict_AOD_ami['999999'] = 0
dict_SUSY_ami['999998'] = 0
dict_SUSY_ami['999999'] = 0

#append raw number in output file
with open(inputFile) as File:
  with open(outputFile, 'w') as outFile:
    for aLine in File:
      #remove and ignore comments in line
      aLine = aLine.split("#")[0]
      if len(aLine)==0: continue

      #find sample id
      elements = aLine.rstrip().split()
      sampleID = elements[0]

      #compare ami and rucio data
      if dict_SUSY_ami[sampleID] != dict_SUSY_rucio[sampleID]:
        print "ami and rucio are different", dict_SUSY_ami[sampleID], dict_SUSY_rucio[sampleID]

      #output file
      outStr = aLine.rstrip() + " " + str(dict_AOD_ami[sampleID]) + " " + str(dict_SUSY_ami[sampleID])
      print outStr
      outFile.write(outStr+"\n")
