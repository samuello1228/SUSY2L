MCSample=[]

#Signal
#input_ami="ami_data/SigSample_ami.txt"
#MCSample.append("../../multiLepSearch/script/MCSig_sample_list.txt")

#MCBG
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
MCSample.append("../../multiLepSearch/script/MCBGVVSherpa_sample_list.txt")
MCSample.append("../../multiLepSearch/script/MCBGVgammaSherpa_sample_list.txt")
MCSample.append("../../multiLepSearch/script/MCBGhiggs_sample_list.txt")

dict_sampleName={}

for sample in MCSample:
  with open(sample) as sampleFiles:
    for aLine in sampleFiles:
      #remove and ignore comments in line
      aLine = aLine.split("#")[0]
      if len(aLine)==0: continue

      #read the line
      name = aLine.rstrip()
      elements = name.split(".")

      #output file
      print elements[1], name
      dict_sampleName[elements[1]] = name

#fill dictionary: raw number for each sample id
dict_AOD={}
dict_SUSY={}

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
    dict_AOD[elements[0]] = elements[1]
    dict_SUSY[elements[0]] = elements[2]

#For Data driven BG
dict_AOD['999998'] = 0
dict_AOD['999999'] = 0
dict_SUSY['999998'] = 0
dict_SUSY['999999'] = 0

#append raw number in output file
#inputFile="SigSample.txt"
#outputFile="SigSample_add.txt"

inputFile="BGSample_XS.txt"
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
