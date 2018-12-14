inputFile="DataSample_lumi.txt"
outputFile="DataSample_add.txt"
input_ami="ami_data/DataSample_ami.txt"

#fill dictionary: sample name for each period
dict_sampleName={}
#fill dictionary: rucio raw number for each period
dict_SUSY_rucio={}

with open("../../multiLepSearch/script/Data_sample_list.txt") as sampleFiles:
  for aLine in sampleFiles:
    #remove and ignore comments in line
    aLine = aLine.split("#")[0]
    if len(aLine)==0: continue

    #read the line
    name = aLine.rstrip()
    elements = name.split(".")

    #read rucio data
    rucioFilePath = "rucio_data/" + name + ".txt"
    keyName = elements[0] + "." + elements[1]
    with open(rucioFilePath) as rucioFile:
      for aLine2 in rucioFile:
          name2 = aLine2.rstrip()

          if "Total events" in name2:
            number_str=""
            for digit in name2:
              if digit.isdigit():
                number_str += digit
            print name2, number_str
            dict_SUSY_rucio[keyName] = number_str

          #if "Total size" in name2:
            #print name2

    #fill dictionary: sample name
    print keyName, name
    dict_sampleName[keyName] = name

#fill dictionary: ami raw number for each period
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

#append raw number in output file
with open(inputFile) as File:
  with open(outputFile, 'w') as outFile:
    for aLine in File:
      #remove and ignore comments in line
      aLine = aLine.split("#")[0]
      if len(aLine)==0: continue

      #find year
      elements = aLine.rstrip().split()
      year = elements[0]

      #Add events
      nAOD = 0
      nSUSY = 0
      for period in dict_SUSY_ami:
        if year in period:
          #compare ami and rucio data
          if dict_SUSY_ami[period] != dict_SUSY_rucio[period]:
            print "ami and rucio are different", dict_SUSY_ami[period], dict_SUSY_rucio[period]

          #add
          nAOD += int(dict_AOD_ami[period])
          nSUSY += int(dict_SUSY_ami[period])

      #output file
      outStr = aLine.rstrip() + " " + str(nAOD) + " " + str(nSUSY)
      print outStr
      outFile.write(outStr+"\n")
