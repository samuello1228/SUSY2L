#fill dictionary: sample name for each sample id
dict_sampleName={}
with open("SigSample_XS.txt") as originalFile:
  for aLine in originalFile:
    #read the line
    elements = aLine.rstrip().split()
    dict_sampleName[elements[0]] = elements[1]

#print dict_sampleName

#output file
with open("MGPy8EG_A14N23LO_C1N2_Wh_2L7_XX_YY.txt") as SUSYFile:
  with open("SigSample_XS.txt", 'w') as outFile:
    isFirstLine = 1
    for aLine in SUSYFile:
      #remove and ignore comments in line
      aLine = aLine.split("#")
      if len(aLine[0])==0: continue
      
      #read the line
      elements = aLine[0].split()
      masses = aLine[1].rstrip()
      
      if isFirstLine == 1:
        kfactor_125 = float(elements[3])
      else:
        kfactor_127 = float(elements[3])
        kfactor_total = kfactor_125 + kfactor_127
        
        masses = masses.replace("(","")
        masses = masses.replace(")","")
        masses = masses.replace(",","")
        
        #output file
        outStr = elements[0] + " " + dict_sampleName[elements[0]] + masses + " " + elements[2] + " " + str(kfactor_total) + " 1. " + str(float(elements[4])*10) + "E-01"
        print outStr, "\n"
        outFile.write(outStr+"\n")
        
      #inter-change isFirstLine
      if isFirstLine == 1:
        isFirstLine = 0
      else:
        isFirstLine = 1
