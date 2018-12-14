inputFile="../BGSample_XS.txt"

input_SUSY="susy_crosssections_13TeV_00-08-60.txt"
outputFile="BGSample_XS_00-08-60.txt"

#input_SUSY="susy_crosssections_13TeV_00-08-71.txt"
#outputFile="BGSample_XS_00-08-71.txt"

#fill dictionary: cross section from the SUSYTools file
dict_SUSY_XS={}

with open(input_SUSY) as SUSYFile:
  for aLine in SUSYFile:
    #remove and ignore comments in line
    aLine = aLine.split("#")[0]
    if len(aLine)==0: continue

    #read the line
    name = aLine.rstrip()

    #remove blank line
    if not name:
      continue
    
    #split by whitespace
    elements = name.split()
    
    #ignore other lines
    if not elements[0].isdigit():
      continue
    
    dict_SUSY_XS[elements[0]] = aLine

#print dict_SUSY_XS

with open(inputFile) as InputFile:
  with open(outputFile, 'w') as OutputFile:
    for aLine in InputFile:
      #read the line
      name = aLine.rstrip()
      elements = name.split()

      if elements[0] in dict_SUSY_XS:
        newLine = dict_SUSY_XS[elements[0]]
        newLine = newLine.replace("MadGraphPythia8EvtGen_ttbarWll","MadGraphPythia8EvtGen_A14NNPDF23LO_ttbarWll")
        newLine = newLine.replace("MadGraph_3topSM","MadGraphPythia8EvtGen_A14NNPDF23_3top_SM")
        newLine = newLine.replace("1 .","1.")
        
        if newLine[-1] != "\n":
          newLine = newLine + "\n"

        print newLine
        OutputFile.write(newLine)
      else:
        print name
        OutputFile.write(aLine)
