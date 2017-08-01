sigFilesTxt="SigSample_ami.txt"
SUSYFilesTxt="SigXS_SUSYTools.txt"
print "reading", sigFilesTxt

dict1={}

with open(SUSYFilesTxt) as SUSYFiles:
  for aLine2 in SUSYFiles:
    name = aLine2.rstrip()
    
    elements = name.split()
    mass1 = elements[7]
    mass1 = mass1.replace("(","")
    mass1 = mass1.replace(",","")

    dict1[(mass1,elements[1])] = elements[2]

dict1[("525.0", "125")] = "0.0266028"
dict1[("550.0", "125")] = "0.021644"
dict1[("575.0", "125")] = "0.0177159"

dict1[("525.0", "127")] = "0.0106671"
dict1[("550.0", "127")] = "0.00852674"
dict1[("575.0", "127")] = "0.00686761"
print dict1

with open(sigFilesTxt) as sigFiles:
  with open("SigSample.txt", 'w') as outFile:
    for aLine1 in sigFiles:
      #remove and ignore comments in line
      aLine1 = aLine1.split("#")[0]
      if len(aLine1)==0: continue

      #read the line
      name = aLine1.rstrip()
      elements = name.split(" ")
      mass1 = elements[2];
      XS = float(dict1[(mass1,"125")]) + float(dict1[(mass1,"127")])

      #output file
      outStr = elements[0] + " " + elements[1] + " " + elements[2] + " " + elements[3] + " " + str(XS) + " 1. " + elements[5]
      print outStr, "\n"
      outFile.write(outStr+"\n")
