from glob import glob
import pdb
import os
thisdir=os.path.dirname(__file__)
output_file=(os.path.join(thisdir,"list_output.txt"))

output_info={}

output=open(output_file).readlines()

for i in output[0:]:
    if i.strip().startswith("#") or i.strip()=="": continue
    ID_position = i.index("skim_")
    #Dot_position = i.index(".root")
    MCID=i[ID_position+5:ID_position+11]
    output_info[i]={"MCID":MCID}
#for add in output_info:
#    print add,output_info[add]["MCID"]
