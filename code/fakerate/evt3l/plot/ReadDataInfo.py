##from ReadDataLuminosity import data_list.txt

from glob import glob
import pdb
import os
thisdir=os.path.dirname(__file__)
data_info_file=(os.path.join(thisdir,"data_list.txt"))

data_info={}

data=open(data_info_file).readlines()

for i in data:
    if i.strip().startswith("#") or i.strip()=="": continue
    dataname,Ldata = i.split()
    
    data_info[dataname]={"name":dataname,
                         "Ldata":float(Ldata)}

# print out all information
#for sample in data_info:
#    print data_info[sample]["Ldata"]
