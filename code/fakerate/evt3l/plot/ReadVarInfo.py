# Read setting for each physics variables


from glob import glob
import pdb
import os
thisdir=os.path.dirname(__file__)
var_info_file=(os.path.join(thisdir,"var_list.txt"))

var_info={}

var=open(var_info_file).readlines()

for i in var:
    if i.strip().startswith("#") or i.strip()=="": continue
    branch,name,unit,nbin,min,max,log = i.split()
    
    var_info[branch]={"name":name,
                      "unit":unit,
                      "nbin":int(nbin),
                      "min":float(min),
                      "max":float(max),
                      "log":int(log)}

# print out all information
#for sample in data_info:
#    print data_info[sample]["Ldata"]
