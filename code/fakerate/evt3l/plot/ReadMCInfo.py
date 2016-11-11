##from ReadMCInfo import MC_info

from glob import glob
import pdb
import os
thisdir=os.path.dirname(__file__)
MC_info_file=(os.path.join(thisdir,"susy_crosssections_13TeV.txt"))

MC_info={}

data=open(MC_info_file).readlines()
header=data[0].split()

for i in data[1:]:
    if i.strip().startswith("#") or i.strip()=="": continue
    MCID, name, xsec, kfac, eff, relunc = i.split()
    
    MC_info[MCID]={"name":name,
                   "xsec":float(xsec),
                   "eff":float(eff),
                   "kfac":float(kfac),
                   "relunc":float(relunc)}

# print out all information
#for sample in MC_info:
#   print sample, MC_info[sample]["k_factor"], MC_info[sample]["FilterEff"], MC_info[sample]["Xsec"], MC_info[sample]["Sample"]
