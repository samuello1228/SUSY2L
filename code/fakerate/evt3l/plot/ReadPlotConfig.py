##from ReadMCInfo import MC_info

from glob import glob
import pdb
import os
thisdir=os.path.dirname(__file__)
plot_config_file=(os.path.join(thisdir,"plot_config.txt"))

plot_config_info={}

config=open(plot_config_file).readlines()

for i in config[1:]:
    if i.strip().startswith("#") or i.strip()=="": continue
    process, name, stack, color, MCIDs = i.split()
    
    for MCID in MCIDs.split(","): 
        plot_config_info[MCID]={"MCID":MCID,
                                "color":int(color),
                                "process":process,
                                "name":name,
                                "stack":stack}
	

# print out all information
#for config_i in plot_config_info:
#   print config_i, plot_config_info[config_i]["scale"]
