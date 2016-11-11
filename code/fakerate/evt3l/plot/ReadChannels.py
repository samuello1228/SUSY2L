##from ReadMCInfo import MC_info

from glob import glob
import pdb
import os
thisdir=os.path.dirname(__file__)
channel_file=(os.path.join(thisdir,"channel_names.txt"))

channel_info=[]

config=open(channel_file).readlines()

for i in config:
    if i.strip().startswith("#") or i.strip()=="": continue
    temp = i.split()
    if(temp[0]=="channel:"): 
        channel_info=channel_info+temp[1:]
	

# print out all information
#for config_i in plot_config_info:
#   print config_i, plot_config_info[config_i]["scale"]
