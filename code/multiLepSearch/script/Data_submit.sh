#!/bin/bash
tag=v8.13
dataPRW=GoodRunsLists/data16_13TeV/20161101/physics_25ns_20.7.lumicalc.OflLumi-13TeV-005.root,GoodRunsLists/data15_13TeV/20160720/physics_25ns_20.7.lumicalc.OflLumi-13TeV-005.root
grl=GoodRunsLists/data16_13TeV/20161101/physics_25ns_20.7.xml,GoodRunsLists/data15_13TeV/20160720/physics_25ns_20.7.xml
#mcPRW=dev/SUSYTools/merged_prw_mc15c_latest.root
mcPRW=dev/SUSYTools/merged_prw_mc15c_Jan31.root
#mcPRW=dev/SUSYTools/merged_prw_mc15c_Feb13.root

file=../multiLepSearch/script/Data_sample_list.txt
k=${tag}.data 
../multiLepSearch/util/run_ss_selection.py --driver grid --inputList ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --grl ${grl} --outputTag ${tag} -o ${k} -w -a 0 --study ss --doSys 0
