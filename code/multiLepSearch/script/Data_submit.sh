#!/bin/bash
tag=v11
grl=GoodRunsLists/data16_13TeV/20170215/physics_25ns_20.7.xml,GoodRunsLists/data15_13TeV/20160720/physics_25ns_20.7.xml
dataPRW=GoodRunsLists/data16_13TeV/20170215/physics_25ns_20.7.lumicalc.OflLumi-13TeV-008.root,GoodRunsLists/data15_13TeV/20160720/physics_25ns_20.7.lumicalc.OflLumi-13TeV-005.root
SUSYconf=multiLepSearch/sel_conf/SUSYTools_multilepAnaMoriond.conf

#mcPRW=dev/SUSYTools/merged_prw_mc15c_latest.root
mcPRW=dev/SUSYTools/merged_prw_mc15c_Jan31.root
#mcPRW=dev/SUSYTools/merged_prw_mc15c_Feb13.root

file=../multiLepSearch/script/Data_sample_list.txt
k=${tag}.data 
../multiLepSearch/util/run_ss_selection.py --driver grid --inputList ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --grl ${grl} --outputTag ${tag} -o ${k} -w -a 0 --study fakes --doSys 0  --conf ${SUSYconf}
