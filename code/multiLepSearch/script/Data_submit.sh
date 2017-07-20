#!/bin/bash
tag=v10.4
grl=GoodRunsLists/data16_13TeV/20170215/physics_25ns_20.7.xml,GoodRunsLists/data15_13TeV/20160720/physics_25ns_20.7.xml
dataPRW=GoodRunsLists/data16_13TeV/20170215/physics_25ns_20.7.lumicalc.OflLumi-13TeV-008.root,GoodRunsLists/data15_13TeV/20160720/physics_25ns_20.7.lumicalc.OflLumi-13TeV-005.root
SUSYconf=SUSYTools/SUSYTools_Default.conf

mcPRW=multiLepSearch/prw_MC/merged_prw_mc15c_Jun15.root

file=../multiLepSearch/script/Data_sample_list.txt
k=${tag}.data 
../multiLepSearch/util/run_ss_selection.py --driver grid --inputList ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --grl ${grl} --outputTag ${tag} -o ${k} -a 0 --study ss --doSys 0  --conf ${SUSYconf} --fast

#file=/afs/cern.ch/work/c/clo/sample/data16_13TeV.00300279.physics_Main.merge.DAOD_SUSY2.f705_m1606_p2950/DAOD_SUSY2.10314997._000013.pool.root.1

#../multiLepSearch/util/run_ss_selection.py --driver direct -f ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --grl ${grl} -o t1 -w -a 0 --study ss --doSys 0  --conf ${SUSYconf}
