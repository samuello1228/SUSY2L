#!/bin/bash
tag=v12.1
grl=GoodRunsLists/data16_13TeV/20170720/physics_25ns_20.7.xml,GoodRunsLists/data15_13TeV/20160720/physics_25ns_20.7.xml
dataPRW=GoodRunsLists/data16_13TeV/20170720/physics_25ns_20.7.lumicalc.OflLumi-13TeV-009.root,GoodRunsLists/data15_13TeV/20160720/physics_25ns_20.7.lumicalc.OflLumi-13TeV-005.root
SUSYconf=multiLepSearch/sel_conf/WHSS_0910.conf

#study=ss
study=fakes

mcPRW=multiLepSearch/prw_MC/merged_prw_mc15c_Sep19.root

file=../multiLepSearch/script/Data_sample_list.txt
k=${tag}.data 
../multiLepSearch/util/run_ss_selection.py --driver grid --inputList ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --grl ${grl} --outputTag ${k} -o ${k} -a 0 --study ${study} --doSys 0  --conf ${SUSYconf} --fast

#file=/afs/cern.ch/work/c/clo/sample/data16_13TeV.00300279.physics_Main.merge.DAOD_SUSY2.f705_m1606_p2950/DAOD_SUSY2.10314997._000013.pool.root.1

#../multiLepSearch/util/run_ss_selection.py --driver direct -f ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --grl ${grl} -o t1 -w -a 0 --study ${study} --doSys 0  --conf ${SUSYconf} # --cutflow | tee cutflow.log
