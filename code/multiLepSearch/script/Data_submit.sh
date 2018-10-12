#!/bin/bash
tag=v13.5

#My GRL
#grl=GoodRunsLists/data16_13TeV/20170720/physics_25ns_20.7.xml,GoodRunsLists/data15_13TeV/20160720/physics_25ns_20.7.xml

#My Data PRW
#dataPRW=GoodRunsLists/data16_13TeV/20170720/physics_25ns_20.7.lumicalc.OflLumi-13TeV-009.root,GoodRunsLists/data15_13TeV/20160720/physics_25ns_20.7.lumicalc.OflLumi-13TeV-005.root

#Peter GRL
grl=multiLepSearch/GRL/data15_13TeV.periodAllYear_DetStatus-v79-repro20-02_DQDefects-00-02-02_PHYS_StandardGRL_All_Good_25ns.xml,multiLepSearch/GRL/data16_13TeV.periodAllYear_DetStatus-v88-pro20-21_DQDefects-00-02-04_PHYS_StandardGRL_All_Good_25ns.xml

#Peter Data PRW
dataPRW=multiLepSearch/prw_Data/ilumicalc_histograms_None_276262-284484.root,multiLepSearch/prw_Data/ilumicalc_histograms_None_297730-311481_OflLumi-13TeV-009.root

#Some dummy MC PRW
mcPRW=multiLepSearch/prw_MC/merged_prw_mc15c_Sep19.root

#SUSYTools Config File
SUSYconf=multiLepSearch/sel_conf/WHSS_0910.conf

study=ss
#study=fakes_Peter
#study=fakes

if [ "$study" = "fakes_Peter" ]
then
    dataPRW=multiLepSearch/prw_Data/ilumicalc_histograms_None_297730-311481_OflLumi-13TeV-009.root
fi

file=../multiLepSearch/script/Data_sample_list.txt
k=${tag}.data 
../multiLepSearch/util/run_ss_selection.py --driver grid --inputList ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --grl ${grl} --outputTag ${k} -o ${k} -a 0 --study ${study} --doSys 0  --conf ${SUSYconf} --fast

#file=/afs/cern.ch/work/c/clo/sample/data16_13TeV.00300279.physics_Main.merge.DAOD_SUSY2.f705_m1606_p2950/DAOD_SUSY2.10314997._000013.pool.root.1

#../multiLepSearch/util/run_ss_selection.py --driver direct -f ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --grl ${grl} -o t1 -w -a 0 --study ${study} --doSys 0  --conf ${SUSYconf} # --cutflow | tee cutflow.log
