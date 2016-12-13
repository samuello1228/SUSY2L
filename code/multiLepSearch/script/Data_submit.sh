#!/bin/bash
tag=v8.0
dataPRW=multiLepSearch/prw_Data/ilumicalc_histograms_None_297730-303560_v80_Nicky.root
grl=multiLepSearch/GRL/data16_13TeV.periodAllYear_DetStatus-v80-pro20-08_DQDefects-00-02-02_PHYS_StandardGRL_All_Good_25ns.xml
mcPRW=multiLepSearch/prw_MC/merged_prw_mc15c_July27_afterFix.root
file=../multiLepSearch/script/Data_sample_list.txt
k=${tag}.data 
../multiLepSearch/util/run_ss_selection.py --driver grid --inputList ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --grl ${grl} --outputTag ${tag} -o ${k} -w -a 0 --study ss --doSys 0
