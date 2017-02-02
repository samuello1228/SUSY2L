#!/bin/bash
tag=v8.8
dataPRW=multiLepSearch/prw_Data/ilumicalc_histograms_None_297730-311481_OflLumi-13TeV-005.root
grl=multiLepSearch/GRL/physics_25ns_20.7.xml
mcPRW=multiLepSearch/prw_MC/merged_prw_mc15c_July27_afterFix.root
file=../multiLepSearch/script/Data_sample_list.txt
k=${tag}.data 
../multiLepSearch/util/run_ss_selection.py --driver grid --inputList ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --grl ${grl} --outputTag ${tag} -o ${k} -w -a 0 --study ss --doSys 0