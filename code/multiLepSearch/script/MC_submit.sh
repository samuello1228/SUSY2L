#!/bin/bash
tag=v9.3
dataPRW=GoodRunsLists/data16_13TeV/20161101/physics_25ns_20.7.lumicalc.OflLumi-13TeV-005.root,GoodRunsLists/data15_13TeV/20160720/physics_25ns_20.7.lumicalc.OflLumi-13TeV-005.root

#For background MC
#mcPRW=dev/SUSYTools/merged_prw_mc15c_latest.root
# mcPRW=dev/SUSYTools/merged_prw_mc15c_Jan31.root
mcPRW=dev/SUSYTools/merged_prw_mc15c_Feb13.root

# file=../multiLepSearch/script/MCBG_sample_list.txt
#file=../multiLepSearch/script/MCBGZjetsSherpa_sample_list.txt
#file=../multiLepSearch/script/MCBGVVSherpa_sample_list.txt
#file=../multiLepSearch/script/MCBGVgammaSherpa_sample_list.txt
# file=../multiLepSearch/script/MCBGDYSherpa_sample_list.txt
file=../multiLepSearch/script/MCBGZeeSherpa.txt
#file=../multiLepSearch/script/MCBG_llll.txt
#file=../multiLepSearch/script/MCZeePowheg.txt
# file=./testZee.txt

# k=${tag}.MCBG
k=${tag}.MCZee

#For signal MC
#mcPRW=multiLepSearch/prw_MC/merged_prw_mc15c_Slep0d95.root,dev/PileupReweighting/mc15c_v2_defaults.NotRecommended.prw.root
#file=../multiLepSearch/script/MCSig_sample_list.txt
#k=${tag}.MCSig

driver=grid
# driver=direct

../multiLepSearch/util/run_ss_selection.py --driver ${driver} --inputList ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o ${k} -w -a 1 --study ss --mcMatch dR --doSys 0 --ChargeID 1

