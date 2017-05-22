#!/bin/bash
tag=v9.3.1
dataPRW=GoodRunsLists/data16_13TeV/20170215/physics_25ns_20.7.lumicalc.OflLumi-13TeV-008.root,GoodRunsLists/data15_13TeV/20160720/physics_25ns_20.7.lumicalc.OflLumi-13TeV-005.root
#mcPRW=dev/SUSYTools/merged_prw_mc15c_Jan31.root,multiLepSearch/prw_MC/merged_prw_mc15c_Slep0d95.root,dev/PileupReweighting/mc15c_v2_defaults.NotRecommended.prw.root
mcPRW=dev/PileupReweighting/mc15c_v2_defaults.NotRecommended.prw.root
SUSYconf=multiLepSearch/sel_conf/SUSYTools_multilepAnaMoriond.conf

#For background MC
#mcPRW=dev/SUSYTools/merged_prw_mc15c_latest.root

#file=../multiLepSearch/script/MCBG_sample_list.txt,../multiLepSearch/script/MCBGZjetsSherpa_sample_list.txt,../multiLepSearch/script/MCBGVVSherpa_sample_list.txt,../multiLepSearch/script/MCBGVgammaSherpa_sample_list.txt,../multiLepSearch/script/MCBGDYSherpa_sample_list.txt
#file=../multiLepSearch/script/MCBG_sample_list.txt
#file=../multiLepSearch/script/MCBGZjetsSherpa_sample_list.txt
file=../multiLepSearch/script/MCBGVVSherpa_sample_list.txt
#file=../multiLepSearch/script/MCBGVgammaSherpa_sample_list.txt
#file=../multiLepSearch/script/MCBGDYSherpa_sample_list.txt

#file=../multiLepSearch/script/MCBGZeeSherpaSelected.txt
#file=../multiLepSearch/script/MCBG_llll.txt

k=${tag}.MCBG

#For signal MC
#file=../multiLepSearch/script/MCSig_sample_list.txt
#file=../multiLepSearch/script/MCSig_sample_list_new.txt
#k=${tag}.MCSig

../multiLepSearch/util/run_ss_selection.py --driver grid --inputList ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o ${k} -w -a 1 --study ss --mcMatch dR --doSys 0 --ChargeID 1 --conf ${SUSYconf}
