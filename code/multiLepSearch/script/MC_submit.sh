#!/bin/bash
tag=v11.3.MCTC
dataPRW=GoodRunsLists/data16_13TeV/20170215/physics_25ns_20.7.lumicalc.OflLumi-13TeV-008.root,GoodRunsLists/data15_13TeV/20160720/physics_25ns_20.7.lumicalc.OflLumi-13TeV-005.root
mcPRW=multiLepSearch/prw_MC/merged_prw_mc15c_Jun15.root,multiLepSearch/prw_MC/merged_prw_mc15c_Slep0d95.root,dev/PileupReweighting/mc15c_v2_defaults.NotRecommended.prw.root
# SUSYconf=multiLepSearch/sel_conf/SUSYTools_multilepAnaMoriond.conf
SUSYconf=multiLepSearch/sel_conf/SUSYTools_Wh.conf

#For background MC
#mcPRW=dev/SUSYTools/merged_prw_mc15c_latest.root

# file=../multiLepSearch/script/MCBG_sample_list.txt,../multiLepSearch/script/MCBGZjetsSherpa_sample_list.txt,../multiLepSearch/script/MCBGVVSherpa_sample_list.txt,../multiLepSearch/script/MCBGVgammaSherpa_sample_list.txt,../multiLepSearch/script/MCBGDYSherpa_sample_list.txt
#file=../multiLepSearch/script/MCBG_sample_list.txt
#file=../multiLepSearch/script/MCBGZjetsSherpa_sample_list.txt
#file=../multiLepSearch/script/MCBGVVSherpa_sample_list.txt
#file=../multiLepSearch/script/MCBGVgammaSherpa_sample_list.txt
#file=../multiLepSearch/script/MCBGDYSherpa_sample_list.txt
file=../multiLepSearch/script/MCFakes_sample_list.txt,../multiLepSearch/script/MCZPowheg.txt
# file=../multiLepSearch/script/MCZPowheg.txt
#file=../multiLepSearch/script/MCBGZeeSherpaSelected.txt
#file=../multiLepSearch/script/MCBG_llll.txt

k=${tag}.MCBG

#For signal MC
#file=../multiLepSearch/script/MCSig_sample_list.txt
#file=../multiLepSearch/script/MCSig_sample_list_new.txt
#k=${tag}.MCSig

../multiLepSearch/util/run_ss_selection.py --driver grid --inputList ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o ${k} -w -a 1 --study fakes --mcMatch MCTC --doSys 0 --ChargeID 1 --conf ${SUSYconf}

