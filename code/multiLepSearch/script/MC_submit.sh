#!/bin/bash
tag=v8.12
dataPRW=multiLepSearch/prw_Data/ilumicalc_histograms_None_297730-311481_OflLumi-13TeV-005.root

#For background MC
#mcPRW=multiLepSearch/prw_MC/merged_prw_mc15c_July27_afterFix.root
#mcPRW=multiLepSearch/prw_MC/merged_prw_mc15c_Jan31.root
mcPRW=dev/SUSYTools/merged_prw_mc15c_latest.root

file=../multiLepSearch/script/MCBG_sample_list.txt
#file=../multiLepSearch/script/MCBGZjetsSherpa_sample_list.txt
#file=../multiLepSearch/script/MCBGVVSherpa_sample_list.txt
#file=../multiLepSearch/script/MCBGVgammaSherpa_sample_list.txt
#file=../multiLepSearch/script/MCBGZeeSherpaSelected.txt
#file=../multiLepSearch/script/MCBG_llll.txt

k=${tag}.MCBG
#k=${tag}.MCZee

#For signal MC
#mcPRW=multiLepSearch/prw_MC/merged_prw_mc15c_Slep0d95.root
#mcPRW="multiLepSearch/prw_MC/merged_prw_mc15c_Slep0d95.root,dev/PileupReweighting/mc15c_v2_defaults.NotRecommended.prw.root"
#file=../multiLepSearch/script/MCSig_sample_list.txt
#k=${tag}.MCSig

../multiLepSearch/util/run_ss_selection.py --driver grid --inputList ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o ${k} -w -a 1 --study ss --mcMatch dR --doSys 0 --ChargeID 1

