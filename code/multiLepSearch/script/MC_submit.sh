#!/bin/bash
tag=v10.1
dataPRW=GoodRunsLists/data16_13TeV/20170215/physics_25ns_20.7.lumicalc.OflLumi-13TeV-008.root,GoodRunsLists/data15_13TeV/20160720/physics_25ns_20.7.lumicalc.OflLumi-13TeV-005.root
mcPRW=multiLepSearch/prw_MC/merged_prw_mc15c_Jun15.root,multiLepSearch/prw_MC/merged_prw_mc15c_Slep0d95.root,dev/PileupReweighting/mc15c_v2_defaults.NotRecommended.prw.root
SUSYconf=SUSYTools/SUSYTools_Default.conf
#SUSYconf=multiLepSearch/sel_conf/SUSYTools_Default_CFT.conf

#For background MC
#mcPRW=dev/SUSYTools/merged_prw_mc15c_latest.root

#Samuel
#file=../multiLepSearch/script/MCBGsingletop_sample_list.txt
#file=../multiLepSearch/script/MCBGttbar_sample_list.txt
#file=../multiLepSearch/script/MCBGttV_sample_list.txt
#file=../multiLepSearch/script/MCBGVVSherpa_sample_list.txt

#Dongliang
file=../multiLepSearch/script/MCBGWjetsSherpa_sample_list.txt,../multiLepSearch/script/MCBGZjetsSherpa_sample_list.txt,../multiLepSearch/script/MCBGDYSherpa_sample_list.txt
#file=../multiLepSearch/script/MCBGWjetsSherpa_sample_list.txt
#file=../multiLepSearch/script/MCBGZjetsSherpa_sample_list.txt
#file=../multiLepSearch/script/MCBGDYSherpa_sample_list.txt

#not ready
#file=../multiLepSearch/script/MCBGmultitop_sample_list.txt
#VVV
#file=../multiLepSearch/script/MCBGVgammaSherpa_sample_list.txt
#file=../multiLepSearch/script/MCBGhiggs_sample_list.txt

k=${tag}.MCBG

#For signal MC
#file=../multiLepSearch/script/MCSig_sample_list.txt
#file=../multiLepSearch/script/MCSig_sample_list_new.txt
#k=${tag}.MCSig

../multiLepSearch/util/run_ss_selection.py --driver grid --inputList ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o ${k} -w -a 1 --study ss --mcMatch dR --doSys 0 --ChargeID 1 --conf ${SUSYconf}

#file=/afs/cern.ch/work/c/clo/sample/mc15_13TeV.392847.MGPy8EG_A14N23LO_C1N2_Slep_500_450_0p95_2L5.merge.DAOD_SUSY2.e5668_a766_a821_r7676_p2949/DAOD_SUSY2.10662422._000002.pool.root.1

#file=/afs/cern.ch/work/c/clo/sample/mc15_13TeV.364104.Sherpa_221_NNPDF30NNLO_Zmumu_MAXHTPTV70_140_CFilterBVeto.merge.DAOD_SUSY2.e5271_s2726_r7772_r7676_p2949/DAOD_SUSY2.11084041._000011.pool.root.1

#file=/eos/atlas/user/d/dzhang/Samples/Wh/mc15_13TeV.993820.MGPy8EG_A14N13LO_C1N2_Wh_2L_175_0.merge.DAOD_SUSY2.e5678_a766_a821_r7676_p2949_p2972_pUM999999/merge.DAOD_SUSY2.993820._00001.pool.root.1

#file=/eos/atlas/user/d/dzhang/Samples/Wh/mc15_13TeV.993821.MGPy8EG_A14N13LO_C1N2_Wh_2L_165_35.merge.DAOD_SUSY2.e5678_a766_a821_r7676_p2949_p2972_pUM999999/merge.DAOD_SUSY2.993821._00001.pool.root.1

#file=/eos/atlas/user/d/dzhang/Samples/Wh/mc15_13TeV.993822.MGPy8EG_A14N13LO_C1N2_Wh_2L_400_0.merge.DAOD_SUSY2.e5678_a766_a821_r7676_p2949_p2972_pUM999999/merge.DAOD_SUSY2.993822._00001.pool.root.1

#../multiLepSearch/util/run_ss_selection.py --driver direct -f ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} -o t1 -w -a 1 --study ss --mcMatch dR --doSys 0 --ChargeID 1 --conf ${SUSYconf}
