#!/bin/bash
tag=v11.0
dataPRW=GoodRunsLists/data16_13TeV/20170720/physics_25ns_20.7.lumicalc.OflLumi-13TeV-009.root,GoodRunsLists/data15_13TeV/20160720/physics_25ns_20.7.lumicalc.OflLumi-13TeV-005.root
mcPRW=multiLepSearch/prw_MC/merged_prw_mc15c_Jun15.root,multiLepSearch/prw_MC/Wh_all_merged_PRW.root,dev/PileupReweighting/mc15c_v2_defaults.NotRecommended.prw.root
#SUSYconf=SUSYTools/SUSYTools_Default.conf
SUSYconf=multiLepSearch/sel_conf/SUSYTools_Wh_update.conf
#SUSYconf=multiLepSearch/sel_conf/SUSYTools_Wh_update_Loose.conf

#For background MC
#mcPRW=dev/SUSYTools/merged_prw_mc15c_latest.root

#Samuel
#file=../multiLepSearch/script/MCBGsingletop_sample_list.txt,../multiLepSearch/script/MCBGttbar_sample_list.txt,../multiLepSearch/script/MCBGttV_sample_list.txt,../multiLepSearch/script/MCBGVVSherpa_sample_list.txt,../multiLepSearch/script/MCBGVgammaSherpa_sample_list.txt,../multiLepSearch/script/MCBGhiggs_sample_list.txt
#file=../multiLepSearch/script/MCBGsingletop_sample_list.txt
#file=../multiLepSearch/script/MCBGttbar_sample_list.txt
#file=../multiLepSearch/script/MCBGttV_sample_list.txt
#file=../multiLepSearch/script/MCBGVVSherpa_sample_list.txt
#file=../multiLepSearch/script/MCBGVgammaSherpa_sample_list.txt
#file=../multiLepSearch/script/MCBGhiggs_sample_list.txt

#Dongliang
file=../multiLepSearch/script/MCBGmultitop_sample_list.txt,../multiLepSearch/script/MCBGVVVSherpa_sample_list.txt,../multiLepSearch/script/MCBGWjetsSherpa_sample_list.txt,../multiLepSearch/script/MCBGZjetsSherpa_sample_list.txt,../multiLepSearch/script/MCBGDYSherpa_sample_list.txt
#file=../multiLepSearch/script/MCBGWjetsSherpa_sample_list.txt
#file=../multiLepSearch/script/MCBGZjetsSherpa_sample_list.txt
#file=../multiLepSearch/script/MCBGDYSherpa_sample_list.txt
#file=../multiLepSearch/script/MCBGmultitop_sample_list.txt
#file=../multiLepSearch/script/MCBGVVVSherpa_sample_list.txt

k=${tag}.MCBG

#For signal MC
#file=../multiLepSearch/script/MCSig_sample_list.txt
#k=${tag}.MCSig

../multiLepSearch/util/run_ss_selection.py --driver grid --inputList ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${k} -o ${k} -a 1 --study ss --mcMatch TruthLink --doSys 0 --ChargeID 1 --conf ${SUSYconf} --fast


# file="mc15_13TeV.*.MGPy8EG_A14N23LO_C1N2_Wh_hall_*_2L7.merge.DAOD_SUSY2.e6153_a766_a821_r7676_p2949"
# ../multiLepSearch/util/run_ss_selection.py --driver grid --inputDS ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${k} -o ${k} -a 1 --study ss --mcMatch TruthLink --doSys 0 --ChargeID 1 --conf ${SUSYconf} --fast

#For running on local files

#file=/afs/cern.ch/work/c/clo/sample/mc15_13TeV.364104.Sherpa_221_NNPDF30NNLO_Zmumu_MAXHTPTV70_140_CFilterBVeto.merge.DAOD_SUSY2.e5271_s2726_r7772_r7676_p2949/DAOD_SUSY2.11084041._000011.pool.root.1
#eventNumber=47

#file=/afs/cern.ch/work/c/clo/sample/mc15_13TeV.410009.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_dil.merge.DAOD_SUSY2.e4511_s2608_s2183_r7725_r7676_p2949/DAOD_SUSY2.10524294._000001.pool.root.1
#eventNumber=5799815

#file=/afs/cern.ch/work/c/clo/sample/mc15_13TeV.363491.Sherpa_221_NNPDF30NNLO_lllv.merge.DAOD_SUSY2.e5332_s2726_r7772_r7676_p2949/DAOD_SUSY2.10506279._000001.pool.root.1
#eventNumber=263

#file=/eos/atlas/user/d/dzhang/Samples/Wh/mc15_13TeV.993820.MGPy8EG_A14N13LO_C1N2_Wh_2L_175_0.merge.DAOD_SUSY2.e5678_a766_a821_r7676_p2949_p2972_pUM999999/merge.DAOD_SUSY2.993820._00001.pool.root.1

#file=/eos/atlas/user/d/dzhang/Samples/Wh/mc15_13TeV.993821.MGPy8EG_A14N13LO_C1N2_Wh_2L_165_35.merge.DAOD_SUSY2.e5678_a766_a821_r7676_p2949_p2972_pUM999999/merge.DAOD_SUSY2.993821._00001.pool.root.1

#file=/eos/atlas/user/d/dzhang/Samples/Wh/mc15_13TeV.993822.MGPy8EG_A14N13LO_C1N2_Wh_2L_400_0.merge.DAOD_SUSY2.e5678_a766_a821_r7676_p2949_p2972_pUM999999/merge.DAOD_SUSY2.993822._00001.pool.root.1
#eventNumber=10566

#file=/afs/cern.ch/user/c/clo/work/sample/mc15_13TeV.364162.Sherpa_221_NNPDF30NNLO_Wmunu_MAXHTPTV140_280_CVetoBVeto.merge.DAOD_SUSY2.e5340_s2726_r7772_r7676_p2949/DAOD_SUSY2.11134137._000011.pool.root.1

#mcPRW=dev/PileupReweighting/mc15c_v2_defaults.NotRecommended.prw.root

#../multiLepSearch/util/run_ss_selection.py --driver direct -f ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} -o t1 -w -a 1 --study ss --mcMatch TruthLink --doSys 0 --ChargeID 1 --conf ${SUSYconf} # --printEvent ${eventNumber} --cutflow | tee cutflow.log

#grep -e '^[0-9]\+$' cutflow.log |tee trigger.txt
