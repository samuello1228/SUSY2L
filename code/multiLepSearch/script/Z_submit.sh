#!/bin/bash
dataPRW=GoodRunsLists/data16_13TeV/20170215/physics_25ns_20.7.lumicalc.OflLumi-13TeV-008.root,GoodRunsLists/data15_13TeV/20160720/physics_25ns_20.7.lumicalc.OflLumi-13TeV-005.root
mcPRW=multiLepSearch/prw_MC/merged_prw_mc15c_Jun15.root,multiLepSearch/prw_MC/merged_prw_mc15c_Slep0d95.root,dev/PileupReweighting/mc15c_v2_defaults.NotRecommended.prw.root
# SUSYconf=multiLepSearch/sel_conf/SUSYTools_multilepAnaMoriond.conf
SUSYconf=multiLepSearch/sel_conf/SUSYTools_Wh.conf

####
tag=v11.4.1.MCTCZ
echo "########## MCTCZ Sherpa #############"
k=${tag}.Sherpa
file=/eos/atlas/user/g/ggallard/xAOD_forTesting/mc15_13TeV.364118.Sherpa_221_NNPDF30NNLO_Zee_MAXHTPTV70_140_CFilterBVeto.merge.DAOD_SUSY2.e5299_s2726_r7772_r7676_p2949/DAOD_SUSY2.11084342._000011.pool.root.1
../multiLepSearch/util/run_ss_selection.py --driver direct -f ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o ${k} -w -a 1 --study fakes --mcMatch MCTCZ --doSys 0 --ChargeID 1 --conf ${SUSYconf}

echo "########## MCTCZ PowhegPythia #############"
k=${tag}.PowhegPythia
file=/eos/atlas/user/g/ggallard/xAOD_forTesting/mc15_13TeV.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee.merge.DAOD_SUSY2.e3601_s2576_s2132_r7725_r7676_p2949/DAOD_SUSY2.11608161._000011.pool.root.1
../multiLepSearch/util/run_ss_selection.py --driver direct -f ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o ${k} -w -a 1 --study fakes --mcMatch MCTCZ --doSys 0 --ChargeID 1 --conf ${SUSYconf}

####
tag=v11.4.1.dRZ
echo "########## dRZ Sherpa #############"
k=${tag}.Sherpa
file=/eos/atlas/user/g/ggallard/xAOD_forTesting/mc15_13TeV.364118.Sherpa_221_NNPDF30NNLO_Zee_MAXHTPTV70_140_CFilterBVeto.merge.DAOD_SUSY2.e5299_s2726_r7772_r7676_p2949/DAOD_SUSY2.11084342._000011.pool.root.1
../multiLepSearch/util/run_ss_selection.py --driver direct -f ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o ${k} -w -a 1 --study fakes --mcMatch dRZ --doSys 0 --ChargeID 1 --conf ${SUSYconf}

echo "########## dRZ PowhegPythia #############"
k=${tag}.PowhegPythia
file=/eos/atlas/user/g/ggallard/xAOD_forTesting/mc15_13TeV.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee.merge.DAOD_SUSY2.e3601_s2576_s2132_r7725_r7676_p2949/DAOD_SUSY2.11608161._000011.pool.root.1
../multiLepSearch/util/run_ss_selection.py --driver direct -f ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o ${k} -w -a 1 --study fakes --mcMatch dRZ --doSys 0 --ChargeID 1 --conf ${SUSYconf}