#!/bin/bash
rc compile

tag=v11.5.1.cf15
grl=GoodRunsLists/data16_13TeV/20170720/physics_25ns_20.7.xml,GoodRunsLists/data15_13TeV/20160720/physics_25ns_20.7.xml
dataPRW=GoodRunsLists/data16_13TeV/20170720/physics_25ns_20.7.lumicalc.OflLumi-13TeV-009.root,GoodRunsLists/data15_13TeV/20160720/physics_25ns_20.7.lumicalc.OflLumi-13TeV-005.root

# dataPRW=GoodRunsLists/data16_13TeV/20170720/physics_25ns_20.7.lumicalc.OflLumi-13TeV-009.root,GoodRunsLists/data15_13TeV/20160720/physics_25ns_20.7.xml
mcPRW=multiLepSearch/prw_MC/merged_prw_mc15c_Jun15.root,multiLepSearch/prw_MC/Wh_all_merged_PRW.root,dev/PileupReweighting/mc15c_v2_defaults.NotRecommended.prw.root
mcPRW=dev/PileupReweighting/mc15c_v2_defaults.NotRecommended.prw.root
SUSYconf=multiLepSearch/sel_conf/SUSYTools_Wh_update.conf

# Use Samuel's or my framework?
st=fakes
# st=ss

if [ -d ${tag} ]; then
    echo "Previous directories/files with tag "${tag}" found."
    echo "Now rm'ing them"
    rm -r ${tag}.*
fi

if false; then
    ###########
    # Cutflow #
    ###########
    echo "### Executing cutflow for signal sample ############################"
    k=${tag}.Sig
    file=/eos/atlas/user/d/dzhang/Samples/Wh/mc15_13TeV.993822.MGPy8EG_A14N13LO_C1N2_Wh_2L_400_0.merge.DAOD_SUSY2.e5678_a766_a821_r7676_p2949_p2972_pUM999999/merge.DAOD_SUSY2.993822._00001.pool.root.1
    ../multiLepSearch/util/run_ss_selection.py --driver direct -f ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o ${k} -w -a 1 --study ${st} --mcMatch MCTC --doSys 0 --ChargeID 1 --conf ${SUSYconf} --printEvent 10566 --cutflow | tee Sig.log  
   
    echo "### Executing cutflow for Zmm sample ############################"
    k=${tag}.Zmm
    file=/eos/atlas/user/g/ggallard/xAOD_forTesting/mc15_13TeV.364104.Sherpa_221_NNPDF30NNLO_Zmumu_MAXHTPTV70_140_CFilterBVeto.merge.DAOD_SUSY2.e5271_s2726_r7772_r7676_p2949/mc15_13TeV/DAOD_SUSY2.11084041._000011.pool.root.1
    ../multiLepSearch/util/run_ss_selection.py --driver direct -f ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o ${k} -w -a 1 --study ${st} --mcMatch dR --doSys 0 --ChargeID 1 --conf ${SUSYconf} --printEvent 47 --cutflow | tee Zmm.log 
    
    echo "### Executing cutflow for ttbar sample ############################"
    k=${tag}.ttbar
    file=/eos/atlas/user/g/ggallard/xAOD_forTesting/mc15_13TeV.410009.PowhegPythiaEvtGen_P2012_ttbar_hdamp172p5_dil.merge.DAOD_SUSY2.e4511_s2608_s2183_r7725_r7676_p2949/DAOD_SUSY2.10524294._000001.pool.root.1
    ../multiLepSearch/util/run_ss_selection.py --driver direct -f ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o ${k} -w -a 1 --study ${st} --mcMatch MCTC --doSys 0 --ChargeID 1 --conf ${SUSYconf} --printEvent 5799815 --cutflow  | tee ttbar.log 

    echo "### Executing cutflow for VV sample ############################"
    k=${tag}.VV
    file=/eos/atlas/user/g/ggallard/xAOD_forTesting/mc15_13TeV.363491.Sherpa_221_NNPDF30NNLO_lllv.merge.DAOD_SUSY2.e5332_s2726_r7772_r7676_p2949/DAOD_SUSY2.10506279._000001.pool.root.1
    ../multiLepSearch/util/run_ss_selection.py --driver direct -f ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o ${k} -w -a 1 --study ${st} --mcMatch MCTC --doSys 0 --ChargeID 1 --conf ${SUSYconf} --printEvent 263 --cutflow | tee VV.log 

    echo "### Executing cutflow for data sample ############################"
    k=${tag}.data
    file=/eos/atlas/user/g/ggallard/xAOD_forTesting/data16_13TeV/DAOD_SUSY2.10314997._000013.pool.root.1
    ../multiLepSearch/util/run_ss_selection.py --driver direct -f ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o ${k} -w -a 0 --study ${st} --mcMatch dR --doSys 0 --ChargeID 1 --conf ${SUSYconf} --cutflow --grl ${grl} | tee data.log  

    # echo "### Executing cutflow for Zee sample ############################"
    # k=${tag}.Zee
    # file=/eos/atlas/user/g/ggallard/xAOD_forTesting/mc15_13TeV.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee.merge.DAOD_SUSY2.e3601_s2576_s2132_r7725_r7676_p2949/DAOD_SUSY2.11608161._000011.pool.root.1
    # ../multiLepSearch/util/run_ss_selection.py --driver direct -f ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o ${k} -w -a 1 --study ${st} --mcMatch dR --doSys 0 --ChargeID 1 --conf ${SUSYconf} --printEvent 47 --cutflow | tee Zmm.log 

    ./printCutflow.py ${tag} > ${tag}.log

else 
    ###################
    # Feature testing #
    ###################
    
    echo "### Executing test for signal sample ############################"
    k=${tag}.Sig
    file=/eos/atlas/user/d/dzhang/Samples/Wh/mc15_13TeV.993822.MGPy8EG_A14N13LO_C1N2_Wh_2L_400_0.merge.DAOD_SUSY2.e5678_a766_a821_r7676_p2949_p2972_pUM999999/merge.DAOD_SUSY2.993822._00001.pool.root.1
    ../multiLepSearch/util/run_ss_selection.py --driver direct -f ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o ${k} -w -a 1 --study ${st} --mcMatch MCTC --doSys 0 --ChargeID 1 --conf ${SUSYconf} --printEvent 10566 --cutflow -n 200| tee Sig.log

    # k=${tag}.ZeeP.dRZ
    # file=/eos/atlas/user/g/ggallard/xAOD_forTesting/mc15_13TeV.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee.merge.DAOD_SUSY2.e3601_s2576_s2132_r7725_r7676_p2949/DAOD_SUSY2.11608161._000011.pool.root.1
    # ../multiLepSearch/util/run_ss_selection.py --driver direct -f ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o ${k} -w -a 1 --study ${st} --mcMatch dRZ --doSys 0 --ChargeID 1 --conf ${SUSYconf} 
fi
