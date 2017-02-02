#!/bin/bash
tag=v9.1

grl=GoodRunsLists/data16_13TeV/20161101/physics_25ns_20.7.xml,GoodRunsLists/data15_13TeV/20160720/physics_25ns_20.7.xml
dataPRW=GoodRunsLists/data16_13TeV/20161101/physics_25ns_20.7.lumicalc.OflLumi-13TeV-005.root,GoodRunsLists/data15_13TeV/20160720/physics_25ns_20.7.lumicalc.OflLumi-13TeV-005.root

#For background MC
#mcPRW=multiLepSearch/prw_MC/merged_prw_mc15c.root
# file=../multiLepSearch/script/MCBG_sample_list.txt
#file=../multiLepSearch/script/MCBGZjetsSherpa_sample_list.txt
#file=../multiLepSearch/script/MCBGVVSherpa_sample_list.txt
#file=../multiLepSearch/script/MCBGVgammaSherpa_sample_list.txt
#file=../multiLepSearch/script/MCBGZeeSherpaSelected.txt

k=${tag}.MCBG
mcPRW=multiLepSearch/prw_MC/merged_prw_mc15c_July27_afterFix.root

# # ../multiLepSearch/util/run_ss_selection.py --driver condor --samplesDir /net/s3_datad/Data15/MC15/SUSY2 --samplePattern ".*/mc15_13TeV.361627.(.*).merge.DAOD_SUSY2.e4093_s2608_s2183_r7772_r7676_p2879" --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o output/${k} -a 1 --study ss --mcMatch TruthLink --doSys 1
# ../multiLepSearch/util/run_ss_selection.py --driver condor --samplesDir /net/s3_datad/Data15/MC15/SUSY2 --samplePattern ".*/mc15_13TeV.(\d*\..*).merge.DAOD_SUSY2.(.*)" --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o output/${k} -a 1 --study ss --mcMatch TruthLink --doSys 1
../multiLepSearch/util/run_ss_selection.py --driver condor --samplesDir /net/s3_datad/Data15/MC15/SUSY2 --samplePattern ".*/mc15_13TeV.(41001\d.PowhegPythiaEvtGen_P2012_Wt_dilepton.*).merge.DAOD_SUSY2.(.*)" --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o output/${k} -a 1 --study ss --mcMatch TruthLink --doSys 0

###############################################
# k=${tag}.MCSig
# mcPRW=multiLepSearch/prw_MC/merged_prw_mc15c_Slep0d95.root,dev/PileupReweighting/mc15c_v2_defaults.NotRecommended.prw.root
# 
# ../multiLepSearch/util/run_ss_selection.py --driver condor --samplesDir /net/s3_datad/Data15/MC15/SUSY2 --samplePattern ".*/mc15_13TeV.(\d*\..*).0p95_2L5.merge.DAOD_SUSY2.(.*)" --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o output/${k} -w -a 1 --study ss --mcMatch TruthLink --doSys 1

## For testing
# ../multiLepSearch/util/run_ss_selection.py -d /net/s3_datad/Data15/MC15/SUSY2/mc15_13TeV.392880.MGPy8EG_A14N23LO_C1N2_Slep_900_0_0p95_2L5.merge.DAOD_SUSY2.e5129_a766_a821_r7676_p2688 --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o ${k} -w -a 1 --study ss --mcMatch TruthLink --doSys 1

# ../multiLepSearch/util/run_ss_selection.py -d /net/s3_datad/Data15/MC15/SUSY2/mc15_13TeV.392883.MGPy8EG_A14N23LO_C1N2_Slep_1000_700_0p95_2L5.merge.DAOD_SUSY2.e5129_a766_a821_r7676_p2688 --outputTag ${tag} -o ${k} -w -a 1 --study ss --mcMatch TruthLink --doSys 1 --dataPRW ${dataPRW} --mcPRW ${mcPRW}
