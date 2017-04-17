#!/bin/bash
tag=v15.0
fix=""
extraOpt=""
grl=GoodRunsLists/data16_13TeV/20170215/physics_25ns_20.7.xml,GoodRunsLists/data15_13TeV/20160720/physics_25ns_20.7.xml
dataPRW=GoodRunsLists/data16_13TeV/20170215/physics_25ns_20.7.lumicalc.OflLumi-13TeV-008.root,GoodRunsLists/data15_13TeV/20160720/physics_25ns_20.7.lumicalc.OflLumi-13TeV-005.root
mcPRW=dev/SUSYTools/merged_prw_mc15c_Jan31.root,multiLepSearch/prw_MC/merged_prw_mc15c_Slep0d95.root,dev/PileupReweighting/mc15c_v2_defaults.NotRecommended.prw.root
k=${tag}.MC

#####################
#For background MC
#file=../multiLepSearch/script/MCBGZjetsSherpa_sample_list.txt
#file=../multiLepSearch/script/MCBGVVSherpa_sample_list.txt
#file=../multiLepSearch/script/MCBGVgammaSherpa_sample_list.txt
#file=../multiLepSearch/script/MCBGDYSherpa_sample_list.txt
#file=../multiLepSearch/script/MCBGZeeSherpaSelected.txt
#file=../multiLepSearch/script/MCBG_llll.txt
# file=../multiLepSearch/script/MCBGZjetsSherpa_sample_list.txt,../multiLepSearch/script/MCBGVVSherpa_sample_list.txt,../multiLepSearch/script/MCBGVgammaSherpa_sample_list.txt,../multiLepSearch/script/MCBGDYSherpa_sample_list.txt
# file=${file},../multiLepSearch/script/MCBG_sample_list.txt
# # file=remaining.txt
# fix=".3"
# 
# ../multiLepSearch/util/run_ss_selection.py --driver grid --inputList ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o ${k}BG${fix} -a 1 --study ss --mcMatch TruthLink --doSys 0 --conf multiLepSearch/sel_conf/SUSYTools_multilepAnaMoriond.conf

# Add extra option for "broken" jobs resubmission
# ../multiLepSearch/util/run_ss_selection.py --driver grid --inputList ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o ${k}BG${fix} -a 1 --study ss --mcMatch TruthLink --doSys 0 --extraOptions "--allowTaskDuplication"

# ../multiLepSearch/util/run_ss_selection.py --driver grid --inputDS mc15_13TeV.*.MGPy8EG_A14N23LO_C1N2_Slep_*_0p95_2L5.merge.DAOD_SUSY2.*_p2949 --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o ${k} -a 1 --study ss --mcMatch TruthLink --doSys 0
################
file=../multiLepSearch/script/MCBGZjetsSherpa_sample_list.txt,../multiLepSearch/script/MCBGVVSherpa_sample_list.txt,../multiLepSearch/script/MCBGVgammaSherpa_sample_list.txt,../multiLepSearch/script/MCBGDYSherpa_sample_list.txt,../multiLepSearch/script/MCSig_sample_p2949.txt,../multiLepSearch/script/MCBG_sample_list.txt

fix=".1"
# extraOpt="--nFilesPerNode 20"
# ../multiLepSearch/util/run_ss_selection.py --driver condor --samplesDir /net/s3_datad/Data15/MC15/SUSY2/ --sampleList ${file} --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o output/${k}${fix} -a 1 --study ss --mcMatch TruthLink --doSys 0 $extraOpt --conf multiLepSearch/sel_conf/SUSYTools_multilepAnaMoriond.conf

../multiLepSearch/util/run_ss_selection.py --driver grid --inputDS mc15_13TeV.*.Sherpa_221_NNPDF30NNLO_W*MAXH*SUSY2*p2879 --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o output/${k}${fix} -a 1 --study ss --mcMatch TruthLink --doSys 0 $extraOpt --conf multiLepSearch/sel_conf/SUSYTools_multilepAnaMoriond.conf
################


# # ../multiLepSearch/util/run_ss_selection.py --driver condor --samplesDir /net/s3_datad/Data15/MC15/SUSY2 --samplePattern ".*/mc15_13TeV.361627.(.*).merge.DAOD_SUSY2.e4093_s2608_s2183_r7772_r7676_p2879" --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o output/${k} -a 1 --study ss --mcMatch TruthLink --doSys 1
# ../multiLepSearch/util/run_ss_selection.py --driver condor --samplesDir /net/s3_datad/Data15/MC15/SUSY2 --samplePattern ".*/mc15_13TeV.(\d*\..*).merge.DAOD_SUSY2.(.*)" --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o output/${k} -a 1 --study ss --mcMatch TruthLink --doSys 1
# ../multiLepSearch/util/run_ss_selection.py --driver condor --samplesDir /net/s3_datad/Data15/MC15/SUSY2 --samplePattern ".*/mc15_13TeV.(41001\d.PowhegPythiaEvtGen_P2012_Wt_dilepton.*).merge.DAOD_SUSY2.(.*)" --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o output/${k} -a 1 --study ss --mcMatch TruthLink --doSys 0

###############################################
# ../multiLepSearch/util/run_ss_selection.py --driver condor --samplesDir /net/s3_datad/Data15/MC15/SUSY2 --samplePattern ".*/mc15_13TeV.(\d*\..*).0p95_2L5.merge.DAOD_SUSY2.(.*)" --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o output/${k} -a 1 --study ss --mcMatch TruthLink --doSys 1

## For testing
# ../multiLepSearch/util/run_ss_selection.py -d /net/s3_datad/Data15/MC15/SUSY2/mc15_13TeV.364131.Sherpa_221_NNPDF30NNLO_Ztautau_MAXHTPTV70_140_CVetoBVeto.merge.DAOD_SUSY2.e5307_s2726_r7772_r7676_p2879/ --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o ${k} -w -a 1 --study ss --mcMatch TruthLink --doSys 0
# ../multiLepSearch/util/run_ss_selection.py -d /net/s3_datad/Data15/MC15/SUSY2/mc15_13TeV.392880.MGPy8EG_A14N23LO_C1N2_Slep_900_0_0p95_2L5.merge.DAOD_SUSY2.e5668_a766_a821_r7676_p2949 --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o ${k} -w -a 1 --study ss --mcMatch TruthLink --doSys 0
