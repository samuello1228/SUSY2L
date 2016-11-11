#!/bin/bash
tag=v5.0.s
dataPRW=multiLepSearch/ilumicalc_histograms_None_276262-280614.root
mcPRW=multiLepSearch/mc15a_defaults.NotRecommended.prw.root
grl=multiLepSearch/data15_13TeV.periodAllYear_DetStatus-v67-pro19-02_DQDefects-00-01-02_PHYS_StandardGRL_All_Good.xml
for i in `cat ../multiLepSearch/script/Data_sample_list.txt`; do

j=${i##*data15_13TeV.}
k=${tag}.${j%%.merge.*}
echo $k
../multiLepSearch/util/run_2l_selection.py --driver grid --inputDS $i --dataPRW ${dataPRW} --mcPRW ${mcPRW} --grl ${grl} --outputTag ${tag} -o ${k} -w -a 0

done
