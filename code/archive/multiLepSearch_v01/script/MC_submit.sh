#!/bin/bash
tag=v5.0.s
dataPRW=multiLepSearch/ilumicalc_histograms_None_276262-280614.root
mcPRW=multiLepSearch/mc15a_defaults.NotRecommended.prw.root
for i in `cat ../multiLepSearch/script/MC_sample_list.txt`; do

j=${i##*mc15_13TeV.}
k=${tag}.${j%%.merge.*}
echo $k
../multiLepSearch/util/run_2l_selection.py --driver grid --inputDS $i --dataPRW ${dataPRW} --mcPRW ${mcPRW} --outputTag ${tag} -o ${k} -w -a 1

done
