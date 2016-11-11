#!/bin/bash
#https://twiki.cern.ch/twiki/bin/view/AtlasProtected/ExtendedPileupReweighting

#setup:
#setupATLAS
#voms-proxy-init -voms atlas
#asetup AthAnalysisBase,2.4.15,here
#athena -c "FILES='/afs/cern.ch/work/c/clo/public/sample/mc15_13TeV:mc15_13TeV.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee.merge.AOD.e3601_s2576_s2132_r7725_r7676/AOD.07981919._000001.pool.root.1';" PileupReweighting/generatePRW_jobOptions.py

#localSetupPandaClient currentJedi
tag=v1.5
#for i in `cat ../multiLepSearch/script/MCBG_sample_list.txt`; do
for i in `cat ../multiLepSearch/script/MCSig_sample_list.txt`; do
a=${i%%.DAOD_SUSY2.*}.AOD.${i##*.DAOD_SUSY2.}
a=${a%%_p*}
echo $a

j=${i##*mc15_13TeV.}
k=user.clo.${tag}.prw.${j%%.merge.*}
echo $k
pathena PileupReweighting/generatePRW_jobOptions.py --inDS ${a} --outDS ${k}

done
