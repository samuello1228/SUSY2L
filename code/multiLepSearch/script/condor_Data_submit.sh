#!/bin/bash
tag=v14.0
# grlDir=
# /cvmfs/atlas.cern.ch/repo/sw/database/GroupData/GoodRunsLists/data16_13TeV/20161101/physics_25ns_20.7.xml

grl=GoodRunsLists/data16_13TeV/20170215/physics_25ns_20.7.xml,GoodRunsLists/data15_13TeV/20160720/physics_25ns_20.7.xml
dataPRW=GoodRunsLists/data16_13TeV/20170215/physics_25ns_20.7.lumicalc.OflLumi-13TeV-008.root,GoodRunsLists/data15_13TeV/20160720/physics_25ns_20.7.lumicalc.OflLumi-13TeV-005.root

k=${tag}.data 
extraOpt='-w'
../multiLepSearch/util/run_ss_selection.py --driver condor --samplesDir /net/s3_datad/Data15/DxAOD/SUSY/ --samplePattern ".*/(data1.)_13TeV.(.*).physics_Main.PhysCont.DAOD_SUSY2.*" --dataPRW ${dataPRW} --grl ${grl} --outputTag ${tag} -o output/${k} -a 0 --study ss --doSys 1 --nFilesPerNode 300 --conf multiLepSearch/sel_conf/SUSYTools_multilepAnaMoriond.conf 

# ../multiLepSearch/util/run_ss_selection.py -f /net/s3_datad/Data15/DxAOD/SUSY/data15_13TeV.periodH.physics_Main.PhysCont.DAOD_SUSY2.grp15_v01_p2880/DAOD_SUSY2.09982394._000146.pool.root.1 --dataPRW ${dataPRW} --grl ${grl} --outputTag ${tag} --conf multiLepSearch/sel_conf/SUSYTools_multilepAnaMoriond.conf -o testD -w -a 0 --study ss --doSys 0

# ../multiLepSearch/util/run_zph_selection.py -f /net/s3_datad/Data15/DxAOD/SUSY/data15_13TeV.periodH.physics_Main.PhysCont.DAOD_SUSY2.grp15_v01_p2880/DAOD_SUSY2.09982394._000146.pool.root.1 --dataPRW ${dataPRW} --grl ${grl} --outputTag ${tag} -o testD -w -a 0 --study ss --doSys 0

# ../multiLepSearch/util/run_zph_selection.py --driver condor -d /net/s3_datad/Data15/DxAOD/SUSY/data15_13TeV.periodH.physics_Main.PhysCont.DAOD_SUSY2.grp15_v01_p2880/ --dataPRW ${dataPRW} --grl ${grl} --outputTag ${tag} -o testP -w -a 0 --study ss --doSys 0

# ../multiLepSearch/util/run_zph_selection.py --driver condor --samplesDir /net/s3_datad/Data15/DxAOD/SUSY/ --samplePattern ".*/(data1.)_13TeV.(.*).physics_Main.PhysCont.DAOD_SUSY2.*" --dataPRW ${dataPRW} --grl ${grl} --outputTag ${tag} -o output/${k}_zph -a 0 --study ss --doSys 1 --nFilesPerNode 300 $extraOpt
 
