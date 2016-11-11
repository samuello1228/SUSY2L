mkdir AnalysisBase-02-04-17
cd AnalysisBase-02-04-17

#https://svnweb.cern.ch/trac/atlasinst/browser/Institutes/Michigan/SUSY/code/multiLepSearch/
svn co -r 414670 svn+ssh://svn.cern.ch/reps/atlasinst/Institutes/Michigan/SUSY/code/multiLepSearch

setupATLAS
#localSetupRucioClients
#voms-proxy-init -voms atlas
#localSetupPandaClient currentJedi

#https://svnweb.cern.ch/trac/atlasoff/browser/PhysicsAnalysis/SUSYPhys/SUSYTools/tags/
rcSetup Base,2.4.17
svn co svn+ssh://svn.cern.ch/reps/atlasoff/PhysicsAnalysis/MuonID/MuonIDAnalysis/MuonEfficiencyCorrections/tags/MuonEfficiencyCorrections-03-03-08 MuonEfficiencyCorrections
svn co svn+ssh://svn.cern.ch/reps/atlasoff/PhysicsAnalysis/ElectronPhotonID/ElectronPhotonFourMomentumCorrection/tags/ElectronPhotonFourMomentumCorrection-00-02-13 ElectronPhotonFourMomentumCorrection

rc find_packages
rc clean
rc compile
mkdir run
cd run

../multiLepSearch/util/run_ss_selection.py -f /afs/cern.ch/work/c/clo/public/sample/mc15_13TeV.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee.merge.DAOD_SUSY2.e3601_s2576_s2132_r7725_r7676_p2666/DAOD_SUSY2.08677613._000001.pool.root.1 --dataPRW multiLepSearch/prw_Data/ilumicalc_histograms_None_297730-300908_v78.root --outputTag test1 -o t1 -w -a 1 --study ss

../multiLepSearch/util/run_ss_selection.py -f /afs/cern.ch/work/c/clo/sample/data16_13TeV.00300784.physics_Main.merge.DAOD_SUSY2.f708_m1606_p2667_tid08616061_00/DAOD_SUSY2.08616061._000002.pool.root.1 --outputTag test1 -o t1 -w -a 0 --study ss

../multiLepSearch/util/run_ss_selection.py -f /afs/cern.ch/work/c/clo/sample/mc15_13TeV.361064.Sherpa_CT10_lllvSFMinus.merge.DAOD_SUSY2.e3836_s2608_s2183_r7725_r7676_p2666_tid08648281_00/DAOD_SUSY2.08648288._000002.pool.root.1 --outputTag test1 -o t1 -w -a 1 --nevents 10000 --study ss

../multiLepSearch/util/run_ss_selection.py -f /afs/cern.ch/work/c/clo/sample/mc15_13TeV.361064.Sherpa_CT10_lllvSFMinus.merge.DAOD_SUSY2.e3836_s2608_s2183_r7725_r7676_p2666_tid08648281_00/DAOD_SUSY2.08648288._000002.pool.root.1 --outputTag test1 -o t1 -w -a 1 --study ss

root -l t1/data-myOutput/test.root

# . ../multiLepSearch/script/MC_submit.sh

