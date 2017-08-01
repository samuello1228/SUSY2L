mkdir AnalysisBase-02-04-31
cd AnalysisBase-02-04-31

#https://gitlab.cern.ch/hku/SUSY2L/tree/master/code/multiLepSearch
git clone https://:@gitlab.cern.ch:8443/hku/SUSY2L.git
cd SUSY2L/code/

setupATLAS
#voms-proxy-init -voms atlas
#lsetup rucio
#lsetup panda

#https://svnweb.cern.ch/trac/atlasoff/browser/PhysicsAnalysis/SUSYPhys/SUSYTools/tags/
rcSetup Base,2.4.31

rc find_packages
rc clean
rc compile
mkdir run
cd run

../multiLepSearch/util/run_ss_selection.py -f /afs/cern.ch/work/c/clo/sample/mc15_13TeV.392847.MGPy8EG_A14N23LO_C1N2_Slep_500_450_0p95_2L5.merge.DAOD_SUSY2.e5668_a766_a821_r7676_p2949/DAOD_SUSY2.10662422._000002.pool.root.1 -o t1 -w -a 1 --study ss --mcPRW "multiLepSearch/prw_MC/merged_prw_mc15c_Slep0d95.root,dev/PileupReweighting/mc15c_v2_defaults.NotRecommended.prw.root" --ChargeID 1 --conf SUSYTools/SUSYTools_Default.conf

root -l t1/data-myOutput/test.root

# . ../multiLepSearch/script/MC_submit.sh

