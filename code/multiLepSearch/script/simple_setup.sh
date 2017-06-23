mkdir AnalysisBase-02-04-30
cd AnalysisBase-02-04-30

#https://gitlab.cern.ch/hku/SUSY2L/tree/master/code/multiLepSearch
git clone https://:@gitlab.cern.ch:8443/hku/SUSY2L.git
cd SUSY2L/code/

setupATLAS
#localSetupRucioClients
#voms-proxy-init -voms atlas
#localSetupPandaClient

#https://svnweb.cern.ch/trac/atlasoff/browser/PhysicsAnalysis/SUSYPhys/SUSYTools/tags/
rcSetup Base,2.4.30

rc find_packages
rc clean
rc compile
mkdir run
cd run

../multiLepSearch/util/run_ss_selection.py -f /afs/cern.ch/work/c/clo/sample/mc15_13TeV.392825.MGPy8EG_A14N23LO_C1N2_Slep_200_180_0p95_2L5.merge.DAOD_SUSY2.e5668_a766_a838_r7676_p2949/DAOD_SUSY2.10748609._000001.pool.root.1 --outputTag test1 -o t1 -w -a 1 --study ss --mcPRW "multiLepSearch/prw_MC/merged_prw_mc15c_Slep0d95.root,dev/PileupReweighting/mc15c_v2_defaults.NotRecommended.prw.root" --ChargeID 1 --conf SUSYTools/SUSYTools_Default.conf

../multiLepSearch/util/run_ss_selection.py -f /afs/cern.ch/work/c/clo/sample/data16_13TeV.00297730.physics_Main.merge.DAOD_SUSY2.f694_m1583_p2880/DAOD_SUSY2.09980028._000001.pool.root.1 --outputTag test1 -o t1 -w -a 0 --study ss

../multiLepSearch/util/run_ss_selection.py -f /afs/cern.ch/work/c/clo/sample/mc15_13TeV.363491.Sherpa_221_NNPDF30NNLO_lllv.merge.DAOD_SUSY2.e5332_s2726_r7772_r7676_p2879/DAOD_SUSY2.10345373._000001.pool.root.1  --outputTag test1 -o t1 -w -a 1 --nevents 10000 --study ss

root -l t1/data-myOutput/test.root

# . ../multiLepSearch/script/MC_submit.sh

