mkdir AnalysisSUSY-02-03-23
cd AnalysisSUSY-02-03-23

#https://svnweb.cern.ch/trac/atlasinst/browser/Institutes/Michigan/SUSY/code/multiLepSearch/
svn co svn+ssh://svn.cern.ch/reps/atlasinst/Institutes/Michigan/SUSY/code/multiLepSearch

setupATLAS
#localSetupDQ2Client
#voms-proxy-init -voms atlas
#localSetupPandaClient --noAthenaCheck

#https://svnweb.cern.ch/trac/atlasoff/browser/AnalysisSUSY/tags/
#https://svnweb.cern.ch/trac/atlasoff/browser/PhysicsAnalysis/SUSYPhys/SUSYTools/tags/
rcSetup SUSY,2.3.23

rc find_packages
rc clean
rc compile
mkdir run
cd run

../multiLepSearch/util/run_2l_selection.py -f /afs/cern.ch/user/c/clo/work/public/sample/mc15_13TeV.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee.merge.AOD.e3601_s2576_s2132_r6765_r6282_tid05568523_00/AOD.05568523._000001.pool.root.1 -t test1 -o t1 -w -a 1

# . ../multiLepSearch/script/MC_submit.sh

