g++ -O3 -Wall -Wextra -std=c++11 -o flip_rates flip_rates.cxx `root-config --cflags --glibs`
./flip_rates ZeeCandidate binsMC.txt `~/Zee/printFiles.sh ~/work/NTuples/user.kristin.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee.e3601_s2576_s2132_r7267_r6282`
python llh.py binsMC.txt ratesMC.txt
