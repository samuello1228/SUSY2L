g++ -O3 -Wall -Wextra -std=c++11 -o flip_rates flip_rates.cxx `root-config --cflags --glibs`
./flip_rates ZeeCandidate binsDATA.txt `~/Zee/printFiles.sh /afs/cern.ch/work/a/anburger/public/ntuples_data_ForPrerecommendations23-05-16`
python llh.py binsDATA.txt ratesDATA.txt
