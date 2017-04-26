# BDT for SUSY2L analysis

Prepared by Gabriel Gallardo Apr 2017

## What is a BDT?
A BDT is a boosted decision tree. It is one of the multi-variable analysis techniques implemented in the [TMVA](http://tmva.sourceforge.net/) package of ROOT. 

For details about TMVA documentation of this package, please refer to its [user guide](http://tmva.sourceforge.net/docu/TMVAUsersGuide.pdf). 

## BDT optimization workflow
The BDT workflow consists of the following steps:
- Generate NTuple from xAOD using the [multiLepSearch](../multiLepSearch/README) package *(documentation may not be up to date)*
- [Train BDTs with different hyper-parameters](#bdt-training)
- [Compare overtraining and performance of trained BDTs](#overtraining-and-sensitivity-checks)
- Apply BDT discriminant to NTuples


## BDT training

The BDT training is done by `optWithTMVA.py`.

### Set up
To execute this script, you need to have RootCore set up in your directory. This is because the script uses scripts from both `multiLepSearch` and `SUSYTools`.

One way to approach this is to first check out the SUSY2L git project to some location, and create soft links from the multiLepSearch package, the CutOpt directory, and the optWithTMVA.py to another directory in which you will execute the script.

```sh
# Checkout git repository
cd ~/Documents/
kinit -f username@CERN.CH
git clone https://:@gitlab.cern.ch:8443/hku/SUSY2L.git

# Set up workspace
mkdir opt
cd opt
ln -s ../SUSY2L/code/multiLepSearch ../SUSY2L/code/CutOpt ../SUSY2L/code/CutOpt .

# Prepare for execution
cd ~/Documents/opt
setupATLAS
rcSetup
rc find_packages
rc compile

```

### Local Execution

```sh
python CutOpt/optWithTMVA.py -b $bkgFiles -s $sigFiles -d $dataFiles -c 0 -o $outFile --NTrees $nTrees --NodeSize $nodeSize --Depth $depth
```

**List of options:**   
- '-b', "--bkgList": list of BKG input ROOT file
- '-s', "--sigList": list of SIG input ROOT file
- '-d', "--dataList": list of DATA input ROOT file
- '-t', "--nominalTreeName": "name of the nominal tree in input file, default="evt2l"
- '-l', "--inputLumi": total Luminosity(fb-1) in dataList samples, default=33257.2
- '-o', "--outFile": name of output file, default="testoutput.root"
- '-c', "--channel": Lepton flavor channels, see [multiLepSearch/python/ssUtil.py](../multiLepSearch/python/ssUtil.py) for documentation
- "--NodeSize": minimum node size (in %) of the trees, default=5
- "--NTrees": number of trees in forest, default=400
- "--Depth": tree depth, default=2

### "Batch" execution on HKU computer

The script `runAll.sh` runs a set of BDT training jobs with different hyper-parameters. This is for BDT optimization.

You can customize your parameter list and input file list in the first few lines of the script:

```sh
channels="1 2 3 4 10 11 12 13 14" # Channels as defined in ssUtil # Note that channel 0 has to be done separately because of some quirk of the script
masses="20 50 100"
nTreesList="100 200 400 600"
nodeSizeList="5 7 10"
depthList="2 3 4 5"

bkgFiles="CutOpt/GabrielFiles/bkgFiles.txt"
dataFiles="CutOpt/GabrielFiles/dataFiles.txt"

# For the signal files, replace the sigFiles variable with the appropriate expression:
sigFiles="CutOpt/GabrielFiles/sigFiles_"$mass".txt"

```


The script executes a small number of jobs at a time, waits for them to complete, then executes another small set of jobs. This is to make use of the RAM and CPU limitations of the HKU computer

### Batch execution on Condor

Condor submission is still a work in progress

## Overtraining and sensitivity checks

Workflow: 
- makeOTPlots.py
- calcCLs.py
- Process output CSV's in spreadsheet application

**makeOTPlots.py**  
Execute `makeOTPlots.py`, which is dependent on 
- `overtraining.py` for the overtraining checks, and
- `mvaeffs.cxx` for the signal/background efficiency and sensitivity checks checks
- [tqdm](https://github.com/tqdm/tqdm) a python package which allows for the use of good-looking, easy to configure progress bars. Install by:
  +  `pip install --user tqdm` or 
  + [Download tqdm](https://github.com/tqdm/tqdm) and `cd tqdm; python setup.py install --user`
  + Alternatively you can comment out the lines with tqdm and comment in the corresponding lines for execution without it

The outputs of this script are:
- CSV spreadsheet files:
  + checksOT.csv: Overtraining checks for each of the trained BDT parameters
  + viableOT.csv: Overtraining checks for only the BDTs which pass the overtraining check
  + checksSig.csv: Significance checks for BDTs which pass overtraining checks
- Plots for overtraining checks and efficiency/sensitivity checks
  +  Saved into the ./bdtDir/plots directory

Note that significance here is calculated as S/sqrt(S+B)

**calcCLs.py**   
This script takes `checksSig.csv` from `makeOTPlots.py` as input, replaces the columns of significance with CLs, and outputs to `checksCLs.csv`

The CSVs can be imported into Excel or Google Spreadsheets to sort the columns based on mass splitting/point and sensitivity to find the best settings for each splitting/mass point. Suggested sorting order:
- dm
- m(C1)
- CLs

The script `makeCLsHist.py` is also provided to plot the maximum/minimum CLs across the mass parameter space.

### Execution

```sh
# Set up root first
setupATLAS
lsetup root
ln -s CutOpt/calcCLs.py CutOpt/makeOTPlots.py . # Optional

# Execution
./makeOTPlots.py ./bdtDir # Here the directory bdtDir* contains the root files and weights from BDT training

# To execute over all the directories output by runAll.sh
./makeOTPlots.py ./Output_*

./calcCLs.py

```

