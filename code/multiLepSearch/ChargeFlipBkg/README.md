# README: ChargeFlipBkg 
Last update: 28 Jul 2016

This code package
- measures electron charge misidentification rates, and 
- applies the rates to predict the number of events in the same-sign dilepton channel arising from opposite sign dilepton events where an electron has had its charge misidentified.

All code is ready for 20.7 derivations.

Two separate programs for measuring chargeMisID are included in this package. See the "[About the chargeMisID code](#about-the-chargemisid-code)" section below

## Contents
- [Documentation of sub-packages](#documentation-of-subpackages)
- [The workflow](#the-workflow)
    + [An example workflow](#an-example-workflow)
- [About the chargeMisID code](#about-the-chargemisid-code)

## Documentation of sub-packages
- [NTuple conversion](convert/README.md)
- [In house chargeMisID](chargeMisID/README.md)
- [EGamma chargeMisID](chargeMisID/EG/README)
- [pt correction](ptcorr/README.md)
- [SS estimation](estimate/README.md)
- [ChargeFlipTool](ChargeFlipTool/README.md)

## The workflow 
1. NTuples are made using the UM/HK SUSY code package `multipLepSelection/ssEvtSelection`
2. Convert NTuples into a format acceptable for the chargeMisID and pt correction tools
3. Find the chargeMisID rates 
    - EGamma code depends on converted NTuples
    - In-house code depends on original NTuples
4. Find the pt correction necessary to correct the predicted SS distributions to match the observed SS distributions
    - Depends on converted NTuples
5. Save rates and pt correction values to a "friend" NTuple of those made in step 1 (see note below)
6. Loop over events in the NTuples made in steps 1&5 to plot histograms

- *Note*: As an alternative to steps 5-6, one can also make histograms directly without making the "friend" NTuple

### An example workflow
```sh
# This is just a sequential list of commands one might use to execute this workflow. 
# It is up to the user to make sure every script/program listed here is configured properly before execution. 
# This section is *not* meant to be copied verbatim.

# STEP 1b: Get the list of NTuples and save it to common/inFileList.txt (with full path) and convert/inFileList.txt (filenames only)
ls /afs/cern.ch/user/g/ggallard/work/NTuples/user.clo.v7.1.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_myOutput.root > convert/inFileList.txt

# STEP 2: Convert the NTuples to a format acceptable by the EGamma chargeMisID tool
cd ChargeFlipBkg/convert
make mcnew

# STEP 3: Find chargeMisID rates (use either EGamma or in-house code)

## EGamma code
cd ChargeFlipBkg/chargeMisID/EG/
./runningCode_SusyMC.sh             # For MC
./runningCode_SusyDATA.sh           # For data
./moveFiles.sh <someDirectory>      # Convenience script to move relevant output files to some location

## In-house code
cd ChargeFlipBkg/chargeMisID
root -l -b -q chargeMisID.C+(...)   # See chargeMisID/README.md for documentation about arguments

## Make plots if necessary
cd ChargeFlipBkg/chargeMisID
root -l -b -q compareDataMC.C(...)  # See chargeMisID/README.md for documentation about arguments

# STEP 4: Find pt corrections
cd ChargeFlipBkg/ptcorr
root -l -b -q makeDEhistos.C(...)   # See ptcorr/README.md for documentation about arguments
./drawDEhistos.py                   # Make plots and saves 2D profiles of 3D histograms to the output root file

# STEP 5: Save misID rates and pt corrections to "friend" NTuple
cd ChargeFlipBkg/estimate
root -l -b -q saveSSprediction.C+(...)      # See estimate/README.md for documentation

# STEP 6: Make histograms from original NTuples and "friend" NTuple
root -l -b -q extractSSprediction.C+(...)   # See estimate/README.md for documentation

# STEP 5 (alternative): Direct to histograms
root -l -b -q getSSPrediction.C(...)        # See estimate/README.md for documentation
 
```

## About the chargeMisID code

### EGamma's chargeMisID code  
*Found in chargeMisID/EG/*
- The chargeMisID code adopted was supplied by the EGamma group. 
- It can be found at `svn+ssh://svn.cern.ch/reps/atlasperf/CombPerf/EGamma/TagAndProbe/Run2/ChargeMisId`.
- The version in this package at the time of writing was based on the trunk which was checked out 16 June 2016 10:44am CERN time. It was then updated and modified to accept a converted version of the UM/HK SUSY NTuples

### In-house chargeMisID code
*Found in chargeMisID/*  
- [chargeMisID.C](chargeMisID/chargeMisID.C) was developed in-house independently from the EGamma code. 
    + chargeMisID.C has the advantage of being (much) faster, (much) easier to automate, and it runs directly on the UM/HK SUSY NTuples, thus eliminating the intermediate step of NTuples conversion
    + The results obtained by this program and EGamma's program are nearly identical
    + The only caveat is that sideband background subtraction, a necessary feature for running on data, is *not* well-implemented.
