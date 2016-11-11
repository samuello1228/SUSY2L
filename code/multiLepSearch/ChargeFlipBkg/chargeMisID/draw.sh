#!/bin/bash

today=0926		# Date/time label

ml=80			# mZ lower limit
mr=100			# mZ upper limit
sl=20			# Left sideband
sr=20			# Right sideband

side=$sl$sr		# Sideband size

dataFiles="../common/inFileList-data.txt"
# mcFiles="inFileList-MC.txt"
mcFilesDR="../common/inFileList-Zee.txt"
# mcFilesTL="inFileList-MCTL.txt"
ttbarFiles="../common/inFileList-MCttbar.txt"

echo "Preparing to produce plots for chargeMisID."
echo "Please check the following parameters before proceeding:"
echo " > Tag: "$today 2>&1 | tee log.txt
echo " > Mass window: "$ml" "$mr 2>&1 | tee -a log.txt
echo " > Sideband size: "$sl" "$sr 2>&1 | tee -a log.txt
read -p "Is this correct? Type 'y' if yes, anything else otherwise: " yn
if [ $yn != 'y' ] 
then 
	echo "Execution canceled"
	exit 1
fi
echo 2>&1 | tee -a log.txt
echo "Producing plots for chargeMisID."
echo 2>&1 | tee -a log.txt

###############
# Label folders
dataLooseNoSub=$today"-data-loose-"$ml$mr"-noSub"
dataSignalNoSub=$today"-data-signal-"$ml$mr"-noSub"
dataLooseSub=$today"-data-loose-"$ml$mr"-Sub-"$sl$sr
dataSignalSub=$today"-data-signal-"$ml$mr"-Sub-"$sl$sr

mcLooseDR=$today"-MC-loose-"$ml$mr
mcSignalDR=$today"-MC-signal-"$ml$mr

ttbarLooseDR=$today"-ttbar-loose-"$ml$mr
ttbarSignalDR=$today"-ttbar-signal-"$ml$mr



###############
# Make plots for each case
##############

# echo "############### Data, Loose, noSub. ###############" 2>&1 | tee -a log.txt
# root -l -b -q "chargeMisID.C+(\"$dataLooseNoSub\", \"$dataFiles\", \"ml=$ml,mr=$mr,sb=00,el=l,mc=n\")" 2>&1 | tee -a log.txt
# echo

echo "############### Data, Signal, noSub. ###############" 2>&1 | tee -a log.txt
root -l -b -q "chargeMisID.C+(\"$dataSignalNoSub\", \"$dataFiles\", \"ml=$ml,mr=$mr,sb=00,el=s,mc=n\")" 2>&1 | tee -a log.txt
echo

# echo "############### Data, Loose, Sub. ###############" 2>&1 | tee -a log.txt
# root -l -b -q "chargeMisID.C+(\"$dataLooseSub\", \"$dataFiles\", \"ml=$ml,mr=$mr,sl=$sl,sr=$sr,el=l,mc=n\")" 2>&1 | tee -a log.txt
# echo

# echo "############### Data, Signal, Sub. ###############" 2>&1 | tee -a log.txt
# root -l -b -q "chargeMisID.C+(\"$dataSignalSub\", \"$dataFiles\", \"ml=$ml,mr=$mr,sl=$sl,sr=$sr,el=s,mc=n\")" 2>&1 | tee -a log.txt
# echo

# echo "############### MC, Loose. (DR) ###############" 2>&1 | tee -a log.txt
# root -l -b -q "chargeMisID.C+(\"$mcLooseDR\", \"$mcFilesDR\", \"ml=$ml,mr=$mr,sb=00,el=l,mc=y,pt=y\")" 2>&1 | tee -a log.txt
# echo "Make dPT histograms" 2>&1 | tee -a log.txt
# cd $mcLooseDR
# ../drawDPThistos.py dPThistos.root -b 2>&1 | tee -a log.txt
# cd ..
# echo

echo "############### MC, Signal. (DR) ###############" 2>&1 | tee -a log.txt
root -l -b -q "chargeMisID.C+(\"$mcSignalDR\", \"$mcFilesDR\", \"ml=$ml,mr=$mr,sb=00,el=s,mc=y,pt=y,tt=n\")" 2>&1 | tee -a log.txt
echo "Make dPT histograms" 2>&1 | tee -a log.txt
cd $mcSignalDR
../drawDPThistos.py dPThistos.root -b 2>&1 | tee -a log.txt
cd ..
echo
echo 

# echo "############### MC, Loose. (ttbar) ###############" 2>&1 | tee -a log.txt
# root -l -b -q "chargeMisID.C+(\"$ttbarLooseDR\", \"$ttbarFiles\", \"mc=y,el=l,tt=y,pt=y\")" 2>&1 | tee -a log.txt
# echo "Make dPT histograms" 2>&1 | tee -a log.txt
# cd $ttbarLooseDR
# ../drawDPThistos.py dPThistos.root -b 2>&1 | tee -a log.txt
# cd ..
# echo

# echo "############### MC, Signal. (ttbar) ###############" 2>&1 | tee -a log.txt
# root -l -b -q "chargeMisID.C+(\"$ttbarSignalDR\", \"$ttbarFiles\", \"mc=y,el=s,tt=y,pt=y\")" 2>&1 | tee -a log.txt
# echo "Make dPT histograms" 2>&1 | tee -a log.txt
# cd $ttbarSignalDR
# ../drawDPThistos.py dPThistos.root -b 2>&1 | tee -a log.txt
# cd ..
# echo
# echo 

############## 
# Make comparison plots
##############
# echo "############### Comparison plots for signal (DR) ###############" 2>&1 | tee -a log.txt
# comparisonSignalDR=$today"-comparison-signal-DR"
# root -l -b -q "compareDataMC.C(\"$dataSignalNoSub/output.root\", \"$dataSignalSub/output.root\", \"0\", \"$comparisonSignalDR\")" 2>&1 | tee -a log.txt
# echo

# echo "############### Comparison plots for signal (DR) ###############" 2>&1 | tee -a log.txt
# comparisonSignalDR=$today"-comparison-signal-DR"
# root -l -b -q "compareDataMC.C(\"$dataSignalNoSub/output.root\", \"$dataSignalSub/output.root\", \"$mcSignalDR/output.root\", \"$comparisonSignalDR\")" 2>&1 | tee -a log.txt
# echo

# echo "############### Comparison plots for looseBaseline (DR) ###############" 2>&1 | tee -a log.txt
# comparisonLooseDR=$today"-comparison-loose-DR"
# root -l -b -q "compareDataMC.C(\"$dataLooseNoSub/output.root\", \"$dataLooseSub/output.root\", \"$mcLooseDR/output.root\", \"$comparisonLooseDR\")" 2>&1 | tee -a log.txt
# echo



############
# TruthLink
###########

# mcLooseTL=$today"-MC-loose-TL-"$ml$mr
# mcSignalTL=$today"-MC-signal-TL-"$ml$mr

# echo "############### MC, Loose. (TL) ###############" 2>&1 | tee -a log.txt
# root -l -b -q "chargeMisID.C+(\"$mcLooseTL\", \"$mcFilesTL\", \"ml=$ml,mr=$mr,sb=00,el=l,mc=y\")" 2>&1 | tee -a log.txt
# echo

# echo "############### MC, Signal. (TL) ###############" 2>&1 | tee -a log.txt
# root -l -b -q "chargeMisID.C+(\"$mcSignalTL\", \"$mcFilesTL\", \"ml=$ml,mr=$mr,sb=00,el=s,mc=y\")" 2>&1 | tee -a log.txt
# echo 

# echo "############### Comparison plots for signal (TL) ###############" 2>&1 | tee -a log.txt
# comparisonSignalTL=$today"-comparison-signal-TL"
# root -l -b -q "compareDataMC.C+(\"$dataSignalNoSub/output.root\", \"$dataSignalSub/output.root\", \"$mcSignalTL/output.root\", \"$comparisonSignalTL\")" 2>&1 | tee -a log.txt
# echo

# echo "############### Comparison plots for looseBaseline (TL) ###############" 2>&1 | tee -a log.txt
# comparisonLooseTL=$today"-comparison-loose-TL"
# root -l -b -q "compareDataMC.C+(\"$dataLooseNoSub/output.root\", \"$dataLooseSub/output.root\", \"$mcLooseTL/output.root\", \"$comparisonLooseTL\")" 2>&1 | tee -a log.txt
# echo


############
# MCTruthClassifier
#############
# mcLooseMCTC=$today"-MC-loose-MCTC-"$ml$mr
# mcSignalMCTC=$today"-MC-signal-MCTC-"$ml$mr

# echo "############### MC, Loose. (MCTC) ###############" 2>&1 | tee -a log.txt
# root -l -b -q "chargeMisID.C+(\"$mcLooseMCTC\", \"$mcFiles\", \"ml=$ml,mr=$mr,sb=00,el=l,mc=y\")" 2>&1 | tee -a log.txt
# echo

# echo "############### MC, Signal. (MCTC) ###############" 2>&1 | tee -a log.txt
# root -l -b -q "chargeMisID.C+(\"$mcSignalMCTC\", \"$mcFiles\", \"ml=$ml,mr=$mr,sb=00,el=s,mc=y\")" 2>&1 | tee -a log.txt
# echo 

# echo "############### Comparison plots for signal (MCTC) ###############" 2>&1 | tee -a log.txt
# comparisonSignal=$today"-comparison-signal-MCTC"
# root -l -b -q "compareDataMC.C+(\"$dataSignalNoSub/output.root\", \"$dataSignalSub/output.root\", \"$mcSignalMCTC/output.root\", \"$comparisonSignal\")" 2>&1 | tee -a log.txt
# echo

# echo "############### Comparison plots for looseBaseline (MCTC) ###############" 2>&1 | tee -a log.txt
# comparisonLoose=$today"-comparison-loose-MCTC"
# root -l -b -q "compareDataMC.C+(\"$dataLooseNoSub/output.root\", \"$dataLooseSub/output.root\", \"$mcLooseMCTC/output.root\", \"$comparisonLoose\")" 2>&1 | tee -a log.txt


