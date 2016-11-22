#!/bin/bash

channels="2 10 11 12"
masses="20 50 100 +"

bkgFiles="CutOpt/GabrielFiles/bkgFiles.txt"
dataFiles="CutOpt/GabrielFiles/dataFiles.txt"

for chan in $channels
do
	for mass in $masses
	do
	    outFile="Channel"$chan"_dm"$mass".root"
	    sigFiles="CutOpt/GabrielFiles/sigFiles_"$mass".txt"
	    python CutOpt/optWithTMVA.py -b $bkgFiles -s $sigFiles -d $dataFiles -c $chan -o $outFile&
	    echo "Submitted jobs for channel "$chan" dm="$mass
	done
	wait
done

for mass in $masses
do
	outFile="Channel"0"_dm"$mass".root"
	sigFiles="CutOpt/GabrielFiles/sigFiles_"$mass".txt"
	python CutOpt/optWithTMVA.py -b $bkgFiles -s $sigFiles -d $dataFiles -c 0 -o $outFile&
done

wait

echo "All done!"