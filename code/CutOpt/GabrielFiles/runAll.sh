#!/bin/bash

channels="0 1 2 10 11 12"
#masses="20 50 100 +"

bkgFiles="CutOpt/GabrielFiles/bkgFiles.txt"
sigFiles="CutOpt/GabrielFiles/sigFiles_"$1".txt"
dataFiles="CutOpt/GabrielFiles/dataFiles.txt"
runJobs(){
	for chan in $channels
	do
        outFile="Channel"$chan"_dm"$1".root"
    	python CutOpt/optWithTMVA.py -b $bkgFiles -s $sigFiles -d $dataFiles -c $chan -o $outFile&
	done
}

runJobs; echo "HEEYYYYYY All jobs finished!"
