#!/bin/bash

channels="0" #1 2 10 11 12"
masses="20" # 50 100 +"
runJobs(){
	for chan in $channels
	do
		for mass in $masses
		do
			python CutOpt/optWithTMVA.py -b CutOpt/GabrielFiles/bkgFiles.txt -s CutOpt/GabrielFiles/sigFiles_$mass.txt -d CutOpt/GabrielFiles/dataFiles.txt -c $chan -o Channel_$chan_dm$mass.root;
		done
	done
}

runJobs; echo "HEEYYYYYY All jobs finished!"
