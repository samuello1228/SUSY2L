#!/bin/bash

channels="0 1 2 10 11 12"
masses="20 50 100 +"
runJobs(){
	for chan in $channels
	do
		for mass in $masses
		{python CutOpt/optWithTMVA.py -b GabrielFiles/bkgFiles.txt -s GabrielFiles/sigFiles_$mass.txt -d GabrielFiles/dataFiles.txt -c $channel -o Channel_$chan_dm$mass.root};
	done
}

runJobs; echo "HEEYYYYYY All jobs finished!"
