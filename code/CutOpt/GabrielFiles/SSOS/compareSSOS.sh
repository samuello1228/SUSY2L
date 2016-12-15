#!/bin/bash

channels="1 2 3 10 11 12 13 20 21 22 23"
masses="all 20 50 100 +"

for mass in $masses
do
	for chan in $channels
	do
		newdir="dm"$mass"_channel"$chan
		mkdir $newdir
		cd $newdir
	    sigFiles="../SUSY2L/code/CutOpt/GabrielFiles/sigFiles_"$mass".txt"
	    root -l -b -q "plot.cxx(\"$sigFiles\", $chan)"
	    cd ..
	done

	newdir="dm"$mass"_channel0"
	mkdir newdir
	cd newdir
    sigFiles="CutOpt/GabrielFiles/sigFiles_"$mass".txt"
    root -l -b -q "plot.cxx(\"$sigFiles\", 0)"
    cd ..

done