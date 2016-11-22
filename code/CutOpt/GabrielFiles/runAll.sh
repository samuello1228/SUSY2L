#!/bin/bash

channels="1 2 10 11 12"
masses="20 50 100 +"
nTreesList="400 600 1000"
nodeSizeList"5 7 10"

bkgFiles="CutOpt/GabrielFiles/bkgFiles.txt"
dataFiles="CutOpt/GabrielFiles/dataFiles.txt"

for nTrees in $nTreesList; do; for nodeSize in $nodeSizeList
do
	for chan in $channels
	do
		for mass in $masses
		do
		    outFile="Channel"$chan"_dm"$mass".root"
		    sigFiles="CutOpt/GabrielFiles/sigFiles_"$mass".txt"
		    python CutOpt/optWithTMVA.py -b $bkgFiles -s $sigFiles -d $dataFiles -c $chan -o $outFile --NTrees $nTrees --NodeSize $nodeSize&
		    echo "Submitted jobs for channel "$chan" dm="$mass
		done
		wait
	done

	# Channel 0 
	for mass in $masses
	do
		outFile="Channel"0"_dm"$mass".root"
		sigFiles="CutOpt/GabrielFiles/sigFiles_"$mass".txt"
		python CutOpt/optWithTMVA.py -b $bkgFiles -s $sigFiles -d $dataFiles -c 0 -o $outFile --NTrees $nTrees --NodeSize $nodeSize&
	done
	wait
	dirName="Output_"$nTrees"_NodeSize"$nodeSize
	mkdir $dirName
	mv Channel*.root $dirName
done; done

echo "All done!"