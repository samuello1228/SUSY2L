#!/bin/bash

channels="1 2 3 10 11 12 13 20 21 22 23"
masses1="all 20"
masses2="50 100 +"
nTreesList="100 200 400 600 800 1000"
nodeSizeList="5 7 10"

bkgFiles="CutOpt/GabrielFiles/bkgFiles.txt"
dataFiles="CutOpt/GabrielFiles/dataFiles.txt"

for nTrees in $nTreesList; do for nodeSize in $nodeSizeList
do
	for chan in $channels
	do
		for mass in $masses1
		do
		    outFile="Channel"$chan"_dm"$mass".root"
		    sigFiles="CutOpt/GabrielFiles/sigFiles_"$mass".txt"
		    python CutOpt/optWithTMVA.py -b $bkgFiles -s $sigFiles -d $dataFiles -c $chan -o $outFile --NTrees $nTrees --NodeSize $nodeSize&
		    echo "Submitted jobs for channel "$chan" dm="$mass" nTrees="$nTrees" nodeSize="$nodeSize
		done
		wait
		for mass in $masses2
		do
		    outFile="Channel"$chan"_dm"$mass".root"
		    sigFiles="CutOpt/GabrielFiles/sigFiles_"$mass".txt"
		    python CutOpt/optWithTMVA.py -b $bkgFiles -s $sigFiles -d $dataFiles -c $chan -o $outFile --NTrees $nTrees --NodeSize $nodeSize&
		    echo "Submitted jobs for channel "$chan" dm="$mass" nTrees="$nTrees" nodeSize="$nodeSize
		done
		wait
	done

    # For executing on selected channels only
	#dirName="Output_"$nTrees"_NodeSize"$nodeSize
	#mv Channel*.root $dirName
    #continue

	# Channel 0 
	for mass in $masses1
	do
		outFile="Channel"0"_dm"$mass".root"
		sigFiles="CutOpt/GabrielFiles/sigFiles_"$mass".txt"
		python CutOpt/optWithTMVA.py -b $bkgFiles -s $sigFiles -d $dataFiles -c 0 -o $outFile --NTrees $nTrees --NodeSize $nodeSize&
		echo "Submitted jobs for channel 0 dm="$mass" nTrees="$nTrees" nodeSize="$nodeSize
	done
	wait

	for mass in $masses2
	do
		outFile="Channel"0"_dm"$mass".root"
		sigFiles="CutOpt/GabrielFiles/sigFiles_"$mass".txt"
		python CutOpt/optWithTMVA.py -b $bkgFiles -s $sigFiles -d $dataFiles -c 0 -o $outFile --NTrees $nTrees --NodeSize $nodeSize&
		echo "Submitted jobs for channel 0 dm="$mass" nTrees="$nTrees" nodeSize="$nodeSize
	done
	wait
	dirName="Output_"$nTrees"_NodeSize"$nodeSize
	mkdir $dirName
	mv Channel*.root $dirName
    mv weights $dirName
done; done

echo "All done!"
