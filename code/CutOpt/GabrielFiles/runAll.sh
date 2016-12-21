#!/bin/bash

channels="1 2 3 4 10 11 12 13 14 20 21 22 23 24"
masses="20 50 100"
# masses1="all 20"
# masses2="50 100 +"
nTreesList="100 200 400 600 800 1000"
nodeSizeList="5 7 10"

bkgFiles="CutOpt/GabrielFiles/bkgFiles.txt"
dataFiles="CutOpt/GabrielFiles/dataFiles.txt"

for nTrees in $nTreesList; do for nodeSize in $nodeSizeList
do
	for chan in $channels
	do
		for mass in $masses
		do
			if [ $mass -eq "100" ]
			then	
				chan1=`expr $chan + 100`
			else 
				chan1=$chan
			fi
		    outFile="Channel"$chan1"_dm"$mass".root"
		    sigFiles="CutOpt/GabrielFiles/sigFiles_"$mass".txt"
		    python CutOpt/optWithTMVA.py -b $bkgFiles -s $sigFiles -d $dataFiles -c $chan1 -o $outFile --NTrees $nTrees --NodeSize $nodeSize&
		    echo "Submitted jobs for channel "$chan1" dm="$mass" nTrees="$nTrees" nodeSize="$nodeSize
		done
		wait
	done

    # For executing on selected channels only
	#dirName="Output_"$nTrees"_NodeSize"$nodeSize
	#mv Channel*.root $dirName
    #continue

	# Channel 0 
	outFile="Channel0_dm20.root"
	sigFiles="CutOpt/GabrielFiles/sigFiles_20.txt"
	python CutOpt/optWithTMVA.py -b $bkgFiles -s $sigFiles -d $dataFiles -c 0 -o $outFile --NTrees $nTrees --NodeSize $nodeSize&
	echo "Submitted jobs for channel 0 dm=20 nTrees="$nTrees" nodeSize="$nodeSize

	outFile="Channel0_dm50.root"
	sigFiles="CutOpt/GabrielFiles/sigFiles_50.txt"
	python CutOpt/optWithTMVA.py -b $bkgFiles -s $sigFiles -d $dataFiles -c 0 -o $outFile --NTrees $nTrees --NodeSize $nodeSize&
	echo "Submitted jobs for channel 0 dm=50 nTrees="$nTrees" nodeSize="$nodeSize

	outFile="Channel0_dm100.root"
	sigFiles="CutOpt/GabrielFiles/sigFiles_100.txt"
	python CutOpt/optWithTMVA.py -b $bkgFiles -s $sigFiles -d $dataFiles -c 100 -o $outFile --NTrees $nTrees --NodeSize $nodeSize&
	echo "Submitted jobs for channel 0 dm=100 nTrees="$nTrees" nodeSize="$nodeSize

	dirName="Output_"$nTrees"_NodeSize"$nodeSize
	mkdir $dirName
	mv Channel*.root $dirName
    mv weights $dirName
done; done

echo "All done!"
