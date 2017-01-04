#!/bin/bash

# channels="1 2 3 4 10 11 12 13 14 20 21 22 23 24"
channels="1 2 3 4 10 11 12 13 14"
masses="20 50 100"
# masses1="all 20"
# masses2="50 100 +"
nTreesList="100 200 400 600"
nodeSizeList="5 7 10"
depthList="2 3 4 5"

bkgFiles="CutOpt/GabrielFiles/bkgFiles.txt"
dataFiles="CutOpt/GabrielFiles/dataFiles.txt"

for nTrees in $nTreesList; do for nodeSize in $nodeSizeList; do for depth in $depthList
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
		    outFile="dm"$mass"_Channel"$chan1".root"
		    sigFiles="CutOpt/GabrielFiles/sigFiles_"$mass".txt"
		    python CutOpt/optWithTMVA.py -b $bkgFiles -s $sigFiles -d $dataFiles -c $chan1 -o $outFile --NTrees $nTrees --NodeSize $nodeSize --Depth $depth&
		    echo "Submitted jobs for channel "$chan1" dm="$mass" nTrees="$nTrees" nodeSize="$nodeSize" depth="$depth
		done
		wait
	done

	# Channel 0 
	outFile="dm20_Channel0.root"
	sigFiles="CutOpt/GabrielFiles/sigFiles_20.txt"
	python CutOpt/optWithTMVA.py -b $bkgFiles -s $sigFiles -d $dataFiles -c 0 -o $outFile --NTrees $nTrees --NodeSize $nodeSize --Depth $depth&
	echo "Submitted jobs for channel 0 dm=20 nTrees="$nTrees" nodeSize="$nodeSize" depth="$depth

	outFile="dm50_Channel0.root"
	sigFiles="CutOpt/GabrielFiles/sigFiles_50.txt"
	python CutOpt/optWithTMVA.py -b $bkgFiles -s $sigFiles -d $dataFiles -c 0 -o $outFile --NTrees $nTrees --NodeSize $nodeSize --Depth $depth&
	echo "Submitted jobs for channel 0 dm=50 nTrees="$nTrees" nodeSize="$nodeSize" depth="$depth

	outFile="dm100_Channel100.root"
	sigFiles="CutOpt/GabrielFiles/sigFiles_100.txt"
	python CutOpt/optWithTMVA.py -b $bkgFiles -s $sigFiles -d $dataFiles -c 100 -o $outFile --NTrees $nTrees --NodeSize $nodeSize --Depth $depth&
	echo "Submitted jobs for channel 0 dm=100 nTrees="$nTrees" nodeSize="$nodeSize" depth="$depth

    wait

	dirName="Output_"$nTrees"_NodeSize"$nodeSize"_Depth"$depth
	mkdir $dirName
	mv dm*.root $dirName
    mv weights $dirName
done; done; done

echo "All done!"
