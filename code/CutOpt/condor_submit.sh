#!/bin/bash

channels="0 1 2 3 4 10 11 12 13 14"
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

	dirName="Output_"$nTrees"_NodeSize"$nodeSize"_Depth"$depth
	if [ -d "$dirName" ]; then 
		rm -rf $dirName
	fi
	mkdir $dirName
	# cd $dirName

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
		    condor_submit MASS=$mass CHAN=$chan1 BKG=$bkgFiles DATA=$dataFiles NTREES=$nTrees NODESIZE=$nodeSize DEPTH=$depth condor.sub
		done
	done

done; done; done

