#!/bin/bash

# channels="1 2 3 4 10 11 12 13 14 20 21 22 23 24"
# channels="1 2 3 4 10 11 12 13 14"
# masses="20 50 100"
# masses1="all 20"
# masses2="50 100 +"
nTreesList="100 200 400 600"
nodeSizeList="5 7 10"
depthList="2 3 4"


channels="3 13"
masses="20 50"

bkgFiles="CutOpt/GabrielFiles/MCbkgFiles.txt"
dataFiles="CutOpt/GabrielFiles/data2016.txt"

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
		    outFile="dm"$mass"_Channel"$chan1
		    sigFiles="CutOpt/GabrielFiles/sigFiles_"$mass".txt"
		    python CutOpt/optWithTMVA.py -b $bkgFiles -s $sigFiles -d $dataFiles -c $chan1 -o $outFile".root" --NTrees $nTrees --NodeSize $nodeSize --Depth $depth > $outFile".log"&
		    echo "Submitted jobs for channel "$chan1" dm="$mass" nTrees="$nTrees" nodeSize="$nodeSize" depth="$depth
		done
		# wait
	done
	wait

	# # Channel 0 
	# for mass in $masses
	# do 
	# 	outFile="dm"$mass"_Channel0"
	# 	sigFiles="CutOpt/GabrielFiles/sigFiles_"$mass".txt"
	# 	python CutOpt/optWithTMVA.py -b $bkgFiles -s $sigFiles -d $dataFiles -c 0 -o $outFile".root" --NTrees $nTrees --NodeSize $nodeSize --Depth $depth > $outFile".log"&
	# 	echo "Submitted jobs for channel 0 dm="$mass" nTrees="$nTrees" nodeSize="$nodeSize" depth="$depth
	# done 

	# wait

	dirName="Output_"$nTrees"_NodeSize"$nodeSize"_Depth"$depth
	mkdir $dirName
	mv dm*.log dm*.root $dirName
    mv weights $dirName
done; done; done

echo "All done!"
