#!/bin/bash

echo "Channel,ISR,Flavor,dm,NTrees,NodeSize,Depth,SigKS,BkgKS,SigChi2,BkgChi2" > checksOT.csv

for dir in $@
do
	if [ ! -d $dir ];then
		continue
	fi
	files=`ls $dir| grep "root$"`

	echo $dir > $dir/cutInfo.txt
	for file in $files
	do
	    echo "## Now processing $dir/$file ##"

	    ISR="";Flavor="";dm="";NTrees="";NodeSize="";Depth="";

	   	dm=`echo $file | grep -o -E "dm[0-9+al]*" | sed 's/[^0-9+al]*//g'`
	    NodeSize=`echo $dir | grep -o -E "NodeSize[0-9]+" | grep -o -E [0-9]+`
	    NTrees=`echo $dir | grep -o -E "_[0-9]*_" | sed 's/[^0-9]*//g'`
	    channel=`echo $file | grep -o -E "Channel[0-9]*" | sed 's/[^0-9]*//g'`
	    Depth=`echo $dir | grep -o -E "Depth[0-9]" | sed 's/[^0-9]*//g'`

	    if [ $((($channel%100)/10)) -eq "1" ]; then
	    	ISR="yes"
	    elif [ $((($channel%100)/10)) -eq "0" ]; then
	    	ISR="no"
	    elif [ $((($channel%100)/10)) -eq "2" ]; then
	    	ISR="comb"
        fi

	    if [ $((channel%10)) -eq "0" ]; then
	    	Flavor="ee"
	    elif [ $((channel%10)) -eq "1" ]; then
	    	Flavor="eu"
	    elif [ $((channel%10)) -eq "2" ]; then
	    	Flavor="uu"
	    elif [ $((channel%10)) -eq "3" ]; then
	    	Flavor="comb"
	    elif [ $((channel%10)) -eq "4" ]; then
	    	Flavor="ll"
        fi

	    # if [ "$dm" != "+" ]; then
	    # 	continue
	    # fi

	   	# echo "dm="$dm
	    # echo "NodeSize="$NodeSize
	    # echo "NTrees="$NTrees
	    # echo "Channel="$channel
	    # echo "ISR="$ISR", Flavor="$Flavor

        # Compares overtraining 
	    ./overtraining.py $dir/$file >> $dir/cutInfo.txt
        
	   	printf "$channel,$ISR,$Flavor,$dm,$NTrees,$NodeSize,$Depth," >> checksOT.csv
	   	printf `cat trainingtest.csv`"\n" >> checksOT.csv
	    # echo >> checksOT.csv

	    name=${file%%.root}
	    if [ ! -d "$dir/plots" ];then
	    	mkdir $dir/plots
	    fi
	    mv plots/overtrain_BDTD.eps $dir/plots/overtrain_$name.eps 
	    touch $dir/plots/overtrain_$name.eps
	done
done
