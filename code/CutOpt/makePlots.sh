#!/bin/bash

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

        # Compares overtraining 
	    ./overtraining.py $dir/$file >> $dir/cutInfo.txt
        
        # Draws signal/background efficiency curves
	    root -l -b -q "mvaeffs.cxx(\"$dir/$file\")" >> $dir/cutInfo.txt
	    name=${file%%.root}
	    if [ ! -d "$dir/plots" ];then
	    	mkdir $dir/plots
	    fi
	    mv plots/overtrain_BDTD.eps $dir/plots/overtrain_$name.eps 
	    #mv plots/overtrain_BDTD.png $dir/plots/overtrain_$name.png
	    mv plots/mvaeffs_BDTD.eps $dir/plots/mvaeffs_$name.eps
	done
done
