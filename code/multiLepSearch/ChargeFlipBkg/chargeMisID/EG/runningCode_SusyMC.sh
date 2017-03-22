#!/bin/sh

listFiles(){
	inFiles=`find $1 | grep -E ".*\.myOutput.root.?[0-9]?$"`
	allFiles=""
	for file in $inFiles
	do
		allFiles=$allFiles,$file
	done
	echo $allFiles
}


echo "Did you check the data/MC tag in flip_rates.cxx? ('y' for yes)"
read input
if [ $input != "y" ] 
then
	echo "Please check it then"
	exit
fi

g++ -O3 -Wall -Wextra -std=c++11 -o flip_rates flip_rates.cxx `root-config --cflags --glibs`
if [ $? -eq 0 ]; then
    echo "Compile successful"
else
    echo "Compile failed"
    exit
fi

./flip_rates ZeeCandidate binsMC.txt `listFiles ~/work/NTuples/convertedZeePowheg`
# ./flip_rates ZeeCandidate binsMC.txt `catFiles MCfiles.txt`
if [ $? -eq 0 ]; then
	echo "flip_rates successful"
else
    echo "flip_rates failed"
    exit
fi
python llh.py binsMC.txt ratesMC.txt

