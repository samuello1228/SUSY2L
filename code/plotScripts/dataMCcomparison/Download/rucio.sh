#!/bin/bash
tag=v7.0

#For Data 
for i in `cat Data_sample_list.txt`; do
j=${i##*data16_13TeV.}

#For MC
#for i in `cat MCBG_sample_list.txt`; do
#for i in `cat MCSig_sample_list.txt`; do
#j=${i##*mc15_13TeV.}

k=${tag}.${j%%.merge.*}
echo $k

#rucio ls --filter type=container  user.clo.${k}_myOutput.root/
rucio download user.clo.${k}_myOutput.root/

done
