#!/bin/bash
filepath=../../../multiLepSearch/script/Data_sample_list.txt
echo opening $filepath

for sample in $(cat $filepath);
do
    echo $sample
    md5="$(md5sum ${sample}.txt)"
    isEmpty=0

    #check whether the output file is empty
    for element in $md5;
    do
        if [ "$element" = "d41d8cd98f00b204e9800998ecf8427e" ];
        then
            echo $element
            isEmpty=1
            echo "empty"
        fi
    done

    if [ $isEmpty -ne 1 ];
    then
        #continue
        x=1
    fi

    rucio list-files ${sample} |tee ${sample}.txt
done
