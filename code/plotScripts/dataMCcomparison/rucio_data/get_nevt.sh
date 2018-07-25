#!/bin/bash
#filelist='MCBGZjetsSherpa_sample_list.txt MCBGWjetsSherpa_sample_list.txt MCBGDYSherpa_sample_list.txt MCBGmultitop_sample_list.txt MCBGmultitop_fast_sample_list.txt MCBGVVVSherpa_sample_list.txt MCBGsingletop_sample_list.txt MCBGttbar_sample_list.txt MCBGttV_sample_list.txt MCBGVVSherpa_sample_list.txt MCBGVgammaSherpa_sample_list.txt MCBGhiggs_sample_list.txt'
filelist='MCSig_sample_list.txt'
for filename in $filelist;
do
    filepath=../../../multiLepSearch/script/$filename
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
        fi

        rucio list-files mc15_13TeV:${sample} |tee ${sample}.txt
    done
done
