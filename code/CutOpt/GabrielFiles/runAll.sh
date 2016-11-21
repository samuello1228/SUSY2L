#!/bin/bash
python CutOpt/optWithTMVA.py -b bkgFiles.txt -s sigFiles_20.txt -d dataFiles.txt -c ISR -o BDT_ISR_dm20.root
python CutOpt/optWithTMVA.py -b bkgFiles.txt -s sigFiles_50.txt -d dataFiles.txt -c ISR -o BDT_ISR_dm50.root
python CutOpt/optWithTMVA.py -b bkgFiles.txt -s sigFiles_100.txt -d dataFiles.txt -c ISR -o BDT_ISR_dm100.root
python CutOpt/optWithTMVA.py -b bkgFiles.txt -s sigFiles_+.txt -d dataFiles.txt -c ISR -o BDT_ISR_dm+.root
python CutOpt/optWithTMVA.py -b bkgFiles.txt -s sigFiles_20.txt -d dataFiles.txt -o BDT_noISR_dm20.root
python CutOpt/optWithTMVA.py -b bkgFiles.txt -s sigFiles_50.txt -d dataFiles.txt -o BDT_noISR_dm50.root
python CutOpt/optWithTMVA.py -b bkgFiles.txt -s sigFiles_100.txt -d dataFiles.txt -o BDT_noISR_dm100.root
python CutOpt/optWithTMVA.py -b bkgFiles.txt -s sigFiles_+.txt -d dataFiles.txt -o BDT_noISR_dm+.root
