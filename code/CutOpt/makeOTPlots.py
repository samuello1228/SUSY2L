#!/usr/bin/env python

import sys, os, subprocess, ROOT, csv, overtraining, re#, ssUtil
import tqdm # Progress bar # Find and comment out all lines with 'tqdm' if not installed
from os import listdir, path, mkdir, rename
from os.path import isdir, isfile
from subprocess import call
from csv import reader
from ROOT import TChain, TH1F, gDirectory, gROOT
from tqdm import tqdm

import multiLepSearch
from multiLepSearch import ssUtil

gROOT.SetBatch(True)
####################
# Config variables #
####################
C1masses = (200, 300, 400, 500, 600, 700, 800, 900, 1000)
luminosity = 33257.2
# channels = (0, 1, 2, 3, 4, 10, 11, 12, 13, 14)
channels = (0, 1, 2, 3, 4, 10, 11, 12, 13, 14, 100, 101, 102, 103, 104, 110, 111, 112, 113, 114)
sigFilesTxt = "CutOpt/GabrielFiles/allSig.txt"
MCbkgFilesTxt = "CutOpt/GabrielFiles/MCbkgFiles.txt"
dataBkgFilesTxt = "CutOpt/GabrielFiles/data2016.txt"

bkgDir = "/srv/SUSY/ntuple/AnalysisBase-02-04-26-da7031fc/"
sigDir = "/srv/SUSY/ntuple/AnalysisBase-02-04-26-4dcc2f47/"

skipDirs = ["Zee_MAX", "Zmumu", "Ztautau", "P2012"]

testRun = False
####################

def makeSampleIdDB():
	sampleIdDB = {}
	fileList = open(sigFilesTxt,"r").readlines()
	# print fileList
	for aLine in fileList:
		if "MGPy8EG" not in aLine or "root/user" not in aLine: continue
		# if (re.search("MCSig.[0-9]{6}", aLine)) is None: continue
		sampleID = int((re.search("12.[0-9]{6}", aLine).group()).split('.')[1])
		masses = re.search("Slep_[0-9]{3,4}_[0-9]{1,4}", aLine).group()
		# print masses
		mC1 = int(masses.split('_')[1])
		mN1 = int(masses.split('_')[2])
		sampleIdDB[(mC1, mN1)] = sampleID
	return sampleIdDB
	# print sampleIdDB

def getTotDict(nDict):
	nTotDict = {}

	for chan in channels:
		nTotDict[chan] = 0

	for sampleID in nDict:
		for chan in channels:
			nTotDict[chan] += nDict[sampleID][chan]

	return nTotDict
	
def getNDict(dictCSV):
	allLines = open(dictCSV).readlines()

	headers = allLines[0].split('\n')[0].split(',')
	nCols = len(headers)

	nDict = {}
	for row in allLines[1:]:
		entries = row.split('\n')[0].split(',')
		sampleDict = {int(headers[i]): float(entries[i]) for i in range(1,nCols)}
		nDict[int(entries[0])] = sampleDict
	return nDict

#### Utility functions
def whichISR(nISR):
	return {
		0 : 'no',
		1 : 'yes',
		2 : 'comb'
	}[nISR]

def whichFlav(nFlav):
	return {
		0 : 'ee',
		1 : 'em',
		2 : 'mm',
		3 : 'comb',
		4 : 'll'
	}[nFlav]

def notOvertrained(overtrainingNums):
	SigKS = overtrainingNums[0]
	BkgKS = overtrainingNums[1]
	SigChi2 = overtrainingNums[2]
	BkgChi2 = overtrainingNums[3]
	return (SigKS > 0.05 and BkgKS > 0.05 and SigKS < 0.95 and BkgKS < 0.95 ) 
			# and (SigChi2 > 0.05 and BkgChi2 > 0.05 and SigChi2 < 0.95 and BkgChi2 < 0.95)

#### Main function
if __name__ == '__main__':
	ROOT.gIgnoreErrorLevel = ROOT.kInfo # Quiet INFO messages

	# Load nSig and nBkg
	sampleIdDB = makeSampleIdDB()
	nBkgDict = getNDict("nBkgDict.csv")
	nSigDict = getNDict("nSigDict.csv")
	nBkgTotDict = getTotDict(nBkgDict)

	# Initialize output .csv's
	outOT = open("checksOT.csv", "w")
	outOT.write("Channel,ISR,Flavor,dm,NTrees,NodeSize,Depth,SigKS,BkgKS,SigChi2,BkgChi2\n")

	outViableOT = open("viableOT.csv", "w")
	outViableOT.write("Channel,ISR,Flavor,dm,NTrees,NodeSize,Depth,SigKS,BkgKS,SigChi2,BkgChi2\n")

	outSig = open("checksSig.csv", "w")
	outSig.write("m(C1),m(N1),Channel,ISR,Flavor,NTrees,NodeSize,Depth,")
	outSig.write("nSig,nBkg,BDTopt,nSig(BDTopt),nBkg(BDTopt),sigma(BDTopt),")
	outSig.write("nSig(0),nBkg(0),sigma(0),")
	outSig.write("nSig(0.1),nBkg(0.1),sigma(0.1),")
	outSig.write("nSig(0.2),nBkg(0.2),sigma(0.2),")
	outSig.write("nSig(0.3),nBkg(0.3),sigma(0.3)\n")

	for directory in tqdm(sys.argv[1:]):
		directory = directory.rstrip('/')
		if not(isdir(directory)): continue
		files = listdir(directory)

		plotDir = "%s/plots" % directory
		if not isdir(plotDir): mkdir(plotDir)

		for file in tqdm(files):
			if not(file.endswith(".root")):	continue

			# print ("## Now processing %s/%s##" % (directory, file))
			tqdm.write("## Now processing %s/%s##" % (directory, file))
			dm = int(re.search("[0-9]+", re.search("dm[0-9]+", file).group()).group())
			NodeSize = int(re.search("[0-9]+", re.search("NodeSize[0-9]+", directory).group()).group())
			NTrees = int(re.search("[0-9]+", re.search("_[0-9]+_", directory).group()).group())
			channel = int(re.search("[0-9]+", re.search("Channel[0-9]+", file).group()).group())
			Depth = int(re.search("[0-9]+", re.search("Depth[0-9]+", directory).group()).group())
			ISR = whichISR((channel%100)/10)
			Flavor = whichFlav(channel%10)

			outOT.write("%s,%s,%s,%s,%s,%s,%s," 
				% (channel, ISR, Flavor, dm, NTrees, NodeSize, Depth))
			overtrainingNums = overtraining.runCheck("%s/%s" % (directory, file))
			outOT.write("%1.3f,%1.3f,%1.3f,%1.3f\n" % overtrainingNums)

			overtrainPlotname = "%s/dm%d_Channel%d.eps" % (plotDir, dm, channel)
			if isfile("plots/overtrain_BDTD.eps"):
				os.rename("plots/overtrain_BDTD.eps", overtrainPlotname)
			else: open(overtrainPlotname, 'a').close()
			
			if not(notOvertrained(overtrainingNums)): continue

			outViableOT.write("%s,%s,%s,%s,%s,%s,%s," 
				% (channel, ISR, Flavor, dm, NTrees, NodeSize, Depth))
			outViableOT.write("%1.3f,%1.3f,%1.3f,%1.3f\n" % overtrainingNums)

			for mass in tqdm(C1masses): 
				nSig = nSigDict.get(sampleIdDB[(mass, mass-dm)],{}).get(channel,0)
				if nSig < 1:
					# print "nSig<1. Skipping."
					# tqdm.write("nSig<1. Skipping.")
					continue
				nBkg = nBkgTotDict[channel]

				# print "## mC1 = %d, mN1 = %d, chan = %d, nSig %f\n" % (mass, mass-dm, channel, nSig)
				# tqdm.write("## mC1 = %d, mN1 = %d, chan = %d, nSig %f\n" % (mass, mass-dm, channel, nSig))

				call('root -l -b -q "mvaeffs.cxx(\\"%s/%s\\", %d, %f)"' % (directory, file, int(nSig), int(nBkg)), shell=True)
				
				outSig.write("%d,%d,%d,%s,%s,%d,%d,%d," 
					% (mass, mass-dm, channel, ISR, Flavor, NTrees, NodeSize, Depth))
				with open("effs.csv") as effFile:
					outSig.write(effFile.readline())
				os.rename("plots/mvaeffs_BDTD.eps", "%s/%d_%d_Channel%d.eps" % (plotDir, mass, mass-dm, channel))

	outOT.close()
	outViableOT.close()
	outSig.close()
	print "\n\n\n\n"
