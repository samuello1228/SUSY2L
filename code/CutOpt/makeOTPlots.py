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

#### Load MC sum-of-weights dictionary ####
sumWdict = {}
sumWdict.update( ssUtil.getSumWfromNTUPList( sigFilesTxt ) )
sumWdict.update( ssUtil.getSumWfromNTUPList( MCbkgFilesTxt ) )
print "sumWdict loaded\n"


#### Load nSig and nBkg. Normalize later. #### 
def getN(fileDir): 
	# import ssUtil
	from multiLepSearch.ssUtil import getCut

	# Initialize nDict
	nDict = {}
	weight = "(isMC? (ElSF * MuSF * BtagSF * weight * pwt) : (fLwt+qFwt))"
	can = ROOT.TCanvas() # for Draw(). Not necessary, but gets rid of some out of place INFO message

	fileList = listdir(fileDir)
	nFiles = len(fileList)
	n = 0
	sw = ROOT.TStopwatch(); sw.Start() # Stopwatch
	for line in tqdm(fileList):
		n += 1
		if testRun and n>1: break
		# print "Folder %d of %d: %s" % (n, nFiles, line)
		tqdm.write("Folder %d of %d: %s" % (n, nFiles, line))

		# Skip specified directories
		if any(substr in line for substr in skipDirs):
			continue

		# Get sampleID for MC, set sampleID of data to 0
		if "physics" not in line:
			match = re.search(".[0-9]{6}.", line) # Find dataset ID
			if not match:
				# print "Cannot infer datasetID from filename %s , skipped" % line
				tqdm.write("Cannot infer datasetID from filename %s , skipped" % line)
				continue

			sampleID = int(match.group()[1:-1])
		else:
			sampleID = 0 # data

		tc = TChain("evt2l"); Add = tc.Add
		for f in listdir("%s/%s" % (fileDir, line)):
			if f.endswith(".root"):
				tc.Add("%s/%s/%s" % (fileDir,line,f))

		# The next two lines are some attempt to speed up the for loop
		GetEntries=tc.GetEntries
		GetHist = gDirectory.Get

		nSample = {}

		# for chan in channels #  If tqdm not available, comment out the next line and comment in the next line
		for chan in tqdm(channels):
			tc.Draw("%s>>hist" % weight, "%s && %s" % (getCut(chan), weight))
			nSample[chan] = GetHist("hist").GetSumOfWeights()

		if nDict.get(sampleID,None) is None: # Save the dictionary if it doesn't exist for the given sampleID
			nDict[sampleID]=nSample
		else: 
			for chan in channels: 			 # Sum the entries in the dictionary if it already exists
				nDict[sampleID][chan] += nSample[chan]
			
	sw.Stop(); sw.Print() # Print stopwatch
	return nDict

# Load nBkg and nSig. Normalize later.
# Cannot be weighted during loading time because of some conflict with sumWdict and xsecDB
nBkgDict = getN(bkgDir)
print "nBkgDict loaded\n"

nSigDict = getN(sigDir)
print "nSigDict loaded\n"

#### Load cross section database #### 
ROOT.gROOT.Macro("$ROOTCOREDIR/scripts/load_packages.C")
from ROOT.SUSY import CrossSectionDB
xsecDB = CrossSectionDB("%s/data/SUSYTools/mc15_13TeV/"%os.environ["ROOTCOREBIN"] )
print "xsecDB loaded\n"

#### Normalize nSig and nBkg to luminosity and save total nBkg to nBkgTotdict ####
nBkgTotDict = {}
# Normalize nBkg and nSig and save total nBkg to nBkgTotDict 
def weightNDict():
	print "Weighting nBkgDict..."

	# for sampleID in nBkgDict:
	for sampleID in tqdm(nBkgDict):
		treeWeight = 1
		if sampleID!=0:
			mcSumW = sumWdict.get(sampleID, -1) # Get sum of weights
			xSECxEff = -1.
			xSECxEff = xsecDB.xsectTimesEff(sampleID) # Get xSec * filterEff

			treeWeight = xSECxEff * luminosity / mcSumW # Weight this tree
			if treeWeight<=0 : 
				# print "Encounter <=0 weight sample %d , skipped" % sampleID
				tqdm.write("Encounter <=0 weight sample %d , skipped" % sampleID)
				continue
		for chan in channels:
			nBkgDict[sampleID][chan] *= treeWeight
			nBkgTotDict[chan] = nBkgTotDict.get(chan, 0) + nBkgDict[sampleID][chan]
	print "nBkgDict weighted, nBkgTotDict loaded"

	print "Weighting nSigDict..."
	# for sampleID in nSigDict:
	for sampleID in tqdm(nSigDict):
		if sampleID==0: continue
		mcSumW = sumWdict.get(sampleID, -1) # Get sum of weights
		xSECxEff = -1.
		xSECxEff = xsecDB.xsectTimesEff(sampleID, 125) # Get xSec * filterEff

		treeWeight = xSECxEff * luminosity / mcSumW # Weight this tree
		if treeWeight<=0 : 
			# print "Encounter <=0 weight sample %d , skipped" % sampleID
			tqdm.write("Encounter <=0 weight sample %d , skipped" % sampleID)
			continue
		for chan in channels:
			nSigDict[sampleID][chan] = nSigDict[sampleID][chan] * treeWeight
	print "nSigDict weighted\n"

weightNDict()

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

sampleIdDB = makeSampleIdDB()

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
