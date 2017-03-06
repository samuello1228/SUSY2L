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
# channels = (0, 1, 2, 3, 4, 10, 11, 12, 13, 14, 100, 101, 102, 103, 104, 110, 111, 112, 113, 114)
channels = (3, 13)

## Updated NTuples
sigFilesTxt = "CutOpt/GabrielFiles/allSig.txt"
MCbkgFilesTxt = "CutOpt/GabrielFiles/MCbkgFiles.txt"
dataBkgFilesTxt = "CutOpt/GabrielFiles/data2016.txt"

bkgDir = "/srv/SUSY/ntuple/AnalysisBase-02-04-26-da7031fc/"
sigDir = "/srv/SUSY/ntuple/AnalysisBase-02-04-26-4dcc2f47/"

## Old NTuples
bkgDir = "/srv/SUSY/ychan/allBkg/"
sigDir = "/srv/SUSY/ychan/sig/"

sigFilesTxt = "oldLists/allSig.txt"
MCbkgFilesTxt = "oldLists/MCbkgFiles.txt"
dataBkgFilesTxt = "oldLists/data10.txt"


skipDirs = ["Zee_MAX", "Zmumu", "Ztautau", "P2012_ttbar", "P2012_Wt" "ttH125", "ZZZ"]

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
	can = ROOT.TCanvas() # for Draw(). Not necessary, but gets rid of some out of place INFO message

	fileList = [d for d in listdir(fileDir) if isdir("%s/%s" % (fileDir,d)) and d.endswith(".root")]
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
			weight = "(ElSF * MuSF * BtagSF * weight * pwt)"
		else:
			sampleID = 0 # data

		tc = TChain("evt2l"); Add = tc.Add
		for f in listdir("%s/%s" % (fileDir, line)):
			if f.endswith(".root"):
				tc.Add("%s/%s/%s" % (fileDir,line,f))

		# The next two lines are some attempt to speed up the for loop
		GetEntries=tc.GetEntries
		GetHist = gDirectory.Get

		nSample = {}; nQF = {}; nCF = {}

		## This part is not very elegant, but at least it works
		if sampleID != 0: # For MC
			# for chan in channels #  If tqdm not available, comment out the next line and comment in the next line
			for chan in tqdm(channels):
				tc.Draw("%s>>hist" % weight, "%s && %s" % (getCut(chan), weight))
				nSample[chan] = GetHist("hist").GetSumOfWeights()
			if nDict.get(sampleID,None) is None: # Save the dictionary if it doesn't exist for the given sampleID
				nDict[sampleID]=nSample
			else: 
				for chan in channels: 			 # Sum the entries in the dictionary if it already exists
					nDict[sampleID][chan] += nSample[chan]
		else: # For data
			tcf = TChain("PP1_evt2l");
			for f in listdir("%s/%s" % (fileDir, line)):
				if f.endswith(".root"):
					tcf.Add("%s/%s/%s" % (fileDir,line,f))
			tc.AddFriend(tcf)
			# for chan in channels #  If tqdm not available, comment out the next line and comment in the next line
			for chan in tqdm(channels):
				# OLDNTUPLES
				weight = "PP1_evt2l.qFwt"
				tc.Draw("%s>>hist1" % weight, "%s && %s" % (getCut(chan), weight))
				nQF[chan] = GetHist("hist1").GetSumOfWeights()

				weight = "PP1_evt2l.fLwt"
				tc.Draw("%s>>hist2" % weight, "%s && %s" % (getCut(chan), weight))
				nCF[chan] = GetHist("hist2").GetSumOfWeights()

			# Assign chargeFlip events sampleID = -1,  fakeLepton sampleID = -2
			if nDict.get(-1,None) is None: # Save the dictionary if it doesn't exist for the given sampleID
				nDict[-1]=nQF
			else: 
				for chan in channels: 			 # Sum the entries in the dictionary if it already exists
					nDict[-1][chan] += nQF[chan]

			if nDict.get(-2,None) is None: # Save the dictionary if it doesn't exist for the given sampleID
				nDict[-2]=nCF
			else: 
				for chan in channels: 			 # Sum the entries in the dictionary if it already exists
					nDict[-2][chan] += nCF[chan]
			
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


#### To find source of MC background from sampleID
MCbkgList = open(MCbkgFilesTxt).readlines()
def sourceBkg(sampleID):
	if sampleID == -1: return "qFlip"
	if sampleID == -2: return "fLep"
	sampleID = "%d" % sampleID; sampleType = None
	for s in MCbkgList:
		if sampleID in s: sampleType = ssUtil.guessSampleType(s)
	return sampleID if sampleType is None else sampleType


#### Normalize nSig and nBkg to luminosity and save total nBkg to nBkgTotdict ####
nBkgTotDict = {}; nBkgSourceDict={}
# Normalize nBkg and nSig and save total nBkg to nBkgTotDict 
def weightNDict():
	print "Weighting nBkgDict..."

	# for sampleID in nBkgDict:
	for sampleID in tqdm(nBkgDict):
		treeWeight = 1
		if sampleID>0:
			mcSumW = sumWdict.get(sampleID, -1) # Get sum of weights
			xSECxEff = -1.
			xSECxEff = xsecDB.xsectTimesEff(sampleID) # Get xSec * filterEff

			treeWeight = xSECxEff * luminosity / mcSumW # Weight this tree
			if treeWeight<=0 : 
				# print "Encounter <=0 weight sample %d , skipped" % sampleID
				tqdm.write("Encounter <=0 weight sample %d , skipped" % sampleID)
				continue
		else:
			treeWeight = 1 * 33257.2 / 10064.3
		for chan in channels:
			nBkgDict[sampleID][chan] *= treeWeight
			nBkgTotDict[chan] = nBkgTotDict.get(chan, 0) + nBkgDict[sampleID][chan]

		# Save to nBkgSourceDict
		bkgSource = sourceBkg(sampleID)
		if nBkgSourceDict.get(bkgSource, None) is None:
			nBkgSourceDict[bkgSource] = nBkgDict[sampleID]
		else:
			for chan in channels:
				nBkgSourceDict[bkgSource][chan] += nBkgDict[sampleID][chan]

	print "nBkgDict weighted, nBkgTotDict loaded"

	print "Weighting nSigDict..."
	# for sampleID in nSigDict:
	for sampleID in tqdm(nSigDict):
		if sampleID<=0: continue
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

# Save bkg counts to CSV
def saveDict():
	with open("nSigDict.csv", "w") as outCSV:

		# Header row
		outCSV.write("sampleID,")
		outCSV.write(",".join(["%d" % int(chan) for chan in channels]))
		outCSV.write("\n")

		for sampleID in nSigDict:
			outCSV.write("%d,"%sampleID)
			outCSV.write(",".join("%f" % nSigDict[sampleID][chan] for chan in channels))
			outCSV.write("\n")

	with open("nBkgDict.csv", "w") as outCSV:

		# Header row
		outCSV.write("sampleID,")
		outCSV.write(",".join(["%d" % int(chan) for chan in channels]))
		outCSV.write("\n")

		for sampleID in nBkgDict:
			outCSV.write("%d,"%sampleID)
			outCSV.write(",".join("%f" % nBkgDict[sampleID][chan] for chan in channels))
			outCSV.write("\n")

	with open("nBkgSourceDict.csv", "w") as outCSV:

		# Header row
		outCSV.write("source,")
		outCSV.write(",".join(["%d" % int(chan) for chan in channels]))
		outCSV.write("\n")

		for sampleID in nBkgSourceDict:
			outCSV.write("%s,"%sampleID)
			outCSV.write(",".join("%f" % nBkgSourceDict[sampleID][chan] for chan in channels))
			outCSV.write("\n")

saveDict()
