#!/usr/bin/env python

import sys, os, subprocess, ssUtil, ROOT
from os import listdir, path, mkdir, rename
from os.path import isdir, isfile
import re
import overtraining
from subprocess import call
ROOT.gROOT.Macro("$ROOTCOREDIR/scripts/load_packages.C")
from ROOT.SUSY import CrossSectionDB
xsecDB = CrossSectionDB("%s/data/SUSYTools/mc15_13TeV/"%os.environ["ROOTCOREBIN"] )

#####
# Config variables
#####
C1masses = (200, 300, 400, 500, 600, 700, 800, 900, 1000)
luminosity = 33257.2
channels = (0, 1, 2, 3, 4, 10, 11, 12, 13, 14)
sigFilesTxt = "CutOpt/GabrielFiles/allSig.txt"
MCbkgFilesTxt = "CutOpt/GabrielFiles/MCbkgFiles.txt"
dataBkgFilesTxt = "CutOpt/GabrielFiles/data2016.txt"
#####

sampleIdDB = {}
sumWdict = {}
sumWdict.update( ssUtil.getSumWfromNTUPList( sigFilesTxt ) )
sumWdict.update( ssUtil.getSumWfromNTUPList( options.MCbkgFilesTxt ) )

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

def loadXSec():
	xSecFile = open("SigSample.txt", "r")
	xSecDict = {}
	for aLine in xSecFile:
		aLine = aLine.split("#")[0]
		if len(aLine)==0: continue

		elements = aLine.split(' ')
		mC1 = elements[2]
		mN1 = elements[3]
		xSec = elements[4]
		xSecEff = elements[5]

		nBefore = elements [-2]
		nAfter = elements[-1].rstrip('\n')
		preSelEff = float(nAfter)/float(nBefore)

		xSecDict[(int(mC1), int(mN1))] = (float(xSec), float(xSecEff), preSelEff)
	xSecFile.close()
	return xSecDict

def makeChanDict(ChanNs):
	n = int(ChanNs[0])
	ChanDict = {
		0: float(ChanNs[1])/n,
		1: float(ChanNs[2])/n,
		2: float(ChanNs[3])/n,
		3: float(ChanNs[4])/n,
		4: float(ChanNs[5])/n,
		10: float(ChanNs[6])/n,
		11: float(ChanNs[7])/n,
		12: float(ChanNs[8])/n,
		13: float(ChanNs[9])/n,
		14: float(ChanNs[10])/n
	}
	return ChanDict

def makeSampleIdDB():
	fileList = open(sigFilesTxt,"r").readlines()
	# print fileList
	for aLine in fileList:
		if "MGPy8EG" not in aLine or "root/user" not in aLine: continue
		# if (re.search("MCSig.[0-9]{6}", aLine)) is None: continue
		sampleID = int((re.search("12.[0-9]{6}", aLine).group()).split('.')[1])
		masses = re.search("Slep_[0-9]{3,4}_[0-9]{1,4}", aLine).group()
		print masses
		mC1 = int(masses.split('_')[1])
		mN1 = int(masses.split('_')[2])
		sampleIdDB[(mC1, mN1)] = sampleID
	print sampleIdDB


def loadEffs():
	effFile = open("CutEffs.csv")
	effDict = {}
	next(effFile)
	for line in effFile:
		elements = line.split(',')
		effDict[(int(elements[0]), int(elements[1]))] = makeChanDict(elements[2:])
		# print (int(elements[0]), int(elements[1]))

	effFile.close()
	return effDict

def getNsig(mC1,mN1,chan):
	xSECxEff = xsecDB.xsectTimesEff(sampleIdDB[(mC1,mN1)], 125)
	mcSumW = sumWdict.get(sampleIdDB[(mC1,mN1)], -1)
	if mcSumW < 0: 
		print "mcSumW <0!"
		return 0

	tc = ROOT.TChain("evt2l")
	for line in open(sigFilesTxt).readlines():
		tc.Add(line)

	treeWeight = xSECxEff * luminosity / mcSumW
	if treeWeight<=0:
		print "Encounter <=0 weight sample %s , skipped" % infname 
		return 0
	weight = "(ElSF * MuSF * BtagSF * weight * pwt)"
	return treeWeight*tc.GetEntries("%s*(%s)"% (weight, ssUtil.getCut(chan)))

def getNbkg(chan):
	nBkg = 0
	weight = "(ElSF * MuSF * BtagSF * weight * pwt)"
	for line in open(MCbkgFilesTxt).readlines():
		match = re.search(".[0-9]{6}.", line)
		if not match:
			print "Cannot infer datasetID from filename %s , skipped" % line
			continue

		sampleID = int(match.group()[1:-1])
		mcSumW = sumWdict.get(sampleID, -1)
		xSECxEff = -1.
		xSECxEff = xsecDB.xsectTimesEff(sampleID)

		treeWeight = xSECxEff * options.inputLumi / mcSumW

		if treeWeight<=0 : 
			print "Encounter <=0 weight sample %s , skipped" % line
			continue
		tc = ROOT.TChain("evt2l")
		tc.Add(line)
		
		nBkg += treeWeight*tc.GetEntries("%s*(%s)"% (weight, ssUtil.getCut(chan)))

	tc = ROOT.TChain("evt2l")
	for line in open (dataBkgFilesTxt).readlines():
		tc.Add(line)
	nBkg += tc.GetEntries("%s * (fLwt+qFwt) * (%s)" % (weight, ssUtil.getCut(chan)))
	return nBkg

xSecDict = loadXSec()
effDict = loadEffs()
nBkgDict = {}
for chan in channels:
	nBkgDict[chan]=getNbkg(chan)

if __name__ == '__main__':
	outOT = open("checksOT.csv", "w")
	outOT.write("Channel,ISR,Flavor,dm,NTrees,NodeSize,Depth,SigKS,BkgKS,SigChi2,BkgChi2\n")

	outViableOT = open("viableOT.csv", "w")
	outViableOT.write("Channel,ISR,Flavor,dm,NTrees,NodeSize,Depth,SigKS,BkgKS,SigChi2,BkgChi2\n")

	outSig = open("checksSig.csv", "w")
	outSig.write("m(C1),m(N1),Channel,ISR,Flavor,NTrees,NodeSize,Depth,xSec,")
	outSig.write("nSig,nBkg,BDTopt,nSig(BDTopt),nBkg(BDTopt),sigma(BDTopt),")
	outSig.write("nSig(0),nBkg(0),sigma(0),")
	outSig.write("nSig(0.1),nBkg(0.1),sigma(0.1),")
	outSig.write("nSig(0.2),nBkg(0.2),sigma(0.2),")
	outSig.write("nSig(0.3),nBkg(0.3),sigma(0.3)\n")

	for directory in sys.argv[1:]:
		directory = directory.rstrip('/')
		if not(isdir(directory)): continue
		files = listdir(directory)

		plotDir = "%s/plots" % directory
		if not isdir(plotDir): mkdir(plotDir)

		for file in files:
			if not(file.endswith(".root")):	continue

			print ("## Now processing %s/%s##" % (directory, file))
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

			for mass in C1masses: 
				xSec = xSecDict[(mass, mass-dm)]
				selEff = effDict[(mass, mass-dm)][channel%100]
				nEvents = xSec[0]*xSec[1]*xSec[2]*luminosity*selEff
				print "## mC1 = %d, mN1 = %d, chan = %d, nEvents %f\n" % (mass, mass-dm, channel, nEvents)
				if nEvents <= 1: continue
				call('root -l -b -q "mvaeffs.cxx(\\"%s/%s\\", %d, %f)"' % (directory, file, int(nEvents), luminosity), shell=True)
				outSig.write("%d,%d,%d,%s,%s,%d,%d,%d,%f," 
					% (mass, mass-dm, channel, ISR, Flavor, NTrees, NodeSize, Depth, xSec[0]))
				with open("effs.csv") as effFile:
					outSig.write(effFile.readline())
				os.rename("plots/mvaeffs_BDTD.eps", "%s/%d_%d_Channel%d.eps" % (plotDir, mass, mass-dm, channel))

	outOT.close()
	outViableOT.close()
	outSig.close()