#!/usr/bin/env python

import sys, os, subprocess
from os import listdir, path
from os.path import isdir
import re
import overtraining
from subprocess import call

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
	return (SigKS > 0.05 and BkgKS > 0.05 and SigChi2 > 0.05 and BkgChi2 > 0.05 
		and SigKS < 0.95 and BkgKS < 0.95 and SigChi2 < 0.95 and BkgChi2 < 0.95)

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
	ChanDict = {
		n: ChanNs[0],
		0: ChanNs[1],
		1: ChanNs[2],
		2: ChanNs[3],
		3: ChanNs[4],
		4: ChanNs[5],
		10: ChanNs[6],
		11: ChanNs[7],
		12: ChanNs[8],
		13: ChanNs[9],
		14: ChanNs[10]
	}
	
def loadEffs():
	effFile = open("CutEffs.csv")
	effDict = {}

	for line in effFile:
		elements = line.split(',')
		efficiency[(int(elements[0]), int(elements[1]))] = makeChanDict(elements[2:])

	effFile.close()
	return effDict

C1masses = (200, 300, 400, 500, 600, 700, 800, 900, 1000)
luminosity = 10064.3 # /1000 pb-1 # TODO: Update for higher luminosity


if __name__ == '__main__':
	outOT = open("checksOT.csv", "w")
	outOT.write("Channel,ISR,Flavor,dm,NTrees,NodeSize,Depth,SigKS,BkgKS,SigChi2,BkgChi2\n")

	outSig = open("checksSig.csv", "w")
	outSig.write("m(C1),m(N1),Channel,ISR,Flavor,NTrees,NodeSize,Depth,xSec,nSig,nBkg,MaxSig,BDTcut\n")

	xSecDict = loadXSec()

	for directory in sys.argv[1:]:
		directory = directory.rstrip('/')
		if not(isdir(directory)): continue
		files = listdir(directory)

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

			if not(notOvertrained(overtrainingNums)): continue

			# TODO: Missing preselection efficiency
			for mass in C1masses:
				xSec = xSecDict[(mass, mass-dm)]
				nEvents = xSec[0]*xSec[1]*xSec[2]*luminosity
				call('root -l -b -q "mvaeffs.cxx(\\"%s/%s\\", %d)"' % (directory, file, int(nEvents)), shell=True)
				outSig.write("%d,%d,%d,%s,%s,%d,%d,%d,%f," 
					% (mass, mass-dm, channel, ISR, Flavor, NTrees, NodeSize, Depth, xSec[0]))
				with open('effs.csv') as effFile:
					outSig.write(effFile.readline())

	outOT.close()