#!/usr/bin/env python

import multiLepSearch, subprocess, ROOT, os
from multiLepSearch import ssUtil
# from ssUtil import getCut
from subprocess import call
from os import listdir, path


channels = (0, 1, 2, 3, 4, 10, 11, 12, 13, 14)

rootDir="/home/Shared/SUSY/ychan/sig"
dirList = listdir(rootDir)

def loadData(directory):
	fileList = listdir("%s/%s" % (rootDir, directory))
	chain = ROOT.TChain("evt2l")
	for f in fileList:
		chain.Add("%s/%s/%s" % (rootDir, directory, f))
	return chain

def getCutEffs():
	outCSV = open("CutEffs.csv", "w")
	outCSV.write("mC1,mN1,nEntries,0,1,2,3,4,10,11,12,13,14\n")
	for directory in dirList:
		if not path.isdir("%s/%s" % (rootDir, directory)): continue
		chain = loadData(directory)

		elements = directory.split("_")
		mC1 = int(elements[4])
		mN1 = int(elements[5])
		print "%s: %d" % (directory, chain.GetEntries())

		outCSV.write("%d,%d,%d" % (mC1, mN1, chain.GetEntries()))
		for chan in channels:
			if (mC1 - mN1) == 100: chan += 100
			outCSV.write(",%d" % chain.GetEntries(ssUtil.getCut(chan)))
		outCSV.write("\n")

	outCSV.close()

if __name__ == "__main__":
	getCutEffs()

