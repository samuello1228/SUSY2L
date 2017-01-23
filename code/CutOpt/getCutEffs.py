#!/usr/bin/env python

import ssUtil, subprocess, ROOT
from ssUtil import getCut
from subprocess import call


channels = (0, 1, 2, 3, 4, 10, 11, 12, 13, 14)

rootDir="/home/Shared/SUSY/ychan/sig"
dirList = call("ls %s" % rootDir, shell=True)

def loadData(directory)
	fileList = call('ls %s/%s | grep -E "root$"' % (rootDir, directory), shell=True)
	chain = ROOT.TChain("evt2l")
	for f in fileList:
		chain.Add(f)
	return chain

def getCutEffs():
	outCSV = open("CutEffs.csv", "w")
	outCSV.write("mC1,mN1,nEntries,0,1,2,3,4,10,11,12,13,14\n")
	for directory in dirList:
		chain = loadData(directory)

		elements = directory.split("_")
		mC1 = int(elements[4])
		mN1 = int(elements[5])

		outCSV.write("%d,%d,%d" % (mC1, mN1, chain.GetEntries()))
		for chan in channels:
			if (mC1 - mN1) == 100: chan += 100
			outCSV.write("%d", chain.GetEntries(getCut(chan)))
		outCSV.write("\n")

	outCSV.close()

if __name__ == "__main__"
	getCutEffs()

