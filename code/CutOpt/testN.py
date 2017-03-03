#!/usr/bin/env python

import makeOTPlots
from makeOTPlots import getNDict, getTotDict, makeSampleIdDB

nBkgDict = getTotDict(getNDict("nBkgDict.csv"))
nSigDict = getNDict("nSigDict.csv")
sampleIdDB = makeSampleIdDB()

print "nSignal:"
for m in range(200,1000, 100):
	print (m, m-20), nSigDict[sampleIdDB[(m,m-20)]]

print "\nnBackground:"
print nBkgDict