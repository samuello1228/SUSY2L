#!/usr/bin/env python

import makeOTPlots
from makeOTPlots import getNDict, getTotDict, makeSampleIdDB
import multiLepSearch
from multiLepSearch.ssUtil import loadSampleDict, guessSampleType

MCbkgFilesTxt = "CutOpt/GabrielFiles/MCbkgFiles.txt"
sampleDictTxt = "multiLepSearch/script/MCBgIDs_2LSS.txt"


nBkgDict = getNDict("nBkgDict.csv")
nBkgSourceDict = {}
nSigDict = getNDict("nSigDict.csv")
sampleIdDB = makeSampleIdDB()

MCbkgList = open(MCbkgFilesTxt).readlines()
sampleDict = loadSampleDict(sampleDictTxt)
def sourceBkg(sampleID):
	if sampleID == 0: return "qF or fL"
	sampleType = sampleDict.get(sampleID, {}).get("phyType", None)
	if sampleType is not None: return sampleType
	sampleID = "%d" % sampleID; 
	for s in MCbkgList:
		if sampleID in s: 
			sampleType = guessSampleType(s)
			break 
	return sampleID if sampleType is None else sampleType

print "nSignal:"
for m in range(200,1000, 100):
	print (m, m-20), nSigDict[sampleIdDB[(m,m-20)]]

print "\nnBackground:"
for sampleID in nBkgDict:
	bkgSource = sourceBkg(sampleID)
	if nBkgSourceDict.get(bkgSource, None) is None:
		nBkgSourceDict[bkgSource] = nBkgDict[sampleID]
	else:
		for chan in (channels for channels in nBkgDict[0]):
			nBkgSourceDict[bkgSource][chan] += nBkgDict[sampleID][chan]
for source in nBkgSourceDict:
	print source, nBkgSourceDict[source]
print "Total", getTotDict(nBkgDict)


