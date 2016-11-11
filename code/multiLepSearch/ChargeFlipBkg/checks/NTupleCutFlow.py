import ROOT
from ROOT import TFile, TH1

oldFiles = []
newFiles = []

oldFilesTxt = "../common/inFileList-MCdR-old.txt"
newFilesTxt = "../common/inFileList-MCdR.txt"
outFilename = "NTupleCutFlow.root"

with open(oldFilesTxt) as f:
	for filename in f:
		oldFiles.append(filename[:-1])

print "Read", oldFilesTxt

with open(newFilesTxt) as f:
	for filename in f:
		newFiles.append(filename[:-1])
print "Read", newFilesTxt

firstLoop = True
for f in oldFiles:
	file = TFile.Open(f) 
	if firstLoop:
		hOldCutFlow = file.Get("hCutFlow")
		hOldCutFlow.SetDirectory(0)
		firstLoop = False
	else:
		hOldCutFlow.Add(file.Get("hCutFlow"))
	file.Close()
print "Added hCutFlow histograms from old files"

firstLoop = True
for f in newFiles:
	file = TFile.Open(f)
	if firstLoop:
		hNewCutFlow = file.Get("hCutFlow")
		hNewCutFlow.SetDirectory(0)
		firstLoop = False
	else:
		hNewCutFlow.Add(file.Get("hCutFlow"))
	file.Close()
print "Added hCutFlow histograms from new files"

outFile = TFile.Open(outFilename, "recreate")
hOldCutFlow.SetDirectory(outFile)
hOldCutFlow.SetName("hOldCutFlow")
hNewCutFlow.SetDirectory(outFile)
hNewCutFlow.SetName("hNewCutFlow")
outFile.Write()
print "Histograms written to", outFilename
outFile.Close()