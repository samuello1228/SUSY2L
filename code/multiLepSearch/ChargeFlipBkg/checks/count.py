import ROOT
from ROOT import TFile

nEventsOld = 0
nEventsNew = 0

with open("../common/inFileList-MCdR-old.txt") as oldFilelist:
	for filename in oldFilelist:
		f = TFile.Open(filename[:-1])
		hCutflow = f.Get("hCutFlow")
		nEvents = int(hCutflow.GetBinContent(1))
		print nEvents, "in", filename
		nEventsOld = nEventsOld + nEvents
print nEventsOld, "in old NTuples"

with open("../common/inFileList-MCdR.txt") as newFilelist:
	for filename in newFilelist:
		f = TFile.Open(filename[:-1])
		hCutflow = f.Get("hCutFlow")
		nEvents = int(hCutflow.GetBinContent(1))
		print nEvents, "in", filename
		nEventsNew = nEventsNew + nEvents

print nEventsOld, "in old NTuples"
print nEventsNew, "in new NTuples"
print "New - old = ", nEventsNew - nEventsOld