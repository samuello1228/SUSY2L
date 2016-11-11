#!/usr/bin/env python

# It is possible to write code to loop over a tree/ntuple in python. 
# This is the proof of concept
# BUT! Be warned. The for loop will take an incredibly long time.
# Use C++ for reading trees instead, the speedup is worth it
# -- Gabriel Gallardo 8 July 2016


import ROOT 
from ROOT import TH1D, TFile, TChain

hOldEta = TH1D("hOldEta", "#eta distribution of electrons;|#eta|", 100, -2.5, 2.5)
hOldLPt = TH1D("hOldLPt", "Leading p_{T} distribution;p_{T} (GeV)", 200, 20, 200)
hOldSLPt = TH1D("hOldSLPt", "Subleading p_{T} distribution;p_{T} (GeV)", 200, 20, 200)

hNewEta = hOldEta.Clone("hNewEta")
hNewLPt = hOldLPt.Clone("hNewLPt")
hNewSLPt = hOldSLPt.Clone("hNewSLPt")

f = TFile.Open("PYdistros.root", "recreate")

hOldEta.SetDirectory(f)
hOldLPt.SetDirectory(f)
hOldSLPt.SetDirectory(f)
hNewEta.SetDirectory(f)
hNewLPt.SetDirectory(f)
hNewSLPt.SetDirectory(f)

oldChain = TChain("ZeeCandidate")

with open("../common/inFileList-MCconverted-old.txt") as oldFilelist:
	for filename in oldFilelist:
		oldChain.Add(filename[:-1])
		print "Added ", filename, " to oldChain"

nEntries = oldChain.GetEntries()
print nEntries, " entries in OldChain"
i = 1
for evt in oldChain:
	hOldEta.Fill(evt.elCand1_cl_eta)
	hOldEta.Fill(evt.elCand2_cl_eta)
	hOldLPt.Fill(max(evt.elCand1_pt, evt.elCand2_pt))
	hOldSLPt.Fill(min(evt.elCand1_pt, evt.elCand2_pt))
	if(i%50000 == 0):
		print "Event #", i
	i+=1


newChain = TChain("ZeeCandidate")

with open("../common/inFileList-MCconverted.txt") as newFilelist:
	for filename in newFilelist:
		newChain.Add(filename)
		print "Added ", filename, " to newChain"

nEntries = newChain.GetEntries()
print nEntries, " entries in NewChain"
i = 1
for evt in newChain:
	hNewEta.Fill(evt.elCand1_cl_eta)
	hNewEta.Fill(evt.elCand2_cl_eta)
	hNewLPt.Fill(max(evt.elCand1_pt, evt.elCand2_pt))
	hNewSLPt.Fill(min(evt.elCand1_pt, evt.elCand2_pt))
	if(i%50000 == 0):
		print "Event #", i
	i+=1

f.Write()