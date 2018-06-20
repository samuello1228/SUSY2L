#!/usr/bin/python

import ROOT, sys, os
import os.path

checkpoints = [
	("AOD", False, False),
	("SUSY2", False, False),
	("GRL", False, False),
	("Trigger", False, False),
	("LArErr", False, False),
	("TileErr", False, False),
	("CoreFlag", False, False),
	("PV", False, False),
	("BadMuon", False, False),
	("cosMuon", False, False),
	("BadJet", False, True),

	(">=1BaseEl", False, False),
	(">=1SigEl", False, False),
	(">=1BaseMu", False, False),
	(">=1SigMu", False, False),
	(">=1BaseLep", True, False),
	(">=1SigLep", True, True),

	(">=2SigLep", True, True),

	("=2SigLep", False, False),
	("=2BaseLep and =2SigLep", True, False),
	("=3SigLep", False, False),
	("=3BaseLep and =3SigLep", False, False),
]

def printCutflow(h):
	for cp in checkpoints:
		binNum = h.GetXaxis().FindBin(cp[0])
		n = h.GetBinContent(binNum)

		print "%25s : %f" % (cp[0], n),

		if cp[1]:
			binNum = h.GetXaxis().FindBin("%s,w" % cp[0])
			n = h.GetBinContent(binNum)
			print "   %f" % n
		else:
			print ""

		if cp[2]: 
			print ""


def printFile(pair):
	fname, eventNum = pair[0], pair[1]
	f = ROOT.TFile(fname)
	hCutFlow = f.Get("hCutFlow")

	printCutflow(hCutFlow)
	
	t = f.Get("evt2l")
	print "Weights for event %d" % eventNum
	t.Scan("evt.weight:evt.pwt:evt.ElSF:evt.MuSF:evt.BtagSF:evt.JvtSF:evt.trigSF", "evt.event==%d"%eventNum)
	

fileList = [
	("./%s.Sig/data-myOutput/test.root" % sys.argv[1], 10566),
	("./%s.Zmm/data-myOutput/test.root" % sys.argv[1], 47),
	("./%s.ttbar/data-myOutput/test.root" % sys.argv[1], 5799815),
	("./%s.VV/data-myOutput/test.root" % sys.argv[1], 263),
	("./%s.data/data-myOutput/test.root" % sys.argv[1], -1)
]

for f in fileList:
	if os.path.isfile(f[0]):
		print "Cutflow for ", f
		printFile(f)
		print '\n\n'