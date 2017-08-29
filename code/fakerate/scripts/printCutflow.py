#!/usr/bin/python

import ROOT, sys, os
import os.path

def printCutflow(pair):
	fname, eventNum = pair[0], pair[1]
	f = ROOT.TFile(fname)
	hCutFlow = f.Get("hCutFlow")

	nBins = hCutFlow.GetXaxis().GetNbins()
	for i in range(1, nBins+1):
		if hCutFlow.GetXaxis().GetBinLabel(i)!=hCutFlow.GetXaxis().GetBinLabel(i+1):
			if hCutFlow.GetXaxis().GetBinLabel(i)=="":
				print ""
			else:
				print "%25s : %f" % (hCutFlow.GetXaxis().GetBinLabel(i), hCutFlow.GetBinContent(i))
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
		printCutflow(f)
		print '\n\n'