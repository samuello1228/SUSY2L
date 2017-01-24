#!/usr/bin/env python

import ROOT,csv

def makeSigHist():
	ROOT.gStyle.SetOptStat(0)
	ROOT.gStyle.SetPalette(ROOT.kAquamarine)
	outFile = ROOT.TFile("significance.root", "recreate")
	hSig = ROOT.TH2D("hSig", "Max significance;mC1 [GeV];mN1 [GeV]", 41, 180, 1000, 46, 80, 1000)
	c = ROOT.TCanvas()

	masses = (200, 300, 400, 500, 600, 700, 800, 900, 1000)
	dms = (20, 50, 100)
	channels = (0, 1, 2, 3, 4, 10, 11, 12, 13, 14, 100, 101, 102, 103, 104, 110, 111, 112, 113, 114)

	with open("checksSig.csv") as sigFile:
		reader = csv.reader(sigFile)
		next(reader)
		sigDict = {(int(row[0]), int(row[1]), int(row[2])): float(row[11]) for row in reader}
		# print sigDict

	for mass in masses:
		for dm in dms:
			maxSig=0
			for channel in channels:
				if dm == 100: channel +=100
				try:
					sig = sigDict[(mass, mass-dm, channel)]
				except KeyError:
					continue
				if sig > maxSig: maxSig = sig
			hSig.Fill(mass, mass-dm, maxSig)
	hSig.Draw("colz text")
	line = ROOT.TLine()
	l0 = line.DrawLine(180, 180, 1000, 1000)
	l0.SetLineStyle(1)
	l0.SetLineWidth(1)

	l20 = line.DrawLine(180, 160, 1000, 980)
	l20.SetLineStyle(3)
	l20.SetLineWidth(1)

	l50 = line.DrawLine(180, 130, 1000, 950)
	l50.SetLineStyle(3)
	l50.SetLineWidth(1)

	l100 = line.DrawLine(180, 80, 1000, 900)
	l100.SetLineStyle(3)
	l100.SetLineWidth(1)

	c.Print("sig.pdf")
	outFile.Write()
	outFile.Close()


if __name__ == '__main__':
	makeSigHist()