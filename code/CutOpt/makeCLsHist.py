#!/usr/bin/env python

import ROOT,csv

col={
	-1 : 14,
	0  : 17,
	0.1: 20, 
	0.2: 23,
	0.3: 26
}

def DrawLines():
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

def makeSigDict(cut):
	with open("checksCLs.csv") as sigFile:
		reader = csv.DictReader(sigFile)
		next(reader)
		cut = "%1.1f" % cut if cut > 0 else "0" if cut == 0 else "BDTopt"
		sigDict = {(int(row["m(C1)"]), int(row["m(N1)"]), int(row["Channel"]), int(row["NTrees"]), int(row["NodeSize"]), int(row["Depth"])): float(row["CLs(%s)"%cut]) if float(row["CLs(%s)"%cut]) >0 else 0 for row in reader}
		# print sigDict
	return sigDict

def makeSigHist():
	ROOT.gStyle.SetOptStat(0)
	ROOT.gStyle.SetPalette(ROOT.kAquamarine)
	outFile = ROOT.TFile("significance.root", "recreate")
	hSigMax = ROOT.TH2D("hSigMax", "Max significance at each point;mC1 [GeV];mN1 [GeV]", 41, 180, 1000, 46, 80, 1000)
	hSig0 = ROOT.TH2D("hSig0", "Max significance at BDT=0;mC1 [GeV];mN1 [GeV]", 41, 180, 1000, 46, 80, 1000)
	hSig1 = ROOT.TH2D("hSig1", "Max significance at BDT=0.1;mC1 [GeV];mN1 [GeV]", 41, 180, 1000, 46, 80, 1000)
	hSig2 = ROOT.TH2D("hSig2", "Max significance at BDT=0.2;mC1 [GeV];mN1 [GeV]", 41, 180, 1000, 46, 80, 1000)
	hSig3 = ROOT.TH2D("hSig3", "Max significance at BDT=0.2;mC1 [GeV];mN1 [GeV]", 41, 180, 1000, 46, 80, 1000)
	c = ROOT.TCanvas()

	masses = (200, 300, 400, 500, 600, 700, 800, 900, 1000)
	dms = (20, 50, 100)
	channels = (0, 1, 2, 3, 4, 10, 11, 12, 13, 14)
	nsizes = (5, 7, 10)
	depths = (2, 3, 4, 5)
	ntrees = (100, 200, 400, 600)

	sigMaxDict = makeSigDict(-1)
	sig0Dict = 	makeSigDict(0)
	sig1Dict = makeSigDict(0.1)
	sig2Dict = makeSigDict(0.2)
	sig3Dict = makeSigDict(0.3)

	for mass in masses:
		for dm in dms:
			minCLs=0.
			sig0 = 0.
			sig1 = 0.
			sig2 = 0.
			sig3 = 0.
			for channel in channels: 
				for nsize in nsizes: 
					for depth in depths: 
						for n in ntrees:
							if dm == 100: channel +=100
							try:
								sig = sigMaxDict[(mass, mass-dm, channel, n, nsize, depth)]
							except KeyError:
								sig = 0
							if sig < minCLs: minCLs = sig

							try:
								sig = sig0Dict[(mass, mass-dm, channel, n, nsize, depth)]
							except KeyError:
								sig = 0
							if sig < sig0: sig0 = sig

							try:
								sig = sig1Dict[(mass, mass-dm, channel, n, nsize, depth)]
							except KeyError:
								sig = 0
							if sig < sig1: sig1 = sig

							try:
								sig = sig2Dict[(mass, mass-dm, channel, n, nsize, depth)]
							except KeyError:
								sig = 0
							if sig < sig2: sig2 = sig

							try:
								sig = sig3Dict[(mass, mass-dm, channel, n, nsize, depth)]
							except KeyError:
								sig = 0
							if sig < sig3: sig3 = sig

			hSigMax.Fill(mass, mass-dm, minCLs)
			hSig0.Fill(mass, mass-dm, sig0)
			hSig1.Fill(mass, mass-dm, sig1)
			hSig2.Fill(mass, mass-dm, sig2)
			hSig3.Fill(mass, mass-dm, sig3)

	hSigMax.Draw("colz text")
	DrawLines()
	c.Print("CLsMin.pdf")

	hSig0.Draw("colz text")
	DrawLines()
	c.Print("CLs0.pdf")

	hSig1.Draw("colz text")
	DrawLines()
	c.Print("CLs1.pdf")

	hSig2.Draw("colz text")
	DrawLines()
	c.Print("CLs2.pdf")

	hSig3.Draw("colz text")
	DrawLines()
	c.Print("CLs3.pdf")

	outFile.Write()
	outFile.Close()


if __name__ == '__main__':
	makeSigHist()