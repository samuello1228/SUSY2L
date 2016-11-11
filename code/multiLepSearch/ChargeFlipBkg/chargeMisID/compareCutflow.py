#!/usr/bin/env python
import ROOT
from ROOT import TFile, TCanvas, TLegend

fGabriel = TFile.Open("0926-data-signal-80100-noSub/checks.root")
hGabriel = fGabriel.Get("hCutflow")

fEG = TFile.Open("EG/checks_80_100_0_0.root")
hEG = fEG.Get("hCutflow")

hGabriel.SetLineColor(ROOT.kRed)

leg = TLegend(0.7, 0.8, 0.9, 0.9)
leg.AddEntry(hGabriel, "Gabriel")
leg.AddEntry(hEG, "EGamma")

c = TCanvas()
hGabriel.Draw()
hEG.Draw("same")
leg.Draw()
c.Print("Cutflow.pdf")

with open("Cutflow.txt", 'w') as fOut:
	fOut.write("Gabriel \t EGamma \t G - EG\n")
	for i in range(1,9):
		fOut.write("%0d" % hGabriel.GetBinContent(i))
		fOut.write("\t")
		fOut.write("%0d" % hEG.GetBinContent(i))
		fOut.write("\t")
		fOut.write("%0d" % (hGabriel.GetBinContent(i)-hEG.GetBinContent(i)))
		fOut.write("\n")