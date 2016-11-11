#!/usr/bin/env python
import sys
import ROOT
from ROOT import TH1, TFile, TCanvas, TStyle, gStyle, TLine, TLegend

def setRed(h):
	h.SetLineColor(ROOT.kRed)

def drawComparisons(hGabriel, hEGamma):
	setRed(hGabriel)
	c = TCanvas("c", "c", 1)
	c.Divide(1,2)
	c.cd(1).SetPad(0,0.2,1,1)
	c.cd(2).SetPad(0,0,1,0.25)

	c.cd(1)
	hEGamma.Draw()
	hGabriel.Draw("same")

	leg = TLegend(0.7,0.9,0.9,0.75)
	leg.AddEntry(hEGamma, "EGamma", "l")
	leg.AddEntry(hGabriel, "Gabriel", "l")
	leg.Draw()

	c.cd(2)
	hRatio = hGabriel.Clone("hRatio")
	hRatio.Divide(hEGamma)

	hRatio.SetTitle("")
	hRatio.GetXaxis().SetLabelOffset(999);
	hRatio.GetXaxis().SetLabelSize(0);
	hRatio.GetXaxis().SetTitle("");
	hRatio.GetYaxis().SetTitle("Mine/EG")
	hRatio.GetYaxis().SetTitleSize(0.1);
	hRatio.GetYaxis().SetTitleOffset(0.4);
	hRatio.GetYaxis().SetLabelSize(0.08);
	hRatio.SetMarkerStyle(20);
	hRatio.GetYaxis().SetRangeUser(0.985,1.015);
	hRatio.Draw()

	line = TLine()
	line1 = line.DrawLine(hRatio.GetXaxis().GetXmin(), 1., hRatio.GetXaxis().GetXmax(), 1.)
	line1.SetLineStyle(2)
	line1.SetLineWidth(2)
	
	name = hEGamma.GetName()
	name = name + ".pdf"
	c.Print(name)

### MAIN FUNCTION ###

if len(sys.argv)!= 3:
	print "Wrong number of arguments"
	print "Usage: ./check.py Gabriel/file.root EGamma/file.root"
	print ""
	sys.exit()

fGabriel = TFile(sys.argv[1])
fEGamma = TFile(sys.argv[2])

hPtGabriel = []
hPtGabriel.append(fGabriel.Get("hPtAll"))
hPtGabriel.append(fGabriel.Get("hPtPreselected"))
hPtGabriel.append(fGabriel.Get("hPtLoose"))
hPtGabriel.append(fGabriel.Get("hPtSignal"))
hPtGabriel.append(fGabriel.Get("hPtSignalZ"))

hPtEGamma = []
hPtEGamma.append(fEGamma.Get("hPtAll"))
hPtEGamma.append(fEGamma.Get("hPtPreselected"))
hPtEGamma.append(fEGamma.Get("hPtLoose"))
hPtEGamma.append(fEGamma.Get("hPtSignal"))
hPtEGamma.append(fEGamma.Get("hPtSignalZ"))

hEtaGabriel = []
hEtaGabriel.append(fGabriel.Get("hEtaAll"))
hEtaGabriel.append(fGabriel.Get("hEtaPreselected"))
hEtaGabriel.append(fGabriel.Get("hEtaLoose"))
hEtaGabriel.append(fGabriel.Get("hEtaSignal"))
hEtaGabriel.append(fGabriel.Get("hEtaSignalZ"))

hEtaEGamma = []
hEtaEGamma.append(fEGamma.Get("hEtaAll"))
hEtaEGamma.append(fEGamma.Get("hEtaPreselected"))
hEtaEGamma.append(fEGamma.Get("hEtaLoose"))
hEtaEGamma.append(fEGamma.Get("hEtaSignal"))
hEtaEGamma.append(fEGamma.Get("hEtaSignalZ"))

hMassGabriel = []
hMassGabriel.append(fGabriel.Get("hMassAll"))
hMassGabriel.append(fGabriel.Get("hMassPreselected"))
hMassGabriel.append(fGabriel.Get("hMassLoose"))
hMassGabriel.append(fGabriel.Get("hMassSignal"))
hMassGabriel.append(fGabriel.Get("hMassSignalZ"))
hMassGabriel[4].GetXaxis().SetRangeUser(80,100)

hMassEGamma = []
hMassEGamma.append(fEGamma.Get("hMassAll"))
hMassEGamma.append(fEGamma.Get("hMassPreselected"))
hMassEGamma.append(fEGamma.Get("hMassLoose"))
hMassEGamma.append(fEGamma.Get("hMassSignal"))
hMassEGamma.append(fEGamma.Get("hMassSignalZ"))
hMassEGamma[4].GetXaxis().SetRangeUser(80,100)

hCutflowGabriel = fGabriel.Get("hCutflow")
hCutflowEGamma = fEGamma.Get("hCutflow")

gStyle.SetOptStat(0)

for i in range(len(hPtGabriel)):
	drawComparisons(hPtGabriel[i], hPtEGamma[i])
	drawComparisons(hEtaGabriel[i], hEtaEGamma[i])
	drawComparisons(hMassGabriel[i], hMassEGamma[i])

drawComparisons(hCutflowGabriel, hCutflowEGamma)

