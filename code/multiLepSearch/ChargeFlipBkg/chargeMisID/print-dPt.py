#!/usr/bin/env python
import ROOT
from ROOT import TH2, TFile, gStyle, TCanvas, TLegend

gStyle.SetOptStat(0)
gStyle.SetPalette(ROOT.kBird)

f = TFile("output.root")

Total2D = f.Get("hDPtPt")
PtFlipped = f.Get("hDPtFlippedPt")
PtNoFlip = f.Get("hDPtNoFlipPt")
EtaFlipped = f.Get("hDPtFlippedEta")
EtaNoFlip = f.Get("hDPtNoFlipEta")

c = TCanvas("c", "c", 1)
c.SetLogz()

EtaFlipped.Draw(" colz")
c.Print("EtaFlipped.pdf")

EtaNoFlip.Draw(" colz")
c.Print("EtaNoFlip.pdf")

c.SetLogx()

PtFlipped.Draw(" colz")
c.Print("PtFlipped.pdf")

PtNoFlip.Draw(" colz")
c.Print("PtNoFlip.pdf")

####
c.SetLogz(False)
c.SetLogx(False)

total = Total2D.ProjectionY()
flipped = PtFlipped.ProjectionY()
noFlip = PtNoFlip.ProjectionY()

c.SetLogy()

total.Draw()

flipped.SetLineColor(ROOT.kGreen+3)
flipped.Draw("same")

noFlip.SetLineColor(ROOT.kRed)
noFlip.Draw("same")

leg = TLegend(0.15, 0.7, 0.4, 0.9)
leg.AddEntry(total, "All electrons", "l")
leg.AddEntry(flipped, "Charge flipped electrons", "l")
leg.AddEntry(noFlip, "Charge ok electrons", "l")
leg.Draw()

c.Print("dPt.pdf")
