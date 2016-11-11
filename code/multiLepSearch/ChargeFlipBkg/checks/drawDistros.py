#!/usr/bin/env python
import ROOT
from ROOT import TH1D, TFile, TCanvas

f = TFile.Open("distros.root")

c = TCanvas()

hOldEta = f.Get("hOldEta")
hOldLPt = f.Get("hOldLPt")
hOldSLPt = f.Get("hOldSLPt")
hNewEta = f.Get("hNewEta")
hNewLPt = f.Get("hNewLPt")
hNewSLPt = f.Get("hNewSLPt")

hOldEta.SetLineColor(ROOT.kRed)
hOldLPt.SetLineColor(ROOT.kRed)
hOldSLPt.SetLineColor(ROOT.kRed)
hNewEta.SetLineColor(ROOT.kBlue)
hNewLPt.SetLineColor(ROOT.kBlue)
hNewSLPt.SetLineColor(ROOT.kBlue)

hOldEta.DrawNormalized()
# hNewEta.DrawNormalized("same")
# c.Print("Eta.pdf")

# hOldLPt.DrawNormalized()
# hNewLPt.DrawNormalized("same")
# c.Print("LPt.pdf")

# hOldSLPt.DrawNormalized()
# hNewSLPt.DrawNormalized("same")
# c.Print("SLPt.pdf")

hNewEta.Draw()
hOldEta.Draw("same")
c.Print("Eta.pdf")

hNewLPt.Draw()
hOldLPt.Draw("same")
c.Print("LPt.pdf")

hNewSLPt.Draw()
hOldSLPt.Draw("same")
c.Print("SLPt.pdf")