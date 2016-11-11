#!/usr/bin/env python
import ROOT, os, sys
from ROOT import TH2, TFile, gStyle, TCanvas, TLegend, TAxis, TLine
from sys import argv

gStyle.SetOptStat(0)
# gStyle.SetPalette(ROOT.kBird)

# os.chdir(argv[1])

f = TFile("mllWeights.root")

SS = f.Get("hMassSS")
ExpSS = f.Get("hMassExpSS")
SF = f.Get("hMassSF")

SS.GetYaxis().SetRangeUser(max(1.0,SS.GetMinimum()*0.5), max(SS.GetMaximum(),ExpSS.GetMaximum())*5)
SS.SetLineColor(ROOT.kAzure+1)
SS.SetFillColor(ROOT.kAzure+1)
SS.SetMarkerColor(ROOT.kAzure+1)

ExpSS.SetLineColor(ROOT.kRed)
ExpSS.SetMarkerColor(ROOT.kRed)

SF.SetLineColor(ROOT.kBlack)
SF.SetMarkerColor(ROOT.kBlack)


c = TCanvas()
c.Divide(1,2)
c.cd(1).SetPad(0, 0.25, 1, 1)
c.cd(2).SetPad(0, 0, 1, 0.25)


c.cd(1)
c.cd(1).SetLogy()
SS.Draw("bar")
ExpSS.Draw("same e")

l = TLegend(0.5, 0.9, 0.9, 0.7)
l.AddEntry(SS, "Observed SS events", "f")
l.AddEntry(ExpSS, "Predicted SS events from OS", "lp")
l.Draw()

c.cd(2)
SF.SetTitle("")
SFx = SF.GetXaxis()
SFx.SetLabelOffset(999)
SFx.SetLabelSize(0)
SFx.SetTitle("")
SFy = SF.GetYaxis()
SFy.SetTitle("Obs/Pred")
SFy.SetTitleSize(0.1)
SFy.SetTitleOffset(0.4)
SFy.SetLabelSize(0.08)
# SF.SetMarkerStyle(20)
SFy.SetRangeUser(0.5,2)
SF.Draw("p e")

l = TLine()
l1 = l.DrawLine(SFx.GetXmin(), 1., SFx.GetXmax(), 1.)
l1.SetLineStyle(2)
l1.SetLineWidth(2)
l1.Draw()

c.Print("mll-unweighted.pdf")


