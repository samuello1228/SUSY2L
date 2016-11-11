#!/usr/bin/env python
import ROOT
from ROOT import TH1, TFile, TLegend, TCanvas, gStyle

gStyle.SetOptStat(0)

f = TFile.Open("ttbarHistos.root")

hLeadingPtOS = f.Get("hLeadingPtOS")
hLeadingPtSS = f.Get("hLeadingPtSS")
hSubleadingPtOS = f.Get("hSubleadingPtOS")
hSubleadingPtSS = f.Get("hSubleadingPtSS")
hInvMassOS = f.Get("hInvMassOS")
hInvMassSS = f.Get("hInvMassSS")
hLeadingEtaOS = f.Get("hLeadingEtaOS")
hLeadingEtaSS = f.Get("hLeadingEtaSS")
hSubleadingEtaOS = f.Get("hSubleadingEtaOS")
hSubleadingEtaSS = f.Get("hSubleadingEtaSS")
hLeadingPhiOS = f.Get("hLeadingPhiOS")
hLeadingPhiSS = f.Get("hLeadingPhiSS")
hSubleadingPhiOS = f.Get("hSubleadingPhiOS")
hSubleadingPhiSS = f.Get("hSubleadingPhiSS")

hLeadingPtSS.SetLineColor(ROOT.kRed)
hSubleadingPtSS.SetLineColor(ROOT.kRed)
hInvMassSS.SetLineColor(ROOT.kRed)
hLeadingEtaSS.SetLineColor(ROOT.kRed)
hSubleadingEtaSS.SetLineColor(ROOT.kRed)
hLeadingPhiSS.SetLineColor(ROOT.kRed)
hSubleadingPhiSS.SetLineColor(ROOT.kRed)

leg = TLegend(0.7, 0.8, 0.9, 0.9)
leg.AddEntry(hLeadingPtOS, "OS e^{#pm} e^{#mp}")
leg.AddEntry(hLeadingPtSS, "SS e^{#mp} e^{#mp}")

can = TCanvas()
can.SetLogy()

hLeadingPtOS.GetYaxis().SetRangeUser(0.5, hLeadingPtOS.GetMaximum()*1.9)
hLeadingPtOS.Draw()
hLeadingPtSS.Draw("same")
leg.Draw()
can.Print("hLeadingPt.pdf")

hSubleadingPtOS.GetYaxis().SetRangeUser(0.5, hLeadingPtOS.GetMaximum()*1.9)
hSubleadingPtOS.Draw()
hSubleadingPtSS.Draw("same")
leg.Draw()
can.Print("hSubleadingPt.pdf")

hInvMassOS.GetYaxis().SetRangeUser(0.5, hLeadingPtOS.GetMaximum()*1.9)
hInvMassOS.Draw()
hInvMassSS.Draw("same")
leg.Draw()
can.Print("hInvMass.pdf")

hLeadingEtaOS.GetYaxis().SetRangeUser(0.5, hLeadingPtOS.GetMaximum()*1.9)
hLeadingEtaOS.Draw()
hLeadingEtaSS.Draw("same")
leg.Draw()
can.Print("hLeadingEta.pdf")

hSubleadingEtaOS.GetYaxis().SetRangeUser(0.5, hLeadingPtOS.GetMaximum()*1.9)
hSubleadingEtaOS.Draw()
hSubleadingEtaSS.Draw("same")
leg.Draw()
can.Print("hSubleadingEta.pdf")

hLeadingPhiOS.GetYaxis().SetRangeUser(0.5, hLeadingPtOS.GetMaximum()*1.9)
hLeadingPhiOS.Draw()
hLeadingPhiSS.Draw("same")
leg.Draw()
can.Print("hLeadingPhi.pdf")

hSubleadingPhiOS.GetYaxis().SetRangeUser(0.5, hLeadingPtOS.GetMaximum()*1.9)
hSubleadingPhiOS.Draw()
hSubleadingPhiSS.Draw("same")
leg.Draw()
can.Print("hSubleadingPhi.pdf")
