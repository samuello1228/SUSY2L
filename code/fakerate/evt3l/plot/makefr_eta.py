import ROOT
import math
from ROOT import TCanvas, THStack, TH1F, gRandom, TPad, TLegend, TLatex
import sys, getopt
import os

c=ROOT.TCanvas()
result = ROOT.TFile("graph.root")
ht = result.Get("t_evt3l_lep_eta_sub")
hl = result.Get("l_evt3l_lep_eta_sub")
houtput  = ht
nbins = ht.GetNbinsX()
for i in range (1,nbins+1):
    numerator = ht.GetBinContent(i)
    denominator = hl.GetBinContent(i)
    if(denominator!=0):
        ratio  = numerator/denominator
        error = math.sqrt(1/denominator * ratio * (1-ratio))
    else:
        ratio = 0
        error = 0
    print numerator, denominator, ratio, error
    houtput.SetBinContent(i, ratio)
    houtput.SetBinError(i, error)
#tight.Divide(tight, loose, 1, 1, "B")
#houtput.SetMarkerStyle(8)
houtput.SetMaximum(1)
houtput.SetMinimum(0)
houtput.SetTitle("Fake Rate vs. Eta")
houtput.GetXaxis().SetTitle("Eta")
houtput.GetYaxis().SetTitle("Fake Rate")
houtput.Draw("pe")

c.SaveAs("FakeRate_Eta.png")
result.Close()

