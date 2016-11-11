#!/usr/bin/env python
import ROOT
from ROOT import TH1I, TChain, TCanvas

chGabriel = TChain("evt2l")
chEGamma = TChain("ZeeCandidate")

c = TCanvas("c", "Canvas", 1)
hGabriel = TH1I("hGabriel", "Compare", 10, 0, 10)
hGabriel.Fill(1)
hGabriel.Draw()
c.Print("try.pdf")