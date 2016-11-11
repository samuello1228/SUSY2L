import ROOT
from ROOT import TH2D, TH1, TCanvas, TLegend

f = ROOT.TFile("output.root")

Total2D = f.Get("hDPtPt")
Flipped2D = f.Get("hDPtFlippedPt")
NoFlip2D = f.Get("hDPtNoFlipPt")

total = Total2D.ProjectionY()
flipped = Flipped2D.ProjectionY()
noFlip = NoFlip2D.ProjectionY()

c = TCanvas("c", "c", 1)
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

