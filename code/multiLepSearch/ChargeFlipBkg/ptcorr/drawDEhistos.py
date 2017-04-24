#!/usr/bin/env python
import ROOT
from ROOT import TH3D, TFile, TCanvas, gStyle, TLegend, TLine, gPad
import sys

ROOT.gROOT.SetBatch(True)
gStyle.SetOptStat(0)
gStyle.SetPalette(ROOT.kLightTerrain)
f = TFile.Open(sys.argv[1], "update")

h = []
h.append(f.Get("hDE"))
h.append(f.Get("hDEflipped"))
h.append(f.Get("hDEok"))

h.append(f.Get("hDPT"))
h.append(f.Get("hDPTflipped"))
h.append(f.Get("hDPTok"))

# h.append(f.Get("hDR"))
# h.append(f.Get("hDRflipped"))
# h.append(f.Get("hDRok"))

c = TCanvas("c", "c", 1000, 500)
c.SetLogx()

# hAvg = []
for hist in h:
	hProj = hist.Project3DProfile("xy")
	hProj.Draw("colz text")
	hProj.SetDirectory(f)
	filename = hist.GetName() + ".pdf"
	c.Print(filename)
	# hAvg.append(hist.Project3DProfile("xy").Clone(hist.GetName()+"Avg"))
f.Write()

h[1].SetLineColor(ROOT.kRed)
h[4].SetLineColor(ROOT.kRed)
# h[7].SetLineColor(ROOT.kRed)

h[2].SetLineColor(ROOT.kGreen+4)
h[5].SetLineColor(ROOT.kGreen+4)
# h[8].SetLineColor(ROOT.kGreen+4)

etaAxis = h[0].GetXaxis()
ptAxis = h[0].GetYaxis()
nEtaBins = etaAxis.GetNbins()
nPtBins = ptAxis.GetNbins()

cDE = TCanvas()
cDE.SetLogy()

cDPT = TCanvas()
cDPT.SetLogy()

cDR = TCanvas()
cDR.SetLogy()

cDE.Print("dE.pdf[")
cDPT.Print("dPT.pdf[")
cDR.Print("dR.pdf[")

leg = TLegend(0.6, 0.9, 0.9, 0.7)
leg.AddEntry(h[0], "All electrons", "l")
leg.AddEntry(h[1], "Flipped electrons", "l")
leg.AddEntry(h[2], "ok electrons", "l")

cDE.Divide(3,2)
cDPT.Divide(3,2)
cDR.Divide(3,2)

l = TLine()
l.SetLineStyle(2)

for i in range(nPtBins):
	ptBinLow = ptAxis.GetBinLowEdge(i+1)
	ptBinUp = ptAxis.GetBinUpEdge(i+1)
	sPDFTitle = str(ptBinLow) + "_p_{T}_" + str(ptBinUp)
	legArr = []
	for j in range(nEtaBins-1):
		k = j
		if(j>2): 
			k += 1

		etaBinLow = etaAxis.GetBinLowEdge(k+1)
		etaBinUp= etaAxis.GetBinUpEdge(k+1)

		sLegendHeader = str(etaBinLow) + "<|#eta|<" + str(etaBinUp) + ", "  + str(ptBinLow) + "<p_{T}<" + str(ptBinUp)
		sName = str(etaBinLow) + "_" + str(etaBinUp) + "__" + str(ptBinLow) + "_" + str(ptBinUp) + "_"
		legArr.append(leg.Clone())
		legArr[j].SetHeader(sLegendHeader)	

		hProj = []
		hProj.append(h[0].ProjectionZ(sName + "DEall", k+1, k+1, i+1, i+1))
		hProj.append(h[1].ProjectionZ(sName + "DEflipped", k+1, k+1, i+1, i+1))
		hProj.append(h[2].ProjectionZ(sName + "DEok", k+1, k+1, i+1, i+1))

		hProj.append(h[3].ProjectionZ(sName + "DPTall", k+1, k+1, i+1, i+1))
		hProj.append(h[4].ProjectionZ(sName + "DPTflipped", k+1, k+1, i+1, i+1))
		hProj.append(h[5].ProjectionZ(sName + "DPTok", k+1, k+1, i+1, i+1))

		# hProj.append(h[6].ProjectionZ(sName + "DRall", k+1, k+1, i+1, i+1))
		# hProj.append(h[7].ProjectionZ(sName + "DRflipped", k+1, k+1, i+1, i+1))
		# hProj.append(h[8].ProjectionZ(sName + "DRok", k+1, k+1, i+1, i+1))

		cDE.cd(j+1).SetLogy()
		hProj[0].Draw()
		hProj[1].Draw("same")
		hProj[2].Draw("same")
		legArr[j].Draw()
		l.DrawLine(0, 0, 0, hProj[0].GetMaximum()*1.9)

		cDPT.cd(j+1).SetLogy()
		hProj[3].Draw()
		hProj[4].Draw("same")
		hProj[5].Draw("same")
		legArr[j].Draw()	
		l.DrawLine(0, 0, 0, hProj[3].GetMaximum()*1.9)

		cDR.cd(j+1).SetLogy()
		# hProj[6].Draw()
		# hProj[7].Draw("same")
		# hProj[8].Draw("same")
		# legArr[j].Draw()	
		# l.DrawLine(0, 0, 0, hProj[6].GetMaximum()*1.9)

		# for m in range(6):
		# 	hAvg[m].SetBinContent(i+1, k+1, hProj[m].GetMean())
		# 	if m==4:
		# 		print ptBinLow, ptBinUp, etaBinLow, etaBinUp, "dPt flipped=", hProj[m].GetMean()
		# 	hAvg[m].SetBinError(i+1, k+1, hProj[m].GetMeanError())

	cDE.Print("dE.pdf", "Title:"+sPDFTitle)
	cDPT.Print("dPT.pdf", "Title:"+sPDFTitle)
	cDR.Print("dR.pdf", "Title:"+sPDFTitle)

cDE.Print("dE.pdf]")
cDPT.Print("dPT.pdf]")
cDR.Print("dR.pdf]")

# c.cd()
# for hist in hAvg:
# 	hist.Draw("colz text")
# 	filename = hist.GetName() + ".pdf"
# 	c.Print(filename)
