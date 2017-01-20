#!/usr/bin/env python

import sys, math, os

import ROOT
from ROOT import TMVA, TMath, TString
from ROOT.TMVA import TMVAGlob
from subprocess import call

ROOT.gROOT.SetBatch(True) 

def setTrainHistStyle(h, col):
    h.SetMarkerColor( col );
    h.SetMarkerSize( 0.7 );
    h.SetMarkerStyle( 20 );
    h.SetLineWidth( 1 );
    h.SetLineColor( col );

def runCheck(filename):
    ## Input
    f = ROOT.TFile.Open(filename)

    if f is None:
    	print ("File %s could not be opened" % (filename))

    if not (ROOT.gDirectory.cd("Method_BDT/BDTD")):
    	print "No directory for BDT method found\n"
    	quit()

    ## Get histograms and set styles
    sig = f.Get("Method_BDT/BDTD/MVA_BDTD_S")
    sigOv = f.Get("Method_BDT/BDTD/MVA_BDTD_Train_S")
    bkg = f.Get("Method_BDT/BDTD/MVA_BDTD_B")
    bkgOv = f.Get("Method_BDT/BDTD/MVA_BDTD_Train_B")

    TMVAGlob.SetSignalAndBackgroundStyle(sig, bkg)
    TMVAGlob.NormalizeHists(sig, bkg)
    TMVAGlob.NormalizeHists(sigOv, bkgOv)
    setTrainHistStyle(sigOv, sig.GetLineColor())
    setTrainHistStyle(bkgOv, bkg.GetLineColor())

    ## Frame
    c = ROOT.TCanvas("canvas", "Overtraining check", 800, 600)

    ROOT.gStyle.SetOptStat(0)

    xmin = TMath.Max(TMath.Min(sig.GetMean() - 10*sig.GetRMS(), 
                                bkg.GetMean() - 10*bkg.GetRMS() ),
                    sig.GetXaxis().GetXmin() )
    xmax = TMath.Min( TMath.Max(sig.GetMean() + 10*sig.GetRMS(), 
                                bkg.GetMean() + 10*bkg.GetRMS() ),
                    sig.GetXaxis().GetXmax() )
    ymin = 0
    ymax = TMath.Max( sig.GetMaximum(), bkg.GetMaximum() )*1.3
    frame = ROOT.TH2F("ovCheck", "Overtraining check", 500, xmin, xmax, 500, ymin, ymax)
    frame.GetXaxis().SetTitle("BDTD response")
    frame.GetYaxis().SetTitle("(1/N) dN^{ }/^{ }dx")

    frame.Draw()

    ## Legends
    legend=ROOT.TLegend( c.GetLeftMargin(), 1 - c.GetTopMargin() - 0.12, 
                           c.GetLeftMargin() + 0.35, 1 - c.GetTopMargin() );
    legend.SetFillStyle( 1 );
    lSigEntry = "Signal"
    lBkgEntry = "Background"
    legend.AddEntry(sig, lSigEntry + " (test sample)", "F");
    legend.AddEntry(bkg, lBkgEntry + " (test sample)", "F");
    legend.SetBorderSize(1);
    legend.SetMargin(0.2);
    legend.Draw("same");

    legend2= ROOT.TLegend( 1 - c.GetRightMargin() - 0.35, 1 - c.GetTopMargin() - 0.12,
                                  1 - c.GetRightMargin(), 1 - c.GetTopMargin() );
    legend2.SetFillStyle( 1 );
    legend2.SetBorderSize(1);
    legend2.AddEntry(sigOv,"Signal (training sample)","P");
    legend2.AddEntry(bkgOv,"Background (training sample)","P");
    legend2.SetMargin( 0.1 );
    legend2.Draw("same");

    sig.Draw("samehist")
    sigOv.Draw("e1 same")
    bkg.Draw("samehist")
    bkgOv.Draw("e1 same")

    ## Tests
    kolS = kolB = -1
    if math.isnan(sig.ComputeIntegral(True)) or math.isnan(sigOv.ComputeIntegral(True)):
        kolS = -1
    else:
        kolS = sig.KolmogorovTest(sigOv)

    if math.isnan(bkg.ComputeIntegral(True)) or math.isnan(bkgOv.ComputeIntegral(True)):
        kolB = -1
    else:
        kolB = bkg.KolmogorovTest(bkgOv)

    print ("Kolmogorov-Smirnov test values for signal (background): %1.3f (%1.3f)" % (kolS, kolB))

    chi2S = sig.Chi2Test( sigOv, "WW");
    chi2B = bkg.Chi2Test( bkgOv, "WW");

    print ("Chi2 test values for signal (background): %1.3f (%1.3f)" % (chi2S, chi2B))

    tHeader = ROOT.TText(0.165, 0.74, "Sig\t\t\t\t\t\t\t\t\tBkg")
    tHeader.SetNDC()
    tHeader.SetTextSize(0.032)
    tHeader.AppendPad()

    tTableText = "#splitline{KS:\t%1.3f \t%1.3f}{#chi^{2}: \t\t%1.3f \t%1.3f}" % (kolS, kolB, chi2S, chi2B)
    tTable = ROOT.TLatex(0.12, 0.69, tTableText)
    tTable.SetNDC()
    tTable.SetTextSize(0.032)
    tTable.AppendPad()

    ## Print to EPS
    c.Print("plots/overtrain_BDTD.eps")
    # c.Print("plots/overtrain_BDTD.pdf")

    ## Save to temp text file
    return (kolS, kolB,chi2S,chi2B)

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print "Put file as argument"
        quit()
    f = open("trainingtest.csv", 'w')
    f.write("%1.3f,%1.3f,%1.3f,%1.3f" % runCheck(sys.argv))
    f.close()