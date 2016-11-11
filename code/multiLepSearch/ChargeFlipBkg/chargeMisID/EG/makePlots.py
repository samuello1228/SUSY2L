#!/bin/env python2
#from QFramework import *
import ROOT
from ROOT import *
from ROOT import Double
from array import array
import time
import os
from os import path

##########################################################################
# Giulia Gonella <giulia.gonella@cern.ch>                                #
#                                                                        #
# Macro for getting the TH2 histo with rates from the fit, change axis   #
# and print it.                                                          #
# Also, separating every pt bin and plotting all together                #
# as a function of eta.                                                  #
# Also, getting the TH2 and every TH1 and save them in a .root file      #
# (input for buildRatesFile.py)                                          #
##########################################################################
                                                                        

#**********************************************************************************************
# WHAT TO CHANGE:

#** Choose the TAG for MC or DATA
# tag = "DATA"
tag = "MC"

#** Choose the mass window
low = "80"
up  = "100"

#** Choose sidebands ranges 
bl = "0"
br = "0"

#** Choose the binning

PT  = [20.0, 60.0, 90.0, 130.0, 1000.0]
ETA = [0, 0.5, 1, 1.37, 1.52, 1.8, 2.0, 2.47]

#**********************************************************************************************



##############################################################################################
##############################################################################################


features = tag+"_"+low+"_"+up+"_"+bl+"_"+br
filename = features


#==== DEFINED FUNCTIONS ===============================================================================================================

#------ To extract 1D plots-------------------------------------------------------

def D1plot(a, binsEta):
    listofhisto = []
    for j in range (a.GetNbinsY()):
        h = TH1F("Rates1D_"+tag+"pt"+str(j+1),"Rates1D_"+tag+"pt"+str(j+1), len(binsEta)-1, binsEta)
        for i in range(a.GetNbinsX()):
            Nbin= a.GetBin(i+1,j+1)
            value = a.GetBinContent(Nbin)
            err = a.GetBinError(Nbin)
            h.SetBinContent(i+1,value)
            h.SetBinError(i+1,err)
        listofhisto.append(h)
    return listofhisto
#----------------------------------------------------------------------------------


#------ To set 1D plots' style-----------------------------------------------------
def setStyle(h, col):
    h.SetLineColor(col)
    h.SetLineWidth(2)
    h.SetMarkerStyle(20)
    h.SetMarkerSize(0.5)
    h.SetMarkerColor(col)
    h.GetYaxis().SetTitleSize(0.05)
    h.GetXaxis().SetTitleSize(0.05)
    return
#----------------------------------------------------------------------------------

#==== END DEFINED FUNCTIONS==============================================================================================================



gStyle.SetPaintTextFormat("4.5f")
gStyle.SetOptStat(000000)


# Axes switching for 2D plot=============================================================================================== 

inputfile = TFile(low+".000000_"+up+".000000_"+bl+".000000_"+br+".000000_"+tag+".root")
originalH = inputfile.Get(low+".0_"+up+".0_"+bl+".0_"+br+".0_"+tag+"_misid")
binsx_pt  = array( 'f', PT)
binsy_eta = array( 'f', ETA)
h = TH2F("Rates_"+tag,"Rates_"+tag, len(PT)-1, binsx_pt, len(ETA)-1, binsy_eta)
#cO = TCanvas("DATA_old", "DATA_old")
c = TCanvas(tag, tag)
    
for i in range(len(binsy_eta)-1):
    for j in range(len(binsx_pt)-1):
        Nbin  = originalH.GetBin(i+1,j+1)
        value = originalH.GetBinContent(Nbin)
        err   = originalH.GetBinError(Nbin)
        Nbinnew = h.GetBin(j+1,i+1)
        h.SetBinContent(Nbinnew, value)
        h.SetBinError(Nbinnew, err)
        #print "From histo ---> value = ", value, " error = ", err

h.GetXaxis().SetTitle("p_{T} GeV")
h.GetYaxis().SetTitle("|#eta|")
h.GetYaxis().SetTitleSize(0.05)
h.GetXaxis().SetTitleSize(0.05)
h.GetXaxis().SetTitleOffset(0.8)
c.cd()
h.SetTitle("Charge mis-ID rates "+features)
h.Draw("COLZTEXTERROR")

#cO.cd()
#a.Draw("COLZTEXT")    
#a.SetTitle("ratesDATA old")


# 1d plotting===============================================================================================================

c1d = TCanvas("1D", "1D")
leg = TLegend(0.2,0.6,0.4,0.8)
leg.SetBorderSize(0)
leg.SetTextSize(0.035)

listOf1Dhistos = D1plot(originalH, binsy_eta)
hpt1 = listOf1Dhistos[0]
hpt2 = listOf1Dhistos[1]
hpt3 = listOf1Dhistos[2]
hpt4 = listOf1Dhistos[3]

setStyle(hpt1, 1)
setStyle(hpt2, 2)
setStyle(hpt3, 3)
setStyle(hpt4, 4)


hpt1.SetMaximum(0.5) #hpt4.GetMaximum()*1.1
hpt1.SetTitle(features)
hpt1.GetXaxis().SetTitle("|#eta|")
hpt1.GetYaxis().SetTitle("#epsilon_{mis-ID}")
hpt1.GetYaxis().SetTitleOffset(0.8)
hpt1.GetXaxis().SetTitleOffset(0.8)


leg.AddEntry(hpt1, "p_{T} #epsilon [25,60] GeV", "l")
leg.AddEntry(hpt2, "p_{T} #epsilon [60,90] GeV", "l")
leg.AddEntry(hpt3, "p_{T} #epsilon [90,130] GeV", "l")
leg.AddEntry(hpt4, "p_{T} #epsilon [130,1000] GeV", "l")


#------ Drawing---------------

c1d.cd()
gPad.SetLogy()
hpt1.Draw("E")
hpt2.Draw("SAMEE")
hpt3.Draw("SAMEE")
hpt4.Draw("SAMEE")
#hpt5.Draw("SAMEE")
leg.Draw("SAME")
c1d.Update()


#----- Saving in the root file---------------

rootfile = ROOT.TFile(filename+".root", "recreate")

h.Write()

hpt1.Write()
hpt2.Write()
hpt3.Write()
hpt4.Write()

rootfile.Close()

#----- Saving as pdf

c.Print("Rates_"+features+".pdf")
c1d.Print("RatesLog_1D_"+features+".pdf")


# time.sleep(5000)
