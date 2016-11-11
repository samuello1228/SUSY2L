#!/bin/env python2
#from QFramework import *
import ROOT
from ROOT import *
from ROOT import Double
from array import array
import time
import os
from os import path
import math

##########################################################################
# Giulia Gonella <giulia.gonella@cern.ch>                                #
#                                                                        #
# Macro for building the input file for ElectronChargeCorrectionTool,    #
# with the two TH2 histograms with charge mis-ID rates for DATA and MC   #
##########################################################################

PT  = [25.0, 60.0, 90.0, 130.0, 150.0]
ETA = [0.0, 0.5, 1, 1.37, 1.52, 1.8, 2.0, 2.47]
bins_pt  = array( 'f', PT)
bins_eta = array( 'f', ETA)                                                                   

eleID = "LooseBaseline"
eleISO = ""


#===== Defined functions===============

def SetHistoNicely(h):

    h.GetXaxis().SetTitle("p_{T} GeV")
    h.GetYaxis().SetTitle("|#eta|")
    h.GetYaxis().SetTitleSize(0.05)
    h.GetXaxis().SetTitleSize(0.05)
    h.GetXaxis().SetTitleOffset(0.8)
    return

#=============================================


datafile = TFile("DATA_80_100_20_20.root") # change here for the file you want to use to build the final root file
mcfile   = TFile("MC_80_100_0_0.root")     # change here for the file you want to use to build the final root file

datarates = datafile.Get("Rates_DATA")
mcrates   = mcfile.Get("Rates_MC")


newfile = ROOT.TFile("ChMisIDRates"+eleID+eleISO+".root", "recreate")

#==== Getting histograms =====

#--------------SF_Flip

SF = datarates.Clone()
SF.Divide(mcrates)


#--------------data_Flip_stat 

hdata_stat = TH2F("data_stat", "data_stat" , len(PT)-1, bins_pt, len(ETA)-1, bins_eta)

for i in range(len(bins_pt)-1):
    for j in range(len(bins_eta)-1):
        Nbin  = datarates.GetBin(i+1,j+1)
        data_errstat = datarates.GetBinError(Nbin)
        hdata_stat.SetBinContent(Nbin, data_errstat)
       

SetHistoNicely(hdata_stat) 


#--------------MC_Flip_stat 

hmc_stat = TH2F("mc_stat", "mc_stat" , len(PT)-1, bins_pt, len(ETA)-1, bins_eta)

for i in range(len(bins_pt)-1):
    for j in range(len(bins_eta)-1):
        Nbin  = mcrates.GetBin(i+1,j+1)
        mc_errstat = mcrates.GetBinError(Nbin)
        hmc_stat.SetBinContent(Nbin, mc_errstat)

SetHistoNicely(hmc_stat)

#--------------data_Flip_sys

### !! here axis of input histos are inverted: pay attention when getting the values and filling the new histo!!

data_Sysmll_file = TFile("Systematics/Sys_mll_DATA.root") 
data_Sysbkg_file = TFile("Systematics/Sys_bkgsub.root") 

hdata_syst_mll = data_Sysmll_file.Get("VariationsHisto_DATA")
hdata_syst_bkg = data_Sysbkg_file.Get("VariationsHisto_DATA")

hdata_syst = TH2F("mc_stat", "mc_stat" , len(PT)-1, bins_pt, len(ETA)-1, bins_eta)

for i in range(len(bins_pt)-1):
    for j in range(len(bins_eta)-1):
        Nbin = hdata_syst_mll.GetBin(j+1,i+1)
        mll_err = hdata_syst_mll.BinContent(Nbin)
        bkg_err = hdata_syst_bkg.BinContent(Nbin)
        total_sys = math.sqrt(pow(mll_err, 2) + pow(bkg_err, 2))
        NbinNew = hdata_syst.GetBin(i+1,j+1)
        hdata_syst.SetBinContent(NbinNew, total_sys)

SetHistoNicely(hdata_syst)


#--------------SF_Flip_sys

SF_Syst_file = TFile("Syst_SF.root")

hSF_syst = SF_Syst_file.Get("VariationsHisto_SF")



#==== Writing to file =====

datarates.Write("data_Flip")
mcrates.Write("MC_Flip")
SF.Write("SF_Flip")
datastat.Write("data_Flip_stat")
mcstat.Write("MC_Flip_stat")
hdata_syst.Write("data_Flip_sys")
hSF_syst.Write("SF_Flip_sys")


newfile.Close()
