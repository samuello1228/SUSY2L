"""
 **********************************************************************************
 * Project: HistFitter - A ROOT-based package for statistical data analysis       *
 * Package: HistFitter                                                            *
 *                                                                                *
 * Description:                                                                   *
 *      Minimal example configuration with two different uncertainties            * 
 *                                                                                *
 * Authors:                                                                       *
 *      HistFitter group, CERN, Geneva                                            *
 *                                                                                *
 * Redistribution and use in source and binary forms, with or without             *
 * modification, are permitted according to the terms listed in the file          *
 * LICENSE.                                                                       *
 **********************************************************************************
"""

################################################################
## In principle all you have to setup is defined in this file ##
################################################################
from configManager import configMgr
from ROOT import kBlack,kWhite,kGray,kRed,kPink,kMagenta,kViolet,kBlue,kAzure,kCyan,kTeal,kGreen,kSpring,kYellow,kOrange
from configWriter import fitConfig,Measurement,Channel,Sample
from systematic import Systematic
from math import sqrt

import os

# Setup for ATLAS plotting
from ROOT import gROOT
#gROOT.LoadMacro("./macros/AtlasStyle.C")
import ROOT
#ROOT.SetAtlasStyle()

##########################

# Set observed and expected number of events in counting experiment
ndata     =  8.3 	# Number of events observed in data

nbkg_VV      = 2.5	# Number of predicted bkg events
nbkg_VV_Err  = 0.6

nbkg_Vg      = 1.8
nbkg_Vg_Err  = 1.0

nbkg_cF      = 0
nbkg_cF_Err  = 0

nbkg_FK      = 4.0
nbkg_FK_Err  = 1.7

nsig      =  10.9  	# Number of predicted signal events
nsigErr   =  1.6  	# (Absolute) Statistical error on signal estimate
lumiError = 0.039 	# Relative luminosity uncertainty

# Set uncorrelated systematics for bkg and signal (1 +- relative uncertainties)
ucb = Systematic("ucb", configMgr.weights, 1.2,0.8, "user","userOverallSys")
ucs = Systematic("ucs", configMgr.weights, 1.1,0.9, "user","userOverallSys")

# correlated systematic between background and signal (1 +- relative uncertainties)
corb = Systematic("cor",configMgr.weights, [1.1],[0.9], "user","userHistoSys")
cors = Systematic("cor",configMgr.weights, [1.15],[0.85], "user","userHistoSys")

##########################

# Setting the parameters of the hypothesis test
configMgr.doExclusion=True # True=exclusion, False=discovery
#configMgr.nTOYs=5000
configMgr.calculatorType=2 # 2=asymptotic calculator, 0=frequentist calculator
configMgr.testStatType=3   # 3=one-sided profile likelihood test statistic (LHC default)
configMgr.nPoints=20       # number of values scanned of signal-strength for upper-limit determination of signal strength.

configMgr.writeXML = True
configMgr.keepSignalRegionType = True

##########################

# Give the analysis a name
configMgr.analysisName = "MyUserTest"
configMgr.outputFileName = "results/%s_Output.root"%configMgr.analysisName

# Define cuts
configMgr.cutsDict["UserRegion"] = "1."

# Define weights
configMgr.weights = "1."

# Define samples
bkgSample_VV = Sample("Bkg_VV",kGreen-9)
bkgSample_VV.setStatConfig(True)
bkgSample_VV.buildHisto([nbkg_VV],"UserRegion","cuts",0.5)
bkgSample_VV.buildStatErrors([nbkg_VV_Err],"UserRegion","cuts")
bkgSample_VV.addSystematic(corb)
bkgSample_VV.addSystematic(ucb)

bkgSample_Vg = Sample("Bkg_Vg",kBlue)
bkgSample_Vg.setStatConfig(True)
bkgSample_Vg.buildHisto([nbkg_Vg],"UserRegion","cuts",0.5)
bkgSample_Vg.buildStatErrors([nbkg_Vg_Err],"UserRegion","cuts")
bkgSample_Vg.addSystematic(corb)
bkgSample_Vg.addSystematic(ucb)

bkgSample_cF = Sample("Bkg_cF",kYellow)
bkgSample_cF.setStatConfig(True)
bkgSample_cF.buildHisto([nbkg_cF],"UserRegion","cuts",0.5)
bkgSample_cF.buildStatErrors([nbkg_cF_Err],"UserRegion","cuts")
bkgSample_cF.addSystematic(corb)
bkgSample_cF.addSystematic(ucb)

bkgSample_FK = Sample("Bkg_FK",kRed)
bkgSample_FK.setStatConfig(True)
bkgSample_FK.buildHisto([nbkg_FK],"UserRegion","cuts",0.5)
bkgSample_FK.buildStatErrors([nbkg_FK_Err],"UserRegion","cuts")
bkgSample_FK.addSystematic(corb)
bkgSample_FK.addSystematic(ucb)

sigSample = Sample("Sig",kPink)
sigSample.setNormFactor("mu_Sig",1.,0.,100.)
sigSample.setStatConfig(True)
sigSample.setNormByTheory()
sigSample.buildHisto([nsig],"UserRegion","cuts",0.5)
sigSample.buildStatErrors([nsigErr],"UserRegion","cuts")
sigSample.addSystematic(cors)
sigSample.addSystematic(ucs)

dataSample = Sample("Data",kBlack)
dataSample.setData()
dataSample.buildHisto([ndata],"UserRegion","cuts",0.5)

# Define top-level
ana = configMgr.addFitConfig("SPlusB")
ana.addSamples([bkgSample_VV,bkgSample_Vg,bkgSample_cF,bkgSample_FK,sigSample,dataSample])
ana.setSignalSample(sigSample)

# Define measurement
meas = ana.addMeasurement(name="NormalMeasurement",lumi=1.0,lumiErr=lumiError)
meas.addPOI("mu_Sig")
#meas.addParamSetting("Lumi",True,1)

# Add the channel
chan = ana.addChannel("cuts",["UserRegion"],1,0.5,1.5)
ana.addSignalChannels([chan])

# These lines are needed for the user analysis to run
# Make sure file is re-made when executing HistFactory
if configMgr.executeHistFactory:
    if os.path.isfile("data/%s.root"%configMgr.analysisName):
        os.remove("data/%s.root"%configMgr.analysisName) 
