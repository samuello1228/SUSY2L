"""
 **********************************************************************************
 Script used to test limit setting using HistFitter.
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

import re

# Setup for ATLAS plotting
from ROOT import gROOT
#gROOT.LoadMacro("./macros/AtlasStyle.C")
import ROOT
#ROOT.SetAtlasStyle()

##########################
class ParmV:
    def __init__(self, text):
        vs = text.split('+/-')
        self.e2 = 0
        try:
            self.val = float(vs[0])
        except TypeError:
            self.val = 0
        if len(vs)>1:
            try:
                self.e2 = pow(float(vs[1]),2)
            except TypeError:
                self.e2 = 0
    def add(self, value2):
        self.val += value2.val
        self.e2 += value2.e2
    def v(self):
        return self.val
    def e(self):
        return sqrt(self.e2)
    def __str__(self):
        return "%.2f+/-%.2f" % (self.val, sqrt(self.e2))

    def __repr__(self):
        return "%.2f+/-%.2f" % (self.val, sqrt(self.e2))


# Set observed and expected number of events in counting experiment
nbkg_tot     =  [] 	# Number of events observed in data
nbkg_tot_Err     =  [] 	# Number of events observed in data
nbkg_VV      = []	# Number of predicted bkg events
nbkg_VV_Err  = []
nbkg_Vg      = []
nbkg_Vg_Err  = []
nbkg_cF      = []
nbkg_cF_Err  = []
nbkg_FK      = []
nbkg_FK_Err  = []
nsig_20      =  []  	# Number of predicted signal events
nsig_20_Err   =  []  	# (Absolute) Statistical error on signal estimate
nsig_50      =  []
nsig_50_Err   =  []
nsig_100      =  []
nsig_100_Err   =  []
nsig_200      =  []
nsig_200_Err   =  []
nsig_300      =  []
nsig_300_Err   =  []
rgs = []

dir1 = os.path.dirname(os.environ['ROOTCOREBIN'])
with open(dir1+"/plotScripts/dataMCcomparison/latex/data/expN/SR.txt") as f1:
    for line in f1.readlines()[1:]:
        line = line.replace("+/- ","+/-")
        fs = line.split()
        # ['ISR_mT_100_inf_ptll_no_cut_MET_150_inf_emu', '3.6+/-0.7', '0.6+/-0.3', '0', '2.2+/-0.7', '6.5+/-1.1', '0.6+/-0.2', '-0.049', '6', '2.2+/-0.5', '0.421', '17', '3.9+/-0.9', '0.871', '18', '4.0+/-1.0', '0.908', '15', '3.9+/-1.1', '0.887', '13']
        pVV = ParmV(fs[1])
        pVg = ParmV(fs[2])
        pcF = ParmV(fs[3])
        pFK = ParmV(fs[4])
        pTot = ParmV(fs[5]) # total as data
        pSig20 = ParmV(fs[6])
        pSig50 = ParmV(fs[9])
        pSig100 = ParmV(fs[12])
        pSig200 = ParmV(fs[15])
        pSig300 = ParmV(fs[18])

        nbkg_tot.append(pTot.v())
        nbkg_tot_Err.append(pTot.e())
        nbkg_VV.append(pVV.v())
        nbkg_VV_Err.append(pVV.e())
        nbkg_Vg.append(pVg.v())
        nbkg_Vg_Err.append(pVg.e())
        nbkg_cF.append(pcF.v())
        nbkg_cF_Err.append(pcF.e())
        nbkg_FK.append(pFK.v())
        nbkg_FK_Err.append(pFK.e())
        nsig_20.append(pSig20.v())
        nsig_20_Err.append(pSig20.e())
        nsig_50.append(pSig50.v())
        nsig_50_Err.append(pSig50.e())
        nsig_100.append(pSig100.v())
        nsig_100_Err.append(pSig100.e())
        nsig_200.append(pSig200.v())
        nsig_200_Err.append(pSig200.e())
        nsig_300.append(pSig300.v())
        nsig_300_Err.append(pSig300.e())

        rgs.append(line[0])

ndata = nbkg_tot
nsig = nsig_100
nsigErr = nsig_100_Err

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
configMgr.analysisName = "MyCombinedTest"
configMgr.outputFileName = "results/%s_Output.root"%configMgr.analysisName

# Define cuts
configMgr.cutsDict["UserRegion"] = "1."

# Define weights
configMgr.weights = "1."

var = 'SRi'

# Define samples
bkgSample_VV = Sample("Bkg_VV",kGreen-9)
bkgSample_VV.setStatConfig(True)
bkgSample_VV.buildHisto(nbkg_VV,"UserRegion",var,0.5)
bkgSample_VV.buildStatErrors(nbkg_VV_Err,"UserRegion",var)
# bkgSample_VV.addSystematic(corb)
# bkgSample_VV.addSystematic(ucb)

bkgSample_Vg = Sample("Bkg_Vg",kBlue)
bkgSample_Vg.setStatConfig(True)
bkgSample_Vg.buildHisto(nbkg_Vg,"UserRegion",var,0.5)
bkgSample_Vg.buildStatErrors(nbkg_Vg_Err,"UserRegion",var)
# bkgSample_Vg.addSystematic(corb)
# bkgSample_Vg.addSystematic(ucb)

bkgSample_cF = Sample("Bkg_cF",kYellow)
bkgSample_cF.setStatConfig(True)
bkgSample_cF.buildHisto(nbkg_cF,"UserRegion",var,0.5)
bkgSample_cF.buildStatErrors(nbkg_cF_Err,"UserRegion",var)
# bkgSample_cF.addSystematic(corb)
# bkgSample_cF.addSystematic(ucb)

bkgSample_FK = Sample("Bkg_FK",kCyan)
bkgSample_FK.setStatConfig(True)
bkgSample_FK.buildHisto(nbkg_FK,"UserRegion",var,0.5)
bkgSample_FK.buildStatErrors(nbkg_FK_Err,"UserRegion",var)
# bkgSample_FK.addSystematic(corb)
# bkgSample_FK.addSystematic(ucb)

sigSample = Sample("Sig",kPink)
sigSample.setNormFactor("mu_Sig",1.,0.,100.)
sigSample.setStatConfig(True)
sigSample.setNormByTheory()
sigSample.buildHisto(nsig,"UserRegion",var,0.5)
sigSample.buildStatErrors(nsigErr,"UserRegion",var)
# sigSample.addSystematic(cors)
# sigSample.addSystematic(ucs)

dataSample = Sample("Data",kBlack)
dataSample.setData()
dataSample.buildHisto(ndata,"UserRegion",var,0.5)

# Define top-level
ana = configMgr.addFitConfig("SPlusB")
ana.addSamples([bkgSample_VV,bkgSample_Vg,bkgSample_cF,bkgSample_FK,sigSample,dataSample])
ana.setSignalSample(sigSample)

# Define measurement
meas = ana.addMeasurement(name="NormalMeasurement",lumi=1.0,lumiErr=lumiError)
meas.addPOI("mu_Sig")
#meas.addParamSetting("Lumi",True,1)

# Add the channel
# chan = ana.addChannel("cuts",["UserRegion"],54,0.5,54.5)
chan = ana.addChannel(var,["UserRegion"],54,0.5,54.5)
# chan = ana.addChannel("cuts",["UserRegion"]*54,54,0.5,54.5)
# chan = ana.addChannel("cuts",["SR"+str(i) for i in range(54)],54,0.5,54.5)
ana.addSignalChannels([chan])

# These lines are needed for the user analysis to run
# Make sure file is re-made when executing HistFactory
if configMgr.executeHistFactory:
    if os.path.isfile("data/%s.root"%configMgr.analysisName):
        os.remove("data/%s.root"%configMgr.analysisName) 
