#!/usr/bin/env python

import os, sys, re, operator, overtraining, tqdm, subprocess, csv
from os import listdir, path, mkdir, rename
from os.path import isdir, isfile
from re import search
from subprocess import call
from tqdm import tqdm
import ROOT
from ROOT.Math import poisson_cdf as PP

OTtolerance = 0.1

C1masses = (200, 300, 400, 500, 600, 700, 800, 900, 1000)
luminosity = 33257.2
sigFilesTxt = "CutOpt/GabrielFiles/sigFiles_all.txt"

headersOT = [
	"dm",
	"channel",
	"NTrees",
	"MNS",
	"Depth",
	"ROC-int",
	"passOT",
	"R@B=0.01",
	"R@B=0.10",
	"R@B=0.30",
	"KSsig",
	"KSbkg",
	"Chi2Sig",
	"Chi2Bkg"
	]

headersSig = [
	"m(C1)","m(N1)","dm", "channel","NTrees","MNS","Depth","ROC-int",
	"nSig","nBkg",
	"BDTopt", "CLs(BDTopt)",
	"nSig(0)","nBkg(0)","CLs(0)",
	"nSig(0.1)","nBkg(0.1)","CLs(0.1)",
	"nSig(0.2)","nBkg(0.2)","CLs(0.2)",
	"nSig(0.3)","nBkg(0.3)","CLs(0.3)"
	]

'''
TODO:
- makeSampleIdDb
- getNDict
- getTotDict
'''

def makeSampleIdDB():
	sampleIdDB = {}
	fileList = open(sigFilesTxt,"r").readlines()

	for aLine in fileList:
		if "MGPy8EG_A14N23LO_C1N2_Slep" not in aLine: continue

		sampleID = int((re.search(".[0-9]{6}.", aLine).group()).split('.')[1])
		masses = re.search("Slep_[0-9]{3,4}_[0-9]{1,4}", aLine).group().split('_')
		mC1 = int(masses[1])
		mN1 = int(masses[2])

		sampleIdDB[(mC1, mN1)] = sampleID
	return sampleIdDB

def getTotDict(nDict):
	nTotDict = {}
	channels = nDict[nDict.keys()[0]].keys()

	for chan in channels: 
		nTotDict[chan] = 0

	for sampleID in nDict:
		for chan in channels:
			nTotDict[chan] += nDict[sampleID][chan]

	return nTotDict
	
def getNDict(dictCSV):
	allLines = open(dictCSV).readlines()

	headers = allLines[0].split('\n')[0].split(',')
	nCols = len(headers)

	nDict = {}
	for row in allLines[1:]:
		entries = row.split('\n')[0].split(',')
		sampleDict = {int(headers[i]): float(entries[i]) for i in range(1,nCols)}
		nDict[int(entries[0])] = sampleDict
	return nDict

#### Main function
if __name__ == '__main__':
	ROOT.gIgnoreErrorLevel = ROOT.kInfo # Quiet INFO messages

	# Load nSig and nBkg
	sampleIdDB = makeSampleIdDB()
	nBkgDict = getNDict("nBkgDict.csv")
	nSigDict = getNDict("nSigDict.csv")
	nBkgTotDict = getTotDict(nBkgDict)

	# Initialize CSVs
	outOT = open("OTnew.csv", "w")
	outOT.write(",".join(headersOT))
	outOT.write("\n")

	outSig = open("SigNew.csv", "w")
	outSig.write(",".join(headersSig))
	outSig.write("\n")

	rowsOT = []
	rowsSig = []

	for directory in sys.argv[1:]:
		for f in listdir(directory):
			if isdir(f): continue
			if not f.endswith("log"): continue
			else: file = f.split(".")[0] + ".root"
			metadata = directory.split("_")
			metadata = metadata + f.split(".")[0].split("_")

			plotDir = "%s/plots" % directory
			if not isdir(plotDir): mkdir(plotDir)

			dm = int(search("[0-9]+", metadata[4]).group(0))
			chan = int(search("[0-9]+", metadata[5]).group(0))
			ntrees = int(metadata[1])
			mns = int(search("[0-9]+", metadata[2]).group(0))
			depth = int(search("[0-9]+", metadata[3]).group(0))
			ratios = []
			roc = 0

			findROCnext = False
			passOT = " "
			for line in reversed(open("%s/%s" % (directory, f)).readlines()):
				if line.startswith("--- Factory                  : BDTD           :"):
					if findROCnext:
						roc = round(int(line.split("|")[0].split(":")[2].split(".")[4])/1000., 3)
						break
					else:
						findROCnext = True
						effs = line.split(":")[2].split("      ")
							
						for n in effs:
							n = n.split(")")[0].split("(")
							r = round(float(n[1]) / float(n[0]), 3) # Training / testing, rounded to 3 places
							ratios.append(r) 
						if ( int(abs(ratios[0]-1)< OTtolerance) + int(abs(ratios[1]-1)< OTtolerance) + int(abs(ratios[2]-1)< OTtolerance) ) == 3:
							passOT = "*"

			## KS test
			overtrainingNums = overtraining.runCheck("%s/%s" % (directory, file))

			# Rename / touch KS plots
			overtrainPlotname = "%s/dm%d_Channel%d.eps" % (plotDir, dm, chan)
			if isfile("plots/overtrain_BDTD.eps"):
				os.rename("plots/overtrain_BDTD.eps", overtrainPlotname)
			else: open(overtrainPlotname, 'a').close()						

			rowsOT.append( (dm, chan, ntrees, mns, depth, roc, passOT) + tuple(["%s" % str(r) for r in ratios] ) + overtrainingNums)

			if passOT==" ": continue

			## CLs checks
			for mass in tqdm(C1masses):
				if mass==1000 and dm%100!=0: continue # Only for mC1=1000, (1000, 100n) samples available
				chan1 = chan%100
				nSig = nSigDict.get(sampleIdDB[(mass, mass-dm)],{}).get(chan1,0)
				if nSig < 1: continue
				nBkg = nBkgTotDict[chan1]

				call('root -l -b -q "mvaeffs.cxx(\\"%s/%s\\", %d, %f)"' % (directory, file, int(nSig), int(nBkg)), shell=True)
				
				minCLs = (1, 1)
				with open("effs.csv") as effFile:
					reader = csv.reader(effFile)
					effList = list(reader)[0][6:] # Only extract nSig, nBkg for BDT = 0, 0.1, 0.2, 0.3
					for i in range(0,12,3):
						# CLs formula is PP(nData|S+B)/PP(nBkg|B)
						nS = float(effList[i])
						nB = float(effList[i+1])
						nD = int(nB) # Take nData = nBkg for this MC study
						CLs = PP(nD, nS+nB) / PP(nD, nB)
						effList[i+2] = CLs
						if minCLs[1]>CLs: minCLs = (round(i/30.,1), CLs) 
				rowsSig.append( (mass, mass-dm, dm, chan, ntrees, mns, depth, roc, nSig, nBkg) + minCLs + tuple(effList[:-1])) 
				os.rename("plots/mvaeffs_BDTD.eps", "%s/%d_%d_Channel%d.eps" % (plotDir, mass, mass-dm, chan))

	rowsOT = sorted(rowsOT, key=lambda x : (x[0], x[1], -x[5]) )
	rowsSig = sorted(rowsSig, key=lambda x: (x[2], x[0], x[-1]))


	# Export to CSV
	for line in rowsOT:
		outOT.write(",".join(["%s" % str(l) for l in line]))
		outOT.write("\n")
	for line in rowsSig:
		outSig.write(",".join(["%s" % str(l) for l in line]))
		outSig.write("\n")

	# Plot graph
	c = ROOT.TCanvas()

	h = ROOT.TH2D("minCLs", "Minimum CLs;m(C1) [GeV];m(N1) [GeV]", 50, 100, 1000, 41, 80, 900)

	currentMass = (0, 0)
	for line in rowsSig:
		if currentMass == (line[0], line[1]): continue
		currentMass = (line[0], line[1])
		h.Fill(line[0], line[1], line[-1])
		print (line[0], line[1], line[-1])
	h.Draw("text")

	line = ROOT.TLine()
	l0 = line.DrawLine(100, 100, 900, 900)
	l0.SetLineStyle(1)
	l0.SetLineWidth(1)

	l20 = line.DrawLine(100, 80, 920, 900)
	l20.SetLineStyle(3)
	l20.SetLineWidth(1)

	l50 = line.DrawLine(130, 80, 950, 900)
	l50.SetLineStyle(3)
	l50.SetLineWidth(1)

	l100 = line.DrawLine(180, 80, 1000, 900)
	l100.SetLineStyle(3)
	l100.SetLineWidth(1)

	c.Print("minCLs.pdf")




	
