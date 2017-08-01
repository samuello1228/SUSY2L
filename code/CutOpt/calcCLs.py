#!/usr/bin/env python

##################
# calcCLs.py
#
# Gabriel Gallardo Mar 2017
# Replace significance with CLs
#    CLs = p(s+b)/[1-p(b)]
#
##################

import ROOT
from ROOT import TMath, Math
import csv
import collections
# from TMath import Poisson

# nSig = 112.83
# n
# nBkg = 113.71
# nData = nBkg

def getCLs(nSig, nBkg, nData):
	nData = float(nData)
	nSig = float(nSig)
	nBkg = float(nBkg)
	# if nBkg>=1: return TMath.PoissonI(nData, nSig+nBkg) / TMath.PoissonI(nData, nBkg) 
	if nBkg>=1: return Math.poisson_cdf(int(nData), nSig+nBkg) / Math.poisson_cdf(int(nData), nBkg)
	else: return 2.

if __name__ == "__main__":
	with open("checksSig.csv","r") as inFile:
		with open("checksCLs.csv","w") as outFile:
			start = inFile.tell()
			outFile.write(inFile.readline().replace("sigma", "CLs"))
			inFile.seek(start)
			r = csv.DictReader(inFile)
			# outFile.write("Channel  ISR  Flavor  NTrees  NodeSize  Depth  xSec      nSig  nBkg     BDTopt  nSig(BDTopt)  nBkg(BDTopt)  sigma(BDTopt)  nSig(0)  nBkg(0) sigma(0)  nSig(0.1)  nBkg(0.1)  sigma(0.1)  nSig(0.2)  nBkg(0.2)  sigma(0.2)  nSig(0.3)  nBkg(0.3)  sigma(0.3)")
			for row in r:
				# outFile.write(",".join("%s" % el for el in ("Romeo", "Juliet")))
				sigOpt = getCLs(row["nSig(BDTopt)"], row["nBkg(BDTopt)"], row["nBkg(BDTopt)"]) 
				sig1 = getCLs(row["nSig(0.1)"], row["nBkg(0.1)"], row["nBkg(0.1)"])
				sig2 = getCLs(row["nSig(0.2)"], row["nBkg(0.2)"], row["nBkg(0.2)"])
				sig3 = getCLs(row["nSig(0.3)"], row["nBkg(0.3)"], row["nBkg(0.3)"])
				sig0 = getCLs(row["nSig(0)"], row["nBkg(0)"], row["nBkg(0)"])

				outFile.write(",".join("%s" % el for el in (row["m(C1)"], row["m(N1)"], row["Channel"], row["ISR"], row["Flavor"], row["NTrees"], row["NodeSize"], row["Depth"], row["nSig"], row["nBkg"], row["BDTopt"], row["nSig(BDTopt)"], row["nBkg(BDTopt)"], sigOpt, row["nSig(0)"], row["nBkg(0)"], sig0, row["nSig(0.1)"], row["nBkg(0.1)"], sig1, row["nSig(0.2)"], row["nBkg(0.2)"] , sig2, row["nSig(0.3)"], row["nBkg(0.3)"] ,sig3)))#
				outFile.write("\n")