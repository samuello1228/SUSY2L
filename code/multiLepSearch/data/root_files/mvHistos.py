#!/usr/bin/env python

import ROOT
from ROOT import TFile

def mvChargeFlipRates():
	outF = TFile("chargeFlipRates_old.root", "recreate")

	fDataSignal = TFile("chargeMisID_Zee_data_signal.root")
	fDataLoose = TFile("chargeMisID_Zee_data_loose.root")
	fMCSignal = TFile("chargeMisID_Zee_MC_signal.root")
	fMCLoose = TFile("chargeMisID_Zee_MC_looseBaseline.root")

	hDataSignal = fDataSignal.Get("hFlipProb")
	hDataSignal.SetName("hFlipProb_data")
	hDataSignal.SetTitle("ChargeMisID from data (signal)")

	hDataLoose = fDataLoose.Get("hFlipProb")
	hDataLoose.SetName("hFlipProb_data_loose")
	hDataLoose.SetTitle("ChargeMisID from data (loose)")


	hMCSignal = fMCSignal.Get("hFlipProb")
	hMCSignal.SetName("hFlipProb_MCLH")
	hMCSignal.SetTitle("ChargeMisID from MCLH (signal)")

	hMCTruthSignal = fMCSignal.Get("hMCFlip")
	hMCTruthSignal.SetName("hFlipProb_MCtruth")
	hMCTruthSignal.SetTitle("ChargeMisID from MC Truth (signal)")

	hMCLoose = fMCLoose.Get("hFlipProb")
	hMCLoose.SetName("hFlipProb_MCLH_loose")
	hMCLoose.SetTitle("ChargeMisID from MCLH (loose)")

	hMCTruthLoose = fMCLoose.Get("hMCFlip")
	hMCTruthLoose.SetName("hFlipProb_MCtruth_loose")
	hMCTruthLoose.SetTitle("ChargeMisID from MC Truth (loose)")

	hDataSignal.SetDirectory(outF)
	hDataLoose.SetDirectory(outF)
	hMCSignal.SetDirectory(outF)
	hMCTruthSignal.SetDirectory(outF)
	hMCLoose.SetDirectory(outF)
	hMCTruthLoose.SetDirectory(outF)

	outF.Write()
	outF.Close()

	fDataSignal.Close()
	fDataLoose.Close()
	fMCSignal.Close()
	fMCLoose.Close()

def mvDPT():
	outF = TFile("dPThistos_old.root", "recreate")
	inF = TFile("dPT_signal.root")
	keys = inF.GetListOfKeys()

	for k in keys:
		h = inF.Get(k.GetName())
		h.SetDirectory(outF)


	inF.Close()
	inF = TFile("dPT_loose.root")
	keys = inF.GetListOfKeys()

	for k in keys:
		h = inF.Get(k.GetName())
		h.SetName(h.GetName()+"_loose")
		h.SetDirectory(outF)

	outF.Write()
	outF.Close()
	inF.Close()