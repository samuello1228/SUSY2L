#include "ChargeFlipTool.h"
#include <iostream>
#include <math.h>
// Implementation of ChargeFlipTool

#define CFERROR(f, m) \
	std::cout << "ERROR in " << f << "of ChargeFlipTool: " << m << std::endl;

ChargeFlipTool::ChargeFlipTool()
{
	hMisIDRates 	= NULL;
	hMCRates 		= NULL;
	hdPTok 			= NULL;
	hdPTflipped 	= NULL;
}

ChargeFlipTool::ChargeFlipTool(const TString pathToMisIDfile, const TString pathToDPTfile, const TString pathToMCMisIDfile)
{
	ChargeFlipTool();
	if(!loadDPTfile(pathToDPTfile)) 
		CFERROR("ChargeFlipTool()", "Initialization failed.");
	if(!loadMisIDfile(pathToMisIDfile)) 
		CFERROR("ChargeFlipTool()", "Initialization failed.");
	if(pathToMCMisIDfile!="" && !loadMCMisIDfile(pathToMCMisIDfile))
		CFERROR("ChargeFlipTool()", "Initialization failed.");
}

ChargeFlipTool::~ChargeFlipTool()
{
	if(hMisIDRates){
		delete hMisIDRates;
		hMisIDRates = NULL;
	}
	if(hMCRates){
		delete hMCRates;
		hMCRates = NULL;
	}
	if(hdPTok){
		delete hdPTok;
		hdPTok = NULL;
	}
	if(hdPTflipped){
		delete hdPTflipped;
		hdPTflipped = NULL;
	}
}

// Load hMisIDRates
bool ChargeFlipTool::loadMisIDfile(const TString pathToFile, const TString histName, const TString histSys)
{
	TFile *f = TFile::Open(pathToFile);
	if(!f){
		CFERROR("loadMisIDfile()", pathToFile + " not found.")
		return false;
	}

	hMisIDRates = (TH2*) f->Get(histName);
	if(!hMisIDRates){
		CFERROR("loadMisIDfile()", histName + " in " + pathToFile + " not found.")
		return false;
	}

	hMisIDRates->SetDirectory(0);
	f->Close();
	return true;
}

// Load hMCRates for Data/MC SF (use on data only)
bool ChargeFlipTool::loadMCMisIDfile(const TString pathToFile, const TString histName, const TString histSys)
{
	TFile *f = TFile::Open(pathToFile);
	if(!f){
		CFERROR("loadMCMisIDfile()", pathToFile + " not found.")
		return false;
	}

	hMCRates = (TH2*) f->Get(histName);
	if(!hMisIDRates){
		CFERROR("loadMCMisIDfile()", histName + " in " + pathToFile + " not found.")
		return false;
	}

	hMCRates->SetDirectory(0);
	f->Close();
	return true;
}

// Load dPT file
bool ChargeFlipTool::loadDPTfile(const TString pathToFile, const TString histOk, const TString histFlipped, const TString histOkSys, const TString histFlippedSys)
{
	TFile *f = TFile::Open(pathToFile);
	if(!f){
		CFERROR("loadDPTfile()", pathToFile + " not found.")
		return false;
	}

	hdPTok = (TH2*) f->Get(histOk);
	if(!hdPTok){
		CFERROR("loadDPTfile()", histOk + " in " + pathToFile + " not found.")
		return false;
	}

	hdPTflipped = (TH2*) f->Get(histFlipped);
	if(!hdPTflipped){
		CFERROR("loadDPTfile()", histFlipped + " in " + pathToFile + " not found.")
		return false;
	}

	hdPTok->SetDirectory(0);
	hdPTflipped->SetDirectory(0);
	f->Close();
	return true;
}

// Get data/MC scale factor [rate(data)/rate(MC)] (use on data only) 
ChargeFlipTool::ResErrorPair ChargeFlipTool::getDataMCSF(const double eta, const double pt) const
{
	if(!hMCRates){
		CFERROR("getDataMCSF()", "MC histogram not loaded! Call loadMCMisIDfile() before getDataMCSF().");
		return errPair;
	}
	if(!hMisIDRates){
		CFERROR("getDataMCSF()", "Rates histogram not loaded! Call loadMisIDfile() before getDataMCSF().");
		return errPair;
	}

	ChargeFlipTool::ResErrorPair dataPair = getChargeFlipRate(eta, pt);
	ChargeFlipTool::ResErrorPair mcPair = getChargeFlipRate(eta, pt, true);

	double sf = dataPair.first/mcPair.first;
	double uncertainty = fabs(sf)*sqrt(pow(dataPair.second/dataPair.first,2) + pow(mcPair.second/mcPair.first,2));

	return std::make_pair(sf, uncertainty);
}

// Get charge flip rate
ChargeFlipTool::ResErrorPair ChargeFlipTool::getChargeFlipRate(const double eta, const double pt, const bool MC) const
{
	TH2* h = 0;
	if(MC){
		if(!hMCRates){
			CFERROR("getChargeFlipRate()", "MC histogram not loaded! Call loadMCMisIDfile() before getChargeFlipRate(eta, pt, true).");
			return errPair;
		}
		h = hMCRates;
	} else {
		if(!hMisIDRates){
			CFERROR("getChargeFlipRate()", "Rates histogram not loaded! Call loadMisIDfile() before getChargeFlipRate()");
			return errPair;
		}
		h = hMisIDRates;
	}

	int i = h->FindBin(fabs(eta), pt);
	return std::make_pair(h->GetBinContent(i), h->GetBinError(i));
}

// Get corrected pt
ChargeFlipTool::ResErrorPair ChargeFlipTool::getCorrectedPt(const double eta, const double pt) const
{
	if(!hdPTflipped){
		CFERROR("getCorrectedPt()", "hdPTflipped histogram not loaded! Call loadDPTfile() before calling getCorrectedPt().");
		return errPair;
	}
	if(!hdPTok){
		CFERROR("getCorrectedPt()", "hdPTok histogram not loaded! Call loadDPTfile() before calling getCorrectedPt(),");
		return errPair;
	}

	ChargeFlipTool::ResErrorPair okPair = getDPT(eta, pt, false);
	ChargeFlipTool::ResErrorPair flippedPair = getDPT(eta, pt, true);

	if(flippedPair == errPair || okPair == errPair)
	{
		CFERROR("getCorrectedPt()", "Could not retrieve pt correction.");
		return errPair;
	}

	double correctedPt = pt + flippedPair.first - okPair.first;
	// double uncertainty = flippedPair.second;
	double uncertainty = sqrt(pow(flippedPair.second,2) + pow(okPair.second,2));

	return std::make_pair(correctedPt, uncertainty);
}

ChargeFlipTool::ResErrorPair ChargeFlipTool::getDPT(const double eta, const double pt, const bool chargeFlipped) const
{
	TH2* h = 0;
	if(chargeFlipped) h = hdPTflipped;
	else h = hdPTok;

	if(!h) return errPair;

	int i = h->FindBin(pt, fabs(eta));
	return std::make_pair(h->GetBinContent(i), h->GetBinError(i));
}
