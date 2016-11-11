#ifndef CHARGEFLIPTOOL
#define CHARGEFLIPTOOL

#include <TH2.h>
#include <TString.h>
#include <utility> // for std::pair

const std::pair<double, double> errPair = std::make_pair(-1,0);

class ChargeFlipTool{
typedef std::pair<double, double> ResErrorPair;

public:
	ChargeFlipTool();
	ChargeFlipTool(const TString pathToMisIDfile, const TString pathToDPTfile, const TString pathToMCMisIDfile="");
	~ChargeFlipTool();

	// Load hMisIDRates
	bool loadMisIDfile(const TString pathToFile, 
					const TString histName="hFlipProb",
					const TString histSys="");

	// Load hMCRates for Data/MC SF (use on data only)
	bool loadMCMisIDfile(const TString pathToFile, 
					const TString histName="hFlipProb",
					const TString histSys="");

	// Load dPT file
	bool loadDPTfile(const TString pathToFile, 
					const TString histOk="hDPTok_pxy", 
					const TString histFlipped="hDPTflipped_pxy",
					const TString histOkSys="",
					const TString histFlippedSys="");

	// Get data/MC scale factor [rate(data)/rate(MC)] (for using MC to predict data) 
	ResErrorPair getDataMCSF(const double eta, const double pt) const;

	// Get charge flip rate
	ResErrorPair getChargeFlipRate(const double eta, const double pt, const bool MC=false) const;

	// Get corrected pt
	ResErrorPair getCorrectedPt(const double eta, const double pt) const;	


private:
	TH2*	hMisIDRates;
	TH2*	hMCRates;

	TH2*	hdPTok;
	TH2*	hdPTflipped;

	ResErrorPair getDPT(const double eta, const double pt, const bool chargeFlipped) const;

};

#endif