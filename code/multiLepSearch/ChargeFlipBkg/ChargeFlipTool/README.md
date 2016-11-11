# README: ChargeFlipTool

ChargeFlipTool
- Written by Gabriel Gallardo Aug 2016 as part of the [ChargeFlipBkg](../) code package.
- Provides interfaces to access misID rates and corrected pT

Interfaces are self-explanatory:

## typedef
### ChargeFlipTool::ResErrorPair
- A `std::vector<double, double>`
- Used to return the result in the first element and the statistical uncertainty in the second element.

## Public functions
### ChargeFlipTool()
Default constructor

### ChargeFlipTool(const TString pathToMisIDfile, const TString pathToDPTfile, const TString pathToMCMisIDfile="")
Initialize ChargeFlipTool with specified files

### ~ChargeFlipTool()
Destructor

### bool loadMisIDfile(const TString pathToFile, const TString histName="hFlipProb", const TString histSys="")
Load misIDfile. Can specify name of histogram in file as well. 

### bool loadMCMisIDfile(const TString pathToFile, const TString histName="hFlipProb", const TString histSys="")
Load MC misID file used for calculation of data/MC scale factors. 

### bool loadDPTfile(const TString pathToFile, const TString histOk="hDPTok_pxy", const TString histFlipped="hDPTflipped_pxy", const TString histOkSys="", const TString histFlippedSys="")
Load dPT file used for pt correction

### ResErrorPair getDataMCSF(const double eta, const double pt) const
Get data/MC scale factor rate_(data)/rate_(MC) (for using MC to predict data) 

### ResErrorPair getChargeFlipRate(const double eta, const double pt, const bool MC=false) const
Get the charge flip rate

### ResErrorPair getCorrectedPt(const double eta, const double pt) const; 
Get the corrected Pt = Pt_reco + dPt_flipped - dPt_ok