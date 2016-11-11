/* saveSSprediction.C
 * Gabriel Gallardo 25 July 2016
 * Rewrite of getSSPrediction.C
 *
 * Saves chargeMisID rates and dPt into a TTree
 * Add the TTree as a friend of the evt2l TTree when making histograms in a separate piece of code
 *
 * Consider: Save in separate files?
 */
#include <iostream>
#include <fstream>
#include <utility> // for std::pair

#include <TString.h>
#include <TH2D.h>
#include <TTree.h>
#include <TFile.h>
#include <TSystem.h>
#include <TChain.h>

#include "../common/evt2l.C"

// ========= CONFIGURATION =========== //
TString defaultOut="outputDir";
TString defaultNTupleList="../common/inFileList-MC.txt";
TString defaultRatesFile="80.000000_100.000000_0.000000_0.000000_MC.root";
TString defaultDPtFile="dEHistos.root";

// TString defaultRatesHistogram="80.0_100.0_0.0_0.0_MC_misid";
TString defaultRatesHistogram="hFlipProb";
TString defaultDPtHistogram="hDPTflipped_pxy";

bool onlySignal = true;

// ========= INFRASTRUCTURE =========== //
struct vars{
   double chargeFlipWeight;

   double e1rate;
   double e1rateErr;
   double e2rate;
   double e2rateErr;
   
   double e1dPt;
   double e1dPtErr;
   double e2dPt;
   double e2dPtErr;
};

TH2 *hChargeMisID=0;
TH2 *hdPt=0;
evt2l *evt2lTree=0;
vars *misIDvars=0;
TTree *misIDtree=0;
TFile *fOut=0;

void saveSSprediction(const TString outputDir, const TString ntupleList, const TString ratesFile, const TString dPtfile, const TString opt);
bool initialize(const TString ouputDir, const TString ntupleList, const TString ratesFile, const TString dPtfile, const TString opt);
TChain* loadData(const TString);
bool passCuts();
std::pair<double, double> getMisIDrate(int i);
std::pair<double, double> getDPt(int i);
void finalize();

// ==================== FUNCTION DEFINITIONS ==================== //
// Progress bar 
inline void loadbar(unsigned int x, unsigned int n, unsigned int w = 50)
{
   if ( (x != n) && (x % (n/100+1) != 0) ) return;

   float ratio    =  x/(float)n;
   unsigned int c =  ratio * w;

   std::cout << setw(3) << (int)(ratio*100) << "% [";
   for (unsigned int x=0; x<c; x++) std::cout << "=";
   for (unsigned int x=c; x<w; x++) std::cout << " ";
   std::cout << "]\r" << flush;
}

// Main
void saveSSprediction(const TString outputDir=defaultOut, 
						const TString ntupleList=defaultNTupleList, 
						const TString ratesFile=defaultRatesFile, 
						const TString dPtfile=defaultDPtFile, 
						const TString opt="")
{
	if(!initialize(outputDir, ntupleList, ratesFile, dPtfile, opt)){
		std::cout << "Exited from initialize()" << endl;
		return;
	}

	std::pair<double, double> res;
	long long nEntries = evt2lTree->fChain->GetEntries();
	for(long long i=0; i<nEntries; i++){
		loadbar(i+1, nEntries);
		evt2lTree->GetEntry(i);

		if(!passCuts()){
			misIDvars->e1rate = -1;
			misIDvars->e1rateErr = -1;
			misIDvars->e2rate = -1;
			misIDvars->e2rateErr = -1;
			misIDvars->e1dPt = -1;
			misIDvars->e1dPtErr = -1;
			misIDvars->e2dPt = -1;
			misIDvars->e2dPtErr = -1;
			misIDvars->chargeFlipWeight = -1;
		} else {
			res = getMisIDrate(1);
			misIDvars->e1rate = res.first;
			misIDvars->e1rateErr = res.second;

			res = getMisIDrate(2);
			misIDvars->e2rate = res.first;
			misIDvars->e2rateErr = res.second;

			res = getDPt(1);
			misIDvars->e1dPt = res.first;
			misIDvars->e1dPtErr = res.second;

			res = getDPt(2);
			misIDvars->e2dPt = res.first;
			misIDvars->e2dPtErr = res.second;

			// pss = e1(1-e2) + e2(1-e1) = e1 + e2 - 2(e1)(e2)
			// flipWeight = pss / (1-pss)
			double pss = misIDvars->e1rate + misIDvars->e2rate - 2*(misIDvars->e1rate)*(misIDvars->e2rate);
			misIDvars->chargeFlipWeight = pss / (1-pss);
		}

		misIDtree->Fill();
	}
	misIDtree->Print();
	misIDtree->Write();

	finalize();
	std::cout << "Exiting saveSSprediction.C" << endl;
}

// Initialize
bool initialize(const TString outputDir, const TString ntupleList, const TString ratesFile, const TString dPtfile, const TString opt)
{
	///////////////////////////////////
	// Get charge misID rate histogram
	///////////////////////////////////
	TFile* fRates = TFile::Open(ratesFile);
	if(!fRates){
		std::cout << "Could not open charge misID file at " << ratesFile << endl;
		return false;
	} 
	hChargeMisID = (TH2*) fRates->Get(defaultRatesHistogram);
	if(!hChargeMisID){
		std::cout << "Could not get charge misID histogram " << defaultRatesHistogram << " from " << ratesFile << endl;
		return false;
	}
	hChargeMisID->SetDirectory(0);
	fRates->Close();

	///////////////////////////////////
	// Get dPt histogram
	///////////////////////////////////
	TFile* fdPt = TFile::Open(dPtfile);
	if(!fdPt){
		std::cout << "Could not open dPt file at " << dPtfile << endl;
		return false;
	}
	hdPt = (TH2*) fdPt->Get(defaultDPtHistogram);
	if(!hdPt){
		std::cout << "Could not get dPt histogram " << defaultDPtHistogram << " from " << dPtfile << endl;
		return false;
	}
	hdPt->SetDirectory(0);
	fdPt->Close();

	///////////////////////////////////
	// Load data
	///////////////////////////////////
	TChain *ch = loadData(ntupleList);
	if(!ch){
		std::cout << "Could not load data from " << ntupleList << endl;
		return false;
	}
	evt2lTree = new evt2l(ch);
	std::cout << evt2lTree->fChain->GetEntries() << " events in input\n";

	///////////////////////////////////
	// Create output .root file
	///////////////////////////////////
	if(!gSystem->cd(outputDir)){
      gSystem->mkdir(outputDir);
      if(!gSystem->cd(outputDir)){
      	std::cout << "Could not create output directory " << outputDir << endl;
      	return false;
      }
	}
	fOut = TFile::Open("misIDdecorations.root", "recreate");
	if(!fOut){
		std::cout << "Could not create output file " << outputDir << "/misIDdecorations.root" << endl;
		return false;
	}

	///////////////////////////////////
	// Setup output TTree
	///////////////////////////////////
	misIDvars = new vars();
	misIDtree = new TTree("misIDvars", "Charge misID variables");
	misIDtree->Branch("e1rate", &(misIDvars->e1rate), "e1rate/D");
	misIDtree->Branch("e1rateErr", &(misIDvars->e1rateErr), "e1rateErr/D");
	misIDtree->Branch("e2rate", &(misIDvars->e2rate), "e2rate/D");
	misIDtree->Branch("e2rateErr", &(misIDvars->e2rateErr), "e2rateErr/D");
	misIDtree->Branch("e1dPt", &(misIDvars->e1dPt), "e1dPt/D");
	misIDtree->Branch("e1dPtErr", &(misIDvars->e1dPtErr), "e1dPtErr/D");
	misIDtree->Branch("e2dPt", &(misIDvars->e2dPt), "e2dPt/D");
	misIDtree->Branch("e2dPtErr", &(misIDvars->e2dPtErr), "e2dPtErr/D");
	misIDtree->Branch("chargeFlipWeight", &(misIDvars->chargeFlipWeight), "chargeFlipWeight/D");
	misIDtree->SetDirectory(fOut);

	///////////////////////////////////
	// Read options
	///////////////////////////////////
	if(opt!="") // Finished if no options to read
	{
		TString key, value;
		int nChar = opt.Length();
		for(int i=0; i<nChar; i++){
			key = value = "";
			while(opt[i]!='=' && opt[i]!='\0'){
				key += opt[i];
				i++;
			} i++;
			while(opt[i]!=',' && opt[i]!='\0'){
				value += opt[i];
				i++;
			}
	
			if(key=="el"){
				if(value=='s') onlySignal = true;
				else if (value=='l') onlySignal = false;
				else {
					std::cout << "Option unrecognized. " << endl;
					return false;
				}
			} else {
				std::cout << "Option unrecognized. " << endl;
				return false;
			}
		}
	}

	if(onlySignal){
		std::cout << "Selecting Signal events" << endl;
	} else {
		std::cout << "Selecting LooseBaseline events" << endl;
	}

	return true;
}

TChain* loadData(TString fileList){
   //return a TChain linked to the data files
   TChain* tc = new TChain("evt2l");

   if (fileList){
      ifstream inFiles(fileList);
      if (inFiles.is_open()){
         for( std::string line; getline( inFiles, line ); ){tc->Add(line.c_str());}
         inFiles.close();
      }
   }
   return tc;
}


/* Three cuts are applied:
 * 1. At least one trigger, exactly two leptons, both leptons are electrons
 * 2. LooseBaseline (all leptons saved are LooseBaseline, no actual cut required here)
 * 3. Signal
 * Function returns false if cuts are not passed.
 */
bool passCuts(){
	bool baseCut = (evt2lTree->sig_trigCode>0 && evt2lTree->leps_ == 2
		     && int(abs(evt2lTree->leps_ID[0]/1000))==11 && int(abs(evt2lTree->leps_ID[1]/1000))==11); 
	if(!baseCut) return false;

   // Select signal events
   if(onlySignal && !(((evt2lTree->leps_lFlag[0] & 2)/2) && ((evt2lTree->leps_lFlag[1] & 2)/2))) return false;

   return true;
}

// Return misID rate from histogram
std::pair<double, double> getMisIDrate(int i)
{
	std::pair<double, double> res;
	int binNum = hChargeMisID->FindBin(fabs(evt2lTree->leps_eta[i-1]), evt2lTree->leps_pt[i-1]);
	res.first = hChargeMisID->GetBinContent(binNum);
	res.second = hChargeMisID->GetBinError(binNum);
	return res;
}

// Return dPt from histogram
std::pair<double, double> getDPt(int i)
{
	std::pair<double, double> res;
	int binNum = hdPt->FindBin(evt2lTree->leps_pt[i-1], fabs(evt2lTree->leps_eta[i-1]));
	res.first = hdPt->GetBinContent(binNum);
	res.second = hdPt->GetBinError(binNum);
	return res;
}

void finalize(){
	delete hChargeMisID;
	delete hdPt;
	delete evt2lTree;
	delete misIDvars;
	delete misIDtree;
	delete fOut;
}