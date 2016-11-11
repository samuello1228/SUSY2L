/* matchEtaShape.C
 * A root script by Gabriel Gallardo
 * 16 July 2016, based on getMllWeights.cpp
 * Input: UM/HK format NTuples
 * 
 * root -l -b -q matchEtaShape.C+
 */

#include "../common/common.C"
#include <iostream>
#include <iomanip>
#include "TH1D.h"
#include "TLine.h"
#include <fstream>
#include <string>
#include "TLegend.h"
#include "TColor.h"
#include "TSystem.h"

using namespace std;

// CONFIG
enum sel{LOOSEBASELINE, SIGNAL};
const sel eventSelected=LOOSEBASELINE;

// ############ INFRASTRUCTURE ############ //
TFile* outFile = 0;

int nPtBins; 
int nEtaBins;
double* PtEdges;
double* EtaEdges;

evt2l *evt2lTree;
TH2D* hChargeMisID;


double flipProb(const int i){
	int binNumber = hChargeMisID->FindBin(fabs(evt2lTree->leps_eta[i]), evt2lTree->leps_pt[i]);
	return hChargeMisID->GetBinContent(binNumber);
}

int matchEtaShape(string outDir="0721b-loose", string misIDfile="../chargeMisID/0721-MC-loose-80100/Gabriel/output.root"){
	
	// INITIALIZE

	TFile* fMisIDfile = TFile::Open(misIDfile.c_str());
	if(!fMisIDfile){
		cout << "Could not open chargeMisID histogram file " << misIDfile <<  endl;
		return -1;
	}
	hChargeMisID = (TH2D*) fMisIDfile->Get("hFlipProb");
	if(!hChargeMisID){
		cout << "Could not extract chargeMisID histogram" << endl;
		return -2;
	}
	hChargeMisID->SetDirectory(0);
	fMisIDfile->Close();

	nPtBins = hChargeMisID->GetYaxis()->GetNbins(); 
	nEtaBins = hChargeMisID->GetXaxis()->GetNbins();
	PtEdges = const_cast<double*>(hChargeMisID->GetYaxis()->GetXbins()->GetArray());
	EtaEdges = const_cast<double*>(hChargeMisID->GetXaxis()->GetXbins()->GetArray());

	TChain *inList = COMMON::loadData("../common/inFileList-MCdR.txt");
	if(!inList){
		cout << "Input NTuples could not be loaded" << endl;
		return -3;
	} 

	if(!(gSystem->cd(outDir.c_str()))){
		gSystem->mkdir(outDir.c_str());
		if(!gSystem->cd(outDir.c_str())){
			cout << "Could not create directory " << outDir;
			return -4;
		}
	}

	outFile = TFile::Open("etaWeights.root", "recreate");

	TH2D *hEtaSS = new TH2D("hEtaSS", "Observed #eta distribution of same sign pairs;Leading #eta;Subleading #eta", 50, 0, 2.47, 50, 0, 2.47);
	hEtaSS->Sumw2();
	hEtaSS->SetDirectory(outFile);

	TH2D *hEtaExpSS = (TH2D*) hEtaSS->Clone("hEtaExpSS");
	hEtaExpSS->SetTitle("Predicted #eta distribution of same sign electron pairs");	
	hEtaExpSS->Sumw2();
	hEtaExpSS->SetDirectory(outFile);

	TH2D *hEtaOS = (TH2D*) hEtaSS->Clone("hEtaOS");
	hEtaOS->SetTitle("Observed #eta distribution of opposite sign electron pairs");
	hEtaOS->Sumw2();
	hEtaOS->SetDirectory(outFile);

	TH2D *hEtaSF = (TH2D*) hEtaSS->Clone("hEtaSF");
	hEtaSF->SetTitle("#eta scale factors");
	hEtaSF->Sumw2();
	hEtaSF->SetDirectory(outFile);

	// LOOP OVER ALL EVENTS

	evt2lTree = new evt2l(inList);
	long long nEntries = evt2lTree->fChain->GetEntries();
	cout << nEntries << " events in TChain" << endl;
	int n = 0;
	for(long long i=0; i<nEntries; i++){
		COMMON::loadbar(i+1, nEntries);
		evt2lTree->GetEntry(i);

		// ----- EVENT SELECTION -------//
		// At least one trigger, passes GRL, and leptons 0 and 1 are both electrons
		bool baseCut = (evt2lTree->sig_trigCode>0
		             && int(abs(evt2lTree->leps_ID[0]/1000))==11 && int(abs(evt2lTree->leps_ID[1]/1000))==11); 
		if (!baseCut) continue;

		// LooseBaseline
		if(!((evt2lTree->leps_lFlag[0] & 1<<0) && (evt2lTree->leps_lFlag[1] & 1<<0))) continue;

		// Signal
		if(eventSelected==SIGNAL && !((evt2lTree->leps_lFlag[0] & 1<<1) && (evt2lTree->leps_lFlag[1] & 1<<1))) continue;

		int charge1 = ((evt2lTree->leps_ID[0]>0) ? 1 : -1);
		int charge2 = ((evt2lTree->leps_ID[1]>0) ? 1 : -1);
		bool isSS = ((charge1 + charge2) != 0);

		double weight = 1;
		weight *= evt2lTree->evt_weight;
		// pileupWeight *= evt_pwt;

		if(isSS){
			hEtaSS->Fill(fabs(evt2lTree->leps_eta[0]), fabs(evt2lTree->leps_eta[1]), weight);
		}
		else {
			double flip1 = flipProb(0);
			double flip2 = flipProb(1);
			double pSS=flip1+flip2-2*flip1*flip2;
			double flipWeight=pSS/(1-pSS);

			hEtaExpSS->Fill(fabs(evt2lTree->leps_eta[0]), fabs(evt2lTree->leps_eta[1]), flipWeight*weight);
			hEtaOS->Fill(fabs(evt2lTree->leps_eta[0]), fabs(evt2lTree->leps_eta[1]), weight);
		}
		n++;
	}
	cout << endl << n << endl;
	// CALCULATE SCALE FACTORS
	hEtaSF->Add(hEtaSS);
	hEtaSF->Divide(hEtaExpSS);

	outFile->Write();

	return 0;
}