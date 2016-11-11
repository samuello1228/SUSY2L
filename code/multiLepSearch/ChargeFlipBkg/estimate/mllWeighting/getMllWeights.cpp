/* getMllWeights.C
 * A root script by Gabriel Gallardo
 * 1 July 2016
 * Input: UM/HK format NTuples
 * 
 * root -l -b -q getMllWeights.C+
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
const sel eventSelected=SIGNAL;

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

int getMllWeights(string outDir="Signal-LeadingPt", string misIDfile="MC-signal-EGamma-noPRW.root"){
	
	// INITIALIZE

	TFile* fMisIDfile = TFile::Open(misIDfile.c_str());
   if(!fMisIDfile){
      cout << "Could not open chargeMisID histogram file " << misIDfile <<  endl;
      return -1;
   }
   hChargeMisID = (TH2D*) fMisIDfile->Get("80.0_100.0_0.0_0.0_MC_misid");
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

	TChain *inList = COMMON::loadData("inFileList-MCdR.txt");
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

	outFile = TFile::Open("mllWeights.root", "recreate");

	TH1D *hMassSS = new TH1D("hMassSS", "Observed invariant mass distribution of same sign electron pairs;m_{ll} (GeV);nEvents/0.5 GeV", 400, 20, 220);
	hMassSS->Sumw2();
	hMassSS->SetDirectory(outFile);

	TH1D *hMassExpSS = (TH1D*) hMassSS->Clone("hMassExpSS");
	hMassExpSS->SetTitle("Predicted invariant mass distribution of same sign electron pairs");	
	hMassExpSS->Sumw2();
	hMassExpSS->SetDirectory(outFile);

	TH1D *hMassOS = (TH1D*) hMassSS->Clone("hMassOS");
	hMassOS->SetTitle("Observed invariant mass distribution of opposite sign electron pairs");
	hMassOS->Sumw2();
	hMassOS->SetDirectory(outFile);

	TH1D *hMassSF = (TH1D*) hMassSS->Clone("hMassSF");
	hMassSF->SetTitle("Invariant mass scale factors;m_{ll} (GeV);SF (Observed/Predicted)");
	hMassSF->Sumw2();
	hMassSF->SetDirectory(outFile);

	// LOOP OVER ALL EVENTS

	evt2lTree = new evt2l(inList);
	long long nEntries = evt2lTree->fChain->GetEntries();
	cout << nEntries << " events in TChain" << endl;
	for(long long i=0; i<nEntries; i++){
		COMMON::loadbar(i+1, nEntries);
		evt2lTree->GetEntry(i);

		// ----- EVENT SELECTION -------//
      // At least one trigger, passes GRL, and leptons 0 and 1 are both electrons
      bool baseCut = (evt2lTree->sig_trigCode>0 && evt2lTree->evt_cuts==1
                     && int(abs(evt2lTree->leps_ID[0]/1000))==11 && int(abs(evt2lTree->leps_ID[1]/1000))==11); 
      if (!baseCut) continue;

      // LooseBaseline
      if(!((evt2lTree->leps_lFlag[0] & 1<<0) && (evt2lTree->leps_lFlag[1] & 1<<0))) continue;

      // Signal
      if(eventSelected==SIGNAL && !((evt2lTree->leps_lFlag[0] & 1<<1) && (evt2lTree->leps_lFlag[1] & 1<<1))) continue;
      
      int charge1 = ((evt2lTree->leps_ID[0]>0) ? 1 : -1);
      int charge2 = ((evt2lTree->leps_ID[1]>0) ? 1 : -1);
      bool isSS = ((charge1 + charge2) != 0);

      double pileupWeight = 1;
      // pileupWeight *= evt_pwt;
      // double mll = max(evt2lTree->leps_pt[0],evt2lTree->leps_pt[1]);
      double mll = evt2lTree->l12_m;


      if(isSS){
      	hMassSS->Fill(mll, pileupWeight);
      }
      else {
      	double flip1 = flipProb(0);
      	double flip2 = flipProb(1);
		   double pSS=flip1+flip2-2*flip1*flip2;
		   double flipWeight=pSS/(1-pSS);

		   hMassExpSS->Fill(mll, flipWeight*pileupWeight);
		   hMassOS->Fill(mll, pileupWeight);
      }
	}
	cout << endl;
	// CALCULATE SCALE FACTORS
	hMassSF->Add(hMassSS);
	hMassSF->Divide(hMassExpSS);

	outFile->Write();

	return 0;
}
