#include "TH1D.h"
#include "TFile.h"
#include "TChain.h"
#include <fstream>
#include <iostream>

inline void loadbar(unsigned int x, unsigned int n, unsigned int w = 50)
{
  if ( (x != n) && (x % (n/100+1) != 0) ) return;

  float ratio    =  x/(float)n;
  unsigned int c =  ratio * w;

  cout << setw(3) << (int)(ratio*100) << "% [";
  for (unsigned int x=0; x<c; x++) cout << "=";
  for (unsigned int x=c; x<w; x++) cout << " ";
  cout << "]\r" << flush;
}

struct event{
	float elCand1_cl_eta;
	float elCand2_cl_eta;

	float elCand1_pt;
	float elCand2_pt;

	int elCand1_ID;
	int elCand2_ID;
};

void checkDistros(){
	TH1D* hOldEta = new TH1D("hOldEta", "#eta distribution of electrons;|#eta|", 100, -2.5, 2.5);
	TH1D* hOldLPt = new TH1D("hOldLPt", "Leading p_{T} distribution;p_{T} (GeV)", 200, 20, 200);
	TH1D* hOldSLPt = new TH1D("hOldSLPt", "Subleading p_{T} distribution;p_{T} (GeV)", 200, 20, 200);

	TH1D* hNewEta = (TH1D*) hOldEta->Clone("hNewEta");
	TH1D* hNewLPt = (TH1D*) hOldLPt->Clone("hNewLPt");
	TH1D* hNewSLPt = (TH1D*) hOldSLPt->Clone("hNewSLPt");

	TFile *f = TFile::Open("distros.root", "recreate");

	hOldEta->SetDirectory(f);
	hOldLPt->SetDirectory(f);
	hOldSLPt->SetDirectory(f);
	hNewEta->SetDirectory(f);
	hNewLPt->SetDirectory(f);
	hNewSLPt->SetDirectory(f);

	TChain* oldChain = new TChain("ZeeCandidate");
	ifstream oldFilelist("../common/inFileList-MCconverted-old.txt");

	string filename;
	while (!oldFilelist.eof()){
		oldFilelist >> filename;
		oldChain->Add(filename.c_str());
		cout << "Added " << filename << " to `oldChain`" << endl;
	}

	event evt;
	oldChain->SetBranchStatus("*",0);
	#define CONNECT(b) evt.b = 0; oldChain->SetBranchStatus(#b,1); oldChain->SetBranchAddress(#b,&(evt.b)); 
	CONNECT(elCand1_ID);
	CONNECT(elCand2_ID);
	CONNECT(elCand1_pt);
	CONNECT(elCand2_pt);
	CONNECT(elCand1_cl_eta);
	CONNECT(elCand2_cl_eta);
	#undef CONNECT

	long long nEvents = oldChain->GetEntries();
	cout << nEvents << " in oldChain" << endl; 
	for(long long i=0; i< nEvents; i++){
		oldChain->GetEntry(i);
		if (int(fabs(evt.elCand1_ID/1000))!=11 || int(fabs(evt.elCand2_ID/1000))!=11) continue;

		hOldEta->Fill(evt.elCand1_cl_eta);
		hOldEta->Fill(evt.elCand2_cl_eta);
		hOldLPt->Fill(max(evt.elCand1_pt, evt.elCand2_pt));
		hOldSLPt->Fill(min(evt.elCand1_pt, evt.elCand2_pt));
		loadbar(i,nEvents);
	}
	cout << endl;
	oldFilelist.close();
	oldChain->SetBranchStatus("*",0);

	TChain* newChain = new TChain("ZeeCandidate");
	ifstream newFilelist("../common/inFileList-MCconverted.txt");
	while (!newFilelist.eof()){
		newFilelist >> filename;
		newChain->Add(filename.c_str());
		cout << "Added " << filename << " to `newChain`" << endl;
	}

	event evtNew;

	newChain->SetBranchStatus("*",0);
	#define CONNECT(b) evtNew.b = 0; newChain->SetBranchStatus(#b,1); newChain->SetBranchAddress(#b,&(evtNew.b)); 
	CONNECT(elCand1_ID);
	CONNECT(elCand2_ID);
	CONNECT(elCand1_pt);
	CONNECT(elCand2_pt);
	CONNECT(elCand1_cl_eta);
	CONNECT(elCand2_cl_eta);
	#undef CONNECT

	nEvents = newChain->GetEntries();
	cout << nEvents << " in newChain" << endl;
	for(long long i=0; i< nEvents; i++){
		newChain->GetEntry(i);
		if (int(fabs(evtNew.elCand1_ID/1000))!=11 || int(fabs(evtNew.elCand2_ID/1000))!=11) continue;

		hNewEta->Fill(evtNew.elCand1_cl_eta);
		hNewEta->Fill(evtNew.elCand2_cl_eta);
		hNewLPt->Fill(max(evtNew.elCand1_pt, evtNew.elCand2_pt));
		hNewSLPt->Fill(min(evtNew.elCand1_pt, evtNew.elCand2_pt));
		loadbar(i,nEvents);
	}
	cout << endl;
	newFilelist.close();
	newChain->SetBranchStatus("*",0);

	f->Write();
}
