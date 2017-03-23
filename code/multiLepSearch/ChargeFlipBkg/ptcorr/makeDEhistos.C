/* makeDEhistos.C
 * Gabriel Gallardo
 * 5 July 2015
 */

#include <fstream>
#include "../common/common.C"

// ======= CONFIGURATION ======== //
bool signalOnly = true;
bool applyPRW = true;
bool passQID = true;

string defaultOutDir="../QiD-on/ptcorr";
string defaultInFileTxt="/afs/cern.ch/user/g/ggallard/private/SUSY2L/code/multiLepSearch/ChargeFlipBkg/common/inFileList-ZeeConverted.txt";

// ========= MAIN FUNCTION ======== //
int makeDEhistos(string outDir=defaultOutDir, string inFileTxt=defaultInFileTxt){

	if(signalOnly){
		cout << "Selecting Signal events" << endl;
	} else {
		cout << "Selecting LooseBaseline events" << endl;
	}

	ifstream fileList(inFileTxt.c_str());

	if (!gSystem->cd(outDir.c_str())){
		gSystem->mkdir(outDir.c_str());
		if(!gSystem->cd(outDir.c_str())){
			cout << "Output directory " << outDir << " could not be created." << endl;
			return -1;
		}
	}

	TFile* outFile = TFile::Open("dEhistos.root", "recreate");

	double etas[] = {0, 0.5, 1, 1.37, 1.52, 1.8, 2.0, 2.47}; // 7 bins
	double pts[] = {20, 30, 40, 50, 60, 80, 120, 1000.0}; // 7 bins
	double dEbins[101];

	for(int i=0; i<101; i++){
		dEbins[i]=i-50;
	}

	TH3D* hDE = new TH3D("hDE", "#Delta E for all electrons;Reco |#eta|;Reco p_{T} (GeV);#Delta E (GeV)", 7, etas, 7, pts, 100, dEbins);

	hDE->SetDirectory(outFile);

	TH3D* hDEflipped = (TH3D*) hDE->Clone("hDEflipped");
	hDEflipped->SetTitle("#Delta E for charge flipped electrons");
	hDEflipped->SetDirectory(outFile);

	TH3D* hDEok = (TH3D*) hDE->Clone("hDEok");
	hDEok->SetTitle("#Delta E for charge ok electrons");
	hDEok->SetDirectory(outFile);

	TH3D* hDPT = (TH3D*) hDE->Clone("hDPT");
	hDPT->SetTitle("#Delta p_{T} for all electrons");
	hDPT->GetZaxis()->SetTitle("#Delta p_{T} (GeV)");
	hDPT->SetDirectory(outFile);

	TH3D* hDPTflipped = (TH3D*) hDPT->Clone("hDPTflipped");
	hDPTflipped->SetTitle("#Delta p_{T} for charge flipped electrons");
	hDPTflipped->SetDirectory(outFile);

	TH3D* hDPTok = (TH3D*) hDPT->Clone("hDPTok");
	hDPTok->SetTitle("#Delta p_{T} for charge ok electrons");
	hDPTok->SetDirectory(outFile);

	double dRbins[101];
	for(int i=0; i<101; i++){
		dRbins[i] = i*(0.04/100);
	}

	TH3D* hDR = new TH3D("hDR", "#Delta R for all electrons;Reco |#eta|;Reco p_{T} (GeV);#Delta R", 7, etas, 7, pts, 100, dRbins);
	hDR->SetDirectory(outFile);

	TH3D* hDRflipped = (TH3D*) hDR->Clone("hDRflipped");
	hDRflipped->SetTitle("#Delta R for charge flipped electrons");
	hDRflipped->SetDirectory(outFile);

	TH3D* hDRok = (TH3D*) hDR->Clone("hDRok");
	hDRok->SetTitle("#Delta R for charge ok electrons");
	hDRok->SetDirectory(outFile);

	TH1D* hTruePt = new TH1D("hTruePt", "True pt", 100, 0, 200);
	hTruePt->SetDirectory(outFile);

	TChain* ch = new TChain("ZeeCandidate");
	string fname;
	while(!fileList.eof()){
		fileList >> fname;
		ch->Add(fname.c_str());
	}

	COMMON::ev event;

	ch->SetBranchStatus("*",0);

	#define CONNECT(b) event.b = 0; ch->SetBranchStatus(#b,1); ch->SetBranchAddress(#b,&(event.b)); 
	CONNECT(MCEvtWeight)
	CONNECT(MCPileupWeight)

	CONNECT(elCand1_cl_eta)
	CONNECT(elCand1_pt)
	CONNECT(elCand1_charge)
	CONNECT(elCand1_E)
	
	CONNECT(elCand1_origE)
	CONNECT(elCand1_origPt)
	CONNECT(elCand1_origCharge)

	CONNECT(elCand2_cl_eta)
	CONNECT(elCand2_pt)
	CONNECT(elCand2_charge)
	CONNECT(elCand2_E)

	CONNECT(elCand2_origE)
	CONNECT(elCand2_origPt)
	CONNECT(elCand2_origCharge)

	CONNECT(elCand1_flag)
	CONNECT(elCand2_flag)

	CONNECT(elCand1_flag)
	CONNECT(elCand2_flag)

	CONNECT(elCand1_dRwOrig)
	CONNECT(elCand2_dRwOrig)

	CONNECT(elCand1_qID);
	CONNECT(elCand2_qID);
	#undef CONNECT


	long long nEvents = ch->GetEntries();
	float deltaE, deltaPt;
	cout << nEvents << " events in chain." << endl;
	for(int i=0; i<nEvents; i++){
		ch->GetEntry(i);
		COMMON::loadbar(i+1,nEvents);

		if(signalOnly){
			if(!((event.elCand1_flag & 2)/2)) continue;
			if(!((event.elCand2_flag & 2)/2)) continue;
		}
		if(passQID){
			if(!event.elCand2_qID || !event.elCand1_qID) continue;
		}

		double w = event.MCEvtWeight;
		if(applyPRW) w*= event.MCPileupWeight;

		if(event.elCand1_origCharge!=0){
			deltaE = event.elCand1_E - event.elCand1_origE;
			deltaPt = event.elCand1_pt - event.elCand1_origPt;

			hDE->Fill(event.elCand1_cl_eta, event.elCand1_pt, deltaE, w);
			hDPT->Fill(event.elCand1_cl_eta, event.elCand1_pt, deltaPt, w);
			hDR->Fill(event.elCand1_cl_eta, event.elCand1_pt, event.elCand1_dRwOrig, w);

			hTruePt->Fill(event.elCand1_origPt);

			if(event.elCand1_origCharge != event.elCand1_charge){
				hDEflipped->Fill(event.elCand1_cl_eta, event.elCand1_pt, deltaE, w);
				hDPTflipped->Fill(event.elCand1_cl_eta, event.elCand1_pt, deltaPt, w);
				hDRflipped->Fill(event.elCand1_cl_eta, event.elCand1_pt, event.elCand1_dRwOrig, w);
			} else {
				hDEok->Fill(event.elCand1_cl_eta, event.elCand1_pt, deltaE, w);
				hDPTok->Fill(event.elCand1_cl_eta, event.elCand1_pt, deltaPt, w);
				hDRok->Fill(event.elCand1_cl_eta, event.elCand1_pt, event.elCand1_dRwOrig, w);
			}
		}

		if(event.elCand2_origCharge!=0){
			deltaE = event.elCand2_E - event.elCand2_origE;
			deltaPt = event.elCand2_pt - event.elCand2_origPt;

			hDE->Fill(event.elCand2_cl_eta, event.elCand2_pt, deltaE, w);
			hDPT->Fill(event.elCand2_cl_eta, event.elCand2_pt, deltaPt, w);
			hDR->Fill(event.elCand2_cl_eta, event.elCand2_pt, event.elCand2_dRwOrig, w);

			hTruePt->Fill(event.elCand2_origPt);

			if(event.elCand2_origCharge != event.elCand2_charge){
				hDEflipped->Fill(event.elCand2_cl_eta, event.elCand2_pt, deltaE, w);
				hDPTflipped->Fill(event.elCand2_cl_eta, event.elCand2_pt, deltaPt, w);
				hDRflipped->Fill(event.elCand2_cl_eta, event.elCand2_pt, event.elCand2_dRwOrig, w);
			} else {
				hDEok->Fill(event.elCand2_cl_eta, event.elCand2_pt, deltaE, w);
				hDPTok->Fill(event.elCand2_cl_eta, event.elCand2_pt, deltaPt, w);
				hDRok->Fill(event.elCand2_cl_eta, event.elCand2_pt, event.elCand2_dRwOrig, w);
			}
		}
	}

	outFile->Write();
	std::cout << "Done" << "\n";
	return 0;
}
