/*

getEffs.C
A ROOT macro to derive efficiencies for the matrix method of
non-prompt lepton estimation for the 2LSS C1N2->Wh analysis.

Written by [Gabriel Gallardo](mailto:gabriel.gallardo@cern.ch)
June-July 2017


Execution:
root -l -b -q getEffs.C+
*/

#include <vector>
#include "support/evt2l.C"
#include "support/obj_def.h"
#include "support/MCTruthClassifierDefs.h"
#include <TChain.h>
#include <TTree.h>
#include <TH2D.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <TProfile.h>
#include <TColor.h>
#include <TMath.h>
#include <TRegexp.h>
#include <TSystem.h>
#include <TLorentzVector.h>
#include "SUSYTools/SUSYCrossSection.h"

using std::cout;
using std::endl;

// ======= CONFIG VARIABLES ======== //
const double ptBins[] = {20, 40, 80, 150};
const int nPtBins = sizeof(ptBins)/sizeof(double)-1;
const double etaBins[] = {0, 1.37, 2.5};
const int nEtaBins = sizeof(etaBins)/sizeof(double)-1;
const double LUMI=36074.56;

// const TString inFileData="support/inFileList-dataTest.txt";
const TString inFileData="support/inFileList-data.txt";

// TString inFileBBCC="support/inFileList-.txt";
const TString inFileTTbar="support/inFileList-ttbar.txt";
const TString inFileZjets="support/inFileList-Zjets.txt";
const TString inFileWjets="support/inFileList-Wjets.txt";
const TString inFileVV="support/inFileList-VV.txt";
const TString inFileVgamma="support/inFileList-Vgamma.txt";
const TString inFileSingleTop="support/inFileList-SingleTop.txt";
const TString inFileMC="support/inFileList-allMC.txt";

const bool DEBUG=true;
const bool VERBOSE=DEBUG || true;

const TString outDir="output/";


bool passRatesCR(evt2l* tree)
{
	// Exactly two leptons
	bool pass = true;
	pass *= tree->leps_==2;
	pass *= tree->jets_>=1;

	return pass;
}

bool passLepCompCR(evt2l* tree)
{
	bool pass = true;
	pass *= tree->leps_==2;
	pass *= tree->leps_ID[0]*tree->leps_ID[1]>0;
	pass *= tree->jets_>=1;

	// // No b-jets
	// for(int i(0); i>tree->jets_;i++) pass*= !(tree->jets_jFlag[i] & JT_BJET);

	return pass;
}

// ======= INFRASTRUCTURE ===== //

SUSY::CrossSectionDB *xsecDB;
TFile* effsFile;
TFile* nEventsFile;
TH2D* hPrototype= new TH2D("hPrototype", ";p_{T} [GeV];|#eta|", nPtBins, ptBins, nEtaBins, etaBins);

enum LEP_SOURCE {
	REAL=-1,
	HEAVY=0,
	LIGHT=1,
	CONV=2,
	OTHER=3, // To check lepton composition
};
const int N_FAKES_SOURCE=4;

const TString LS_TOSTRING(LEP_SOURCE s)
{
	if(s==REAL) return "Real";
	if(s==HEAVY) return "Heavy";
	if(s==LIGHT) return "Light";
	if(s==CONV) return "Conv";
	if(s==OTHER) return "Other";
	return "Unknown";
}

enum LEP_PROC
{
	TTBAR=0,
	// BBCC,
	ZJETS,
	WJETS,
	VV,
	VGAMMA,
	SINGLETOP,
	// ALL_MC,
};
const int N_PROC=6;

const TString LP_TOSTRING(LEP_PROC p)
{
	if(p==TTBAR) return "TTBar";
	// if(p==BBCC) return "BBCC";
	if(p==ZJETS) return "Zjets";
	if(p==WJETS) return "Wjets";
	if(p==VV) return "VV";
	if(p==VGAMMA) return "Vgamma";
	if(p==SINGLETOP) return "SingleTop";
	return "Unknown";
}

enum DATA_TYPE {
	DATA,
	MC,
};

struct lepInfo
{
	double pt;
	double eta;
	double w;
	bool isTag;
	bool isTight;
};
class Histos;
std::vector<Histos*> allHistos;

class Histos
{
private:
	TH2D* NTight;
	TH2D* NLoose;
	TH2D* Eff;
	bool ignoreTagTag=false;

public:
	Histos(): Histos("Test") {}
	Histos(TString name)
	{
		NTight = (TH2D*) hPrototype->Clone(name+"_NTight");
		NLoose = (TH2D*) hPrototype->Clone(name+"_NLoose");
		Eff = (TH2D*) hPrototype->Clone(name+"_Eff"); Eff->SetTitle(name+"_Eff");

		NTight->SetDirectory(nEventsFile);
		NLoose->SetDirectory(nEventsFile);
		Eff->SetDirectory(effsFile);

		allHistos.push_back(this);
	}
	~Histos(){}

	bool Fill(lepInfo &lep)
	{
		NLoose->Fill(lep.pt, fabs(lep.eta), lep.w);
		if(lep.isTight) NTight->Fill(lep.pt, fabs(lep.eta), lep.w);
		return true;
	}

	bool Fill(lepInfo &lep1, lepInfo &lep2)
	{
		if (ignoreTagTag && (lep1.isTag && lep2.isTag)) return false;
		if (lep1.isTag)
		{
			NLoose->Fill(lep2.pt, fabs(lep2.eta), lep2.w);
			if(lep2.isTight) NTight->Fill(lep2.pt, fabs(lep2.eta), lep2.w);
			// return true;
		}
		if (lep2.isTag)
		{
			NLoose->Fill(lep1.pt, fabs(lep1.eta), lep1.w);
			if(lep1.isTight) NTight->Fill(lep1.pt, fabs(lep1.eta), lep1.w);
			// return true;
		}
		return true;
	}

	bool calcEff() {Eff->Add(NTight); Eff->Divide(NLoose); return true;}

	void Draw()
	{
		TCanvas can;

		Eff->Draw("text colz");
		can.Print(outDir+TString(Eff->GetName())+".pdf");

		// Eff->ProfileX()->Draw();
		// can.Print(outDir+TString(Eff->GetName())+"Pt.pdf");

		// Eff->ProfileY()->Draw();
		// can.Print(outDir+TString(Eff->GetName())+"Eta.pdf");

		NTight->Draw("text");
		can.Print(outDir+TString(NTight->GetName())+".pdf");

		NLoose->Draw("text");
		can.Print(outDir+TString(NLoose->GetName())+".pdf");

	}

	void vetoTagTag(bool b=true){ignoreTagTag=b;}
};

class LepProp;
class LepProp
{
private:
	bool DEBUG=false;
	bool normalized=false;
	typedef std::vector<double> binsOfEta;
	typedef std::vector<binsOfEta> binsOfPt;
	typedef std::vector<binsOfPt> binsOfProc;
	typedef std::vector<binsOfProc> binsOfSource;

	binsOfSource fakes;
	binsOfSource reals;

	TString name;

	bool AddSourceToVec(binsOfSource &v)
	{
		binsOfProc procRow;
		procRow.reserve(N_PROC);
		for(int nProc=0; nProc<N_PROC; nProc++) // Add all processes
		{
			binsOfPt ptRow;
			ptRow.reserve(nPtBins);
			for(int ptBin=0; ptBin < nPtBins; ptBin++) // Add pt bins
			{
				binsOfEta etaRow(nEtaBins, 0.0); // Reserve bins of eta
				ptRow.push_back(etaRow);
			}
			procRow.push_back(ptRow);
		}
		v.push_back(procRow);
		return true;
	}

	bool NormalizeVec(binsOfSource &v)
	{
		for(int ptBin=0; ptBin<nPtBins; ptBin++)
		for(int etaBin=0; etaBin<nEtaBins; etaBin++)
		{
			double sumW = 0;
			for(uint s=0; s<v.size(); s++)
			for(int p=0; p<N_PROC; p++)
				sumW +=v[s][p][ptBin][etaBin];

			for(uint s=0; s<v.size(); s++)
			for(int p=0; p<N_PROC; p++)
				v[s][p][ptBin][etaBin] /= sumW;
		}
		
		return true;
	}

	int GetPtBin(double pt)
	{
		if (pt<ptBins[0]) return -1; // Underflow
		for(int i=0; i<nPtBins; i++) if (pt<ptBins[i+1]) return i;
		return nPtBins-1; // Overflow, place in largest ptBin.
	}

	int GetEtaBin(double eta)
	{
		eta = TMath::Abs(eta);
		if (eta<etaBins[0]) return -1; // Underflow, not really possible for 
		for(int i=0; i<nEtaBins; i++) if (eta<etaBins[i+1]) return i;
		return nEtaBins-1; // Overflow, place in largest etaBin.
	}

public:
	LepProp():LepProp("test"){};
	LepProp(TString name){
		this->name = name;
		AddSourceToVec(reals);
		for(int i=0; i<N_FAKES_SOURCE; i++) AddSourceToVec(fakes);

		if (DEBUG){
			cout << "LepProp " << name << "initialized." << endl;
			cout << "The reals vector has size " << reals.size() << endl;
			cout << "The fakes vector has size " << fakes.size() << endl;
		}
	}
	~LepProp(){}

	bool Fill(LEP_SOURCE s, LEP_PROC p, double pt, double eta, double w)
	{
		if (DEBUG) cout << "### Enter " << name << "::Fill()" << endl;
		int ptBin = GetPtBin(pt); int etaBin = GetEtaBin(eta);
		if (DEBUG) cout << "### Pt/Eta Bin acquired" << endl;
		if(ptBin<0 || etaBin<0) return false;

		if(s==REAL) {
			if (DEBUG) cout << "### Filling real for " << LP_TOSTRING(p) << " " << LS_TOSTRING(s) << endl;
			reals[0][p][ptBin][etaBin] +=w;
		}
		else {
			if (DEBUG) cout << "### Filling fakes for " << LP_TOSTRING(p) << " " << LS_TOSTRING(s) << endl;
			fakes[s][p][ptBin][etaBin]+=w;
		}
		if (DEBUG) cout << "### LepProp filled. Return now." << endl;

		return true;
	}

	bool Normalize()
	{
		if(normalized) return true;
		NormalizeVec(reals);
		NormalizeVec(fakes);
		return true;
	};

	bool Write(){

		Normalize(); // Normalize histograms first

		// PROTOTYPE HISTOGRAMS
		TFile *lepCompFile = new TFile(outDir+"lepComp.root", "recreate");

		TH2D* hFakesPrototype2D = new TH2D("hFakesPrototype2D", "Fake lepton composition fractions;Process;Source", N_PROC, 0, N_PROC, N_FAKES_SOURCE, 0, N_FAKES_SOURCE);
		TH1D* hRealsPrototype = new TH1D("hRealsPrototype", "Real lepton composition fractions;Process;Fraction", N_PROC, 0, N_PROC);

		hFakesPrototype2D->SetDirectory(0);
		hRealsPrototype->SetDirectory(0);

		// Label bins
		for(int p=0; p<N_PROC; p++)
		{
			hRealsPrototype->GetXaxis()->SetBinLabel(p+1, LP_TOSTRING((LEP_PROC) p));
			hFakesPrototype2D->GetXaxis()->SetBinLabel(p+1, LP_TOSTRING((LEP_PROC) p));
		}
		for(int s=0;s<N_FAKES_SOURCE; s++)
		{
			hFakesPrototype2D->GetYaxis()->SetBinLabel(s+1, LS_TOSTRING((LEP_SOURCE) s));
		}


		// INTIALIZE HISTOGRAMS
		struct lpHistos
		{
			TH2* hFake;
			TH1* hReal;
			TH1* hFakeProc;
			TH1* hFakeSource;
		};

		std::vector< std::vector<lpHistos> > histVector;
		histVector.reserve(nPtBins);

		for(int ptBin=0; ptBin<nPtBins; ptBin++)
		{
			std::vector<lpHistos> v;
			v.reserve(nEtaBins);
			TString ptStr = "PT_";
			ptStr+=int(ptBins[ptBin]);
			ptStr+="_";
			ptStr+=int(ptBins[ptBin+1]);

			TString ptTitle = " p_{T} #in["; ptTitle+=ptBins[ptBin]; ptTitle+=",";
			ptTitle += ptBins[ptBin+1]; ptTitle += ")";

			for(int etaBin=0; etaBin<nEtaBins; etaBin++)
			{
				lpHistos tmp;

				TString ptEtaStr = ptStr + "_ETA_";
				ptEtaStr+=int(etaBins[etaBin]*100);
				ptEtaStr+="_";
				ptEtaStr+=int(etaBins[etaBin+1]*100);

				char tmpStr[4];
				sprintf(tmpStr, "%4.2f", etaBins[etaBin]);
				TString ptEtaTitle = ptTitle + ", |#eta| #in["; ptEtaTitle+=tmpStr; ptEtaTitle+=",";
				sprintf(tmpStr, "%4.2f", etaBins[etaBin+1]);
				ptEtaTitle += tmpStr; ptEtaTitle += ")";

				TString fakesName = (TString) "h" + name + (TString) + "_FakeComp_" + ptEtaStr;
				tmp.hFake = (TH2*) hFakesPrototype2D->Clone(fakesName);
				tmp.hFake->SetTitle(ptEtaTitle.Prepend(tmp.hFake->GetTitle()));
				tmp.hFake->SetDirectory(lepCompFile);
				tmp.hFake->GetZaxis()->SetRangeUser(0,1);

				TString realsName = (TString) "h" + name + (TString) + "_RealComp_" + ptEtaStr;
				tmp.hReal = (TH1*) hRealsPrototype->Clone(realsName);
				tmp.hReal->SetTitle(ptEtaTitle.Prepend(tmp.hReal->GetTitle()));
				tmp.hReal->SetDirectory(lepCompFile);
				tmp.hReal->GetYaxis()->SetRangeUser(0,1);

				v.push_back(tmp);
			}
			histVector.push_back(v);
		}


		// FILL HISTOGRAMS
		for(int p=0; p<N_PROC; p++)
		for(int ptBin=0; ptBin<nPtBins; ptBin++)
		for(int etaBin=0; etaBin<nEtaBins; etaBin++)
		{
			histVector[ptBin][etaBin].hReal->Fill(p, reals[0][p][ptBin][etaBin]);
			for(int s=0;s<N_FAKES_SOURCE; s++)
			{
				histVector[ptBin][etaBin].hFake->Fill(p, s, fakes[s][p][ptBin][etaBin]);
			}
		}

		// OUTPUT HISTOGRAMS
		TCanvas *c = new TCanvas();
		for(int ptBin=0; ptBin<nPtBins; ptBin++)
		for(int etaBin=0; etaBin<nEtaBins; etaBin++)
		{
			// Real composition 
			histVector[ptBin][etaBin].hReal->SetFillColor(kAzure+10);
			histVector[ptBin][etaBin].hReal->Draw("hist bar text");
			c->Print(outDir+((TString)".pdf").Prepend(histVector[ptBin][etaBin].hReal->GetName()));

			// Fake composition 2D
			histVector[ptBin][etaBin].hFake->Draw("colz text");
			c->Print(outDir+((TString)".pdf").Prepend(histVector[ptBin][etaBin].hFake->GetName()));

			// Fake composition by process
			histVector[ptBin][etaBin].hFakeProc = histVector[ptBin][etaBin].hFake->ProjectionX();
			histVector[ptBin][etaBin].hFakeProc->SetDirectory(lepCompFile);
			histVector[ptBin][etaBin].hFakeProc->SetName(((TString)"_byProcess").Prepend(histVector[ptBin][etaBin].hFake->GetName()));
			histVector[ptBin][etaBin].hFakeProc->SetTitle(((TString)" by process").Prepend(histVector[ptBin][etaBin].hFake->GetTitle()));
			histVector[ptBin][etaBin].hFakeProc->GetXaxis()->SetRangeUser(0, N_PROC);
			histVector[ptBin][etaBin].hFakeProc->GetYaxis()->SetRangeUser(0, 1);
			histVector[ptBin][etaBin].hFakeProc->SetFillColor(kAzure+10);
			histVector[ptBin][etaBin].hFakeProc->Draw("hist bar text");
			c->Print(outDir+((TString)".pdf").Prepend(histVector[ptBin][etaBin].hFakeProc->GetName()));

			// Fake composition by source
			histVector[ptBin][etaBin].hFakeSource = histVector[ptBin][etaBin].hFake->ProjectionY();
			histVector[ptBin][etaBin].hFakeSource->SetDirectory(lepCompFile);
			histVector[ptBin][etaBin].hFakeSource->SetName(((TString)"_bySource").Prepend(histVector[ptBin][etaBin].hFake->GetName()));
			histVector[ptBin][etaBin].hFakeSource->SetTitle(((TString)" by source").Prepend(histVector[ptBin][etaBin].hFake->GetTitle()));
			histVector[ptBin][etaBin].hFakeSource->GetXaxis()->SetRangeUser(0, N_FAKES_SOURCE);
			histVector[ptBin][etaBin].hFakeSource->GetYaxis()->SetRangeUser(0, 1);
			histVector[ptBin][etaBin].hFakeSource->SetFillColor(kAzure+10);
			histVector[ptBin][etaBin].hFakeSource->Draw("hist bar text");
			c->Print(outDir+((TString)".pdf").Prepend(histVector[ptBin][etaBin].hFakeSource->GetName()));

		}
		lepCompFile->Write();
		lepCompFile->Close();


		return true;
	} // Should write this

	double GetWeight(LEP_SOURCE s, LEP_PROC p, double pt, double eta)
	{
		int ptBin = GetPtBin(pt); int etaBin = GetEtaBin(eta);
		if(ptBin<0 || etaBin<0) return 0;
		if(s==REAL) return reals[0][p][ptBin][etaBin];
		else return fakes[s][p][ptBin][etaBin];
	}
};
LepProp* lp_Mu;
LepProp* lp_El;

TH2* hMu_MCTC;
TH2* hEl_MCTC;

Histos* hMu_Real_MCTP;
Histos* hEl_Real_MCTP;
Histos* hMu_Heavy_MCTP;
Histos* hEl_Heavy_MCTP;
Histos* hMu_Light_MCTP;
Histos* hEl_Light_MCTP;
Histos* hEl_Conv_MCTP;

// ======= FUNCTION PROTOTYPES ===== //

int getEffs();
bool initialize();
TChain* loadData(TString, bool isMC=false);
bool finalize();
bool doTP(evt2l* tree, Histos* hMu_Real, Histos* hEl_Real, Histos* hMu_Heavy, Histos* hEl_Heavy, Histos* hEl_Conv);
bool loopMC(evt2l* tree, Histos *hMu_Real, Histos *hEl_Real, Histos *hMu_Heavy, Histos *hEl_Heavy, Histos *hMu_Light, Histos *hEl_Light, Histos *hEl_Conv, LEP_PROC p);
LEP_SOURCE castSource(MCTC::ParticleOrigin);
void labelMCTChist(TH2* h);

// ======== FUNCTION DEFINITIONS ======= //
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

void labelMCTChist(TH2* h)
{
	MCTC::ParticleDef pd;

	for(int i=1; i<=PARTICLEORIGIN; i++)
		h->GetYaxis()->SetBinLabel(i,pd.sParticleOrigin[i-1].c_str());

	for(int i=1; i<=PARTICLETYPES; i++)
		h->GetXaxis()->SetBinLabel(i,pd.sParticleType[i-1].c_str());

	return;
}

bool initialize()
{
	if(DEBUG) cout << "Begin initialize()" <<endl;
	gStyle->SetOptStat(0); 
	gStyle->SetPalette(kLightTerrain);

	effsFile = new TFile(outDir+"efficiencies.root", "recreate");
	nEventsFile = new TFile(outDir+"nEvents.root", "recreate");

	gROOT->Macro("$ROOTCOREDIR/scripts/load_packages.C");
	xsecDB = new SUSY::CrossSectionDB(gSystem->ExpandPathName("$ROOTCOREBIN/data/SUSYTools/mc15_13TeV/"));
	if (!xsecDB){
	   cout << "CrossSectionDB could not be loaded" << endl;
	   return false;
	}

	hMu_MCTC = new TH2D("hMu_MCTC", "MCTC Mu;Type;Origin", PARTICLETYPES, 0, PARTICLETYPES, PARTICLEORIGIN, 0, PARTICLEORIGIN);
	hEl_MCTC = new TH2D("hEl_MCTC", "MCTC El;Type;Origin", PARTICLETYPES, 0, PARTICLETYPES, PARTICLEORIGIN, 0, PARTICLEORIGIN);

	hMu_MCTC->SetDirectory(effsFile);
	hEl_MCTC->SetDirectory(effsFile);

	lp_El = new LepProp("El");
	lp_Mu = new LepProp("Mu");

	hMu_Real_MCTP = new Histos("hMu_Real_MCTP");
	hEl_Real_MCTP = new Histos("hEl_Real_MCTP");
	hMu_Heavy_MCTP = new Histos("hMu_Heavy_MCTP");
	hEl_Heavy_MCTP = new Histos("hEl_Heavy_MCTP");
	hMu_Light_MCTP = new Histos("hMu_Light_MCTP");
	hEl_Light_MCTP = new Histos("hEl_Light_MCTP");
	hEl_Conv_MCTP = new Histos("hEl_Conv_MCTP");


	if(VERBOSE) cout << "initialize() complete" << endl;
	return true;
}

int getEffs()
{
	initialize();
	TChain* tcMC = new TChain("evt2l");

	// lp_El->Fill(REAL, TTBAR, 25, 1.4, 1);
	// lp_Mu->Fill(REAL, TTBAR, 25, 1.4, 1);

	// TTBAR
	evt2l* ttEvts = new evt2l(loadData(inFileTTbar, true));  tcMC->Add((TChain*) ttEvts->fChain);
	Histos* hMu_Real_TT = new Histos("hMu_Real_TT");
	Histos* hEl_Real_TT = new Histos("hEl_Real_TT");
	Histos* hMu_Heavy_TT = new Histos("hMu_Heavy_TT");
	Histos* hEl_Heavy_TT = new Histos("hEl_Heavy_TT");
	Histos* hMu_Light_TT = new Histos("hMu_Light_TT");
	Histos* hEl_Light_TT = new Histos("hEl_Light_TT");
	Histos* hEl_Conv_TT = new Histos("hEl_Conv_TT");
	cout << "Begin loop over ttBar: " << ttEvts->fChain->GetEntries() << " events" << endl;
	loopMC(ttEvts, hMu_Real_TT, hEl_Real_TT, hMu_Heavy_TT, hEl_Heavy_TT, hMu_Light_TT, hEl_Light_TT, hEl_Conv_TT, TTBAR);
	cout << "ttBar loop finished" << endl << endl;

	// BBBAR/CCBAR
	// evt2l* bbccEvts = new evt2l(loadData(inFileBBCC, true)); // tcMC->Add((TChain*) bbccEvts->fChain);
	// Histos* = new Histos hMu_Real_BBCC("hMu_Real_BBCC");
	// Histos* = new Histos hEl_Real_BBCC("hEl_Real_BBCC");
	// Histos* = new Histos hMu_Heavy_BBCC("hMu_Heavy_BBCC");
	// Histos* = new Histos hEl_Heavy_BBCC("hEl_Heavy_BBCC");
	// Histos* = new Histos hMu_Light_BBCC("hMu_Light_BBCC");
	// Histos* = new Histos hEl_Light_BBCC("hEl_Light_BBCC");
	// Histos* = new Histos hEl_Conv_BBCC("hEl_Conv_BBCC");
	// cout << "Begin loop over bb/cc: " << bbccEvts->fChain->GetEntries() << " events" << endl;
	// loopMC(bbccEvts, hMu_Real_BBCC, hEl_Real_BBCC, hMu_Heavy_BBCC, hEl_Heavy_BBCC, hMu_Light_BBCC, hEl_Light_BBCC, hEl_Conv_BBCC);
	// cout << "bb/cc loop finished" << endl << endl;

	// // Z+jets
	evt2l* zJetsEvts = new evt2l(loadData(inFileZjets, true));  tcMC->Add((TChain*) zJetsEvts->fChain);
	Histos* hMu_Real_Zjets = new Histos("hMu_Real_Zjets");
	Histos* hEl_Real_Zjets = new Histos("hEl_Real_Zjets");
	Histos* hMu_Heavy_Zjets = new Histos("hMu_Heavy_Zjets");
	Histos* hEl_Heavy_Zjets = new Histos("hEl_Heavy_Zjets");
	Histos* hMu_Light_Zjets = new Histos("hMu_Light_Zjets");
	Histos* hEl_Light_Zjets = new Histos("hEl_Light_Zjets");
	Histos* hEl_Conv_Zjets = new Histos("hEl_Conv_Zjets");
	cout << "Begin loop over Z+jets: " << zJetsEvts->fChain->GetEntries() << " events" << endl;
	loopMC(zJetsEvts, hMu_Real_Zjets, hEl_Real_Zjets, hMu_Heavy_Zjets, hEl_Heavy_Zjets, hMu_Light_Zjets, hEl_Light_Zjets, hEl_Conv_Zjets, ZJETS);
	cout << "Z+jets loop finished" << endl << endl;

	// // W+jets
	evt2l* wJetsEvts = new evt2l(loadData(inFileWjets, true));  tcMC->Add((TChain*) wJetsEvts->fChain);
	Histos* hMu_Real_Wjets = new Histos("hMu_Real_Wjets");
	Histos* hEl_Real_Wjets = new Histos("hEl_Real_Wjets");
	Histos* hMu_Heavy_Wjets = new Histos("hMu_Heavy_Wjets");
	Histos* hEl_Heavy_Wjets = new Histos("hEl_Heavy_Wjets");
	Histos* hMu_Light_Wjets = new Histos("hMu_Light_Wjets");
	Histos* hEl_Light_Wjets = new Histos("hEl_Light_Wjets");
	Histos* hEl_Conv_Wjets = new Histos("hEl_Conv_Wjets");
	cout << "Begin loop over W+jets: " << wJetsEvts->fChain->GetEntries() << " events" << endl;
	loopMC(wJetsEvts, hMu_Real_Wjets, hEl_Real_Wjets, hMu_Heavy_Wjets, hEl_Heavy_Wjets, hMu_Light_Wjets, hEl_Light_Wjets, hEl_Conv_Wjets, WJETS);
	cout << "W+jets loop finished" << endl << endl;

	// DIBOSON
	evt2l* vvEvts = new evt2l(loadData(inFileVV, true));  tcMC->Add((TChain*) vvEvts->fChain);
	Histos* hMu_Real_VV = new Histos ("hMu_Real_VV");
	Histos* hEl_Real_VV = new Histos ("hEl_Real_VV");
	Histos* hMu_Heavy_VV = new Histos ("hMu_Heavy_VV");
	Histos* hEl_Heavy_VV = new Histos ("hEl_Heavy_VV");
	Histos* hMu_Light_VV = new Histos ("hMu_Light_VV");
	Histos* hEl_Light_VV = new Histos ("hEl_Light_VV");
	Histos* hEl_Conv_VV = new Histos ("hEl_Conv_VV");
	cout << "Begin loop over VV: " << vvEvts->fChain->GetEntries() << " events" << endl;
	loopMC(vvEvts, hMu_Real_VV, hEl_Real_VV, hMu_Heavy_VV, hEl_Heavy_VV, hMu_Light_VV, hEl_Light_VV, hEl_Conv_VV, VV);
	cout << "VV loop finished" << endl << endl;

	// VGAMMA
	evt2l* vgammaEvts = new evt2l(loadData(inFileVgamma, true));  tcMC->Add((TChain*) vgammaEvts->fChain);
	Histos* hMu_Real_Vgamma = new Histos("hMu_Real_Vgamma");
	Histos* hEl_Real_Vgamma = new Histos("hEl_Real_Vgamma");
	Histos* hMu_Heavy_Vgamma = new Histos("hMu_Heavy_Vgamma");
	Histos* hEl_Heavy_Vgamma = new Histos("hEl_Heavy_Vgamma");
	Histos* hMu_Light_Vgamma = new Histos("hMu_Light_Vgamma");
	Histos* hEl_Light_Vgamma = new Histos("hEl_Light_Vgamma");
	Histos* hEl_Conv_Vgamma = new Histos("hEl_Conv_Vgamma");
	cout << "Begin loop over Vgamma: " << vvEvts->fChain->GetEntries() << " events" << endl;
	loopMC(vgammaEvts, hMu_Real_Vgamma, hEl_Real_Vgamma, hMu_Heavy_Vgamma, hEl_Heavy_Vgamma, hMu_Light_Vgamma, hEl_Light_Vgamma, hEl_Conv_Vgamma, VGAMMA);
	cout << "Vgamma loop finished" << endl << endl;

	// SingleTop
	evt2l* singleTopEvts = new evt2l(loadData(inFileSingleTop, true));  tcMC->Add((TChain*) vgammaEvts->fChain);
	Histos* hMu_Real_SingleTop = new Histos("hMu_Real_SingleTop");
	Histos* hEl_Real_SingleTop = new Histos("hEl_Real_SingleTop");
	Histos* hMu_Heavy_SingleTop = new Histos("hMu_Heavy_SingleTop");
	Histos* hEl_Heavy_SingleTop = new Histos("hEl_Heavy_SingleTop");
	Histos* hMu_Light_SingleTop = new Histos("hMu_Light_SingleTop");
	Histos* hEl_Light_SingleTop = new Histos("hEl_Light_SingleTop");
	Histos* hEl_Conv_SingleTop = new Histos("hEl_Conv_SingleTop");
	cout << "Begin loop over single top: " << vvEvts->fChain->GetEntries() << " events" << endl;
	loopMC(singleTopEvts, hMu_Real_SingleTop, hEl_Real_SingleTop, hMu_Heavy_SingleTop, hEl_Heavy_SingleTop, hMu_Light_SingleTop, hEl_Light_SingleTop, hEl_Conv_SingleTop, SINGLETOP);
	cout << "single top loop finished" << endl << endl;

	// ALL MC
	// evt2l* mcEvts = new evt2l(tcMC); 
	// Histos* hMu_Real_MC = new Histos("hMu_Real_MC");
	// Histos* hEl_Real_MC = new Histos("hEl_Real_MC");
	// Histos* hMu_Heavy_MC = new Histos("hMu_Heavy_MC");
	// Histos* hEl_Heavy_MC = new Histos("hEl_Heavy_MC");
	// Histos* hMu_Light_MC = new Histos("hMu_Light_MC");
	// Histos* hEl_Light_MC = new Histos("hEl_Light_MC");
	// Histos* hEl_Conv_MC = new Histos("hEl_Conv_MC");
	// cout << "Begin loop over all MC: " << mcEvts->fChain->GetEntries() << " events" << endl;
	// doTP(mcEvts, hMu_Real_MC, hEl_Real_MC, hMu_Heavy_MC, hEl_Heavy_MC, hEl_Conv_MC);
	// cout << "All MC loop finished" << endl;

	// DATA TAG-AND-PROBE
	evt2l* dataEvts = new evt2l(loadData(inFileData));
	cout << "Executing tag-and-probe on data" << endl;
	Histos* hMu_Real_Data = new Histos("hMu_Real_Data");
	Histos* hEl_Real_Data = new Histos("hEl_Real_Data");
	Histos* hMu_Heavy_Data = new Histos("hMu_Heavy_Data");
	Histos* hEl_Heavy_Data = new Histos("hEl_Heavy_Data");
	Histos* hEl_Conv_Data = new Histos("hEl_Conv_Data");
	int nDataEntries = dataEvts->fChain->GetEntries();
	cout << "Begin loop over data" << endl;
	for(int i=0; i<nDataEntries; i++)
	{
		loadbar(i+1, nDataEntries);
		dataEvts->GetEntry(i);
		doTP(dataEvts, hMu_Real_Data, hEl_Real_Data, hMu_Heavy_Data, hEl_Heavy_Data, hEl_Conv_Data);
	}
	cout << "Data tag-and-probe finished" << endl;

	finalize();

	return 0;
}

bool doTP(evt2l* tree, Histos* hMu_Real, Histos* hEl_Real, Histos* hMu_Heavy, Histos* hEl_Heavy, Histos* hEl_Conv)
{
	if (tree->sig_trigCode==0) return false;
	if (tree->leps_==2)
	{			

		// mm and ee
		bool mumu = TMath::Abs(int(tree->leps_ID[0]/1000)*int(tree->leps_ID[0]/1000))==169;
		bool elel = TMath::Abs(int(tree->leps_ID[0]/1000)*int(tree->leps_ID[0]/1000))==121;
		bool isOS = (tree->leps_ID[0]*tree->leps_ID[1])<0;

		// "Tight" defined as signal electron
		bool lep1isTight = (tree->leps_lFlag[0] & IS_SIGNAL);
		bool lep2isTight = (tree->leps_lFlag[1] & IS_SIGNAL);

		/** UNWEIGHTED REAL EFFICIENCIES MEASURED BY Z_ll TAG-AND-PROBE IN DATA **/
		if(fabs(tree->l12_m - 91)<10 && isOS) // Select Z mass pair
		{
			// Tag requirement: pT >25, passes signal ("tight") cut, passes trigger, and is trigger matched
			bool lep1isTag = (tree->leps_pt[0]>25 && lep1isTight && (tree->leps_ID[0] & 1) && (tree->leps_ID[0] & 2));
			bool lep2isTag = (tree->leps_pt[1]>25 && lep2isTight && (tree->leps_ID[1] & 1) && (tree->leps_ID[1] & 2));

			lepInfo l1info = {tree->leps_pt[0], tree->leps_eta[0], 1, lep1isTag, lep1isTight};
			lepInfo l2info = {tree->leps_pt[1], tree->leps_eta[1], 1, lep2isTag, lep2isTight};
			if (mumu) hMu_Real->Fill(l1info, l2info);
			else if (elel) hEl_Real->Fill(l1info, l2info);	
		}

		/* UNWEIGHTED HEAVY EFFICIENCIES MEASURED BY ll TAG-AND-PROBE IN DATA */
		/** in W/signal suppressed region **/
		else if (tree->sig_Met<40)
		{
			// To suppress W background and signal
			// if (tree->sig_Met>40) continue; TMath::Sqrt(2*leps[0].pt*sig.Met*(1-TMath::Cos(leps[0].MET_dPhi)))<40

			std::pair<double, double> mT_lep_MET;
			mT_lep_MET.first  = TMath::Sqrt(2*tree->leps_pt[0]*tree->sig_Met*(1-TMath::Cos(tree->leps_MET_dPhi[0])));
			mT_lep_MET.second = TMath::Sqrt(2*tree->leps_pt[1]*tree->sig_Met*(1-TMath::Cos(tree->leps_MET_dPhi[1])));
			if(mT_lep_MET.first>40 || mT_lep_MET.second>40) return false;

			bool lep1isTag = (tree->leps_pt[0]>20 && lep1isTight && TMath::Abs(int(tree->leps_ID[0]/1000))==13 && (tree->leps_ID[0] & 1) && (tree->leps_ID[0] & 2));
			bool lep2isTag = (tree->leps_pt[1]>20 && lep2isTight && TMath::Abs(int(tree->leps_ID[1]/1000))==13 && (tree->leps_ID[1] & 1) && (tree->leps_ID[1] & 2));

			lepInfo l1info = {tree->leps_pt[0], tree->leps_eta[0], 1, lep1isTag, lep1isTight};
			lepInfo l2info = {tree->leps_pt[1], tree->leps_eta[1], 1, lep2isTag, lep2isTight};

			// Select mumu with mll cut, or emu
			if (mumu && tree->l12_m>40) hMu_Heavy->Fill(l1info, l2info);
			else hEl_Heavy->Fill(l1info, l2info);
		}		
	}

	/** UNWEIGHTED PHOTON CONVERSTION EFFICIENCIES MEASURED BY mumu-e TAG-AND-PROBE IN DATA **/ 
	else if (tree->leps_==3)
	{
		// Select Z_mumu pair and electron 
		if ((int(tree->leps_ID[0]/1000)*int(tree->leps_ID[1]/1000) != -169)	|| TMath::Abs(int(tree->leps_ID[2])/1000)!= 11)
			return false;

		// Dimuon trigger match
		if (!(tree->sig_trigMask && tree->sig_trigCode)) return false;

		TLorentzVector p1, p2, p3;
		p1.SetPtEtaPhiM(tree->leps_pt[0], tree->leps_eta[0], tree->leps_phi[0], 0.105658);
		p2.SetPtEtaPhiM(tree->leps_pt[1], tree->leps_eta[1], tree->leps_phi[1], 0.105658);
		p3.SetPtEtaPhiM(tree->leps_pt[2], tree->leps_eta[2], tree->leps_phi[2], 0.000511);

		if (TMath::Abs((p1+p2+p3).M()-91)>10) return false;

		bool elIsTight = (tree->leps_lFlag[2] & IS_SIGNAL);

		lepInfo elInfo = {tree->leps_pt[2], tree->leps_eta[2], 1, false, elIsTight};
		hEl_Conv->Fill(elInfo);
	}
	return true;
}

bool loopMC(evt2l* tree, Histos *hMu_Real, Histos *hEl_Real, Histos *hMu_Heavy,  Histos *hEl_Heavy, Histos *hMu_Light,  Histos *hEl_Light, Histos *hEl_Conv, LEP_PROC p)
{
	bool DEBUG=false;
	using namespace MCTC;

	// if(DEBUG){
		// if(lp_Mu) cout << "lp_Mu initialized" << endl;
		// if(lp_El) cout << "lp_El initialized" << endl;
	// }
	if(tree && DEBUG) cout << "tree exists" << endl;
	if(!tree) return false;

	int nEntries = tree->fChain->GetEntries();
	if (DEBUG) cout << "Begin loop MC" << endl;
	for(int i=0; i<nEntries; i++)
	{
		loadbar(i+1, nEntries);
		tree->GetEntry(i); 
		double w = tree->evt_weight*tree->evt_pwt*tree->fChain->GetTree()->GetWeight()*tree->evt_ElSF*tree->evt_MuSF;

		// Measure unweighted efficiencies
		if(passRatesCR(tree)){ for(int j(0); j<2; j++)
		{
			if(tree->leps_pt[j]<20) continue;
			bool lepIsTight = tree->leps_lFlag[j] & IS_SIGNAL;
			bool isRecoEl = TMath::Abs(int(tree->leps_ID[j]/1000)) == 11;
			bool isRecoMu = TMath::Abs(int(tree->leps_ID[j]/1000)) == 13;

			ParticleType   type = static_cast<ParticleType>(tree->leps_truthType[j]);
			ParticleOrigin orig = static_cast<ParticleOrigin>(tree->leps_truthOrig[j]);
			LEP_SOURCE 	   source = castSource(orig);

			lepInfo l = {tree->leps_pt[j], tree->leps_eta[j], w, true, lepIsTight};

			if (isRecoEl && type==IsoElectron) hEl_Real->Fill(l); 
			else if (isRecoMu && type==IsoMuon) hMu_Real->Fill(l); 
			else if (isRecoEl && (type==UnknownElectron || type==NonIsoElectron || type==BkgElectron))
			{
				if (source==HEAVY) hEl_Heavy->Fill(l); 
				else if (source==CONV) hEl_Conv->Fill(l);
				else if (source==LIGHT) hEl_Light->Fill(l); 
			}
			else if (isRecoMu && (type==UnknownMuon || type==NonIsoMuon || type==BkgMuon))
			{
				if (source==HEAVY) hMu_Heavy->Fill(l);
				else if (source==LIGHT) hMu_Light->Fill(l);
			}
			if (DEBUG) cout << "# Filled rates histograms for lepton " << j << endl;
		}}

		// Measure composition fractions
		if (passLepCompCR(tree)){ for(int j(0); j<2; j++)
		{
			if(tree->leps_pt[j]<20) continue;

			bool lepIsTight = tree->leps_lFlag[j] & IS_SIGNAL;
			bool isRecoEl = TMath::Abs(int(tree->leps_ID[j]/1000)) == 11;
			bool isRecoMu = TMath::Abs(int(tree->leps_ID[j]/1000)) == 13;

			ParticleType   type = static_cast<ParticleType>(tree->leps_truthType[j]);
			ParticleOrigin orig = static_cast<ParticleOrigin>(tree->leps_truthOrig[j]);
			LEP_SOURCE 	   source = castSource(orig);

			if (isRecoEl)
			{
				// if (DEBUG && type!=IsoElectron) cout << LS_TOSTRING(source) << "\t" << int(orig) << " electron" << endl;
				if (type!=IsoElectron) lp_El->Fill(source, p, tree->leps_pt[j], tree->leps_eta[j], w);
				else lp_El->Fill(REAL, p, tree->leps_pt[j], tree->leps_eta[j], w);
				hEl_MCTC->Fill(type, orig);
			}
			else if (isRecoMu)
			{
				// if (DEBUG && type!=IsoMuon) cout << LS_TOSTRING(source) << "\t" << int(orig) << " muon" << endl;
				if (type!=IsoMuon) lp_Mu->Fill(source, p, tree->leps_pt[j], tree->leps_eta[j], w);
				else lp_Mu->Fill(REAL, p, tree->leps_pt[j], tree->leps_eta[j], w);
				hMu_MCTC->Fill(type, orig);
			}
			if (DEBUG) cout << "# Filled lepton composition histograms for lepton " << j << endl;
		}}

		// Tag and probe for data/MC scale factors
		doTP(tree, hMu_Real_MCTP, hEl_Real_MCTP, hMu_Heavy_MCTP, hEl_Heavy_MCTP, hEl_Conv_MCTP);

		if (DEBUG) cout << "# After lepton loop" << endl;
	}
	cout << endl;
	return true;
}


LEP_SOURCE castSource(MCTC::ParticleOrigin orig)
{
	using namespace MCTC;
	if (orig==CharmedMeson || orig==BottomMeson || orig==CCbarMeson || orig==BBbarMeson 
		|| orig==CharmedBaryon || orig==BottomBaryon)
		return HEAVY;
	else if (orig==LightMeson || orig==StrangeMeson || orig==LightBaryon || orig==StrangeBaryon 
		|| orig==PionDecay || orig==KaonDecay || orig==PiZero|| orig==QCD || orig==NonDefined)
		return LIGHT;
	else if (orig==PhotonConv || orig==BremPhot || orig==PromptPhot || orig==UndrPhot 
		|| orig==ISRPhot || orig==FSRPhot)
		return CONV;
	else if (orig==NonDefined) return LIGHT;

	return OTHER;
}


bool finalize()
{
	if (VERBOSE) cout << "Cleaning up..." << endl;

	if (VERBOSE) cout << "Histograms" << endl;
	for (auto h : allHistos){h->calcEff(); h->Draw();}

	if (VERBOSE) cout << "Writing lepton composition" << endl;
	lp_El->Normalize(); lp_El->Write();
	lp_Mu->Normalize(); lp_Mu->Write();

	if (VERBOSE) cout << "Draw MCTC " << endl;
	TCanvas *can = new TCanvas("cMCTC", "cMCTC", 1280,800); can->SetGrid();
	if (hMu_MCTC) {labelMCTChist(hMu_MCTC); hMu_MCTC->Draw("text"); can->Print(outDir+"hMu_MCTC.pdf");}
	if (hEl_MCTC) {labelMCTChist(hEl_MCTC); hEl_MCTC->Draw("text");can->Print(outDir+"hEl_MCTC.pdf");}
	delete can; can=0;

	if (VERBOSE) cout << "Write to file" << endl;
	effsFile->Write(); effsFile->Close();
	nEventsFile->Write(); nEventsFile->Close();
	for (auto h: allHistos) if (h!=0) {delete h; h=0;}

	if (xsecDB!=0){delete xsecDB; xsecDB=0;}
	if (lp_Mu!=0){delete lp_Mu; lp_Mu=0;}
	if (lp_El!=0){delete lp_El; lp_El=0;}

	cout << "All done!!" << endl << endl;
	return true;
}

TChain* loadData(TString fileList, bool isMC){
	//return a TChain linked to the data files
	bool DEBUG = false;
	TChain* tc = new TChain("evt2l");

	ifstream inF(fileList.Data());
	vector<TString> allFiles; double sumW = 0;
	for (std::string line; getline( inF, line ); ){
		if(DEBUG) cout << line << endl;
		if (line[0]=='#') continue;
	   	allFiles.push_back(line);
	}

 	if (false && isMC) // Weight MC trees
 	{
 		cout << "Weighting MC trees from " << fileList << endl;
 		for(auto it=allFiles.begin(); it!=allFiles.end(); )
 		{
 			std::vector<TString> shortList;

 			// Get sampleID
 			auto s = *it;
 			int sampleID = TString(s(s.Index("TeV")+4,6)).Atoi();
 			if(DEBUG) cout << "sampleID: " << sampleID << '\n';
 			int s1 = sampleID;
 			while(sampleID == s1)
 			{
				shortList.push_back(s);
				it++; if (it==allFiles.end()) break;
				s = *it;
				s1 = TString(s(s.Index("TeV")+4,6)).Atoi();
 			} 

	 		double sumW(0);
	 		for(auto s : shortList)
	 		{
				TFile f(s, "read");
				if(!f.IsOpen()){ cout << "Could not open " << s << endl; continue; }

				TH1* h = (TH1*) f.Get("hCutFlow");
				if (!h){ cout << "Could not get hCutFlow histogram for " << s << endl; }
				else sumW += h->GetBinContent(2);
				f.Close(0);
	 		}
	 		double xSecxEff = 1;
	 		if (DEBUG) cout << "Getting cross section now..." << endl;
	 		if (DEBUG && xsecDB) cout << "database is initialized" << endl; 
	 		xSecxEff *= xsecDB->xsectTimesEff(sampleID);
	 		if (DEBUG) cout << "xSecxEff = " << xSecxEff << endl;
	 		double mcEvtW = xSecxEff * LUMI / sumW;

	 		for(auto s : shortList)
	 		{
	 			TFile f(s, "read");
	 			if(!f.IsOpen()){ cout << "Could not open " << s << endl; continue; }

	 			TTree *t = (TTree*) f.Get("evt2l");

	 			TString fname= TString(s(0,s.Length()-5));
	 			fname+="_WEIGHTED.root";

	 			TFile *newFile = new TFile(fname, "recreate");
	 			TTree *newTree = t->CloneTree(-1, "fast");

	 			newTree->SetWeight(mcEvtW);

	 			newTree->SetDirectory(newFile);
	 			newFile->Write();
	 			newFile->Close();
	 		}
 		} 
 		cout << "Trees weighted" << endl;		
 	}

 	cout << "Loading trees from " << fileList << endl;
 	for(auto f : allFiles)
 	{
 		if(isMC) {
 			TString fname= TString(f(0,f.Length()-5));
 			fname+="_WEIGHTED.root";
 			tc->Add(fname);
 		} else tc->Add(f);
 	}

	return tc;
}
