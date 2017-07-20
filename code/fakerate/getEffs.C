/*

getEffs.C
A ROOT macro to derive efficiencies for the matrix method of
non-prompt lepton estimation for the 2LSS C1N2->Wh analysis.

Written by [Gabriel Gallardo](mailto:gabriel.gallardo@cern.ch)
June-July 2017

Works on the evt2l trees produced by the multiLepSearch package found at
https://gitlab.cern.ch/hku/SUSY2L

Execution:

	```
	setupATLAS
	rcSetup
	rc compile # Make sure SUSYTools is there
	root -l -b -q getEffs.C+
	```

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
#include <TDirectory.h>

using std::cout;
using std::endl;

// ======= CONFIG VARIABLES ======== //
const double ptBins[] = {20, 40, 80, 150};
const int nPtBins = sizeof(ptBins)/sizeof(double)-1;
const double etaBins[] = {0, 1.37, 2.5};
const int nEtaBins = sizeof(etaBins)/sizeof(double)-1;
const double LUMI=36074.56;

// const TString inFileData="../support/inFileList-dataTest.txt";
const TString inFileData="../support/inFileList-data.txt";

// TString inFileBBCC="../support/inFileList-.txt";
const TString inFileTTbar="../support/inFileList-ttbar.txt";
const TString inFileZjets="../support/inFileList-Zjets.txt";
// const TString inFileZjets="../support/inFileList-ZPowheg.txt";
const TString inFileWjets="../support/inFileList-Wjets.txt";
const TString inFileVV="../support/inFileList-VV.txt";
const TString inFileVgamma="../support/inFileList-Vgamma.txt";
const TString inFileSingleTop="../support/inFileList-SingleTop.txt";
const TString inFileMC="../support/inFileList-allMC.txt";

const bool DEBUG=false;
const bool VERBOSE=DEBUG || true;

const TString outDir="output/";

const bool measureUnweightedRates = true;
const bool measureLepComp = true;
const bool measureSF = true;
const bool doDataTP = true;

bool passRatesCR(evt2l* tree){
	/* RATES MEASUREMENT CONTROL REGION

	- Exactly two baseline leptons
		- May or may not be signal
		- No sign requirement
	- No more than three jets
	- No b-jets
	- Relative missing ET < 70?
	*/

	if (!measureUnweightedRates) return false;
	// Exactly two leptons
	bool pass = true;
	pass *= tree->leps_==2;
	pass *= tree->jets_<=3;
	for(int i=0; i<tree->jets_; i++) pass *= !(tree->jets_jFlag[i] & JT_BJET);

	return pass;
}

bool passLepCompCR(evt2l* tree)
{
	/* LEPTON COMPOSOTION CONTROL REGION
	- Exactly two leptons
		- Same sign
		- Signal
	- No more than three jets
	- No b-jets
	- Relative missing ET < 70?

	Recall f= N(tight)/N(loose), should require signal leptons
	*/

	if (!measureLepComp) return false;
	bool pass = true;
	pass *= tree->leps_==2;
	pass *= (tree->leps_lFlag[0] & IS_SIGNAL) && (tree->leps_lFlag[1] & IS_SIGNAL);
	pass *= tree->leps_ID[0]*tree->leps_ID[1]>0;
	pass *= tree->jets_<=3;

	// // No b-jets
	for(int i(0); i>tree->jets_;i++) pass*= !(tree->jets_jFlag[i] & JT_BJET);

	return pass;
}

// ======= INFRASTRUCTURE ===== //

SUSY::CrossSectionDB *xsecDB;

// Output histogram file structure
TFile* outFile;
TDirectory* dir_weightedRates;
TDirectory* dir_unweightedRatesMeasurement;
TDirectory* dir_unweightedRatesMeasurement_rates;
TDirectory* dir_unweightedRatesMeasurement_nEvents;
TDirectory* dir_sfMeasurement;
TDirectory* dir_sfMeasurement_scaleFactors;
TDirectory* dir_sfMeasurement_effs;
TDirectory* dir_sfMeasurement_nEvents;
TDirectory* dir_lepCompMeasurement;
TDirectory* dir_lepCompMeasurement_lepComp;
TDirectory* dir_lepCompMeasurement_nEvents;

TH2D* hPrototype= new TH2D("hPrototype", ";p_{T} [GeV];|#eta|", nPtBins, ptBins, nEtaBins, etaBins);

enum LEP_SOURCE {
	REAL=-1,
	HEAVY=0,
	LIGHT=1,
	CONV=2,
	// OTHER=3, // To check lepton composition
};
const int N_FAKES_SOURCE=3;

const TString LS_TOSTRING(LEP_SOURCE s)
{
	if(s==REAL) return "Real";
	if(s==HEAVY) return "Heavy";
	if(s==LIGHT) return "Light";
	if(s==CONV) return "Conv";
	// if(s==OTHER) return "Other";
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

enum LEP_TYPE {
	ELEC,
	MUON,
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
	binsOfSource fakeSumW2; 
	binsOfSource nFakeEntries;

	binsOfSource reals; 
	binsOfSource realSumW2; 
	binsOfSource nRealEntries;

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

	bool NormalizeVec(binsOfSource &v, binsOfSource &vSumW2)
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
				{v[s][p][ptBin][etaBin] /= sumW; vSumW2[s][p][ptBin][etaBin] /= (sumW*sumW);}
		}
		
		return true;
	}

public:
	LepProp():LepProp("test"){};
	LepProp(TString name){
		this->name = name;
		AddSourceToVec(reals);
		for(int i=0; i<N_FAKES_SOURCE; i++) AddSourceToVec(fakes);

		// Initialize sumW vectors
		realSumW2=reals;
		fakeSumW2=fakes;

		// Initialize counters for uncertainties
		nRealEntries = reals;
		nFakeEntries = fakes;

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
			realSumW2[0][p][ptBin][etaBin] += w*w;
			nRealEntries[0][p][ptBin][etaBin] += 1;
		}
		else {
			if (DEBUG) cout << "### Filling fakes for " << LP_TOSTRING(p) << " " << LS_TOSTRING(s) << endl;
			fakes[s][p][ptBin][etaBin]+=w;
			fakeSumW2[s][p][ptBin][etaBin] += w*w;
			nFakeEntries[s][p][ptBin][etaBin] += 1;
		}
		if (DEBUG) cout << "### LepProp filled. Return now." << endl;

		return true;
	}

	bool Normalize()
	{
		if(normalized) return true;
		NormalizeVec(reals, realSumW2);
		NormalizeVec(fakes, fakeSumW2);
		return true;
	};

	bool Draw()
	{

		Normalize(); // Normalize histograms first

		// PROTOTYPE HISTOGRAMS
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

				TString fakesName = (TString) "h" + name + (TString) "_FakeComp_" + ptEtaStr;
				tmp.hFake = (TH2*) hFakesPrototype2D->Clone(fakesName);
				tmp.hFake->SetTitle(ptEtaTitle.Prepend(tmp.hFake->GetTitle()));
				tmp.hFake->SetDirectory(dir_lepCompMeasurement_lepComp);
				tmp.hFake->GetZaxis()->SetRangeUser(0,1);

				TString realsName = (TString) "h" + name + (TString) "_RealComp_" + ptEtaStr;
				tmp.hReal = (TH1*) hRealsPrototype->Clone(realsName);
				tmp.hReal->SetTitle(ptEtaTitle.Prepend(tmp.hReal->GetTitle()));
				tmp.hReal->SetDirectory(dir_lepCompMeasurement_lepComp);
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
			double statUnc = TMath::Sqrt(realSumW2[0][p][ptBin][etaBin]);
			histVector[ptBin][etaBin].hReal->SetBinError(p+1, statUnc);

			for(int s=0;s<N_FAKES_SOURCE; s++)
			{
				histVector[ptBin][etaBin].hFake->Fill(p, s, fakes[s][p][ptBin][etaBin]);
				double statUnc = TMath::Sqrt(fakeSumW2[s][p][ptBin][etaBin]);
				histVector[ptBin][etaBin].hFake->SetBinError(p+1, s+1, statUnc);
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
			c->Print(((TString)".pdf").Prepend(histVector[ptBin][etaBin].hReal->GetName()));

			// Fake composition 2D
			histVector[ptBin][etaBin].hFake->Draw("colz text");
			c->Print(((TString)".pdf").Prepend(histVector[ptBin][etaBin].hFake->GetName()));

			// Fake composition by process
			histVector[ptBin][etaBin].hFakeProc = histVector[ptBin][etaBin].hFake->ProjectionX();
			histVector[ptBin][etaBin].hFakeProc->SetDirectory(dir_lepCompMeasurement_lepComp);
			histVector[ptBin][etaBin].hFakeProc->SetName(((TString)"_byProcess").Prepend(histVector[ptBin][etaBin].hFake->GetName()));
			histVector[ptBin][etaBin].hFakeProc->SetTitle(((TString)" by process").Prepend(histVector[ptBin][etaBin].hFake->GetTitle()));
			histVector[ptBin][etaBin].hFakeProc->GetXaxis()->SetRangeUser(0, N_PROC);
			histVector[ptBin][etaBin].hFakeProc->GetYaxis()->SetRangeUser(0, 1);
			histVector[ptBin][etaBin].hFakeProc->SetFillColor(kAzure+10);
			histVector[ptBin][etaBin].hFakeProc->Draw("bar text");
			c->Print(((TString)".pdf").Prepend(histVector[ptBin][etaBin].hFakeProc->GetName()));

			// Fake composition by source
			histVector[ptBin][etaBin].hFakeSource = histVector[ptBin][etaBin].hFake->ProjectionY();
			histVector[ptBin][etaBin].hFakeSource->SetDirectory(dir_lepCompMeasurement_lepComp);
			histVector[ptBin][etaBin].hFakeSource->SetName(((TString)"_bySource").Prepend(histVector[ptBin][etaBin].hFake->GetName()));
			histVector[ptBin][etaBin].hFakeSource->SetTitle(((TString)" by source").Prepend(histVector[ptBin][etaBin].hFake->GetTitle()));
			histVector[ptBin][etaBin].hFakeSource->GetXaxis()->SetRangeUser(0, N_FAKES_SOURCE);
			histVector[ptBin][etaBin].hFakeSource->GetYaxis()->SetRangeUser(0, 1);
			histVector[ptBin][etaBin].hFakeSource->SetFillColor(kAzure+10);
			histVector[ptBin][etaBin].hFakeSource->Draw("bar text");
			c->Print(((TString)".pdf").Prepend(histVector[ptBin][etaBin].hFakeSource->GetName()));
		}

		return true;
	} 

	std::pair<double, double> GetWeight(LEP_SOURCE s, LEP_PROC p, double pt, double eta)
	{
		int ptBin = GetPtBin(pt); int etaBin = GetEtaBin(eta);
		return GetWeight(s, p, ptBin, etaBin);
	}

	std::pair<double, double> GetWeight(LEP_SOURCE s, LEP_PROC p, int ptBin, int etaBin)
	{
		if(ptBin<0 || etaBin<0) return std::make_pair(0,0);
		if(s==REAL) return std::make_pair(reals[0][p][ptBin][etaBin], TMath::Sqrt(realSumW2[0][p][ptBin][etaBin]));
		else return std::make_pair(fakes[s][p][ptBin][etaBin], TMath::Sqrt(fakeSumW2[s][p][ptBin][etaBin]));
	}

	static int GetPtBin(double pt)
	{
		if (pt<ptBins[0]) return -1; // Underflow
		for(int i=0; i<nPtBins; i++) if (pt<ptBins[i+1]) return i;
		return nPtBins-1; // Overflow, place in largest ptBin.
	}

	static int GetEtaBin(double eta)
	{
		eta = TMath::Abs(eta);
		if (eta<etaBins[0]) return -1; // Underflow, not really possible for eta
		for(int i=0; i<nEtaBins; i++) if (eta<etaBins[i+1]) return i;
		return nEtaBins-1; // Overflow, place in largest etaBin.
	}
};

std::vector<Histos*> allHistos;

class Histos
{
private:
	TH2D* NTight;
	TH2D* NLoose;
	TH2D* Eff;
	bool ignoreTagTag=false;
	bool EffCalced=false;

public:
	Histos(): Histos("Test") {}
	Histos(TString name)
	{
		NTight = (TH2D*) hPrototype->Clone(name+"_NTight"); NTight->SetTitle(name+"_NTight"); NTight->Sumw2();
		NLoose = (TH2D*) hPrototype->Clone(name+"_NLoose"); NLoose->SetTitle(name+"_NLoose"); NLoose->Sumw2();
		Eff = (TH2D*) hPrototype->Clone(name+"_Eff"); Eff->SetTitle(name+"_Eff"); Eff->Sumw2();

		NTight->SetDirectory(dir_unweightedRatesMeasurement_nEvents);
		NLoose->SetDirectory(dir_unweightedRatesMeasurement_nEvents);
		Eff->SetDirectory(dir_unweightedRatesMeasurement_rates);

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
			double pt = lep2.pt > ptBins[nPtBins-1]? ptBins[nPtBins-1]-1 : lep2.pt;
			NLoose->Fill(pt, fabs(lep2.eta), lep2.w);
			if(lep2.isTight) NTight->Fill(pt, fabs(lep2.eta), lep2.w);
			// return true;
		}
		if (lep2.isTag)
		{
			double pt = lep1.pt > ptBins[nPtBins-1]? ptBins[nPtBins-1]-1 : lep1.pt;
			NLoose->Fill(pt, fabs(lep1.eta), lep1.w);
			if(lep1.isTight) NTight->Fill(pt, fabs(lep1.eta), lep1.w);
			// return true;
		}
		return true;
	}

	bool calcEff() 
	{
		if (EffCalced) return true;
		Eff->Add(NTight); Eff->Divide(NLoose); 
		EffCalced=true; 
		return true;
	}

	void Draw()
	{
		TCanvas can;

		Eff->Draw("text colz");
		can.Print(TString(Eff->GetName())+".pdf");

		// Eff->ProfileX()->Draw();
		// can.Print(TString(Eff->GetName())+"Pt.pdf");

		// Eff->ProfileY()->Draw();
		// can.Print(TString(Eff->GetName())+"Eta.pdf");

		NTight->Draw("text");
		can.Print(TString(NTight->GetName())+".pdf");

		NLoose->Draw("text");
		can.Print(TString(NLoose->GetName())+".pdf");

	}

	void SetHistDir(TDirectory* nDir, TDirectory* effDir)
	{
		NTight->SetDirectory(nDir);
		NLoose->SetDirectory(nDir);
		Eff->SetDirectory(effDir);
	}

	void vetoTagTag(bool b=true){ignoreTagTag=b;}

	TH2* GetEffHist() {
		calcEff();
		return Eff;
	}

	std::pair<double, double> GetEff(double pt, double eta)
	{
		int ptBin = LepProp::GetPtBin(pt); 
		int etaBin = LepProp::GetPtBin(eta);
		return GetEff(ptBin, etaBin);
	}

	std::pair<double, double> GetEff(int ptBin, int etaBin)
	{
		return std::make_pair(
			Eff->GetBinContent(ptBin, etaBin),
			Eff->GetBinError(ptBin, etaBin)
		);
	}

};

class SFhist
{
private:
	TString name;

	TH2 const* hData;
	TH2 const* hMC;
	TH2* hRatio;

	// SFhist(): SFhist("hist"){};
	// SFhist(TString name) : SFhist(name, 0, 0) {};
public:
	
	SFhist(TString n, Histos* d, Histos* m):
		name(n), hData(d->GetEffHist()), hMC(m->GetEffHist())
	{
		int nameLength = ((TString)hData->GetName()).Length();
		TString histName = (TString)(((TString) hData->GetName())(0, nameLength-8));
		histName += "SF";

		hRatio = (TH2*) hData->Clone(histName); 
		hRatio->SetTitle(histName);

		hRatio->Divide(hMC);

		hRatio->SetDirectory(dir_sfMeasurement_scaleFactors);

        Draw();
	}

	void Draw()
	{
		TCanvas *c = new TCanvas();
		hRatio->Draw("colz text");
		c->Print(hRatio->GetName()+(TString)".pdf");
	}

	std::pair<double, double> GetSF(double pt, double eta)
	{
		int ptBin = LepProp::GetPtBin(pt);
		int etaBin = LepProp::GetEtaBin(eta);
		return GetSF(ptBin, etaBin);
	}

	std::pair<double, double> GetSF(int ptBin, int etaBin)
	{
		return std::make_pair(
			hRatio->GetBinContent(ptBin, etaBin),
			hRatio->GetBinError(ptBin, etaBin)
		);
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

Histos* hMu_Real_Data;
Histos* hEl_Real_Data;
Histos* hMu_Heavy_Data;
Histos* hEl_Heavy_Data;
Histos* hEl_Conv_Data;


TH2* hEl_final_fakeRate;
TH2* hMu_final_fakeRate;
TH2* hEl_final_realRate;
TH2* hMu_final_realRate;

std::vector< std::vector<Histos*> > hElReal;
std::vector< std::vector<Histos*> > hElFake;
std::vector< std::vector<Histos*> > hMuReal;
std::vector< std::vector<Histos*> > hMuFake;

std::vector<SFhist*> sfhEl_FakeScaleFactors;
std::vector<SFhist*> sfhMu_FakeScaleFactors;
SFhist* hEl_RealScaleFactors;
SFhist* hMu_RealScaleFactors;
// ======= FUNCTION PROTOTYPES ===== //

int getEffs();
bool initialize();
TChain* loadData(TString, bool isMC=false);
bool finalize();
bool doTP(evt2l* tree, Histos* hMu_Real, Histos* hEl_Real, Histos* hMu_Heavy, Histos* hEl_Heavy, Histos* hEl_Conv, double w=1);
bool loopMC(evt2l* tree, Histos *hMu_Real, Histos *hEl_Real, Histos *hMu_Heavy, Histos *hEl_Heavy, Histos *hMu_Light, Histos *hEl_Light, Histos *hEl_Conv, LEP_PROC p);
LEP_SOURCE castSource(MCTC::ParticleType type, MCTC::ParticleOrigin orig, LEP_TYPE l);
void labelMCTChist(TH2* h);
void calcFinalEffs();

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
	cout << "Begin initialize()" <<endl;

	// Set output directory
	if(!gSystem->cd(outDir)){
		gSystem->mkdir(outDir);
		gSystem->cd(outDir);
	}

	// Drawing options
	gStyle->SetOptStat(0); 
	gStyle->SetPalette(kLightTerrain);

	// Define output file and structure
	outFile = new TFile("output.root", "recreate");
	// {
		dir_weightedRates = outFile->mkdir("weightedRates");

		dir_unweightedRatesMeasurement = outFile->mkdir("unweightedRatesMeasurement");
		dir_unweightedRatesMeasurement_rates = dir_unweightedRatesMeasurement->mkdir("rates");
		dir_unweightedRatesMeasurement_nEvents = dir_unweightedRatesMeasurement->mkdir("nEvents");

		dir_sfMeasurement = outFile->mkdir("sfMeasurement");
		dir_sfMeasurement_scaleFactors = dir_sfMeasurement->mkdir("scaleFactors");
		dir_sfMeasurement_effs = dir_sfMeasurement->mkdir("effs");
		dir_sfMeasurement_nEvents = dir_sfMeasurement->mkdir("nEvents");

		dir_lepCompMeasurement = outFile->mkdir("lepCompMeasurement");
		dir_lepCompMeasurement_lepComp = dir_lepCompMeasurement->mkdir("lepComp");
		dir_lepCompMeasurement_nEvents = dir_lepCompMeasurement->mkdir("nEvents");
	// }

	// Load cross section database
	gROOT->Macro("$ROOTCOREDIR/scripts/load_packages.C");
	xsecDB = new SUSY::CrossSectionDB(gSystem->ExpandPathName("$ROOTCOREBIN/data/SUSYTools/mc15_13TeV/"));
	if (!xsecDB){
	   cout << "CrossSectionDB could not be loaded" << endl;
	   return false;
	} else cout << "CrossSectionDB loaded" << endl;

	// // //  Set up some histograms // // // 

	// Unweighted efficiencies 
	std::vector<Histos*> v(N_PROC, 0);
	hElFake.reserve(N_FAKES_SOURCE);
	for(int i=0; i<N_FAKES_SOURCE; i++) {hElFake.push_back(v);}
	hElReal.reserve(1); hElReal.push_back(v);

	hMuFake.reserve(N_FAKES_SOURCE);
	for(int i=0; i<N_FAKES_SOURCE; i++) {hMuFake.push_back(v);}
	hMuReal.reserve(1); hMuReal.push_back(v);

	// Lepton composition
	hMu_MCTC = new TH2D("hMu_MCTC", "MCTC Mu;Type;Origin", PARTICLETYPES, 0, PARTICLETYPES, PARTICLEORIGIN, 0, PARTICLEORIGIN);
	hEl_MCTC = new TH2D("hEl_MCTC", "MCTC El;Type;Origin", PARTICLETYPES, 0, PARTICLETYPES, PARTICLEORIGIN, 0, PARTICLEORIGIN);
	hMu_MCTC->SetDirectory(dir_lepCompMeasurement_nEvents);
	hEl_MCTC->SetDirectory(dir_lepCompMeasurement_nEvents);

	lp_El = new LepProp("El");
	lp_Mu = new LepProp("Mu");

	// Tag-and-probe histograms for MC
	hMu_Real_MCTP = new Histos("hMu_Real_MCTP"); hMu_Real_MCTP->SetHistDir(dir_sfMeasurement_nEvents, dir_sfMeasurement_effs);
	hEl_Real_MCTP = new Histos("hEl_Real_MCTP"); hEl_Real_MCTP->SetHistDir(dir_sfMeasurement_nEvents, dir_sfMeasurement_effs);
	hMu_Heavy_MCTP = new Histos("hMu_Heavy_MCTP"); hMu_Heavy_MCTP->SetHistDir(dir_sfMeasurement_nEvents, dir_sfMeasurement_effs);
	hEl_Heavy_MCTP = new Histos("hEl_Heavy_MCTP"); hEl_Heavy_MCTP->SetHistDir(dir_sfMeasurement_nEvents, dir_sfMeasurement_effs);
	hMu_Light_MCTP = new Histos("hMu_Light_MCTP"); hMu_Light_MCTP->SetHistDir(dir_sfMeasurement_nEvents, dir_sfMeasurement_effs);
	hEl_Light_MCTP = new Histos("hEl_Light_MCTP"); hEl_Light_MCTP->SetHistDir(dir_sfMeasurement_nEvents, dir_sfMeasurement_effs);
	hEl_Conv_MCTP = new Histos("hEl_Conv_MCTP"); hEl_Conv_MCTP->SetHistDir(dir_sfMeasurement_nEvents, dir_sfMeasurement_effs);

	// Tag-and-probe histograms for Data
	hMu_Real_Data = new Histos("hMu_Real_Data"); hMu_Real_Data->SetHistDir(dir_sfMeasurement_nEvents, dir_sfMeasurement_effs);
	hEl_Real_Data = new Histos("hEl_Real_Data"); hEl_Real_Data->SetHistDir(dir_sfMeasurement_nEvents, dir_sfMeasurement_effs);
	hMu_Heavy_Data = new Histos("hMu_Heavy_Data"); hMu_Heavy_Data->SetHistDir(dir_sfMeasurement_nEvents, dir_sfMeasurement_effs);
	hEl_Heavy_Data = new Histos("hEl_Heavy_Data"); hEl_Heavy_Data->SetHistDir(dir_sfMeasurement_nEvents, dir_sfMeasurement_effs);
	hEl_Conv_Data = new Histos("hEl_Conv_Data"); hEl_Conv_Data->SetHistDir(dir_sfMeasurement_nEvents, dir_sfMeasurement_effs);

	// SF histograms
	sfhEl_FakeScaleFactors.reserve(N_FAKES_SOURCE);
	sfhMu_FakeScaleFactors.reserve(N_FAKES_SOURCE);
	for(int i=0; i<N_FAKES_SOURCE; i++) {sfhEl_FakeScaleFactors[i]=0; sfhMu_FakeScaleFactors[i]=0;}
	hEl_RealScaleFactors = 0;
	hMu_RealScaleFactors = 0;

	cout << "initialize() complete" << endl;
	return true;
}

int getEffs()
{
	initialize(); 
	TChain* tcMC = new TChain("evt2l");


	// BBBAR/CCBAR
	// evt2l* bbccEvts = new evt2l(loadData(inFileBBCC, true)); // tcMC->Add((TChain*) bbccEvts->fChain);
	// hMuReal[0][BBCC]= new Histos hMu_Real_BBCC("hMu_Real_BBCC");
	// hElReal[0][BBCC]= new Histos hEl_Real_BBCC("hEl_Real_BBCC");
	// hMuFake[HEAVY][BBCC]= new Histos hMu_Heavy_BBCC("hMu_Heavy_BBCC");
	// hElFake[HEAVY][BBCC]= new Histos hEl_Heavy_BBCC("hEl_Heavy_BBCC");
	// hMuFake[LIGHT][BBCC]= new Histos hMu_Light_BBCC("hMu_Light_BBCC");
	// hElFake[LIGHT][BBCC]= new Histos hEl_Light_BBCC("hEl_Light_BBCC");
	// hElFake[CONV][BBCC]= new Histos hEl_Conv_BBCC("hEl_Conv_BBCC");
	// cout << "Begin loop over bb/cc: " << bbccEvts->fChain->GetEntries() << " events" << endl;
	// loopMC(bbccEvts, hMuReal[0][BBCC], hElReal[0][BBCC], hMuFake[HEAVY][BBCC], hElFake[HEAVY][BBCC], hMuFake[LIGHT][BBCC], hElFake[LIGHT][BBCC], hElFake[CONV][BBCC], BBCC);
	// cout << "bb/cc loop finished" << endl << endl;
	//
	
	// TTBAR
	evt2l* ttEvts = new evt2l(loadData(inFileTTbar, true));  tcMC->Add((TChain*) ttEvts->fChain);
	hMuReal[0][TTBAR] = new Histos("hMu_Real_TT");
	hElReal[0][TTBAR] = new Histos("hEl_Real_TT");
	hMuFake[HEAVY][TTBAR] = new Histos("hMu_Heavy_TT");
	hElFake[HEAVY][TTBAR] = new Histos("hEl_Heavy_TT");
	hMuFake[LIGHT][TTBAR] = new Histos("hMu_Light_TT");
	hElFake[LIGHT][TTBAR] = new Histos("hEl_Light_TT");
	hElFake[CONV][TTBAR] = new Histos("hEl_Conv_TT");
	cout << "Begin loop over ttBar: " << ttEvts->fChain->GetEntries() << " events" << endl;
	loopMC(ttEvts, hMuReal[0][TTBAR], hElReal[0][TTBAR], hMuFake[HEAVY][TTBAR], hElFake[HEAVY][TTBAR], hMuFake[LIGHT][TTBAR], hElFake[LIGHT][TTBAR], hElFake[CONV][TTBAR], TTBAR);
	cout << "ttBar loop finished" << endl << endl;


	// // Z+jets
	evt2l* zJetsEvts = new evt2l(loadData(inFileZjets, true));  tcMC->Add((TChain*) zJetsEvts->fChain);
	hMuReal[0][ZJETS] = new Histos("hMu_Real_Zjets");
	hElReal[0][ZJETS] = new Histos("hEl_Real_Zjets");
	hMuFake[HEAVY][ZJETS] = new Histos("hMu_Heavy_Zjets");
	hElFake[HEAVY][ZJETS] = new Histos("hEl_Heavy_Zjets");
	hMuFake[LIGHT][ZJETS] = new Histos("hMu_Light_Zjets");
	hElFake[LIGHT][ZJETS] = new Histos("hEl_Light_Zjets");
	hElFake[CONV][ZJETS] = new Histos("hEl_Conv_Zjets");
	cout << "Begin loop over Z+jets: " << zJetsEvts->fChain->GetEntries() << " events" << endl;
	loopMC(zJetsEvts, hMuReal[0][ZJETS], hElReal[0][ZJETS], hMuFake[HEAVY][ZJETS], hElFake[HEAVY][ZJETS], hMuFake[LIGHT][ZJETS], hElFake[LIGHT][ZJETS], hElFake[CONV][ZJETS], ZJETS);
	cout << "Z+jets loop finished" << endl << endl;

	// // W+jets
	evt2l* wJetsEvts = new evt2l(loadData(inFileWjets, true));  tcMC->Add((TChain*) wJetsEvts->fChain);
	hMuReal[0][WJETS] = new Histos("hMu_Real_Wjets");
	hElReal[0][WJETS] = new Histos("hEl_Real_Wjets");
	hMuFake[HEAVY][WJETS] = new Histos("hMu_Heavy_Wjets");
	hElFake[HEAVY][WJETS] = new Histos("hEl_Heavy_Wjets");
	hMuFake[LIGHT][WJETS] = new Histos("hMu_Light_Wjets");
	hElFake[LIGHT][WJETS] = new Histos("hEl_Light_Wjets");
	hElFake[CONV][WJETS] = new Histos("hEl_Conv_Wjets");
	cout << "Begin loop over W+jets: " << wJetsEvts->fChain->GetEntries() << " events" << endl;
	loopMC(wJetsEvts, hMuReal[0][WJETS], hElReal[0][WJETS], hMuFake[HEAVY][WJETS], hElFake[HEAVY][WJETS], hMuFake[LIGHT][WJETS], hElFake[LIGHT][WJETS], hElFake[CONV][WJETS], WJETS);
	cout << "W+jets loop finished" << endl << endl;

	// DIBOSON
	evt2l* vvEvts = new evt2l(loadData(inFileVV, true));  tcMC->Add((TChain*) vvEvts->fChain);
	hMuReal[0][VV] = new Histos ("hMu_Real_VV");
	hElReal[0][VV] = new Histos ("hEl_Real_VV");
	hMuFake[HEAVY][VV] = new Histos ("hMu_Heavy_VV");
	hElFake[HEAVY][VV] = new Histos ("hEl_Heavy_VV");
	hMuFake[LIGHT][VV] = new Histos ("hMu_Light_VV");
	hElFake[LIGHT][VV] = new Histos ("hEl_Light_VV");
	hElFake[CONV][VV] = new Histos ("hEl_Conv_VV");
	cout << "Begin loop over VV: " << vvEvts->fChain->GetEntries() << " events" << endl;
	loopMC(vvEvts, hMuReal[0][VV], hElReal[0][VV], hMuFake[HEAVY][VV], hElFake[HEAVY][VV], hMuFake[LIGHT][VV], hElFake[LIGHT][VV], hElFake[CONV][VV], VV);
	cout << "VV loop finished" << endl << endl;

	// VGAMMA
	evt2l* vgammaEvts = new evt2l(loadData(inFileVgamma, true));  tcMC->Add((TChain*) vgammaEvts->fChain);
	hMuReal[0][VGAMMA] = new Histos("hMu_Real_Vgamma");
	hElReal[0][VGAMMA] = new Histos("hEl_Real_Vgamma");
	hMuFake[HEAVY][VGAMMA] = new Histos("hMu_Heavy_Vgamma");
	hElFake[HEAVY][VGAMMA] = new Histos("hEl_Heavy_Vgamma");
	hMuFake[LIGHT][VGAMMA] = new Histos("hMu_Light_Vgamma");
	hElFake[LIGHT][VGAMMA] = new Histos("hEl_Light_Vgamma");
	hElFake[CONV][VGAMMA] = new Histos("hEl_Conv_Vgamma");
	cout << "Begin loop over Vgamma: " << vvEvts->fChain->GetEntries() << " events" << endl;
	loopMC(vgammaEvts, hMuReal[0][VGAMMA], hElReal[0][VGAMMA], hMuFake[HEAVY][VGAMMA], hElFake[HEAVY][VGAMMA], hMuFake[LIGHT][VGAMMA], hElFake[LIGHT][VGAMMA], hElFake[CONV][VGAMMA], VGAMMA);
	cout << "Vgamma loop finished" << endl << endl;

	// SingleTop
	evt2l* singleTopEvts = new evt2l(loadData(inFileSingleTop, true));  tcMC->Add((TChain*) singleTopEvts->fChain);
	hMuReal[0][SINGLETOP] = new Histos("hMu_Real_SingleTop");
	hElReal[0][SINGLETOP] = new Histos("hEl_Real_SingleTop");
	hMuFake[HEAVY][SINGLETOP] = new Histos("hMu_Heavy_SingleTop");
	hElFake[HEAVY][SINGLETOP] = new Histos("hEl_Heavy_SingleTop");
	hMuFake[LIGHT][SINGLETOP] = new Histos("hMu_Light_SingleTop");
	hElFake[LIGHT][SINGLETOP] = new Histos("hEl_Light_SingleTop");
	hElFake[CONV][SINGLETOP] = new Histos("hEl_Conv_SingleTop");
	cout << "Begin loop over single top: " << vvEvts->fChain->GetEntries() << " events" << endl;
	loopMC(singleTopEvts, hMuReal[0][SINGLETOP], hElReal[0][SINGLETOP], hMuFake[HEAVY][SINGLETOP], hElFake[HEAVY][SINGLETOP], hMuFake[LIGHT][SINGLETOP], hElFake[LIGHT][SINGLETOP], hElFake[CONV][SINGLETOP], SINGLETOP);
	cout << "Single top loop finished" << endl << endl;


	// DATA TAG-AND-PROBE
	if(doDataTP){
		evt2l* dataEvts = new evt2l(loadData(inFileData));
		cout << "Executing tag-and-probe on data" << endl;
		int nDataEntries = dataEvts->fChain->GetEntries();
		cout << "Begin loop over data:" << nDataEntries << " events" << endl;
		for(int i=0; i<nDataEntries; i++)
		{
			loadbar(i+1, nDataEntries);
			dataEvts->GetEntry(i);
			if (DEBUG) cout << "Event no " << i << endl;
			doTP(dataEvts, hMu_Real_Data, hEl_Real_Data, hMu_Heavy_Data, hEl_Heavy_Data, hEl_Conv_Data);
		}
		cout << endl;
		cout << "Data tag-and-probe finished" << endl;
	}

	if (measureSF) 
	{
		sfhEl_FakeScaleFactors[HEAVY] = new SFhist("El_Heavy", hEl_Heavy_Data, hEl_Heavy_MCTP);
		// sfhEl_FakeScaleFactors[LIGHT] = new SFhist("El_Light", hEl_Light_Data, hEl_Light_MCTP);
		sfhEl_FakeScaleFactors[CONV] = new SFhist("El_Conv", hEl_Conv_Data, hEl_Conv_MCTP);

		sfhMu_FakeScaleFactors[HEAVY] = new SFhist("Mu_Heavy", hMu_Heavy_Data, hMu_Heavy_MCTP);
		// sfhMu_FakeScaleFactors[LIGHT] = new SFhist("Mu_Light", hMu_Light_Data, hMu_Light_MCTP);

		hEl_RealScaleFactors = new SFhist("El_Real", hEl_Real_Data, hEl_Real_MCTP);
		hMu_RealScaleFactors = new SFhist("Mu_Real", hMu_Real_Data, hMu_Real_MCTP);
	}

	if (measureSF && measureLepComp && measureUnweightedRates) calcFinalEffs();
	// */
	finalize();

	return 0;
}

bool doTP(evt2l* tree, Histos* hMu_Real, Histos* hEl_Real, Histos* hMu_Heavy, Histos* hEl_Heavy, Histos* hEl_Conv, double w)
{
	// const bool DEBUG = true;

	if (DEBUG) cout << "Enter doTP()" << endl;
	// if (tree->sig_trigCode==0) return false;
	if (tree->leps_==2)
	{			
		// mm and ee
		bool mumu = TMath::Abs(int(tree->leps_ID[0]/1000)*int(tree->leps_ID[1]/1000))==169;
		bool elel = TMath::Abs(int(tree->leps_ID[0]/1000)*int(tree->leps_ID[1]/1000))==121;
		bool isOS = (tree->leps_ID[0]*tree->leps_ID[1])<0;

		// "Tight" defined as signal electron
		bool lep1isTight = (tree->leps_lFlag[0] & IS_SIGNAL);
		bool lep2isTight = (tree->leps_lFlag[1] & IS_SIGNAL);

		/** UNWEIGHTED REAL EFFICIENCIES MEASURED BY Z_ll TAG-AND-PROBE IN DATA **/
		if(fabs(tree->l12_m - 91)<10 && isOS) // Select Z mass pair
		{
			if (DEBUG) cout << "Zll\t";
			// Tag requirement: pT >25, passes signal ("tight") cut, passes trigger, and is trigger matched
			bool lep1isTag = (tree->leps_pt[0]>25 && lep1isTight && (TMath::Abs(tree->leps_ID[0]) & 1) && (TMath::Abs(tree->leps_ID[0]) & 2));
			bool lep2isTag = (tree->leps_pt[1]>25 && lep2isTight && (TMath::Abs(tree->leps_ID[1]) & 1) && (TMath::Abs(tree->leps_ID[1]) & 2));

			if (DEBUG) cout << "Lep 1: " << (lep1isTag?"*":" ") << '\t'
							<< "Lep 2: " << (lep2isTag?"*":" ") << '\t';

			lepInfo l1info = {tree->leps_pt[0], tree->leps_eta[0], w, lep1isTag, lep1isTight};
			lepInfo l2info = {tree->leps_pt[1], tree->leps_eta[1], w, lep2isTag, lep2isTight};
			if (mumu) hMu_Real->Fill(l1info, l2info);
			else if (elel) hEl_Real->Fill(l1info, l2info);	
			// if (DEBUG) cout << "Histograms filled" << endl;
		}

		/* UNWEIGHTED HEAVY EFFICIENCIES MEASURED BY ll TAG-AND-PROBE IN DATA */
		/** in W/signal suppressed region **/
		else if (tree->sig_Met<40)
		{
			// To suppress W background and signal
			// if (tree->sig_Met>40) continue; TMath::Sqrt(2*leps[0].pt*sig.Met*(1-TMath::Cos(leps[0].MET_dPhi)))<40

			if (DEBUG) cout << "Heavy lep tag-and-probe" << endl;
			std::pair<double, double> mT_lep_MET;
			mT_lep_MET.first  = TMath::Sqrt(2*tree->leps_pt[0]*tree->sig_Met*(1-TMath::Cos(tree->leps_MET_dPhi[0])));
			mT_lep_MET.second = TMath::Sqrt(2*tree->leps_pt[1]*tree->sig_Met*(1-TMath::Cos(tree->leps_MET_dPhi[1])));
			if(mT_lep_MET.first>40 || mT_lep_MET.second>40) return false;

			bool lep1isTag = (tree->leps_pt[0]>20 && lep1isTight && TMath::Abs(int(tree->leps_ID[0]/1000))==13 && (TMath::Abs(tree->leps_ID[0]) & 1) && (TMath::Abs(tree->leps_ID[0]) & 2));
			bool lep2isTag = (tree->leps_pt[1]>20 && lep2isTight && TMath::Abs(int(tree->leps_ID[1]/1000))==13 && (TMath::Abs(tree->leps_ID[1]) & 1) && (TMath::Abs(tree->leps_ID[1]) & 2));

			if (DEBUG) cout << "Lep 1: " << (lep1isTag?"*":" ") << '\t'
							<< "Lep 2: " << (lep2isTag?"*":" ") << '\t';

			lepInfo l1info = {tree->leps_pt[0], tree->leps_eta[0], w, lep1isTag, lep1isTight};
			lepInfo l2info = {tree->leps_pt[1], tree->leps_eta[1], w, lep2isTag, lep2isTight};

			// Select mumu with mll cut, or emu
			if (mumu && tree->l12_m>40) hMu_Heavy->Fill(l1info, l2info);
			else hEl_Heavy->Fill(l1info, l2info);
			if (DEBUG) cout << "Histograms filled" << endl;

		}		
	}

	/** UNWEIGHTED PHOTON CONVERSTION EFFICIENCIES MEASURED BY mumu-e TAG-AND-PROBE IN DATA **/ 
	else if (false && tree->leps_==3)
	{
		if (DEBUG) cout << "PhotonConv lep tag-and-probe" << endl;

		// Select Z_mumu pair and electron 
		if ((int(tree->leps_ID[0]/1000)*int(tree->leps_ID[1]/1000) != -169)	|| TMath::Abs(int(tree->leps_ID[2])/1000)!= 11)
			return false;

		// Dimuon trigger match
		if (!(tree->sig_trigMask & tree->sig_trigCode)) return false;

		TLorentzVector p1, p2, p3;
		p1.SetPtEtaPhiM(tree->leps_pt[0], tree->leps_eta[0], tree->leps_phi[0], 0.105658);
		p2.SetPtEtaPhiM(tree->leps_pt[1], tree->leps_eta[1], tree->leps_phi[1], 0.105658);
		p3.SetPtEtaPhiM(tree->leps_pt[2], tree->leps_eta[2], tree->leps_phi[2], 0.000511);

		if (TMath::Abs((p1+p2+p3).M()-91)>10) return false;

		bool elIsTight = (tree->leps_lFlag[2] & IS_SIGNAL);

		lepInfo elInfo = {tree->leps_pt[2], tree->leps_eta[2], w, false, elIsTight};
		hEl_Conv->Fill(elInfo);
		if (DEBUG) cout << "Histograms filled" << endl;

	}

	if (DEBUG) cout << "End of doTP()" << endl;
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

			LEP_TYPE recoLepType;
			if (TMath::Abs(int(tree->leps_ID[j]/1000)) == 11) recoLepType=ELEC;
			else if (TMath::Abs(int(tree->leps_ID[j]/1000)) == 13) recoLepType=MUON;
			else continue;

			ParticleType   type = static_cast<ParticleType>(tree->leps_truthType[j]);
			ParticleOrigin orig = static_cast<ParticleOrigin>(tree->leps_truthOrig[j]);
			LEP_SOURCE 	   source = castSource(type, orig, recoLepType);

			// Charge flip. (+lep_ID) * (-pdgId) is no flip
			// if (source==CONV && recoLepType==ELEC && tree->leps_ID[j]*tree->leps_firstEgMotherPdgId[j]>0)
			// 	continue;

			lepInfo l = {tree->leps_pt[j], tree->leps_eta[j], w, true, lepIsTight};

			if (recoLepType==ELEC){ switch (source)
			{
				case REAL: hEl_Real->Fill(l); break;
				case HEAVY: hEl_Heavy->Fill(l); break;
				case LIGHT: hEl_Light->Fill(l); break;
				case CONV: hEl_Conv->Fill(l); break;
				default: break;
			}}

			else if (recoLepType==MUON){ switch (source)
			{
				case REAL: hMu_Real->Fill(l); break;
				case HEAVY: hMu_Heavy->Fill(l); break;
				case LIGHT: hMu_Light->Fill(l); break;
				// case CONV: hMu_Conv->Fill(l); break;
				default: break;
			}}
			if (DEBUG) cout << "# Filled rates histograms for lepton " << j << endl;
		}}

		// Measure composition fractions
		if (passLepCompCR(tree)){ for(int j(0); j<2; j++)
		{
			if(tree->leps_pt[j]<20) continue;

			bool lepIsTight = tree->leps_lFlag[j] & IS_SIGNAL;
			LEP_TYPE recoLepType;
			if (TMath::Abs(int(tree->leps_ID[j]/1000)) == 11) recoLepType=ELEC;
			else if (TMath::Abs(int(tree->leps_ID[j]/1000)) == 13) recoLepType=MUON;
			else continue;

			ParticleType   type = static_cast<ParticleType>(tree->leps_truthType[j]);
			ParticleOrigin orig = static_cast<ParticleOrigin>(tree->leps_truthOrig[j]);
			LEP_SOURCE 	   source = castSource(type, orig, recoLepType);

			if (recoLepType==ELEC)
			{
				lp_El->Fill(source, p, tree->leps_pt[j], tree->leps_eta[j], w);
				hEl_MCTC->Fill(type, orig);
			}
			else if (recoLepType==MUON)
			{
				lp_Mu->Fill(source, p, tree->leps_pt[j], tree->leps_eta[j], w);
				hMu_MCTC->Fill(type, orig);
			}
			if (DEBUG) cout << "# Filled lepton composition histograms for lepton " << j << endl;
		}}

		// Tag and probe for data/MC scale factors
		if(measureSF) doTP(tree, hMu_Real_MCTP, hEl_Real_MCTP, hMu_Heavy_MCTP, hEl_Heavy_MCTP, hEl_Conv_MCTP, w);

		if (DEBUG) cout << "# After lepton loop" << endl;
	}
	cout << endl;
	return true;
}

// /*
void calcFinalEffs()
{
	for(auto h: allHistos) h->calcEff(); 
    lp_El->Normalize();
    lp_Mu->Normalize();

	hEl_final_fakeRate = new TH2D("hEl_final_fakeRate", "El Fake rate;p_{T} [GeV];|#eta|", nPtBins, ptBins, nEtaBins, etaBins);
	hMu_final_fakeRate = new TH2D("hMu_final_fakeRate", "Mu Fake rate;p_{T} [GeV];|#eta|", nPtBins, ptBins, nEtaBins, etaBins);
	hEl_final_realRate = new TH2D("hEl_final_realRate", "El Real rate;p_{T} [GeV];|#eta|", nPtBins, ptBins, nEtaBins, etaBins);
	hMu_final_realRate = new TH2D("hMu_final_realRate", "Mu Real rate;p_{T} [GeV];|#eta|", nPtBins, ptBins, nEtaBins, etaBins);

	hEl_final_fakeRate->SetDirectory(dir_weightedRates);
	hMu_final_fakeRate->SetDirectory(dir_weightedRates);
	hEl_final_realRate->SetDirectory(dir_weightedRates);
	hMu_final_realRate->SetDirectory(dir_weightedRates);

	for(int etaBin=0; etaBin<nEtaBins; etaBin++)
	for(int  ptBin=0;  ptBin< nPtBins;  ptBin++)
	{
		// Initialize
		double elRealRate, elRealUnc, muRealRate, muRealUnc;
		elRealRate = elRealUnc = muRealRate = muRealUnc = 0;
		double elFakeRate, elFakeUnc, muFakeRate, muFakeUnc;
		elFakeRate = elFakeUnc = muFakeRate = muFakeUnc = 0;

		// Each term is Rates are Sum(e). e = (unweighted eff) * (SF) * (LepComp fraction)
		// Uncertainty Sqrt(Sum((Unc_i)^2)), (Unc_i)^2 = e^2 ( Sum(Sq(sigma/val)) ), assuming uncorrelated uncertainties
		for(int p=0; p<N_PROC; p++)
		{
			auto elEff 	= hElReal[0][p]->GetEff(ptBin, etaBin);
			auto elSF 	= std::make_pair(1.,0.); // hEl_RealScaleFactors->GetSF(ptBin, etaBin);
			auto elLP 	= lp_El->GetWeight(REAL, (LEP_PROC) p, ptBin, etaBin);

			elRealRate += elEff.first*elSF.first*elLP.first;
			elRealUnc  += elRealRate*elRealRate * (pow(elEff.second/elEff.first,2)+pow(elSF.second/elSF.first,2)+pow(elLP.second/elLP.first,2));

			auto muEff 	= hMuReal[0][p]->GetEff(ptBin, etaBin);
			auto muSF 	= std::make_pair(1.,0.); // hMu_RealScaleFactors->GetSF(ptBin, etaBin);
			auto muLP 	= lp_Mu->GetWeight(REAL, (LEP_PROC) p, ptBin, etaBin);

			muRealRate += muEff.first*muSF.first*muLP.first;
			muRealUnc  += muRealRate*muRealRate * (pow(muEff.second/muEff.first,2)+pow(muSF.second/muSF.first,2)+pow(muLP.second/muLP.first,2));

			for(int s=0; s<N_FAKES_SOURCE; s++)
			{

				auto elEff 	= hElFake[s][p]->GetEff(ptBin, etaBin);
				auto elSF 	= std::make_pair(1.,0.); // ((LEP_SOURCE) s)==LIGHT? std::make_pair(1.,0.0) : sfhEl_FakeScaleFactors[s]->GetSF(ptBin, etaBin);
				auto elLP 	= lp_El->GetWeight((LEP_SOURCE) s, (LEP_PROC) p, ptBin, etaBin);

				elFakeRate += elEff.first*elSF.first*elLP.first;
				elFakeUnc  += elFakeRate*elFakeRate * (pow(elEff.second/elEff.first,2)+pow(elSF.second/elSF.first,2)+pow(elLP.second/elLP.first,2));

				if(((LEP_SOURCE) s)==CONV) continue;

				auto muEff 	= hMuFake[s][p]->GetEff(ptBin, etaBin);
				auto muSF 	= std::make_pair(1.,0.); // ((LEP_SOURCE) s)==LIGHT? std::make_pair(1.,0.0) : sfhMu_FakeScaleFactors[s]->GetSF(ptBin, etaBin);
				auto muLP 	= lp_Mu->GetWeight((LEP_SOURCE) s, (LEP_PROC) p, ptBin, etaBin);

				muFakeRate += muEff.first*muSF.first*muLP.first;
				muFakeUnc  += muFakeRate*muFakeRate * (pow(muEff.second/muEff.first,2)+pow(muSF.second/muSF.first,2)+pow(muLP.second/muLP.first,2));
			}
		}

		elRealUnc = TMath::Sqrt(elRealUnc);
		muRealUnc = TMath::Sqrt(muRealUnc);

		elFakeUnc = TMath::Sqrt(elFakeUnc);
		muFakeUnc = TMath::Sqrt(muFakeUnc);

		// Fill histograms
		hEl_final_fakeRate->SetBinContent(ptBin, etaBin, elFakeRate); hEl_final_fakeRate->SetBinError(ptBin, etaBin, elFakeUnc);
		hMu_final_fakeRate->SetBinContent(ptBin, etaBin, muFakeRate); hMu_final_fakeRate->SetBinError(ptBin, etaBin, muFakeUnc);
		hEl_final_realRate->SetBinContent(ptBin, etaBin, elRealRate); hEl_final_realRate->SetBinError(ptBin, etaBin, elRealUnc);
		hMu_final_realRate->SetBinContent(ptBin, etaBin, muRealRate); hMu_final_realRate->SetBinError(ptBin, etaBin, muRealUnc);
	}
    TCanvas c;
	hEl_final_fakeRate->Draw("colz text"); c.Print(TString(hEl_final_fakeRate->GetName())+".pdf");
	hMu_final_fakeRate->Draw("colz text"); c.Print(TString(hMu_final_fakeRate->GetName())+".pdf");
	hEl_final_realRate->Draw("colz text"); c.Print(TString(hEl_final_realRate->GetName())+".pdf");
	hMu_final_realRate->Draw("colz text"); c.Print(TString(hMu_final_realRate->GetName())+".pdf");
}

// // */

LEP_SOURCE castSource(MCTC::ParticleType type, MCTC::ParticleOrigin orig, LEP_TYPE l)
{
	/* cf https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/FakeObjectBgEstimation */
	using namespace MCTC;

	if (l==MUON && type==IsoMuon) return REAL;
	if (l==ELEC && type==IsoElectron) return REAL;

	if (orig==CharmedMeson || orig==BottomMeson || orig==CCbarMeson || orig==BBbarMeson 
		|| orig==CharmedBaryon || orig==BottomBaryon || orig==JPsi)
		return HEAVY;

	if (orig==LightMeson || orig==StrangeMeson || orig==LightBaryon || orig==StrangeBaryon 
		|| orig==PionDecay || orig==KaonDecay || orig==PiZero|| orig==NonDefined || orig==DalitzDec)
		return LIGHT;

	if (l==ELEC && orig==Mu && type==NonIsoElectron) return LIGHT;

	if (orig==PhotonConv || orig==BremPhot || orig==PromptPhot || orig==UndrPhot 
		|| orig==ISRPhot || orig==FSRPhot)
		return CONV;

	 return LIGHT;
}

bool finalize()
{
	if (VERBOSE) cout << "Cleaning up..." << endl;

	if (VERBOSE) cout << "Histograms" << endl;
	if (sfhEl_FakeScaleFactors[HEAVY]) sfhEl_FakeScaleFactors[HEAVY]->Draw();
	if (sfhEl_FakeScaleFactors[LIGHT]) sfhEl_FakeScaleFactors[LIGHT]->Draw();
	if (sfhEl_FakeScaleFactors[CONV]) sfhEl_FakeScaleFactors[CONV]->Draw();
	if (sfhMu_FakeScaleFactors[HEAVY]) sfhMu_FakeScaleFactors[HEAVY]->Draw();
	if (sfhMu_FakeScaleFactors[LIGHT]) sfhMu_FakeScaleFactors[LIGHT]->Draw();
	if (hEl_RealScaleFactors) hEl_RealScaleFactors->Draw();
	if (hMu_RealScaleFactors) hMu_RealScaleFactors->Draw();

	for(auto h : allHistos) h->Draw();

	if (VERBOSE) cout << "Drawing lepton composition" << endl;
	lp_El->Normalize(); lp_El->Draw();
	lp_Mu->Normalize(); lp_Mu->Draw();

	if (VERBOSE) cout << "Draw MCTC " << endl;
	TCanvas *can = new TCanvas("cMCTC", "cMCTC", 1280,800); can->SetGrid();
	if (hMu_MCTC) {labelMCTChist(hMu_MCTC); hMu_MCTC->Draw("text"); can->Print("hMu_MCTC.pdf");}
	if (hEl_MCTC) {labelMCTChist(hEl_MCTC); hEl_MCTC->Draw("text");can->Print("hEl_MCTC.pdf");}
	delete can; can=0;

	if (VERBOSE) cout << "Write to file" << endl;
	outFile->Write();
	for (auto h: allHistos) if (h!=0) {delete h; h=0;}

	if (xsecDB!=0){delete xsecDB; xsecDB=0;}
	if (lp_Mu!=0){delete lp_Mu; lp_Mu=0;}
	if (lp_El!=0){delete lp_El; lp_El=0;}

	cout << "All done!!" << endl << endl;
	return true;
}

TChain* loadData(TString fileList, bool isMC){
	//return a TChain linked to the data files
	// bool DEBUG = true;
	TChain* tc = new TChain("evt2l");

	ifstream inF(fileList.Data());
	vector<TString> allFiles; double sumW = 0;
	for (std::string line; getline( inF, line ); ){
		if(DEBUG) cout << line << endl;
		if (line[0]=='#') continue;
	   	allFiles.push_back(line);
	}

 	if (fileList==inFileZjets && isMC) // Weight MC trees
 	{
 		cout << "Weighting MC trees from " << fileList << endl;
 		for(auto it=allFiles.begin(); it!=allFiles.end(); )
 		{
 			std::vector<TString> shortList;

 			// Get sampleID
 			auto s = *it;
 			int sampleID = TString(s(s.Index("TeV")+4,6)).Atoi();
 			if(DEBUG) cout << "sampleID: " << sampleID << '\t';
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
	 		// if (DEBUG) cout << "Getting cross section now..." << endl;
	 		// if (DEBUG && xsecDB) cout << "database is initialized" << endl; 
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
 	if(DEBUG) cout << "Trees loaded " << endl;

	return tc;
}
