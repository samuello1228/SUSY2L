#include "support/evt2l.C"
#include "support/obj_def.h"
#include <TH1D.h>
#include <TString.h>
#include <iostream>
#include <fstream>
#include <TFile.h>
#include <vector>
#include <TChain.h>
#include <TSystem.h>

using std::cout; using std::endl;

// ============ CONFIG VARIABLES ================= //
const bool DEBUG = true;

const TString inFileData="../support/inFileList-dataTest.txt";
// const TString inFileData="../support/inFileList-data.txt";

// TString inFileBBCC="../support/inFileList-.txt";
const TString inFileTTbar="../support/inFileList-ttbar.txt";
const TString inFileZjets="../support/inFileList-Zjets.txt";
const TString inFileWjets="../support/inFileList-Wjets.txt";
const TString inFileVV="../support/inFileList-VV.txt";
const TString inFileVgamma="../support/inFileList-Vgamma.txt";
const TString inFileSingleTop="../support/inFileList-SingleTop.txt";
const TString inFileMC="../support/inFileList-allMC.txt";

const TString outDir="validation/";

// ============ INFRASTRUCTURE ================= //
TFile* outFile;
enum DATA_TYPE{
	DATA,
	MC,
};
const int N_DATA_TYPE=2;

TString DataTypeToString(DATA_TYPE d)
{
	switch(d){
		case (DATA): return "DATA";
		case (MC): return "MC";
		default: return "";
	};
}
enum HIST_TYPE
{
	mll,
	LeadingPt,
	SubleadingPt,
	LeadingEta,
	SubleadingEta,
	// meff,
	mt2,
	met,
}; const int N_HIST_TYPE=7;


enum LEP_TYPE {ELEC, MUON};

TString HistTypeToString(HIST_TYPE h)
{
	switch(h){
		case (mll): return "mll";
		case (LeadingPt): return "LeadingPt";
		case (SubleadingPt): return "SubleadingPt";
		case (LeadingEta): return "LeadingEta";
		case (SubleadingEta): return "SubleadingEta";
		// case (meff): return "meff";
		case (mt2): return "mt2";
		case (met): return "met";
		default: return "";
	};
}
struct histType{
	TString name;
	TString title;
	TString xLabel;
	TString unit;

	int nBins;
	double xLow;
	double xHigh;
};

std::vector<histType> allHistTypes;

class SetOfHist{
public:
	TString regionName;
	typedef std::vector<TH1*> dataRow;
	std::vector<dataRow> hist;

	SetOfHist(): SetOfHist("test"){};
	SetOfHist(TString regionName)
	{
		this->regionName = regionName;
		hist.reserve(2);
		for(int i=0; i<2; i++){
			dataRow v; v.reserve(allHistTypes.size());
			for(uint j=0; j<allHistTypes.size(); j++)
			{
				TString hName = regionName+"_"; hName += allHistTypes[j].name; hName += i==0? "_DATA": "_MC";
				TString hTitle = allHistTypes[j].title + " in "; hTitle += regionName;
						hTitle += ";"+ allHistTypes[j].xLabel; hTitle += " " + allHistTypes[j].unit;
						hTitle += ";"; hTitle += "Events/"; hTitle += (allHistTypes[j].xHigh-allHistTypes[j].xLow)/allHistTypes[j].nBins; hTitle += allHistTypes[j].unit;
				TH1* h = new TH1D(hName, hTitle, allHistTypes[j].nBins, allHistTypes[j].xLow, allHistTypes[j].xHigh);
				h->SetDirectory(outFile);
				v.push_back(h);
			}
			hist.push_back(v);
		}
	}
	void Fill(evt2l* tree, DATA_TYPE d)
	{
		double w = tree->evt_weight*tree->evt_pwt*tree->fChain->GetTree()->GetWeight()*tree->evt_ElSF*tree->evt_MuSF;
		hist[d][mll]->Fill(tree->l12_m, w);
		hist[d][LeadingPt]->Fill(tree->leps_pt[0], w);
		hist[d][SubleadingPt]->Fill(tree->leps_pt[1], w);
		hist[d][LeadingEta]->Fill(tree->leps_eta[0], w);
		hist[d][SubleadingEta]->Fill(tree->leps_eta[1], w);
		hist[d][mt2]->Fill(tree->sig_mT2, w);
		hist[d][met]->Fill(tree->sig_Met, w);
	}

	void Draw()
	{
		TCanvas *c = new TCanvas();
		for(int j=0; j<N_HIST_TYPE; j++)
		{
			hist[MC][j]->SetFillColor(kAzure+4);
			hist[MC][j]->SetLineColor(kAzure+4);
			hist[MC][j]->SetMarkerColor(kAzure+4);

			hist[DATA][j]->SetLineColor(kBlack); 
			hist[DATA][j]->SetMarkerColor(kBlack); 
			hist[DATA][j]->SetMarkerStyle(20);

			hist[MC][j]->Draw("bar e");
			hist[DATA][j]->Draw("same");

			TString printName = TString(hist[MC][j]->GetName())(0,TString(hist[MC][j]->GetName()).Length()-3);
			printName+=".pdf";
			c->Print(printName);

		}

	}

};

SetOfHist *hHeavy;
SetOfHist *hLight;
SetOfHist *hConv;
SetOfHist *hReal;
SetOfHist *hAll;


// ============ FUNCTION DELCARATIONS ================= //
void initialize();
void makeValidationPlots();
void loop(evt2l* tree, DATA_TYPE d);
TChain* loadData(TString fileList, bool isMC);
int castLepPdgId(int lepID);
TChain* loadData(TString, DATA_TYPE);

// ============ FUNCTION DELCARATIONS ================= //
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

void makeValidationPlots()
{
	initialize();

	TChain* tMC = new TChain("evt2l");
	tMC->Add(loadData(inFileTTbar, MC));
	tMC->Add(loadData(inFileZjets, MC));
	tMC->Add(loadData(inFileWjets, MC));
	tMC->Add(loadData(inFileVV, MC));
	tMC->Add(loadData(inFileVgamma, MC));
	tMC->Add(loadData(inFileSingleTop, MC));
	evt2l* mcTree = new evt2l(tMC);
	loop(mcTree, MC);

	evt2l* dataTree = new evt2l(loadData(inFileData, DATA));
	loop(dataTree, DATA);

	hAll->Draw();
	outFile->Write();
	outFile->Close();
}

void loop(evt2l* tree, DATA_TYPE d)
{
	cout << "Beginning loop of type " << DataTypeToString(d) << endl;
	// if(!tree) {cout << "Tree is empty!!" << endl; return;}
	// else cout << "Tree is not empty" << endl;

	// if(!(tree->fChain)) {cout << "no chain!!" << endl; return;}
	// else cout << "Chain is there" << endl;

	int nEntries = tree->fChain->GetEntries();
	for(int i=0; i<nEntries; i++)
	{
		loadbar(i+1, nEntries);
		tree->GetEntry(i);
		bool pass = true;
		pass *= tree->leps_==2;
		pass *= (tree->leps_lFlag[0] & IS_SIGNAL) && (tree->leps_lFlag[1] & IS_SIGNAL);
		pass *= tree->leps_ID[0]*tree->leps_ID[1]>0;
		pass *= tree->jets_<=3;

		hAll->Fill(tree, d);
	}

}

int castLepPdgId(int lepID){return int(lepID/1000);}

void initialize()
{
	if(!gSystem->cd(outDir)){
		gSystem->mkdir(outDir);
		gSystem->cd(outDir);
	}

	outFile = new TFile("validation.root", "recreate");

	allHistTypes.push_back({"mll", "m_{ll}", "m_{ll}", "[GeV]", 10, 20, 220});
	allHistTypes.push_back({"LeadingPt", "leading p_{T}", "p_{T}^{l1}", "[GeV]", 10, 20, 220});
	allHistTypes.push_back({"SubleadingPt", "subleading p_{T}", "p_{T}^{l2}", "[GeV]", 10, 20, 220});
	allHistTypes.push_back({"LeadingEta", "leading |#eta|", "|#eta^{l1}|", "", 4, 0, 2.4});
	allHistTypes.push_back({"SubleadingEta", "subleading #eta", "|#eta^{l2}|", "", 4, 0, 2.4});
	allHistTypes.push_back({"mt2", "m_{T2}", "m_{T2}", "[GeV]", 10, 20, 220});
	allHistTypes.push_back({"met", "MET", "MET", "[GeV]", 10, 20, 220});

	// hHeavy = new SetOfHist("hHeavy");
	// hLight = new SetOfHist("hLight");
	// hConv = new SetOfHist("hConv");
	// hReal = new SetOfHist("hReal");
	hAll = new SetOfHist("hAll");
}

TChain* loadData(TString fileList, DATA_TYPE d){
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

 	cout << "Loading trees from " << fileList << endl;
 	for(auto f : allFiles)
 	{
 		if(d==MC) {
 			TString fname= TString(f(0,f.Length()-5));
 			fname+="_WEIGHTED.root";
 			tc->Add(fname);
 		} else tc->Add(f);
 	}

	return tc;
}