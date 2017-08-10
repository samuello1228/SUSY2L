#include "../support/evt2l.C"
#include "../support/obj_def.h"
#include "../support/MCTruthClassifierDefs.h"
#include <TString.h>
#include <TChain.h>
#include <fstream>
#include <iostream>

using namespace MCTC;
using std::cout; 
using std::endl;

// const TString inFileZjets="../support/inFileList-Zjets.txt";
const TString inFileZjets="../support/inFileList-ZPowheg.txt";
int nEvents = -1;


enum LEP_SOURCE {
	REAL=-1,
	HEAVY=0,
	LIGHT=1,
	CONV=2,
	// OTHER=3, // To check lepton composition
};
enum LEP_TYPE {
	ELEC,
	MUON,
};

TChain* loadData(TString);
void printDecay();
LEP_SOURCE castSource(MCTC::ParticleType type, MCTC::ParticleOrigin orig, LEP_TYPE l);
bool passLepCompCR(evt2l* tree);
const TString LS_TOSTRING(LEP_SOURCE s);


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

void printDecay()
{
	using namespace MCTC;
	ofstream outTxt("decays.txt");
	ofstream outLightText("decays-Light.txt");

	// lightOriginEl = TH1D("lightOriginEl", "Origin of light fakes (El)", 2, 0, 2);
	// lightOriginEl.SetBinLabel(1, "UnknownEl+NonDefined origin");
	// lightOriginEl.SetBinLabel(2, "Other");

	// lightOriginMu = TH1D("lightOriginMu", "Origin of light fakes (Mu)", 2, 0, 2);
	// lightOriginMu.SetBinLabel(1, "UnknownMu+NonDefined origin");
	// lightOriginMu.SetBinLabel(2, "Other");

	evt2l* zTree = new evt2l(loadData(inFileZjets));

	int nTotalEvents = zTree->fChain->GetEntries();
	if (nEvents<0 || nEvents>nTotalEvents) nEvents = nTotalEvents;

	outTxt << "There are " << nTotalEvents << " events in " << inFileZjets
			<< ", we will look at " << nEvents << endl;

	int nReal, nHeavy, nLight, nConv;
	nReal = nHeavy = nLight = nConv = 0;

	for(int i=0; i<nEvents; i++)
	{
		zTree->GetEntry(i);
		loadbar(i+1, nEvents);
		if(!passLepCompCR(zTree)) continue;
		outTxt << "Event #" << i << " ================ "<< endl;
		double w = zTree->evt_weight*zTree->evt_pwt*zTree->fChain->GetTree()->GetWeight()*zTree->evt_ElSF*zTree->evt_MuSF;
		for(int j=0; j<2; j++)
		{
			outTxt << "Lepton " << j << ", ";
			if(zTree->leps_pt[j]<20) {outTxt << "Pt <20" << endl; continue;}

			bool lepIsTight = zTree->leps_lFlag[j] & IS_SIGNAL;
			outTxt << (lepIsTight? "Tight":"Loose");
			LEP_TYPE recoLepType;
			if (TMath::Abs(int(zTree->leps_ID[j]/1000)) == 11) {recoLepType=ELEC; outTxt << " Reco el, ";}
			else if (TMath::Abs(int(zTree->leps_ID[j]/1000)) == 13) {recoLepType=MUON; outTxt << " Reco mu, ";}
			else continue;

			ParticleType   type = static_cast<ParticleType>(zTree->leps_truthType[j]);
			ParticleOrigin orig = static_cast<ParticleOrigin>(zTree->leps_truthOrig[j]);
			LEP_SOURCE 	   source = castSource(type, orig, recoLepType);

			outTxt << "Type " << (int) type << ", Origin " << (int) orig << ", ";
			outTxt << LS_TOSTRING(source) << endl << "  ==> ";

			switch(source){
				case (REAL):   nReal++; break;
				case (LIGHT): nLight++;
							outLightText << "Reco " << (recoLepType==ELEC? "el":"mu") << " from event #" << i
										<< "\t| Type " << (int) type << ", Origin " << (int) orig << endl;
								break;
				case (HEAVY): nHeavy++; break;
				case (CONV):   nConv++; break;
			}


			int truthI = zTree->leps_truthI[j];
			while(truthI >=0)
			{
				outTxt << zTree->truths_pdgId[truthI] << " ";
				if (source==LIGHT) outLightText << zTree->truths_pdgId[truthI] << " ";
				truthI = zTree->truths_motherI[truthI];
			}
			outTxt << endl;
			if(source==LIGHT) outLightText << endl << endl;


		}
		outTxt << endl << endl;
	}
	outTxt << "Total number of electrons: " << endl
			<< "REAL  : " << nReal << endl
			<< "LIGHT : " << nLight << endl
			<< "HEAVY : " << nHeavy << endl
			<< "CONV  : " << nConv  << endl;
	cout << endl << endl;
}

bool passLepCompCR(evt2l* tree)
{
	/* LEPTON COMPOSOTION CONTROL REGION
	- Exactly two leptons
		- Same sign
		- Signal
	- No more than three jets
	- Z mass veto

	Recall f= N(tight)/N(loose), should require signal leptons
	*/
	bool pass = true;
	pass *= tree->leps_==2;
	pass *= (tree->leps_lFlag[0] & IS_SIGNAL) && (tree->leps_lFlag[1] & IS_SIGNAL);
	pass *= tree->leps_ID[0]*tree->leps_ID[1]>0;
	pass *= tree->jets_<=3;
	pass *= fabs(tree->l12_m-91.2)>10;  

	// // No b-jets
	// for(int i(0); i>tree->jets_;i++) pass*= !(tree->jets_jFlag[i] & JT_BJET);

	return pass;
}

const TString LS_TOSTRING(LEP_SOURCE s)
{
	if(s==REAL) return "Real";
	if(s==HEAVY) return "Heavy";
	if(s==LIGHT) return "Light";
	if(s==CONV) return "Conv";
	// if(s==OTHER) return "Other";
	return "Unknown";
}

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


TChain* loadData(TString fileList){
	//return a TChain linked to the data files
	// bool DEBUG = true;
	TChain* tc = new TChain("evt2l");

	ifstream inF(fileList.Data());
	vector<TString> allFiles; double sumW = 0;
	for (std::string line; getline( inF, line ); ){
		if (line[0]=='#') continue;
	   	allFiles.push_back(line);
	}

 	cout << "Loading trees from " << fileList << endl;
 	for(auto f : allFiles)
 	{
		TString fname= TString(f(0,f.Length()-5));
		fname+="_WEIGHTED.root";
		tc->Add(fname);
 	}
	return tc;
}