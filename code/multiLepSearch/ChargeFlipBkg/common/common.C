/* Collection of functions common to these scripts
 * Not designed to run in isolation. Include in header section.
 * Call function by COMMON::fcn()
 * For root macros only, not to be compiled
 */

#ifndef COMMON_C
#define COMMON_C

#include "evt2l.C"
#include "obj_def.h"
#include "TLorentzVector.h"
#include <iostream>
#include <fstream>

namespace COMMON{

	struct ev{
      float MCEvtWeight;
      float MCPileupWeight;

		float elCand1_cl_eta;
		float elCand1_pt;
		int elCand1_charge;
		float elCand1_E;

		float elCand1_origE;
		float elCand1_origPt;
		int elCand1_origCharge;

		float elCand2_cl_eta;
		float elCand2_pt;
		int elCand2_charge;
		float elCand2_E;

		float elCand2_origE;
		float elCand2_origPt;
		int elCand2_origCharge;

		int elCand1_flag;
		int elCand2_flag;

		float elCand1_dRwOrig;
		float elCand2_dRwOrig;
	};

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

	//#define CONNECT(b) event.b = 0; ch->SetBranchStatus(#b,1); ch->SetBranchAddress(#b,&(event.b));

	inline void loadbar(unsigned int x, unsigned int n, unsigned int w);
	int getTruthElecI(evt2l& evt2lTree, int i);
	TChain* loadData(TString fileList);
	double getMll(const evt2l& evt2lTree, int i, int j);
};


inline void COMMON::loadbar(unsigned int x, unsigned int n, unsigned int w = 50)
{
   if ( (x != n) && (x % (n/100+1) != 0) ) return;

   float ratio    =  x/(float)n;
   unsigned int c =  ratio * w;

   std::cout << setw(3) << (int)(ratio*100) << "% [";
   for (unsigned int x=0; x<c; x++) std::cout << "=";
   for (unsigned int x=c; x<w; x++) std::cout << " ";
   std::cout << "]\r" << flush;
}

int COMMON::getTruthElecI(evt2l& evt2lTree, int i){
   // Find matching truth particle
   int tp = evt2lTree.leps_truthI[i];
   if (tp < 0) return -1; // tp = -1 if no truth particle;
   int matchedPDGID = evt2lTree.truths_pdgId[tp];
   if (abs(matchedPDGID) != 11) return -1;

   // Find original electron from Z
   int mother = evt2lTree.truths_motherI[tp];

   while(mother>=0 && tp >= 0 && evt2lTree.truths_pdgId[mother] != 23 && abs(evt2lTree.truths_pdgId[mother]) < 100){ // while mother exists and is not Z
      tp = mother;
      mother = evt2lTree.truths_motherI[tp];
   }

   // Reject if the particle found is not from Z
   if (tp <0 || mother <0 || evt2lTree.truths_pdgId[mother]!=23){ 
      return -1;
   }

   // Reject if particle found is not electron
   if (abs(evt2lTree.truths_pdgId[tp])==11) return tp;
   else {
      return -1;
   }
}

TChain* COMMON::loadData(TString fileList){
   //return a TChain linked to the data files
   TChain* tc = new TChain("evt2l");

   if (fileList){
      std::ifstream inFiles(fileList);
      if (inFiles.is_open()){
         for( std::string line; getline( inFiles, line ); ){tc->Add(line.c_str());}
         inFiles.close();
      }
   }
   return tc;
}

double COMMON::getMll(const evt2l& evt2lTree, int i=0, int j=1){
	TLorentzVector p1, p2;
	p1.SetPtEtaPhiM(evt2lTree.leps_pt[i], evt2lTree.leps_eta[i], evt2lTree.leps_phi[0], 0.000511);
	p2.SetPtEtaPhiM(evt2lTree.leps_pt[j], evt2lTree.leps_eta[j], evt2lTree.leps_phi[1], 0.000511);

	return (p1+p2).M();
}

#endif