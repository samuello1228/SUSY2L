/* ptcorr.C
 * Find pt diff for MC. Reco vs true
 *
 */

#include "../common/evt2l.C"
#include "../common/obj_def.h"
#include <iostream>
#include <iomanip>
#include "TH3D.h"
#include "TLine.h"
#include <fstream>
#include <string>
#include "TLegend.h"
#include "TColor.h"
#include "TSystem.h"
#include "TMath.h"

inline void loadbar(unsigned int x, unsigned int n, unsigned int w = 50){
   if ( (x != n) && (x % (n/100+1) != 0) ) return;

   float ratio    =  x/(float)n;
   unsigned int c =  ratio * w;

   std::cout << std::setw(3) << (int)(ratio*100) << "% [";
   for (unsigned int x=0; x<c; x++) cout << "=";
   for (unsigned int x=c; x<w; x++) cout << " ";
   cout << "]\r" << flush;
}

TChain* loadData(TString fileList){
   //return a TChain linked to the data files
   TChain* tc = new TChain("evt2l");

   if (fileList){
      ifstream inFiles;
      inFiles.open(fileList);
      // cout << "Analysing NTuples:" << endl;
      if (inFiles.is_open()){
         for( std::string line; getline( inFiles, line ); ){
            tc->Add(line.c_str());
            // cout << line << endl;
         }
         inFiles.close();
      } else {
         delete tc;
         tc = 0;
      }
   } else {
      delete tc;
      tc = 0;
   }
   return tc;
}

int getTruthElecI(evt2l& evt2lTree, int i){
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
      // while (evt2lTree.truths_pdgId[tp] == evt2lTree.truths_pdgId[mother]){
      //    tp = mother;
      //    mother = evt2lTree.truths_motherI[tp];
      // }
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

double pt2E(double pt, double eta){
   double ret = TMath::Exp(-eta);
   ret = TMath::Sin(2*ret);
   ret = pt / ret;
   return ret;
}

int EnergyCorr(){
   TChain *inList = loadData("/afs/cern.ch/user/g/ggallard/Zee/common/inFileList-MCdR.txt");
   if(!inList){
      cout << "Input NTuples could not be loaded" << endl;
      return -1;
   } 
   evt2l evt2lTree(inList);

   TFile* outFile = TFile::Open("output.root", "recreate");

	const Double_t PtEdges[]  = {20, 30, 40, 50, 60, 80, 120, 1000.0};
	const Int_t      nPtBins   = 7;
	const Double_t EtaEdges[] = {0, 0.5, 1, 1.37, 1.52, 1.8, 2.0, 2.47}; 
	const Int_t      nEtaBins = 7;

   TH3D* hDE = new TH3D("hDE", "E difference (reco - true);Reconstructed #eta;Reconstructed p_{T} (GeV);#Delta E (GeV)", nEtaBins, EtaEdges, nPtBins, PtEdges);
   TH3D* hDEflipped = (TH3D*) hDE->Clone("hDEflipped");
   TH3D* hDEnoFlip = (TH3D*) hDE->Clone("hDEnoFlip");
   hDE->SetDirectory(outFile);
   hDEflipped->SetDirectory(outFile);

   long long nEntries = evt2lTree.fChain->GetEntries();
   cout << nEntries << " events in TChain" << endl;
   for (long long i =0; i<nEntries; i++){
      loadbar(i+1, nEntries); // loadbar
      evt2lTree.GetEntry(i);

      for(int i=0; i<2; i++){
         if(getTruthElecI(evt2lTree, i)>=0){
            double eta = evt2lTree.leps_eta[i]; double recoPt = evt2lTree.leps_pt[i];
            double recoE = pt2E(recoPt, eta);
            double matchedE = pt2E(evt2lTree.truths_pt[i], eta);
            double dE = recoE - matchedE;
            bool hasFlipped = ((charge*evt2lTree.leps_ID[i]) < 0);

            hDE->Fill(eta, recoPt, dE)
            if(hasFlipped) hDEflipped->Fill(eta, recoPt, dE);
            else hDEnoFlip->Fill(eta, recoPt, dE);
         }
      }
   }
   cout << endl;
   hDPt->Write();

   return 0;
}