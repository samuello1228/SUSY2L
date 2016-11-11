/* ConvertNTuples.cxx
 *
 * Takes text file of UM/HK root files and converts them to input NTuples
 * of the EGamma ChargeMisID tool. 
 *
 * ./ConvertNTuples /inDir/ inFileList.txt /outDir/
 * inFileList.txt has NO PATH. in directory specified in first argument
 * inDir, outDir has to be FULL PATHS
 *
 * Compile with:
 * g++ -O3 -Wall -Wextra -std=c++11 -o ConvertNTuples ConvertNTuples.cxx `root-config --cflags --glibs`
 */

#include <vector>
#include "TSystem.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include "susyEvts.h"
#include "TH1D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TPRegexp.h"

using namespace std;

TString inDir;
TString outDir;
bool isMC = true;
// double PTSCALE = 1.0; // 1000 for GeV to MeV. 

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

double pt2E(double pt, double eta){
  double ret = TMath::Exp(-eta);
  ret = 2*TMath::ATan(ret);
  ret = pt / TMath::Sin(ret);
  return ret;
}

int getOrigElecI(const susyEvts* mEvts, int i){
  stringstream tempOut;
  // Find matching truth particle
  int tp = mEvts->leps[i].truthI;
  if (tp < 0) return -1; // tp = -1 if no truth particle;
  int matchedPDGID = mEvts->truths[tp].pdgId;
  tempOut << "PDGIDs: " << matchedPDGID << '\t';
  if (abs(matchedPDGID) != 11) {
    tempOut << " Not an electron" << endl;
    return -1; // Reject if not matched to an electron
  }

  // Find original electron from Z
  int mother = mEvts->truths[tp].motherI;

  // while mother exists and is not photon or electron
  while(mother>=0 && tp >= 0 && (abs(mEvts->truths[mother].pdgId) == 11 || mEvts->truths[mother].pdgId == 22)  && abs(mEvts->truths[mother].pdgId) < 100){ 
    tp = mother;
    mother = mEvts->truths[tp].motherI;
    tempOut << mEvts->truths[tp].pdgId << '\t';
  }
  tempOut << mEvts->truths[mother].pdgId << '\t';  

  // Reject if the particle found is invalid
  if (tp <0 || mother <0){ 
    tempOut << "Particle invalid" << endl;
    return -1; 
  }

  // Reject if particle found is not electron
  if (abs(mEvts->truths[tp].pdgId)==11){
    // cout << mEvts->truths[tp].pdgId << endl;
    return tp;
  } else if (mEvts->truths[tp].pdgId==22){
    tempOut << "Photon found." << endl;
    cout << flush << '\r' << tempOut.str() << flush;
    return -2;
  }
  else{ 
    tempOut << "Not electron" << endl;
    cout << flush << '\r' << tempOut.str() << flush;
    return -1;
  }
}

// int getOrigElecI(const susyEvts* mEvts, int i){
//   stringstream tempOut;
//   int tp = mEvts->leps[i].truthI;
//   if (tp < 0) return -1;

//   if(mEvts->truths[tp].particleType==2)
//   {
//     // tempOut << " ";
//     if(mEvts->truths[tp].particleOrigin==13){
//       // tempOut << " from Z";
//     } else {
//       tempOut << "Isolated electron not from Z";
//       tp = -1;
//     }
//   } else if(mEvts->truths[tp].particleType==4)
//   {
//     char temp[100];
//     sprintf(temp, "Bkg electron %d/%d/%d\t", mEvts->truths[tp].pdgId, mEvts->truths[tp].particleType, mEvts->truths[tp].particleOrigin);
//     tempOut << temp;
//     int mother = mEvts->truths[tp].motherI;
//     if(mEvts->truths[tp].particleOrigin!=5)
//     {
//       tempOut << "Not from photon conversion";
//       tp = -1;
//     } else 
//     {
//       while(mother>=0 && tp>=0 && abs(mEvts->truths[mother].pdgId) == 11)
//       {
//         sprintf(temp, "%d/%d/%d\t", mEvts->truths[tp].pdgId, mEvts->truths[tp].particleType, mEvts->truths[tp].particleOrigin);
//         tempOut << temp;
//         tp = mother;
//         mother = mEvts->truths[tp].motherI;
//       }
//       sprintf(temp, "%d/%d/%d\t", mEvts->truths[mother].pdgId, mEvts->truths[mother].particleType, mEvts->truths[mother].particleOrigin);
//       tempOut << temp;

//       if(mEvts->truths[tp].particleType==2 && mEvts->truths[tp].particleOrigin==13)
//       {
//         tempOut << "2/13 electron found";
//       } else tp = -1;
//     }
  
//   }
//   if(tempOut.str()!="") cout << tempOut.str() << endl;
//   return tp;
// }

bool convert(TString file){

  TPRegexp reDir("user.clo.v.*myOutput.root/.{0}");
  TPRegexp reFile("user.clo.[0-9]{7}._[0-9]{6}.myOutput.root$");
  // TPRegexp reDir("user.ggallard.v.*myOutput.root/.{0}");
  // TPRegexp reFile("user.ggallard.[0-9]{7}._[0-9]{6}.myOutput.root.?[0-9]?$");
  outDir = file(reDir).Data();
  outDir = outDir(0,outDir.First('/'));
  TString outFilename = file(reFile).Data();
  // outDir = "./";
  // outFilename = "Zee.root";
  cout << "Outdir  = " << outDir << endl;
  cout << "Outfile = " << outFilename << endl;

  // Opening out file
  if(!gSystem->cd(outDir)){
    gSystem->mkdir(outDir);
    if(!gSystem->cd(outDir)){
      cout << "Could not create directory " << outDir << endl;
      return false;
    }
  }
  gSystem->cd("..");

  TFile *outFile = TFile::Open(outDir+"/"+outFilename, "RECREATE");
  if (!outFile){
    cout << "Output file " << file << " could not be created." << endl;
    return false;
  }

  // ===== INITIALIZE OUT TREE ====== //
  float MCEvtWeight;
  float MCPileupWeight; 
  float Zcand_M; 
  unsigned long int trigCode;
  unsigned int cuts;

  // Electron 1
  int elCand1_charge; 
  float elCand1_pt; 
  float elCand1_cl_eta; 
  float elCand1_phi;
  float elCand1_E;
  int elCand1_ID;

  // Electron 2
  int elCand2_charge; 
  float elCand2_pt; 
  float elCand2_cl_eta; 
  float elCand2_phi;
  float elCand2_E;
  int elCand2_ID;

  // Variables for cuts (SUSY)
  int elCand1_flag;
  int elCand2_flag; 

  // Variables for cuts (both)
  float elCand1_d0significance;  // = d0 / sig_d0 
  float elCand2_d0significance; 

  // Truth pt
  float elCand1_truthPt;
  float elCand1_truthE;
  float elCand2_truthPt;
  float elCand2_truthE;

  // Original electron pt
  float elCand1_origPt;
  float elCand1_origE;
  int elCand1_origCharge;
  float elCand2_origPt;
  float elCand2_origE;
  int elCand2_origCharge;

  float elCand1_dRwOrig;
  float elCand2_dRwOrig;

  TTree *outTree = new TTree("ZeeCandidate", "ZeeCandidate");

  string var;
  #define NEWBRANCH(b,v)                  \
    var = #b; var += "/"; var +=  #v;     \
    outTree->Branch(#b, &b, var.c_str());

  NEWBRANCH(MCEvtWeight,F);
  NEWBRANCH(MCPileupWeight, F);
  NEWBRANCH(Zcand_M, F);
  NEWBRANCH(trigCode,l);
  NEWBRANCH(cuts,i);

  NEWBRANCH(elCand1_charge, I);
  NEWBRANCH(elCand1_pt, F);
  NEWBRANCH(elCand1_cl_eta, F);
  NEWBRANCH(elCand1_phi, F);
  NEWBRANCH(elCand1_E, F);
  NEWBRANCH(elCand1_ID,I);

  NEWBRANCH(elCand2_charge, I);
  NEWBRANCH(elCand2_pt, F);
  NEWBRANCH(elCand2_cl_eta, F);
  NEWBRANCH(elCand2_phi, F);
  NEWBRANCH(elCand2_E, F);
  NEWBRANCH(elCand2_ID,I);

  NEWBRANCH(elCand1_flag,I);
  NEWBRANCH(elCand2_flag,I);
  NEWBRANCH(elCand1_d0significance, F);
  NEWBRANCH(elCand2_d0significance, F);

  
  NEWBRANCH(elCand1_truthPt,F);
  NEWBRANCH(elCand2_truthPt,F);
  NEWBRANCH(elCand1_truthE,F);
  NEWBRANCH(elCand2_truthE,F);

  NEWBRANCH(elCand1_origPt,F);
  NEWBRANCH(elCand2_origPt,F);
  NEWBRANCH(elCand1_origE,F);
  NEWBRANCH(elCand2_origE,F);

  NEWBRANCH(elCand1_origCharge,I);
  NEWBRANCH(elCand2_origCharge,I);

  NEWBRANCH(elCand1_dRwOrig, F);
  NEWBRANCH(elCand2_dRwOrig, F);



  #undef NEWBRANCH

  // === INITIALIZE IN TREE ==== //
  TChain *inChain = new TChain("evt2l");
  inChain->Add(file);
  susyEvts* mEvts = new susyEvts(inChain);

  // == FILL TREE
  long int nEntries = inChain->GetEntries();
  for(long int i=0; i<nEntries; i++){
    loadbar(i+1,nEntries);
    mEvts->GetEntry(i);

    // Require exactly two leptons, and the leptons are electrons
    if(mEvts->leps.size()!=2) continue;
    if (int(fabs(mEvts->leps[0].ID/1000))!=11 || int(fabs(mEvts->leps[1].ID/1000))!=11) continue;

    MCEvtWeight = mEvts->evt.weight*mEvts->evt.ElSF*mEvts->evt.MuSF;
    MCPileupWeight = mEvts->evt.pwt;
    trigCode = mEvts->sig.trigCode;
    cuts = mEvts->evt.cuts;

    elCand1_pt = mEvts->leps[0].pt;
    elCand1_cl_eta = mEvts->leps[0].eta;
    elCand1_charge = (mEvts->leps[0].ID > 0) - (mEvts->leps[0].ID < 0); // Charge = sign(ID)
    elCand1_flag = mEvts->leps[0].lFlag;
    elCand1_d0significance = mEvts->leps[0].d0/mEvts->leps[0].d0Err;
    elCand1_phi = mEvts->leps[0].phi;
    elCand1_ID = mEvts->leps[0].ID;
    elCand1_E = pt2E(elCand1_pt, elCand1_cl_eta);

    TLorentzVector p1;
    p1.SetPtEtaPhiM(elCand1_pt, elCand1_cl_eta, elCand1_phi, 0.000511);

    if(isMC && mEvts->leps[0].truthI>=0){
      elCand1_truthPt = mEvts->truths[mEvts->leps[0].truthI].pt;
      elCand1_truthE = pt2E(elCand1_truthPt, elCand1_cl_eta);

      int origElecI = getOrigElecI(mEvts, 0);
      if(origElecI>=0){
        elCand1_origPt = mEvts->truths[origElecI].pt;
        elCand1_origE = pt2E(elCand1_origPt, elCand1_cl_eta);
        elCand1_origCharge = (mEvts->truths[origElecI].pdgId < 0) - (mEvts->truths[origElecI].pdgId > 0);

        TLorentzVector ptruth;
        ptruth.SetPtEtaPhiM(elCand1_origPt, mEvts->truths[origElecI].eta, mEvts->truths[origElecI].phi, 0.000511);
        if(elCand1_origPt==0) cout << mEvts->truths[origElecI].eta << endl;
        elCand1_dRwOrig = p1.DeltaR(ptruth);

      } else {elCand1_origPt = elCand1_origE = -1; elCand1_origCharge = elCand1_dRwOrig = 0; }
    } 
    else {
      elCand1_origPt = elCand1_origE = elCand1_truthE = elCand1_truthPt = -1;
      elCand1_dRwOrig = elCand1_origCharge = 0;
    }

    elCand2_pt = mEvts->leps[1].pt;
    elCand2_cl_eta = mEvts->leps[1].eta;
    elCand2_charge = (mEvts->leps[1].ID > 0) - (mEvts->leps[1].ID < 0); // Charge = sign(ID)
    elCand2_flag = mEvts->leps[1].lFlag;
    elCand2_d0significance = mEvts->leps[1].d0/mEvts->leps[1].d0Err;
    elCand2_phi = mEvts->leps[1].phi;
    elCand2_ID = mEvts->leps[1].ID;
    elCand2_E = pt2E(elCand2_pt, elCand2_cl_eta);

    TLorentzVector p2;
    p2.SetPtEtaPhiM(elCand2_pt, elCand2_cl_eta, elCand2_phi, 0.000511);

    if(isMC && mEvts->leps[1].truthI>=0){
      elCand2_truthPt = mEvts->truths[mEvts->leps[1].truthI].pt;
      elCand2_truthE = pt2E(elCand2_truthPt, elCand2_cl_eta);

      int origElecI = getOrigElecI(mEvts, 1);
      if(origElecI>=0){
        elCand2_origPt = mEvts->truths[origElecI].pt;
        elCand2_origE = pt2E(elCand2_origPt, elCand2_cl_eta);
        elCand2_origCharge = (mEvts->truths[origElecI].pdgId < 0) - (mEvts->truths[origElecI].pdgId > 0);

        TLorentzVector ptruth;
        ptruth.SetPtEtaPhiM(elCand2_origPt, mEvts->truths[origElecI].eta, mEvts->truths[origElecI].phi, 0.000511);
        if(elCand2_origPt==0) cout << mEvts->truths[origElecI].eta << endl;
        elCand2_dRwOrig = p2.DeltaR(ptruth);

      } else {elCand2_origPt = elCand2_origE = -1; elCand2_dRwOrig = elCand2_origCharge = 0;}
    } 
    else {
      elCand2_origPt = elCand2_origE = elCand2_truthE = elCand2_truthPt = -1;
      elCand2_dRwOrig = elCand2_origCharge = 0;
    }  

    Zcand_M = (p1+p2).M();

    outTree->Fill();
  }
  outFile->Write();
  return true;
}

int main(int argc, char *argv[]){
	if (argc!=3){
		cout << "Wrong number of arguments. Two required. " << endl
			<< "./ConvertNTuples inFileList.txt /FullPath/outDir/" << endl;
		return -1;
	}

	ifstream inFiles(argv[1]);

  outDir = argv[2];
	if(!gSystem->cd(outDir)){
		gSystem->mkdir(outDir);
		if(!gSystem->cd(outDir)){
      cout << "Output directory could not be created" << endl;
      return -3;
    }
	}
  
  vector<TString> allFiles;
  string file;
  while(getline(inFiles, file)){
    allFiles.push_back(file);
  }
  int nFiles = allFiles.size();
  int n = 1;
	for (TString file : allFiles){
    cout << "File " << n << " of " << nFiles << ": " << file << endl; 
    // convert1(file);
		if (!convert(file)) return -1;
    cout << endl << "-" << endl;
    n++;
    // if(n>1) break;
	}

}
