//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jun 23 16:57:17 2015 by ROOT version 6.02/05
// from TTree data/Electron event data
// found on file: hist-sample.root
//////////////////////////////////////////////////////////


// TODO

#ifndef TagProbePlot_h
#define TagProbePlot_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TH1.h>
//#include <TProfile.h>
#include <iostream>
#include <TH2.h>
#include <string>
#include <iomanip>
#include <sstream>
#include <fstream>

// Header file for the classes stored in the TTree if any.
#include "TObject.h"

#define Z_MASS_UPPER_LIMIT 100
#define Z_MASS_LOWER_LIMIT 80

#define PRINT_STAT(MSG, VALUE)                                       \
{                                                                    \
   std::cout << MSG << ": " << VALUE << std::endl;                   \
}                                                                    \


class TagProbePlot : public TSelector {
private:
// "enum" for accessing elements in histogram arrays
   static const int ERROR = 0;
   static const int N_PROBES = 1;
   static const int N_SSPROBES = 2;
   static const int N_OSPROBES = 3;
   
   static const int DEFAULT = 0;
   static const int LOOSE = 1;
   static const int MEDIUM = 2;
   static const int TIGHT = 3;
   static const int MC = 10;
   
   static const int AUTO = 0;
   static const int MAN = 1;
   
   static const bool ABS_ETA = true;
   
   //TFile *outputFile = 0;

public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   
   // Bin edges. Defined in accompanying C file. Remember to update nBins accordingly
   static const Double_t   *PtEdges;
   static const Int_t      nPtBins = 9;

   static const Double_t   *EtaEdges; 
   static const Int_t      nEtaBins = 7;
   
   static const Bool_t onlyZMCevents = false;
   
   static const Bool_t  print = false; // Set to true to output to pdf
   static const Bool_t  toFile = true; // Set to true to output to output.root
   static const Bool_t  drawEta = true;
   static const Bool_t  drawPt  = true;
   static const Bool_t  drawEtaPt = false; // Draw eta dist in diff pt bins
   static const Bool_t  drawPtEta = false; // Draw pt dist in diff eta bins
   
   
   // string submitDir;
   
   // For error(Pt)
   TH1            *hPt[4][4];
   TH1            *hPt_MC[4][4];
   
   //For error(eta)
   TH1            *hEta[4][4];
   TH1            *hEta_MC[4][4];
   
   // For error(InvMass)
   TH1            *hMass[4][4];
   TH1            *hMass_MC[4][4];
   
   // For error(eta, pt)
   TH2            *hEP[4][4];
   TH2            *hEP_MC[4][4];
   
   // Distributions of basic properties of ee pair
   TH1            *hMassDist[4][4];
   TH1            *hEtaDist[4][4];
   TH1            *hPtDist[4][4];
   TH1            *hMassDist_MC[4][4];
   TH1            *hEtaDist_MC[4][4];
   TH1            *hPtDist_MC[4][4];
   
   TH2            *hEtaMass[4][4];
   TH2            *hPtMass[4][4];
   TH2            *hEtaMassNoZ[4][4];
   TH2            *hPtMassNoZ[4][4];

   
   // To check for particle origin 
   TH1            *hParticleOrigin = 0;
   TH1            *hParticleOriginSS = 0;
   TH2            *hPOvsD0 = 0;
   TH1            *hParents = 0;
   TH1            *hPairOrigin = 0;
   
   // Statistics
   Long64_t          nRealE = 0;
   Long64_t          nTotalE = 0;
   Long64_t          nZ = 0;
   Long64_t          nRealZ = 0;
   Long64_t          nMCevents = 0;
   Long64_t          nTotEvents = 0;
   Long64_t          nEvents = 0;
   Long64_t          nPassedEvents = 0;
   Long64_t          nSSEvents = 0;
   Long64_t          nSSEvents_MC = 0;
   Long64_t          nNonZProbes = 0;
   Long64_t          nZMCandSS = 0;
   Long64_t          nTpPairs = 0;

   Long64_t          nFromZ = 0;
   
   TH1               *hPtFromPhotonOnly = 0;
   TH1               *hPtFromPEZ = 0;
   TH1               *hTagCuts = 0;
   
   // output 
   std::ofstream     file_eParents;
   std::ofstream     file_nProbes;
   std::ofstream     file_misID;
   TFile             *outputFile = 0;

   //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//
   // AUTO ADDED AFTER THIS                                            //
   //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   UInt_t          fUniqueID;
   UInt_t          fBits;
   UInt_t          fE1_fUniqueID;
   UInt_t          fE1_fBits;
   Double_t        fE1_fCharge;
   Double_t        fE1_fEta;
   Double_t        fE1_fPt;
   Double_t        fE1_fD0;
   Double_t        fE1_fD0Error;
   Int_t           fE1_fSiHits;
   Int_t           fE1_fBHits;
   Bool_t          fE1_fTriggerMatched;
   Bool_t          fE1_fTag;
   Int_t           fE1_fPassCut;
   UInt_t          fE2_fUniqueID;
   UInt_t          fE2_fBits;
   Double_t        fE2_fCharge;
   Double_t        fE2_fEta;
   Double_t        fE2_fPt;
   Double_t        fE2_fD0;
   Double_t        fE2_fD0Error;
   Int_t           fE2_fSiHits;
   Int_t           fE2_fBHits;
   Bool_t          fE2_fTriggerMatched;
   Bool_t          fE2_fTag;
   Int_t           fE2_fPassCut;
   Double_t        fInvMass;
   Double_t        fdPhi;
   Bool_t          fSSevent;
   Bool_t          fMCevent;
   UInt_t          ftE1_fUniqueID;
   UInt_t          ftE1_fBits;
   Double_t        ftE1_fCharge;
   Double_t        ftE1_fEta;
   Double_t        ftE1_fPt;
   Int_t           ftE1_fParticleType;
   Int_t           ftE1_fParticleOrigin;
   Int_t           ftE1_fTruthType;
   Int_t           ftE1_fTruthOrigin;
   Int_t           ftE1_fTruthBkgType;
   Int_t           ftE1_fTruthBkgOrigin;
   Double_t        ftE1_fTruthBkgCharge;
   Bool_t          ftE1_fIsElectron;
   Bool_t          ftE1_fFromZ;
   Double_t        ftE1_fDeltaR;
   UInt_t          ftE2_fUniqueID;
   UInt_t          ftE2_fBits;
   Double_t        ftE2_fCharge;
   Double_t        ftE2_fEta;
   Double_t        ftE2_fPt;
   Int_t           ftE2_fParticleType;
   Int_t           ftE2_fParticleOrigin;
   Int_t           ftE2_fTruthType;
   Int_t           ftE2_fTruthOrigin;
   Int_t           ftE2_fTruthBkgType;
   Int_t           ftE2_fTruthBkgOrigin;
   Double_t        ftE2_fTruthBkgCharge;
   Bool_t          ftE2_fIsElectron;
   Bool_t          ftE2_fFromZ;
   Double_t        ftE2_fDeltaR;
   Double_t        ftInvMass;
   Double_t        ftdPhi;
   Bool_t          ftSSevent;
   Bool_t          ftZEEevent;

   // List of branches
   TBranch        *b_Event_fUniqueID;   //!
   TBranch        *b_Event_fBits;   //!
   TBranch        *b_Event_fE1_fUniqueID;   //!
   TBranch        *b_Event_fE1_fBits;   //!
   TBranch        *b_Event_fE1_fCharge;   //!
   TBranch        *b_Event_fE1_fEta;   //!
   TBranch        *b_Event_fE1_fPt;   //!
   TBranch        *b_Event_fE1_fD0;   //!
   TBranch        *b_Event_fE1_fD0Error;   //!
   TBranch        *b_Event_fE1_fSiHits;   //!
   TBranch        *b_Event_fE1_fBHits;   //!
   TBranch        *b_Event_fE1_fTriggerMatched;   //!
   TBranch        *b_Event_fE1_fTag;   //!
   TBranch        *b_Event_fE1_fPassCut;   //!
   TBranch        *b_Event_fE2_fUniqueID;   //!
   TBranch        *b_Event_fE2_fBits;   //!
   TBranch        *b_Event_fE2_fCharge;   //!
   TBranch        *b_Event_fE2_fEta;   //!
   TBranch        *b_Event_fE2_fPt;   //!
   TBranch        *b_Event_fE2_fD0;   //!
   TBranch        *b_Event_fE2_fD0Error;   //!
   TBranch        *b_Event_fE2_fSiHits;   //!
   TBranch        *b_Event_fE2_fBHits;   //!
   TBranch        *b_Event_fE2_fTriggerMatched;   //!
   TBranch        *b_Event_fE2_fTag;   //!
   TBranch        *b_Event_fE2_fPassCut;   //!
   TBranch        *b_Event_fInvMass;   //!
   TBranch        *b_Event_fdPhi;   //!
   TBranch        *b_Event_fSSevent;   //!
   TBranch        *b_Event_fMCevent;   //!
   TBranch        *b_Event_ftE1_fUniqueID;   //!
   TBranch        *b_Event_ftE1_fBits;   //!
   TBranch        *b_Event_ftE1_fCharge;   //!
   TBranch        *b_Event_ftE1_fEta;   //!
   TBranch        *b_Event_ftE1_fPt;   //!
   TBranch        *b_Event_ftE1_fParticleType;   //!
   TBranch        *b_Event_ftE1_fParticleOrigin;   //!
   TBranch        *b_Event_ftE1_fTruthType;   //!
   TBranch        *b_Event_ftE1_fTruthOrigin;   //!
   TBranch        *b_Event_ftE1_fTruthBkgType;   //!
   TBranch        *b_Event_ftE1_fTruthBkgOrigin;   //!
   TBranch        *b_Event_ftE1_fTruthBkgCharge;   //!
   TBranch        *b_Event_ftE1_fIsElectron;   //!
   TBranch        *b_Event_ftE1_fFromZ;   //!
   TBranch        *b_Event_ftE1_fDeltaR;   //!
   TBranch        *b_Event_ftE2_fUniqueID;   //!
   TBranch        *b_Event_ftE2_fBits;   //!
   TBranch        *b_Event_ftE2_fCharge;   //!
   TBranch        *b_Event_ftE2_fEta;   //!
   TBranch        *b_Event_ftE2_fPt;   //!
   TBranch        *b_Event_ftE2_fParticleType;   //!
   TBranch        *b_Event_ftE2_fParticleOrigin;   //!
   TBranch        *b_Event_ftE2_fTruthType;   //!
   TBranch        *b_Event_ftE2_fTruthOrigin;   //!
   TBranch        *b_Event_ftE2_fTruthBkgType;   //!
   TBranch        *b_Event_ftE2_fTruthBkgOrigin;   //!
   TBranch        *b_Event_ftE2_fTruthBkgCharge;   //!
   TBranch        *b_Event_ftE2_fIsElectron;   //!
   TBranch        *b_Event_ftE2_fFromZ;   //!
   TBranch        *b_Event_ftE2_fDeltaR;   //!
   TBranch        *b_Event_ftInvMass;   //!
   TBranch        *b_Event_ftdPhi;   //!
   TBranch        *b_Event_ftSSevent;   //!
   TBranch        *b_Event_ftZEEevent;   //!

   TagProbePlot(TTree * /*tree*/ =0) : fChain(0) { }
   virtual ~TagProbePlot() { }
   virtual Int_t   Version() const { return 2; }
   virtual void    Begin(TTree *tree);
   virtual void    SlaveBegin(TTree *tree);
   virtual void    Init(TTree *tree);
   virtual Bool_t  Notify();
   virtual Bool_t  Process(Long64_t entry);
   virtual Int_t   GetEntry(Long64_t entry, Int_t getall = 0) { return fChain ? fChain->GetTree()->GetEntry(entry, getall) : 0; }
   virtual void    SetOption(const char *option) { fOption = option; }
   virtual void    SetObject(TObject *obj) { fObject = obj; }
   virtual void    SetInputList(TList *input) { fInput = input; }
   virtual TList  *GetOutputList() const { return fOutput; }
   virtual void    SlaveTerminate();
   virtual void    Terminate();
   
   ClassDef(TagProbePlot,0);

   //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//
   // USER ADDED AFTER THIS                                            //
   //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//
   
   inline bool MisId(int code);
   inline void ProcessData();
   void ProcessMC();   
   void FillHist(Double_t var, TH1* h[], int code, bool absolute = false);
   void FillEPHist(Double_t eta, Double_t pt, TH2* h[], int code, bool absEta = false);
   //void FillEtaPtHist(Double_t eta, Double_t pt, TH1* h[][3], int code, bool absEta = false);
   void FillArrays();
   void DivideHist(TH1* h[]);
   inline Bool_t PassedEventCut();
   inline Bool_t EventIsZMC();
   bool ElectronIsTag(int elec, bool isData=true, int code=MAN);
   void PrintParticleOrigin();
   void ToCSV();
   void ToCSV(int);
   void ToFile();
   
   // elec = electron number
   // cut: 1 = light cut; 2 = medium cut; 3 = tight cut;
   Bool_t PassedLHCut(int elec, int cut);

};

#endif

#ifdef TagProbePlot_cxx
void TagProbePlot::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fChain->SetMakeClass(1);
   
   nTotEvents = fChain->GetEntries();

    fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_Event_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_Event_fBits);
   fChain->SetBranchAddress("fE1.fUniqueID", &fE1_fUniqueID, &b_Event_fE1_fUniqueID);
   fChain->SetBranchAddress("fE1.fBits", &fE1_fBits, &b_Event_fE1_fBits);
   fChain->SetBranchAddress("fE1.fCharge", &fE1_fCharge, &b_Event_fE1_fCharge);
   fChain->SetBranchAddress("fE1.fEta", &fE1_fEta, &b_Event_fE1_fEta);
   fChain->SetBranchAddress("fE1.fPt", &fE1_fPt, &b_Event_fE1_fPt);
   fChain->SetBranchAddress("fE1.fD0", &fE1_fD0, &b_Event_fE1_fD0);
   fChain->SetBranchAddress("fE1.fD0Error", &fE1_fD0Error, &b_Event_fE1_fD0Error);
   fChain->SetBranchAddress("fE1.fSiHits", &fE1_fSiHits, &b_Event_fE1_fSiHits);
   fChain->SetBranchAddress("fE1.fBHits", &fE1_fBHits, &b_Event_fE1_fBHits);
   fChain->SetBranchAddress("fE1.fTriggerMatched", &fE1_fTriggerMatched, &b_Event_fE1_fTriggerMatched);
   fChain->SetBranchAddress("fE1.fTag", &fE1_fTag, &b_Event_fE1_fTag);
   fChain->SetBranchAddress("fE1.fPassCut", &fE1_fPassCut, &b_Event_fE1_fPassCut);
   fChain->SetBranchAddress("fE2.fUniqueID", &fE2_fUniqueID, &b_Event_fE2_fUniqueID);
   fChain->SetBranchAddress("fE2.fBits", &fE2_fBits, &b_Event_fE2_fBits);
   fChain->SetBranchAddress("fE2.fCharge", &fE2_fCharge, &b_Event_fE2_fCharge);
   fChain->SetBranchAddress("fE2.fEta", &fE2_fEta, &b_Event_fE2_fEta);
   fChain->SetBranchAddress("fE2.fPt", &fE2_fPt, &b_Event_fE2_fPt);
   fChain->SetBranchAddress("fE2.fD0", &fE2_fD0, &b_Event_fE2_fD0);
   fChain->SetBranchAddress("fE2.fD0Error", &fE2_fD0Error, &b_Event_fE2_fD0Error);
   fChain->SetBranchAddress("fE2.fSiHits", &fE2_fSiHits, &b_Event_fE2_fSiHits);
   fChain->SetBranchAddress("fE2.fBHits", &fE2_fBHits, &b_Event_fE2_fBHits);
   fChain->SetBranchAddress("fE2.fTriggerMatched", &fE2_fTriggerMatched, &b_Event_fE2_fTriggerMatched);
   fChain->SetBranchAddress("fE2.fTag", &fE2_fTag, &b_Event_fE2_fTag);
   fChain->SetBranchAddress("fE2.fPassCut", &fE2_fPassCut, &b_Event_fE2_fPassCut);
   fChain->SetBranchAddress("fInvMass", &fInvMass, &b_Event_fInvMass);
   fChain->SetBranchAddress("fdPhi", &fdPhi, &b_Event_fdPhi);
   fChain->SetBranchAddress("fSSevent", &fSSevent, &b_Event_fSSevent);
   fChain->SetBranchAddress("fMCevent", &fMCevent, &b_Event_fMCevent);
   fChain->SetBranchAddress("ftE1.fUniqueID", &ftE1_fUniqueID, &b_Event_ftE1_fUniqueID);
   fChain->SetBranchAddress("ftE1.fBits", &ftE1_fBits, &b_Event_ftE1_fBits);
   fChain->SetBranchAddress("ftE1.fCharge", &ftE1_fCharge, &b_Event_ftE1_fCharge);
   fChain->SetBranchAddress("ftE1.fEta", &ftE1_fEta, &b_Event_ftE1_fEta);
   fChain->SetBranchAddress("ftE1.fPt", &ftE1_fPt, &b_Event_ftE1_fPt);
   fChain->SetBranchAddress("ftE1.fParticleType", &ftE1_fParticleType, &b_Event_ftE1_fParticleType);
   fChain->SetBranchAddress("ftE1.fParticleOrigin", &ftE1_fParticleOrigin, &b_Event_ftE1_fParticleOrigin);
   fChain->SetBranchAddress("ftE1.fTruthType", &ftE1_fTruthType, &b_Event_ftE1_fTruthType);
   fChain->SetBranchAddress("ftE1.fTruthOrigin", &ftE1_fTruthOrigin, &b_Event_ftE1_fTruthOrigin);
   fChain->SetBranchAddress("ftE1.fTruthBkgType", &ftE1_fTruthBkgType, &b_Event_ftE1_fTruthBkgType);
   fChain->SetBranchAddress("ftE1.fTruthBkgOrigin", &ftE1_fTruthBkgOrigin, &b_Event_ftE1_fTruthBkgOrigin);
   fChain->SetBranchAddress("ftE1.fTruthBkgCharge", &ftE1_fTruthBkgCharge, &b_Event_ftE1_fTruthBkgCharge);
   fChain->SetBranchAddress("ftE1.fIsElectron", &ftE1_fIsElectron, &b_Event_ftE1_fIsElectron);
   fChain->SetBranchAddress("ftE1.fFromZ", &ftE1_fFromZ, &b_Event_ftE1_fFromZ);
   fChain->SetBranchAddress("ftE1.fDeltaR", &ftE1_fDeltaR, &b_Event_ftE1_fDeltaR);
   fChain->SetBranchAddress("ftE2.fUniqueID", &ftE2_fUniqueID, &b_Event_ftE2_fUniqueID);
   fChain->SetBranchAddress("ftE2.fBits", &ftE2_fBits, &b_Event_ftE2_fBits);
   fChain->SetBranchAddress("ftE2.fCharge", &ftE2_fCharge, &b_Event_ftE2_fCharge);
   fChain->SetBranchAddress("ftE2.fEta", &ftE2_fEta, &b_Event_ftE2_fEta);
   fChain->SetBranchAddress("ftE2.fPt", &ftE2_fPt, &b_Event_ftE2_fPt);
   fChain->SetBranchAddress("ftE2.fParticleType", &ftE2_fParticleType, &b_Event_ftE2_fParticleType);
   fChain->SetBranchAddress("ftE2.fParticleOrigin", &ftE2_fParticleOrigin, &b_Event_ftE2_fParticleOrigin);
   fChain->SetBranchAddress("ftE2.fTruthType", &ftE2_fTruthType, &b_Event_ftE2_fTruthType);
   fChain->SetBranchAddress("ftE2.fTruthOrigin", &ftE2_fTruthOrigin, &b_Event_ftE2_fTruthOrigin);
   fChain->SetBranchAddress("ftE2.fTruthBkgType", &ftE2_fTruthBkgType, &b_Event_ftE2_fTruthBkgType);
   fChain->SetBranchAddress("ftE2.fTruthBkgOrigin", &ftE2_fTruthBkgOrigin, &b_Event_ftE2_fTruthBkgOrigin);
   fChain->SetBranchAddress("ftE2.fTruthBkgCharge", &ftE2_fTruthBkgCharge, &b_Event_ftE2_fTruthBkgCharge);
   fChain->SetBranchAddress("ftE2.fIsElectron", &ftE2_fIsElectron, &b_Event_ftE2_fIsElectron);
   fChain->SetBranchAddress("ftE2.fFromZ", &ftE2_fFromZ, &b_Event_ftE2_fFromZ);
   fChain->SetBranchAddress("ftE2.fDeltaR", &ftE2_fDeltaR, &b_Event_ftE2_fDeltaR);
   fChain->SetBranchAddress("ftInvMass", &ftInvMass, &b_Event_ftInvMass);
   fChain->SetBranchAddress("ftdPhi", &ftdPhi, &b_Event_ftdPhi);
   fChain->SetBranchAddress("ftSSevent", &ftSSevent, &b_Event_ftSSevent);
   fChain->SetBranchAddress("ftZEEevent", &ftZEEevent, &b_Event_ftZEEevent);

   //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//
   // USER ADDED AFTER THIS                                            //
   //\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//\\//
   
   //============================//
   // Initialize all histograms //
   //===========================//
   // Initialize error() histograms to 0
   for(int i = 0; i<4; i++){
      for(int j = 0; j<4; j++){
         hPt[i][j] = 0;
         hPt_MC[i][j] = 0;
         
         hEta[i][j] = 0;
         hEta_MC[i][j] = 0;
         
         hMass[i][j] = 0;
         hMass_MC[i][j] = 0;
         
         hEP[i][j] = 0;
         hEP_MC[i][j] = 0;
         
         hMassDist[i][j] = 0;
         hEtaDist[i][j] = 0;
         hPtDist[i][j] = 0;
         
         hMassDist_MC[i][j] = 0;
         hEtaDist_MC[i][j] = 0;
         hPtDist_MC[i][j] = 0;
         
         hEtaMass[i][j] = 0;
         /*
         for(int k = 0; k<nPtBins; k++){
            hEtaPt[k][i][j] = 0;
            hEtaPt_MC[k][i][j] = 0;
         }
         */
      }
      
   }
   
   /*
   // -------- InvMassDist HISTOGRAMS ------ //
   
   for(int i = 1; i<=3; i++){
      string sMassDistName = "hMassDist";
      string sMassDistTitle = "Invariant mass distribution of "
      switch(i){
         case LOOSE:    sMassDistName += "Loose"; break;
         case MEDIUM:   sMassDistName += "Medium"; break;
         case TIGHT:    sMassDistName += "Tight"; break;
         default continue;
      }
      for (int j = 1; j <=3; j++){
         switch(j){
            case N_PROBES:    sMassDistName += "NProbes"; sMassDistTitle += "all probes "; break;
            case N_SSPROBES:  sMassDistName += "NSSProbes"; sMassDistTitle += "SS probes "; break;
            case N_OSPROBES:  sMassDistName += "NOSProbes"; sMassDistTitle += "OS probes "; break;
            default continue;
         }
         switch(i){
            case LOOSE:    sMassDistTitle += "(Loose probe)"; break;
            case MEDIUM:   sMassDistTitle += "(Medium probe)"; break;
            case TIGHT:    sMassDistTitle += "(Tight probe)"; break;
            default continue;
         }
         string sMassDistTitleMC = sMassDistTitle;
         sMassDistTitle += ";GeV;Number of probes";
         hMassDist[i][j] = new TH1(sMassDistName.c_str(), sMassDistTitle.c_str(); Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
         
         sMassDistTitleMC += " (MC);GeV;Number of probes";
         sMassDistName += "_MC";
         hMassDist_MC[i][j] = new TH1(sMassDistName.c_str(), sMassDistTitleMC.c_str(); Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
         }
   }
   
   // ---------- EtaDist HISTOGRAMS ------------// 
   for(int i = 1; i<=3; i++){
      string sEtaDistName = "hEtaDist";
      string sEtaDistTitle = "Eta distribution of "
      switch(i){
         case LOOSE:    sEtaDistName += "Loose"; break;
         case MEDIUM:   sEtaDistName += "Medium"; break;
         case TIGHT:    sEtaDistName += "Tight"; break;
         default continue;
      }
      for (int j = 1; j <=3; j++){
         switch(j){
            case N_PROBES:    sEtaDistName += "NProbes"; sEtaDistTitle += "all probes "; break;
            case N_SSPROBES:  sEtaDistName += "NSSProbes"; sEtaDistTitle += "SS probes "; break;
            case N_OSPROBES:  sEtaDistName += "NOSProbes"; sEtaDistTitle += "OS probes "; break;
            default continue;
         }
         switch(i){
            case LOOSE:    sEtaDistTitle += "(Loose probe)"; break;
            case MEDIUM:   sEtaDistTitle += "(Medium probe)"; break;
            case TIGHT:    sEtaDistTitle += "(Tight probe)"; break;
            default continue;
         }
         string sEtaDistTitleMC = sEtaDistTitle;
         if(ABS_ETA) {sEtaDistTitle +=";|#eta|;Number of probes"; sEtaDistTitleMC += " (MC);|#eta|;Number of probes";}
         else {sEtaDistTitle +=";#eta;Number of probes"; sEtaDistTitleMC += " (MC);#eta;Number of probes";}
         
         hEtaDist[i][j] = new TH1(sEtaDistName.c_str(), sEtaDistTitle.c_str(); nEtaBins, EtaEdges);
         sEtaDistName += "_MC";
         hEtaDist_MC[i][j] = new TH1(sEtaDistName.c_str(), sEtaDistTitleMC.c_str(); nEtaBins, EtaEdges);
   }
   }
   // ------------- PtDist HISTOGRAMS --------------//
    for(int i = 1; i<=3; i++){
      string sPtDistName = "hPtDist";
      string sPtDistTitle = "p_{T} distribution of "
      switch(i){
         case LOOSE:    sPtDistName += "Loose"; break;
         case MEDIUM:   sPtDistName += "Medium"; break;
         case TIGHT:    sPtDistName += "Tight"; break;
         default continue;
      }
      for (int j = 1; j <=3; j++){
         switch(j){
            case N_PROBES:    sPtDistName += "NProbes"; sPtDistTitle += "all probes "; break;
            case N_SSPROBES:  sPtDistName += "NSSProbes"; sPtDistTitle += "SS probes "; break;
            case N_OSPROBES:  sPtDistName += "NOSProbes"; sPtDistTitle += "OS probes "; break;
            default continue;
         }
         switch(i){
            case LOOSE:    sPtDistTitle += "(Loose probe)"; break;
            case MEDIUM:   sPtDistTitle += "(Medium probe)"; break;
            case TIGHT:    sPtDistTitle += "(Tight probe)"; break;
            default continue;
         }
         string sPtDistTitleMC = sPtDistTitle;
         sPtDistTitle += ";p_{T} (GeV);Number of probes";
         sPtDistTitleMC += " (MC);p_{T} (GeV);Number of probes";
         hPtDist[i][j] = new TH1(sPtDistName.c_str(), sPtDistTitle.c_str(); nPtBins, PtEdges);
         sPtDistName += "_MC";
         hPtDist_MC[i][j] = new TH1(sPtDistName.c_str(), sPtDistTitleMC.c_str(); nPtBins, PtEdges);
      }
   }
   */
   
   // -------- MISID(PT) HISTOGRAMS ----------//
   hPt[DEFAULT][ERROR] = new TH1D("hPtDefaultError", "Charge misidentification as a function of p_{T};p_{T} (GeV);Misidentification rate", nPtBins, PtEdges);
   
   hPt[LOOSE][ERROR] = new TH1D("hPtLooseError", "Charge misidentification as a function of p_{T} (on loose cut)", nPtBins, PtEdges);
   hPt[LOOSE][N_PROBES] = new TH1D("hPtLooseNProbes", "", nPtBins, PtEdges);
   hPt[LOOSE][N_SSPROBES] = new TH1D("hPtLooseNSSProbes", "", nPtBins, PtEdges);
   hPt[LOOSE][N_OSPROBES] = new TH1D("hPtLooseNOSProbes", "", nPtBins, PtEdges);
   
   hPt[MEDIUM][ERROR] = new TH1D("hPtMediumError", "Charge misidentification as a function of p_{T} (on medium cut)", nPtBins, PtEdges);
   hPt[MEDIUM][N_PROBES] = new TH1D("hPtMediumNProbes", "", nPtBins, PtEdges);
   hPt[MEDIUM][N_SSPROBES] = new TH1D("hPtMediumNSSProbes", "", nPtBins, PtEdges);
   hPt[MEDIUM][N_OSPROBES] = new TH1D("hPtMediumNOSProbes", "", nPtBins, PtEdges);
   
   hPt[TIGHT][ERROR] = new TH1D("hPtTightError", "Charge misidentification as a function of p_{T} (on tight cut)", nPtBins, PtEdges);
   hPt[TIGHT][N_PROBES] = new TH1D("hPtTightNProbes", "", nPtBins, PtEdges);
   hPt[TIGHT][N_SSPROBES] = new TH1D("hPtTightNSSProbes", "", nPtBins, PtEdges);
   hPt[TIGHT][N_OSPROBES] = new TH1D("hPtTightNOSProbes", "", nPtBins, PtEdges);
   
   hPt_MC[DEFAULT][ERROR] = new TH1D("hPtDefaultError_MC", "Charge misidentification as a function of p_{T} (MC);p_{T} (GeV);Misidentification rate", nPtBins, PtEdges);
   hPt_MC[DEFAULT][N_PROBES] = new TH1D("hPtDefaultNProbes_MC", "", nPtBins, PtEdges);
   hPt_MC[DEFAULT][N_SSPROBES] = new TH1D("hPtDefaultNSSProbes_MC", "", nPtBins, PtEdges);
      
   hPt_MC[LOOSE][ERROR] = new TH1D("hPtLooseError_MC", "Charge misidentification as a function of p_{T} (MC loose cut)", nPtBins, PtEdges);
   hPt_MC[LOOSE][N_PROBES] = new TH1D("hPtLooseNProbes_MC", "", nPtBins, PtEdges);
   hPt_MC[LOOSE][N_SSPROBES] = new TH1D("hPtLooseNSSProbes_MC", "", nPtBins, PtEdges);
   hPt_MC[LOOSE][N_OSPROBES] = new TH1D("hPtLooseNOSProbes_MC", "", nPtBins, PtEdges);
   
   hPt_MC[MEDIUM][ERROR] = new TH1D("hPtMediumError_MC", "Charge misidentification as a function of p_{T} (MC medium cut)", nPtBins, PtEdges);
   hPt_MC[MEDIUM][N_PROBES] = new TH1D("hPtMediumNPRobes_MC", "", nPtBins, PtEdges);
   hPt_MC[MEDIUM][N_SSPROBES] = new TH1D("hPtMediumNSSProbes_MC", "", nPtBins, PtEdges);
   hPt_MC[MEDIUM][N_OSPROBES] = new TH1D("hPtMediumNOSProbes_MC", "", nPtBins, PtEdges);
   
   hPt_MC[TIGHT][ERROR] = new TH1D("hPtTightError_MC", "Charge misidentification as a function of p_{T} (MC tight cut)", nPtBins, PtEdges);
   hPt_MC[TIGHT][N_PROBES] = new TH1D("hPtTightNProbes_MC", "", nPtBins, PtEdges);
   hPt_MC[TIGHT][N_SSPROBES] = new TH1D("hPtTightNSSProbes_MC", "", nPtBins, PtEdges);
   hPt_MC[TIGHT][N_OSPROBES] = new TH1D("hPtTightNOSProbes_MC", "", nPtBins, PtEdges);
   

   // ---------- MISID(ETA) HISTOGRAMS -----------//
   hEta[DEFAULT][ERROR] = new TH1D("hEtaDefaultError", "Charge misidentification as a function of #eta;#eta;Misidentification rate", nEtaBins, EtaEdges);
   
   hEta[LOOSE][ERROR] = new TH1D("hEtaLooseError", "Charge misidentification as a function #eta (on loose cut)", nEtaBins, EtaEdges);
   hEta[LOOSE][N_PROBES] = new TH1D("hEtaLooseNProbes", "", nEtaBins, EtaEdges);
   hEta[LOOSE][N_SSPROBES] = new TH1D("hEtaLooseNSSProbes", "", nEtaBins, EtaEdges);
   hEta[LOOSE][N_OSPROBES] = new TH1D("hEtaLooseNOSProbes", "", nEtaBins, EtaEdges);
   
   hEta[MEDIUM][ERROR] = new TH1D("hEtaMediumError", "Charge misidentification as a function #eta (on medium cut)", nEtaBins, EtaEdges);
   hEta[MEDIUM][N_PROBES] = new TH1D("hEtaMediumNProbes", "", nEtaBins, EtaEdges);
   hEta[MEDIUM][N_SSPROBES] = new TH1D("hEtaMediumNSSProbes", "", nEtaBins, EtaEdges);
   hEta[MEDIUM][N_OSPROBES] = new TH1D("hEtaMediumNOSProbes", "", nEtaBins, EtaEdges);
   
   hEta[TIGHT][ERROR] = new TH1D("hEtaTightError", "Charge misidentification as a function #eta (on tight cut)", nEtaBins, EtaEdges);
   hEta[TIGHT][N_PROBES] = new TH1D("hEtaTightNProbes", "", nEtaBins, EtaEdges);
   hEta[TIGHT][N_SSPROBES] = new TH1D("hEtaTightNSSProbes", "", nEtaBins, EtaEdges);
   hEta[TIGHT][N_OSPROBES] = new TH1D("hEtaTightNOSProbes", "", nEtaBins, EtaEdges);
   
   hEta_MC[DEFAULT][ERROR] = new TH1D("hEtaDefaultError_MC", "Charge misidentification as a function #eta (MC)", nEtaBins, EtaEdges);
   hEta_MC[DEFAULT][N_PROBES] = new TH1D("hEtaDefaultNProbes_MC", "", nEtaBins, EtaEdges);
   hEta_MC[DEFAULT][N_SSPROBES] = new TH1D("hEtaDefaultNSSProbes_MC", "", nEtaBins, EtaEdges);
   
   hEta_MC[LOOSE][ERROR] = new TH1D("hEtaLooseError_MC", "Charge misidentification as a function #eta (MC loose cut)", nEtaBins, EtaEdges);
   hEta_MC[LOOSE][N_PROBES] = new TH1D("hEtaLooseNProbes_MC", "", nEtaBins, EtaEdges);
   hEta_MC[LOOSE][N_SSPROBES] = new TH1D("hEtaLooseNSSProbes_MC", "", nEtaBins, EtaEdges);
   hEta_MC[LOOSE][N_OSPROBES] = new TH1D("hEtaLooseNOSProbes_MC", "", nEtaBins, EtaEdges);
   
   hEta_MC[MEDIUM][ERROR] = new TH1D("hEtaMediumError_MC", "Charge misidentification as a function #eta (MC medium cut)", nEtaBins, EtaEdges);
   hEta_MC[MEDIUM][N_PROBES] = new TH1D("hEtaMediumNProbes_MC", "", nEtaBins, EtaEdges);
   hEta_MC[MEDIUM][N_SSPROBES] = new TH1D("hEtaMediumNSSProbes_MC", "", nEtaBins, EtaEdges);
   hEta_MC[MEDIUM][N_OSPROBES] = new TH1D("hEtaMediumNOSProbes_MC", "", nEtaBins, EtaEdges);
   
   hEta_MC[TIGHT][ERROR] = new TH1D("hEtaTightError_MC", "Charge misidentification as a function #eta (MC tight cut)", nEtaBins, EtaEdges);
   hEta_MC[TIGHT][N_PROBES] = new TH1D("hEtaTightNProbes_MC", "", nEtaBins, EtaEdges);
   hEta_MC[TIGHT][N_SSPROBES] = new TH1D("hEtaTightNSSProbes_MC", "", nEtaBins, EtaEdges);
   hEta_MC[TIGHT][N_OSPROBES] = new TH1D("hEtaTightNOSProbes_MC", "", nEtaBins, EtaEdges);

   // -------- MISID(MASS) HISTOGRAMS ----------//
   hMass[DEFAULT][ERROR] = new TH1D("hMassDefaultError", "Charge misidentification as a function of reconstructed invariant mass", Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   
   hMass[LOOSE][ERROR] = new TH1D("hMassLooseError", "Charge misidentification as a function of reconstructed invariant mass (on loose cut)", Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   hMass[LOOSE][N_PROBES] = new TH1D("hMassLooseNProbes", "", Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   hMass[LOOSE][N_SSPROBES] = new TH1D("hMassLooseNSSProbes", "", Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   hMass[LOOSE][N_OSPROBES] = new TH1D("hMassLooseNOSProbes", "", Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);

   hMass[MEDIUM][ERROR] = new TH1D("hMassMediumError", "Charge misidentification as a function of reconstructed invariant mass (on medium cut)", Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   hMass[MEDIUM][N_PROBES] = new TH1D("hMassMediumNProbes", "", Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   hMass[MEDIUM][N_SSPROBES] = new TH1D("hMassMediumNSSProbes", "", Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   hMass[MEDIUM][N_OSPROBES] = new TH1D("hMassMediumNOSProbes", "", Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   
   hMass[TIGHT][ERROR] = new TH1D("hMassTightError", "Charge misidentification as a function of reconstructed invariant mass (on tight cut)", Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   hMass[TIGHT][N_PROBES] = new TH1D("hMassTightNProbes", "", Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   hMass[TIGHT][N_SSPROBES] = new TH1D("hMassTightNSSProbes", "", Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   hMass[TIGHT][N_OSPROBES] = new TH1D("hMassTightNOSProbes", "", Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   
   hMass_MC[DEFAULT][ERROR] = new TH1D("hMassDefaultError_MC", "Charge misidentification as a function of reconstructed invariant mass (MC)", Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   hMass_MC[DEFAULT][N_PROBES] = new TH1D("hMassDefaultNProbes_MC", "", Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   hMass_MC[DEFAULT][N_SSPROBES] = new TH1D("hMassDefaultNSSProbes_MC", "", Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   
   hMass_MC[LOOSE][ERROR] = new TH1D("hMassLooseError_MC", "Charge misidentification as a function of reconstructed invariant mass (MC loose cut)", Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   hMass_MC[LOOSE][N_PROBES] = new TH1D("hMassLooseNProbes_MC", "", Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   hMass_MC[LOOSE][N_SSPROBES] = new TH1D("hMassLooseNSSProbes_MC", "", Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   hMass_MC[LOOSE][N_OSPROBES] = new TH1D("hMassLooseNOSProbes_MC", "", Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);


   hMass_MC[MEDIUM][ERROR] = new TH1D("hMassMediumError_MC", "Charge misidentification as a function of reconstructed invariant mass (MC medium cut)", Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   hMass_MC[MEDIUM][N_PROBES] = new TH1D("hMassMediumNProbes_MC", "", Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   hMass_MC[MEDIUM][N_SSPROBES] = new TH1D("hMassMediumNSSProbes_MC", "", Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   hMass_MC[MEDIUM][N_OSPROBES] = new TH1D("hMassMediumNOSProbes_MC", "", Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   
   hMass_MC[TIGHT][ERROR] = new TH1D("hMassTightError_MC", "Charge misidentification as a function of reconstructed invariant mass (MC tight cut)", Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   hMass_MC[TIGHT][N_PROBES] = new TH1D("hMassTightNProbes_MC", "", Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   hMass_MC[TIGHT][N_SSPROBES] = new TH1D("hMassTightNSSProbes_MC", "", Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   hMass_MC[TIGHT][N_OSPROBES] = new TH1D("hMassTightNOSProbes_MC", "", Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
 

   // -------- PARTICLE ORIGIN HISTOGRAMS ----------//
   
   // Implement hParticleOrigin with array notation? To take advantage of FillHist()
   hParticleOrigin = new TH1F("partOrig", "Particle Origin", 200, 0, 200);
   hParticleOriginSS = new TH1F("partOrigSS", "Origin of SS events", 200, 0, 200);
   
   // Pair origin
   hPairOrigin = new TH2D("hPairOrigin", "Origin of SS pairs", 200, 0, 200, 200, 0, 200);
   
   //hPOvsD0 = new TH2F("POvsD0", "Particle Origin vs d0;D0/D0error;SSevent", 100, -20, 20, 5,0,2);
   //hPtFromPhotonOnly = new TH1F("PTfromP", "Pt from photons only", 50, -1, 1);
   //hPtFromPEZ = new TH1F("PTfromPEZ", "PT from PEZ", 50, -1, 1);
   //hParents = new TH1D("parents", "Photon parents", 100, 20, 120);

   
   
   // ----------- MISID(PT, ETA) HISTOGRAMS ---------------- //
   hEP[LOOSE][ERROR] = new TH2D("hEPLooseError", "Charge misID as a function of #eta and p_{T} (on loose cuts);#eta;;p_{T}", nEtaBins, EtaEdges, nPtBins, PtEdges);
   hEP[LOOSE][N_PROBES] = new TH2D("hEPLooseNProbes", "", nEtaBins, EtaEdges, nPtBins, PtEdges);
   hEP[LOOSE][N_SSPROBES] = new TH2D("hEPLooseNSSProbes", "", nEtaBins, EtaEdges, nPtBins, PtEdges);
   hEP[LOOSE][N_OSPROBES] = new TH2D("hEPLooseNOSProbes", "", nEtaBins, EtaEdges, nPtBins, PtEdges);

   
   hEP[MEDIUM][ERROR] = new TH2D("hEPMediumError", "Charge misID as a function of #eta and p_{T} (on medium cuts);#eta;p_{T} (GeV)", nEtaBins, EtaEdges, nPtBins, PtEdges);
   hEP[MEDIUM][N_PROBES] = new TH2D("hEPMediumNProbes", "", nEtaBins, EtaEdges, nPtBins, PtEdges);
   hEP[MEDIUM][N_SSPROBES] = new TH2D("hEPMediumNSSProbes", "", nEtaBins, EtaEdges, nPtBins, PtEdges);
   hEP[MEDIUM][N_OSPROBES] = new TH2D("hEPMediumNOSProbes", "", nEtaBins, EtaEdges, nPtBins, PtEdges);
   
   hEP[TIGHT][ERROR] = new TH2D("hEPTightError", "Charge misID as a function of #eta and p_{T} (on tight cuts);#eta;p_{T} (GeV)", nEtaBins, EtaEdges, nPtBins, PtEdges);
   hEP[TIGHT][N_PROBES] = new TH2D("hEPTightNProbes", "", nEtaBins, EtaEdges, nPtBins, PtEdges);
   hEP[TIGHT][N_SSPROBES] = new TH2D("hEPTightNSSProbes", "", nEtaBins, EtaEdges, nPtBins, PtEdges);
   hEP[TIGHT][N_OSPROBES] = new TH2D("hEPTightNOSProbes", "", nEtaBins, EtaEdges, nPtBins, PtEdges);
   
   hEP_MC[LOOSE][ERROR] = new TH2D("hEPLooseError_MC", "Charge misID as a function of #eta and p_{T} (MC loose cuts);#eta;p_{T} (GeV)", nEtaBins, EtaEdges, nPtBins, PtEdges);
   hEP_MC[LOOSE][N_PROBES] = new TH2D("hEPLooseNProbes_MC", "", nEtaBins, EtaEdges, nPtBins, PtEdges);
   hEP_MC[LOOSE][N_SSPROBES] = new TH2D("hEPLooseNSSProbes_MC", "", nEtaBins, EtaEdges, nPtBins, PtEdges);
   hEP_MC[LOOSE][N_OSPROBES] = new TH2D("hEPLooseNOSProbes_MC", "", nEtaBins, EtaEdges, nPtBins, PtEdges);
   
   hEP_MC[MEDIUM][ERROR] = new TH2D("hEPMediumError_MC", "Charge misID as a function of #eta and p_{T} (MC medium cuts);#eta;p_{T} (GeV)", nEtaBins, EtaEdges, nPtBins, PtEdges);
   hEP_MC[MEDIUM][N_PROBES] = new TH2D("hEPMediumNProbes_MC", "", nEtaBins, EtaEdges, nPtBins, PtEdges);
   hEP_MC[MEDIUM][N_SSPROBES] = new TH2D("hEPMediumNSSProbes_MC", "", nEtaBins, EtaEdges, nPtBins, PtEdges);
   hEP_MC[MEDIUM][N_OSPROBES] = new TH2D("hEPMediumNOSProbes_MC", "", nEtaBins, EtaEdges, nPtBins, PtEdges);
   
   hEP_MC[TIGHT][ERROR] = new TH2D("hEPTightError_MC", "Charge misID as a function of #eta and p_{T} (MC tight cuts);#eta;p_{T} (GeV)", nEtaBins, EtaEdges, nPtBins, PtEdges);
   hEP_MC[TIGHT][N_PROBES] = new TH2D("hEPTightNProbes_MC", "", nEtaBins, EtaEdges, nPtBins, PtEdges);
   hEP_MC[TIGHT][N_SSPROBES] = new TH2D("hEPTightNSSProbes_MC", "", nEtaBins, EtaEdges, nPtBins, PtEdges);
   hEP_MC[TIGHT][N_OSPROBES] = new TH2D("hEPTightNOSProbes_MC", "", nEtaBins, EtaEdges, nPtBins, PtEdges);
   
   
   // =============== FINDING NOISE =================== //
   hEtaMass[TIGHT][ERROR] = new TH2D("hEtaMassTightError", "EtaMass MisID;#eta;mll (GeV)", nEtaBins, EtaEdges, Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   hEtaMass[TIGHT][N_PROBES] = new TH2D("hEtaMassTightNProbes", "EtaMass NProbes;#eta;mll (GeV)", nEtaBins, EtaEdges, Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   hEtaMass[TIGHT][N_SSPROBES] = new TH2D("hEtaMassTightNSSProbes", "EtaMass NSSProbes;#eta;mll (GeV)",  nEtaBins, EtaEdges, Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   hEtaMass[TIGHT][N_OSPROBES] = new TH2D("hEtaMassTightNOSProbes", "EtaMass NOSProbes;#eta;mll (GeV)", nEtaBins, EtaEdges, Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   
   hPtMass[TIGHT][ERROR] = new TH2D("hPtMassTightError", "PtMass MisID;p_{T} (GeV);mll (GeV)", nPtBins, PtEdges, Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   hPtMass[TIGHT][N_PROBES] = new TH2D("hPtMassTightNProbes", "PtMass NProbes;p_{T} (GeV);mll (GeV)", nPtBins, PtEdges, Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   hPtMass[TIGHT][N_SSPROBES] = new TH2D("hPtMassTightNSSProbes", "PtMass NSSProbes;p_{T} (GeV);mll (GeV)",  nPtBins, PtEdges, Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   hPtMass[TIGHT][N_OSPROBES] = new TH2D("hPtMassTightNOSProbes", "PtMass NOSProbes;p_{T} (GeV);mll (GeV)", nPtBins, PtEdges, Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   
   hEtaMassNoZ[TIGHT][ERROR] = new TH2D("hEtaMassNoZTightError", "EtaMass NoZ MisID;#eta;mll (GeV)", nEtaBins, EtaEdges, Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   hEtaMassNoZ[TIGHT][N_PROBES] = new TH2D("hEtaMassNoZTightNProbes", "EtaMass NoZ NProbes;#eta;mll (GeV)", nEtaBins, EtaEdges, Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   hEtaMassNoZ[TIGHT][N_SSPROBES] = new TH2D("hEtaMassNoZTightNSSProbes", "EtaMass NoZ NSSProbes;#eta;mll (GeV)",  nEtaBins, EtaEdges, Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   hEtaMassNoZ[TIGHT][N_OSPROBES] = new TH2D("hEtaMassNoZTightNOSProbes", "EtaMass NoZ NOSProbes;#eta;mll (GeV)", nEtaBins, EtaEdges, Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   
   
   hPtMassNoZ[TIGHT][ERROR] = new TH2D("hPtMassNoZTightError", "PtMass NoZ MisID;p_{T} (GeV);mll (GeV)", nPtBins, PtEdges, Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   hPtMassNoZ[TIGHT][N_PROBES] = new TH2D("hPtMassNoZTightNProbes", "PtMass NoZ NProbes;p_{T} (GeV);mll (GeV)", nPtBins, PtEdges, Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   hPtMassNoZ[TIGHT][N_SSPROBES] = new TH2D("hPtMassNoZTightNSSProbes", "PtMass NoZ NSSProbes;p_{T} (GeV);mll (GeV)",  nPtBins, PtEdges, Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   hPtMassNoZ[TIGHT][N_OSPROBES] = new TH2D("hPtMassNoZTightNOSProbes", "PtMass NoZ NOSProbes;p_{T} (GeV);mll (GeV)", nPtBins, PtEdges, Z_MASS_UPPER_LIMIT - Z_MASS_LOWER_LIMIT, Z_MASS_LOWER_LIMIT, Z_MASS_UPPER_LIMIT);
   /* OLD */
   /*
   for(int i = 0; i<nPtBins; i++){
      std::ostringstream lowPt;
      lowPt << std::fixed << std::setprecision(0) << PtEdges[i];
      
      std::ostringstream highPt;
      highPt << std::fixed << std::setprecision(0) << PtEdges[i+1];
      
      
      string sHistName = "Charge misID as a function of |#eta| in p_{T} range " + lowPt.str() + "-" + highPt.str() + " GeV;#eta;Error rate";
      string sPadName = "ErrorEtaPt_" + lowPt.str();
      
      hEtaPt[i][ERROR] = new TH1D(sPadName.c_str(), sHistName.c_str(), nEtaBins, EtaEdges);
      hEtaPt[i][N_PROBES] = 0;
      hEtaPt[i][N_SSPROBES] = 0;
            
      hEtaPtLoose[i][ERROR] = new TH1D((sPadName+"L").c_str(), (sHistName + " (loose cuts)").c_str(), nEtaBins, EtaEdges);
      hEtaPtMedium[i][ERROR] = new TH1D((sPadName+"M").c_str(), (sHistName + " (medium cuts)").c_str(), nEtaBins, EtaEdges);
      hEtaPtTight[i][ERROR] = new TH1D((sPadName+"T").c_str(), (sHistName + " (tight cuts)").c_str(), nEtaBins, EtaEdges);
      
      hEtaPt_MC[i][ERROR] = new TH1D((sPadName+"MC").c_str(), (sHistName + " (MC)").c_str(), nEtaBins, EtaEdges);
      
      hEtaPtLoose[i][N_PROBES] = new TH1D((sPadName+"nL").c_str(), "", nEtaBins, EtaEdges);
      hEtaPtMedium[i][N_PROBES] = new TH1D((sPadName+"nM").c_str(), "", nEtaBins, EtaEdges);
      hEtaPtTight[i][N_PROBES] = new TH1D((sPadName+"nT").c_str(), "", nEtaBins, EtaEdges);
      hEtaPt_MC[i][N_PROBES] = new TH1D((sPadName+"nMC").c_str(), "", nEtaBins, EtaEdges);
      
      hEtaPtLoose[i][N_SSPROBES] = new TH1D((sPadName+"sL").c_str(), "", nEtaBins, EtaEdges);
      hEtaPtMedium[i][N_SSPROBES] = new TH1D((sPadName+"sM").c_str(), "", nEtaBins, EtaEdges);
      hEtaPtTight[i][N_SSPROBES] = new TH1D((sPadName+"sT").c_str(), "", nEtaBins, EtaEdges);
      hEtaPt_MC[i][N_SSPROBES] = new TH1D((sPadName+"sMC").c_str(), "", nEtaBins, EtaEdges);
   }
   //*/
   
   //======== Statistics =============//
   hTagCuts = new TH1I("hTagCuts", "Number of electrons passing tag cuts", 6, 0, 6);
   hTagCuts->GetXaxis()->SetBinLabel(1, "Passed event selection");
   hTagCuts->GetXaxis()->SetBinLabel(2, "Tight electron");
   hTagCuts->GetXaxis()->SetBinLabel(3, "|#eta| < 1.37");
   hTagCuts->GetXaxis()->SetBinLabel(4, "d0/#sigma_{d0} < 1.5");
   hTagCuts->GetXaxis()->SetBinLabel(5, "At least 1 B-layer hit");
   hTagCuts->GetXaxis()->SetBinLabel(6, "At least 9 Si-layer hits");   
   
   //===============================//
   // Initialize output file stream //
   //===============================//
   file_eParents.open("parents.txt");
   file_nProbes.open("nProbes.csv");
   file_misID.open("misID.csv");
   
   

   PRINT_STAT("Initialization", "Complete");
}


// Doesn't really do anything
Bool_t TagProbePlot::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

#endif // corresponds to #ifdef TagProbePlot_cxx
