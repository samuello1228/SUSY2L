//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jul  7 14:50:19 2016 by ROOT version 6.04/14
// from TTree evt2l/a angles tree
// found on file: /afs/cern.ch/user/g/ggallard/work/NTuples/user.clo.v7.1.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_myOutput.root/user.clo.8907010._000001.myOutput.root
//////////////////////////////////////////////////////////

#ifndef evt2l_h
#define evt2l_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class evt2l {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
   static const Int_t kMaxleps = 4;
   static const Int_t kMaxjets = 20;
   static const Int_t kMaxtruths = 13;

   // Declaration of leaf types
   ULong64_t       evt_run;
   ULong64_t       evt_event;
   ULong64_t       evt_lumiBlock;
   ULong64_t       evt_actualMu;
   UInt_t          evt_index;
   UInt_t          evt_cuts;
   UInt_t          evt_trig;
   Int_t           evt_flag;
   Int_t           evt_pass;
   Float_t         evt_averageMu;
   Float_t         evt_Xsec;
   Float_t         evt_weight;
   Float_t         evt_pwt;
   Float_t         evt_ElSF;
   Float_t         evt_MuSF;
   Float_t         evt_BtagSF;
   Float_t         evt_qFwt;
   Int_t           leps_;
   Float_t         leps_pt[kMaxleps];   //[leps_]
   Float_t         leps_eta[kMaxleps];   //[leps_]
   Float_t         leps_phi[kMaxleps];   //[leps_]
   Float_t         leps_MET_dPhi[kMaxleps];   //[leps_]
   Float_t         leps_topoetcone20[kMaxleps];   //[leps_]
   Float_t         leps_topoetcone30[kMaxleps];   //[leps_]
   Float_t         leps_topoetcone40[kMaxleps];   //[leps_]
   Float_t         leps_ptcone20[kMaxleps];   //[leps_]
   Float_t         leps_ptcone30[kMaxleps];   //[leps_]
   Float_t         leps_ptcone40[kMaxleps];   //[leps_]
   Float_t         leps_mT[kMaxleps];   //[leps_]
   Float_t         leps_d0[kMaxleps];   //[leps_]
   Float_t         leps_d0Err[kMaxleps];   //[leps_]
   Float_t         leps_z0[kMaxleps];   //[leps_]
   Float_t         leps_z0Err[kMaxleps];   //[leps_]
   Float_t         leps_z0sinTheta[kMaxleps];   //[leps_]
   Float_t         leps_d0sig[kMaxleps];   //[leps_]
   Float_t         leps_truthProb[kMaxleps];   //[leps_]
   Float_t         leps_SF_Loose[kMaxleps];   //[leps_]
   Float_t         leps_SF_Medium[kMaxleps];   //[leps_]
   Float_t         leps_SF_Tight[kMaxleps];   //[leps_]
   Float_t         leps_wt1[kMaxleps];   //[leps_]
   Float_t         leps_wt2[kMaxleps];   //[leps_]
   Float_t         leps_wt3[kMaxleps];   //[leps_]
   Int_t           leps_ID[kMaxleps];   //[leps_]
   Int_t           leps_author[kMaxleps];   //[leps_]
   Int_t           leps_truthI[kMaxleps];   //[leps_]
   Int_t           leps_truthType[kMaxleps];   //[leps_]
   Int_t           leps_truthOrig[kMaxleps];   //[leps_]
   Int_t           leps_isTight[kMaxleps];   //[leps_]
   UInt_t          leps_Q[kMaxleps];   //[leps_]
   UInt_t          leps_lFlag[kMaxleps];   //[leps_]
   UInt_t          leps_nBHits[kMaxleps];   //[leps_]
   UInt_t          leps_nPixHits[kMaxleps];   //[leps_]
   UInt_t          leps_nSCTHits[kMaxleps];   //[leps_]
   UInt_t          leps_nPixHoles[kMaxleps];   //[leps_]
   UInt_t          leps_nSCTHoles[kMaxleps];   //[leps_]
   UInt_t          leps_nTRTHits[kMaxleps];   //[leps_]
   UInt_t          leps_nTRTOutliers[kMaxleps];   //[leps_]
   Float_t         l12_pt;
   Float_t         l12_eta;
   Float_t         l12_phi;
   Float_t         l12_MET_dPhi;
   Float_t         l12_m;
   Float_t         l12_dPhi;
   Float_t         l12_dR;
   Int_t           jets_;
   Float_t         jets_pt[kMaxjets];   //[jets_]
   Float_t         jets_eta[kMaxjets];   //[jets_]
   Float_t         jets_phi[kMaxjets];   //[jets_]
   Float_t         jets_MET_dPhi[kMaxjets];   //[jets_]
   UInt_t          jets_jFlag[kMaxjets];   //[jets_]
   Float_t         jets_e[kMaxjets];   //[jets_]
   Int_t           truths_;
   Float_t         truths_pt[kMaxtruths];   //[truths_]
   Float_t         truths_eta[kMaxtruths];   //[truths_]
   Float_t         truths_phi[kMaxtruths];   //[truths_]
   Int_t           truths_pdgId[kMaxtruths];   //[truths_]
   Int_t           truths_barcode[kMaxtruths];   //[truths_]
   Int_t           truths_motherI[kMaxtruths];   //[truths_]
   Int_t           truths_matchI[kMaxtruths];   //[truths_]
   ULong64_t       sig_trigCode;
   Float_t         sig_Met;
   Float_t         sig_MetRel;
   Float_t         sig_MetX;
   Float_t         sig_MetY;
   Float_t         sig_mT2;
   Float_t         sig_HT;
   UInt_t          sig_nEl;
   UInt_t          sig_nMu;
   UInt_t          sig_nTau;
   UInt_t          sig_nJet;
   UInt_t          sig_nPV;
   UInt_t          sig_nVtx;

   // List of branches
   TBranch        *b_evt;   //!
   TBranch        *b_leps_;   //!
   TBranch        *b_leps_pt;   //!
   TBranch        *b_leps_eta;   //!
   TBranch        *b_leps_phi;   //!
   TBranch        *b_leps_MET_dPhi;   //!
   TBranch        *b_leps_topoetcone20;   //!
   TBranch        *b_leps_topoetcone30;   //!
   TBranch        *b_leps_topoetcone40;   //!
   TBranch        *b_leps_ptcone20;   //!
   TBranch        *b_leps_ptcone30;   //!
   TBranch        *b_leps_ptcone40;   //!
   TBranch        *b_leps_mT;   //!
   TBranch        *b_leps_d0;   //!
   TBranch        *b_leps_d0Err;   //!
   TBranch        *b_leps_z0;   //!
   TBranch        *b_leps_z0Err;   //!
   TBranch        *b_leps_z0sinTheta;   //!
   TBranch        *b_leps_d0sig;   //!
   TBranch        *b_leps_truthProb;   //!
   TBranch        *b_leps_SF_Loose;   //!
   TBranch        *b_leps_SF_Medium;   //!
   TBranch        *b_leps_SF_Tight;   //!
   TBranch        *b_leps_wt1;   //!
   TBranch        *b_leps_wt2;   //!
   TBranch        *b_leps_wt3;   //!
   TBranch        *b_leps_ID;   //!
   TBranch        *b_leps_author;   //!
   TBranch        *b_leps_truthI;   //!
   TBranch        *b_leps_truthType;   //!
   TBranch        *b_leps_truthOrig;   //!
   TBranch        *b_leps_isTight;   //!
   TBranch        *b_leps_Q;   //!
   TBranch        *b_leps_lFlag;   //!
   TBranch        *b_leps_nBHits;   //!
   TBranch        *b_leps_nPixHits;   //!
   TBranch        *b_leps_nSCTHits;   //!
   TBranch        *b_leps_nPixHoles;   //!
   TBranch        *b_leps_nSCTHoles;   //!
   TBranch        *b_leps_nTRTHits;   //!
   TBranch        *b_leps_nTRTOutliers;   //!
   TBranch        *b_l12;   //!
   TBranch        *b_jets_;   //!
   TBranch        *b_jets_pt;   //!
   TBranch        *b_jets_eta;   //!
   TBranch        *b_jets_phi;   //!
   TBranch        *b_jets_MET_dPhi;   //!
   TBranch        *b_jets_jFlag;   //!
   TBranch        *b_jets_e;   //!
   TBranch        *b_truths_;   //!
   TBranch        *b_truths_pt;   //!
   TBranch        *b_truths_eta;   //!
   TBranch        *b_truths_phi;   //!
   TBranch        *b_truths_pdgId;   //!
   TBranch        *b_truths_barcode;   //!
   TBranch        *b_truths_motherI;   //!
   TBranch        *b_truths_matchI;   //!
   TBranch        *b_sig;   //!

   evt2l(TTree *tree=0);
   virtual ~evt2l();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef evt2l_cxx
evt2l::evt2l(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/afs/cern.ch/user/g/ggallard/work/NTuples/user.clo.v7.1.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_myOutput.root/user.clo.8907010._000001.myOutput.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/afs/cern.ch/user/g/ggallard/work/NTuples/user.clo.v7.1.361106.PowhegPythia8EvtGen_AZNLOCTEQ6L1_Zee_myOutput.root/user.clo.8907010._000001.myOutput.root");
      }
      f->GetObject("evt2l",tree);

   }
   Init(tree);
}

evt2l::~evt2l()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t evt2l::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t evt2l::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void evt2l::Init(TTree *tree)
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
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("evt", &evt_run, &b_evt);
   fChain->SetBranchAddress("leps", &leps_, &b_leps_);
   fChain->SetBranchAddress("leps.pt", leps_pt, &b_leps_pt);
   fChain->SetBranchAddress("leps.eta", leps_eta, &b_leps_eta);
   fChain->SetBranchAddress("leps.phi", leps_phi, &b_leps_phi);
   fChain->SetBranchAddress("leps.MET_dPhi", leps_MET_dPhi, &b_leps_MET_dPhi);
   fChain->SetBranchAddress("leps.topoetcone20", leps_topoetcone20, &b_leps_topoetcone20);
   fChain->SetBranchAddress("leps.topoetcone30", leps_topoetcone30, &b_leps_topoetcone30);
   fChain->SetBranchAddress("leps.topoetcone40", leps_topoetcone40, &b_leps_topoetcone40);
   fChain->SetBranchAddress("leps.ptcone20", leps_ptcone20, &b_leps_ptcone20);
   fChain->SetBranchAddress("leps.ptcone30", leps_ptcone30, &b_leps_ptcone30);
   fChain->SetBranchAddress("leps.ptcone40", leps_ptcone40, &b_leps_ptcone40);
   fChain->SetBranchAddress("leps.mT", leps_mT, &b_leps_mT);
   fChain->SetBranchAddress("leps.d0", leps_d0, &b_leps_d0);
   fChain->SetBranchAddress("leps.d0Err", leps_d0Err, &b_leps_d0Err);
   fChain->SetBranchAddress("leps.z0", leps_z0, &b_leps_z0);
   fChain->SetBranchAddress("leps.z0Err", leps_z0Err, &b_leps_z0Err);
   fChain->SetBranchAddress("leps.z0sinTheta", leps_z0sinTheta, &b_leps_z0sinTheta);
   fChain->SetBranchAddress("leps.d0sig", leps_d0sig, &b_leps_d0sig);
   fChain->SetBranchAddress("leps.truthProb", leps_truthProb, &b_leps_truthProb);
   fChain->SetBranchAddress("leps.SF_Loose", leps_SF_Loose, &b_leps_SF_Loose);
   fChain->SetBranchAddress("leps.SF_Medium", leps_SF_Medium, &b_leps_SF_Medium);
   fChain->SetBranchAddress("leps.SF_Tight", leps_SF_Tight, &b_leps_SF_Tight);
   fChain->SetBranchAddress("leps.wt1", leps_wt1, &b_leps_wt1);
   fChain->SetBranchAddress("leps.wt2", leps_wt2, &b_leps_wt2);
   fChain->SetBranchAddress("leps.wt3", leps_wt3, &b_leps_wt3);
   fChain->SetBranchAddress("leps.ID", leps_ID, &b_leps_ID);
   fChain->SetBranchAddress("leps.author", leps_author, &b_leps_author);
   fChain->SetBranchAddress("leps.truthI", leps_truthI, &b_leps_truthI);
   fChain->SetBranchAddress("leps.truthType", leps_truthType, &b_leps_truthType);
   fChain->SetBranchAddress("leps.truthOrig", leps_truthOrig, &b_leps_truthOrig);
   fChain->SetBranchAddress("leps.isTight", leps_isTight, &b_leps_isTight);
   fChain->SetBranchAddress("leps.Q", leps_Q, &b_leps_Q);
   fChain->SetBranchAddress("leps.lFlag", leps_lFlag, &b_leps_lFlag);
   fChain->SetBranchAddress("leps.nBHits", leps_nBHits, &b_leps_nBHits);
   fChain->SetBranchAddress("leps.nPixHits", leps_nPixHits, &b_leps_nPixHits);
   fChain->SetBranchAddress("leps.nSCTHits", leps_nSCTHits, &b_leps_nSCTHits);
   fChain->SetBranchAddress("leps.nPixHoles", leps_nPixHoles, &b_leps_nPixHoles);
   fChain->SetBranchAddress("leps.nSCTHoles", leps_nSCTHoles, &b_leps_nSCTHoles);
   fChain->SetBranchAddress("leps.nTRTHits", leps_nTRTHits, &b_leps_nTRTHits);
   fChain->SetBranchAddress("leps.nTRTOutliers", leps_nTRTOutliers, &b_leps_nTRTOutliers);
   fChain->SetBranchAddress("l12", &l12_pt, &b_l12);
   fChain->SetBranchAddress("jets", &jets_, &b_jets_);
   fChain->SetBranchAddress("jets.pt", jets_pt, &b_jets_pt);
   fChain->SetBranchAddress("jets.eta", jets_eta, &b_jets_eta);
   fChain->SetBranchAddress("jets.phi", jets_phi, &b_jets_phi);
   fChain->SetBranchAddress("jets.MET_dPhi", jets_MET_dPhi, &b_jets_MET_dPhi);
   fChain->SetBranchAddress("jets.jFlag", jets_jFlag, &b_jets_jFlag);
   fChain->SetBranchAddress("jets.e", jets_e, &b_jets_e);
   fChain->SetBranchAddress("truths", &truths_, &b_truths_);
   fChain->SetBranchAddress("truths.pt", truths_pt, &b_truths_pt);
   fChain->SetBranchAddress("truths.eta", truths_eta, &b_truths_eta);
   fChain->SetBranchAddress("truths.phi", truths_phi, &b_truths_phi);
   fChain->SetBranchAddress("truths.pdgId", truths_pdgId, &b_truths_pdgId);
   fChain->SetBranchAddress("truths.barcode", truths_barcode, &b_truths_barcode);
   fChain->SetBranchAddress("truths.motherI", truths_motherI, &b_truths_motherI);
   fChain->SetBranchAddress("truths.matchI", truths_matchI, &b_truths_matchI);
   fChain->SetBranchAddress("sig", &sig_trigCode, &b_sig);
   Notify();
}

Bool_t evt2l::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void evt2l::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t evt2l::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef evt2l_cxx
