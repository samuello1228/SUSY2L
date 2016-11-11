//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jul 27 12:07:10 2016 by ROOT version 6.04/14
// from TTree misIDvars/Charge misID variables
// found on file: 0727-MC-loose/misIDdecorations.root
//////////////////////////////////////////////////////////

#ifndef misIDvars_h
#define misIDvars_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class misIDvars {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Double_t        e1rate;
   Double_t        e1rateErr;
   Double_t        e2rate;
   Double_t        e2rateErr;
   Double_t        e1dPt;
   Double_t        e1dPtErr;
   Double_t        e2dPt;
   Double_t        e2dPtErr;
   Double_t        chargeFlipWeight;

   // List of branches
   TBranch        *b_e1rate;   //!
   TBranch        *b_e1rateErr;   //!
   TBranch        *b_e2rate;   //!
   TBranch        *b_e2rateErr;   //!
   TBranch        *b_e1dPt;   //!
   TBranch        *b_e1dPtErr;   //!
   TBranch        *b_e2dPt;   //!
   TBranch        *b_e2dPtErr;   //!
   TBranch        *b_chargeFlipWeight;   //!

   misIDvars(TTree *tree=0);
   virtual ~misIDvars();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef misIDvars_cxx
misIDvars::misIDvars(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("0727-MC-loose/misIDdecorations.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("0727-MC-loose/misIDdecorations.root");
      }
      f->GetObject("misIDvars",tree);

   }
   Init(tree);
}

misIDvars::~misIDvars()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t misIDvars::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t misIDvars::LoadTree(Long64_t entry)
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

void misIDvars::Init(TTree *tree)
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

   fChain->SetBranchAddress("e1rate", &e1rate, &b_e1rate);
   fChain->SetBranchAddress("e1rateErr", &e1rateErr, &b_e1rateErr);
   fChain->SetBranchAddress("e2rate", &e2rate, &b_e2rate);
   fChain->SetBranchAddress("e2rateErr", &e2rateErr, &b_e2rateErr);
   fChain->SetBranchAddress("e1dPt", &e1dPt, &b_e1dPt);
   fChain->SetBranchAddress("e1dPtErr", &e1dPtErr, &b_e1dPtErr);
   fChain->SetBranchAddress("e2dPt", &e2dPt, &b_e2dPt);
   fChain->SetBranchAddress("e2dPtErr", &e2dPtErr, &b_e2dPtErr);
   fChain->SetBranchAddress("chargeFlipWeight", &chargeFlipWeight, &b_chargeFlipWeight);
   Notify();
}

Bool_t misIDvars::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void misIDvars::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t misIDvars::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef misIDvars_cxx
