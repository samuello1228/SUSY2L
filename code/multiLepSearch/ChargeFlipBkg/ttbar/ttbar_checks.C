#define ttbar_checks_cxx
#include "ttbar_checks.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

// Progress bar 
void loadbar(unsigned int x, unsigned int n, unsigned int w = 50)
{
   if ( (x != n) && (x % (n/100+1) != 0) ) return;

   float ratio    =  x/(float)n;
   unsigned int c =  ratio * w;

   std::cout << setw(3) << (int)(ratio*100) << "% [";
   for (unsigned int x=0; x<c; x++) std::cout << "=";
   for (unsigned int x=c; x<w; x++) std::cout << " ";
   std::cout << "]\r" << flush;
}


void ttbar_checks::Loop()
{
//   In a ROOT session, you can do:
//      root> .L ttbar_checks.C
//      root> ttbar_checks t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   bool onlySignal = false;

   TFile* fOut = TFile::Open("ttbarHistos.root", "recreate");
   TH1D* hLeadingPtOS = new TH1D("hLeadingPtOS", "Leading p_{T};leading p_{T} (GeV);Number of weighted events", 200, 20, 220);
   TH1D* hLeadingPtSS = new TH1D("hLeadingPtSS", "Leading p_{T};leading p_{T} (GeV);Number of weighted events", 200, 20, 220);
   TH1D* hSubleadingPtOS = new TH1D("hSubleadingPtOS", "Subleading p_{T};subleading p_{T} (GeV);Number of weighted events", 200, 20, 220);
   TH1D* hSubleadingPtSS = new TH1D("hSubleadingPtSS", "Subleading p_{T};subleading p_{T} (GeV);Number of weighted events", 200, 20, 220);
   TH1D* hInvMassOS = new TH1D("hInvMassOS", "Invariant mass;m_{ee} (GeV);Number of weighted events", 200, 20, 220);
   TH1D* hInvMassSS = new TH1D("hInvMassSS", "Invariant mass;m_{ee} (GeV);Number of weighted events", 200, 20, 220);
   TH1D* hLeadingEtaOS = new TH1D("hLeadingEtaOS", "Leading |#eta|;leading |#eta|;Number of weighted events", 100, 0, 2.47);
   TH1D* hLeadingEtaSS = new TH1D("hLeadingEtaSS", "Leading |#eta|;leading |#eta|;Number of weighted events", 100, 0, 2.47);
   TH1D* hSubleadingEtaOS = new TH1D("hSubleadingEtaOS", "Subeading |#eta|;subleading |#eta|;Number of weighted events", 100, 0, 2.47);
   TH1D* hSubleadingEtaSS = new TH1D("hSubleadingEtaSS", "Subeading |#eta|;subleading |#eta|;Number of weighted events", 100, 0, 2.47);
   TH1D* hLeadingPhiOS = new TH1D("hLeadingPhiOS", "Leading #phi;leading #phi;Number of weighted events", 100, -3.15, 3.15);
   TH1D* hLeadingPhiSS = new TH1D("hLeadingPhiSS", "Leading #phi;leading #phi;Number of weighted events", 100, -3.15, 3.15);
   TH1D* hSubleadingPhiOS = new TH1D("hSubleadingPhiOS", "Subleading #phi;subleading #phi;Number of weighted events", 100, -3.15, 3.15);
   TH1D* hSubleadingPhiSS = new TH1D("hSubleadingPhiSS", "Subleading #phi;subleading #phi;Number of weighted events", 100, -3.15, 3.15);
   hLeadingPtOS->SetDirectory(fOut);
   hLeadingPtSS->SetDirectory(fOut);
   hSubleadingPtOS->SetDirectory(fOut);
   hSubleadingPtSS->SetDirectory(fOut);
   hInvMassOS->SetDirectory(fOut);
   hInvMassSS->SetDirectory(fOut);
   hLeadingEtaOS->SetDirectory(fOut);
   hLeadingEtaSS->SetDirectory(fOut);
   hSubleadingEtaOS->SetDirectory(fOut);
   hSubleadingEtaSS->SetDirectory(fOut);
   hLeadingPhiOS->SetDirectory(fOut);
   hLeadingPhiSS->SetDirectory(fOut);
   hSubleadingPhiOS->SetDirectory(fOut);
   hSubleadingPhiSS->SetDirectory(fOut);

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
   	loadbar(jentry, nentries);
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

		bool baseCut = sig_trigCode>0 && leps_ == 2 && int(abs(leps_ID[0]/1000)) == 11 && int(abs(leps_ID[1]/1000)) == 11;
		if (!baseCut) continue;

      int nJets(0), nBJets(0);
      for(int i=0; i<jets_; i++){
         if(jets_jFlag[i] & 2) nJets++;
         if(jets_jFlag[i] & 32) nBJets++;
      }
      if(nJets<2) continue;
      if(nBJets!=2) continue;

      if(onlySignal && !(((leps_lFlag[0] & 2)/2) && ((leps_lFlag[1] & 2)/2))) continue;

      bool SSevent = ((leps_ID[0]>0) == (leps_ID[1]>0));
      double w = evt_weight*evt_pwt;
		
		if(SSevent){
		   hLeadingPtSS->Fill(leps_pt[0], w);
			hSubleadingPtSS->Fill(leps_pt[1], w);
			hInvMassSS->Fill(l12_m, w);
			hLeadingEtaSS->Fill(leps_eta[0], w);
			hSubleadingEtaSS->Fill(leps_eta[1], w);
			hLeadingPhiSS->Fill(leps_phi[0], w);
			hSubleadingPhiSS->Fill(leps_phi[1], w);
		}
		else {
		   hLeadingPtOS->Fill(leps_pt[0], w);
			hSubleadingPtOS->Fill(leps_pt[1], w);
			hInvMassOS->Fill(l12_m, w);
			hLeadingEtaOS->Fill(leps_eta[0], w);
			hSubleadingEtaOS->Fill(leps_eta[1], w);
			hLeadingPhiOS->Fill(leps_phi[0], w);
			hSubleadingPhiOS->Fill(leps_phi[1], w);
		}
   }
   std::cout << std::endl;
   fOut->Write();
}
