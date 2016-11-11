#define analysis_cxx
// The class definition in analysis.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// root> T->Process("analysis.C")
// root> T->Process("analysis.C","some options")
// root> T->Process("analysis.C+")


#include "analysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <iostream>
#include <AnglesUtil.hpp>
#include <TMath.h>
#include <TLorentzVector.h>

using namespace std;

void analysis::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
   //_file = new TFile("output/output.root","recreate");
   _ttree = new TTree("t_evt3l","t_evt3l");
   _ltree = new TTree("l_evt3l","l_evt3l");
   //_mtree = new TTree("m_evt3l","m_evt3l");
   _hcutflow2 = new TH1F("hcutflow2", "cut flow 2", 10, 0, 10);
   _hcutflow2->SetDirectory(0);


   _ttree->Branch("lep_pt", &lep_pt, "lep_pt/F");
   _ttree->Branch("lep_eta", &lep_eta, "lep_eta/F");
   _ttree->Branch("lep_SF",&lep_SF,"lep_SF/F");
   _ttree->Branch("z_m",&z_m, "z_m/F");
   _ttree->Branch("z_pt",&z_pt, "z_pt/F");
   _ttree->Branch("z_eta",&z_eta, "z_eta/F");
   _ttree->Branch("lep_z_dphi",&lep_z_dphi, "lep_z_dphi/F");
   _ttree->Branch("mu_mu_dphi",&mu_mu_dphi, "mu_mu_dphi/F");
   _ttree->Branch("MET",&MET, "MET/F");
   _ttree->Branch("averageMu",&averageMu,"averageMu/F");

   // _ttree->Branch("dphi", &dphi, "dphi/F");
   _ttree->Branch("Mt",&Mt,"Mt/F");
   // _ttree->Branch("lep_ratioet20", &lep_ratioet20, "lep_ratioet20/F");
   // _ttree->Branch("lep_ratioet30", &lep_ratioet30, "lep_ratioet30/F");
   // _ttree->Branch("lep_ratioet40", &lep_ratioet40, "lep_ratioet40/F");

   
   _ltree->Branch("lep_pt", &lep_pt, "lep_pt/F");
   _ltree->Branch("lep_eta", &lep_eta, "lep_eta/F");
   _ltree->Branch("lep_SF",&lep_SF,"lep_SF/F");
   _ltree->Branch("z_m",&z_m, "z_m/F");
   _ltree->Branch("z_pt",&z_pt, "z_pt/F");
   _ltree->Branch("z_eta",&z_eta, "z_eta/F");
   _ltree->Branch("lep_z_dphi",&lep_z_dphi, "lep_z_dphi/F");
   _ltree->Branch("mu_mu_dphi",&mu_mu_dphi, "mu_mu_dphi/F");
   _ltree->Branch("MET",&MET, "MET/F");
   _ltree->Branch("averageMu",&averageMu,"averageMu/F");
   // _ltree->Branch("dphi", &dphi, "dphi/F");
   _ltree->Branch("Mt",&Mt,"Mt/F");
   // _ltree->Branch("lep_ratioet20", &lep_ratioet20, "lep_ratioet20/F");
   // _ltree->Branch("lep_ratioet30", &lep_ratioet30, "lep_ratioet30/F");
   // _ltree->Branch("lep_ratioet40", &lep_ratioet40, "lep_ratioet40/F");

  
   // _mtree->Branch("lep_pt", &lep_pt, "lep_pt/F");
   // _mtree->Branch("lep_eta", &lep_eta, "lep_eta/F");
   // _mtree->Branch("lep_SF",&lep_SF,"lep_SF/F");
   // _mtree->Branch("z_m",&z_m, "z_m/F");
   // _mtree->Branch("z_pt",&z_pt, "z_pt/F");
   // _mtree->Branch("z_eta",&z_eta, "z_eta/F");
   // _mtree->Branch("lep_z_dphi",&lep_z_dphi, "lep_z_dphi/F");
   // _mtree->Branch("mu_mu_dphi",&mu_mu_dphi, "mu_mu_dphi/F");
   // _mtree->Branch("MET",&MET, "MET/F");
   // _mtree->Branch("averageMu",&averageMu,"averageMu/F");
   // _mtree->Branch("dphi", &dphi, "dphi/F");
   //_mtree->Branch("Mt",&Mt,"Mt/F");
   // _mtree->Branch("lep_ratioet20", &lep_ratioet20, "lep_ratioet20/F");
   // _mtree->Branch("lep_ratioet30", &lep_ratioet30, "lep_ratioet30/F");
   // _mtree->Branch("lep_ratioet40", &lep_ratioet40, "lep_ratioet40/F");

   tight = 0;
   loose = 0;
   medium = 0;
   eee=eemu=mumue=mumumu = 0;


   //_trigs = new TH1F("trigs", "n pass trigger", 15, 0, 15);

}

void analysis::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t analysis::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either analysis::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.
   fChain->GetTree()->GetEntry(entry);
   if(entry%1000 == 0) cout<<"Processing event : "<<entry<<endl;
   _hcutflow2->Fill("ntuple",1);


   //trigger trig[3] HLT_2mu10 
   //trigger trig[0] HLT_mu20_iloose_L1MU15
   if (!(sig_trigCode & (1<<3))) return kTRUE;
   _hcutflow2->Fill("passtrigger",1);


   //grl
   if(!evt_cuts) return kTRUE;
   _hcutflow2->Fill("passgrl",1);


   //select electron/muon channel
   int count = 0;
   int el = -1;
   int mu[2] = {-1,-1};
   for ( int i=0 ; i<leps_ ; i++)
   {
    
    if(TMath::Abs(leps_ID[i]/1000)<12.) el = i;
    else
    {
      mu[count] = i;
      count++;
    }

   }
   //if(count == 0) eee++;
   //if(count == 1) eemu++;
   //if(count == 2) mumue++;
   //if(count == 3) mumumu++;
   //if(count > 3 ) cout<<"interesting"<<endl;

   //select mumue channel
   if(!(count==2)) return kTRUE;
   _hcutflow2->Fill("selectmumue",1);

   if( (TMath::Abs(leps_eta[el]) < 1.52) && (TMath::Abs(leps_eta[el]) > 1.37)) return kTRUE;
   //if(1.37<TMath::Abs(leps_eta[el])<1.52) return kTRUE;
   _hcutflow2->Fill("removecrackregion",1);


   // decode quality
   bool IsTightEl, IsLooseEl;
   IsLooseEl = true;
   IsTightEl = false;
   if(leps_isTight[el]) IsTightEl = true;

   


   //calculation
   TLorentzVector mu1, mu2, z, electron, vl, w;
   mu1.SetPtEtaPhiM(leps_pt[mu[0]],leps_eta[mu[0]],leps_phi[mu[0]],0.105658);
   mu2.SetPtEtaPhiM(leps_pt[mu[1]],leps_eta[mu[1]],leps_phi[mu[1]],0.105658);
   z = mu1 + mu2;


   electron.SetPtEtaPhiM(leps_pt[el],leps_eta[el],leps_phi[el],0.0005);
   vl.SetPxPyPzE(sig_MetX,sig_MetY,0,sig_MetRel);
   w = electron + vl;
   Mt = TMath::Sqrt((electron.Pt()+vl.Pt())*(electron.Pt()+vl.Pt())-(w.Pt())*(w.Pt())); 

   lep_z_dphi = kinem::delta_phi(leps_phi[el], z.Phi());
   mu_mu_dphi = kinem::delta_phi(mu1.Phi(), mu2.Phi());


   

   // //Cuts
   bool wzselection = true;
   if(wzselection)
   {
      //Mt requirement
      if(Mt < 40) return kTRUE;
      _hcutflow2->Fill("mt>40",1);
      //z mass
      if(z.M()< 80) return kTRUE;
      if(z.M() > 100) return kTRUE;
      _hcutflow2->Fill("z mass",1);
      //select w events for testing
      if(sig_MetRel < 25) return kTRUE;
      _hcutflow2->Fill("met>25GeV",1);

   }

   bool fakeratecal = false;
   if(fakeratecal)
   {
      //back to back requirement for sig_lep and sig_jet
      if(!(kinem::delta_phi(leps_phi[el], z.Phi()) > 2.6)) return kTRUE;
      _hcutflow2->Fill("pass backtoback",1);
      //z mass
      if(z.M()< 80) return kTRUE;
      if(z.M() > 100) return kTRUE;
      _hcutflow2->Fill("z mass",1);
      //met
      if(sig_MetRel >  25) return kTRUE;
      _hcutflow2->Fill("met<=25GeV",1);


   }

   
   
   //Write Tree
   z_pt = z.Pt();
   z_eta = z.Eta();
   z_m= z.M();

   lep_pt = leps_pt[el];
   lep_eta = leps_eta[el];
   averageMu = evt_averageMu;

   //lep_ratioet20 = leps_topoetcone20[0] / lep_pt;
   //lep_ratioet30 = leps_topoetcone30[0] / lep_pt;
   //lep_ratioet40 = leps_topoetcone40[0] / lep_pt;
   MET = sig_MetRel;

   if(IsTightEl) 
   {
     lep_SF = evt_ElSF * evt_weight;
     _ttree->Fill();
     tight++;
   }

   if(IsLooseEl) 
   {
     lep_SF = evt_ElSF * evt_weight;
     _ltree->Fill();
     loose++;
   }

   return kTRUE;
}

void analysis::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void analysis::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.
   _ttree->Write();
   _ltree->Write();
   _hcutflow2->Write();
   //_file->Write("All");
   //_file->Close();
   cout<<"tight:"<<tight<<" ;loose:"<<loose<<"   ;medium"<<medium<<endl;
   //cout<<"eee:"<<eee<<";  eemu"<<eemu<<";  mumue"<<mumue<<";mumumu"<<mumumu<<endl;
}
