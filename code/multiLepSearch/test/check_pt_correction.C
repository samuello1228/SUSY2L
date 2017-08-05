#include "../multiLepSearch/susyEvts.h"
#include <TChain.h>
#include <iostream>
#include <memory>
#include <TLorentzVector.h>
using namespace std;


std::unique_ptr<TLorentzVector> getLepton(const L_PAR& l){
  std::unique_ptr<TLorentzVector> ptr(new TLorentzVector());
  float m = abs(l.ID)==11000?0.511:105.6583745;

  ptr->SetPtEtaPhiM(l.pt,l.eta,l.phi,m*0.001);

  return ptr;
}



int main()
{
  cout << "testing " << endl;
  TChain* ch = new TChain("evt2l");
  ch->Add("~/links/SAMPLES/R20/susyNtuple/user.dzhang.v10.5.MCSig_myOutput.root/user.dzhang.mc15_13TeV.393873.MGPy8EG_A14N23LO_C1N2_Wh_hall_325p0_75p0_2L7.11771539._000035.myOutput.root");
  susyEvts* mEvts = new susyEvts(ch);

  for(long int i=0; i<ch->GetEntries(); i++){
    mEvts->GetEntry(i);
    auto kp=mEvts->vleps;

    if(abs((*kp)[0].ID)!=11000) continue;
    if(abs((*kp)[1].ID)!=11000) continue;
    if((*kp)[1].ID * (*kp)[1].ID < 0) continue;


    auto l1 = getLepton((*kp)[0]);
    auto l2 = getLepton((*kp)[1]);
    auto l12 = *l1 + *l2;

    cout << l1->Pt() << " " << l2->Pt() << " " << l12.M() << " - " << mEvts->l12.m << " = " << l12.M() - mEvts->l12.m << endl;

    cout << kp->size() << endl;
    for(auto x: *kp){cout << x.pt << endl;}
  }

  return 0;
}
