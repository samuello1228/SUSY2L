#include "../multiLepSearch/susyEvts.h"
#include <TChain.h>
#include <iostream>
using namespace std;

int main()
{
  cout << "testing " << endl;
  TChain* ch = new TChain("evt2l");
  ch->Add("~/links/SAMPLES/DC14/evt2l/user.clo.test1.204549_Herwigpp_UEEE4_CTEQ6L1_sM_wC_WW_C1_240_N1_0_myOutput.root.29953328/user.clo.5579861._000001.myOutput.root");
  susyEvts* mEvts = new susyEvts(ch);

  for(long int i=0; i<ch->GetEntries(); i++){
    mEvts->GetEntry(i);
    cout << mEvts->sig.mT2 << endl;
    auto kp=mEvts->vleps;
    cout << kp->size() << endl;
    for(auto x: *kp){cout << x.pt << endl;}
  }

  return 0;
}
