#include <../multiLepSearch/obj_def.h>
#include <iostream>
#include <vector>
#include <TChain.h>
#include <TInterpreter.h>

using namespace std;
int test2(){

  gInterpreter->GenerateDictionary("vector<L_PAR>","../multiLepSearch/obj_def.h;vector"); 
  TChain* ch = new TChain("evt2l");
  ch->Add("/eos/atlas/user/c/clo/ntuple/AnalysisBase-02-04-31-ccd99030/mc15_13TeV.993821.MGPy8EG_A14N13LO_C1N2_Wh_2L_165_35.merge.DAOD_SUSY2.e5678_a766_a821_r7676_p2949_p2972.root");
  ch->Show(0);
  
  auto vleps = new vector< L_PAR >();
  ch->SetBranchAddress("leps", &vleps);

  int Z; 
//   ch->SetBranchAddress("leps", &Z);

  ch->GetEntry(0);
  cout << "the size of the leptons is " << vleps->size() << " " << Z << endl;
  cout << "the lepton 0 pt= " << vleps->at(0).pt << endl;

  cout << "testing " << endl;
  return 0;
}

