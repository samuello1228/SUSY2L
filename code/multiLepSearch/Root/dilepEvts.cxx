#include "multiLepSearch/dilepEvts.h"

dilepEvts::dilepEvts(TString name, TFile* outFile)
{
  this->name = name;
  this->outFile = outFile;
}

int dilepEvts::fill(TString streamName, bool saveSysWeights){
  auto result = mTrees.find(streamName);
  TTree* outTree = (result==mTrees.end())? getNewTree(streamName, saveSysWeights) : result->second;
  return outTree->Fill();
}

void dilepEvts::setWeight(TString wName, float val){
  mWeights[wName] = val;
}

TTree* dilepEvts::getNewTree(TString streamName, bool saveSysWeights){
  TTree* myTree = new TTree(streamName, "a dilep NTUP tree");
  myTree->Branch("evtInfo" , &evtInfo);
  myTree->Branch("l12" , &l12);
  myTree->Branch("sig" , &sig);
  myTree->Branch("leps", &leps);
  myTree->Branch("jets", &jets);

  if (saveSysWeights){
    for (auto& aPair : mWeights){ myTree->Branch(aPair.first, &aPair.second);}
  }

  if (outFile) myTree->SetDirectory(outFile);
  mTrees[streamName] = myTree;
  return myTree;
}
