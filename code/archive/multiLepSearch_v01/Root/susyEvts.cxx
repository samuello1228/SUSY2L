#include "multiLepSearch/susyEvts.h"

TTree* susyEvts::makeTree(TString treename)
{
  tree2 = new TTree(treename, "a angles tree");
  tree2->Branch("evt", &evt, EVT_s.c_str());
  tree2->Branch("leps", &vleps);
  tree2->Branch("l12", &l12, R_PAR_s.c_str());
  tree2->Branch("jets", &vjets);
  tree2->Branch("truths", &vtruths);
  tree2->Branch("sig", &sig, SIGNATURE_s.c_str());
  return tree2;
}

void susyEvts::getTree(TChain* tr)
{
  tree1 = tr;
  tree1->SetBranchAddress("evt", &evt);
  tree1->SetBranchAddress("leps", &vleps);
  tree1->SetBranchAddress("l12", &l12);
  tree1->SetBranchAddress("jets", &vjets);
  tree1->SetBranchAddress("truths", &vtruths);
  tree1->SetBranchAddress("sig", &sig);
}

Int_t susyEvts::Next(){
  if(m_el == m_elist->GetN()) return -1;
  Int_t treenum=0;
  Long64_t treeEntry = m_elist->GetEntryAndTree(m_el,treenum);
  Long64_t chainEntry = treeEntry+tree1->GetTreeOffset()[treenum];
  m_el++;
  return tree1->GetEntry(chainEntry);
}
