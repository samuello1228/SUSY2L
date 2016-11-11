#ifndef dilepEvts_H
#define dilepEvts_H
#include <TTree.h>
#include <TChain.h>
#include <TFile.h>
#include <TEntryList.h>
#include "multiLepSearch/dilep_objDef.h"
// #include "obj_def.h"
#include <string>
#include <map>
class dilepEvts{

 public:
   dilepEvts(TString name="dilep", TFile* outFile=NULL);
   ~dilepEvts();

   int fill(TString streamName = "evt2l", bool saveSysWeights = false);
   void setWeight(TString wName, float val);

   ntupEvtInfo evtInfo;
   ntupSig sig;
   ntupP12 l12;
   std::vector<ntupLep> leps;
   std::vector<ntupJet> jets;

 private:
   TTree* getNewTree(TString streamName, bool saveSysWeights);

   TString name;
   TFile* outFile;
   std::map<TString, TTree*> mTrees;
   std::map<TString, float> mWeights;

};
#endif //dilepEvts_H
