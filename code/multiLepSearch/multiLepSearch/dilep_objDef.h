#ifndef dilep_objDef_H
#define dilep_objDef_H

#include <string>
#include <vector>
using std::string;

class ntupP12{
public:
  float m;
  float pt;
  float eta;
  float phi;
  float dPhi;
  float dR;
  float MET_dPhi;
};

class ntupLep{
public:
  float pt;
  float eta;
  float phi;
  float topoetcone20;
  float topoetcone30;
  float topoetcone40;
  float ptcone20;
  float ptcone30;
  float ptcone40;
  float mT;
  float d0;
  float d0Err;
  float z0;
  float z0Err;
  float truthProb;
  float MET_dPhi;
  int ID;
  int author;
  bool isBad;
  bool isCosmic;
};

class ntupJet{
public:
  float pt;
  float eta;
  float phi;
  float MET_dPhi;
  bool passOR;     
  bool isBaseline;
  bool isBad;      
  bool isSig;      
  bool isBJetLoose;
  bool isBJet;     

};

class ntupSig{
public:
  unsigned long int trigCode;//trigger info
  float Met;
  float MetRel;
  float MetX;
  float MetY;
  float mT2;
  unsigned nEl;
  unsigned nMu;
  unsigned nTau;
  unsigned nJet;
  unsigned nPV;
  unsigned nVtx;
  unsigned nCosmic;
  unsigned nBJet;
  float ElSF;
  float MuSF;
  float FakeLepWeight0;
  float CFlipWeight0;
};

class ntupEvtInfo{
public:
  unsigned long int run;
  unsigned long int event;
  unsigned long int lumiBlock;
  unsigned long int actualMu;
  //unsigned int index;
  //unsigned int cuts;
  unsigned int trig;
  //int flag;
  //int pass;
  float averageMu;
  float weight;
  bool passGRL;
};
#endif
