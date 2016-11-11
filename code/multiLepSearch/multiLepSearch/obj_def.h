#ifndef OBJECT_H
#define OBJECT_H

#include <string>
#include <vector>
using std::string;

struct PAR0{
  float pt;
  float eta;
  float phi;
};
const string PAR0_s = "pt/F:eta/F:phi/F";

struct PAR:PAR0{
  float MET_dPhi;
};
const string PAR_s = PAR0_s+":MET_dPhi/F";

struct R_PAR:PAR{
  float m;
  float dPhi;
  float dR;
};
const string R_PAR_s = PAR_s+":m/F:dPhi:dR";

struct L_PAR:PAR{
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
  float z0sinTheta;
  float d0sig;
  float truthProb;
  float SF_Loose;
  float SF_Medium;
  float SF_Tight;
  float wt1;
  float wt2;
  float wt3;
  int ID;
  int author;
  int truthI;
  int truthType;
  int truthOrig;
  int isTight;
  unsigned int Q;
  unsigned int lFlag;
  unsigned int nBHits;
  unsigned int nPixHits;
  unsigned int nSCTHits;
  unsigned int nPixHoles;
  unsigned int nSCTHoles;
  unsigned int nTRTHits;
  unsigned int nTRTOutliers;
};
const string L_PAR_s = PAR_s+":topoetcone20/F:topoetcone30:topoetcone40:ptcone20:ptcone30:ptcone40:mT:d0:d0Err:z0:z0Err:z0sinTheta:d0sig:truthProb:SF_Loose:SF_Medium:SF_Tight:wt1:wt2:wt3:ID/I:isTight/I:author/I:truthI/I:truthType:truthOrig:Q/i:lFlag/i:nBHits/i:nPixHits:nSCTHits:nPixHoles:nSCTHoles:nTRTHits:nTRTOutliers";

struct EL_Par:L_PAR{
  int elID;
};
const string EL_PAR_s = L_PAR_s+":elID/i";

struct MU_Par:L_PAR{
  int muID;
};
const string MU_PAR_s = L_PAR_s+":muID/i";

struct J_PAR:PAR{
  unsigned int jFlag ;
  float e;
};
const string J_PAR_s = PAR_s+":jFlag/I:e/F";

struct EVT{
  unsigned long int run;
  unsigned long int event;
  unsigned long int lumiBlock;
  unsigned long int actualMu;
  unsigned int index;
  unsigned int cuts;
  unsigned int trig;
  unsigned int isMC;
  int flag;
  int pass;
  float averageMu;
  float Xsec;
  float weight;
  float pwt;
  float ElSF;
  float MuSF;
  float BtagSF;
  float qFwt;
  float fLwt;
 };
const string EVT_s = "run/l:event/l:lumiBlock/l:actualMu/l:index/i:cuts/i:trig/i:isMC/i:flag/I:pass/I:averageMu/F:Xsec/F:weight/F:pwt/F:ElSF/F:MuSF/F:BtagSF/F:qFwt/F:fLwt/F";
struct SIGNATURE{
  unsigned long int trigCode;//trigger info
  float Met;
  float MetRel;
  float MetX;
  float MetY;
  float mT2;
  float HT;
  unsigned nEl;
  unsigned nMu;
  unsigned nTau;
  unsigned nJet;
  unsigned nPV;
  unsigned nVtx;
};
const string SIGNATURE_s = "trigCode/l:Met/F:MetRel/F:MetX/F:MetY/F:mT2/F:HT:nEl/i:nMu/i:nTau/i:nJet/i:nPV/i:nVtx/i";

struct TR_PAR:PAR0{
  int pdgId; 
  int barcode;
  int motherI;
  int matchI;
  int particleType;
  int particleOrigin;
};
const string TR_PAR_s = PAR0_s+"pdgId/I:barcode/I:motherI/I:matchI/I:particleType/I:particleOrigin/I";

enum FLAGS{
  PASS_GRL = 1<<0,
};

enum OFLAGS{
  IS_BASELINE = 1<<0, //1
  IS_SIGNAL = 1<<1, //2
  IS_PASSOR = 1<<2, //4
  IS_BAD = 1<<3, //8
  IS_LOOSEBASELINE = 1<<4, //16 
  ISO_LOOSE = 1<<8, /// isolation
  ISO_TIGHT = 1<<9,
  ISO_GRAD = 1<<10,
};

enum EFLAGS{
  EL_OQ = 1<<4, //16
};

enum MFLAGS{
  MU_COSMIC = 1<<4, //16
};

enum JFALGS{
  JT_BJET_LOOSE = 1<<4, //16
  JT_BJET = 1<<5, //32
};

typedef std::vector< L_PAR > LEPTONS;
typedef std::vector< PAR > PARTICLES;
typedef std::vector< J_PAR > JETS;
typedef std::vector< TR_PAR > TRUTHS;

#endif
