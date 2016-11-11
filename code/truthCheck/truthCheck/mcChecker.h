#ifndef multiLepSearch_mcChecker_H
#define multiLepSearch_mcChecker_H

#include <EventLoop/Algorithm.h>
#include <TTree.h>
#include <xAODTruth/TruthParticle.h>
#include <vector>

struct kPar{
  float pt;
  float eta;
  float phi;
  int id;
};
const std::string kPar_s = "pt/F:eta:phi:id/I";

struct mPar:kPar{
  float mll;
};
const std::string mPar_s = kPar_s+":mll/F";

struct xPar:kPar{
  float MET_dPhi;
};
const std::string xPar_s = kPar_s+":MET_dPhi/F";

class mcChecker : public EL::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  // float cutValue;



  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:
  // Tree *myTree; //!
  // TH1 *myHist; //!
  std::string outputName;

  // this is a standard constructor
  mcChecker ();

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();
  void setKinPar(kPar& p, const xAOD::TruthParticle* t=nullptr);
  void setKinPar(kPar& p, TLorentzVector& t, int pdg=-999);
  void setKinPar(mPar& p, const xAOD::TruthParticle* t=nullptr);
  void setKinPar(mPar& p, TLorentzVector& t, int pdg=-999);
  void showDecay(const xAOD::TruthParticle* t, std::string x="");
  size_t getDaughters(const xAOD::TruthParticle* t, std::vector< const xAOD::TruthParticle* >& ds);

private:
  TTree *tree; //!
  kPar pC1; //!
  kPar pN2; //!
  kPar pN1w; //!
  kPar pN1z; //!
  mPar pW; //!
  mPar pZ; //!
  xPar pl1; //!
  xPar pl2; //!
  xPar pl3; //!
  xPar pl4; //!
  float mSC; //!
  float mSS; //!
  float MET; //!
  float MET_phi; //!
  float SS_dPhi; //!
  float SS_dR; //!
  std::vector<kPar> pOther; //!
  std::vector<xPar> pJets; //!

public:
  // this is needed to distribute the algorithm to the workers
  ClassDef(mcChecker, 1);
};

#endif
