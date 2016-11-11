#ifndef multiLepSearch_dilepSelection_H
#define multiLepSearch_dilepSelection_H

#include <EventLoop/Algorithm.h>
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "ElectronPhotonSelectorTools/AsgElectronLikelihoodTool.h"
// #include "ElectronIsolationSelection/IsolationSelectionTool.h"
#include "IsolationSelection/IsolationSelectionTool.h"
#include "MuonSelectorTools/MuonSelectionTool.h"

#include "xAODMissingET/MissingETContainer.h"
#include "xAODMissingET/MissingETAuxContainer.h"

// #include "CPAnalysisExamples/errorcheck.h"
#include "SUSYTools/SUSYObjDef_xAOD.h"
#include "dilepEvts.h"

//trigger
#include <TrigConfxAOD/xAODConfigTool.h>
#include <TrigDecisionTool/TrigDecisionTool.h>

//pileup reweighting
#include "PileupReweighting/PileupReweightingTool.h"
#include "AsgTools/ToolHandle.h"

//bkg weight tools
#include "multiLepSearch/IChargeFlipBkgTool.h"
#include "multiLepSearch/ChargeFlipBkgTool.h"
#include "multiLepSearch/IFakeLepBkgTool.h"
#include "multiLepSearch/FakeLepBkgTool.h"

//systematics
#include "PATInterfaces/SystematicRegistry.h"
///////////////////////

class GoodRunsListSelectionTool;

class dilepSelection : public EL::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  // float cutValue;



  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:
  std::string outputName;
  std::string outputTreeName;
  //int nLepCutExactly;
  //int nLepCutMin;
  //int nJetCutExactly;
  //int nJetCutMin;

  int m_dataType;
  int m_doSys;
  std::string m_grlFile;
  std::string m_susyToolCfgFile;
  std::vector< std::string > trigNames;
  std::vector< std::string > PRW_confFiles;
  std::vector< std::string > PRW_lcalcFiles;

  // this is a standard constructor
  dilepSelection (std::string name="dilepSelection");

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
  void fillLepton(xAOD::IParticle* p, ntupLep& l);

private:
  std::string m_name; //!

  TH1 *m_hCutFlow; //!
  TH1 *m_hTrigs; //!

  xAOD::TEvent *m_event;  //!
  xAOD::TStore *m_store;  //!

  ST::SUSYObjDef_xAOD* m_objTool; //!
  GoodRunsListSelectionTool* m_grl; //!
  IChargeFlipBkgTool* mChargeFlipBkgTool; //!
  IFakeLepBkgTool*    mFakeLepBkgTool; //!

  bool firstLoop;//!
  int m_eventCounter; //!

  int m_vxTrkNMin; //!
  float m_vxTrkPtMin; //!

  xAOD::ElectronContainer* electrons_copy; //!
  xAOD::ShallowAuxContainer* electrons_copyaux; //!

  xAOD::MuonContainer* muons_copy; //!
  xAOD::ShallowAuxContainer* muons_copyaux; //!

  xAOD::JetContainer* jets_copy; //!
  xAOD::ShallowAuxContainer* jets_copyaux; //!

  xAOD::TauJetContainer* taus_copy; //!
  xAOD::ShallowAuxContainer* taus_copyaux; //!

  xAOD::MissingETContainer* metcst; //!
  xAOD::MissingETAuxContainer* metcst_aux; //!

  dilepEvts* m_dilepEvt; //!
  std::vector<CP::SystematicSet> m_sysList; //!

  // this is needed to distribute the algorithm to the workers
  ClassDef(dilepSelection, 1);

};

#endif
