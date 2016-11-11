#ifndef multiLepSearch_ljetSelection_H
#define multiLepSearch_ljetSelection_H

#include <EventLoop/Algorithm.h>
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "ElectronPhotonSelectorTools/AsgElectronLikelihoodTool.h"
// #include "ElectronIsolationSelection/IsolationSelectionTool.h"
#include "IsolationSelection/IsolationSelectionTool.h"
#include "MuonSelectorTools/MuonSelectionTool.h"

#include "CPAnalysisExamples/errorcheck.h"
#include "SUSYTools/SUSYObjDef_xAOD.h"
#include "ljetEvts.h"

//trigger
#include <TrigConfxAOD/xAODConfigTool.h>
#include <TrigDecisionTool/TrigDecisionTool.h>

//pileup reweighting
#include "PileupReweighting/PileupReweightingTool.h"
#include "AsgTools/ToolHandle.h"

// Electron ID scale factors
#include "ElectronEfficiencyCorrection/AsgElectronEfficiencyCorrectionTool.h"

///////////////////////

class GoodRunsListSelectionTool;

class ljetSelection : public EL::Algorithm
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
  int nLepCutExactly;
  int nLepCutMin;
  int nJetCutExactly;
  int nJetCutMin;

  float ElPtCut;
  float Eld0SigCut;
  float Elz0Cut;

  float MuPtCut;
  float Mud0SigCut;
  float Muz0Cut;

  float JetPtCut;
  float JetEtaCut;
  float JetJvtCut;

  int m_isMC;
  std::string m_grlFile;
  std::vector< std::string > trigNames;

  // this is a standard constructor
  ljetSelection (std::string name="ljetSelection");

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
  void fillLepton(xAOD::IParticle* p, L_PAR& l, unsigned int i);
  int addTruthPar(const xAOD::TruthParticle* p, TRUTHS& v, int saveParent=0);

private:
  // Tree *myTree; //!
  TH1 *m_hCutFlow; //!
  TH1 *m_hTrigs; //!
  xAOD::TEvent *m_event;  //!
  xAOD::TStore *m_store;  //!
  std::string m_name; //!
  ST::SUSYObjDef_xAOD* m_objTool; //!
  CP::MuonSelectionTool* m_muonSelTool; //!

//////// trigger
  TrigConf::xAODConfigTool *configTool;//!
  ToolHandle<TrigConf::ITrigConfigTool> *configHandle;//!
  Trig::TrigDecisionTool *trigDecTool;//!
  bool start;//!

//   static const int nTrig = 15;
//   static int trigCount[nTrig];
//   static const string trigNames[nTrig];
///////////////////////////////////////

  ljetEvts* m_susyEvt; //!
  int m_eventCounter; //!
  int m_vxTrkNMin; //!
  float m_vxTrkPtMin; //!
  AsgElectronLikelihoodTool* m_LHToolTight2015; //!
  AsgElectronLikelihoodTool* m_LHToolMedium2015; //!
  AsgElectronLikelihoodTool* m_LHToolLoose2015; //!
  GoodRunsListSelectionTool *m_grl; //!
//   std::map<OFLAGS, CP::IsolationSelectionTool* > m_isoTools; //!
  CP::IsolationSelectionTool* m_isoTool; //!

  // electron ID SF tools
  AsgElectronEfficiencyCorrectionTool* m_elecEfficiencySFTool_loose; //!
  AsgElectronEfficiencyCorrectionTool* m_elecEfficiencySFTool_tight; //!

  // this is needed to distribute the algorithm to the workers
  ClassDef(ljetSelection, 1);

  //pileup reweighting
  ToolHandle<CP::IPileupReweightingTool> m_prwTool;
};

#endif
