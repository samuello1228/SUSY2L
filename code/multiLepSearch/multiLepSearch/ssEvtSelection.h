#ifndef multiLepSearch_ssEvtSelection_H
#define multiLepSearch_ssEvtSelection_H

#include <EventLoop/Algorithm.h>
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
// #include "ElectronPhotonSelectorTools/AsgElectronLikelihoodTool.h"
// #include "ElectronIsolationSelection/IsolationSelectionTool.h"
#include "IsolationSelection/IsolationSelectionTool.h"
// #include "MuonSelectorTools/MuonSelectionTool.h"

// #include "CPAnalysisExamples/errorcheck.h"
#include "SUSYTools/SUSYObjDef_xAOD.h"
#include "SUSYTools/SUSYCrossSection.h"
#include "susyEvts.h"

//trigger
#include <TrigConfxAOD/xAODConfigTool.h>
#include <TrigDecisionTool/TrigDecisionTool.h>

//pileup reweighting
#include "PileupReweighting/PileupReweightingTool.h"
#include "AsgTools/ToolHandle.h"

//systematics
// #include "PATInterfaces/SystematicRegistry.h"

// MCTruthClassifier
#include "MCTruthClassifier/MCTruthClassifier.h"

class GoodRunsListSelectionTool;
class TH1;

class ChargeFlipBkgTool;
class FakeLepBkgTool;

struct TRIGCONF{
   uint32_t runStart;
   uint32_t runEnd;
   std::vector<std::string> eeTrig;
   std::vector<std::string> emTrig;
   std::vector<std::string> mmTrig;
   unsigned long int ee_mask;
   unsigned long int em_mask;
   unsigned long int mm_mask;
};

class ssEvtSelection : public EL::Algorithm
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
  std::string CF_outputName;
  std::string CF_outputTreeName;
  std::string CF_derivationName;
  std::string study;
  std::string mcTruthMatch;
  bool useChargeIDSelector;
  int doSys;
  //int CF_nLepCutExactly;
  //int CF_nLepCutMin;
  //int CF_nJetCutExactly;
  //int CF_nJetCutMin;

  //float CF_ElPtCut;
  //float CF_Eld0SigCut;
  //float CF_Elz0Cut;

  //float CF_MuPtCut;
  //float CF_Mud0SigCut;
  //float CF_Muz0Cut;

  //float CF_JetPtCut;
  //float CF_JetEtaCut;
  //float CF_JetJvtCut;
  int CF_vxTrkNMin;
  float CF_vxTrkPtMin;
  float CF_mT2_m0;

  int CF_isMC;
  std::vector< std::string > CF_grlFiles;
  std::string CF_ConfigFile;
  std::vector< std::string > CF_trigNames;

  /// PRW files
  std::vector< std::string > CF_PRW_confFiles;
  std::vector< std::string > CF_PRW_lcalcFiles;

  // MCTruthClassifier
  MCTruthClassifier* m_truthClassifier; //!


  double ECIDS_OP; //!
  std::string ECIDS_trainingFile; //!

  // this is a standard constructor
  ssEvtSelection (std::string name="ss2lSelection");

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
  EL::StatusCode fillLepton(xAOD::IParticle* p, L_PAR& l, unsigned int i);
  EL::StatusCode fillLepton(xAOD::Electron* p, L_PAR& l, unsigned int i);
  EL::StatusCode fillLepton(xAOD::Muon* p, L_PAR& l, unsigned int i);
  void fillLeptonCommon(xAOD::IParticle* p, L_PAR& l);
  int addTruthPar(const xAOD::TruthParticle* p, TRUTHS& v, int saveParent=0);
  bool isSS(xAOD::IParticle* p1, xAOD::IParticle* p2){
    static SG::AuxElement::Accessor<float> charge("charge");
    return charge(*p1)*charge(*p2)>0;
  }
  void setupTriggers();
  TRIGCONF* getTriggerConf(uint32_t run);

 protected:
  std::string m_name; //!
  // Tree *myTree; //!
  TH1 *m_hCutFlow; //!
  TH1 *m_hCutFlowNominal; //!
  TH1 *m_hCutFlowDummy; //!
  TH1 *m_hTrigs; //!
//   xAOD::TEvent *m_event;  //!
//   xAOD::TStore *m_store;  //!
  std::vector<ST::SystInfo> m_systInfoList; //!
  std::vector<susyEvts*   > m_susyEvtList; //!
  susyEvts* m_susyEvt; //!
  ST::SUSYObjDef_xAOD* m_objTool; //!
//   CP::MuonSelectionTool* m_muonSelTool; //!

  ChargeFlipBkgTool* mChargeFlipBkgTool; //!
  FakeLepBkgTool*    mFakeLepBkgTool;  //!

  //////// trigger
  //TrigConf::xAODConfigTool *configTool;//!
  //ToolHandle<TrigConf::ITrigConfigTool> *configHandle;//!
  //Trig::TrigDecisionTool *trigDecTool;//!
  bool start;//!

  //static const int nTrig = 15;
  //static int trigCount[nTrig];
  //static const string trigNames[nTrig];
  ///////////////////////////////////////

//   std::vector<susyEvts*> m_susyEvtList; //!

  int m_eventCounter; //!
//   AsgElectronLikelihoodTool* m_LHToolTight2015; //!
//   AsgElectronLikelihoodTool* m_LHToolMedium2015; //!
//   AsgElectronLikelihoodTool* m_LHToolLoose2015; //!
  GoodRunsListSelectionTool *m_grl; //!
  //std::map<OFLAGS, CP::IsolationSelectionTool* > m_isoTools; //!
  CP::IsolationSelectionTool* m_isoTool; //!
  TH1* mh_ElChargeFlip; //!
  std::vector< TRIGCONF* > m_trigSel; //!
  TRIGCONF* m_nowTrigSel{0}; //!
  std::string m_em_eKey;
  std::string m_em_mKey;
  std::string m_ee_Key;

 private:
  // this is needed to distribute the algorithm to the workers
  ClassDef(ssEvtSelection, 1);
};

#endif
