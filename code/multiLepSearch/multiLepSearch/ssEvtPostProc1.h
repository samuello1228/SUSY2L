#ifndef multiLepSearch_ssEvtPostProc1_H
#define multiLepSearch_ssEvtPostProc1_H

#include "TObject.h"
#include "SUSYTools/SUSYCrossSection.h"
#include "susyEvts.h"

#include "AsgTools/ToolHandle.h"
#include <vector>

class ChargeFlipBkgTool;
class FakeLepBkgTool;

class ssEvtPostProc1 : public TObject
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.

  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:
  // Tree *myTree; //!
  // TH1 *myHist; //!
  std::string inputFile;
  std::string norminalTreeName;
  bool isData;

  ChargeFlipBkgTool* mChargeFlipBkgTool;
  FakeLepBkgTool*    mFakeLepBkgTool;

  // this is a standard constructor
  ssEvtPostProc1 (std::string name="ss2lPostProc1");
  ~ssEvtPostProc1();

  int execute();

 protected:
  string m_name;//!

  SUSY::CrossSectionDB* m_XsecDB;  //!
  TFile* inF; //!
  TH1* mh_ElChargeFlip; //!
  bool isFirstInit; //!


  std::vector<susyEvts*> rawSusyEvtsList                  ; //!
  std::vector<susyEvts*> outSusyEvtsList                  ; //!
  std::vector<susyEvts*> outSusyEvtsList_CFLIP_SYS__1up   ; //!
  std::vector<susyEvts*> outSusyEvtsList_CFLIP_SYS__1dn   ; //!
  std::vector<susyEvts*> outSusyEvtsList_FAKE_LEP_E_SYS__1up; //!
  std::vector<susyEvts*> outSusyEvtsList_FAKE_LEP_E_SYS__1dn; //!
  std::vector<susyEvts*> outSusyEvtsList_FAKE_LEP_U_SYS__1up; //!
  std::vector<susyEvts*> outSusyEvtsList_FAKE_LEP_U_SYS__1dn; //!






  int initialize ();
  int runLoop();
  int finalize ();

  int recalPtRelatedVar (susyEvts* inTree);

 private:

  // this is needed to distribute the algorithm to the workers
  ClassDef(ssEvtPostProc1, 1);
};

#endif
