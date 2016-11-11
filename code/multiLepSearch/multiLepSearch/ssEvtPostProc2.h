#ifndef multiLepSearch_ssEvtPostProc2_H
#define multiLepSearch_ssEvtPostProc2_H

#include "TObject.h"
#include "SUSYTools/SUSYCrossSection.h"
#include "susyEvts.h"

#include "AsgTools/ToolHandle.h"
#include <vector>

#include "TMVA/Reader.h"

using namespace std;

class ssEvtPostProc2 : public TObject
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.

  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:
  std::string inputFile;
  std::string norminalTreeName;
  std::string outTreePrefix;
  bool isData;

  double mcEvtW;

  std::vector<string> BDTNameList   ;
  std::vector<string> varSetList    ;
  std::vector<string> weightFileList;

  std::vector<string> basicBDTVarNameList;
  std::vector<string>   isrBDTVarNameList;

  std::vector<string> basicBDTVarFormulaList;
  std::vector<string>   isrBDTVarFormulaList;

  // this is a standard constructor
  ssEvtPostProc2 (std::string name="ss2lPostProc2");
  ~ssEvtPostProc2();

  int execute();

 protected:
  string m_name;//!

  SUSY::CrossSectionDB* m_XsecDB;  //!
  TFile* inF; //!

  TMVA::Reader *reader   ;    //!
  TMVA::Reader *readerISR;    //!

  vector<Float_t> basicBDTVarValues; //!
  vector<Float_t>   isrBDTVarValues; //!




  int initialize ();
  int runLoop();
  int finalize ();

 private:

  // this is needed to distribute the algorithm to the workers
  ClassDef(ssEvtPostProc2, 1);
};

#endif
