#ifndef multiLepSearch_aDAODMaker_H
#define multiLepSearch_aDAODMaker_H

#include <EventLoop/Algorithm.h>
#include "SUSYTools/SUSYObjDef_xAOD.h"

//trigger
// #include <TrigConfxAOD/xAODConfigTool.h>
// #include <TrigDecisionTool/TrigDecisionTool.h>


class aDAODMaker : public EL::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  // float cutValue;
  int CF_isMC;
  std::vector< std::string > CF_grlFiles;
  std::string CF_ConfigFile;
  std::vector< std::string > CF_trigNames;




  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:
  // Tree *myTree; //!
  // TH1 *myHist; //!



  // this is a standard constructor
  aDAODMaker ();

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

 private:
  std::string m_name; //!

  ST::SUSYObjDef_xAOD* m_objTool; //!

  // this is needed to distribute the algorithm to the workers
  ClassDef(aDAODMaker, 1);
};

#endif