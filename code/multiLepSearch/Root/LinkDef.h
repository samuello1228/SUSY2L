#include <multiLepSearch/dibosonSelection.h>

#include <multiLepSearch/aDAODMaker.h>

#include <multiLepSearch/ssEvtSelection.h>
#include <multiLepSearch/ssEvtPostProc1.h>
#include <multiLepSearch/ssEvtPostProc2.h>

#include <multiLepSearch/ChargeFlipBkgTool.h>
#include <multiLepSearch/FakeLepBkgTool.h>

#include <multiLepSearch/obj_def.h>
#include <multiLepSearch/susyEvts.h>
#include <multiLepSearch/dilepSelection.h>
#include <multiLepSearch/dilep_objDef.h>
#include <vector>

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;

#endif

#ifdef __CINT__
#pragma link C++ class dilepSelection+;
// #pragma link C++ class anaHelper+;
#endif

#ifdef __MAKECINT__
// #pragma link C++ all struct+;
// #pragma link C++ all typedef+;
#pragma link C++ struct SIGNATURE+;
#pragma link C++ struct EVT+;
#pragma link C++ struct PAR0+;
#pragma link C++ struct PAR+;
#pragma link C++ struct R_PAR+;
#pragma link C++ struct L_PAR+;
#pragma link C++ struct J_PAR+;
#pragma link C++ struct TR_PAR+;
#pragma link C++ class vector<PAR>+;
#pragma link C++ class vector<L_PAR>+;
#pragma link C++ class vector<J_PAR>+;
#pragma link C++ class vector<TR_PAR>+;
#pragma link C++ class ntupP12+;
#pragma link C++ class ntupLep+;
#pragma link C++ class ntupJet+;
#pragma link C++ class ntupSig+;
#pragma link C++ class ntupEvtInfo+;
#pragma link C++ class vector<ntupLep>+;
#pragma link C++ class vector<ntupJet>+;
#endif

#ifdef __CINT__
#pragma link C++ class ssEvtSelection+;
#pragma link C++ class ssEvtPostProc1+;
#pragma link C++ class ssEvtPostProc2+;

#pragma link C++ class ChargeFlipBkgTool+;
#pragma link C++ class FakeLepBkgTool+;
#endif

#ifdef __CINT__
#pragma link C++ class aDAODMaker+;
#endif

#ifdef __CINT__
#pragma link C++ class dibosonSelection+;
#endif
