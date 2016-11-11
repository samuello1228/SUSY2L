#include <multiLepSearch/obj_def.h>
#include <multiLepSearch/susyEvts.h>
#include <multiLepSearch/ljetEvts.h>
// #include <multiLepSearch/anaHelper.h>
#include <multiLepSearch/ljetSelection.h>
#include <multiLepSearch/evtSelection.h>
#include <vector>

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;

#endif

#ifdef __CINT__
#pragma link C++ class evtSelection+;
#pragma link C++ class ljetSelection+;
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
#endif
