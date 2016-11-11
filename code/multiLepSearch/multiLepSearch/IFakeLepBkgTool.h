#ifndef IFakeLepBkgTool_H
#define IFakeLepBkgTool_H

#include "AsgTools/IAsgTool.h"

#include "xAODEgamma/Electron.h"
#include "xAODMuon/Muon.h"

#include <vector>

class IFakeLepBkgTool : virtual public asg::IAsgTool {
  ASG_TOOL_INTERFACE(IFakeLepBkgTool)

public:

  virtual double GetWeight(std::vector<xAOD::IParticle*> pList, int sigma=0, int nthSys=0) const = 0;

  virtual double GetWeight(std::vector<double> &eta, std::vector<double> &pt, std::vector<int> &pdgID, std::vector<bool> &isSig,
		           int sigma=0, int nthSys=0) const =0;

};

#endif
