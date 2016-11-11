#ifndef IChargeFlipBkgTool_H
#define IChargeFlipBkgTool_H

#include "AsgTools/IAsgTool.h"

#include "xAODEgamma/Electron.h"
#include "xAODMuon/Muon.h"

#include <vector>

class IChargeFlipBkgTool : virtual public asg::IAsgTool {
  ASG_TOOL_INTERFACE(IChargeFlipBkgTool)

public:

  virtual double GetWeight(std::vector<xAOD::IParticle*> &pList, int sigma=0, int nthSys=0) const = 0;

  virtual double GetWeight(std::vector<double> &eta, std::vector<double> &pt, int sigma=0, int nthSys=0) const = 0;

  virtual std::vector<double> GetCorrectedPt(std::vector<xAOD::IParticle*> &pList, int sigma=0, int nthSys=0) const = 0;

  virtual std::vector<double> GetCorrectedPt(std::vector<double> &eta, std::vector<double> &pt, int sigma=0, int nthSys=0) const = 0;

  virtual double GetCorrectedPt(const double &eta, const double &pt, const int sigma, const int nthSys) const = 0;

};

#endif
