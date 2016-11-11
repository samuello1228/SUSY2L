#ifndef FakeLepBkgTool_H
#define FakeLepBkgTool_H

// David Adams.
// January 2014
//
// This is a simple ASG dual-use tool intended as an
// example for tool developers.

#include "multiLepSearch/IFakeLepBkgTool.h"
#include "AsgTools/AsgTool.h"

#include "TFile.h"
#include "TH2D.h"

class FakeLepBkgTool
: public asg::AsgTool,
  virtual public IFakeLepBkgTool {
ASG_TOOL_CLASS(FakeLepBkgTool, IFakeLepBkgTool)

public:

  FakeLepBkgTool( const std::string& myname );

  virtual StatusCode initialize();

  virtual StatusCode finalize();

  virtual double GetWeight(std::vector<xAOD::IParticle*> pList, int sigma=0, int nthSys=0) const;

  virtual double GetWeight(std::vector<double> &eta, std::vector<double> &pt, std::vector<int> &pdgID, std::vector<bool> &isSig,
		           int sigma=0, int nthSys=0) const;


private:
  enum method_t {kMATRIX, kFAKE_FACTOR};
  
  std::string m_inputFileName;
  std::string m_method;
  method_t m_methodEnum;
  

  std::string m_RealeEffHistoName;
  std::string m_RealuEffHistoName;
  std::string m_FakeeEffHistoName;
  std::string m_FakeuEffHistoName;

  std::string m_eFakeFactorHistoName;
  std::string m_uFakeFactorHistoName;

  TFile* inputFile;
  TH2D* hRealeEff;
  TH2D* hRealuEff;
  TH2D* hFakeeEff;
  TH2D* hFakeuEff;
  TH2D* heFakeFactor;
  TH2D* huFakeFactor;
};

#endif
