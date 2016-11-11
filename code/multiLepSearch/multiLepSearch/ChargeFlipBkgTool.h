#ifndef ChargeFlipBkgTool_H
#define ChargeFlipBkgTool_H

// David Adams.
// January 2014
//
// This is a simple ASG dual-use tool intended as an
// example for tool developers.

#include "multiLepSearch/IChargeFlipBkgTool.h"
#include "AsgTools/AsgTool.h"

#include "TFile.h"
#include "TH2D.h"

class ChargeFlipBkgTool
: public asg::AsgTool,
  virtual public IChargeFlipBkgTool {
ASG_TOOL_CLASS(ChargeFlipBkgTool, IChargeFlipBkgTool)

public:

  ChargeFlipBkgTool( const std::string& myname );

  virtual StatusCode initialize();

  virtual StatusCode finalize();

  virtual double GetWeight(std::vector<xAOD::IParticle*> &pList, int sigma=0, int nthSys=0) const;

  virtual double GetWeight(std::vector<double> &eta, std::vector<double> &pt, int sigma=0, int nthSys=0) const;

  virtual std::vector<double> GetCorrectedPt(std::vector<xAOD::IParticle*> &pList, int sigma=0, int nthSys=0) const;

  virtual std::vector<double> GetCorrectedPt(std::vector<double> &eta, std::vector<double> &pt, int sigma=0, int nthSys=0) const;

  virtual double GetCorrectedPt(const double &eta, const double &pt, const int sigma, const int nthSys) const;
  
private:

  std::string m_inputRatesFileName;
  std::string m_inputRatesHistoName;
  std::string m_inputRatesSysHistoName;

  std::string m_inputDPtFileName;
  std::string m_inputDPtOkHistoName;
  std::string m_inputDPtFlippedHistoName;


  TFile* ratesFile;
  TH2D* hFlipProb;
  TH2D* hFlipProb_sys;

  TFile* dPtFile;
  TH2D* hdPTok;
  TH2D* hdPTflipped;

};

#endif
