#include "multiLepSearch/ChargeFlipBkgTool.h"

#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODEgamma/Electron.h"

using namespace std;

static SG::AuxElement::Decorator<char> dec_signal("signal");

//**********************************************************************

ChargeFlipBkgTool::
ChargeFlipBkgTool( const string& myname )
: AsgTool(myname) {
  declareProperty("InputRatesFileName" , m_inputRatesFileName = "$ROOTCOREBIN/data/multiLepSearch/root_files/chargeMisID_Zee_MC_looseBaseline.root");
  declareProperty("InputRatesHistoName", m_inputRatesHistoName= "hFlipProb_stat"      );
  declareProperty("InputRatesSysHistoName", m_inputRatesSysHistoName = "hFlipProb_sys");
  declareProperty("InputDPtFileName", m_inputDPtFileName      = "$ROOTCOREBIN/data/multiLepSearch/root_files/dPT_loose.root");
  declareProperty("InputDPtOkHistoName", m_inputDPtOkHistoName= "hDPTok_pxy" );
  declareProperty("InputDPtFlippedHistoName", m_inputDPtFlippedHistoName = "hDPTflipped_pxy");

  ratesFile   = NULL;
  hFlipProb   = NULL;
  dPtFile     = NULL;
  hdPTok      = NULL;
  hdPTflipped = NULL;
}

//**********************************************************************

StatusCode ChargeFlipBkgTool::initialize() {

  ATH_MSG_VERBOSE( "Initialising tool " << name() );

  ratesFile = new TFile(m_inputRatesFileName.c_str(), "READ");
  if (!ratesFile) return  StatusCode::FAILURE;

  hFlipProb = (TH2D*) ratesFile->FindObjectAny(m_inputRatesHistoName.c_str());
  if (!hFlipProb) return  StatusCode::FAILURE;

  hFlipProb_sys = (TH2D*) ratesFile->FindObjectAny(m_inputRatesSysHistoName.c_str());
  if(!hFlipProb_sys) return StatusCode::FAILURE;

  dPtFile = new TFile(m_inputDPtFileName.c_str(), "READ");
  if (!dPtFile) return StatusCode::FAILURE;

  hdPTok = (TH2D*) dPtFile->FindObjectAny(m_inputDPtOkHistoName.c_str());
  if (!hdPTok) return StatusCode::FAILURE;

  hdPTflipped = (TH2D*) dPtFile->FindObjectAny(m_inputDPtFlippedHistoName.c_str());
  if (!hdPTflipped) return StatusCode::FAILURE;

  return StatusCode::SUCCESS;
}


StatusCode ChargeFlipBkgTool::finalize() {
  if (ratesFile) ratesFile->Close();
  if (dPtFile) dPtFile->Close();
  return StatusCode::SUCCESS;
}

//**********************************************************************
//
double ChargeFlipBkgTool::GetWeight(vector<xAOD::IParticle*> &pList, int sigma, int nthSys) const{
  ATH_MSG_DEBUG("GetWeight -- with IParticle");

  if (!dec_signal(*pList[0]) || !dec_signal(*pList[1])) return 0.0; //must be 2 signal leps


  vector<double> eta ; eta.reserve(2);
  vector<double> pt  ; pt.reserve(2);

  xAOD::Electron* ele;
  for (auto aParticle : pList){
    ele = dynamic_cast<xAOD::Electron*>(aParticle);
    if (ele && dec_signal(*aParticle) ){ //require signal electron
      eta.push_back( ele->eta() );
      pt.push_back( ele->pt()*0.001 ); //*0.001 to GeV
    }
  }

  return GetWeight(eta, pt, sigma, nthSys);
}
 
//inputs are the eta and pt(GeV) of signal electrons
double ChargeFlipBkgTool::GetWeight(vector<double> &eta, vector<double> &pt, int sigma, int nthSys) const{
  ATH_MSG_DEBUG("GetWeight");

  int bin;
  vector<double > prob; prob.reserve(2);

  for (unsigned int i=0; i<eta.size(); i++){
      bin = hFlipProb->FindBin( fabs(eta[i]), pt[i]);
      prob.push_back( bin==-1 ? 0.0 : hFlipProb->GetBinContent(bin) );
      if (bin==-1) ATH_MSG_ERROR("No flipProb defined for this e: " << fabs(eta[i]) << " " << pt[i] );
  }
  
  double weight = 0.0, ptot = 0.0;
  if (prob.size() >2){
    ATH_MSG_ERROR("More then 2 signal e in input, use only first 2");
    ptot = prob[0]+prob[1]-2*prob[0]*prob[1];
  }else if (prob.size()==2){
    ATH_MSG_DEBUG("2P: " << prob[0] << " " << prob[1]);
    ptot = prob[0]+prob[1]-2*prob[0]*prob[1];
  }else if (prob.size()==1){
    ATH_MSG_DEBUG("1P: " << prob[0] );
    ptot = prob[0];
  }
  weight = ptot/(1.-ptot);

  if (sigma==0) return weight;

  //dummy systematics for now
  double uncertainty = 0;
  vector<double > unc; unc.reserve(2);
  if (nthSys==0){
    for (unsigned int i=0; i<eta.size(); i++){
      bin = hFlipProb->FindBin( fabs(eta[i]), pt[i]);
      unc.push_back( bin==-1 ? 0.0 : hFlipProb->GetBinError(bin) );
      if (bin==-1) ATH_MSG_ERROR("No flipProb defined for this e: " << fabs(eta[i]) << " " << pt[i] );
    }
  } else if (nthSys==1) {
    for (unsigned int i=0; i<eta.size(); i++){
      bin = hFlipProb_sys->FindBin( fabs(eta[i]), pt[i]);
      unc.push_back( bin==-1 ? 0.0 : hFlipProb_sys->GetBinError(bin) );
      if (bin==-1) ATH_MSG_ERROR("No flipProb defined for this e: " << fabs(eta[i]) << " " << pt[i] );
    }
  }

  if(prob.size()>=2){
    // Uncertainty calc based on w = p0 + p1 - 2*p0*p1
    uncertainty = 1./pow(1.-ptot,2) * sqrt( pow((1.-2*prob[1])*unc[0],2) + pow((1.-2*prob[0])*unc[1],2) );
    //ATH_MSG_ERROR("prob: " << prob[0] << " " << prob[1]);
    //ATH_MSG_ERROR("unc: " << unc[0] << " " << unc[1] << " " << uncertainty);
  } else if(prob.size()==1) {
    uncertainty = 1./pow(1.-ptot,2) * unc[0];
  }

  return weight + sigma*uncertainty;
}
//**********************************************************************

//this do not modify the input pList
std::vector<double> ChargeFlipBkgTool::GetCorrectedPt(std::vector<xAOD::IParticle*> &pList, int sigma, int nthSys) const
{
  ATH_MSG_DEBUG("GetCorrectedPt -- with IParticle");

  std::vector<double> correctedPt;
  for(auto aParticle : pList){
    correctedPt.push_back(aParticle->pt() *0.001); // *0.001 to GeV
  }

  std::vector<int> lepCharge(2,0);
  std::vector<int> lepPdgID(2,0);
  xAOD::Electron* tmpE;  xAOD::Muon* tmpMu;

  tmpE = NULL; tmpMu = NULL;
  tmpE = dynamic_cast<xAOD::Electron*>(pList[0]);
  if (tmpE){ lepCharge[0] = tmpE->charge(); lepPdgID[0] = 11 * tmpE->charge();}
  else{ 
    tmpMu = dynamic_cast<xAOD::Muon*>(pList[0]);
    if (tmpMu){ lepCharge[0] = tmpMu->charge(); lepPdgID[0] = 13 * tmpMu->charge();}
  }
  
  tmpE = NULL; tmpMu = NULL;
  tmpE = dynamic_cast<xAOD::Electron*>(pList[1]);
  if (tmpE){ lepCharge[1] = tmpE->charge(); lepPdgID[1] = 11 * tmpE->charge();}
  else{ 
    tmpMu = dynamic_cast<xAOD::Muon*>(pList[1]);
    if (tmpMu){ lepCharge[1] = tmpMu->charge(); lepPdgID[1] = 13 * tmpMu->charge();}
  }

  if (!dec_signal(*pList[0]) || !dec_signal(*pList[1])  //must be 2 signal leps
      || lepCharge[0]*lepCharge[1]==-1 )     //must be OS events
  {
    ATH_MSG_DEBUG("Not an OS 2 signal lepton event. No pt correction applied.");
    return correctedPt;
  }

  if (abs(lepPdgID[0])!=11 && abs(lepPdgID[1])!=11){
    ATH_MSG_DEBUG("Neither signal lepton in OS pair is an electron. No pt correction applied");
    return correctedPt;
  }

  std::vector<double> eta; eta.push_back(pList[0]->eta()); eta.push_back(pList[1]->eta());

  // If both are electrons, apply correction with electron most likely to flip
  if(abs(lepPdgID[0])==11 && abs(lepPdgID[1]==11))
    return GetCorrectedPt(eta, correctedPt, sigma, nthSys);

  // Find the signal lepton which is an electron and apply pt correction
  for(int i=0; i<2; i++){
    // continue if not an electron
    if(!dynamic_cast<xAOD::Electron*>(pList[i])) continue;
    correctedPt[i] = GetCorrectedPt(pList[i]->eta(), pList[i]->pt(), sigma, nthSys);
  }

  return correctedPt;
}

//this modify the input pt
std::vector<double> ChargeFlipBkgTool::GetCorrectedPt(std::vector<double> &eta, std::vector<double> &pt, int sigma, int nthSys) const
{
  ATH_MSG_DEBUG("GetCorrectedPt -- with vector eta & pt");

  double flipProb0 = hFlipProb->GetBinContent(hFlipProb->FindBin(fabs(eta[0]), pt[0]));
  double flipProb1 = hFlipProb->GetBinContent(hFlipProb->FindBin(fabs(eta[1]), pt[1]));

  int i = ((flipProb0>flipProb1) ? 0 : 1);
  pt[i] = GetCorrectedPt(eta[i], pt[i], sigma, nthSys);

  return pt;
}

//this do not modify the input pt
double ChargeFlipBkgTool::GetCorrectedPt(const double &eta, const double &pt, const int sigma, const int nthSys) const
{
  ATH_MSG_DEBUG("GetCorrectedPt");
  int bin = hdPTflipped->FindBin(pt, fabs(eta));
  double correctedPt = pt + hdPTflipped->GetBinContent(bin) - hdPTok->GetBinContent(bin);

  if(nthSys==0 && sigma==1) correctedPt += sqrt(pow(hdPTflipped->GetBinErrorUp(bin),2) + pow(hdPTok->GetBinErrorUp(bin),2));
  if(nthSys==0 && sigma==-1) correctedPt -= sqrt(pow(hdPTflipped->GetBinErrorLow(bin),2) + pow(hdPTok->GetBinErrorLow(bin),2));

  return correctedPt;
}
