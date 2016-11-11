#include "multiLepSearch/FakeLepBkgTool.h"

#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODEgamma/Electron.h"
#include "xAODMuon/Muon.h"

#include "Math/SMatrix.h"

using namespace std;

static SG::AuxElement::Decorator<char> dec_signal("signal");

//**********************************************************************

FakeLepBkgTool::
FakeLepBkgTool( const string& myname )
: AsgTool(myname) {
  declareProperty("InputFileName" , m_inputFileName = "$ROOTCOREBIN/data/multiLepSearch/root_files/RealFakeLepEff_dummy.root");
  declareProperty("Method", m_method= "FakeFactor");  //choice: FakeFactor / Matrix

  declareProperty("RealeEffHistoName", m_RealeEffHistoName= "RealeEff");
  declareProperty("RealuEffHistoName", m_RealuEffHistoName= "RealuEff");
  declareProperty("FakeeEffHistoName", m_FakeeEffHistoName= "FakeeEff");
  declareProperty("FakeuEffHistoName", m_FakeuEffHistoName= "FakeuEff");

  declareProperty("eFakeFactorHistoName", m_eFakeFactorHistoName= "h_ff_ele_v2");
  declareProperty("uFakeFactorHistoName", m_uFakeFactorHistoName= "h_ff_mu");


  inputFile = NULL;
  hRealeEff = NULL; hFakeeEff = NULL;
  hRealuEff = NULL; hFakeuEff = NULL;
}

//**********************************************************************

StatusCode FakeLepBkgTool::initialize() {

   ATH_MSG_VERBOSE( "Initialising tool " << name() );

   inputFile = new TFile(m_inputFileName.c_str(), "READ");
   if (!inputFile) return  StatusCode::FAILURE;

   if      (m_method=="Matrix"    ){ m_methodEnum = kMATRIX; }
   else if (m_method=="FakeFactor"){ m_methodEnum = kFAKE_FACTOR; }
   else { 
     ATH_MSG_ERROR("Unrecognized method choice for FakeLep background estimation: " << m_method);
     return StatusCode::FAILURE;
   }

   if ( m_methodEnum == kMATRIX){
     hRealeEff = (TH2D*) inputFile->FindObjectAny(m_RealeEffHistoName.c_str());
     if (!hRealeEff) return  StatusCode::FAILURE;

     hRealuEff = (TH2D*) inputFile->FindObjectAny(m_RealuEffHistoName.c_str());
     if (!hRealuEff) return  StatusCode::FAILURE;

     hFakeeEff = (TH2D*) inputFile->FindObjectAny(m_FakeeEffHistoName.c_str());
     if (!hFakeeEff) return  StatusCode::FAILURE;

     hFakeuEff = (TH2D*) inputFile->FindObjectAny(m_FakeuEffHistoName.c_str());
     if (!hFakeuEff) return  StatusCode::FAILURE;
   }
   else if ( m_methodEnum == kFAKE_FACTOR){
     heFakeFactor = (TH2D*) inputFile->FindObjectAny(m_eFakeFactorHistoName.c_str());
     if (!heFakeFactor) return  StatusCode::FAILURE;

     huFakeFactor = (TH2D*) inputFile->FindObjectAny(m_uFakeFactorHistoName.c_str());
     if (!huFakeFactor) return  StatusCode::FAILURE;
   }

   return StatusCode::SUCCESS;
}

StatusCode FakeLepBkgTool::finalize() {
   if (inputFile) inputFile->Close();
   return StatusCode::SUCCESS;
}

//**********************************************************************
 
double FakeLepBkgTool::GetWeight(vector<xAOD::IParticle*> pList, int sigma, int nthSys) const{
  ATH_MSG_DEBUG("GetWeight -- with IParticle");

  if (pList.size()<2) return 0.0; //safty check

  //stupid slow way to check particle type
  xAOD::Electron* e1 = dynamic_cast<xAOD::Electron*>(pList[0]);
  xAOD::Electron* e2 = dynamic_cast<xAOD::Electron*>(pList[1]);
  xAOD::Muon*     u1 = dynamic_cast<xAOD::Muon*    >(pList[0]);
  xAOD::Muon*     u2 = dynamic_cast<xAOD::Muon*    >(pList[1]);

  vector<double> pt, eta;
  vector<int> pdgID;
  vector<bool> isSig;

  //particle 1
  if (e1){
    pt.push_back( e1->pt()*0.001 ); eta.push_back(e1->eta()); pdgID.push_back(11);
  }else if (u1){
    pt.push_back( u1->pt()*0.001 ); eta.push_back(u1->eta()); pdgID.push_back(13);
  }else{
    return 0.0;
  }
  isSig.push_back( dec_signal(*pList[0]) );

  //particle 2
  if (e2){
    pt.push_back( e2->pt()*0.001 ); eta.push_back(e2->eta()); pdgID.push_back(11);
  }else if (u2){
    pt.push_back( u2->pt()*0.001 ); eta.push_back(u2->eta()); pdgID.push_back(13);
  }else{
    return 0.0;
  }
  isSig.push_back( dec_signal(*pList[1]) );

  return GetWeight(eta, pt, pdgID, isSig, sigma, nthSys);

}

double FakeLepBkgTool::GetWeight(vector<double> &eta, vector<double> &pt, vector<int> &pdgID, vector<bool> &isSig, int sigma, int nthSys) const{
  ATH_MSG_DEBUG("GetWeight");

  if (pt.size()<2) return 0.0; //safty check

  double weight = 0.0;

  if ( m_methodEnum == kMATRIX){

    int binR1, binR2;
    int binF1, binF2;
    double realEff1, realEff2, fakeEff1, fakeEff2;
    TH2D *hRealEff, *hFakeEff;

    //particle 1
    if      (abs(pdgID[0])==11){ hRealEff = hRealeEff; hFakeEff = hFakeeEff;}
    else if (abs(pdgID[0])==13){ hRealEff = hRealuEff; hFakeEff = hFakeuEff;}
    else { return 0;}

    binR1    = hRealEff->Fill( pt[0], fabs(eta[0]), 0. ); if (binR1==-1) return 0.0;
    realEff1 = hRealEff->GetBinContent(binR1);
    binF1    = hFakeEff->Fill( pt[0], fabs(eta[0]), 0. ); if (binF1==-1) return 0.0;
    fakeEff1 = hFakeEff->GetBinContent(binF1);

    //particle 2
    if      (abs(pdgID[1])==11){ hRealEff = hRealeEff; hFakeEff = hFakeeEff;}
    else if (abs(pdgID[1])==13){ hRealEff = hRealuEff; hFakeEff = hFakeuEff;}
    else { return 0;}

    binR2    = hRealEff->Fill( pt[1], fabs(eta[1]), 0. ); if (binR2==-1) return 0.0;
    realEff2 = hRealEff->GetBinContent(binR2);
    binF2    = hFakeEff->Fill( pt[1], fabs(eta[1]), 0. ); if (binF2==-1) return 0.0;
    fakeEff2 = hFakeEff->GetBinContent(binF2);

    //--------------------------------------------------------
    //Use matrix method, via inverting the matrix numerically
    //--------------------------------------------------------
    //double elements[16] = 
    //  { 
    //    (   realEff1)*(   realEff2), (   realEff1)*(   fakeEff2), (   fakeEff1)*(   realEff2), (   fakeEff1)*(   fakeEff2) ,
    //    (   realEff1)*(1.-realEff2), (   realEff1)*(1.-fakeEff2), (   fakeEff1)*(1.-realEff2), (   fakeEff1)*(1.-fakeEff2) ,
    //    (1.-realEff1)*(   realEff2), (1.-realEff1)*(   fakeEff2), (1.-fakeEff1)*(   realEff2), (1.-fakeEff1)*(   fakeEff2) ,
    //    (1.-realEff1)*(1.-realEff2), (1.-realEff1)*(1.-fakeEff2), (1.-fakeEff1)*(1.-realEff2), (1.-fakeEff1)*(1.-fakeEff2)
    //  };

    //ROOT::Math::SMatrix<double,4> mat(elements,16);
    //if (!mat.Invert()){
    //  ATH_MSG_ERROR("eff " << realEff1 << " " << realEff2 << " " << fakeEff1 << " " << fakeEff2);
    //  ATH_MSG_ERROR("bin " << binR1 << " " << binR2 << " " << binF1 << " " << binF2);
    //  return 0.0;
    //}
    //
    //int tarCol = ((!dec_signal(*pList[0]))<<1) + ((!dec_signal(*pList[1]))<<0);
    //double weight = mat(1, tarCol)* realEff1 * fakeEff2 + 
    //                mat(2, tarCol)* fakeEff1 * realEff2 + 
    //                mat(3, tarCol)* fakeEff1 * fakeEff2 ;
    
    //--------------------------------------------------------
    //Use matrix method, via inverting the matrix analyticaly
    //For the formula see:
    //http://live.sympy.org/?evaluate=e1%2Ce2%2Cf1%2Cf2%20%3D%20symbols%28%22e1%20e2%20f1%20f2%22%29%0A%23--%0Ae1%0A%23--%0Am%20%3D%20Matrix%28[[e1%2Ce2]%2C[f1%2Cf2]]%29%0A%23--%0Am%0A%23--%0Am.inv%28Abs%28%29%0A%23--%0Am%0A%23--%0Am.inv%0A%23--%0Am.inv%28%29%0A%23--%0Am.inv%28%29*m%0A%23--%0Asimplify%28m.inv%28%29*m%29%0A%23--%0Asimplify%28m.inv%28%29%29%0A%23--%0Am%20%3D%20Matrix%28[[e1*e2%2Ce1*f2%2Cf1*e2%2Cf1*f2]%2C[e1*%281-e2%29%2Ce1*%281-f2%29%2Cf1*%281-e2%29%2Cf1*%281-f2%29]%2C[%281-e1%29*e2%2C%281-e1%29*f2%2C%281-f1%29*e2%2C%281-f1%29*f2]%2C[%281-e1%29*%281-e2%29%2C%281-e1%29*%281-f2%29%2C%281-f1%29*%281-e2%29%2C%281-f1%29*%281-f2%29]]%29%0A%23--%0Am%0A%23--%0Am2%20%3D%20m.inv%28%29%0A%23--%0Asimplify%28m2%29%0A%23--%0Am2s%20%3D%20simplify%28m2%29%0A%23--%0A
    //--------------------------------------------------------
    double d = 1./(realEff1*realEff2 + fakeEff1*fakeEff2 - realEff1*fakeEff2 - fakeEff1*realEff2);
    double Nrf, Nfr, Nff;

    if       ( isSig[0] &&  isSig[1]){
      //Ntt
      Nrf = d * ( -realEff2*fakeEff1 + realEff2 + fakeEff1 -1. );
      Nfr = d * ( -realEff1*fakeEff2 + realEff1 + fakeEff2 -1. );
      Nff = d * (  realEff1*realEff2 - realEff1 - realEff2 +1. );
    }else if ( isSig[0] && !isSig[1]){
      //Ntl
      Nrf = d * ( -realEff2*(fakeEff1-1.) );
      Nfr = d * ( -fakeEff2*(realEff1-1.) );
      Nff = d * (  realEff2*(realEff1-1.) );
    }else if (!isSig[0] &&  isSig[1]){
      //Nlt
      Nrf = d * ( -fakeEff1*(realEff2-1.) );
      Nfr = d * ( -realEff1*(fakeEff2-1.) );
      Nff = d * (  realEff1*(realEff2-1.) );
    }else if (!isSig[0] && !isSig[1]){
      //Nll
      Nrf = d * ( -realEff2*fakeEff1 );
      Nfr = d * ( -realEff1*fakeEff2 );
      Nff = d * (  realEff1*realEff2 );
    }else{
      ATH_MSG_ERROR("Impossible case");
      Nrf = 0.;
      Nfr = 0.;
      Nff = 0.;
    }

    weight = Nrf * realEff1 * fakeEff2 + 
             Nfr * fakeEff1 * realEff2 + 
             Nff * fakeEff1 * fakeEff2 ;
  }

  else if ( m_methodEnum == kFAKE_FACTOR){

    int binF = -1;
    TH2D* hFakeFactor = NULL;

    weight = 1.0;

    //particle 1
    if (!isSig[0]){
      if      (abs(pdgID[0])==11){ hFakeFactor = heFakeFactor;}
      else if (abs(pdgID[0])==13){ hFakeFactor = huFakeFactor;}
      else                       { hFakeFactor = NULL;}

      if (hFakeFactor){
        binF    = hFakeFactor->Fill( pt[0], fabs(eta[0]), 0. );
        weight *= (binF==-1? 0.0 : hFakeFactor->GetBinContent(binF));
      }
    }

    //particle 2
    if (!isSig[1]){
      if      (abs(pdgID[1])==11){ hFakeFactor = heFakeFactor;}
      else if (abs(pdgID[1])==13){ hFakeFactor = huFakeFactor;}
      else { hFakeFactor = NULL;}

      if (hFakeFactor){
        binF    = hFakeFactor->Fill( pt[1], fabs(eta[1]), 0. );
        weight *= (binF==-1? 0.0 : hFakeFactor->GetBinContent(binF));
      }
    }

    if (isSig[0] && isSig[1]){weight=0.0;}

  }

  if (std::isnan(weight)){
    ATH_MSG_WARNING("Get NaN weight, now set weight to 0 instead");
    ATH_MSG_WARNING("eta   " <<   eta[0] << " " <<   eta[1]);
    ATH_MSG_WARNING("pt    " <<    pt[0] << " " <<    pt[1]);
    ATH_MSG_WARNING("pdgID " << pdgID[0] << " " << pdgID[1]);
    ATH_MSG_WARNING("isSig " << isSig[0] << " " << isSig[1]);
    weight = 0.0;
  }
  
  //dummy systematics for now
  if (nthSys==0 && sigma== 1){ weight*=1.3; }
  if (nthSys==0 && sigma==-1){ weight*=0.7; }

  return weight;
}
//**********************************************************************
