#include "Analysis/TruthElectron.h"

TruthElectron::TruthElectron()
: BaseElectron(){}

TruthElectron::~TruthElectron(){}

void TruthElectron::Set(){
	BaseElectron::Set();
	SetTruth(0, 0, 0);
   SetMCTruthClassifier(0,0,0,0,0);
}

void TruthElectron::Set(Double_t charge, Double_t eta, Double_t pt){
	BaseElectron::Set(charge, eta, pt);
	SetTruth(0, 0, 0);
}

void TruthElectron::SetTruth(Bool_t isE, Int_t partType, Int_t partOri){
   fIsElectron = isE;   
   fParticleType = partType;
   fParticleOrigin = partOri;

   fFromZ = 0;
   fDeltaR = 0;
}

void TruthElectron::SetMCTruthClassifier(Int_t type, Int_t origin, Int_t bkgType, Int_t bkgOrigin, Double_t bkgCharge){
   fTruthType = type;
   fTruthOrigin = origin;
   fTruthBkgType = bkgType;
   fTruthBkgOrigin = bkgOrigin;
   fTruthBkgCharge = bkgCharge;
}

void TruthElectron::PrintTruth(std::ostream& stream){
   stream << "Truth:: Type: " ;
   if (fCharge == 1) stream << "+";
   else if (fCharge == -1) stream << "-";
   stream << fTruthType << " Origin: " << fTruthOrigin;
   
   if(fTruthType == 4){
      stream << " Bkg:: Type: ";
      if (fTruthBkgCharge == 1) stream << "+";
      else if (fTruthBkgCharge == -1) stream << "-";
      stream << fTruthBkgType << " Origin: " << fTruthBkgOrigin;
   }
   stream << std::endl;
}