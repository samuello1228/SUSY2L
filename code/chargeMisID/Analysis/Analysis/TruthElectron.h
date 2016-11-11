#ifndef TRUTH_ELECTRON_H
#define TRUTH_ELECTRON_H

/* TruthElectron.h
 *
 * Gabriel Gallardo
 * Version 150804.16
 * Defines truth electron objects 
 */

#include "BaseElectron.h"

class TruthElectron : public BaseElectron{
public:
	// For MC
   Int_t fParticleType; // pdgId of TruthParticle
   Int_t fParticleOrigin; // pdgId of 0-th parent of TruthParticle

   Int_t fTruthType;
   Int_t fTruthOrigin;
   Int_t fTruthBkgType;
   Int_t fTruthBkgOrigin;
   Double_t fTruthBkgCharge;

	Bool_t fIsElectron;
   Bool_t fFromZ;
   
   Double_t fDeltaR;

	TruthElectron();
	~TruthElectron();
	void Set(); // Make virtual?
	void Set(Double_t charge, Double_t eta, Double_t pt); // Make virtual?
   void SetTruth(Bool_t isE, Int_t partType, Int_t partOrig);
   void SetMCTruthClassifier(Int_t type, Int_t origin, Int_t bkgType, Int_t bkgOrigin, Double_t bkgCharge);
   void PrintTruth(std::ostream&);

	ClassDef(TruthElectron, 1);
};
#endif