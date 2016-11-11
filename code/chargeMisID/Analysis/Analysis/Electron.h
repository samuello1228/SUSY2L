#ifndef ELECTRON_H
#define ELECTRON_H

/* Electron.h
 *
 * Gabriel Gallardo
 * Version 150804.16
 * Defines reconstructed electron objects 
 */
#include "BaseElectron.h"

class Electron : public BaseElectron{
public:
   // -----------------------------------------
   // Instance variables
   // -----------------------------------------
   
   // Detector information
   Double_t fD0;
   Double_t fD0Error;
   Int_t fSiHits;
   Int_t fBHits;
   Bool_t fTriggerMatched;
   
   
   // Tag-and-probe information
   Bool_t fTag;
   
   // Pass cut indicator (binary)
   // (fPassCut & 0x2) => passed light cut;
   // (fPassCut & 0x4) => passed medium cut;
   // (fPassCut & 0x8) => passed tight cut;
   Int_t fPassCut;
 

  
   // -----------------------------------------
   // Public methods
   // -----------------------------------------
   
   // Default constructor. Does nothing.
   Electron();
   
   // Constructor sets basic information of electron, sets tag-and-probe analysis info to 0.
   Electron(Double_t charge, Double_t eta, Double_t pt);
   virtual ~Electron();
   
   // Empty setter, sets all values to 0
   void Set();
   
   // Set basic information of electron
   void Set(Double_t charge, Double_t eta, Double_t pt);
   
   // Set additional information used for tag-and-probe analysis, calls setTag()
   void SetInfo(Double_t d0, Double_t sigma, Int_t SiHits, Int_t BHits, Bool_t TriggerMatched);
   
   // Set whether the electron is eligible to be a tag electron
   void setTag();
   
   void PrintExtra(std::ostream&);

   void SetPassCut(Bool_t light, Bool_t medium, Bool_t tight);

   ClassDef(Electron, 1);

};
#endif