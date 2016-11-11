/* Electron.cxx
 *
 * Gabriel Gallardo
 * Version 150804.16
 * Defines reconstructed electron objects 
 */

#include "Analysis/Electron.h"
#include <iostream>

// ELECTRON ---------------------- //

// Default constructor. Does nothing
Electron::Electron(){}
   
// Constructor sets basic information of electron, sets tag-and-probe analysis info to 0.
Electron::Electron(Double_t charge, Double_t eta, Double_t pt)
   : BaseElectron(charge, eta, pt)
{
   Set(charge, eta, pt);
   SetInfo(0, 0, 0, 0, false);    
}

// Destructor. Does nothing.
Electron::~Electron(){}
   
//Empty setter. Sets all values to 0
void Electron::Set(){
   Set(0, 0, 0);
}

// Set basic information of electron
void Electron::Set(Double_t charge, Double_t eta, Double_t pt){
   BaseElectron::Set(charge, eta, pt);
   fPassCut = 0;
   SetInfo(0, 0, 0, 0, false);
}
   
// Set additional information used for tag-and-probe analysis
void Electron::SetInfo(Double_t d0, Double_t sigma, Int_t SiHits, Int_t BHits, Bool_t TriggerMatched){
   fD0 = d0;
   fD0Error = sigma;
   fSiHits = SiHits;
   fBHits = BHits;
   fTriggerMatched = TriggerMatched;
   
   setTag();      
}

// Set whether the electron is eligible to be a tag electron
void Electron::setTag(){
   fTag = true;
   if (!fTriggerMatched) fTag = false;
   else if((fD0Error == 0) || (fD0 / fD0Error >= 1.5)) fTag = false;
   else if (abs(fEta) >= 2.0) fTag = false;
   else if (fSiHits <= 9) fTag = false;
   else if (fBHits < 1) fTag = false;

}

void Electron::PrintExtra(std::ostream& stream){
   stream <<  "d0: " << fD0 << " | d0error: " << fD0Error
            << " | Si: " << fSiHits << " | B: " << fBHits << " | Tag : ";
   if(fTag) stream << "Yes";
   else stream << "No";
   stream <<std::endl;
}

void Electron::SetPassCut(Bool_t light, Bool_t medium, Bool_t tight){
   if(light) fPassCut = fPassCut | 0x2;
   if(medium) fPassCut = fPassCut | 0x4;
   if(tight) fPassCut = fPassCut | 0x8;
}