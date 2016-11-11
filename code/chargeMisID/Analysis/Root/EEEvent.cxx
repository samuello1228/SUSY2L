/* EEEvent.cxx
 *
 * Gabriel Gallardo
 * Version 150622.12
 * Defines Z->ee event objects 
 */

#include "Analysis/EEEvent.h"
#include <iostream>
#include "TMath.h"
#include "xAODEgamma/Egamma.h"

// Default constructor. Does nothing.
EEEvent::EEEvent(){}

// Saves basic information about the event
EEEvent::EEEvent(Electron e1, Electron e2, Double_t mass, Double_t dPhi, Bool_t SSevent){
   Set(e1, e2, mass, dPhi, SSevent);
}

// Destructor. Does nothing.
EEEvent::~EEEvent(){}

// Saves basic information about the event
void EEEvent::Set(Electron e1, Electron e2, Double_t mass, Double_t dPhi, Bool_t SSevent){
   fE1 = e1; 
   fE2 = e2; 
   fInvMass = mass;
   fdPhi = dPhi;
   fSSevent = SSevent;
}


void EEEvent::PrintInfo(std::ostream& stream){
      
      stream << "> ------------ ELECTRON 1 ------------" << std::endl;   
      stream << "> "; fE1.PrintBasic(stream);
      stream << "> "; fE1.PrintExtra(stream);
      stream << "> "; ftE1.PrintTruth(stream);
      stream << "> ------------ ELECTRON 2 ------------" << std::endl;
      stream << "> "; fE2.PrintBasic(stream);
      stream << "> "; fE2.PrintExtra(stream);
      stream << "> "; ftE2.PrintTruth(stream);
      stream << "> ------------ EVENT INFO ------------" << std::endl;
   
      stream << "> Invariant mass of electron pair = " << fInvMass << " GeV" << std::endl;
      stream << "> dPhi of electron pair = " << fdPhi / TMath::Pi() *180 << " degrees" << std::endl; 
      stream << "========================================" << std::endl;
      stream << std::endl;

}

// Saves truth information about the event
void EEEvent::SetTruth(TruthElectron e1, TruthElectron e2, Double_t mass, Double_t dPhi, Bool_t SSevent, Bool_t isZevent){
   fMCevent = true;
   ftE1 = e1; 
   ftE2 = e2; 
   ftInvMass = mass;
   ftdPhi = dPhi;
   ftSSevent = SSevent;
   ftZEEevent = isZevent;

}

// Use this for non MC events
void EEEvent::SetTruth(Bool_t isMC){
   fMCevent = isMC;
   ftE1.Set(); 
   ftE2.Set();
   ftInvMass = ftdPhi = 0;
   ftSSevent = 0;
}
