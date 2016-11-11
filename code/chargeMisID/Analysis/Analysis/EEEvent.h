#ifndef EEEVENT_H
#define EEEVENT_H

/* EEEvent.h
 *
 * Gabriel Gallardo
 * Version 150804.16
 * Defines Z->ee event objects 
 */

#include "TObject.h"
#include "Electron.h"
#include "TruthElectron.h"

// Z->ee event object
class EEEvent : public TObject{
public:
   // -----------------------------------------
   // Instance variables
   // -----------------------------------------
   
   // Detected electrons
   Electron fE1;
   Electron fE2; 
   
   // Event information
   Double_t fInvMass;
   Double_t fdPhi;
   Bool_t fSSevent; 
   
   // true if is MC event
   Bool_t fMCevent;
   
   // MC truth information
   TruthElectron ftE1;
   TruthElectron ftE2;
   Double_t ftInvMass;
   Double_t ftdPhi;
   Bool_t ftSSevent;
   Bool_t ftZEEevent;

   
public:
   // -----------------------------------------
   // Public methods
   // -----------------------------------------
   
   // Default constructor. Does nothing.
   EEEvent();
   
   // Saves basic information about the event
   EEEvent(Electron e1, Electron e2, Double_t mass, Double_t dPhi, Bool_t SSevent);
   virtual ~EEEvent();
   
   // Saves basic information about the event
   void Set(Electron e1, Electron e2, Double_t mass, Double_t dPhi, Bool_t SSevent);
   
   // Saves truth information about the event (use only on MC data)
   void SetTruth(TruthElectron e1, TruthElectron e2, Double_t mass, Double_t dPhi, Bool_t SSevent, Bool_t ZEEvent);
   
   // Use on non-MC events to set event as non-MC
   void SetTruth(Bool_t isMC);
   
   ClassDef(EEEvent, 1);
   
   void PrintInfo(std::ostream&);
};

#endif