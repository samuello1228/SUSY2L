#ifndef BASE_ELECTRON_H
#define BASE_ELECTRON_H

#include "TObject.h"
#include <iostream>

/* BaseElectron.h
 * 
 * Gabriel Gallardo
 * Version 150804.16
 * Defines basic properties of electron
 * For inheriting
 */

class BaseElectron : public TObject{
public:
   // Basic properties of electron
   
   Double_t fCharge;
   Double_t fEta;
   Double_t fPt;
   
   // Methods
   BaseElectron();
   
   BaseElectron(Double_t charge, Double_t eta, Double_t pt);
   
   ~BaseElectron();
   
   virtual void Set();
   
   void Set(Double_t charge, Double_t eta, Double_t pt);
   
   void PrintBasic(std::ostream& stream);
   
   ClassDef(BaseElectron, 1);
};

#endif