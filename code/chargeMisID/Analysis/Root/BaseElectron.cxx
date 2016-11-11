/* BaseElectron.cxx
 * 
 * Gabriel Gallardo
 * Version 150805.16
 * Defines basic properties of electron
 */

#include "Analysis/BaseElectron.h"
#include <iostream>


BaseElectron::BaseElectron(){
   Set(0, 0, 0);
}

BaseElectron::BaseElectron(Double_t charge, Double_t eta, Double_t pt){
   Set(charge, eta, pt);
}

BaseElectron::~BaseElectron(){}

void BaseElectron::Set(){
   Set(0, 0, 0);
}

void BaseElectron::Set(Double_t charge, Double_t eta, Double_t pt){
   fCharge = charge;
   fEta = eta;
   fPt = pt;
}
   
void BaseElectron::PrintBasic(std::ostream& stream){
   if (fCharge == 1.0) stream << "Charge: +e" << "\t";
   else if (fCharge == -1.0) stream << "Charge: -e" << "\t";
   else stream << "Charge :   " << "\t";
   stream << "Pt: " << fPt << "\t" << "Eta: " << fEta << std::endl;
}