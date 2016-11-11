#ifndef Analysis_TwoElectrons_H
#define Analysis_TwoElectrons_H

//#define GRL // Uncomment to include GRL code
#define DECAY 

#include <EventLoop/Algorithm.h>
#include "TH1.h"
#include "TTree.h"
#include "xAODEgamma/ElectronContainer.h"
#include "EEEvent.h"
#include "ElectronPhotonSelectorTools/AsgElectronLikelihoodTool.h"
#include "MCTruthClassifier/MCTruthClassifier.h"
#include "xAODTruth/TruthParticleContainer.h"
#include <fstream>


#ifdef GRL
#include "GoodRunsLists/GoodRunsListSelectionTool.h"
#endif

// Analysis algorithm class declaration
class TwoElectrons : public EL::Algorithm
{
   // put your configuration variables here as public variables.
   // that way they can be set directly from CINT and python.

public:
   enum Matching { 
      manualDR = 0,
      autoOld = 1,
      MCTruthClass = 2
    };

   // CONFIG variables
   std::string submitDir; 
   Matching matchMode = MCTruthClass;
   Bool_t decay = true;
   Bool_t testRunDirect = false;
   std::string outputName;
   
   const Double_t ZMASS = 91.2;


   
private:  
   // TTree to store data
   TTree *t_data; //!
   EEEvent *ev; //!
   Electron *elec1; //!
   Electron *elec2; //!

   // Data for truth info
   TruthElectron *tElec1; //!
   TruthElectron *tElec2; //!

   // Statistics
   Int_t n_TotalEventCount; //!
   Int_t n_InterestingEventCount; //!
   Int_t n_OSEventCount; //!   // Opposite sign event count
   Int_t n_SSEventCount; //!   // Same sign event count
   Double_t d_HighestPt; //!
   Double_t d_HighestMass; //!
   
   std::ofstream f_decayTree; //!
   
   Int_t electronCount;
   
   // Good runs list
#ifdef GRL

   GoodRunsListSelectionTool *m_grl; //!
#endif 

   // Electron selection
   AsgElectronLikelihoodTool* m_LHToolTight2015; //!
   AsgElectronLikelihoodTool* m_LHToolMedium2015; //!
   AsgElectronLikelihoodTool* m_LHToolLoose2015; //!
   
   // MCTruthClassifier
   MCTruthClassifier* m_truthClassifier; //!
   

private:
   // float cutValue;
   void SaveTagProbeInfo(xAOD::Electron*, Electron*);
   
   // Saves all the relevant truth information of an event
   void SaveMCInfo(const xAOD::Electron* e1, const xAOD::Electron* e2);
   
   // Saves all the relevant truth information of an electron
   void SaveMCInfo(const xAOD::TruthParticle* tp, TruthElectron *tElec);

   // Saves all truth information of event using MCTruthClassifier
   void SaveMCTruthClassifier(const xAOD::Electron* e1, const xAOD::Electron* e2);
   
   const xAOD::TruthParticle* SaveMCTruthClassifier(const xAOD::Electron* e1, TruthElectron *tElec);
   
   void PrintResults(std::ostream& );
   
   bool PassedEventCut();
   
   EL::StatusCode GetTruthParticle(const xAOD::Electron*, xAOD::TruthParticle*&);

   // Outputs decay tree to file
   void PrintDecayTree(const xAOD::Electron* e, std::ostream& );

   // variables that don't get filled at submission time should be
   // protected from being send from the submission node to the worker
   // node (done by the //!)
public:
  // Tree *myTree; //!




  // this is a standard constructor
  TwoElectrons ();

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();

  // this is needed to distribute the algorithm to the workers
  ClassDef(TwoElectrons, 1);
};

#endif
