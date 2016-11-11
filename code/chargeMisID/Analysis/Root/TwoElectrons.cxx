/////// About TwoElectrons.cxx ///////
/* Written by Gabriel Gallardo, The University of Hong Kong
 *
 * August 2015
 * 
 * This algorithm extracts electrons from a given set of events and places relevant info
 * in EEEvent and Electron classes.
 *
 * Information from AsgElectronLikelihoodTool are saved in the electron class
 *
 *
 *
 */

#include <Analysis/TwoElectrons.h>
#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/tools/Message.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/Egamma.h"
#include <Analysis/EEEvent.h>
#include <iostream>
#include "xAODEventInfo/EventInfo.h"
#include <TSystem.h>
#include <queue>
#include "xAODTruth/TruthParticle.h"
#include "xAODTruth/xAODTruthHelpers.h"
#include <TFile.h>
#include <utility> // For std::pair< , >

#define EL_RETURN_CHECK( CONTEXT, EXP )                     \
   do {                                                     \
      if( ! EXP.isSuccess() ) {                             \
         Error( CONTEXT,                                    \
                XAOD_MESSAGE( "Failed to execute: %s" ),    \
                #EXP );                                     \
         return EL::StatusCode::FAILURE;                    \
      }                                                     \
   } while( false )


// this is needed to distribute the algorithm to the workers
ClassImp(TwoElectrons)




typedef ElementLink< xAOD::TruthParticleContainer > TruthLink;
static SG::AuxElement::Accessor< TruthLink > acc_truthLink("truthParticleLink"); // ID track, electron
static SG::AuxElement::Accessor< int > acc_truthType("truthType");
static SG::AuxElement::Accessor< int > acc_truthOrig("truthOrigin");


// Initialize all pointers to electron and EEevent objects, set stat counters to 0
TwoElectrons :: TwoElectrons ()
{
   elec1 = new Electron();
   elec2 = new Electron();
   ev = new EEEvent();
   
   tElec1 = new TruthElectron();
   tElec2 = new TruthElectron();
   
   n_TotalEventCount = 0;
   d_HighestPt =0;
   n_InterestingEventCount = 0;
   n_SSEventCount = 0;
   n_OSEventCount = 0;
}


// job.useXAOD()
EL::StatusCode TwoElectrons :: setupJob (EL::Job& job)
{
  job.useXAOD();
  EL_RETURN_CHECK( "setupJob()", xAOD::Init() );
  
  return EL::StatusCode::SUCCESS;
}


// Initialize tree for taking data
EL::StatusCode TwoElectrons :: histInitialize ()
{
   // Here you do everything that needs to be done at the very
   // beginning on each worker node, e.g. create histograms and output
   // trees.  This method gets called before any input files are
   // connected.

   if(testRunDirect){
      t_data = new TTree("data", "Electron event data");
      t_data->Branch("Event", "EEEvent", ev);
      wk()->addOutput(t_data);
      return EL::StatusCode::SUCCESS;
   }

   // Initialize TTree
   TFile *outputFile = wk()->getOutputFile (outputName);
   t_data = new TTree("data", "Electron event data");
   t_data->SetDirectory (outputFile);

   t_data->Branch("Event", "EEEvent", ev);
   //wk()->addOutput(t_data);

   return EL::StatusCode::SUCCESS;
}


// Does nothing
EL::StatusCode TwoElectrons :: fileExecute ()
{
   // Here you do everything that needs to be done exactly once for every
   // single file, e.g. collect a list of all lumi-blocks processed
   return EL::StatusCode::SUCCESS;
}


// Does nothing
EL::StatusCode TwoElectrons :: changeInput (bool /*firstFile*/)
{
   // Here you do everything you need to do when we change input files,
   // e.g. resetting branch addresses on trees.  If you are using
   // D3PDReader or a similar service this method is not needed.
   return EL::StatusCode::SUCCESS;
}


// Initialize GRL (currently inactive), electron selection tools, MCTruthClassifier, and output stream for DecayTree.txt
EL::StatusCode TwoElectrons :: initialize ()
{
   //xAOD::TEvent* event = wk()->xaodEvent(); // Uncomment if necessary

#ifdef GRL 
   // === GOOD RUNS LIST INITIALIZATION ===
   // Update GRL file path
   m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
   const char* GRLFilePath = "$ALRB_TutorialData/data12_8TeV.periodAllYear_DetStatus-v61-pro14-02_DQDefects-00-01-00_PHYS_StandardGRL_All_Good.xml";
   const char* fullGRLFilePath = gSystem->ExpandPathName (GRLFilePath);
   std::vector<std::string> vecStringGRL;
   vecStringGRL.push_back(fullGRLFilePath);
   EL_RETURN_CHECK("initialize()",m_grl->setProperty( "GoodRunsListVec", vecStringGRL));
   EL_RETURN_CHECK("initialize()",m_grl->setProperty("PassThrough", false)); // if true (default) will ignore result of GRL and will just pass all events
   EL_RETURN_CHECK("initialize()",m_grl->initialize());
#endif

   // === ELECTRON SELECTION TOOLS INITIALIZATION ===
   m_LHToolTight2015  = new AsgElectronLikelihoodTool ("m_LHToolTight2015");
   m_LHToolMedium2015 = new AsgElectronLikelihoodTool ("m_LHToolMedium2015"); 
   m_LHToolLoose2015  = new AsgElectronLikelihoodTool ("m_LHToolLoose2015");

   // initialize the primary vertex container for the tool to have access to the number of vertices used to adapt cuts based on the pileup
   EL_RETURN_CHECK("initialize()", m_LHToolTight2015->setProperty("primaryVertexContainer","PrimaryVertices"));
   EL_RETURN_CHECK("initialize()", m_LHToolMedium2015->setProperty("primaryVertexContainer","PrimaryVertices"));
   EL_RETURN_CHECK("initialize()", m_LHToolLoose2015->setProperty("primaryVertexContainer","PrimaryVertices"));
   EL_RETURN_CHECK("initialize()", m_LHToolLoose2015->setProperty("primaryVertexContainer","PrimaryVertices"));

   // define the config files
   std::string confDir = "ElectronPhotonSelectorTools/offline/mc15_20150408/";
   EL_RETURN_CHECK("initialize()", m_LHToolTight2015->setProperty("ConfigFile",confDir+"ElectronLikelihoodTightOfflineConfig2015.conf"));
   EL_RETURN_CHECK("initialize()", m_LHToolMedium2015->setProperty("ConfigFile",confDir+"ElectronLikelihoodMediumOfflineConfig2015.conf"));
   EL_RETURN_CHECK("initialize()", m_LHToolLoose2015->setProperty("ConfigFile",confDir+"ElectronLikelihoodLooseOfflineConfig2015.conf"));

   EL_RETURN_CHECK("initialize()", m_LHToolTight2015->initialize());
   EL_RETURN_CHECK("initialize()", m_LHToolMedium2015->initialize());
   EL_RETURN_CHECK("initialize()", m_LHToolLoose2015->initialize());

   // === MC-TRUTH-CLASSIFIER INITIALIZATION ===
   m_truthClassifier = new MCTruthClassifier("m_truthClassifier"); // I think this is the proper call.

	if(decay) f_decayTree.open(submitDir+"/DecayTree.txt");
	
   return EL::StatusCode::SUCCESS;
}


// Data analysis and saving
EL::StatusCode TwoElectrons :: execute ()
{
   xAOD::TEvent* event = wk()->xaodEvent();
   n_TotalEventCount++;
   
   // Setup pointer to electron container
   const xAOD::ElectronContainer* electrons = 0;
   EL_RETURN_CHECK( "execute()", event->retrieve( electrons, "Electrons" ));

   // Check if event is MC data
   const xAOD::EventInfo* eventInfo = 0;
   EL_RETURN_CHECK("execute()",event->retrieve( eventInfo, "EventInfo"));  
   bool isMC = false;
   if(eventInfo->eventType( xAOD::EventInfo::IS_SIMULATION ) ) isMC = true; 
      
   // GRL: Reject if event does not pass GRL
#ifdef GRL
   if(!isMC){ 
      if(!m_grl->passRunLB(*eventInfo)){
         return EL::StatusCode::SUCCESS; 
      }
   } 
#endif

   // Temporary electron pointers
   xAOD::Electron* e1 = 0;
   xAOD::Electron* e2 = 0;
   
   // Initialize temporary storage objects;
   elec1->Set();
   elec2->Set();
   tElec1->Set();
   tElec2->Set();

   // Count of electrons in this event meeting the selection criteria
   electronCount = 0;

   //=============================//
   // Preliminary Event Selection //
   //=============================//
   for(auto e : *electrons){
      // Electron must pass loose criteria of ElectronLikelihoodTool
      if(!m_LHToolLoose2015->accept(e)) continue;

      Double_t d_PtElectron = e->pt()*0.001;  // Divide by 1000 to obtain pt in GeV
      Double_t d_EtaElectron = e->eta();
       
		// Preliminary selection criteria: Pt > 20GeV  && |eta| < 2.5
      if(!(
         d_PtElectron >= 20 &&
         fabs(d_EtaElectron) < 2.5
         )) continue;
      electronCount++;
      if(d_PtElectron > d_HighestPt) d_HighestPt = d_PtElectron;  
    
      // Save pointers. For cases where more than 2 meet the criteria, the two with mass closest to ZMASS are selected
      if(electronCount == 1) e1 = const_cast<xAOD::Electron*>(e); 
      else if(electronCount == 2) e2 = const_cast<xAOD::Electron*>(e);
      else {
      	/* // Choose highest Pt
      	if (e1->pt() < e2->pt() ){
      	   if (e1->pt() < e->pt()){
      		   e1 = const_cast<xAOD::Electron*>(e);
      		}
      	} else if (e2->pt() < e->pt()){
      	   e2 = const_cast<xAOD::Electron*>(e);
      	}
      	//*/
      	
      	//* Choose mass closest to 91.2 GeV
			// Find absolute difference between invariant mass of pair and mass of Z boson
			double dInvMass12 = fabs( ( e1->p4() + e2->p4() ).M() / 1000 - ZMASS );
			double dInvMass13 = fabs( ( e1->p4() + e->p4() ).M() / 1000 - ZMASS );
			double dInvMass32 = fabs( ( e->p4() + e2->p4() ).M() / 1000 - ZMASS );   
			
			if(dInvMass13 < dInvMass32){
				if(dInvMass13 < dInvMass12){
					e2 = const_cast<xAOD::Electron*>(e);
				}
			}
			else if(dInvMass32 < dInvMass12){
				e1 = const_cast<xAOD::Electron*>(e);
			}
			//*/
      }
      
   }
   // End loop over each electron
   
   //==================================================================//
   // Retrieving information from electron objects and saving to TTree //
   //==================================================================//
   if(electronCount < 2){ 
   	// Nothing to do if less than 2 eligible electrons were found
      return EL::StatusCode::SUCCESS;
   }
   else{
      // Is it a same sign event?
      Bool_t SSevent = ( e1->charge() == e2->charge() );

		// Calculate invariant mass of the two electrons (in GeV)
		auto v_pe1 = e1->p4();
		auto v_pe2 = e2->p4();
      double InvMass = ( v_pe1 + v_pe2 ).M() / 1000;                        

		// Calculate dPhi
      Double_t dPhi = v_pe1.DeltaPhi(v_pe2);

      // Electron objects to store charge, eta, pt
      elec1->Set(e1->charge(), e1->eta(), e1->pt()*0.001);
      elec2->Set(e2->charge(), e2->eta(), e2->pt()*0.001);
      
      // Save electron pass-cut info
      elec1->SetPassCut(m_LHToolLoose2015->accept(e1), m_LHToolMedium2015->accept(e1), m_LHToolTight2015->accept(e1));
      elec2->SetPassCut(m_LHToolLoose2015->accept(e2), m_LHToolMedium2015->accept(e2), m_LHToolTight2015->accept(e2));
      
      // Save Tag-and-probe information (detector info)
      SaveTagProbeInfo(e1, elec1);
      SaveTagProbeInfo(e2, elec2);

		// Save event info
      ev->Set(*elec1, *elec2, InvMass, dPhi, SSevent); 
      
      // Do not save event into tree if event cut was not passed
      if(!PassedEventCut()) return EL::StatusCode::SUCCESS;

      if(matchMode != MCTruthClass){
         // Save MC truth info
         if (isMC) SaveMCInfo(e1, e2);
         else ev->SetTruth(false);
      } else {
         SaveMCTruthClassifier(const_cast<xAOD::Electron*>(e1), const_cast<xAOD::Electron*>(e2));
      }

      // Fill tree
      t_data->Fill();
      
      // Update statistics
      if (InvMass > d_HighestMass) d_HighestMass = InvMass;
      n_InterestingEventCount++;
      if(SSevent) n_SSEventCount++;
      else n_OSEventCount++;

      // Print results
      PrintResults(std::cout); // Print to std::cout


      // Print decay tree for events which are not Z->ee  
      if(isMC && (ev->ftSSevent || !(ev->ftZEEevent))){
         if(tElec1->fParticleOrigin != 23 || tElec2->fParticleOrigin != 23){
            if(ev->ftSSevent) std::cout << "SS " ;
            if(!ev->ftZEEevent) std::cout << "NotZ";
            std::cout << std::endl;
            PrintDecayTree(e1, std::cout);
            PrintDecayTree(e2, std::cout);
            std::cout << std::endl;
            
            if(decay){
               PrintResults(f_decayTree);
               if(ev->ftSSevent) f_decayTree << "SS " ;
               if(!ev->ftZEEevent) f_decayTree << "NotZ";
               std::cout << std::endl;
               PrintDecayTree(e1, f_decayTree);
               PrintDecayTree(e2, f_decayTree);
               f_decayTree << std::endl;
            }
         }
      }


   }
   return EL::StatusCode::SUCCESS;
}


// Does nothing
EL::StatusCode TwoElectrons :: postExecute ()
{
   return EL::StatusCode::SUCCESS;
}


// Output basic statistics to terminal, clean up pointers.
EL::StatusCode TwoElectrons :: finalize ()
{
   //xAOD::TEvent* event = wk()->xaodEvent(); // (Un)comment if necessary
   
   std::cout << "Highest recorded pt = " << d_HighestPt << "GeV" << std::endl 
            << "Highest recorded mass = " << d_HighestMass << "GeV" << std::endl
            << "No. of events printed = " << n_InterestingEventCount <<" / " << n_TotalEventCount << std::endl
            << "SS events: " << n_SSEventCount << "\tOSevents: " << n_OSEventCount << std::endl
            <<"MisID rate: " << n_SSEventCount/ (double) n_InterestingEventCount * 100.0 << "% " << std::endl << std::endl;
            
   
   // Cleanup 
   if (elec1) { delete elec1; elec1 = 0; }
   if (elec2) { delete elec2; elec2 = 0; }
   if (ev) { delete ev; ev = 0; }
   if (tElec1) { delete tElec1; tElec1 = 0; }
   if (tElec2) { delete tElec2; tElec2 = 0; }
   if (m_truthClassifier) { delete m_truthClassifier; m_truthClassifier = 0; }
#ifdef GRL
   if (m_grl) { delete m_grl;  m_grl = 0; }
#endif
	
   if(decay) f_decayTree.close();

   return EL::StatusCode::SUCCESS;
}


// Does nothing
EL::StatusCode TwoElectrons :: histFinalize ()
{

   return EL::StatusCode::SUCCESS;
}


// Saves d0, d0error, nBHits, nSiHits, trigger matched condition (the latter to be implemented)
void TwoElectrons :: SaveTagProbeInfo(xAOD::Electron* e, Electron *elec)
{
   // Calculate additional values of electron for Tag/Probe analysis
   const xAOD::TrackParticle* trk = e->trackParticle(); // TriggerMatched condition?
   if (trk){
   	// Save d0
      Double_t d0 = trk->d0();
      
   	// Save d0error
      xAOD::ParametersCovMatrix_t cov = trk->definingParametersCovMatrix();
      Double_t d0error = sqrt(cov(0,0));
      
      // Save nBhits
      uint8_t BHits = 0;
      trk->summaryValue(BHits,xAOD::numberOfBLayerHits);
      
      // Save nSiHits
      uint8_t SiHits = 0;
      trk->summaryValue(SiHits, xAOD::numberOfSCTHits);
      
      // Save isTriggerMatched?
      bool triggerMatched = true;

      // Insert trigger matched checking code here
      
      elec->SetInfo(d0, d0error, SiHits, BHits, triggerMatched); 
   }
   else elec->SetInfo(0,0,0,0,false);
}


// Uses MCTruthClassifier to save truth information
void TwoElectrons :: SaveMCTruthClassifier(const xAOD::Electron* e1, const xAOD::Electron* e2)
{
   auto tE1 = SaveMCTruthClassifier(e1, tElec1);
   auto tE2 = SaveMCTruthClassifier(e2, tElec2);

   Bool_t isZevent = false;
   Bool_t SSevent = false;
   Double_t InvMass = 0;
   Double_t dPhi = 0;

   if(tE1 && tE2){
      if(tElec1->fIsElectron && tElec1->fFromZ && tElec2->fIsElectron && tElec2->fFromZ) isZevent = true;
      SSevent = (tElec1->fCharge == tElec2->fCharge);

      auto v_pe1 = tE1->p4();
      auto v_pe2 = tE2->p4();
      InvMass = ( v_pe1 + v_pe2 ).M() / 1000;                         

      dPhi = v_pe1.DeltaPhi(v_pe2);
   }

   ev->SetTruth(*tElec1, *tElec2, InvMass, dPhi, SSevent, isZevent);
}

const xAOD::TruthParticle* TwoElectrons :: SaveMCTruthClassifier(const xAOD::Electron* elec, TruthElectron* tElec)
{
   tElec->Set();

   std::pair<MCTruthPartClassifier::ParticleType, MCTruthPartClassifier::ParticleOrigin> res;
   res = m_truthClassifier->particleTruthClassifier(elec);

   tElec->fTruthType = res.first;
   tElec->fTruthOrigin = res.second;

   auto truthP = m_truthClassifier->getGenPart();
   if(truthP == 0) return truthP; // guard in case truth particle is not found

   tElec->Set(truthP->charge(), truthP->eta(), truthP->pt()/1000.0);

   // if it is a background electron, check its origin
   if(tElec->fTruthType == 4){
      res = m_truthClassifier->checkOrigOfBkgElec(truthP);
      tElec->fTruthBkgType = res.first;
      tElec->fTruthBkgOrigin = res.second;
      
      auto mother = m_truthClassifier->getBkgElecMother();
      if (mother == 0) ;
      else tElec->fTruthBkgCharge = mother->charge();
   }

   if (tElec->fTruthType >= 1 && tElec->fTruthType <=4) tElec->fIsElectron = true;
   if (tElec->fTruthOrigin == 13 || tElec->fTruthBkgOrigin == 13) tElec->fFromZ = true;

   tElec->fDeltaR = m_truthClassifier->getdeltaRMatch();

   return truthP;
}

// Saves all the relevant truth information of an event
void TwoElectrons :: SaveMCInfo(const xAOD::Electron* e1, const xAOD::Electron* e2)
{
   Bool_t isZevent = true;
   
   // Save for electron 1
	xAOD::TruthParticle *tE1 = 0;
	GetTruthParticle(e1, tE1);
	if(tE1){
	   SaveMCInfo(tE1, tElec1); 
	   tElec1->fDeltaR = (tE1->p4()).DeltaR(e1->p4());
	}
	else {
      tElec1->Set();
	   isZevent = false;	
	}
   
   // Save for electron 2
	xAOD::TruthParticle *tE2 = 0;
	GetTruthParticle(e2, tE2);
	if(tE2){
	   SaveMCInfo(tE2, tElec2);
	   tElec2->fDeltaR = (tE2->p4()).DeltaR(e2->p4());
	}
	else {
	   tElec2->Set();
	   isZevent = false;
	         
	}
   
   Bool_t SSevent = 0;
   Double_t InvMass = 0;
   Double_t dPhi = 0;
   
   if(tE1 && tE2){
      // if the parents are not the same or if the parent is not a Z => not Z event
      if(!(tE1->parent(0) == tE2->parent(0)) || !tE1->parent(0)->isZ()){
         isZevent = false;
      }
      
      SSevent = (tElec1->fCharge == tElec2->fCharge);
		
		auto v_pe1 = tE1->p4();
		auto v_pe2 = tE2->p4();
      InvMass = ( v_pe1 + v_pe2 ).M() / 1000;                         

      dPhi = v_pe1.DeltaPhi(v_pe2);
   }
   
   ev->SetTruth(*tElec1, *tElec2, InvMass, dPhi, SSevent, isZevent);
}


// Saves all the relevant truth information of an electron
void TwoElectrons :: SaveMCInfo(const xAOD::TruthParticle* tp, TruthElectron *tElec){
	tElec->Set();
	if (!tp) return;

	Bool_t isElectron = tp->isElectron();
	Int_t partType = tp->pdgId();
	Int_t partOrig = tp->parent(0)->pdgId();
	
	tElec->Set(tp->charge(), tp->eta(), tp->pt() / 1000);
	tElec->SetTruth(isElectron, partType, partOrig);

	if(tp->nParents()==1 && tp->parent(0)->isZ()) tElec->fFromZ = true;
}


/** PassedEventCut()
 * Apply secondary event cuts here
 * Returns whether the event passed the cut
 */
bool TwoElectrons :: PassedEventCut(){
   const double d0cut = 3;
   if(
      // Invariant mass cut: 70 < InvMass < 110 (90 ± 20 GeV)
      ev->fInvMass >= 70 &&
      ev->fInvMass <= 110 &&
      
      // d0/d0error cut
      fabs(elec1->fD0/elec1->fD0Error) <= d0cut &&
      fabs(elec2->fD0/elec2->fD0Error) <= d0cut
   ) return true;
   else return false;
}


/** PrintResults()
 * Prints results of a single event to specified ostream stream
 * Call on std::cout for terminal
 * Call on filestream for output to file
 * Used by execute() to output to terminal, and output to f_decayTree
 */
void TwoElectrons :: PrintResults(std::ostream& stream){
   stream << "#" << n_TotalEventCount << ": " << "\t"
           << "Electron count: " << electronCount << std::endl;

   ev->PrintInfo(stream);
}


/** GetTruthParticle()
 * Assigns truth particle of e to truth
 * Modify here to choose manual matching (dR matching) or auto matching (truthLink).
 *
 * This function returns a StatusCode and not a pointer because it requires use of the 
 * event->retrieve() method, which necessarily is enclosed inside the EL_RETURN_CHECK()
 * macro, which itself returns an EL::StatusCode. The workaround is to pass the pointer
 * to be set as a reference.
 */
EL::StatusCode TwoElectrons :: GetTruthParticle(const xAOD::Electron* e, xAOD::TruthParticle* &truth){
   xAOD::TEvent* event = wk()->xaodEvent();
   
   xAOD::TruthParticle *tp = 0;
   
   //===============================//
   // Manual match                  //
   //===============================//
   
   //*
   if (matchMode == manualDR){
      // Set up truth particle container
      const xAOD::TruthParticleContainer* tContainer = 0;
      EL_RETURN_CHECK("GetTruthParticle", event->retrieve( tContainer, "TruthParticles" ));
   
      // Momentum vectors
      TLorentzVector p_tp(0,0,0,0);
      TLorentzVector p_e = e->p4();
      bool firstTry = true;
   
      // Find electron with smallest dR < 0.1
      for(auto particle : *tContainer){
         if(!(particle->isElectron())) continue;
   
         TLorentzVector p_particle = particle->p4();
         Double_t dR = p_e.DeltaR(p_particle);
         if (dR > 0.1) continue;
         
         if (firstTry || dR < p_e.DeltaR(p_tp)){
            tp = const_cast<xAOD::TruthParticle*>(particle);
            p_tp = p_particle;
            firstTry = false;
         }
         else continue;
      }
   }   
   
   //===============================//
   // Auto match                    //
   //===============================//
   else if (matchMode == autoOld){
      auto tLink = acc_truthLink(*e);
      if(!tLink.isValid()) return EL::StatusCode::SUCCESS;
      tp = const_cast<xAOD::TruthParticle*>(*tLink);
   }

   // Find most parent electron from Z
   if (!tp) return EL::StatusCode::SUCCESS; // return if not matched
   while(tp->nParents() == 1){
      if (!(tp->parent(0)->isElectron() || tp->parent(0)->isPhoton())) break; // break if parent is neither electron nor photon
      if(tp->parent(0)->nParents() == 0) break;
      tp = const_cast<xAOD::TruthParticle*>(tp->parent(0));
   }
   
   truth = tp;
   return EL::StatusCode::SUCCESS;
}


//===================================================//
// Prints decay tree to file attached at f_decayTree //
//===================================================//



/** PrintDecayTree()
 * Prints decay path of specified electron e
 * Marks the truth particle matched by GetTruthParticle() with '*'
 * Marks the truth particle matched by truthLink with 'm'
 */
void TwoElectrons :: PrintDecayTree(const xAOD::Electron* e, std::ostream& stream){
   // struct used for implementation of PrintDecayTree()
   struct decayParticle {
      xAOD::TruthParticle* p;
      Int_t particleNumber;
      Int_t motherNumber;
   };

	stream << "Decay tree: " << std::endl;
	stream << "Orig\tNo.\tPdgID\tStatus\tBarcode\tMother\tDaughters\n";
	std::queue<decayParticle> decays; // Queue (FIFO) for particles to be examined 
	Int_t nParticle = 1; // Particle counter
	
   /*
	auto tLink = acc_truthLink(*e); // for matching later
	xAOD::TruthParticle* tE = 0;
	GetTruthParticle(e, tE);
   */
   auto res = m_truthClassifier->particleTruthClassifier(e);
   auto tE = m_truthClassifier->getGenPart();


	if (tE == 0){
		stream << "No decay tree. " << std::endl;
		return;
	}

   // Find ancestors of particle (decay path)
	auto mother = m_truthClassifier->getMother();
   if (mother == 0) {
      stream << "No mother found. " << std::endl;
      return;
   }

   if(res.second == 4) {
      mother = m_truthClassifier->getBkgElecMother();
      if (mother->pdgId() == 22){
         m_truthClassifier->particleTruthClassifier(mother);
         mother = m_truthClassifier->getMother();
      }
   }
	
	// Save particulars of mother particle
	decayParticle tempParticle;
	tempParticle.p = const_cast<xAOD::TruthParticle*>(mother);
	tempParticle.particleNumber = nParticle;
	tempParticle.motherNumber = 0;
	decays.push(tempParticle);
	nParticle++;
	
	// Go through queue and output properties of particles
	while( !decays.empty() ){
		auto particle = decays.front(); // Access first element in queue
		
		if(particle.p == tE) stream << "*";
		stream<<" \t";
		
		stream << particle.particleNumber << "\t";
		
		// Output names of selected particles
		int pdgId = particle.p->pdgId();
		if (pdgId == 11){ stream << "e-" << "\t";} 
		else if (pdgId == -11){ stream << "e+" << "\t";}
		else if (pdgId == 23){ stream << "Z" << "\t";}
		else if (pdgId == 22){ stream << "photon" << "\t";}
		else if (pdgId == 24){ stream << "W+" << "\t";}
		else if (pdgId == -24){ stream << "W-" << "\t";}
		else if (pdgId == 111){ stream << "pion-0" << "\t";}
		else if (pdgId == 2212){ stream << "proton" << "\t";}
		else { stream << pdgId << "\t";}
		   
		stream << particle.p->status() << "\t";
		stream << particle.p->barcode() << "\t";
		stream << particle.motherNumber << "\t";
		
		// Save daughter particles into queue for output later
		for(unsigned int i = 0; i < particle.p->nChildren(); i++){
			tempParticle.p = const_cast<xAOD::TruthParticle*>(particle.p->child(i));
			tempParticle.particleNumber = nParticle;
			tempParticle.motherNumber = particle.particleNumber;
			decays.push(tempParticle);
			stream << nParticle << "\t";
			nParticle++;
		}
		
		decays.pop(); // Pop first element in queue
		stream << std::endl;
		
	}

	
}
