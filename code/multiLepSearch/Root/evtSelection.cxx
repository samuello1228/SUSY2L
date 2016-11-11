#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <multiLepSearch/evtSelection.h>
// #include <multiLepSearch/MT2.h>
#include <multiLepSearch/anaHelper.h>
#include "xAODMissingET/MissingETContainer.h"
#include "xAODMissingET/MissingETAuxContainer.h"
// #include "xAODTruth/TruthParticleContainer.h"
// #include "xAODTruth/TruthParticle.h"
// #include "xAODMuon/Muon.h"
#include "GoodRunsLists/GoodRunsListSelectionTool.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODCutFlow/CutBookkeeperContainer.h"
#include "xAODCutFlow/CutBookkeeper.h"
#include "PathResolver/PathResolver.h"
#include "xAODRootAccess/tools/Message.h"
// #include "MuonSelectorTools/MuonSelectionTool.h"

#include "PATInterfaces/SystematicVariation.h" 
#include "PATInterfaces/SystematicsUtil.h"

#include <TError.h>
#include <algorithm>
#include <vector>
#include <multiLepSearch/mCHECK.h>


using namespace xAOD;

// static SG::AuxElement::Decorator<char> dec_quality("quality");
static SG::AuxElement::Decorator<char> dec_baseline("baseline");
static SG::AuxElement::Decorator<char> dec_signal("signal");
static SG::AuxElement::Decorator<char> dec_passOR("passOR");
static SG::AuxElement::Decorator<char> dec_bad("bad");
static SG::AuxElement::Decorator<char> dec_bjet_loose("bjet_loose");
static SG::AuxElement::Decorator<char> dec_bjet("bjet");
static SG::AuxElement::Decorator<char> dec_cosmic("cosmic");
typedef ElementLink< xAOD::TruthParticleContainer > TruthLink;
static SG::AuxElement::Accessor< int > acc_truthType("truthType");
static SG::AuxElement::Accessor< int > acc_truthOrig("truthOrigin");
// static SG::AuxElement::Accessor< float > acc_truthProb("truthMatchProbability"); // only ID track
static SG::AuxElement::Accessor< TruthLink > acc_truthLink("truthParticleLink"); // ID track, electron
// this is needed to distribute the algorithm to the workers
ClassImp(evtSelection)


const float iGeV = 0.001;
const float GeV = 1000;

evtSelection :: evtSelection (std::string name)
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
  m_name = name;
  //m_susyEvt = new susyEvts();
  m_vxTrkNMin = 2;
  m_vxTrkPtMin = 2*GeV;
  outputName = "test.root";
  outputTreeName = "evt2l";
  m_grl = 0;
  nLepCutExactly = 2;
  nLepCutMin = 2;
  nJetCutExactly = 2;
  nJetCutMin = 2;

  ElPtCut = 25000;
  Eld0SigCut = 5;
  Elz0Cut = 0.5;

  MuPtCut = 25000;
  Mud0SigCut = 3;
  Muz0Cut = 0.5;

  JetPtCut = 20000;
  JetEtaCut = 2.8;
  JetJvtCut = 0.64;

  m_isMC = 0;
  m_doSys = 0;
}

EL::StatusCode evtSelection :: setupJob (EL::Job& job)
{
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.
  job.useXAOD ();
  xAOD::Init( m_name.c_str() ).ignore(); // call before opening first file

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode evtSelection :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.
  TFile *outputFile = wk()->getOutputFile(outputName);

  m_hCutFlow = new TH1F("hCutFlow", "cut flow", 20, 0, 20);
  m_hCutFlow->SetDirectory(outputFile);

  m_hTrigs = new TH1F("hTrigs", "n pass trigger", trigNames.size(), 0, trigNames.size());
  for(unsigned int i=0; i<trigNames.size(); i++)
  {
    m_hTrigs->GetXaxis()->SetBinLabel(i+1,trigNames[i].c_str());
  }
  m_hTrigs->SetDirectory(outputFile);
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode evtSelection :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  m_event = wk()->xaodEvent();

  //https://twiki.cern.ch/twiki/bin/view/AtlasProtected/AnalysisMetadata
  //get the MetaData tree once a new file is opened, with
  TTree *MetaData = dynamic_cast<TTree*>(wk()->inputFile()->Get("MetaData"));
  if (!MetaData) {
    Error("fileExecute()", "MetaData not found! Exiting.");
    return EL::StatusCode::FAILURE;
  }
  MetaData->LoadTree(0);

  //check if file is from a DxAOD
  bool m_isDerivation = !MetaData->GetBranch("StreamAOD");

  if(m_isDerivation){
   const xAOD::CutBookkeeperContainer* incompleteCBC = nullptr;
   if(!m_event->retrieveMetaInput(incompleteCBC, "IncompleteCutBookkeepers").isSuccess()){
      Error("fileExecute()","Failed to retrieve IncompleteCutBookkeepers from MetaData! Exiting.");
      return EL::StatusCode::FAILURE;
    }
    if ( incompleteCBC->size() != 0 ) {
      Warning("fileExecute()","Found incomplete Bookkeepers! Check file for corruption.");
      //return EL::StatusCode::FAILURE;
    }
    // Now, let's find the actual information
    const xAOD::CutBookkeeperContainer* completeCBC = 0;
    if(!m_event->retrieveMetaInput(completeCBC, "CutBookkeepers").isSuccess()){
      Error("fileExecute()","Failed to retrieve CutBookkeepers from MetaData! Exiting.");
      return EL::StatusCode::FAILURE;
    }

    // First, let's find the smallest cycle number,
    // i.e., the original first processing step/cycle
    int minCycle = 10000;
    for ( auto cbk : *completeCBC ) {
      if ( minCycle > cbk->cycle() ){ minCycle = cbk->cycle(); }
    }
    // Now, let's actually find the right one that contains all the needed info...
    const xAOD::CutBookkeeper* allEventsCBK=0;
    for ( auto cbk :  *completeCBC ) {
      if ( minCycle == cbk->cycle() && cbk->name() == "AllExecutedEvents" ){
      allEventsCBK = cbk;
      break;
      }
    }
    uint64_t nEventsProcessed  = allEventsCBK->nAcceptedEvents();
    //double sumOfWeights        = allEventsCBK->sumOfEventWeights();
    //double sumOfWeightsSquared = allEventsCBK->sumOfEventWeightsSquared();
    m_hCutFlow->Fill("AOD", nEventsProcessed);
  }
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode evtSelection :: changeInput (bool /*firstFile*/)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode evtSelection :: initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.
  m_event = wk()->xaodEvent();
  m_store = wk()->xaodStore();

  // as a check, let's see the number of events in our xAOD
  Info("initialize()", "Number of events = %lli", m_event->getEntries() ); // print long long int

  // count number of events
  m_eventCounter = 0;

  m_objTool = new ST::SUSYObjDef_xAOD("SUSYObjDef_xAOD");
  ///////////////////////////////////////////////////////////////////////////////////////////
  // Configure the SUSYObjDef instance
  // CHECK(m_objTool->setProperty("IsData",1) ) ;
  // CHECK(m_objTool->setProperty("IsAtlfast",0) );
  // Info("initialize()", "DataSource = %d", m_isMC);

  ST::ISUSYObjDef_xAODTool::DataSource ds = static_cast<ST::ISUSYObjDef_xAODTool::DataSource>(m_isMC); 
  CHECK(m_objTool->setProperty("DataSource",ds));
  //CHECK(m_objTool->setProperty("DataSource",m_isMC));
//   CHECK(m_objTool->setProperty("EleId","TightLH") );
//   CHECK(m_objTool->setProperty("EleIdBaseline","MediumLH") );
//   //CHECK(m_objTool->setProperty("TauId","Tight") );
//   CHECK(m_objTool->setProperty("PhotonIsoWP","Cone20") );
// 
  //CHECK(m_objTool->setProperty("IsoWP","GradientLoose") );
  
  CHECK(m_objTool->setProperty( "PRWConfigFiles", PRW_confFiles));
  CHECK(m_objTool->setProperty( "PRWLumiCalcFiles", PRW_lcalcFiles));
  //CHECK(m_prwTool->setProperty( "DataScaleFactor", 1.0/1.16)); //already done by SUSYTool

  m_objTool->msg().setLevel( MSG::ERROR );
  //m_objTool->msg().setLevel( MSG::WARNING );
  //m_objTool->msg().setLevel( MSG::INFO );
  //m_objTool->msg().setLevel( MSG::VERBOSE );
  
  CHECK(m_objTool->initialize().isSuccess());
  
  m_LHToolTight2015    = new AsgElectronLikelihoodTool("m_LHToolTight2015");
  m_LHToolMedium2015   = new AsgElectronLikelihoodTool("m_LHToolMedium2015"); 
  m_LHToolLoose2015    = new AsgElectronLikelihoodTool("m_LHToolLoose2015");

  CHECK(m_LHToolTight2015->setProperty("primaryVertexContainer","PrimaryVertices"));
  CHECK(m_LHToolMedium2015->setProperty("primaryVertexContainer","PrimaryVertices"));
  CHECK(m_LHToolLoose2015->setProperty("primaryVertexContainer","PrimaryVertices"));

  std::string confDir = "ElectronPhotonSelectorTools/offline/mc15_20150224/";
  CHECK(m_LHToolTight2015->setProperty("ConfigFile",confDir+"ElectronLikelihoodTightOfflineConfig2015.conf"));
  CHECK(m_LHToolMedium2015->setProperty("ConfigFile",confDir+"ElectronLikelihoodMediumOfflineConfig2015.conf"));
  CHECK(m_LHToolLoose2015->setProperty("ConfigFile",confDir+"ElectronLikelihoodLooseOfflineConfig2015.conf"));

  CHECK(m_LHToolTight2015->initialize());
  CHECK(m_LHToolMedium2015->initialize());
  CHECK(m_LHToolLoose2015->initialize());

  if(!m_isMC){
    m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
    std::vector<std::string> vecStringGRL;
    vecStringGRL.push_back(PathResolverFindCalibFile(m_grlFile));
    CHECK(m_grl->setProperty( "GoodRunsListVec", vecStringGRL));
    CHECK(m_grl->setProperty("PassThrough", false)); // if true (default) will ignore result of GRL and will just pass all events
    CHECK(m_grl->initialize());
  }

  /// isolation tools
  m_isoTool = new CP::IsolationSelectionTool("isoTool");
  CHECK(m_isoTool->setProperty("MuonWP","Loose"));
  CHECK(m_isoTool->setProperty("ElectronWP","Loose"));
  CHECK(m_isoTool->initialize());

  CHECK(m_isoTool->addWP("Tight", xAOD::Type::Muon));
  CHECK(m_isoTool->addWP("Tight", xAOD::Type::Electron));

  CHECK(m_isoTool->addWP("Gradient", xAOD::Type::Muon));
  CHECK(m_isoTool->addWP("Gradient", xAOD::Type::Electron));

  CHECK(m_isoTool->addWP("GradientLoose", xAOD::Type::Muon));
  CHECK(m_isoTool->addWP("GradientLoose", xAOD::Type::Electron));

  CHECK(m_isoTool->addWP("LooseTrackOnly", xAOD::Type::Muon));
  CHECK(m_isoTool->addWP("LooseTrackOnly", xAOD::Type::Electron));

  CHECK(m_isoTool->addWP("FixedCutLoose", xAOD::Type::Muon));
  CHECK(m_isoTool->addWP("FixedCutLoose", xAOD::Type::Electron));

  CHECK(m_isoTool->addWP("FixedCutTightTrackOnly", xAOD::Type::Muon));
  CHECK(m_isoTool->addWP("FixedCutTightTrackOnly", xAOD::Type::Electron));

  CHECK(m_isoTool->addWP("FixedCutTight", xAOD::Type::Electron));


  /// muon selector tool
  m_muonSelTool = new CP::MuonSelectionTool("MuonSelector");
  CHECK(m_muonSelTool->setProperty( "MaxEta", 2.7 ));
  CHECK(m_muonSelTool->setProperty( "MuQuality", 3));
  CHECK(m_muonSelTool->initialize().isSuccess());

  //trigger /////////////////////
  //Use the one in SUSYTool
  //configTool = new TrigConf::xAODConfigTool("xAODConfigTool2");
  //configHandle = new ToolHandle<TrigConf::ITrigConfigTool>(configTool);
  //CHECK((*configHandle)->initialize());

  //trigDecTool = new Trig::TrigDecisionTool("TrigDecTool2");
  //CHECK(trigDecTool->setProperty("ConfigTool",*configHandle));
  //CHECK(trigDecTool->setProperty("TrigDecisionKey","xTrigDecision"));
  //CHECK(trigDecTool->initialize());
  start = true; //indicates the 1st loop
  ///////////////////////////////////////////

  const CP::SystematicRegistry& registry = CP::SystematicRegistry::getInstance();
  const CP::SystematicSet& recommendedSystematics = registry.recommendedSystematics(); // get list of recommended systematics
  m_sysList = CP::make_systematics_vector(recommendedSystematics); 

  TFile *outputFile = wk()->getOutputFile(outputName);
  susyEvts* m_susyEvt;

  // print list of recommended systematics and prepare output tree for them
  std::vector<CP::SystematicSet>::const_iterator sysListItr;
  std::cout << "======  Systematic list  ======" << std::endl;
  for (sysListItr = m_sysList.begin(); sysListItr != m_sysList.end(); ++sysListItr){
    if((*sysListItr).name()!="" && !m_doSys) continue;

    if((*sysListItr).name()=="") std::cout << "Nominal (no syst) "  << std::endl;
    else std::cout << "Systematic: " << (*sysListItr).name() << std::endl;

    m_susyEvt = new susyEvts();
    m_susyEvt->makeTree(outputTreeName+(*sysListItr).name());
    m_susyEvt->tree2->SetDirectory(outputFile);
    m_susyEvtList.push_back(m_susyEvt);
  }
  std::cout << "===============================" << std::endl;

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode evtSelection :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.

  m_hCutFlow->Fill("All", 1);
  // print every 100 events, so we know where we are:
  if((m_eventCounter % 1000)==0) Info("execute()", "Event number = %i", m_eventCounter );
  m_eventCounter++;

  // loop over recommended systematics
  int sysIdx = 0;
  std::vector<CP::SystematicSet>::const_iterator sysListItr;
  for (sysListItr = m_sysList.begin(); sysListItr != m_sysList.end(); ++sysListItr){
    if((*sysListItr).name()!="" && !m_doSys) continue;
  
    susyEvts* m_susyEvt = m_susyEvtList[sysIdx];
    sysIdx++;
    //---------------------------------------------------------------------------------------
    // config tools to apply the current sys
    CP::SystematicCode ret(CP::SystematicCode::Ok);

    if (m_objTool) {
      ret  = m_objTool->applySystematicVariation(*sysListItr);
      if ( ret != CP::SystematicCode::Ok) {
        ATH_MSG_ERROR("Cannot configure SUSYObjDefTool for systematic var. " << sysListItr->name() );
        continue;
      } else {
        ATH_MSG_VERBOSE("SUSYObjDef configured for systematic var. " << sysListItr->name() );
      }
    }
    //---------------------------------------------------------------------------------------

    //----------------------------
    // Event information
    //--------------------------- 
    const xAOD::EventInfo* eventInfo = 0;
    if( ! m_event->retrieve(eventInfo, "EventInfo").isSuccess()){
      Error("execute()", "Failed to retrieve event info collection. Exiting.");
      return EL::StatusCode::FAILURE;
     }
    m_susyEvt->evt.run = eventInfo->runNumber();
    m_susyEvt->evt.event = eventInfo->eventNumber();
    m_susyEvt->evt.lumiBlock = eventInfo->lumiBlock();
    m_susyEvt->evt.actualMu = eventInfo->actualInteractionsPerCrossing();
    m_susyEvt->evt.cuts = 0;
    m_susyEvt->truths.clear(); 

    if(m_isMC)
    { 
      m_susyEvt->evt.weight = m_objTool->GetPileupWeight();
      m_susyEvt->evt.averageMu = eventInfo->averageInteractionsPerCrossing();
    }
    else
    {
      m_susyEvt->evt.weight = 1;
      m_susyEvt->evt.averageMu = m_objTool->GetCorrectedAverageInteractionsPerCrossing();
    }

    /// grl
    if(m_isMC || m_grl->passRunLB(*eventInfo)) m_susyEvt->evt.cuts |= PASS_GRL;

    // Get the Electrons from the event
    xAOD::ElectronContainer* electrons_copy(0);
    xAOD::ShallowAuxContainer* electrons_copyaux(0);
    CHECK(m_objTool->GetElectrons(electrons_copy,electrons_copyaux));

    // Get the Muons from the event
    xAOD::MuonContainer* muons_copy(0);
    xAOD::ShallowAuxContainer* muons_copyaux(0);
    CHECK( m_objTool->GetMuons(muons_copy,muons_copyaux) );

    ///jet
    xAOD::JetContainer* jets_copy(0);
    xAOD::ShallowAuxContainer* jets_copyaux(0);
    CHECK( m_objTool->GetJets(jets_copy,jets_copyaux) );

    /// overlap removal
    if(m_objTool->GetPrimVtx() == nullptr) continue; 
    CHECK( m_objTool->OverlapRemoval(electrons_copy, muons_copy, jets_copy) );
 
    /// dumping objects
    std::vector< xAOD::IParticle* > sel_Ls;
    sel_Ls.reserve(10);

    // Electrons
    int nEL = 0;
    for(auto el: *electrons_copy){
      if(!m_objTool->IsSignalElectron(*el, ElPtCut, Eld0SigCut, Elz0Cut) || !dec_passOR(*el)) continue;
      sel_Ls.push_back(el);
      nEL++;
     }
 
    // Muons
    int nMU = 0;
    for(auto mu: *muons_copy){
      if(!m_objTool->IsSignalMuon(*mu, MuPtCut, Mud0SigCut, Muz0Cut) || !dec_passOR(*mu)) continue;
      sel_Ls.push_back(mu);
      nMU++;
    }
    m_susyEvt->sig.nEl = nEL;
    m_susyEvt->sig.nMu = nMU;
    // jets
    int nJET = 0;
    std::vector< xAOD::IParticle* > jet_Ls;
    jet_Ls.reserve(10);
    for(auto jet: *jets_copy){
      if(!m_objTool->IsSignalJet(*jet, JetPtCut, JetEtaCut) || !dec_passOR(*jet)) continue;
      nJET++;
      jet_Ls.push_back(jet);
    }



    if((nEL+nMU == nLepCutExactly || nEL+nMU >= nLepCutMin) && (nJET == nJetCutExactly || nJET >= nJetCutMin)){
      m_hCutFlow->Fill("NLepton", 1);
      /////////trigger
      if(start){
        auto chainGroup=m_objTool->GetTrigChainGroup(".*");
        Info("execute ()", "----------Available Triggers ---------" );
        for(auto &trig : chainGroup->getListOfTriggers()) {Info("execute ()", "%s",trig.c_str());}
        Info("execute ()", "----------Using Triggers ---------" );
        for(auto &trig : trigNames) {Info("execute ()", "%s",trig.c_str());};

        start=false;
      }

      //Assign trigCode
      m_susyEvt->sig.trigCode = 0;
      unsigned long int ADD = 1;
      for(unsigned int i=0; i<trigNames.size(); i++){
        if(m_objTool->IsTrigPassed(trigNames[i])) 
        {
          m_susyEvt->sig.trigCode += ADD;
          m_hTrigs->Fill(trigNames[i].c_str(), 1);
        }
        ADD *= 2;
      }
 
    //Scale factor
    if(m_isMC)
    {
      m_susyEvt->evt.ElSF = m_objTool->GetTotalElectronSF(*electrons_copy);
      m_susyEvt->evt.MuSF = m_objTool->GetTotalMuonSF(*muons_copy,true,true,"HLT_mu26_imedium");
    }
    else
    {
      m_susyEvt->evt.ElSF = 1;
      m_susyEvt->evt.MuSF = 1;
    }
//  /////////////
      /// met
      xAOD::MissingETContainer* metcst = new xAOD::MissingETContainer();
      xAOD::MissingETAuxContainer* metcst_aux = new xAOD::MissingETAuxContainer();
      metcst->setStore(metcst_aux);
      CHECK( m_objTool->GetMET(*metcst,jets_copy,electrons_copy,muons_copy) );
      auto metFinal = (*metcst)["Final"];
      m_susyEvt->sig.MetRel = metFinal->met()*iGeV;
      m_susyEvt->sig.MetX = metFinal->mpx()*iGeV;
      m_susyEvt->sig.MetY = metFinal->mpy()*iGeV;
      TLorentzVector metV(metFinal->mpx(), metFinal->mpy(), 0, metFinal->met());

      /// save leptons
      m_susyEvt->leps.resize(nEL+nMU);
      /// sort the vector
      std::sort(sel_Ls.begin(), sel_Ls.end(), [](xAOD::IParticle* a, xAOD::IParticle* b)->bool{return a->pt()>b->pt();});
      for(unsigned int i=0;i<sel_Ls.size(); i++){
        fillLepton(m_susyEvt, sel_Ls[i], m_susyEvt->leps[i], i);
        m_susyEvt->leps[i].MET_dPhi = metV.DeltaPhi(sel_Ls[i]->p4());
      }

      /// jets
      m_susyEvt->sig.nJet = jet_Ls.size();
      m_susyEvt->jets.resize(jet_Ls.size());
      int i=0;
       for(auto j0: jet_Ls){
         auto j = dynamic_cast<xAOD::Jet*>(j0);
        m_susyEvt->jets[i].pt = j->pt()*iGeV; 
        m_susyEvt->jets[i].eta = j->eta(); 
        m_susyEvt->jets[i].phi = j->phi();
        m_susyEvt->jets[i].MET_dPhi = metV.DeltaPhi(j->p4());
        unsigned int& flag = m_susyEvt->jets[i].jFlag;
        flag = 0;
        if(dec_baseline(*j)) flag |= IS_BASELINE;
        if(dec_passOR(*j)) flag |= IS_PASSOR;
        if(m_objTool->IsBadJet(*j)) flag |= IS_BAD;
        if(dec_signal(*j)) flag |= IS_SIGNAL;
        if(dec_bjet_loose(*j)) flag |= JT_BJET_LOOSE;
        if(m_objTool->IsBJet(*j)) flag |= JT_BJET;
        i++;
      }

      /// two leading leptons
      if(nEL+nMU>1){
        xAOD::IParticle* l1 = sel_Ls[0];
        xAOD::IParticle* l2 = sel_Ls[1];
        TLorentzVector ll(l1->p4()+l2->p4());
        m_susyEvt->l12.m = ll.M()*iGeV;
        m_susyEvt->l12.pt = ll.Pt()*iGeV; 
        m_susyEvt->l12.eta = ll.Eta(); 
        m_susyEvt->l12.phi = ll.Phi(); 
        m_susyEvt->l12.dPhi = l1->p4().DeltaPhi(l2->p4()); 
        m_susyEvt->l12.dR = l1->p4().DeltaR(l2->p4()); 
        m_susyEvt->l12.MET_dPhi = metV.DeltaPhi(ll);

        /// mT2
        auto tl1 = l1->p4();
        auto tl2 = l2->p4();
        //m_susyEvt->sig.mT2 =  asymm_mt2_lester_bisect::get_mT2(tl1.M(), tl1.Px(), tl1.Py(), tl2.M(), tl2.Px(), tl2.Py(), metFinal->mpx(), metFinal->mpy(), 1, 1)*iGeV;
        m_susyEvt->sig.mT2 =  anaHelper::get_mT2(tl1.M(), tl1.Px(), tl1.Py(), tl2.M(), tl2.Px(), tl2.Py(), metFinal->mpx(), metFinal->mpy(), 1, 1)*iGeV;
      }

      // Primary vertex
      m_susyEvt->sig.nPV=0;
      m_susyEvt->sig.nVtx=0;
      const xAOD::VertexContainer* vertices(0);
      if( m_event->retrieve( vertices, "PrimaryVertices" ).isSuccess() ) {
        for( const auto& vx : *vertices ) {
          m_susyEvt->sig.nVtx++;
          if(vx->vertexType() != xAOD::VxType::PriVtx) continue;
          int nTrk = 0;
          for(size_t i=0; i<vx->nTrackParticles(); i++){
            const xAOD::TrackParticle* trk = vx->trackParticle(i);
            if(trk && trk->pt()>m_vxTrkPtMin) nTrk++;
           }
          if(nTrk<m_vxTrkNMin) continue;
          m_susyEvt->sig.nPV++;
        }
      }else{
        Error("Failed to retrieve VertexContainer %s, returning NULL", "PrimaryVertices");
      }

      /// tau
      //xAOD::TauJetContainer* taus_copy(0);
      //xAOD::ShallowAuxContainer* taus_copyaux(0);
      //CHECK( m_objTool->GetTaus(taus_copy,taus_copyaux) );
      m_susyEvt->sig.nTau = 0;
      //for(auto tau: *taus_copy){
        //if(!m_objTool->IsSignalTau(*tau)) continue;
        //m_susyEvt->sig.nTau++;
      //}

      /// fill events
      m_susyEvt->fill();

      //m_store->record(taus_copy, "ShallowCopiedTaus");
      //m_store->record(taus_copyaux, "ShallowCopiedTausAux.");
      delete metcst;
      delete metcst_aux;
    }

    m_store->record(muons_copy       , sysListItr->name()+"ShallowCopiedMuons");
    m_store->record(muons_copyaux    , sysListItr->name()+"ShallowCopiedMuonsAux.");
    m_store->record(electrons_copy   , sysListItr->name()+"ShallowCopiedElectrons");
    m_store->record(electrons_copyaux, sysListItr->name()+"ShallowCopiedElectronsAux.");
    m_store->record(jets_copy        , sysListItr->name()+"ShallowCopiedJets");
    m_store->record(jets_copyaux     , sysListItr->name()+"ShallowCopiedJetsAux.");
  }

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode evtSelection :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode evtSelection :: finalize ()
{ 
  // This method is the mirror image of initialize(), meaning it gets
  // called after the last event has been processed on the worker node
  // and allows you to finish up any objects you created in
  // initialize() before they are written to disk.  This is actually
  // fairly rare, since this happens separately for each worker node.
  // Most of the time you want to do your post-processing on the
  // submission node after all your histogram outputs have been
  // merged.  This is different from histFinalize() in that it only
  // gets called on worker nodes that processed input events.
  delete m_objTool;
  //for(int i=0; i<nTrig; i++) std::cout << i+1<<" "<<trigNames[i]<<":   " <<trigCount[i] <<std::endl;
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode evtSelection :: histFinalize ()
{
  // This method is the mirror image of histInitialize(), meaning it
  // gets called after the last event has been processed on the worker
  // node and allows you to finish up any objects you created in
  // histInitialize() before they are written to disk.  This is
  // actually fairly rare, since this happens separately for each
  // worker node.  Most of the time you want to do your
  // post-processing on the submission node after all your histogram
  // outputs have been merged.  This is different from finalize() in
  // that it gets called on all worker nodes regardless of whether
  // they processed input events.
  
  //m_susyEvt->writeTree("2lEvt");

  return EL::StatusCode::SUCCESS;
}

void evtSelection :: fillLepton(susyEvts* m_susyEvt, xAOD::IParticle* p, L_PAR& l, unsigned int index)
{
    l.pt = p->pt()*iGeV;
    l.eta = p->eta();
    l.phi = p->phi();
    l.Q = 0;
    l.lFlag = 0;

    xAOD::Electron* el = 0; 
    xAOD::Muon* mu = dynamic_cast<xAOD::Muon*>(p);
    if(mu)
    {
      l.ID = 13000;
      //l.author = mu->allAuthors();
      l.author = mu->author();
      //l.ID += mu->quality();
      l.ID += m_muonSelTool->getQuality(*mu);
      l.ID *= mu->charge();
      if(!mu->isolation(l.topoetcone20, xAOD::Iso::topoetcone20)) Error("fillLepton(mu)", "topoetcone20 failed."); 
      if(!mu->isolation(l.topoetcone30, xAOD::Iso::topoetcone30)) Error("fillLepton(mu)", "topoetcone30 failed."); 
      if(!mu->isolation(l.topoetcone40, xAOD::Iso::topoetcone40)) Error("fillLepton(mu)", "topoetcone40 failed."); 
      if(!mu->isolation(l.ptcone20, xAOD::Iso::ptcone20)) Error("fillLepton(mu)", "ptcone20 failed."); 
      if(!mu->isolation(l.ptcone30, xAOD::Iso::ptcone30)) Error("fillLepton(mu)", "ptcone30 failed."); 
      if(!mu->isolation(l.ptcone40, xAOD::Iso::ptcone40)) Error("fillLepton(mu)", "ptcone40 failed."); 

      
      if(m_objTool->IsBadMuon(*mu, 0.1)) l.lFlag |= IS_BAD;
      if(m_objTool->IsCosmicMuon(*mu, 3, 0.5)) l.lFlag |= MU_COSMIC;
      const xAOD::TrackParticle* trk = mu->primaryTrackParticle();
      if(trk)
      {
        l.d0 = trk->d0();
        l.z0 = trk->z0();
        const xAOD::ParametersCovMatrix_t cov = trk->definingParametersCovMatrix();
        l.d0Err = sqrt(cov(0, 0));
        l.z0Err = sqrt(cov(1, 1));
      }
      //if(mu->passesIDCuts()) l.Q |= 1<<0;
      //if(mu->passesHighPtCuts()) l.Q |= 1<<1;

      uint8_t n;
      l.nBHits      = mu->summaryValue(n, numberOfBLayerHits) ?n:0;
      l.nPixHits    = mu->summaryValue(n, numberOfPixelHits)  ?n:0;
      l.nSCTHits    = mu->summaryValue(n, numberOfSCTHits)    ?n:0;
      l.nPixHoles   = mu->summaryValue(n, numberOfPixelHoles) ?n:0;
      l.nSCTHoles   = mu->summaryValue(n, numberOfSCTHoles)   ?n:0;
      l.nTRTHits    = mu->summaryValue(n, numberOfTRTHits)    ?n:0;
      l.nTRTOutliers= mu->summaryValue(n, numberOfTRTOutliers)?n:0;
      
      /// truth
      if(m_isMC)
       {
         //l.SF_iso = m_objTool->GetSignalMuonSF(*mu);
         //l.SF = m_objTool->GetSignalMuonSF(*mu, true, false);
         
         if(trk){
          l.truthType = acc_truthType(*trk);
          l.truthOrig = acc_truthOrig(*trk);
         }

        /// save truth match and parents if exist, otherwise save -1.
        auto tl = acc_truthLink(*mu);
        //l.truthI = tl.isValid()? addTruthPar(*tl, m_susyEvt->truths, -1):-1;
        //if(l.truthI>=0) m_susyEvt->truths[l.truthI].matchI = index;

	  
        if(tl.isValid())
         {
          l.truthI = addTruthPar(*tl, m_susyEvt->truths, -1);
          m_susyEvt->truths[l.truthI].matchI = index;
         }else l.truthI = -1;
       }
      unsigned long iso = m_isoTool->accept(*mu).getCutResultBitSet().to_ulong();
      l.lFlag |= (iso<<20); // shift 20 bits

     }else{
       el = dynamic_cast<xAOD::Electron*>(p);
       if(el)
       {
         l.ID = 11000;
         l.author = el->author();
         /// eletron quality cut not sure.
         if(m_LHToolTight2015->accept(el)) l.ID += 15;
         else if(m_LHToolMedium2015->accept(el)) l.ID += 7;
         else if(m_LHToolLoose2015->accept(el)) l.ID += 2;

         //if(el->isGoodOQ(el->OQ())) l.lFlag |= EL_OQ;

         l.ID *= el->charge();
         if(!el->isolationValue(l.topoetcone20, xAOD::Iso::topoetcone20)) Error("fillLepton(el)", "topoetcone20 failed."); 
         if(!el->isolationValue(l.topoetcone30, xAOD::Iso::topoetcone30)) Error("fillLepton(el)", "topoetcone30 failed."); 
         if(!el->isolationValue(l.topoetcone40, xAOD::Iso::topoetcone40)) Error("fillLepton(el)", "topoetcone40 failed."); 
         if(!el->isolationValue(l.ptcone20, xAOD::Iso::ptcone20)) Error("fillLepton(el)", "ptcone20 failed."); 
         if(!el->isolationValue(l.ptcone30, xAOD::Iso::ptcone30)) Error("fillLepton(el)", "ptcone30 failed."); 
         if(!el->isolationValue(l.ptcone40, xAOD::Iso::ptcone40)) Error("fillLepton(el)", "ptcone40 failed."); 
         const xAOD::TrackParticle* trk = el->trackParticle(0);
         if(trk)
         {
           l.d0 = trk->d0();
           l.z0 = trk->z0();
           const xAOD::ParametersCovMatrix_t cov = trk->definingParametersCovMatrix();
           l.d0Err = sqrt(cov(0, 0));
           l.z0Err = sqrt(cov(1, 1));
         }
         l.nBHits       =el->trackParticleSummaryIntValue(numberOfBLayerHits);
         l.nPixHits     =el->trackParticleSummaryIntValue(numberOfPixelHits);
         l.nSCTHits     =el->trackParticleSummaryIntValue(numberOfSCTHits);
         //l.nPixHoles    =el->trackParticleSummaryIntValue(numberOfPixelHoles);
         //l.nSCTHoles    =el->trackParticleSummaryIntValue(numberOfSCTHoles);
         l.nTRTHits     =el->trackParticleSummaryIntValue(numberOfTRTHits);
         l.nTRTOutliers =el->trackParticleSummaryIntValue(numberOfTRTOutliers);

         if(m_isMC)
         {
           //l.SF = m_objTool->GetSignalElecSF(*el, true, true,false, false);
           //l.SF_iso = m_objTool->GetSignalElecSF(*el,true, true, false, true);
           l.truthType = acc_truthType(*el);
           l.truthOrig = acc_truthOrig(*el);

           /// save truth match and parents if exist, otherwise save -1.
           auto tl = acc_truthLink(*el);
           //if(tl.isValid()) l.truthI = tl.isValid()? addTruthPar(*tl, m_susyEvt->truths, -1):-1;
           
           if(tl.isValid())
           {
             l.truthI = addTruthPar(*tl, m_susyEvt->truths, -1);
             m_susyEvt->truths[l.truthI].matchI = index;
           }
           else l.truthI = -1;
           
           //auto trk = el->trackParticle();
           //l.truthProb = trk?acc_truthProb(*trk):-1;
         }
         unsigned long iso = m_isoTool->accept(*el).getCutResultBitSet().to_ulong();
         l.lFlag |= (iso<<20); // shift 20 bits
       }
     }

    l.topoetcone20 *= iGeV;
    l.topoetcone30 *= iGeV;
    l.topoetcone40 *= iGeV;
    l.ptcone20 *= iGeV;
    l.ptcone30 *= iGeV;
    l.ptcone40 *= iGeV;

    ///
    unsigned int& flag = l.lFlag;
    if(dec_baseline(*p)) flag |= IS_BASELINE;
    if(dec_signal(*p)) flag |= IS_SIGNAL;
    if(dec_passOR(*p)) flag |= IS_PASSOR;

    /// check iso


    return;
}
int evtSelection::addTruthPar(const xAOD::TruthParticle* p, TRUTHS& v, int pLevel){
    /// check if already exist
    const int bcode = p->barcode();
    int nv = v.size();
    for(int i=0;i<nv;i++){if(v[i].barcode == bcode) {return i;}}

    v.emplace_back();
    TR_PAR& t = v.back();
    t.barcode = bcode;
    t.pt = p->pt()*iGeV;
    t.eta = t.pt>1?p->eta():-999;
    t.phi = p->phi();
    t.pdgId = p->pdgId(); 
    t.matchI = -1;
    t.motherI = -1;

    /// add parents if exist
    if(pLevel){
      auto m = p->parent();
      if(m){
        while(m->parent() && m->pdgId() == m->parent()->pdgId()) m=m->parent();
        int id = addTruthPar(m, v, pLevel-1);
        v[nv].motherI = id;
      }
    }
    return nv;
}
