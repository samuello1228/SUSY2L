#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <multiLepSearch/dilepSelection.h>
// #include <multiLepSearch/MT2.h>
#include <multiLepSearch/anaHelper.h>
// #include "xAODTruth/TruthParticleContainer.h"
// #include "xAODTruth/TruthParticle.h"
// #include "xAODMuon/Muon.h"
#include "GoodRunsLists/GoodRunsListSelectionTool.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODCutFlow/CutBookkeeperContainer.h"
#include "xAODCutFlow/CutBookkeeper.h"
#include "PathResolver/PathResolver.h"
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
static SG::AuxElement::Decorator<float> dec_d0sig("d0sig");
typedef ElementLink< xAOD::TruthParticleContainer > TruthLink;
static SG::AuxElement::Accessor< int > acc_truthType("truthType");
static SG::AuxElement::Accessor< int > acc_truthOrig("truthOrigin");
// static SG::AuxElement::Accessor< float > acc_truthProb("truthMatchProbability"); // only ID track
static SG::AuxElement::Accessor< TruthLink > acc_truthLink("truthParticleLink"); // ID track, electron
// this is needed to distribute the algorithm to the workers
ClassImp(dilepSelection)


const float iGeV = 0.001;
const float GeV = 1000;

dilepSelection :: dilepSelection (std::string name)
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
  m_name = name;

  m_vxTrkNMin  = 2;
  m_vxTrkPtMin = 2*GeV;

  outputName = "test.root";
  m_grl = 0;

  //nLepCutExactly = 2;
  //nLepCutMin = 2;
  //nJetCutExactly = 2;
  //nJetCutMin = 2;

  m_dataType = 0;
  m_doSys = 0;
  m_grlFile = "multiLepSearch/data15_13TeV.periodAllYear_DetStatus-v63-pro18-01_DQDefects-00-01-02_PHYS_StandardGRL_All_Good.xml";
  m_susyToolCfgFile = "multiLepSearch/sel_conf/SS_selection.conf";

  electrons_copy    = NULL;
  electrons_copyaux = NULL;
  muons_copy        = NULL;
  muons_copyaux     = NULL;
  jets_copy         = NULL;
  jets_copyaux      = NULL;
  taus_copy         = NULL;
  taus_copyaux      = NULL;
  metcst            = NULL;
  metcst_aux        = NULL;
}

EL::StatusCode dilepSelection :: setupJob (EL::Job& job)
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



EL::StatusCode dilepSelection :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.
  TFile *outputFile = wk()->getOutputFile(outputName);

  m_hCutFlow = new TH1F("hCutFlow", "cut flow", 20, 0, 20);
  m_hCutFlow->SetDirectory(outputFile);

  m_hTrigs = new TH1F("hTrigs", "n pass trigger", trigNames.size(), 0, trigNames.size());
  for(unsigned int i=0; i<trigNames.size(); i++){
    m_hTrigs->GetXaxis()->SetBinLabel(i+1,trigNames[i].c_str());
  }
  m_hTrigs->SetDirectory(outputFile);
  
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode dilepSelection :: fileExecute ()
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



EL::StatusCode dilepSelection :: changeInput (bool /*firstFile*/)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode dilepSelection :: initialize ()
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

  Info("initialize()", "Number of events = %lli", m_event->getEntries() ); // print long long int

  // count number of events
  m_eventCounter = 0;

  // SUSYTool
  m_objTool = new ST::SUSYObjDef_xAOD("SUSYObjDef_xAOD");

  CHECK(m_objTool->setProperty("DataSource",static_cast<ST::SUSYObjDef_xAOD::DataSource>(m_dataType) ));
  CHECK(m_objTool->setProperty("ConfigFile", m_susyToolCfgFile));
  CHECK(m_objTool->setProperty("PRWConfigFiles", PRW_confFiles));
  CHECK(m_objTool->setProperty("PRWLumiCalcFiles", PRW_lcalcFiles));

  //m_objTool->msg().setLevel( MSG::VERBOSE );
  m_objTool->msg().setLevel( MSG::ERROR );
  //m_objTool->msg().setLevel( MSG::INFO );

  CHECK(m_objTool->initialize().isSuccess());

  // GRLTool
  if(m_objTool->isData()){
    m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
    std::vector<std::string> vecStringGRL;
    vecStringGRL.push_back(PathResolverFindCalibFile(m_grlFile));
    CHECK(m_grl->setProperty( "GoodRunsListVec", vecStringGRL));
    CHECK(m_grl->setProperty("PassThrough", false)); // if true (default) will ignore result of GRL and will just pass all events
    CHECK(m_grl->initialize());
  }

  // ChargeFlipBkgTool
  mChargeFlipBkgTool = new ChargeFlipBkgTool("ChargeFlipBkgTool");
  CHECK(mChargeFlipBkgTool->initialize());

  // FakeLepBkgTool
  mFakeLepBkgTool    = new FakeLepBkgTool("FakeLepBkgTool");
  CHECK(mFakeLepBkgTool->initialize());

  // Systematics
  const CP::SystematicRegistry& registry = CP::SystematicRegistry::getInstance();
  const CP::SystematicSet& recommendedSystematics = registry.recommendedSystematics(); // get list of recommended systematics
  m_sysList = CP::make_systematics_vector(recommendedSystematics); 

  // Init output trees and conenct it to outFile
  TFile *outputFile = wk()->getOutputFile(outputName);
  m_dilepEvt = new dilepEvts("dilep", outputFile);

  // print list of recommended systematics
  std::vector<ST::SystInfo> systInfoList = m_objTool->getSystInfoList();
  std::cout << "======  Systematic list  ======" << std::endl;
  for (auto& aSystInfo : systInfoList){
    if(aSystInfo.systset.name()=="") std::cout << "Nominal (no syst) "  << std::endl;
    else if (m_doSys) std::cout << "Systematic: " << aSystInfo.systset.name() << std::endl;
  }
  std::cout << "===============================" << std::endl;

  firstLoop = true; //indicates the 1st loop

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode dilepSelection :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.

  // print every 100 events, so we know where we are:
  if((m_eventCounter % 1000)==0) Info("execute()", "Event number = %i", m_eventCounter );
  m_eventCounter++;

  // loop over recommended systematics
  std::vector<ST::SystInfo> systInfoList = m_objTool->getSystInfoList();
  for (const auto& aSystInfo : systInfoList){
    //cases:
    //if the sys affect kinematics, rerun everything and dedicate a tree to it in output
    //if the sys affect weight only, give a weight field for it in nominal tree
    //if the sys affect both, rerun everything and dedicate a tree to it in output
    //if the sys affect none, skip (except if its the nominal "sys")
    
    bool isNominal = (aSystInfo.systset.name()=="") ;
    if(!isNominal){
      if(!m_doSys) continue;  
      if(m_objTool->isData()) continue; // sys in systInfoList only applies on MC
      if(!aSystInfo.affectsKinematics) continue;
    }

    //---------------------------------------------------------------------------------------
    // config tools to apply the current sys
    //---------------------------------------------------------------------------------------
    CP::SystematicCode ret(CP::SystematicCode::Ok);

    if (m_objTool) {
      ret  = m_objTool->applySystematicVariation(aSystInfo.systset);
      if ( ret != CP::SystematicCode::Ok) {
        ATH_MSG_ERROR("Cannot configure SUSYObjDefTool for systematic var. " << aSystInfo.systset.name() );
        continue;
      } else {
        ATH_MSG_VERBOSE("SUSYObjDef configured for systematic var. " << aSystInfo.systset.name() );
      }
    }

    if(isNominal) m_hCutFlow->Fill("All", 1);// only fill cutflow for nomianl sys

    //----------------------------
    // Cut on Primary vertex
    //----------------------------
    m_dilepEvt->sig.nPV=0;
    m_dilepEvt->sig.nVtx=0;
    const xAOD::VertexContainer* vertices(0);
    if( m_event->retrieve( vertices, "PrimaryVertices" ).isSuccess() ) {
      for( const auto& vx : *vertices ) {
        m_dilepEvt->sig.nVtx++;
        if(vx->vertexType() != xAOD::VxType::PriVtx) continue;
        int nTrk = 0;
        for(size_t i=0; i<vx->nTrackParticles(); i++){
          const xAOD::TrackParticle* trk = vx->trackParticle(i);
          if(trk && trk->pt()>m_vxTrkPtMin) nTrk++;
         }
        if(nTrk<m_vxTrkNMin) continue;
        m_dilepEvt->sig.nPV++;
      }
    }else{
      Error("Failed to retrieve VertexContainer %s, returning NULL", "PrimaryVertices");
    }
    if(m_dilepEvt->sig.nPV<1) return EL::StatusCode::SUCCESS;
    if(isNominal) m_hCutFlow->Fill("nPV", 1);// only fill for nomianl sys

    //----------------------------
    // Cut on GRL, Event information
    //--------------------------- 
    const xAOD::EventInfo* eventInfo = 0;
    if( ! m_event->retrieve(eventInfo, "EventInfo").isSuccess()){
      Error("execute()", "Failed to retrieve event info collection. Exiting.");
      return EL::StatusCode::FAILURE;
     }
    m_dilepEvt->evtInfo.run = eventInfo->runNumber();
    m_dilepEvt->evtInfo.event = eventInfo->eventNumber();
    m_dilepEvt->evtInfo.lumiBlock = eventInfo->lumiBlock();
    m_dilepEvt->evtInfo.actualMu = eventInfo->actualInteractionsPerCrossing();

    if(!m_objTool->isData()){ 
      m_dilepEvt->evtInfo.weight = m_objTool->GetPileupWeight();
      m_dilepEvt->evtInfo.averageMu = eventInfo->averageInteractionsPerCrossing();
    }else{
      m_dilepEvt->evtInfo.weight = 1.0;
      m_dilepEvt->evtInfo.averageMu = m_objTool->GetCorrectedAverageInteractionsPerCrossing();
    }

    /// grl
    if(m_objTool->isData() && !m_grl->passRunLB(*eventInfo)) continue;
    m_dilepEvt->evtInfo.passGRL = true;

    if(isNominal) m_hCutFlow->Fill("GRL", 1);// only fill for nomianl sys
    
    //---------------------------------------------------------------------------------------
    //Obtain "raw" objects
    //---------------------------------------------------------------------------------------
    //Quote dongliang:
    //"The store takes care of the deletion of the shadow copies. 
    //So all the shadow copies should be added to store"

    // Get the Electrons from the event
    CHECK(m_objTool->GetElectrons(electrons_copy,electrons_copyaux));
    m_store->record(electrons_copy   , aSystInfo.systset.name()+"ShallowCopiedElectrons");
    m_store->record(electrons_copyaux, aSystInfo.systset.name()+"ShallowCopiedElectronsAux.");

    // Get the Muons from the event
    CHECK( m_objTool->GetMuons(muons_copy,muons_copyaux) );
    m_store->record(muons_copy       , aSystInfo.systset.name()+"ShallowCopiedMuons");
    m_store->record(muons_copyaux    , aSystInfo.systset.name()+"ShallowCopiedMuonsAux.");

    ///jet
    CHECK( m_objTool->GetJets(jets_copy,jets_copyaux) );
    m_store->record(jets_copy        , aSystInfo.systset.name()+"ShallowCopiedJets");
    m_store->record(jets_copyaux     , aSystInfo.systset.name()+"ShallowCopiedJetsAux.");

    ///tau
    //FIXME: Skip tau for MC for now, as GetTau crash with
    //TauAnalysisTools::CommonSmearingTool::checkTruthMatch(const TauJet&) const): 
    //No truth match information available. Please run TauTruthMatchingTool first.
    //if (m_objTool->isData()){
      CHECK( m_objTool->GetTaus(taus_copy,taus_copyaux) );
      m_store->record(taus_copy        , aSystInfo.systset.name()+"ShallowCopiedTaus");
      m_store->record(taus_copyaux     , aSystInfo.systset.name()+"ShallowCopiedTausAux.");
    //}

    //---------------------------------------------------------------------------------------

    //----------------------
    /// overlap removal
    //----------------------
    CHECK( m_objTool->OverlapRemoval(electrons_copy, muons_copy, jets_copy) );
 
    //----------------------
    /// pick out baseline objects
    //----------------------
    std::vector< xAOD::IParticle* > baseLeps;
    baseLeps.reserve(10);

    // Tau (only for veto)
    int nTau = 0;
    //FIXME: Skip tau for MC for now, as GetTau crash with
    //TauAnalysisTools::CommonSmearingTool::checkTruthMatch(const TauJet&) const): 
    //No truth match information available. Please run TauTruthMatchingTool first.
    //if (m_objTool->isData()){
      for (auto tau: *taus_copy) {
        if (dec_signal(*tau)) {nTau+=1;}
      }
    //}
    //if (hasTau) continue;    //veto events with tau
    //if((*sysListItr).name()=="") m_hCutFlow->Fill("noTau", 1);// only fill for nomianl sys

    // Electrons
    int nEL = 0;
    for(auto el: *electrons_copy){
      //if(!m_objTool->IsSignalElectron(*el, ElPtCut, Eld0SigCut, Elz0Cut) || !dec_passOR(*el)) continue;
      if(!dec_passOR(*el)) continue;
      baseLeps.push_back(el);
      nEL++;
    }
 
    // Muons
    int nMU = 0;
    int nCosmic = 0;
    for(auto mu: *muons_copy){
      //if(!m_objTool->IsSignalMuon(*mu, MuPtCut, Mud0SigCut, Muz0Cut) || !dec_passOR(*mu)) continue;
      if(!dec_passOR(*mu)) continue;
      if( dec_cosmic(*mu)) nCosmic += 1;
      baseLeps.push_back(mu);
      nMU++;
    }

    // Jets
    int nJET = 0;
    int nBJet = 0;
    std::vector< xAOD::IParticle* > baseJets;
    baseJets.reserve(10);
    for(auto jet: *jets_copy){
      if(!dec_passOR(*jet)) continue;
      if(m_objTool->IsBJet(*jet)) nBJet += 1; //FIXME:: not sure if it should go before or after overlap removal
      nJET++;
      baseJets.push_back(jet);
    }

    //std::cout << "N "  <<  nEL << " " <<  nMU << " " <<  nJET << std::endl;    



    //----------------------
    //Cut For 2LSS: veto evt with cosmics
    //----------------------
    //if (hasCosmic) continue;
    //if((*sysListItr).name()=="") m_hCutFlow->Fill("noCosmicMu", 1);// only fill for nomianl sys

    //----------------------
    //Cut For 2LSS: veto evt with B jets
    //----------------------
    //if(nBJet) continue;
    //if((*sysListItr).name()=="") m_hCutFlow->Fill("noBJets", 1);// only fill for nomianl sys

    //----------------------
    //Cut For 2LSS: exactly 2Lep, nJet>=0
    //----------------------
    //if((nEL+nMU == nLepCutExactly || nEL+nMU >= nLepCutMin) && (nJET == nJetCutExactly || nJET >= nJetCutMin))
    if(nEL+nMU != 2) continue;
    if(isNominal) m_hCutFlow->Fill("NLepton", 1);// only fill for nomianl sys



    //--------------------------------------------------------------------------------------------
    //Fill up data valid for all output stream
    //--------------------------------------------------------------------------------------------
    
    //Assign trigCode
    if(firstLoop){
      auto chainGroup=m_objTool->GetTrigChainGroup(".*");
      Info("execute ()", "----------Available Triggers ---------" );
      for(auto &trig : chainGroup->getListOfTriggers()) {Info("execute ()", "%s",trig.c_str());}
      Info("execute ()", "----------Using Triggers ---------" );
      for(auto &trig : trigNames) {Info("execute ()", "%s",trig.c_str());};
      if (trigNames.size()>64) Error("execute()", "Cannot treat >64 trigger types.");
    }

    m_dilepEvt->sig.trigCode = 0;
    uint64_t ADD = 1; //at most hold 64 type of triggers
    for(unsigned int i=0; i<trigNames.size(); i++){
      if(m_objTool->IsTrigPassed(trigNames[i])){
        m_dilepEvt->sig.trigCode += ADD;
        m_hTrigs->Fill(trigNames[i].c_str(), 1);
      }
      ADD *= 2;
    }

    //if (m_dilepEvt->sig.trigCode == 0) {
    //  Info("execute ()", "Interesting Trigger pass cut but not used:" );
    //  auto chainGroup=m_objTool->GetTrigChainGroup(".*");
    //  for(auto &trig : chainGroup->getListOfTriggers()) {
    //    if(m_objTool->IsTrigPassed(trig.c_str())) Info("execute ()", "%s",trig.c_str());
    //  }
    //  Info("execute ()", "---End Trig---" );
    //} 

    // met
    metcst     = new xAOD::MissingETContainer();
    metcst_aux = new xAOD::MissingETAuxContainer();
    metcst->setStore(metcst_aux);
    m_store->record(metcst    , aSystInfo.systset.name()+"ShallowCopiedMetcst");
    m_store->record(metcst_aux, aSystInfo.systset.name()+"ShallowCopiedMetcstAux.");

    CHECK( m_objTool->GetMET(*metcst,jets_copy,electrons_copy,muons_copy) );
    auto metFinal = (*metcst)["Final"];
    m_dilepEvt->sig.Met  = metFinal->met()*iGeV;
    m_dilepEvt->sig.MetX = metFinal->mpx()*iGeV;
    m_dilepEvt->sig.MetY = metFinal->mpy()*iGeV;
    TLorentzVector metV(metFinal->mpx(), metFinal->mpy(), 0, metFinal->met());
    float minMetdPhi = FLT_MAX; //Met's dPhi from nearest obj (e/mu/jet)
    
    // save leptons
    m_dilepEvt->sig.nTau = nTau;
    m_dilepEvt->sig.nCosmic = nCosmic;
    m_dilepEvt->sig.nEl = nEL;
    m_dilepEvt->sig.nMu = nMU;
    m_dilepEvt->leps.resize(nEL+nMU);
    std::sort(baseLeps.begin(), baseLeps.end(), [](xAOD::IParticle* a, xAOD::IParticle* b)->bool{return a->pt()>b->pt();});
    for(unsigned int i=0;i<baseLeps.size(); i++){
      fillLepton(baseLeps[i], m_dilepEvt->leps[i]);
      m_dilepEvt->leps[i].mT = sqrt(2.*(metV.Pt()*baseLeps[i]->pt()-metV.Px()*baseLeps[i]->p4().Px()-metV.Py()*baseLeps[i]->p4().Py()))*iGeV;
      m_dilepEvt->leps[i].MET_dPhi = metV.DeltaPhi(baseLeps[i]->p4());
      minMetdPhi = std::min(minMetdPhi, std::abs(m_dilepEvt->leps[i].MET_dPhi));
    }

    // save jets
    m_dilepEvt->sig.nBJet = nBJet;
    m_dilepEvt->sig.nJet = baseJets.size();
    m_dilepEvt->jets.resize(baseJets.size());
    std::sort(baseJets.begin(), baseJets.end(), [](xAOD::IParticle* a, xAOD::IParticle* b)->bool{return a->pt()>b->pt();});
    int i=0;
    for(auto j0: baseJets){
      auto j = dynamic_cast<xAOD::Jet*>(j0);
      m_dilepEvt->jets[i].pt  = j->pt()*iGeV; 
      m_dilepEvt->jets[i].eta = j->eta(); 
      m_dilepEvt->jets[i].phi = j->phi();
      m_dilepEvt->jets[i].MET_dPhi = metV.DeltaPhi(j->p4());
      minMetdPhi = std::min(minMetdPhi, std::abs(m_dilepEvt->jets[i].MET_dPhi));
      m_dilepEvt->jets[i].passOR     = dec_passOR(    *j)? true : false;
      m_dilepEvt->jets[i].isBaseline = dec_baseline(  *j)? true : false;
      m_dilepEvt->jets[i].isBad      = dec_bad(       *j)? true : false;
      m_dilepEvt->jets[i].isSig      = dec_signal(    *j)? true : false;
      m_dilepEvt->jets[i].isBJetLoose= dec_bjet_loose(*j)? true : false;
      m_dilepEvt->jets[i].isBJet     = dec_bjet(      *j)? true : false;
      i++;
    }

    m_dilepEvt->sig.MetRel = m_dilepEvt->sig.Met;
    if (minMetdPhi<1.570796327) m_dilepEvt->sig.MetRel *= sin(minMetdPhi);

    // two leading leptons
    if(nEL+nMU>1){
      xAOD::IParticle* l1 = baseLeps[0];
      xAOD::IParticle* l2 = baseLeps[1];
      TLorentzVector ll(l1->p4()+l2->p4());
      m_dilepEvt->l12.m    = ll.M()*iGeV;
      m_dilepEvt->l12.pt   = ll.Pt()*iGeV; 
      m_dilepEvt->l12.eta  = ll.Eta(); 
      m_dilepEvt->l12.phi  = ll.Phi(); 
      m_dilepEvt->l12.dPhi = l1->p4().DeltaPhi(l2->p4()); 
      m_dilepEvt->l12.dR   = l1->p4().DeltaR(l2->p4()); 
      m_dilepEvt->l12.MET_dPhi = metV.DeltaPhi(ll);

      /// mT2
      auto tl1 = l1->p4();
      auto tl2 = l2->p4();
      m_dilepEvt->sig.mT2 = anaHelper::get_mT2(tl1.M(), tl1.Px(), tl1.Py(), tl2.M(), tl2.Px(), tl2.Py(), metFinal->mpx(), metFinal->mpy(), 1, 1)*iGeV;
    }

    // Scale factor
    m_dilepEvt->sig.ElSF = m_objTool->GetTotalElectronSF(*electrons_copy);
    m_dilepEvt->sig.MuSF = m_objTool->GetTotalMuonSF(*muons_copy,true,true,"HLT_mu26_imedium");
    //FIXME: need to consider tau SF or not if we just veto events with it?
    //m_dilepEvt->sig.TauSF = ??? 
    //FIXME: need to consider Btag SF or not if we just veto events with it?
    //       also, need to get only passOR jets into BTagSF calculation as I only run m_objtool->IsBJet on those
    //m_dilepEvt->sig.BTagSF = m_objTool->BTagSF(???)

    // Weight based systematics for MC, only do on nominal tree
    if(m_doSys && isNominal && !m_objTool->isData()){
      for (const auto& bSystInfo : systInfoList){
        if(bSystInfo.affectsWeights && !bSystInfo.affectsKinematics){

          // config tools to apply the sys
          if (m_objTool->applySystematicVariation(bSystInfo.systset) != CP::SystematicCode::Ok) {
            ATH_MSG_ERROR("Cannot configure SUSYObjDefTool for systematic var. " << bSystInfo.systset.name() ); continue;
          }else { 
            ATH_MSG_VERBOSE("SUSYObjDef configured for systematic var. " << bSystInfo.systset.name() );
	  }

          float sysWeight = 1.0;
          bool  hasEffect = true;
	  switch (bSystInfo.affectsType) {
            //case ST::SystObjType::Unknown     : break;
	    //case ST::SystObjType::Jet         : break;
	    case ST::SystObjType::Egamma      : sysWeight = m_objTool->GetTotalElectronSF(*electrons_copy); break;
	    case ST::SystObjType::Electron    : sysWeight = m_objTool->GetTotalElectronSF(*electrons_copy); break;
	    //case ST::SystObjType::Photon      : break;
	    case ST::SystObjType::Muon        : sysWeight = m_objTool->GetTotalMuonSF(*muons_copy,true,true,"HLT_mu26_imedium"); break;
	    //case ST::SystObjType::Tau         : break;
	    //case ST::SystObjType::BTag        : break;
	    //case ST::SystObjType::MET_TST     : break;
	    //case ST::SystObjType::MET_CST     : break;
	    case ST::SystObjType::EventWeight : sysWeight = m_objTool->isData()? 1.0 : m_objTool->GetPileupWeight(); break;
	    default                           : hasEffect = false;
	  }

	  if (hasEffect)  m_dilepEvt->setWeight( bSystInfo.systset.name(), sysWeight);
	}
      }
    }

    // Weight based systematics for data, only do on nominal tree
    if(m_doSys && isNominal && m_objTool->isData()){
      m_dilepEvt->setWeight( "CFLIP_0__1up"    , 1.0);
      m_dilepEvt->setWeight( "CFLIP_0__1down"  , 1.0);
      m_dilepEvt->setWeight( "FAKELEP_0__1up"   , 1.0);
      m_dilepEvt->setWeight( "FAKELEP_0__1down" , 1.0);
    }

    //--------------------------------------------------------------------------------------------
    //Decide stream and fill up data valid for that stream
    //--------------------------------------------------------------------------------------------
    
    bool is2SigLeps = dec_signal(*baseLeps[0]) && dec_signal(*baseLeps[1]);
    bool isSS       = ((m_dilepEvt->leps[0].ID>0)==(m_dilepEvt->leps[1].ID>0)); //because we filled that info in fillLepton

    m_dilepEvt->sig.CFlipWeight0   = 1.0;
    m_dilepEvt->sig.FakeLepWeight0 = 1.0;

    //2L OS, charge flip bkg, assign flip prob weight and save as charge flip tree
    if (is2SigLeps && !isSS && m_objTool->isData() && m_dilepEvt->sig.nEl>0){
      m_dilepEvt->sig.CFlipWeight0 = mChargeFlipBkgTool->GetWeight( baseLeps, 0, 0);
      if (m_doSys){
        m_dilepEvt->setWeight( "CFLIP_0__1up"  , mChargeFlipBkgTool->GetWeight( baseLeps, 1, 0) );
        m_dilepEvt->setWeight( "CFLIP_0__1down", mChargeFlipBkgTool->GetWeight( baseLeps,-1, 0) );
      }
      m_dilepEvt->fill("CFlip_"+ aSystInfo.systset.name(), (isNominal && m_doSys) );
    }

    //2LSS, both signal leptons -> target signal
    if (is2SigLeps && isSS){
      m_dilepEvt->fill("Data_"+ aSystInfo.systset.name(), (isNominal && m_doSys) );
    } 
    
    //2LSS, no check on leptons loose/tight to get fake bkg weight
    if (isSS && m_objTool->isData()){
      m_dilepEvt->sig.FakeLepWeight0 = mFakeLepBkgTool->GetWeight( baseLeps, 0, 0); //sigma, n-th sys
      if (m_doSys){
        m_dilepEvt->setWeight( "FAKELEP_0__1up"  , mFakeLepBkgTool->GetWeight( baseLeps, 1, 0) );
        m_dilepEvt->setWeight( "FAKELEP_0__1down", mFakeLepBkgTool->GetWeight( baseLeps,-1, 0) );
      }
      m_dilepEvt->fill("FakeLep_"+ aSystInfo.systset.name(), (isNominal && m_doSys) );
    }

    firstLoop = false;
    
  }

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode dilepSelection :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode dilepSelection :: finalize ()
{
  //const char* APP_NAME = m_name.c_str();
  
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



EL::StatusCode dilepSelection :: histFinalize ()
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
  
  //m_dilepEvt->writeTree("2lEvt");

  return EL::StatusCode::SUCCESS;
}

void dilepSelection :: fillLepton(xAOD::IParticle* p, ntupLep& l)
{
    l.pt  = p->pt()*iGeV;
    l.eta = p->eta();
    l.phi = p->phi();

    xAOD::Muon* mu = dynamic_cast<xAOD::Muon*>(p);
    if(mu){
      l.ID = 13; l.ID *= mu->charge();
      //l.author = mu->allAuthors();
      l.author = mu->author();

      if(!mu->isolation(l.topoetcone20, xAOD::Iso::topoetcone20)) Error("fillLepton(mu)", "topoetcone20 failed."); 
      if(!mu->isolation(l.topoetcone30, xAOD::Iso::topoetcone30)) Error("fillLepton(mu)", "topoetcone30 failed."); 
      if(!mu->isolation(l.topoetcone40, xAOD::Iso::topoetcone40)) Error("fillLepton(mu)", "topoetcone40 failed."); 
      if(!mu->isolation(l.ptcone20    , xAOD::Iso::ptcone20    )) Error("fillLepton(mu)", "ptcone20 failed."); 
      if(!mu->isolation(l.ptcone30    , xAOD::Iso::ptcone30    )) Error("fillLepton(mu)", "ptcone30 failed."); 
      if(!mu->isolation(l.ptcone40    , xAOD::Iso::ptcone40    )) Error("fillLepton(mu)", "ptcone40 failed."); 

      l.isBad    = dec_bad(   *mu)? true : false;
      l.isCosmic = dec_cosmic(*mu)? true : false;

      const xAOD::TrackParticle* trk = mu->primaryTrackParticle();
      if(trk){
        l.d0 = trk->d0();
        l.z0 = trk->z0();
        const xAOD::ParametersCovMatrix_t cov = trk->definingParametersCovMatrix();
        l.d0Err = sqrt(cov(0, 0));
        l.z0Err = sqrt(cov(1, 1));
      }
      
    }

    xAOD::Electron* el = dynamic_cast<xAOD::Electron*>(p);
    if(el){
      l.ID = 11; l.ID *= el->charge();
      l.author = el->author();

      //if(el->isGoodOQ(el->OQ())) l.lFlag |= EL_OQ;

      if(!el->isolationValue(l.topoetcone20, xAOD::Iso::topoetcone20)) Error("fillLepton(el)", "topoetcone20 failed."); 
      if(!el->isolationValue(l.topoetcone30, xAOD::Iso::topoetcone30)) Error("fillLepton(el)", "topoetcone30 failed."); 
      if(!el->isolationValue(l.topoetcone40, xAOD::Iso::topoetcone40)) Error("fillLepton(el)", "topoetcone40 failed."); 
      if(!el->isolationValue(l.ptcone20    , xAOD::Iso::ptcone20    )) Error("fillLepton(el)", "ptcone20 failed."); 
      if(!el->isolationValue(l.ptcone30    , xAOD::Iso::ptcone30    )) Error("fillLepton(el)", "ptcone30 failed."); 
      if(!el->isolationValue(l.ptcone40    , xAOD::Iso::ptcone40    )) Error("fillLepton(el)", "ptcone40 failed."); 

      const xAOD::TrackParticle* trk = el->trackParticle(0);
      if(trk){
        l.d0 = trk->d0();
        l.z0 = trk->z0();
        const xAOD::ParametersCovMatrix_t cov = trk->definingParametersCovMatrix();
        l.d0Err = sqrt(cov(0, 0));
        l.z0Err = sqrt(cov(1, 1));
      }
    }

    l.topoetcone20 *= iGeV;
    l.topoetcone30 *= iGeV;
    l.topoetcone40 *= iGeV;
    l.ptcone20 *= iGeV;
    l.ptcone30 *= iGeV;
    l.ptcone40 *= iGeV;

    return;
}

