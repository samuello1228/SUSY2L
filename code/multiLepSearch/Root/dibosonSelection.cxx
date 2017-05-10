#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <EventLoop/OutputStream.h>
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODCutFlow/CutBookkeeperContainer.h"
#include "xAODCutFlow/CutBookkeeper.h"
#include "PathResolver/PathResolver.h"
#include "GoodRunsLists/GoodRunsListSelectionTool.h"
#include "xAODMuon/Muon.h"
#include "xAODEventInfo/EventInfo.h"
#include "multiLepSearch/mCHECK.h"
#include "multiLepSearch/dibosonSelection.h"

#include <TError.h>
#include <vector>

#include <TH1D.h>
#include <TFile.h>

// this is needed to distribute the algorithm to the workers
ClassImp(dibosonSelection)

using namespace xAOD;
using namespace std;

static SG::AuxElement::Accessor<char> dec_baseline("baseline");
static SG::AuxElement::Accessor<char> dec_signal("signal");
static SG::AuxElement::Accessor<char> dec_passOR("passOR");
static SG::AuxElement::Accessor<char> dec_bad("bad");
static SG::AuxElement::Accessor<char> dec_bjet_loose("bjet_loose");
static SG::AuxElement::Accessor<char> dec_bjet("bjet");
static SG::AuxElement::Accessor<char> dec_cosmic("cosmic");
typedef ElementLink< xAOD::TruthParticleContainer > TruthLink;
static SG::AuxElement::Accessor< TruthLink > acc_truthLink("truthParticleLink"); // ID track, electron

const float iGeV = 0.001;
const float GeV = 1000;

dibosonSelection :: dibosonSelection ()
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
}



EL::StatusCode dibosonSelection :: setupJob (EL::Job& job)
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

  EL::OutputStream out(CF_outputName);
  job.outputAdd(out);

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode dibosonSelection :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.
  TFile *outputFile = wk()->getOutputFile(CF_outputName);

  m_hCutFlow = new TH1D("hCutFlow", "cut flow", 60, 0, 60);
  m_hCutFlow->SetDirectory(outputFile);

  m_zphEvt = new ZgammaEvt();
  m_zphEvt->makeTree(CF_outputTreeName);
  m_zphEvt->tree2->SetDirectory(outputFile);

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode dibosonSelection :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed

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
    if(!wk()->xaodEvent()->retrieveMetaInput(incompleteCBC, "IncompleteCutBookkeepers").isSuccess()){
      Error("fileExecute()","Failed to retrieve IncompleteCutBookkeepers from MetaData! Exiting.");
      return EL::StatusCode::FAILURE;
    }
    if( incompleteCBC->size() != 0 ) {
      Warning("fileExecute()","Found incomplete Bookkeepers! Check file for corruption.");
      //return EL::StatusCode::FAILURE;
    }
    // Now, let's find the actual information
    const xAOD::CutBookkeeperContainer* completeCBC = 0;
    if(!wk()->xaodEvent()->retrieveMetaInput(completeCBC, "CutBookkeepers").isSuccess()){
    Error("fileExecute()","Failed to retrieve CutBookkeepers from MetaData! Exiting.");
    return EL::StatusCode::FAILURE;
    }

    // First, let's find the smallest cycle number,
    // i.e., the original first processing step/cycle
    const xAOD::CutBookkeeper* allEventsCBK=0;
    const xAOD::CutBookkeeper* DxAODEventsCBK=0;
    int maxCycle = -1;
    for (const auto& cbk: *completeCBC) {
      if (cbk->cycle() > maxCycle && cbk->name() == "AllExecutedEvents" && cbk->inputStream() == "StreamAOD") {
       allEventsCBK = cbk;
       maxCycle = cbk->cycle();
      }
      if ( cbk->name() == CF_derivationName+"KernelSkim"){
        DxAODEventsCBK = cbk;
      }
    }

//     m_hCutFlow = m_hCutFlowNominal;

    if(allEventsCBK){
      m_hCutFlow->Fill("AOD", allEventsCBK->nAcceptedEvents());
      m_hCutFlow->Fill("aSumW", allEventsCBK->sumOfEventWeights());
      m_hCutFlow->Fill("aSumW2", allEventsCBK->sumOfEventWeightsSquared());
    }

    if(DxAODEventsCBK){
      m_hCutFlow->Fill(CF_derivationName.c_str(), DxAODEventsCBK->nAcceptedEvents());
      m_hCutFlow->Fill("dSumW", DxAODEventsCBK->sumOfEventWeights());
      m_hCutFlow->Fill("dSumW2", DxAODEventsCBK->sumOfEventWeightsSquared());
    }
  }
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode dibosonSelection :: changeInput (bool /*firstFile*/)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode dibosonSelection :: initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.

  /// SUSYTools
  for(auto x: CF_PRW_confFiles){Info("CF_PRW_confFiles", "%s", x.c_str());}
  m_objTool = new ST::SUSYObjDef_xAOD("SUSYObjDef_xAOD");
  ST::ISUSYObjDef_xAODTool::DataSource ds = static_cast<ST::ISUSYObjDef_xAODTool::DataSource>(CF_isMC); 
  CHECK(m_objTool->setProperty("DataSource",ds));
  CHECK(m_objTool->setProperty("ConfigFile", CF_ConfigFile));
  CHECK(m_objTool->setProperty("PRWConfigFiles", CF_PRW_confFiles));
  CHECK(m_objTool->setProperty("PRWLumiCalcFiles", CF_PRW_lcalcFiles));
  CHECK(m_objTool->initialize().isSuccess());


  /// GRL
  if(!CF_isMC){
    vector<string> vecStringGRL;
    for(auto& g: CF_grlFiles) vecStringGRL.push_back(PathResolverFindCalibFile(g));
    m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
    CHECK(m_grl->setProperty("GoodRunsListVec", vecStringGRL));
    CHECK(m_grl->setProperty("PassThrough", false)); // if true (default) will ignore result of GRL and will just pass all events
    CHECK(m_grl->initialize());
  }

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode dibosonSelection :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.

  auto sc = EL::StatusCode::SUCCESS;
  // Primary vertex
  unsigned int nPV=0;
  //m_susyEvt->sig.nVtx=0;
  const xAOD::VertexContainer* vertices(0);
  if( wk()->xaodEvent()->retrieve( vertices, "PrimaryVertices" ).isSuccess() ) {
    for(const auto& vx : *vertices) {
      if(vx->vertexType() != xAOD::VxType::PriVtx) continue;
      nPV++;
    }
  }else{
    Error("Failed to retrieve VertexContainer %s, returning NULL", "PrimaryVertices");
  }

  if(nPV<1) return sc;
  m_hCutFlow->Fill("PV", 1);


  CHECK(m_objTool->ApplyPRWTool());

  /// photon Check
  m_photons = nullptr;
  ShallowAuxContainer* photons_copyaux(0);
  CHECK(m_objTool->GetPhotons(m_photons,photons_copyaux, true));
  if(m_photons->size()==0) return EL::StatusCode::SUCCESS; /// no need to go further if no photon

  // Get the Electrons from the event
  m_electrons = nullptr;
  ShallowAuxContainer* electrons_copyaux(0);
  CHECK(m_objTool->GetElectrons(m_electrons,electrons_copyaux, true));

  /// muon Check
  m_muons = nullptr;
  ShallowAuxContainer* muons_copyaux(0);
  CHECK(m_objTool->GetMuons(m_muons,muons_copyaux, true));

  /// at lest two Same Flavour leptons
  if(m_electrons->size()<2 && m_muons->size()<2) return EL::StatusCode::SUCCESS; /// no need to go further if no photon

  ///jet
  m_jets = nullptr;
  ShallowAuxContainer* jets_copyaux(0);
  CHECK(m_objTool->GetJets(m_jets,jets_copyaux, true));

  /// overlap removal
  CHECK(m_objTool->OverlapRemoval(m_electrons, m_muons, m_jets, m_photons));

  /// get the photons
  vector< Photon* > sel_phos;
  for(auto el: *m_photons){if(dec_passOR(*el)) sel_phos.push_back(el);}
//   for(auto el: *m_photons){
//     sel_phos.push_back(el);
//    }
  if(sel_phos.size()==0) return EL::StatusCode::SUCCESS;

  /// get the two electrons
  vector< Electron* > dilep_els;
  for(auto el: *m_electrons){if(dec_passOR(*el) && dec_signal(*el)) dilep_els.push_back(el);}

  /// get the two muons
  vector< xAOD::Muon* > dilep_mus;
  for(auto l: *m_muons){if(dec_passOR(*l) && dec_signal(*l)) dilep_mus.push_back(l);}

  /// save photons
  m_zphEvt->phos.clear();
  sortByPt(sel_phos);
  for(auto ph: sel_phos){
//     Info("photons", "pt=%.2f, eta=%.2f", ph->pt()*iGeV, ph->eta());
    /// fill photons
    m_zphEvt->phos.emplace_back();
    auto& gam = m_zphEvt->phos.back();
    gam.pt = ph->pt()*iGeV;
    gam.eta = ph->eta();
    gam.phi = ph->phi();
    gam.ID = ph->conversionType();
    gam.mT = ph->conversionRadius();
   }
  auto& ph0 = sel_phos[0];

  /// save the event
  const xAOD::EventInfo* eventInfo(nullptr);
  if(! wk()->xaodEvent()->retrieve(eventInfo, "EventInfo").isSuccess()){
    Error("execute()", "Failed to retrieve event info collection. Exiting.");
    return EL::StatusCode::SUCCESS;
   }

  m_zphEvt->evt.event = eventInfo->eventNumber();
  m_zphEvt->evt.weight = CF_isMC?eventInfo->mcEventWeight():1;
  m_zphEvt->evt.isMC = CF_isMC? 1:0;
  m_hCutFlow->Fill("nSumW", m_zphEvt->evt.weight);
  m_zphEvt->truths.clear();

  /// ee gamma channel
  if(dilep_els.size()>1){
    sortByPt(dilep_els);
    if(dilep_els[0]->charge() * dilep_els[1]->charge() < 0){ /// OS, fill the events
      m_zphEvt->leps.resize(dilep_els.size());
      for(size_t i=0; i<dilep_els.size(); i++){fillLepton(dilep_els[i], m_zphEvt->leps[i], i);} 

      auto& v0 = ph0->p4();
      auto& v1 = dilep_els[0]->p4();
      auto& v2 = dilep_els[1]->p4();

      TLorentzVector ll(v1+v2);
      m_zphEvt->l12.m = ll.M()*iGeV;
      m_zphEvt->l12.pt = ll.Pt()*iGeV; 
      m_zphEvt->l12.eta = ll.Eta(); 
      m_zphEvt->l12.phi = ll.Phi(); 
      m_zphEvt->l12.dPhi = v1.DeltaPhi(v2); 
      m_zphEvt->l12.dR = v1.DeltaR(v2);

      TLorentzVector vllp(ll+v0);
      m_zphEvt->llp.m = vllp.M()*iGeV;
      m_zphEvt->llp.pt = vllp.Pt()*iGeV; 
      m_zphEvt->llp.eta = vllp.Eta(); 
      m_zphEvt->llp.phi = vllp.Phi(); 
      m_zphEvt->llp.dPhi = ll.DeltaPhi(v0); 
      m_zphEvt->llp.dR = ll.DeltaR(v0);

      m_zphEvt->evt.flag = 1;
      m_zphEvt->fill();
     }
  }

  /// mm gamma channel
  if(dilep_mus.size()>1){
    sortByPt(dilep_mus);
    if(dilep_mus[0]->charge() * dilep_mus[1]->charge() < 0){ /// OS, fill the events
      m_zphEvt->leps.resize(dilep_mus.size());
      for(size_t i=0; i<dilep_mus.size(); i++){fillLepton(dilep_mus[i], m_zphEvt->leps[i], i);} 

      auto& v0 = ph0->p4();
      auto& v1 = dilep_mus[0]->p4();
      auto& v2 = dilep_mus[1]->p4();

      TLorentzVector ll(v1+v2);
      m_zphEvt->l12.m = ll.M()*iGeV;
      m_zphEvt->l12.pt = ll.Pt()*iGeV; 
      m_zphEvt->l12.eta = ll.Eta(); 
      m_zphEvt->l12.phi = ll.Phi(); 
      m_zphEvt->l12.dPhi = v1.DeltaPhi(v2); 
      m_zphEvt->l12.dR = v1.DeltaR(v2);

      TLorentzVector vllp(ll+v0);
      m_zphEvt->llp.m = vllp.M()*iGeV;
      m_zphEvt->llp.pt = vllp.Pt()*iGeV; 
      m_zphEvt->llp.eta = vllp.Eta(); 
      m_zphEvt->llp.phi = vllp.Phi(); 
      m_zphEvt->llp.dPhi = ll.DeltaPhi(v0); 
      m_zphEvt->llp.dR = ll.DeltaR(v0);

      m_zphEvt->evt.flag = 2;
      m_zphEvt->fill();
     }
  }

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode dibosonSelection :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode dibosonSelection :: finalize ()
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
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode dibosonSelection :: histFinalize ()
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
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode dibosonSelection :: fillLepton(xAOD::Electron* el, L_PAR& l, unsigned int index)
{
  l.ID = 11000;
  l.ID *= el->charge();
  //l.author = el->author();

  //isolation
  //if(!el->isolationValue(l.topoetcone20, xAOD::Iso::topoetcone20)) Error("fillLepton(el)", "topoetcone20 failed."); 
  //if(!el->isolationValue(l.topoetcone30, xAOD::Iso::topoetcone30)) Error("fillLepton(el)", "topoetcone30 failed."); 
  //if(!el->isolationValue(l.topoetcone40, xAOD::Iso::topoetcone40)) Error("fillLepton(el)", "topoetcone40 failed."); 
  //if(!el->isolationValue(l.ptcone20, xAOD::Iso::ptcone20)) Error("fillLepton(el)", "ptcone20 failed."); 
  //if(!el->isolationValue(l.ptcone30, xAOD::Iso::ptcone30)) Error("fillLepton(el)", "ptcone30 failed."); 
  //if(!el->isolationValue(l.ptcone40, xAOD::Iso::ptcone40)) Error("fillLepton(el)", "ptcone40 failed."); 

  /*
  //track
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
  */

  if(CF_isMC)
  {
      //l.truthType = acc_truthType(*el);
      //l.truthOrig = acc_truthOrig(*el);

      /// save truth match and parents if exist, otherwise save -1.
      auto tl = acc_truthLink(*el);
      //if(tl.isValid()) l.truthI = tl.isValid()? addTruthPar(*tl, m_zphEvt->truths, -1):-1;
      
      if(tl.isValid()){
        l.truthI = addTruthPar(*tl, m_zphEvt->truths, -1);
        m_zphEvt->truths[l.truthI].matchI = index;
       }
      else l.truthI = -1;
      
      //auto trk = el->trackParticle();
      //l.truthProb = trk?acc_truthProb(*trk):-1;
  }

  fillLeptonCommon(el, l);
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode dibosonSelection :: fillLepton(xAOD::Muon* mu, L_PAR& l, unsigned int index)
{
  l.ID = 13000;
  l.ID *= mu->charge();
  //l.author = mu->author();

  //isolation
  //if(!mu->isolation(l.topoetcone20, xAOD::Iso::topoetcone20)) Error("fillLepton(mu)", "topoetcone20 failed."); 
  //if(!mu->isolation(l.topoetcone30, xAOD::Iso::topoetcone30)) Error("fillLepton(mu)", "topoetcone30 failed."); 
  //if(!mu->isolation(l.topoetcone40, xAOD::Iso::topoetcone40)) Error("fillLepton(mu)", "topoetcone40 failed."); 
  //if(!mu->isolation(l.ptcone20, xAOD::Iso::ptcone20)) Error("fillLepton(mu)", "ptcone20 failed."); 
  //if(!mu->isolation(l.ptcone30, xAOD::Iso::ptcone30)) Error("fillLepton(mu)", "ptcone30 failed."); 
  //if(!mu->isolation(l.ptcone40, xAOD::Iso::ptcone40)) Error("fillLepton(mu)", "ptcone40 failed."); 

  /*
  //track
  const xAOD::TrackParticle* trk = mu->primaryTrackParticle();
  if(trk)
  {
    l.d0 = trk->d0();
    l.z0 = trk->z0();
    const xAOD::ParametersCovMatrix_t cov = trk->definingParametersCovMatrix();
    l.d0Err = sqrt(cov(0, 0));
    l.z0Err = sqrt(cov(1, 1));
  }

  uint8_t n;
  l.nBHits      = mu->summaryValue(n, numberOfBLayerHits) ?n:0;
  l.nPixHits    = mu->summaryValue(n, numberOfPixelHits)  ?n:0;
  l.nSCTHits    = mu->summaryValue(n, numberOfSCTHits)    ?n:0;
  l.nPixHoles   = mu->summaryValue(n, numberOfPixelHoles) ?n:0;
  l.nSCTHoles   = mu->summaryValue(n, numberOfSCTHoles)   ?n:0;
  l.nTRTHits    = mu->summaryValue(n, numberOfTRTHits)    ?n:0;
  l.nTRTOutliers= mu->summaryValue(n, numberOfTRTOutliers)?n:0;
  */

  /// truth
  if(CF_isMC)
  {
    /*
    if(trk){
      l.truthType = acc_truthType(*trk);
      l.truthOrig = acc_truthOrig(*trk);
    }
    */

    /// save truth match and parents if exist, otherwise save -1.
    auto tl = acc_truthLink(*mu);
    //l.truthI = tl.isValid()? addTruthPar(*tl, m_zphEvt->truths, -1):-1;
    //if(l.truthI>=0) m_zphEvt->truths[l.truthI].matchI = index;

    if(tl.isValid())
    {
      l.truthI = addTruthPar(*tl, m_zphEvt->truths, -1);
      m_zphEvt->truths[l.truthI].matchI = index;
      }else l.truthI = -1;
  }
  fillLeptonCommon(mu, l);
  return EL::StatusCode::SUCCESS;
}

void dibosonSelection :: fillLeptonCommon(xAOD::IParticle* p, L_PAR& l)
{
  //l.topoetcone20 *= iGeV;
  //l.topoetcone30 *= iGeV;
  //l.topoetcone40 *= iGeV;
  //l.ptcone20 *= iGeV;
  //l.ptcone30 *= iGeV;
  //l.ptcone40 *= iGeV;

  // 
  l.pt = p->pt()*iGeV;
  l.eta = p->eta();
  l.phi = p->phi(); 

  //flag 
  l.lFlag = 0;
  unsigned int& flag = l.lFlag;
  if(dec_baseline(*p)) flag |= IS_BASELINE;
  if(dec_signal(*p)) flag |= IS_SIGNAL;
  if(dec_passOR(*p)) flag |= IS_PASSOR;

  //l.z0sinTheta = dec_z0sinTheta(*p);
  //l.d0sig = dec_d0sig(*p);
}

EL::StatusCode dibosonSelection :: fillLepton(xAOD::IParticle* p, L_PAR& l, unsigned int index)
{
  xAOD::Muon* mu = dynamic_cast<xAOD::Muon*>(p);
  if(mu) fillLepton(mu, l, index);
  else{
    xAOD::Electron* el = dynamic_cast<xAOD::Electron*>(p);
    if(el) fillLepton(el, l, index);
  }
  return EL::StatusCode::SUCCESS;
}
int dibosonSelection::addTruthPar(const xAOD::TruthParticle* p, TRUTHS& v, int pLevel){
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

  /*
  std::pair<MCTruthPartClassifier::ParticleType, MCTruthPartClassifier::ParticleOrigin> res;
  res = m_truthClassifier->particleTruthClassifier(p);
  t.particleType = res.first;
  t.particleOrigin = res.second;
  */

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


