#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <multiLepSearch/aDAODMaker.h>
#include "EventLoop/OutputStream.h"
//Infrastructure include(s):
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include <multiLepSearch/mCHECK.h>
//
// this is needed to distribute the algorithm to the workers
ClassImp(aDAODMaker)



aDAODMaker :: aDAODMaker ():m_name("aDAODMaker")
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
}



EL::StatusCode aDAODMaker :: setupJob (EL::Job& job)
{
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.

  job.useXAOD ();
  CHECK(xAOD::Init()); // call before opening first file

  EL::OutputStream out("outputLabel", "xAOD");
  job.outputAdd(out);

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode aDAODMaker :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode aDAODMaker :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode aDAODMaker :: changeInput (bool /*firstFile*/)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode aDAODMaker :: initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.
  //

  xAOD::TEvent* event = wk()->xaodEvent();

  TFile *file = wk()->getOutputFile("outputLabel");
  CHECK(event->writeTo(file));

//   event->setAuxItemList( "MuonsAux.", "truthType" );

  /// SUSYTools
  m_objTool = new ST::SUSYObjDef_xAOD("SUSYObjDef_xAOD");
  //m_objTool->msg().setLevel( MSG::ERROR );
  //m_objTool->msg().setLevel( MSG::WARNING );

  ST::ISUSYObjDef_xAODTool::DataSource ds = static_cast<ST::ISUSYObjDef_xAODTool::DataSource>(CF_isMC); 
  CHECK(m_objTool->setProperty("DataSource",ds));
  CHECK(m_objTool->setProperty("ConfigFile", CF_ConfigFile));
//   CHECK(m_objTool->setProperty("PRWConfigFiles", CF_PRW_confFiles));
//   CHECK(m_objTool->setProperty("PRWLumiCalcFiles", CF_PRW_lcalcFiles));
  CHECK(m_objTool->initialize().isSuccess());



  return EL::StatusCode::SUCCESS;
}



EL::StatusCode aDAODMaker :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.
  xAOD::TEvent* event = wk()->xaodEvent();
//   CHECK(event->copy("Muons"));

//     CHECK(m_objTool->ApplyPRWTool());

    // Get the Electrons from the event
    xAOD::ElectronContainer* electrons_copy(0);
    xAOD::ShallowAuxContainer* electrons_copyaux(0);
    electrons_copyaux->setShallowIO( false );
    CHECK(m_objTool->GetElectrons(electrons_copy,electrons_copyaux));
    CHECK(event->record(electrons_copy   , "ShallowCopiedElectrons"));
    CHECK(event->record(electrons_copyaux, "ShallowCopiedElectronsAux."));

//     CHECK(wk()->xaodStore()->record(electrons_copy   , "ShallowCopiedElectrons"));
//     CHECK(wk()->xaodStore()->record(electrons_copyaux, "ShallowCopiedElectronsAux."));

    /// photon Check
//     xAOD::PhotonContainer* photons_copy(0);
//     xAOD::ShallowAuxContainer* photons_copyaux(0);
//     CHECK( m_objTool->GetPhotons(photons_copy,photons_copyaux) );
//     CHECK(wk()->xaodStore()->record(photons_copy       , "ShallowCopiedPhotons"));
//     CHECK(wk()->xaodStore()->record(photons_copyaux    , "ShallowCopiedPhotonsAux."));
// 
//     /// muon Check
//     xAOD::MuonContainer* muons_copy(0);
//     xAOD::ShallowAuxContainer* muons_copyaux(0);
//     CHECK( m_objTool->GetMuons(muons_copy,muons_copyaux) );
//     CHECK(wk()->xaodStore()->record(muons_copy       , "ShallowCopiedMuons"));
//     CHECK(wk()->xaodStore()->record(muons_copyaux    , "ShallowCopiedMuonsAux."));
// 
//     ///jet
//     xAOD::JetContainer* jets_copy(0);
//     xAOD::ShallowAuxContainer* jets_copyaux(0);
//     CHECK( m_objTool->GetJets(jets_copy,jets_copyaux) );
//     CHECK(wk()->xaodStore()->record(jets_copy        , "ShallowCopiedJets"));
//     CHECK(wk()->xaodStore()->record(jets_copyaux     , "ShallowCopiedJetsAux."));
// 
//     /// overlap removal
//     CHECK( m_objTool->OverlapRemoval(electrons_copy, muons_copy, jets_copy) );



  /// Need electron, muon, tau, photon, jets, MET, truth?

  /// do selection
  

  event->fill();

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode aDAODMaker :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode aDAODMaker :: finalize ()
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

  // finalize and close our output xAOD file:
  xAOD::TEvent* event = wk()->xaodEvent();
  TFile *file = wk()->getOutputFile("outputLabel");
  CHECK(event->finishWritingTo(file));

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode aDAODMaker :: histFinalize ()
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
