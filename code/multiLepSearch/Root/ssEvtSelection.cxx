#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <EventLoop/OutputStream.h>
#include <multiLepSearch/ssEvtSelection.h>
#include <multiLepSearch/anaHelper.h>
#include "xAODMissingET/MissingETContainer.h"
#include "xAODMissingET/MissingETAuxContainer.h"
// #include "xAODTruth/TruthParticleContainer.h"
// #include "xAODTruth/TruthParticle.h"
#include "xAODMuon/Muon.h"
#include "GoodRunsLists/GoodRunsListSelectionTool.h"
#include "xAODEventInfo/EventInfo.h"
#include "xAODCutFlow/CutBookkeeperContainer.h"
#include "xAODCutFlow/CutBookkeeper.h"
#include "PathResolver/PathResolver.h"
#include "xAODRootAccess/tools/Message.h"
// #include "MuonSelectorTools/MuonSelectionTool.h"
#include "xAODTruth/TruthParticle.h"

// #include "PATInterfaces/SystematicVariation.h" 
// #include "PATInterfaces/SystematicsUtil.h"

#include "xAODCore/tools/IOStats.h"
#include "xAODCore/tools/ReadStats.h"

#include <multiLepSearch/ChargeFlipBkgTool.h>
#include <multiLepSearch/FakeLepBkgTool.h>

#include <TError.h>
#include <algorithm>
#include <vector>
#include <multiLepSearch/mCHECK.h>

#include <TH1D.h>
#include <TFile.h>

using namespace xAOD;
using namespace std;

// static SG::AuxElement::Decorator<char> dec_quality("quality");
static SG::AuxElement::Accessor<char> dec_baseline("baseline");
static SG::AuxElement::Accessor<char> dec_signal("signal");
static SG::AuxElement::Accessor<char> dec_passOR("passOR");
static SG::AuxElement::Accessor<char> dec_bad("bad");
static SG::AuxElement::Accessor<char> dec_bjet_loose("bjet_loose");
static SG::AuxElement::Accessor<char> dec_bjet("bjet");
static SG::AuxElement::Accessor<char> dec_cosmic("cosmic");
typedef ElementLink< xAOD::TruthParticleContainer > TruthLink;
//static SG::AuxElement::Accessor< int > acc_truthType("truthType");
//static SG::AuxElement::Accessor< int > acc_truthOrig("truthOrigin");
// static SG::AuxElement::Accessor< float > acc_truthProb("truthMatchProbability"); // only ID track
static SG::AuxElement::Accessor< TruthLink > acc_truthLink("truthParticleLink"); // ID track, electron
//static SG::AuxElement::Accessor<float> dec_z0sinTheta("z0sinTheta");
//static SG::AuxElement::Accessor<float> dec_d0sig("d0sig");

const float iGeV = 0.001;
const float GeV = 1000;


// this is needed to distribute the algorithm to the workers
ClassImp(ssEvtSelection)



ssEvtSelection :: ssEvtSelection(string name):m_name(name),m_susyEvt(0),m_grl(0){
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().

  //m_susyEvt = new susyEvts();

  CF_vxTrkNMin = 0;
  CF_vxTrkPtMin = -1*GeV;
  CF_outputName = "test.root";
  CF_outputTreeName = "evt3l";
  CF_derivationName = "SUSY2";
  //CF_nLepCutExactly = 2;
  //CF_nLepCutMin = 2;
  //CF_nJetCutExactly = 2;
  //CF_nJetCutMin = 2;

  //CF_ElPtCut = 20000;
  //CF_Eld0SigCut = 5;
  //CF_Elz0Cut = 0.5;

  //CF_MuPtCut = 20000;
  //CF_Mud0SigCut = 3;
  //CF_Muz0Cut = 0.5;

  //CF_JetPtCut = 20000;
  //CF_JetEtaCut = 2.8;
  //CF_JetJvtCut = 0.64;

  CF_mT2_m0 = 1;

  CF_isMC = 0;
  study = "3l";
  doSys = 0;

  // MCTC, TruthLink, dR
  mcTruthMatch = "MCTC"; 

  // ElectronChargeIDSelector working points
  // Defined here: https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/ElectronChargeFlipTaggerTool
  ECIDS_OP=-0.28087;
  ECIDS_trainingFile="ElectronPhotonSelectorTools/ChargeID/ECIDS_20161125for2017Moriond.root";
}



EL::StatusCode ssEvtSelection :: setupJob (EL::Job& job)
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



EL::StatusCode ssEvtSelection :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.
  TFile *outputFile = wk()->getOutputFile(CF_outputName);
  //m_susyEvt->makeTree(CF_outputTreeName);     //moved to initialize()
  //m_susyEvt->tree2->SetDirectory(outputFile);

  m_hCutFlow = new TH1D("hCutFlow", "cut flow", 60, 0, 60);
  m_hCutFlow->SetDirectory(outputFile);
  // m_hCutFlow->SetDirectory(0);
  // m_susyEvt->tree2->GetUserInfo()->Add(m_hCutFlow);
  m_hCutFlow->GetXaxis()->SetBinLabel(1,"AOD");
  m_hCutFlow->GetXaxis()->SetBinLabel(2,"aSumW");
  m_hCutFlow->GetXaxis()->SetBinLabel(3,"aSumW2");
  m_hCutFlow->GetXaxis()->SetBinLabel(4,"SUSY2");
  m_hCutFlow->GetXaxis()->SetBinLabel(5,"dSumW");
  m_hCutFlow->GetXaxis()->SetBinLabel(6,"dSumW2");
  m_hCutFlow->GetXaxis()->SetBinLabel(7,"All");
  m_hCutFlow->GetXaxis()->SetBinLabel(8,"GRL");
  m_hCutFlow->GetXaxis()->SetBinLabel(9,"LArErr");
  m_hCutFlow->GetXaxis()->SetBinLabel(10,"TileErr");
  m_hCutFlow->GetXaxis()->SetBinLabel(11,"CoreFlag");
  m_hCutFlow->GetXaxis()->SetBinLabel(12,"PV");
  m_hCutFlow->GetXaxis()->SetBinLabel(13,"BadMuon");
  m_hCutFlow->GetXaxis()->SetBinLabel(14,"cosMuon");
  m_hCutFlow->GetXaxis()->SetBinLabel(15,"BadJet");
  m_hCutFlow->GetXaxis()->SetBinLabel(16,"nSumW");
  m_hCutFlow->GetXaxis()->SetBinLabel(20,"=2BaseLep");
  m_hCutFlow->GetXaxis()->SetBinLabel(21,"=2SigLep");
  m_hCutFlow->GetXaxis()->SetBinLabel(30,">=1BaseEl");
  m_hCutFlow->GetXaxis()->SetBinLabel(31,">=1SigEl");
  m_hCutFlow->GetXaxis()->SetBinLabel(32,">=1BaseMu");
  m_hCutFlow->GetXaxis()->SetBinLabel(33,">=1SigMu");
  m_hCutFlow->GetXaxis()->SetBinLabel(34,">=1BaseLep");
  m_hCutFlow->GetXaxis()->SetBinLabel(35,">=1SigLep");
  m_hCutFlow->GetXaxis()->SetBinLabel(36,">=2BaseEl");
  m_hCutFlow->GetXaxis()->SetBinLabel(37,">=2SigEl");
  m_hCutFlow->GetXaxis()->SetBinLabel(38,">=2BaseMu");
  m_hCutFlow->GetXaxis()->SetBinLabel(39,">=2SigMu");
  m_hCutFlow->GetXaxis()->SetBinLabel(40,">=2BaseLep");
  m_hCutFlow->GetXaxis()->SetBinLabel(41,">=2SigLep");
  m_hCutFlow->GetXaxis()->SetBinLabel(42,">=1BaseJet");
  m_hCutFlow->GetXaxis()->SetBinLabel(43,">=1SigJet");
  m_hCutFlow->GetXaxis()->SetBinLabel(44,">=1BJet");
  m_hCutFlow->GetXaxis()->SetBinLabel(45,"=3BaseEl");
  m_hCutFlow->GetXaxis()->SetBinLabel(46,"=3SigEl");
  m_hCutFlow->GetXaxis()->SetBinLabel(47,"=3BaseMu");
  m_hCutFlow->GetXaxis()->SetBinLabel(48,"=3SigMu");
  m_hCutFlow->GetXaxis()->SetBinLabel(49,"=3BaseLep");
  m_hCutFlow->GetXaxis()->SetBinLabel(50,"=3SigLep");
  m_hCutFlow->GetXaxis()->SetBinLabel(51,"=3BaseLep and =3SigLep");

  m_hCutFlowNominal = m_hCutFlow;
  m_hCutFlowDummy = new TH1D( *(TH1D*)m_hCutFlow );
  m_hCutFlowDummy->SetDirectory(0); //prevent saving this hist, set it not belonging to current opend file/dir

  m_hTrigs = new TH1D("hTrigs", "n pass trigger", CF_trigNames.size(), 0, CF_trigNames.size());
  for(unsigned int i=0; i<CF_trigNames.size(); i++){
    m_hTrigs->GetXaxis()->SetBinLabel(i+1,CF_trigNames[i].c_str());
    std::cout<<CF_trigNames[i].c_str()<<std::endl;
  }
  m_hTrigs->SetDirectory(outputFile);
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ssEvtSelection :: fileExecute ()
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

    m_hCutFlow = m_hCutFlowNominal;

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



EL::StatusCode ssEvtSelection :: changeInput (bool /*firstFile*/)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ssEvtSelection :: initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.

  // as a check, let's see the number of events in our xAOD
  Info("initialize()", "Number of events = %lli", wk()->xaodEvent()->getEntries() ); // print long long int

  // count number of events
  m_eventCounter = 0;

  /// IsolationSelectionTool
  m_isoTool = new CP::IsolationSelectionTool("my_isoTest");
  CHECK(m_isoTool->setProperty("MuonWP","Loose"));
  CHECK(m_isoTool->setProperty("ElectronWP","Loose"));
  CHECK(m_isoTool->initialize());

  CHECK(m_isoTool->addWP("Gradient", xAOD::Type::Muon)); //2
  CHECK(m_isoTool->addWP("GradientLoose", xAOD::Type::Muon)); //4
  CHECK(m_isoTool->addWP("LooseTrackOnly", xAOD::Type::Muon)); //8
  CHECK(m_isoTool->addWP("FixedCutTight", xAOD::Type::Muon)); //10
  CHECK(m_isoTool->addWP("FixedCutTightTrackOnly", xAOD::Type::Muon)); //20
  CHECK(m_isoTool->addWP("FixedCutLoose", xAOD::Type::Muon)); //40
  CHECK(m_isoTool->addWP("FixedCutHighPtTrackOnly", xAOD::Type::Muon));//80
  CHECK(m_isoTool->addWP("Gradient", xAOD::Type::Electron)); //2
  CHECK(m_isoTool->addWP("GradientLoose", xAOD::Type::Electron)); //4
  CHECK(m_isoTool->addWP("LooseTrackOnly", xAOD::Type::Electron)); //8
  CHECK(m_isoTool->addWP("FixedCutTight", xAOD::Type::Electron)); //10
  CHECK(m_isoTool->addWP("FixedCutTightTrackOnly", xAOD::Type::Electron)); //20
  CHECK(m_isoTool->addWP("FixedCutLoose", xAOD::Type::Electron)); //40
  CHECK(m_isoTool->addWP("FixedCutHighPtCaloOnly", xAOD::Type::Electron)); //80
  CHECK(m_isoTool->addWP("FixedCutTrackCone40", xAOD::Type::Electron)); //100

  /// SUSYTools
  m_objTool = new ST::SUSYObjDef_xAOD("SUSYObjDef_xAOD");
  //m_objTool->msg().setLevel( MSG::ERROR );
  //m_objTool->msg().setLevel( MSG::WARNING );
//   m_objTool->msg().setLevel( MSG::VERBOSE );

  for(auto x: CF_PRW_confFiles){Info("CF_PRW_confFiles", "%s", x.c_str());}

  ST::ISUSYObjDef_xAODTool::DataSource ds = static_cast<ST::ISUSYObjDef_xAODTool::DataSource>(CF_isMC); 
  CHECK(m_objTool->setProperty("DataSource",ds));
  CHECK(m_objTool->setProperty("ConfigFile", CF_ConfigFile));
  CHECK(m_objTool->setProperty("PRWConfigFiles", CF_PRW_confFiles));
  CHECK(m_objTool->setProperty("PRWLumiCalcFiles", CF_PRW_lcalcFiles));
  CHECK(m_objTool->initialize().isSuccess());

  mChargeFlipBkgTool = new ChargeFlipBkgTool("MyQFlipTool");
  CHECK(mChargeFlipBkgTool->setProperty("InputRatesFileName" , "$ROOTCOREBIN/data/multiLepSearch/root_files/chargeFlipRates.root"));
  CHECK(mChargeFlipBkgTool->setProperty("InputRatesHistoName", "hFlipProb_data"));
  CHECK(mChargeFlipBkgTool->initialize());

  mFakeLepBkgTool = new FakeLepBkgTool("MyFLepTool");
  //CHECK(mFakeLepBkgTool->setProperty("Method", "Matrix"));
  //CHECK(mFakeLepBkgTool->setProperty("InputFileName"    , "$ROOTCOREBIN/data/multiLepSearch/root_files/RealFakeLepEff_dummy.root"));
  //CHECK(mFakeLepBkgTool->setProperty("RealeEffHistoName", "RealeEff"));
  //CHECK(mFakeLepBkgTool->setProperty("RealuEffHistoName", "RealuEff"));
  //CHECK(mFakeLepBkgTool->setProperty("FakeeEffHistoName", "FakeeEff"));
  //CHECK(mFakeLepBkgTool->setProperty("FakeuEffHistoName", "FakeuEff"));

  CHECK(mFakeLepBkgTool->setProperty("Method", "FakeFactor"));
  CHECK(mFakeLepBkgTool->setProperty("InputFileName"    , "$ROOTCOREBIN/data/multiLepSearch/root_files/fakefactor_2D_Data16.root"));
  CHECK(mFakeLepBkgTool->setProperty("eFakeFactorHistoName", "h_ff_ele"));
  //CHECK(mFakeLepBkgTool->setProperty("eFakeFactorHistoName", "h_ff_ele_v2")); //this histo has problem of bin error being just the sqrt of bin content
  CHECK(mFakeLepBkgTool->setProperty("uFakeFactorHistoName", "h_ff_mu"));

  CHECK(mFakeLepBkgTool->initialize());

  /// GRL
  if(!CF_isMC){
    m_grl = new GoodRunsListSelectionTool("GoodRunsListSelectionTool");
    vector<string> vecStringGRL;
    for(unsigned int i=0;i<CF_grlFiles.size();i++)
    {
      vecStringGRL.push_back(PathResolverFindCalibFile(CF_grlFiles[i]));
    }
    CHECK(m_grl->setProperty("GoodRunsListVec", vecStringGRL));
    CHECK(m_grl->setProperty("PassThrough", false)); // if true (default) will ignore result of GRL and will just pass all events
    CHECK(m_grl->initialize());
  }

  m_truthClassifier = new MCTruthClassifier("m_truthClassifier");

  //prepare list of systematics to do
  TFile *outputFile = wk()->getOutputFile(CF_outputName);
  std::vector<ST::SystInfo> fullSystInfoList = m_objTool->getSystInfoList();
  m_systInfoList.clear();

  // --- find Nominal (no syst) config and place it at position 0
  for (auto& aSystInfo : fullSystInfoList){
    if(aSystInfo.systset.name()==""){
      m_susyEvt = new susyEvts();
      m_susyEvt->makeTree(CF_outputTreeName);
      m_susyEvt->tree2->SetDirectory(outputFile);
      m_systInfoList.push_back( aSystInfo );
      m_susyEvtList.push_back( m_susyEvt );
      break;
    }
  }

  // --- append other sys if doSys is true, also sys in systInfoList only applies on MC
  if (doSys && CF_isMC){
    for (auto& aSystInfo : fullSystInfoList){
      if(aSystInfo.systset.name()==""){ continue; }
      //TODO : reject some sys if it is in a user provided blacklist
      if(aSystInfo.systset.name().find("TAU")!=std::string::npos ){ continue; }
      if(aSystInfo.systset.name().find("JET")!=std::string::npos ){ continue; }

      //skip these as right now this don't work with any mu Trigger
      //missing histo in /cvmfs/atlas.cern.ch/repo/sw/database/GroupData/MuonEfficiencyCorrections/160624_ICHEP/muontrigger_sf_2015_mc15c_v01.root
      if(aSystInfo.systset.name()=="MUON_EFF_TrigSystUncertainty__1down"){ continue; } 
      if(aSystInfo.systset.name()=="MUON_EFF_TrigSystUncertainty__1up"  ){ continue; }
      
      m_susyEvt = new susyEvts();
      if(aSystInfo.affectsKinematics){
        m_susyEvt->makeKinematicsSysTree(CF_outputTreeName + "_" + aSystInfo.systset.name(), m_susyEvtList[0]);
      }else{
        m_susyEvt->makeWeightOnlyTree(CF_outputTreeName + "_" + aSystInfo.systset.name(), m_susyEvtList[0]);
      }
      //m_susyEvt->tree2->SetDirectory(outputFile); //done in makeXXXXTree
      m_systInfoList.push_back( aSystInfo );
      m_susyEvtList.push_back( m_susyEvt );
    }
  }

  // Todo : add sys other than those m_objTool(SUSYTool) recommends, eg. theoretical sys, data bkg sys?

  // --- print out accepted sys
  std::cout << "======  Systematic list  ======" << std::endl;
  for (auto& aSystInfo : m_systInfoList){
    if(aSystInfo.systset.name()==""){
      std::cout << "Nominal (no syst) "  << std::endl;
    }else{
      std::cout << "Systematic: " << aSystInfo.systset.name() << std::endl;
    }
  }
  std::cout << "===============================" << std::endl;


  setupTriggers();


  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ssEvtSelection :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.
  auto sc = EL::StatusCode::SUCCESS;

  m_hCutFlow = m_hCutFlowNominal;
  m_hCutFlow->Fill("All", 1);

  // print every 100 events, so we know where we are:
  if((m_eventCounter % 1000)==0) Info("execute()", "Event number = %i", m_eventCounter );
  m_eventCounter++;

  // loop over systematics
  for (unsigned int iSyst=0; iSyst<m_systInfoList.size();iSyst++){
    
    if(iSyst>0 && !m_systInfoList[iSyst].affectsKinematics) continue; //weight only syst are treated elsewhere, see end of this iSyst loop

    //---------------------------------------------------------------------------------------
    // config tools to apply the current sys
    //---------------------------------------------------------------------------------------
    CP::SystematicCode ret(CP::SystematicCode::Ok);

    if (m_objTool) {
      ret  = m_objTool->applySystematicVariation(m_systInfoList[iSyst].systset);
      if ( ret != CP::SystematicCode::Ok) {
        ATH_MSG_ERROR("Cannot configure SUSYObjDefTool for systematic var. " << m_systInfoList[iSyst].systset.name() );
        continue;
      } else {
        ATH_MSG_VERBOSE("SUSYObjDef configured for systematic var. " << m_systInfoList[iSyst].systset.name() );
      }
    }

    // set the pointer to the correct object
    m_susyEvt = m_susyEvtList[iSyst];
    if (iSyst!=0) m_hCutFlow = m_hCutFlowDummy; //to avoid repeated filling of cutflow 



    //----------------------------
    // Event information
    //--------------------------- 
    const xAOD::EventInfo* eventInfo(nullptr);
    if( ! wk()->xaodEvent()->retrieve(eventInfo, "EventInfo").isSuccess()){
      Error("execute()", "Failed to retrieve event info collection. Exiting.");
      return sc;
    }

    /// GRL
    if(CF_isMC || m_grl->passRunLB(*eventInfo)){
      m_hCutFlow->Fill("GRL", 1);
    }else return sc;

    /// LArError, tileError, CoreFlags
    if(eventInfo->errorState(xAOD::EventInfo::LAr)!=xAOD::EventInfo::Error){
      m_hCutFlow->Fill("LArErr", 1);
    }
    else return sc;

    if(eventInfo->errorState(xAOD::EventInfo::Tile)!=xAOD::EventInfo::Error){
      m_hCutFlow->Fill("TileErr", 1);
    }
    else return sc;

    if(eventInfo->isEventFlagBitSet(xAOD::EventInfo::Core, 18)) {
    //if(eventInfo->eventFlags(EventInfo::Core) & 0x40000) 
      ATH_MSG_WARNING("This event is incompletely built. Skipping.");
      return sc;
    }
    m_hCutFlow->Fill("CoreFlag", 1);

    // Primary vertex
    unsigned int nPV=0;
    //m_susyEvt->sig.nVtx=0;
    const xAOD::VertexContainer* vertices(0);
    if( wk()->xaodEvent()->retrieve( vertices, "PrimaryVertices" ).isSuccess() ) {
      for( const auto& vx : *vertices ) {
        //m_susyEvt->sig.nVtx++;
        if(vx->vertexType() != xAOD::VxType::PriVtx) continue;
        int nTrk = 0;
        for(size_t i=0; i<vx->nTrackParticles(); i++){
          const xAOD::TrackParticle* trk = vx->trackParticle(i);
          if(trk && trk->pt()>CF_vxTrkPtMin) nTrk++;
        }
        if(nTrk<CF_vxTrkNMin) continue;
        nPV++;
      }
    }else{
      Error("Failed to retrieve VertexContainer %s, returning NULL", "PrimaryVertices");
    }

    if(nPV<1) return sc;
    m_hCutFlow->Fill("PV", 1);

    CHECK(m_objTool->ApplyPRWTool());

    // Get the Electrons from the event
    xAOD::ElectronContainer* electrons_copy(0);
    xAOD::ShallowAuxContainer* electrons_copyaux(0);
    CHECK(m_objTool->GetElectrons(electrons_copy,electrons_copyaux));
    CHECK(wk()->xaodStore()->record(electrons_copy   , m_systInfoList[iSyst].systset.name()+"ShallowCopiedElectrons"));
    CHECK(wk()->xaodStore()->record(electrons_copyaux, m_systInfoList[iSyst].systset.name()+"ShallowCopiedElectronsAux."));

    /// photon Check
    xAOD::PhotonContainer* photons_copy(0);
    xAOD::ShallowAuxContainer* photons_copyaux(0);
    CHECK( m_objTool->GetPhotons(photons_copy,photons_copyaux) );
    CHECK(wk()->xaodStore()->record(photons_copy       , m_systInfoList[iSyst].systset.name()+"ShallowCopiedPhotons"));
    CHECK(wk()->xaodStore()->record(photons_copyaux    , m_systInfoList[iSyst].systset.name()+"ShallowCopiedPhotonsAux."));

    /// muon Check
    xAOD::MuonContainer* muons_copy(0);
    xAOD::ShallowAuxContainer* muons_copyaux(0);
    CHECK( m_objTool->GetMuons(muons_copy,muons_copyaux) );
    CHECK(wk()->xaodStore()->record(muons_copy       , m_systInfoList[iSyst].systset.name()+"ShallowCopiedMuons"));
    CHECK(wk()->xaodStore()->record(muons_copyaux    , m_systInfoList[iSyst].systset.name()+"ShallowCopiedMuonsAux."));

    ///jet
    xAOD::JetContainer* jets_copy(0);
    xAOD::ShallowAuxContainer* jets_copyaux(0);
    CHECK( m_objTool->GetJets(jets_copy,jets_copyaux) );
    CHECK(wk()->xaodStore()->record(jets_copy        , m_systInfoList[iSyst].systset.name()+"ShallowCopiedJets"));
    CHECK(wk()->xaodStore()->record(jets_copyaux     , m_systInfoList[iSyst].systset.name()+"ShallowCopiedJetsAux."));

    /// overlap removal
    CHECK( m_objTool->OverlapRemoval(electrons_copy, muons_copy, jets_copy) );

    // Bad muon
    bool hasBadMuon = false;
    bool isCosmicMuon = false;
    for(auto mu: *muons_copy){
      if(dec_baseline(*mu) && dec_bad(*mu)){
        hasBadMuon = true;
        break;
      }

      if(dec_baseline(*mu) && dec_passOR(*mu) && dec_cosmic(*mu)){
        isCosmicMuon = true;
      }
    }
    if(hasBadMuon) continue;
    m_hCutFlow->Fill("BadMuon", 1);

    // Cosmic muon
    if(isCosmicMuon) continue;
    m_hCutFlow->Fill("cosMuon", 1);

    /// Bad jet
    bool hasBadJet = false;
    for(auto jet: *jets_copy){
      if(dec_baseline(*jet) && dec_passOR(*jet) && dec_bad(*jet)){
        hasBadJet = true;
        break;
      }
    }
    if(hasBadJet) continue;
    m_hCutFlow->Fill("BadJet", 1);

    //// let's select the two leptons
    vector< IParticle* > sel_Ls;
    vector< IParticle* > sig_Ls;
    sel_Ls.reserve(3);
    sig_Ls.reserve(3);

    int nBaseEl = 0;
    int nSigEl = 0;
    // Electrons
    for(auto el: *electrons_copy){
      //if(el->pt()<CF_ElPtCut) continue;
      if(dec_passOR(*el)){
        if     (dec_signal  (*el)){ sig_Ls.push_back(el);nBaseEl++;nSigEl++;}
        else if(dec_baseline(*el)){ sel_Ls.push_back(el);nBaseEl++;}
      }
    }
   
    int nBaseMu = 0;
    int nSigMu = 0;
    // Muons
    for(auto mu: *muons_copy){
      //if(mu->pt()<CF_MuPtCut) continue;
      if(dec_passOR(*mu)){
        if     (dec_signal  (*mu)){ sig_Ls.push_back(mu);nBaseMu++;nSigMu++;}
        else if(dec_baseline(*mu)){ sel_Ls.push_back(mu);nBaseMu++;}
      }
    }

    int totLs = sel_Ls.size() + sig_Ls.size();

    /// sort the two leptons for signal pairs and double fake
    sort(sel_Ls.begin(), sel_Ls.end(), [](xAOD::IParticle* a, xAOD::IParticle* b)->bool{return a->pt()>b->pt();});
    sort(sig_Ls.begin(), sig_Ls.end(), [](xAOD::IParticle* a, xAOD::IParticle* b)->bool{return a->pt()>b->pt();});

    bool cutflow = false;
    //cutflow = true;

    vector< IParticle* > dilepPair(2, nullptr);
    if(study == "ss"){
      bool keep = false;

      //this catches 2SigLepSS and 2SigLepOS(i.e. charge flip)
      if (sig_Ls.size() == 2) {
	dilepPair[0] = sig_Ls[0];
	dilepPair[1] = sig_Ls[1];
        keep = true;
        m_susyEvt->evt.flag = 1;
      } 

      //this catches 1SigLep1FakeLepSS -.-
      if (sig_Ls.size() == 1) {
	int sigLepSign = 0;
        xAOD::Muon* mu = dynamic_cast<xAOD::Muon*>(sig_Ls[0]);
        if(mu) sigLepSign = mu->charge();
        else{
          xAOD::Electron* el = dynamic_cast<xAOD::Electron*>(sig_Ls[0]);
          if(el) sigLepSign = el->charge();
        }
	dilepPair[0] = sig_Ls[0];

	for (auto p : sel_Ls){
	  int baseLepSign = -999;
          xAOD::Muon* mu = dynamic_cast<xAOD::Muon*>(p);
          if(mu) baseLepSign = mu->charge();
          else{
            xAOD::Electron* el = dynamic_cast<xAOD::Electron*>(p);
            if(el) baseLepSign = el->charge();
          }
          if (baseLepSign==sigLepSign){
	    dilepPair[1] = p;
            keep = true;
            m_susyEvt->evt.flag = 2;
	    break;
	  }
	}
      }

      //this catches 2FakeLepSS
      if (sig_Ls.size() == 0) {
        if (sel_Ls.size() >=2) { 
	  int baseLep0Sign = 0;
          xAOD::Muon* mu = dynamic_cast<xAOD::Muon*>(sel_Ls[0]);
          if(mu) baseLep0Sign = mu->charge();
          else{
            xAOD::Electron* el = dynamic_cast<xAOD::Electron*>(sel_Ls[0]);
           if(el) baseLep0Sign = el->charge();
          }

	  for (unsigned int i = 1; i<sel_Ls.size(); i++){
	    int baseLep1Sign = -999;
            xAOD::Muon* mu = dynamic_cast<xAOD::Muon*>(sel_Ls[1]);
            if(mu) baseLep1Sign = mu->charge();
            else{
              xAOD::Electron* el = dynamic_cast<xAOD::Electron*>(sel_Ls[1]);
              if(el) baseLep1Sign = el->charge();
            }
            if (baseLep0Sign==baseLep1Sign){
	      dilepPair[0] = sel_Ls[0];
	      dilepPair[1] = sel_Ls[1];
  	      keep = true;
              m_susyEvt->evt.flag = 3;
	      break;
	    }
	  }

	  if (dilepPair[0]==nullptr && sel_Ls.size()>2){
	    //no one has same sign as the first lep
	    //=> everyone except the first lep is Same sign among themselves
	    dilepPair[0] = sel_Ls[1];
	    dilepPair[1] = sel_Ls[2];
  	    keep = true;
            m_susyEvt->evt.flag = 3;
	  }
	}
      } 

      if (!cutflow && !keep) {continue;}
    }
    if(study == "3l" && totLs != 3) continue;

    sel_Ls.insert( sel_Ls.begin(), sig_Ls.begin(), sig_Ls.end());
    sort(sel_Ls.begin(), sel_Ls.end(), [](xAOD::IParticle* a, xAOD::IParticle* b)->bool{return a->pt()>b->pt();});

    if(!cutflow && study == "ss"){
      //for ss study we put the dilepPair at front, then others sorted by pT at tail
      sort(dilepPair.begin(), dilepPair.end(), [](xAOD::IParticle* a, xAOD::IParticle* b)->bool{return a->pt()>b->pt();});
      for (auto p : sel_Ls){
        if (p==dilepPair[0] || p==dilepPair[1])continue;
	dilepPair.push_back(p);
      }
      sel_Ls = dilepPair;
    }

    if(study == "ss")
    {
      if(totLs == 2) m_hCutFlow->Fill("=2BaseLep", 1);
      else if(totLs == 3) m_hCutFlow->Fill("=3BaseLep", 1);
      if(sig_Ls.size() == 2) m_hCutFlow->Fill("=2SigLep", 1);
      else if(sig_Ls.size() == 3) m_hCutFlow->Fill("=3SigLep", 1);

      if(nBaseEl >= 1) m_hCutFlow->Fill(">=1BaseEl", 1);
      if(nSigEl >= 1) m_hCutFlow->Fill(">=1SigEl", 1);
      if(nBaseMu >= 1) m_hCutFlow->Fill(">=1BaseMu", 1);
      if(nSigMu >= 1) m_hCutFlow->Fill(">=1SigMu", 1);
      if(totLs >= 1) m_hCutFlow->Fill(">=1BaseLep", 1);
      if(sig_Ls.size() >= 1) m_hCutFlow->Fill(">=1SigLep", 1);

      if(nBaseEl >= 2) m_hCutFlow->Fill(">=2BaseEl", 1);
      if(nSigEl >= 2) m_hCutFlow->Fill(">=2SigEl", 1);
      if(nBaseMu >= 2) m_hCutFlow->Fill(">=2BaseMu", 1);
      if(nSigMu >= 2) m_hCutFlow->Fill(">=2SigMu", 1);
      if(totLs >= 2) m_hCutFlow->Fill(">=2BaseLep", 1);
      if(sig_Ls.size() >= 2) m_hCutFlow->Fill(">=2SigLep", 1);

      if(nBaseEl == 3) m_hCutFlow->Fill("=3BaseEl", 1);
      if(nSigEl == 3) m_hCutFlow->Fill("=3SigEl", 1);
      if(nBaseMu == 3) m_hCutFlow->Fill("=3BaseMu", 1);
      if(nSigMu == 3) m_hCutFlow->Fill("=3SigMu", 1);

      if(totLs == 3 && sig_Ls.size() == 3) m_hCutFlow->Fill("=3BaseLep and =3SigLep", 1);
    }
    /////////////////////////
    // Save informations
    ////////////////////////
    
    //unused variables
    m_susyEvt->evt.cut = 0;

    //event information
    m_susyEvt->evt.event = eventInfo->eventNumber();
    //m_susyEvt->evt.lumiBlock = eventInfo->lumiBlock();
    //m_susyEvt->evt.actualMu = eventInfo->actualInteractionsPerCrossing();
    
    if(!CF_isMC) m_susyEvt->evt.weight = 1;
    else if(fabs(eventInfo->mcEventWeight())<100) m_susyEvt->evt.weight = eventInfo->mcEventWeight();
    else m_susyEvt->evt.weight = 1;
    
    m_susyEvt->evt.isMC = CF_isMC? 1:0;
    m_hCutFlow->Fill("nSumW", m_susyEvt->evt.weight);

    m_susyEvt->truths.clear();

    // charge flip weight and pT correction (have to apply before other calculation that use lep pT)
    m_susyEvt->evt.qFwt = 0.0;
    m_susyEvt->evt.qFwt_sys_1up = 0.0;
    m_susyEvt->evt.qFwt_sys_1dn = 0.0;
    if ((sel_Ls.size()==2)&&(m_susyEvt->evt.flag==1)){
        //ugly code to get lep0 type and charge :(
        xAOD::Electron* tmpE0 = NULL;  xAOD::Muon* tmpMu0 = NULL;
        int sigLepSign0 = 0;
        tmpMu0 = dynamic_cast<xAOD::Muon*>(sel_Ls[0]);
        if(tmpMu0) sigLepSign0 = tmpMu0->charge();
        else{
          tmpE0 = dynamic_cast<xAOD::Electron*>(sel_Ls[0]);
          if(tmpE0) sigLepSign0 = tmpE0->charge();
        }
        //get lep1 type and charge
        xAOD::Electron* tmpE1 = NULL;  xAOD::Muon* tmpMu1 = NULL;
        int sigLepSign1 = 0;
        tmpMu1 = dynamic_cast<xAOD::Muon*>(sel_Ls[1]);
        if(tmpMu1) sigLepSign1 = tmpMu1->charge();
        else{
          tmpE1 = dynamic_cast<xAOD::Electron*>(sel_Ls[1]);
          if(tmpE1) sigLepSign1 = tmpE1->charge();
        }
        if (sigLepSign0!=sigLepSign1){
          m_susyEvt->evt.qFwt = mChargeFlipBkgTool->GetWeight( sel_Ls ,0,0);
          m_susyEvt->evt.qFwt_sys_1up = 0; // mChargeFlipBkgTool->GetWeight( sel_Ls , 1,-1);
          m_susyEvt->evt.qFwt_sys_1dn = 0; // mChargeFlipBkgTool->GetWeight( sel_Ls ,-1,-1);
          auto tmpPt = mChargeFlipBkgTool->GetCorrectedPt( sel_Ls ,0,0);
          if(tmpPt.size()==2){
            if (tmpE0) tmpE0->setP4( tmpPt[0]*1000., tmpE0->eta(), tmpE0->phi(), tmpE0->m());
            if (tmpE1) tmpE1->setP4( tmpPt[1]*1000., tmpE1->eta(), tmpE1->phi(), tmpE1->m());
            //ATH_MSG_ERROR("E0" << tmpPt[0] << " " << sel_Ls[0]->pt());
            //ATH_MSG_ERROR("E1" << tmpPt[1] << " " << sel_Ls[1]->pt());
          }
        }
    }
    //fake lep weight
    m_susyEvt->evt.fLwt = 0.0;
    m_susyEvt->evt.fLwt_e_sys_1up = 0.0; 
    m_susyEvt->evt.fLwt_e_sys_1dn = 0.0; 
    m_susyEvt->evt.fLwt_u_sys_1up = 0.0; 
    m_susyEvt->evt.fLwt_u_sys_1dn = 0.0; 
    if ((sel_Ls.size()==2)&&( (m_susyEvt->evt.flag==2) || (m_susyEvt->evt.flag==3) )){
      m_susyEvt->evt.fLwt = mFakeLepBkgTool->GetWeight(sel_Ls, 0,0);
      m_susyEvt->evt.fLwt_e_sys_1up = mFakeLepBkgTool->GetWeight(sel_Ls, 1,0);
      m_susyEvt->evt.fLwt_e_sys_1dn = mFakeLepBkgTool->GetWeight(sel_Ls,-1,0);
      m_susyEvt->evt.fLwt_u_sys_1up = mFakeLepBkgTool->GetWeight(sel_Ls, 1,1);
      m_susyEvt->evt.fLwt_u_sys_1dn = mFakeLepBkgTool->GetWeight(sel_Ls,-1,1);
      //ATH_MSG_ERROR("FW " << mFakeLepBkgTool->GetWeight(sel_Ls, 0,0));
    }

    /// met
    xAOD::MissingETContainer* metcst = new xAOD::MissingETContainer();
    xAOD::MissingETAuxContainer* metcst_aux = new xAOD::MissingETAuxContainer();
    metcst->setStore(metcst_aux);
    CHECK( m_objTool->GetMET(*metcst,jets_copy,electrons_copy,muons_copy,photons_copy) );

    auto metFinal = (*metcst)["Final"];
    m_susyEvt->sig.Met  = metFinal->met()*iGeV;
    m_susyEvt->sig.MetX = metFinal->mpx()*iGeV;
    m_susyEvt->sig.MetY = metFinal->mpy()*iGeV;
    TLorentzVector metV(metFinal->mpx(), metFinal->mpy(), 0, metFinal->met());

    // jets
    vector< xAOD::IParticle* > jet_Ls;
    jet_Ls.reserve(10);
    for(auto jet: *jets_copy){
      if(cutflow && dec_baseline(*jet) && dec_passOR(*jet)) jet_Ls.push_back(jet);
      if(!cutflow && dec_signal(*jet) && dec_passOR(*jet)) jet_Ls.push_back(jet);
    }
    sort(jet_Ls.begin(), jet_Ls.end(), [](xAOD::IParticle* a, xAOD::IParticle* b)->bool{return a->pt()>b->pt();});
    if(jet_Ls.size() >= 1) m_hCutFlow->Fill(">=1BaseJet", 1);

    m_susyEvt->sig.HT = 0;
    //m_susyEvt->sig.nJet = jet_Ls.size();
    m_susyEvt->jets.resize(jet_Ls.size());
    int nSigJet = 0;
    int nBJet = 0;
    int nISR = 0;
    vector< xAOD::Jet* > cjet_Ls;
    cjet_Ls.reserve(10);
    for(unsigned int i=0;i<jet_Ls.size(); i++){
      auto j = dynamic_cast<xAOD::Jet*>(jet_Ls[i]);
      
      m_susyEvt->jets[i].pt = j->pt()*iGeV; 
      m_susyEvt->jets[i].eta = j->eta(); 
      m_susyEvt->jets[i].phi = j->phi();
      m_susyEvt->jets[i].MET_dPhi = metV.DeltaPhi(j->p4());
      unsigned int& flag = m_susyEvt->jets[i].jFlag;
      flag = 0;
      if(dec_signal(*j)) {flag |= IS_SIGNAL;nSigJet++;}
      if(dec_bjet_loose(*j)) flag |= JT_BJET_LOOSE;
      if(m_objTool->IsBJet(*j)) {flag |= JT_BJET;nBJet++;}
      if(m_susyEvt->jets[i].pt > 40 && fabs(m_susyEvt->jets[i].eta) < 2.4) nISR++;

      //Central jets
      if(fabs(j->eta())<2.4 && j->pt()>20e3 && !(flag & JT_BJET))
      {
        flag |= JT_CJET;
        cjet_Ls.push_back(j);
        m_susyEvt->sig.HT += j->pt()*iGeV;
      }
    }

    if(nSigJet >= 1) m_hCutFlow->Fill(">=1SigJet", 1);
    if(nBJet >= 1) m_hCutFlow->Fill(">=1BJet", 1);
    if(cutflow) continue;

    if(jet_Ls.size()>=2)  m_susyEvt->sig.mjj = (jet_Ls[0]->p4()+jet_Ls[1]->p4()).M() *iGeV;
    else m_susyEvt->sig.mjj = -1;

    //The two leading leptons
    xAOD::IParticle* l1 = dilepPair[0];
    xAOD::IParticle* l2 = dilepPair[1];
    // cout << l1->pt() << " " << l2->pt() << endl;

    m_susyEvt->sig.mlj = -1;
    if(cjet_Ls.size()>=1 && cjet_Ls.size()<=3)
    {
      TLorentzVector JetSystem;
      if(cjet_Ls.size()==1)
      {
        JetSystem = cjet_Ls[0]->p4();
      }
      else if(cjet_Ls.size()==2 || cjet_Ls.size()==3)
      {
        JetSystem = cjet_Ls[0]->p4()+cjet_Ls[1]->p4();
      }

      double dR1 = JetSystem.DeltaR(l1->p4());
      double dR2 = JetSystem.DeltaR(l2->p4());
      xAOD::IParticle* lmindR = nullptr;
      if(dR1<dR2) lmindR = l1;
      else lmindR = l2;
      m_susyEvt->sig.mlj = (JetSystem+lmindR->p4()).M() *iGeV;
    }

    /// save leptons
    m_susyEvt->leps.resize(sel_Ls.size());
    for(unsigned int i=0;i<sel_Ls.size(); i++){
      auto& l = m_susyEvt->leps[i];
      fillLepton(sel_Ls[i], l, i);
      l.MET_dPhi = metV.DeltaPhi(sel_Ls[i]->p4());
      //if(study == "3l") m_susyEvt->leps[i].isTight = (m_susyEvt->leps[i].lFlag & 2)/2;
      /// mT
      auto xt(sel_Ls[i]->p4()+metV); /// sqrt((E_1+E_2)^2-(p_T1+pT_2))
      l.mT = sqrt(xt.Et2()-xt.Perp2())*iGeV;
      
      //// dR with leading jet, and the smallest dR
      l.jet0_dR = jet_Ls.size()>0?jet_Ls[0]->p4().DeltaR(sel_Ls[i]->p4()):-1;
      l.jet_dRm = l.jet0_dR;
      for(size_t j=1; j<jet_Ls.size(); j++){auto dRn=jet_Ls[j]->p4().DeltaR(sel_Ls[i]->p4()); if(dRn<l.jet_dRm) l.jet_dRm = dRn;}
    }

    ///MetRel correction factor (See Eq6, 2LSS note: https://cds.cern.ch/record/1747285/files/ATL-COM-PHYS-2014-954.pdf)
    float minMetdPhi = FLT_MAX; //Met's dPhi from nearest obj (e/mu/jet)
    for(unsigned int i=0;i<sel_Ls.size(); i++){
      minMetdPhi = std::min(minMetdPhi, std::abs(m_susyEvt->leps[i].MET_dPhi));
    }
    for(unsigned int i=0;i<jet_Ls.size(); i++){
      minMetdPhi = std::min(minMetdPhi, std::abs(m_susyEvt->jets[i].MET_dPhi));
    }
    m_susyEvt->sig.MetRel = m_susyEvt->sig.Met;
    if (minMetdPhi<1.570796327) m_susyEvt->sig.MetRel *= sin(minMetdPhi);

    /*
    /// tau
    xAOD::TauJetContainer* taus_copy(0);
    xAOD::ShallowAuxContainer* taus_copyaux(0);
    wk()->xaodStore()->record(taus_copy, "ShallowCopiedTaus");
    wk()->xaodStore()->record(taus_copyaux, "ShallowCopiedTausAux.");
    CHECK( m_objTool->GetTaus(taus_copy,taus_copyaux) );
    m_susyEvt->sig.nTau = 0;
    for(auto tau: *taus_copy){
      if(dec_signal(*tau)) m_susyEvt->sig.nTau++;
    }
    */

    //Assign trigCode
    m_susyEvt->sig.trigCode = 0;
    unsigned long int ADD = 1;
    for(unsigned int i=0; i<CF_trigNames.size(); i++){
      if(m_objTool->IsTrigPassed(CF_trigNames[i])){
        m_susyEvt->sig.trigCode += ADD;
        m_hTrigs->Fill(CF_trigNames[i].c_str(), 1);
      }
      ADD *= 2;
    }

    //l12
    TLorentzVector ll(l1->p4()+l2->p4());
    m_susyEvt->l12.m = ll.M()*iGeV;
    m_susyEvt->l12.pt = ll.Pt()*iGeV; 
    m_susyEvt->l12.eta = ll.Eta(); 
    m_susyEvt->l12.phi = ll.Phi(); 
    m_susyEvt->l12.dPhi = l1->p4().DeltaPhi(l2->p4()); 
    m_susyEvt->l12.dR = l1->p4().DeltaR(l2->p4()); 
    m_susyEvt->l12.MET_dPhi = metV.DeltaPhi(ll);
    if(jet_Ls.size()>=1) m_susyEvt->l12.jet0_dPhi = ll.DeltaPhi(jet_Ls[0]->p4());
    else m_susyEvt->l12.jet0_dPhi = -100;

    m_susyEvt->sig.HT += (l1->pt()+l2->pt())*iGeV;

    /// mT2
    auto tl1 = l1->p4();
    auto tl2 = l2->p4();
    m_susyEvt->sig.mT2 =  anaHelper::get_mT2(tl1.M(), tl1.Px(), tl1.Py(), tl2.M(), tl2.Px(), tl2.Py(), metFinal->mpx(), metFinal->mpy(), CF_mT2_m0*GeV, CF_mT2_m0*GeV)*iGeV;

    //loop over systematics and vary scale factor for weight based systematics
    for (unsigned int jSyst=0; jSyst<m_systInfoList.size();jSyst++){
      
      if (jSyst>0){
          if (iSyst!=0) break; //no need to treat weight based sys for non-nominal tree
          if (!m_systInfoList[jSyst].affectsKinematics){
	    // apply variation for weight based syst.
            if (m_objTool) {
              ret  = m_objTool->applySystematicVariation(m_systInfoList[jSyst].systset);
              if ( ret != CP::SystematicCode::Ok) {
                ATH_MSG_ERROR("Cannot configure SUSYObjDefTool for systematic var. " << m_systInfoList[jSyst].systset.name() );
                continue;
              } else {
                ATH_MSG_VERBOSE("SUSYObjDef configured for systematic var. " << m_systInfoList[jSyst].systset.name() );
              }
            }

            // set the pointer to the correct object
            m_susyEvt = m_susyEvtList[jSyst];
	    m_susyEvt->evt  = m_susyEvtList[0]->evt;  // init value to that of the nominal susyEvt
	    m_susyEvt->leps = m_susyEvtList[0]->leps; // needed by qFwt calculation
	    //m_susyEvt->l12  = m_susyEvtList[0]->l12;
	    //m_susyEvt->jets = m_susyEvtList[0]->jets;
	    //m_susyEvt->truths = m_susyEvtList[0]->truths;
	    //m_susyEvt->sig = m_susyEvtList[0]->sig;
            if (jSyst!=0) m_hCutFlow = m_hCutFlowDummy; //to avoid repeated filling of cutflow 
	  
          }else{
	    //skip non-weight syst as they are treated in the iSyst for loop
            continue;
	  }
      }

      //pileup weights
      if(CF_isMC){ 
        m_susyEvt->evt.pwt = m_objTool->GetPileupWeight();
        m_susyEvt->evt.averageMu = eventInfo->averageInteractionsPerCrossing();
      }else{
        m_susyEvt->evt.pwt = 1;
        m_susyEvt->evt.averageMu = m_objTool->GetCorrectedAverageInteractionsPerCrossing();
      }

      //// get trigger info, if passed, and the SF
      auto& sEvt = m_susyEvt->evt;
      //sEvt.run = eventInfo->runNumber();
      sEvt.run = m_objTool->GetRunNumber();
      auto trigCut = getTriggerConf(sEvt.run);

      //12 channel
      m_susyEvt->evt.flag = 0;
      if(totLs == 2)
      {
        if(TMath::Abs(m_susyEvt->leps[0].ID) == 11000 &&
           TMath::Abs(m_susyEvt->leps[1].ID) == 11000 ){
          m_susyEvt->evt.flag += 1;
          m_susyEvt->sig.trigMask = trigCut->ee_mask;
        }

        else if(TMath::Abs(m_susyEvt->leps[0].ID) == 11000 &&
                TMath::Abs(m_susyEvt->leps[1].ID) == 13000 ){
          m_susyEvt->evt.flag += 3;
          m_susyEvt->sig.trigMask = trigCut->em_mask;
        }

        else if(TMath::Abs(m_susyEvt->leps[0].ID) == 13000 &&
                TMath::Abs(m_susyEvt->leps[1].ID) == 11000 ){
          m_susyEvt->evt.flag += 3;
          m_susyEvt->sig.trigMask = trigCut->em_mask;
        }
        else if(TMath::Abs(m_susyEvt->leps[0].ID) == 13000 &&
                TMath::Abs(m_susyEvt->leps[1].ID) == 13000 ){
          m_susyEvt->evt.flag += 2;
          m_susyEvt->sig.trigMask = trigCut->mm_mask;
        }

        if((m_susyEvt->leps[0].ID > 0 && m_susyEvt->leps[1].ID > 0) ||
           (m_susyEvt->leps[0].ID < 0 && m_susyEvt->leps[1].ID < 0) )
          m_susyEvt->evt.flag += 3;

        if(nISR==1) m_susyEvt->evt.flag += 6;
        else if(nISR!=0) m_susyEvt->evt.flag = 0;
      }

      //Scale factor
      if(CF_isMC){
        if( study == "ss" )
        {
          if(sEvt.flag%3 == 1){ // ee
            //sEvt.ElSF = m_objTool->GetTotalElectronSF(*electrons_copy, true, true, true, true, m_ee_Key, true);
            sEvt.ElSF = m_objTool->GetTotalElectronSF(*electrons_copy, true, true, true, true, m_ee_Key, false);
            sEvt.MuSF = 1;
          }else if(sEvt.flag%3 == 2){ // mumu
            sEvt.ElSF = 1;
            sEvt.MuSF = m_objTool->GetTotalMuonSF(*muons_copy, true, true, trigCut->mmTrig[0]);
          }else if(sEvt.flag%3 == 0){ // emu
            sEvt.ElSF = m_objTool->GetTotalElectronSF(*electrons_copy, true, true, true, true, m_em_eKey);
            sEvt.MuSF = m_objTool->GetTotalMuonSF(*muons_copy, true, true, m_em_mKey);
          }
        }
        else if(study == "3l")
        {
          sEvt.ElSF = m_objTool->GetTotalElectronSF(*electrons_copy,true,true,false,true,"HLT_mu24_iloose_L1MU15");
          sEvt.MuSF = m_objTool->GetTotalMuonSF(*muons_copy,true,true,"HLT_mu24_iloose_L1MU15");
        }
        sEvt.BtagSF = m_objTool->BtagSF(jets_copy);
      }else{
        sEvt.ElSF = 1;
        sEvt.MuSF = 1;
        sEvt.BtagSF = 1;
      }


      /// fill events
      ATH_MSG_VERBOSE("Fill " << iSyst << " " << jSyst );
      m_susyEvt->fill();
    }
    delete metcst;
    delete metcst_aux;

  } //end iSyst for loop
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ssEvtSelection :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc. 

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode ssEvtSelection :: finalize ()
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

//  xAOD::IOStats::instance().stats().printSmartSlimmingBranchList();

 for(auto x: m_trigSel) delete x;
 m_trigSel.clear();

 return EL::StatusCode::SUCCESS;
}



EL::StatusCode ssEvtSelection :: histFinalize ()
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

EL::StatusCode ssEvtSelection :: fillLepton(xAOD::Electron* el, L_PAR& l, unsigned int index)
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
  l.isoPass = m_isoTool->accept(*el).getCutResultBitSet().to_ulong();

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
    if (mcTruthMatch == "TruthLink") {
      //l.truthType = acc_truthType(*el);
      //l.truthOrig = acc_truthOrig(*el);

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

    else if (mcTruthMatch == "MCTC"){
      //std::pair<MCTruthPartClassifier::ParticleType, MCTruthPartClassifier::ParticleOrigin> res;
      //res = m_truthClassifier->particleTruthClassifier(el);

      //l.truthType = res.first;
      //l.truthOrig = res.second;

      auto truthP = m_truthClassifier->getGenPart();
      if (truthP) {
        l.truthI = addTruthPar(truthP, m_susyEvt->truths, -1);
        m_susyEvt->truths[l.truthI].matchI = index;
      } else l.truthI = -1;
    }

    else if (mcTruthMatch == "dR"){
      xAOD::TruthParticle *tp = 0;

      // Set up truth particle container
      const xAOD::TruthParticleContainer* tContainer = 0;
      CHECK(wk()->xaodEvent()->retrieve( tContainer, "TruthParticles" ));
   
      // Momentum vectors
      TLorentzVector p_tp(0,0,0,0);
      TLorentzVector p_e = el->p4();
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
      if (tp){
        l.truthI = addTruthPar(tp, m_susyEvt->truths, -1);
        m_susyEvt->truths[l.truthI].matchI = index;
      } else l.truthI = -1;
    }
  }

  return EL::StatusCode::SUCCESS;
}

EL::StatusCode ssEvtSelection :: fillLepton(xAOD::Muon* mu, L_PAR& l, unsigned int index)
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
  l.isoPass = m_isoTool->accept(*mu).getCutResultBitSet().to_ulong();

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
    //l.truthI = tl.isValid()? addTruthPar(*tl, m_susyEvt->truths, -1):-1;
    //if(l.truthI>=0) m_susyEvt->truths[l.truthI].matchI = index;

    if(tl.isValid())
    {
      l.truthI = addTruthPar(*tl, m_susyEvt->truths, -1);
      m_susyEvt->truths[l.truthI].matchI = index;
      }else l.truthI = -1;
  }
  fillLeptonCommon(mu, l);
  return EL::StatusCode::SUCCESS;
}

void ssEvtSelection :: fillLeptonCommon(xAOD::IParticle* p, L_PAR& l)
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

EL::StatusCode ssEvtSelection :: fillLepton(xAOD::IParticle* p, L_PAR& l, unsigned int index)
{
  xAOD::Muon* mu = dynamic_cast<xAOD::Muon*>(p);
  if(mu) fillLepton(mu, l, index);
  else{
    xAOD::Electron* el = dynamic_cast<xAOD::Electron*>(p);
    if(el) fillLepton(el, l, index);
  }
  return EL::StatusCode::SUCCESS;
}
int ssEvtSelection::addTruthPar(const xAOD::TruthParticle* p, TRUTHS& v, int pLevel){
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

void ssEvtSelection::setupTriggers(){
  /// try di lepton trigger for now
  /// 2015
  m_trigSel.push_back(new TRIGCONF{276073,284484,{"HLT_2e12_lhloose_L12EM10VH"},{"HLT_e17_lhloose_mu14"},{"HLT_mu18_mu8noL1"},0,0,0}); /// 2015 data
  /// 2016: A-D3
//   m_trigSel.push_back(new TRIGCONF{296939,302872,{"HLT_2e15_lhvloose_nod0_L12EM13VH"},{"HLT_e17_lhloose_nod0_mu14"},{"HLT_mu20_mu8noL1"},0,0,0});
  m_trigSel.push_back(new TRIGCONF{296939,302872,{"HLT_2e17_lhvloose_nod0"},{"HLT_e17_lhloose_nod0_mu14"},{"HLT_mu20_mu8noL1"},0,0,0});
  /// 2016: D4-
  m_trigSel.push_back(new TRIGCONF{302919,311481 ,{"HLT_2e17_lhvloose_nod0"},{"HLT_e17_lhloose_nod0_mu14"},{"HLT_mu22_mu8noL1"},0,0,0});

  m_ee_Key = "DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0";
  m_em_eKey = "MULTI_L_2015_e17_lhloose_2016_e17_lhloose_nod0";
  m_em_mKey = "HLT_mu14";

//   /// 2016: A
//   m_trigSel.push_back(new TRIGCONF{296939,300287,{"HLT_2e15_lhvloose_nod0_L12EM13VH"},{"HLT_e17_lhloose_nod0_mu14"},{"HLT_mu20_mu8noL1"},0,0,0});
// 
//   /// 2016: B-D3
//   m_trigSel.push_back(new TRIGCONF{300345,302872,{"HLT_2e15_lhvloose_nod0_L12EM13VH"},{"HLT_e17_lhloose_nod0_mu14"},{"HLT_mu20_mu8noL1"},0,0,0});
// 
//   /// 2016: D4-E
//   m_trigSel.push_back(new TRIGCONF{302919,303892 ,{"HLT_2e17_lhvloose_nod0"},{"HLT_e17_lhloose_nod0_mu14"},{"HLT_mu22_mu8noL1"},0,0,0});
// 
//   /// 2016: F
//   m_trigSel.push_back(new TRIGCONF{303943,304494 ,{"HLT_2e17_lhvloose_nod0"},{"HLT_e17_lhloose_nod0_mu14"},{"HLT_mu22_mu8noL1"},0,0,0});
// 
//   /// 2016: G1-G2
//   m_trigSel.push_back(new TRIGCONF{305291,305293 ,{"HLT_2e17_lhvloose_nod0"},{"HLT_e17_lhloose_nod0_mu14"},{"HLT_mu22_mu8noL1"},0,0,0});
// 
//   /// 2016: G3-I3
//   m_trigSel.push_back(new TRIGCONF{305380,307601 ,{"HLT_2e17_lhvloose_nod0"},{"HLT_e17_lhloose_nod0_mu14"},{"HLT_mu22_mu8noL1"},0,0,0});
// 
//   /// 2016: I4-
//   m_trigSel.push_back(new TRIGCONF{307619,311481 ,{"HLT_2e17_lhvloose_nod0"},{"HLT_e17_lhloose_nod0_mu14"},{"HLT_mu22_mu8noL1"},0,0,0});


//   /// 2016: temp
//   m_trigSel.push_back(new TRIGCONF{-1,-1,{"", ""},{"",""},{"",""}});

  for(auto& t: m_trigSel){
    unsigned long int m1(1);
    for(auto& x: CF_trigNames){
      if(find(t->eeTrig.begin(),t->eeTrig.end(),x)!=t->eeTrig.end()) t->ee_mask |= m1;
      if(find(t->emTrig.begin(),t->emTrig.end(),x)!=t->emTrig.end()) t->em_mask |= m1;
      if(find(t->mmTrig.begin(),t->mmTrig.end(),x)!=t->mmTrig.end()) t->mm_mask |= m1;
      m1 *= 2;
    }
  }

  return;
}

TRIGCONF* ssEvtSelection::getTriggerConf(uint32_t run){
  /// use the cached one if possible to save time
  if(m_nowTrigSel && run>=m_nowTrigSel->runStart && run<=m_nowTrigSel->runEnd) return m_nowTrigSel;

  /// if not cached, fine the new one
  m_nowTrigSel = nullptr;
  for(auto& x: m_trigSel){
    if(run>=x->runStart && run<=x->runEnd){
      m_nowTrigSel = x;
      break;
    }
  }

  return m_nowTrigSel; 
}
