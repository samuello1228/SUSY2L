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
static SG::AuxElement::Accessor<int> acc_truthType("truthType");
static SG::AuxElement::Accessor<int> acc_truthOrig("truthOrigin");
// static SG::AuxElement::Accessor< float > acc_truthProb("truthMatchProbability"); // only ID track
static SG::AuxElement::Accessor< TruthLink > acc_truthLink("truthParticleLink"); // ID track, electron
typedef ElementLink< xAOD::ElectronContainer > RecoElLink;
static SG::AuxElement::Accessor< RecoElLink > acc_recoElLink("recoElectronLink");
typedef ElementLink< xAOD::MuonContainer > RecoMuLink;
static SG::AuxElement::Accessor< RecoMuLink > acc_recoMuLink("recoMuonLink");
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

  CF_outputName = "test.root";
  CF_outputTreeName = "evt2l";
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
  study = "ss";
  SampleName = "";
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

  //list of cutflow
  cutflowList.push_back("AOD");
  cutflowList.push_back("aSumW");
  cutflowList.push_back("aSumW2");
  cutflowList.push_back("SUSY2");
  cutflowList.push_back("dSumW");
  cutflowList.push_back("dSumW2");
  cutflowList.push_back("All");
  cutflowList.push_back("GRL");
  cutflowList.push_back("Trigger");
  cutflowList.push_back("LAr+Tile+SCT+CoreFlag");
  cutflowList.push_back("PV");
  cutflowList.push_back("cosMuon");
  cutflowList.push_back("BadMuon");
  cutflowList.push_back("BadJet");

  cutflowList.push_back("=0BaseLep");
  cutflowList.push_back("=1BaseLep");
  cutflowList.push_back("=2BaseLep");
  cutflowList.push_back("=3BaseLep");
  cutflowList.push_back("=4BaseLep");
  cutflowList.push_back(">=5BaseLep");

  cutflowList.push_back("=0BaseEl");
  cutflowList.push_back("=1BaseEl");
  cutflowList.push_back("=2BaseEl");
  cutflowList.push_back("=3BaseEl");
  cutflowList.push_back("=4BaseEl");
  cutflowList.push_back(">=5BaseEl");

  cutflowList.push_back("=0BaseMu");
  cutflowList.push_back("=1BaseMu");
  cutflowList.push_back("=2BaseMu");
  cutflowList.push_back("=3BaseMu");
  cutflowList.push_back("=4BaseMu");
  cutflowList.push_back(">=5BaseMu");

  cutflowList.push_back("=0SigLep");
  cutflowList.push_back("=1SigLep");
  cutflowList.push_back("=2SigLep");
  cutflowList.push_back("=3SigLep");
  cutflowList.push_back("=4SigLep");
  cutflowList.push_back(">=5SigLep");

  cutflowList.push_back("=0SigEl");
  cutflowList.push_back("=1SigEl");
  cutflowList.push_back("=2SigEl");
  cutflowList.push_back("=3SigEl");
  cutflowList.push_back("=4SigEl");
  cutflowList.push_back(">=5SigEl");

  cutflowList.push_back("=0SigMu");
  cutflowList.push_back("=1SigMu");
  cutflowList.push_back("=2SigMu");
  cutflowList.push_back("=3SigMu");
  cutflowList.push_back("=4SigMu");
  cutflowList.push_back(">=5SigMu");

  cutflowList.push_back("=2BaseLep and =2SigLep");
  cutflowList.push_back("=2BaseLep and =2SigLep,w");

  cutflowList.push_back("=0SigJet");
  cutflowList.push_back("=1SigJet");
  cutflowList.push_back("=2SigJet");
  cutflowList.push_back("=3SigJet");
  cutflowList.push_back("=4SigJet");
  cutflowList.push_back("=5SigJet");
  cutflowList.push_back("=6SigJet");
  cutflowList.push_back("=7SigJet");
  cutflowList.push_back("=8SigJet");
  cutflowList.push_back("=9SigJet");
  cutflowList.push_back(">=10SigJet");

  cutflowList.push_back("=0BJet");
  cutflowList.push_back("=1BJet");
  cutflowList.push_back("=2BJet");
  cutflowList.push_back("=3BJet");
  cutflowList.push_back("=4BJet");
  cutflowList.push_back(">=5BJet");

  cutflowList.push_back(">=2BaseLep");
  cutflowList.push_back(">=2BaseLep,w");
  cutflowList.push_back(">=2SigLep");
  cutflowList.push_back(">=2SigLep,w");
  cutflowList.push_back("SS leptons");
  cutflowList.push_back(">=1 passOR jet");
  cutflowList.push_back(">=1 signal jet");
  cutflowList.push_back("Z veto");

  m_hCutFlow = new TH1D("hCutFlow", "cut flow", cutflowList.size(), 0, cutflowList.size());
  for(unsigned int i=0;i<cutflowList.size();i++) m_hCutFlow->GetXaxis()->SetBinLabel(i+1,cutflowList[i].Data());
  m_hCutFlow->SetDirectory(outputFile);
  // m_hCutFlow->SetDirectory(0);
  // m_susyEvt->tree2->GetUserInfo()->Add(m_hCutFlow);

  m_hCutFlowNominal = m_hCutFlow;
  m_hCutFlowDummy = new TH1D( *(TH1D*)m_hCutFlow );
  m_hCutFlowDummy->SetDirectory(0); //prevent saving this hist, set it not belonging to current opend file/dir

  m_hTrigs = new TH1D("hTrigs", "n pass trigger", CF_trigNames_2015.size()+CF_trigNames_2016.size(), 0, CF_trigNames_2015.size()+CF_trigNames_2016.size());
  for(unsigned int i=0; i<CF_trigNames_2015.size(); i++){
    m_hTrigs->GetXaxis()->SetBinLabel(i+1,CF_trigNames_2015[i].c_str());
    std::cout<<CF_trigNames_2015[i].c_str()<<std::endl;
  }
  for(unsigned int i=0; i<CF_trigNames_2016.size(); i++){
    m_hTrigs->GetXaxis()->SetBinLabel(i+1+CF_trigNames_2015.size(),CF_trigNames_2016[i].c_str());
    std::cout<<CF_trigNames_2016[i].c_str()<<std::endl;
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
  // m_objTool->msg().setLevel( MSG::VERBOSE );

  for(auto x: CF_PRW_confFiles){Info("CF_PRW_confFiles", "%s", x.c_str());}

  ST::ISUSYObjDef_xAODTool::DataSource ds = static_cast<ST::ISUSYObjDef_xAODTool::DataSource>(CF_isMC); 

  //Get shower type
  int ShowerType = CF_isMC ? m_objTool->getMCShowerType(SampleName) : 0;

  CHECK(m_objTool->setProperty("DataSource",ds));
  CHECK(m_objTool->setProperty("ConfigFile", CF_ConfigFile));
  CHECK(m_objTool->setProperty("ShowerType", ShowerType));
  CHECK(m_objTool->setProperty("PRWConfigFiles", CF_PRW_confFiles));
  CHECK(m_objTool->setProperty("PRWLumiCalcFiles", CF_PRW_lcalcFiles));
  CHECK(m_objTool->initialize().isSuccess());

  mChargeFlipBkgTool = new ChargeFlipBkgTool("MyQFlipTool");
  CHECK(mChargeFlipBkgTool->setProperty("InputRatesFileName" , "$ROOTCOREBIN/data/multiLepSearch/root_files/chargeFlipRates.root"));
  CHECK(mChargeFlipBkgTool->setProperty("InputRatesHistoName", "hFlipProb_data"));
  CHECK(mChargeFlipBkgTool->initialize());

  mFakeLepBkgTool = new FakeLepBkgTool("MyFLepTool");
  CHECK(mFakeLepBkgTool->setProperty("Method", "Matrix"));
  CHECK(mFakeLepBkgTool->setProperty("InputFileName"    , "$ROOTCOREBIN/data/multiLepSearch/root_files/RealFakeLepEff_WhSS_Aug2017.root"));
  CHECK(mFakeLepBkgTool->setProperty("RealeEffHistoName", "RealeEff"));
  CHECK(mFakeLepBkgTool->setProperty("RealuEffHistoName", "RealuEff"));
  CHECK(mFakeLepBkgTool->setProperty("FakeeEffHistoName", "FakeeEff"));
  CHECK(mFakeLepBkgTool->setProperty("FakeuEffHistoName", "FakeuEff"));

  /*
  CHECK(mFakeLepBkgTool->setProperty("Method", "FakeFactor"));
  CHECK(mFakeLepBkgTool->setProperty("InputFileName"    , "$ROOTCOREBIN/data/multiLepSearch/root_files/fakefactor_2D_Data16.root"));
  CHECK(mFakeLepBkgTool->setProperty("eFakeFactorHistoName", "h_ff_ele"));
  //CHECK(mFakeLepBkgTool->setProperty("eFakeFactorHistoName", "h_ff_ele_v2")); //this histo has problem of bin error being just the sqrt of bin content
  CHECK(mFakeLepBkgTool->setProperty("uFakeFactorHistoName", "h_ff_mu"));
  */

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

  //Initialization of the Trigger SF Tool
  Info("initialize()", "Initializing Tool: \t %s", "TrigGlobalEfficiencyCorrectionTool");
  if( !this->InitializeTriggerTools() ){
    Error("initialize()", "Cannot initialize trigger SF Tools. Exiting." ); return EL::StatusCode::FAILURE; }
  TriggerSFTool = new TrigGlobalEfficiencyCorrectionTool("TrigGlobalEfficiencyCorrectionTool");
  ANA_CHECK( TriggerSFTool->setProperty("ElectronEfficiencyTools",electronEffTools) );
  ANA_CHECK( TriggerSFTool->setProperty("ElectronScaleFactorTools",electronSFTools) );
  ANA_CHECK( TriggerSFTool->setProperty("MuonTools",muonTools) );
  ANA_CHECK( TriggerSFTool->setProperty("ListOfLegsPerTool",LegsPerTool) );
  ANA_CHECK( TriggerSFTool->setProperty("TriggerCombination2015", "e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose || mu20_iloose_L1MU15_OR_mu40 || 2e12_lhloose_L12EM10VH || mu18_mu8noL1 || e17_lhloose_mu14 || e7_lhmedium_mu24") );
  ANA_CHECK( TriggerSFTool->setProperty("TriggerCombination2016","e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0 || mu26_imedium_OR_mu50 || 2e17_lhvloose_nod0 || e17_lhloose_nod0_mu14 || mu22_mu8noL1 || e7_lhmedium_nod0_mu24") );
  
  // ANA_CHECK( TriggerSFTool->setProperty("TriggerCombination2015", "2e12_lhloose_L12EM10VH || mu18_mu8noL1 || e17_lhloose_mu14") );
  // ANA_CHECK( TriggerSFTool->setProperty("TriggerCombination2016","2e17_lhvloose_nod0 || e17_lhloose_nod0_mu14 || mu22_mu8noL1") );

  ANA_CHECK( TriggerSFTool->initialize() );
  TriggerSFTool->msg().setLevel( MSG::VERBOSE );
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

    CHECK(m_objTool->ApplyPRWTool());
    bool passTrig=false;
    if (m_objTool->treatAsYear()==2015) for(auto trig : CF_trigNames_2015){
        if(m_objTool->IsTrigPassed(trig)) {passTrig=true; break;}
    }
    else if (m_objTool->treatAsYear()==2016) for(auto trig : CF_trigNames_2016){
        if(m_objTool->IsTrigPassed(trig)) {passTrig=true; break;}
    }
    if (passTrig) m_hCutFlow->Fill("Trigger", 1);
    else return sc;

    /// LArError, tileError, SCTError, CoreFlags
    bool isLArError = (eventInfo->errorState(xAOD::EventInfo::LAr)  == xAOD::EventInfo::Error ||
                       eventInfo->errorState(xAOD::EventInfo::Tile) == xAOD::EventInfo::Error ||
                       eventInfo->errorState(xAOD::EventInfo::SCT)  == xAOD::EventInfo::Error ||
                       eventInfo->isEventFlagBitSet(xAOD::EventInfo::Core, 18)                );

    if(!CF_isMC && isLArError) return sc;
    m_hCutFlow->Fill("LAr+Tile+SCT+CoreFlag", 1);

    // Primary vertex
    const xAOD::Vertex* PV = m_objTool->GetPrimVtx();
    int nTrk = PV ? PV->nTrackParticles() : 0;
    if( !nTrk ) return sc;
    m_hCutFlow->Fill("PV", 1);

    // Get the Electrons from the event
    xAOD::ElectronContainer* electrons_copy(0);
    xAOD::ShallowAuxContainer* electrons_copyaux(0);
    CHECK(m_objTool->GetElectrons(electrons_copy,electrons_copyaux));
    CHECK(wk()->xaodStore()->record(electrons_copy   , m_systInfoList[iSyst].systset.name()+"ShallowCopiedElectrons"));
    CHECK(wk()->xaodStore()->record(electrons_copyaux, m_systInfoList[iSyst].systset.name()+"ShallowCopiedElectronsAux."));

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

    xAOD::JetContainer* jets_copy_signal = new xAOD::JetContainer(SG::VIEW_ELEMENTS);
    CHECK(wk()->xaodStore()->record(jets_copy_signal        , m_systInfoList[iSyst].systset.name()+"ShallowCopiedJets_signal"));

    /// overlap removal
    CHECK( m_objTool->OverlapRemoval(electrons_copy, muons_copy, jets_copy) );

    // Bad muon and Cosmic muon
    bool hasBadMuon = false;
    bool isCosmicMuon = false;
    for(auto mu: *muons_copy){
      if(dec_baseline(*mu) && dec_bad(*mu)){
        hasBadMuon = true;
      }

      if(dec_baseline(*mu) && dec_passOR(*mu) && dec_cosmic(*mu)){
        isCosmicMuon = true;
      }
    }

    if(isCosmicMuon) continue;
    m_hCutFlow->Fill("cosMuon", 1);

    if(hasBadMuon) continue;
    m_hCutFlow->Fill("BadMuon", 1);

    /// Bad jet
    bool hasBadJet = false;
    bool hasBaselineJet = false;
    bool hasSignalJet = false;
    for(auto jet: *jets_copy){
      if(dec_baseline(*jet) && dec_passOR(*jet)){
        hasBaselineJet = true;
        if(dec_bad(*jet)) hasBadJet = true;
      }

      if(dec_signal(*jet)) hasSignalJet = true;
    }
    if(hasBadJet) continue;
    m_hCutFlow->Fill("BadJet", 1);

    //// count the leptons
    vector< IParticle* > sel_Ls;
    vector< IParticle* > sig_Ls;
    sel_Ls.reserve(3);
    sig_Ls.reserve(3);
    vector< IParticle* > base_el;
    vector< IParticle* > base_mu;

    int nBaseEl = 0;
    int nSigEl = 0;
    // Electrons

    trigElec.clear(); trigMuon.clear();

    for(auto el: *electrons_copy){
      //if(el->pt()<CF_ElPtCut) continue;
      if(dec_passOR(*el)){
        if     (dec_signal  (*el)){
           sig_Ls.push_back(el);nBaseEl++;nSigEl++;
           base_el.push_back(el);
           trigElec.push_back(el);
        }
        else if(dec_baseline(*el)){ sel_Ls.push_back(el);base_el.push_back(el);nBaseEl++;}

      }
      // trigElec.push_back(el);
    }
   
    int nBaseMu = 0;
    int nSigMu = 0;
    // Muons
    for(auto mu: *muons_copy){
      //if(mu->pt()<CF_MuPtCut) continue;
      if(dec_passOR(*mu)){
        if     (dec_signal  (*mu)){ 
          sig_Ls.push_back(mu);nBaseMu++;nSigMu++;
          base_mu.push_back(mu);
          trigMuon.push_back(mu);
        }
        else if(dec_baseline(*mu)){ sel_Ls.push_back(mu);base_mu.push_back(mu);nBaseMu++;}
      }
      // trigMuon.push_back(mu);
    }

    int& nBaseLep = m_susyEvt->sig.nBaseLep;
    int& nSigLep = m_susyEvt->sig.nSigLep;
    nBaseLep = sel_Ls.size() + sig_Ls.size();
    nSigLep  = sig_Ls.size();

    //cutflow for number of leptons
    {
      //Baseline leptons
      if     (nBaseLep==0) m_hCutFlow->Fill("=0BaseLep", 1);
      else if(nBaseLep==1) m_hCutFlow->Fill("=1BaseLep", 1);
      else if(nBaseLep==2) m_hCutFlow->Fill("=2BaseLep", 1);
      else if(nBaseLep==3) m_hCutFlow->Fill("=3BaseLep", 1);
      else if(nBaseLep==4) m_hCutFlow->Fill("=4BaseLep", 1);
      else                 m_hCutFlow->Fill(">=5BaseLep", 1);

      if     (nBaseEl==0) m_hCutFlow->Fill("=0BaseEl", 1);
      else if(nBaseEl==1) m_hCutFlow->Fill("=1BaseEl", 1);
      else if(nBaseEl==2) m_hCutFlow->Fill("=2BaseEl", 1);
      else if(nBaseEl==3) m_hCutFlow->Fill("=3BaseEl", 1);
      else if(nBaseEl==4) m_hCutFlow->Fill("=4BaseEl", 1);
      else                m_hCutFlow->Fill(">=5BaseEl", 1);

      if     (nBaseMu==0) m_hCutFlow->Fill("=0BaseMu", 1);
      else if(nBaseMu==1) m_hCutFlow->Fill("=1BaseMu", 1);
      else if(nBaseMu==2) m_hCutFlow->Fill("=2BaseMu", 1);
      else if(nBaseMu==3) m_hCutFlow->Fill("=3BaseMu", 1);
      else if(nBaseMu==4) m_hCutFlow->Fill("=4BaseMu", 1);
      else                m_hCutFlow->Fill(">=5BaseMu", 1);

      //Signal leptons
      if     (nSigLep==0) m_hCutFlow->Fill("=0SigLep", 1);
      else if(nSigLep==1) m_hCutFlow->Fill("=1SigLep", 1);
      else if(nSigLep==2) m_hCutFlow->Fill("=2SigLep", 1);
      else if(nSigLep==3) m_hCutFlow->Fill("=3SigLep", 1);
      else if(nSigLep==4) m_hCutFlow->Fill("=4SigLep", 1);
      else                m_hCutFlow->Fill(">=5SigLep", 1);

      if     (nSigEl==0) m_hCutFlow->Fill("=0SigEl", 1);
      else if(nSigEl==1) m_hCutFlow->Fill("=1SigEl", 1);
      else if(nSigEl==2) m_hCutFlow->Fill("=2SigEl", 1);
      else if(nSigEl==3) m_hCutFlow->Fill("=3SigEl", 1);
      else if(nSigEl==4) m_hCutFlow->Fill("=4SigEl", 1);
      else               m_hCutFlow->Fill(">=5SigEl", 1);

      if     (nSigMu==0) m_hCutFlow->Fill("=0SigMu", 1);
      else if(nSigMu==1) m_hCutFlow->Fill("=1SigMu", 1);
      else if(nSigMu==2) m_hCutFlow->Fill("=2SigMu", 1);
      else if(nSigMu==3) m_hCutFlow->Fill("=3SigMu", 1);
      else if(nSigMu==4) m_hCutFlow->Fill("=4SigMu", 1);
      else               m_hCutFlow->Fill(">=5SigMu", 1);

      if(nBaseLep == 2 && nSigLep == 2) m_hCutFlow->Fill("=2BaseLep and =2SigLep", 1);
    }

    // jets
    vector< xAOD::IParticle* > jet_Ls;
    jet_Ls.reserve(10);
    int& nBJet = m_susyEvt->sig.nBJet;
    nBJet = 0;
    for(auto jet: *jets_copy){
      if(!dec_passOR(*jet)) continue;
      if(!dec_signal(*jet)) continue;

      float jvt = 0;
      if(jet->isAvailable<float>("Jvt")) {
        //don't know the reason why accessing jvt will change the information of jet, and MET will be changed. use a new copy.
        xAOD::Jet jet_copy = (*jet);
        jvt = jet_copy.auxdata<float>("Jvt");
      }
      if(jet->pt()<60000 && fabs(jet->eta()) < 2.4 && jvt < 0.59) continue;
      bool isbjet = m_objTool->IsBJet(*jet);
      if(study != "fakes_Peter" && isbjet && fabs(jet->eta()) >= 2.4) continue;

      jet_Ls.push_back(jet);
      jets_copy_signal->push_back(jet);
      if(isbjet) nBJet++;
    }
    int& nSigJet = m_susyEvt->sig.nJet;
    nSigJet = jet_Ls.size();

    //Cutflow for signal jets
    if     (nSigJet==0) m_hCutFlow->Fill("=0SigJet", 1);
    else if(nSigJet==1) m_hCutFlow->Fill("=1SigJet", 1);
    else if(nSigJet==2) m_hCutFlow->Fill("=2SigJet", 1);
    else if(nSigJet==3) m_hCutFlow->Fill("=3SigJet", 1);
    else if(nSigJet==4) m_hCutFlow->Fill("=4SigJet", 1);
    else if(nSigJet==5) m_hCutFlow->Fill("=5SigJet", 1);
    else if(nSigJet==6) m_hCutFlow->Fill("=6SigJet", 1);
    else if(nSigJet==7) m_hCutFlow->Fill("=7SigJet", 1);
    else if(nSigJet==8) m_hCutFlow->Fill("=8SigJet", 1);
    else if(nSigJet==9) m_hCutFlow->Fill("=9SigJet", 1);
    else                m_hCutFlow->Fill(">=10SigJet", 1);

    //Cutflow for b-jets
    if     (nBJet==0) m_hCutFlow->Fill("=0BJet", 1);
    else if(nBJet==1) m_hCutFlow->Fill("=1BJet", 1);
    else if(nBJet==2) m_hCutFlow->Fill("=2BJet", 1);
    else if(nBJet==3) m_hCutFlow->Fill("=3BJet", 1);
    else if(nBJet==4) m_hCutFlow->Fill("=4BJet", 1);
    else              m_hCutFlow->Fill(">=5BJet", 1);

    //at least 2 baseline leptons
    if(nBaseLep < 2) continue;
    m_hCutFlow->Fill(">=2BaseLep", 1);

    /// sort leptons
    sort(sel_Ls.begin(), sel_Ls.end(), [](xAOD::IParticle* a, xAOD::IParticle* b)->bool{return a->pt()>b->pt();});
    sort(sig_Ls.begin(), sig_Ls.end(), [](xAOD::IParticle* a, xAOD::IParticle* b)->bool{return a->pt()>b->pt();});
    sel_Ls.insert( sel_Ls.begin(), sig_Ls.begin(), sig_Ls.end());

    /// sort the two leptons
    if(sel_Ls[0]->pt() < sel_Ls[1]->pt()) {
      IParticle* temp = sel_Ls[0];
      sel_Ls[0] = sel_Ls[1];
      sel_Ls[1] = temp;
    }

    //The two leptons
    xAOD::IParticle* l1 = sel_Ls[0];
    xAOD::IParticle* l2 = sel_Ls[1];

    /////////////////////////
    // Save informations
    ////////////////////////
    
    //unused variables
    m_susyEvt->evt.flag = 0;

    //event information
    m_susyEvt->evt.event = eventInfo->eventNumber();
    //m_susyEvt->evt.lumiBlock = eventInfo->lumiBlock();
    //m_susyEvt->evt.actualMu = eventInfo->actualInteractionsPerCrossing();
    
    double TotalWeight = 1;
    if(CF_isMC) TotalWeight *= eventInfo->mcEventWeight();
    if(!CF_isMC) m_susyEvt->evt.weight = 1;
    else m_susyEvt->evt.weight = eventInfo->mcEventWeight();
    
    m_susyEvt->evt.isMC = CF_isMC? 1:0;

    m_susyEvt->truths.clear();

    // charge flip weight and pT correction (have to apply before other calculation that use lep pT)
    m_susyEvt->evt.qFwt = 0.0;
    //m_susyEvt->evt.qFwt_sys_1up = 0.0;
    //m_susyEvt->evt.qFwt_sys_1dn = 0.0;
    if (nBaseLep == 2 && nSigLep == 2){
      if (!isSS(sel_Ls[0], sel_Ls[1])){
        m_susyEvt->evt.qFwt = mChargeFlipBkgTool->GetWeight( sel_Ls ,0,0);
        //m_susyEvt->evt.qFwt_sys_1up = mChargeFlipBkgTool->GetWeight( sel_Ls , 1,0);
        //m_susyEvt->evt.qFwt_sys_1dn = mChargeFlipBkgTool->GetWeight( sel_Ls ,-1,0);
        // auto tmpPt = mChargeFlipBkgTool->GetCorrectedPt( sel_Ls ,0,0);
        // if(tmpPt.size()==2){
        //   if (tmpE0) tmpE0->setP4( tmpPt[0]*1000., tmpE0->eta(), tmpE0->phi(), tmpE0->m());
        //   if (tmpE1) tmpE1->setP4( tmpPt[1]*1000., tmpE1->eta(), tmpE1->phi(), tmpE1->m());
          //ATH_MSG_ERROR("E0" << tmpPt[0] << " " << sel_Ls[0]->pt());
          //ATH_MSG_ERROR("E1" << tmpPt[1] << " " << sel_Ls[1]->pt());
        // }
      }
    }
    //fake lep weight
    m_susyEvt->evt.fLwt = 0.0;
    //m_susyEvt->evt.fLwt_e_sys_1up = 0.0; 
    //m_susyEvt->evt.fLwt_e_sys_1dn = 0.0; 
    //m_susyEvt->evt.fLwt_u_sys_1up = 0.0; 
    //m_susyEvt->evt.fLwt_u_sys_1dn = 0.0; 
    //if (nBaseLep == 2){
      //m_susyEvt->evt.fLwt = mFakeLepBkgTool->GetWeight(sel_Ls, 0,0);
      //m_susyEvt->evt.fLwt_e_sys_1up = mFakeLepBkgTool->GetWeight(sel_Ls, 1,0);
      //m_susyEvt->evt.fLwt_e_sys_1dn = mFakeLepBkgTool->GetWeight(sel_Ls,-1,0);
      //m_susyEvt->evt.fLwt_u_sys_1up = mFakeLepBkgTool->GetWeight(sel_Ls, 1,1);
      //m_susyEvt->evt.fLwt_u_sys_1dn = mFakeLepBkgTool->GetWeight(sel_Ls,-1,1);
      //ATH_MSG_ERROR("FW " << mFakeLepBkgTool->GetWeight(sel_Ls, 0,0));
    //}

    /// photon Check
    xAOD::PhotonContainer* photons_copy(0);
    xAOD::ShallowAuxContainer* photons_copyaux(0);
    CHECK( m_objTool->GetPhotons(photons_copy,photons_copyaux) );
    CHECK(wk()->xaodStore()->record(photons_copy       , m_systInfoList[iSyst].systset.name()+"ShallowCopiedPhotons"));
    CHECK(wk()->xaodStore()->record(photons_copyaux    , m_systInfoList[iSyst].systset.name()+"ShallowCopiedPhotonsAux."));

    /// met
    xAOD::MissingETContainer* metcst = new xAOD::MissingETContainer();
    xAOD::MissingETAuxContainer* metcst_aux = new xAOD::MissingETAuxContainer();
    metcst->setStore(metcst_aux);
    CHECK( m_objTool->GetMET(*metcst,jets_copy,electrons_copy,muons_copy,photons_copy) );

    auto metFinal = (*metcst)["Final"];
    m_susyEvt->sig.Met  = metFinal->met()*iGeV;
    m_susyEvt->sig.MetX = metFinal->mpx()*iGeV;
    m_susyEvt->sig.MetY = metFinal->mpy()*iGeV;
    m_susyEvt->sig.MetPhi = metFinal->phi();
    //TLorentzVector metV(metFinal->mpx(), metFinal->mpy(), 0, metFinal->met());

    //sort jets
    sort(jet_Ls.begin(), jet_Ls.end(), [](xAOD::IParticle* a, xAOD::IParticle* b)->bool{return a->pt()>b->pt();});

    m_susyEvt->sig.HT = 0;
    m_susyEvt->jets.resize(jet_Ls.size());
    for(unsigned int i=0;i<jet_Ls.size(); i++){
      auto j = dynamic_cast<xAOD::Jet*>(jet_Ls[i]);
      
      m_susyEvt->jets[i].pt = j->pt()*iGeV; 
      m_susyEvt->jets[i].eta = j->eta(); 
      m_susyEvt->jets[i].phi = j->phi();
      m_susyEvt->jets[i].m = j->m()*iGeV;
      unsigned int& flag = m_susyEvt->jets[i].jFlag;
      flag = 0;
      if(dec_signal(*j)) flag |= IS_SIGNAL;
      //if(dec_bjet_loose(*j)) flag |= JT_BJET_LOOSE;
      if(m_objTool->IsBJet(*j)) flag |= JT_BJET;

      m_susyEvt->sig.HT += j->pt()*iGeV;
    }

    // Save Z decay chain
    // Info("execute()", "Save Z decay chain");
    if(CF_isMC && mcTruthMatch.back()=='Z')
    {
      const xAOD::TruthParticle* truthZ=0;
      const xAOD::TruthParticleContainer* tContainer = 0;
      CHECK(wk()->xaodEvent()->retrieve( tContainer, "TruthParticles" ));
      for(auto particle : *tContainer){
        if(!particle->isZ()) continue;
        truthZ = particle;
        while(truthZ->parent() && truthZ->parent()->isZ()){truthZ = truthZ->parent();}
        break;
      }

      addTruthPar(truthZ, m_susyEvt->truths, 1);
    }
    // Info("execute()", "Z decay saved");

    if(jet_Ls.size()>=2)  m_susyEvt->sig.mjj = (jet_Ls[0]->p4()+jet_Ls[1]->p4()).M() *iGeV;
    else m_susyEvt->sig.mjj = -1;

    if(jet_Ls.size() == 0) m_susyEvt->sig.mlj = -1;
    else
    {
      TLorentzVector JetSystem;
      if(jet_Ls.size()==1)
      {
        JetSystem = jet_Ls[0]->p4();
      }
      else
      {
        JetSystem = jet_Ls[0]->p4()+jet_Ls[1]->p4();
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
    //trigElec.clear(); trigMuon.clear(); 
    // Info("execute()", "There are %lu leptons in this event", sel_Ls.size());
    for(unsigned int i=0;i<sel_Ls.size(); i++){
      auto& l = m_susyEvt->leps[i];
      fillLepton(sel_Ls[i], l, i);
      //l.MET_dPhi = metV.DeltaPhi(sel_Ls[i]->p4());

      {
        /// mT
        TVector3 lepVec;
        lepVec.SetPtEtaPhi(l.pt, 0, l.phi);
        TLorentzVector Vt;
        Vt.SetPxPyPzE( m_susyEvt->sig.Met * TMath::Cos(m_susyEvt->sig.MetPhi) + lepVec.Px(),
                       m_susyEvt->sig.Met * TMath::Sin(m_susyEvt->sig.MetPhi) + lepVec.Py(),
                       0,
                       m_susyEvt->sig.Met + lepVec.Pt() );
        l.mT = Vt.M();
        //l.mT = sqrt(2*l.pt*m_susyEvt->sig.Met*(1-cos(l.MET_dPhi)));
      }
    }
    // Info("execute()", "After fillLepton");

    /*
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
    */
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
    unsigned long int ADD;
    // There are 8 triggers for each year. 
    // Use the first byte to store 2015
    // Use the second byte to store 2016
    if (m_objTool->treatAsYear()==2015){ 
      ADD = 1; 
      for(auto trig : CF_trigNames_2015){
        if(m_objTool->IsTrigPassed(trig)){
          m_susyEvt->sig.trigCode += ADD;
          m_hTrigs->Fill(trig.c_str(), 1);
        }
        ADD *= 2;
      }
    }
    else if (m_objTool->treatAsYear()==2016){ 
      ADD = 1 << CF_trigNames_2015.size(); 
      for(auto trig : CF_trigNames_2016){
        if(m_objTool->IsTrigPassed(trig)){
          m_susyEvt->sig.trigCode += ADD;
          m_hTrigs->Fill(trig.c_str(), 1);
        }
        ADD *= 2;
      }
    }

    //l12
    TLorentzVector ll(l1->p4()+l2->p4());
    m_susyEvt->l12.m = ll.M()*iGeV;
    m_susyEvt->l12.pt = ll.Pt()*iGeV; 
    m_susyEvt->l12.eta = ll.Eta(); 
    m_susyEvt->l12.phi = ll.Phi(); 
    m_susyEvt->l12.dPhi = l1->p4().DeltaPhi(l2->p4()); 
    m_susyEvt->l12.dR = l1->p4().DeltaR(l2->p4()); 
    if(jet_Ls.size()>=1) m_susyEvt->l12.jet0_dPhi = ll.DeltaPhi(jet_Ls[0]->p4());
    else m_susyEvt->l12.jet0_dPhi = -100;

    m_susyEvt->sig.HT += (l1->pt()+l2->pt())*iGeV;

    /// mT2
    auto tl1 = l1->p4();
    auto tl2 = l2->p4();
    //https://svnweb.cern.ch/trac/atlasphys-susy/browser/Physics/SUSY/Tools/CalcGenericMT2
    m_susyEvt->sig.mT2 = anaHelper::get_mT2(tl1.M(), tl1.Px(), tl1.Py(),
                                            tl2.M(), tl2.Px(), tl2.Py(),
                                            //metFinal->mpx(),
                                            //metFinal->mpy(),
                                            metFinal->met() * TMath::Cos(metFinal->phi()),
                                            metFinal->met() * TMath::Sin(metFinal->phi()),
                                            //CF_mT2_m0*GeV, CF_mT2_m0*GeV) *iGeV;
                                            tl1.M(), tl2.M()) *iGeV;

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
        TotalWeight *= m_susyEvt->evt.pwt;
      }else{
        m_susyEvt->evt.pwt = 1;
        m_susyEvt->evt.averageMu = m_objTool->GetCorrectedAverageInteractionsPerCrossing();
      }

      //// get trigger info, if passed, and the SF
      auto& sEvt = m_susyEvt->evt;
      //sEvt.run = eventInfo->runNumber();
      sEvt.run = m_objTool->GetRunNumber();

      //Scale factor
      if(CF_isMC){
        // Info("execute()", "Calculate scale factors");
        // Trigger SF from SUSYTools turned off as of Aug 3 2017. Replaced by TrigGlobalEfficiencyTool 11 Aug 2017.
        sEvt.ElSF = m_objTool->GetTotalElectronSF(*electrons_copy, true, true, false /*trigSF*/, true, "", true);
        sEvt.MuSF = m_objTool->GetTotalMuonSF(*muons_copy, true, true, "");
        sEvt.BtagSF = m_objTool->BtagSF(jets_copy_signal);
        sEvt.JvtSF = m_objTool->JVT_SF(jets_copy_signal);

        double trigSF=1;
        auto trigRet = TriggerSFTool->getEfficiencyScaleFactor(sEvt.run,trigElec, trigMuon, trigSF); //sEvt.run

        sEvt.trigSF = trigSF;

        if (eventInfo->eventNumber()==4576526 || eventInfo->eventNumber()==5093803 ||
	    eventInfo->eventNumber()==5095103 || eventInfo->eventNumber()==4660617 ||
		eventInfo->eventNumber()==5094844 || eventInfo->eventNumber()==4660570 ||
		eventInfo->eventNumber()==5096519){

          std::cout << "################################################### DEBUGGING #####################################" << std::endl;
          std::cout << "-- Trigger SF for event : " << eventInfo->eventNumber() << ",   trigger SF :" << trigSF << std::endl;
          std::cout << "---- Run : " <<  sEvt.run << std::endl;
          //std::cout << "--- Trigger Electron : " << trigElec << std::endl;
          //std::cout << "--- Trigger Muon: " << trigMuon << std::endl;
        }

        // DEBUG MESSAGES. COMMENT OUT TO READ
        // if (trigRet==CP::CorrectionCode::OutOfValidityRange){
        //   cout << "Event " << m_eventCounter << ": Year " << m_objTool->treatAsYear() << ", trigSF = " << sEvt.trigSF << endl;
        //   cout << sel_Ls.size() << " selected leptons" << endl;
        //   for (uint i=0; i<sel_Ls.size(); i++)
        //   {
        //     cout << "Lepton " << i;
        //     xAOD::Muon* mu = dynamic_cast<xAOD::Muon*>(sel_Ls[i]);
        //     xAOD::Electron* el = dynamic_cast<xAOD::Electron*>(sel_Ls[i]);
        //     if(mu) cout << " muon" ;
        //     else if (el) cout << " electron";
        //     else cout << " unknown flavor";

        //     cout << endl << "PT: " << sel_Ls[i]->pt() << "\t" << "ETA: " << sel_Ls[i]->eta() << endl;      
        //   }
        //   cout << endl;
        // }
        
        TotalWeight *= sEvt.ElSF*sEvt.MuSF*sEvt.BtagSF*sEvt.JvtSF*sEvt.trigSF;

      }else{
        sEvt.ElSF = 1;
        sEvt.MuSF = 1;
        sEvt.BtagSF = 1;
        sEvt.JvtSF = 1;
      }

      if(nBaseLep >= 2) m_hCutFlow->Fill(">=2BaseLep,w", TotalWeight);
      if(nSigLep >= 2) m_hCutFlow->Fill(">=2SigLep,w", TotalWeight);
      if(nBaseLep == 2 && nSigLep == 2) m_hCutFlow->Fill("=2BaseLep and =2SigLep,w", TotalWeight);

      m_susyEvt->sig.isSS = 0;
      m_susyEvt->sig.JetCut = (hasBaselineJet && hasSignalJet);
      m_susyEvt->sig.isZ = isZ(base_el) || isZ(base_mu);
      //cutflow for Dani ntuple
      if(nSigLep >= 2)
      {
        m_hCutFlow->Fill(">=2SigLep", 1);
        if(isSS(sig_Ls[0], sig_Ls[1]) || nSigLep >= 3)
        {
          m_susyEvt->sig.isSS = 1;
          m_hCutFlow->Fill("SS leptons", 1);
          if(hasBaselineJet)
          {
            m_hCutFlow->Fill(">=1 passOR jet", 1);
            if(hasSignalJet)
            {
              m_hCutFlow->Fill(">=1 signal jet", 1);
              if(!m_susyEvt->sig.isZ) m_hCutFlow->Fill("Z veto", 1);
            }
          }
        }
      }

      if(CF_isMC && printEvent>=0 && m_susyEvt->evt.event == (uint) printEvent)
      {
        cout<<"event: "<<m_susyEvt->evt.event<<endl;
        cout<<"mcweight: "<<eventInfo->mcEventWeight()<<endl;
        cout<<"pileup: "<<m_susyEvt->evt.pwt<<endl;
        cout<<"ele: "<<sEvt.ElSF<<endl;
        cout<<"mu: "<<sEvt.MuSF<<endl;
        cout<<"bjet: "<<sEvt.BtagSF<<endl;
        cout<<"jvt: "<<sEvt.JvtSF<<endl;
        cout<<"trigSF: "<<sEvt.trigSF<<endl;
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

  //Print Cutflow
  std::cout << " Cutflow summary: " << std::endl;
  std::cout << "============================" << std::endl;
  for(unsigned int i=1;i<=cutflowList.size();i++)
  {
    TString LabelName = m_hCutFlow->GetXaxis()->GetBinLabel(i);
    std::cout<<LabelName<<": ";
    if(LabelName.Contains(",w")) std::cout<< m_hCutFlow->GetBinContent(i) << std::endl;
    else std::cout<< int(m_hCutFlow->GetBinContent(i)) << std::endl;
  }
  std::cout << "============================" << std::endl; 

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

  // Info("fillLepton(el)", "Before trigger matching");
  // Save single lepton trigger matching info for fakes
  if (study=="fakes")
  {
    if (m_objTool->treatAsYear() == 2015) // 2015 data
      {
        if (m_objTool->IsTrigPassed("HLT_e24_lhmedium_L1EM20VH")
          || m_objTool->IsTrigPassed("HLT_e60_lhmedium") 
          || m_objTool->IsTrigPassed("HLT_e120_lhloose") ) l.ID +=1;
        if (m_objTool->IsTrigMatched(el, "HLT_e24_lhmedium_L1EM20VH")
          || m_objTool->IsTrigMatched(el, "HLT_e60_lhmedium") 
          || m_objTool->IsTrigMatched(el, "HLT_e120_lhloose") ) l.ID +=2;
      }
      else // 2016 data
      {
        if (m_objTool->IsTrigPassed("HLT_e26_lhtight_nod0_ivarloose")
          || m_objTool->IsTrigPassed("HLT_e60_lhmedium_nod0") 
          || m_objTool->IsTrigPassed("HLT_e140_lhloose_nod0") ) l.ID +=1;
        if (m_objTool->IsTrigMatched(el, "HLT_e26_lhtight_nod0_ivarloose")
          || m_objTool->IsTrigMatched(el, "HLT_e60_lhmedium_nod0") 
          || m_objTool->IsTrigMatched(el, "HLT_e140_lhloose_nod0") ) l.ID +=2;
      }
  }
  // Info("fillLepton(el)", "After trigger matching");

 // trigElec.push_back(el);

  l.ID *= el->charge();
  //l.author = el->author();

  //isolation
  //if(!el->isolationValue(l.topoetcone20, xAOD::Iso::topoetcone20)) Error("fillLepton(el)", "topoetcone20 failed."); 
  //if(!el->isolationValue(l.topoetcone30, xAOD::Iso::topoetcone30)) Error("fillLepton(el)", "topoetcone30 failed."); 
  //if(!el->isolationValue(l.topoetcone40, xAOD::Iso::topoetcone40)) Error("fillLepton(el)", "topoetcone40 failed."); 
  //if(!el->isolationValue(l.ptcone20, xAOD::Iso::ptcone20)) Error("fillLepton(el)", "ptcone20 failed."); 
  //if(!el->isolationValue(l.ptcone30, xAOD::Iso::ptcone30)) Error("fillLepton(el)", "ptcone30 failed."); 
  //if(!el->isolationValue(l.ptcone40, xAOD::Iso::ptcone40)) Error("fillLepton(el)", "ptcone40 failed."); 
  //l.isoPass = m_isoTool->accept(*el).getCutResultBitSet().to_ulong();

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

  // Info("fillLepton(el)", "Before CF_isMC");
  if(CF_isMC)
  {
    std::pair<MCTruthPartClassifier::ParticleType, MCTruthPartClassifier::ParticleOrigin> res;
    // Info("fillLepton(el)", "About to run particleTruthClassifier()");
    res = m_truthClassifier->particleTruthClassifier(el);
    // Info("fillLepton(el)", "After particleTruthClassifier()");

    //l.truthType = res.first;
    //l.truthOrig = res.second;
    l.truthType = acc_truthType(*el);
    l.truthOrig = acc_truthOrig(*el);
    l.firstEgMotherPdgId = el->auxdata<int>("firstEgMotherPdgId");
    l.lepTruth = (l.truthType==2) || //IsoElectron
                 (l.truthOrig==5 && abs(l.firstEgMotherPdgId)==11 && el->charge() * l.firstEgMotherPdgId < 0); //PhotonConv 

    const xAOD::TruthParticle *tp = 0;
    if (mcTruthMatch == "TruthLink" || (mcTruthMatch=="tryAll" && !tp)) {
      auto tl = acc_truthLink(*el);
      if (tl.isValid()) tp = *tl;
    }

    else if (mcTruthMatch == "MCTC" || (mcTruthMatch=="tryAll" && !tp)){
      tp = m_truthClassifier->getGenPart();
    }

    else if (mcTruthMatch == "reverseTruthLink" || (mcTruthMatch=="tryAll" && !tp)){
      const xAOD::TruthParticleContainer* tContainer = 0;
      CHECK(wk()->xaodEvent()->retrieve( tContainer, "TruthParticles" ));
      for(auto particle : *tContainer){
        if(!(particle->isElectron())) continue;
        auto matchedPart = acc_recoElLink(*particle);
        if(matchedPart.isValid() && const_cast<xAOD::Electron*>(*matchedPart)==el) {
          tp = particle; break;
        }
      }
    }

    else if (mcTruthMatch == "dR" || (mcTruthMatch=="tryAll" && !tp)){
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
          tp = particle;
          p_tp = p_particle;
          firstTry = false;
        }
        else continue;
      }
    }
    // Info("fillLepton(el)", "After trying to find a match");

    /// save truth match and parents if exist, otherwise save -1.
    if (tp!=0) {
      l.truthI = addTruthPar(tp, m_susyEvt->truths, -1);
      m_susyEvt->truths[l.truthI].matchI = index;
    } else l.truthI = -1;
    // Info("fillLepton(el)", "truthI = %d", l.truthI);

  }
  fillLeptonCommon(el, l);
  return EL::StatusCode::SUCCESS;
}

EL::StatusCode ssEvtSelection :: fillLepton(xAOD::Muon* mu, L_PAR& l, unsigned int index)
{
  l.ID = 13000;

  // Save trigger information for fakes study
  if (study=="fakes")
  {
    if (m_objTool->treatAsYear() == 2015)
    {
      if (m_objTool->IsTrigPassed("HLT_mu20_iloose_L1MU15")
        || m_objTool->IsTrigPassed("HLT_mu40")) l.ID +=1;
      if (m_objTool->IsTrigMatched(mu, "HLT_mu20_iloose_L1MU15")
        || m_objTool->IsTrigMatched(mu, "HLT_mu40")) l.ID +=2;
    }
    else
    {
      if (m_objTool->IsTrigPassed("HLT_mu25_imedium")
        || m_objTool->IsTrigPassed("HLT_mu50")) l.ID +=1;
      if (m_objTool->IsTrigMatched(mu, "HLT_mu25_imedium")
        || m_objTool->IsTrigMatched(mu, "HLT_mu50")) l.ID +=2;
    }
  }
  
  //trigMuon.push_back(mu);

  l.ID *= mu->charge();
  //l.author = mu->author();

  //isolation
  //if(!mu->isolation(l.topoetcone20, xAOD::Iso::topoetcone20)) Error("fillLepton(mu)", "topoetcone20 failed."); 
  //if(!mu->isolation(l.topoetcone30, xAOD::Iso::topoetcone30)) Error("fillLepton(mu)", "topoetcone30 failed."); 
  //if(!mu->isolation(l.topoetcone40, xAOD::Iso::topoetcone40)) Error("fillLepton(mu)", "topoetcone40 failed."); 
  //if(!mu->isolation(l.ptcone20, xAOD::Iso::ptcone20)) Error("fillLepton(mu)", "ptcone20 failed."); 
  //if(!mu->isolation(l.ptcone30, xAOD::Iso::ptcone30)) Error("fillLepton(mu)", "ptcone30 failed."); 
  //if(!mu->isolation(l.ptcone40, xAOD::Iso::ptcone40)) Error("fillLepton(mu)", "ptcone40 failed."); 
  //l.isoPass = m_isoTool->accept(*mu).getCutResultBitSet().to_ulong();

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
    // Info("fillLepton(mu)", "Before MCTruthClassifier");
    std::pair<MCTruthPartClassifier::ParticleType, MCTruthPartClassifier::ParticleOrigin> res;
    res = m_truthClassifier->particleTruthClassifier(mu);
    l.firstEgMotherPdgId = 0;

    //l.truthType = res.first;
    //l.truthOrig = res.second;
    l.truthType = acc_truthType(*mu);
    l.truthOrig = acc_truthOrig(*mu);
    l.lepTruth = (l.truthType==6); //IsoMuon

    const xAOD::TruthParticle *tp = 0;
    if(mcTruthMatch=="TruthLink" || (mcTruthMatch=="tryAll" && !tp))
    {
      auto tl = acc_truthLink(*mu);
      if (tl.isValid()) tp = *tl;
    }

    // MCTruthClassifier
    else if (mcTruthMatch=="MCTC" || (mcTruthMatch=="tryAll" && !tp)){
      tp = m_truthClassifier->getGenPart();
    }

    // Reverse truthLink
    else if (mcTruthMatch=="reverseTruthLink" || (mcTruthMatch=="tryAll" && !tp)){
      const xAOD::TruthParticleContainer* tContainer = 0;
      CHECK(wk()->xaodEvent()->retrieve( tContainer, "TruthParticles" ));
      for(auto particle : *tContainer){
        if(!(particle->isMuon())) continue;
        auto matchedPart = acc_recoMuLink(*particle);
        if(matchedPart.isValid() && (const_cast<xAOD::Muon*>(*matchedPart))==mu) {
          tp = particle; break;
        }
      }
    }

    // dR
    else if (mcTruthMatch=="dR" || (mcTruthMatch=="tryAll" && !tp)){
      // Set up truth particle container
      const xAOD::TruthParticleContainer* tContainer = 0;
      CHECK(wk()->xaodEvent()->retrieve( tContainer, "TruthParticles" ));
   
      // Momentum vectors
      TLorentzVector p_tp(0,0,0,0);
      TLorentzVector p_m = mu->p4();
      bool firstTry = true;
   
      // Find electron with smallest dR < 0.1
      for(auto particle : *tContainer){
        if(!(particle->isElectron())) continue;
   
        TLorentzVector p_particle = particle->p4();
        Double_t dR = p_m.DeltaR(p_particle);
        if (dR > 0.1) continue;
         
        if (firstTry || dR < p_m.DeltaR(p_tp)){
          tp = particle;
          p_tp = p_particle;
          firstTry = false;
        }
        else continue;
      }
    }
    
    /// save truth match and parents if exist, otherwise save -1.
    if (tp!=NULL) {
      l.truthI = addTruthPar(tp, m_susyEvt->truths, -1);
      m_susyEvt->truths[l.truthI].matchI = index;
    } else l.truthI = -1;
    // Info("fillLepton(mu)", "l.truthI = %d", l.truthI);

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
  l.truthType = 0;
  l.truthOrig = 0;
  l.firstEgMotherPdgId = 0;
  // Info("fillLepton(p)", "Casting leptons");
  xAOD::Muon* mu = dynamic_cast<xAOD::Muon*>(p);
  if(mu) {
    fillLepton(mu, l, index);
    if(CF_isMC && mcTruthMatch=="MCTCZ") m_truthClassifier->particleTruthClassifier(mu);
  }
  else{
    xAOD::Electron* el = dynamic_cast<xAOD::Electron*>(p);
    if(el) fillLepton(el, l, index);
    if(CF_isMC && mcTruthMatch=="MCTCZ") m_truthClassifier->particleTruthClassifier(el);
  }

  // Info("fillLepton(p)", "Before Z fill");
  if(CF_isMC && mcTruthMatch.back()=='Z')
  {
    const xAOD::TruthParticle* tp=0;
    if(mcTruthMatch=="dRZ")
    {
      // Info("fillLepton(p)", "MCTruthMatch is dRZ");
      TLorentzVector p_tp(0,0,0,0);
      TLorentzVector p_lep = p->p4();
      bool firstTry = true;
    
      // Find particle with smallest dR < 0.1
      const xAOD::TruthParticleContainer* tContainer = 0;
      CHECK(wk()->xaodEvent()->retrieve( tContainer, "TruthParticles" ));
      for(auto particle : *tContainer){ 
        if (particle->p4().Pt()==0) continue;
        TLorentzVector p_particle = particle->p4();
        Double_t dR = p_lep.DeltaR(p_particle);
        if (dR > 0.1) continue;
         
        if (firstTry || dR < p_lep.DeltaR(p_tp)){
          tp = particle;
          p_tp = p_particle;
          firstTry = false;
        }
        else continue;
      }
    } else if (mcTruthMatch=="MCTCZ")
    {
      tp = m_truthClassifier->getGenPart();
    }

    if(tp!=NULL) {
      // Info("fillLepton(p)", "Found particle");
      l.truthI=addTruthPar(tp, m_susyEvt->truths, 0);
      m_susyEvt->truths[l.truthI].matchI = index;
    }
    else l.truthI=-1;
  }
  return EL::StatusCode::SUCCESS;
}
int ssEvtSelection::addTruthPar(const xAOD::TruthParticle* p, TRUTHS& v, int pLevel){
  /// check if already exist
  if(!p) return -1;
  if(p->status()!=1 && p->status()!=2 && p->status()!=4) return -1;
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
  t.status = p->status();

  /*
  std::pair<MCTruthPartClassifier::ParticleType, MCTruthPartClassifier::ParticleOrigin> res;
  res = m_truthClassifier->particleTruthClassifier(p);
  t.particleType = res.first;
  t.particleOrigin = res.second;
  */

  /// add first parent and children if exist
  if(pLevel){
    auto m = p->parent();
    if(m){
      while(m->parent() && m->pdgId() == m->parent()->pdgId()) m=m->parent();
      int id = addTruthPar(m, v, pLevel-1);
      v[nv].motherI = id;

      for(uint i=0; i<p->nChildren(); i++) addTruthPar(p->child(i), v, pLevel+1);
    }
  }

  return nv;
}

bool ssEvtSelection :: InitializeTriggerTools(std::string var){

  const char* APP_NAME = "InitializeTriggerTools()";
  LegsPerTool.clear();
  Info( APP_NAME, "Starting to initialize single lepton efficiency tools for GlobalEfficiencyCorrectionTool[%s]", var.c_str() );

  electronEffTools.clear();
  electronSFTools.clear();
  muonTools.clear();

  //Trigger efficiency tool for single lepton trigger
  AsgElectronEfficiencyCorrectionTool* elEffTool1 = new AsgElectronEfficiencyCorrectionTool("ElTrigEff-1"+var);
  elEffTool1->msg().setLevel( MSG::FATAL );
  elEffTool1->setProperty("MapFilePath","/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/ElectronEfficiencyCorrection/2015_2016/rel20.7/Moriond_February2017_v1/map0.txt").ignore();
  elEffTool1->setProperty("TriggerKey","Eff_SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0").ignore();
  elEffTool1->setProperty("IdKey","Medium").ignore();
  elEffTool1->setProperty("IsoKey","FixedCutTight").ignore();
  elEffTool1->setProperty("CorrelationModel","TOTAL").ignore();
  if (CF_isMC){
       if (CF_isMC==2) elEffTool1->setProperty("ForceDataType", (int) PATCore::ParticleDataType::Fast).ignore();
       else elEffTool1->setProperty("ForceDataType", (int) PATCore::ParticleDataType::Full).ignore();
  }

  if(elEffTool1->initialize() != StatusCode::SUCCESS){ Error(APP_NAME, "Unable to initialize ElectronTriggerEfficiencyTool-1"); return false;}
  else{ Info(APP_NAME, "Initialized ElectronTriggerEfficiencyTool-1");}
  LegsPerTool["ElTrigEff-1"+var] = "e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose,e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0";

  //Trigger SF tool for single lepton trigger
  AsgElectronEfficiencyCorrectionTool* elSFTool1 = new AsgElectronEfficiencyCorrectionTool("ElTrigSF-1"+var);
  elSFTool1->msg().setLevel( MSG::FATAL );
  elSFTool1->setProperty("MapFilePath","/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/ElectronEfficiencyCorrection/2015_2016/rel20.7/Moriond_February2017_v1/map0.txt").ignore();
  elSFTool1->setProperty("TriggerKey","SINGLE_E_2015_e24_lhmedium_L1EM20VH_OR_e60_lhmedium_OR_e120_lhloose_2016_e26_lhtight_nod0_ivarloose_OR_e60_lhmedium_nod0_OR_e140_lhloose_nod0").ignore();
  elSFTool1->setProperty("IdKey","Medium").ignore();
  elSFTool1->setProperty("IsoKey","FixedCutTight").ignore();
  elSFTool1->setProperty("CorrelationModel","TOTAL").ignore();

  if (CF_isMC){
     if (CF_isMC==2) elSFTool1->setProperty("ForceDataType", (int) PATCore::ParticleDataType::Fast).ignore(); 
     else elSFTool1->setProperty("ForceDataType", (int) PATCore::ParticleDataType::Full).ignore();
  }
  if(elSFTool1->initialize() != StatusCode::SUCCESS){ Error(APP_NAME, "Unable to initialize ElectronTriggerSFTool-1"); return false;}
  else{ Info(APP_NAME, "Initialized ElectronTriggerSFTool-1");}
  LegsPerTool["ElTrigSF-1"+var] = LegsPerTool["ElTrigEff-1"+var];

  //Trigger efficiency tool for dielectron trigger 
  AsgElectronEfficiencyCorrectionTool* elEffTool2 = new AsgElectronEfficiencyCorrectionTool("ElTrigEff-2"+var);
  elEffTool2->msg().setLevel( MSG::FATAL );
  elEffTool2->setProperty("MapFilePath","/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/ElectronEfficiencyCorrection/2015_2016/rel20.7/Moriond_February2017_v1/map0.txt").ignore();
  elEffTool2->setProperty("TriggerKey","Eff_DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0").ignore();
  elEffTool2->setProperty("IdKey","Medium").ignore();
  elEffTool2->setProperty("IsoKey","FixedCutTight").ignore();
  elEffTool2->setProperty("CorrelationModel","TOTAL").ignore();
  if (CF_isMC){
     if (CF_isMC==2) elEffTool2->setProperty("ForceDataType", (int) PATCore::ParticleDataType::Fast).ignore();
     else elEffTool2->setProperty("ForceDataType", (int) PATCore::ParticleDataType::Full).ignore();
   }
  if(elEffTool2->initialize() != StatusCode::SUCCESS){ Error(APP_NAME, "Unable to initialize ElectronTriggerEfficiencyTool-2"); return false;}
  else{ Info(APP_NAME, "Initialized ElectronTriggerEfficiencyTool-2");}
  LegsPerTool["ElTrigEff-2"+var] = "e12_lhloose_L1EM10VH,e17_lhvloose_nod0";

  //Trigger SF tool for dilelectron trigger 
  AsgElectronEfficiencyCorrectionTool* elSFTool2 = new AsgElectronEfficiencyCorrectionTool("ElTrigSF-2"+var);
  elSFTool2->msg().setLevel( MSG::FATAL );
  elSFTool2->setProperty("MapFilePath","/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/ElectronEfficiencyCorrection/2015_2016/rel20.7/Moriond_February2017_v1/map0.txt").ignore();
  elSFTool2->setProperty("TriggerKey","DI_E_2015_e12_lhloose_L1EM10VH_2016_e17_lhvloose_nod0").ignore();
  elSFTool2->setProperty("IdKey","Medium").ignore();
  elSFTool2->setProperty("IsoKey","FixedCutTight").ignore();
  elSFTool2->setProperty("CorrelationModel","TOTAL").ignore();
  if (CF_isMC){
     if (CF_isMC==2)  elSFTool2->setProperty("ForceDataType", (int) PATCore::ParticleDataType::Fast).ignore();
     else elSFTool2->setProperty("ForceDataType", (int) PATCore::ParticleDataType::Full).ignore();
  }
  if(elSFTool2->initialize() != StatusCode::SUCCESS){ Error(APP_NAME, "Unable to initialize ElectronTriggerSFTool-2"); return false;}
  else{ Info(APP_NAME, "Initialized ElectronTriggerSFTool-2");}
  LegsPerTool["ElTrigSF-2"+var] = LegsPerTool["ElTrigEff-2"+var];

  //Trigger efficiency tool for dilepton trigger (e17_mu14)
  AsgElectronEfficiencyCorrectionTool* elEffTool3 = new AsgElectronEfficiencyCorrectionTool("ElTrigEff-3"+var);
  elEffTool3->msg().setLevel( MSG::FATAL );
  elEffTool3->setProperty("MapFilePath","/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/ElectronEfficiencyCorrection/2015_2016/rel20.7/Moriond_February2017_v1/map0.txt").ignore();
  elEffTool3->setProperty("TriggerKey","Eff_MULTI_L_2015_e17_lhloose_2016_e17_lhloose_nod0").ignore();
  elEffTool3->setProperty("IdKey","Medium").ignore();
  elEffTool3->setProperty("IsoKey","FixedCutTight").ignore();
  elEffTool3->setProperty("CorrelationModel","TOTAL").ignore();
  if (CF_isMC){
      if (CF_isMC==1) elEffTool3->setProperty("ForceDataType", (int) PATCore::ParticleDataType::Full).ignore();
      else elEffTool3->setProperty("ForceDataType", (int) PATCore::ParticleDataType::Fast).ignore();
   }
  if(elEffTool3->initialize() != StatusCode::SUCCESS){ Error(APP_NAME, "Unable to initialize ElectronTriggerEfficiencyTool-3"); return false;}
  else{ Info(APP_NAME, "Initialized ElectronTriggerEfficiencyTool-3");}
  LegsPerTool["ElTrigEff-3"+var] = "e17_lhloose,e17_lhloose_nod0";

  //Trigger SF tool for dilepton trigger (e17_mu14)
  AsgElectronEfficiencyCorrectionTool* elSFTool3 = new AsgElectronEfficiencyCorrectionTool("ElTrigSF-3"+var);
  elSFTool3->msg().setLevel( MSG::FATAL );
  elSFTool3->setProperty("MapFilePath","/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/ElectronEfficiencyCorrection/2015_2016/rel20.7/Moriond_February2017_v1/map0.txt").ignore();
  elSFTool3->setProperty("TriggerKey","MULTI_L_2015_e17_lhloose_2016_e17_lhloose_nod0").ignore();
  elSFTool3->setProperty("IdKey","Medium").ignore();
  elSFTool3->setProperty("IsoKey","FixedCutTight").ignore();
  elSFTool3->setProperty("CorrelationModel","TOTAL").ignore();

  if (CF_isMC){
     if (CF_isMC==2) elSFTool3->setProperty("ForceDataType", (int) PATCore::ParticleDataType::Fast).ignore();
     else elSFTool3->setProperty("ForceDataType", (int) PATCore::ParticleDataType::Full).ignore();
  }
  if(elSFTool3->initialize() != StatusCode::SUCCESS){ Error(APP_NAME, "Unable to initialize ElectronTriggerSFTool-3"); return false;}
  else{ Info(APP_NAME, "Initialized ElectronTriggerSFTool-3");}
  LegsPerTool["ElTrigSF-3"+var] = LegsPerTool["ElTrigEff-3"+var];

  //Trigger efficiency tool for dilepton trigger (e7_mu24)
  AsgElectronEfficiencyCorrectionTool* elEffTool4 = new AsgElectronEfficiencyCorrectionTool("ElTrigEff-4"+var);
  elEffTool4->msg().setLevel( MSG::FATAL );
  elEffTool4->setProperty("MapFilePath","/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/ElectronEfficiencyCorrection/2015_2016/rel20.7/Moriond_February2017_v1/map0.txt").ignore();
  elEffTool4->setProperty("TriggerKey","Eff_MULTI_L_2015_e7_lhmedium_2016_e7_lhmedium_nod0").ignore();
  elEffTool4->setProperty("IdKey","Medium").ignore();
  elEffTool4->setProperty("IsoKey","FixedCutTight").ignore();
  elEffTool4->setProperty("CorrelationModel","TOTAL").ignore();

  if (CF_isMC){
      if (CF_isMC==2) elEffTool4->setProperty("ForceDataType", (int) PATCore::ParticleDataType::Fast).ignore();
      else elEffTool4->setProperty("ForceDataType", (int) PATCore::ParticleDataType::Full).ignore();
  }
  if(elEffTool4->initialize() != StatusCode::SUCCESS){ Error(APP_NAME, "Unable to initialize ElectronTriggerEfficiencyTool-4"); return false;}
  else{ Info(APP_NAME, "Initialized ElectronTriggerEfficiencyTool-4");}
  LegsPerTool["ElTrigEff-4"+var] = "e7_lhmedium,e7_lhmedium_nod0";

  //Trigger SF tool for dilepton trigger (e7_mu24)
  AsgElectronEfficiencyCorrectionTool* elSFTool4 = new AsgElectronEfficiencyCorrectionTool("ElTrigSF-4"+var);
  elSFTool4->msg().setLevel( MSG::FATAL );
  elSFTool4->setProperty("MapFilePath","/cvmfs/atlas.cern.ch/repo/sw/database/GroupData/ElectronEfficiencyCorrection/2015_2016/rel20.7/Moriond_February2017_v1/map0.txt").ignore();
  elSFTool4->setProperty("TriggerKey","MULTI_L_2015_e7_lhmedium_2016_e7_lhmedium_nod0").ignore();
  elSFTool4->setProperty("IdKey","Medium").ignore();
  elSFTool4->setProperty("IsoKey","FixedCutTight").ignore();
  elSFTool4->setProperty("CorrelationModel","TOTAL").ignore();

  if (CF_isMC){
     if (CF_isMC==2) elSFTool4->setProperty("ForceDataType", (int) PATCore::ParticleDataType::Fast).ignore();
     else elSFTool4->setProperty("ForceDataType", (int) PATCore::ParticleDataType::Full).ignore();
  }
  if(elSFTool4->initialize() != StatusCode::SUCCESS){ Error(APP_NAME, "Unable to initialize ElectronTriggerSFTool-4"); return false;}
  else{ Info(APP_NAME, "Initialized ElectronTriggerSFTool-4");}
  LegsPerTool["ElTrigSF-4"+var] = LegsPerTool["ElTrigEff-4"+var];

  // Trigger Efficiency tool for muon (2015)
  CP::MuonTriggerScaleFactors* toolMuons1 =  new CP::MuonTriggerScaleFactors("MuonTrigEff-2015"+var);
  toolMuons1->msg().setLevel( MSG::FATAL );
  toolMuons1->setProperty("CalibrationRelease", "170209_Moriond").ignore();
  toolMuons1->setProperty("MuonQuality", "Tight").ignore();
  toolMuons1->setProperty("Isolation", "GradientLoose").ignore();
  toolMuons1->setProperty("Year", "2015").ignore();
  if(toolMuons1->initialize() != StatusCode::SUCCESS){ Error(APP_NAME, "Unable to initialize MuonTriggerScaleFactorsTool"); return false;}
  else{ Info(APP_NAME, "Initialized MuonTriggerScaleFactorsTool for 2015");}

  // Trigger Efficiency tool for muon (2016)
  CP::MuonTriggerScaleFactors* toolMuons2 =  new CP::MuonTriggerScaleFactors("MuonTrigEff-2016"+var);
  toolMuons2->msg().setLevel( MSG::FATAL );
  toolMuons2->setProperty("CalibrationRelease", "170209_Moriond").ignore();
  toolMuons2->setProperty("MuonQuality", "Tight").ignore();
  toolMuons2->setProperty("Isolation", "GradientLoose").ignore();
  toolMuons2->setProperty("Year", "2016").ignore();
  if(toolMuons2->initialize() != StatusCode::SUCCESS){ Error(APP_NAME, "Unable to initialize MuonTriggerScaleFactorsTool"); return false;}

  electronEffTools.push_back(elEffTool1);
  electronEffTools.push_back(elEffTool2);
  electronEffTools.push_back(elEffTool3);
  electronEffTools.push_back(elEffTool4);
  electronSFTools.push_back(elSFTool1);
  electronSFTools.push_back(elSFTool2);
  electronSFTools.push_back(elSFTool3);
  electronSFTools.push_back(elSFTool4);
  muonTools.push_back(toolMuons1);
  muonTools.push_back(toolMuons2);

  Info( APP_NAME, "Initialized single lepton efficiency tools for GlobalEfficiencyCorrectionTool[%s]", var.c_str() );
  return true;
}
