#include <TError.h>
#include <algorithm>
#include <vector>
#include <EventLoop/StatusCode.h>

#include <multiLepSearch/mCHECK.h>
#include <multiLepSearch/ssEvtPostProc1.h>
#include <multiLepSearch/ChargeFlipBkgTool.h>
#include <multiLepSearch/FakeLepBkgTool.h>
#include <multiLepSearch/obj_def.h>

#include <TH1D.h>
#include <TFile.h>
#include <TRegexp.h>
#include <TKey.h>

using namespace std;

// this is needed to distribute the algorithm to the workers
ClassImp(ssEvtPostProc1)

ssEvtPostProc1 :: ssEvtPostProc1(string name):m_name(name){
  inputFile="dummy.root";
  norminalTreeName="evt2l";
  isData=true;

  //m_susyEvt = new susyEvts();
}

ssEvtPostProc1 :: ~ssEvtPostProc1(){
}


int ssEvtPostProc1 :: initialize ()
{
  // ChargeFlipBkgTool
  CHECK(mChargeFlipBkgTool->initialize());
  //mChargeFlipBkgTool->msg().setLevel( MSG::VERBOSE );

  // FakeLepBkgTool
  CHECK(mFakeLepBkgTool->initialize());
  //mFakeLepBkgTool->msg().setLevel( MSG::VERBOSE );

  inF = new TFile(inputFile.c_str(), "UPDATE");

  //Loop over objects in the file
  TIter next(inF->GetListOfKeys());
  TKey* key;
  int nKeys = inF->GetNkeys();
  for (int iKey=0; iKey<nKeys; iKey++){
    key = (TKey *) next();
    if (!(TString(key->GetClassName())=="TTree")) continue;

    TTree* rawTree = (TTree*)key->ReadObj();
    TString treeName = key->GetName();
    if (treeName.Index(norminalTreeName)!=0) continue; //only look at norminal tree and its systematics

    TString outTreeName = "PP1_"+treeName;

    //Check if the tree needs update
    bool needUpdate = false;
    if (treeName==norminalTreeName && isData ) needUpdate = true;

    if (!needUpdate){
      //no updated need: make an empty tree and addfriend the raw tree
      TTree* outTree = new TTree( outTreeName, "Added PostProcess1 variables");
      outTree->SetDirectory(inF);
      outTree->AddFriend(rawTree);
      int nEntries = rawTree->GetEntries();
      for (int i=0;i<nEntries;i++){outTree->Fill();}
      outTree->Write();
      delete outTree;
    }else{
      //need update: make a new tree and addfriend the raw tree
      susyEvts* rawSusyEvts = new susyEvts(rawTree);
      rawSusyEvtsList.push_back(rawSusyEvts);

      susyEvts* outSusyEvts                   = new susyEvts(); //new nominal tree with qFlip and fakeLep weight updated
      susyEvts* outSusyEvts_CFLIP_SYS__1up    = new susyEvts();
      susyEvts* outSusyEvts_CFLIP_SYS__1dn    = new susyEvts();
      susyEvts* outSusyEvts_FAKE_LEP_SYS__1up = new susyEvts();
      susyEvts* outSusyEvts_FAKE_LEP_SYS__1dn = new susyEvts();

      outSusyEvts                  ->makeWeightOnlyTree(outTreeName , rawSusyEvts);
      outSusyEvts_CFLIP_SYS__1up   ->makeWeightOnlyTree(outTreeName + "_CFLIP_SYS__1up", rawSusyEvts);
      outSusyEvts_CFLIP_SYS__1dn   ->makeWeightOnlyTree(outTreeName + "_CFLIP_SYS__1down", rawSusyEvts);
      outSusyEvts_FAKE_LEP_SYS__1up->makeWeightOnlyTree(outTreeName + "_FAKE_LEP_SYS__1up", rawSusyEvts);
      outSusyEvts_FAKE_LEP_SYS__1dn->makeWeightOnlyTree(outTreeName + "_FAKE_LEP_SYS__1down", rawSusyEvts);

      outSusyEvtsList                  .push_back(outSusyEvts                  );
      outSusyEvtsList_CFLIP_SYS__1up   .push_back(outSusyEvts_CFLIP_SYS__1up   );
      outSusyEvtsList_CFLIP_SYS__1dn   .push_back(outSusyEvts_CFLIP_SYS__1dn   );
      outSusyEvtsList_FAKE_LEP_SYS__1up.push_back(outSusyEvts_FAKE_LEP_SYS__1up);
      outSusyEvtsList_FAKE_LEP_SYS__1dn.push_back(outSusyEvts_FAKE_LEP_SYS__1dn);

      outSusyEvts                  ->tree2->SetDirectory(inF);
      outSusyEvts_CFLIP_SYS__1up   ->tree2->SetDirectory(inF);
      outSusyEvts_CFLIP_SYS__1dn   ->tree2->SetDirectory(inF);
      outSusyEvts_FAKE_LEP_SYS__1up->tree2->SetDirectory(inF);
      outSusyEvts_FAKE_LEP_SYS__1dn->tree2->SetDirectory(inF);
    }
  }

  return 0;
}


int ssEvtPostProc1 :: runLoop ()
{
  int nTrees = rawSusyEvtsList.size();
  for (int i=0; i<nTrees; i++){
    susyEvts* rawSusyEvts = rawSusyEvtsList[i];
    susyEvts* outSusyEvts                   = outSusyEvtsList                  [i];
    susyEvts* outSusyEvts_CFLIP_SYS__1up    = outSusyEvtsList_CFLIP_SYS__1up   [i];
    susyEvts* outSusyEvts_CFLIP_SYS__1dn    = outSusyEvtsList_CFLIP_SYS__1dn   [i];
    susyEvts* outSusyEvts_FAKE_LEP_SYS__1up = outSusyEvtsList_FAKE_LEP_SYS__1up[i];
    susyEvts* outSusyEvts_FAKE_LEP_SYS__1dn = outSusyEvtsList_FAKE_LEP_SYS__1dn[i];

    cout << "PostProcessing " << rawSusyEvts->tree1->GetName() << endl;
    vector<double> eta, pt;
    vector<int> pdgID;
    vector<bool> isSig;

    while (rawSusyEvts->Next()>0){
      outSusyEvts                  ->evt = rawSusyEvts->evt;
      outSusyEvts_CFLIP_SYS__1up   ->evt = rawSusyEvts->evt;
      outSusyEvts_CFLIP_SYS__1dn   ->evt = rawSusyEvts->evt;
      outSusyEvts_FAKE_LEP_SYS__1up->evt = rawSusyEvts->evt;
      outSusyEvts_FAKE_LEP_SYS__1dn->evt = rawSusyEvts->evt;

      //look for OS 2 signal leptons to get charge flip weight
      if ( 
             rawSusyEvts->leps.size()>=2 && 
            (rawSusyEvts->leps[0].lFlag & IS_SIGNAL) &&
	    (rawSusyEvts->leps[1].lFlag & IS_SIGNAL) && 
           ((rawSusyEvts->leps[0].ID>0) != (rawSusyEvts->leps[1].ID>0))
      ){
        eta.clear();
	pt.clear();

        for (int idx=0;idx<1;idx++){
          if (abs(rawSusyEvts->leps[idx].ID)==11000){
	    eta.push_back( rawSusyEvts->leps[idx].eta );
	    pt.push_back( rawSusyEvts->leps[idx].pt );
	  }
        }

        //for (auto aLep: rawSusyEvts->leps){
        //  if (abs(aLep.ID)==11000){
	//    eta.push_back( aLep.eta );
	//    pt.push_back( aLep.pt );
	//  }
        //}

        outSusyEvts                  ->evt.qFwt = mChargeFlipBkgTool->GetWeight(eta, pt, 0,0);
        outSusyEvts_CFLIP_SYS__1up   ->evt.qFwt = mChargeFlipBkgTool->GetWeight(eta, pt, 1,0);
        outSusyEvts_CFLIP_SYS__1dn   ->evt.qFwt = mChargeFlipBkgTool->GetWeight(eta, pt,-1,0);
      }else{
        outSusyEvts                  ->evt.qFwt = 0.0;
        outSusyEvts_CFLIP_SYS__1up   ->evt.qFwt = 0.0;
        outSusyEvts_CFLIP_SYS__1dn   ->evt.qFwt = 0.0;
      }

      //look for 2LSS, no check on leptons loose/tight to get fake bkg weight
      if ( 
            rawSusyEvts->leps.size()>=2 && 
          ((rawSusyEvts->leps[0].ID>0) == (rawSusyEvts->leps[1].ID>0))
      ){
        eta.clear();
	pt.clear();
	pdgID.clear();
	isSig.clear();
        for (auto aLep: rawSusyEvts->leps){
          if      (abs(aLep.ID)==11000){ pdgID.push_back(11);}
	  else if (abs(aLep.ID)==13000){ pdgID.push_back(13);}
	  else {continue;}
	  eta.push_back( aLep.eta );
	  pt.push_back( aLep.pt );
	  isSig.push_back(aLep.lFlag & IS_SIGNAL);
        }
        outSusyEvts                  ->evt.fLwt = mFakeLepBkgTool->GetWeight(eta, pt, pdgID, isSig, 0,0);
        outSusyEvts_FAKE_LEP_SYS__1up->evt.fLwt = mFakeLepBkgTool->GetWeight(eta, pt, pdgID, isSig, 1,0);
        outSusyEvts_FAKE_LEP_SYS__1dn->evt.fLwt = mFakeLepBkgTool->GetWeight(eta, pt, pdgID, isSig,-1,0);
      }else{
        outSusyEvts                  ->evt.fLwt = 0.0;
        outSusyEvts_FAKE_LEP_SYS__1up->evt.fLwt = 0.0;
        outSusyEvts_FAKE_LEP_SYS__1dn->evt.fLwt = 0.0;
      }

      outSusyEvts                  ->fill();
      outSusyEvts_CFLIP_SYS__1up   ->fill();
      outSusyEvts_CFLIP_SYS__1dn   ->fill();
      outSusyEvts_FAKE_LEP_SYS__1up->fill();
      outSusyEvts_FAKE_LEP_SYS__1dn->fill();
    }

    outSusyEvts                  ->tree2->Write();
    outSusyEvts_CFLIP_SYS__1up   ->tree2->Write();
    outSusyEvts_CFLIP_SYS__1dn   ->tree2->Write();
    outSusyEvts_FAKE_LEP_SYS__1up->tree2->Write();
    outSusyEvts_FAKE_LEP_SYS__1dn->tree2->Write();
  }
  
  return 0;
}


int ssEvtPostProc1 :: finalize ()
{
  for (unsigned int i=0; i<rawSusyEvtsList.size();i++){
    delete rawSusyEvtsList                  [i];
    delete outSusyEvtsList                  [i];
    delete outSusyEvtsList_CFLIP_SYS__1up   [i];
    delete outSusyEvtsList_CFLIP_SYS__1dn   [i];
    delete outSusyEvtsList_FAKE_LEP_SYS__1up[i];
    delete outSusyEvtsList_FAKE_LEP_SYS__1dn[i];
  }
  rawSusyEvtsList                  .clear();
  outSusyEvtsList                  .clear();
  outSusyEvtsList_CFLIP_SYS__1up   .clear();
  outSusyEvtsList_CFLIP_SYS__1dn   .clear();
  outSusyEvtsList_FAKE_LEP_SYS__1up.clear();
  outSusyEvtsList_FAKE_LEP_SYS__1dn.clear();
  inF->Flush();
  inF->Close();
  return 0;
}


int ssEvtPostProc1 :: execute ()
{
  initialize();
  runLoop();
  finalize();
  std::cout << "done" << endl;
  return 0;
}



