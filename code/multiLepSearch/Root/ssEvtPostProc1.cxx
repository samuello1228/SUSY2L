#include <TError.h>
#include <algorithm>
#include <vector>
#include <EventLoop/StatusCode.h>

#include <multiLepSearch/mCHECK.h>
#include <multiLepSearch/ssEvtPostProc1.h>
#include <multiLepSearch/ChargeFlipBkgTool.h>
#include <multiLepSearch/FakeLepBkgTool.h>
#include <multiLepSearch/obj_def.h>
#include <multiLepSearch/anaHelper.h>

#include <TH1D.h>
#include <TFile.h>
#include <TRegexp.h>
#include <TKey.h>
#include <TLorentzVector.h>

using namespace std;

// this is needed to distribute the algorithm to the workers
ClassImp(ssEvtPostProc1)

ssEvtPostProc1 :: ssEvtPostProc1(string name):m_name(name){
  inputFile="dummy.root";
  norminalTreeName="evt2l";
  isData=true;
  isFirstInit=true;

  //m_susyEvt = new susyEvts();
}

ssEvtPostProc1 :: ~ssEvtPostProc1(){
}


int ssEvtPostProc1 :: initialize ()
{
  if (isFirstInit){
    // ChargeFlipBkgTool
    CHECK(mChargeFlipBkgTool->initialize());
    //mChargeFlipBkgTool->msg().setLevel( MSG::VERBOSE );

    // FakeLepBkgTool
    CHECK(mFakeLepBkgTool->initialize());
    //mFakeLepBkgTool->msg().setLevel( MSG::VERBOSE );

    isFirstInit = false;
  }

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
      susyEvts* outSusyEvts_FAKE_LEP_E_SYS__1up = new susyEvts();
      susyEvts* outSusyEvts_FAKE_LEP_E_SYS__1dn = new susyEvts();
      susyEvts* outSusyEvts_FAKE_LEP_U_SYS__1up = new susyEvts();
      susyEvts* outSusyEvts_FAKE_LEP_U_SYS__1dn = new susyEvts();

      outSusyEvts                  ->makePtCorrTree(outTreeName , rawSusyEvts);
      outSusyEvts_CFLIP_SYS__1up   ->makePtCorrTree(outTreeName + "_CFLIP_SYS__1up", rawSusyEvts);
      outSusyEvts_CFLIP_SYS__1dn   ->makePtCorrTree(outTreeName + "_CFLIP_SYS__1down", rawSusyEvts);
      outSusyEvts_FAKE_LEP_E_SYS__1up->makeWeightOnlyTree(outTreeName + "_FAKE_LEP_E_SYS__1up", rawSusyEvts);
      outSusyEvts_FAKE_LEP_E_SYS__1dn->makeWeightOnlyTree(outTreeName + "_FAKE_LEP_E_SYS__1down", rawSusyEvts);
      outSusyEvts_FAKE_LEP_U_SYS__1up->makeWeightOnlyTree(outTreeName + "_FAKE_LEP_U_SYS__1up", rawSusyEvts);
      outSusyEvts_FAKE_LEP_U_SYS__1dn->makeWeightOnlyTree(outTreeName + "_FAKE_LEP_U_SYS__1down", rawSusyEvts);

      outSusyEvtsList                  .push_back(outSusyEvts                  );
      outSusyEvtsList_CFLIP_SYS__1up   .push_back(outSusyEvts_CFLIP_SYS__1up   );
      outSusyEvtsList_CFLIP_SYS__1dn   .push_back(outSusyEvts_CFLIP_SYS__1dn   );
      outSusyEvtsList_FAKE_LEP_E_SYS__1up.push_back(outSusyEvts_FAKE_LEP_E_SYS__1up);
      outSusyEvtsList_FAKE_LEP_E_SYS__1dn.push_back(outSusyEvts_FAKE_LEP_E_SYS__1dn);
      outSusyEvtsList_FAKE_LEP_U_SYS__1up.push_back(outSusyEvts_FAKE_LEP_U_SYS__1up);
      outSusyEvtsList_FAKE_LEP_U_SYS__1dn.push_back(outSusyEvts_FAKE_LEP_U_SYS__1dn);

      outSusyEvts                  ->tree2->SetDirectory(inF);
      outSusyEvts_CFLIP_SYS__1up   ->tree2->SetDirectory(inF);
      outSusyEvts_CFLIP_SYS__1dn   ->tree2->SetDirectory(inF);
      outSusyEvts_FAKE_LEP_E_SYS__1up->tree2->SetDirectory(inF);
      outSusyEvts_FAKE_LEP_E_SYS__1dn->tree2->SetDirectory(inF);
      outSusyEvts_FAKE_LEP_U_SYS__1up->tree2->SetDirectory(inF);
      outSusyEvts_FAKE_LEP_U_SYS__1dn->tree2->SetDirectory(inF);
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
    susyEvts* outSusyEvts_FAKE_LEP_E_SYS__1up = outSusyEvtsList_FAKE_LEP_E_SYS__1up[i];
    susyEvts* outSusyEvts_FAKE_LEP_E_SYS__1dn = outSusyEvtsList_FAKE_LEP_E_SYS__1dn[i];
    susyEvts* outSusyEvts_FAKE_LEP_U_SYS__1up = outSusyEvtsList_FAKE_LEP_U_SYS__1up[i];
    susyEvts* outSusyEvts_FAKE_LEP_U_SYS__1dn = outSusyEvtsList_FAKE_LEP_U_SYS__1dn[i];

    cout << "PostProcessing " << rawSusyEvts->tree1->GetName() << endl;
    vector<double> eta, pt;
    vector<int> pdgID;
    vector<bool> isSig;

    while (rawSusyEvts->Next()>0){
      outSusyEvts                  ->evt   = rawSusyEvts->evt;
      outSusyEvts_CFLIP_SYS__1up   ->evt   = rawSusyEvts->evt;
      outSusyEvts_CFLIP_SYS__1dn   ->evt   = rawSusyEvts->evt;
      outSusyEvts_FAKE_LEP_E_SYS__1up->evt = rawSusyEvts->evt;
      outSusyEvts_FAKE_LEP_E_SYS__1dn->evt = rawSusyEvts->evt;
      outSusyEvts_FAKE_LEP_U_SYS__1up->evt = rawSusyEvts->evt;
      outSusyEvts_FAKE_LEP_U_SYS__1dn->evt = rawSusyEvts->evt;

      outSusyEvts                  ->leps = rawSusyEvts->leps;
      outSusyEvts_CFLIP_SYS__1up   ->leps = rawSusyEvts->leps;
      outSusyEvts_CFLIP_SYS__1dn   ->leps = rawSusyEvts->leps;

      outSusyEvts                  ->l12 = rawSusyEvts->l12;
      outSusyEvts_CFLIP_SYS__1up   ->l12 = rawSusyEvts->l12;
      outSusyEvts_CFLIP_SYS__1dn   ->l12 = rawSusyEvts->l12;

      outSusyEvts                  ->sig = rawSusyEvts->sig;
      outSusyEvts_CFLIP_SYS__1up   ->sig = rawSusyEvts->sig;
      outSusyEvts_CFLIP_SYS__1dn   ->sig = rawSusyEvts->sig;

      outSusyEvts                  ->jets = rawSusyEvts->jets; //not updated but needed for HT recalulation
      outSusyEvts_CFLIP_SYS__1up   ->jets = rawSusyEvts->jets;
      outSusyEvts_CFLIP_SYS__1dn   ->jets = rawSusyEvts->jets;

      //look for OS 2 signal leptons to get charge flip weight
      if ( 
             rawSusyEvts->leps.size()==2 && 
            (rawSusyEvts->leps[0].lFlag & IS_SIGNAL) &&
	    (rawSusyEvts->leps[1].lFlag & IS_SIGNAL) && 
           ((rawSusyEvts->leps[0].ID>0) != (rawSusyEvts->leps[1].ID>0))
      ){
        eta.clear();
	pt.clear();

        for (uint idx=0;idx<rawSusyEvts->leps.size();idx++){
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

        //apply pT corr from change flip
        if (abs(rawSusyEvts->leps[0].ID)==11000 && abs(rawSusyEvts->leps[1].ID)==11000){
          vector<double> ptCorr    = pt;
          vector<double> ptCorr1up = pt;
          vector<double> ptCorr1dn = pt;

          mChargeFlipBkgTool->GetCorrectedPt(eta, ptCorr   , 0,0); //this modifies pt in place
          mChargeFlipBkgTool->GetCorrectedPt(eta, ptCorr1up, 1,0); //this modifies pt in place
          mChargeFlipBkgTool->GetCorrectedPt(eta, ptCorr1dn,-1,0); //this modifies pt in place

          //cout << "pT corr " << pt[0] << " " << pt[1] << " " << ptCorr[0] << " " << ptCorr[1] << endl;

          for (uint idx=0;idx<pt.size();idx++){
            outSusyEvts                  ->leps[idx].pt = ptCorr   [idx];
            outSusyEvts_CFLIP_SYS__1up   ->leps[idx].pt = ptCorr1up[idx];
            outSusyEvts_CFLIP_SYS__1dn   ->leps[idx].pt = ptCorr1dn[idx];
          }

        }else if (abs(rawSusyEvts->leps[0].ID)==11000 || abs(rawSusyEvts->leps[1].ID)==11000){
          int eleIdx = abs(rawSusyEvts->leps[0].ID)==11000 ? 0 : 1;
          outSusyEvts                  ->leps[eleIdx].pt = mChargeFlipBkgTool->GetCorrectedPt(eta[eleIdx], pt[eleIdx] , 0,0);
          outSusyEvts_CFLIP_SYS__1up   ->leps[eleIdx].pt = mChargeFlipBkgTool->GetCorrectedPt(eta[eleIdx], pt[eleIdx] , 1,0);
          outSusyEvts_CFLIP_SYS__1dn   ->leps[eleIdx].pt = mChargeFlipBkgTool->GetCorrectedPt(eta[eleIdx], pt[eleIdx] ,-1,0);

          //cout << "pT corr " << pt[eleIdx] << " " << outSusyEvts->leps[eleIdx].pt << " " << endl;
        }

        recalPtRelatedVar (outSusyEvts                );
        recalPtRelatedVar (outSusyEvts_CFLIP_SYS__1up );
        recalPtRelatedVar (outSusyEvts_CFLIP_SYS__1dn );


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
        outSusyEvts                    ->evt.fLwt = mFakeLepBkgTool->GetWeight(eta, pt, pdgID, isSig, 0,0);
        outSusyEvts_FAKE_LEP_E_SYS__1up->evt.fLwt = mFakeLepBkgTool->GetWeight(eta, pt, pdgID, isSig, 1,0);
        outSusyEvts_FAKE_LEP_E_SYS__1dn->evt.fLwt = mFakeLepBkgTool->GetWeight(eta, pt, pdgID, isSig,-1,0);
        outSusyEvts_FAKE_LEP_U_SYS__1up->evt.fLwt = mFakeLepBkgTool->GetWeight(eta, pt, pdgID, isSig, 1,1);
        outSusyEvts_FAKE_LEP_U_SYS__1dn->evt.fLwt = mFakeLepBkgTool->GetWeight(eta, pt, pdgID, isSig,-1,1);
      }else{
        outSusyEvts                    ->evt.fLwt = 0.0;
        outSusyEvts_FAKE_LEP_E_SYS__1up->evt.fLwt = 0.0;
        outSusyEvts_FAKE_LEP_E_SYS__1dn->evt.fLwt = 0.0;
        outSusyEvts_FAKE_LEP_U_SYS__1up->evt.fLwt = 0.0;
        outSusyEvts_FAKE_LEP_U_SYS__1dn->evt.fLwt = 0.0;
      }

      outSusyEvts                  ->fill();
      outSusyEvts_CFLIP_SYS__1up   ->fill();
      outSusyEvts_CFLIP_SYS__1dn   ->fill();
      outSusyEvts_FAKE_LEP_E_SYS__1up->fill();
      outSusyEvts_FAKE_LEP_E_SYS__1dn->fill();
      outSusyEvts_FAKE_LEP_U_SYS__1up->fill();
      outSusyEvts_FAKE_LEP_U_SYS__1dn->fill();
    }
    outSusyEvts                  ->tree2->Write();
    outSusyEvts_CFLIP_SYS__1up   ->tree2->Write();
    outSusyEvts_CFLIP_SYS__1dn   ->tree2->Write();
    outSusyEvts_FAKE_LEP_E_SYS__1up->tree2->Write();
    outSusyEvts_FAKE_LEP_E_SYS__1dn->tree2->Write();
    outSusyEvts_FAKE_LEP_U_SYS__1up->tree2->Write();
    outSusyEvts_FAKE_LEP_U_SYS__1dn->tree2->Write();
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
    delete outSusyEvtsList_FAKE_LEP_E_SYS__1up[i];
    delete outSusyEvtsList_FAKE_LEP_E_SYS__1dn[i];
    delete outSusyEvtsList_FAKE_LEP_U_SYS__1up[i];
    delete outSusyEvtsList_FAKE_LEP_U_SYS__1dn[i];
  }
  rawSusyEvtsList                  .clear();
  outSusyEvtsList                  .clear();
  outSusyEvtsList_CFLIP_SYS__1up   .clear();
  outSusyEvtsList_CFLIP_SYS__1dn   .clear();
  outSusyEvtsList_FAKE_LEP_E_SYS__1up.clear();
  outSusyEvtsList_FAKE_LEP_E_SYS__1dn.clear();
  outSusyEvtsList_FAKE_LEP_U_SYS__1up.clear();
  outSusyEvtsList_FAKE_LEP_U_SYS__1dn.clear();
  inF->Flush();
  inF->Close();
  return 0;
}


int ssEvtPostProc1 :: execute ()
{
  initialize();
  runLoop();
  finalize();
  std::cout << "done " << inputFile << endl;
  return 0;
}

int ssEvtPostProc1 :: recalPtRelatedVar (susyEvts* inTree)
{
  TLorentzVector l1, l2;
  L_PAR *tmpLep = NULL;
  double tmpM;

  tmpLep = &inTree->leps[0]; tmpM = tmpLep->ID==11000? 0.000511 : 0.105658;
  l1.SetPtEtaPhiM(tmpLep->pt, tmpLep->eta, tmpLep->phi, tmpM);
  tmpLep = &inTree->leps[1]; tmpM = tmpLep->ID==11000? 0.000511 : 0.105658;
  l2.SetPtEtaPhiM(tmpLep->pt, tmpLep->eta, tmpLep->phi, tmpM);

  TLorentzVector ll(l1+l2);
  TLorentzVector metV(inTree->sig.MetX, inTree->sig.MetY, 0, inTree->sig.Met);

  inTree->l12.m = ll.M();
  inTree->l12.pt = ll.Pt(); 
  inTree->l12.eta = ll.Eta(); 
  inTree->l12.phi = ll.Phi(); 
  //inTree->l12.dPhi = l1.DeltaPhi(l2); 
  //inTree->l12.dR = l1.DeltaR(l2); 
  //inTree->l12.MET_dPhi = metV.DeltaPhi(ll);

  inTree->sig.HT = 0.0;
  inTree->sig.HT += (l1.Pt()+l2.Pt());
  for (auto j : inTree->jets){ inTree->sig.HT += j.pt;} //FIXME somehow this don't add jet pt
  
  double CF_mT2_m0 = 1.;
  inTree->sig.mT2 =  anaHelper::get_mT2(l1.M()*1000. , l1.Px()*1000. , l1.Py()*1000. ,
                                        l2.M()*1000. , l2.Px()*1000. , l2.Py()*1000. ,
                                        metV.Px()*1000. , metV.Py()*1000., CF_mT2_m0*1000., CF_mT2_m0*1000.)*0.001;
  return 0;
}


