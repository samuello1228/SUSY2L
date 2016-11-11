#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <truthCheck/mcChecker.h>

#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include <xAODTruth/TruthEventContainer.h>
#include <TFile.h>
#include <TTree.h>
#include "xAODJet/Jet.h"
#include "xAODJet/JetContainer.h"
// #include <AsgTools/MessageCheck.h>


const float iGeV = 0.001;
const float GeV = 1000.;
// this is needed to distribute the algorithm to the workers
ClassImp(mcChecker)

mcChecker :: mcChecker ()
{
  // Here you put any code for the base initialization of variables,
  // e.g. initialize all pointers to 0.  Note that you should only put
  // the most basic initialization here, since this method will be
  // called on both the submission and the worker node.  Most of your
  // initialization code will go into histInitialize() and
  // initialize().
}



EL::StatusCode mcChecker :: setupJob (EL::Job& job)
{
  // Here you put code that sets up the job on the submission object
  // so that it is ready to work with your algorithm, e.g. you can
  // request the D3PDReader service or add output files.  Any code you
  // put here could instead also go into the submission script.  The
  // sole advantage of putting it here is that it gets automatically
  // activated/deactivated when you add/remove the algorithm from your
  // job, which may or may not be of value to you.

  // let's initialize the algorithm to use the xAODRootAccess package
  job.useXAOD ();
  if(!xAOD::Init().isSuccess()) Error("setupJob", "xAOD::Init() not isSuccess()"); // call before opening first file

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode mcChecker :: histInitialize ()
{
  // Here you do everything that needs to be done at the very
  // beginning on each worker node, e.g. create histograms and output
  // trees.  This method gets called before any input files are
  // connected.
  TFile *outputFile = wk()->getOutputFile(outputName);
  tree = new TTree ("tree", "tree");
  tree->SetDirectory(outputFile);
  tree->Branch("C1", &pC1, kPar_s.c_str());
  tree->Branch("N2", &pN2, kPar_s.c_str());
  tree->Branch("N1w", &pN1w, kPar_s.c_str());
  tree->Branch("N1z", &pN1z, kPar_s.c_str());
  tree->Branch("W", &pW, mPar_s.c_str());
  tree->Branch("Z", &pZ, mPar_s.c_str());
  tree->Branch("l1", &pl1, xPar_s.c_str());
  tree->Branch("l2", &pl2, xPar_s.c_str());
  tree->Branch("l3", &pl3, xPar_s.c_str());
  tree->Branch("l4", &pl4, xPar_s.c_str()); 
  tree->Branch("Other", &pOther); 
  tree->Branch("Jets",  &pJets); 
  tree->Branch("mSC", &mSC);
  tree->Branch("mSS", &mSS);
  tree->Branch("MET", &MET);
  tree->Branch("MET_phi", &MET_phi);
  tree->Branch("SS_dR", &SS_dR);
  tree->Branch("SS_dPhi", &SS_dPhi);

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode mcChecker :: fileExecute ()
{
  // Here you do everything that needs to be done exactly once for every
  // single file, e.g. collect a list of all lumi-blocks processed
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode mcChecker :: changeInput (bool /*firstFile*/)
{
  // Here you do everything you need to do when we change input files,
  // e.g. resetting branch addresses on trees.  If you are using
  // D3PDReader or a similar service this method is not needed.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode mcChecker :: initialize ()
{
  // Here you do everything that you need to do after the first input
  // file has been connected and before the first event is processed,
  // e.g. create additional histograms based on which variables are
  // available in the input files.  You can also create all of your
  // histograms and trees in here, but be aware that this method
  // doesn't get called if no events are processed.  So any objects
  // you create here won't be available in the output if you have no
  // input events.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode mcChecker :: execute ()
{
  // Here you do everything that needs to be done on every single
  // events, e.g. read input variables, apply cuts, and fill
  // histograms and trees.  This is where most of your actual analysis
  // code will go.

  const xAOD::TruthEventContainer* mcEvents = 0;
  if(! wk()->xaodEvent()->retrieve(mcEvents, "TruthEvents").isSuccess()){
    Error("execute()", "Failed to retrieve TruthEvents collection. Exiting.");
    return EL::StatusCode::FAILURE;
   }

  for(auto e: *mcEvents){
    const xAOD::TruthParticle *tp1(0), *tp2(0), *tp3(0), *tp4(0), *tC1(0), *tN2(0), *tN1w(0), *tN1z(0), *tW(0), *tZ(0);
    TLorentzVector metV;
    pJets.clear();

    for(auto p: e->truthParticleLinks()){
      if(p.isValid() && (*p)->status()==62){
        if((*p)->pdgId()==1000023){
//           Info("execute()", "par %d is %d (%d), m=%.2f", (*p)->barcode(), (*p)->pdgId(), (*p)->status(), (*p)->m());
//           showDecay(*p, "->");

          tN2 = *p;
          /// find ll N1z
          std::vector< const xAOD::TruthParticle* > ds(3,nullptr);
          auto k = (*p);
          while(getDaughters(k, ds)==1 && k->pdgId()==ds[0]->pdgId()) k=ds[0];
          for(auto l: ds){
            if(l->isLepton()&&abs(l->pdgId())<1000000) {tp3 = l;}
            else if(l->pdgId() == 1000022) {tN1z = l;}
            else{
              std::vector< const xAOD::TruthParticle* > ds2(3, nullptr);
              while(getDaughters(l, ds2)==1 && l->pdgId()==ds2[0]->pdgId()) l=ds2[0];
              if(l->pdgId()==23) {
                tZ = l;
                for(auto m: ds2){
                  if(m->pdgId()==-11 || m->pdgId()==-13 || m->pdgId()==-15) tp3 = m;
                  else tp4 = m;
                 }
               }else{
                 for(auto m: ds2){
                   if(m->isLepton()&&abs(m->pdgId())<1000000){
                     tp4 = m;
                    }else{
                     tN1z = m;
                    }
                 }
               }
             }
           }
         }else if(abs((*p)->pdgId())==1000024){
//           Info("execute()", "par %d is %d (%d)", (*p)->barcode(), (*p)->pdgId(), (*p)->status());
//           showDecay(*p, "->");

          tC1 = *p;
          /// find ll N1w
          std::vector< const xAOD::TruthParticle* > ds(3,nullptr);
          auto k = (*p);
          while(getDaughters(k, ds)==1 && k->pdgId()==ds[0]->pdgId()) k=ds[0];
          for(auto l: ds){
            if(l->isLepton()&&abs(l->pdgId())<1000000) {tp1 = l;}
            else if(l->pdgId() == 1000022) {tN1w = l;}
            else{
              std::vector< const xAOD::TruthParticle* > ds2(3, nullptr);
              while(getDaughters(l, ds2)==1 && l->pdgId()==ds2[0]->pdgId()) l=ds2[0];
              if(abs(l->pdgId())==24) {
                tW = l;
                for(auto m: ds2){
                  if(m->isCharged()) tp2 = m;
                  else tp1 = m;
                 }
               }else{
                 for(auto m: ds2){
                   if(m->isLepton()&&abs(m->pdgId())<1000000){
                     tp2 = m;
                    }else{
                     tN1w = m;
                    }
                 }
               }
             }
           }
         }else continue;
       }
     }
//     Info("execute()", "pC1: %d pN2: %d pN1w: %d pN1z: %d", tC1?tC1->barcode():-999, tN2?tN2->barcode():-999, tN1w?tN1w->barcode():-999, tN1z?tN1z->barcode():-999);
    setKinPar(pC1, tC1); 
    setKinPar(pN2, tN2);
    setKinPar(pN1w, tN1w); 
    setKinPar(pN1z, tN1z);

//     Info("execute()", "p1: %d p2: %d p3: %d p4: %d", tp1?tp1->barcode():-999, tp2?tp2->barcode():-999, tp3?tp3->barcode():-999, tp4?tp4->barcode():-999);
    setKinPar(pl1, tp1); 
    setKinPar(pl2, tp2); 
    setKinPar(pl3, tp3); 
    setKinPar(pl4, tp4); 
 
    if(tN1w) metV += tN1w->p4();
    if(tN1z) metV += tN1z->p4();
    if(tp1&&tp1->isNeutral()) metV += tp1->p4();
    else if(tp2&&tp2->isNeutral()) metV += tp2->p4();

    MET = metV.Pt()*iGeV;
    MET_phi = metV.Phi();

    if(tp1) pl1.MET_dPhi = metV.DeltaPhi(tp1->p4()); 
    if(tp2) pl2.MET_dPhi = metV.DeltaPhi(tp2->p4()); 
    if(tp3) pl3.MET_dPhi = metV.DeltaPhi(tp3->p4()); 
    if(tp4) pl4.MET_dPhi = metV.DeltaPhi(tp4->p4());

    if(tW || !tp1 || !tp2) setKinPar(pW, tW);
    else{
      auto x1 = tp1->p4()+tp2->p4();
      setKinPar(pW, x1, 999);
    }

    if(tZ || !tp3 || !tp4) setKinPar(pZ, tZ);
    else{
      auto x1 = tp3->p4()+tp4->p4();
      setKinPar(pZ, x1, 999);
    }

    auto tpl = (tp1&&tp1->isCharged())?tp1:tp2;

    if(tpl){
      if(tpl->pdgId()>0){
        mSS = tp3?(tp3->p4()+tpl->p4()).M()*iGeV:-999;
        mSC = tp4?(tp4->p4()+tpl->p4()).M()*iGeV:-999;
        SS_dR = tp3?tpl->p4().DeltaR(tp3->p4()):-999;
        SS_dPhi = tp3?tpl->p4().DeltaPhi(tp3->p4()):-999;
      }else{
        mSC = tp3?(tp3->p4()+tpl->p4()).M()*iGeV:-999;
        mSS = tp4?(tp4->p4()+tpl->p4()).M()*iGeV:-999;
        SS_dR = tp4?tpl->p4().DeltaR(tp4->p4()):-999;
        SS_dPhi = tp4?tpl->p4().DeltaPhi(tp4->p4()):-999;
      }
    }

    /// save jets
    // Retrieve the truth jets
    const xAOD::JetContainer* truthJets = 0;
    if(! event->retrieve(truthJets, "AntiKt4TruthJets").isSuccess()){
      Error("execute()", "Failed to retrieve event info collection. Exiting.");
     return EL::StatusCode::FAILURE;
     }

    for(const auto& truthJet : *truthJets) {
      if(truthJet->pt() < 20.*GeV) continue; // pT > 20 GeV 
      if(fabs(truthJet->eta()) > 4.5) continue; // |eta| < 4.5

      pJets.emplace_back();
      auto& jet1 = pJets.back();
      jet1.pt = truthJet->pt()*iGeV;
      jet1.eta = truthJet->eta();
      jet1.phi = truthJet->phi();
      jet1.MET_dPhi = metV.DeltaPhi(truthJet->p4());;
      jet1.id = truthJet->auxdata<int>("PartonTruthLabelID") 
    } // end loop over truth jets

    // Sort by Pt
    std::sort(jets->begin(),jets->end(),SortByPt());

    tree->Fill();
  }

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode mcChecker :: postExecute ()
{
  // Here you do everything that needs to be done after the main event
  // processing.  This is typically very rare, particularly in user
  // code.  It is mainly used in implementing the NTupleSvc.
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode mcChecker :: finalize ()
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



EL::StatusCode mcChecker :: histFinalize ()
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

void mcChecker :: setKinPar(mPar& p, const xAOD::TruthParticle* t)
{
  if(t){
    p.pt = t->pt()*iGeV; p.eta = t->eta(); p.phi = t->phi(); p.id = t->pdgId(); p.mll = t->m()*iGeV;
  }else{
    p.pt = -999; p.eta = -999; p.phi = -999; p.id = -999; p.mll = -999;
  }
}

void mcChecker :: setKinPar(mPar& p, TLorentzVector& t, int pdgid)
{
  p.pt = t.Pt()*iGeV; p.eta = t.Eta(); p.phi = t.Phi(); p.id = pdgid; p.mll = t.M()*iGeV;
}

void mcChecker :: setKinPar(kPar& p, const xAOD::TruthParticle* t)
{
  if(t){
    p.pt = t->pt()*iGeV; p.eta = t->eta(); p.phi = t->phi(); p.id = t->pdgId();
  }else{
    p.pt = -999; p.eta = -999; p.phi = -999; p.id = -999;
  }
}

void mcChecker :: setKinPar(kPar& p, TLorentzVector& t, int pdgid)
{
  p.pt = t.Pt()*iGeV; p.eta = t.Eta(); p.phi = t.Phi(); p.id = pdgid;
}

size_t mcChecker :: getDaughters(const xAOD::TruthParticle* t, std::vector< const xAOD::TruthParticle* >& ds)
{
   ds.clear();
   auto vx = t->decayVtx();
   if(vx){
     for(auto l: vx->outgoingParticleLinks()){
       if(l.isValid()) ds.push_back(*l);
      }
    }
   return ds.size();
}
void mcChecker :: showDecay(const xAOD::TruthParticle* t, std::string x)
{
   auto vx = t->decayVtx();
   if(vx){
     for(auto l: vx->outgoingParticleLinks()){
       if(l.isValid()){
         Info("execute()", "%s %d is %d (%d), m=%.2f, pT,eta,phi(%.2f,%.2f,%.2f)", x.c_str(), (*l)->barcode(), (*l)->pdgId(), (*l)->status(), (*l)->m(), (*l)->pt()*iGeV, (*l)->eta(), (*l)->phi());
         if((*l)->status()!=1) showDecay(*l, "--"+x);
       }
     }
   }
}
