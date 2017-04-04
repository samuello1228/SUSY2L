#ifndef multiLepSearch_dibosonSelection_H
#define multiLepSearch_dibosonSelection_H

#include <EventLoop/Algorithm.h>
#include <SUSYTools/SUSYObjDef_xAOD.h>
#include "multiLepSearch/obj_def.h"
#include <TTree.h>
#include <algorithm>

class ZgammaEvt{
 public:
   ZgammaEvt(){tree2=0;}
//    ZgammaEvt(TTree* tr){getTree(tr); tree2=0; resetIter();}
   ~ZgammaEvt(){tree1=0; if(tree2){delete tree2; tree2=0;}}

   TString version;
   TTree* makeTree(TString treename="zphEvt"){
     tree2 = new TTree(treename, "a angles tree");
     tree2->Branch("evt", &evt, EVT_s.c_str());
     tree2->Branch("leps", &vleps);
     tree2->Branch("phos", &vphos);
     tree2->Branch("l12", &l12, R_PAR_s.c_str());
     tree2->Branch("llp", &llp, R_PAR_s.c_str());
     tree2->Branch("jets", &vjets);
     tree2->Branch("truths", &vtruths);
     tree2->Branch("sig", &sig, SIGNATURE_s.c_str());
     return tree2;
   }
//    TTree* makeWeightOnlyTree(TString treename, ZgammaEvt* parent);
//    TTree* makePtCorrTree(TString treename, ZgammaEvt* parent);
//    TTree* makeKinematicsSysTree(TString treename, ZgammaEvt* parent);
//    Int_t GetEntry(Long64_t entry = 0, Int_t getall = 0){return tree1->GetEntry(entry, getall);}
   int fill(){return tree2->Fill();}
   void writeTree(TString treename="zphEvt"){tree2->Write(treename);}
//    Int_t Next();
//    void resetIter(){m_elist=tree1->GetEntryList(); m_el = 0;}

   TTree* tree1;
   TTree* tree2;

   EVT evt;
   SIGNATURE sig;
   LEPTONS leps;
   LEPTONS* vleps = &leps;
   LEPTONS phos;
   LEPTONS* vphos = &phos;
   R_PAR l12;
   R_PAR llp;
   JETS jets;
   JETS* vjets = &jets;
   TRUTHS truths;
   TRUTHS* vtruths = &truths;

 private:
//    void getTree(TTree* t){};
//    TEntryList* m_elist;
//    Long64_t m_el;
};

class GoodRunsListSelectionTool;
class TH1;

// namespace xAOD{
//   class ElectronContainer;
//   class PhotonContainer;
//   class MuonContainer;
//   class JetContainer; 
// }


class dibosonSelection : public EL::Algorithm
{
  // put your configuration variables here as public variables.
  // that way they can be set directly from CINT and python.
public:
  // float cutValue;



  // variables that don't get filled at submission time should be
  // protected from being send from the submission node to the worker
  // node (done by the //!)
public:
  std::string CF_outputName;
  std::string CF_outputTreeName;
  std::string CF_derivationName;

  int CF_isMC{0};
  std::vector< std::string > CF_grlFiles;
  std::string CF_ConfigFile;
  std::vector< std::string > CF_trigNames;

  /// PRW files
  std::vector< std::string > CF_PRW_confFiles;
  std::vector< std::string > CF_PRW_lcalcFiles;

  // this is a standard constructor
  dibosonSelection ();

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();
  EL::StatusCode fillLepton(xAOD::IParticle* p, L_PAR& l, unsigned int i);
  EL::StatusCode fillLepton(xAOD::Electron* p, L_PAR& l, unsigned int i);
  EL::StatusCode fillLepton(xAOD::Muon* p, L_PAR& l, unsigned int i);
  void fillLeptonCommon(xAOD::IParticle* p, L_PAR& l);
  int addTruthPar(const xAOD::TruthParticle* p, TRUTHS& v, int saveParent=0);

  template <class T>
    void sortByPt(std::vector< T* >& Ls){sort(Ls.begin(), Ls.end(), [](T* a, T* b)->bool{return a->pt()>b->pt();});}

//   void sortByPt(std::vector< xAOD::IParticle* >& Ls){
//     sort(Ls.begin(), Ls.end(), [](xAOD::IParticle* a, xAOD::IParticle* b)->bool{return a->pt()>b->pt();});
//    }

private:
  std::string m_name; //!
  TH1 *m_hCutFlow{0}; //!
  xAOD::ElectronContainer* m_electrons{0};//!
  xAOD::PhotonContainer* m_photons{0}; //!
  xAOD::MuonContainer* m_muons{0}; //!
  xAOD::JetContainer* m_jets{0}; //!
  ST::SUSYObjDef_xAOD* m_objTool{0}; //!
  GoodRunsListSelectionTool *m_grl; //!

  ZgammaEvt* m_zphEvt{0}; //!
  // this is needed to distribute the algorithm to the workers
  ClassDef(dibosonSelection, 1);
};

#endif
