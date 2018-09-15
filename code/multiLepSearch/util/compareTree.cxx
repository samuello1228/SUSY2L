#include <iostream>
#include <vector>
#include <algorithm>
using namespace std;

#include "TChain.h"
#include "TFile.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TMath.h"

//#include "MT2_ROOT.h"
#include "../multiLepSearch/susyEvts.h"

//variable for Peter ntuple
Long64_t EventNumber;
Long64_t RunNumber;
Long64_t PileUpHash;
float xAODWeightSum;
int randomRN;
int configYear;
int MCId;
float EventWeight;
float PileUpWeight;
float ObjectSF;
float triggerSF;
float chargeFlipSF;
float weight;
float EtMiss;
float EtMissPhi;
float meff;
float minv;
float mt;
float ht;
float mj;
float mjj;
float mjjdR;
float mjRec;
float truthHt;
float truthMET;
int isZevent;
int nW;
int NlepBL;
int NlepSig;
int NlepBLEta;
int Njet25;
int Njet35;
int Njet40;
int Njet50;
int Nbjet;
std::vector<float>* jetBtag;
std::vector<TLorentzVector>* jetTLV;
std::vector<TLorentzVector>* lepTLV;
std::vector<TLorentzVector>* lepTLV_BL;
std::vector<float>* lepZ0;
std::vector<float>* lepCharges;
std::vector<float>* lepCFmedium;
std::vector<float>* lepCFloose;
std::vector<float>* lepTruth;
int SSChannel;
bool SameSign;
std::vector<std::string>* SysNames;
std::vector<float>* SysWeights;
std::vector<float>* PDFWeights;
std::vector<float>* VarWeights;
std::vector<float>* CFTWeights;
std::vector<float>* TrigWeights;
std::vector<std::string>* triggerInfo;
int PDGId1;
int PDGId2;
int GLDec1;
int GLDec2;
int IntPerX;

void checkError(float x, float y, TString name, float limit)
{
    float error = fabs((y-x)/x);
    if(error > limit)
    {
        cout<<name.Data()<<" are different: "<<x<<", "<<y<<", "<<error<<endl;
    }
}

struct pt_index
{
    float pt;
    int index;
};

int main()
{
    //Peter
    TString path1 = "/afs/cern.ch/user/c/clo/AnalysisBase-02-04-31-Peter1/Output_old/data-output/mc15_13TeV.361071.Sherpa_CT10_lllvjj_EW6.merge.DAOD_SUSY2.e3836_s2726_r7772_r7676_p2949.root";
    const TString treeName1 = "DiLeptonTree_SRall";
    TChain* tree1 = new TChain(treeName1);
    tree1->Add(path1.Data());
    cout<<"tree1 has events: "<<tree1->GetEntries()<<endl;

    //configure input branches
    tree1->SetBranchAddress("EventNumber",  &EventNumber);
    tree1->SetBranchAddress("RunNumber",    &RunNumber);
    tree1->SetBranchAddress("PileUpHash",   &PileUpHash);
    tree1->SetBranchAddress("xAODWeightSum", &xAODWeightSum);
    tree1->SetBranchAddress("RandomRN",     &randomRN);
    tree1->SetBranchAddress("ConfigYear",   &configYear);
    tree1->SetBranchAddress("MCId",         &MCId);
    tree1->SetBranchAddress("MCWeight",     &EventWeight);
    tree1->SetBranchAddress("PileUpWeight", &PileUpWeight);
    tree1->SetBranchAddress("ObjectSF",     &ObjectSF);
    tree1->SetBranchAddress("TriggerSF",    &triggerSF);
    tree1->SetBranchAddress("ChargeFlipSF", &chargeFlipSF);
    tree1->SetBranchAddress("TotalWeight",  &weight);
    tree1->SetBranchAddress("EtMiss",       &EtMiss);
    tree1->SetBranchAddress("EtMissPhi",    &EtMissPhi);
    tree1->SetBranchAddress("Meff",         &meff);
    tree1->SetBranchAddress("Minv",         &minv);
    tree1->SetBranchAddress("Mt",           &mt);
    tree1->SetBranchAddress("Ht",           &ht);
    tree1->SetBranchAddress("Mj",           &mj);
    tree1->SetBranchAddress("Mjj",          &mjj);
    tree1->SetBranchAddress("MjjdR",        &mjjdR);
    tree1->SetBranchAddress("MjRec",        &mjRec);
    tree1->SetBranchAddress("TruthHt",      &truthHt);
    tree1->SetBranchAddress("TruthMET",     &truthMET);
    tree1->SetBranchAddress("Zevent",       &isZevent);
    tree1->SetBranchAddress("nW",           &nW);
    tree1->SetBranchAddress("NlepBL",       &NlepBL);
    tree1->SetBranchAddress("NlepSig",      &NlepSig);
    tree1->SetBranchAddress("NlepBLEta",    &NlepBLEta);
    tree1->SetBranchAddress("Njet25",       &Njet25);
    tree1->SetBranchAddress("Njet35",       &Njet35);
    tree1->SetBranchAddress("Njet40",       &Njet40);
    tree1->SetBranchAddress("Njet50",       &Njet50);
    tree1->SetBranchAddress("Nbjet",        &Nbjet);
    tree1->SetBranchAddress("jetBtag",      &jetBtag);
    tree1->SetBranchAddress("jetTLV",       &jetTLV);
    tree1->SetBranchAddress("lepTLV",       &lepTLV);
    tree1->SetBranchAddress("lepTLV_BL",    &lepTLV_BL);
    tree1->SetBranchAddress("lepZ0",        &lepZ0);
    tree1->SetBranchAddress("lepCharges",   &lepCharges);
    tree1->SetBranchAddress("lepCFmedium",  &lepCFmedium);
    tree1->SetBranchAddress("lepCFloose",   &lepCFloose);
    tree1->SetBranchAddress("lepTruth",     &lepTruth);
    tree1->SetBranchAddress("SSChannel",    &SSChannel);
    tree1->SetBranchAddress("isSameSign",   &SameSign);
    tree1->SetBranchAddress("SysNames",     &SysNames);
    tree1->SetBranchAddress("SysWeights",   &SysWeights);
    tree1->SetBranchAddress("PDFWeights",   &PDFWeights);
    tree1->SetBranchAddress("VarWeights",   &VarWeights);
    tree1->SetBranchAddress("CFTWeights",   &CFTWeights);
    tree1->SetBranchAddress("TrigWeights",  &TrigWeights);
    tree1->SetBranchAddress("triggerInfo",  &triggerInfo);
    tree1->SetBranchAddress("PDGId1",       &PDGId1);
    tree1->SetBranchAddress("PDGId2",       &PDGId2);
    tree1->SetBranchAddress("GLDec1",       &GLDec1);
    tree1->SetBranchAddress("GLDec2",       &GLDec2);
    tree1->SetBranchAddress("IntPerX",      &IntPerX);

    //Our ntuple
    TChain* tree2 = new TChain("evt2l");
    tree2->Add("/afs/cern.ch/user/c/clo/AnalysisBase-02-04-31-test/SUSY2L/code/run/t1_old/data-myOutput/test.root");
    //tree2->Add("/afs/cern.ch/user/c/clo/AnalysisBase-02-04-31-test/SUSY2L/code/run/t1/data-myOutput/test.root");
    susyEvts* evt = new susyEvts(tree2);
    cout<<"tree2 has events: "<<tree2->GetEntries()<<endl;

    //check the total number of events
    if(tree1->GetEntries() != tree2->GetEntries())
    {
        cout<<"Error: Number of events are different."<<endl;
        return 0;
    }

    //loop over events
    for (long i=0;i<tree1->GetEntries();i++)
    {
        tree1->GetEntry(i);
        evt->GetEntry(i);

        //event number
        if(EventNumber != long(evt->evt.event) )
        {
            cout<<"Error: Event number are different."<<endl;
            return 0;
        }
        //cout<<"Event number: "<<EventNumber<<", index: "<<i<<endl;

        //leptons
        //find sorting index
        vector<pt_index> BaseLep;
        for(unsigned int j=0;j<evt->leps.size();j++)
        {
            pt_index temp;
            temp.pt = evt->leps[j].pt;
            temp.index = j;
            BaseLep.push_back(temp);
        }
        sort(BaseLep.begin(), BaseLep.end(), [](pt_index a, pt_index b)->bool{return a.pt > b.pt;});

        //check baseline leptons
        if(NlepBL != int(lepTLV_BL->size())) cout<<"nBaseLep are different."<<endl;
        if(NlepBL != evt->sig.nBaseLep) cout<<"nBaseLep are different."<<endl;
        if(NlepBL != int(evt->leps.size())) cout<<"nBaseLep are different."<<endl;
        for(unsigned int j=0;j<lepTLV_BL->size();j++)
        {
            int index = BaseLep[j].index;
            checkError(lepTLV_BL->at(j).Pt(), evt->leps[index].pt * 1000, "pt", 2e-7);
            checkError(lepTLV_BL->at(j).Eta(), evt->leps[index].eta, "eta", 2e-7);
            checkError(lepTLV_BL->at(j).Phi(), evt->leps[index].phi, "phi", 2e-7);
        }

        //check signal leptons
        if(NlepSig != int(lepTLV->size())) cout<<"nSigLep are different."<<endl;
        if(NlepSig != int(lepCharges->size())) cout<<"nSigLep are different."<<endl;
        if(NlepSig != evt->sig.nSigLep) cout<<"nSigLep are different."<<endl;
        for(unsigned int j=0;j<lepTLV->size();j++)
        {
            checkError(lepTLV->at(j).Pt(), evt->leps[j].pt * 1000, "pt", 2e-7);
            checkError(lepTLV->at(j).Eta(), evt->leps[j].eta, "eta", 2e-7);
            checkError(lepTLV->at(j).Phi(), evt->leps[j].phi, "phi", 2e-7);
            int charge1 = lepCharges->at(j);
            int charge2 = evt->leps[j].ID;
            if(charge2 > 0) charge2 = 1;
            else charge2 = -1;
            if(charge1 != charge2) cout<<"charge are different."<<endl;

            //lepTruth
        }

        //check jets
        if(int(jetTLV->size()) != evt->sig.nJet) cout<<"nJet are different."<<endl;
        if(jetTLV->size() != evt->jets.size()) cout<<"nJet are different."<<endl;
        if(Nbjet != evt->sig.nBJet) cout<<"nBJet are different."<<endl;
        for(unsigned int j=0;j<jetTLV->size();j++)
        {
            checkError(jetTLV->at(j).Pt(), evt->jets[j].pt * 1000, "pt", 2e-7);
            checkError(jetTLV->at(j).Eta(), evt->jets[j].eta, "eta", 2e-7);
            checkError(jetTLV->at(j).Phi(), evt->jets[j].phi, "phi", 2e-7);
            bool isbjet1 = jetBtag->at(j);
            bool isbjet2 = evt->jets[j].jFlag & JT_BJET;
            if(isbjet1 != isbjet2) cout<<"btag are different."<<endl;
        }

        //check weight
        checkError(EventWeight, evt->evt.weight, "MC weight", 2e-7);
        checkError(PileUpWeight, evt->evt.pwt, "pile up weight", 2e-7);
        checkError(triggerSF, evt->evt.trigSF, "trigger weight", 2e-7);
        checkError(ObjectSF * chargeFlipSF, evt->evt.ElSF * evt->evt.MuSF * evt->evt.BtagSF * evt->evt.JvtSF, "object weight", 4e-7);
        checkError(weight, evt->evt.weight * evt->evt.pwt * evt->evt.trigSF * evt->evt.ElSF * evt->evt.MuSF * evt->evt.BtagSF * evt->evt.JvtSF, "Total weight", 5e-7);

        //MET
        checkError(EtMiss, evt->sig.Met * 1000, "MET", 2e-7);
        checkError(EtMissPhi, evt->sig.MetPhi, "MetPhi", 2e-7);

        //HT
        if(NlepSig == 2)
        {
            checkError(ht, evt->sig.HT * 1000, "HT", 3e-7);
            checkError(meff, (evt->sig.HT + evt->sig.Met) * 1000, "meff", 3e-7);
            checkError(minv, evt->l12.m * 1000, "mll", 2e-7);
        }

        //leading mT
        checkError(mt, evt->leps[0].mT * 1000, "mT", 3e-7);
    }

    //delete objects
    delete tree1;
    delete evt;
    delete tree2;
    return 0;
}

