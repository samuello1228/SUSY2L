// C++
#include <iostream>
using namespace std;
// ROOT
#include "TChain.h"
#include "TString.h"
#include "TFile.h"
#include "TH2D.h"
// My packages
#include "../multiLepSearch/susyEvts.h"

int MCId;
Long64_t evn;

int NlepBL;
int nSigLep;
int isTruthLep1;
int isTruthLep2;
int nJets20;
int nBJets20;
Float_t Pt_l;
Float_t Pt_subl;
Float_t DeltaEtaLep;
Float_t met;
Float_t mt;
Float_t meff;
Float_t mljj_comb;
Float_t MT2;

//weight
Float_t mcweight;
Float_t puweight;
Float_t totweight;
Float_t lumiScaling;

Float_t wmu_nom;
Float_t wel_nom;
Float_t wjet_nom;
Float_t wpu_nom_bkg;
Float_t wtrig_nom;
Float_t wchflip_nom;
Float_t wttV_nom;

void checkError(float x, float y, TString name, float limit)
{
    float error = fabs((y-x)/x);
    if(error > limit)
    {
        cout<<"Event number: "<<evn<<": "<<name.Data()<<" are different: "<<x<<", "<<y<<", "<<error<<endl;
    }
}

int main()
{
    const TString CHAIN_NAME = "WZ_nom";
    
    //read old tree
    TString path = "/eos/user/d/dkoeck/WHSS/HistFitterTrees/Trees/IncludingLepTruth/Backgrounds_inclTruth_tWZ_tH.36100.root";
    TChain* tree1 = new TChain(CHAIN_NAME);
    tree1->Add(path.Data());
    cout << "Dani ntuple has events: " <<tree1->GetEntries()<<endl;

    tree1->SetBranchAddress("mcweight", &mcweight);
    tree1->SetBranchAddress("puweight", &puweight);
    tree1->SetBranchAddress("totweight", &totweight);
    tree1->SetBranchAddress("lumiScaling", &lumiScaling);

    tree1->SetBranchAddress("wmu_nom", &wmu_nom);
    tree1->SetBranchAddress("wel_nom", &wel_nom);
    tree1->SetBranchAddress("wjet_nom", &wjet_nom);
    tree1->SetBranchAddress("wpu_nom_bkg", &wpu_nom_bkg);
    tree1->SetBranchAddress("wtrig_nom", &wtrig_nom);
    tree1->SetBranchAddress("wchflip_nom", &wchflip_nom);
    tree1->SetBranchAddress("wttV_nom", &wttV_nom);

    tree1->SetBranchAddress("MCId", &MCId);
    tree1->SetBranchAddress("evn", &evn);

    tree1->SetBranchAddress("NlepBL", &NlepBL);
    tree1->SetBranchAddress("nSigLep", &nSigLep);
    tree1->SetBranchAddress("isTruthLep1", &isTruthLep1);
    tree1->SetBranchAddress("isTruthLep2", &isTruthLep2);
    tree1->SetBranchAddress("nJets20", &nJets20);
    tree1->SetBranchAddress("nBJets20", &nBJets20);
    tree1->SetBranchAddress("Pt_l", &Pt_l);
    tree1->SetBranchAddress("Pt_subl", &Pt_subl);
    tree1->SetBranchAddress("DeltaEtaLep", &DeltaEtaLep);
    tree1->SetBranchAddress("met", &met);
    tree1->SetBranchAddress("mt", &mt);
    tree1->SetBranchAddress("meff", &meff);
    tree1->SetBranchAddress("mljj_comb", &mljj_comb);
    tree1->SetBranchAddress("MT2", &MT2);
    //tree1->SetBranchAddress("", &);

    //Our ntuple
    TChain* tree2 = new TChain("evt2l");
    //tree2->Add("/afs/cern.ch/user/c/clo/AnalysisBase-02-04-31-test/SUSY2L/code/run/t1_old/data-myOutput/test.root");
    //tree2->Add("/afs/cern.ch/user/c/clo/AnalysisBase-02-04-31-test/SUSY2L/code/run/t1/data-myOutput/test.root");
    //TString path2 = "/eos/user/c/clo/ntuple/AnalysisBase-02-04-31-c1d5464d/user.clo.v13.0.";
    TString path2 = "/eos/user/c/clo/ntuple/AnalysisBase-02-04-31-new/user.clo.v13.2.";
    tree2->Add(path2 + "VV_CT10_myOutput.root/user.clo.mc15_13TeV.361071.*.myOutput.root*"); int sampleID = 361071; double XS = 0.042287*0.91;
    //tree2->Add(path2 + "VV_221_myOutput.root/user.clo.mc15_13TeV.363491.*.myOutput.root*"); int sampleID = 363491;
    susyEvts* evt = new susyEvts(tree2);
    cout<<"Our ntuple has events: "<<tree2->GetEntries()<<endl;

    std::vector<TString> cutflowList;
    cutflowList.push_back("=2BaseLep && =2SigLep");
    cutflowList.push_back("=2BaseLep && =2SigLep,w");
    cutflowList.push_back("pt1");
    cutflowList.push_back("pt1,w");
    cutflowList.push_back("pt2");
    cutflowList.push_back("pt2,w");
    cutflowList.push_back("isTruthLep1");
    cutflowList.push_back("isTruthLep1,w");
    cutflowList.push_back("isTruthLep2");
    cutflowList.push_back("isTruthLep2,w");
    cutflowList.push_back("nBJets20");
    cutflowList.push_back("nBJets20,w");
    cutflowList.push_back("nJet");
    cutflowList.push_back("nJet,w");
    cutflowList.push_back("DeltaEtaLep");
    cutflowList.push_back("DeltaEtaLep,w");
    cutflowList.push_back("met");
    cutflowList.push_back("met,w");
    cutflowList.push_back("mt");
    cutflowList.push_back("mt,w");
    cutflowList.push_back("meff");
    cutflowList.push_back("meff,w");
    cutflowList.push_back("mljj");
    cutflowList.push_back("mljj,w");
    cutflowList.push_back("MT2");
    cutflowList.push_back("MT2,w");
    //cutflowList.push_back("");
    //cutflowList.push_back(",w");

    TH1D* hCutflow_Dani = new TH1D("cutflow_Dani", "cutflow_Dani", cutflowList.size(), 0, cutflowList.size());
    TH1D* hCutflow_Samuel = new TH1D("cutflow_Samuel", "cutflow_Samuel", cutflowList.size(), 0, cutflowList.size());
    for(unsigned int i=0;i<cutflowList.size();i++)    
    {
        hCutflow_Dani->GetXaxis()->SetBinLabel(i+1,cutflowList[i].Data());
        hCutflow_Samuel->GetXaxis()->SetBinLabel(i+1,cutflowList[i].Data());
    }

    int SR = 1;
    //SR = 2;

    //cutflow for Dani ntuple
    //for (long i = 0; i < 100000; i++)
    for (long i = 0; i < tree1->GetEntries(); i++)
    {
        tree1->GetEntry(i);

        if(MCId != sampleID) continue;
        if(NlepBL !=2) continue;
        if(nSigLep !=2) continue;
        double weight = totweight*lumiScaling;
        hCutflow_Dani->Fill("=2BaseLep && =2SigLep",1);
        hCutflow_Dani->Fill("=2BaseLep && =2SigLep,w",weight);

        if(Pt_l < 25000) continue;
        hCutflow_Dani->Fill("pt1",1);
        hCutflow_Dani->Fill("pt1,w",weight);

        if(Pt_subl < 25000) continue;
        hCutflow_Dani->Fill("pt2",1);
        hCutflow_Dani->Fill("pt2,w",weight);

        if(isTruthLep1 != 1) continue;
        hCutflow_Dani->Fill("isTruthLep1",1);
        hCutflow_Dani->Fill("isTruthLep1,w",weight);

        if(isTruthLep2 != 1) continue;
        hCutflow_Dani->Fill("isTruthLep2",1);
        hCutflow_Dani->Fill("isTruthLep2,w",weight);

        if(nBJets20 !=0) continue;
        hCutflow_Dani->Fill("nBJets20",1);
        hCutflow_Dani->Fill("nBJets20,w",weight);

        if(SR == 1)
        {
            if(nJets20 !=1) continue;
            hCutflow_Dani->Fill("nJet",1);
            hCutflow_Dani->Fill("nJet,w",weight);

            if(fabs(DeltaEtaLep)>1.5) continue;
            hCutflow_Dani->Fill("DeltaEtaLep",1);
            hCutflow_Dani->Fill("DeltaEtaLep,w",weight);

            if(met<100000) continue;
            hCutflow_Dani->Fill("met",1);
            hCutflow_Dani->Fill("met,w",weight);

            if(mt<140000) continue;
            hCutflow_Dani->Fill("mt",1);
            hCutflow_Dani->Fill("mt,w",weight);

            if(meff<260000) continue;
            hCutflow_Dani->Fill("meff",1);
            hCutflow_Dani->Fill("meff,w",weight);

            if(mljj_comb>=180000) continue;
            hCutflow_Dani->Fill("mljj",1);
            hCutflow_Dani->Fill("mljj,w",weight);

            if(MT2<80000) continue;
            hCutflow_Dani->Fill("MT2",1);
            hCutflow_Dani->Fill("MT2,w",weight);
        }
        else if(SR == 2)
        {
            if(nJets20 !=2 && nJets20 !=3) continue;
            hCutflow_Dani->Fill("nJet",1);
            hCutflow_Dani->Fill("nJet,w",weight);

            hCutflow_Dani->Fill("DeltaEtaLep",1);
            hCutflow_Dani->Fill("DeltaEtaLep,w",weight);

            if(met<100000) continue;
            hCutflow_Dani->Fill("met",1);
            hCutflow_Dani->Fill("met,w",weight);

            if(mt<120000) continue;
            hCutflow_Dani->Fill("mt",1);
            hCutflow_Dani->Fill("mt,w",weight);

            if(meff<240000) continue;
            hCutflow_Dani->Fill("meff",1);
            hCutflow_Dani->Fill("meff,w",weight);

            if(mljj_comb>=130000) continue;
            hCutflow_Dani->Fill("mljj",1);
            hCutflow_Dani->Fill("mljj,w",weight);

            if(MT2<70000) continue;
            hCutflow_Dani->Fill("MT2",1);
            hCutflow_Dani->Fill("MT2,w",weight);
        }
    }

    for(unsigned int j=1;j<=cutflowList.size();j++)
    {
        cout<<hCutflow_Dani->GetXaxis()->GetBinLabel(j)<<": ";
        if(j%2 == 1) cout<<int(hCutflow_Dani->GetBinContent(j))<<endl;
        if(j%2 == 0) cout<<std::fixed<<std::setprecision(10)<<hCutflow_Dani->GetBinContent(j)<<endl;
    }

    double nwAOD = 0;
    TObjArray* fileList = tree2->GetListOfFiles();
    for(int i=0;i<fileList->GetEntries();i++)
    {
        cout<<fileList->At(i)->GetTitle()<<endl;
        TFile *f1 = TFile::Open(fileList->At(i)->GetTitle());
        TH1D *h1 = (TH1D*) f1->Get("hCutFlow");
        //for(unsigned int j=1;j<30;j++) icout<<j<<" "<<h1->GetBinContent(j)<<endl;;
        nwAOD += h1->GetBinContent(2);
        
        delete f1;
    }

    double lumi = 36100;
    double commonWeight = XS *lumi /nwAOD;
    //cutflow for My ntuple
    for (long i = 0; i < tree2->GetEntries(); i++)
    {
        tree2->GetEntry(i);

        if(evt->sig.nSigLep < 2) continue;
        if(!evt->sig.isSS) continue;
        if(!evt->sig.JetCut) continue;
        if(evt->sig.isZ) continue;

        if(evt->sig.nBaseLep !=2) continue;
        if(evt->sig.nSigLep !=2) continue;
        double weight = evt->evt.weight * evt->evt.pwt * evt->evt.trigSF * evt->evt.ElSF * evt->evt.MuSF * evt->evt.BtagSF * evt->evt.JvtSF * commonWeight;
        hCutflow_Samuel->Fill("=2BaseLep && =2SigLep",1);
        hCutflow_Samuel->Fill("=2BaseLep && =2SigLep,w",weight);

        if(evt->leps[0].pt < 25) continue;
        hCutflow_Samuel->Fill("pt1",1);
        hCutflow_Samuel->Fill("pt1,w",weight);

        if(evt->leps[1].pt < 25) continue;
        hCutflow_Samuel->Fill("pt2",1);
        hCutflow_Samuel->Fill("pt2,w",weight);

        if(evt->leps[0].lepTruth != 1) continue;
        hCutflow_Samuel->Fill("isTruthLep1",1);
        hCutflow_Samuel->Fill("isTruthLep1,w",weight);

        if(evt->leps[1].lepTruth != 1) continue;
        hCutflow_Samuel->Fill("isTruthLep2",1);
        hCutflow_Samuel->Fill("isTruthLep2,w",weight);

        if(evt->sig.nBJet !=0) continue;
        hCutflow_Samuel->Fill("nBJets20",1);
        hCutflow_Samuel->Fill("nBJets20,w",weight);

        if(SR == 1)
        {
            if(evt->sig.nJet !=1) continue;
            hCutflow_Samuel->Fill("nJet",1);
            hCutflow_Samuel->Fill("nJet,w",weight);

            if(fabs(evt->leps[0].eta - evt->leps[1].eta) > 1.5) continue;
            hCutflow_Samuel->Fill("DeltaEtaLep",1);
            hCutflow_Samuel->Fill("DeltaEtaLep,w",weight);

            if(evt->sig.Met<100) continue;
            hCutflow_Samuel->Fill("met",1);
            hCutflow_Samuel->Fill("met,w",weight);

            if(evt->leps[0].mT<140) continue;
            hCutflow_Samuel->Fill("mt",1);
            hCutflow_Samuel->Fill("mt,w",weight);

            if(evt->sig.HT + evt->sig.Met < 260) continue;
            hCutflow_Samuel->Fill("meff",1);
            hCutflow_Samuel->Fill("meff,w",weight);

            if(evt->sig.mlj>=180) continue;
            hCutflow_Samuel->Fill("mljj",1);
            hCutflow_Samuel->Fill("mljj,w",weight);

            if(evt->sig.mT2<80) continue;
            hCutflow_Samuel->Fill("MT2",1);
            hCutflow_Samuel->Fill("MT2,w",weight);
        }
        else if(SR == 2)
        {
            if(evt->sig.nJet !=2 && evt->sig.nJet !=3) continue;
            hCutflow_Samuel->Fill("nJet",1);
            hCutflow_Samuel->Fill("nJet,w",weight);

            hCutflow_Samuel->Fill("DeltaEtaLep",1);
            hCutflow_Samuel->Fill("DeltaEtaLep,w",weight);

            if(evt->sig.Met<100) continue;
            hCutflow_Samuel->Fill("met",1);
            hCutflow_Samuel->Fill("met,w",weight);

            if(evt->leps[0].mT<120) continue;
            hCutflow_Samuel->Fill("mt",1);
            hCutflow_Samuel->Fill("mt,w",weight);

            if(evt->sig.HT + evt->sig.Met < 240) continue;
            hCutflow_Samuel->Fill("meff",1);
            hCutflow_Samuel->Fill("meff,w",weight);

            if(evt->sig.mlj>=130) continue;
            hCutflow_Samuel->Fill("mljj",1);
            hCutflow_Samuel->Fill("mljj,w",weight);

            if(evt->sig.mT2<70) continue;
            hCutflow_Samuel->Fill("MT2",1);
            hCutflow_Samuel->Fill("MT2,w",weight);
        }
    }

    for(unsigned int j=1;j<=cutflowList.size();j++)
    {
        cout<<hCutflow_Samuel->GetXaxis()->GetBinLabel(j)<<": ";
        if(j%2 == 1) cout<<int(hCutflow_Samuel->GetBinContent(j))<<endl;
        if(j%2 == 0) cout<<std::fixed<<std::setprecision(10)<<hCutflow_Samuel->GetBinContent(j)<<endl;
    }

    int j = 0;
    int count = 0;
    for (long i = 0; i < tree2->GetEntries(); i++)
    {
        tree2->GetEntry(i);

        if(evt->sig.nSigLep < 2) continue;
        if(!evt->sig.isSS) continue;
        if(!evt->sig.JetCut) continue;
        if(evt->sig.isZ) continue;

        tree1->GetEntry(j);

        //check event number
        if(evn != long(evt->evt.event) )
        {
            //search for event number
            bool isFound = false;
            cout<<"Searching event number: "<<evt->evt.event<<endl;
            for (long k = 0; k < tree1->GetEntries(); k++)
            {
                tree1->GetEntry(k);
                if(evn == long(evt->evt.event))
                {
                    isFound = true;
                    j=k;
                    cout<<"Event number: "<<evt->evt.event<<" is found. Dani index is "<<k<<endl;
                    break;
                }
            }

            if(!isFound)
            {
                cout<<"Event number: "<<evt->evt.event<<" cannot be found. End. "<<endl;
                return 0;
            }
        }
        //cout<<"Event number: "<<evn<<", index: "<<j<<endl;

        //wight
        checkError(mcweight, evt->evt.weight, "mc gen weight", 2e-7);
        checkError(puweight, evt->evt.pwt, "pileup weight", 2e-7);
        checkError(puweight, wpu_nom_bkg, "pileup weight", 2e-7);
        checkError(wel_nom * wchflip_nom, evt->evt.ElSF, "electron weight", 4e-7);
        checkError(wmu_nom, evt->evt.MuSF, "muon weight", 2e-7);
        checkError(wjet_nom, evt->evt.BtagSF * evt->evt.JvtSF, "jet weight", 2e-7);
        checkError(wtrig_nom, evt->evt.trigSF, "trigger weight", 2e-7);
        //checkError(wttV_nom, evt->evt.ttVSF, "ttV weight", 2e-7);
        checkError(totweight, evt->evt.weight * evt->evt.pwt * evt->evt.trigSF * evt->evt.ElSF * evt->evt.MuSF * evt->evt.BtagSF * evt->evt.JvtSF, "Total weight", 4e-7);

        //number of leptons
        if(NlepBL != evt->sig.nBaseLep) cout<<"nBaseLep are different."<<endl;
        if(NlepBL != int(evt->leps.size())) cout<<"nBaseLep are different."<<endl;
        if(nSigLep != evt->sig.nSigLep) cout<<"nSigLep are different."<<endl;

        //lepTruth
        if(isTruthLep1 != evt->leps[0].lepTruth) cout<<"lepTruth1 are different."<<endl;
        if(isTruthLep2 != evt->leps[1].lepTruth) cout<<"lepTruth2 are different."<<endl;

        //number of jets
        if(nJets20 != evt->sig.nJet) cout<<"nJet are different."<<endl;
        if(nBJets20 != evt->sig.nBJet) cout<<"nBJet are different."<<endl;

        //variables
        checkError(Pt_l, evt->leps[0].pt * 1000, "pt1", 2e-7);
        checkError(Pt_subl, evt->leps[1].pt * 1000, "pt2", 2e-7);
        checkError(DeltaEtaLep, fabs(evt->leps[0].eta - evt->leps[1].eta), "dEta", 2e-7);
        checkError(met, evt->sig.Met * 1000, "MET", 2e-7);
        checkError(mt, evt->leps[0].mT * 1000, "mT", 4e-7);

        //meff
        float pt_sum = 0;
        for(int k=2;k<nSigLep;k++) pt_sum += evt->leps[k].pt;
        checkError(meff, (evt->sig.HT + evt->sig.Met + pt_sum) * 1000, "meff", 4e-7);

        //mlj
        if(nSigLep == 2 && nJets20 > 0) checkError(mljj_comb, evt->sig.mlj * 1000, "mlj", 2e-7);

        //MT2
        if(MT2>10000) checkError(MT2, evt->sig.mT2 * 1000, "MT2", 3e-6);

        j++;
        count++;
    }
    cout<<"Total number of event in Dani selection: "<<count<<endl;

    delete hCutflow_Dani;
    delete hCutflow_Samuel;
    delete tree1;
}

